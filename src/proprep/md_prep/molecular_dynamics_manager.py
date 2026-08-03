"""
Molecular Dynamics Manager - Unified MD Setup and Execution

Combines template/workflow setup with simulation execution in a transparent,
user-controlled interface. No black-box operations or hidden presets.
"""

import os
import time
import json
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
)
from proprep.utils.file_browser import (
    remap_recorded_index, annotate_selected_path,
    remap_recorded_index_by_key, annotate_recorded_key,
)
from .user_data_manager import UserDataManager
from .amber_workflow_components import AMBERWorkflowCore, AMBERMonitor
from .workflow_centric_step1 import WorkflowCentricStep1Manager
from .molecular_dynamics_commands import (
    AnalyzeSimulationsCommand,
    SetupSingleSimulationCommand,
    PerformSimulationsCommand,
    ImportTemplateCommand,
    ImportWorkflowCommand,
)
from .layout_helpers import (
    WorkflowLayoutFormatter,
    TemplatePreviewFormatter,
    format_simulation_name,
    format_file_path
)


@dataclass
class SimulationConfig:
    """Configuration for a single simulation."""
    name: str
    template_id: str
    mdin_path: str
    engine: str
    prmtop: Optional[str] = None
    rst7: Optional[str] = None
    mpi_tasks: Optional[int] = None
    gpu_ids: Optional[str] = None
    hardware_config: Optional[Dict] = None
    # Workflow-related fields
    workflow_id: Optional[str] = None  # ID of parent workflow if part of one
    workflow_step: Optional[int] = None  # Step number within workflow
    depends_on: Optional[str] = None  # Previous step dependency
    # Execution status
    status: str = "active"  # "active" or "hold"
    # Restraint configuration from Step 2
    restraints: Optional[Dict] = None
    # Display fields (kept separate from the combined 'name' used as queue key)
    step_name: Optional[str] = None  # e.g., "Initial Energy Minimization"
    structure_label: Optional[str] = None  # e.g., "structure_with_md_names"
    # Simulation type from protocol step (minimization, heating, equilibration, production)
    simulation_type: Optional[str] = None
    # Protocol-centric fields (edits live on protocol steps, not custom template files)
    parameter_overrides: Optional[Dict] = None  # Structured parameter edits
    mdin_content_override: Optional[str] = None  # Full MDIN text if user edited via nano


@dataclass
class WorkflowConfig:
    """Configuration for a multi-step workflow on the same system."""
    workflow_id: str
    name: str
    description: str
    system_prmtop: str  # Single topology for all steps
    initial_rst7: str   # Starting coordinates
    steps: List[SimulationConfig]  # Ordered list of simulations
    hardware_config: Optional[Dict] = None  # Default hardware for all steps
    preset_id: Optional[str] = None  # Reference to workflow preset if used
    # Constant pH / Redox MD settings (applied to production steps only)
    cpin_file: Optional[str] = None  # Path to CPIN file from cpinutil
    cpin_config: Optional[Dict] = None  # Full config dict from workspace
    cpmd_settings: Optional[Dict] = None  # e.g. {'icnstph': 2, 'solvph': 7.0, 'ntcnstph': 100, 'ntrelax': 200}


class SimulationQueue:
    """Manage ordered queue of simulations and workflows."""
    
    def __init__(self, workspace=None):
        self.workspace = workspace
        self._local_queue: List[SimulationConfig] = []
        self._workflows: Dict[str, WorkflowConfig] = {}  # Store workflow configurations
    
    @property
    def queue(self):
        """Get queue from workspace or local storage."""
        # Only reload from workspace if local queue is empty
        if not self._local_queue and self.workspace:
            queue_data = self.workspace.get('md_simulation_queue', [])
            workflow_data = self.workspace.get('md_workflows', {})
            
            # Load workflows first
            for wf_id, wf_dict in workflow_data.items():
                if isinstance(wf_dict, dict):
                    self._workflows[wf_id] = WorkflowConfig(
                        workflow_id=wf_dict.get('workflow_id', wf_id),
                        name=wf_dict.get('name', ''),
                        description=wf_dict.get('description', ''),
                        system_prmtop=wf_dict.get('system_prmtop', ''),
                        initial_rst7=wf_dict.get('initial_rst7', ''),
                        steps=[],  # Will be populated from queue
                        hardware_config=wf_dict.get('hardware_config'),
                        preset_id=wf_dict.get('preset_id'),
                        cpin_file=wf_dict.get('cpin_file'),
                        cpin_config=wf_dict.get('cpin_config'),
                        cpmd_settings=wf_dict.get('cpmd_settings'),
                    )
            
            # Convert dicts back to SimulationConfig objects
            for item in queue_data:
                if isinstance(item, dict):
                    config = SimulationConfig(
                        name=item.get('name', ''),
                        template_id=item.get('template_id', ''),
                        mdin_path=item.get('mdin_path', ''),
                        engine=item.get('engine', ''),
                        prmtop=item.get('prmtop'),
                        rst7=item.get('rst7'),
                        mpi_tasks=item.get('mpi_tasks'),
                        gpu_ids=item.get('gpu_ids'),
                        workflow_id=item.get('workflow_id'),
                        workflow_step=item.get('workflow_step'),
                        depends_on=item.get('depends_on'),
                        status=item.get('status', 'active'),  # Default to active for old configs
                        restraints=item.get('restraints'),
                        step_name=item.get('step_name'),
                        structure_label=item.get('structure_label'),
                        simulation_type=item.get('simulation_type'),
                        parameter_overrides=item.get('parameter_overrides'),
                        mdin_content_override=item.get('mdin_content_override'),
                    )
                    if 'hardware_config' in item:
                        config.hardware_config = item.get('hardware_config', {})
                    self._local_queue.append(config)
                elif isinstance(item, SimulationConfig):
                    self._local_queue.append(item)
        return self._local_queue
    
    def _sync_to_workspace(self):
        """Sync local queue and workflows to workspace."""
        if self.workspace:
            # Sync simulation queue
            queue_data = []
            for config in self._local_queue:
                item = {
                    'name': config.name,
                    'template_id': config.template_id,
                    'mdin_path': config.mdin_path,
                    'engine': config.engine,
                    'prmtop': config.prmtop,
                    'rst7': config.rst7,
                    'mpi_tasks': config.mpi_tasks,
                    'gpu_ids': config.gpu_ids,
                    'workflow_id': config.workflow_id,
                    'workflow_step': config.workflow_step,
                    'depends_on': config.depends_on
                }
                if config.hardware_config:
                    item['hardware_config'] = config.hardware_config
                if config.restraints:
                    item['restraints'] = config.restraints
                if config.step_name:
                    item['step_name'] = config.step_name
                if config.structure_label:
                    item['structure_label'] = config.structure_label
                if config.simulation_type:
                    item['simulation_type'] = config.simulation_type
                if config.parameter_overrides:
                    item['parameter_overrides'] = config.parameter_overrides
                if config.mdin_content_override:
                    item['mdin_content_override'] = config.mdin_content_override
                queue_data.append(item)
            self.workspace.set('md_simulation_queue', queue_data)
            
            # Sync workflows
            workflow_data = {}
            for wf_id, workflow in self._workflows.items():
                wf_entry = {
                    'workflow_id': workflow.workflow_id,
                    'name': workflow.name,
                    'description': workflow.description,
                    'system_prmtop': workflow.system_prmtop,
                    'initial_rst7': workflow.initial_rst7,
                    'hardware_config': workflow.hardware_config,
                    'preset_id': workflow.preset_id,
                }
                if workflow.cpin_file:
                    wf_entry['cpin_file'] = workflow.cpin_file
                if workflow.cpin_config:
                    wf_entry['cpin_config'] = workflow.cpin_config
                if workflow.cpmd_settings:
                    wf_entry['cpmd_settings'] = workflow.cpmd_settings
                workflow_data[wf_id] = wf_entry
            self.workspace.set('md_workflows', workflow_data)
        
    def add_simulation(self, config: SimulationConfig, index: Optional[int] = None):
        """Add simulation at index (or end if None)."""
        queue = self.queue  # Get current queue
        if index is None:
            self._local_queue.append(config)
        else:
            self._local_queue.insert(index, config)
        self._sync_to_workspace()
            
    def remove_simulation(self, index: int) -> bool:
        """Remove simulation at index."""
        queue = self.queue  # Get current queue
        if 0 <= index < len(self._local_queue):
            self._local_queue.pop(index)
            self._sync_to_workspace()
            return True
        return False
        
    def move_simulation(self, from_idx: int, to_idx: int) -> bool:
        """Move simulation from one position to another."""
        queue = self.queue  # Get current queue
        if 0 <= from_idx < len(self._local_queue) and 0 <= to_idx < len(self._local_queue):
            item = self._local_queue.pop(from_idx)
            self._local_queue.insert(to_idx, item)
            self._sync_to_workspace()
            return True
        return False
    
    def clear(self):
        """Clear all simulations."""
        self._local_queue.clear()
        self._sync_to_workspace()
    
    def __len__(self) -> int:
        """Return number of simulations in queue."""
        return len(self.queue)
    
    def __iter__(self):
        """Allow iteration over simulation configs."""
        return iter(self.queue)
    
    def __getitem__(self, index):
        """Allow indexing into the queue."""
        return self.queue[index]
        
    def add_workflow(self, workflow: WorkflowConfig):
        """Add a workflow and its simulations to the queue."""
        self._workflows[workflow.workflow_id] = workflow
        
        # Add all workflow steps to the queue
        for step in workflow.steps:
            self._local_queue.append(step)
        
        self._sync_to_workspace()
    
    def remove_workflow(self, workflow_id: str) -> bool:
        """Remove a workflow and all its simulations from the queue."""
        if workflow_id not in self._workflows:
            return False
            
        # Remove all simulations belonging to this workflow
        self._local_queue = [
            sim for sim in self._local_queue 
            if sim.workflow_id != workflow_id
        ]
        
        # Remove workflow config
        del self._workflows[workflow_id]
        self._sync_to_workspace()
        return True
    
    def get_workflow_simulations(self, workflow_id: str) -> List[SimulationConfig]:
        """Get all simulations belonging to a workflow."""
        return [
            sim for sim in self._local_queue 
            if sim.workflow_id == workflow_id
        ]
    
    def display_queue(self, console: Console):
        """Display current queue with workflow grouping."""
        if not self.queue:
            console.print("[yellow]No simulations queued[/yellow]")
            return

        # Helper to extract structure name and step name from simulation name
        def parse_sim_name(sim_name):
            """Extract structure and step names from format: {structure}_{step_name}"""
            parts = sim_name.split('_')
            if len(parts) >= 2:
                step_keywords = ['Energy', 'System', 'NPT', 'NVT', 'Production', 'Equilibration', 'Heating', 'Minimization']
                for i, part in enumerate(parts):
                    if any(keyword in part for keyword in step_keywords):
                        structure = '_'.join(parts[:i])
                        step = '_'.join(parts[i:])
                        return structure, step
            return sim_name, ""

        console.print(f"\n[bold]Simulation Queue ({len(self.queue)} simulations):[/bold]")

        # Group simulations by workflow
        standalone = []
        workflows_sims = {}

        for config in self.queue:
            if config.workflow_id:
                if config.workflow_id not in workflows_sims:
                    workflows_sims[config.workflow_id] = []
                workflows_sims[config.workflow_id].append(config)
            else:
                standalone.append(config)

        # Display standalone simulations first
        if standalone:
            console.print("\n  [cyan]Standalone Simulations:[/cyan]")
            for i, config in enumerate(standalone, 1):
                status_color = "green" if config.status == "active" else "yellow"
                status_display = f"[{status_color}]{config.status.upper()}[/{status_color}]"
                structure, step = parse_sim_name(config.name)
                console.print(f"    {i}. {step} ({status_display})")

        # Display workflows
        for wf_id, sims in workflows_sims.items():
            # Get structure name from first simulation
            structure_name = "Unknown"
            if sims:
                structure_name, _ = parse_sim_name(sims[0].name)

            if wf_id in self._workflows:
                wf = self._workflows[wf_id]
                console.print(f"\n  [cyan]{structure_name} - {wf.name}[/cyan]")
                console.print(f"    [grey50]{wf.description}[/grey50]")
            else:
                console.print(f"\n  [cyan]{structure_name} - Protocol[/cyan]")

            for step, config in enumerate(sorted(sims, key=lambda x: x.workflow_step or 0), 1):
                # Extract just the step name
                _, step_name = parse_sim_name(config.name)

                status_color = "green" if config.status == "active" else "yellow"
                status_display = f"[{status_color}]{config.status.upper()}[/{status_color}]"
                console.print(f"    Step {step}: {step_name} ({status_display})")
            
    def get_all(self) -> List[Tuple[int, str, str, str]]:
        """Get all simulations as (index, name, template_id, mdin_path)."""
        return [
            (i+1, config.name, config.template_id, config.mdin_path)
            for i, config in enumerate(self.queue)
        ]
        
    def __len__(self):
        return len(self.queue)
        
    def __getitem__(self, index):
        return self.queue[index]


@register_module
class MolecularDynamicsManager(ProcessingModule):
    """Unified MD setup and execution manager."""
    
    NAME = "Molecular Dynamics Manager"
    CATEGORY = "preparation"
    DESCRIPTION = "Configure, execute, monitor, and analyze MD simulations"
    VERSION = "3.0.0"
    REQUIRES = ["Topology Generator"]  # Need topology/coords
    PRIORITY = 11  # Run after TLEaP
    
    def __init__(self):
        super().__init__()
        self.user_data_manager = None
        self.workflow_core = AMBERWorkflowCore()
        self._simulation_queue = None  # Will be initialized when processor is set
        self._workflow_outputs = {}  # Track output files from completed workflow steps
        # Layout formatters for consistent UX
        self.layout = WorkflowLayoutFormatter(self.console)
        self.template_preview = TemplatePreviewFormatter(self.console)
        # Cluster profile + run plan, set during Step 4 when the user loads them
        self._cluster_profile = None
        self._run_plan = None
        # Don't set initial working dir here - will get it from processor/workspace
        
    def set_processor(self, processor):
        """Set the processor reference."""
        self.processor = processor
        # Propagate to the core so its interactive prompts record context
        if hasattr(self, "workflow_core") and self.workflow_core is not None:
            self.workflow_core.processor = processor
        if hasattr(processor, 'console'):
            self.user_data_manager = UserDataManager(console=processor.console)
        else:
            self.user_data_manager = UserDataManager(console=Console())

    @property
    def console(self):
        """Get console from processor if available."""
        if (hasattr(self, "processor") and self.processor and 
            hasattr(self.processor, "console")):
            return self.processor.console
        else:
            return Console()
    
    @property
    def workspace(self):
        """Get workspace from processor if available."""
        if (hasattr(self, "processor") and self.processor and 
            hasattr(self.processor, "_get_workspace")):
            return self.processor._get_workspace()
        return None
    
    @property
    def simulation_queue(self):
        """Get simulation queue, creating with workspace if needed."""
        if self._simulation_queue is None:
            self._simulation_queue = SimulationQueue(workspace=self.workspace)
        # Inject template resolver for better display
        if hasattr(self, 'user_data_manager') and self.user_data_manager:
            self._simulation_queue._template_resolver = self.user_data_manager.load_custom_template
        return self._simulation_queue
    
    @property
    def working_directory(self):
        """Get the proper working directory for file operations."""
        import os
        
        # First, try to use the actual current working directory
        # This should be where the user launched ProPrep from
        cwd = Path(os.getcwd())
        
        # Check if we're in a package installation directory (not user's test directory)
        cwd_str = str(cwd)
        in_package = ('site-packages' in cwd_str or 
                     cwd_str.endswith('/src/proprep') or
                     cwd_str.endswith('/src/proprep/tests') or
                     cwd_str.endswith('/proprep/tests'))
        
        if in_package:
            # We're in package installation directory, try alternatives

            # Try to get project directory from workspace (set at startup)
            if self.workspace and self.workspace.has('project_directory'):
                project_dir = self.workspace.get('project_directory')
                if project_dir and Path(project_dir).exists():
                    return Path(project_dir)

            # Try to get from workspace (legacy)
            if self.workspace and self.workspace.has('working_directory'):
                workspace_dir = self.workspace.get('working_directory')
                if workspace_dir and Path(workspace_dir).exists():
                    return Path(workspace_dir)

            # Try environment variable
            if 'PROPREP_WORKING_DIR' in os.environ:
                env_dir = Path(os.environ['PROPREP_WORKING_DIR'])
                if env_dir.exists():
                    return env_dir

            # Final fallback
            return Path.home()
        else:
            # We're in a normal directory, use it
            return cwd

    def get_menu_options(self) -> Dict[str, str]:
        """Get available menu options."""
        return {
            "setup": "Setup and configure simulations",
            "execute": "Execute simulations",
            "monitor": "Monitor simulation progress",
            "analyze": "Analyze completed simulations",
            "library": "Manage Templates & Protocols",
            "clusters": "Manage cluster profiles (HPC configurations)",
            "plans": "Manage run plans (per-protocol resource assignments)",
        }

    def get_enhanced_menu_options(self, workspace):
        """
        Get menu options with enhanced status information.

        Args:
            workspace: Current workspace

        Returns:
            List of MenuOption objects with status
        """
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus

        options = []

        # Check workspace state
        simulation_queue = workspace.get('md_simulation_queue', [])
        has_setup = len(simulation_queue) > 0

        # Check for running simulations
        has_running = False
        if hasattr(self, 'running_simulations'):
            has_running = len(self.running_simulations) > 0

        # Option 1: Setup - always available (can use existing prmtop/rst7)
        if has_setup:
            status = OptionStatus.COMPLETED
        else:
            status = OptionStatus.AVAILABLE

        options.append(MenuOption(
            key="1",
            description="Setup and configure simulations",
            status=status
        ))

        # Option 2: Execute - requires setup
        if has_setup:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to setup simulations first] ○"

        options.append(MenuOption(
            key="2",
            description="Execute simulations",
            status=status,
            dependency_text=dep_text
        ))

        # Option 3: Monitor - requires running simulations
        if has_running:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[No simulations currently running] ○"

        options.append(MenuOption(
            key="3",
            description="Monitor simulation progress",
            status=status,
            dependency_text=dep_text
        ))

        # Option 4: Analyze - always available
        options.append(MenuOption(
            key="4",
            description="Analyze completed simulations",
            status=OptionStatus.AVAILABLE
        ))

        # Option 5: Manage Templates & Protocols - always available
        options.append(MenuOption(
            key="5",
            description="Manage Templates & Protocols",
            status=OptionStatus.AVAILABLE
        ))

        # Option 6: Manage cluster profiles (HPC configurations) — always available
        options.append(MenuOption(
            key="6",
            description="Manage cluster profiles (HPC configurations)",
            status=OptionStatus.AVAILABLE
        ))

        # Option 7: Manage run plans — available once a profile exists and a
        # protocol is set up. Keep it listed either way so the entry point is
        # discoverable; handle_menu_option gates the actual action.
        options.append(MenuOption(
            key="7",
            description="Manage run plans (per-protocol resource assignments)",
            status=OptionStatus.AVAILABLE
        ))

        return options

    def get_menu_suggestion(self, workspace):
        """
        Get a suggestion for the next recommended action.

        Args:
            workspace: Current workspace

        Returns:
            Suggestion text or None
        """
        simulation_queue = workspace.get('md_simulation_queue', [])
        has_setup = len(simulation_queue) > 0

        has_running = False
        if hasattr(self, 'running_simulations'):
            has_running = len(self.running_simulations) > 0

        if not has_setup:
            return "Setup simulations (option 1). Import templates/protocols via Manage Templates & Protocols (option 5) if needed"
        elif has_running:
            num_running = len(self.running_simulations) if hasattr(self, 'running_simulations') else 0
            return f"✓ {num_running} simulation(s) running. Monitor progress (option 3) or press [m] to return to the main menu"
        else:
            num_queued = len(simulation_queue)
            return f"✓ {num_queued} simulation(s) configured. Execute them (option 2) or press [m] to return to the main menu"

    def get_workspace_requirements(self) -> List[str]:
        """List what this module needs from workspace.

        The module is usable once a topology is loaded together with some
        form of coordinates — an rst7 for simulation setup, an md_structure_pairs
        entry, or a trajectory for post-hoc analysis. can_process() enforces
        this OR relationship; the listed keys are the minimum set for display.
        """
        return ["parm7_file"]

    def availability_note(self, workspace):
        """Menu note when unavailable (○) — distinguishes a missing
        topology from missing coordinates."""
        if self.can_process(workspace):
            return None
        if workspace.get("parm7_file") is None:
            return "Needs a topology — run Topology Generator first"
        return "Needs coordinates (rst7) or a trajectory"

    def can_process(self, workspace) -> bool:
        """Accept any usable combination of topology + coordinates/trajectory."""
        if workspace.get("parm7_file") is None:
            return False
        has_rst7 = workspace.get("rst7_file") is not None
        has_pairs = bool(workspace.get("md_structure_pairs") or [])
        has_traj = bool(workspace.get("trajectory_files") or [])
        return has_rst7 or has_pairs or has_traj

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return [
            "md_simulation_queue",
            "md_workflows",
            "md_template_assignments",
            "restraint_integration_config",
            "topology_extracted_pdb",
            "md_structure_pairs",
            "preferred_amber_engine",
            "mpi_tasks",
            "gpu_ids",
        ]

    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection using command pattern."""
        if option == "setup":
            command = SetupSingleSimulationCommand(self.processor)
            return command.execute()

        elif option == "execute":
            command = PerformSimulationsCommand(self.processor)
            return command.execute()

        elif option == "monitor":
            self._monitor_simulation()
            return True

        elif option == "analyze":
            command = AnalyzeSimulationsCommand(self.processor)
            return command.execute()

        elif option == "library":
            self._manage_library()
            return True

        elif option == "clusters":
            self._manage_clusters()
            return True

        elif option == "plans":
            self._manage_plans()
            return True

        elif option == "view":
            command = WorkspaceOverviewCommand(self.processor)
            return command.execute()

        return False

    def _show_system_resources(self):
        """Display available system resources."""
        try:
            # CPU info
            total_cpus = os.cpu_count()
            try:
                load_avg = os.getloadavg()[0]  # 1-minute average
                idle_cpus = max(1, int(total_cpus - load_avg))
            except (OSError, AttributeError):
                # getloadavg not available on all platforms
                idle_cpus = total_cpus // 2  # Conservative estimate
                
            # GPU info
            gpu_info = self._get_gpu_info()
            
            resource_text = f"System Resources: {total_cpus} CPU cores ({idle_cpus} appear idle)"
            if gpu_info and gpu_info.get('available', 0) > 0:
                gpu_devices = gpu_info.get('devices', [])
                gpu_ids = [gpu['id'] for gpu in gpu_devices]
                resource_text += f", {len(gpu_devices)} GPUs available (IDs: {', '.join(gpu_ids)})"
            else:
                resource_text += ", No GPUs detected"
                
            self.console.print(f"[grey50]{resource_text}[/grey50]")
            
        except Exception as e:
            self.console.print(f"[grey50]System Resources: Unable to detect ({e})[/grey50]")

    def _is_process_running(self, pid: int) -> bool:
        """Check if a process with given PID is still running."""
        try:
            import os
            import signal
            # Send signal 0 to check if process exists
            os.kill(pid, 0)
            return True
        except (OSError, ProcessLookupError):
            return False
    
    def _get_cpu_info(self) -> Dict[str, int]:
        """Get CPU information."""
        try:
            cpu_count = os.cpu_count() or 1
            return {
                'available': cpu_count,
                'total': cpu_count
            }
        except Exception:
            return {'available': 1, 'total': 1}

    def _get_gpu_info(self) -> Dict[str, Any]:
        """Get GPU information if available."""
        gpu_devices = []
        try:
            result = subprocess.run(
                ["nvidia-smi", "--query-gpu=index,name,memory.free", 
                 "--format=csv,noheader,nounits"],
                capture_output=True, text=True, timeout=5
            )
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        parts = [p.strip() for p in line.split(',')]
                        if len(parts) >= 3:
                            gpu_devices.append({
                                'id': parts[0],
                                'name': parts[1],
                                'memory_free': f"{parts[2]} MB"
                            })
        except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
            pass  # nvidia-smi not available or failed
            
        return {
            'available': len(gpu_devices),
            'devices': gpu_devices,
            'raw_devices': gpu_devices  # For backward compatibility
        }

    def _setup_single_simulations(self) -> bool:
        """Manage single simulation queue with 6-step workflow.

        Uses a step-navigation loop so that 'Previous' goes back exactly
        one step instead of restarting from the beginning.

        Step return conventions:
            True / "EXECUTE_NOW" — advance to next step (or execute)
            False                — go back to previous step
            None                 — exit workflow entirely
        """
        from .amber_controller import AmberController

        # Create controller instance for template management
        controller = AmberController(processor=self.processor)

        self.console.print(f"\n[bold cyan]===== Simulation Setup =====[/bold cyan]")

        # Load existing queue from workspace if available
        self._load_queue_from_workspace()

        # State shared across steps
        structure_pairs = None
        current_step = 0
        self.multi_structure_mode = None
        self._active_structure_label = None

        while current_step <= 5:
            if current_step == 0:
                # Step 0: Structure Files Selection
                structure_pairs = self._step0_structure_file_selection()
                if not structure_pairs:
                    return True  # Nothing selected — back to main menu

                # Prompt for multi-structure mode if multiple pairs selected
                if len(structure_pairs) > 1:
                    self.multi_structure_mode = self._prompt_multi_structure_mode(len(structure_pairs))
                else:
                    self.multi_structure_mode = None

                current_step += 1

            elif current_step == 1:
                # Step 1: Workflow-Centric Template Configuration
                result = self._step1_workflow_centric_configuration(structure_pairs)
                if result is None:
                    return True
                elif result is False:
                    current_step -= 1  # Back to Step 0
                else:
                    current_step += 1

            elif self.multi_structure_mode == "separate" and current_step == 2:
                # Per-structure loop: restraints (step 2) and hardware (step 4) per structure,
                # then queue sequencing (step 3) and review (step 5) globally
                structure_labels = self._get_structure_labels()
                all_done = True

                # Per-structure steps: restraints then hardware
                per_struct_steps = [2, 4]

                self.console.print(
                    "\n[grey50]  Restraints and hardware will be configured for each structure individually.[/grey50]"
                )
                self.console.print(
                    "[grey50]  Queue management and review will follow once all structures are configured.[/grey50]\n"
                )

                for struct_idx, label in enumerate(structure_labels):
                    self._active_structure_label = label
                    self.console.print(
                        f"\n[bold cyan]━━━ Configuring: {label} "
                        f"({struct_idx + 1} of {len(structure_labels)}) ━━━[/bold cyan]"
                    )

                    step_pos = 0  # index into per_struct_steps
                    while step_pos < len(per_struct_steps):
                        ps = per_struct_steps[step_pos]
                        result = self._run_wizard_step(ps)
                        if result is None:
                            self._active_structure_label = None
                            return True  # Exit
                        elif result is False:
                            if step_pos == 0:
                                # Back from restraints of first structure → back to step 1
                                if struct_idx == 0:
                                    all_done = False
                                    break
                                else:
                                    continue
                            else:
                                step_pos -= 1
                        else:
                            step_pos += 1

                    if not all_done:
                        break

                self._active_structure_label = None

                if not all_done:
                    current_step -= 1  # Back to Step 1
                else:
                    current_step = 3  # Continue to global queue sequencing

            else:
                # Standard flow (single structure or "all together" mode)
                result = self._run_wizard_step(current_step)
                if result is None:
                    return True
                elif result is False:
                    if self.multi_structure_mode == "separate" and current_step == 3:
                        # Back from queue sequencing → re-enter per-structure loop
                        current_step = 2
                    else:
                        current_step -= 1
                elif result == "EXECUTE_NOW":
                    self._perform_simulations()
                    return True
                else:
                    if self.multi_structure_mode == "separate" and current_step == 3:
                        # After queue sequencing, skip hardware (already done per-structure)
                        current_step = 5
                    else:
                        current_step += 1

        return True

    def _run_wizard_step(self, step_num):
        """Execute a single wizard step by number. Returns True/False/None/'EXECUTE_NOW'."""
        if step_num == 2:
            return self._step2_restraint_integration()
        elif step_num == 3:
            return self._step3_queue_sequencing()
        elif step_num == 4:
            return self._step4_engine_configuration()
        elif step_num == 5:
            result = self._step5_review_and_save()
            if result == "EXECUTE_NOW":
                return "EXECUTE_NOW"
            return result
        return True

    def _prompt_multi_structure_mode(self, num_structures):
        """Prompt user to choose how to configure multiple structures."""
        self.console.print(f"\n[bold]Multiple structures selected ({num_structures}).[/bold]")
        self.console.print("How would you like to configure them?\n")
        self.console.print("   [cyan](a)[/cyan] Apply same settings to all structures")
        self.console.print("   [cyan](s)[/cyan] Configure each structure separately\n")

        choice = prompt_with_context(
            self.processor,
            "Configuration mode",
            choices=["a", "s"],
            default="a",
            module="MD Manager - Multi-Structure",
            description="Choose multi-structure configuration mode",
            options_map={
                "a": "Apply same settings to all structures",
                "s": "Configure each structure separately"
            }
        )

        mode = "all" if choice == "a" else "separate"
        if mode == "all":
            self.console.print("[grey50]Settings will be applied to all structures together.[/grey50]\n")
        else:
            self.console.print("[grey50]You will configure each structure individually.[/grey50]\n")
        return mode

    def _get_structure_labels(self):
        """Return unique structure labels from the simulation queue, preserving order."""
        seen = set()
        labels = []
        for sim in self.simulation_queue.queue:
            label = sim.structure_label or sim.name
            if label not in seen:
                seen.add(label)
                labels.append(label)
        return labels

    def _get_active_queue(self):
        """Return simulations for the active structure, or all if no filter is set."""
        if not getattr(self, '_active_structure_label', None):
            return list(self.simulation_queue.queue)
        return [s for s in self.simulation_queue.queue
                if s.structure_label == self._active_structure_label]

    def _display_queue_status(self):
        """Display compact queue status."""
        if not self.simulation_queue:
            self.console.print("[grey50]Queue empty[/grey50]")
        else:
            count = len(self.simulation_queue)
            self.console.print(f"[grey50]Queue: {count} simulation{'s' if count != 1 else ''}[/grey50]")

    def _show_command_help(self):
        """Show command help."""
        self.console.print(f"\n[bold]Available Commands:[/bold]")
        self.console.print("  [cyan]view[/cyan]                    - Show full queue")
        self.console.print("  [cyan]view N[/cyan]                  - Show details of simulation N")
        self.console.print("  [cyan]view templates[/cyan]          - Show available templates")
        self.console.print("  [cyan]add N[/cyan]                   - Add template N to end of queue")
        self.console.print("  [cyan]add N before/after/at M[/cyan] - Add template N at specific position")
        self.console.print("  [cyan]create[/cyan]                  - Create new template with wizard")
        self.console.print("  [cyan]modify N[/cyan]                - Modify template (metadata/content)")
        self.console.print("  [cyan]remove N[/cyan]                - Remove simulation N")
        self.console.print("  [cyan]move X to Y[/cyan]             - Move simulation from X to Y")
        self.console.print("  [cyan]swap X Y[/cyan]                - Swap simulations X and Y")
        self.console.print("  [cyan]import[/cyan]                  - Browse and import .mdin file")
        self.console.print("  [cyan]import /path/file.mdin[/cyan]  - Import specific file")
        self.console.print("  [cyan]execute[/cyan]                 - Run the queue")
        self.console.print("  [cyan]help[/cyan]                    - Show this help")
        self.console.print("  [cyan]done[/cyan]                    - Exit setup")

    def _view_specific_simulation(self, identifier):
        """View details of a specific simulation."""
        if not self.simulation_queue:
            self.console.print("[yellow]No simulations in queue[/yellow]")
            return
            
        try:
            # Try to parse as number first
            index = int(identifier) - 1
            if 0 <= index < len(self.simulation_queue):
                config = self.simulation_queue.queue[index]
            else:
                self.console.print(f"[red]Position {identifier} out of range (1-{len(self.simulation_queue)})[/red]")
                return
        except ValueError:
            # Try to find by name
            config = None
            for c in self.simulation_queue.queue:
                if c.name.lower() == identifier.lower():
                    config = c
                    break
            if not config:
                self.console.print(f"[red]Simulation '{identifier}' not found[/red]")
                return
                
        # Display detailed information
        self.console.print(f"\n[bold]Simulation Details:[/bold]")
        self.console.print(f"  Name: {config.name}")
        self.console.print(f"  Template: {config.template_id}")
        self.console.print(f"  Engine: {config.engine if config.engine else '[not set]'}")
        self.console.print(f"  MDIN Path: {config.mdin_path}")
        if hasattr(config, 'parameters') and config.parameters:
            self.console.print(f"  Parameters:")
            for key, value in config.parameters.items():
                self.console.print(f"    {key} = {value}")

    def _view_available_templates(self, controller):
        """Show available templates with numbers."""
        templates = self.user_data_manager.list_templates()
        
        if not templates:
            self.console.print("[yellow]No .mdin templates found[/yellow]")
            return
            
        # Use the existing categorized display with proper priority ordering
        template_choices = controller._display_categorized_templates(templates, show_numbers=True)
        
        # Store the template mapping for add commands
        self._template_choices = template_choices

    def _handle_add_command(self, command, controller):
        """Handle add commands: add N [before/after/at M]"""
        parts = command.split()
        if len(parts) < 2:
            self.console.print("[red]Usage: add N [before/after/at M][/red]")
            return
            
        try:
            template_num = int(parts[1])
        except ValueError:
            self.console.print("[red]Template number must be an integer[/red]")
            return
            
        # Get available templates and create the mapping like view templates does
        templates = self.user_data_manager.list_templates()
        if not templates:
            self.console.print("[yellow]No templates available[/yellow]")
            return
            
        # Use same display logic as view templates to get consistent numbering
        template_choices = controller._display_categorized_templates(templates, show_numbers=True)
        
        if template_num < 1 or template_num > len(template_choices):
            self.console.print(f"[red]Template {template_num} not found (1-{len(template_choices)})[/red]")
            return
            
        # Get template by number from the choices
        template_id = template_choices[str(template_num)]
        template_metadata = templates[template_id]
        
        # Parse position
        if len(parts) >= 3:
            if len(parts) == 4 and parts[2] in ["before", "after", "at"]:
                try:
                    pos_ref = int(parts[3])
                except ValueError:
                    self.console.print("[red]Position must be an integer[/red]")
                    return
                    
                if parts[2] == "before":
                    pos = pos_ref - 1
                elif parts[2] == "after":
                    pos = pos_ref
                else:  # at
                    pos = pos_ref - 1
                    
                if pos < 0:
                    pos = 0
                elif pos > len(self.simulation_queue):
                    pos = len(self.simulation_queue)
            else:
                self.console.print("[red]Usage: add N [before/after/at M][/red]")
                return
        else:
            # Just "add N" - append at end
            pos = len(self.simulation_queue)
            
        # Get simulation name
        default_name = f"{template_metadata['name'].replace(' ', '_')}"
        sim_name = prompt_with_context(
            self.processor,
            "Simulation name",
            default=default_name,
            module="MD Manager - Template Management",
            description="Enter simulation name"
        )
        
        # Create and add simulation
        config = SimulationConfig(
            name=sim_name,
            template_id=template_id,
            mdin_path=template_metadata.get("template_path", template_id),
            engine=""
        )
        
        self.simulation_queue.add_simulation(config, pos if pos < len(self.simulation_queue) else None)
        self.console.print(f"[green]✓ Added '{sim_name}' at position {pos+1}[/green]")

    def _handle_create_command(self, controller):
        """Handle create command: create new template with wizard and add to queue"""
        from .amber_annotated_templates import SimulationType
        
        self.console.print("\n[bold cyan]Create New Template with Wizard[/bold cyan]")
        
        # Let the wizard handle simulation type selection entirely
        try:
            # Start with a default type - the wizard will let user change it
            base_template = controller.template_system.load_template(SimulationType.MINIMIZATION)
            configured_template = controller._configure_with_wizard(base_template)
            
            if not configured_template:
                self.console.print("[yellow]Template creation cancelled[/yellow]")
                return
                
            # Display the generated template content
            self.console.print(f"\n[bold cyan]Generated Template Content:[/bold cyan]")
            template_content = configured_template.generate_mdin_content()
            
            # Display with line numbers for easy reference
            lines = template_content.split('\n')
            for i, line in enumerate(lines, 1):
                self.console.print(f"[grey50]{i:3}│[/grey50] {line}")
                
            self.console.print(f"\n[bold]Template Actions:[/bold]")
            self.console.print("  1. Use as-is and add to queue")
            self.console.print("  2. Edit template content")
            self.console.print("  3. Cancel")

            action = prompt_with_context(
                self.processor,
                "Select action",
                choices=["1","2","3"],
                default="1",
                module="MD Manager - Template Creation",
                description="Template action",
                options_map={"1": "Use as-is and add to queue", "2": "Edit template content", "3": "Cancel"}
            )
            
            if action == "3":
                self.console.print("[yellow]Template creation cancelled[/yellow]")
                return
            elif action == "2":
                # Launch editor to modify the content
                edited_content = self._edit_template_content(template_content)
                if edited_content:
                    # Update the template with edited content
                    # This is a bit hacky but works for now
                    configured_template._mdin_content = edited_content
                    template_content = edited_content
                else:
                    self.console.print("[yellow]Edit cancelled, using original content[/yellow]")
                    
            # Get template name from user
            default_name = f"{configured_template.simulation_type.value}_custom"
            template_name = prompt_with_context(
                self.processor,
                "Template name",
                default=default_name,
                module="MD Manager - Template Creation",
                description="Enter template name"
            )
            
            # Save as custom template
            template_id = controller._save_as_custom_template(configured_template)
            
            if template_id:
                # Get the template metadata to populate mdin_path
                template_metadata = self.user_data_manager.list_templates()[template_id]
                
                # Create simulation config and add to queue
                config = SimulationConfig(
                    name=template_name,
                    template_id=template_id,
                    mdin_path=template_metadata.get('template_path', ''),  # Populate immediately
                    engine=""  # Will be configured during execution
                )
                
                # Add to queue
                self.simulation_queue.add_simulation(config)
                self.console.print(f"[green]✓ Created template '{template_name}' and added to queue[/green]")
            else:
                self.console.print("[red]Failed to save custom template[/red]")
                
        except Exception as e:
            self.console.print(f"[red]Error creating template: {e}[/red]")
            
    def _edit_template_content(self, content):
        """Launch editor to modify template content."""
        import tempfile
        import os
        import subprocess
        
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mdin', delete=False) as f:
            f.write(content)
            temp_file = f.name
            
        try:
            # Launch editor
            editor = os.environ.get('EDITOR', 'nano')
            self.console.print(f"\n[cyan]Launching {editor} to edit template...[/cyan]")
            
            result = subprocess.run([editor, temp_file], check=True)
            
            # Read back the edited content
            with open(temp_file, 'r') as f:
                edited_content = f.read()
                
            self.console.print("[green]✓ Template edited successfully[/green]")
            return edited_content
            
        except subprocess.CalledProcessError:
            self.console.print("[red]Error during editing[/red]")
            return None
        except FileNotFoundError:
            # Try vi if nano not found
            try:
                subprocess.run(['vi', temp_file], check=True)
                with open(temp_file, 'r') as f:
                    edited_content = f.read()
                self.console.print("[green]✓ Template edited successfully[/green]")
                return edited_content
            except:
                self.console.print("[red]No suitable editor found. Please set EDITOR environment variable[/red]")
                return None
        finally:
            # Clean up temporary file
            try:
                os.unlink(temp_file)
            except:
                pass

    def _handle_modify_command(self, command, controller):
        """Handle modify commands: modify N - modify the template itself"""
        parts = command.split()
        
        if len(parts) < 2:
            self.console.print("[red]Usage: modify N[/red]")
            return
            
        try:
            index = int(parts[1]) - 1
        except ValueError:
            self.console.print("[red]Simulation number must be an integer[/red]")
            return
            
        if not (0 <= index < len(self.simulation_queue)):
            self.console.print(f"[red]Position {parts[1]} out of range (1-{len(self.simulation_queue)})[/red]")
            return
            
        config = self.simulation_queue.queue[index]
        
        # Get the actual template data
        templates = self.user_data_manager.list_templates()
        if config.template_id not in templates:
            self.console.print(f"[red]Template '{config.template_id}' not found[/red]")
            return
            
        template_metadata = templates[config.template_id]
        
        # Display current template information
        self._display_template_details(template_metadata, config, controller)
        
        # Offer modification options
        self._modify_template_interactive(config, template_metadata, controller)
    
    def _display_template_details(self, template_metadata, config, controller=None):
        """Display template metadata and content."""
        self.console.print(f"\n[bold cyan]Template Details for: {config.name}[/bold cyan]\n")
        
        # Display metadata
        self.console.print("[bold]Template Metadata:[/bold]")
        self.console.print(f"  Name: {template_metadata.get('name', 'Unknown')}")
        self.console.print(f"  Description: {template_metadata.get('description', 'No description')}")
        self.console.print(f"  Type: {template_metadata.get('type', 'Unknown')}")
        self.console.print(f"  Priority: {template_metadata.get('priority', 'Not set')}")
        
        # Get and display template content
        template_path = template_metadata.get('template_path', config.mdin_path)
        if template_path:
            full_template_path = self.user_data_manager.template_base_dir / template_path
            if full_template_path.exists():
                self.console.print("\n[bold]Template Content (.mdin):[/bold]")
                self.console.print("[grey50]" + "─" * 60 + "[/grey50]")
                
                try:
                    with open(full_template_path, 'r') as f:
                        content = f.read()
                    
                    # Format content with aligned comments (use AmberController's formatting method)
                    if controller:
                        formatted_content = controller._format_mdin_content(content)
                    else:
                        # Fallback to unformatted content if no controller provided
                        formatted_content = content
                    
                    # Display full content with syntax highlighting (consistent with ProPrep standards)
                    from rich.syntax import Syntax
                    syntax = Syntax(formatted_content, "fortran", theme="monokai", line_numbers=True)
                    self.console.print(syntax)
                    
                except Exception as e:
                    self.console.print(f"  [red]Error reading template: {e}[/red]")
                    
                self.console.print("[grey50]" + "─" * 60 + "[/grey50]")
            else:
                self.console.print(f"[yellow]Template file not found at: {full_template_path}[/yellow]")
        else:
            self.console.print("[yellow]No template path specified[/yellow]")
            
    def _modify_template_interactive(self, config, template_metadata, controller):
        """Interactive template modification options."""
        self.console.print("\n[bold]Template Modification Options:[/bold]")
        self.console.print("  1. Edit metadata (name, description, priority)")
        self.console.print("  2. Edit template content (launch editor)")
        self.console.print("  3. Change to different template")
        self.console.print("  4. Cancel")

        choice = prompt_with_context(
            self.processor,
            "Select option",
            choices=["1","2","3","4"],
            default="4",
            module="MD Manager - Template Modification",
            description="Template modification action",
            options_map={
                "1": "Edit metadata (name, description, priority)",
                "2": "Edit template content (launch editor)",
                "3": "Change to different template",
                "4": "Cancel"
            }
        )
        
        if choice == "1":
            self._modify_template_metadata(config, template_metadata)
        elif choice == "2":
            self._edit_template_directly(config, template_metadata)
        elif choice == "3":
            self._change_to_different_template(config, controller)
        # Choice 4 = Cancel, do nothing
            
    def _modify_template_metadata(self, config, template_metadata):
        """Modify template metadata fields."""
        self.console.print("\n[bold]Modify Metadata:[/bold]")
        self.console.print("[grey50]Type keyword to modify, or press Enter to skip[/grey50]\n")
        
        # Name
        response = prompt_with_context(
            self.processor,
            "Modify 'name'? (Enter new value or press Enter to skip)",
            default="",
            module="MD Manager - Metadata Edit",
            description="Modify template name"
        )
        if response:
            old_name = template_metadata.get('name', '')
            template_metadata['name'] = response
            # Update config name too
            config.name = response.replace(' ', '_')
            self.console.print(f"  [green]✓ Updated name from '{old_name}' to '{response}'[/green]")

        # Description
        response = prompt_with_context(
            self.processor,
            "Modify 'description'? (Enter new value or press Enter to skip)",
            default="",
            module="MD Manager - Metadata Edit",
            description="Modify template description"
        )
        if response:
            template_metadata['description'] = response
            self.console.print(f"  [green]✓ Updated description[/green]")

        # Priority (with validation)
        response = prompt_with_context(
            self.processor,
            "Modify 'priority'? (Enter number or press Enter to skip)",
            default="",
            module="MD Manager - Metadata Edit",
            description="Modify template priority"
        )
        if response:
            try:
                priority = int(response)
                template_metadata['priority'] = priority
                self.console.print(f"  [green]✓ Updated priority to {priority}[/green]")
            except ValueError:
                self.console.print(f"  [red]Invalid priority value (must be integer)[/red]")
                
        # Save metadata changes - create custom template if modifying builtin
        self._save_template_metadata_changes(config, template_metadata)
        
    def _edit_template_directly(self, config, template_metadata):
        """Launch editor to directly edit template file."""
        template_path = template_metadata.get('template_path', config.mdin_path)
        if not template_path or not Path(template_path).exists():
            self.console.print("[red]Template file not found[/red]")
            return
            
        # Launch editor (using system default or vi/nano)
        import os
        import subprocess
        
        editor = os.environ.get('EDITOR', 'nano')  # Default to nano if no EDITOR set
        self.console.print(f"\n[cyan]Launching {editor} to edit template...[/cyan]")
        
        try:
            subprocess.run([editor, str(template_path)], check=True)
            self.console.print("[green]✓ Template edited successfully[/green]")
        except subprocess.CalledProcessError:
            self.console.print("[red]Error editing template[/red]")
        except FileNotFoundError:
            # Try vi if nano not found
            try:
                subprocess.run(['vi', str(template_path)], check=True)
                self.console.print("[green]✓ Template edited successfully[/green]")
            except:
                self.console.print("[red]No suitable editor found. Please set EDITOR environment variable[/red]")
                
    def _change_to_different_template(self, config, controller):
        """Change to a different existing template."""
        templates = self.user_data_manager.list_templates()
        
        if not templates:
            self.console.print("[yellow]No other templates available[/yellow]")
            return
            
        # Display available templates
        template_choices = controller._display_categorized_templates(templates, show_numbers=True)

        choice = prompt_with_context(
            self.processor,
            "Select new template number",
            default="",
            module="MD Manager - Template Selection",
            description="Select different template"
        )
        
        if choice in template_choices:
            new_template_id = template_choices[choice]
            new_metadata = templates[new_template_id]
            
            config.template_id = new_template_id
            config.mdin_path = new_metadata.get("template_path", new_template_id)
            
            self.console.print(f"[green]✓ Changed to template: {new_metadata.get('name', new_template_id)}[/green]")
        else:
            self.console.print("[yellow]Template change cancelled[/yellow]")
            
    def _save_template_metadata_changes(self, config, template_metadata):
        """Save template metadata changes, creating custom template if needed."""
        template_path = template_metadata.get('template_path', config.mdin_path)
        full_template_path = self.user_data_manager.template_base_dir / template_path
        
        # Check if this is a builtin template
        is_builtin = config.template_id.startswith("builtin_")
        
        if is_builtin:
            # Create a custom copy of the builtin template with updated metadata
            self.console.print("\n[cyan]Creating custom copy of builtin template...[/cyan]")
            
            # Read the original template content
            try:
                with open(full_template_path, 'r') as f:
                    original_content = f.read()
            except Exception as e:
                self.console.print(f"[red]Error reading template: {e}[/red]")
                return
                
            # Update the header comments with new metadata
            updated_content = self._update_template_header(original_content, template_metadata)
            
            # Create custom template using UserDataManager
            custom_name = template_metadata.get('name', config.name)
            custom_description = template_metadata.get('description', 'Modified builtin template')
            sim_type = template_metadata.get('simulation_type', 'minimization')
            
            try:
                new_template_id = self.user_data_manager.create_template(
                    name=custom_name,
                    description=custom_description,
                    simulation_type=sim_type,
                    content=updated_content,
                    author="User"
                )
                
                # Update the config to use the new custom template
                config.template_id = new_template_id
                new_metadata = self.user_data_manager.list_templates()[new_template_id]
                config.mdin_path = new_metadata.get('template_path', '')
                
                self.console.print(f"[green]✓ Created custom template: {custom_name}[/green]")
                self.console.print(f"[green]✓ Simulation now uses custom template ID: {new_template_id}[/green]")
                
            except Exception as e:
                self.console.print(f"[red]Error creating custom template: {e}[/red]")
                
        else:
            # This is a custom template, update both the file headers and metadata
            self.console.print("\n[cyan]Updating custom template...[/cyan]")
            
            # Update the template file with new header metadata
            try:
                with open(full_template_path, 'r') as f:
                    current_content = f.read()
                    
                updated_content = self._update_template_header(current_content, template_metadata)
                
                with open(full_template_path, 'w') as f:
                    f.write(updated_content)
                    
                self.console.print("[green]✓ Template file updated with new metadata[/green]")
                
            except Exception as e:
                self.console.print(f"[red]Error updating template file: {e}[/red]")
                return
                
            # Also update the metadata.json file
            try:
                metadata_file = self.user_data_manager.user_template_dir / "metadata.json"
                all_metadata = {}
                
                if metadata_file.exists():
                    with open(metadata_file, 'r') as f:
                        all_metadata = json.load(f)
                        
                if "custom_templates" not in all_metadata:
                    all_metadata["custom_templates"] = {}
                    
                all_metadata["custom_templates"][config.template_id] = template_metadata
                
                with open(metadata_file, 'w') as f:
                    json.dump(all_metadata, f, indent=2)
                    
                self.console.print("[green]✓ Metadata file updated successfully[/green]")
                
            except Exception as e:
                self.console.print(f"[red]Error updating metadata file: {e}[/red]")
                
    def _update_template_header(self, content, metadata):
        """Update template header comments with new metadata."""
        lines = content.split('\n')
        
        # Find where the header ends (first non-comment line or &cntrl)
        header_end = 0
        for i, line in enumerate(lines):
            if not line.strip().startswith('!') and line.strip() and not line.strip().startswith('#'):
                header_end = i
                break
                
        # Build new header
        new_header = []
        if 'name' in metadata:
            new_header.append(f"! TEMPLATE: {metadata['name']}")
        if 'description' in metadata:
            new_header.append(f"! DESCRIPTION: {metadata['description']}")
        if 'simulation_type' in metadata:
            new_header.append(f"! TYPE: {metadata['simulation_type']}")
        if 'priority' in metadata:
            new_header.append(f"! PRIORITY: {metadata['priority']}")
        if 'author' in metadata:
            new_header.append(f"! AUTHOR: {metadata['author']}")
        
        # Keep any original header lines not covered above
        for line in lines[:header_end]:
            if (line.strip().startswith('! VERSION:') or 
                line.strip().startswith('! SOURCE:') or
                line.strip().startswith('! DATE:')):
                new_header.append(line)
                
        new_header.append("")  # Blank line after header
        
        # Reconstruct content
        return '\n'.join(new_header + lines[header_end:])

    def _load_queue_from_workspace(self):
        """Load existing queue configuration from workspace."""
        workspace_queue_file = Path.cwd() / "simulation_queue.json"
        
        if workspace_queue_file.exists():
            try:
                with open(workspace_queue_file, 'r') as f:
                    queue_data = json.load(f)
                    
                # Reconstruct queue from saved data
                for sim_data in queue_data.get('simulations', []):
                    config = SimulationConfig(
                        name=sim_data['name'],
                        template_id=sim_data['template_id'],
                        mdin_path=sim_data['mdin_path'],
                        engine=sim_data.get('engine', 'pmemd'),  # Default to pmemd instead of empty
                        prmtop=sim_data.get('prmtop'),
                        rst7=sim_data.get('rst7'),
                        mpi_tasks=sim_data.get('mpi_tasks'),
                        gpu_ids=sim_data.get('gpu_ids'),
                        hardware_config=sim_data.get('hardware_config'),
                        workflow_id=sim_data.get('workflow_id'),
                        workflow_step=sim_data.get('workflow_step'),
                        depends_on=sim_data.get('depends_on')
                    )
                    # Restore any additional parameters
                    if 'parameters' in sim_data:
                        config.parameters = sim_data['parameters']
                    self.simulation_queue.add_simulation(config)
                    
                self.console.print(f"[green]Loaded {len(queue_data.get('simulations', []))} simulations from workspace[/green]")
                
            except Exception as e:
                self.console.print(f"[yellow]Could not load previous queue: {e}[/yellow]")

    def _step1_workflow_centric_configuration(self, structure_pairs):
        """Step 1: Workflow-centric template configuration with multi-structure support."""
        try:
            # Initialize workflow-centric Step 1 manager
            step1_manager = WorkflowCentricStep1Manager(
                console=self.console,
                user_data_manager=self.user_data_manager,
                workspace=self.workspace
            )
            
            # Execute the new workflow-centric Step 1.
            # Passing multi_structure_mode lets the protocol step honor the
            # wizard-level apply-to-all / configure-separately choice instead
            # of re-asking the user.
            structure_assignments = step1_manager.execute_step1(
                structure_pairs,
                multi_structure_mode=self.multi_structure_mode,
            )
            
            if not structure_assignments:
                return False  # Go back to Step 0
                
            # Convert assignments to simulation queue entries
            self.simulation_queue.clear()  # Clear any existing simulations
            
            for assignment in structure_assignments:
                # Load the custom workflow
                custom_workflow = self.user_data_manager.load_custom_workflow(assignment.custom_workflow_id)
                if not custom_workflow:
                    self.console.print(f"[red]Error: Could not load protocol {assignment.custom_workflow_id}[/red]")
                    continue
                    
                # Create simulation configs for each step in the workflow
                previous_step_name = None
                for i, step in enumerate(custom_workflow.get('steps', [])):
                    step_name = f"{assignment.structure_name}_{step['name']}"

                    config = SimulationConfig(
                        name=step_name,
                        template_id=step.get('template_ref', step.get('custom_template_id', '')),
                        mdin_path='',  # Will be generated later
                        engine='pmemd',  # Default
                        prmtop=str(assignment.structure_pair.get('prmtop', '')),
                        rst7=str(assignment.structure_pair.get('rst7', '')),
                        workflow_id=assignment.custom_workflow_id,
                        workflow_step=i + 1,
                        depends_on=previous_step_name,  # Depends on previous step in workflow
                        step_name=step['name'],
                        structure_label=assignment.structure_name,
                        simulation_type=step.get('type'),
                        parameter_overrides=step.get('parameter_overrides'),
                        mdin_content_override=step.get('mdin_content_override'),
                    )
                    self.simulation_queue.add_simulation(config)
                    previous_step_name = step_name  # Save for next iteration

            # Offer constant pH MD (if a CPIN was generated in the workspace)
            # and register a WorkflowConfig so the downstream CpHMD machinery
            # engages. The live path otherwise never registers a WorkflowConfig,
            # so _get_workflow_for_step always returned None and CpHMD silently
            # never applied — see _maybe_enable_cphmd_for_assignments.
            self._maybe_enable_cphmd_for_assignments(structure_assignments)

            self.console.print(f"[green]✓ Added {len(self.simulation_queue)} simulations to queue[/green]")
            return True
            
        except KeyboardInterrupt:
            return None  # Exit setup
        except Exception as e:
            self.console.print(f"[red]Error in Step 1: {e}[/red]")
            return False

    def _maybe_enable_cphmd_for_assignments(self, structure_assignments) -> None:
        """Offer constant pH MD (if a CPIN exists in the workspace) and register
        a metadata-only WorkflowConfig for each structure whose topology matches
        the CPIN, so the production-step CpHMD machinery engages
        (``_get_workflow_for_step`` -> CPIN staging / icnstph+solvph injection /
        ``-cpin -cpout -cprestrt`` flags).

        The live setup path builds SimulationConfigs but never registered a
        WorkflowConfig, so every downstream CpHMD guard resolved the workflow to
        ``None`` and the (otherwise complete) machinery no-opped. This is the
        single hook that turns it on.

        No-op when there is no CPIN in the workspace, no production step is
        queued, or the user declines — keeping non-CpHMD runs unchanged.
        """
        # Only meaningful if at least one queued step is a production step.
        if not any(self._is_production_step(c) for c in self.simulation_queue.queue):
            return

        # _check_and_offer_cpmd reads cpin_config/cpin_file from the workspace;
        # the prmtop arg is informational only.
        rep_prmtop = ""
        if structure_assignments:
            rep_prmtop = str(structure_assignments[0].structure_pair.get('prmtop', ''))
        cpmd_info = self._check_and_offer_cpmd(rep_prmtop)
        if not cpmd_info:
            return

        cpin_config = cpmd_info.get('cpin_config') or {}
        # Topologies the CPIN is valid for: the cpinutil input and, for explicit
        # solvent, the radii-corrected *_cpin.prmtop that MD actually runs on.
        cpin_prmtops = {
            os.path.basename(p) for p in (
                cpin_config.get('prmtop_file'),
                cpin_config.get('modified_prmtop'),
            ) if p
        }

        registered = 0
        for assignment in structure_assignments:
            prmtop = str(assignment.structure_pair.get('prmtop', ''))
            rst7 = str(assignment.structure_pair.get('rst7', ''))
            # Don't attach a CPIN generated for one structure to a different
            # protein queued in the same wizard run.
            if cpin_prmtops and os.path.basename(prmtop) not in cpin_prmtops:
                continue
            wid = assignment.custom_workflow_id
            # Register metadata only (steps=[]): the queue already holds the
            # SimulationConfigs, and add_workflow() would re-append them. The
            # persistence model keeps workflow metadata and queue sims separate.
            self.simulation_queue._workflows[wid] = WorkflowConfig(
                workflow_id=wid,
                name=assignment.structure_name,
                description="Constant pH MD",
                system_prmtop=prmtop,
                initial_rst7=rst7,
                steps=[],
                cpin_file=cpmd_info.get('cpin_file'),
                cpin_config=cpin_config,
                cpmd_settings=cpmd_info.get('cpmd_settings'),
            )
            registered += 1

        if registered:
            self.simulation_queue._sync_to_workspace()
            self.console.print(
                f"[grey50]  Constant pH MD will be applied to production steps "
                f"of {registered} structure(s).[/grey50]"
            )
        else:
            self.console.print(
                "[yellow]  Note: the generated CPIN does not match the selected "
                "structure topology, so constant pH MD was not enabled.[/yellow]"
            )

    def _step1_template_selection_legacy(self, controller, structure_pairs):
        """Legacy Step 1: Select from existing templates or create new templates and assign to structure pairs."""
        self.console.print(f"\n[bold cyan]===== Step 1: Template Selection & Assignment =====[/bold cyan]")
        
        # Track template-structure assignments
        template_assignments = []  # List of {template_id, structure_pair, name}
        
        while True:
            # Display structure pairs
            self.console.print(f"\n[bold]Available Structure Pairs:[/bold]")
            for i, pair in enumerate(structure_pairs, 1):
                self.console.print(f"  {i}. {pair['name']}")
                
            # Display current assignments
            if template_assignments:
                self.console.print(f"\n[bold cyan]Template-Structure Assignments:[/bold cyan]")
                for i, assignment in enumerate(template_assignments, 1):
                    templates = self.user_data_manager.list_templates()
                    template_name = templates.get(assignment['template_id'], {}).get('name', assignment['template_id'])
                    self.console.print(f"  [{i}] {template_name} → {assignment['structure_pair']['name']}")
            else:
                self.console.print("\n[grey50]No template assignments yet[/grey50]")
                
            self.console.print(f"\n[bold]Template Selection Options:[/bold]")
            self.console.print("  1. Manage template library")
            self.console.print("  2. Assign template from library to structure pair")
            
            # Dynamic numbering for conditional options
            next_num = 3
            option_map = {}
            
            if template_assignments:
                self.console.print(f"  {next_num}. Remove assignment")
                option_map[str(next_num)] = "remove"
                next_num += 1
                
            self.console.print(f"\n[bold]Navigation:[/bold]")
            if template_assignments:
                self.console.print("  n. Next step (Queue Sequencing)")
            else:
                self.console.print("  [grey50]n. Next step (requires at least 1 assignment)[/grey50]")
            self.console.print("  b. ← Back: Structure Files")
            self.console.print("  x. Exit setup")
            
            # Build valid choices dynamically
            valid_choices = ["1", "2"] + list(option_map.keys())
            if template_assignments:
                valid_choices.append("n")
            valid_choices.extend(["b", "x"])

            # Build options map for context
            context_options_map = {
                "1": "Manage template library",
                "2": "Assign template from library to structure pair"
            }
            if template_assignments:
                if "3" in option_map:
                    context_options_map["3"] = "Remove assignment"
                context_options_map["n"] = "Next step (Queue Sequencing)"
            context_options_map["b"] = "← Back: Structure Files"
            context_options_map["x"] = "Exit setup"

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=valid_choices,
                default="1",
                module="MD Manager - Step 0",
                description="Template selection",
                options_map=context_options_map
            )
            
            if choice == "1":
                # Manage template library
                self._manage_template_library_extended(controller)
                    
            elif choice == "2":
                # Assign template from library to structure pair
                result = self._select_existing_template(controller)
                if result:
                    template_ids = result if isinstance(result, list) else [result]
                    for template_id in template_ids:
                        pair = self._select_structure_for_template(structure_pairs, template_id)
                        if pair:
                            template_assignments.append({
                                'template_id': template_id,
                                'structure_pair': pair,
                                'name': f"{template_id}_{pair['name']}"
                            })
                            templates = self.user_data_manager.list_templates()
                            template_name = templates.get(template_id, {}).get('name', template_id)
                            self.console.print(f"[green]✓ Template '{template_name}' assigned to {pair['name']}[/green]")
                    
            elif choice in option_map:
                # Handle dynamic options
                action = option_map[choice]
                if action == "remove":
                    self._remove_template_assignment(template_assignments)
                
            elif choice == "n" and template_assignments:
                # Save template assignments to workspace
                if self.workspace:
                    workspace_assignments = []
                    for assignment in template_assignments:
                        workspace_assignments.append({
                            'template_id': assignment['template_id'],
                            'structure_pair_name': assignment['structure_pair']['name'],
                            'name': assignment['name']
                        })
                    self.workspace.set('md_template_assignments', workspace_assignments)
                
                # Next step - clear existing queue and add all assignments with structure files
                self.simulation_queue.clear()  # Clear any existing simulations
                for assignment in template_assignments:
                    templates = self.user_data_manager.list_templates()
                    template_metadata = templates[assignment['template_id']]
                    config = SimulationConfig(
                        name=f"{template_metadata.get('name', assignment['template_id'])}_{assignment['structure_pair']['name']}",
                        template_id=assignment['template_id'],
                        mdin_path=template_metadata.get('template_path', ''),
                        engine="",
                        prmtop=str(assignment['structure_pair']['prmtop']),
                        rst7=str(assignment['structure_pair']['rst7'])
                    )

                    self.simulation_queue.add_simulation(config)
                        
                self.console.print(f"[green]✓ Added {len(template_assignments)} simulations to queue[/green]")
                return True
                
            elif choice == "b":
                # Previous step - go back to structure selection
                return False
                    
            elif choice == "x":
                # Exit setup
                return None

        return True

    def _step2_restraint_integration(self):
        """
        Step 3 (internal step 2): Restraint Integration.

        Shows a table of current restraint state from protocol templates, with
        options to edit positional restraints, generate masks, configure DISANG,
        or set up GROUP restraints.

        Returns:
            True: Continue to next step
            False: Go back to previous step
            None: Exit setup
        """
        if not self.simulation_queue or len(self.simulation_queue) == 0:
            self.console.print(self.layout.step_header(step_num=2))
            self.console.print("\n[yellow]No simulations in queue. Go back and select a protocol first.[/yellow]")
            choice = prompt_with_context(
                self.processor, "\nNavigation", choices=["b", "x"], default="b",
                module="MD Manager - Step 3 Restraints", description="Empty queue navigation",
                options_map={"b": "← Back: Protocol Selection", "x": "Exit setup"}
            )
            return False if choice == "b" else None

        while True:
            state_list = self._build_restraint_state()

            self.console.print(self.layout.step_header(step_num=2))
            if self._active_structure_label:
                self.console.print(f"  [bold yellow]Structure: {self._active_structure_label}[/bold yellow]")
            self.console.print("[grey50]  Positional restraints (restraintmask) hold atoms near reference coordinates.[/grey50]")
            self.console.print("[grey50]  DISANG restraints enforce distance, angle, or torsion targets.[/grey50]")
            self.console.print()

            self._display_restraint_table(state_list)

            self.console.print("[bold]Options:[/bold]")
            self.console.print("   [cyan](e)[/cyan] Edit positional restraints")
            self.console.print("   [cyan](g)[/cyan] Generate new restraint mask [grey50](atom type, chain, redox site selection)[/grey50]")
            self.console.print("   [cyan](d)[/cyan] Configure DISANG restraints [grey50](distance/angle/torsion)[/grey50]")
            self.console.print("   [cyan](r)[/cyan] Configure GROUP restraints [grey50](advanced: multi-group with FIND criteria)[/grey50]")
            self.console.print("   [cyan](s)[/cyan] Continue with current restraints")
            self.console.print("\n[bold]Navigation:[/bold]")
            self.console.print("   [yellow](b)[/yellow] ← Back: Protocol Selection")
            self.console.print("   [grey50](x)[/grey50] Exit setup")
            self.console.print()

            choice = prompt_with_context(
                self.processor, "Choose option",
                choices=["e", "g", "d", "r", "s", "b", "x"], default="s",
                module="MD Manager - Step 3 Restraints",
                description="Restraint integration options",
                options_map={
                    "e": "Edit positional restraints",
                    "g": "Generate new restraint mask",
                    "d": "Configure DISANG restraints",
                    "r": "Configure GROUP restraints",
                    "s": "Continue with current restraints",
                    "b": "← Back: Protocol Selection",
                    "x": "Exit setup"
                }
            )

            if choice == "e":
                self._edit_positional_restraints(state_list)
            elif choice == "g":
                self._generate_and_apply_mask(state_list)
            elif choice == "d":
                self._configure_disang_integration(state_list)
            elif choice == "r":
                self._configure_group_integration(state_list)
            elif choice == "s":
                return True
            elif choice == "b":
                return False
            elif choice == "x":
                return None

    def _open_restraint_manager(self):
        """
        Open MD Restraint Manager module and return to MD Manager after completion.

        Returns:
            True if user completed restraint management
            False if user cancelled or exited
        """
        try:
            from proprep.structure_prep.md_restraint_manager import MDRestraintModule

            self.console.print("\n[cyan]═══ Opening MD Restraint Manager ═══[/cyan]")
            self.console.print("[grey50]Configure restraints and return here when done[/grey50]\n")

            # Initialize the module
            restraint_module = MDRestraintModule()
            restraint_module.processor = self.processor

            # Initialize to set up the restraint manager
            restraint_module.initialize()

            # Always extract PDB from topology — ensures we use the prmtop/rst7
            # selected in Step 1 with correct AMBER-consecutive numbering
            extracted_pdb = self._extract_pdb_from_topology()
            if extracted_pdb:
                self.workspace.set('topology_extracted_pdb', extracted_pdb)
            else:
                self.console.print("[red]✗ Could not extract PDB from topology.[/red]")
                self.console.print("[grey50]Ensure prmtop and rst7 files are available.[/grey50]")
                return False

            # Process the workspace (opens interactive menu)
            result = restraint_module.process(self.workspace)

            # Return to MD Manager
            self.console.print("\n[cyan]═══ Returning to MD Manager ═══[/cyan]\n")

            # Return True if processing completed successfully
            return result is not None

        except ImportError as e:
            self.console.print(f"[red]Error: Could not import MD Restraint Manager: {e}[/red]")
            return False
        except Exception as e:
            self.console.print(f"[red]Error opening MD Restraint Manager: {e}[/red]")
            import traceback
            traceback.print_exc()
            return False

    def _build_restraint_state(self):
        """
        Build list of effective restraint parameters for each simulation in the queue.

        Resolution order per step (highest priority first):
        1. sim_config.restraints dict (set by Step 3 edits)
        2. sim_config.parameter_overrides dict (set by Step 2 protocol editing)
        3. Template defaults (parsed from resolved MDIN content)

        When _active_structure_label is set, only includes sims for that structure.
        """
        state_list = []
        active_queue = self._get_active_queue()

        for i, sim_config in enumerate(active_queue, 1):
            name = sim_config.step_name or sim_config.name or sim_config.template_id
            restraints = sim_config.restraints or {}
            overrides = sim_config.parameter_overrides or {}

            # Determine simulation type from template metadata
            sim_type = ""
            if hasattr(self, 'user_data_manager') and self.user_data_manager:
                try:
                    tdata = self.user_data_manager.load_custom_template(sim_config.template_id)
                    if tdata:
                        sim_type = tdata.get('simulation_type', tdata.get('type', ''))
                except Exception:
                    pass

            # --- Positional restraints ---
            mask = None
            weight = None
            group_active = 'group' in restraints

            if 'restraintmask' in restraints:
                # Step 3 override (highest priority)
                rc = restraints['restraintmask']
                mask = rc.get('mask') if isinstance(rc, dict) else rc
                weight = rc.get('weight', 10.0) if isinstance(rc, dict) else 10.0
            elif 'restraintmask' in overrides:
                # Step 2 parameter override
                mask = str(overrides['restraintmask']).strip("'\"")
                weight = overrides.get('restraint_wt')
                if weight is not None:
                    weight = float(weight)
            else:
                # Fall back to parsing template MDIN content
                try:
                    mdin_text = self._resolve_mdin_content(sim_config)
                    if mdin_text:
                        params = self._parse_mdin_params(mdin_text)
                        ntr_val = params.get('ntr')
                        if ntr_val and int(ntr_val) == 1:
                            raw_mask = params.get('restraintmask')
                            if raw_mask:
                                mask = str(raw_mask).strip("'\"")
                            wt = params.get('restraint_wt')
                            if wt is not None:
                                weight = float(wt)
                except Exception:
                    pass

            # --- DISANG ---
            disang_file = None
            if 'disang' in restraints:
                disang_file = restraints['disang'].get('file')

            state_list.append({
                'index': i,
                'name': name,
                'type': sim_type,
                'sim_config': sim_config,
                'mask': mask,
                'weight': weight,
                'disang_file': disang_file,
                'group': group_active,
            })

        return state_list

    def _display_restraint_table(self, state_list):
        """Render a Rich table of current restraint state per protocol step.

        When the state list spans multiple structures (i.e. the wizard is in
        apply-to-all mode with no active structure filter), each structure is
        rendered as its own sub-table with a header, mirroring the grouping
        style of the queue step. Row indices stay global so selection commands
        like '2-4' still reference the same entries.
        """
        from rich.table import Table
        from rich.markup import escape

        # Group entries by structure_label while preserving first-seen order.
        groups: "dict[str, list]" = {}
        for entry in state_list:
            label = getattr(entry['sim_config'], 'structure_label', None) or ""
            groups.setdefault(label, []).append(entry)

        multi_group = len(groups) > 1 and any(label for label in groups)

        def render_rows(table, entries):
            for entry in entries:
                if entry['group']:
                    mask_str = "[yellow](GROUP)[/yellow]"
                elif entry['mask']:
                    # Show the full mask — these strings carry the semantics the
                    # user needs to verify (which residues/atoms are held), so
                    # truncating defeats the purpose of the table. Rich wraps the
                    # column when the terminal is narrow rather than clipping.
                    mask_str = escape(entry['mask'])
                else:
                    mask_str = "[grey50](none)[/grey50]"

                if entry['group']:
                    group_data = entry['sim_config'].restraints.get('group', [])
                    if isinstance(group_data, list) and group_data:
                        fcs = [str(g.get('force_constant', '?')) for g in group_data]
                        weight_str = ",".join(fcs)
                    else:
                        weight_str = "[grey50]—[/grey50]"
                elif entry['weight'] is not None:
                    weight_str = f"{entry['weight']:.1f}"
                else:
                    weight_str = "[grey50]—[/grey50]"

                disang_str = Path(entry['disang_file']).name if entry['disang_file'] else "[grey50]—[/grey50]"

                table.add_row(
                    str(entry['index']),
                    entry['name'],
                    mask_str,
                    weight_str,
                    disang_str,
                )

        def new_table():
            t = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
            t.add_column("#", style="grey50", width=3, justify="right")
            t.add_column("Step", min_width=16)
            t.add_column("Mask", min_width=10)
            t.add_column("Weight", min_width=8, justify="right")
            t.add_column("DISANG", min_width=8)
            return t

        self.console.print("[bold]Protocol Restraints:[/bold]")

        if multi_group:
            for label, entries in groups.items():
                display_label = label or "(unlabeled)"
                self.console.print(f"\n  [cyan]{display_label}[/cyan]")
                table = new_table()
                render_rows(table, entries)
                self.console.print(table)
        else:
            self.console.print()
            table = new_table()
            render_rows(table, state_list)
            self.console.print(table)

        # GROUP footnote
        has_group = any(e['group'] for e in state_list)
        if has_group:
            self.console.print("\n[grey50]  * GROUP replaces restraintmask on those steps[/grey50]")

        self.console.print()

    def _edit_positional_restraints(self, state_list, prefill_mask=None):
        """
        Batch-edit positional restraints: select steps, set mask, per-step weights.

        Args:
            state_list: Current restraint state from _build_restraint_state()
            prefill_mask: Optional mask string to pre-fill the mask prompt
        """
        from rich.markup import escape

        while True:
            self.console.print("\n[bold cyan]Edit Positional Restraints[/bold cyan]\n")
            self._display_restraint_table(state_list)

            # Step selection
            selection = prompt_with_context(
                self.processor,
                "Select steps to modify (e.g., 2-4, all, b to go back)",
                default="all",
                module="MD Manager - Step 3 Restraints",
                description="Select steps for positional restraint editing"
            )

            if selection.strip().lower() == "b":
                return

            try:
                selected_indices = self._parse_template_selection(selection, len(state_list))
            except ValueError as e:
                self.console.print(f"[red]Error: {e}[/red]")
                continue

            # Check GROUP conflicts
            group_conflicts = [state_list[idx - 1] for idx in selected_indices if state_list[idx - 1]['group']]
            if group_conflicts:
                names = ", ".join(e['name'] for e in group_conflicts)
                self.console.print(f"\n[yellow]Warning: {names} have GROUP restraints active.[/yellow]")
                self.console.print("[yellow]Positional restraints cannot coexist with GROUP.[/yellow]")
                remove_group = confirm_with_context(
                    self.processor, "Remove GROUP from these steps?", default=False,
                    module="MD Manager - Step 3 Restraints",
                    description="Resolve GROUP/restraintmask conflict"
                )
                if not remove_group:
                    # Remove conflicting indices from selection
                    group_idx = {e['index'] for e in group_conflicts}
                    selected_indices = [idx for idx in selected_indices if idx not in group_idx]
                    if not selected_indices:
                        self.console.print("[yellow]No steps remaining to edit.[/yellow]")
                        continue

            # Determine default mask
            if prefill_mask:
                default_mask = prefill_mask
                prefill_mask = None  # Only use once
            else:
                selected_masks = [state_list[idx - 1]['mask'] for idx in selected_indices if state_list[idx - 1]['mask']]
                if selected_masks and len(set(selected_masks)) == 1:
                    default_mask = selected_masks[0]
                elif selected_masks:
                    default_mask = selected_masks[0]
                else:
                    default_mask = "!@H="

            self.console.print(f"\n[grey50]Common mask presets: !@H= (heavy atoms), @CA,C,O,N (backbone), @CA (alpha carbons)[/grey50]")
            self.console.print(f"[grey50]Enter 'none' to clear restraints for selected steps[/grey50]")
            mask_input = prompt_with_context(
                self.processor,
                f"Restraint mask for selected steps",
                default=default_mask,
                module="MD Manager - Step 3 Restraints",
                description="Enter restraint mask for selected steps"
            )
            user_mask = mask_input.strip()

            # Handle "none" as clearing restraints
            if user_mask.lower() == "none":
                user_mask = None

            # Per-step weight prompts
            step_edits = []
            if user_mask is None:
                # Clearing restraints — no weight needed
                for idx in selected_indices:
                    step_edits.append((idx, None, None))
            else:
                self.console.print("\n[bold]Set weight per step[/bold] [grey50](Enter to keep current)[/grey50]")
                for idx in selected_indices:
                    entry = state_list[idx - 1]
                    current_wt = entry['weight'] if entry['weight'] is not None else 10.0
                    wt_input = prompt_with_context(
                        self.processor,
                        f"  Step {idx} — {entry['name']} [{current_wt}]",
                        default=str(current_wt),
                        module="MD Manager - Step 3 Restraints",
                        description=f"Restraint weight for {entry['name']}"
                    )
                    try:
                        user_weight = float(wt_input)
                    except ValueError:
                        user_weight = current_wt
                    step_edits.append((idx, user_mask, user_weight))

            # Preview: build temporary updated state for display
            preview_state = []
            edits_by_idx = {idx: (m, w) for idx, m, w in step_edits}
            for entry in state_list:
                if entry['index'] in edits_by_idx:
                    m, w = edits_by_idx[entry['index']]
                    preview_entry = dict(entry)
                    preview_entry['mask'] = m
                    preview_entry['weight'] = w
                    preview_entry['group'] = False  # GROUP will be cleared
                    preview_state.append(preview_entry)
                else:
                    preview_state.append(entry)

            self.console.print("\n[bold]Updated restraints:[/bold]\n")
            self._display_restraint_table(preview_state)

            self.console.print("   (e) Edit more steps")
            self.console.print("   (s) Accept and continue")
            self.console.print("   (b) Discard changes\n")

            action = prompt_with_context(
                self.processor, "Action",
                choices=["e", "s", "b"], default="s",
                module="MD Manager - Step 3 Restraints",
                description="Accept, edit more, or discard restraint changes",
                options_map={
                    "e": "Edit more steps",
                    "s": "Accept and continue",
                    "b": "Discard changes"
                }
            )

            if action in ("s", "e"):
                # Apply edits to sim_config.restraints
                for idx, mask, weight in step_edits:
                    sc = state_list[idx - 1]['sim_config']
                    if sc.restraints is None:
                        sc.restraints = {}
                    if mask is None:
                        # Clear restraints for this step
                        sc.restraints.pop('restraintmask', None)
                    else:
                        sc.restraints['restraintmask'] = {'mask': mask, 'weight': weight}
                    # Clear GROUP conflict
                    if 'group' in sc.restraints and group_conflicts:
                        sc.restraints.pop('group', None)
                # Rebuild state_list to reflect applied edits
                state_list[:] = self._build_restraint_state()
                if action == "s":
                    return
                # action == "e": loop continues with updated state
            elif action == "b":
                return

    def _generate_and_apply_mask(self, state_list):
        """
        Generate a restraint mask via the Restraint Manager's structure-aware
        mask builder, then feed the result into the positional restraint editor.
        """
        try:
            from proprep.structure_prep.md_restraint_manager import MDRestraintModule

            # Initialize restraint module
            restraint_module = MDRestraintModule()
            restraint_module.processor = self.processor
            restraint_module.initialize()

            # Always extract PDB from topology for correct AMBER numbering
            extracted_pdb = self._extract_pdb_from_topology()
            if extracted_pdb:
                self.workspace.set('topology_extracted_pdb', extracted_pdb)
            else:
                self.console.print("[red]Could not extract PDB from topology.[/red]")
                return

            self.console.print("\n[cyan]═══ Opening Restraint Mask Generator ═══[/cyan]\n")

            result = restraint_module.restraint_manager.generate_restraint_mask(
                workspace=self.workspace, interactive=True
            )

            self.console.print("\n[cyan]═══ Returning to MD Manager ═══[/cyan]\n")

            if result is None:
                self.console.print("[yellow]Mask generation cancelled.[/yellow]")
                return

            # Read generated mask from workspace
            generated_mask = self.workspace.get('redox_restraint_mask')
            if generated_mask:
                self.console.print(f"[green]Generated mask:[/green] {generated_mask}")
                # Feed into the batch edit flow with the mask pre-filled
                updated_state = self._build_restraint_state()
                self._edit_positional_restraints(updated_state, prefill_mask=generated_mask)
            else:
                self.console.print("[yellow]No mask was generated.[/yellow]")

        except ImportError as e:
            self.console.print(f"[red]Error: Could not import MD Restraint Manager: {e}[/red]")
        except Exception as e:
            self.console.print(f"[red]Error generating restraint mask: {e}[/red]")
            import traceback
            traceback.print_exc()

    def _configure_disang_integration(self, state_list):
        """
        Configure DISANG restraints: open Restraint Manager, then assign
        the exported DISANG file to selected protocol steps.
        """
        # Open the full Restraint Manager (user can configure/import/export DISANG)
        if not self._open_restraint_manager():
            return

        # Check for DISANG file after user returns
        disang_file = self.workspace.get('disang_file') if self.workspace else None
        if not disang_file or not Path(disang_file).exists():
            self.console.print("[yellow]No DISANG file found in workspace. Configure DISANG restraints in the Restraint Manager first.[/yellow]")
            return

        # Show summary
        md_restraints = self.workspace.get('md_restraints', []) if self.workspace else []
        type_counts = {'distance': 0, 'angle': 0, 'torsion': 0}
        for r in md_restraints:
            rtype = (r.get('restraint_type', '') if isinstance(r, dict) else getattr(r, 'restraint_type', '')).lower()
            if rtype in type_counts:
                type_counts[rtype] += 1
        type_str = ", ".join(f"{c} {t}" for t, c in type_counts.items() if c > 0)

        self.console.print(f"\n[bold]DISANG file:[/bold] {Path(disang_file).name}")
        if type_str:
            self.console.print(f"[bold]Restraints:[/bold] {len(md_restraints)} ({type_str})")
        self.console.print()

        # Step selection
        self.console.print("[bold]Select steps for DISANG integration:[/bold]")
        for entry in state_list:
            self.console.print(f"  [{entry['index']}] {entry['name']}")
        self.console.print()

        selection = prompt_with_context(
            self.processor,
            "Apply DISANG to which steps? (e.g., 2-5, all, none)",
            default="all",
            module="MD Manager - Step 3 Restraints",
            description="Select steps for DISANG integration"
        )

        if selection.strip().lower() == "none":
            return

        try:
            selected_indices = self._parse_template_selection(selection, len(state_list))
        except ValueError as e:
            self.console.print(f"[red]Error: {e}[/red]")
            return

        # Dump frequency
        dump_freq_input = prompt_with_context(
            self.processor,
            "Restraint dump frequency (steps)",
            default="500",
            module="MD Manager - Step 3 Restraints",
            description="DISANG dump frequency"
        )
        try:
            dump_freq = int(dump_freq_input)
        except ValueError:
            dump_freq = 500

        # Apply to selected steps
        disang_stem = Path(disang_file).stem
        for idx in selected_indices:
            sc = state_list[idx - 1]['sim_config']
            if sc.restraints is None:
                sc.restraints = {}
            sc.restraints['disang'] = {
                'file': disang_file,
                'dumpave_file': f"{disang_stem}_dump.txt",
                'listout_file': f"{disang_stem}_violations.txt",
                'dump_freq': dump_freq,
            }

        self.console.print(f"\n[green]✓ DISANG applied to {len(selected_indices)} step(s)[/green]")

    def _configure_group_integration(self, state_list):
        """
        Configure GROUP restraints: open GROUP config wizard, then assign
        to selected protocol steps (with conflict warnings for restraintmask).
        """
        try:
            from proprep.structure_prep.md_restraint_manager import MDRestraintModule

            restraint_module = MDRestraintModule()
            restraint_module.processor = self.processor
            restraint_module.initialize()

            self.console.print("\n[cyan]═══ Opening GROUP Configuration ═══[/cyan]\n")
            result = restraint_module._configure_group_restraints()
            self.console.print("\n[cyan]═══ Returning to MD Manager ═══[/cyan]\n")

            if not result:
                return

        except ImportError as e:
            self.console.print(f"[red]Error: Could not import MD Restraint Manager: {e}[/red]")
            return
        except Exception as e:
            self.console.print(f"[red]Error configuring GROUP restraints: {e}[/red]")
            return

        # Read GROUP data from workspace
        group_data = self.workspace.get('group_restraints') if self.workspace else None
        if not group_data:
            self.console.print("[yellow]No GROUP specification configured.[/yellow]")
            return

        # Show summary
        self.console.print(f"\n[bold]GROUP specification:[/bold] {len(group_data)} group(s)")
        for g in group_data:
            self.console.print(f"  • {g.get('title', 'Unnamed')} — {g.get('force_constant', '?')} kcal/mol/Å²")
        self.console.print()

        # Step selection
        selection = prompt_with_context(
            self.processor,
            "Apply GROUP to which steps? (e.g., 2-5, all, none)",
            default="all",
            module="MD Manager - Step 3 Restraints",
            description="Select steps for GROUP integration"
        )

        if selection.strip().lower() == "none":
            return

        try:
            selected_indices = self._parse_template_selection(selection, len(state_list))
        except ValueError as e:
            self.console.print(f"[red]Error: {e}[/red]")
            return

        # Check for restraintmask conflicts
        mask_conflicts = [state_list[idx - 1] for idx in selected_indices if state_list[idx - 1]['mask']]
        if mask_conflicts:
            names = ", ".join(e['name'] for e in mask_conflicts)
            self.console.print(f"\n[yellow]Warning: {names} have positional restraints (restraintmask).[/yellow]")
            self.console.print("[yellow]GROUP replaces restraintmask — these masks will be removed.[/yellow]")
            proceed = confirm_with_context(
                self.processor, "Continue?", default=True,
                module="MD Manager - Step 3 Restraints",
                description="Confirm GROUP replaces restraintmask"
            )
            if not proceed:
                return

        # Apply GROUP to selected steps
        for idx in selected_indices:
            sc = state_list[idx - 1]['sim_config']
            if sc.restraints is None:
                sc.restraints = {}
            sc.restraints['group'] = group_data
            # Clear conflicting restraintmask
            sc.restraints.pop('restraintmask', None)

        self.console.print(f"\n[green]✓ GROUP applied to {len(selected_indices)} step(s)[/green]")


    def _parse_template_selection(self, selection, max_index):
        """
        Parse user's template selection string.

        Args:
            selection: Selection string (e.g., "1,2,5" or "1-3" or "all")
            max_index: Maximum valid index

        Returns:
            list: List of selected indices (1-based)

        Raises:
            ValueError: If selection is invalid
        """
        selection = selection.strip().lower()

        if selection == "all":
            return list(range(1, max_index + 1))

        if selection == "none":
            return []

        # Parse comma-separated and ranges
        indices = []
        parts = selection.split(',')

        for part in parts:
            part = part.strip()

            if '-' in part:
                # Range like "1-3"
                try:
                    start, end = part.split('-')
                    start_idx = int(start.strip())
                    end_idx = int(end.strip())

                    if start_idx < 1 or end_idx > max_index:
                        raise ValueError(f"Range {part} out of bounds (1-{max_index})")

                    indices.extend(range(start_idx, end_idx + 1))
                except ValueError as e:
                    raise ValueError(f"Invalid range '{part}': {e}")
            else:
                # Single number
                try:
                    idx = int(part)
                    if idx < 1 or idx > max_index:
                        raise ValueError(f"Index {idx} out of bounds (1-{max_index})")
                    indices.append(idx)
                except ValueError:
                    raise ValueError(f"Invalid index '{part}'")

        # Remove duplicates and sort
        return sorted(set(indices))

    def _get_template_display_name(self, template_id):
        """Get display name for a template."""
        if hasattr(self, 'user_data_manager') and self.user_data_manager:
            try:
                template_data = self.user_data_manager.load_custom_template(template_id)
                if template_data and 'name' in template_data:
                    return template_data['name']
            except:
                pass
        return template_id

    def _get_structure_display_name(self, full_name):
        """Get display name for a structure."""
        if '_' in full_name:
            parts = full_name.split('_')
            return parts[0] if parts else full_name
        return full_name

    def _apply_restraints_to_queue(self, config):
        """
        Apply restraint configuration to simulation queue.

        Args:
            config: Restraint configuration dict with 'positional' and 'disang' keys
        """
        if not self.simulation_queue or len(self.simulation_queue) == 0:
            return

        # Iterate through simulation queue and apply restraints
        for sim_config in self.simulation_queue.queue:
            template_id = sim_config.template_id

            # Apply positional restraints if configured for this template
            if template_id in config['positional']:
                restraint_config = config['positional'][template_id]
                # Store in simulation config for later template rendering
                if sim_config.restraints is None:
                    sim_config.restraints = {}
                sim_config.restraints['restraintmask'] = restraint_config

            # Apply DISANG restraints if configured for this template
            if template_id in config['disang']:
                disang_config = config['disang'][template_id]
                # Store complete DISANG configuration
                if sim_config.restraints is None:
                    sim_config.restraints = {}
                sim_config.restraints['disang'] = {
                    'file': disang_config['file'],
                    'dumpave_file': disang_config['dumpave_file'],
                    'listout_file': disang_config['listout_file'],
                    'dump_freq': disang_config['dump_freq']
                }

    def _step3_queue_sequencing(self):
        """Step 3: Allow user to define the queuing sequence."""
        # Context: simulations execute sequentially, each step's output coordinates
        # feed into the next step as input.

        # Handle single simulation case
        if len(self.simulation_queue) <= 1:
            sim = self.simulation_queue.queue[0] if len(self.simulation_queue) > 0 else None

            if sim:
                status_lines = [
                    f"Simulation Queue: 1 simulation",
                    f"  1. {format_simulation_name(sim.name, max_length=60)}"
                ]

                if sim.workflow_id:
                    status_lines.append(f"     Part of protocol (Step {sim.workflow_step})")

                status_lines.append("")
                status_lines.append("✓ Only one simulation - proceeding to hardware configuration")
            else:
                status_lines = ["Simulation Queue: Empty"]

            actions = [
                ("n", "Continue to Hardware Configuration"),
                ("b", "← Back: Protocol Selection"),
                ("x", "Exit setup")
            ]

            self.layout.simple_prompt(
                step_num=3,
                status_lines=status_lines,
                actions=actions
            )

            choice = prompt_with_context(
                self.processor,
                "Enter choice",
                choices=["n", "b", "x"],
                default="n",
                module="MD Manager - Step 3 Queue",
                description="Queue sequencing (single simulation)",
                options_map={
                    "n": "Continue to Hardware Configuration",
                    "b": "← Back: Protocol Selection",
                    "x": "Exit setup"
                }
            )

            if choice == "x":
                return None
            elif choice == "b":
                return False
            else:
                return True

        # Multiple simulations - show queue management
        while True:
            # Helper function to extract structure name from simulation name
            def extract_structure_name(sim_name):
                parts = sim_name.split('_')
                if len(parts) >= 2:
                    step_keywords = ['Energy', 'System', 'NPT', 'NVT', 'Production', 'Equilibration', 'Heating', 'Minimization']
                    for i, part in enumerate(parts):
                        if any(keyword in part for keyword in step_keywords):
                            return '_'.join(parts[:i])
                    return '_'.join(parts[:-1]) if len(parts) > 1 else parts[0]
                return sim_name

            # Group simulations by structure
            structure_groups = {}
            for i, sim in enumerate(self.simulation_queue.queue, 1):
                structure_name = extract_structure_name(sim.name)
                if structure_name not in structure_groups:
                    structure_groups[structure_name] = []
                structure_groups[structure_name].append((i, sim))

            # Build table data with structure grouping
            headers = ["#", "Structure", "Simulation", "Step", "Status"]
            rows = []

            for struct_idx, (structure_name, sims) in enumerate(structure_groups.items()):
                for sim_idx, (global_idx, sim) in enumerate(sims):
                    # Extract step name (remove structure prefix)
                    sim_display_name = sim.name
                    if sim.name.startswith(structure_name + '_'):
                        sim_display_name = sim.name[len(structure_name) + 1:]

                    # Get step number
                    step_info = str(sim.workflow_step) if sim.workflow_step else "—"

                    # Status display
                    if sim.status == "active":
                        status_display = "[green]ACTIVE[/green]"
                    else:
                        status_display = "[yellow]HOLD[/yellow]"

                    rows.append([
                        str(global_idx),
                        structure_name if sim_idx == 0 else "",  # Only show structure name on first row
                        sim_display_name,
                        step_info,
                        status_display
                    ])

                # Add separator row between structures (except after last structure)
                if struct_idx < len(structure_groups) - 1:
                    rows.append(["---", "---", "---", "---", "---"])

            # Count total structures
            num_structures = len(structure_groups)
            title_suffix = f"{num_structures} structure(s), {len(self.simulation_queue.queue)} simulation(s)"

            actions = [
                ("k", "Keep current order and continue"),
                ("r", "Reorder simulations"),
                ("t", "Toggle status (ACTIVE/HOLD) for a simulation"),
                ("b", "← Back: Protocol Selection"),
                ("x", "Exit setup")
            ]

            self.layout.config_table(
                step_num=3,
                headers=headers,
                rows=rows,
                actions=actions,
                title=f"Simulation Queue - {title_suffix}"
            )

            choice = prompt_with_context(
                self.processor,
                "Enter choice",
                choices=["k", "r", "t", "b", "x"],
                default="k",
                module="MD Manager - Step 3 Queue",
                description="Queue sequencing (multiple simulations)",
                options_map={
                    "k": "Keep current order (continue to next step)",
                    "r": "Reorder simulations",
                    "t": "Toggle status (ACTIVE/HOLD)",
                    "b": "← Back: Protocol Selection",
                    "x": "Exit setup"
                }
            )

            if choice == "k":
                return True
            elif choice == "r":
                self._reorder_simulations()
            elif choice == "t":
                self._toggle_simulation_status()
            elif choice == "b":
                return False
            elif choice == "x":
                return None

    def _step4_engine_configuration(self):
        """Step 4: Set engine and hardware for each queued simulation."""

        # Get system info
        cpu_info = self._get_cpu_info()
        gpu_info = self._get_gpu_info()

        while True:
            # Build status lines
            status_lines = self._build_step3_status_lines(cpu_info, gpu_info)

            # Build actions
            actions = self._build_step3_actions()

            # Display step header and status
            self.console.print(self.layout.step_header(step_num=4))
            if self._active_structure_label:
                self.console.print(f"  [bold yellow]Structure: {self._active_structure_label}[/bold yellow]")
            self.console.print("[grey50]  Minimization on CPU: GPU fixed-point precision (SPFP) can silently clip[/grey50]")
            self.console.print("[grey50]  large forces from atomic clashes in unrelaxed starting structures.[/grey50]")
            self.console.print("[grey50]  NPT equilibration on CPU: the GPU PME grid is fixed at run start and[/grey50]")
            self.console.print("[grey50]  cannot reorganize as box dimensions change during density equilibration.[/grey50]")
            self.console.print("[grey50]  NVT/production on GPU: once the structure is relaxed and density has[/grey50]")
            self.console.print("[grey50]  converged, GPU acceleration is safe and dramatically faster.[/grey50]")
            self.console.print()

            for line in status_lines:
                self.console.print(f"  {line}")

            if status_lines:
                self.console.print()

            # Display options and navigation with section headers
            self._display_step3_menu(actions)

            # Build valid choices from actions
            valid_choices = [k for k, _ in actions]

            # Build options map from actions
            options_map = {k: d for k, d in actions}

            choice = prompt_with_context(
                self.processor,
                "Enter choice",
                choices=valid_choices,
                default="1",
                module="MD Manager - Step 4 Hardware",
                description="Hardware configuration",
                options_map=options_map
            )

            if choice == "1":
                self._apply_recommended_engines(gpu_info)
            elif choice == "2":
                if len(self.simulation_queue.queue) == 1:
                    self._configure_individual_engines(gpu_info['available'] > 0)
                else:
                    self._configure_individual_engines(gpu_info['available'] > 0)
            elif choice == "3":
                if len(self.simulation_queue.queue) == 1:
                    # Single sim: option 3 is unassigned; fall through to re-prompt.
                    self.console.print("[yellow]Option not available in single-simulation mode[/yellow]")
                else:
                    self._configure_bulk_engines(gpu_info['available'] > 0)
            elif choice == "p":
                self._action_load_cluster_profile()
            elif choice == "l":
                if self._cluster_profile is None:
                    self.console.print("[yellow]Load a cluster profile first (option 'p')[/yellow]")
                else:
                    self._action_load_run_plan()
            elif choice == "s":
                if self._cluster_profile is None:
                    self.console.print("[yellow]Load a cluster profile first (option 'p')[/yellow]")
                elif not any(bool(sim.engine) for sim in self.simulation_queue.queue):
                    self.console.print("[yellow]Assign resources first (option 1 or 2) before saving a run plan[/yellow]")
                else:
                    self._action_save_run_plan()
            elif choice == "n":
                # Check if all engines are set
                unset_count = sum(1 for config in self.simulation_queue.queue if not config.engine)
                if unset_count > 0:
                    self.console.print(f"[yellow]⚠ Warning: {unset_count} simulations have no engine set[/yellow]")
                    continue_str = prompt_with_context(
                        self.processor,
                        "Continue anyway?",
                        choices=["y", "n"],
                        default="n",
                        module="MD Manager - Step 3: Queue Management",
                        description="Continue with unset engines",
                        options_map={"y": "Yes, continue anyway", "n": "No, go back"}
                    )
                    if not (continue_str.lower() == "y"):
                        continue
                return True
            elif choice == "h":
                self._display_engine_help()
            elif choice == "b":
                return False
            elif choice == "x":
                return None

        return True

    def _build_step3_status_lines(self, cpu_info, gpu_info):
        """Build status lines for Step 3 display."""
        status_lines = []
        active_queue = self._get_active_queue()

        # Hardware availability
        gpu_str = f"{gpu_info['available']} GPUs" if gpu_info['available'] > 0 else "0 GPUs"
        status_lines.append(f"Available: {cpu_info['available']} CPUs, {gpu_str}")
        status_lines.append("")

        # Current configuration (compact)
        if len(active_queue) == 1:
            sim = active_queue[0]
            sim_name = format_simulation_name(sim.name, max_length=50)
            status_lines.append(f"Simulation: {sim_name}")
            engine_display = sim.engine if sim.engine else "[not set]"
            status_lines.append(f"Current engine: {engine_display}")
        else:
            status_lines.append(f"Simulations: {len(active_queue)} configured")
            # Show summary
            engine_counts = {}
            for sim in active_queue:
                engine = sim.engine if sim.engine else "[not set]"
                engine_counts[engine] = engine_counts.get(engine, 0) + 1
            for engine, count in engine_counts.items():
                status_lines.append(f"  {count}x {engine}")

        return status_lines

    def _display_step3_menu(self, actions):
        """Display Step 4 menu with section headers matching other steps."""
        # Separate options from navigation
        navigation_keys = {'n', 'h', 'b', 'x'}
        options = [(k, d) for k, d in actions if k not in navigation_keys]
        navigation = [(k, d) for k, d in actions if k in navigation_keys]

        # Display options section
        if options:
            self.console.print("What would you like to do?", style="bold")
            for key, description in options:
                self.console.print(f"   {key}. {description}")
            self.console.print()

        # Display navigation section with cyan header
        self.console.print("[bold cyan]Navigation:[/bold cyan]")
        for key, description in navigation:
            # Determine style based on key
            if key == 'x':
                style = "dim"
            elif key == 'n':
                style = "green"
            elif key == 'b':
                style = "yellow"
            elif key == 'h':
                style = "cyan"
            else:
                style = "white"

            key_str = f"[{style}]{key}[/{style}]"
            self.console.print(f"   {key_str}. {description}")
        self.console.print()

    def _build_step3_actions(self):
        """Build actions menu for Step 4."""
        actions = []

        actions.append(("1", "Apply recommended assignments (CPU for min/NPT, GPU for rest)"))

        if len(self.simulation_queue.queue) == 1:
            # Single simulation - simplified workflow
            actions.append(("2", "Change engine for this simulation"))
            actions.append(("n", "Keep current and continue"))
        else:
            # Multiple simulations
            actions.append(("2", "Configure individual simulation engines"))
            actions.append(("3", "Set all simulations to same engine"))
            actions.append(("n", "Next step (Review & Save)"))

        # Cluster profile + run plan actions. Session-replay-safe: always
        # visible in the same positions; unavailable items render dim.
        profile_loaded = self._cluster_profile is not None
        any_engine_assigned = any(
            bool(sim.engine) for sim in self.simulation_queue.queue
        )
        if profile_loaded:
            p_label = f"Switch / unload cluster profile (active: {self._cluster_profile.name})"
        else:
            p_label = "Load cluster profile"
        actions.append(("p", p_label))
        if profile_loaded:
            actions.append(("l", "Load run plan"))
        else:
            actions.append(("l", "[grey50]Load run plan[/grey50] [grey50 italic](requires cluster profile)[/grey50 italic]"))
        if profile_loaded and any_engine_assigned:
            actions.append(("s", "Save current assignments as run plan"))
        elif not profile_loaded:
            actions.append(("s", "[grey50]Save current assignments as run plan[/grey50] [grey50 italic](requires cluster profile)[/grey50 italic]"))
        else:
            actions.append(("s", "[grey50]Save current assignments as run plan[/grey50] [grey50 italic](no assignments yet)[/grey50 italic]"))

        actions.append(("h", "Why CPU vs. GPU? (detailed explanation)"))
        actions.append(("b", "← Back: Queue Management"))
        actions.append(("x", "Exit setup"))

        return actions

    def _display_engine_help(self):
        """Display detailed explanation of CPU vs GPU engine selection."""
        from rich.panel import Panel

        help_text = """\
[bold cyan]1. Minimizations: Run on CPU for Numerical Precision[/bold cyan]

The GPU code (pmemd.cuda) uses a hybrid precision scheme called SPFP
(Single Precision Fixed Point), designed to balance speed and accuracy
for normal MD simulations. In SPFP, individual force calculations are
done in single precision, but accumulated using a fixed-point integer
representation.

One limitation of SPFP is that forces can be truncated if they overflow
the fixed-precision representation. This should never be a problem during
MD simulations for any well-behaved system. However, for minimization or
very early in the heating phase it can present a problem — especially if
two atoms are close to each other and thus have large VDW repulsions.

The fixed-precision number used on the GPU has a finite dynamic range,
about 100x that ever experienced in MD at 300 K for a reasonable system.
That sounds like a lot, but a bad starting structure can easily exceed
even that. The insidious part is that the GPU minimization may not crash
or throw an error — it just produces a subtly wrong minimized structure.
When you start heating from it, the simulation can violently blow up
because the GPU-minimized geometry has subtle unresolved strain that
the CPU would have properly minimized away.

[bold cyan]2. NPT Equilibration: Run on CPU (PME Grid Limitation)[/bold cyan]

Particle Mesh Ewald (PME) splits long-range electrostatics into a
real-space part (computed directly) and a reciprocal-space part computed
on a 3D grid. The grid spacing and dimensions are set at the start of
the run based on the initial box size. In NVT or NVE simulations, the
box doesn't change, so this is fine. But in NPT, the barostat rescales
the box volume every step as the system equilibrates toward the correct
density.

The GPU code does not automatically reorganize grid cells. If the
periodic box dimensions change too much from their initial values —
which happens readily when the system density is far from equilibrium —
the GPU code will halt with an error: "Periodic box dimensions have
changed too much from their initial values."

The CPU code handles this gracefully because it can reorganize the
neighbor list and PME grid on the fly as the box changes. The GPU code,
by contrast, sets up the grid in GPU memory at the start of the run and
keeps it fixed for performance — reallocating and reorganizing that data
structure on the GPU mid-run is expensive and architecturally complex.

The practical consequence is that early NPT equilibration — when density
is changing rapidly because a fresh solvated system starts at the wrong
density — is precisely the worst case for the GPU code. Running that
phase on CPU sidesteps the problem entirely, and you can switch to the
GPU for production once the box dimensions have stabilized.

[bold cyan]Summary[/bold cyan]

[bold]Minimization[/bold]  → CPU: SPFP can silently clip large forces from clashing atoms
[bold]NPT equil.[/bold]    → CPU: GPU PME grid is fixed, can't handle volume changes
[bold]NVT/production[/bold] → GPU: structure is relaxed, density stable — safe and fast

Both issues are architectural trade-offs in the GPU code to achieve high
performance for production MD — they are not bugs, just limitations that
are well-known and documented by the AMBER developers."""

        self.console.print()
        self.console.print(Panel(help_text, title="CPU vs. GPU Engine Selection",
                                 border_style="bright_blue", expand=False, padding=(1, 2)))

    def _step5_review_and_save(self):
        """Step 5: Generate summary and save queue data to workspace."""

        profile_loaded = self._cluster_profile is not None
        plan_loaded = self._run_plan is not None

        # Binding header: protocol / cluster / run plan.
        # Protocol names live in the queue's WorkflowConfig map, keyed by
        # workflow_id on each SimulationConfig.
        protocol_names = []
        try:
            wf_map = getattr(self.simulation_queue, '_workflows', {}) or {}
            seen = set()
            for sim in self.simulation_queue.queue:
                wf_id = getattr(sim, 'workflow_id', None)
                if not wf_id or wf_id in seen:
                    continue
                seen.add(wf_id)
                wf = wf_map.get(wf_id)
                if wf and getattr(wf, 'name', ''):
                    protocol_names.append(wf.name)
        except Exception:
            pass
        protocol_label = ", ".join(protocol_names) if protocol_names else "[grey50](unnamed)[/grey50]"

        header_lines = []
        header_lines.append(f"[bold]Protocol:[/bold] {protocol_label}")
        if profile_loaded:
            p = self._cluster_profile
            header_lines.append(
                f"[bold]Cluster profile:[/bold] {p.display_name or p.name} "
                f"[grey50]({p.name}, {p.source})[/grey50]"
            )
        else:
            header_lines.append("[bold]Cluster profile:[/bold] [grey50](none — manual mode)[/grey50]")
        if plan_loaded:
            rp = self._run_plan
            header_lines.append(
                f"[bold]Run plan:[/bold] {rp.name} "
                f"[grey50](cluster={rp.cluster_name}, {rp.source})[/grey50]"
            )
        elif profile_loaded:
            header_lines.append(
                "[bold]Run plan:[/bold] [grey50]not saved — using current Step 5 assignments "
                "(hit 's' in Step 5 to save as a named plan)[/grey50]"
            )
        else:
            header_lines.append("[bold]Run plan:[/bold] [grey50]none[/grey50]")
        self.console.print()
        self.console.print(Panel(
            "\n".join(header_lines),
            title="Configuration",
            border_style="cyan",
            expand=False,
            padding=(0, 1),
        ))

        # Build table data. When a profile is loaded, show the resource
        # class and wall time so the user can audit per-step binding.
        if profile_loaded:
            headers = ["#", "System", "Simulation", "Ensemble", "Sim Time",
                       "Engine", "Resource", "Wall time"]
        else:
            headers = ["#", "System", "Simulation", "Ensemble", "Sim Time", "Engine"]
        rows = []

        # Reuse the same structure name extraction as Step 4
        def extract_structure_name(sim_name):
            parts = sim_name.split('_')
            if len(parts) >= 2:
                step_keywords = ['Energy', 'System', 'NPT', 'NVT', 'Production', 'Equilibration', 'Heating', 'Minimization']
                for idx, part in enumerate(parts):
                    if any(keyword in part for keyword in step_keywords):
                        return '_'.join(parts[:idx])
                return '_'.join(parts[:-1]) if len(parts) > 1 else parts[0]
            return sim_name

        # Group by structure for display (show structure name only on first row)
        structure_groups = {}
        for i, config in enumerate(self.simulation_queue.queue, 1):
            structure_name = extract_structure_name(config.name)
            if structure_name not in structure_groups:
                structure_groups[structure_name] = []
            structure_groups[structure_name].append((i, config))

        for structure_name, sims in structure_groups.items():
            for sim_idx, (global_idx, config) in enumerate(sims):
                # Extract simulation display name (remove structure prefix)
                sim_display_name = config.name
                if config.name.startswith(structure_name + '_'):
                    sim_display_name = config.name[len(structure_name) + 1:]

                # Extract ensemble and sim time from template parameters
                ensemble = ""
                sim_time = ""
                t_config = {}

                # Parse parameters from MDIN template content
                if hasattr(self, 'user_data_manager') and self.user_data_manager and config.template_id:
                    try:
                        content, _ = self.user_data_manager.get_template_content(config.template_id)
                        if content:
                            t_config = self._parse_mdin_params(content)
                    except Exception:
                        pass

                # Apply parameter overrides on top of template config
                if config.parameter_overrides:
                    t_config.update(config.parameter_overrides)

                # Determine ensemble and sim time from config
                try:
                    imin = int(t_config.get('imin', 0))
                except (ValueError, TypeError):
                    imin = 0

                if imin == 1:
                    ensemble = "—"
                    maxcyc = t_config.get('maxcyc', '')
                    sim_time = f"{maxcyc} cycles" if maxcyc else ""
                else:
                    try:
                        ntp = int(t_config.get('ntp', 0))
                    except (ValueError, TypeError):
                        ntp = 0
                    ensemble = "NPT" if ntp > 0 else "NVT"

                    nstlim = t_config.get('nstlim')
                    dt = t_config.get('dt')
                    if nstlim and dt:
                        try:
                            total_ps = float(nstlim) * float(dt)
                            if total_ps >= 1000:
                                sim_time = f"{total_ps / 1000:.1f} ns"
                            else:
                                sim_time = f"{total_ps:.0f} ps"
                        except (ValueError, TypeError):
                            pass

                engine_display = config.engine if config.engine else "[not set]"

                base_row = [
                    str(global_idx),
                    structure_name if sim_idx == 0 else "",
                    sim_display_name,
                    ensemble,
                    sim_time,
                    engine_display
                ]
                if profile_loaded:
                    hw = config.hardware_config or {}
                    resource = hw.get('resource_class', '') or "[grey50][not assigned][/grey50]"
                    wall = hw.get('time_limit', '')
                    if not wall and resource and hw.get('resource_class'):
                        try:
                            cls_def = self._cluster_profile.resolve_class(hw['resource_class'])
                            wall = cls_def.get('default_time', '')
                        except KeyError:
                            wall = ""
                    wall_display = wall or "[grey50](class default)[/grey50]"
                    rows.append(base_row + [resource, wall_display])
                else:
                    rows.append(base_row)

        # Build actions. Options are always in the same order with stable
        # keys; unavailable ones render dim with a reason tag (session-
        # replay-safe: selection is handled, not hidden).
        sbatch_available = shutil.which("sbatch") is not None

        action_4_label = "Write MDIN + SLURM scripts (no launch)"
        if not profile_loaded:
            action_4_label = f"{action_4_label} [grey50 italic](manual wizard — no profile loaded)[/grey50 italic]"
        action_5_label = "Write files and submit via sbatch now"
        if not profile_loaded:
            action_5_label = (
                f"[grey50]{action_5_label}[/grey50] "
                f"[grey50 italic](requires cluster profile)[/grey50 italic]"
            )
        elif not sbatch_available:
            action_5_label = (
                f"[grey50]{action_5_label}[/grey50] "
                f"[grey50 italic](sbatch not on PATH)[/grey50 italic]"
            )

        actions = [
            ("1", "Save setup to workspace (execute later with 'run' command)"),
            ("2", "Write MDIN only (no SLURM, no launch)"),
            ("3", "Write files and run locally now"),
            ("4", action_4_label),
            ("5", action_5_label),
            ("b", "← Back: Hardware Configuration"),
            ("x", "Exit without saving"),
        ]

        # Display using table layout
        self.layout.config_table(
            step_num=5,
            headers=headers,
            rows=rows,
            actions=actions,
            title="Simulation Configuration Summary"
        )

        choice = prompt_with_context(
            self.processor,
            "Enter choice",
            choices=["1", "2", "3", "4", "5", "b", "x"],
            default="1",
            module="MD Manager - Step 5 Review",
            description="Review and save configuration",
            options_map={
                "1": "Save setup for later execution",
                "2": "Write MDIN files only",
                "3": "Write files and run locally now",
                "4": "Write MDIN + SLURM scripts (no launch)",
                "5": "Write files and submit via sbatch now",
                "b": "← Back: Hardware Configuration",
                "x": "Exit without saving",
            }
        )

        if choice == "1":
            self._save_queue_to_workspace()
            self.console.print("[green]✓ Configuration saved successfully[/green]")
            return True
        elif choice == "2":
            self._save_queue_to_workspace()
            self._prepare_simulation_files()
            self.console.print("[green]✓ MDIN files written successfully[/green]")
            return True
        elif choice == "3":
            self._save_queue_to_workspace()
            self.console.print("[green]✓ Configuration saved, switching to execution...[/green]")
            return "EXECUTE_NOW"
        elif choice == "4":
            self._save_queue_to_workspace()
            self._prepare_simulation_files()
            if profile_loaded:
                self._write_slurm_scripts_from_profile(submit=False)
            else:
                # Fall back to the manual SLURM wizard for users who haven't
                # set up a cluster profile yet.
                self._configure_slurm_mode()
            return True
        elif choice == "5":
            if not profile_loaded:
                self.console.print(
                    "[yellow]SLURM generation requires a cluster profile. "
                    "Go back to Step 4 and use option 'p' to load one.[/yellow]"
                )
                return True
            if not sbatch_available:
                self.console.print(
                    "[yellow]sbatch not found on PATH — use option 4 to "
                    "write scripts offline, then submit from a login node.[/yellow]"
                )
                return True
            self._save_queue_to_workspace()
            self._prepare_simulation_files()
            self._write_slurm_scripts_from_profile(submit=True)
            return True
        elif choice == "b":
            return False
        elif choice == "x":
            self.console.print("[yellow]Configuration cancelled[/yellow]")
            return True

        return True

    @staticmethod
    def _parse_mdin_params(mdin_content: str) -> dict:
        """Parse key AMBER parameters from MDIN file content."""
        import re
        params = {}
        # Match quoted values (which may contain ! and ,) or unquoted values
        for match in re.finditer(r"(\w+)\s*=\s*('[^']*'|\"[^\"]*\"|[^,!\n]+)", mdin_content):
            key = match.group(1).strip()
            value = match.group(2).strip().rstrip(',')
            try:
                if '.' in value:
                    params[key] = float(value)
                else:
                    params[key] = int(value)
            except ValueError:
                params[key] = value
        return params

    def _save_queue_to_workspace(self):
        """Save queue configuration to workspace for execution."""
        # Queue is already saved to workspace automatically via SimulationQueue._sync_to_workspace()
        # This method is kept for compatibility but now just ensures sync
        if self.simulation_queue:
            # Force a sync to make sure workspace is up to date
            self.simulation_queue._sync_to_workspace()
            
        # Optionally also save to JSON file for external tools/debugging
        if self.workspace and self.workspace.get('md_save_json_backup', False):
            workspace_queue_file = Path.cwd() / "simulation_queue.json"
            
            queue_data = {
                "simulations": [],
                "created": datetime.now().isoformat(),
                "total_simulations": len(self.simulation_queue)
            }
            
            for config in self.simulation_queue.queue:
                sim_data = {
                    "name": config.name,
                    "template_id": config.template_id,
                    "mdin_path": config.mdin_path,
                    "engine": config.engine,
                    "prmtop": config.prmtop,
                    "rst7": config.rst7
                }
                # Save hardware configuration if available
                if hasattr(config, 'hardware_config') and config.hardware_config:
                    sim_data['hardware_config'] = config.hardware_config
                # Save any additional parameters
                if hasattr(config, 'parameters'):
                    sim_data['parameters'] = config.parameters
                queue_data["simulations"].append(sim_data)
                
            try:
                with open(workspace_queue_file, 'w') as f:
                    json.dump(queue_data, f, indent=2)
            except Exception as e:
                self.console.print(f"[red]Error saving JSON backup: {e}[/red]")

    def _prepare_simulation_files(self):
        """Write all simulation files without launching simulations."""
        if not self.simulation_queue:
            self.console.print("[yellow]No simulations in queue to prepare[/yellow]")
            return False

        # Use unified status checker
        md_status = self._get_md_ready_status()
        if not md_status["ready"]:
            self.console.print("[red]Cannot prepare files: Missing requirements[/red]")
            for missing_item in md_status["missing"]:
                self.console.print(f"  • {missing_item}")
            return False

        self.console.print(f"\n[bold]Preparing Simulation Files[/bold]")
        self.simulation_queue.display_queue(self.console)

        # Count active simulations
        active_count = sum(1 for config in self.simulation_queue.queue if config.status == "active")
        hold_count = len(self.simulation_queue.queue) - active_count

        if active_count == 0:
            self.console.print(f"\n[yellow]No active simulations to prepare (all are on hold)[/yellow]")
            return False
        elif hold_count > 0:
            self.console.print(f"\n[grey50]Will prepare {active_count} active simulations, skip {hold_count} on hold[/grey50]")

        # Group simulations by workflow_id
        workflows = {}  # workflow_id -> list of sim_configs
        standalone = []  # standalone simulations (no workflow_id)

        for sim_config in self.simulation_queue.queue:
            if sim_config.status == "hold":
                continue
            if sim_config.workflow_id:
                if sim_config.workflow_id not in workflows:
                    workflows[sim_config.workflow_id] = []
                workflows[sim_config.workflow_id].append(sim_config)
            else:
                standalone.append(sim_config)

        # Sort workflow steps by workflow_step number
        for workflow_id in workflows:
            workflows[workflow_id].sort(key=lambda x: x.workflow_step or 0)

        # Create run directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = Path.cwd() / "simulations" / f"batch_{timestamp}"
        run_dir.mkdir(parents=True, exist_ok=True)

        self.console.print(f"\n[bold cyan]Created run directory:[/bold cyan] {run_dir}")

        original_cwd = Path.cwd()

        import os
        import shutil

        prepared_count = 0
        skipped_count = 0

        # For single-structure workflows, write directly into batch dir (no extra nesting)
        num_groups = len(workflows) + len(standalone)

        try:
            # Prepare workflows first
            for workflow_id, workflow_steps in workflows.items():
                # Get structure files from this workflow's first simulation
                wf_first = workflow_steps[0]
                prmtop_file = Path(wf_first.prmtop).resolve()
                coord_file = Path(wf_first.rst7).resolve()

                # Use structure label for directory naming when available
                structure_label = wf_first.structure_label
                wf_config = self.simulation_queue._workflows.get(workflow_id)
                protocol_name = structure_label or (wf_config.name if wf_config and wf_config.name else None)

                if num_groups == 1:
                    # Single workflow — write directly into batch dir
                    self._prepare_workflow_files(run_dir, workflow_steps, prmtop_file, coord_file, original_cwd, nest=False, protocol_name=protocol_name)
                else:
                    self._prepare_workflow_files(run_dir, workflow_steps, prmtop_file, coord_file, original_cwd, nest=True, protocol_name=protocol_name)
                prepared_count += len(workflow_steps)

            # Prepare standalone simulations
            for i, sim_config in enumerate(standalone):
                sa_prmtop = Path(sim_config.prmtop).resolve()
                sa_coord = Path(sim_config.rst7).resolve()
                self._prepare_standalone_simulation(run_dir, sim_config, i, sa_prmtop, sa_coord, original_cwd)
                prepared_count += 1

        finally:
            os.chdir(original_cwd)

        # Summary
        self.console.print(f"\n[bold green]File Preparation Complete[/bold green]")
        self.console.print(f"  ✓ Prepared: {prepared_count} simulation(s)")
        if skipped_count > 0:
            self.console.print(f"  ⊘ Skipped: {skipped_count} simulation(s) (on hold)")
        self.console.print(f"\n[cyan]Files written to:[/cyan] {run_dir.relative_to(original_cwd)}")

        if workflows:
            self.console.print(f"\n[bold]Protocol directories contain:[/bold]")
            self.console.print(f"  • stepN_*.mdin - AMBER input files for each step")
            self.console.print(f"  • *.prmtop - Topology file (shared by all steps)")
            self.console.print(f"  • *.rst7 - Initial coordinates (for step 1)")
            self.console.print(f"  • run_workflow.sh - Master script that chains all steps")
            self.console.print(f"  • README.txt - Protocol details")

        if standalone:
            self.console.print(f"\n[bold]Standalone simulation directories contain:[/bold]")
            self.console.print(f"  • simulation.mdin - AMBER input file")
            self.console.print(f"  • run_simulation.sh - Executable script")
            self.console.print(f"  • README.txt - Simulation details")

        self.console.print(f"\nTo execute: cd into each directory and run the .sh script")
        self.console.print(f"Or use 'Execute simulations' from the main menu")

        return True

    def _prepare_workflow_files(self, run_dir: Path, workflow_steps: List, prmtop_file: Path,
                                coord_file: Path, original_cwd: Path, nest: bool = True,
                                protocol_name: str = None):
        """Prepare files for a multi-step workflow in a single directory."""
        import shutil

        # Use user-assigned protocol name, fall back to prmtop stem
        workflow_name = protocol_name or prmtop_file.stem

        # Create workflow directory (or use run_dir directly for single-structure case)
        if nest:
            workflow_dir = run_dir / f"workflow_{workflow_name.replace(' ', '_')}"
            workflow_dir.mkdir(exist_ok=True)
        else:
            workflow_dir = run_dir

        self.console.print(f"\n[bold blue]Preparing protocol: {workflow_name} ({len(workflow_steps)} steps)[/bold blue]")

        # Copy topology and initial coordinates once
        prmtop_name = prmtop_file.name
        coord_name = coord_file.name

        shutil.copy2(prmtop_file, workflow_dir / prmtop_name)
        shutil.copy2(coord_file, workflow_dir / coord_name)
        self.console.print(f"[grey50]  ✓ Copied structure files[/grey50]")

        # CpHMD: we don't stage cpin/modified-prmtop at the workflow level any
        # more; they're copied into each production step's subdir below, which
        # matches the SLURM path layout.
        wf_config = None
        if workflow_steps and workflow_steps[0].workflow_id:
            wf_config = self.simulation_queue._workflows.get(workflow_steps[0].workflow_id)

        # Prepare each step's subdir + mdin file. Layout mirrors SLURM:
        #   workflow_dir/
        #     <prmtop>, <coord>.rst7, run_workflow.sh, README.txt
        #     step1/simulation.mdin
        #     step2/simulation.mdin
        #     ...
        # mdin_files stores each step's relative-from-workflow-dir path so the
        # README and run-script builders can reference them consistently.
        mdin_files = []
        step_names = []

        for i, sim_config in enumerate(workflow_steps, 1):
            step_key = f"step{i}"
            step_dir = workflow_dir / step_key
            step_dir.mkdir(exist_ok=True)

            mdin_file = step_dir / "simulation.mdin"
            rel_mdin = f"{step_key}/simulation.mdin"

            template_content = self._resolve_mdin_content(sim_config)

            if template_content:
                # Per-step sim_dir so DISANG/DUMPAVE/LISTOUT paths written by
                # _apply_configured_restraints land next to the simulation.
                template_content = self._apply_configured_restraints(
                    template_content, sim_config, step_dir
                )
                with open(mdin_file, 'w') as f:
                    f.write(template_content)
                mdin_files.append(rel_mdin)
                step_names.append(sim_config.name)
                self.console.print(f"[grey50]  ✓ Created {rel_mdin}[/grey50]")
            else:
                self._create_basic_mdin_file(mdin_file)
                mdin_files.append(rel_mdin)
                step_names.append(sim_config.name)
                self.console.print(f"[yellow]  ⚠ Created basic {rel_mdin}[/yellow]")

            # Stage CpHMD files per production step (mirrors SLURM).
            if wf_config and wf_config.cpin_file and self._is_production_step(sim_config):
                cpin_path = Path(wf_config.cpin_file)
                if cpin_path.exists():
                    shutil.copy2(cpin_path, step_dir / cpin_path.name)
                    self.console.print(f"[grey50]  ✓ Staged CPIN into {step_key}/ ({cpin_path.name})[/grey50]")
                if wf_config.cpin_config:
                    mod_prmtop = wf_config.cpin_config.get('modified_prmtop')
                    if mod_prmtop and Path(mod_prmtop).exists():
                        shutil.copy2(mod_prmtop, step_dir / Path(mod_prmtop).name)
                        self.console.print(f"[grey50]  ✓ Staged modified topology into {step_key}/[/grey50]")

        # Check if last step is production - offer extended production cycles
        extended_production_cycles = 0
        last_step_type = self._get_template_type(workflow_steps[-1].template_id) if workflow_steps else ""

        if last_step_type == "production":
            from rich.prompt import Confirm, IntPrompt
            if confirm_with_context(
                self.processor,
                "\n[cyan]Add extended production run cycles?[/cyan]",
                default=False,
                module="MD Manager - Workflow",
                description="Add extended production run cycles",
            ):
                extended_production_cycles = int_prompt_with_context(
                    self.processor,
                    "Number of additional production cycles",
                    default=10,
                    module="MD Manager - Workflow",
                    description="Number of additional production cycles",
                )
                self.console.print(f"[green]Will add {extended_production_cycles} extended production cycles[/green]")

        # Create master run script that chains all steps
        self._create_workflow_run_script(workflow_dir, workflow_steps, mdin_files,
                                        prmtop_name, coord_name, extended_production_cycles)
        self.console.print(f"[grey50]  ✓ Created master run script (run_workflow.sh)[/grey50]")

        # Create README
        readme_file = workflow_dir / "README.txt"
        with open(readme_file, 'w') as f:
            f.write(f"Protocol: {workflow_name}\n")
            f.write(f"Structure: {prmtop_file.stem}\n")
            f.write(f"Steps: {len(workflow_steps)}\n")
            f.write(f"Topology: {prmtop_name}\n")
            f.write(f"Initial Coordinates: {coord_name}\n")
            f.write(f"\nProtocol Steps:\n")
            for i, (step_name, mdin_file) in enumerate(zip(step_names, mdin_files), 1):
                f.write(f"  {i}. {step_name} ({mdin_file})\n")
            # Collect hardware summary across all steps
            engine_counts = {}
            gpu_ids_set = set()
            mpi_tasks_val = None
            for step in workflow_steps:
                hw = step.hardware_config or {}
                engine = hw.get('engine', step.engine)
                engine_counts[engine] = engine_counts.get(engine, 0) + 1
                if hw.get('gpu_ids') is not None:
                    gpu_ids_set.add(hw['gpu_ids'])
                if hw.get('mpi_tasks') is not None:
                    mpi_tasks_val = hw['mpi_tasks']
            f.write(f"\nHardware Configuration:\n")
            for engine, count in engine_counts.items():
                f.write(f"  {count}x {engine}\n")
            if mpi_tasks_val is not None:
                f.write(f"  mpi_tasks: {mpi_tasks_val}\n")
            if gpu_ids_set:
                f.write(f"  gpu_ids: {', '.join(sorted(gpu_ids_set))}\n")
        self.console.print(f"[grey50]  ✓ Created README[/grey50]")

        self.console.print(f"[green]✓ Protocol prepared in: {workflow_dir.relative_to(original_cwd)}[/green]")

    def _prepare_standalone_simulation(self, run_dir: Path, sim_config, index: int,
                                      prmtop_file: Path, coord_file: Path, original_cwd: Path):
        """Prepare files for a standalone (non-workflow) simulation."""
        import shutil

        self.console.print(f"\n[bold blue]Preparing standalone: {sim_config.name}[/bold blue]")

        # Create simulation directory
        sim_dir = run_dir / f"{index+1:02d}_{sim_config.name.replace(' ', '_')}"
        sim_dir.mkdir(exist_ok=True)

        prmtop_name = prmtop_file.name
        coord_name = coord_file.name

        try:
            # Copy topology and coordinates
            shutil.copy2(prmtop_file, sim_dir / prmtop_name)
            shutil.copy2(coord_file, sim_dir / coord_name)

            # Create MDIN file
            mdin_file = sim_dir / "simulation.mdin"
            template_content = self._resolve_mdin_content(sim_config)

            if template_content:
                # Apply restraints if configured
                if sim_config.restraints:
                    if 'group' in sim_config.restraints:
                        template_content = self._apply_group_to_template(
                            template_content,
                            sim_config.restraints['group']
                        )
                        self.console.print(f"[grey50]  ✓ Applied GROUP restraints[/grey50]")
                    elif 'restraintmask' in sim_config.restraints:
                        template_content = self._apply_restraintmask_to_template(
                            template_content,
                            sim_config.restraints['restraintmask']
                        )
                        self.console.print(f"[grey50]  ✓ Applied positional restraints[/grey50]")
                    if 'disang' in sim_config.restraints:
                        template_content = self._apply_disang_to_template(
                            template_content,
                            sim_config.restraints['disang'],
                            sim_dir
                        )
                        self.console.print(f"[grey50]  ✓ Applied DISANG restraints[/grey50]")

                with open(mdin_file, 'w') as f:
                    f.write(template_content)
                self.console.print(f"[grey50]  ✓ Created MDIN file[/grey50]")
            else:
                self._create_basic_mdin_file(mdin_file)
                self.console.print(f"[yellow]  ⚠ Created basic MDIN file[/yellow]")

            # Create README
            readme_file = sim_dir / "README.txt"
            hw_config = sim_config.hardware_config or {}
            with open(readme_file, 'w') as f:
                f.write(f"Simulation: {sim_config.name}\n")
                f.write(f"Template: {sim_config.template_id}\n")
                f.write(f"Topology: {prmtop_name}\n")
                f.write(f"Coordinates: {coord_name}\n")
                f.write(f"\nHardware Configuration:\n")
                f.write(f"  engine: {hw_config.get('engine', sim_config.engine)}\n")
                if hw_config.get('gpu_ids') is not None:
                    f.write(f"  gpu_ids: {hw_config['gpu_ids']}\n")
                if hw_config.get('mpi_tasks') is not None:
                    f.write(f"  mpi_tasks: {hw_config['mpi_tasks']}\n")
            self.console.print(f"[grey50]  ✓ Created README[/grey50]")

            # Generate run script
            output_prefix = sim_config.name.replace(' ', '_').lower()
            mdout_file = sim_dir / f"{output_prefix}.mdout"
            restart_file = sim_dir / f"{output_prefix}.rst7"
            trajectory_file = sim_dir / f"{output_prefix}.nc"
            reference_file = sim_dir / coord_name

            cmd = self._build_amber_command_line(
                engine=sim_config.engine,
                mdin_file=mdin_file,
                topology_file=sim_dir / prmtop_name,
                coordinate_file=sim_dir / coord_name,
                output_file=mdout_file,
                restart_file=restart_file,
                trajectory_file=trajectory_file,
                reference_file=reference_file,
                hardware_config=hw_config
            )

            self._create_run_script(sim_dir, cmd, sim_config, hw_config)
            self.console.print(f"[grey50]  ✓ Created run script (run_simulation.sh)[/grey50]")

            self.console.print(f"[green]✓ Files prepared in: {sim_dir.relative_to(original_cwd)}[/green]")

        except Exception as e:
            self.console.print(f"[red]✗ Error preparing {sim_config.name}: {e}[/red]")

    def _step0_structure_file_selection(self):
        """Step 0: Select prmtop/rst7 structure file pairs for simulations."""

        # Load existing structure pairs from workspace if available
        structure_pairs = []
        if self.workspace and self.workspace.has('md_structure_pairs'):
            saved_pairs = self.workspace.get('md_structure_pairs', [])
            for pair_data in saved_pairs:
                # Verify files still exist
                prmtop_path = Path(pair_data['prmtop'])
                rst7_path = Path(pair_data['rst7'])
                if prmtop_path.exists() and rst7_path.exists():
                    structure_pairs.append({
                        'name': pair_data['name'],
                        'prmtop': prmtop_path,
                        'rst7': rst7_path
                    })

        while True:
            # Scan for available structure files and detect smart pairs
            workspace_files, directory_files = self._scan_structure_files_enhanced()
            suggested_pairs = self._detect_smart_structure_pairs(workspace_files, directory_files)

            # Build status lines for display
            status_lines = self._build_step0_status_lines(structure_pairs, suggested_pairs)

            # Build actions menu (combines options + navigation)
            actions, option_map = self._build_step0_actions(structure_pairs, suggested_pairs, workspace_files, directory_files)

            # Display step header and status
            self.console.print(self.layout.step_header(step_num=0))
            self.console.print("[grey50]  Each simulation requires a topology (.prmtop) defining the system's atoms,[/grey50]")
            self.console.print("[grey50]  bonds, and force field parameters, and a coordinate file (.rst7/.inpcrd)[/grey50]")
            self.console.print("[grey50]  defining atomic positions.[/grey50]")
            self.console.print()

            for line in status_lines:
                self.console.print(f"  {line}")

            if status_lines:
                self.console.print()

            # Display options and navigation with section headers
            self._display_step0_menu(actions, option_map, structure_pairs)

            # Get user choice
            valid_choices = list(option_map.keys())

            # Build options_map for context (convert tuples to strings)
            context_options_map = {}
            for key, value in actions:
                if isinstance(option_map.get(key), tuple):
                    # Extract readable description from tuple
                    context_options_map[key] = value
                else:
                    context_options_map[key] = value

            choice = prompt_with_context(
                self.processor,
                "Enter choice",
                choices=valid_choices,
                module="MD Manager - Step 0 Structure Files",
                description="Structure file selection",
                options_map=context_options_map
            )
            
            if choice == "x":
                return None
            elif choice == "n":
                if structure_pairs:
                    # Save structure pairs to workspace before returning
                    if self.workspace:
                        workspace_pairs_data = []
                        for pair in structure_pairs:
                            workspace_pairs_data.append({
                                'name': pair['name'],
                                'prmtop': str(pair['prmtop']),
                                'rst7': str(pair['rst7'])
                            })
                        self.workspace.set('md_structure_pairs', workspace_pairs_data)
                    return structure_pairs
                else:
                    self.console.print("[yellow]Please add at least one structure file pair first[/yellow]")
                    continue
            
            # Handle action-based choices
            action = option_map[choice]
            if isinstance(action, tuple) and action[0] == "use_suggested":
                # Use specific suggested pair
                pair = suggested_pairs[action[1]]
                structure_pairs.append({
                    'name': pair['name'],
                    'prmtop': pair['prmtop'],
                    'rst7': pair['rst7']
                })
                self.console.print(f"[green]✓ Added {pair['name']} ({pair['match_type']})[/green]")
                
            elif action == "use_all_suggested":
                # Use all suggested pairs
                for pair in suggested_pairs:
                    structure_pairs.append({
                        'name': pair['name'],
                        'prmtop': pair['prmtop'],
                        'rst7': pair['rst7']
                    })
                self.console.print(f"[green]✓ Added {len(suggested_pairs)} suggested pairs[/green]")
                
            elif action == "create_custom":
                # Manual pair creation from available files
                pair = self._create_custom_structure_pair(workspace_files, directory_files)
                if pair:
                    structure_pairs.append(pair)
                    self.console.print(f"[green]✓ Added custom pair: {pair['name']}[/green]")
                    
            elif action == "browse_files":
                # Browse for files not in current locations
                pair = self._browse_for_structure_pair()
                if isinstance(pair, list):
                    # Multiple pairs from "find pairs"
                    structure_pairs.extend(pair)
                    self.console.print(f"[green]✓ Added {len(pair)} matching pairs[/green]")
                elif pair:
                    # Single pair from manual selection
                    structure_pairs.append(pair)
                    self.console.print(f"[green]✓ Added structure pair: {pair['name']}[/green]")
                    
            elif action == "remove":
                self._remove_structure_pair(structure_pairs)
            elif action == "help":
                self._show_step0_help()

    def _build_step0_status_lines(self, structure_pairs, suggested_pairs):
        """Build status lines for Step 0 display."""
        status_lines = []

        # Current pairs
        if structure_pairs:
            status_lines.append(f"Structure Pairs: {len(structure_pairs)} selected")
            for i, pair in enumerate(structure_pairs, 1):
                prmtop_name = pair['prmtop'].name
                rst7_name = pair['rst7'].name
                status_lines.append(f"  {i}. {prmtop_name} + {rst7_name}")
        else:
            status_lines.append("Structure Pairs: None selected")

        # Suggested pairs
        if suggested_pairs:
            status_lines.append("")
            status_lines.append(f"Auto-detected: {len(suggested_pairs)} pairs available")

        return status_lines

    def _display_step0_menu(self, actions, option_map, structure_pairs):
        """Display Step 0 menu with section headers matching Step 2 style."""
        # Separate options from navigation
        navigation_keys = {'n', 'h', 'x'}
        options = [(k, d) for k, d in actions if k not in navigation_keys]
        navigation = [(k, d) for k, d in actions if k in navigation_keys]

        # Display options section
        if options:
            self.console.print("What would you like to do?", style="bold")
            for key, description in options:
                self.console.print(f"  {key:>2}. {description}")
            self.console.print()

        # Display navigation section with cyan header
        self.console.print("[bold cyan]Navigation:[/bold cyan]")
        for key, description in navigation:
            # Determine style based on key
            if key == 'x':
                style = "dim"
            elif key == 'n':
                style = "green"
            elif key == 'h':
                style = "cyan"
            else:
                style = "white"

            key_str = f"[{style}]{key}[/{style}]"
            self.console.print(f"   {key_str}. {description}")
        self.console.print()

    def _build_step0_actions(self, structure_pairs, suggested_pairs, workspace_files, directory_files):
        """Build actions menu with option mapping for Step 0."""
        actions = []
        option_map = {}

        # Suggested pair options
        for i, pair in enumerate(suggested_pairs):
            key = str(len(actions) + 1)
            match_type = pair.get('match_type', 'exact match')
            actions.append((key, f"Use {pair['name']} ({match_type})"))
            option_map[key] = ("use_suggested", i)

        # Use all suggested
        if len(suggested_pairs) > 1:
            key = str(len(actions) + 1)
            actions.append((key, f"Use all suggested pairs ({len(suggested_pairs)})"))
            option_map[key] = "use_all_suggested"

        # Manual options
        has_topology = any(workspace_files.get('prmtop', [])) or any(directory_files.get('prmtop', []))
        has_coordinates = any(workspace_files.get('rst7', [])) or any(directory_files.get('rst7', []))

        if has_topology and has_coordinates:
            key = str(len(actions) + 1)
            actions.append((key, "Create custom pair (manual selection)"))
            option_map[key] = "create_custom"

        # Browse for files
        key = str(len(actions) + 1)
        actions.append((key, "Browse for structure files (topology + coordinates)"))
        option_map[key] = "browse_files"

        # Remove pair option
        if structure_pairs:
            key = str(len(actions) + 1)
            actions.append((key, "Remove structure pair"))
            option_map[key] = "remove"

        # Navigation
        if structure_pairs:
            actions.append(("n", "Next step (Protocol Selection)"))
            option_map["n"] = "next"
        else:
            # Still allow 'n' but it will show warning
            option_map["n"] = "next"

        actions.append(("h", "Help"))
        option_map["h"] = "help"

        actions.append(("x", "Exit setup"))
        option_map["x"] = "exit"

        return actions, option_map

    def _show_step0_help(self):
        """Show help for Step 0."""
        help_text = """
[bold cyan]Structure File Selection Help[/bold cyan]

[bold]What are structure files?[/bold]
MD simulations require TWO files per structure:

1. [cyan]Topology file[/cyan] (.prmtop, .parm7)
   → Contains system parameters, atom types, and force field information

2. [cyan]Coordinate file[/cyan] (.rst7, .inpcrd)
   → Contains 3D positions of all atoms in the system

[bold]Selection Process:[/bold]
• First, select a topology file (.prmtop/.parm7)
• Then, ProPrep will search for matching coordinate files (.rst7/.inpcrd)

[bold]Auto-detection:[/bold]
• ProPrep automatically finds matching pairs in your workspace and working directory
• Match types:
  - [green]exact match[/green]: Same basename (e.g., system.prmtop + system.rst7)
  - [yellow]modified[/yellow]: Modified topology for constant pH MD (system_modified.prmtop + system.rst7)

[bold]Tips:[/bold]
• Use suggested pairs for quick setup
• Browse files if your structure files are in a different location
• You can select multiple pairs for batch simulations

[grey50]Press Enter to continue...[/grey50]
"""
        self.console.print(Panel(help_text, border_style="bright_blue", padding=(1, 2)))
        input()

    def _scan_structure_files_enhanced(self):
        """
        Scan for structure files in workspace and working directory separately.
        Supports both .prmtop/.parm7 and .rst7/.inpcrd extensions.
        """
        from pathlib import Path
        
        current_dir = Path.cwd()
        
        # Scan workspace (files from ProPrep modules)
        workspace_files = {
            'prmtop': [],
            'rst7': []
        }
        
        if self.workspace:
            # Get workspace files from various ProPrep modules. Modules may
            # store the same physical file under multiple keys (canonical
            # `parm7_file`/`rst7_file` plus descriptive aliases like
            # `prmtop_titrated`/`rst7_titrated` from pb_titrate). Dedup by
            # resolved path so the same file isn't surfaced twice to the
            # pair-builder.
            seen_prmtop_resolved = set()
            seen_rst7_resolved = set()
            for key, paths in self.workspace.items():
                if isinstance(paths, str) and (paths.endswith('.prmtop') or paths.endswith('.parm7')):
                    p = Path(paths)
                    if not p.exists():
                        continue
                    resolved = p.resolve()
                    if resolved in seen_prmtop_resolved:
                        continue
                    seen_prmtop_resolved.add(resolved)
                    workspace_files['prmtop'].append(p)
                elif isinstance(paths, str) and (paths.endswith('.rst7') or paths.endswith('.inpcrd')):
                    p = Path(paths)
                    if not p.exists():
                        continue
                    resolved = p.resolve()
                    if resolved in seen_rst7_resolved:
                        continue
                    seen_rst7_resolved.add(resolved)
                    workspace_files['rst7'].append(p)

        # Scan current directory
        directory_files = {
            'prmtop': list(current_dir.glob("*.prmtop")) + list(current_dir.glob("*.parm7")),
            'rst7': list(current_dir.glob("*.rst7")) + list(current_dir.glob("*.inpcrd"))
        }

        return workspace_files, directory_files
    
    def _detect_smart_structure_pairs(self, workspace_files, directory_files):
        """
        Detect structure pairs using smart matching:
        1. Exact basename matches
        2. Modified topology matches (base_name_* + base_name)

        Deduplicates files that appear in both workspace and directory.
        Workspace files take priority over directory files.
        """
        pairs = []

        # Deduplicate files: workspace takes priority over directory
        # Use resolved paths for comparison
        workspace_prmtop_paths = {f.resolve() for f in workspace_files['prmtop']}
        workspace_rst7_paths = {f.resolve() for f in workspace_files['rst7']}

        # Build deduplicated lists with source tracking
        all_prmtop = []
        prmtop_sources = {}
        for f in workspace_files['prmtop']:
            resolved = f.resolve()
            all_prmtop.append(f)
            prmtop_sources[resolved] = "workspace"
        for f in directory_files['prmtop']:
            resolved = f.resolve()
            if resolved not in workspace_prmtop_paths:
                all_prmtop.append(f)
                prmtop_sources[resolved] = "directory"

        all_rst7 = []
        rst7_sources = {}
        for f in workspace_files['rst7']:
            resolved = f.resolve()
            all_rst7.append(f)
            rst7_sources[resolved] = "workspace"
        for f in directory_files['rst7']:
            resolved = f.resolve()
            if resolved not in workspace_rst7_paths:
                all_rst7.append(f)
                rst7_sources[resolved] = "directory"

        # For each rst7, find the best matching prmtop
        # Priority: _cpin version > exact match > other modified versions
        for rst7 in all_rst7:
            rst7_base = rst7.stem
            rst7_resolved = rst7.resolve()
            rst7_source = rst7_sources[rst7_resolved]

            # Find all potential matches for this rst7
            exact_match = None
            cpin_match = None
            other_modified = []

            for prmtop in all_prmtop:
                prmtop_base = prmtop.stem
                prmtop_resolved = prmtop.resolve()
                prmtop_source = prmtop_sources[prmtop_resolved]

                # Determine combined source
                if prmtop_source == "workspace" and rst7_source == "workspace":
                    source = "workspace"
                elif prmtop_source == "workspace" or rst7_source == "workspace":
                    source = "mixed"
                else:
                    source = "directory"

                # Exact match
                if prmtop_base == rst7_base:
                    exact_match = {
                        'name': rst7_base,
                        'prmtop': prmtop,
                        'rst7': rst7,
                        'match_type': "exact match",
                        'source': source
                    }

                # Modified topology match (topology starts with rst7_base + "_")
                elif prmtop_base.startswith(f"{rst7_base}_"):
                    suffix = prmtop_base[len(rst7_base)+1:]
                    match_info = {
                        'name': rst7_base,  # Use base name for cpin matches
                        'prmtop': prmtop,
                        'rst7': rst7,
                        'match_type': f"modified ({suffix})",
                        'source': source
                    }
                    if suffix == "cpin":
                        cpin_match = match_info
                        cpin_match['match_type'] = "cpin-ready"  # Mark as constant pH ready
                    else:
                        other_modified.append(match_info)

            # Add the best match: cpin > exact > other modified
            if cpin_match:
                pairs.append(cpin_match)
            elif exact_match:
                pairs.append(exact_match)
            elif other_modified:
                # Add first other modified match
                pairs.append(other_modified[0])

        # Sort pairs: cpin-ready first, then exact matches, then by source preference
        def sort_key(pair):
            if pair['match_type'] == "cpin-ready":
                match_priority = 0
            elif pair['match_type'] == "exact match":
                match_priority = 1
            else:
                match_priority = 2
            source_priority = {"workspace": 0, "mixed": 1, "directory": 2}[pair['source']]
            return (match_priority, source_priority, pair['name'])

        pairs.sort(key=sort_key)
        return pairs
    
    def _display_structure_files_by_source(self, workspace_files, directory_files, suggested_pairs):
        """Display available structure files organized by source."""
        
        # Always show what files were found, even if empty
        total_files = (len(workspace_files['prmtop']) + len(workspace_files['rst7']) + 
                      len(directory_files['prmtop']) + len(directory_files['rst7']))
        
        if total_files == 0:
            self.console.print(f"\n[yellow]No structure files found in workspace or working directory[/yellow]")
            self.console.print(f"[grey50]Working directory: {self.working_directory}[/grey50]")
            self.console.print(f"[grey50]Looking for: *.prmtop, *.parm7, *.rst7, *.inpcrd[/grey50]")
            return
        
        # Show workspace files
        if workspace_files['prmtop'] or workspace_files['rst7']:
            self.console.print(f"\n📁 [bold]Workspace Files:[/bold]")
            if workspace_files['prmtop']:
                topo_names = [f.name for f in workspace_files['prmtop']]
                self.console.print(f"  TOPOLOGY: {', '.join(topo_names)}")
            if workspace_files['rst7']:
                coord_names = [f.name for f in workspace_files['rst7']]
                self.console.print(f"  COORDS: {', '.join(coord_names)}")
        
        # Show suggested pairs
        if suggested_pairs:
            self.console.print(f"\n✅ [bold]Detected Structure Pairs:[/bold]")
            for i, pair in enumerate(suggested_pairs, 1):
                prmtop_name = pair['prmtop'].name
                rst7_name = pair['rst7'].name
                match_type = pair['match_type']
                source = pair['source']
                self.console.print(f"  [{i}] {prmtop_name} + {rst7_name} ({match_type}, {source})")
        
        # Show unpaired files if any
        unpaired_prmtop = []
        unpaired_rst7 = []
        all_prmtop = workspace_files['prmtop'] + directory_files['prmtop']
        all_rst7 = workspace_files['rst7'] + directory_files['rst7']
        
        for prmtop in all_prmtop:
            if not any(pair['prmtop'] == prmtop for pair in suggested_pairs):
                unpaired_prmtop.append(prmtop)
                
        for rst7 in all_rst7:
            if not any(pair['rst7'] == rst7 for pair in suggested_pairs):
                unpaired_rst7.append(rst7)
        
        if unpaired_prmtop or unpaired_rst7:
            self.console.print(f"\n⚠️  [bold]Unpaired Files:[/bold]")
            if unpaired_prmtop:
                names = [f.name for f in unpaired_prmtop]
                self.console.print(f"  TOPOLOGY: {', '.join(names)}")
            if unpaired_rst7:
                names = [f.name for f in unpaired_rst7]
                self.console.print(f"  COORDS: {', '.join(names)}")
    
    def _create_custom_structure_pair(self, workspace_files, directory_files):
        """Create a custom structure pair by manually selecting files."""
        all_prmtop = workspace_files['prmtop'] + directory_files['prmtop']
        all_rst7 = workspace_files['rst7'] + directory_files['rst7']
        
        if not all_prmtop or not all_rst7:
            self.console.print("[red]Need both topology and coordinate files to create a pair[/red]")
            return None
        
        from rich.prompt import IntPrompt
        
        # Select topology file
        self.console.print(f"\n[bold]Select Topology File:[/bold]")
        for i, prmtop in enumerate(all_prmtop, 1):
            source = "workspace" if prmtop in workspace_files['prmtop'] else "directory"
            self.console.print(f"  {i}. {prmtop.name} ({source})")

        try:
            topo_choice_str = prompt_with_context(
                self.processor,
                "Topology file",
                choices=[str(i) for i in range(1, len(all_prmtop)+1)],
                module="MD Manager - Step 0 Structure Files",
                description="Select topology file"
            )
            topo_choice = int(topo_choice_str)
            selected_prmtop = all_prmtop[topo_choice - 1]
        except:
            self.console.print("[yellow]Selection cancelled[/yellow]")
            return None

        # Select coordinate file
        self.console.print(f"\n[bold]Select Coordinate File:[/bold]")
        for i, rst7 in enumerate(all_rst7, 1):
            source = "workspace" if rst7 in workspace_files['rst7'] else "directory"
            self.console.print(f"  {i}. {rst7.name} ({source})")

        try:
            coord_choice_str = prompt_with_context(
                self.processor,
                "Coordinate file",
                choices=[str(i) for i in range(1, len(all_rst7)+1)],
                module="MD Manager - Step 0 Structure Files",
                description="Select coordinate file"
            )
            coord_choice = int(coord_choice_str)
            selected_rst7 = all_rst7[coord_choice - 1]
        except:
            self.console.print("[yellow]Selection cancelled[/yellow]")
            return None
        
        # Create pair
        pair_name = f"{selected_prmtop.stem}_{selected_rst7.stem}"
        return {
            'name': pair_name,
            'prmtop': selected_prmtop,
            'rst7': selected_rst7
        }
    
    def _browse_for_structure_pair(self):
        """Browse for structure files not in current locations."""
        # Reuse existing browse functionality
        return self._select_structure_pair()

    def _find_matching_pairs_in_dir(self, search_dir):
        """Find matching prmtop/rst7 pairs from a specific directory."""
        self.console.print(f"\n[bold cyan]Finding Matching Structure Pairs[/bold cyan]")
        self.console.print(f"[grey50]Searching recursively in {search_dir}...[/grey50]")
        
        # Find all prmtop and rst7 files recursively
        prmtop_files = list(search_dir.rglob("*.prmtop"))
        rst7_files = list(search_dir.rglob("*.rst7"))
        
        if not prmtop_files or not rst7_files:
            self.console.print("[yellow]No prmtop or rst7 files found in this directory tree[/yellow]")
            if prmtop_files:
                self.console.print(f"[grey50]Found {len(prmtop_files)} prmtop files but no rst7 files[/grey50]")
            elif rst7_files:
                self.console.print(f"[grey50]Found {len(rst7_files)} rst7 files but no prmtop files[/grey50]")
            input("Press Enter to continue browsing...")
            return None
        
        # Find and display matching pairs
        matching_pairs = self._match_structure_files(prmtop_files, rst7_files, search_dir)
        
        if not matching_pairs:
            self.console.print("[yellow]No matching pairs found (files with same base name)[/yellow]")
            self.console.print(f"[grey50]Found {len(prmtop_files)} prmtop and {len(rst7_files)} rst7 files, but no matching base names[/grey50]")
            input("Press Enter to continue browsing...")
            return None
        
        # Display and let user select
        selected = self._select_from_found_pairs(matching_pairs, search_dir)
        return selected
    
    def _find_matching_structure_pairs(self):
        """Find prmtop/rst7 files with matching base names."""
        from pathlib import Path
        
        self.console.print(f"\n[bold cyan]Finding Matching Structure Pairs[/bold cyan]")
        
        # Ask for directory to search
        search_dir_input = prompt_with_context(
            self.processor, "Enter directory to search (or press Enter for current directory)",
            default="",
            module="MD Manager",
            description="Directory to search for structure pairs",
        ).strip()
        if search_dir_input:
            search_dir = Path(search_dir_input).expanduser()
            if not search_dir.exists() or not search_dir.is_dir():
                self.console.print(f"[red]Directory not found: {search_dir}[/red]")
                return []
        else:
            search_dir = Path.cwd()
        
        self.console.print(f"[grey50]Searching in {search_dir}...[/grey50]")
        
        # Find all prmtop and rst7 files recursively
        prmtop_files = list(search_dir.rglob("*.prmtop"))
        rst7_files = list(search_dir.rglob("*.rst7"))
        
        if not prmtop_files or not rst7_files:
            self.console.print("[yellow]No prmtop or rst7 files found[/yellow]")
            return []
        
        # Find matching pairs
        matching_pairs = self._match_structure_files(prmtop_files, rst7_files, search_dir)
        
        if not matching_pairs:
            self.console.print("[yellow]No matching pairs found (files with same base name)[/yellow]")
            return []
        
        # Display and let user select
        return self._select_from_found_pairs(matching_pairs, search_dir)
    
    def _match_structure_files(self, prmtop_files, rst7_files, base_dir):
        """Match prmtop and rst7 files by base name."""
        # Create dictionaries indexed by base name
        prmtop_by_base = {}
        for prmtop in prmtop_files:
            base_name = prmtop.stem
            if base_name not in prmtop_by_base:
                prmtop_by_base[base_name] = []
            prmtop_by_base[base_name].append(prmtop)
        
        rst7_by_base = {}
        for rst7 in rst7_files:
            base_name = rst7.stem
            if base_name not in rst7_by_base:
                rst7_by_base[base_name] = []
            rst7_by_base[base_name].append(rst7)
        
        # Find matching pairs
        matching_pairs = []
        for base_name in prmtop_by_base:
            if base_name in rst7_by_base:
                # For each combination of matching files
                for prmtop in prmtop_by_base[base_name]:
                    for rst7 in rst7_by_base[base_name]:
                        # Prefer pairs in the same directory
                        if prmtop.parent == rst7.parent:
                            matching_pairs.append({
                                'base_name': base_name,
                                'prmtop': prmtop,
                                'rst7': rst7,
                                'name': base_name,
                                'directory': prmtop.parent
                            })
                            break  # Only take first matching rst7 in same dir
        
        return matching_pairs
    
    def _select_from_found_pairs(self, matching_pairs, search_dir):
        """Display found pairs and let user select."""
        from pathlib import Path
        
        # Display found pairs
        self.console.print(f"\n[green]Found {len(matching_pairs)} matching pair(s):[/green]")
        for i, pair in enumerate(matching_pairs, 1):
            relative_dir = pair['directory'].relative_to(search_dir) if pair['directory'] != search_dir else Path(".")
            self.console.print(f"  {i:2}. {pair['base_name']} in {relative_dir}")
            
        # Let user select which pairs to add
        self.console.print("\n[bold]Select pairs to add:[/bold]")
        self.console.print("  Enter numbers separated by commas (e.g., 1,3,5)")
        self.console.print("  Enter 'all' to add all pairs")
        self.console.print("  Enter 'cancel' to cancel")
        
        _pair_map = {str(i): pair['base_name'] for i, pair in enumerate(matching_pairs, 1)}
        _pair_map["all"] = "Add all pairs"
        _pair_map["cancel"] = "Cancel"
        selection = prompt_with_context(
            self.processor,
            "Selection",
            module="MD Manager - Structure Pairs",
            description="Select prmtop/rst7 pairs to add",
            options_map=_pair_map,
        ).strip().lower()
        
        if selection == 'cancel':
            return []
        elif selection == 'all':
            return matching_pairs
        else:
            # Parse comma-separated numbers
            selected_pairs = []
            try:
                indices = [int(x.strip()) - 1 for x in selection.split(',')]
                for idx in indices:
                    if 0 <= idx < len(matching_pairs):
                        selected_pairs.append(matching_pairs[idx])
                    else:
                        self.console.print(f"[yellow]Skipping invalid index: {idx + 1}[/yellow]")
                return selected_pairs
            except ValueError:
                self.console.print("[red]Invalid selection format[/red]")
                return []
    
    def _scan_workspace_for_structure_files(self):
        """Scan current workspace for prmtop and rst7 files (original format for single sim setup)."""
        from pathlib import Path
        
        current_dir = Path.cwd()
        workspace_files = {
            'prmtop': list(current_dir.glob("*.prmtop")),
            'rst7': list(current_dir.glob("*.rst7"))
        }
        
        return workspace_files
    
    def _scan_workspace_for_structure_pairs(self):
        """Scan current workspace and return structure file pairs for workflow setup."""
        from pathlib import Path
        
        current_dir = Path.cwd()
        prmtop_files = list(current_dir.glob("*.prmtop"))
        rst7_files = list(current_dir.glob("*.rst7"))
        rst7_files.extend(list(current_dir.glob("*.inpcrd")))  # Also include inpcrd files
        
        # Try to create sensible pairs based on filename matching
        pairs = []
        
        for prmtop in prmtop_files:
            prmtop_base = prmtop.stem
            
            # Look for matching rst7/inpcrd file
            best_match = None
            for rst7 in rst7_files:
                rst7_base = rst7.stem
                
                # Perfect match
                if prmtop_base == rst7_base:
                    best_match = rst7
                    break
                # Partial match (e.g., "system" matches "system_solv")
                elif prmtop_base in rst7_base or rst7_base in prmtop_base:
                    if not best_match:  # Take first partial match if no perfect match
                        best_match = rst7
            
            if best_match:
                pairs.append((str(prmtop), str(best_match)))
                rst7_files.remove(best_match)  # Remove to avoid duplicate pairing
            else:
                # No match found, pair with first available rst7 if any
                if rst7_files:
                    pairs.append((str(prmtop), str(rst7_files.pop(0))))
                    
        return pairs
        
    def _select_structure_pair(self):
        """Interactive selection of prmtop/rst7 file pair."""
        self.console.print(f"\n[bold]Select Structure File Pair[/bold]")
        
        # Select prmtop file or find pairs
        prmtop_result = self._browse_for_structure_file("prmtop")
        if not prmtop_result:
            return None
            
        # Check if user selected "find pairs"
        if isinstance(prmtop_result, tuple) and prmtop_result[0] == "FIND_PAIRS":
            # Return the found pairs directly
            return prmtop_result[1]
            
        prmtop_file = prmtop_result
            
        # After selecting prmtop from find, automatically search for rst7 in same location
        prmtop_dir = prmtop_file.parent
        
        # Automatically search for rst7 files in the prmtop directory
        self.console.print(f"\n[bold]Now searching for RST7 files in the same location...[/bold]")
        rst7_result = self._find_structure_files(prmtop_dir, "rst7", original_start_dir=prmtop_dir)
        
        if rst7_result == "BROWSE":
            # User wants to browse for rst7 elsewhere
            rst7_file = self._browse_for_structure_file("rst7", start_dir=prmtop_dir)
        else:
            rst7_file = rst7_result
            
        if not rst7_file:
            return None
            
        # Use prmtop stem as base name; only add rst7 stem if different
        if prmtop_file.stem == rst7_file.stem:
            name = prmtop_file.stem
        else:
            name = f"{prmtop_file.stem}_{rst7_file.stem}"

        return {
            'prmtop': prmtop_file,
            'rst7': rst7_file,
            'name': name
        }
        
    def _browse_for_structure_file(self, file_type, start_dir=None):
        """Browse for prmtop or rst7 files (similar to .mdin browser).

        Thin wrapper over the shared file browser: unified bare-N / q UX,
        filename-based session replay, and the historical `find` / `find pairs`
        commands preserved via extra_commands. Returns a Path, None on cancel,
        or the ("FIND_PAIRS", pairs) tuple from `find pairs`.
        """
        from pathlib import Path
        from datetime import datetime
        from proprep.utils.file_browser import file_browser

        start = Path(start_dir) if start_dir else self.working_directory

        def _ctime_detail(p):
            try:
                return "created " + datetime.fromtimestamp(os.path.getctime(p)).strftime("%m/%d/%Y")
            except OSError:
                return ""

        extra = {
            "find": (f"Search recursively for .{file_type} files",
                     lambda cur: self._find_structure_files(Path(cur), file_type, start_dir)),
        }
        if file_type == "prmtop":
            def _find_pairs(cur):
                pairs = self._find_matching_pairs_in_dir(Path(cur))
                return ("FIND_PAIRS", pairs) if pairs else None
            extra["find pairs"] = ("Find matching prmtop/rst7 pairs", _find_pairs)

        return file_browser(
            directory=str(start),
            extensions=[f".{file_type}"],
            console=self.console,
            processor=self.processor,
            label=f"{file_type} file",
            entry_detail=_ctime_detail,
            path_factory=Path,
            extra_commands=extra,
            module="MD Manager - Structure Browser",
        )



    def _find_structure_files(self, current_dir, file_type, original_start_dir=None):
        """Search recursively for prmtop or rst7 files."""
        from pathlib import Path
        
        self.console.print(f"\n[bold cyan]Searching for .{file_type} files in {current_dir}...[/bold cyan]")
        
        # Find all files of the target type recursively
        target_files = list(current_dir.rglob(f"*.{file_type}"))
        
        if not target_files:
            self.console.print(f"[grey50]No .{file_type} files found[/grey50]")
            return None
            
        self.console.print(f"[green]Found {len(target_files)} .{file_type} files:[/green]")
        
        # Display found files
        for i, file_path in enumerate(target_files, 1):
            relative_path = file_path.relative_to(current_dir)
            creation_time = datetime.fromtimestamp(file_path.stat().st_ctime)
            date_str = creation_time.strftime("%m/%d/%Y")
            self.console.print(f"  {i:2}. {relative_path} (created {date_str})")
            
        # Let user select
        while True:
            # For prmtop, just show select/cancel
            # For rst7, also show option to browse elsewhere
            if file_type == "rst7" and original_start_dir:
                prompt = f"\nSelect file (1-{len(target_files)}), 'browse' to return to browser, or 'cancel': "
            else:
                prompt = f"\nSelect file (1-{len(target_files)}) or 'cancel': "

            choice = prompt_with_context(
                self.processor, prompt.strip().rstrip(':').strip(),
                default="cancel", module="MD Manager - File Search",
                description="Select from recursively-found structure files",
            ).strip().lower()
            
            if choice == 'cancel':
                return None
            elif choice == 'browse' and file_type == "rst7":
                # Return special marker to go back to browser
                return "BROWSE"
                
            try:
                file_num = int(choice)
                if 1 <= file_num <= len(target_files):
                    selected_file = target_files[file_num - 1]
                    self.console.print(f"[green]Selected: {selected_file}[/green]")
                    return selected_file
                else:
                    self.console.print(f"[red]Please enter a number between 1 and {len(target_files)}[/red]")
            except ValueError:
                if file_type == "rst7" and original_start_dir:
                    self.console.print("[red]Please enter a valid number, 'browse', or 'cancel'[/red]")
                else:
                    self.console.print("[red]Please enter a valid number or 'cancel'[/red]")

    def _remove_structure_pair(self, structure_pairs):
        """Remove a structure pair from the list."""
        if not structure_pairs:
            return
            
        self.console.print(f"\n[bold]Remove Structure Pair[/bold]")
        for i, pair in enumerate(structure_pairs, 1):
            self.console.print(f"  {i}. {pair['name']}")

        try:
            # Build options map
            options_map = {}
            for i, pair in enumerate(structure_pairs, 1):
                options_map[str(i)] = pair['name']

            choice_str = prompt_with_context(
                self.processor,
                "Select pair to remove (number)",
                module="MD Manager - Structure Management",
                description="Select structure pair to remove"
            )
            choice = int(choice_str)
            if 1 <= choice <= len(structure_pairs):
                removed = structure_pairs.pop(choice - 1)
                self.console.print(f"[green]✓ Removed pair: {removed['name']}[/green]")
            else:
                self.console.print("[red]Invalid selection[/red]")
        except ValueError:
            self.console.print("[red]Please enter a valid number[/red]")
    
    def _select_structure_for_template(self, structure_pairs, template_id):
        """Select which structure pair to assign to a template."""
        templates = self.user_data_manager.list_templates()
        template_name = templates.get(template_id, {}).get('name', template_id)
        
        self.console.print(f"\n[bold]Assign template '{template_name}' to structure pair:[/bold]")
        for i, pair in enumerate(structure_pairs, 1):
            self.console.print(f"  {i}. {pair['name']}")
        self.console.print("  0. Skip assignment")

        try:
            # Build options map
            options_map = {}
            for i, pair in enumerate(structure_pairs, 1):
                options_map[str(i)] = pair['name']
            options_map["0"] = "Skip assignment"

            choice_str = prompt_with_context(
                self.processor,
                "Select structure pair",
                default="1",
                module="MD Manager - Structure Management",
                description="Assign template to structure pair",
                options_map=options_map
            )
            choice = int(choice_str)
            if choice == 0:
                return None
            elif 1 <= choice <= len(structure_pairs):
                return structure_pairs[choice - 1]
            else:
                self.console.print("[red]Invalid selection[/red]")
                return None
        except ValueError:
            self.console.print("[red]Please enter a valid number[/red]")
            return None
    
    def _remove_template_assignment(self, template_assignments):
        """Remove a template-structure assignment."""
        if not template_assignments:
            return
            
        self.console.print(f"\n[bold]Remove Template Assignment[/bold]")
        templates = self.user_data_manager.list_templates()
        for i, assignment in enumerate(template_assignments, 1):
            template_name = templates.get(assignment['template_id'], {}).get('name', assignment['template_id'])
            self.console.print(f"  {i}. {template_name} → {assignment['structure_pair']['name']}")

        try:
            # Build options map
            options_map = {}
            for i, assignment in enumerate(template_assignments, 1):
                template_name = templates.get(assignment['template_id'], {}).get('name', assignment['template_id'])
                options_map[str(i)] = f"{template_name} → {assignment['structure_pair']['name']}"

            choice_str = prompt_with_context(
                self.processor,
                "Select assignment to remove (number)",
                module="MD Manager - Structure Management",
                description="Select template assignment to remove",
                options_map=options_map
            )
            choice = int(choice_str)
            if 1 <= choice <= len(template_assignments):
                removed = template_assignments.pop(choice - 1)
                self.console.print(f"[green]✓ Removed assignment[/green]")
            else:
                self.console.print("[red]Invalid selection[/red]")
        except ValueError:
            self.console.print("[red]Please enter a valid number[/red]")
    
    def _preview_library_templates(self):
        """Preview templates from the library with categorized display."""
        templates = self.user_data_manager.list_templates()
        if not templates:
            self.console.print("[yellow]No templates available to preview[/yellow]")
            return
            
        # Get the controller for proper categorized display
        from .amber_controller import AmberController
        controller = AmberController(processor=self.processor)
        
        self.console.print(f"\n[bold]Select Template to Preview:[/bold]")
        
        # Use the same categorized display as template selection
        template_choices = controller._display_categorized_templates(templates, show_numbers=True)
        self.console.print(f"\n[grey50]0. Cancel preview[/grey50]")
        
        # Get valid choices
        valid_choices = list(template_choices.keys()) + ["0"]

        # Build options map for context
        context_options = {"0": "Cancel preview"}
        for key, template_id in template_choices.items():
            template_name = templates[template_id].get('name', template_id)
            context_options[key] = f"Preview: {template_name}"

        choice = prompt_with_context(
            self.processor,
            "Select template to preview",
            choices=valid_choices,
            default="0",
            module="MD Manager - Template Library",
            description="Select template to preview",
            options_map=context_options
        )
        
        if choice == "0":
            return
            
        # Get the selected template
        if choice in template_choices:
            template_id = template_choices[choice]
            template_metadata = templates[template_id]
            
            # Create a temporary config for preview (same pattern as existing _preview_template_content)
            temp_config = SimulationConfig(
                name=template_metadata.get('name', template_id),
                template_id=template_id,
                mdin_path=template_metadata.get('template_path', ''),
                engine="preview"
            )
            
            # Use the existing display method for consistency
            self._display_template_details(template_metadata, temp_config, controller)
        else:
            self.console.print("[red]Invalid selection[/red]")

    def _preview_assignment_templates(self, template_assignments):
        """Preview templates from assignments."""
        if not template_assignments:
            return
            
        self.console.print(f"\n[bold]Preview Template Content[/bold]")
        templates = self.user_data_manager.list_templates()
        for i, assignment in enumerate(template_assignments, 1):
            template_name = templates.get(assignment['template_id'], {}).get('name', assignment['template_id'])
            self.console.print(f"  {i}. {template_name}")

        try:
            # Build options map
            options_map = {}
            for i, assignment in enumerate(template_assignments, 1):
                template_name = templates.get(assignment['template_id'], {}).get('name', assignment['template_id'])
                options_map[str(i)] = template_name

            choice_str = prompt_with_context(
                self.processor,
                "Select template to preview (number)",
                module="MD Manager - Template Preview",
                description="Select template to preview content",
                options_map=options_map
            )
            choice = int(choice_str)
            if 1 <= choice <= len(template_assignments):
                template_id = template_assignments[choice - 1]['template_id']
                # Use existing preview method
                self._display_template_details(template_id)
            else:
                self.console.print("[red]Invalid selection[/red]")
        except ValueError:
            self.console.print("[red]Please enter a valid number[/red]")

    def _handle_remove_command(self, command):
        """Handle remove commands: remove N"""
        parts = command.split()
        if len(parts) != 2:
            self.console.print("[red]Usage: remove N[/red]")
            return
            
        try:
            index = int(parts[1]) - 1
        except ValueError:
            self.console.print("[red]Position must be an integer[/red]")
            return
            
        if not (0 <= index < len(self.simulation_queue)):
            self.console.print(f"[red]Position {parts[1]} out of range (1-{len(self.simulation_queue)})[/red]")
            return
            
        config = self.simulation_queue.queue[index]
        self.simulation_queue.remove_simulation(index)
        self.console.print(f"[green]✓ Removed '{config.name}'[/green]")

    def _interactive_file_import(self):
        """Enhanced interactive file browser for importing .mdin files."""
        current_dir = Path.cwd()
        recent_dirs = [Path.home(), Path.cwd()]  # Start with common directories
        bookmarks = {"home": Path.home(), "cwd": Path.cwd()}
        
        self.console.print("\n[bold cyan]Enhanced File Browser[/bold cyan]")
        self.console.print("[grey50]Commands: number, path, 'ls', 'find pattern', 'cd path', 'bookmark name', 'recent', 'help', 'cancel'[/grey50]")
        
        while True:
            self.console.print(f"\n[bold]Current:[/bold] {current_dir}")
            
            # Get directory contents
            dirs, mdin_files, other_files = self._scan_directory(current_dir)
            
            if dirs is None:  # Permission error
                self.console.print("[red]Permission denied[/red]")
                if recent_dirs:
                    current_dir = recent_dirs[-1]
                    continue
                else:
                    return
                    
            # Display directory contents
            items = []
            item_count = 0
            
            # Parent directory
            if current_dir.parent != current_dir:
                item_count += 1
                items.append((f"[{item_count:2}] [grey50].. (parent)[/grey50]", current_dir.parent, "parent"))
                
            # Subdirectories
            for d in sorted(dirs):
                item_count += 1
                items.append((f"[{item_count:2}] [blue]{d.name}/[/blue]", d, "dir"))
                
            # .mdin files
            for f in sorted(mdin_files):
                item_count += 1
                size = f.stat().st_size if f.exists() else 0
                size_str = self._format_file_size(size)
                items.append((f"[{item_count:2}] [green]{f.name}[/green] [grey50]({size_str})[/grey50]", f, "file"))
            
            # Display items
            for display_text, path, item_type in items:
                self.console.print(f"  {display_text}")
                
            if not items:
                self.console.print("  [grey50](empty directory)[/grey50]")
                
            # Show other files count if any
            if other_files:
                self.console.print(f"  [grey50]({len(other_files)} other files not shown)[/grey50]")
                
            # Get user input
            choice = prompt_with_context(
                self.processor,
                "\n[bold]>[/bold]",
                default="cancel",
                module="MD Manager - File Browser",
                description="File browser command"
            ).strip()
            
            if choice.lower() in ["cancel", "exit", "quit"]:
                return
                
            elif choice.lower() == "help":
                self._show_browser_help()
                continue
                
            elif choice.lower().startswith("ls"):
                # Handle ls with optional pattern: ls, ls *.mdin, ls */*.txt, etc.
                parts = choice.split(maxsplit=1)
                pattern = parts[1] if len(parts) > 1 else "*"
                self._list_files(current_dir, pattern)
                continue
                
            elif choice.lower().startswith("find "):
                # Enhanced find: find pattern, find *.mdin, find **/equilibration*
                pattern = choice[5:].strip()
                if self._enhanced_find(current_dir, pattern):
                    return  # File was imported
                continue
                
            elif choice.lower().startswith("cd "):
                new_path = choice[3:].strip()
                result_dir = self._change_directory(current_dir, new_path, bookmarks)
                if result_dir:
                    if result_dir not in recent_dirs:
                        recent_dirs.append(result_dir)
                        if len(recent_dirs) > 10:  # Keep last 10
                            recent_dirs = recent_dirs[-10:]
                    current_dir = result_dir
                continue
                
            elif choice.lower().startswith("bookmark "):
                name = choice[9:].strip()
                if name:
                    bookmarks[name] = current_dir
                    self.console.print(f"[green]Bookmarked '{current_dir}' as '{name}'[/green]")
                continue
                
            elif choice.lower() == "recent":
                self._show_recent_dirs(recent_dirs)
                continue
                
            # Try numeric choice
            try:
                choice_num = int(choice)
                if 1 <= choice_num <= len(items):
                    selected_item = items[choice_num - 1]
                    if selected_item[2] == "file":
                        self._import_mdin_file(selected_item[1])
                        return
                    else:
                        new_dir = selected_item[1]
                        if new_dir not in recent_dirs:
                            recent_dirs.append(new_dir)
                        current_dir = new_dir
                        continue
                else:
                    self.console.print(f"[red]Invalid selection (1-{len(items)})[/red]")
                    continue
            except ValueError:
                pass
                
            # Try as direct path
            result = self._handle_path_input(choice, current_dir, bookmarks)
            if result == "imported":
                return
            elif result:
                if result not in recent_dirs:
                    recent_dirs.append(result)
                current_dir = result
            else:
                self.console.print("[red]Invalid input. Type 'help' for commands.[/red]")
                
    def _scan_directory(self, directory):
        """Scan directory and categorize contents."""
        try:
            all_items = list(directory.iterdir())
            dirs = [item for item in all_items if item.is_dir()]
            mdin_files = [item for item in all_items if item.is_file() and item.suffix.lower() == ".mdin"]
            other_files = [item for item in all_items if item.is_file() and item.suffix.lower() != ".mdin"]
            return dirs, mdin_files, other_files
        except PermissionError:
            return None, None, None
            
    def _format_file_size(self, size_bytes):
        """Format file size in human readable format."""
        if size_bytes < 1024:
            return f"{size_bytes}B"
        elif size_bytes < 1024**2:
            return f"{size_bytes/1024:.1f}KB"
        elif size_bytes < 1024**3:
            return f"{size_bytes/1024**2:.1f}MB"
        else:
            return f"{size_bytes/1024**3:.1f}GB"
            
    def _show_browser_help(self):
        """Show file browser help."""
        self.console.print("\n[bold]File Browser Commands:[/bold]")
        self.console.print("  [cyan]number[/cyan]           - Select numbered item")
        self.console.print("  [cyan]/path/to/dir[/cyan]    - Navigate to absolute path")
        self.console.print("  [cyan]../subdir[/cyan]       - Navigate to relative path")
        self.console.print("  [cyan]~[/cyan]               - Go to home directory")
        self.console.print("  [cyan]ls[/cyan]              - List all files")
        self.console.print("  [cyan]ls *.mdin[/cyan]       - List files matching glob pattern")
        self.console.print("  [cyan]ls */*.mdin[/cyan]     - List files in subdirectories")
        self.console.print("  [cyan]find pattern[/cyan]    - Search recursively for pattern")
        self.console.print("  [cyan]find *.mdin[/cyan]     - Search using glob patterns")
        self.console.print("  [cyan]find **/equil*[/cyan]  - Deep search with patterns")
        self.console.print("  [cyan]cd path[/cyan]         - Change to directory")
        self.console.print("  [cyan]bookmark name[/cyan]   - Bookmark current directory")
        self.console.print("  [cyan]recent[/cyan]          - Show recently visited directories")
        self.console.print("  [cyan]cancel[/cyan]          - Exit browser")
        
    def _list_files(self, current_dir, pattern):
        """List files with glob pattern support."""
        import glob
        from pathlib import Path
        
        try:
            self.console.print(f"\n[bold]Listing: {pattern} in {current_dir}[/bold]")
            
            # Handle glob patterns
            if "*" in pattern or "?" in pattern or "[" in pattern:
                # Use glob for pattern matching
                search_path = current_dir / pattern
                matches = glob.glob(str(search_path), recursive=True)
                
                if not matches:
                    self.console.print("  [grey50]No files match pattern[/grey50]")
                    return
                    
                # Separate files and directories
                files = []
                dirs = []
                
                for match_str in sorted(matches):
                    match_path = Path(match_str)
                    if match_path.is_file():
                        files.append(match_path)
                    elif match_path.is_dir():
                        dirs.append(match_path)
                
                # Display directories first
                for directory in dirs:
                    rel_path = directory.relative_to(current_dir) if directory.is_relative_to(current_dir) else directory
                    self.console.print(f"  [blue]{rel_path}/[/blue]")
                    
                # Display files with details
                for i, file_path in enumerate(files, 1):
                    rel_path = file_path.relative_to(current_dir) if file_path.is_relative_to(current_dir) else file_path
                    size = file_path.stat().st_size if file_path.exists() else 0
                    
                    # Highlight .mdin files and add selection numbers
                    if file_path.suffix.lower() == ".mdin":
                        self.console.print(f"  [{i:2}] [green]{rel_path}[/green] [grey50]({self._format_file_size(size)})[/grey50]")
                    else:
                        self.console.print(f"  [{i:2}] {rel_path} [grey50]({self._format_file_size(size)})[/grey50]")
                
                # Allow importing .mdin files from results
                if any(f.suffix.lower() == ".mdin" for f in files):
                    choice = prompt_with_context(
                        self.processor,
                        "Select .mdin file number to import, or press Enter to continue",
                        default="",
                        module="MD Manager - File Browser",
                        description="Import .mdin file from search results"
                    )
                    if choice.isdigit():
                        choice_num = int(choice)
                        if 1 <= choice_num <= len(files) and files[choice_num-1].suffix.lower() == ".mdin":
                            self._import_mdin_file(files[choice_num-1])
                            return True
            else:
                # Simple listing without patterns (fallback to original behavior)
                all_items = sorted(current_dir.iterdir())
                for item in all_items:
                    if item.is_dir():
                        self.console.print(f"  [blue]{item.name}/[/blue]")
                    else:
                        size = item.stat().st_size if item.exists() else 0
                        self.console.print(f"  {item.name} [grey50]({self._format_file_size(size)})[/grey50]")
                        
        except Exception as e:
            self.console.print(f"[red]Error listing files: {e}[/red]")
        
        return False
        
    def _enhanced_find(self, current_dir, pattern):
        """Enhanced find with glob pattern support."""
        import glob
        from pathlib import Path
        
        try:
            # If pattern contains glob characters, use glob search
            if "*" in pattern or "?" in pattern or "[" in pattern:
                self.console.print(f"\n[bold]Finding: {pattern} (recursive glob search)[/bold]")
                
                # Use ** for recursive search if not already present
                if "**" not in pattern and "/" not in pattern:
                    search_pattern = f"**/*{pattern}*" if not pattern.startswith("*") else f"**/{pattern}"
                else:
                    search_pattern = pattern
                    
                search_path = current_dir / search_pattern
                matches = glob.glob(str(search_path), recursive=True)
                
                # Filter for .mdin files only
                mdin_matches = []
                for match_str in matches:
                    match_path = Path(match_str)
                    if match_path.is_file() and match_path.suffix.lower() == ".mdin":
                        mdin_matches.append(match_path)
                        
            else:
                # Simple name-based search (original behavior)
                self.console.print(f"\n[bold]Searching for '*{pattern}*' in {current_dir} (recursive)[/bold]")
                mdin_matches = []
                for item in current_dir.rglob(f"*{pattern}*"):
                    if item.is_file() and item.suffix.lower() == ".mdin":
                        mdin_matches.append(item)
            
            if not mdin_matches:
                self.console.print("  [grey50]No .mdin files found[/grey50]")
                return False
                
            # Sort and display results
            sorted_matches = sorted(mdin_matches)
            for i, file_path in enumerate(sorted_matches[:30], 1):  # Show max 30
                try:
                    rel_path = file_path.relative_to(current_dir)
                    size = file_path.stat().st_size
                    self.console.print(f"  [{i:2}] [green]{rel_path}[/green] [grey50]({self._format_file_size(size)})[/grey50]")
                except ValueError:
                    # File is outside current directory
                    size = file_path.stat().st_size
                    self.console.print(f"  [{i:2}] [green]{file_path}[/green] [grey50]({self._format_file_size(size)})[/grey50]")
                    
            if len(sorted_matches) > 30:
                self.console.print(f"  [grey50]... and {len(sorted_matches)-30} more matches[/grey50]")
                
            # Allow selection from results
            choice = prompt_with_context(
                self.processor,
                "Select file number to import, or press Enter to continue",
                default="",
                module="MD Manager - File Browser",
                description="Import file from search results (list command)"
            )
            if choice.isdigit():
                choice_num = int(choice)
                if 1 <= choice_num <= min(len(sorted_matches), 30):
                    self._import_mdin_file(sorted_matches[choice_num-1])
                    return True
                    
        except Exception as e:
            self.console.print(f"[red]Search error: {e}[/red]")
            
        return False
        
    def _find_files(self, current_dir, pattern):
        """Find files matching pattern."""
        self.console.print(f"\n[bold]Searching for '*{pattern}*' in {current_dir}:[/bold]")
        matches = []
        try:
            for item in current_dir.rglob(f"*{pattern}*"):
                if item.is_file() and item.suffix.lower() == ".mdin":
                    rel_path = item.relative_to(current_dir)
                    matches.append((str(rel_path), item))
                    
            if matches:
                for i, (rel_path, full_path) in enumerate(matches[:20], 1):  # Show max 20
                    size = full_path.stat().st_size
                    self.console.print(f"  [{i:2}] [green]{rel_path}[/green] [grey50]({self._format_file_size(size)})[/grey50]")
                if len(matches) > 20:
                    self.console.print(f"  [grey50]... and {len(matches)-20} more[/grey50]")
                    
                # Allow selection from search results
                choice = prompt_with_context(
                    self.processor,
                    "Select file number to import, or press Enter to continue",
                    default="",
                    module="MD Manager - File Browser",
                    description="Import file from search results (find command)"
                )
                if choice.isdigit():
                    choice_num = int(choice)
                    if 1 <= choice_num <= min(len(matches), 20):
                        self._import_mdin_file(matches[choice_num-1][1])
                        return True
            else:
                self.console.print("  [grey50]No .mdin files found matching pattern[/grey50]")
        except Exception as e:
            self.console.print(f"[red]Search error: {e}[/red]")
        return False
        
    def _change_directory(self, current_dir, path_str, bookmarks):
        """Change directory with various path formats."""
        try:
            # Check bookmarks first
            if path_str in bookmarks:
                return bookmarks[path_str]
                
            # Handle path
            if path_str.startswith("~"):
                new_path = Path(path_str).expanduser().resolve()
            elif path_str.startswith("/"):
                new_path = Path(path_str).resolve()
            else:
                new_path = (current_dir / path_str).resolve()
                
            if new_path.is_dir():
                return new_path
            else:
                self.console.print(f"[red]Directory not found: {path_str}[/red]")
                return None
        except Exception as e:
            self.console.print(f"[red]Invalid path: {e}[/red]")
            return None
            
    def _show_recent_dirs(self, recent_dirs):
        """Show recently visited directories."""
        if len(recent_dirs) <= 1:
            self.console.print("[grey50]No recent directories[/grey50]")
            return
            
        self.console.print("\n[bold]Recent directories:[/bold]")
        for i, directory in enumerate(recent_dirs[-10:], 1):
            self.console.print(f"  [{i:2}] {directory}")

        choice = prompt_with_context(
            self.processor,
            "Select directory number, or press Enter to continue",
            default="",
            module="MD Manager - File Browser",
            description="Select from recent directories"
        )
        if choice.isdigit():
            choice_num = int(choice)
            if 1 <= choice_num <= len(recent_dirs[-10:]):
                return recent_dirs[-10:][choice_num-1]
        return None
        
    def _handle_path_input(self, input_str, current_dir, bookmarks):
        """Handle various path input formats."""
        try:
            # Check if it's a bookmark
            if input_str in bookmarks:
                return bookmarks[input_str]
                
            # Try to resolve as path
            if input_str.startswith("~"):
                path = Path(input_str).expanduser().resolve()
            elif input_str.startswith("/"):
                path = Path(input_str).resolve()
            else:
                path = (current_dir / input_str).resolve()
                
            if path.is_file() and path.suffix.lower() == ".mdin":
                self._import_mdin_file(path)
                return "imported"
            elif path.is_dir():
                return path
            else:
                return None
        except Exception:
            return None

    def _direct_file_import(self, command):
        """Handle direct file import: import /path/file.mdin"""
        parts = command.split(maxsplit=1)
        if len(parts) != 2:
            self.console.print("[red]Usage: import /path/file.mdin[/red]")
            return
            
        file_path = Path(parts[1]).expanduser().resolve()
        self._import_mdin_file(file_path)

    def _import_mdin_file(self, file_path):
        """Import a specific .mdin file."""
        if not file_path.exists():
            self.console.print(f"[red]File not found: {file_path}[/red]")
            return
            
        if file_path.suffix.lower() != ".mdin":
            self.console.print("[red]File must have .mdin extension[/red]")
            return

        sim_name = prompt_with_context(
            self.processor,
            "Simulation name",
            default=file_path.stem,
            module="MD Manager - Import",
            description="Enter name for imported simulation"
        )
        
        config = SimulationConfig(
            name=sim_name,
            template_id=str(file_path),
            mdin_path=str(file_path),
            engine=""
        )
        
        self.simulation_queue.add_simulation(config)
        self.console.print(f"[green]✓ Imported '{sim_name}' from {file_path}[/green]")

    def _check_and_offer_cpmd(self, prmtop_file: str) -> Optional[Dict]:
        """Check workspace for constant pH config and offer to enable CpHMD.

        Returns a dict with cpin_file, cpin_config, and cpmd_settings if the
        user accepts, or None if no CPIN file exists or the user declines.
        """
        workspace = self.workspace
        if not workspace:
            return None

        cpin_config = workspace.get('cpin_config')
        cpin_file = workspace.get('cpin_file')
        if not cpin_config or not cpin_file:
            return None

        # Verify the CPIN file still exists on disk
        from pathlib import Path
        if not Path(cpin_file).exists():
            return None

        # Show what was configured
        sim_type = cpin_config.get('simulation_type', 'unknown')
        num_res = cpin_config.get('num_residues', 0)
        self.console.print(f"\n[bold cyan]Constant pH MD Available[/bold cyan]")
        self.console.print(f"  CPIN file: {Path(cpin_file).name}")
        self.console.print(f"  Solvent model: {sim_type}")
        self.console.print(f"  Titratable residues: {num_res}")

        # Check for modified prmtop (explicit solvent)
        modified_prmtop = cpin_config.get('modified_prmtop')
        if modified_prmtop:
            self.console.print(f"  Modified topology: {Path(modified_prmtop).name}")

        enable = confirm_with_context(
            self.processor,
            "Enable constant pH MD for production steps?",
            default=True,
            module="MD Manager - CpHMD",
            description="Enable constant pH molecular dynamics"
        )
        if not enable:
            return None

        # Determine icnstph value from simulation type
        icnstph = 2 if sim_type == 'explicit' else 1

        # Get target pH
        default_ph = 7.0
        ph_str = prompt_with_context(
            self.processor,
            "Target pH",
            default=str(default_ph),
            module="MD Manager - CpHMD",
            description="Solvent pH for constant pH MD"
        )
        try:
            solvph = float(ph_str)
        except ValueError:
            solvph = default_ph

        # Get protonation state change frequency
        ntcnstph_str = prompt_with_context(
            self.processor,
            "Steps between protonation state attempts (ntcnstph)",
            default="100",
            module="MD Manager - CpHMD",
            description="MC protonation state change frequency"
        )
        try:
            ntcnstph = int(ntcnstph_str)
        except ValueError:
            ntcnstph = 100

        # Salt concentration (reference energies were derived with saltcon=0.1)
        saltcon_str = prompt_with_context(
            self.processor,
            "Salt concentration in M (saltcon)",
            default="0.1",
            module="MD Manager - CpHMD",
            description="Salt concentration for CpHMD reference energies"
        )
        try:
            saltcon = float(saltcon_str)
        except ValueError:
            saltcon = 0.1

        # Build settings dict
        cpmd_settings = {
            'icnstph': icnstph,
            'solvph': solvph,
            'ntcnstph': ntcnstph,
            'saltcon': saltcon,
        }

        # Implicit solvent: igb must match what was used to derive reference energies
        if icnstph == 1:
            igb = cpin_config.get('igb', 2)
            cpmd_settings['igb'] = igb

        # Explicit solvent needs ntrelax
        if icnstph == 2:
            ntrelax_str = prompt_with_context(
                self.processor,
                "Solvent relaxation steps after state change (ntrelax)",
                default="200",
                module="MD Manager - CpHMD",
                description="Solvent relaxation dynamics steps (explicit solvent only)"
            )
            try:
                cpmd_settings['ntrelax'] = int(ntrelax_str)
            except ValueError:
                cpmd_settings['ntrelax'] = 200

        self.console.print(f"\n[green]✓ Constant pH MD enabled[/green]")
        self.console.print(f"  icnstph={cpmd_settings['icnstph']}, solvph={solvph}, "
                         f"ntcnstph={ntcnstph}, saltcon={saltcon}")
        if 'igb' in cpmd_settings:
            self.console.print(f"  igb={cpmd_settings['igb']} (from CPIN generation)")
        if 'ntrelax' in cpmd_settings:
            self.console.print(f"  ntrelax={cpmd_settings['ntrelax']}")
        self.console.print(f"  [grey50]CpHMD flags will be applied to production steps only[/grey50]")

        return {
            'cpin_file': cpin_file,
            'cpin_config': cpin_config,
            'cpmd_settings': cpmd_settings,
        }

    def _perform_simulations(self) -> bool:
        """Execute simulations with context-aware smart redirect."""
        # Restore any running simulations from previous session
        self._restore_running_simulations()

        # Get MD readiness status from workspace (single source of truth)
        md_status = self._get_md_ready_status()

        # Check if we have any simulations configured
        has_simulations = len(self.simulation_queue) > 0
        
        if not has_simulations:
            # Context-aware redirect: No simulations configured
            self.console.print(f"\n[bold cyan]===== Setup Required =====[/bold cyan]")
            self.console.print(f"[yellow]No simulations are currently configured.[/yellow]")
            self.console.print(f"\nWould you like to set up simulations now?")
            self.console.print(f"1. Setup simulations for submission")
            self.console.print(f"2. ← Return to main menu")

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2"],
                default="2",
                module="MD Manager - Execution",
                description="HPC submission setup",
                options_map={
                    "1": "Setup simulations for submission",
                    "2": "← Return to main menu"
                }
            )
            if choice == "1":
                command = SetupSingleSimulationCommand(self.processor)
                return command.execute()
            else:
                return True  # Return to main menu
        
        # We have simulations, show the execution interface
        while True:
            self.console.print(f"\n[bold cyan]===== Execute Simulations =====[/bold cyan]")

            # Display workflow status if there are running simulations or pending workflow steps
            self._display_workflow_status()

            # Check if we can proceed with simulations using unified status
            if not md_status["ready"]:
                self.console.print(f"\n[bold red]Cannot proceed: Missing requirements[/bold red]")
                for missing_item in md_status["missing"]:
                    self.console.print(f"  • {missing_item}")
                return_str = prompt_with_context(
                    self.processor,
                    "Return to menu?",
                    choices=["y", "n"],
                    default="y",
                    module="MD Manager - Step 5: Execute",
                    description="Return to main menu due to missing requirements",
                    options_map={"y": "Yes, return to menu", "n": "No, stay here"}
                )
                if return_str.lower() == "y":
                    return True
            
            # Display available actions
            if len(self.simulation_queue) > 0:
                self.console.print(f"\n[bold]Execute Options:[/bold]")
                self.console.print(f"  1. Execute simulation queue ({len(self.simulation_queue)} simulations)")
                self.console.print(f"  2. Manage queue (reorder, hold, remove)")
            else:
                self.console.print(f"\n[bold]Execute Options:[/bold]")
                self.console.print("  1. [grey50]Execute simulation queue (none queued)[/grey50]")
                self.console.print("  2. [grey50]Manage queue (none queued)[/grey50]")
                
            self.console.print(f"\n[bold]Other Options:[/bold]")
            self.console.print("  3. Monitor running simulation") 
            self.console.print("  4. Analyze completed simulations")
            
            self.console.print(f"\n[bold]Navigation:[/bold]")
            self.console.print("  b. ← Back")
            
            valid_choices = ["3", "4", "b"]
            if len(self.simulation_queue) > 0:
                valid_choices = ["1", "2", "3", "4", "b"]

            # Build options map
            context_options = {
                "3": "Monitor running simulation",
                "4": "Analyze completed simulations",
                "b": "← Back"
            }
            if len(self.simulation_queue) > 0:
                context_options["1"] = f"Execute simulation queue ({len(self.simulation_queue)} simulations)"
                context_options["2"] = "Manage queue (reorder, hold, remove)"

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=valid_choices,
                default="b",
                module="MD Manager - Execution",
                description="Execution menu",
                options_map=context_options
            )
            
            if choice == "1" and len(self.simulation_queue) > 0:
                self._execute_simulation_queue()
            elif choice == "2" and len(self.simulation_queue) > 0:
                self._manage_simulation_queue()
            elif choice == "3":
                self._monitor_simulation()
            elif choice == "4":
                self._analyze_simulations()
            elif choice == "b":
                return True
    
    # _display_step_progress method removed - now using WorkflowLayoutFormatter.step_header()

    def _manage_simulation_queue(self):
        """Manage simulation queue: reorder, hold, remove simulations."""
        while True:
            self.console.print(f"\n[bold cyan]===== Queue Management =====[/bold cyan]")
            
            # Display current queue with status
            self.simulation_queue.display_queue(self.console)
            
            self.console.print(f"\n[bold]Management Options:[/bold]")
            self.console.print("  1. Manage protocols (reorder, remove entire protocols)")
            self.console.print("  2. Manage protocol steps (reorder, remove steps within protocols)")
            self.console.print("  3. Set simulation status (active/hold)")
            self.console.print("  4. Clear entire queue")
            
            self.console.print(f"\n[bold]Navigation:[/bold]")
            self.console.print("  b. ← Back to execution menu")

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3", "4", "b"],
                default="b",
                module="MD Manager - Queue Management",
                description="Queue management options",
                options_map={
                    "1": "Manage protocols (reorder, remove entire protocols)",
                    "2": "Manage protocol steps (reorder, remove steps within protocols)",
                    "3": "Set simulation status (active/hold)",
                    "4": "Clear entire queue",
                    "b": "← Back to execution menu"
                }
            )

            if choice == "1":
                self._manage_workflows()
            elif choice == "2":
                self._manage_workflow_steps()
            elif choice == "3":
                self._set_simulation_status()
            elif choice == "4":
                confirm_str = prompt_with_context(
                    self.processor,
                    "Clear entire simulation queue?",
                    choices=["y", "n"],
                    default="n",
                    module="MD Manager - Queue Management",
                    description="Confirm clear queue",
                    options_map={"y": "Yes, clear queue", "n": "No, cancel"}
                )
                if confirm_str.lower() == "y":
                    self.simulation_queue.clear()
                    self.console.print("[green]✓ Queue cleared[/green]")
            elif choice == "b":
                break
    
    def _manage_workflows(self):
        """Manage entire workflows: reorder, remove."""
        # Group simulations by workflow
        standalone_sims = []
        workflow_groups = {}
        
        for config in self.simulation_queue.queue:
            if config.workflow_id:
                if config.workflow_id not in workflow_groups:
                    workflow_groups[config.workflow_id] = []
                workflow_groups[config.workflow_id].append(config)
            else:
                standalone_sims.append(config)
        
        if not workflow_groups and not standalone_sims:
            self.console.print("[yellow]No protocols to manage[/yellow]")
            return
        
        while True:
            self.console.print(f"\n[bold cyan]Protocol Management[/bold cyan]")
            
            # Display workflows
            workflow_list = []
            display_index = 1
            
            # Show standalone simulations as a group
            if standalone_sims:
                self.console.print(f"\n[cyan]Standalone Simulations ({len(standalone_sims)} simulations):[/cyan]")
                for sim in standalone_sims:
                    self.console.print(f"    • {sim.name}")
                workflow_list.append(("standalone", standalone_sims))
                display_index += 1
            
            # Show workflows
            for wf_id, sims in workflow_groups.items():
                wf_name = wf_id
                if wf_id in self.simulation_queue._workflows:
                    wf = self.simulation_queue._workflows[wf_id]
                    wf_name = f"{wf.name}"
                else:
                    # Try to get a more readable name from workspace metadata
                    try:
                        md_metadata = self.workspace.get('md_metadata', {})
                        custom_workflows = md_metadata.get('custom_workflows', {})
                        for workflow_name, workflow_info in custom_workflows.items():
                            if workflow_info.get('id') == wf_id:
                                wf_name = workflow_name
                                break
                        else:
                            # Fallback: Use first simulation's structure name as prefix
                            if sims and hasattr(sims[0], 'name'):
                                structure_part = sims[0].name.split('_')[0]
                                wf_name = f"Protocol for {structure_part}"
                            else:
                                wf_name = f"Protocol {wf_id[:8]}..."
                    except:
                        wf_name = f"Protocol {wf_id[:8]}..."
                
                sorted_sims = sorted(sims, key=lambda x: x.workflow_step or 0)
                self.console.print(f"\n[cyan]{display_index}. Protocol: {wf_name} ({len(sorted_sims)} steps):[/cyan]")
                for sim in sorted_sims:
                    status_color = "green" if sim.status == "active" else "yellow"
                    status_display = f"[{status_color}]{sim.status.upper()}[/{status_color}]"
                    self.console.print(f"    Step {sim.workflow_step}: {sim.name} ({status_display})")
                
                workflow_list.append((wf_id, sorted_sims))
                display_index += 1
            
            self.console.print(f"\n[bold]Protocol Management Options:[/bold]")
            self.console.print("  1. Reorder protocols")
            self.console.print("  2. Remove entire protocol")
            self.console.print("  b. Back to queue management")

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "b"],
                default="b",
                module="MD Manager - Protocol Management",
                description="Protocol management actions",
                options_map={
                    "1": "Reorder protocols",
                    "2": "Remove entire protocol",
                    "b": "Back to queue management"
                }
            )
            
            if choice == "b":
                break
            elif choice == "1":
                self._reorder_workflows(workflow_list)
            elif choice == "2":
                self._remove_workflows(workflow_list)
    
    def _reorder_workflows(self, workflow_list):
        """Reorder entire workflows."""
        if len(workflow_list) < 2:
            self.console.print("[yellow]Need at least 2 protocols to reorder[/yellow]")
            return
        
        self.console.print("\n[bold]Current Protocol Order:[/bold]")
        for i, (wf_id, sims) in enumerate(workflow_list, 1):
            if wf_id == "standalone":
                self.console.print(f"  {i}. Standalone Simulations")
            else:
                wf_name = wf_id
                if wf_id in self.simulation_queue._workflows:
                    wf_name = self.simulation_queue._workflows[wf_id].name
                self.console.print(f"  {i}. {wf_name}")
        
        self.console.print("\n[bold]Reorder Instructions:[/bold]")
        self.console.print("Enter new order as comma-separated numbers (e.g., '2,1,3' to move first protocol to second position)")

        try:
            new_order = prompt_with_context(
                self.processor,
                "New order",
                module="MD Manager - Queue Management",
                description="Enter new protocol order (comma-separated)"
            ).strip()
            indices = [int(x.strip()) - 1 for x in new_order.split(',')]
            
            if len(indices) != len(workflow_list):
                self.console.print(f"[red]Must specify {len(workflow_list)} positions[/red]")
                return
                
            if set(indices) != set(range(len(workflow_list))):
                self.console.print("[red]Invalid order - must use each position exactly once[/red]")
                return
            
            # Reorder workflows
            reordered_workflows = [workflow_list[i] for i in indices]
            
            # Rebuild the entire queue in new order
            new_queue = []
            for wf_id, sims in reordered_workflows:
                new_queue.extend(sims)
            
            # Replace queue contents
            self.simulation_queue.clear()
            for config in new_queue:
                self.simulation_queue.add_simulation(config)
            
            self.console.print("[green]✓ Protocols reordered successfully[/green]")
            
        except (ValueError, IndexError):
            self.console.print("[red]Invalid input format[/red]")
    
    def _remove_workflows(self, workflow_list):
        """Remove entire workflows."""
        self.console.print("\n[bold]Select protocols to remove:[/bold]")
        for i, (wf_id, sims) in enumerate(workflow_list, 1):
            if wf_id == "standalone":
                self.console.print(f"  {i}. Standalone Simulations ({len(sims)} simulations)")
            else:
                wf_name = wf_id
                if wf_id in self.simulation_queue._workflows:
                    wf_name = self.simulation_queue._workflows[wf_id].name
                self.console.print(f"  {i}. {wf_name} ({len(sims)} steps)")
        
        try:
            selections = prompt_with_context(
                self.processor,
                "Enter protocol numbers to remove (comma-separated)",
                module="MD Manager - Queue Management",
                description="Select protocols to remove"
            ).strip()
            if not selections:
                return
                
            indices = [int(x.strip()) - 1 for x in selections.split(',')]
            valid_indices = [i for i in indices if 0 <= i < len(workflow_list)]
            
            if not valid_indices:
                self.console.print("[red]No valid protocol numbers provided[/red]")
                return
            
            # Get workflows to remove
            workflows_to_remove = [workflow_list[i] for i in valid_indices]
            
            # Confirm removal
            total_sims = sum(len(sims) for _, sims in workflows_to_remove)
            confirm_str = prompt_with_context(
                self.processor,
                f"Remove {len(workflows_to_remove)} protocol(s) containing {total_sims} simulations?",
                choices=["y", "n"],
                default="n",
                module="MD Manager - Queue Management",
                description="Confirm protocol removal from queue",
                options_map={"y": "Yes, remove protocols", "n": "No, cancel"}
            )
            if not (confirm_str.lower() == "y"):
                return
            
            # Remove simulations
            sims_to_remove = []
            for wf_id, sims in workflows_to_remove:
                sims_to_remove.extend(sims)
                if wf_id == "standalone":
                    self.console.print(f"[yellow]Removed: Standalone Simulations ({len(sims)} simulations)[/yellow]")
                else:
                    wf_name = wf_id
                    if wf_id in self.simulation_queue._workflows:
                        wf_name = self.simulation_queue._workflows[wf_id].name
                    self.console.print(f"[yellow]Removed: {wf_name} ({len(sims)} steps)[/yellow]")
            
            # Remove from queue
            for sim_to_remove in sims_to_remove:
                for i, config in enumerate(list(self.simulation_queue.queue)):
                    if config.name == sim_to_remove.name and config.template_id == sim_to_remove.template_id:
                        self.simulation_queue.remove_simulation(i)
                        break
            
            self.console.print(f"[green]✓ Removed {len(workflows_to_remove)} protocol(s)[/green]")
            
        except (ValueError, IndexError):
            self.console.print("[red]Invalid input format[/red]")
    
    def _manage_workflow_steps(self):
        """Manage steps within workflows: reorder, remove steps."""
        # Group simulations by workflow
        workflow_groups = {}
        
        for config in self.simulation_queue.queue:
            if config.workflow_id:
                if config.workflow_id not in workflow_groups:
                    workflow_groups[config.workflow_id] = []
                workflow_groups[config.workflow_id].append(config)
        
        if not workflow_groups:
            self.console.print("[yellow]No protocols found. Use standalone simulation management instead.[/yellow]")
            return
        
        # Let user select which workflow to manage
        workflow_list = list(workflow_groups.items())
        
        self.console.print(f"\n[bold cyan]Select Protocol to Manage:[/bold cyan]")
        for i, (wf_id, sims) in enumerate(workflow_list, 1):
            wf_name = wf_id
            if wf_id in self.simulation_queue._workflows:
                wf = self.simulation_queue._workflows[wf_id]
                wf_name = f"{wf.name}"
            else:
                # Try to get a more readable name from workspace metadata
                try:
                    md_metadata = self.workspace.get('md_metadata', {})
                    custom_workflows = md_metadata.get('custom_workflows', {})
                    for workflow_name, workflow_info in custom_workflows.items():
                        if workflow_info.get('id') == wf_id:
                            wf_name = workflow_name
                            break
                    else:
                        # Fallback: Use first simulation's structure name as prefix
                        if sims and hasattr(sims[0], 'name'):
                            structure_part = sims[0].name.split('_')[0]
                            wf_name = f"Protocol for {structure_part}"
                        else:
                            wf_name = f"Protocol {wf_id[:8]}..."
                except:
                    wf_name = f"Protocol {wf_id[:8]}..."
            
            sorted_sims = sorted(sims, key=lambda x: x.workflow_step or 0)
            self.console.print(f"  {i}. {wf_name} ({len(sorted_sims)} steps)")
        
        try:
            while True:
                choice_str = prompt_with_context(
                    self.processor,
                    "Select protocol number",
                    default="1",
                    module="MD Manager - Queue Management",
                    description="Select protocol to manage"
                )
                try:
                    choice = int(choice_str)
                    if 1 <= choice <= len(workflow_list):
                        break
                    else:
                        self.console.print(f"[red]Please enter a number between 1 and {len(workflow_list)}[/red]")
                except ValueError:
                    self.console.print("[red]Please enter a valid number[/red]")
            
            selected_wf_id, selected_sims = workflow_list[choice - 1]
            self._manage_single_workflow_steps(selected_wf_id, selected_sims)
            
        except KeyboardInterrupt:
            self.console.print("[yellow]Cancelled[/yellow]")
    
    def _manage_single_workflow_steps(self, workflow_id: str, workflow_sims: list):
        """Manage steps within a specific workflow."""
        while True:
            # Sort steps by workflow_step number
            sorted_sims = sorted(workflow_sims, key=lambda x: x.workflow_step or 0)
            
            wf_name = workflow_id
            if workflow_id in self.simulation_queue._workflows:
                wf = self.simulation_queue._workflows[workflow_id]
                wf_name = f"{wf.name}"
            else:
                # Try to get a more readable name from workspace metadata
                try:
                    md_metadata = self.workspace.get('md_metadata', {})
                    custom_workflows = md_metadata.get('custom_workflows', {})
                    for workflow_name, workflow_info in custom_workflows.items():
                        if workflow_info.get('id') == workflow_id:
                            wf_name = workflow_name
                            break
                    else:
                        # Fallback: Use first simulation's structure name as prefix
                        if workflow_sims and hasattr(workflow_sims[0], 'name'):
                            structure_part = workflow_sims[0].name.split('_')[0]
                            wf_name = f"Protocol for {structure_part}"
                        else:
                            wf_name = f"Protocol {workflow_id[:8]}..."
                except:
                    wf_name = f"Protocol {workflow_id[:8]}..."
            
            self.console.print(f"\n[bold cyan]Managing Steps in: {wf_name}[/bold cyan]")
            
            # Display workflow steps
            for i, sim in enumerate(sorted_sims, 1):
                status_color = "green" if sim.status == "active" else "yellow"
                status_display = f"[{status_color}]{sim.status.upper()}[/{status_color}]"
                self.console.print(f"  Step {sim.workflow_step}: {sim.name} ({status_display})")
            
            self.console.print(f"\n[bold]Step Management Options:[/bold]")
            self.console.print("  1. Reorder steps")
            self.console.print("  2. Remove steps")

            self.console.print(f"\n[bold]Navigation:[/bold]")
            self.console.print("  b. ← Back to protocol selection")

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "b"],
                default="b",
                module="MD Manager - Step Management",
                description="Protocol step management actions",
                options_map={
                    "1": "Reorder steps",
                    "2": "Remove steps",
                    "b": "← Back to protocol selection"
                }
            )
            
            if choice == "b":
                break
            elif choice == "1":
                self._reorder_workflow_steps(workflow_id, sorted_sims)
                # Refresh the sims list after reordering
                workflow_sims[:] = [s for s in self.simulation_queue.queue if s.workflow_id == workflow_id]
            elif choice == "2":
                self._remove_workflow_steps(workflow_id, sorted_sims)
                # Refresh the sims list after removal
                workflow_sims[:] = [s for s in self.simulation_queue.queue if s.workflow_id == workflow_id]
                if not workflow_sims:  # If all steps removed, exit
                    self.console.print("[yellow]All steps removed from protocol[/yellow]")
                    break
    
    def _reorder_workflow_steps(self, workflow_id: str, workflow_sims: list):
        """Reorder steps within a workflow."""
        if len(workflow_sims) < 2:
            self.console.print("[yellow]Need at least 2 steps to reorder[/yellow]")
            return
        
        self.console.print("\n[bold]Current Step Order:[/bold]")
        for i, sim in enumerate(workflow_sims, 1):
            self.console.print(f"  {i}. Step {sim.workflow_step}: {sim.name}")
        
        self.console.print("\n[bold]Reorder Instructions:[/bold]")
        self.console.print("Enter new order as comma-separated numbers (e.g., '1,3,2' to swap steps 2 and 3)")

        try:
            new_order = prompt_with_context(
                self.processor,
                "New order",
                module="MD Manager - Queue Management",
                description="Enter new step order (comma-separated)"
            ).strip()
            indices = [int(x.strip()) - 1 for x in new_order.split(',')]
            
            if len(indices) != len(workflow_sims):
                self.console.print(f"[red]Must specify {len(workflow_sims)} positions[/red]")
                return
                
            if set(indices) != set(range(len(workflow_sims))):
                self.console.print("[red]Invalid order - must use each position exactly once[/red]")
                return
            
            # Reorder steps and update step numbers
            reordered_sims = [workflow_sims[i] for i in indices]
            
            # Update workflow_step numbers to match new order
            for step_num, sim in enumerate(reordered_sims, 1):
                sim.workflow_step = step_num
            
            # Update the original list
            workflow_sims[:] = reordered_sims
            
            # Sync to workspace
            self.simulation_queue._sync_to_workspace()
            
            self.console.print("[green]✓ Protocol steps reordered successfully[/green]")
            
        except (ValueError, IndexError):
            self.console.print("[red]Invalid input format[/red]")
    
    def _remove_workflow_steps(self, workflow_id: str, workflow_sims: list):
        """Remove steps from a workflow."""
        self.console.print("\n[bold]Select steps to remove:[/bold]")
        for i, sim in enumerate(workflow_sims, 1):
            status_color = "green" if sim.status == "active" else "yellow"
            status_display = f"[{status_color}]{sim.status.upper()}[/{status_color}]"
            self.console.print(f"  {i}. Step {sim.workflow_step}: {sim.name} ({status_display})")
        
        try:
            selections = prompt_with_context(
                self.processor,
                "Enter step numbers to remove (comma-separated)",
                module="MD Manager - Queue Management",
                description="Select steps to remove from protocol"
            ).strip()
            if not selections:
                return
                
            indices = [int(x.strip()) - 1 for x in selections.split(',')]
            valid_indices = [i for i in indices if 0 <= i < len(workflow_sims)]
            
            if not valid_indices:
                self.console.print("[red]No valid step numbers provided[/red]")
                return
            
            # Get steps to remove
            steps_to_remove = [workflow_sims[i] for i in valid_indices]

            # Confirm removal
            confirm_str = prompt_with_context(
                self.processor,
                f"Remove {len(steps_to_remove)} step(s) from protocol?",
                choices=["y", "n"],
                default="n",
                module="MD Manager - Queue Management",
                description="Confirm step removal from protocol",
                options_map={"y": "Yes, remove steps", "n": "No, cancel"}
            )
            if not (confirm_str.lower() == "y"):
                return
            
            # Remove from queue
            for step_to_remove in steps_to_remove:
                for i, config in enumerate(list(self.simulation_queue.queue)):
                    if (config.name == step_to_remove.name and 
                        config.template_id == step_to_remove.template_id and
                        config.workflow_id == workflow_id):
                        self.simulation_queue.remove_simulation(i)
                        self.console.print(f"[yellow]Removed: Step {step_to_remove.workflow_step}: {step_to_remove.name}[/yellow]")
                        break
            
            # Renumber remaining steps
            remaining_sims = [s for s in self.simulation_queue.queue if s.workflow_id == workflow_id]
            remaining_sims.sort(key=lambda x: x.workflow_step or 0)
            for step_num, sim in enumerate(remaining_sims, 1):
                sim.workflow_step = step_num
            
            self.simulation_queue._sync_to_workspace()
            self.console.print(f"[green]✓ Removed {len(steps_to_remove)} step(s) and renumbered remaining steps[/green]")
            
        except (ValueError, IndexError):
            self.console.print("[red]Invalid input format[/red]")
                
    def _toggle_simulation_status(self):
        """Set ACTIVE/HOLD status for simulations using concise syntax."""
        queue = self.simulation_queue.queue
        if not queue:
            self.console.print("[yellow]Queue is empty[/yellow]")
            return

        self.console.print("\n[bold]Simulation Status:[/bold]")
        for i, sim in enumerate(queue, 1):
            status_color = "green" if sim.status == "active" else "yellow"
            self.console.print(f"  {i:>2}. {sim.name}  [{status_color}]{sim.status.upper()}[/{status_color}]")

        self.console.print()
        self.console.print("[bold]Commands:[/bold]")
        self.console.print("  hold \\[n]                 Put on hold (e.g., hold 3)")
        self.console.print("  active \\[n]               Set active (e.g., active 1)")
        self.console.print("  hold \\[n],\\[n]             Multiple (e.g., hold 1,2,5)")
        self.console.print("  active \\[n]-\\[n]           Range (e.g., active 1-4)")
        self.console.print("  hold all                 Hold all simulations")
        self.console.print("  active all               Activate all simulations")
        self.console.print("  b                        Back")

        command = prompt_with_context(
            self.processor,
            "Command",
            default="b",
            module="MD Manager - Step 3 Queue",
            description="Set simulation status (e.g. 'hold 1-3', 'active 2,5')"
        )

        command = command.strip().lower()
        if command == "b":
            return

        # Parse: <status> <range>
        parts = command.split(None, 1)
        if len(parts) != 2 or parts[0] not in ("hold", "active"):
            self.console.print("[red]Invalid syntax. Examples: hold 1-3, active 2,5, hold all[/red]")
            return

        new_status = parts[0]
        range_str = parts[1].strip()

        # Parse indices
        indices = self._parse_index_range(range_str, len(queue))
        if indices is None:
            return

        for idx in indices:
            queue[idx].status = new_status

        count = len(indices)
        status_color = "green" if new_status == "active" else "yellow"
        self.console.print(f"[{status_color}]{count} simulation(s) set to {new_status.upper()}[/{status_color}]")

    def _parse_index_range(self, range_str: str, queue_len: int):
        """Parse index expressions like '1,2,5', '1-4', 'all' into 0-based indices."""
        if range_str == "all":
            return list(range(queue_len))

        indices = set()
        for part in range_str.split(','):
            part = part.strip()
            if '-' in part:
                bounds = part.split('-', 1)
                try:
                    start, end = int(bounds[0]), int(bounds[1])
                    if start < 1 or end > queue_len or start > end:
                        self.console.print(f"[red]Range {part} out of bounds (1-{queue_len})[/red]")
                        return None
                    indices.update(range(start - 1, end))
                except ValueError:
                    self.console.print(f"[red]Invalid range: {part}[/red]")
                    return None
            else:
                try:
                    n = int(part)
                    if n < 1 or n > queue_len:
                        self.console.print(f"[red]{n} out of bounds (1-{queue_len})[/red]")
                        return None
                    indices.add(n - 1)
                except ValueError:
                    self.console.print(f"[red]Invalid number: {part}[/red]")
                    return None

        return sorted(indices)

    def _expand_reorder_input(self, text: str) -> list:
        """Expand a reorder input string supporting ranges and individual numbers.

        Examples:
            '6-19,1-5'   → [6,7,...,19,1,2,...,5]
            '1,3,2,4'    → [1,3,2,4]
            '6-19 1-5'   → [6,7,...,19,1,2,...,5]
        """
        indices = []
        for token in text.replace(',', ' ').split():
            token = token.strip()
            if '-' in token:
                start, end = token.split('-', 1)
                s, e = int(start.strip()), int(end.strip())
                if s <= e:
                    indices.extend(range(s, e + 1))
                else:
                    indices.extend(range(s, e - 1, -1))
            else:
                indices.append(int(token))
        return indices

    def _reorder_simulations(self):
        """Reorder simulations in queue."""
        if len(self.simulation_queue) < 2:
            self.console.print("[yellow]Need at least 2 simulations to reorder[/yellow]")
            return

        self.console.print("\n[bold]Current Order:[/bold]")
        for i, config in enumerate(self.simulation_queue.queue, 1):
            self.console.print(f"  {i}. {config.name}")

        self.console.print("\n[bold]Reorder Instructions:[/bold]")
        self.console.print("Enter new order using numbers and/or ranges")
        self.console.print("[grey50]  Examples: '1,3,2,4' or '6-19,1-5' or '6-19 1-5'[/grey50]")

        try:
            new_order = prompt_with_context(
                self.processor,
                "New order",
                module="MD Manager - Queue Management",
                description="Enter new simulation order (numbers/ranges)"
            ).strip()
            indices = [x - 1 for x in self._expand_reorder_input(new_order)]
            
            if len(indices) != len(self.simulation_queue):
                self.console.print(f"[red]Must specify {len(self.simulation_queue)} positions[/red]")
                return
                
            if set(indices) != set(range(len(self.simulation_queue))):
                self.console.print("[red]Invalid order - must use each position exactly once[/red]")
                return
                
            # Reorder the queue using proper SimulationQueue methods
            old_queue = list(self.simulation_queue.queue)  # Create a copy
            new_queue = [old_queue[i] for i in indices]
            
            # Clear current queue and rebuild with new order
            self.simulation_queue.clear()
            for config in new_queue:
                # Update workflow step to match new order
                config.workflow_step = len(self.simulation_queue._local_queue) + 1
                self.simulation_queue.add_simulation(config)
            
            self.console.print("[green]✓ Queue reordered successfully[/green]")
            
        except (ValueError, IndexError):
            self.console.print("[red]Invalid input format[/red]")
    
    def _set_simulation_status(self):
        """Set simulation status (active/hold)."""
        if not self.simulation_queue:
            self.console.print("[yellow]No simulations to manage[/yellow]")
            return
            
        while True:
            # Display current queue with status
            self.simulation_queue.display_queue(self.console)
            
            self.console.print("\n[bold]Status Management Options:[/bold]")
            self.console.print("  1. Set individual simulation status")
            self.console.print("  2. Set all simulations to active")
            self.console.print("  3. Set all simulations to hold")

            self.console.print(f"\n[bold]Navigation:[/bold]")
            self.console.print("  b. ← Back to queue management")

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1","2","3","b"],
                default="1",
                module="MD Manager - Status Management",
                description="Simulation status management",
                options_map={
                    "1": "Set individual simulation status",
                    "2": "Set all simulations to active",
                    "3": "Set all simulations to hold",
                    "b": "← Back to queue management"
                }
            )
            
            if choice == "b":
                break
            elif choice == "1":
                self._set_individual_status()
            elif choice == "2":
                self._set_all_status("active")
            elif choice == "3":
                self._set_all_status("hold")
    
    def _set_individual_status(self):
        """Set status for individual simulations."""
        self.console.print("\n[bold]Select simulations to change status:[/bold]")
        
        # Create numbered list for selection
        sim_list = list(self.simulation_queue.queue)
        for i, config in enumerate(sim_list, 1):
            status_color = "green" if config.status == "active" else "yellow"
            status_display = f"[{status_color}]{config.status.upper()}[/{status_color}]"
            self.console.print(f"  {i}. {config.name} - {status_display}")
        
        try:
            selections = prompt_with_context(
                self.processor,
                "Enter simulation numbers to change (comma-separated)",
                module="MD Manager - Queue Management",
                description="Select simulations to change status"
            ).strip()
            if not selections:
                return

            indices = [int(x.strip()) - 1 for x in selections.split(',')]
            valid_indices = [i for i in indices if 0 <= i < len(sim_list)]

            if not valid_indices:
                self.console.print("[red]No valid simulation numbers provided[/red]")
                return

            # Ask for new status
            new_status = prompt_with_context(
                self.processor,
                "Set status to",
                choices=["active", "hold"],
                default="active",
                module="MD Manager - Queue Management",
                description="Select new status",
                options_map={"active": "Active (will run)", "hold": "Hold (skip)"}
            )
            
            # Update status for selected simulations
            changed_count = 0
            for i in valid_indices:
                if sim_list[i].status != new_status:
                    sim_list[i].status = new_status
                    changed_count += 1
                    self.console.print(f"[green]✓ Set {sim_list[i].name} to {new_status.upper()}[/green]")
            
            if changed_count > 0:
                self.simulation_queue._sync_to_workspace()
                self.console.print(f"[green]✓ Updated {changed_count} simulation(s)[/green]")
            else:
                self.console.print("[yellow]No status changes made[/yellow]")
                
        except (ValueError, IndexError):
            self.console.print("[red]Invalid input format[/red]")
    
    def _set_all_status(self, status: str):
        """Set all simulations to the same status."""
        changed_count = 0
        for config in self.simulation_queue.queue:
            if config.status != status:
                config.status = status
                changed_count += 1
        
        if changed_count > 0:
            self.simulation_queue._sync_to_workspace()
            self.console.print(f"[green]✓ Set all {len(self.simulation_queue)} simulation(s) to {status.upper()}[/green]")
        else:
            self.console.print(f"[yellow]All simulations already {status.upper()}[/yellow]")
        
    def _remove_simulations(self):
        """Remove selected simulations from queue."""
        if not self.simulation_queue:
            self.console.print("[yellow]No simulations to remove[/yellow]")
            return
            
        self.console.print("\n[bold]Select simulations to remove:[/bold]")
        for i, config in enumerate(self.simulation_queue.queue, 1):
            self.console.print(f"  {i}. {config.name}")
            
        try:
            selections = prompt_with_context(
                self.processor,
                "Enter simulation numbers to remove (comma-separated)",
                module="MD Manager - Queue Management",
                description="Select simulations to remove"
            ).strip()
            indices = [int(x.strip()) - 1 for x in selections.split(',')]
            
            # Validate indices
            valid_indices = [i for i in indices if 0 <= i < len(self.simulation_queue)]
            if not valid_indices:
                self.console.print("[red]No valid simulation numbers provided[/red]")
                return
                
            # Remove in reverse order to maintain correct indices
            for i in sorted(valid_indices, reverse=True):
                removed_sim = self.simulation_queue.queue[i]
                self.simulation_queue.queue.pop(i)
                self.console.print(f"[yellow]Removed: {removed_sim.name}[/yellow]")
                
            # Update workflow steps after removal
            for step_num, config in enumerate(self.simulation_queue.queue, 1):
                if config.workflow_step:
                    config.workflow_step = step_num
                    
            self.simulation_queue._sync_to_workspace()
            self.console.print(f"[green]✓ Removed {len(valid_indices)} simulation(s)[/green]")
            
        except (ValueError, IndexError):
            self.console.print("[red]Invalid input format[/red]")

    def _show_settings_help(self) -> bool:
        """Show settings and help."""
        self.console.print(f"\n[bold cyan]===== Settings & Help =====[/bold cyan]")
        
        self.console.print(f"\n[bold]Available AMBER Engines:[/bold]")
        self.console.print("  • sander      — Standard CPU engine with broadest feature support")
        self.console.print("  • pmemd       — Optimized CPU engine (single core)")
        self.console.print("  • pmemd.MPI   — Optimized CPU engine parallelized across multiple cores")
        self.console.print("  • pmemd.cuda  — GPU-accelerated engine, requires NVIDIA GPU")
        
        self.console.print(f"\n[bold]Hardware Configuration:[/bold]")
        self.console.print("  • For CPU: Specify number of MPI tasks")
        self.console.print("  • For GPU: Specify device ID(s), comma-separated for multi-GPU")
        
        self.console.print(f"\n[bold]File Requirements:[/bold]")
        self.console.print("  • Topology file (.prmtop) from Topology Generator")
        self.console.print("  • Coordinate file (.rst7/.inpcrd) from Topology Generator")
        self.console.print("  • Input templates (.mdin) from template setup")
        
        input("Press Enter to continue...")
        return True

    def _view_simulation_queue(self):
        """View detailed simulation queue."""
        if not self.simulation_queue:
            self.console.print("[yellow]No simulations in queue[/yellow]")
            return
            
        self.console.print(f"\n[bold]Detailed Queue View:[/bold]")
        
        for i, config in enumerate(self.simulation_queue.queue, 1):
            self.console.print(f"\n[bold cyan]{i}. {config.name}[/bold cyan]")
            self.console.print(f"   Template ID: {config.template_id}")
            self.console.print(f"   Input file: {config.mdin_path}")
            if hasattr(config, 'engine') and config.engine:
                self.console.print(f"   Engine: {config.engine}")
                if config.mpi_tasks:
                    self.console.print(f"   MPI tasks: {config.mpi_tasks}")
                if config.gpu_ids:
                    self.console.print(f"   GPU IDs: {config.gpu_ids}")
        
        input("\nPress Enter to continue...")

    def workspace_overview(self) -> bool:
        """Display comprehensive overview of all available MD files and resources."""
        self._display_workspace_overview()
        
        self.console.print(f"\n[bold]Navigation:[/bold]")
        self.console.print("  b. ← Back to main menu")

        prompt_with_context(
            self.processor,
            "Press Enter to continue",
            choices=["b"],
            default="b",
            module="MD Manager - Workspace Overview",
            description="Return to main menu"
        )
        return True

    def _add_simulation_to_queue(self, controller):
        """Add new simulation to queue."""
        self.console.print(f"\n[bold]Add Simulation to Queue[/bold]")
        
        # Get available templates
        templates = self.user_data_manager.list_templates()
        if not templates:
            self.console.print("[yellow]No templates available. Create templates first.[/yellow]")
            return
            
        # Use the template selection interface from the controller
        self.console.print("\n[bold]Available Templates:[/bold]")
        
        # Display categorized templates
        template_choices = controller._display_categorized_templates(templates, show_numbers=True)
        
        # Add back option
        back_num = len(templates) + 1
        template_choices[str(back_num)] = None
        self.console.print(f"\n[bold]Options:[/bold]")
        self.console.print(f"  {back_num}. ← Back")

        # Build options map for context
        context_options = {str(back_num): "← Back"}
        for key, template_id in template_choices.items():
            if template_id:  # Skip the back option
                template_name = templates[template_id].get('name', template_id)
                context_options[key] = f"Select: {template_name}"

        choice = prompt_with_context(
            self.processor,
            "Select template",
            choices=list(template_choices.keys()),
            default=str(back_num),
            module="MD Manager - Single Simulation Setup",
            description="Select template for simulation",
            options_map=context_options
        )
        selected_template = template_choices[choice]

        if not selected_template:
            return

        # Get simulation details
        template_metadata = templates[selected_template]
        default_name = f"{template_metadata['name'].replace(' ', '_')}"

        sim_name = prompt_with_context(
            self.processor,
            "Simulation name",
            default=default_name,
            module="MD Manager - Single Simulation Setup",
            description="Enter simulation name"
        )
        
        # Ask for insertion position
        if self.simulation_queue:
            self.console.print(f"\nCurrent queue has {len(self.simulation_queue)} simulations")
            position = prompt_with_context(
                self.processor,
                f"Insert at position (1-{len(self.simulation_queue)+1})",
                default=str(len(self.simulation_queue)+1),
                module="MD Manager - Queue Management",
                description="Select insertion position in queue"
            )
            try:
                pos = int(position) - 1  # Convert to 0-based
                if pos < 0 or pos > len(self.simulation_queue):
                    pos = len(self.simulation_queue)  # Append at end
            except ValueError:
                pos = len(self.simulation_queue)  # Append at end
        else:
            pos = 0
            
        # Create simulation config
        # For now, we'll configure hardware later during execution
        config = SimulationConfig(
            name=sim_name,
            template_id=selected_template,
            mdin_path=template_metadata.get("template_path", selected_template),
            engine=""  # Will be configured during execution
        )
        
        self.simulation_queue.add_simulation(config, pos if pos < len(self.simulation_queue) else None)
        
        self.console.print(f"[green]✓ Added '{sim_name}' to queue at position {pos+1}[/green]")

    def _modify_queue_simulation(self):
        """Modify simulation in queue."""
        if not self.simulation_queue:
            self.console.print("[yellow]No simulations in queue to modify[/yellow]")
            return
            
        self.console.print(f"\n[bold cyan]===== Modify Queue Simulation =====[/bold cyan]")
        self.simulation_queue.display_queue(self.console)

        try:
            choice_str = prompt_with_context(
                self.processor,
                f"Select simulation to modify (1-{len(self.simulation_queue)})",
                default="1",
                module="MD Manager - Queue Modification",
                description="Select simulation to modify"
            )
            choice = int(choice_str)
            
            if 1 <= choice <= len(self.simulation_queue):
                idx = choice - 1
                config = self.simulation_queue.queue[idx]
                
                self.console.print(f"\n[bold]Modifying: {config.name}[/bold]")
                self.console.print(f"Current settings:")
                self.console.print(f"  Name: [cyan]{config.name}[/cyan]")
                self.console.print(f"  Template: [cyan]{config.template_id}[/cyan]")  
                self.console.print(f"  MDIN path: [cyan]{config.mdin_path}[/cyan]")
                self.console.print(f"  Engine: [cyan]{config.engine}[/cyan]")
                if config.mpi_tasks:
                    self.console.print(f"  MPI tasks: [cyan]{config.mpi_tasks}[/cyan]")
                if config.gpu_ids:
                    self.console.print(f"  GPU IDs: [cyan]{config.gpu_ids}[/cyan]")
                
                self.console.print(f"\n[bold]What would you like to modify?[/bold]")
                self.console.print("1. Simulation name", highlight=False)
                self.console.print("2. Template ID", highlight=False)
                self.console.print("3. MDIN file path", highlight=False)
                self.console.print("4. Execution engine", highlight=False)
                self.console.print("5. Hardware settings (MPI/GPU)", highlight=False)
                self.console.print("6. Cancel", highlight=False)

                modify_choice = prompt_with_context(
                    self.processor,
                    "Select option",
                    choices=["1","2","3","4","5","6"],
                    default="6",
                    module="MD Manager - Queue Modification",
                    description="Select modification option",
                    options_map={
                        "1": "Simulation name",
                        "2": "Template ID",
                        "3": "MDIN file path",
                        "4": "Execution engine",
                        "5": "Hardware settings (MPI/GPU)",
                        "6": "Cancel"
                    }
                )

                if modify_choice == "1":
                    new_name = prompt_with_context(
                        self.processor,
                        "New simulation name",
                        default=config.name,
                        module="MD Manager - Queue Modification",
                        description="Enter new simulation name"
                    )
                    config.name = new_name
                    self.console.print(f"[green]✓ Updated name to: {new_name}[/green]")

                elif modify_choice == "2":
                    new_template = prompt_with_context(
                        self.processor,
                        "New template ID",
                        default=config.template_id,
                        module="MD Manager - Queue Modification",
                        description="Enter new template ID"
                    )
                    config.template_id = new_template
                    self.console.print(f"[green]✓ Updated template ID to: {new_template}[/green]")

                elif modify_choice == "3":
                    new_mdin = prompt_with_context(
                        self.processor,
                        "New MDIN file path",
                        default=config.mdin_path,
                        module="MD Manager - Queue Modification",
                        description="Enter new MDIN file path"
                    )
                    confirm_str = "y"
                    if not Path(new_mdin).exists():
                        confirm_str = prompt_with_context(
                            self.processor,
                            f"File {new_mdin} doesn't exist. Use anyway?",
                            choices=["y", "n"],
                            default="n",
                            module="MD Manager - Queue Modification",
                            description="Confirm non-existent file",
                            options_map={"y": "Yes, use anyway", "n": "No, cancel"}
                        )
                    if Path(new_mdin).exists() or confirm_str.lower() == "y":
                        config.mdin_path = new_mdin
                        self.console.print(f"[green]✓ Updated MDIN path to: {new_mdin}[/green]")
                    else:
                        self.console.print("[yellow]MDIN path not changed[/yellow]")
                        
                elif modify_choice == "4":
                    self.console.print("Available engines:")
                    self.console.print("1. sander      [grey50]— standard CPU, broadest feature support[/grey50]", highlight=False)
                    self.console.print("2. pmemd       [grey50]— optimized CPU, single core[/grey50]", highlight=False)
                    self.console.print("3. pmemd.MPI   [grey50]— optimized CPU, multi-core parallel[/grey50]", highlight=False)
                    self.console.print("4. pmemd.cuda  [grey50]— GPU-accelerated, requires NVIDIA GPU[/grey50]", highlight=False)

                    engine_choice = prompt_with_context(
                        self.processor,
                        "Select engine",
                        choices=["1","2","3","4"],
                        default="2",
                        module="MD Manager - Queue Modification",
                        description="Select MD engine",
                        options_map={
                            "1": "sander (single CPU)",
                            "2": "pmemd (single CPU, optimized)",
                            "3": "pmemd.MPI (multi-CPU)",
                            "4": "pmemd.cuda (GPU acceleration)"
                        }
                    )
                    engine_map = {
                        "1": "sander",
                        "2": "pmemd",
                        "3": "pmemd.MPI",
                        "4": "pmemd.cuda"
                    }

                    new_engine = engine_map[engine_choice]
                    config.engine = new_engine

                    # Clear previous hardware settings
                    config.mpi_tasks = None
                    config.gpu_ids = None

                    # Set new hardware settings if needed
                    if new_engine == "pmemd.MPI":
                        cpu_info = self._get_cpu_info()
                        tasks_str = prompt_with_context(
                            self.processor,
                            "Number of MPI tasks",
                            default=str(min(16, cpu_info['available'])),
                            module="MD Manager - Queue Modification",
                            description="Enter number of MPI tasks"
                        )
                        tasks = int(tasks_str)
                        config.mpi_tasks = tasks
                    elif new_engine == "pmemd.cuda":
                        gpu_ids = prompt_with_context(
                            self.processor,
                            "GPU IDs (comma-separated)",
                            default="0",
                            module="MD Manager - Queue Modification",
                            description="Enter GPU device IDs"
                        )
                        config.gpu_ids = gpu_ids
                        
                    self.console.print(f"[green]✓ Updated engine to: {new_engine}[/green]")
                    
                elif modify_choice == "5":
                    if config.engine == "pmemd.MPI":
                        new_tasks_str = prompt_with_context(
                            self.processor,
                            "Number of MPI tasks",
                            default=str(config.mpi_tasks or 4),
                            module="MD Manager - Queue Modification",
                            description="Enter number of MPI tasks"
                        )
                        new_tasks = int(new_tasks_str)
                        config.mpi_tasks = new_tasks
                        self.console.print(f"[green]✓ Updated MPI tasks to: {new_tasks}[/green]")
                    elif config.engine == "pmemd.cuda":
                        new_gpu_ids = prompt_with_context(
                            self.processor,
                            "GPU IDs (comma-separated)",
                            default=config.gpu_ids or "0",
                            module="MD Manager - Queue Modification",
                            description="Enter GPU device IDs"
                        )
                        config.gpu_ids = new_gpu_ids
                        self.console.print(f"[green]✓ Updated GPU IDs to: {new_gpu_ids}[/green]")
                    else:
                        self.console.print(f"[yellow]No hardware settings available for {config.engine}[/yellow]")
                        
                elif modify_choice == "6":
                    self.console.print("[yellow]Modification cancelled[/yellow]")
                    return
                    
                self.console.print(f"\n[green]Simulation successfully modified[/green]")
                
            else:
                self.console.print("[red]Invalid selection[/red]")
                
        except KeyboardInterrupt:
            self.console.print("[yellow]Modification cancelled[/yellow]")

    def _remove_from_queue(self):
        """Remove simulation from queue."""
        if not self.simulation_queue:
            self.console.print("[yellow]No simulations in queue to remove[/yellow]")
            return
            
        self.simulation_queue.display_queue(self.console)

        try:
            choice_str = prompt_with_context(
                self.processor,
                f"Select simulation to remove (1-{len(self.simulation_queue)})",
                default="1",
                module="MD Manager - Queue Management",
                description="Select simulation to remove"
            )
            choice = int(choice_str)

            if 1 <= choice <= len(self.simulation_queue):
                config = self.simulation_queue.queue[choice-1]

                confirm_str = prompt_with_context(
                    self.processor,
                    f"Remove '{config.name}' from queue?",
                    choices=["y", "n"],
                    default="n",
                    module="MD Manager - Queue Management",
                    description="Confirm removal",
                    options_map={"y": "Yes, remove", "n": "No, cancel"}
                )
                if confirm_str.lower() == "y":
                    self.simulation_queue.remove_simulation(choice-1)
                    self.console.print(f"[green]✓ Removed '{config.name}' from queue[/green]")
            else:
                self.console.print("[red]Invalid selection[/red]")
                
        except KeyboardInterrupt:
            self.console.print("[yellow]Removal cancelled[/yellow]")

    def _reorder_queue(self):
        """Reorder simulation queue with command-line interface."""
        if len(self.simulation_queue) < 2:
            self.console.print("[yellow]Need at least 2 simulations to reorder[/yellow]")
            return
            
        self.console.print(f"\n[bold cyan]===== Reorder Simulation Queue =====[/bold cyan]")
        self.console.print("[grey50]Commands: 'move X to Y', 'swap X Y', 'help', or 'done'[/grey50]")
        
        while True:
            # Show current queue with indices
            self.simulation_queue.display_queue(self.console)
            
            try:
                command = prompt_with_context(
                    self.processor,
                    f"\n[bold blue]queue>[/bold blue]",
                    default="done",
                    module="MD Manager - Queue Reordering",
                    description="Enter queue command (move, swap, help, or done)"
                ).strip().lower()
                
                if command == "done":
                    self.console.print("[green]✓ Queue reordering complete[/green]")
                    break
                elif command == "help":
                    self._show_reorder_help()
                    continue
                elif command.startswith("move"):
                    self._handle_move_command(command)
                elif command.startswith("swap"):
                    self._handle_swap_command(command)
                else:
                    self.console.print("[red]Unknown command. Type 'help' for available commands.[/red]")
                    
            except KeyboardInterrupt:
                self.console.print("\n[yellow]Queue reordering cancelled[/yellow]")
                break
            except Exception as e:
                self.console.print(f"[red]Error: {e}[/red]")

    def _show_reorder_help(self):
        """Show help for queue reordering commands."""
        self.console.print(f"\n[bold]Queue Reordering Commands:[/bold]")
        self.console.print("  [cyan]move X to Y[/cyan]  - Move simulation at position X to position Y")
        self.console.print("  [cyan]swap X Y[/cyan]     - Swap simulations at positions X and Y")
        self.console.print("  [cyan]help[/cyan]         - Show this help message")
        self.console.print("  [cyan]done[/cyan]         - Finish reordering and return to menu")
        self.console.print("\nExamples:")
        self.console.print("  [grey50]move 3 to 1[/grey50]   - Move simulation #3 to position #1")
        self.console.print("  [grey50]swap 2 4[/grey50]      - Swap simulations at positions #2 and #4")

    def _handle_move_command(self, command: str):
        """Handle 'move X to Y' command."""
        try:
            # Parse "move X to Y"
            parts = command.split()
            if len(parts) != 4 or parts[2] != "to":
                self.console.print("[red]Invalid syntax. Use: move X to Y[/red]")
                return
                
            from_pos = int(parts[1])
            to_pos = int(parts[3])
            
            # Convert to 0-based indices
            from_idx = from_pos - 1
            to_idx = to_pos - 1
            
            # Validate indices
            queue_len = len(self.simulation_queue)
            if not (0 <= from_idx < queue_len):
                self.console.print(f"[red]Invalid source position {from_pos}. Must be 1-{queue_len}[/red]")
                return
            if not (0 <= to_idx < queue_len):
                self.console.print(f"[red]Invalid destination position {to_pos}. Must be 1-{queue_len}[/red]")
                return
                
            if from_idx == to_idx:
                self.console.print("[yellow]Source and destination are the same[/yellow]")
                return
                
            # Get simulation name for feedback
            sim_name = self.simulation_queue.queue[from_idx].name
            
            # Perform the move
            if self.simulation_queue.move_simulation(from_idx, to_idx):
                self.console.print(f"[green]✓ Moved '{sim_name}' from position {from_pos} to {to_pos}[/green]")
            else:
                self.console.print(f"[red]Failed to move simulation[/red]")
                
        except ValueError:
            self.console.print("[red]Invalid position numbers. Use integers only.[/red]")
        except Exception as e:
            self.console.print(f"[red]Error parsing move command: {e}[/red]")

    def _handle_swap_command(self, command: str):
        """Handle 'swap X Y' command."""
        try:
            # Parse "swap X Y"
            parts = command.split()
            if len(parts) != 3:
                self.console.print("[red]Invalid syntax. Use: swap X Y[/red]")
                return
                
            pos1 = int(parts[1])
            pos2 = int(parts[2])
            
            # Convert to 0-based indices  
            idx1 = pos1 - 1
            idx2 = pos2 - 1
            
            # Validate indices
            queue_len = len(self.simulation_queue)
            if not (0 <= idx1 < queue_len):
                self.console.print(f"[red]Invalid position {pos1}. Must be 1-{queue_len}[/red]")
                return
            if not (0 <= idx2 < queue_len):
                self.console.print(f"[red]Invalid position {pos2}. Must be 1-{queue_len}[/red]")
                return
                
            if idx1 == idx2:
                self.console.print("[yellow]Cannot swap a simulation with itself[/yellow]")
                return
                
            # Get simulation names for feedback
            sim1_name = self.simulation_queue.queue[idx1].name
            sim2_name = self.simulation_queue.queue[idx2].name
            
            # Perform the swap using SimulationQueue's queue attribute directly
            # (since move_simulation doesn't support swapping)
            self.simulation_queue.queue[idx1], self.simulation_queue.queue[idx2] = \
                self.simulation_queue.queue[idx2], self.simulation_queue.queue[idx1]
                
            self.console.print(f"[green]✓ Swapped '{sim1_name}' and '{sim2_name}' (positions {pos1} ↔ {pos2})[/green]")
            
        except ValueError:
            self.console.print("[red]Invalid position numbers. Use integers only.[/red]")
        except Exception as e:
            self.console.print(f"[red]Error parsing swap command: {e}[/red]")

    def _import_mdin_to_queue(self):
        """Import .mdin file to queue."""
        self.console.print(f"\n[bold]Import .mdin File[/bold]")

        filepath = prompt_with_context(
            self.processor,
            "Path to .mdin file",
            module="MD Manager - Import",
            description="Enter path to .mdin file"
        )
        
        if not Path(filepath).exists():
            self.console.print("[red]File not found[/red]")
            return
            
        try:
            # Read file content
            with open(filepath, 'r') as f:
                content = f.read()
                
            # Try to parse header for metadata
            metadata = self.user_data_manager._parse_template_header(content)
            
            # Get simulation details
            if metadata.get('name'):
                default_name = metadata['name'].replace(' ', '_')
            else:
                default_name = Path(filepath).stem
                
            sim_name = prompt_with_context(
                self.processor,
                "Simulation name",
                default=default_name,
                module="MD Manager - Import",
                description="Enter simulation name"
            )
            
            # Create temporary template ID
            import uuid
            temp_template_id = f"imported_{uuid.uuid4().hex[:8]}"
            
            # Save as temporary file in workspace for execution
            # (This is a simplified approach - in full implementation, 
            # we'd integrate with workspace file management)
            temp_path = Path(filepath).resolve()
            
            # Create simulation config
            config = SimulationConfig(
                name=sim_name,
                template_id=temp_template_id,
                mdin_path=str(temp_path),
                engine=""
            )
            
            # Ask for insertion position
            if self.simulation_queue:
                position = prompt_with_context(
                    self.processor,
                    f"Insert at position (1-{len(self.simulation_queue)+1})",
                    default=str(len(self.simulation_queue)+1),
                    module="MD Manager - Import",
                    description="Select insertion position in queue"
                )
                try:
                    pos = int(position) - 1
                    if pos < 0 or pos > len(self.simulation_queue):
                        pos = None  # Append at end
                except ValueError:
                    pos = None
            else:
                pos = None
                
            self.simulation_queue.add_simulation(config, pos)
            
            self.console.print(f"[green]✓ Imported '{sim_name}' to queue[/green]")
            
        except Exception as e:
            self.console.print(f"[red]Error importing file: {e}[/red]")

    def _show_simulation_status(self):
        """Show simulation readiness status."""
        self.console.print(f"\n[bold]Simulation Readiness Status[/bold]")
        
        # Create status table
        table = Table()
        table.add_column("Component", style="bright_blue", width=20)
        table.add_column("Status", style="white", width=12)
        table.add_column("Details", style="grey50", width=40)
        
        # Check for topology file
        prmtop_files = self._find_workspace_files("*.prmtop")
        if prmtop_files:
            prmtop_file = prmtop_files[0]  # Use first one found
            table.add_row(
                "Topology",
                "✓ Ready",
                f"{prmtop_file.name}"
            )
        else:
            table.add_row(
                "Topology",
                "✗ Missing", 
                "Run Topology Generator first"
            )
            
        # Check for coordinate files
        coord_files = (self._find_workspace_files("*.rst7") + 
                      self._find_workspace_files("*.inpcrd"))
        if coord_files:
            coord_file = coord_files[0]  # Use first one found
            table.add_row(
                "Coordinates",
                "✓ Ready",
                f"{coord_file.name}"
            )
        else:
            table.add_row(
                "Coordinates",
                "✗ Missing",
                "Run Topology Generator first"
            )
            
        # Check single simulation queue
        if self.simulation_queue:
            table.add_row(
                "Single Simulations",
                f"✓ {len(self.simulation_queue)} queued",
                ", ".join([config.name for config in self.simulation_queue.queue[:3]]) + 
                ("..." if len(self.simulation_queue) > 3 else "")
            )
        else:
            table.add_row(
                "Single Simulations",
                "- None",
                "Setup simulations first"
            )
            
        # Check for workflow files
        workflow_files = self._find_workspace_files("*_workflow.json")
        if workflow_files:
            table.add_row(
                "Protocols",
                f"✓ {len(workflow_files)} available",
                ", ".join([f.stem.replace('_workflow', '') for f in workflow_files[:2]]) +
                ("..." if len(workflow_files) > 2 else "")
            )
        else:
            table.add_row(
                "Protocols",
                "- None",
                "Create protocols first"
            )
            
        # Check for running simulations
        if self.workflow_core.checkpoint_manager.checkpoint_exists():
            table.add_row(
                "Active Simulation",
                "🔄 Checkpoint exists",
                "Simulation may be paused or failed"
            )
            
        self.console.print(table)
        
        # Summary message
        prerequisites_ready = bool(prmtop_files and coord_files)
        simulations_ready = bool(self.simulation_queue or workflow_files)
        
        if not prerequisites_ready:
            self.console.print("\n[red]⚠ Missing required files from TLEaP[/red]")
        elif not simulations_ready:
            self.console.print("\n[yellow]⚠ No simulations configured[/yellow]")
        else:
            ready_count = len(self.simulation_queue) + len(workflow_files)
            self.console.print(f"\n[green]✓ {ready_count} simulation(s) ready to execute[/green]")

    def _find_workspace_files(self, pattern: str) -> List[Path]:
        """Find files in workspace matching pattern."""
        try:
            return list(Path.cwd().glob(pattern))
        except Exception:
            return []

    def _extract_pdb_from_topology(self) -> Optional[str]:
        """
        Extract a PDB file from available topology (prmtop/rst7) using cpptraj.

        Returns:
            Path to extracted PDB file, or None if extraction failed.
        """
        import subprocess
        import tempfile

        # Find prmtop and rst7 files
        prmtop_file = None
        rst7_file = None

        # Check simulation queue first
        if self.simulation_queue and len(self.simulation_queue) > 0:
            first_config = self.simulation_queue.queue[0]
            if first_config.prmtop and Path(first_config.prmtop).exists():
                prmtop_file = first_config.prmtop
            if first_config.rst7 and Path(first_config.rst7).exists():
                rst7_file = first_config.rst7

        # Check workspace for structure pairs
        if not prmtop_file and self.workspace:
            structure_pairs = self.workspace.get('md_structure_pairs', [])
            if structure_pairs:
                pair = structure_pairs[0]
                if pair.get('prmtop') and Path(pair['prmtop']).exists():
                    prmtop_file = pair['prmtop']
                if pair.get('rst7') and Path(pair['rst7']).exists():
                    rst7_file = pair['rst7']

        # Fall back to scanning directory
        if not prmtop_file:
            prmtop_files = self._find_workspace_files("*.prmtop")
            if prmtop_files:
                prmtop_file = str(prmtop_files[0])

        if not rst7_file:
            rst7_files = self._find_workspace_files("*.rst7") + self._find_workspace_files("*.inpcrd")
            if rst7_files:
                rst7_file = str(rst7_files[0])

        if not prmtop_file or not rst7_file:
            return None

        # Generate output PDB path - use project directory from workspace
        prmtop_path = Path(prmtop_file).resolve()

        # Use project_directory from workspace (set by main.py at session start)
        output_dir = None
        if self.workspace:
            output_dir = self.workspace.get('project_directory')

        if output_dir and Path(output_dir).is_dir():
            output_pdb = Path(output_dir) / f"{prmtop_path.stem}_extracted.pdb"
        else:
            # Fall back to same directory as prmtop
            output_pdb = prmtop_path.parent / f"{prmtop_path.stem}_extracted.pdb"

        # Create cpptraj input
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.in', delete=False) as f:
                f.write(f"parm {prmtop_file}\n")
                f.write(f"trajin {rst7_file}\n")
                f.write(f"trajout {output_pdb} pdb\n")
                f.write("go\n")
                cpptraj_input = f.name

            self.console.print(f"[grey50]Extracting PDB from topology using cpptraj...[/grey50]")

            result = subprocess.run(
                ['cpptraj', '-i', cpptraj_input],
                capture_output=True,
                text=True,
                timeout=60
            )

            # Clean up input file
            Path(cpptraj_input).unlink(missing_ok=True)

            if result.returncode == 0 and output_pdb.exists():
                self.console.print(f"[green]✓ Extracted PDB: {output_pdb.name}[/green]")
                return str(output_pdb)
            else:
                self.console.print(f"[red]✗ cpptraj failed: {result.stderr[:200]}[/red]")
                return None

        except subprocess.TimeoutExpired:
            self.console.print("[red]✗ cpptraj timed out[/red]")
            return None
        except FileNotFoundError:
            self.console.print("[red]✗ cpptraj not found in PATH[/red]")
            return None
        except Exception as e:
            self.console.print(f"[red]✗ Error extracting PDB: {e}[/red]")
            return None

    def _get_md_ready_status(self) -> Dict[str, Any]:
        """Get MD readiness status from workspace - single source of truth."""
        status = {
            'ready': False,
            'structure_pairs': [],
            'template_assignments': [],
            'simulation_queue': [],
            'missing': [],
            'has_topology': False,
            'has_coordinates': False,
            'amber_ready': False
        }
        
        # Get structure pairs from workspace
        if self.workspace:
            structure_pairs = self.workspace.get('md_structure_pairs', [])
            status['structure_pairs'] = structure_pairs
            
            # Verify files exist
            valid_pairs = []
            for pair in structure_pairs:
                prmtop_path = Path(pair['prmtop'])
                rst7_path = Path(pair['rst7'])
                if prmtop_path.exists() and rst7_path.exists():
                    valid_pairs.append(pair)
                    status['has_topology'] = True
                    status['has_coordinates'] = True
            
            # Get template assignments
            status['template_assignments'] = self.workspace.get('md_template_assignments', [])
            
            # Get simulation queue
            queue_data = self.workspace.get('md_simulation_queue', [])
            status['simulation_queue'] = queue_data
            
            # Check structure files from both single sim pairs AND workflow simulations
            has_valid_structures = bool(valid_pairs)
            
            # Also check for structure files in simulation queue (from workflows)
            if not has_valid_structures and queue_data:
                for sim in queue_data:
                    prmtop_path = sim.get('prmtop')
                    rst7_path = sim.get('rst7')
                    if prmtop_path and rst7_path:
                        prmtop_file = Path(prmtop_path)
                        rst7_file = Path(rst7_path)
                        if prmtop_file.exists() and rst7_file.exists():
                            has_valid_structures = True
                            status['has_topology'] = True
                            status['has_coordinates'] = True
                            break
            
            # Determine readiness
            if has_valid_structures and queue_data:
                status['amber_ready'] = True
                status['ready'] = True
            else:
                if not has_valid_structures:
                    status['missing'].append('Valid structure file pairs (prmtop/rst7)')
                if not queue_data:
                    status['missing'].append('Configured simulations')
        
        return status
    
    def _detect_workspace_status(self) -> Dict[str, Any]:
        """Detect and analyze workspace file status for MD simulations."""
        status = {
            "topology_files": [],
            "coordinate_files": [],
            "template_files": [],
            "workflow_files": [],
            "simulation_directories": [],
            "amber_ready": False,
            "missing_requirements": [],
            "warnings": []
        }
        
        # Look for topology files (.prmtop, .parm7)
        topology_patterns = ["*.prmtop", "*.parm7"]
        for pattern in topology_patterns:
            status["topology_files"].extend(self._find_workspace_files(pattern))
        
        # Look for coordinate files (.rst7, .inpcrd, .crd)
        coordinate_patterns = ["*.rst7", "*.inpcrd", "*.crd"]
        for pattern in coordinate_patterns:
            status["coordinate_files"].extend(self._find_workspace_files(pattern))
            
        # Also check files referenced in simulation queue configurations
        queue_topology_files = []
        queue_coordinate_files = []
        
        # Check current simulation queue
        if hasattr(self, 'simulation_queue') and self.simulation_queue:
            for config in self.simulation_queue.queue:
                if config.prmtop and Path(config.prmtop).exists():
                    queue_topology_files.append(Path(config.prmtop))
                if config.rst7 and Path(config.rst7).exists():
                    queue_coordinate_files.append(Path(config.rst7))
        
        # Check saved simulation queue files
        workspace_queue_file = Path.cwd() / "simulation_queue.json"
        if workspace_queue_file.exists():
            try:
                with open(workspace_queue_file, 'r') as f:
                    queue_data = json.load(f)
                for sim in queue_data.get('simulations', []):
                    prmtop_path = sim.get('prmtop')
                    rst7_path = sim.get('rst7')
                    if prmtop_path and Path(prmtop_path).exists():
                        queue_topology_files.append(Path(prmtop_path))
                    if rst7_path and Path(rst7_path).exists():
                        queue_coordinate_files.append(Path(rst7_path))
            except Exception:
                pass
        
        # Add queue files to status (avoid duplicates)
        for f in queue_topology_files:
            if f not in status["topology_files"]:
                status["topology_files"].append(f)
        for f in queue_coordinate_files:
            if f not in status["coordinate_files"]:
                status["coordinate_files"].append(f)
            
        # Look for template files (.mdin)
        status["template_files"] = self._find_workspace_files("*.mdin")
        
        # Look for workflow files (.json)
        json_files = self._find_workspace_files("*.json")
        # Filter for workflow files (contain workflow-like structure)
        for json_file in json_files:
            try:
                import json
                with open(json_file, 'r') as f:
                    data = json.load(f)
                # Simple check for workflow structure
                if isinstance(data, dict) and ("steps" in data or "workflow" in data or "simulations" in data):
                    status["workflow_files"].append(json_file)
            except:
                continue
                
        # Also check for custom workflows and active workflows in simulation queue
        if hasattr(self, 'user_data_manager') and self.user_data_manager:
            try:
                # Check custom workflow directory
                custom_workflow_dir = Path(self.user_data_manager.user_workflow_dir) / "custom"
                if custom_workflow_dir.exists():
                    custom_workflows = list(custom_workflow_dir.glob("*.json"))
                    status["workflow_files"].extend(custom_workflows)
                    
                # Check for active workflows in simulation queue
                if hasattr(self, 'simulation_queue') and self.simulation_queue:
                    active_workflows = set()
                    for config in self.simulation_queue.queue:
                        if hasattr(config, 'workflow_id') and config.workflow_id:
                            active_workflows.add(config.workflow_id)
                    # Add count of active workflows if not already counted
                    if active_workflows and not status["workflow_files"]:
                        status["workflow_files"] = list(active_workflows)  # Just store workflow IDs
            except Exception:
                pass  # Fail gracefully if workflow detection fails
                
        # Look for existing simulation directories
        simulations_path = Path.cwd() / "simulations"
        if simulations_path.exists():
            status["simulation_directories"] = [
                d for d in simulations_path.iterdir() 
                if d.is_dir() and not d.name.startswith('.')
            ]
            
        # Determine if AMBER-ready
        has_topology = bool(status["topology_files"])
        has_coordinates = bool(status["coordinate_files"])
        status["amber_ready"] = has_topology and has_coordinates
        
        # Check for missing requirements
        if not has_topology:
            status["missing_requirements"].append("Topology file (.prmtop/.parm7) - Run TLEaP first")
        if not has_coordinates:
            status["missing_requirements"].append("Coordinate file (.rst7/.inpcrd) - Run TLEaP first")
            
        # Add warnings for common issues
        if len(status["topology_files"]) > 1:
            status["warnings"].append(f"Multiple topology files found ({len(status['topology_files'])})")
        if len(status["coordinate_files"]) > 1:
            status["warnings"].append(f"Multiple coordinate files found ({len(status['coordinate_files'])})")
            
        # Check if files are paired correctly (same base name)
        if has_topology and has_coordinates:
            topo_bases = {f.stem for f in status["topology_files"]}
            coord_bases = {f.stem for f in status["coordinate_files"]}
            if not topo_bases.intersection(coord_bases):
                status["warnings"].append("Topology and coordinate files may not be paired (different base names)")
                
        return status

    def _display_workspace_overview(self):
        """Display comprehensive overview of all available MD files and resources."""
        self.console.print(f"\n[bold cyan]===== Available Files & Workspace Overview =====[/bold cyan]")
        
        status = self._detect_workspace_status()
        
        # Create status table  
        from rich.table import Table
        table = Table(title="Available Resources Overview", show_header=True, header_style="bold blue")
        table.add_column("Component", style="bright_blue", width=20)
        table.add_column("Status", width=15)
        table.add_column("Details", style="grey50")
        
        # Topology files
        if status["topology_files"]:
            file_names = [f.name for f in status["topology_files"]]
            table.add_row(
                "Topology Files",
                f"[green]✓ {len(status['topology_files'])} found[/green]",
                ", ".join(file_names)
            )
        else:
            table.add_row(
                "Topology Files", 
                "[red]✗ Missing[/red]",
                "Run Topology Generator first"
            )
            
        # Coordinate files  
        if status["coordinate_files"]:
            file_names = [f.name for f in status["coordinate_files"]]
            table.add_row(
                "Coordinate Files",
                f"[green]✓ {len(status['coordinate_files'])} found[/green]",
                ", ".join(file_names)
            )
        else:
            table.add_row(
                "Coordinate Files",
                "[red]✗ Missing[/red]", 
                "Run Topology Generator first"
            )
            
        # Simulation queue
        if len(self.simulation_queue) > 0:
            queue_names = [config.name for config in self.simulation_queue.queue[:3]]
            details = ", ".join(queue_names)
            if len(self.simulation_queue) > 3:
                details += f" (and {len(self.simulation_queue) - 3} more)"
            table.add_row(
                "Queued Simulations",
                f"[green]✓ {len(self.simulation_queue)} ready[/green]",
                details
            )
        else:
            table.add_row(
                "Queued Simulations",
                "[yellow]○ None[/yellow]",
                "Use 'Setup single simulations' to add"
            )
            
        # Workflow files - check both JSON files and queued workflows
        json_workflows = status["workflow_files"]
        queued_workflows = []
        
        # Check for workflows in the simulation queue
        if self.workspace:
            workflows_data = self.workspace.get('md_workflows', {})
            queued_workflows = list(workflows_data.keys())
            
        total_workflows = len(json_workflows) + len(queued_workflows)
        
        if total_workflows > 0:
            # Build details string
            details_parts = []
            if json_workflows:
                json_names = [f.stem for f in json_workflows[:2]]
                details_parts.extend(json_names)
            if queued_workflows:
                # Get workflow names if available
                workflow_names = []
                for wf_id in queued_workflows[:2]:
                    wf_data = workflows_data.get(wf_id, {})
                    wf_name = wf_data.get('name', wf_id)
                    workflow_names.append(wf_name)
                details_parts.extend(workflow_names)
                
            details = ", ".join(details_parts[:2])
            if total_workflows > 2:
                details += f" (and {total_workflows - 2} more)"
                
            table.add_row(
                "Protocol Files",
                f"[green]✓ {total_workflows} available[/green]",
                details
            )
        else:
            table.add_row(
                "Protocol Files",
                "[yellow]○ None[/yellow]",
                "Use 'Setup protocols' or 'Import .json' to add"
            )
            
        # Simulation directories (previous runs)
        if status["simulation_directories"]:
            dir_names = [d.name for d in status["simulation_directories"][:3]]
            details = ", ".join(dir_names)
            if len(status["simulation_directories"]) > 3:
                details += f" (and {len(status['simulation_directories']) - 3} more)"
            table.add_row(
                "Previous Runs",
                f"[blue]📁 {len(status['simulation_directories'])} found[/blue]",
                details
            )
        else:
            table.add_row(
                "Previous Runs",
                "[grey50]○ None[/grey50]",
                "No previous simulations found"
            )
            
        self.console.print(table)
        
        # Overall status summary
        if status["amber_ready"]:
            # Count unique simulations - avoid double counting current queue if it's saved to workflow file
            queue_sims = len(self.simulation_queue)
            workflow_sims = 0
            
            # Only count workflow file simulations if no current queue (to avoid double counting)
            if queue_sims == 0:
                for workflow_file in status["workflow_files"]:
                    try:
                        with open(workflow_file, 'r') as f:
                            data = json.load(f)
                        workflow_sims += len(data.get('simulations', []))
                    except:
                        continue
            
            ready_items = queue_sims + workflow_sims
            if ready_items > 0:
                self.console.print(f"\n[bold green]✅ Ready for MD simulations ({ready_items} simulation(s) configured)[/bold green]")
            else:
                self.console.print(f"\n[bold yellow]⚠️  Files ready, but no simulations configured[/bold yellow]")
        else:
            self.console.print(f"\n[bold red]❌ Not ready for MD simulations[/bold red]")
            
        # Show missing requirements
        if status["missing_requirements"]:
            self.console.print(f"\n[bold red]Missing Requirements:[/bold red]")
            for req in status["missing_requirements"]:
                self.console.print(f"  • {req}")
                
        # Show warnings
        if status["warnings"]:
            self.console.print(f"\n[bold yellow]Warnings:[/bold yellow]") 
            for warning in status["warnings"]:
                self.console.print(f"  ⚠️  {warning}")
                
        return status
        # This is a simplified implementation
        # In full version, this would integrate with the workspace system
        try:
            if hasattr(self.processor, 'workspace') and self.processor.workspace:
                # Use workspace API if available
                workspace_dir = Path(self.processor.workspace.workspace_dir)
                return list(workspace_dir.glob(pattern))
            else:
                # Fallback: search current directory
                return list(Path.cwd().glob(pattern))
        except Exception:
            return []

    def _get_hardware_config(self, sim_name: str) -> Dict[str, Any]:
        """Get hardware configuration for simulation."""
        self.console.print(f"\n[bold]Configure hardware for: {sim_name}[/bold]")
        
        # Show available options
        engines = ["sander", "pmemd", "pmemd.MPI", "pmemd.cuda"]
        engine = prompt_with_context(
            self.processor,
            "Select AMBER engine",
            choices=engines,
            default="pmemd.cuda",
            module="MD Manager - Hardware Configuration",
            description="Select MD engine",
            options_map={
                "sander": "Single CPU",
                "pmemd": "Single CPU (optimized)",
                "pmemd.MPI": "Multi-CPU (MPI)",
                "pmemd.cuda": "GPU acceleration"
            }
        )

        config = {"engine": engine}

        if engine == "pmemd.MPI":
            # Get system info for default
            total_cpus = os.cpu_count()
            try:
                load_avg = os.getloadavg()[0] if hasattr(os, 'getloadavg') else 0
                idle_cpus = max(1, int(total_cpus - load_avg))
            except:
                idle_cpus = total_cpus // 2

            self.console.print(f"Available CPU cores: {total_cpus} ({idle_cpus} appear idle)")
            mpi_tasks_str = prompt_with_context(
                self.processor,
                "Number of MPI tasks",
                default=str(min(16, idle_cpus)),
                module="MD Manager - Hardware Configuration",
                description="Enter number of MPI tasks"
            )
            config["mpi_tasks"] = int(mpi_tasks_str)

        elif engine == "pmemd.cuda":
            gpu_info = self._get_gpu_info()
            if gpu_info['available'] > 0:
                self.console.print("Available GPUs:")
                for gpu in gpu_info['devices']:
                    self.console.print(f"  [{gpu['id']}] {gpu['name']} - {gpu['memory_free']} free")

                default_gpu = gpu_info['devices'][0]['id']
                config["gpu_ids"] = prompt_with_context(
                    self.processor,
                    "GPU device ID(s) (comma-separated for multi-GPU)",
                    default=default_gpu,
                    module="MD Manager - Hardware Configuration",
                    description="Enter GPU device ID(s)"
                )
            else:
                self.console.print("[yellow]No GPUs detected with nvidia-smi[/yellow]")
                config["gpu_ids"] = prompt_with_context(
                    self.processor,
                    "GPU device ID",
                    default="0",
                    module="MD Manager - Hardware Configuration",
                    description="Enter GPU device ID"
                )
                
        return config

    def _execute_simulation_queue(self):
        """Execute all simulations in queue with workflow dependency handling."""
        if not self.simulation_queue:
            self.console.print("[yellow]No simulations in queue to execute[/yellow]")
            return
            
        # Clear workflow output tracking for fresh execution
        self._workflow_outputs.clear()
            
        # Use unified status checker instead of redundant file checks
        md_status = self._get_md_ready_status()
        if not md_status["ready"]:
            self.console.print("[red]Cannot execute: Missing requirements[/red]")
            for missing_item in md_status["missing"]:
                self.console.print(f"  • {missing_item}")
            return
            
        self.console.print(f"\n[bold]Execute Simulation Queue[/bold]")
        self.simulation_queue.display_queue(self.console)
        
        # Count active simulations for confirmation
        active_count = sum(1 for config in self.simulation_queue.queue if config.status == "active")
        hold_count = len(self.simulation_queue.queue) - active_count
        
        # Check if there are any active simulations
        if active_count == 0:
            self.console.print(f"\n[yellow]No active simulations to execute (all are on hold)[/yellow]")
            return
        elif hold_count > 0:
            self.console.print(f"\n[grey50]Will execute {active_count} active simulations, skip {hold_count} on hold[/grey50]")

        confirm_str = prompt_with_context(
            self.processor,
            f"\nProceed with execution of {active_count} active simulations?",
            choices=["y", "n"],
            default="y",
            module="MD Manager - Execution",
            description="Confirm simulation batch execution",
            options_map={"y": "Yes, proceed with execution", "n": "No, cancel"}
        )
        if not (confirm_str.lower() == "y"):
            return
            
        # Store execution order for actual execution (use queue order)
        self.execution_order = list(self.simulation_queue.queue)
            
        # Create run directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = Path.cwd() / "simulations" / f"batch_{timestamp}"
        run_dir.mkdir(parents=True, exist_ok=True)
        
        self.console.print(f"\nCreated run directory: {run_dir}")
        
        # Copy shared files from the first simulation's structure files
        shared_dir = run_dir / "shared"
        shared_dir.mkdir(exist_ok=True)
        
        # Get structure files from the first simulation config
        first_config = self.simulation_queue.queue[0]
        prmtop_file = Path(first_config.prmtop)
        coord_file = Path(first_config.rst7)
        
        # Change working directory to where the structure files are located
        # This ensures all file operations happen in the correct context
        structure_dir = prmtop_file.parent.resolve()
        original_cwd = Path.cwd()
        
        self.console.print(f"[grey50]Changing to structure directory: {structure_dir}[/grey50]")
        import os
        os.chdir(structure_dir)
        
        # Now resolve file paths relative to the correct directory
        prmtop_file = prmtop_file.resolve()
        coord_file = coord_file.resolve()
        
        # Keep original filenames
        prmtop_name = prmtop_file.name
        coord_name = coord_file.name
        
        import shutil
        try:
            shutil.copy2(prmtop_file, shared_dir / prmtop_name)
            shutil.copy2(coord_file, shared_dir / coord_name)
            self.console.print(f"[grey50]Copied structure files to shared directory[/grey50]")
        except Exception as e:
            self.console.print(f"[red]Error copying files: {e}[/red]")
            raise
        
        # Execute each simulation in the queue  
        successful_runs = 0
        failed_runs = 0
        
        try:
            for i, sim_config in enumerate(self.simulation_queue.queue):
                # Check simulation status - skip if on hold
                if sim_config.status == "hold":
                    self.console.print(f"\n[bold yellow]Skipping simulation {i+1}/{len(self.simulation_queue)}: {sim_config.name} (ON HOLD)[/bold yellow]")
                    continue
                    
                self.console.print(f"\n[bold blue]Executing simulation {i+1}/{len(self.simulation_queue)}: {sim_config.name}[/bold blue]")
                
                # Create individual simulation directory
                sim_dir = run_dir / f"{i+1:02d}_{sim_config.name.replace(' ', '_')}"
                sim_dir.mkdir(exist_ok=True)
                
                # Get hardware configuration from simulation config
                hw_config = sim_config.hardware_config or {}
                
                # Execute the simulation
                try:
                    # Check if this simulation has dependencies (workflow step)
                    if sim_config.workflow_id and sim_config.depends_on:
                        # This is a dependent workflow step - mark as pending, don't start yet
                        self.console.print(f"[yellow]⏸️  Queued (depends on: {sim_config.depends_on})[/yellow]")
                        self.console.print(f"[grey50]Will auto-start when dependency completes[/grey50]")
                        # Store pending workflow info for background checker
                        if not hasattr(self, 'pending_workflow_steps'):
                            self.pending_workflow_steps = []
                        pending_info = {
                            'config': sim_config,
                            'sim_dir': sim_dir,
                            'shared_dir': shared_dir,
                            'hardware_config': hw_config
                        }
                        self.pending_workflow_steps.append(pending_info)
                        # Save to file for restoration
                        self._save_pending_workflow_step(pending_info)
                        successful_runs += 1
                        continue

                    # No dependencies - start immediately
                    success = self._run_amber_simulation(
                        sim_config=sim_config,
                        sim_dir=sim_dir,
                        shared_dir=shared_dir,
                        hardware_config=hw_config
                    )

                    if success:
                        successful_runs += 1
                        if sim_config.workflow_id:
                            self.console.print(f"[grey50]First step of protocol started.[/grey50]")
                    else:
                        failed_runs += 1
                        self.console.print(f"[red]✗ Failed to start: {sim_config.name}[/red]")

                        # Ask if we should continue with remaining simulations
                        if i < len(self.simulation_queue) - 1:  # Not the last simulation
                            continue_str = prompt_with_context(
                                self.processor,
                                "Continue with remaining simulations?",
                                choices=["y", "n"],
                                default="y",
                                module="MD Manager - Execution",
                                description="Continue after simulation failure",
                                options_map={"y": "Yes, continue", "n": "No, stop batch execution"}
                            )
                            if not (continue_str.lower() == "y"):
                                self.console.print("[yellow]Batch execution stopped by user[/yellow]")
                                break
                                
                except KeyboardInterrupt:
                    self.console.print("\n[yellow]Execution interrupted by user[/yellow]")
                    break
                except Exception as e:
                    failed_runs += 1
                    self.console.print(f"[red]✗ Error in {sim_config.name}: {e}[/red]")
                    continue
                    
            # Final summary
            self.console.print(f"\n[bold cyan]===== Batch Execution Started =====[/bold cyan]")
            self.console.print(f"[green]✓ Simulations started: {successful_runs}[/green]")
            if failed_runs > 0:
                self.console.print(f"[red]✗ Failed to start: {failed_runs}[/red]")

            # Show pending workflow steps
            if hasattr(self, 'pending_workflow_steps') and self.pending_workflow_steps:
                pending_count = len(self.pending_workflow_steps)
                self.console.print(f"[yellow]⏸️  Pending protocol steps: {pending_count}[/yellow]")
                self.console.print(f"[grey50]These will auto-start when dependencies complete[/grey50]")

            self.console.print(f"[blue]📁 Results directory: {run_dir}[/blue]")
            self.console.print(f"[yellow]💡 Use monitoring features to track simulation progress[/yellow]")

            # Start background workflow manager if there are pending steps
            if hasattr(self, 'pending_workflow_steps') and self.pending_workflow_steps:
                self.console.print(f"[yellow]💡 Protocol steps will auto-start in the background[/yellow]")
                self._start_workflow_background_manager()

        finally:
            # Always restore original working directory
            os.chdir(original_cwd)
            self.console.print(f"[grey50]Restored working directory: {original_cwd}[/grey50]")

        # Return to menu immediately for monitoring
        return True

    def _reap_finished_simulations(self):
        """Move finished simulations out of ``running_simulations``.

        A tracked simulation is "finished" once its process has exited — this
        includes zombie/dead states, since these detached children are never
        ``wait()``-ed on, so a bare ``os.kill(pid, 0)`` would keep reporting a
        finished-but-unreaped process as alive. Finished entries are moved into
        ``self.completed_simulations`` (tagged completed vs. failed by whether
        the restart file was written) so the running banner and monitor reflect
        reality while dependency auto-start and analysis can still see results.
        """
        if not getattr(self, 'running_simulations', None):
            return

        import psutil

        if not hasattr(self, 'completed_simulations'):
            self.completed_simulations = {}

        for name, info in list(self.running_simulations.items()):
            pid = info.get('pid')
            finished = False
            if pid is None:
                finished = True
            else:
                try:
                    proc = psutil.Process(pid)
                    status = proc.status()
                    if status in (psutil.STATUS_ZOMBIE, psutil.STATUS_DEAD) or not proc.is_running():
                        finished = True
                        try:
                            proc.wait(timeout=0)  # reap zombie so it stops lingering
                        except Exception:
                            pass
                except psutil.NoSuchProcess:
                    finished = True
                except Exception:
                    # Can't determine state — leave it tracked rather than guess
                    continue

            if finished:
                self.running_simulations.pop(name, None)
                restart_file = info.get('restart_file')
                info['ended_at'] = datetime.now()
                info['completed'] = bool(restart_file and Path(restart_file).exists())
                self.completed_simulations[name] = info

    def _display_workflow_status(self):
        """Display status bar showing running, finished, and pending simulations."""
        # Reap finished processes first so the banner reflects reality instead of
        # showing every launched simulation as perpetually "running".
        self._reap_finished_simulations()

        has_running = bool(getattr(self, 'running_simulations', None))
        has_pending = bool(getattr(self, 'pending_workflow_steps', None))
        has_completed = bool(getattr(self, 'completed_simulations', None))

        if not (has_running or has_pending or has_completed):
            return

        running_count = len(self.running_simulations) if has_running else 0
        pending_count = len(self.pending_workflow_steps) if has_pending else 0

        self.console.print()

        # Show running simulations
        if running_count > 0:
            for sim_name in self.running_simulations.keys():
                # Extract just the step name (remove structure prefix)
                step_name = sim_name.split('_', 3)[-1] if '_' in sim_name else sim_name
                self.console.print(f"[green]▶ Running:[/green] {step_name}")

        # Show finished simulations once nothing is actively running, so the user
        # gets explicit confirmation a step completed (and can go analyze it).
        if has_completed and running_count == 0:
            for sim_name, info in self.completed_simulations.items():
                step_name = sim_name.split('_', 3)[-1] if '_' in sim_name else sim_name
                if info.get('completed'):
                    self.console.print(f"[blue]✓ Finished:[/blue] {step_name}")
                else:
                    self.console.print(f"[red]✗ Ended (no restart file):[/red] {step_name}")

        # Show pending count
        if pending_count > 0:
            next_pending = self.pending_workflow_steps[0]['config'].name if self.pending_workflow_steps else ""
            step_name = next_pending.split('_', 3)[-1] if '_' in next_pending else next_pending
            self.console.print(f"[yellow]⏸  Next:[/yellow] {step_name} [grey50](+{pending_count - 1} more pending)[/grey50]" if pending_count > 1 else f"[yellow]⏸  Next:[/yellow] {step_name}")

    def _start_workflow_background_manager(self):
        """Start a background thread to monitor and auto-start workflow steps."""
        import threading

        # Check if background manager is already running
        if hasattr(self, '_workflow_manager_thread') and self._workflow_manager_thread and self._workflow_manager_thread.is_alive():
            return  # Already running

        # Create and start background thread
        self._workflow_manager_stop = threading.Event()
        self._workflow_manager_thread = threading.Thread(
            target=self._workflow_background_loop,
            daemon=True  # Daemon thread will exit when main program exits
        )
        self._workflow_manager_thread.start()
        self.console.print(f"[grey50]Started background protocol manager[/grey50]")

    def _workflow_background_loop(self):
        """Background loop that checks for completed dependencies and starts pending steps."""
        import time

        while not self._workflow_manager_stop.is_set():
            try:
                # Check if there are pending workflow steps
                if hasattr(self, 'pending_workflow_steps') and self.pending_workflow_steps:
                    self._check_and_start_pending_workflow_steps(silent=True)

                    # If no more pending steps, exit the loop
                    if not self.pending_workflow_steps:
                        break
                else:
                    # No pending steps, exit
                    break

                # Wait before next check (check every 10 seconds)
                time.sleep(10)

            except Exception as e:
                # Silently log errors - don't interrupt user's console
                break

    def _check_and_start_pending_workflow_steps(self, silent=False):
        """Check if any pending workflow steps can be started (dependencies completed).

        Args:
            silent: If True, suppress console output (for background thread)
        """
        if not hasattr(self, 'pending_workflow_steps') or not self.pending_workflow_steps:
            return

        import psutil

        # Reap finished processes first so a dependency that just completed is
        # recognized here (and surfaced in the running banner) consistently.
        self._reap_finished_simulations()

        still_pending = []

        for pending in self.pending_workflow_steps:
            sim_config = pending['config']
            depends_on = sim_config.depends_on

            # Check if dependency has completed
            dependency_completed = False

            # Fast path: reaping already moved a finished dependency here.
            completed = getattr(self, 'completed_simulations', {})
            if depends_on in completed:
                if completed[depends_on].get('completed'):
                    dependency_completed = True
                elif not silent:
                    restart_file = completed[depends_on].get('restart_file')
                    self.console.print(f"[red]✗ Dependency ended without restart file: {restart_file}[/red]")

            # Look for the dependency in running_simulations
            if not dependency_completed and hasattr(self, 'running_simulations'):
                if depends_on in self.running_simulations:
                    dep_info = self.running_simulations[depends_on]
                    pid = dep_info['pid']

                    try:
                        process = psutil.Process(pid)
                        # Check if process is zombie (finished but not reaped) or truly not running
                        process_status = process.status()
                        is_zombie = process_status == psutil.STATUS_ZOMBIE
                        is_dead = process_status == psutil.STATUS_DEAD

                        if is_zombie or is_dead or not process.is_running():
                            # Process completed (or zombie) - check for restart file
                            restart_file = dep_info['restart_file']
                            if restart_file.exists():
                                dependency_completed = True
                                # Remove from running_simulations
                                del self.running_simulations[depends_on]
                            else:
                                if not silent:
                                    self.console.print(f"[red]✗ Process completed but restart file not found: {restart_file}[/red]")
                    except psutil.NoSuchProcess:
                        # Process no longer exists - check for restart file
                        restart_file = dep_info['restart_file']
                        if restart_file.exists():
                            dependency_completed = True
                            # Remove from running_simulations
                            del self.running_simulations[depends_on]
                        else:
                            if not silent:
                                self.console.print(f"[red]✗ Restart file not found: {restart_file}[/red]")

            if dependency_completed:
                # Start this workflow step (silently in background)
                if not silent:
                    self.console.print(f"\n[bold blue]Auto-starting protocol step: {sim_config.name}[/bold blue]")
                try:
                    success = self._run_amber_simulation(
                        sim_config=sim_config,
                        sim_dir=pending['sim_dir'],
                        shared_dir=pending['shared_dir'],
                        hardware_config=pending['hardware_config'],
                        silent=silent
                    )
                    if success:
                        # Remove from persistence file since it's now running
                        self._remove_pending_workflow_step(sim_config.name)
                    elif not silent:
                        self.console.print(f"[red]✗ Failed to auto-start: {sim_config.name}[/red]")
                except Exception as e:
                    if not silent:
                        self.console.print(f"[red]✗ Error auto-starting {sim_config.name}: {e}[/red]")
            else:
                # Still pending
                still_pending.append(pending)

        # Update pending list
        self.pending_workflow_steps = still_pending

    def _wait_for_simulation_completion(self, sim_config: SimulationConfig, sim_dir: Path) -> bool:
        """
        Wait for a simulation to complete by monitoring the process ID.

        Returns:
            True if simulation completed successfully, False otherwise
        """
        import time
        import psutil

        # Get the PID from running_simulations
        if not hasattr(self, 'running_simulations') or sim_config.name not in self.running_simulations:
            self.console.print(f"[red]Cannot monitor: simulation not found in running list[/red]")
            return False

        sim_info = self.running_simulations[sim_config.name]
        pid = sim_info['pid']

        output_prefix = sim_config.name.replace(' ', '_').lower()
        restart_file = sim_dir / f"{output_prefix}.rst7"
        mdout_file = sim_dir / f"{output_prefix}.mdout"

        self.console.print(f"[grey50]Monitoring process PID: {pid}[/grey50]")

        try:
            process = psutil.Process(pid)
        except psutil.NoSuchProcess:
            self.console.print(f"[red]Process {pid} not found - may have exited immediately[/red]")
            # Check if output files exist to determine success
            if restart_file.exists():
                return True
            else:
                return False

        # Monitor process until it completes
        check_count = 0
        while True:
            try:
                # Check if process is still running
                if not process.is_running():
                    self.console.print(f"[grey50]Process {pid} has completed[/grey50]")

                    # Check exit code if available
                    try:
                        exit_code = process.wait(timeout=1)
                        if exit_code == 0:
                            self.console.print(f"[green]Process exited successfully (exit code: 0)[/green]")
                        else:
                            self.console.print(f"[yellow]Process exited with code: {exit_code}[/yellow]")
                    except:
                        pass

                    # Verify restart file was created
                    if restart_file.exists():
                        return True
                    else:
                        self.console.print(f"[yellow]Warning: Restart file not found after completion[/yellow]")
                        # Give it a moment in case of file system delay
                        time.sleep(2)
                        if restart_file.exists():
                            return True
                        else:
                            return False

                # Process still running - show periodic status
                check_count += 1
                if check_count % 10 == 0:  # Every 30 seconds (10 checks * 3 sec)
                    # Show file size as progress indicator
                    if mdout_file.exists():
                        size_mb = mdout_file.stat().st_size / (1024 * 1024)
                        self.console.print(f"[grey50]Still running... (output file: {size_mb:.1f} MB)[/grey50]")
                    else:
                        self.console.print(f"[grey50]Still running... (waiting for output)[/grey50]")

                # Wait before next check
                time.sleep(3)

            except psutil.NoSuchProcess:
                # Process disappeared - check for successful completion
                self.console.print(f"[grey50]Process completed[/grey50]")
                if restart_file.exists():
                    return True
                else:
                    return False

            except KeyboardInterrupt:
                self.console.print(f"\n[yellow]Monitoring interrupted by user[/yellow]")
                self.console.print(f"[yellow]Note: Process {pid} is still running in background[/yellow]")
                return False

            except Exception as e:
                self.console.print(f"[red]Error monitoring process: {e}[/red]")
                return False

    def _find_workflow_dependency_file(self, workflow_id: str, depends_on: str, run_dir: Path, silent: bool = False) -> Optional[Path]:
        """Find the restart file from a previous workflow step.

        Args:
            silent: If True, suppress console output
        """
        if not depends_on or not workflow_id:
            return None

        # Look for the previous step's output in _workflow_outputs
        dependency_key = f"{workflow_id}:{depends_on}"
        if dependency_key in self._workflow_outputs:
            output_file = self._workflow_outputs[dependency_key]
            if output_file.exists():
                return output_file

        # Fallback: Search for restart files in run directory
        # Look for directories that might contain the previous step
        depends_on_normalized = depends_on.replace(' ', '_').lower()
        
        # Try multiple matching strategies
        potential_dirs = []
        for sim_dir_path in run_dir.glob("*"):
            if sim_dir_path.is_dir():
                dir_name_lower = sim_dir_path.name.lower()
                # Strategy 1: Exact substring match
                if depends_on.lower() in dir_name_lower:
                    potential_dirs.append(sim_dir_path)
                # Strategy 2: Normalized name match (spaces → underscores)
                elif depends_on_normalized in dir_name_lower:
                    potential_dirs.append(sim_dir_path)
                # Strategy 3: Try matching key words from the step name
                elif any(word.lower() in dir_name_lower for word in depends_on.split() if len(word) > 3):
                    potential_dirs.append(sim_dir_path)
        
        # Look for restart files in potential directories
        for sim_dir_path in potential_dirs:
            # Look for output restart files (with lowercase underscored names from AMBER)
            # These have pattern like: stepname.rst7 or stepname.rst
            restart_files = list(sim_dir_path.glob("*.rst*"))

            # Sort to prioritize output files over input files
            # Output files typically have longer names with underscores
            restart_files.sort(key=lambda f: len(f.stem), reverse=True)

            for restart_file in restart_files:
                if restart_file.is_file():
                    # Only accept files that look like AMBER output restart files
                    # These should contain lowercase words from the step name
                    file_stem_lower = restart_file.stem.lower()
                    if '_' in file_stem_lower and any(word in file_stem_lower for word in ['minimization', 'heating', 'equilibration', 'production', 'md', 'npt', 'nvt']):
                        # This looks like an output restart file
                        if not silent:
                            self.console.print(f"[grey50]Found dependency file: {restart_file.name} in {sim_dir_path.name}[/grey50]")
                        return restart_file

        if not silent:
            self.console.print(f"[yellow]Warning: Could not find output from previous step '{depends_on}' for protocol '{workflow_id}'[/yellow]")
            self.console.print(f"[grey50]Searched in: {[d.name for d in potential_dirs]}[/grey50]")
        return None
    
    def _run_amber_simulation(self, sim_config: SimulationConfig, sim_dir: Path,
                            shared_dir: Path, hardware_config: Dict[str, Any], silent: bool = False) -> bool:
        """Run a single AMBER simulation with the specified configuration.

        Args:
            silent: If True, suppress console output (for background auto-start)
        """
        try:
            # Get original filenames from config
            prmtop_path = Path(sim_config.prmtop)
            coord_path = Path(sim_config.rst7)
            prmtop_name = prmtop_path.name
            coord_name = coord_path.name
            
            # Copy topology file to simulation directory
            import shutil
            shutil.copy2(shared_dir / prmtop_name, sim_dir / prmtop_name)
            
            # Handle coordinate file with workflow dependency logic
            actual_coord_file = None
            if sim_config.workflow_id and sim_config.depends_on:
                # This is a workflow step with dependencies - look for previous step's output
                run_dir = sim_dir.parent  # Go up one level to the main run directory
                dependency_file = self._find_workflow_dependency_file(
                    sim_config.workflow_id,
                    sim_config.depends_on,
                    run_dir,
                    silent=silent
                )
                
                if dependency_file:
                    # Use the previous step's restart file as input coordinates
                    actual_coord_file = sim_dir / coord_name
                    shutil.copy2(dependency_file, actual_coord_file)
                    if not silent:
                        self.console.print(f"[grey50]Using coordinates from previous step: {dependency_file.name}[/grey50]")
                else:
                    # Fallback to original coordinates with warning
                    actual_coord_file = sim_dir / coord_name
                    shutil.copy2(shared_dir / coord_name, actual_coord_file)
                    if not silent:
                        self.console.print(f"[yellow]Warning: Using initial coordinates (dependency not found)[/yellow]")
            else:
                # Standalone simulation or first step of workflow - use original coordinates
                actual_coord_file = sim_dir / coord_name
                shutil.copy2(shared_dir / coord_name, actual_coord_file)
                if not silent:
                    if sim_config.workflow_id:
                        self.console.print(f"[grey50]Using initial coordinates (first step of protocol)[/grey50]")
                    else:
                        self.console.print(f"[grey50]Using coordinates: {coord_name}[/grey50]")

            # Get or create mdin file
            mdin_file = sim_dir / "simulation.mdin"

            # Source the mdin text, then run EVERY source through the same
            # restraint pass below. An imported .mdin used to be copy2'd verbatim
            # here, which silently discarded restraints configured in Step 3 --
            # the restraint manager is the single authority, so reading the file
            # rather than copying it is what makes that true for custom inputs.
            source_label = None
            if sim_config.mdin_path and Path(sim_config.mdin_path).exists():
                try:
                    template_content = Path(sim_config.mdin_path).read_text()
                    source_label = f"mdin file: {sim_config.mdin_path}"
                except OSError as e:
                    self.console.print(
                        f"[yellow]Warning: could not read {sim_config.mdin_path} ({e}); "
                        f"copying it unmodified. Restraints will NOT be applied.[/yellow]")
                    shutil.copy2(sim_config.mdin_path, mdin_file)
                    template_content = None
                    source_label = "__copied__"
            else:
                # Get template content from template_id
                template_content = self._resolve_mdin_content(sim_config)
                source_label = f"template: {sim_config.template_id}"

            if source_label != "__copied__":
                if template_content:
                    self._warn_if_restraints_overwrite_input(
                        template_content, sim_config, silent)
                    # Apply restraints if configured (new unified restraint system)
                    if sim_config.restraints:
                        # Apply GROUP restraints (mutually exclusive with restraintmask)
                        if 'group' in sim_config.restraints:
                            template_content = self._apply_group_to_template(
                                template_content,
                                sim_config.restraints['group']
                            )
                            if not silent:
                                self.console.print(f"[grey50]Applied GROUP restraints to template[/grey50]")
                        # Apply positional restraints (restraintmask)
                        elif 'restraintmask' in sim_config.restraints:
                            template_content = self._apply_restraintmask_to_template(
                                template_content,
                                sim_config.restraints['restraintmask']
                            )
                            if not silent:
                                self.console.print(f"[grey50]Applied positional restraints to template[/grey50]")

                        # Apply DISANG restraints
                        if 'disang' in sim_config.restraints:
                            template_content = self._apply_disang_to_template(
                                template_content,
                                sim_config.restraints['disang'],
                                sim_dir
                            )
                            if not silent:
                                self.console.print(f"[grey50]Applied DISANG restraints to template[/grey50]")

                    # Inject CpHMD namelist variables for production steps
                    template_content = self._maybe_inject_cpmd_params(
                        template_content, sim_config, silent)

                    with open(mdin_file, 'w') as f:
                        f.write(template_content)
                    if not silent:
                        self.console.print(f"[grey50]Using {source_label}[/grey50]")
                else:
                    # Fallback to basic mdin file
                    self._create_basic_mdin_file(mdin_file)
                    if not silent:
                        self.console.print(f"[yellow]Warning: Created basic mdin file for {sim_config.name}[/yellow]")

            # Set up output file names based on simulation name
            output_prefix = sim_config.name.replace(' ', '_').lower()
            mdout_file = sim_dir / f"{output_prefix}.mdout"
            restart_file = sim_dir / f"{output_prefix}.rst7"
            trajectory_file = sim_dir / f"{output_prefix}.nc"

            # Build AMBER command - engine comes from sim_config, not hardware_config
            engine = sim_config.engine

            # Use input coordinates as reference for restraints
            reference_file = sim_dir / coord_name

            # Determine topology to use (may be swapped for explicit-solvent CpHMD)
            topology_to_use = sim_dir / prmtop_name

            # Handle CpHMD file setup for production steps
            cpmd_extra_flags = self._setup_cpmd_files(
                sim_config, sim_dir, shared_dir, prmtop_name, output_prefix, silent)
            if cpmd_extra_flags and cpmd_extra_flags.get('modified_prmtop'):
                topology_to_use = sim_dir / cpmd_extra_flags['modified_prmtop']

            cmd = self._build_amber_command_line(
                engine=engine,
                mdin_file=mdin_file,
                topology_file=topology_to_use,
                coordinate_file=sim_dir / coord_name,
                output_file=mdout_file,
                restart_file=restart_file,
                trajectory_file=trajectory_file,
                reference_file=reference_file,
                hardware_config=hardware_config,
                extra_flags=cpmd_extra_flags.get('flags') if cpmd_extra_flags else None,
            )

            # Create run script for documentation
            self._create_run_script(sim_dir, cmd, sim_config, hardware_config)

            # Show command being executed
            if not silent:
                self.console.print(f"\n[bold]Starting simulation in background:[/bold]")
                self.console.print(f"[cyan]{' '.join(cmd)}[/cyan]")
            
            # Execute AMBER asynchronously with proper detachment
            # Redirect output to log files to prevent pipe buffer issues
            stdout_log = sim_dir / "simulation.stdout"
            stderr_log = sim_dir / "simulation.stderr"
            
            with open(stdout_log, 'w') as stdout_file, open(stderr_log, 'w') as stderr_file:
                process = subprocess.Popen(
                    cmd,
                    cwd=sim_dir,
                    stdout=stdout_file,
                    stderr=stderr_file,
                    start_new_session=True  # This handles session detachment safely
                )
            
            # Store process info for monitoring
            if not hasattr(self, 'running_simulations'):
                self.running_simulations = {}
            
            self.running_simulations[sim_config.name] = {
                'sim_dir': sim_dir,
                'sim_config': sim_config,
                'hardware_config': hardware_config,
                'mdout_file': mdout_file,
                'restart_file': restart_file,
                'started_at': datetime.now(),
                'pid': process.pid,
                'stdout_log': stdout_log,
                'stderr_log': stderr_log
            }
            
            # Also save process info and config to file for persistence
            self._save_process_info(sim_config.name, process.pid, str(sim_dir), str(mdout_file))
            self._save_sim_config(sim_config, sim_dir)

            if not silent:
                self.console.print(f"[green]✓ Simulation '{sim_config.name}' started in background (PID: {process.pid})[/green]")
                self.console.print(f"[grey50]Simulation directory: {sim_dir}[/grey50]")
                self.console.print(f"[cyan]Tracking PID: {process.pid}[/cyan]")
                self.console.print(f"[grey50]Restart file will be: {restart_file}[/grey50]")
            
            # For now, return True since we started successfully
            # Actual completion will be checked by monitoring functions
            return True
                
        except Exception as e:
            self.console.print(f"[red]Error starting simulation: {e}[/red]")
            return False

    def _save_process_info(self, sim_name: str, pid: int, sim_dir: str, mdout_file: str):
        """Save process information to file for persistence across restarts."""
        try:
            import json
            process_file = Path.cwd() / "simulations" / ".running_processes.json"
            
            # Load existing processes
            processes = {}
            if process_file.exists():
                with open(process_file, 'r') as f:
                    processes = json.load(f)
            
            # Add new process
            processes[sim_name] = {
                'pid': pid,
                'sim_dir': sim_dir,
                'mdout_file': mdout_file,
                'started_at': datetime.now().isoformat()
            }
            
            # Save back to file
            process_file.parent.mkdir(exist_ok=True)
            with open(process_file, 'w') as f:
                json.dump(processes, f, indent=2)
                
        except Exception as e:
            # Don't fail simulation if we can't save process info
            self.console.print(f"[dark_orange3]Warning: Could not save process info: {e}[/dark_orange3]")

    def _save_sim_config(self, sim_config: SimulationConfig, sim_dir: Path):
        """Save simulation config to the simulation directory for restoration."""
        try:
            import json
            config_file = sim_dir / "config.json"

            config_data = {
                'name': sim_config.name,
                'prmtop': sim_config.prmtop,
                'rst7': sim_config.rst7,
                'template_id': sim_config.template_id,
                'mdin_path': sim_config.mdin_path,
                'engine': sim_config.engine,
                'workflow_id': getattr(sim_config, 'workflow_id', None),
                'workflow_step': getattr(sim_config, 'workflow_step', None),
                'depends_on': getattr(sim_config, 'depends_on', None)
            }

            with open(config_file, 'w') as f:
                json.dump(config_data, f, indent=2)

        except Exception as e:
            # Don't fail simulation if we can't save config
            pass

    def _save_pending_workflow_step(self, pending_info: Dict):
        """Save pending workflow step to file for restoration."""
        try:
            import json
            pending_file = Path.cwd() / "simulations" / ".pending_workflows.json"

            # Load existing pending steps
            pending_steps = []
            if pending_file.exists():
                with open(pending_file, 'r') as f:
                    pending_steps = json.load(f)

            # Add new pending step
            sim_config = pending_info['config']
            pending_steps.append({
                'name': sim_config.name,
                'sim_dir': str(pending_info['sim_dir']),
                'shared_dir': str(pending_info['shared_dir']),
                'hardware_config': pending_info['hardware_config'],
                'workflow_id': sim_config.workflow_id,
                'depends_on': sim_config.depends_on
            })

            # Save back to file
            pending_file.parent.mkdir(exist_ok=True)
            with open(pending_file, 'w') as f:
                json.dump(pending_steps, f, indent=2)

        except Exception as e:
            # Don't fail if we can't save
            pass

    def _remove_pending_workflow_step(self, sim_name: str):
        """Remove a pending workflow step from persistence file."""
        try:
            import json
            pending_file = Path.cwd() / "simulations" / ".pending_workflows.json"

            if not pending_file.exists():
                return

            with open(pending_file, 'r') as f:
                pending_steps = json.load(f)

            # Filter out the completed step
            pending_steps = [step for step in pending_steps if step['name'] != sim_name]

            # Save back
            with open(pending_file, 'w') as f:
                json.dump(pending_steps, f, indent=2)

        except Exception as e:
            # Don't fail if we can't remove
            pass

    def _restore_running_simulations(self):
        """Restore running simulations from .running_processes.json file."""
        try:
            import json
            import psutil

            process_file = Path.cwd() / "simulations" / ".running_processes.json"

            if not process_file.exists():
                return  # No saved processes

            with open(process_file, 'r') as f:
                saved_processes = json.load(f)

            if not saved_processes:
                return

            # Initialize running_simulations if needed
            if not hasattr(self, 'running_simulations'):
                self.running_simulations = {}

            restored_count = 0
            cleaned_count = 0
            processes_to_keep = {}

            for sim_name, process_info in saved_processes.items():
                pid = process_info.get('pid')
                sim_dir = Path(process_info.get('sim_dir'))
                mdout_file = Path(process_info.get('mdout_file'))

                # Check if process is still running
                try:
                    process = psutil.Process(pid)

                    # Check if it's a zombie or dead
                    if process.status() in [psutil.STATUS_ZOMBIE, psutil.STATUS_DEAD]:
                        cleaned_count += 1
                        continue

                    # Process is running - restore it
                    restart_file = sim_dir / f"{sim_name.replace(' ', '_').lower()}.rst7"
                    stdout_log = sim_dir / "simulation.stdout"
                    stderr_log = sim_dir / "simulation.stderr"

                    # Try to load sim config from the simulation directory
                    sim_config = self._load_sim_config_from_dir(sim_dir)

                    if sim_config:
                        self.running_simulations[sim_name] = {
                            'sim_dir': sim_dir,
                            'sim_config': sim_config,
                            'hardware_config': {},  # Will be unknown after restart
                            'mdout_file': mdout_file,
                            'restart_file': restart_file,
                            'started_at': datetime.fromisoformat(process_info.get('started_at', datetime.now().isoformat())),
                            'pid': pid,
                            'stdout_log': stdout_log,
                            'stderr_log': stderr_log,
                            'restored': True
                        }

                        # Add to pending workflow steps if it's part of a workflow
                        if hasattr(sim_config, 'workflow_id') and sim_config.workflow_id:
                            self._restore_workflow_pending_steps(sim_config)

                        processes_to_keep[sim_name] = process_info
                        restored_count += 1
                    else:
                        # Can't restore without config - keep in file for manual cleanup
                        processes_to_keep[sim_name] = process_info

                except psutil.NoSuchProcess:
                    # Process no longer exists - clean it up
                    cleaned_count += 1
                    continue
                except Exception as e:
                    # Any other error - skip this process
                    self.console.print(f"[dark_orange3]Warning: Could not restore process {sim_name}: {e}[/dark_orange3]")
                    continue

            # Update the process file with only the running processes
            if processes_to_keep != saved_processes:
                with open(process_file, 'w') as f:
                    json.dump(processes_to_keep, f, indent=2)

            # Restore pending workflow steps
            pending_count = self._restore_pending_workflow_steps()

            # Report restoration results
            if restored_count > 0:
                self.console.print(f"[green]✓ Restored {restored_count} running simulation(s)[/green]")
            if pending_count > 0:
                self.console.print(f"[yellow]⏸  Restored {pending_count} pending protocol step(s)[/yellow]")
            if cleaned_count > 0:
                self.console.print(f"[grey50]Cleaned up {cleaned_count} completed process(es)[/grey50]")

        except Exception as e:
            self.console.print(f"[dark_orange3]Warning: Could not restore running simulations: {e}[/dark_orange3]")

    def _load_sim_config_from_dir(self, sim_dir: Path) -> Optional[SimulationConfig]:
        """Load simulation config from a simulation directory."""
        try:
            # Look for config.json in simulation directory
            config_file = sim_dir / "config.json"
            if config_file.exists():
                import json
                with open(config_file, 'r') as f:
                    config_data = json.load(f)

                # Reconstruct SimulationConfig
                return SimulationConfig(
                    name=config_data.get('name', sim_dir.name),
                    prmtop=config_data.get('prmtop', ''),
                    rst7=config_data.get('rst7', ''),
                    template_id=config_data.get('template_id', ''),
                    mdin_path=config_data.get('mdin_path'),
                    engine=config_data.get('engine', 'sander'),
                    workflow_id=config_data.get('workflow_id'),
                    workflow_step=config_data.get('workflow_step'),
                    depends_on=config_data.get('depends_on')
                )

            # Fallback: Try to reconstruct from directory contents
            mdin_file = sim_dir / "simulation.mdin"
            if mdin_file.exists():
                # Look for topology and coordinate files
                prmtop_files = list(sim_dir.glob("*.prmtop"))
                coord_files = list(sim_dir.glob("*.rst7")) + list(sim_dir.glob("*.inpcrd"))

                if prmtop_files and coord_files:
                    return SimulationConfig(
                        name=sim_dir.name,
                        prmtop=str(prmtop_files[0]),
                        rst7=str(coord_files[0]),
                        template_id="unknown",
                        mdin_path=str(mdin_file),
                        engine='sander'
                    )

            return None

        except Exception as e:
            self.console.print(f"[dark_orange3]Warning: Could not load config from {sim_dir}: {e}[/dark_orange3]")
            return None

    def _restore_pending_workflow_steps(self) -> int:
        """Restore pending workflow steps from .pending_workflows.json file."""
        try:
            import json
            pending_file = Path.cwd() / "simulations" / ".pending_workflows.json"

            if not pending_file.exists():
                return 0

            with open(pending_file, 'r') as f:
                saved_pending = json.load(f)

            if not saved_pending:
                return 0

            # Initialize pending_workflow_steps if needed
            if not hasattr(self, 'pending_workflow_steps'):
                self.pending_workflow_steps = []

            restored_count = 0

            for pending_data in saved_pending:
                sim_name = pending_data['name']
                sim_dir = Path(pending_data['sim_dir'])
                shared_dir = Path(pending_data['shared_dir'])
                hardware_config = pending_data['hardware_config']
                workflow_id = pending_data['workflow_id']
                depends_on = pending_data['depends_on']

                # Load the simulation config
                sim_config = self._load_sim_config_from_dir(sim_dir)

                if sim_config:
                    # Check if already in pending list
                    already_pending = any(p['config'].name == sim_name for p in self.pending_workflow_steps)

                    if not already_pending:
                        self.pending_workflow_steps.append({
                            'config': sim_config,
                            'sim_dir': sim_dir,
                            'shared_dir': shared_dir,
                            'hardware_config': hardware_config
                        })
                        restored_count += 1

            return restored_count

        except Exception as e:
            # Silent failure - pending workflow restoration is best-effort
            return 0

    def _get_workflow_for_step(self, sim_config: SimulationConfig) -> Optional[WorkflowConfig]:
        """Look up the parent WorkflowConfig for a simulation step."""
        if not sim_config.workflow_id:
            return None
        return self.simulation_queue._workflows.get(sim_config.workflow_id)

    def _is_production_step(self, sim_config: SimulationConfig) -> bool:
        """Return True if this step is a production simulation."""
        if sim_config.simulation_type:
            return sim_config.simulation_type == 'production'
        # Fall back to template_id heuristic
        return self._get_template_type(sim_config.template_id) == 'production'

    def _maybe_inject_cpmd_params(self, mdin_content: str, sim_config: SimulationConfig,
                                   silent: bool = False) -> str:
        """Inject CpHMD namelist variables into mdin content for production steps.

        Only modifies the content if the parent workflow has cpmd_settings
        and this step is a production step.
        """
        if not self._is_production_step(sim_config):
            return mdin_content

        workflow = self._get_workflow_for_step(sim_config)
        if not workflow or not workflow.cpmd_settings:
            return mdin_content

        # Inject each CpHMD parameter as a parameter override
        overrides = {k: v for k, v in workflow.cpmd_settings.items()}
        mdin_content = self._apply_parameter_overrides_to_content(mdin_content, overrides)

        if not silent:
            self.console.print(f"[grey50]Injected CpHMD parameters: "
                             f"icnstph={overrides.get('icnstph')}, "
                             f"solvph={overrides.get('solvph')}[/grey50]")
        return mdin_content

    def _setup_cpmd_files(self, sim_config: SimulationConfig, sim_dir: Path,
                          shared_dir: Path, prmtop_name: str,
                          output_prefix: str, silent: bool = False) -> Optional[Dict]:
        """Set up CpHMD files for a production step.

        Copies the CPIN file, resolves restart chaining between consecutive
        production steps, and swaps the topology for explicit-solvent CpHMD.

        Returns a dict with 'flags' (list of command-line args) and optionally
        'modified_prmtop' (filename to use instead of the original), or None
        if CpHMD is not active for this step.
        """
        import shutil

        if not self._is_production_step(sim_config):
            return None

        workflow = self._get_workflow_for_step(sim_config)
        if not workflow or not workflow.cpin_file:
            return None

        cpin_source = Path(workflow.cpin_file)
        if not cpin_source.exists():
            if not silent:
                self.console.print(f"[yellow]Warning: CPIN file not found: {cpin_source}[/yellow]")
            return None

        # Determine CPIN input: original file or previous production step's cprestrt
        cpin_input_name = cpin_source.name
        cprestrt_from_prev = self._find_previous_cprestrt(sim_config, sim_dir)
        if cprestrt_from_prev and cprestrt_from_prev.exists():
            # Use previous step's cprestrt as our cpin
            shutil.copy2(cprestrt_from_prev, sim_dir / cprestrt_from_prev.name)
            cpin_input_name = cprestrt_from_prev.name
            if not silent:
                self.console.print(f"[grey50]CpHMD restart: using {cprestrt_from_prev.name} as cpin[/grey50]")
        else:
            # First production step: copy original CPIN file
            shutil.copy2(cpin_source, sim_dir / cpin_input_name)
            if not silent:
                self.console.print(f"[grey50]CpHMD: copied {cpin_input_name} to simulation directory[/grey50]")

        # Output files
        cpout_name = f"{output_prefix}.cpout"
        cprestrt_name = f"{output_prefix}.cprestrt"

        flags = [
            "-cpin", cpin_input_name,
            "-cpout", cpout_name,
            "-cprestrt", cprestrt_name,
        ]

        result = {'flags': flags}

        # Handle modified topology for explicit solvent CpHMD
        cpin_config = workflow.cpin_config or {}
        modified_prmtop_path = cpin_config.get('modified_prmtop')
        if modified_prmtop_path and Path(modified_prmtop_path).exists():
            mod_prmtop_name = Path(modified_prmtop_path).name
            shutil.copy2(modified_prmtop_path, sim_dir / mod_prmtop_name)
            result['modified_prmtop'] = mod_prmtop_name
            if not silent:
                self.console.print(f"[grey50]CpHMD: using modified topology {mod_prmtop_name}[/grey50]")

        return result

    def _find_previous_cprestrt(self, sim_config: SimulationConfig,
                                 sim_dir: Path) -> Optional[Path]:
        """Find the cprestrt file from the previous production step in this workflow."""
        workflow = self._get_workflow_for_step(sim_config)
        if not workflow or sim_config.workflow_step is None:
            return None

        # Walk backward through workflow steps to find the most recent production step
        current_step_idx = sim_config.workflow_step
        for prev_step in reversed(workflow.steps[:current_step_idx]):
            if self._is_production_step(prev_step):
                # Look for cprestrt in the previous step's directory
                prev_prefix = prev_step.name.replace(' ', '_').lower()
                # The previous step's directory is a sibling of our sim_dir
                run_dir = sim_dir.parent
                prev_dir = run_dir / prev_step.name.replace(' ', '_').lower()
                cprestrt = prev_dir / f"{prev_prefix}.cprestrt"
                if cprestrt.exists():
                    return cprestrt
                # Also check the pattern used by _run_amber_simulation
                for candidate in prev_dir.glob("*.cprestrt"):
                    return candidate
                break  # Only check the immediately preceding production step

        return None

    def _build_amber_command_line(self, engine: str, mdin_file: Path, topology_file: Path,
                                coordinate_file: Path, output_file: Path, restart_file: Path,
                                trajectory_file: Path, reference_file: Path,
                                hardware_config: Dict[str, Any],
                                extra_flags: Optional[List[str]] = None) -> List[str]:
        """Build AMBER command line based on engine and hardware configuration.

        extra_flags: additional command-line flags (e.g. CpHMD -cpin/-cpout/-cprestrt).
        """

        # Add info file for pmemd
        info_file = output_file.with_suffix('.mdinfo')

        # Use relative paths (just filenames) since all files are in the same directory
        base_args = [
            "-O",  # Overwrite output files
            "-i", mdin_file.name,
            "-p", topology_file.name,
            "-c", coordinate_file.name,
            "-o", output_file.name,
            "-r", restart_file.name,
            "-x", trajectory_file.name,
            "-ref", reference_file.name  # Reference coordinates for restraints
        ]

        # Add -inf flag for pmemd variants
        if engine.startswith("pmemd"):
            base_args.extend(["-inf", info_file.name])

        # Append any extra flags (e.g. CpHMD files)
        if extra_flags:
            base_args.extend(extra_flags)
        
        if engine == "sander":
            return ["sander"] + base_args
        elif engine == "pmemd":
            return ["pmemd"] + base_args
        elif engine == "pmemd.MPI":
            mpi_tasks = hardware_config.get('mpi_tasks', 4)
            return ["mpirun", "-np", str(mpi_tasks), "pmemd.MPI"] + base_args
        elif engine == "pmemd.cuda":
            # Set GPU environment variable if specified
            gpu_ids = hardware_config.get('gpu_ids')
            if gpu_ids:
                import os
                os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_ids)
            return ["pmemd.cuda"] + base_args
        else:
            # Default to sander
            return ["sander"] + base_args

    def _resolve_mdin_content(self, sim_config) -> Optional[str]:
        """Resolve MDIN content for a simulation config using the protocol-centric chain.

        Resolution order:
        1. mdin_content_override — user edited via nano, stored on protocol step
        2. template_id (builtin path) + parameter_overrides — read builtin, apply overrides
        3. parameter_overrides alone — wizard-created step, generate from config
        4. Legacy fallback — follow custom template UUID chain
        """
        # 1. Direct MDIN content override (user edited via nano)
        if sim_config.mdin_content_override:
            return sim_config.mdin_content_override

        template_id = sim_config.template_id or ''

        # 2. Builtin template path + optional parameter overrides
        if template_id.startswith('builtin/') or template_id.startswith('builtin\\'):
            try:
                from proprep.md_prep.user_data_manager import UserDataManager
                user_data_manager = UserDataManager(console=self.console)
                content, _ = user_data_manager.get_template_content(template_id)
                if content and sim_config.parameter_overrides:
                    content = self._apply_parameter_overrides_to_content(
                        content, sim_config.parameter_overrides
                    )
                return content
            except Exception as e:
                self.console.print(f"[yellow]Warning: Could not load builtin template {template_id}: {e}[/yellow]")

        # 3. No template ref but has parameter_overrides — wizard-created step
        if not template_id and sim_config.parameter_overrides:
            return self._generate_mdin_from_overrides(sim_config.parameter_overrides)

        # 4. Legacy fallback — existing chain-following logic
        if template_id:
            return self._get_template_content_from_id(template_id)

        return None

    def _apply_parameter_overrides_to_content(self, content: str, overrides: Dict) -> str:
        """Apply parameter overrides to MDIN content by modifying &cntrl values."""
        if not overrides or not content:
            return content

        lines = content.split('\n')
        new_lines = []
        in_cntrl = False
        applied = set()

        for line in lines:
            stripped = line.strip()
            if stripped.startswith('&cntrl'):
                in_cntrl = True
                new_lines.append(line)
                continue
            if in_cntrl and stripped == '/':
                # Insert any remaining overrides before closing
                for key, value in overrides.items():
                    if key not in applied:
                        if isinstance(value, str) and ' ' in value:
                            new_lines.append(f"  {key}='{value}',")
                        else:
                            new_lines.append(f"  {key}={value},")
                        applied.add(key)
                in_cntrl = False
                new_lines.append(line)
                continue
            if in_cntrl:
                # Check if this line contains a parameter we need to override
                modified = False
                for key, value in overrides.items():
                    if key in applied:
                        continue
                    # Match parameter assignment patterns like "  imin=1,"
                    import re
                    pattern = rf'(\s*){key}\s*=\s*[^,/\n]+(,?)'
                    match = re.search(pattern, stripped)
                    if match:
                        if isinstance(value, str) and ' ' in value:
                            new_val = f"  {key}='{value}',"
                        else:
                            new_val = f"  {key}={value},"
                        new_lines.append(new_val)
                        applied.add(key)
                        modified = True
                        break
                if not modified:
                    new_lines.append(line)
            else:
                new_lines.append(line)

        return '\n'.join(new_lines)

    def _generate_mdin_from_overrides(self, overrides: Dict) -> str:
        """Generate MDIN content from parameter overrides (wizard-created steps)."""
        description = overrides.pop('_description', 'Custom simulation step') if '_description' in overrides else 'Custom simulation step'
        nmr_section = overrides.pop('_nmr_section', '') if '_nmr_section' in overrides else ''

        lines = [description, "&cntrl"]
        for key in sorted(overrides.keys()):
            if key.startswith('_'):
                continue
            value = overrides[key]
            if isinstance(value, str) and ' ' in value:
                lines.append(f"  {key}='{value}',")
            else:
                lines.append(f"  {key}={value},")
        lines.append("/")

        if nmr_section:
            lines.append("")
            lines.extend(nmr_section.split('\n'))

        return '\n'.join(lines) + '\n'

    def _get_template_content_from_id(self, template_id: str) -> Optional[str]:
        """Get template content from template ID, following chain of template references."""
        try:
            from proprep.md_prep.user_data_manager import UserDataManager

            user_data_manager = UserDataManager(console=self.console)

            # Handle different template ID formats
            if template_id.startswith("builtin"):
                # Builtin template path
                content, metadata = user_data_manager.get_template_content(template_id)
                return content
            else:
                # Custom template ID (UUID) - follow the chain
                visited = set()  # Prevent infinite loops
                current_id = template_id

                while current_id and current_id not in visited:
                    visited.add(current_id)

                    template_data = user_data_manager.load_custom_template(current_id)
                    if not template_data:
                        break

                    # Check if we have a template reference
                    if 'template' in template_data:
                        template_ref = template_data['template']

                        # Is it a builtin template path?
                        if template_ref.startswith("builtin"):
                            try:
                                content, _ = user_data_manager.get_template_content(template_ref)
                                return content
                            except:
                                pass
                        else:
                            # It's another custom template UUID - follow the chain
                            current_id = template_ref
                            continue

                    # No template reference or reached end of chain - generate from config
                    return self._generate_mdin_from_template_data(template_data)

        except Exception as e:
            self.console.print(f"[yellow]Warning: Could not load template {template_id}: {e}[/yellow]")
            return None

    def _get_template_type(self, template_id: str) -> str:
        """Get template type (minimization, heating, equilibration, production) from template ID."""
        try:
            from proprep.md_prep.user_data_manager import UserDataManager

            user_data_manager = UserDataManager(console=self.console)

            # Handle different template ID formats
            if template_id.startswith("builtin"):
                # Determine type from filename keywords
                if "minimization" in template_id:
                    return "minimization"
                elif "heating" in template_id:
                    return "heating"
                elif "equilibration" in template_id:
                    return "equilibration"
                elif "production" in template_id:
                    return "production"
                return "step"
            else:
                # Custom template ID (UUID) - load and get type
                template_data = user_data_manager.load_custom_template(template_id)
                if template_data and 'type' in template_data:
                    return template_data['type']
                return "step"

        except Exception as e:
            return "step"

    def _generate_mdin_from_template_data(self, template_data: Dict) -> str:
        """Generate mdin content from template data."""
        config = template_data.get('config', {})
        description = template_data.get('description', 'Custom template')
        nmr_section = template_data.get('nmr_section', '')
        
        lines = [
            description,
            "&cntrl"
        ]
        
        # Add parameters in sorted order for consistency
        for key in sorted(config.keys()):
            value = config[key]
            if isinstance(value, str) and " " in value:
                lines.append(f"  {key}='{value}',")
            else:
                lines.append(f"  {key}={value},")
                
        lines.append("/")
        
        # Add NMR section if present
        if nmr_section:
            lines.append("")
            lines.extend(nmr_section.split('\n'))
            
        return "\n".join(lines) + "\n"
    
    def _apply_redox_restraint_mask(self, template_content: str) -> str:
        """
        Apply redox site restraint mask to template content if available.
        
        Args:
            template_content: Original template content
            
        Returns:
            Modified template content with redox restraint mask applied
        """
        if not self.workspace:
            return template_content
            
        # Check for redox restraint mask in workspace
        redox_mask = self.workspace.get("redox_restraint_mask")
        if not redox_mask:
            return template_content
            
        # Apply the redox restraint mask to restraintmask parameter
        lines = template_content.split('\n')
        modified_lines = []
        
        for line in lines:
            # Look for restraintmask parameter
            if 'restraintmask' in line.lower() and '=' in line:
                # Replace the restraint mask value
                parts = line.split('=', 1)
                if len(parts) == 2:
                    # Keep the parameter name and comments, replace the value
                    prefix = parts[0] + '='
                    # Extract any inline comments
                    if '!' in parts[1]:
                        comment_parts = parts[1].split('!', 1)
                        comment = '!' + comment_parts[1]
                    else:
                        comment = ''
                    
                    # Create new line with redox restraint mask
                    new_line = f"{prefix}'{redox_mask}',  {comment}".rstrip()
                    modified_lines.append(new_line)
                    from rich.markup import escape
                    self.console.print(f"[grey50]Applied redox restraint mask: {escape(redox_mask)}[/grey50]")
                else:
                    modified_lines.append(line)
            else:
                modified_lines.append(line)
        
        return '\n'.join(modified_lines)

    # Restraint keywords the restraint manager owns. If an incoming mdin already
    # sets one of these AND restraints are configured, the manager's value wins --
    # which is the point of centralising, but it means editing input the user
    # hand-wrote, so it gets said out loud rather than done quietly.
    _RESTRAINT_OWNED_KEYWORDS = ("ntr", "restraintmask", "restraint_wt")

    def _warn_if_restraints_overwrite_input(self, mdin_text: str, sim_config,
                                            silent: bool = False) -> list:
        """Report restraint keywords in `mdin_text` that the restraint manager will replace.

        Returns the list of overridden keyword names (empty when nothing clashes),
        so callers and tests can act on it without re-parsing.
        """
        import re

        restraints = getattr(sim_config, 'restraints', None)
        if not restraints or not mdin_text:
            return []
        # DISANG layers on top via nmropt/DISANG files; it does not replace these.
        if not ({'restraintmask', 'group'} & set(restraints)):
            return []

        present = []
        for line in mdin_text.splitlines():
            body = line.split('!')[0]
            for kw in self._RESTRAINT_OWNED_KEYWORDS:
                if kw in present:
                    continue
                if re.search(rf"(?<![A-Za-z0-9_]){kw}\s*=", body, re.IGNORECASE):
                    present.append(kw)
        if present and not silent:
            self.console.print(
                f"[yellow]Note: this input already sets {', '.join(present)}; "
                f"the restraint manager's settings replace them.[/yellow]"
            )
        return present

    def _apply_configured_restraints(self, template_content: str, sim_config, working_dir: Path) -> str:
        """Apply any restraints configured on sim_config to the template.

        GROUP and restraintmask are mutually exclusive; DISANG layers on top
        of either (it writes its own auxiliary file into working_dir). This
        mirrors the logic the batch path has used from day one; kept here
        so the SLURM paths can share it without drifting.
        """
        if not getattr(sim_config, 'restraints', None):
            return template_content
        restraints = sim_config.restraints
        if 'group' in restraints:
            template_content = self._apply_group_to_template(
                template_content, restraints['group']
            )
        elif 'restraintmask' in restraints:
            template_content = self._apply_restraintmask_to_template(
                template_content, restraints['restraintmask']
            )
        if 'disang' in restraints:
            template_content = self._apply_disang_to_template(
                template_content, restraints['disang'], working_dir
            )
        return template_content

    def _apply_restraintmask_to_template(self, template_content: str, restraint_config: dict) -> str:
        """
        Apply restraintmask and restraint_wt to template content.

        Adds or updates:
        1. ntr=1 (enable positional restraints)
        2. restraint_wt (force constant)
        3. restraintmask (atom selection)

        Args:
            template_content: Original template content
            restraint_config: Dict with 'mask' and 'weight' keys

        Returns:
            Modified template content with restraintmask and restraint_wt applied
        """
        # Handle both old string format and new dict format for backward compatibility
        if isinstance(restraint_config, str):
            # Old format - just a mask string, use default weight
            restraintmask = restraint_config
            restraint_wt = 10.0
        else:
            # New format - dict with mask and weight
            restraintmask = restraint_config['mask']
            restraint_wt = restraint_config['weight']

        lines = template_content.split('\n')
        modified_lines = []
        in_cntrl_section = False
        ntr_found = False
        restraint_wt_found = False
        restraintmask_found = False

        for line in lines:
            stripped = line.strip().lower()

            # Skip misleading "no restraint" comments when we're applying restraints
            if in_cntrl_section and stripped.startswith('!') and 'no restraint' in stripped:
                continue  # Skip this line

            # Track if we're in &cntrl section
            if stripped.startswith('&cntrl'):
                in_cntrl_section = True
                modified_lines.append(line)
            elif stripped == '/' and in_cntrl_section:
                # End of &cntrl section - add missing parameters before closing
                if not ntr_found:
                    modified_lines.append(f"  ntr=1,  ! Enable positional restraints")
                if not restraint_wt_found:
                    modified_lines.append(f"  restraint_wt={restraint_wt},  ! Restraint force constant (kcal/mol/Å²)")
                if not restraintmask_found:
                    modified_lines.append(f"  restraintmask='{restraintmask}',  ! Atoms to restrain")
                in_cntrl_section = False
                modified_lines.append(line)
            # Look for ntr parameter (must be the actual parameter name, not substring)
            elif in_cntrl_section and '=' in stripped and stripped.split('=')[0].strip() == 'ntr':
                # Replace with ntr=1
                parts = line.split('=', 1)
                if len(parts) == 2:
                    prefix = parts[0] + '='
                    if '!' in parts[1]:
                        comment_parts = parts[1].split('!', 1)
                        comment = '!' + comment_parts[1]
                    else:
                        comment = ''
                    new_line = f"{prefix}1,  {comment}".rstrip()
                    modified_lines.append(new_line)
                    ntr_found = True
                else:
                    modified_lines.append(line)
            # Look for restraint_wt parameter
            elif 'restraint_wt' in stripped and '=' in stripped and in_cntrl_section:
                # Replace the restraint weight value
                parts = line.split('=', 1)
                if len(parts) == 2:
                    prefix = parts[0] + '='
                    if '!' in parts[1]:
                        comment_parts = parts[1].split('!', 1)
                        comment = '!' + comment_parts[1]
                    else:
                        comment = ''
                    new_line = f"{prefix}{restraint_wt},  {comment}".rstrip()
                    modified_lines.append(new_line)
                    restraint_wt_found = True
                else:
                    modified_lines.append(line)
            # Look for restraintmask parameter
            elif 'restraintmask' in stripped and '=' in stripped and in_cntrl_section:
                # Replace the restraint mask value
                parts = line.split('=', 1)
                if len(parts) == 2:
                    prefix = parts[0] + '='
                    # Find comment after the value (look for ! after closing quote)
                    # The value is typically 'value', or "value",
                    rest = parts[1]
                    comment = ''

                    # Find the closing quote (either ' or ")
                    if "'" in rest:
                        # Find position after closing quote
                        first_quote = rest.index("'")
                        after_first = rest[first_quote + 1:]
                        if "'" in after_first:
                            closing_quote_pos = after_first.index("'") + first_quote + 2
                            after_value = rest[closing_quote_pos:]
                            if '!' in after_value:
                                comment = '!' + after_value.split('!', 1)[1]
                    elif '"' in rest:
                        first_quote = rest.index('"')
                        after_first = rest[first_quote + 1:]
                        if '"' in after_first:
                            closing_quote_pos = after_first.index('"') + first_quote + 2
                            after_value = rest[closing_quote_pos:]
                            if '!' in after_value:
                                comment = '!' + after_value.split('!', 1)[1]

                    new_line = f"{prefix}'{restraintmask}',  {comment}".rstrip()
                    modified_lines.append(new_line)
                    restraintmask_found = True
                else:
                    modified_lines.append(line)
            else:
                modified_lines.append(line)

        return '\n'.join(modified_lines)

    def _format_group_specification(self, group_data: list) -> str:
        """
        Convert GROUP restraint dicts to mdin GROUP text block.

        Args:
            group_data: List of group dicts, each with 'title', 'force_constant',
                        'find_criteria' (list of strings), and 'residue_ranges'
                        (list of (start, end) tuples).

        Returns:
            Formatted GROUP specification text to append after the &cntrl namelist.
        """
        lines = []
        for group in group_data:
            lines.append(group.get('title', 'Restraint group'))
            lines.append(str(group.get('force_constant', 10.0)))
            lines.append("FIND")
            for criterion in group.get('find_criteria', []):
                lines.append(criterion)
            lines.append("SEARCH")
            # RES line: pairs of start end values
            res_values = []
            for start, end in group.get('residue_ranges', []):
                res_values.extend([str(start), str(end)])
            lines.append("RES " + " ".join(res_values))
            lines.append("END")
        lines.append("END")
        return "\n".join(lines)

    def _apply_group_to_template(self, template_content: str, group_data: list) -> str:
        """
        Apply GROUP restraint specification to template content.

        Replaces restraintmask/restraint_wt with ntr=1 and appends the GROUP
        specification block after the &cntrl namelist.

        Args:
            template_content: Original template content
            group_data: List of group dicts (see _format_group_specification)

        Returns:
            Modified template content with GROUP specification appended.
        """
        lines = template_content.split('\n')
        modified_lines = []
        in_cntrl_section = False
        ntr_found = False

        for line in lines:
            stripped = line.strip().lower()

            if stripped.startswith('&cntrl'):
                in_cntrl_section = True
                modified_lines.append(line)
            elif stripped == '/' and in_cntrl_section:
                # End of &cntrl — ensure ntr=1 is present, then close
                if not ntr_found:
                    modified_lines.append("  ntr=1,  ! Enable positional restraints (GROUP)")
                in_cntrl_section = False
                modified_lines.append(line)
                # Append GROUP specification immediately after &cntrl closing
                modified_lines.append(self._format_group_specification(group_data))
            elif in_cntrl_section:
                param_name = stripped.split('=')[0].strip() if '=' in stripped else ''
                # Remove restraintmask and restraint_wt — GROUP replaces them
                if param_name in ('restraintmask', 'restraint_wt'):
                    continue
                if param_name == 'ntr':
                    # Ensure ntr=1
                    parts = line.split('=', 1)
                    if len(parts) == 2:
                        prefix = parts[0] + '='
                        comment = ''
                        if '!' in parts[1]:
                            comment = '!' + parts[1].split('!', 1)[1]
                        modified_lines.append(f"{prefix}1,  {comment}".rstrip())
                    else:
                        modified_lines.append(line)
                    ntr_found = True
                else:
                    modified_lines.append(line)
            else:
                modified_lines.append(line)

        return '\n'.join(modified_lines)

    def _apply_disang_to_template(self, template_content: str, disang_config: dict, sim_dir: Path) -> str:
        """
        Apply DISANG restraint integration to template content.

        Adds:
        1. nmropt=1 to &cntrl section
        2. DUMPFREQ &wt blocks to nmr_section
        3. File redirection lines (DISANG, DUMPAVE, LISTOUT)

        Args:
            template_content: Original template content
            disang_config: DISANG configuration dictionary (contains file, dumpave_file, listout_file, dump_freq)
            sim_dir: Simulation directory for output files

        Returns:
            Modified template content with DISANG integration
        """
        lines = template_content.split('\n')
        modified_lines = []
        in_cntrl_section = False
        nmropt_added = False
        nmr_section_found = False

        # Step 1: Process existing content and add nmropt=1
        for line in lines:
            stripped = line.strip().lower()

            # Track if we're in &cntrl section
            if stripped.startswith('&cntrl'):
                in_cntrl_section = True
                modified_lines.append(line)
            elif stripped == '/' and in_cntrl_section:
                # End of &cntrl section - add nmropt=1 if not already present
                if not nmropt_added:
                    modified_lines.append(f"  nmropt=1,  ! Enable NMR restraints")
                    nmropt_added = True
                in_cntrl_section = False
                modified_lines.append(line)
            elif 'nmropt' in stripped and '=' in stripped and in_cntrl_section:
                # Replace existing nmropt with nmropt=1
                parts = line.split('=', 1)
                if len(parts) == 2:
                    prefix = parts[0]
                    comment = '! Enable NMR restraints'
                    if '!' in parts[1]:
                        comment_parts = parts[1].split('!', 1)
                        comment = '!' + comment_parts[1]
                    new_line = f"{prefix}=1,  {comment}".rstrip()
                    modified_lines.append(new_line)
                    nmropt_added = True
                else:
                    modified_lines.append(line)
            elif stripped.startswith('disang=') or stripped.startswith('dumpave=') or \
                 stripped.startswith('listout=') or stripped.startswith('listin='):
                # Remove existing file redirections — will be re-added below
                continue
            elif stripped.startswith('! disang file redirection'):
                # Remove old section header too
                continue
            elif stripped.startswith('! varying conditions') or '&wt' in stripped:
                # Found existing NMR section
                nmr_section_found = True
                modified_lines.append(line)
            else:
                modified_lines.append(line)

        # Step 2: Add DUMPFREQ &wt blocks and file redirections
        dump_freq = disang_config.get('dump_freq', 500)

        # Create file paths relative to sim_dir
        disang_file = disang_config['file']
        dumpave_file = sim_dir / disang_config['dumpave_file']
        listout_file = sim_dir / disang_config['listout_file']

        # nmropt reads each redirection line into a FIXED 80-char buffer, so
        # 'KEYWORD=' + path must fit in 80 chars or the filename is silently
        # truncated (a deep sim_dir clips DUMPAVE/LISTOUT to the same stem, and
        # they then collide). sander is always launched with cwd=sim_dir (both
        # the batch run_workflow.sh and the SLURM wrappers cd into stepN/), so
        # emit these RELATIVE to sim_dir to keep the lines short. DUMPAVE/LISTOUT
        # collapse to their bare basenames; DISANG becomes a short '../…' hop.
        disang_file = os.path.relpath(str(disang_file), str(sim_dir))
        dumpave_file = os.path.relpath(str(dumpave_file), str(sim_dir))
        listout_file = os.path.relpath(str(listout_file), str(sim_dir))

        # Check if DUMPFREQ already exists in the content
        has_dumpfreq = any('dumpfreq' in line.strip().lower() for line in modified_lines)

        # Add &wt blocks for DUMPFREQ
        if nmr_section_found:
            if not has_dumpfreq:
                # Insert DUMPFREQ before existing &wt type='END'
                final_lines = []
                for line in modified_lines:
                    if "&wt type='END'" in line or '&wt TYPE="END"' in line:
                        # Add DUMPFREQ before END
                        final_lines.append(f"&wt type='DUMPFREQ', istep1={dump_freq}, &end")
                        final_lines.append(line)
                    else:
                        final_lines.append(line)
                modified_lines = final_lines
        else:
            # Add complete nmr section
            modified_lines.append("")
            modified_lines.append("! Varying conditions (&wt blocks)")
            modified_lines.append(f"&wt type='DUMPFREQ', istep1={dump_freq}, &end")
            modified_lines.append("&wt type='END', &end")

        # Step 3: Add file redirection section.
        # NOTE: the restraint-file redirection block (after &wt type='END') is
        # parsed line-by-line by nmropt's own reader, which — unlike the &cntrl
        # namelist — has NO comment support: a leading '!' line is read as a
        # malformed KEYWORD=file redirection and sander aborts with
        # 'Missing "=" and/or filename after keyword'. So emit ONLY the
        # KEYWORD=path lines here (a blank separator line is tolerated).
        modified_lines.append("")
        modified_lines.append(f"DISANG={disang_file}")
        modified_lines.append(f"DUMPAVE={dumpave_file}")
        modified_lines.append(f"LISTOUT={listout_file}")

        return '\n'.join(modified_lines)

    def _create_workflow_run_script(self, workflow_dir: Path, workflow_steps: List,
                                    mdin_files: List[str], prmtop_name: str, coord_name: str,
                                    extended_production_cycles: int = 0):
        """Create a master bash script that chains all workflow steps together.

        If the parent workflow has CpHMD enabled, production steps get
        -cpin/-cpout/-cprestrt flags with proper restart chaining.
        """
        # Look up CpHMD configuration from the parent workflow
        wf_config = None
        if workflow_steps and workflow_steps[0].workflow_id:
            wf_config = self.simulation_queue._workflows.get(workflow_steps[0].workflow_id)
        cpmd_active = wf_config and wf_config.cpin_file
        cpin_filename = Path(wf_config.cpin_file).name if cpmd_active else None
        # Use modified topology for explicit-solvent CpHMD
        cpmd_prmtop = None
        if cpmd_active and wf_config.cpin_config:
            mod = wf_config.cpin_config.get('modified_prmtop')
            if mod and Path(mod).exists():
                cpmd_prmtop = Path(mod).name
        script_file = workflow_dir / "run_workflow.sh"

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(f"# AMBER Protocol Run Script\n")
            f.write(f"# Generated by ProPrep Molecular Dynamics Manager\n")
            f.write(f"# Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"#\n")
            f.write(f"# Protocol: {len(workflow_steps)} sequential steps\n")

            if extended_production_cycles > 0:
                f.write(f"# Extended Production: {extended_production_cycles} additional cycles\n")

            # List engines used (may vary by step)
            engines_used = set(step.engine for step in workflow_steps)
            f.write(f"# Engines: {', '.join(sorted(engines_used))}\n")

            f.write("\n")
            f.write("# Exit on error\n")
            f.write("set -e\n\n")

            # Track the last production step's key so later production steps
            # can chain CpHMD restarts from ../step{prev}/step{prev}.cprestrt.
            prev_prod_step_key: Optional[str] = None

            # Write each step. Each step cds into its own subdirectory
            # (step1/, step2/, ...) and references the shared prmtop + initial
            # coords via '../'. Output artifacts land in that subdir.
            for i, (sim_config, mdin_file) in enumerate(zip(workflow_steps, mdin_files), 1):
                engine = sim_config.engine
                hw_config = sim_config.hardware_config or {}
                step_key = f"step{i}"

                f.write(f"# ========================================\n")
                f.write(f"# Step {i}: {sim_config.name}\n")
                f.write(f"# Engine: {engine}\n")
                f.write(f"# ========================================\n")

                # Paths as seen from inside step_key/. -c/-ref are either the
                # initial coords (step 1) or the previous step's restart.
                if i == 1:
                    input_coord = f"../{coord_name}"
                else:
                    input_coord = f"../step{i-1}/step{i-1}.rst7"

                topology_arg = f"../{prmtop_name}"
                mdin_arg = "simulation.mdin"
                output_file = f"{step_key}.mdout"
                restart_file = f"{step_key}.rst7"
                trajectory_file = f"{step_key}.nc"
                info_file = f"{step_key}.mdinfo"

                # Add GPU environment variable if needed for this step
                if engine == "pmemd.cuda" and hw_config.get('gpu_ids'):
                    f.write(f"export CUDA_VISIBLE_DEVICES={hw_config['gpu_ids']}\n")

                # Checkpoint detection from workflow_dir (outer cwd).
                f.write(f"\n# Check if step already completed\n")
                f.write(f"if [ -f \"{step_key}/{restart_file}\" ]; then\n")
                f.write(f"    echo \"✓ Step {i} already completed ({step_key}/{restart_file} exists), skipping...\"\n")
                f.write(f"else\n")
                f.write(f"    # Run step {i}\n")
                f.write(f"    cd {step_key}\n")

                indent = "    "
                if engine == "sander":
                    cmd_parts = [
                        "sander",
                        "-O",
                        "-i", mdin_arg,
                        "-p", topology_arg,
                        "-c", input_coord,
                        "-o", output_file,
                        "-r", restart_file,
                        "-x", trajectory_file,
                        "-ref", input_coord,
                    ]
                elif engine == "pmemd":
                    cmd_parts = [
                        "pmemd",
                        "-O",
                        "-i", mdin_arg,
                        "-p", topology_arg,
                        "-c", input_coord,
                        "-o", output_file,
                        "-r", restart_file,
                        "-x", trajectory_file,
                        "-ref", input_coord,
                        "-inf", info_file,
                    ]
                elif engine == "pmemd.MPI":
                    mpi_tasks = hw_config.get('mpi_tasks', 4)
                    cmd_parts = [
                        "mpirun", "-np", str(mpi_tasks), "pmemd.MPI",
                        "-O",
                        "-i", mdin_arg,
                        "-p", topology_arg,
                        "-c", input_coord,
                        "-o", output_file,
                        "-r", restart_file,
                        "-x", trajectory_file,
                        "-ref", input_coord,
                        "-inf", info_file,
                    ]
                elif engine == "pmemd.cuda":
                    cmd_parts = [
                        "pmemd.cuda",
                        "-O",
                        "-i", mdin_arg,
                        "-p", topology_arg,
                        "-c", input_coord,
                        "-o", output_file,
                        "-r", restart_file,
                        "-x", trajectory_file,
                        "-ref", input_coord,
                        "-inf", info_file,
                    ]
                else:
                    cmd_parts = [
                        engine,
                        "-O",
                        "-i", mdin_arg,
                        "-p", topology_arg,
                        "-c", input_coord,
                        "-o", output_file,
                        "-r", restart_file,
                        "-x", trajectory_file,
                        "-ref", input_coord,
                    ]

                # CpHMD flags for production steps
                is_prod = self._is_production_step(sim_config)
                if cpmd_active and is_prod:
                    # Swap topology if explicit-solvent CpHMD. cpmd_prmtop is
                    # staged inside this step_dir (see _prepare_workflow_files),
                    # so reference it by basename.
                    if cpmd_prmtop:
                        cmd_parts = [cpmd_prmtop if x == topology_arg else x for x in cmd_parts]

                    cpout_file = f"{step_key}.cpout"
                    cprestrt_file = f"{step_key}.cprestrt"
                    if prev_prod_step_key:
                        cpin_input = f"../{prev_prod_step_key}/{prev_prod_step_key}.cprestrt"
                    else:
                        # First production step — cpin staged in this subdir
                        cpin_input = cpin_filename

                    cmd_parts.extend([
                        "-cpin", cpin_input,
                        "-cpout", cpout_file,
                        "-cprestrt", cprestrt_file,
                    ])
                    prev_prod_step_key = step_key

                f.write(f"{indent}echo \"Starting Step {i}: {sim_config.name}...\"\n")
                f.write(f"{indent}{cmd_parts[0]} \\\n")
                for arg in cmd_parts[1:-1]:
                    f.write(f"{indent}  {arg} \\\n")
                f.write(f"{indent}  {cmd_parts[-1]} < /dev/null\n")
                f.write(f"{indent}EXIT_CODE=$?\n")
                f.write(f"{indent}cd ..\n")
                f.write(f"{indent}\n")
                f.write(f"{indent}if [ $EXIT_CODE -eq 0 ]; then\n")
                f.write(f"{indent}    echo \"✓ Step {i} completed successfully\"\n")
                f.write(f"{indent}else\n")
                f.write(f"{indent}    echo \"✗ Step {i} failed with exit code $EXIT_CODE\"\n")
                f.write(f"{indent}    exit 1\n")
                f.write(f"{indent}fi\n")
                f.write(f"fi\n\n")

            # Add extended production cycles if requested. These run inside
            # the last step's subdirectory — the .rst7/.cprestrt chain stays
            # local, the shared prmtop is referenced via '../'.
            if extended_production_cycles > 0:
                last_step = workflow_steps[-1]
                last_engine = last_step.engine
                last_hw_config = last_step.hardware_config or {}
                last_step_num = len(workflow_steps)
                last_step_key = f"step{last_step_num}"

                f.write("# ========================================\n")
                f.write("# Extended Production Cycles\n")
                f.write(f"# Running {extended_production_cycles} additional production simulations\n")
                f.write(f"# in {last_step_key}/, chained from {last_step_key}.rst7.\n")
                f.write("# ========================================\n\n")

                if last_engine == "pmemd.cuda" and last_hw_config.get('gpu_ids'):
                    f.write(f"export CUDA_VISIBLE_DEVICES={last_hw_config['gpu_ids']}\n\n")

                f.write(f"cd {last_step_key}\n")
                f.write(f"for cycle in $(seq 1 {extended_production_cycles}); do\n")
                f.write(f"    echo \"Starting Extended Production Cycle $cycle/{extended_production_cycles}...\"\n")
                f.write(f"    \n")

                # Input-coord / CPIN selection for cycle 1 vs cycle N+1.
                # All referenced files live in this (last_step_key) subdir.
                f.write(f"    # Cycle 1 chains from {last_step_key}.rst7; later cycles from prior extended restart.\n")
                f.write(f"    if [ $cycle -eq 1 ]; then\n")
                f.write(f"        INPUT_COORD=\"{last_step_key}.rst7\"\n")
                if cpmd_active:
                    # Extended cycles only run when the last protocol step is
                    # production — at that point prev_prod_step_key is this
                    # very step, so its cprestrt sits right here. Fall back to
                    # the staged cpin if somehow no prior cprestrt exists.
                    if prev_prod_step_key:
                        f.write(f"        CPIN_INPUT=\"{last_step_key}.cprestrt\"\n")
                    else:
                        f.write(f"        CPIN_INPUT=\"{cpin_filename}\"\n")
                f.write(f"    else\n")
                f.write(f"        prev_cycle=$((cycle - 1))\n")
                f.write(f"        INPUT_COORD=\"{last_step_key}_extended_${{prev_cycle}}.rst7\"\n")
                if cpmd_active:
                    f.write(f"        CPIN_INPUT=\"{last_step_key}_extended_${{prev_cycle}}.cprestrt\"\n")
                f.write(f"    fi\n")
                f.write(f"    \n")

                ext_prefix = f"{last_step_key}_extended_${{cycle}}"
                # cpmd_prmtop is staged per production step (basename-accessible).
                # Otherwise the shared prmtop is one level up.
                ext_prmtop = cpmd_prmtop if (cpmd_active and cpmd_prmtop) else f"../{prmtop_name}"

                if last_engine == "pmemd.MPI":
                    mpi_tasks = last_hw_config.get('mpi_tasks', 4)
                    engine_prefix = f"mpirun -np {mpi_tasks} pmemd.MPI"
                else:
                    engine_prefix = last_engine

                ext_args = [
                    "-O",
                    "-i simulation.mdin",
                    f"-p {ext_prmtop}",
                    f"-c \"$INPUT_COORD\"",
                    f"-o {ext_prefix}.mdout",
                    f"-r {ext_prefix}.rst7",
                    f"-x {ext_prefix}.nc",
                    f"-ref \"$INPUT_COORD\"",
                ]
                if last_engine.startswith("pmemd"):
                    ext_args.append(f"-inf {ext_prefix}.mdinfo")

                if cpmd_active:
                    ext_args.extend([
                        f"-cpin \"$CPIN_INPUT\"",
                        f"-cpout {ext_prefix}.cpout",
                        f"-cprestrt {ext_prefix}.cprestrt",
                    ])

                f.write(f"    {engine_prefix} \\\n")
                for arg in ext_args[:-1]:
                    f.write(f"      {arg} \\\n")
                f.write(f"      {ext_args[-1]} < /dev/null\n")

                f.write(f"    EXIT_CODE=$?\n")
                f.write(f"    \n")
                f.write(f"    if [ $EXIT_CODE -eq 0 ]; then\n")
                f.write(f"        echo \"✓ Extended Production Cycle $cycle completed successfully\"\n")
                f.write(f"    else\n")
                f.write(f"        echo \"✗ Extended Production Cycle $cycle failed with exit code $EXIT_CODE\"\n")
                f.write(f"        exit 1\n")
                f.write(f"    fi\n")
                f.write(f"done\n")
                f.write(f"cd ..\n\n")

            f.write("# ========================================\n")
            f.write("# Protocol Complete\n")
            f.write("# ========================================\n")
            if extended_production_cycles > 0:
                f.write(f"echo \"All protocol steps and {extended_production_cycles} extended production cycles completed successfully!\"\n")
            else:
                f.write("echo \"All steps completed successfully!\"\n")

        # Make script executable
        import os
        os.chmod(script_file, 0o755)

    def _create_run_script(self, sim_dir: Path, cmd: List[str], sim_config: SimulationConfig,
                          hardware_config: Dict[str, Any]):
        """Create a bash script documenting the exact command used."""
        script_file = sim_dir / "run_simulation.sh"

        with open(script_file, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(f"# AMBER Simulation Run Script\n")
            f.write(f"# Generated by ProPrep Molecular Dynamics Manager\n")
            f.write(f"# Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"#\n")
            f.write(f"# Simulation: {sim_config.name}\n")
            f.write(f"# Template: {sim_config.template_id}\n")
            f.write(f"# Engine: {sim_config.engine}\n")

            if hardware_config:
                if hardware_config.get('mpi_tasks'):
                    f.write(f"# MPI Tasks: {hardware_config['mpi_tasks']}\n")
                if hardware_config.get('gpu_ids'):
                    f.write(f"# GPU IDs: {hardware_config['gpu_ids']}\n")

            f.write("\n")

            # Add environment variables if needed
            if sim_config.engine == "pmemd.cuda" and hardware_config.get('gpu_ids'):
                f.write(f"export CUDA_VISIBLE_DEVICES={hardware_config['gpu_ids']}\n\n")

            # Write the command with line continuations for readability
            f.write("# Execute AMBER simulation\n")
            if len(cmd) > 1:
                # Multi-line format for better readability
                f.write(f"{cmd[0]} \\\n")
                for arg in cmd[1:-1]:
                    f.write(f"  {arg} \\\n")
                f.write(f"  {cmd[-1]} < /dev/null\n")
            else:
                f.write(f"{' '.join(cmd)} < /dev/null\n")

            f.write("\n# End of run script\n")

        # Make script executable
        import os
        os.chmod(script_file, 0o755)
    
    def _create_basic_mdin_file(self, mdin_file: Path):
        """Create a basic minimization mdin file."""
        content = """Basic minimization and short MD
 &cntrl
  imin=1, maxcyc=500, ncyc=250,
  ntb=1, ntp=0, ntf=2, ntc=2,
  cut=10.0, ntpr=50, ntwx=50,
 /
"""
        with open(mdin_file, 'w') as f:
            f.write(content)

    def _create_simulation_log(self, sim_dir: Path, sim_config: SimulationConfig, 
                             hardware_config: Dict[str, Any], result: subprocess.CompletedProcess):
        """Create detailed simulation execution log."""
        log_file = sim_dir / "execution.log"
        
        with open(log_file, 'w') as f:
            f.write(f"AMBER Simulation Execution Log\n")
            f.write(f"=" * 50 + "\n\n")
            f.write(f"Simulation: {sim_config.name}\n")
            f.write(f"Template ID: {sim_config.template_id}\n")
            f.write(f"Engine: {sim_config.engine}\n")
            f.write(f"Return Code: {result.returncode}\n")
            f.write(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            if hardware_config.get('mpi_tasks'):
                f.write(f"MPI Tasks: {hardware_config['mpi_tasks']}\n")
            if hardware_config.get('gpu_ids'):
                f.write(f"GPU IDs: {hardware_config['gpu_ids']}\n")
                
            f.write(f"\nStandard Output:\n")
            f.write(f"-" * 30 + "\n")
            f.write(result.stdout or "No output captured\n")
            
            f.write(f"\nStandard Error:\n")
            f.write(f"-" * 30 + "\n") 
            f.write(result.stderr or "No errors captured\n")

    def _show_batch_results(self, run_dir: Path):
        """Show summary of batch execution results."""
        self.console.print(f"\n[bold cyan]===== Batch Results Summary =====[/bold cyan]")
        self.console.print(f"Results directory: [blue]{run_dir}[/blue]\n")
        
        # List simulation directories
        sim_dirs = [d for d in run_dir.iterdir() if d.is_dir() and not d.name == "shared"]
        sim_dirs.sort()
        
        from rich.table import Table
        table = Table(title="Simulation Results", show_header=True)
        table.add_column("#", width=3)
        table.add_column("Simulation", style="bright_blue")
        table.add_column("Status", width=10)
        table.add_column("Files", style="grey50")
        
        for i, sim_dir in enumerate(sim_dirs, 1):
            sim_name = sim_dir.name.replace(f"{i:02d}_", "").replace("_", " ")
            
            # Check for output files
            mdout_files = list(sim_dir.glob("*.mdout"))
            nc_files = list(sim_dir.glob("*.nc"))
            rst_files = list(sim_dir.glob("*.rst7"))
            
            if mdout_files:
                # Check if simulation completed
                with open(mdout_files[0], 'r') as f:
                    content = f.read()
                if "FINAL RESULTS" in content or "Total wall time:" in content:
                    status = "[green]✓ Complete[/green]"
                else:
                    status = "[yellow]⚠ Partial[/yellow]"
            else:
                status = "[red]✗ Failed[/red]"
                
            file_count = len(mdout_files) + len(nc_files) + len(rst_files)
            files_info = f"{file_count} files"
            
            table.add_row(str(i), sim_name, status, files_info)
            
        self.console.print(table)

    def _execute_workflow(self):
        """Execute workflow."""
        # Check for available workflows
        workflow_files = self._find_workspace_files("*_workflow.json")
        
        if not workflow_files:
            self.console.print("[yellow]No protocols available. Create protocols first.[/yellow]")
            return
            
        # Check prerequisites (same as single simulation)
        prmtop_files = self._find_workspace_files("*.prmtop")
        coord_files = (self._find_workspace_files("*.rst7") + 
                      self._find_workspace_files("*.inpcrd"))
                      
        if not prmtop_files or not coord_files:
            self.console.print("[red]Missing topology or coordinate files. Run TLEaP first.[/red]")
            return
            
        # Select workflow
        if len(workflow_files) == 1:
            selected_workflow = workflow_files[0]
        else:
            self.console.print("\n[bold]Available Protocols:[/bold]")
            for i, wf in enumerate(workflow_files, 1):
                name = wf.stem.replace('_workflow', '')
                self.console.print(f"  {i}. {name}")

            choice_str = prompt_with_context(
                self.processor,
                f"Select protocol (1-{len(workflow_files)})",
                default="1",
                module="MD Manager - Execution",
                description="Select protocol to execute"
            )
            choice = int(choice_str)
            
            if 1 <= choice <= len(workflow_files):
                selected_workflow = workflow_files[choice-1]
            else:
                self.console.print("[red]Invalid selection[/red]")
                return
                
        self.console.print(f"[green]Selected protocol: {selected_workflow.stem}[/green]")
        
        # Load and execute workflow
        try:
            with open(selected_workflow, 'r') as f:
                workflow = json.load(f)
        except Exception as e:
            self.console.print(f"[red]Error loading protocol: {e}[/red]")
            return
            
        # Parse workflow structure
        workflow_name = workflow.get('name', selected_workflow.stem)
        steps = workflow.get('steps', [])
        
        if not steps:
            self.console.print("[red]Protocol contains no steps[/red]")
            return
            
        self.console.print(f"\n[bold]Executing protocol: {workflow_name}[/bold]")
        self.console.print(f"Steps to execute: {len(steps)}")
        
        # Prepare workflow execution directory
        workflow_dir = Path.cwd() / "simulations" / f"workflow_{workflow_name.lower().replace(' ', '_')}"
        workflow_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy prerequisite files to workflow directory
        for prmtop in prmtop_files:
            shutil.copy2(prmtop, workflow_dir)
        for coord in coord_files:
            shutil.copy2(coord, workflow_dir)
            
        # Execute workflow steps sequentially
        results = {}
        failed_step = None
        
        for i, step in enumerate(steps, 1):
            step_name = step.get('name', f'Step_{i}')
            step_type = step.get('type', 'simulation')
            
            self.console.print(f"\n[cyan]--- Executing Step {i}/{len(steps)}: {step_name} ---[/cyan]")
            
            if step_type == 'simulation':
                # Execute simulation step
                success = self._execute_workflow_simulation_step(step, workflow_dir, results)
                if not success:
                    failed_step = step_name
                    break
                    
            elif step_type == 'analysis':
                # Execute analysis step
                success = self._execute_workflow_analysis_step(step, workflow_dir, results)
                if not success:
                    failed_step = step_name
                    break
                    
            elif step_type == 'preparation':
                # Execute preparation step
                success = self._execute_workflow_preparation_step(step, workflow_dir, results)
                if not success:
                    failed_step = step_name
                    break
                    
            else:
                self.console.print(f"[yellow]Unknown step type: {step_type}. Skipping...[/yellow]")
                
            results[step_name] = {'status': 'completed', 'step_type': step_type}
            
        # Report workflow completion
        if failed_step:
            self.console.print(f"\n[red]Protocol failed at step: {failed_step}[/red]")
            self.console.print(f"[yellow]Results available in: {workflow_dir}[/yellow]")
        else:
            self.console.print(f"\n[green]✓ Protocol completed successfully![/green]")
            self.console.print(f"[green]All results saved to: {workflow_dir}[/green]")
            
        # Save workflow execution log
        log_file = workflow_dir / "workflow_execution.log"
        with open(log_file, 'w') as f:
            json.dump({
                'workflow_name': workflow_name,
                'execution_time': datetime.now().isoformat(),
                'status': 'failed' if failed_step else 'completed',
                'failed_step': failed_step,
                'results': results
            }, f, indent=2)

    def _execute_workflow_simulation_step(self, step, workflow_dir, results):
        """Execute a simulation step in workflow."""
        step_name = step.get('name', 'simulation')
        mdin_template = step.get('mdin_template')
        parameters = step.get('parameters', {})
        dependencies = step.get('depends_on', [])
        
        # Check dependencies
        for dep in dependencies:
            if dep not in results or results[dep]['status'] != 'completed':
                self.console.print(f"[red]Dependency not met: {dep}[/red]")
                return False
                
        # Create step directory
        step_dir = workflow_dir / step_name.lower().replace(' ', '_')
        step_dir.mkdir(exist_ok=True)
        
        # Generate mdin file from template and parameters
        if mdin_template:
            mdin_content = mdin_template
            for key, value in parameters.items():
                mdin_content = mdin_content.replace(f'{{{key}}}', str(value))
                
            mdin_file = step_dir / f"{step_name.lower().replace(' ', '_')}.mdin"
            with open(mdin_file, 'w') as f:
                f.write(mdin_content)
        else:
            self.console.print("[yellow]No mdin template provided, using default minimal MD[/yellow]")
            mdin_file = step_dir / "default.mdin"
            with open(mdin_file, 'w') as f:
                f.write("Minimal MD\n &cntrl\n  imin=0, nstlim=1000, dt=0.002\n /\n")
                
        # Copy input files to step directory
        for prmtop in workflow_dir.glob("*.prmtop"):
            shutil.copy2(prmtop, step_dir)
        for coord in list(workflow_dir.glob("*.rst7")) + list(workflow_dir.glob("*.inpcrd")):
            shutil.copy2(coord, step_dir)
            
        # Execute simulation using AmberController
        engine = step.get('engine', 'sander')
        original_cwd = os.getcwd()
        
        try:
            os.chdir(step_dir)
            
            # Use AmberController for execution
            success = self.amber_controller.run_simulation(
                str(mdin_file.name),
                list(step_dir.glob("*.prmtop"))[0].name,
                list(step_dir.glob("*.rst7") or step_dir.glob("*.inpcrd"))[0].name,
                f"{step_name.lower().replace(' ', '_')}.out",
                f"{step_name.lower().replace(' ', '_')}.rst",
                engine=engine
            )
            
            if success:
                self.console.print(f"[green]✓ Step '{step_name}' completed successfully[/green]")
                return True
            else:
                self.console.print(f"[red]✗ Step '{step_name}' failed[/red]")
                return False
                
        finally:
            os.chdir(original_cwd)
            
    def _execute_workflow_analysis_step(self, step, workflow_dir, results):
        """Execute an analysis step in workflow."""
        step_name = step.get('name', 'analysis')
        analysis_type = step.get('analysis_type', 'energy')
        dependencies = step.get('depends_on', [])
        
        # Check dependencies
        for dep in dependencies:
            if dep not in results or results[dep]['status'] != 'completed':
                self.console.print(f"[red]Dependency not met: {dep}[/red]")
                return False
                
        self.console.print(f"[cyan]Running {analysis_type} analysis...[/cyan]")
        
        # Create analysis directory
        analysis_dir = workflow_dir / f"analysis_{step_name.lower().replace(' ', '_')}"
        analysis_dir.mkdir(exist_ok=True)
        
        # Perform basic analysis based on type
        if analysis_type == 'energy':
            # Extract energy information from output files
            output_files = list(workflow_dir.rglob("*.out"))
            if output_files:
                analysis_file = analysis_dir / "energy_analysis.txt"
                with open(analysis_file, 'w') as f:
                    f.write(f"Energy Analysis for {step_name}\n")
                    f.write("=" * 40 + "\n")
                    for out_file in output_files:
                        f.write(f"\nFile: {out_file.name}\n")
                        f.write("Energy data extraction would be implemented here\n")
                        
        elif analysis_type == 'trajectory':
            # Trajectory analysis placeholder
            analysis_file = analysis_dir / "trajectory_analysis.txt"
            with open(analysis_file, 'w') as f:
                f.write(f"Trajectory Analysis for {step_name}\n")
                f.write("=" * 40 + "\n")
                f.write("Trajectory analysis would be implemented here\n")
                
        else:
            self.console.print(f"[yellow]Unknown analysis type: {analysis_type}[/yellow]")
            return False
            
        self.console.print(f"[green]✓ Analysis '{step_name}' completed[/green]")
        return True
        
    def _execute_workflow_preparation_step(self, step, workflow_dir, results):
        """Execute a preparation step in workflow."""
        step_name = step.get('name', 'preparation')
        prep_type = step.get('prep_type', 'copy_files')
        parameters = step.get('parameters', {})
        dependencies = step.get('depends_on', [])
        
        # Check dependencies
        for dep in dependencies:
            if dep not in results or results[dep]['status'] != 'completed':
                self.console.print(f"[red]Dependency not met: {dep}[/red]")
                return False
                
        self.console.print(f"[cyan]Running {prep_type} preparation...[/cyan]")
        
        if prep_type == 'copy_files':
            # Copy specified files
            source_files = parameters.get('source_files', [])
            for source_file in source_files:
                source_path = Path(source_file)
                if source_path.exists():
                    shutil.copy2(source_path, workflow_dir)
                    self.console.print(f"[green]Copied: {source_file}[/green]")
                else:
                    self.console.print(f"[yellow]File not found: {source_file}[/yellow]")
                    
        elif prep_type == 'modify_coordinates':
            # Coordinate modification placeholder
            self.console.print("[yellow]Coordinate modification would be implemented here[/yellow]")
            
        else:
            self.console.print(f"[yellow]Unknown preparation type: {prep_type}[/yellow]")
            return False
            
        self.console.print(f"[green]✓ Preparation '{step_name}' completed[/green]")
        return True

    def _monitor_simulation(self):
        """Monitor running simulation using integrated AMBERMonitor."""
        while True:
            self.console.print(f"\n[bold cyan]===== AMBER Simulation Monitoring =====[/bold cyan]")
            
            # First check for actively running simulations
            active_sims = []
            if hasattr(self, 'running_simulations'):
                for name, sim_info in self.running_simulations.items():
                    pid = sim_info['pid']
                    # Check if process is still running using PID
                    if self._is_process_running(pid):  # Still running
                        active_sims.append((name, sim_info, 'running'))
                    else:  # Process finished
                        active_sims.append((name, sim_info, 'finished'))
            
            # Also look for recent simulation directories (fallback for older simulations)
            recent_sims = []
            simulations_dir = Path.cwd() / "simulations"
            if simulations_dir.exists():
                for sim_dir in simulations_dir.iterdir():
                    if sim_dir.is_dir() and not sim_dir.name.startswith('.'):
                        # Look for .mdout files
                        mdout_files = list(sim_dir.rglob("*.mdout"))
                        if mdout_files:
                            # Get most recent mdout file
                            mdout_file = max(mdout_files, key=os.path.getmtime)
                            # Check if this simulation is not already in active_sims
                            # Compare by both directory path and mdout file path to avoid duplicates
                            already_tracked = False
                            for _, sim_info, _ in active_sims:
                                if (sim_info['sim_dir'] == sim_dir or 
                                    str(sim_info['mdout_file']) == str(mdout_file)):
                                    already_tracked = True
                                    break
                            
                            if not already_tracked:
                                recent_sims.append((sim_dir, mdout_file))
            
            # Combine and display options
            all_options = []
            
            # Add active simulations first
            if active_sims:
                self.console.print("[bold]Currently tracked simulations:[/bold]")
                for name, sim_info, status in active_sims:
                    status_color = "green" if status == "running" else "blue"
                    runtime = datetime.now() - sim_info['started_at']
                    runtime_str = str(runtime).split('.')[0]  # Remove microseconds
                    batch_dir = sim_info['sim_dir'].parent.name  # Get batch directory name
                    self.console.print(f"  • {name}")
                    self.console.print(f"    ([{status_color}]{status}[/{status_color}], runtime: {runtime_str}, batch: {batch_dir})")
                    all_options.append((f"{name} ({batch_dir})", sim_info['mdout_file'], sim_info))
            
            # Add recent simulations
            if recent_sims:
                if active_sims:
                    self.console.print("\n[grey50]Other recent simulations:[/grey50]")
                for sim_dir, mdout_file in recent_sims:
                    mod_time = datetime.fromtimestamp(os.path.getmtime(mdout_file))
                    self.console.print(f"  • {sim_dir.name} (modified: {mod_time.strftime('%H:%M:%S')})")
                    all_options.append((sim_dir.name, mdout_file, None))
            
            if not all_options:
                self.console.print("[yellow]No simulations found for monitoring[/yellow]")
                return
                
            # Select simulation to monitor
            if len(all_options) == 1:
                sim_name, mdout_file, sim_info = all_options[0]
                self.console.print(f"[grey50]Monitoring: {sim_name}[/grey50]")
            else:
                self.console.print(f"\n[bold]Select simulation to monitor:[/bold]")
                for i, (sim_name, mdout_file, sim_info) in enumerate(all_options, 1):
                    self.console.print(f"  {i}. {sim_name}")
                self.console.print(f"  b. ← Back to Molecular Dynamics Manager")
                
                try:
                    valid_choices = [str(i) for i in range(1, len(all_options) + 1)] + ["b"]

                    # Build options map
                    options_map = {}
                    for i, (sim_name, _, _) in enumerate(all_options, 1):
                        options_map[str(i)] = sim_name
                    options_map["b"] = "← Back to MD Manager"

                    choice = prompt_with_context(
                        self.processor,
                        "Select simulation",
                        choices=valid_choices,
                        default="1",
                        module="MD Manager - Monitoring",
                        description="Select simulation to monitor",
                        options_map=options_map
                    )
                    
                    if choice == "b":
                        return  # Back to MD Manager
                    elif choice.isdigit():
                        choice_idx = int(choice) - 1
                        if 0 <= choice_idx < len(all_options):
                            sim_name, mdout_file, sim_info = all_options[choice_idx]
                        else:
                            self.console.print("[red]Invalid selection[/red]")
                            return
                    else:
                        self.console.print("[red]Invalid selection[/red]")
                        return
                except:
                    return
                    
            # Create monitor for the selected simulation
            try:
                monitor = AMBERMonitor(str(mdout_file))
                monitor.parse_amber_output()
                
                if not any(monitor.data.values()):
                    self.console.print("[yellow]No monitoring data found in output file[/yellow]")
                    return
                    
                self.console.print(f"\n[bold]Monitoring: {sim_dir.name}[/bold]")
                self.console.print(f"Output file: [blue]{mdout_file.name}[/blue]")
                
                # Show monitoring interface
                self._show_monitoring_interface(monitor, str(mdout_file))
                # After monitoring interface exits, return to simulation list
                break
                
            except Exception as e:
                self.console.print(f"[red]Error creating monitor: {e}[/red]")

    def _show_monitoring_interface(self, monitor: AMBERMonitor, output_file: str):
        """Show interactive monitoring interface with ASCII plots."""
        while True:
            self.console.print(f"\n[bold cyan]🔍 AMBER Simulation Monitor[/bold cyan]")
            self.console.print("📊 [bold]Data & Status:[/bold]")
            self.console.print(" 1. Show current status")
            self.console.print(" 2. Display energy summary")
            self.console.print(" 3. Display temperature summary")
            self.console.print(" 4. Display pressure summary")
            self.console.print(" 5. Show timing information")
            self.console.print("")
            self.console.print("📈 [bold]ASCII Plots:[/bold]")
            self.console.print(" 6. Plot energy (total)")
            self.console.print(" 7. Plot temperature")
            self.console.print(" 8. Plot pressure")
            self.console.print(" 9. Plot kinetic energy")
            self.console.print("10. Plot potential energy", highlight=False)
            self.console.print("11. Plot RMS gradient", highlight=False)
            self.console.print("")
            self.console.print("🔄 [bold]Control:[/bold]")
            self.console.print(" r. Refresh data")
            self.console.print(" b. ← Back")

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1","2","3","4","5","6","7","8","9","10","11","r","b"],
                default="1",
                module="MD Manager - Monitoring",
                description="Simulation monitoring options",
                options_map={
                    "1": "Show current status",
                    "2": "Display energy summary",
                    "3": "Display temperature summary",
                    "4": "Display pressure summary",
                    "5": "Show timing information",
                    "6": "Plot energy (total)",
                    "7": "Plot temperature",
                    "8": "Plot pressure",
                    "9": "Plot kinetic energy",
                    "10": "Plot potential energy",
                    "11": "Plot RMS gradient",
                    "r": "Refresh data",
                    "b": "← Back"
                }
            )
            
            # Always refresh data before showing anything
            monitor.parse_amber_output()
            
            if choice == "1":
                self._show_simulation_status(monitor)
            elif choice == "2":
                self._show_energy_summary(monitor)
            elif choice == "3":
                self._show_temperature_summary(monitor)
            elif choice == "4":
                self._show_pressure_summary(monitor)
            elif choice == "5":
                self._show_timing_information(monitor)
            elif choice == "6":
                self._show_ascii_plot(monitor, 'total_energy', 'Total Energy', '(kcal/mol)')
            elif choice == "7":
                self._show_ascii_plot(monitor, 'temperature', 'Temperature', '(K)')
            elif choice == "8":
                self._show_ascii_plot(monitor, 'pressure', 'Pressure', '(bar)')
            elif choice == "9":
                self._show_ascii_plot(monitor, 'kinetic_energy', 'Kinetic Energy', '(kcal/mol)')
            elif choice == "10":
                self._show_ascii_plot(monitor, 'potential_energy', 'Potential Energy', '(kcal/mol)')
            elif choice == "11":
                self._show_ascii_plot(monitor, 'rms_gradient', 'RMS Gradient', '(kcal/mol/Å)')
            elif choice == "r":
                self.console.print("[green]🔄 Data refreshed[/green]")
            elif choice == "b":
                break

    def _show_ascii_plot(self, monitor: AMBERMonitor, data_key: str, title: str, units: str):
        """Display ASCII plot using the monitor's built-in plotting."""
        self.console.print(f"\n[bold green]📈 {title} {units}[/bold green]")
        
        # Get the ASCII plot from the monitor
        plot_output = monitor.ascii_plot(data_key, title, units, height=15, width=60)
        
        # Display the plot (convert to Rich console output)
        self.console.print(f"[bright_blue]{plot_output}[/bright_blue]")
        
        # Add some additional info
        if data_key in monitor.data and monitor.data[data_key]:
            values = list(monitor.data[data_key])
            if values:
                latest = values[-1]
                avg = sum(values) / len(values)
                self.console.print(f"\n[grey50]Latest: {latest:.2f} {units} | Average: {avg:.2f} {units} | Points: {len(values)}[/grey50]")

    def _show_simulation_status(self, monitor: AMBERMonitor):
        """Show current simulation status."""
        self.console.print(f"\n[bold blue]Current Simulation Status[/bold blue]")
        
        if not monitor.data.get('step'):
            self.console.print("[bright_yellow]No step data available[/bright_yellow]")
            return
            
        # Get latest values
        latest_step = monitor.data['step'][-1] if monitor.data['step'] else 0
        latest_time = monitor.data['time'][-1] if monitor.data['time'] else 0
        latest_temp = monitor.data['temperature'][-1] if monitor.data['temperature'] else None
        latest_energy = monitor.data['total_energy'][-1] if monitor.data['total_energy'] else None
        
        from rich.table import Table
        status_table = Table(title="Current Status")
        status_table.add_column("Property", style="bright_blue")
        status_table.add_column("Value", style="green")
        
        status_table.add_row("Step", f"{latest_step:,}")
        status_table.add_row("Time (ps)", f"{latest_time:.3f}")
        if latest_temp is not None:
            status_table.add_row("Temperature (K)", f"{latest_temp:.2f}")
        if latest_energy is not None:
            status_table.add_row("Total Energy (kcal/mol)", f"{latest_energy:.2f}")
            
        # Calculate progress if we can estimate
        if latest_time > 0:
            status_table.add_row("Progress", f"Step {latest_step:,} at {latest_time:.1f} ps")
            
        self.console.print(status_table)

    def _show_energy_summary(self, monitor: AMBERMonitor):
        """Show energy analysis summary."""
        self.console.print(f"\n[bold blue]Energy Summary[/bold blue]")
        
        if not monitor.data.get('total_energy'):
            self.console.print("[yellow]No energy data available[/yellow]")
            return
            
        # Calculate statistics
        energies = monitor.data['total_energy']
        if len(energies) > 1:
            avg_energy = sum(energies) / len(energies)
            min_energy = min(energies)
            max_energy = max(energies)
            
            from rich.table import Table
            energy_table = Table(title="Energy Statistics")
            energy_table.add_column("Metric", style="bright_blue")
            energy_table.add_column("Value (kcal/mol)", style="green")
            
            energy_table.add_row("Average", f"{avg_energy:.2f}")
            energy_table.add_row("Minimum", f"{min_energy:.2f}")
            energy_table.add_row("Maximum", f"{max_energy:.2f}")
            energy_table.add_row("Range", f"{max_energy - min_energy:.2f}")
            energy_table.add_row("Current", f"{energies[-1]:.2f}")
            
            self.console.print(energy_table)
            
            # Show trend
            if len(energies) > 10:
                recent_avg = sum(energies[-10:]) / 10
                trend = "stable"
                if recent_avg > avg_energy + 50:
                    trend = "increasing"
                elif recent_avg < avg_energy - 50:
                    trend = "decreasing"
                self.console.print(f"Recent trend: [cyan]{trend}[/cyan]")

    def _show_temperature_summary(self, monitor: AMBERMonitor):
        """Show temperature analysis summary."""
        self.console.print(f"\n[bold blue]Temperature Summary[/bold blue]")
        
        if not monitor.data.get('temperature'):
            self.console.print("[yellow]No temperature data available[/yellow]")
            return
            
        temps = monitor.data['temperature']
        if len(temps) > 1:
            avg_temp = sum(temps) / len(temps)
            min_temp = min(temps)
            max_temp = max(temps)
            
            from rich.table import Table
            temp_table = Table(title="Temperature Statistics")
            temp_table.add_column("Metric", style="bright_blue")
            temp_table.add_column("Value (K)", style="green")
            
            temp_table.add_row("Average", f"{avg_temp:.2f}")
            temp_table.add_row("Minimum", f"{min_temp:.2f}")
            temp_table.add_row("Maximum", f"{max_temp:.2f}")
            temp_table.add_row("Range", f"{max_temp - min_temp:.2f}")
            temp_table.add_row("Current", f"{temps[-1]:.2f}")
            
            self.console.print(temp_table)

    def _show_pressure_summary(self, monitor: AMBERMonitor):
        """Show pressure analysis summary."""
        self.console.print(f"\n[bold blue]Pressure Summary[/bold blue]")
        
        if not monitor.data.get('pressure'):
            self.console.print("[yellow]No pressure data available[/yellow]")
            return
            
        pressures = monitor.data['pressure']
        if len(pressures) > 1:
            avg_pressure = sum(pressures) / len(pressures)
            
            from rich.table import Table
            pressure_table = Table(title="Pressure Statistics")
            pressure_table.add_column("Metric", style="bright_blue")
            pressure_table.add_column("Value (bar)", style="green")
            
            pressure_table.add_row("Average", f"{avg_pressure:.2f}")
            pressure_table.add_row("Current", f"{pressures[-1]:.2f}")
            
            self.console.print(pressure_table)

    def _show_timing_information(self, monitor: AMBERMonitor):
        """Show simulation timing information."""
        self.console.print(f"\n[bold blue]Timing Information[/bold blue]")
        
        if monitor.timing_data:
            from rich.table import Table
            timing_table = Table(title="Performance Metrics")
            timing_table.add_column("Metric", style="bright_blue")
            timing_table.add_column("Value", style="green")
            
            for key, value in monitor.timing_data.items():
                timing_table.add_row(key, str(value))
                
            self.console.print(timing_table)
        else:
            self.console.print("[yellow]No timing information available[/yellow]")

    def _analyze_simulations(self):
        """Analyze completed simulations using integrated analysis tools."""
        self.console.print(f"\n[bold cyan]===== Simulation Analysis =====[/bold cyan]")

        # Offer a shortcut for trajectory analysis using topology + trajectory
        # already loaded into the workspace (e.g. via the Structure Loader).
        if self.workspace:
            ws_prmtop = self.workspace.get("parm7_file")
            ws_trajs = self.workspace.get("trajectory_files") or []
            ws_trajs = [p for p in ws_trajs if os.path.exists(p)]
            if ws_prmtop and os.path.exists(ws_prmtop) and ws_trajs:
                self.console.print(
                    f"\n[bold]Workspace-loaded files available:[/bold]"
                )
                self.console.print(
                    f"  Topology:   [cyan]{os.path.basename(ws_prmtop)}[/cyan]"
                )
                self.console.print(
                    f"  Trajectory: [cyan]{len(ws_trajs)} segment(s)[/cyan]"
                )
                if confirm_with_context(
                    self.processor,
                    "Run trajectory analysis on workspace-loaded files?",
                    default=True,
                    module="MD Manager - Analysis",
                    description="Use workspace-loaded topology + trajectory",
                ):
                    self._run_trajectory_analysis(
                        [Path(p) for p in ws_trajs], Path(ws_prmtop)
                    )
                    return True

        # Look for simulation output directories
        simulations_path = Path.cwd() / "simulations"

        # Find completed simulations recursively: any directory that directly
        # contains a .mdout file. Batch runs nest output as
        # simulations/batch_*/NN_<step>/<name>.mdout, so a top-level glob misses
        # them — collect the leaf step directories instead (downstream code globs
        # each leaf non-recursively for its mdout/nc/prmtop files).
        sim_dirs = []
        if simulations_path.exists():
            seen = set()
            for mdout in sorted(simulations_path.rglob("*.mdout")):
                leaf_dir = mdout.parent
                if leaf_dir in seen:
                    continue
                # Skip hidden directories anywhere along the path
                rel_parts = leaf_dir.relative_to(simulations_path).parts
                if any(part.startswith('.') for part in rel_parts):
                    continue
                seen.add(leaf_dir)
                sim_dirs.append(leaf_dir)
        else:
            self.console.print("[yellow]No simulations directory found[/yellow]")

        if not sim_dirs:
            self.console.print("[yellow]No completed simulations found in simulations/ directory[/yellow]")
            self.console.print("[grey50]Browse manually to analyze output files from other locations[/grey50]")
            sim_dirs = []  # Empty list to trigger manual browse option
            
        # Sort by modification time (most recent first)
        if sim_dirs:
            sim_dirs.sort(key=lambda d: max(f.stat().st_mtime for f in d.glob("*.mdout")), reverse=True)

        selected_item = None

        if len(sim_dirs) == 1:
            # Single simulation found - offer to use it or browse manually
            self.console.print(f"\nFound simulation: [cyan]{sim_dirs[0].name}[/cyan]")

            choice_str = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2"],
                default="1",
                module="MD Manager - Analysis",
                description="Analyze found simulation or browse manually",
                options_map={
                    "1": "Analyze this simulation",
                    "2": "Browse for output files manually"
                }
            )

            if choice_str == "1":
                selected_item = ("dir", sim_dirs[0])
            else:
                selected_item = ("browse", None)

        elif len(sim_dirs) > 1:
            # Multiple simulations found
            self.console.print("\nAvailable simulations for analysis:")
            for i, sim_dir in enumerate(sim_dirs[:10], 1):
                # Get simulation info
                mdout_files = list(sim_dir.glob("*.mdout"))
                latest_mdout = max(mdout_files, key=os.path.getmtime) if mdout_files else None
                if latest_mdout:
                    mod_time = datetime.fromtimestamp(os.path.getmtime(latest_mdout))
                    self.console.print(f"  {i}. {sim_dir.name} (completed: {mod_time.strftime('%Y-%m-%d %H:%M')})")
                else:
                    self.console.print(f"  {i}. {sim_dir.name}")

            # Build options map
            options_map = {}
            for i, sim_dir in enumerate(sim_dirs[:10], 1):
                options_map[str(i)] = sim_dir.name
            options_map["browse"] = "Browse for output files manually"

            try:
                choice_str = prompt_with_context(
                    self.processor,
                    f"Select simulation (1-{min(10, len(sim_dirs))}) or 'browse'",
                    choices=[str(i) for i in range(1, min(10, len(sim_dirs))+1)] + ["browse"],
                    default="1",
                    module="MD Manager - Analysis",
                    description="Select simulation to analyze",
                    options_map=options_map
                )

                if choice_str == "browse":
                    selected_item = ("browse", None)
                else:
                    choice = int(choice_str)
                    if 1 <= choice <= len(sim_dirs):
                        selected_item = ("dir", sim_dirs[choice-1])
                    else:
                        self.console.print("[red]Invalid selection[/red]")
                        return
            except ValueError:
                return
        else:
            # No simulations found - go directly to browse
            selected_item = ("browse", None)

        # Handle the selected option
        if selected_item:
            item_type, item_value = selected_item

            if item_type == "browse":
                # Browse for analysis files (mdout and/or nc)
                selected_files = self._browse_for_analysis_files()
                if selected_files:
                    # Route based on what was selected
                    self._route_analysis_by_file_type(selected_files)
            elif item_type == "dir":
                # Analyze simulation directory
                self._analyze_single_simulation(item_value)

        return True  # Indicate successful execution

    def _analyze_single_simulation(self, sim_dir: Path):
        """Analyze a single completed simulation (energetics and/or trajectory)."""
        self.console.print(f"\n[bold]Analyzing Simulation: {sim_dir.name}[/bold]")

        # Find output files
        mdout_files = list(sim_dir.glob("*.mdout"))
        nc_files = list(sim_dir.glob("*.nc"))
        prmtop_files = list(sim_dir.glob("*.prmtop"))

        # Show what's available
        self.console.print("\n[bold]Available Analysis Files:[/bold]")
        if mdout_files:
            self.console.print(f"  • {len(mdout_files)} mdout file(s) (energetics)")
        if nc_files:
            total_size_mb = sum(f.stat().st_size for f in nc_files) / (1024**2)
            self.console.print(f"  • {len(nc_files)} trajectory file(s) ({total_size_mb:.1f} MB total)")
        if prmtop_files:
            self.console.print(f"  • {len(prmtop_files)} topology file(s)")

        if not mdout_files and not nc_files:
            self.console.print("[red]No analysis files found (need .mdout or .nc files)[/red]")
            return

        # Determine analysis options based on available files
        analysis_options = {}
        if mdout_files:
            analysis_options["1"] = "Energetics analysis (mdout)"
        if nc_files and prmtop_files:
            analysis_options["2"] = "Trajectory analysis (structural)"
        if mdout_files and nc_files and prmtop_files:
            analysis_options["3"] = "Combined analysis (recommended)"

        # If only one option, do it automatically
        if len(analysis_options) == 1:
            choice = list(analysis_options.keys())[0]
            self.console.print(f"\n[grey50]Running {analysis_options[choice]}...[/grey50]")
        else:
            # Ask user to choose
            self.console.print("\n[bold]Select Analysis Type:[/bold]")
            for key, desc in analysis_options.items():
                self.console.print(f"  {key}. {desc}")

            choice = prompt_with_context(
                self.processor,
                "Select analysis type",
                choices=list(analysis_options.keys()),
                default="3" if "3" in analysis_options else "1",
                module="MD Manager - Analysis",
                description="Choose analysis type",
                options_map=analysis_options
            )

        # Perform selected analysis
        try:
            if choice == "1":
                # Energetics only
                mdout_file = max(mdout_files, key=os.path.getmtime)
                monitor = AMBERMonitor.create_historical_monitor(str(mdout_file))

                if not any(monitor.data.values()):
                    self.console.print("[yellow]No analysis data found in output file[/yellow]")
                    return

                self._confirm_simulation_stage(monitor)
                self._run_quality_assessment(monitor, sim_dir.name)

            elif choice == "2":
                # Trajectory only
                # Browse for trajectory files (allows multi-file selection)
                self.console.print("\n[grey50]Select trajectory file(s) to analyze[/grey50]")
                selected_nc_files = self._browse_for_trajectory_files(start_dir=sim_dir)

                if not selected_nc_files:
                    self.console.print("[yellow]No trajectory files selected[/yellow]")
                    return

                # Get topology file
                if len(prmtop_files) == 1:
                    prmtop = prmtop_files[0]
                else:
                    self.console.print("\n[bold]Multiple topology files found:[/bold]")
                    for i, f in enumerate(prmtop_files, 1):
                        self.console.print(f"  {i}. {f.name}")

                    prmtop_choice = prompt_with_context(
                        self.processor,
                        f"Select topology file (1-{len(prmtop_files)})",
                        default="1",
                        module="MD Manager - Trajectory Analysis",
                        description="Select topology"
                    )
                    # Replay by basename so a changed topology list can't mis-pick.
                    prmtop_choice = remap_recorded_index(self.processor, prmtop_files, str(prmtop_choice))
                    prmtop = prmtop_files[int(prmtop_choice) - 1]
                    annotate_selected_path(self.processor, prmtop)

                # Run trajectory analysis
                self._analyze_trajectory(selected_nc_files, prmtop, sim_dir.name)

            elif choice == "3":
                # Combined analysis
                mdout_file = max(mdout_files, key=os.path.getmtime)

                # Run energetics analysis first
                self.console.print("\n" + "="*70)
                self.console.print("ENERGETICS ANALYSIS")
                self.console.print("="*70)

                monitor = AMBERMonitor.create_historical_monitor(str(mdout_file))

                if any(monitor.data.values()):
                    self._confirm_simulation_stage(monitor)
                    self._run_quality_assessment(monitor, sim_dir.name)
                else:
                    self.console.print("[yellow]No energetics data found in mdout file[/yellow]")

                # Then run trajectory analysis
                self.console.print("\n" + "="*70)
                self.console.print("TRAJECTORY ANALYSIS")
                self.console.print("="*70)

                # Browse for trajectory files
                selected_nc_files = self._browse_for_trajectory_files(start_dir=sim_dir)

                if selected_nc_files:
                    # Get topology
                    if len(prmtop_files) == 1:
                        prmtop = prmtop_files[0]
                    else:
                        self.console.print("\n[bold]Multiple topology files found:[/bold]")
                        for i, f in enumerate(prmtop_files, 1):
                            self.console.print(f"  {i}. {f.name}")

                        prmtop_choice = prompt_with_context(
                            self.processor,
                            f"Select topology file (1-{len(prmtop_files)})",
                            default="1",
                            module="MD Manager - Trajectory Analysis",
                            description="Select topology"
                        )
                        # Replay by basename so a changed topology list can't mis-pick.
                        prmtop_choice = remap_recorded_index(self.processor, prmtop_files, str(prmtop_choice))
                        prmtop = prmtop_files[int(prmtop_choice) - 1]
                        annotate_selected_path(self.processor, prmtop)

                    # Run trajectory analysis
                    self._analyze_trajectory(selected_nc_files, prmtop, sim_dir.name)
                else:
                    self.console.print("[yellow]No trajectory files selected - skipping trajectory analysis[/yellow]")

        except Exception as e:
            self.console.print(f"[red]Error analyzing simulation: {e}[/red]")
            import traceback
            traceback.print_exc()

    def _confirm_simulation_stage(self, monitor: AMBERMonitor):
        """Infer simulation stage and ask user to confirm or correct."""
        from rich.panel import Panel

        # Infer stage from data
        inferred_stage, confidence, evidence = monitor.infer_simulation_stage()

        # Get system info
        natom = monitor.get_natom()

        # Display inference
        self.console.print(f"\n[bold cyan]Simulation Type Detection[/bold cyan]")
        self.console.print(f"System size: {natom:,} atoms")

        # Format stage name for display
        stage_display_map = {
            'minimization': 'Minimization',
            'nvt_heating': 'Heating (NVT)',
            'npt_heating': 'Heating (NPT)',
            'nvt_equilibration': 'NVT Equilibration/Production',
            'npt_equilibration': 'NPT Equilibration/Production',
        }

        inferred_display = stage_display_map.get(inferred_stage, inferred_stage)

        # Show evidence
        self.console.print(f"\n[bold]Detected: [green]{inferred_display}[/green][/bold] [grey50](confidence: {confidence})[/grey50]")
        if evidence:
            self.console.print("[grey50]Evidence:[/grey50]")
            for ev in evidence:
                self.console.print(f"  [grey50]• {ev}[/grey50]")

        # Stage options
        stage_options = [
            ('1', 'minimization', 'Minimization'),
            ('2', 'nvt_heating', 'Heating (NVT)'),
            ('3', 'npt_heating', 'Heating (NPT)'),
            ('4', 'nvt_equilibration', 'NVT Equilibration/Production'),
            ('5', 'npt_equilibration', 'NPT Equilibration/Production')
        ]

        # Find the detected stage number for default
        default_choice = '1'
        for choice, stage_key, _ in stage_options:
            if stage_key == inferred_stage:
                default_choice = choice
                break

        # Build options map with detected stage highlighted
        options_map = {}
        for choice, stage_key, stage_label in stage_options:
            if stage_key == inferred_stage:
                options_map[choice] = f"{stage_label} (detected)"
            else:
                options_map[choice] = stage_label

        # Display options
        self.console.print("\n[bold]Simulation Types:[/bold]")
        for choice, stage_key, stage_label in stage_options:
            if stage_key == inferred_stage:
                self.console.print(f"  {choice}. {stage_label} [green](detected)[/green]")
            else:
                self.console.print(f"  {choice}. {stage_label}")

        # Ask user to confirm
        choice = prompt_with_context(
            self.processor,
            "Confirm simulation type",
            choices=[c for c, _, _ in stage_options],
            default=default_choice,
            module="MD Manager - Analysis",
            description="Select simulation type",
            options_map=options_map
        )

        # Map choice to stage
        for ch, stage_key, _ in stage_options:
            if ch == choice:
                monitor.simulation_stage = stage_key
                break

        # Show confirmation
        confirmed_display = stage_display_map.get(monitor.simulation_stage, monitor.simulation_stage)
        if monitor.simulation_stage == inferred_stage:
            self.console.print(f"[green]✓ Confirmed: {confirmed_display}[/green]")
        else:
            self.console.print(f"[yellow]Updated to: {confirmed_display}[/yellow]")

    def _run_quality_assessment(self, monitor: AMBERMonitor, sim_name: str):
        """Run comprehensive quality assessment based on simulation stage."""
        from rich.panel import Panel

        stage = monitor.simulation_stage
        natom = monitor.get_natom()

        # Format stage name for display
        stage_display_map = {
            'minimization': 'Minimization',
            'nvt_heating': 'Heating (NVT)',
            'npt_heating': 'Heating (NPT)',
            'nvt_equilibration': 'NVT Equilibration/Production',
            'npt_equilibration': 'NPT Equilibration/Production',
        }
        stage_display = stage_display_map.get(stage, stage.replace('_', ' ').title())

        self.console.print(f"\n[bold cyan]===== Simulation Analysis =====[/bold cyan]")
        self.console.print(f"Simulation: [cyan]{sim_name}[/cyan]")
        self.console.print(f"Stage: [yellow]{stage_display}[/yellow]")
        self.console.print(f"System: {natom:,} atoms\n")

        # Dispatch to stage-specific assessment
        if stage == 'minimization':
            self._assess_minimization(monitor)
        elif stage in ['nvt_heating', 'npt_heating']:
            self._assess_heating(monitor, stage)
        elif stage == 'npt_equilibration':
            self._assess_npt_equilibration(monitor)
        elif stage == 'nvt_equilibration':
            self._assess_nvt_equilibration_production(monitor, stage)
        else:
            self.console.print(f"[yellow]No specific quality assessment for stage: {stage}[/yellow]")
            self.console.print("[grey50]Showing generic analysis...[/grey50]")
            self._show_generic_analysis(monitor)

    def _assess_minimization(self, monitor: AMBERMonitor):
        """Quality assessment for minimization. Reports raw metrics without
        arbitrary convergence thresholds; the user judges adequacy."""
        from rich.table import Table
        import numpy as np

        self.console.print("[bold]Minimization Metrics[/bold]\n")

        if not monitor.data.get('step'):
            self.console.print("[red]No minimization data available[/red]")
            return

        steps = list(monitor.data['step'])
        energies = list(monitor.data['total_energy'])
        rms_values = list(monitor.data.get('rms_gradient', []))
        gmax_values = list(monitor.data.get('gmax', []))

        # Calculate metrics
        results = {}
        results['total_steps'] = steps[-1] if steps else 0
        results['initial_energy'] = energies[0] if energies else None
        results['final_energy'] = energies[-1] if energies else None
        results['energy_change'] = energies[-1] - energies[0] if len(energies) > 0 else None

        if rms_values:
            results['final_rms'] = rms_values[-1]

        if gmax_values:
            results['final_gmax'] = gmax_values[-1]

        # Check monotonic decrease
        if len(energies) > 1:
            energy_array = np.array(energies)
            results['monotonic_decrease'] = np.all(np.diff(energy_array) <= 0)
        else:
            results['monotonic_decrease'] = None

        # Display results table
        table = Table(title="Minimization Results", show_header=True, header_style="bold bright_blue")
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        table.add_row("Total Steps", f"{results['total_steps']}")
        table.add_row("Initial Energy", f"{results['initial_energy']:.2f} kcal/mol")
        table.add_row("Final Energy", f"{results['final_energy']:.2f} kcal/mol")
        table.add_row("Energy Change", f"{results['energy_change']:.2f} kcal/mol")

        if 'final_rms' in results:
            table.add_row(
                "Final RMS Gradient",
                f"{results['final_rms']:.4f}"
            )

        if 'final_gmax' in results:
            table.add_row(
                "Final GMAX",
                f"{results['final_gmax']:.2f} kcal/mol/Å"
            )

        if results['monotonic_decrease'] is not None:
            table.add_row(
                "Monotonic Decrease",
                "Yes" if results['monotonic_decrease'] else "No"
            )

        self.console.print(table)

        # ASCII Plots
        self.console.print("\n" + "="*70)
        self.console.print("MINIMIZATION - ASCII VISUALIZATION")
        self.console.print("="*70)

        # Energy plot
        if energies:
            # Check if first point is outlier
            skip_first = False
            if len(energies) > 1:
                energy_array = np.array(energies)
                median_rest = np.median(energy_array[1:])
                if abs(energy_array[0]) > 1000 * abs(median_rest):
                    skip_first = True

            if skip_first:
                note = f"\nNote: First point excluded for clarity (E₀={energies[0]:.2e} kcal/mol)"
                self.console.print(note, highlight=False)
                plot_energies = energies[1:]
                plot_steps = steps[1:]
            else:
                plot_energies = energies
                plot_steps = steps

            energy_plot = self._create_ascii_plot(
                plot_energies,
                title="Potential Energy vs Step (Convergence)",
                ylabel="Energy (kcal/mol)",
                x_values=plot_steps
            )
            self.console.print(energy_plot, highlight=False)

        # RMS gradient plot
        if rms_values:
            # Check if first RMS point is an outlier (independently from energy)
            skip_first_rms = False
            if len(rms_values) > 1:
                rms_array = np.array(rms_values)
                median_rest = np.median(rms_array[1:])
                threshold = 10.0  # Use 10x threshold for RMS outlier detection

                if median_rest > 0 and abs(rms_array[0]) > threshold * abs(median_rest):
                    skip_first_rms = True

            if skip_first_rms:
                note = f"\nNote: First point excluded for clarity (RMS₀={rms_values[0]:.2e})"
                self.console.print(note, highlight=False)
                plot_rms = rms_values[1:]
                plot_steps_rms = steps[1:]
            else:
                plot_rms = rms_values
                plot_steps_rms = steps

            self.console.print("\n")
            rms_plot = self._create_ascii_plot(
                plot_rms,
                title="RMS Gradient vs Step (lower is better)",
                ylabel="RMS Gradient",
                x_values=plot_steps_rms
            )
            self.console.print(rms_plot, highlight=False)

    def _assess_heating(self, monitor: AMBERMonitor, stage: str):
        """Quality assessment for heating. Reports raw metrics without
        arbitrary convergence thresholds; the user judges adequacy."""
        from rich.table import Table
        import numpy as np

        # Extract data
        temps = list(monitor.data.get('temperature', []))
        times = list(monitor.data.get('time', []))
        etot = list(monitor.data.get('total_energy', []))
        volumes = list(monitor.data.get('volume', []))

        if not temps or not times:
            self.console.print("[red]Error: Insufficient temperature data for heating analysis[/red]")
            return

        # Get target temperature from mdin parameters
        target_temp = monitor.mdin_params.get('temp0', 300.0)
        is_nvt = 'nvt' in stage.lower()

        # Calculate metrics
        results = {}
        results['total_steps'] = int(times[-1] / 0.002) if times else 0  # Assume 2fs timestep
        results['total_time_ps'] = times[-1] if times else 0
        results['initial_temp'] = temps[0]
        results['final_temp'] = temps[-1]

        # Analyze last 20% of simulation for equilibration check
        temp_array = np.array(temps)
        last_20_idx = max(1, int(0.8 * len(temp_array)))
        final_temps = temp_array[last_20_idx:]

        results['avg_final_temp'] = np.mean(final_temps)
        results['std_final_temp'] = np.std(final_temps)

        # Energy drift analysis
        if len(etot) > 1:
            etot_array = np.array(etot)
            results['initial_energy'] = etot[0]
            results['final_energy'] = etot[-1]
            results['energy_drift'] = etot[-1] - etot[0]

        # Volume variation check (should be constant for NVT)
        if is_nvt and len(volumes) > 1:
            vol_array = np.array(volumes)
            results['volume_variation_pct'] = (np.std(vol_array) / np.mean(vol_array)) * 100

        # Display results table
        table = Table(title="Heating Metrics", show_header=True)
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        # Simulation info
        table.add_row("Total Time", f"{results['total_time_ps']:.2f} ps")

        # Temperature ramping
        table.add_row(
            "Initial Temperature",
            f"{results['initial_temp']:.1f} K"
        )
        table.add_row(
            "Final Temperature",
            f"{results['final_temp']:.1f} K"
        )
        table.add_row(
            "Avg Final Temp (last 20%)",
            f"{results['avg_final_temp']:.1f} ± {results['std_final_temp']:.1f} K (target: {target_temp:.1f} K)"
        )

        # Energy
        if 'energy_drift' in results:
            table.add_row(
                "Energy Drift",
                f"{results['energy_drift']:+.2f} kcal/mol"
            )

        # Volume check for NVT
        if is_nvt and 'volume_variation_pct' in results:
            table.add_row(
                "Volume Variation (NVT)",
                f"{results['volume_variation_pct']:.3f}%"
            )

        self.console.print(table)

        # ASCII Plots
        self.console.print("\n" + "="*70)
        self.console.print(f"{stage.replace('_', ' ').upper()} - ASCII VISUALIZATION")
        self.console.print("="*70)

        # Temperature vs Time plot
        temp_plot = self._create_ascii_plot(
            temps,
            title="Temperature vs Time",
            ylabel="Temperature (K)",
            x_values=times
        )
        self.console.print(temp_plot, highlight=False)

        # Total Energy vs Time plot
        if etot:
            self.console.print("\n")
            energy_plot = self._create_ascii_plot(
                etot,
                title="Total Energy vs Time",
                ylabel="Energy (kcal/mol)",
                x_values=times
            )
            self.console.print(energy_plot, highlight=False)

        # Temperature distribution histogram for last 20%
        if len(final_temps) > 5:
            self.console.print("\n")
            hist_plot = self._create_ascii_histogram(
                final_temps.tolist(),
                title="Temperature Distribution (final 20%)",
                xlabel="Temperature (K)"
            )
            self.console.print(hist_plot, highlight=False)

    def _create_ascii_plot(self, values: list, title: str = "", ylabel: str = "",
                           xlabel: str = "", x_values: list = None, width: int = 70, height: int = 20) -> str:
        """
        Create an ASCII line plot using the QC script's superior design.

        Args:
            values: Data values to plot
            title: Plot title
            ylabel: Y-axis label
            xlabel: X-axis label
            x_values: Optional x-axis values (e.g., steps, time)
            width: Plot width in characters
            height: Plot height in characters

        Returns:
            ASCII plot as string
        """
        import numpy as np

        if not values or len(values) < 2:
            return "Insufficient data to plot"

        original_n = len(values)
        x_min = x_values[0] if x_values else 0
        x_max = x_values[-1] if x_values else original_n - 1

        # Downsample if needed
        if len(values) > width:
            step = len(values) // width
            values = [values[i] for i in range(0, len(values), step)]

        # Normalize data to plot height
        min_val = min(values)
        max_val = max(values)

        if max_val == min_val:
            return f"{title}: Constant value {max_val:.2f}"

        # Normalize to height
        normalized = []
        for v in values:
            norm_val = int(((v - min_val) / (max_val - min_val)) * (height - 1))
            normalized.append(norm_val)

        lines = []

        # Title
        if title:
            lines.append(f"\n{title}")
            lines.append("=" * len(title))

        # Plot rows (top to bottom)
        for row in range(height - 1, -1, -1):
            # Y-axis label
            if row == height - 1:
                line = f"{max_val:10.2e} ┤"
            elif row == 0:
                line = f"{min_val:10.2e} ┤"
            elif row == height // 2:
                line = f"{(min_val + max_val) / 2:10.2e} ┤"
            else:
                line = " " * 10 + " │"

            # Plot points
            for col_val in normalized:
                if col_val == row:
                    line += "●"
                elif col_val > row:
                    line += "│"
                else:
                    line += " "

            lines.append(line)

        # X-axis
        lines.append(" " * 11 + "└" + "─" * len(normalized))

        # X-axis labels
        if x_values:
            x_min_str = f"{int(x_min)}" if x_min >= 1 else f"{x_min:.1f}"
            x_max_str = f"{int(x_max)}" if x_max >= 1 else f"{x_max:.1f}"

            label_len = len(x_min_str) + len(x_max_str)
            available_space = len(normalized)

            if available_space > label_len + 2:
                spacing = available_space - label_len
                x_label = f"{x_min_str}{' ' * spacing}{x_max_str}"
            else:
                x_label = f"{x_min_str}...{x_max_str}"

            lines.append(" " * 11 + x_label)
            lines.append(" " * 11 + f"[{original_n} data points]")
        else:
            lines.append(" " * 11 + f"0{' ' * (len(normalized) - 10)}N={len(values)}")

        # X-label
        if xlabel:
            lines.append(f"\n{xlabel}")

        # Y-label
        if ylabel:
            lines.append(f"\n{ylabel}")

        # Statistics
        mean = np.mean(values)
        std = np.std(values)
        lines.append(f"\nStats: Mean={mean:.2e}, Std={std:.2e}, Min={min_val:.2e}, Max={max_val:.2e}")

        return "\n".join(lines)

    def _create_ascii_histogram(self, values: list, bins: int = 20, width: int = 50,
                                title: str = "", xlabel: str = "") -> str:
        """
        Create an ASCII histogram for distribution visualization.

        Args:
            values: Data values to plot
            bins: Number of bins
            width: Maximum bar width in characters
            title: Plot title
            xlabel: X-axis label

        Returns:
            ASCII histogram as string
        """
        import numpy as np

        if len(values) < 2:
            return "Insufficient data for histogram"

        hist, edges = np.histogram(values, bins=bins)
        max_count = max(hist)

        lines = []

        # Title
        if title:
            lines.append(f"\n{title}")
            lines.append("=" * len(title))

        # Plot bars
        for i, count in enumerate(hist):
            bar_len = int((count / max_count) * width) if max_count > 0 else 0
            bar = "█" * bar_len
            lines.append(f"{edges[i]:8.2f} │{bar} {count}")

        # X-axis label
        if xlabel:
            lines.append(f"\n{xlabel}")

        # Statistics
        mean = np.mean(values)
        std = np.std(values)
        lines.append(f"\nStats: Mean={mean:.2f}, Std={std:.2f}, N={len(values)}")

        return "\n".join(lines)

    def _assess_npt_equilibration(self, monitor: AMBERMonitor):
        """Quality assessment for NPT equilibration/production. Reports raw
        metrics without arbitrary convergence thresholds; the user judges adequacy."""
        from rich.table import Table
        import numpy as np

        # Extract data
        densities = list(monitor.data.get('density', []))
        temps = list(monitor.data.get('temperature', []))
        times = list(monitor.data.get('time', []))
        volumes = list(monitor.data.get('volume', []))
        etot = list(monitor.data.get('total_energy', []))
        pressures = list(monitor.data.get('pressure', []))
        ke = list(monitor.data.get('kinetic_energy', []))
        pe = list(monitor.data.get('potential_energy', []))

        if not densities or not times:
            self.console.print("[red]Error: Insufficient density data for NPT equilibration analysis[/red]")
            return

        # Get target values
        target_temp = monitor.mdin_params.get('temp0', 300.0)
        target_pressure = monitor.mdin_params.get('pres0', 1.0)  # bar

        # Calculate metrics
        results = {}
        results['total_time_ps'] = times[-1] if times else 0

        # Time-based window analysis (first 1ns vs last 1ns)
        time_array = np.array(times)
        density_array = np.array(densities)
        start_time = times[0]
        end_time = times[-1]
        total_duration = end_time - start_time

        if total_duration >= 2000:  # At least 2ns duration for meaningful comparison
            # First 1ns window (relative to start)
            first_1ns_mask = (time_array >= start_time) & (time_array <= start_time + 1000)
            first_1ns_densities = density_array[first_1ns_mask]

            # Last 1ns window (relative to end)
            last_1ns_mask = time_array >= (end_time - 1000)
            last_1ns_densities = density_array[last_1ns_mask]

            if len(first_1ns_densities) > 0 and len(last_1ns_densities) > 0:
                results['first_1ns_density_mean'] = np.mean(first_1ns_densities)
                results['first_1ns_density_std'] = np.std(first_1ns_densities)
                results['last_1ns_density_mean'] = np.mean(last_1ns_densities)
                results['last_1ns_density_std'] = np.std(last_1ns_densities)
                results['density_change'] = abs(results['first_1ns_density_mean'] - results['last_1ns_density_mean'])

                # Use last 1ns as target density and for histogram
                target_density = results['last_1ns_density_mean']
                final_densities = last_1ns_densities
            else:
                # Fall back if masks produce empty arrays
                last_50_idx = max(1, int(0.5 * len(density_array)))
                final_densities = density_array[last_50_idx:]
                target_density = np.mean(final_densities)
                results['avg_density'] = target_density
                results['std_density'] = np.std(final_densities)
        else:
            # Fall back to last 50% for short simulations
            last_50_idx = max(1, int(0.5 * len(density_array)))
            final_densities = density_array[last_50_idx:]
            target_density = np.mean(final_densities)
            results['avg_density'] = target_density
            results['std_density'] = np.std(final_densities)

        # Pressure analysis (first 1ns vs last 1ns, same as density)
        if pressures:
            pressure_array = np.array(pressures)
            if total_duration >= 2000 and len(pressure_array) > 0:
                # Use same time masks as density
                first_1ns_pressures = pressure_array[first_1ns_mask] if len(pressure_array) == len(time_array) else pressure_array[:len(first_1ns_mask)][first_1ns_mask[:len(pressure_array)]]
                last_1ns_pressures = pressure_array[last_1ns_mask] if len(pressure_array) == len(time_array) else pressure_array[last_1ns_mask[:len(pressure_array)]]

                if len(first_1ns_pressures) > 0 and len(last_1ns_pressures) > 0:
                    results['first_1ns_pressure_mean'] = np.mean(first_1ns_pressures)
                    results['first_1ns_pressure_std'] = np.std(first_1ns_pressures)
                    results['last_1ns_pressure_mean'] = np.mean(last_1ns_pressures)
                    results['last_1ns_pressure_std'] = np.std(last_1ns_pressures)
                    results['pressure_change'] = abs(results['first_1ns_pressure_mean'] - results['last_1ns_pressure_mean'])
                    final_pressures = last_1ns_pressures
                else:
                    # Fallback to last 50%
                    last_50_idx_press = max(1, int(0.5 * len(pressure_array)))
                    final_pressures = pressure_array[last_50_idx_press:]
                    results['avg_pressure'] = np.mean(final_pressures)
                    results['std_pressure'] = np.std(final_pressures)
            else:
                # Short simulation: use last 50%
                last_50_idx_press = max(1, int(0.5 * len(pressure_array)))
                final_pressures = pressure_array[last_50_idx_press:]
                results['avg_pressure'] = np.mean(final_pressures)
                results['std_pressure'] = np.std(final_pressures)

        # Temperature analysis (last 50%)
        if temps:
            temp_array = np.array(temps)
            last_50_idx_temp = max(1, int(0.5 * len(temp_array)))
            final_temps = temp_array[last_50_idx_temp:]

            results['avg_temp'] = np.mean(final_temps)
            results['std_temp'] = np.std(final_temps)

        # Volume analysis (last 50%)
        if volumes:
            vol_array = np.array(volumes)
            last_50_idx_vol = max(1, int(0.5 * len(vol_array)))
            final_volumes = vol_array[last_50_idx_vol:]

            results['avg_volume'] = np.mean(final_volumes)
            results['std_volume'] = np.std(final_volumes)
            results['volume_variation_pct'] = (np.std(final_volumes) / np.mean(final_volumes)) * 100

        # Energy analysis (last 50%)
        if etot:
            etot_array = np.array(etot)
            last_50_idx_etot = max(1, int(0.5 * len(etot_array)))
            final_energies = etot_array[last_50_idx_etot:]

            results['avg_total_energy'] = np.mean(final_energies)
            results['std_total_energy'] = np.std(final_energies)
            results['energy_drift'] = final_energies[-1] - final_energies[0]

        # Kinetic and Potential Energy components
        if ke:
            ke_array = np.array(ke)
            last_50_idx_ke = max(1, int(0.5 * len(ke_array)))
            final_ke = ke_array[last_50_idx_ke:]
            results['avg_kinetic_energy'] = np.mean(final_ke)
            results['std_kinetic_energy'] = np.std(final_ke)

        if pe:
            pe_array = np.array(pe)
            last_50_idx_pe = max(1, int(0.5 * len(pe_array)))
            final_pe = pe_array[last_50_idx_pe:]
            results['avg_potential_energy'] = np.mean(final_pe)
            results['std_potential_energy'] = np.std(final_pe)

        # Sampling information
        results['n_snapshots'] = len(times)
        if len(times) > 1:
            time_diffs = np.diff(time_array)
            results['snapshot_interval'] = np.median(time_diffs)  # Use median to handle irregularities

        # Display results table
        table = Table(title="NPT Equilibration/Production Metrics", show_header=True)
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        # Simulation info
        table.add_row("Total Time", f"{end_time:.2f} ps (duration: {total_duration:.2f} ps)")

        # Sampling information
        if 'n_snapshots' in results:
            table.add_row("Snapshots", f"{results['n_snapshots']}")
        if 'snapshot_interval' in results:
            table.add_row("Snapshot Interval", f"{results['snapshot_interval']:.1f} ps")

        # Density comparison (first 1ns vs last 1ns)
        if 'first_1ns_density_mean' in results:
            # Successfully analyzed with 1ns windows
            table.add_row("Analysis Window", "First 1ns vs Last 1ns")
            table.add_row(
                "Density (first 1ns)",
                f"{results['first_1ns_density_mean']:.4f} ± {results['first_1ns_density_std']:.4f} g/cm³"
            )
            table.add_row(
                "Density (last 1ns)",
                f"{results['last_1ns_density_mean']:.4f} ± {results['last_1ns_density_std']:.4f} g/cm³"
            )
            table.add_row(
                "Density Change (|Δρ|)",
                f"{results['density_change']:.4f} g/cm³"
            )
        else:
            # Short simulation or continuation - use last 50%
            reason = "simulation < 2ns" if total_duration < 2000 else "continuation run"
            table.add_row("Analysis Window", f"Last 50% ({reason})")
            table.add_row(
                "Avg Density (final 50%)",
                f"{results['avg_density']:.4f} ± {results['std_density']:.4f} g/cm³"
            )

        # Pressure comparison (first 1ns vs last 1ns)
        if 'first_1ns_pressure_mean' in results:
            table.add_row(
                "Pressure (first 1ns)",
                f"{results['first_1ns_pressure_mean']:.2f} ± {results['first_1ns_pressure_std']:.2f} bar"
            )
            table.add_row(
                "Pressure (last 1ns)",
                f"{results['last_1ns_pressure_mean']:.2f} ± {results['last_1ns_pressure_std']:.2f} bar (target: {target_pressure:.1f} bar)"
            )
            table.add_row(
                "Pressure Change (|ΔP|)",
                f"{results['pressure_change']:.2f} bar"
            )
        elif 'avg_pressure' in results:
            table.add_row(
                "Pressure (final 50%)",
                f"{results['avg_pressure']:.2f} ± {results['std_pressure']:.2f} bar (target: {target_pressure:.1f} bar)"
            )

        # Temperature
        if 'avg_temp' in results:
            table.add_row(
                "Temperature (final 50%)",
                f"{results['avg_temp']:.1f} ± {results['std_temp']:.1f} K (target: {target_temp:.1f} K)"
            )

        # Volume
        if 'avg_volume' in results:
            table.add_row(
                "Volume (final 50%)",
                f"{results['avg_volume']:.1f} ± {results['std_volume']:.1f} A^3 (σ/μ = {results['volume_variation_pct']:.2f}%)"
            )

        # Energy
        if 'energy_drift' in results:
            table.add_row(
                "Total Energy (final 50%)",
                f"{results['avg_total_energy']:.2f} ± {results['std_total_energy']:.2f} kcal/mol"
            )
            table.add_row(
                "Energy Drift (final 50%)",
                f"{results['energy_drift']:+.2f} kcal/mol"
            )

        # Energy components (KE and PE)
        if 'avg_kinetic_energy' in results:
            table.add_row(
                "Kinetic Energy (final 50%)",
                f"{results['avg_kinetic_energy']:.2f} ± {results['std_kinetic_energy']:.2f} kcal/mol"
            )
        if 'avg_potential_energy' in results:
            table.add_row(
                "Potential Energy (final 50%)",
                f"{results['avg_potential_energy']:.2f} ± {results['std_potential_energy']:.2f} kcal/mol"
            )

        self.console.print(table)

        # ASCII Plots
        self.console.print("\n" + "="*70)
        self.console.print("NPT EQUILIBRATION - ASCII VISUALIZATION")
        self.console.print("="*70)

        # Density vs Time plot
        density_plot = self._create_ascii_plot(
            densities,
            title="Density vs Time",
            ylabel="Density (g/cm³)",
            x_values=times
        )
        self.console.print(density_plot, highlight=False)

        # Temperature vs Time plot
        if temps:
            self.console.print("\n")
            temp_plot = self._create_ascii_plot(
                temps,
                title="Temperature vs Time",
                ylabel="Temperature (K)",
                x_values=times
            )
            self.console.print(temp_plot, highlight=False)

        # Volume vs Time plot
        if volumes:
            self.console.print("\n")
            volume_plot = self._create_ascii_plot(
                volumes,
                title="System Volume vs Time",
                ylabel="Volume (A^3)",
                x_values=times
            )
            self.console.print(volume_plot, highlight=False)

        # Pressure vs Time plot
        if pressures:
            self.console.print("\n")
            pressure_plot = self._create_ascii_plot(
                pressures,
                title="Pressure vs Time",
                ylabel="Pressure (bar)",
                x_values=times[:len(pressures)]  # Ensure same length
            )
            self.console.print(pressure_plot, highlight=False)

        # Density distribution histogram for last 50%
        if len(final_densities) > 5:
            self.console.print("\n")
            hist_plot = self._create_ascii_histogram(
                final_densities.tolist(),
                title="Density Distribution (final 50%)",
                xlabel="Density (g/cm³)"
            )
            self.console.print(hist_plot, highlight=False)

    def _assess_nvt_equilibration_production(self, monitor: AMBERMonitor, stage: str):
        """Quality assessment for NVT equilibration/production. Reports raw
        metrics without arbitrary convergence thresholds; the user judges adequacy."""
        from rich.table import Table
        import numpy as np

        self.console.print("[bold]NVT Equilibration/Production Metrics[/bold]\n")

        # Get data
        times = list(monitor.data.get('time', []))
        temperatures = list(monitor.data.get('temperature', []))
        energies = list(monitor.data.get('total_energy', []))
        volumes = list(monitor.data.get('volume', []))
        ke = list(monitor.data.get('kinetic_energy', []))
        pe = list(monitor.data.get('potential_energy', []))

        if not times or not temperatures:
            self.console.print("[red]Insufficient data for NVT assessment[/red]")
            return

        results = {}

        # Basic metrics
        time_array = np.array(times)
        temp_array = np.array(temperatures)
        energy_array = np.array(energies) if energies else None
        volume_array = np.array(volumes) if volumes else None
        ke_array = np.array(ke) if ke else None
        pe_array = np.array(pe) if pe else None

        start_time = times[0]
        end_time = times[-1]
        total_duration = end_time - start_time

        results['start_time'] = start_time
        results['end_time'] = end_time
        results['total_duration'] = total_duration

        # Time-based window analysis (first 1ns vs last 1ns)
        if total_duration >= 2000:  # At least 2ns for meaningful comparison
            # First 1ns window
            first_1ns_mask = (time_array >= start_time) & (time_array <= start_time + 1000)
            first_1ns_temps = temp_array[first_1ns_mask]

            # Last 1ns window
            last_1ns_mask = time_array >= (end_time - 1000)
            last_1ns_temps = temp_array[last_1ns_mask]

            if len(first_1ns_temps) > 0 and len(last_1ns_temps) > 0:
                results['first_1ns_temp_mean'] = np.mean(first_1ns_temps)
                results['first_1ns_temp_std'] = np.std(first_1ns_temps)
                results['last_1ns_temp_mean'] = np.mean(last_1ns_temps)
                results['last_1ns_temp_std'] = np.std(last_1ns_temps)
                results['temp_change'] = abs(results['first_1ns_temp_mean'] - results['last_1ns_temp_mean'])

                # Use last 1ns for final statistics
                target_temp = results['last_1ns_temp_mean']
                final_temps = last_1ns_temps
            else:
                # Fallback: not enough data for 1ns windows
                results['insufficient_data'] = True
                target_temp = np.mean(temp_array)
                final_temps = temp_array
        else:
            # Short simulation: use all data
            results['short_simulation'] = True
            target_temp = np.mean(temp_array)
            final_temps = temp_array

        # Final temperature statistics (last 1ns or all data)
        results['avg_temp'] = np.mean(final_temps)
        results['std_temp'] = np.std(final_temps)
        results['min_temp'] = np.min(final_temps)
        results['max_temp'] = np.max(final_temps)

        # Energy statistics (if available)
        if energy_array is not None and len(energies) > 0:
            # Use same time windows for energy
            if total_duration >= 2000 and 'first_1ns_temp_mean' in results:
                # Get corresponding energies for last 1ns
                final_energies = energy_array[last_1ns_mask]
                results['avg_energy'] = np.mean(final_energies)
                results['std_energy'] = np.std(final_energies)
                results['energy_drift'] = final_energies[-1] - final_energies[0] if len(final_energies) > 1 else 0.0
            else:
                # Use all energy data
                results['avg_energy'] = np.mean(energy_array)
                results['std_energy'] = np.std(energy_array)
                results['energy_drift'] = energy_array[-1] - energy_array[0] if len(energy_array) > 1 else 0.0

        # Volume constancy check (NVT should have constant volume)
        if volume_array is not None and len(volumes) > 0:
            results['avg_volume'] = np.mean(volume_array)
            results['std_volume'] = np.std(volume_array)
            results['min_volume'] = np.min(volume_array)
            results['max_volume'] = np.max(volume_array)
            # Volume variation coefficient (σ/μ)
            results['volume_variation_pct'] = (results['std_volume'] / results['avg_volume'] * 100) if results['avg_volume'] > 0 else 0.0

        # Kinetic and Potential Energy components
        if ke_array is not None and len(ke) > 0:
            # Use same time windows as temperature
            if total_duration >= 2000 and 'first_1ns_temp_mean' in results:
                final_ke = ke_array[last_1ns_mask] if len(ke_array) == len(time_array) else ke_array
                results['avg_kinetic_energy'] = np.mean(final_ke)
                results['std_kinetic_energy'] = np.std(final_ke)
            else:
                results['avg_kinetic_energy'] = np.mean(ke_array)
                results['std_kinetic_energy'] = np.std(ke_array)

        if pe_array is not None and len(pe) > 0:
            # Use same time windows as temperature
            if total_duration >= 2000 and 'first_1ns_temp_mean' in results:
                final_pe = pe_array[last_1ns_mask] if len(pe_array) == len(time_array) else pe_array
                results['avg_potential_energy'] = np.mean(final_pe)
                results['std_potential_energy'] = np.std(final_pe)
            else:
                results['avg_potential_energy'] = np.mean(pe_array)
                results['std_potential_energy'] = np.std(pe_array)

        # Sampling information
        results['n_snapshots'] = len(times)
        if len(times) > 1:
            time_diffs = np.diff(time_array)
            results['snapshot_interval'] = np.median(time_diffs)

        # Display results table
        table = Table(title="NVT Equilibration/Production Metrics", show_header=True)
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        # Simulation info
        table.add_row("Total Time", f"{end_time:.2f} ps (duration: {total_duration:.2f} ps)")

        # Sampling information
        if 'n_snapshots' in results:
            table.add_row("Snapshots", f"{results['n_snapshots']}")
        if 'snapshot_interval' in results:
            table.add_row("Snapshot Interval", f"{results['snapshot_interval']:.1f} ps")

        # Temperature comparison (first 1ns vs last 1ns)
        if 'first_1ns_temp_mean' in results:
            table.add_row("Analysis Window", "First 1ns vs Last 1ns")
            table.add_row(
                "Temperature (first 1ns)",
                f"{results['first_1ns_temp_mean']:.2f} ± {results['first_1ns_temp_std']:.2f} K"
            )
            table.add_row(
                "Temperature (last 1ns)",
                f"{results['last_1ns_temp_mean']:.2f} ± {results['last_1ns_temp_std']:.2f} K"
            )
            table.add_row(
                "Temperature Change (|ΔT|)",
                f"{results['temp_change']:.2f} K"
            )
        elif 'short_simulation' in results:
            table.add_row("Analysis Window", "All data (simulation < 2ns)")
        elif 'insufficient_data' in results:
            table.add_row("Analysis Window", "All data (insufficient points for 1ns windows)")

        # Temperature statistics
        table.add_row(
            "Temperature (final)",
            f"{results['avg_temp']:.2f} ± {results['std_temp']:.2f} K (range: {results['min_temp']:.1f} - {results['max_temp']:.1f} K)"
        )

        # Energy statistics
        if 'avg_energy' in results:
            table.add_row(
                "Total Energy (final)",
                f"{results['avg_energy']:.2f} ± {results['std_energy']:.2f} kcal/mol"
            )
            table.add_row(
                "Energy Drift",
                f"{results['energy_drift']:.2f} kcal/mol"
            )

        # Energy components (KE and PE)
        if 'avg_kinetic_energy' in results:
            table.add_row(
                "Kinetic Energy (final)",
                f"{results['avg_kinetic_energy']:.2f} ± {results['std_kinetic_energy']:.2f} kcal/mol"
            )
        if 'avg_potential_energy' in results:
            table.add_row(
                "Potential Energy (final)",
                f"{results['avg_potential_energy']:.2f} ± {results['std_potential_energy']:.2f} kcal/mol"
            )

        # Volume constancy (NVT should have constant volume)
        if 'avg_volume' in results:
            table.add_row(
                "Volume",
                f"{results['avg_volume']:.1f} ± {results['std_volume']:.1f} A^3 (σ/μ = {results['volume_variation_pct']:.3f}%)"
            )
            table.add_row(
                "Volume Range",
                f"{results['min_volume']:.1f} - {results['max_volume']:.1f} A^3"
            )

        self.console.print(table)

        # Generate plots
        self.console.print("\n[bold]Trajectory Plots:[/bold]\n")

        # Plot 1: Temperature vs Time
        if len(temperatures) > 1:
            plot_data = self._create_ascii_plot(
                temperatures,
                title="Temperature vs Time",
                ylabel="Temperature (K)",
                x_values=times,
                width=80,
                height=15
            )
            self.console.print(plot_data)

        # Plot 2: Energy vs Time
        if len(energies) > 1:
            plot_data = self._create_ascii_plot(
                energies,
                title="Total Energy vs Time",
                ylabel="Energy (kcal/mol)",
                x_values=times[:len(energies)],  # Ensure same length
                width=80,
                height=15
            )
            self.console.print(plot_data)

        # Plot 3: Volume vs Time (should be constant for NVT)
        if len(volumes) > 1:
            plot_data = self._create_ascii_plot(
                volumes,
                title="Volume vs Time (should be constant for NVT)",
                ylabel="Volume (A^3)",
                x_values=times[:len(volumes)],  # Ensure same length
                width=80,
                height=15
            )
            self.console.print(plot_data)

        # Plot 4: Temperature distribution
        if len(final_temps) > 10:
            hist_plot = self._create_ascii_histogram(
                final_temps,
                title=f"Temperature Distribution (final data, mean={results['avg_temp']:.2f} K)",
                xlabel="Temperature (K)",
                width=80
            )
            self.console.print(hist_plot, highlight=False)

    def _show_generic_analysis(self, monitor: AMBERMonitor):
        """Show generic analysis when stage-specific assessment is not available."""
        from rich.table import Table

        self.console.print("\n[bold]Generic Statistics:[/bold]")

        # Energy statistics
        if monitor.data.get('total_energy'):
            energies = list(monitor.data['total_energy'])
            table = Table(show_header=True, header_style="bold bright_blue")
            table.add_column("Property", style="bright_blue")
            table.add_column("Value", style="white")

            table.add_row("Mean Energy", f"{sum(energies)/len(energies):.2f} kcal/mol")
            table.add_row("Min Energy", f"{min(energies):.2f} kcal/mol")
            table.add_row("Max Energy", f"{max(energies):.2f} kcal/mol")

            if monitor.data.get('temperature'):
                temps = list(monitor.data['temperature'])
                table.add_row("Mean Temperature", f"{sum(temps)/len(temps):.2f} K")

            self.console.print(table)

    def _show_analysis_interface(self, monitor: AMBERMonitor, sim_dir: Path, mdout_file: Path):
        """Show interactive analysis interface."""
        while True:
            self.console.print(f"\n[bold]Analysis Options for {sim_dir.name}:[/bold]")
            self.console.print("1. Simulation overview", highlight=False)
            self.console.print("2. Energy analysis", highlight=False)
            self.console.print("3. Temperature analysis", highlight=False)
            self.console.print("4. Pressure analysis", highlight=False)
            self.console.print("5. Trajectory information", highlight=False)
            self.console.print("6. Performance summary", highlight=False)
            self.console.print("7. Export analysis data", highlight=False)
            self.console.print("8. ← Back", highlight=False)
            
            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1","2","3","4","5","6","7","8"],
                default="1",
                module="MD Manager - Analysis",
                description="Select analysis option",
                options_map={
                    "1": "Complete simulation overview",
                    "2": "Energy analysis",
                    "3": "Temperature analysis",
                    "4": "Pressure analysis",
                    "5": "Trajectory information",
                    "6": "Performance summary",
                    "7": "Export analysis data",
                    "8": "← Back"
                }
            )
            
            if choice == "1":
                self._show_simulation_overview(monitor, sim_dir, mdout_file)
            elif choice == "2":
                self._show_detailed_energy_analysis(monitor)
            elif choice == "3":
                self._show_detailed_temperature_analysis(monitor)
            elif choice == "4":
                self._show_detailed_pressure_analysis(monitor)
            elif choice == "5":
                self._show_trajectory_information(sim_dir)
            elif choice == "6":
                self._show_performance_summary(monitor, mdout_file)
            elif choice == "7":
                self._export_analysis_data(monitor, sim_dir)
            elif choice == "8":
                break

    def _show_simulation_overview(self, monitor: AMBERMonitor, sim_dir: Path, mdout_file: Path):
        """Show comprehensive simulation overview."""
        self.console.print(f"\n[bold blue]Simulation Overview[/bold blue]")
        
        from rich.table import Table
        overview_table = Table(title=f"Analysis: {sim_dir.name}")
        overview_table.add_column("Property", style="bright_blue")
        overview_table.add_column("Value", style="green")
        overview_table.add_column("Details", style="grey50")
        
        # Basic information
        overview_table.add_row("Output File", mdout_file.name, f"Size: {mdout_file.stat().st_size / 1024:.1f} KB")
        
        # Simulation progress
        if monitor.data.get('step'):
            total_steps = len(monitor.data['step'])
            final_step = monitor.data['step'][-1] if monitor.data['step'] else 0
            final_time = monitor.data['time'][-1] if monitor.data['time'] else 0
            overview_table.add_row("Total Steps", f"{final_step:,}", f"{total_steps} data points")
            overview_table.add_row("Final Time", f"{final_time:.3f} ps", f"{final_time/1000:.3f} ns")
        
        # Energy information
        if monitor.data.get('total_energy'):
            energies = monitor.data['total_energy']
            avg_energy = sum(energies) / len(energies)
            energy_stability = max(energies) - min(energies)
            overview_table.add_row("Average Energy", f"{avg_energy:.2f} kcal/mol", 
                                 f"Range: {energy_stability:.2f}")
        
        # Temperature information
        if monitor.data.get('temperature'):
            temps = monitor.data['temperature']
            avg_temp = sum(temps) / len(temps)
            temp_stability = max(temps) - min(temps)
            overview_table.add_row("Average Temperature", f"{avg_temp:.2f} K", 
                                 f"Range: {temp_stability:.2f}")
        
        # File information
        nc_files = list(sim_dir.glob("*.nc"))
        rst_files = list(sim_dir.glob("*.rst7"))
        if nc_files:
            trajectory_size = sum(f.stat().st_size for f in nc_files) / (1024**2)
            overview_table.add_row("Trajectory Files", f"{len(nc_files)} files", 
                                 f"Total: {trajectory_size:.1f} MB")
        if rst_files:
            overview_table.add_row("Restart Files", f"{len(rst_files)} files", "Available for continuation")
        
        self.console.print(overview_table)

    def _show_detailed_energy_analysis(self, monitor: AMBERMonitor):
        """Show detailed energy analysis with trends."""
        self.console.print(f"\n[bold blue]Detailed Energy Analysis[/bold blue]")
        
        if not monitor.data.get('total_energy'):
            self.console.print("[yellow]No energy data available[/yellow]")
            return
        
        energies = monitor.data['total_energy']
        steps = monitor.data.get('step', list(range(len(energies))))
        
        # Statistical analysis
        avg_energy = sum(energies) / len(energies)
        min_energy = min(energies)
        max_energy = max(energies)
        std_dev = (sum((e - avg_energy) ** 2 for e in energies) / len(energies)) ** 0.5
        
        from rich.table import Table
        stats_table = Table(title="Energy Statistics")
        stats_table.add_column("Statistic", style="bright_blue")
        stats_table.add_column("Value (kcal/mol)", style="green")
        
        stats_table.add_row("Mean", f"{avg_energy:.2f}")
        stats_table.add_row("Standard Deviation", f"{std_dev:.2f}")
        stats_table.add_row("Minimum", f"{min_energy:.2f}")
        stats_table.add_row("Maximum", f"{max_energy:.2f}")
        stats_table.add_row("Range", f"{max_energy - min_energy:.2f}")
        
        self.console.print(stats_table)
        
        # Energy stability analysis
        if len(energies) > 100:
            # Analyze last 25% for equilibration
            quarter_point = len(energies) * 3 // 4
            recent_energies = energies[quarter_point:]
            recent_avg = sum(recent_energies) / len(recent_energies)
            recent_std = (sum((e - recent_avg) ** 2 for e in recent_energies) / len(recent_energies)) ** 0.5
            
            equilibration_table = Table(title="Equilibration Analysis")
            equilibration_table.add_column("Period", style="bright_blue")
            equilibration_table.add_column("Mean Energy", style="green")
            equilibration_table.add_column("Std Dev", style="yellow")
            
            equilibration_table.add_row("Full Simulation", f"{avg_energy:.2f}", f"{std_dev:.2f}")
            equilibration_table.add_row("Final Quarter", f"{recent_avg:.2f}", f"{recent_std:.2f}")
            
            # Assess equilibration
            if recent_std < std_dev * 0.8:
                equilibration_table.add_row("Assessment", "[green]Well Equilibrated[/green]", "✓")
            else:
                equilibration_table.add_row("Assessment", "[yellow]May Need More Time[/yellow]", "⚠")
                
            self.console.print(equilibration_table)

    def _show_detailed_temperature_analysis(self, monitor: AMBERMonitor):
        """Show detailed temperature analysis."""
        self.console.print(f"\n[bold blue]Detailed Temperature Analysis[/bold blue]")
        
        if not monitor.data.get('temperature'):
            self.console.print("[yellow]No temperature data available[/yellow]")
            return
            
        temps = monitor.data['temperature']
        
        # Temperature statistics
        avg_temp = sum(temps) / len(temps)
        min_temp = min(temps)
        max_temp = max(temps)
        temp_std = (sum((t - avg_temp) ** 2 for t in temps) / len(temps)) ** 0.5
        
        from rich.table import Table
        temp_table = Table(title="Temperature Statistics")
        temp_table.add_column("Statistic", style="bright_blue")
        temp_table.add_column("Value (K)", style="green")
        
        temp_table.add_row("Mean", f"{avg_temp:.2f}")
        temp_table.add_row("Standard Deviation", f"{temp_std:.2f}")
        temp_table.add_row("Minimum", f"{min_temp:.2f}")
        temp_table.add_row("Maximum", f"{max_temp:.2f}")
        temp_table.add_row("Range", f"{max_temp - min_temp:.2f}")
        
        self.console.print(temp_table)
        
        # Temperature control assessment
        if temp_std < 5:
            self.console.print("[green]✓ Excellent temperature control[/green]")
        elif temp_std < 15:
            self.console.print("[yellow]○ Good temperature control[/yellow]")
        else:
            self.console.print("[red]⚠ Poor temperature control - check thermostat settings[/red]")

    def _show_detailed_pressure_analysis(self, monitor: AMBERMonitor):
        """Show detailed pressure analysis."""
        self.console.print(f"\n[bold blue]Detailed Pressure Analysis[/bold blue]")
        
        if not monitor.data.get('pressure'):
            self.console.print("[yellow]No pressure data available[/yellow]")
            return
            
        pressures = monitor.data['pressure']
        avg_pressure = sum(pressures) / len(pressures)
        pressure_std = (sum((p - avg_pressure) ** 2 for p in pressures) / len(pressures)) ** 0.5
        
        from rich.table import Table
        pressure_table = Table(title="Pressure Statistics")
        pressure_table.add_column("Statistic", style="bright_blue")
        pressure_table.add_column("Value (bar)", style="green")
        
        pressure_table.add_row("Mean", f"{avg_pressure:.2f}")
        pressure_table.add_row("Standard Deviation", f"{pressure_std:.2f}")
        pressure_table.add_row("Minimum", f"{min(pressures):.2f}")
        pressure_table.add_row("Maximum", f"{max(pressures):.2f}")
        
        self.console.print(pressure_table)

    def _show_trajectory_information(self, sim_dir: Path):
        """Show trajectory file information."""
        self.console.print(f"\n[bold blue]Trajectory Information[/bold blue]")
        
        nc_files = list(sim_dir.glob("*.nc"))
        if not nc_files:
            self.console.print("[yellow]No trajectory files (.nc) found[/yellow]")
            return
            
        from rich.table import Table
        traj_table = Table(title="Trajectory Files")
        traj_table.add_column("File", style="bright_blue")
        traj_table.add_column("Size (MB)", style="green")
        traj_table.add_column("Modified", style="grey50")
        
        for nc_file in nc_files:
            size_mb = nc_file.stat().st_size / (1024**2)
            mod_time = datetime.fromtimestamp(nc_file.stat().st_mtime)
            traj_table.add_row(nc_file.name, f"{size_mb:.1f}", mod_time.strftime("%Y-%m-%d %H:%M"))
            
        self.console.print(traj_table)

    def _show_performance_summary(self, monitor: AMBERMonitor, mdout_file: Path):
        """Show simulation performance summary."""
        self.console.print(f"\n[bold blue]Performance Summary[/bold blue]")
        
        # Try to extract timing information from mdout file
        try:
            with open(mdout_file, 'r') as f:
                content = f.read()
                
            # Look for timing information in AMBER output
            timing_info = {}
            
            # Look for wall clock time
            wall_time_match = re.search(r'Total wall time:\s+([\d.]+)', content)
            if wall_time_match:
                timing_info["Wall Time"] = f"{float(wall_time_match.group(1)):.1f} seconds"
                
            # Look for performance metrics
            ns_per_day_match = re.search(r'([\d.]+)\s+ns/day', content)
            if ns_per_day_match:
                timing_info["Performance"] = f"{float(ns_per_day_match.group(1)):.2f} ns/day"
                
            if timing_info:
                from rich.table import Table
                perf_table = Table(title="Performance Metrics")
                perf_table.add_column("Metric", style="bright_blue")
                perf_table.add_column("Value", style="green")
                
                for metric, value in timing_info.items():
                    perf_table.add_row(metric, value)
                    
                self.console.print(perf_table)
            else:
                self.console.print("[yellow]No performance metrics found in output file[/yellow]")
                
        except Exception as e:
            self.console.print(f"[red]Error reading performance data: {e}[/red]")

    def _export_analysis_data(self, monitor: AMBERMonitor, sim_dir: Path):
        """Export analysis data to files."""
        self.console.print(f"\n[bold blue]Export Analysis Data[/bold blue]")
        
        if not any(monitor.data.values()):
            self.console.print("[yellow]No data available to export[/yellow]")
            return
            
        export_dir = sim_dir / "analysis"
        export_dir.mkdir(exist_ok=True)
        
        try:
            # Export energy data
            if monitor.data.get('total_energy'):
                energy_file = export_dir / "energy_data.csv"
                with open(energy_file, 'w') as f:
                    f.write("Step,Time_ps,Total_Energy_kcal_mol\n")
                    steps = monitor.data.get('step', range(len(monitor.data['total_energy'])))
                    times = monitor.data.get('time', [0] * len(monitor.data['total_energy']))
                    
                    for step, time_ps, energy in zip(steps, times, monitor.data['total_energy']):
                        f.write(f"{step},{time_ps},{energy}\n")
                self.console.print(f"[green]✓ Energy data exported to {energy_file.name}[/green]")
                
            # Export temperature data
            if monitor.data.get('temperature'):
                temp_file = export_dir / "temperature_data.csv"
                with open(temp_file, 'w') as f:
                    f.write("Step,Time_ps,Temperature_K\n")
                    steps = monitor.data.get('step', range(len(monitor.data['temperature'])))
                    times = monitor.data.get('time', [0] * len(monitor.data['temperature']))
                    
                    for step, time_ps, temp in zip(steps, times, monitor.data['temperature']):
                        f.write(f"{step},{time_ps},{temp}\n")
                self.console.print(f"[green]✓ Temperature data exported to {temp_file.name}[/green]")
                
            # Export summary statistics
            summary_file = export_dir / "analysis_summary.txt"
            with open(summary_file, 'w') as f:
                f.write(f"AMBER Simulation Analysis Summary\n")
                f.write(f"=" * 40 + "\n\n")
                f.write(f"Simulation: {sim_dir.name}\n")
                f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                if monitor.data.get('total_energy'):
                    energies = monitor.data['total_energy']
                    f.write(f"Energy Statistics:\n")
                    f.write(f"  Average: {sum(energies)/len(energies):.2f} kcal/mol\n")
                    f.write(f"  Min: {min(energies):.2f} kcal/mol\n")
                    f.write(f"  Max: {max(energies):.2f} kcal/mol\n\n")
                    
                if monitor.data.get('temperature'):
                    temps = monitor.data['temperature']
                    f.write(f"Temperature Statistics:\n")
                    f.write(f"  Average: {sum(temps)/len(temps):.2f} K\n")
                    f.write(f"  Min: {min(temps):.2f} K\n")
                    f.write(f"  Max: {max(temps):.2f} K\n\n")
                    
            self.console.print(f"[green]✓ Analysis summary exported to {summary_file.name}[/green]")
            self.console.print(f"[grey50]All files saved in: {export_dir}[/grey50]")
            
        except Exception as e:
            self.console.print(f"[red]Error exporting data: {e}[/red]")

    def _browse_for_mdout_file(self, start_dir=None):
        """Browse for .mdout/.out output files (multi-select for batch).

        Thin wrapper over the shared file browser. Returns a list of Paths
        (single picks come back as a one-element list), None on cancel, or the
        single Path returned by the recursive `find` command.
        """
        from pathlib import Path
        from datetime import datetime
        from proprep.utils.file_browser import file_browser

        start = Path(start_dir) if start_dir else self.working_directory

        def _mtime_detail(p):
            try:
                return "modified " + datetime.fromtimestamp(os.path.getmtime(p)).strftime("%m/%d/%Y %H:%M")
            except OSError:
                return ""

        extra = {
            "find": ("Search recursively for .mdout/.out files",
                     lambda cur: self._find_mdout_files(Path(cur))),
        }
        return file_browser(
            directory=str(start),
            extensions=[".mdout", ".out"],
            console=self.console,
            processor=self.processor,
            multi=True,
            label="mdout file",
            entry_detail=_mtime_detail,
            path_factory=Path,
            extra_commands=extra,
            module="MD Manager - MDOUT Browser",
        )



    def _find_mdout_files(self, current_dir):
        """Search recursively for .mdout/.out files."""
        from pathlib import Path

        self.console.print(f"\n[bold cyan]Searching for .mdout/.out files in {current_dir}...[/bold cyan]")

        # Find all mdout files recursively (both .mdout and .out extensions)
        mdout_files = list(current_dir.rglob("*.mdout")) + list(current_dir.rglob("*.out"))

        if not mdout_files:
            self.console.print(f"[grey50]No .mdout/.out files found[/grey50]")
            return None

        self.console.print(f"[green]Found {len(mdout_files)} .mdout/.out files:[/green]")

        # Display found files
        for i, file_path in enumerate(mdout_files, 1):
            relative_path = file_path.relative_to(current_dir)
            mod_time = datetime.fromtimestamp(file_path.stat().st_mtime)
            date_str = mod_time.strftime("%m/%d/%Y %H:%M")
            self.console.print(f"  {i:2}. {relative_path} (modified {date_str})")

        # Let user select
        while True:
            prompt = f"\nSelect file (1-{len(mdout_files)}) or 'cancel': "
            choice = prompt_with_context(
                self.processor, prompt.strip().rstrip(':').strip(),
                default="cancel", module="MD Manager - File Search",
                description="Select from recursively-found mdout files",
            ).strip().lower()

            if choice == 'cancel':
                return None

            try:
                file_num = int(choice)
                if 1 <= file_num <= len(mdout_files):
                    selected_file = mdout_files[file_num - 1]
                    self.console.print(f"[green]Selected: {selected_file}[/green]")
                    return selected_file
                else:
                    self.console.print(f"[red]Invalid selection. Choose 1-{len(mdout_files)}[/red]")
            except ValueError:
                self.console.print("[red]Invalid input. Enter a number or 'cancel'[/red]")

    def _browse_for_analysis_files(self, start_dir=None):
        """Browse for analysis files (.mdout/.out, .nc, .prmtop/.parm7).

        Returns:
            dict with keys:
                'mdout': list of .mdout/.out files
                'nc': list of .nc files
                'prmtop': list of .prmtop/.parm7 files
            or None if cancelled
        """
        from pathlib import Path
        import os

        self.console.print(f"\n[bold cyan]═══════════════════════════════════════════════════════════════[/bold cyan]")
        self.console.print(f"[bold cyan]                    MD ANALYSIS FILE BROWSER                   [/bold cyan]")
        self.console.print(f"[bold cyan]═══════════════════════════════════════════════════════════════[/bold cyan]")
        self.console.print()
        self.console.print(f"[bold yellow]WHAT TO SELECT:[/bold yellow]")
        self.console.print(f"  [cyan]For Energetics Analysis:[/cyan]")
        self.console.print(f"    → Select .mdout/.out file(s) only")
        self.console.print()
        self.console.print(f"  [cyan]For Trajectory Analysis:[/cyan]")
        self.console.print(f"    → Select .nc file(s) AND .prmtop/.parm7 topology file")
        self.console.print(f"    → Or select .nc file(s) only (will prompt for topology)")
        self.console.print()
        self.console.print(f"  [cyan]For Combined Analysis:[/cyan]")
        self.console.print(f"    → Select .mdout/.out + .nc + .prmtop/.parm7 files together")
        self.console.print()
        self.console.print(
            "  [grey50]Files in different directories? Use 'add N' to stage a file, 'cd' to\n"
            "  another directory, 'add N' again, then 'done'. (Don't pick a topology\n"
            "  on its own — pick the trajectory and you'll be prompted for topology.)[/grey50]"
        )
        self.console.print(f"[bold cyan]═══════════════════════════════════════════════════════════════[/bold cyan]")
        self.console.print()

        # Start from provided directory or working directory
        if start_dir:
            current_dir = Path(start_dir)
        else:
            current_dir = self.working_directory

        # Staged selections persist across directory navigations so the user
        # can pick a trajectory in one directory and a topology in another.
        staged_mdout: list = []
        staged_nc: list = []
        staged_prmtop: list = []

        while True:
            self.console.print(f"\n[bold]Current Directory:[/bold] [cyan]{current_dir}[/cyan]")

            # List directories, mdout files, nc files, and topology files in current directory
            try:
                dirs = [d for d in current_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
                mdout_files = [f for f in current_dir.iterdir() if f.is_file() and f.suffix in ['.mdout', '.out']]
                nc_files = [f for f in current_dir.iterdir() if f.is_file() and f.suffix == '.nc']
                # Support both .prmtop and .parm7 (newer AMBER format)
                prmtop_files = [f for f in current_dir.iterdir() if f.is_file() and f.suffix in ['.prmtop', '.parm7']]

                all_items = []

                # Add parent directory if not at root
                if current_dir.parent != current_dir:
                    all_items.append(("📁 .. (parent)", current_dir.parent, "parent"))

                # Add subdirectories
                for d in sorted(dirs):
                    all_items.append((f"📁 {d.name}", d, "dir"))

                # Add mdout files
                for f in sorted(mdout_files):
                    mod_time = datetime.fromtimestamp(f.stat().st_mtime)
                    date_str = mod_time.strftime("%m/%d/%Y %H:%M")
                    all_items.append((f"📄 [mdout] {f.name} (modified {date_str})", f, "mdout"))

                # Add nc files
                for f in sorted(nc_files):
                    mod_time = datetime.fromtimestamp(f.stat().st_mtime)
                    date_str = mod_time.strftime("%m/%d/%Y %H:%M")
                    size_mb = f.stat().st_size / (1024 * 1024)
                    all_items.append((f"📄 [nc] {f.name} ({size_mb:.1f} MB, modified {date_str})", f, "nc"))

                # Add topology files (.prmtop or .parm7)
                for f in sorted(prmtop_files):
                    mod_time = datetime.fromtimestamp(f.stat().st_mtime)
                    date_str = mod_time.strftime("%m/%d/%Y %H:%M")
                    size_kb = f.stat().st_size / 1024
                    # Show the actual extension in the label
                    ext = f.suffix[1:]  # Remove the dot
                    all_items.append((f"📄 [{ext}] {f.name} ({size_kb:.1f} KB, modified {date_str})", f, "prmtop"))

                if all_items:
                    # Display grouped by type with blank lines between groups
                    current_index = 1

                    # First: directories (parent + subdirs)
                    for (display, path, item_type) in all_items:
                        if item_type in ["parent", "dir"]:
                            self.console.print(f"  [{current_index:2}] {display}")
                            current_index += 1

                    # Second: mdout files
                    mdout_items = [(d, p, t) for d, p, t in all_items if t == "mdout"]
                    if mdout_items:
                        self.console.print()  # Blank line before mdout group
                        for (display, path, item_type) in mdout_items:
                            self.console.print(f"  [{current_index:2}] {display}")
                            current_index += 1

                    # Third: nc files
                    nc_items = [(d, p, t) for d, p, t in all_items if t == "nc"]
                    if nc_items:
                        self.console.print()  # Blank line before nc group
                        for (display, path, item_type) in nc_items:
                            self.console.print(f"  [{current_index:2}] {display}")
                            current_index += 1

                    # Fourth: topology files
                    prmtop_items = [(d, p, t) for d, p, t in all_items if t == "prmtop"]
                    if prmtop_items:
                        self.console.print()  # Blank line before topology group
                        for (display, path, item_type) in prmtop_items:
                            self.console.print(f"  [{current_index:2}] {display}")
                            current_index += 1
                else:
                    self.console.print(f"  [grey50]No directories or analysis files found[/grey50]")

            except PermissionError:
                self.console.print("[red]Permission denied accessing directory[/red]")
                current_dir = current_dir.parent
                continue

            # Show staged selections so the user knows what they have so far
            if staged_mdout or staged_nc or staged_prmtop:
                total_staged = len(staged_mdout) + len(staged_nc) + len(staged_prmtop)
                self.console.print(f"\n[bold green]Staged selections ({total_staged}):[/bold green]")
                for f in staged_mdout:
                    self.console.print(f"  • [mdout] {f}")
                for f in staged_nc:
                    self.console.print(f"  • [nc] {f}")
                for f in staged_prmtop:
                    self.console.print(f"  • [prmtop] {f}")

            self.console.print("\n[bold]Commands:[/bold]")
            self.console.print("  [cyan]select N[/cyan]   - Select file by number (submit immediately)")
            self.console.print("  [cyan]select N,M,P[/cyan] - Select multiple files by numbers (submit immediately)")
            self.console.print("  [cyan]add N[/cyan]      - Stage file(s) (e.g. add 3) — keep browsing to add more from other directories")
            self.console.print("  [cyan]add N,M,P[/cyan]  - Stage multiple files at once")
            self.console.print("  [cyan]done[/cyan]       - Submit staged selections")
            self.console.print("  [cyan]clear[/cyan]      - Clear staged selections")
            self.console.print("  [cyan]cd N[/cyan]       - Change to directory by number")
            self.console.print("  [cyan]cd path[/cyan]    - Change to specific directory")
            self.console.print("  [cyan]cd ..[/cyan]      - Go up one level")
            self.console.print(f"  [cyan]find[/cyan]       - Search recursively for analysis files")
            self.console.print("  [cyan]q[/cyan]          - Cancel and return")

            # prompt_with_context (not raw input) so this staging browser is
            # captured by the session recorder like every other browser.
            command = prompt_with_context(
                self.processor,
                "Enter command",
                default="q",
                module="MD Manager - Analysis Browser",
                description="Stage or select MD analysis files",
            ).strip()

            def _parse_nums(s: str):
                """Parse '3' or '2,5,7' into a list of ints. Returns None on bad input."""
                try:
                    return [int(n.strip()) for n in s.split(',')] if s else None
                except ValueError:
                    return None

            def _collect_from_indices(nums):
                """Map item numbers to (mdout, nc, prmtop) lists, skipping dirs/invalid."""
                m, n, p = [], [], []
                for file_num in nums:
                    if 1 <= file_num <= len(all_items):
                        _, path, item_type = all_items[file_num - 1]
                        if item_type == "mdout":
                            m.append(path)
                        elif item_type == "nc":
                            n.append(path)
                        elif item_type == "prmtop":
                            p.append(path)
                        else:
                            self.console.print(f"[red]Item {file_num} is a directory - skipping[/red]")
                    else:
                        self.console.print(f"[red]Invalid item number {file_num} - skipping[/red]")
                return m, n, p

            if command in ("q", "exit"):  # 'exit' kept as a back-compat alias
                return None
            elif command == "done":
                if staged_mdout or staged_nc or staged_prmtop:
                    total = len(staged_mdout) + len(staged_nc) + len(staged_prmtop)
                    self.console.print(f"[green]Submitting {total} staged file(s)[/green]")
                    return {'mdout': staged_mdout, 'nc': staged_nc, 'prmtop': staged_prmtop}
                self.console.print("[yellow]Nothing staged. Use 'add N' to stage files first, or 'select N' to pick a single file.[/yellow]")
            elif command == "clear":
                staged_mdout.clear()
                staged_nc.clear()
                staged_prmtop.clear()
                self.console.print("[green]Cleared staged selections[/green]")
            elif command.startswith("add "):
                nums = _parse_nums(command[4:].strip())
                if nums is None:
                    self.console.print("[red]Usage: add N or add N,M,P[/red]")
                else:
                    m, n, p = _collect_from_indices(nums)
                    if m or n or p:
                        # De-dup against already-staged files
                        for f in m:
                            if f not in staged_mdout:
                                staged_mdout.append(f)
                        for f in n:
                            if f not in staged_nc:
                                staged_nc.append(f)
                        for f in p:
                            if f not in staged_prmtop:
                                staged_prmtop.append(f)
                        added = len(m) + len(n) + len(p)
                        self.console.print(f"[green]Staged {added} file(s). Use 'cd' to browse other directories, then 'done' to submit.[/green]")
                    else:
                        self.console.print("[red]No valid files staged[/red]")
            elif command.startswith("select "):
                nums = _parse_nums(command[7:].strip())
                if nums is None:
                    self.console.print("[red]Usage: select N or select N,M,P[/red]")
                else:
                    m, n, p = _collect_from_indices(nums)
                    if m or n or p:
                        # Merge with anything already staged so we don't lose prior picks
                        out_m = staged_mdout + [f for f in m if f not in staged_mdout]
                        out_n = staged_nc + [f for f in n if f not in staged_nc]
                        out_p = staged_prmtop + [f for f in p if f not in staged_prmtop]
                        total = len(out_m) + len(out_n) + len(out_p)
                        self.console.print(f"[green]Selected {total} file(s):[/green]")
                        for f in out_m:
                            self.console.print(f"  • [mdout] {f.name}")
                        for f in out_n:
                            self.console.print(f"  • [nc] {f.name}")
                        for f in out_p:
                            self.console.print(f"  • [prmtop] {f.name}")
                        return {'mdout': out_m, 'nc': out_n, 'prmtop': out_p}
                    else:
                        self.console.print("[red]No valid files selected[/red]")
            elif command.startswith("cd "):
                path = command[3:].strip()
                if path == "..":
                    current_dir = current_dir.parent
                elif path.isdigit():
                    dir_num = int(path)
                    dir_items = [item for item in all_items if item[2] in ["parent", "dir"]]
                    if 1 <= dir_num <= len(dir_items):
                        current_dir = dir_items[dir_num - 1][1]
                        self.console.print(f"[green]Changed to: {current_dir}[/green]")
                    else:
                        self.console.print(f"[red]Invalid directory number. Choose 1-{len(dir_items)}[/red]")
                else:
                    new_dir = current_dir / path
                    if new_dir.exists() and new_dir.is_dir():
                        current_dir = new_dir
                        self.console.print(f"[green]Changed to: {current_dir}[/green]")
                    else:
                        self.console.print(f"[red]Directory not found: {path}[/red]")
            elif command == "find":
                found_files = self._find_analysis_files(current_dir)
                if found_files:
                    # Merge staged so users don't lose prior picks when they fall back to find
                    out_m = staged_mdout + [f for f in found_files.get('mdout', []) if f not in staged_mdout]
                    out_n = staged_nc + [f for f in found_files.get('nc', []) if f not in staged_nc]
                    out_p = staged_prmtop + [f for f in found_files.get('prmtop', []) if f not in staged_prmtop]
                    return {'mdout': out_m, 'nc': out_n, 'prmtop': out_p}
            else:
                self.console.print("[red]Unknown command[/red]")

    def _find_analysis_files(self, current_dir):
        """Search recursively for .mdout/.out, .nc, .prmtop, and .parm7 files.

        Returns:
            dict with keys 'mdout', 'nc', and 'prmtop' containing lists of files,
            or None if cancelled
        """
        from pathlib import Path

        self.console.print(f"\n[bold cyan]Searching for analysis files in {current_dir}...[/bold cyan]")

        # Find all analysis files recursively (both .mdout and .out extensions)
        mdout_files = list(current_dir.rglob("*.mdout")) + list(current_dir.rglob("*.out"))
        nc_files = list(current_dir.rglob("*.nc"))
        # Support both .prmtop and .parm7 (newer AMBER format)
        prmtop_files = list(current_dir.rglob("*.prmtop")) + list(current_dir.rglob("*.parm7"))

        if not mdout_files and not nc_files and not prmtop_files:
            self.console.print(f"[grey50]No analysis files found[/grey50]")
            return None

        # Combine and display all files
        all_files = []
        for f in mdout_files:
            all_files.append((f, "mdout"))
        for f in nc_files:
            all_files.append((f, "nc"))
        for f in prmtop_files:
            all_files.append((f, "prmtop"))

        self.console.print(f"[green]Found {len(all_files)} analysis file(s):[/green]")

        # Display found files
        for i, (file_path, file_type) in enumerate(all_files, 1):
            relative_path = file_path.relative_to(current_dir)
            mod_time = datetime.fromtimestamp(file_path.stat().st_mtime)
            date_str = mod_time.strftime("%m/%d/%Y %H:%M")
            if file_type == "nc":
                size_mb = file_path.stat().st_size / (1024 * 1024)
                self.console.print(f"  {i:2}. [{file_type}] {relative_path} ({size_mb:.1f} MB, modified {date_str})")
            elif file_type == "prmtop":
                size_kb = file_path.stat().st_size / 1024
                # Show actual extension (prmtop or parm7)
                ext = file_path.suffix[1:]
                self.console.print(f"  {i:2}. [{ext}] {relative_path} ({size_kb:.1f} KB, modified {date_str})")
            else:
                self.console.print(f"  {i:2}. [{file_type}] {relative_path} (modified {date_str})")

        # Let user select (support multiple selection)
        while True:
            choice = prompt_with_context(
                self.processor,
                f"Select file(s) (1-{len(all_files)}, comma-separated for multiple) or 'cancel'",
                default="cancel",
                module="MD Manager - Analysis Browser",
                description="Select from recursively-found analysis files",
            ).strip().lower()

            if choice == 'cancel':
                return None

            try:
                if ',' in choice:
                    # Multiple files
                    file_nums = [int(n.strip()) for n in choice.split(',')]
                    selected_mdout = []
                    selected_nc = []
                    selected_prmtop = []

                    for file_num in file_nums:
                        if 1 <= file_num <= len(all_files):
                            file_path, file_type = all_files[file_num - 1]
                            if file_type == "mdout":
                                selected_mdout.append(file_path)
                            elif file_type == "nc":
                                selected_nc.append(file_path)
                            else:
                                selected_prmtop.append(file_path)
                        else:
                            self.console.print(f"[red]Invalid file number {file_num} - skipping[/red]")

                    if selected_mdout or selected_nc or selected_prmtop:
                        total = len(selected_mdout) + len(selected_nc) + len(selected_prmtop)
                        self.console.print(f"[green]Selected {total} file(s)[/green]")
                        return {'mdout': selected_mdout, 'nc': selected_nc, 'prmtop': selected_prmtop}
                    else:
                        self.console.print("[red]No valid files selected[/red]")
                else:
                    # Single file
                    file_num = int(choice)
                    if 1 <= file_num <= len(all_files):
                        file_path, file_type = all_files[file_num - 1]
                        self.console.print(f"[green]Selected: {file_path}[/green]")
                        if file_type == "mdout":
                            return {'mdout': [file_path], 'nc': [], 'prmtop': []}
                        elif file_type == "nc":
                            return {'mdout': [], 'nc': [file_path], 'prmtop': []}
                        else:
                            return {'mdout': [], 'nc': [], 'prmtop': [file_path]}
                    else:
                        self.console.print(f"[red]Invalid selection. Choose 1-{len(all_files)}[/red]")
            except ValueError:
                self.console.print("[red]Invalid input. Enter number(s) or 'cancel'[/red]")

    def _route_analysis_by_file_type(self, selected_files: dict):
        """Route to appropriate analysis based on selected file types.

        Args:
            selected_files: dict with keys 'mdout', 'nc', and 'prmtop' containing file lists
        """
        mdout_files = selected_files.get('mdout', [])
        nc_files = selected_files.get('nc', [])
        prmtop_files = selected_files.get('prmtop', [])

        # If nc files selected but no prmtop, prompt for topology immediately
        if nc_files and not prmtop_files:
            self.console.print("\n[yellow]Trajectory files selected but no topology file (.prmtop/.parm7)[/yellow]")
            self.console.print("[yellow]Topology file is required for trajectory analysis[/yellow]")

            # Look for topology files in the same directory (both .prmtop and .parm7)
            parent_dir = nc_files[0].parent
            found_prmtop = sorted(list(parent_dir.glob("*.prmtop")) + list(parent_dir.glob("*.parm7")))

            if found_prmtop:
                self.console.print(f"\n[bold]Found {len(found_prmtop)} topology file(s) in same directory:[/bold]")
                for i, f in enumerate(found_prmtop, 1):
                    size_kb = f.stat().st_size / 1024
                    self.console.print(f"  {i}. {f.name} ({size_kb:.1f} KB)")

                if len(found_prmtop) == 1:
                    use_it = prompt_with_context(
                        self.processor,
                        f"Use {found_prmtop[0].name} as topology? (y/n)",
                        choices=["y", "n"],
                        default="y",
                        module="MD Manager - Trajectory Analysis",
                        description="Use found topology file",
                        options_map={"y": "Yes", "n": "No"}
                    )
                    if use_it == "y":
                        prmtop_files = [found_prmtop[0]]
                else:
                    choice = prompt_with_context(
                        self.processor,
                        f"Select topology file (1-{len(found_prmtop)}) or 'n' to browse elsewhere",
                        default="1",
                        module="MD Manager - Trajectory Analysis",
                        description="Select topology file"
                    )
                    # Replay by basename so a changed topology list can't mis-pick.
                    choice = remap_recorded_index(self.processor, found_prmtop, choice)
                    if choice.lower() != 'n':
                        try:
                            idx = int(choice) - 1
                            if 0 <= idx < len(found_prmtop):
                                prmtop_files = [found_prmtop[idx]]
                                annotate_selected_path(self.processor, found_prmtop[idx])
                        except ValueError:
                            pass

            # If still no prmtop, search recursively (down) AND walk up the tree.
            # Common layout: simulations/step6/step6.nc with topology one level up.
            if not prmtop_files:
                self.console.print("\n[bold]Searching for topology files...[/bold]")

                # Downward recursive search from the .nc directory
                found = list(parent_dir.rglob("*.prmtop")) + list(parent_dir.rglob("*.parm7"))

                # Upward search: scan each ancestor directory (non-recursive at each level)
                # up to 4 levels, so we don't scan unrelated sibling subtrees.
                ancestor = parent_dir.parent
                seen = {p.resolve() for p in found}
                for _ in range(4):
                    if ancestor == ancestor.parent:  # filesystem root
                        break
                    try:
                        ancestor_hits = list(ancestor.glob("*.prmtop")) + list(ancestor.glob("*.parm7"))
                    except (PermissionError, OSError):
                        ancestor_hits = []
                    for f in ancestor_hits:
                        if f.resolve() not in seen:
                            found.append(f)
                            seen.add(f.resolve())
                    ancestor = ancestor.parent

                all_prmtop = sorted(found)

                if all_prmtop:
                    self.console.print(f"[green]Found {len(all_prmtop)} topology file(s):[/green]")
                    for i, f in enumerate(all_prmtop, 1):
                        try:
                            relative = f.relative_to(parent_dir)
                        except ValueError:
                            relative = f  # ancestor dir — show absolute path
                        size_kb = f.stat().st_size / 1024
                        self.console.print(f"  {i}. {relative} ({size_kb:.1f} KB)")

                    while True:
                        choice = prompt_with_context(
                            self.processor,
                            f"Select topology file (1-{len(all_prmtop)}), 'browse' to pick manually, or 'cancel'",
                            default="cancel",
                            module="MD Manager - Trajectory Analysis",
                            description="Select topology file for trajectory analysis",
                        ).strip()
                        # Replay by basename so a changed topology list can't mis-pick.
                        choice = remap_recorded_index(self.processor, all_prmtop, choice)
                        if choice.lower() == 'cancel':
                            self.console.print("[yellow]Cannot analyze trajectory without topology - skipping trajectory analysis[/yellow]")
                            nc_files = []  # Clear nc files so we only do mdout analysis
                            break
                        if choice.lower() == 'browse':
                            browsed = self._browse_for_analysis_files(start_dir=parent_dir.parent)
                            if browsed and browsed.get('prmtop'):
                                prmtop_files = [browsed['prmtop'][0]]
                                self.console.print(f"[green]Selected: {prmtop_files[0].name}[/green]")
                                break
                            self.console.print("[yellow]No topology selected via browse[/yellow]")
                            continue
                        try:
                            idx = int(choice) - 1
                            if 0 <= idx < len(all_prmtop):
                                prmtop_files = [all_prmtop[idx]]
                                annotate_selected_path(self.processor, all_prmtop[idx])
                                self.console.print(f"[green]Selected: {prmtop_files[0].name}[/green]")
                                break
                            else:
                                self.console.print(f"[red]Invalid selection. Choose 1-{len(all_prmtop)}[/red]")
                        except ValueError:
                            self.console.print("[red]Invalid input[/red]")
                else:
                    # Nothing found automatically — let the user browse manually before giving up.
                    self.console.print("[yellow]No topology files found automatically.[/yellow]")
                    browse_choice = prompt_with_context(
                        self.processor,
                        "Browse manually for topology file? (y/n)",
                        choices=["y", "n"],
                        default="y",
                        module="MD Manager - Trajectory Analysis",
                        description="Browse manually for topology file",
                        options_map={"y": "Yes", "n": "No"}
                    )
                    if browse_choice == "y":
                        browsed = self._browse_for_analysis_files(start_dir=parent_dir.parent)
                        if browsed and browsed.get('prmtop'):
                            prmtop_files = [browsed['prmtop'][0]]
                            self.console.print(f"[green]Selected: {prmtop_files[0].name}[/green]")
                        else:
                            self.console.print("[red]No topology selected - skipping trajectory analysis[/red]")
                            nc_files = []
                    else:
                        self.console.print("[red]Skipping trajectory analysis[/red]")
                        nc_files = []

        # Determine what was selected
        has_mdout = len(mdout_files) > 0
        has_nc = len(nc_files) > 0

        if has_mdout and has_nc:
            # Both types selected - offer choice
            self.console.print("\n[bold]Selected Files:[/bold]")
            self.console.print(f"  • {len(mdout_files)} mdout file(s) (energetics)")
            self.console.print(f"  • {len(nc_files)} nc file(s) (trajectory)")

            # Display analysis options
            self.console.print("\n[bold]Select Analysis Type:[/bold]")
            self.console.print("  1. Energetics only (mdout)")
            self.console.print("  2. Trajectory only (nc)")
            self.console.print("  3. Both (recommended)")

            analysis_choice = prompt_with_context(
                self.processor,
                "Select analysis type",
                choices=["1", "2", "3"],
                default="3",
                module="MD Manager - Analysis",
                description="Choose analysis type for selected files",
                options_map={
                    "1": "Energetics only (mdout)",
                    "2": "Trajectory only (nc)",
                    "3": "Both (recommended)"
                }
            )

            if analysis_choice == "1":
                # Energetics only
                self._run_energetics_analysis(mdout_files)
            elif analysis_choice == "2":
                # Trajectory only
                self._run_trajectory_analysis(nc_files, prmtop_files[0] if prmtop_files else None)
            else:
                # Both
                self._run_energetics_analysis(mdout_files)
                self._run_trajectory_analysis(nc_files, prmtop_files[0] if prmtop_files else None)

        elif has_mdout:
            # Only mdout files - energetics analysis
            self._run_energetics_analysis(mdout_files)

        elif has_nc:
            # Only nc files - trajectory analysis
            self._run_trajectory_analysis(nc_files, prmtop_files[0] if prmtop_files else None)

        elif prmtop_files:
            # Only a topology file was selected. A topology can't be analyzed on
            # its own — it's a support file that must be paired with a trajectory
            # (.nc) for structural analysis. This commonly happens when the
            # trajectory and topology live in different directories and the user
            # picks the topology first. Mirror the trajectory-first flow: keep
            # the topology and prompt to browse for the trajectory it pairs with.
            self.console.print(f"\n[yellow]Only a topology file was selected ({prmtop_files[0].name}).[/yellow]")
            self.console.print("[yellow]A topology can't be analyzed alone — it must be paired with a trajectory (.nc).[/yellow]")
            browse_choice = prompt_with_context(
                self.processor,
                "Browse for the trajectory (.nc) to pair with this topology? (y/n)",
                choices=["y", "n"],
                default="y",
                module="MD Manager - Trajectory Analysis",
                description="Browse for trajectory to pair with selected topology",
                options_map={"y": "Yes", "n": "No"}
            )
            if browse_choice == "y":
                browsed = self._browse_for_analysis_files(start_dir=prmtop_files[0].parent)
                if browsed and browsed.get('nc'):
                    nc_files = browsed['nc']
                    # Prefer a topology found alongside the trajectory; otherwise
                    # keep the one originally selected.
                    if browsed.get('prmtop'):
                        prmtop_files = browsed['prmtop']
                    self.console.print(f"[green]Pairing {len(nc_files)} trajectory file(s) with {prmtop_files[0].name}[/green]")
                    self._run_trajectory_analysis(nc_files, prmtop_files[0])
                else:
                    self.console.print("[yellow]No trajectory selected - nothing to analyze[/yellow]")
            else:
                self.console.print("[yellow]Nothing to analyze (topology alone)[/yellow]")

        else:
            self.console.print("[yellow]No files selected for analysis[/yellow]")
            self.console.print("[grey50]Tip: select the trajectory (.nc) or energetics (.mdout/.out) file — not just the topology. "
                               "If the trajectory and topology are in different directories, you'll be prompted for the topology after picking the trajectory.[/grey50]")

    def _run_energetics_analysis(self, mdout_files: list):
        """Run energetics analysis on selected mdout files.

        Args:
            mdout_files: list of .mdout file paths
        """
        if len(mdout_files) == 1:
            self._analyze_single_mdout(mdout_files[0])
        else:
            self._analyze_concatenated_mdouts(mdout_files)

    def _run_trajectory_analysis(self, nc_files: list, prmtop: Path = None):
        """Run trajectory analysis on selected nc files.

        Args:
            nc_files: list of .nc file paths
            prmtop: Path to topology file (if already selected)
        """
        # If prmtop not provided, need to find or select topology file
        if not prmtop:
            parent_dir = nc_files[0].parent
            # Support both .prmtop and .parm7
            prmtop_files = sorted(list(parent_dir.glob("*.prmtop")) + list(parent_dir.glob("*.parm7")))

        if prmtop:
            self.console.print(f"\n[green]Using topology: {prmtop.name}[/green]")
        elif prmtop_files:
            if len(prmtop_files) == 1:
                prmtop = prmtop_files[0]
                self.console.print(f"\n[green]Using topology: {prmtop.name}[/green]")
            else:
                # Multiple topology files - let user choose
                self.console.print("\n[bold]Multiple topology files found:[/bold]")
                for i, f in enumerate(prmtop_files, 1):
                    self.console.print(f"  {i}. {f.name}")

                prmtop_choice = prompt_with_context(
                    self.processor,
                    f"Select topology (1-{len(prmtop_files)})",
                    choices=[str(i) for i in range(1, len(prmtop_files) + 1)],
                    default="1",
                    module="MD Manager - Trajectory Analysis",
                    description="Select topology file"
                )
                # Replay by basename so a changed topology list can't mis-pick.
                prmtop_choice = remap_recorded_index(self.processor, prmtop_files, str(prmtop_choice))
                prmtop = prmtop_files[int(prmtop_choice) - 1]
                annotate_selected_path(self.processor, prmtop)
        else:
            # No topology file found - browse for it
            self.console.print("\n[yellow]No topology file (.prmtop/.parm7) found in same directory[/yellow]")

            browse_for_prmtop = prompt_with_context(
                self.processor,
                "Browse for topology file? (y/n)",
                choices=["y", "n"],
                default="y",
                module="MD Manager - Trajectory Analysis",
                description="Browse for topology file",
                options_map={"y": "Yes", "n": "No"}
            )

            if browse_for_prmtop == "y":
                # Search for topology files (both .prmtop and .parm7)
                self.console.print("\n[bold]Searching for topology files...[/bold]")

                # Search for both .prmtop and .parm7 files
                prmtop_found = list(parent_dir.rglob("*.prmtop")) + list(parent_dir.rglob("*.parm7"))
                if prmtop_found:
                    self.console.print(f"[green]Found {len(prmtop_found)} topology file(s):[/green]")
                    for i, f in enumerate(prmtop_found, 1):
                        relative = f.relative_to(parent_dir)
                        self.console.print(f"  {i}. {relative}")

                    while True:
                        choice = prompt_with_context(
                            self.processor,
                            f"Select topology file (1-{len(prmtop_found)}) or 'cancel'",
                            default="cancel",
                            module="MD Manager - Trajectory Analysis",
                            description="Select topology file for trajectory analysis",
                        ).strip()
                        # Replay by basename so a changed topology list can't mis-pick.
                        choice = remap_recorded_index(self.processor, prmtop_found, choice)
                        if choice.lower() == 'cancel':
                            self.console.print("[yellow]Cannot analyze trajectory without topology file[/yellow]")
                            return
                        try:
                            idx = int(choice) - 1
                            if 0 <= idx < len(prmtop_found):
                                prmtop = prmtop_found[idx]
                                annotate_selected_path(self.processor, prmtop_found[idx])
                                break
                            else:
                                self.console.print(f"[red]Invalid selection. Choose 1-{len(prmtop_found)}[/red]")
                        except ValueError:
                            self.console.print("[red]Invalid input[/red]")
                else:
                    self.console.print("[red]No topology files (.prmtop/.parm7) found - cannot analyze trajectory[/red]")
                    return
            else:
                self.console.print("[yellow]Cannot analyze trajectory without topology file[/yellow]")
                return

        if not prmtop:
            self.console.print("[red]No topology file available - cannot analyze trajectory[/red]")
            return

        # Run trajectory analysis
        sim_name = " + ".join([f.stem for f in nc_files])
        self._analyze_trajectory(nc_files, prmtop, sim_name)

    def _analyze_single_mdout(self, mdout_file: Path):
        """Analyze a single .mdout file directly."""
        self.console.print(f"\n[bold]Analyzing Output File: {mdout_file.name}[/bold]")
        self.console.print(f"[grey50]Location: {mdout_file.parent}[/grey50]")

        if not mdout_file.exists():
            self.console.print("[red].mdout file not found[/red]")
            return

        try:
            # Create historical monitor for complete analysis
            monitor = AMBERMonitor.create_historical_monitor(str(mdout_file))

            if not any(monitor.data.values()):
                self.console.print("[yellow]No analysis data found in output file[/yellow]")
                return

            # Confirm simulation stage with user
            self._confirm_simulation_stage(monitor)

            # Run comprehensive quality assessment for the confirmed stage
            self._run_quality_assessment(monitor, mdout_file.name)

        except Exception as e:
            self.console.print(f"[red]Error analyzing output file: {e}[/red]")

    def _analyze_concatenated_mdouts(self, mdout_files: list):
        """Analyze multiple .mdout files as concatenated segments."""
        from collections import deque

        self.console.print(f"\n[bold cyan]Concatenating {len(mdout_files)} Simulation Segments[/bold cyan]")
        for i, f in enumerate(mdout_files, 1):
            self.console.print(f"  {i}. {f.name}")

        try:
            # Parse all files
            monitors = []
            for mdout_file in mdout_files:
                if not mdout_file.exists():
                    self.console.print(f"[red]File not found: {mdout_file}[/red]")
                    return

                self.console.print(f"[grey50]Parsing {mdout_file.name}...[/grey50]")
                monitor = AMBERMonitor.create_historical_monitor(str(mdout_file))

                if not any(monitor.data.values()):
                    self.console.print(f"[yellow]No data found in {mdout_file.name}[/yellow]")
                    return

                monitors.append(monitor)

            # Concatenate data from all monitors
            self.console.print(f"[grey50]Concatenating data from {len(monitors)} files...[/grey50]")
            concatenated_monitor = self._concatenate_monitors(monitors, mdout_files)

            if not any(concatenated_monitor.data.values()):
                self.console.print("[yellow]No data after concatenation[/yellow]")
                return

            # Confirm simulation stage with user
            self._confirm_simulation_stage(concatenated_monitor)

            # Run comprehensive quality assessment on concatenated data
            file_names = " + ".join([f.name for f in mdout_files])
            self._run_quality_assessment(concatenated_monitor, f"Concatenated: {file_names}")

        except Exception as e:
            self.console.print(f"[red]Error analyzing concatenated files: {e}[/red]")
            import traceback
            traceback.print_exc()

    def _concatenate_monitors(self, monitors: list, mdout_files: list):
        """Concatenate data from multiple AMBER monitors."""
        from collections import deque
        import numpy as np

        # Create a new monitor with concatenated data
        concatenated = AMBERMonitor(output_file="concatenated", max_points=100000)

        # Use system info from first monitor
        concatenated.system_info = monitors[0].system_info.copy()
        concatenated.mdin_params = monitors[0].mdin_params.copy()
        concatenated.simulation_stage = monitors[0].simulation_stage

        # Initialize data deques
        for key in monitors[0].data.keys():
            concatenated.data[key] = deque(maxlen=100000)

        # Concatenate all data
        for monitor in monitors:
            for key in monitor.data.keys():
                concatenated.data[key].extend(monitor.data[key])

        self.console.print(f"[green]✓ Concatenated {len(concatenated.data.get('time', []))} data points[/green]")

        return concatenated

    def _show_analysis_interface_for_file(self, monitor: AMBERMonitor, mdout_file: Path):
        """Show interactive analysis interface for a single file."""
        while True:
            self.console.print(f"\n[bold]Analysis Options for {mdout_file.name}:[/bold]")
            self.console.print("1. Simulation overview", highlight=False)
            self.console.print("2. Energy analysis", highlight=False)
            self.console.print("3. Temperature analysis", highlight=False)
            self.console.print("4. Pressure analysis", highlight=False)
            self.console.print("5. Performance summary", highlight=False)
            self.console.print("6. Export analysis data", highlight=False)
            self.console.print("7. ← Back", highlight=False)

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1","2","3","4","5","6","7"],
                default="1",
                module="MD Manager - Analysis",
                description="Select analysis option",
                options_map={
                    "1": "Complete simulation overview",
                    "2": "Energy analysis",
                    "3": "Temperature analysis",
                    "4": "Pressure analysis",
                    "5": "Performance summary",
                    "6": "Export analysis data",
                    "7": "← Back"
                }
            )

            if choice == "1":
                self._show_simulation_overview_for_file(monitor, mdout_file)
            elif choice == "2":
                self._show_detailed_energy_analysis(monitor)
            elif choice == "3":
                self._show_detailed_temperature_analysis(monitor)
            elif choice == "4":
                self._show_detailed_pressure_analysis(monitor)
            elif choice == "5":
                self._show_performance_summary(monitor, mdout_file)
            elif choice == "6":
                self._export_analysis_data_for_file(monitor, mdout_file)
            elif choice == "7":
                break

    def _show_simulation_overview_for_file(self, monitor: AMBERMonitor, mdout_file: Path):
        """Show overview for a single file (without trajectory info)."""
        from rich.panel import Panel
        from rich.table import Table

        self.console.print(f"\n[bold blue]Simulation Overview[/bold blue]")
        self.console.print(f"File: [cyan]{mdout_file.name}[/cyan]")
        self.console.print(f"Location: [grey50]{mdout_file.parent}[/grey50]\n")

        # Create statistics table
        table = Table(show_header=True, header_style="bold bright_blue", box=None)
        table.add_column("Property", style="bright_blue")
        table.add_column("Value", style="white")

        # Add data rows
        if monitor.data.get('step'):
            table.add_row("Total Steps", str(max(monitor.data['step'])))

        if monitor.data.get('time'):
            table.add_row("Simulation Time", f"{max(monitor.data['time']):.2f} ps")

        if monitor.data.get('total_energy'):
            energies = monitor.data['total_energy']
            avg_energy = sum(energies) / len(energies)
            table.add_row("Average Energy", f"{avg_energy:.2f} kcal/mol")

        if monitor.data.get('temperature'):
            temps = monitor.data['temperature']
            avg_temp = sum(temps) / len(temps)
            table.add_row("Average Temperature", f"{avg_temp:.2f} K")

        if monitor.data.get('pressure'):
            pressures = monitor.data['pressure']
            avg_pressure = sum(pressures) / len(pressures)
            table.add_row("Average Pressure", f"{avg_pressure:.2f} bar")

        self.console.print(table)

    def _export_analysis_data_for_file(self, monitor: AMBERMonitor, mdout_file: Path):
        """Export analysis data for a single file."""
        self.console.print(f"\n[bold blue]Export Analysis Data[/bold blue]")

        if not any(monitor.data.values()):
            self.console.print("[yellow]No data available to export[/yellow]")
            return

        # Create export directory next to the mdout file
        export_dir = mdout_file.parent / f"{mdout_file.stem}_analysis"
        export_dir.mkdir(exist_ok=True)

        try:
            # Export energy data
            if monitor.data.get('total_energy'):
                energy_file = export_dir / "energy_data.csv"
                with open(energy_file, 'w') as f:
                    f.write("Step,Time_ps,Total_Energy_kcal_mol\n")
                    steps = monitor.data.get('step', range(len(monitor.data['total_energy'])))
                    times = monitor.data.get('time', [0] * len(monitor.data['total_energy']))

                    for step, time_ps, energy in zip(steps, times, monitor.data['total_energy']):
                        f.write(f"{step},{time_ps},{energy}\n")
                self.console.print(f"[green]✓ Energy data exported to {energy_file.name}[/green]")

            # Export temperature data
            if monitor.data.get('temperature'):
                temp_file = export_dir / "temperature_data.csv"
                with open(temp_file, 'w') as f:
                    f.write("Step,Time_ps,Temperature_K\n")
                    steps = monitor.data.get('step', range(len(monitor.data['temperature'])))
                    times = monitor.data.get('time', [0] * len(monitor.data['temperature']))

                    for step, time_ps, temp in zip(steps, times, monitor.data['temperature']):
                        f.write(f"{step},{time_ps},{temp}\n")
                self.console.print(f"[green]✓ Temperature data exported to {temp_file.name}[/green]")

            # Export summary statistics
            summary_file = export_dir / "analysis_summary.txt"
            with open(summary_file, 'w') as f:
                f.write(f"AMBER Simulation Analysis Summary\n")
                f.write(f"=" * 40 + "\n\n")
                f.write(f"Output File: {mdout_file.name}\n")
                f.write(f"Location: {mdout_file.parent}\n")
                f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

                if monitor.data.get('total_energy'):
                    energies = monitor.data['total_energy']
                    f.write(f"Energy Statistics:\n")
                    f.write(f"  Average: {sum(energies)/len(energies):.2f} kcal/mol\n")
                    f.write(f"  Min: {min(energies):.2f} kcal/mol\n")
                    f.write(f"  Max: {max(energies):.2f} kcal/mol\n\n")

                if monitor.data.get('temperature'):
                    temps = monitor.data['temperature']
                    f.write(f"Temperature Statistics:\n")
                    f.write(f"  Average: {sum(temps)/len(temps):.2f} K\n")
                    f.write(f"  Min: {min(temps):.2f} K\n")
                    f.write(f"  Max: {max(temps):.2f} K\n\n")

            self.console.print(f"[green]✓ Analysis summary exported to {summary_file.name}[/green]")
            self.console.print(f"[grey50]All files saved in: {export_dir}[/grey50]")

        except Exception as e:
            self.console.print(f"[red]Error exporting data: {e}[/red]")

    def _configure_hardware(self) -> bool:
        """Configure hardware settings."""
        self._show_system_resources()
        
        self.console.print(f"\n[bold cyan]===== Hardware Configuration =====[/bold cyan]")
        self.console.print("\nSelect AMBER engine:")
        self.console.print("1. sander (single CPU)", highlight=False)
        self.console.print("2. pmemd (single CPU, optimized)", highlight=False)
        self.console.print("3. pmemd.MPI (multi-CPU)", highlight=False)
        self.console.print("4. pmemd.cuda (GPU acceleration)", highlight=False)
        
        choice = prompt_with_context(
            self.processor,
            "Select engine",
            choices=["1","2","3","4"],
            default="2",
            module="MD Manager - Hardware Configuration",
            description="Select AMBER engine",
            options_map={
                "1": "sander (single CPU)",
                "2": "pmemd (single CPU, optimized)",
                "3": "pmemd.MPI (multi-CPU)",
                "4": "pmemd.cuda (GPU acceleration)"
            }
        )

        engine_map = {
            "1": "sander",
            "2": "pmemd",
            "3": "pmemd.MPI",
            "4": "pmemd.cuda"
        }

        selected_engine = engine_map[choice]
        self.console.print(f"[green]Selected: {selected_engine}[/green]")

        # Store in workspace
        if hasattr(self, 'processor') and self.processor:
            workspace = self.processor.workspace
            workspace.set("preferred_amber_engine", selected_engine)

            if selected_engine == "pmemd.MPI":
                cpu_info = self._get_cpu_info()
                cores_str = prompt_with_context(
                    self.processor,
                    "Number of MPI tasks",
                    default=str(min(16, cpu_info['available'])),
                    module="MD Manager - Hardware Configuration",
                    description="Enter number of MPI tasks"
                )
                cores = int(cores_str)
                workspace.set("mpi_tasks", cores)
            elif selected_engine == "pmemd.cuda":
                gpu_ids = prompt_with_context(
                    self.processor,
                    "GPU IDs (comma-separated)",
                    default="0",
                    module="MD Manager - Hardware Configuration",
                    description="Enter GPU device IDs"
                )
                workspace.set("gpu_ids", gpu_ids)
        
        return True

    def _import_template(self) -> bool:
        """Import .mdin template file and integrate with AmberController."""
        self.console.print(f"\n[bold cyan]===== Import .mdin Template =====[/bold cyan]")
        
        # Find .mdin files in current directory
        mdin_files = list(Path.cwd().glob("*.mdin"))
        
        if not mdin_files:
            self.console.print("[yellow]No .mdin files found in current directory[/yellow]")
            file_path = prompt_with_context(
                self.processor,
                "Enter path to .mdin file",
                module="MD Manager - Import",
                description="Enter .mdin file path"
            )
            mdin_file = Path(file_path)
            if not mdin_file.exists():
                self.console.print("[red]File not found[/red]")
                return False
        else:
            self.console.print("Available .mdin files:")
            for i, mdin_file in enumerate(mdin_files, 1):
                self.console.print(f"  {i}. {mdin_file.name}")

            choice_str = prompt_with_context(
                self.processor,
                f"Select file (1-{len(mdin_files)})",
                default="1",
                module="MD Manager - Import",
                description="Select .mdin file to import"
            )
            choice = int(choice_str)
            if 1 <= choice <= len(mdin_files):
                mdin_file = mdin_files[choice-1]
            else:
                self.console.print("[red]Invalid selection[/red]")
                return False

        try:
            # Read and analyze the mdin file
            with open(mdin_file, 'r') as f:
                mdin_content = f.read()

            # Create template name and description
            template_name = prompt_with_context(
                self.processor,
                "Template name",
                default=mdin_file.stem,
                module="MD Manager - Import",
                description="Enter template name"
            )
            template_desc = prompt_with_context(
                self.processor,
                "Template description (optional)",
                default=f"Imported from {mdin_file.name}",
                module="MD Manager - Import",
                description="Enter template description"
            )
            
            # Determine simulation type based on mdin content
            sim_type = self._analyze_mdin_simulation_type(mdin_content)
            self.console.print(f"[grey50]Detected simulation type: {sim_type}[/grey50]")
            
            # Create custom template directory if it doesn't exist
            custom_templates_dir = Path.cwd() / "md_templates" / "user" / "custom"
            custom_templates_dir.mkdir(parents=True, exist_ok=True)
            
            # Save template file
            template_file = custom_templates_dir / f"{template_name}.mdin"
            with open(template_file, 'w') as f:
                f.write(f"! Template: {template_name}\n")
                f.write(f"! Description: {template_desc}\n") 
                f.write(f"! Imported: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"! Original file: {mdin_file.name}\n")
                f.write("!\n")
                f.write(mdin_content)
            
            # Update user metadata
            metadata_file = custom_templates_dir.parent / "metadata.json"
            metadata = {}
            if metadata_file.exists():
                try:
                    with open(metadata_file, 'r') as f:
                        metadata = json.load(f)
                except:
                    metadata = {}
            
            if "custom_templates" not in metadata:
                metadata["custom_templates"] = {}
                
            metadata["custom_templates"][template_name] = {
                "file": f"{template_name}.mdin",
                "description": template_desc,
                "simulation_type": sim_type,
                "imported_from": str(mdin_file),
                "created": time.strftime('%Y-%m-%d %H:%M:%S')
            }
            
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            self.console.print(f"[green]✓ Template '{template_name}' successfully imported[/green]")
            self.console.print(f"[grey50]Saved to: {template_file}[/grey50]")
            self.console.print(f"[grey50]Available in template system as: {sim_type}/{template_name}[/grey50]")
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error importing template: {e}[/red]")
            return False

    def _analyze_mdin_simulation_type(self, mdin_content: str) -> str:
        """Analyze mdin content to determine simulation type."""
        content_lower = mdin_content.lower()
        
        # Check for minimization
        if "imin=1" in content_lower or "imin = 1" in content_lower:
            return "minimization"
        
        # Check for heating (usually has temperature control)
        if "tempi=" in content_lower and "temp0=" in content_lower:
            tempi_match = re.search(r'tempi\s*=\s*(\d+)', content_lower)
            temp0_match = re.search(r'temp0\s*=\s*(\d+)', content_lower)
            if tempi_match and temp0_match:
                tempi = float(tempi_match.group(1))
                temp0 = float(temp0_match.group(1))
                if tempi < temp0:
                    return "heating"
        
        # Check for production (usually long nstlim)
        nstlim_match = re.search(r'nstlim\s*=\s*(\d+)', content_lower)
        if nstlim_match:
            nstlim = int(nstlim_match.group(1))
            if nstlim > 100000:  # Arbitrary threshold for production
                return "production"
        
        # Check for equilibration keywords
        if "equilibr" in content_lower or "equil" in content_lower:
            return "equilibration"
            
        # Default classification
        if "nstlim=" in content_lower:
            return "equilibration"  # Assume equilibration if not clearly production
        else:
            return "minimization"

    def _import_workflow(self) -> bool:
        """Import .json protocol file and integrate with workflow system."""
        self.console.print(f"\n[bold cyan]===== Import .json Protocol =====[/bold cyan]")
        
        # Find .json files in current directory
        json_files = list(Path.cwd().glob("*.json"))
        workflow_files = []
        
        # Filter for potential workflow files
        for json_file in json_files:
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                # Basic check for workflow structure
                if isinstance(data, dict) and any(key in data for key in ["steps", "workflow", "simulations", "phases"]):
                    workflow_files.append(json_file)
            except:
                continue
        
        if not workflow_files:
            self.console.print("[yellow]No protocol .json files found in current directory[/yellow]")
            file_path = prompt_with_context(
                self.processor,
                "Enter path to .json protocol file",
                module="MD Manager - Import",
                description="Enter .json protocol file path"
            )
            json_file = Path(file_path)
            if not json_file.exists():
                self.console.print("[red]File not found[/red]")
                return False
        else:
            self.console.print("Available protocol .json files:")
            for i, json_file in enumerate(workflow_files, 1):
                self.console.print(f"  {i}. {json_file.name}")

            choice_str = prompt_with_context(
                self.processor,
                f"Select file (1-{len(workflow_files)})",
                default="1",
                module="MD Manager - Import",
                description="Select .json protocol file to import"
            )
            choice = int(choice_str)
            if 1 <= choice <= len(workflow_files):
                json_file = workflow_files[choice-1]
            else:
                self.console.print("[red]Invalid selection[/red]")
                return False

        try:
            # Load and validate the workflow
            with open(json_file, 'r') as f:
                workflow_data = json.load(f)

            # Analyze workflow structure
            workflow_info = self._analyze_workflow_structure(workflow_data)
            self.console.print(f"[grey50]Detected protocol structure: {workflow_info['type']} with {workflow_info['step_count']} steps[/grey50]")

            # Get workflow details
            workflow_name = prompt_with_context(
                self.processor,
                "Protocol name",
                default=json_file.stem,
                module="MD Manager - Import",
                description="Enter protocol name"
            )
            workflow_desc = prompt_with_context(
                self.processor,
                "Protocol description (optional)",
                default=workflow_data.get('description', f"Imported from {json_file.name}"),
                module="MD Manager - Import",
                description="Enter protocol description"
            )
            
            # Create workflows directory structure
            workflows_dir = Path.cwd() / "md_workflows" / "user" / "custom"
            workflows_dir.mkdir(parents=True, exist_ok=True)
            
            # Save workflow file
            workflow_file = workflows_dir / f"{workflow_name}_workflow.json"
            
            # Enhance workflow data with metadata
            enhanced_workflow = {
                "name": workflow_name,
                "description": workflow_desc,
                "imported_from": str(json_file),
                "created": time.strftime('%Y-%m-%d %H:%M:%S'),
                "version": "1.0",
                "type": workflow_info['type']
            }
            
            # Merge with original data
            enhanced_workflow.update(workflow_data)
            
            with open(workflow_file, 'w') as f:
                json.dump(enhanced_workflow, f, indent=2)
            
            # Update workflows metadata
            metadata_file = workflows_dir.parent / "metadata.json"
            metadata = {}
            if metadata_file.exists():
                try:
                    with open(metadata_file, 'r') as f:
                        metadata = json.load(f)
                except:
                    metadata = {}
            
            if "custom_workflows" not in metadata:
                metadata["custom_workflows"] = {}
                
            metadata["custom_workflows"][workflow_name] = {
                "file": f"{workflow_name}_workflow.json",
                "description": workflow_desc,
                "type": workflow_info['type'],
                "step_count": workflow_info['step_count'],
                "imported_from": str(json_file),
                "created": time.strftime('%Y-%m-%d %H:%M:%S')
            }
            
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            self.console.print(f"[green]✓ Protocol '{workflow_name}' successfully imported[/green]")
            self.console.print(f"[grey50]Saved to: {workflow_file}[/grey50]")
            self.console.print(f"[grey50]Type: {workflow_info['type']} protocol with {workflow_info['step_count']} steps[/grey50]")
            
            return True
            
        except json.JSONDecodeError:
            self.console.print("[red]Invalid JSON format[/red]")
            return False
        except Exception as e:
            self.console.print(f"[red]Error importing protocol: {e}[/red]")
            return False

    def _analyze_workflow_structure(self, workflow_data: Dict[str, Any]) -> Dict[str, Any]:
        """Analyze workflow JSON to determine structure and properties."""
        info = {
            "type": "unknown",
            "step_count": 0,
            "has_dependencies": False
        }
        
        # Check for different workflow structures
        if "steps" in workflow_data:
            steps = workflow_data["steps"]
            if isinstance(steps, list):
                info["step_count"] = len(steps)
                info["type"] = "sequential"
                # Check for dependencies
                for step in steps:
                    if isinstance(step, dict) and "depends_on" in step:
                        info["has_dependencies"] = True
                        info["type"] = "dependency-based"
                        break
            elif isinstance(steps, dict):
                info["step_count"] = len(steps)
                info["type"] = "named_steps"
                
        elif "workflow" in workflow_data:
            workflow = workflow_data["workflow"]
            if isinstance(workflow, list):
                info["step_count"] = len(workflow)
                info["type"] = "workflow_list"
                
        elif "simulations" in workflow_data:
            simulations = workflow_data["simulations"]
            if isinstance(simulations, list):
                info["step_count"] = len(simulations)
                info["type"] = "simulation_batch"
                
        elif "phases" in workflow_data:
            phases = workflow_data["phases"]
            if isinstance(phases, (list, dict)):
                info["step_count"] = len(phases)
                info["type"] = "phase_based"
                
        return info
        
    # ====================
    # 4-STEP WORKFLOW HELPER METHODS
    # ====================
    
    def _remove_from_selection(self, selected_templates):
        """Remove templates from the current selection."""
        if not selected_templates:
            self.console.print("[yellow]No templates selected to remove[/yellow]")
            return
            
        self.console.print("\n[bold]Remove Templates from Selection:[/bold]")
        
        # Display current templates with numbers
        templates = self.user_data_manager.list_templates()
        for i, template_id in enumerate(selected_templates, 1):
            template_name = templates.get(template_id, {}).get('name', template_id)
            self.console.print(f"  {i:2}. {template_name}")
            
        self.console.print(f"\n[grey50]0. Cancel removal[/grey50]")
        self.console.print("[grey50]Tip: Enter comma-separated numbers for multiple removals (e.g., '1,3')[/grey50]")
        
        # Get user selection
        valid_choices = [str(i) for i in range(1, len(selected_templates)+1)] + ["0"]
        user_input = prompt_with_context(
            self.processor,
            "Select templates to remove",
            default="0",
            module="MD Manager - Template Library",
            description="Select templates to remove (comma-separated, 0 to cancel)"
        )
        
        if user_input == "0":
            return
            
        # Parse comma-separated selections
        try:
            choices = [choice.strip() for choice in user_input.split(',')]
            indices_to_remove = []
            
            for choice in choices:
                if choice not in valid_choices:
                    self.console.print(f"[red]Invalid choice: {choice}[/red]")
                    return
                if choice != "0":
                    indices_to_remove.append(int(choice) - 1)  # Convert to 0-based
                    
            if not indices_to_remove:
                return
                
            # Sort in reverse order to avoid index shifting issues
            indices_to_remove.sort(reverse=True)
            
            # Remove templates and show what was removed
            removed_names = []
            for index in indices_to_remove:
                if 0 <= index < len(selected_templates):
                    template_id = selected_templates[index]
                    template_name = templates.get(template_id, {}).get('name', template_id)
                    removed_names.append(template_name)
                    selected_templates.pop(index)
                    
            if removed_names:
                self.console.print(f"[green]✓ Removed {len(removed_names)} template{'s' if len(removed_names) != 1 else ''}:[/green]")
                for name in removed_names:
                    self.console.print(f"  - {name}")
                    
        except (ValueError, IndexError) as e:
            self.console.print(f"[red]Error parsing selection: {e}[/red]")
            
    def _manage_template_library_extended(self, controller):
        """Extended template library management including create/modify and preview."""
        self.console.print("\n[bold cyan]===== Manage Template Library =====[/bold cyan]")
        
        while True:
            self.console.print("\n[bold]Template Library Management:[/bold]")
            self.console.print("  1. View library statistics")
            self.console.print("  2. Preview template content")
            self.console.print("  3. Import template from file")
            self.console.print("  4. Create or modify template")
            self.console.print("  5. Apply redox site restraint to template")
            self.console.print("  6. Export template to file") 
            self.console.print("  7. Delete user template")
            self.console.print("  8. Back to template selection")
            
            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3", "4", "5", "6", "7", "8"],
                default="8",
                module="MD Manager - Template Library",
                description="Template library management",
                options_map={
                    "1": "View library statistics",
                    "2": "Preview template content",
                    "3": "Import template from file",
                    "4": "Create or modify template",
                    "5": "Apply redox site restraint to template",
                    "6": "Export template to file",
                    "7": "Delete user template",
                    "8": "Back to template selection"
                }
            )
            
            if choice == "1":
                # View statistics
                self._view_library_statistics()
            elif choice == "2":
                # Preview template content
                self._preview_library_templates()
            elif choice == "3":
                # Import template
                self._import_template_to_library()
            elif choice == "4":
                # Create or modify template
                template_id = self._create_or_modify_template_with_wizard(controller)
                if template_id:
                    templates = self.user_data_manager.list_templates()
                    template_name = templates.get(template_id, {}).get('name', template_id)
                    self.console.print(f"[green]✓ Template '{template_name}' ready for use[/green]")
            elif choice == "5":
                # Apply redox site restraint to template
                self._apply_redox_restraint_to_template()
            elif choice == "6":
                # Export template
                self._export_template_from_library()
            elif choice == "7":
                # Delete template
                self._delete_template_from_library()
            elif choice == "8":
                # Back to template selection
                break
                
    def _apply_redox_restraint_to_template(self):
        """Apply redox site restraint mask to a selected template."""
        # Check if redox restraint mask exists in workspace
        if not self.workspace:
            self.console.print("[yellow]Workspace not available[/yellow]")
            return
            
        redox_mask = self.workspace.get("redox_restraint_mask")
        if not redox_mask:
            self.console.print("[yellow]No redox restraint mask found in workspace[/yellow]")
            self.console.print("[grey50]Generate a redox restraint mask using the Metallo Preparation module first[/grey50]")
            return
            
        from rich.markup import escape
        self.console.print(f"\n[bold]Redox restraint mask found:[/bold] [cyan]{escape(redox_mask)}[/cyan]")
        
        # Select template to apply to
        self.console.print("\n[bold]Select template to apply restraint to:[/bold]")
        templates = self.user_data_manager.list_templates()
        
        if not templates:
            self.console.print("[yellow]No templates available[/yellow]")
            return
            
        # Display templates with numbering
        template_list = []
        choice_num = 1
        for template_id, template_info in templates.items():
            source_tag = f"[{template_info['source']}]"
            if template_info['source'] == 'builtin':
                source_tag = "[cyan][builtin][/cyan]"
            elif template_info['source'] == 'custom':
                source_tag = "[green][custom][/green]"
            elif 'modified' in template_info['source']:
                source_tag = "[yellow][modified][/yellow]"
                
            self.console.print(f"  {choice_num}. {template_info['name']} {source_tag}")
            template_list.append((template_id, template_info))
            choice_num += 1
            
        self.console.print(f"  {choice_num}. Cancel")
        
        choices = [str(i) for i in range(1, choice_num + 1)]

        # Build options map
        options_map = {}
        for i, (template_id, template_info) in enumerate(template_list, 1):
            options_map[str(i)] = template_info['name']
        options_map[str(choice_num)] = "Cancel"

        choice = prompt_with_context(
            self.processor,
            "Select template",
            choices=choices,
            default=str(choice_num),
            module="MD Manager - Template Library",
            description="Select template for restraint application",
            options_map=options_map
        )
        
        if choice == str(choice_num):  # Cancel
            return
            
        selected_template_id, selected_template_info = template_list[int(choice) - 1]
        
        # Get template content
        template_content, template_metadata = self.user_data_manager.get_template_content(selected_template_id)
        if not template_content:
            self.console.print(f"[red]Could not load template content[/red]")
            return
            
        # Show before/after preview
        self.console.print(f"\n[bold]Template:[/bold] {selected_template_info['name']}")
        
        # Find current restraintmask line
        current_mask_line = None
        for line in template_content.split('\n'):
            if 'restraintmask' in line.lower() and '=' in line:
                current_mask_line = line.strip()
                break
                
        if current_mask_line:
            # Escape the mask line to prevent Rich from interpreting it as markup
            from rich.markup import escape
            self.console.print(f"[bold]Current restraint:[/bold] [grey50]{escape(current_mask_line)}[/grey50]")
        else:
            self.console.print("[yellow]No restraintmask parameter found in template[/yellow]")
            add_mask_str = prompt_with_context(
                self.processor,
                "Add restraintmask parameter to template?",
                choices=["y", "n"],
                default="y",
                module="MD Manager - Template Library",
                description="Add missing restraintmask parameter",
                options_map={"y": "Yes, add parameter", "n": "No, cancel"}
            )
            if not (add_mask_str.lower() == "y"):
                return
                
        # Apply the redox restraint mask
        modified_content = self._apply_redox_restraint_mask(template_content)
        
        # Show what will be changed
        new_mask_line = None
        for line in modified_content.split('\n'):
            if 'restraintmask' in line.lower() and '=' in line:
                new_mask_line = line.strip()
                break
                
        if new_mask_line:
            from rich.markup import escape
            self.console.print(f"[bold]New restraint:[/bold] [green]{escape(new_mask_line)}[/green]")
            
        # Save options based on template source
        self.console.print("\n[bold]Save options:[/bold]")
        
        if selected_template_info['source'] == 'builtin':
            # Cannot overwrite builtin templates
            self.console.print("[yellow]Built-in templates cannot be modified directly[/yellow]")
            self.console.print("  1. Save as new custom template")
            self.console.print("  2. Cancel")
            
            save_choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2"],
                default="2",
                module="MD Manager - Template Library",
                description="Save builtin template with restraints",
                options_map={"1": "Save as new custom template", "2": "Cancel"}
            )

            if save_choice == "1":
                # Save as new custom template
                new_name = prompt_with_context(
                    self.processor,
                    "Enter name for new template",
                    default=f"{selected_template_info['name']}_redox",
                    module="MD Manager - Template Library",
                    description="Enter new template name"
                )

                new_description = prompt_with_context(
                    self.processor,
                    "Enter description",
                    default=f"{selected_template_info.get('description', '')} (with redox restraints)",
                    module="MD Manager - Template Library",
                    description="Enter template description"
                )
                
                # Save the new template
                new_template_id = self.user_data_manager.save_custom_template(
                    content=modified_content,
                    name=new_name,
                    description=new_description,
                    template_type=selected_template_info.get('type', 'unknown'),
                    based_on=selected_template_id
                )
                
                if new_template_id:
                    self.console.print(f"[green]✓ Saved as new template: {new_name}[/green]")
                    self.console.print(f"[grey50]Template ID: {new_template_id}[/grey50]")
                else:
                    self.console.print(f"[red]Failed to save template[/red]")
                    
        else:
            # Custom or modified template - can overwrite or save as new
            self.console.print("  1. Overwrite existing template")
            self.console.print("  2. Save as new template")
            self.console.print("  3. Cancel")
            
            save_choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3"],
                default="3",
                module="MD Manager - Template Library",
                description="Save custom template with restraints",
                options_map={
                    "1": "Overwrite existing template",
                    "2": "Save as new template",
                    "3": "Cancel"
                }
            )

            if save_choice == "1":
                # Overwrite existing
                if self.user_data_manager.update_custom_template(
                    template_id=selected_template_id,
                    content=modified_content
                ):
                    self.console.print(f"[green]✓ Updated template: {selected_template_info['name']}[/green]")
                else:
                    self.console.print(f"[red]Failed to update template[/red]")

            elif save_choice == "2":
                # Save as new
                new_name = prompt_with_context(
                    self.processor,
                    "Enter name for new template",
                    default=f"{selected_template_info['name']}_redox",
                    module="MD Manager - Template Library",
                    description="Enter new template name"
                )

                new_description = prompt_with_context(
                    self.processor,
                    "Enter description",
                    default=f"{selected_template_info.get('description', '')} (with redox restraints)",
                    module="MD Manager - Template Library",
                    description="Enter template description"
                )
                
                new_template_id = self.user_data_manager.save_custom_template(
                    content=modified_content,
                    name=new_name,
                    description=new_description,
                    template_type=selected_template_info.get('type', 'unknown'),
                    based_on=selected_template_id
                )
                
                if new_template_id:
                    self.console.print(f"[green]✓ Saved as new template: {new_name}[/green]")
                    self.console.print(f"[grey50]Template ID: {new_template_id}[/grey50]")
                else:
                    self.console.print(f"[red]Failed to save template[/red]")
    
    def _manage_template_library(self):
        """Original template library management - import, export, delete templates."""
        self.console.print("\n[bold cyan]===== Manage Template Library =====[/bold cyan]")
        
        while True:
            self.console.print("\n[bold]Library Management Options:[/bold]")
            self.console.print("  1. Import template from file")
            self.console.print("  2. Export template to file")
            self.console.print("  3. Delete user template")
            self.console.print("  4. View library statistics")
            self.console.print("  5. Back to template selection")
            
            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3", "4", "5"],
                default="5",
                module="MD Manager - Template Library",
                description="Library management options",
                options_map={
                    "1": "Import template from file",
                    "2": "Export template to file",
                    "3": "Delete user template",
                    "4": "View library statistics",
                    "5": "Back to template selection"
                }
            )
            
            if choice == "1":
                # Import template
                self._import_template_to_library()
            elif choice == "2":
                # Export template
                self._export_template_from_library()
            elif choice == "3":
                # Delete template
                self._delete_template_from_library()
            elif choice == "4":
                # View statistics
                self._view_library_statistics()
            elif choice == "5":
                # Back to template selection
                break
                
    def _import_template_to_library(self):
        """Import a template file into the user template library."""
        self.console.print("\n[bold]Import Template to Library[/bold]")
        
        # Enhanced import with directory navigation
        self.console.print("\n[bold]Import Options:[/bold]")
        self.console.print("  1. Enter file path directly")
        self.console.print("  2. Browse and search for .mdin files")
        self.console.print("  3. Back to library management")
        
        import_choice = prompt_with_context(
            self.processor,
            "Select import method",
            choices=["1", "2", "3"],
            default="1",
            module="MD Manager - Template Import",
            description="Select import method",
            options_map={
                "1": "Enter file path directly",
                "2": "Browse and search for .mdin files",
                "3": "Back to library management"
            }
        )

        if import_choice == "3":
            return

        file_path = None
        if import_choice == "1":
            # Direct path entry
            path_input = prompt_with_context(
                self.processor,
                "Enter path to .mdin template file",
                module="MD Manager - Template Import",
                description="Enter file path"
            )
            file_path = Path(path_input).expanduser()
            
            if not file_path.exists():
                self.console.print(f"[red]File not found: {file_path}[/red]")
                return
                
            if file_path.suffix != '.mdin':
                self.console.print("[red]File must have .mdin extension[/red]")
                return
        else:
            # Enhanced file browser
            file_path = self._browse_for_template_import()
            if not file_path:
                self.console.print("[yellow]Import cancelled[/yellow]")
                return
            
        try:
            # Read the template file
            with open(file_path, 'r') as f:
                content = f.read()
                
            # Parse header metadata if present
            lines = content.split('\n')
            metadata = {}
            for line in lines[:10]:  # Check first 10 lines for metadata
                if line.startswith('! TEMPLATE:'):
                    metadata['name'] = line.replace('! TEMPLATE:', '').strip()
                elif line.startswith('! DESCRIPTION:'):
                    metadata['description'] = line.replace('! DESCRIPTION:', '').strip()
                elif line.startswith('! TYPE:'):
                    metadata['simulation_type'] = line.replace('! TYPE:', '').strip()
                elif line.startswith('! PRIORITY:'):
                    try:
                        metadata['priority'] = int(line.replace('! PRIORITY:', '').strip())
                    except:
                        pass
                elif line.startswith('! AUTHOR:'):
                    metadata['author'] = line.replace('! AUTHOR:', '').strip()
                elif line.startswith('! VERSION:'):
                    metadata['version'] = line.replace('! VERSION:', '').strip()
                elif line.startswith('! SOURCE:'):
                    metadata['source'] = line.replace('! SOURCE:', '').strip()
                        
            # Prompt for missing metadata
            self.console.print(f"\n[bold]Template Metadata[/bold]")
            if 'name' not in metadata:
                metadata['name'] = prompt_with_context(
                    self.processor,
                    "Template name",
                    default=file_path.stem,
                    module="MD Manager - Template Import",
                    description="Enter template name"
                )
            if 'description' not in metadata:
                metadata['description'] = prompt_with_context(
                    self.processor,
                    "Template description",
                    module="MD Manager - Template Import",
                    description="Enter template description"
                )
            if 'simulation_type' not in metadata:
                sim_types = ["minimization", "heating", "equilibration", "production", "other"]
                self.console.print("\nSimulation types:")
                for i, st in enumerate(sim_types, 1):
                    self.console.print(f"  {i}. {st}")
                type_choice = prompt_with_context(
                    self.processor,
                    "Select type",
                    choices=[str(i) for i in range(1, 6)],
                    default="1",
                    module="MD Manager - Template Import",
                    description="Select simulation type",
                    options_map={str(i): st for i, st in enumerate(sim_types, 1)}
                )
                metadata['simulation_type'] = sim_types[int(type_choice)-1]
            if 'priority' not in metadata:
                try:
                    priority_input = prompt_with_context(
                        self.processor,
                        "Priority (lower = higher priority)",
                        default="10",
                        module="MD Manager - Template Import",
                        description="Enter template priority"
                    )
                    metadata['priority'] = int(priority_input)
                except:
                    metadata['priority'] = 10
            if 'author' not in metadata:
                metadata['author'] = prompt_with_context(
                    self.processor,
                    "Author name",
                    default="User",
                    module="MD Manager - Template Import",
                    description="Enter author name"
                )
            if 'version' not in metadata:
                metadata['version'] = prompt_with_context(
                    self.processor,
                    "Version",
                    default="1.0",
                    module="MD Manager - Template Import",
                    description="Enter version"
                )
            if 'source' not in metadata:
                metadata['source'] = prompt_with_context(
                    self.processor,
                    "Source (URL or reference)",
                    default="imported",
                    module="MD Manager - Template Import",
                    description="Enter source reference"
                )
                
            # Add metadata header to content if not present
            if not any(line.startswith('! TEMPLATE:') for line in content.split('\n')[:10]):
                header = f"""! TEMPLATE: {metadata['name']}
! DESCRIPTION: {metadata['description']}
! TYPE: {metadata['simulation_type']}
! PRIORITY: {metadata['priority']}
! AUTHOR: {metadata['author']}
! VERSION: {metadata['version']}
! SOURCE: {metadata['source']}

"""
                content = header + content
                
            # Create template using UserDataManager
            template_id = self.user_data_manager.create_template(
                metadata['name'],
                metadata['description'],
                metadata['simulation_type'],
                content
            )
            
            if template_id:
                self.console.print(f"[green]✓ Template '{metadata['name']}' imported successfully[/green]")
                self.console.print(f"[grey50]Template ID: {template_id}[/grey50]")
            else:
                self.console.print("[red]Failed to import template[/red]")
                
        except Exception as e:
            self.console.print(f"[red]Error importing template: {e}[/red]")
            
    def _export_template_from_library(self):
        """Export a template from the library to a file."""
        self.console.print("\n[bold]Export Template from Library[/bold]")
        
        # Get available templates
        templates = self.user_data_manager.list_templates()
        if not templates:
            self.console.print("[yellow]No templates available to export[/yellow]")
            return
            
        # Sort templates by priority (lower number = higher priority)
        sorted_templates = sorted(templates.items(), key=lambda x: x[1].get('priority', 999))
        
        # Display templates for selection
        self.console.print("\n[bold]Select template to export:[/bold]")
        template_choices = {}
        choice_num = 1
        
        for template_id, metadata in sorted_templates:
            template_choices[str(choice_num)] = template_id
            name = metadata.get('name', template_id)
            source = metadata.get('source', 'unknown')
            priority = metadata.get('priority', 'N/A')
            self.console.print(f"  {choice_num}. {name} [{source}] (priority: {priority})")
            choice_num += 1
            
        self.console.print(f"  0. Cancel")
        self.console.print("[grey50]Tip: Enter comma-separated numbers for multiple exports (e.g., '1,3,5')[/grey50]")

        valid_choices = list(template_choices.keys()) + ["0"]
        user_input = prompt_with_context(
            self.processor,
            "Select template(s)",
            default="0",
            module="MD Manager - Template Export",
            description="Select templates to export (comma-separated, 0 to cancel)"
        )
        
        if user_input == "0":
            return
            
        # Parse comma-separated selections
        try:
            choices = [choice.strip() for choice in user_input.split(',')]
            selected_template_ids = []
            
            for choice in choices:
                if choice not in valid_choices:
                    self.console.print(f"[red]Invalid choice: {choice}[/red]")
                    return
                if choice != "0":
                    template_id = template_choices[choice]
                    selected_template_ids.append(template_id)
                    
            if not selected_template_ids:
                return
                
            # Export each selected template
            export_count = 0
            for template_id in selected_template_ids:
                template_metadata = templates[template_id]
                
                if len(selected_template_ids) == 1:
                    # Single template - allow custom filename
                    default_name = f"{template_metadata.get('name', template_id)}.mdin"
                    export_path_str = prompt_with_context(
                        self.processor,
                        "Export filename",
                        default=default_name,
                        module="MD Manager - Template Export",
                        description="Enter export filename"
                    )
                    export_path = Path(export_path_str)
                else:
                    # Multiple templates - use default naming
                    filename = f"{template_metadata.get('name', template_id)}.mdin"
                    export_path = Path(filename)
                    
                # Export using UserDataManager
                if self.user_data_manager.export_content(template_id, export_path):
                    self.console.print(f"[green]✓ Template '{template_metadata.get('name', template_id)}' exported to {export_path}[/green]")
                    export_count += 1
                else:
                    self.console.print(f"[red]Failed to export template '{template_metadata.get('name', template_id)}'[/red]")
                    
            if export_count > 1:
                self.console.print(f"[green]✓ Successfully exported {export_count} templates[/green]")
                
        except Exception as e:
            self.console.print(f"[red]Error parsing selection: {e}[/red]")
            
    def _parse_numeric_selection(self, user_input: str, valid_choices: List[str]) -> Optional[List[str]]:
        """Parse user selection supporting comma-separated values and ranges (e.g., '1,3-5,7').

        Returns an ordered list of selected choices (deduplicated, excluding '0'),
        or None if the input contains an invalid token (error already printed).
        """
        tokens = [tok.strip() for tok in user_input.split(',') if tok.strip()]
        selected: List[str] = []
        seen = set()

        for tok in tokens:
            if '-' in tok and tok.count('-') == 1 and not tok.startswith('-') and not tok.endswith('-'):
                start_s, end_s = tok.split('-', 1)
                start_s, end_s = start_s.strip(), end_s.strip()
                if not (start_s.isdigit() and end_s.isdigit()):
                    self.console.print(f"[red]Invalid range: {tok}[/red]")
                    return None
                start, end = int(start_s), int(end_s)
                if start > end:
                    start, end = end, start
                for n in range(start, end + 1):
                    choice = str(n)
                    if choice not in valid_choices:
                        self.console.print(f"[red]Invalid choice in range {tok}: {choice}[/red]")
                        return None
                    if choice != "0" and choice not in seen:
                        seen.add(choice)
                        selected.append(choice)
            else:
                if tok not in valid_choices:
                    self.console.print(f"[red]Invalid choice: {tok}[/red]")
                    return None
                if tok != "0" and tok not in seen:
                    seen.add(tok)
                    selected.append(tok)

        return selected

    def _delete_template_from_library(self):
        """Delete a user template from the library."""
        self.console.print("\n[bold]Delete Template from Library[/bold]")
        
        # Get only user templates (not builtin)
        all_templates = self.user_data_manager.list_templates()
        user_templates = {k: v for k, v in all_templates.items() if v.get('source') != 'builtin'}
        
        if not user_templates:
            self.console.print("[yellow]No user templates available to delete[/yellow]")
            return
            
        # Display user templates for selection
        self.console.print("\n[bold]Select template to delete:[/bold]")
        template_choices = {}
        choice_num = 1
        
        for template_id, metadata in user_templates.items():
            template_choices[str(choice_num)] = template_id
            name = metadata.get('name', template_id)
            self.console.print(f"  {choice_num}. {name}")
            choice_num += 1
            
        self.console.print(f"  0. Cancel")
        self.console.print("[grey50]Tip: Enter comma-separated numbers and/or ranges (e.g., '1,3-5,7') or 'all' to delete all[/grey50]")

        valid_choices = list(template_choices.keys()) + ["0", "all"]
        user_input = prompt_with_context(
            self.processor,
            "Select template(s) to delete",
            default="0",
            module="MD Manager - Template Library",
            description="Select templates to delete (comma-separated, ranges, 'all', or 0 to cancel)"
        )

        if user_input == "0":
            return

        # Parse comma-separated selections, ranges, or 'all'
        try:
            # Handle 'all' option
            if user_input.lower() == "all":
                selected_template_ids = list(template_choices.values())
            else:
                selected_choices = self._parse_numeric_selection(user_input, valid_choices)
                if selected_choices is None:
                    return
                selected_template_ids = [template_choices[c] for c in selected_choices]
                    
            if not selected_template_ids:
                return
                
            # Show what will be deleted and confirm
            self.console.print(f"\n[bold]Templates to delete:[/bold]")
            for template_id in selected_template_ids:
                template_metadata = user_templates[template_id]
                self.console.print(f"  - {template_metadata.get('name', template_id)}")
                
            # Confirm deletion
            confirm_str = prompt_with_context(
                self.processor,
                f"[red]Delete {len(selected_template_ids)} template{'s' if len(selected_template_ids) != 1 else ''}?[/red]",
                choices=["y", "n"],
                default="n",
                module="MD Manager - Template Library",
                description="Confirm template deletion",
                options_map={"y": "Yes, delete", "n": "No, cancel"}
            )

            if confirm_str.lower() == "y":
                delete_count = 0
                for template_id in selected_template_ids:
                    template_meta = user_templates[template_id]
                    try:
                        deleted_file = False

                        # Delete template file if path exists in metadata
                        template_path_str = template_meta.get('template_path')
                        if template_path_str:
                            template_path = self.user_data_manager.template_base_dir / template_path_str
                            if template_path.exists():
                                template_path.unlink()
                                deleted_file = True

                        # Also check for custom JSON template file (not tracked in metadata)
                        custom_dir = self.user_data_manager.user_template_dir / "custom"
                        if custom_dir.exists():
                            for json_file in custom_dir.glob(f"{template_id}_*.json"):
                                json_file.unlink()
                                deleted_file = True

                        # Remove from metadata if present
                        if template_id in self.user_data_manager.template_metadata:
                            del self.user_data_manager.template_metadata[template_id]

                        if deleted_file or template_id in user_templates:
                            self.console.print(f"[green]✓ Deleted template '{template_meta.get('name', template_id)}'[/green]")
                            delete_count += 1
                    except Exception as e:
                        self.console.print(f"[red]Error deleting template '{template_meta.get('name', template_id)}': {e}[/red]")
                        
                # Save metadata changes
                if delete_count > 0:
                    try:
                        self.user_data_manager._save_metadata(
                            self.user_data_manager.template_metadata,
                            self.user_data_manager.template_metadata_file
                        )
                        self.console.print(f"[green]✓ Successfully deleted {delete_count} template{'s' if delete_count != 1 else ''}[/green]")
                    except Exception as e:
                        self.console.print(f"[red]Error saving metadata: {e}[/red]")
            else:
                self.console.print("[grey50]Deletion cancelled[/grey50]")
                
        except Exception as e:
            self.console.print(f"[red]Error parsing selection: {e}[/red]")

    def _manage_clusters(self):
        """Manage cluster profiles: list, show, create, export, import, delete."""
        from proprep.md_prep import cluster_profile as cp

        def _strip_notes(obj):
            if isinstance(obj, dict):
                return {k: _strip_notes(v) for k, v in obj.items() if not k.startswith("_notes")}
            if isinstance(obj, list):
                return [_strip_notes(i) for i in obj]
            return obj
        # Expose for nested helper that's called via self._clusters_create_from_template.
        self._strip_notes = _strip_notes

        while True:
            self.console.print("\n[bold cyan]===== Cluster Profile Management =====[/bold cyan]")
            entries = cp.list_profiles()
            visible = [e for e in entries if not e['name'].startswith('_')]
            self.console.print("\n[bold]Available profiles:[/bold]")
            if visible:
                for e in visible:
                    marker = f"[grey50]({e['source']})[/grey50]"
                    self.console.print(f"  • {e['display_name']} [grey50]{e['name']}[/grey50] {marker}")
            else:
                self.console.print("  [grey50](none — use option 3 to create from template)[/grey50]")

            self.console.print("\n[bold]Actions:[/bold]")
            self.console.print("  1. Show profile details")
            self.console.print("  2. Edit/fill required fields (quick)")
            self.console.print("  3. Edit full profile JSON in $EDITOR")
            self.console.print("  4. Create new profile from template")
            self.console.print("  5. Export profile to file")
            self.console.print("  6. Import profile from file")
            self.console.print("  7. Delete user profile")
            self.console.print("\n[bold]Navigation:[/bold]")
            self.console.print("  b. ← Back")

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3", "4", "5", "6", "7", "b"],
                default="b",
                module="MD Manager - Cluster Profiles",
                description="Cluster profiles action",
                options_map={
                    "1": "Show profile details",
                    "2": "Edit/fill required fields (quick)",
                    "3": "Edit full profile JSON in $EDITOR",
                    "4": "Create new profile from template",
                    "5": "Export profile to file",
                    "6": "Import profile from file",
                    "7": "Delete user profile",
                    "b": "Back",
                },
            )
            if choice == "b":
                return
            if choice == "1":
                self._clusters_show_details(visible)
            elif choice == "2":
                self._clusters_edit_profile(visible)
            elif choice == "3":
                self._clusters_edit_in_editor(visible)
            elif choice == "4":
                self._clusters_create_from_template()
            elif choice == "5":
                self._clusters_export(visible)
            elif choice == "6":
                self._clusters_import()
            elif choice == "7":
                self._clusters_delete(visible)

    def _clusters_show_details(self, entries):
        from proprep.md_prep import cluster_profile as cp
        if not entries:
            return
        self.console.print("\nAvailable profiles:")
        for i, e in enumerate(entries, 1):
            self.console.print(f"  {i}. {e['name']}")
        raw = prompt_with_context(
            self.processor,
            "Pick one (or 'c' to cancel)",
            default="1",
            module="MD Manager - Cluster Profiles",
            description="Pick profile from list",
        )
        if raw.strip().lower() == 'c':
            return
        try:
            idx = int(raw) - 1
            name = entries[idx]['name']
        except (ValueError, IndexError):
            self.console.print("[red]Invalid selection[/red]")
            return
        profile = cp.load_profile(name)
        self.console.print(f"\n[bold]{profile.display_name}[/bold] [grey50]({profile.name})[/grey50]")
        self.console.print(f"  {profile.description}")
        self.console.print(f"  Account: {profile.cluster.get('account', '') or '[not set]'}")
        self.console.print(f"  Modules: {', '.join(profile.cluster.get('modules', []))}")
        self.console.print(f"  Resource classes: {', '.join(profile.resource_classes)}")
        self.console.print(f"  Defaults: cpu={profile.defaults.get('cpu_class')}, gpu={profile.defaults.get('gpu_class')}")
        errors = cp.validate(profile)
        if errors:
            self.console.print("[yellow]Validation warnings:[/yellow]")
            for e in errors:
                self.console.print(f"  • {e}")

    def _clusters_edit_profile(self, entries):
        """Prompt for required fields of the chosen profile and save to user scope."""
        from proprep.md_prep import cluster_profile as cp
        if not entries:
            return
        self.console.print("\nProfiles with missing fields:")
        candidates = []
        for e in entries:
            p = cp.load_profile(e['name'])
            if p.missing_required_fields() or p.cluster.get('account') == "":
                candidates.append(p)
        if not candidates:
            self.console.print("[green]All profiles have their required fields filled in[/green]")
            return
        for i, p in enumerate(candidates, 1):
            missing = ", ".join(p.missing_required_fields()) or "(account empty)"
            self.console.print(f"  {i}. {p.name} [grey50]— {missing}[/grey50]")
        raw = prompt_with_context(
            self.processor,
            "Pick one (or 'c' to cancel)",
            default="1",
            module="MD Manager - Cluster Profiles",
            description="Pick profile from list",
        )
        if raw.strip().lower() == 'c':
            return
        raw = remap_recorded_index_by_key(self.processor, candidates, lambda p: p.name, str(raw))
        try:
            profile = candidates[int(raw) - 1]
            annotate_recorded_key(self.processor, profile.name)
        except (ValueError, IndexError):
            self.console.print("[red]Invalid selection[/red]")
            return
        for path in profile.missing_required_fields() or ["account"]:
            current = profile.cluster.get(path, "")
            value = prompt_with_context(
                self.processor,
                f"  {path}",
                default=current,
                module="MD Manager - Cluster Profiles",
                description=f"Fill cluster profile field '{path}'",
            )
            profile.cluster[path] = value
        cp.save_profile(profile, scope="user")
        self.console.print(f"[green]✓ Saved user-scope override for '{profile.name}'[/green]")

    def _clusters_edit_in_editor(self, entries):
        """Open the chosen profile's JSON in $EDITOR. Creates a user-scope
        override from the bundled version if one doesn't exist yet, so the
        bundled file is never modified directly.
        """
        from proprep.md_prep import cluster_profile as cp
        if not entries:
            self.console.print("[yellow]No profiles available to edit[/yellow]")
            return
        for i, e in enumerate(entries, 1):
            source_tag = f"[grey50]({e['source']})[/grey50]"
            self.console.print(f"  {i}. {e['name']} {source_tag}")
        raw = prompt_with_context(
            self.processor,
            "Pick one (or 'c' to cancel)",
            default="1",
            module="MD Manager - Cluster Profiles",
            description="Select cluster profile to edit",
        )
        if raw.strip().lower() == 'c':
            return
        raw = remap_recorded_index_by_key(self.processor, entries, lambda e: e['name'], str(raw))
        try:
            name = entries[int(raw) - 1]['name']
            source = entries[int(raw) - 1]['source']
            annotate_recorded_key(self.processor, name)
        except (ValueError, IndexError):
            self.console.print("[red]Invalid selection[/red]")
            return

        user_path = cp.user_dir() / f"{name}.json"
        if source == 'bundled' and not user_path.exists():
            # Materialize a user-scope copy so edits don't go into the
            # shipped package directory.
            profile = cp.load_profile(name)
            cp.save_profile(profile, scope="user")
            self.console.print(
                f"[grey50]Created user override at {user_path} "
                f"(bundled copy untouched).[/grey50]"
            )
        target = user_path if user_path.exists() else Path(entries[int(raw) - 1]['path'])

        editor = os.environ.get("EDITOR") or os.environ.get("VISUAL") or "vi"
        self.console.print(f"[grey50]Opening {target} in {editor}...[/grey50]")
        try:
            result = subprocess.run([editor, str(target)])
            if result.returncode != 0:
                self.console.print(f"[yellow]Editor exited with code {result.returncode}[/yellow]")
        except FileNotFoundError:
            self.console.print(
                f"[red]Editor '{editor}' not found. Set $EDITOR env var "
                f"or edit {target} directly.[/red]"
            )
            return
        except Exception as e:
            self.console.print(f"[red]Error launching editor: {e}[/red]")
            return

        # Validate the edited profile; warn but don't delete on failure.
        try:
            reloaded = cp.load_profile(name)
            errors = cp.validate(reloaded)
            if errors:
                self.console.print("[yellow]Validation warnings after edit:[/yellow]")
                for e in errors:
                    self.console.print(f"  • {e}")
            else:
                self.console.print(f"[green]✓ '{name}' validated OK[/green]")
        except Exception as e:
            self.console.print(f"[red]Profile became unparseable: {e}[/red]")

    def _clusters_create_from_template(self):
        from proprep.md_prep import cluster_profile as cp
        template_path = cp.bundled_dir() / "_template.json"
        if not template_path.exists():
            self.console.print(f"[red]Template not found: {template_path}[/red]")
            return
        name = prompt_with_context(
            self.processor,
            "Name for new profile (e.g., 'my-cluster')",
            default="",
            module="MD Manager - Cluster Profiles",
            description="New cluster profile name",
        ).strip()
        if not name or name.startswith("_"):
            self.console.print("[red]Invalid name (cannot be empty or start with '_')[/red]")
            return
        with open(template_path, "r") as fh:
            data = json.load(fh)
        # Strip all _notes_* keys recursively from the saved copy
        strip = getattr(self, "_strip_notes", None)
        if strip is None:
            def strip(obj):
                if isinstance(obj, dict):
                    return {k: strip(v) for k, v in obj.items() if not k.startswith("_notes")}
                if isinstance(obj, list):
                    return [strip(i) for i in obj]
                return obj
        clean = strip(data)
        clean["name"] = name
        profile = cp.ClusterProfile.from_json(clean)
        path = cp.save_profile(profile, scope="user")
        self.console.print(f"[green]✓ Created user-scope profile at {path}[/green]")
        self.console.print(
            "[grey50]Edit the JSON directly or use 'clusters → edit' to fill in "
            "site-specific values.[/grey50]"
        )

    def _clusters_export(self, entries):
        from proprep.md_prep import cluster_profile as cp
        if not entries:
            return
        for i, e in enumerate(entries, 1):
            self.console.print(f"  {i}. {e['name']}")
        raw = prompt_with_context(
            self.processor,
            "Pick one",
            default="1",
            module="MD Manager - Cluster Profiles",
            description="Select cluster profile to export",
        )
        raw = remap_recorded_index_by_key(self.processor, entries, lambda e: e['name'], str(raw))
        try:
            name = entries[int(raw) - 1]['name']
            annotate_recorded_key(self.processor, name)
        except (ValueError, IndexError):
            self.console.print("[red]Invalid selection[/red]")
            return
        dest_str = prompt_with_context(
            self.processor,
            "Export path",
            default=f"./{name}.json",
            module="MD Manager - Cluster Profiles",
            description="Export file path",
        )
        path = cp.export_profile(name, Path(dest_str))
        self.console.print(f"[green]✓ Exported to {path}[/green]")

    def _clusters_import(self):
        from proprep.md_prep import cluster_profile as cp
        src_str = prompt_with_context(
            self.processor,
            "Path to profile JSON",
            default="",
            module="MD Manager - Cluster Profiles",
            description="Cluster profile JSON path to import",
        )
        src = Path(src_str.strip())
        if not src.exists():
            self.console.print(f"[red]File not found: {src}[/red]")
            return
        rename = prompt_with_context(
            self.processor,
            "Rename on import? (blank = keep name)",
            default="",
            module="MD Manager - Cluster Profiles",
            description="Optional rename on import",
        ).strip()
        path = cp.import_profile(src, scope="user", rename=rename or None)
        self.console.print(f"[green]✓ Imported to {path}[/green]")

    def _clusters_delete(self, entries):
        from proprep.md_prep import cluster_profile as cp
        user_entries = [e for e in entries if e['source'] == 'user']
        if not user_entries:
            self.console.print("[yellow]No user-scope profiles to delete[/yellow]")
            return
        for i, e in enumerate(user_entries, 1):
            self.console.print(f"  {i}. {e['name']}")
        raw = prompt_with_context(
            self.processor,
            "Pick one (or 'c' to cancel)",
            default="c",
            module="MD Manager - Cluster Profiles",
            description="Select user-scope profile to delete",
        )
        if raw.strip().lower() == 'c':
            return
        raw = remap_recorded_index_by_key(self.processor, user_entries, lambda e: e['name'], str(raw))
        try:
            name = user_entries[int(raw) - 1]['name']
            annotate_recorded_key(self.processor, name)
        except (ValueError, IndexError):
            self.console.print("[red]Invalid selection[/red]")
            return
        if confirm_with_context(
            self.processor,
            f"Really delete '{name}'?",
            default=False,
            module="MD Manager - Cluster Profiles",
            description=f"Confirm delete cluster profile {name}",
        ):
            cp.delete_profile(name, scope="user")
            self.console.print(f"[green]✓ Deleted '{name}'[/green]")

    def _manage_plans(self):
        """Manage run plans: list, show, export, import, delete, remap."""
        from proprep.md_prep import run_plan as rp

        while True:
            self.console.print("\n[bold cyan]===== Run Plan Management =====[/bold cyan]")
            plans = rp.list_plans()
            self.console.print("\n[bold]Available run plans:[/bold]")
            if plans:
                for e in plans:
                    self.console.print(
                        f"  • {e['name']} [grey50]({e['protocol_name']}, "
                        f"cluster={e['cluster_name']}, {e['source']})[/grey50]"
                    )
            else:
                self.console.print("  [grey50](none — save one from Step 4)[/grey50]")

            self.console.print("\n[bold]Actions:[/bold]")
            self.console.print("  1. Show plan details")
            self.console.print("  2. Export plan to file")
            self.console.print("  3. Import plan from file")
            self.console.print("  4. Delete user plan")
            self.console.print("  5. Remap plan to another cluster")
            self.console.print("\n[bold]Navigation:[/bold]")
            self.console.print("  b. ← Back")

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3", "4", "5", "b"],
                default="b",
                module="MD Manager - Run Plans",
                description="Run-plan management action",
                options_map={
                    "1": "Show plan details",
                    "2": "Export plan to file",
                    "3": "Import plan from file",
                    "4": "Delete user plan",
                    "5": "Remap plan to another cluster",
                    "b": "Back",
                },
            )
            if choice == "b":
                return
            if choice == "1":
                self._plans_show_details(plans)
            elif choice == "2":
                self._plans_export(plans)
            elif choice == "3":
                self._plans_import()
            elif choice == "4":
                self._plans_delete(plans)
            elif choice == "5":
                self._plans_remap(plans)

    def _plans_show_details(self, plans):
        from proprep.md_prep import run_plan as rp
        if not plans:
            return
        for i, e in enumerate(plans, 1):
            self.console.print(f"  {i}. {e['name']}")
        raw = prompt_with_context(
            self.processor,
            "Pick one",
            default="1",
            module="MD Manager - Run Plans",
            description="Select run plan from list",
        )
        raw = remap_recorded_index_by_key(self.processor, plans, lambda e: e['name'], str(raw))
        try:
            plan = rp.load_plan(plans[int(raw) - 1]['name'])
            annotate_recorded_key(self.processor, plans[int(raw) - 1]['name'])
        except (ValueError, IndexError):
            self.console.print("[red]Invalid selection[/red]")
            return
        self.console.print(f"\n[bold]{plan.name}[/bold]")
        self.console.print(f"  Protocol: {plan.protocol_name} {plan.protocol_version}")
        self.console.print(f"  Cluster:  {plan.cluster_name}")
        self.console.print(f"  Steps:")
        for step_id, sr in plan.step_resources.items():
            overrides_str = f", overrides={sr.overrides}" if sr.overrides else ""
            self.console.print(f"    {step_id}: {sr.class_name} @ {sr.time_limit or '(default)'}{overrides_str}")

    def _plans_export(self, plans):
        from proprep.md_prep import run_plan as rp
        if not plans:
            return
        for i, e in enumerate(plans, 1):
            self.console.print(f"  {i}. {e['name']}")
        raw = prompt_with_context(
            self.processor,
            "Pick one",
            default="1",
            module="MD Manager - Run Plans",
            description="Select run plan from list",
        )
        raw = remap_recorded_index_by_key(self.processor, plans, lambda e: e['name'], str(raw))
        try:
            name = plans[int(raw) - 1]['name']
            annotate_recorded_key(self.processor, name)
        except (ValueError, IndexError):
            self.console.print("[red]Invalid selection[/red]")
            return
        dest_str = prompt_with_context(
            self.processor,
            "Export path",
            default=f"./{name}.json",
            module="MD Manager - Run Plans",
            description="Run plan export file path",
        )
        path = rp.export_plan(name, Path(dest_str))
        self.console.print(f"[green]✓ Exported to {path}[/green]")

    def _plans_import(self):
        from proprep.md_prep import run_plan as rp
        src_str = prompt_with_context(
            self.processor,
            "Path to plan JSON",
            default="",
            module="MD Manager - Run Plans",
            description="Run plan JSON path to import",
        )
        src = Path(src_str.strip())
        if not src.exists():
            self.console.print(f"[red]File not found: {src}[/red]")
            return
        rename = prompt_with_context(
            self.processor,
            "Rename on import? (blank = keep name)",
            default="",
            module="MD Manager - Run Plans",
            description="Optional rename on plan import",
        ).strip()
        path = rp.import_plan(src, scope="user", rename=rename or None)
        self.console.print(f"[green]✓ Imported to {path}[/green]")

    def _plans_delete(self, plans):
        from proprep.md_prep import run_plan as rp
        user_plans = [e for e in plans if e['source'] == 'user']
        if not user_plans:
            self.console.print("[yellow]No user-scope plans to delete[/yellow]")
            return
        for i, e in enumerate(user_plans, 1):
            self.console.print(f"  {i}. {e['name']}")
        raw = prompt_with_context(
            self.processor,
            "Pick one (or 'c' to cancel)",
            default="c",
            module="MD Manager - Run Plans",
            description="Select user-scope plan to delete",
        )
        if raw.strip().lower() == 'c':
            return
        raw = remap_recorded_index_by_key(self.processor, user_plans, lambda e: e['name'], str(raw))
        try:
            name = user_plans[int(raw) - 1]['name']
            annotate_recorded_key(self.processor, name)
        except (ValueError, IndexError):
            self.console.print("[red]Invalid selection[/red]")
            return
        if confirm_with_context(
            self.processor,
            f"Really delete '{name}'?",
            default=False,
            module="MD Manager - Run Plans",
            description=f"Confirm delete run plan {name}",
        ):
            rp.delete_plan(name, scope="user")
            self.console.print(f"[green]✓ Deleted '{name}'[/green]")

    def _plans_remap(self, plans):
        from proprep.md_prep import run_plan as rp
        from proprep.md_prep import cluster_profile as cp
        if not plans:
            return
        for i, e in enumerate(plans, 1):
            self.console.print(f"  {i}. {e['name']} (cluster={e['cluster_name']})")
        raw = prompt_with_context(
            self.processor,
            "Pick a plan to remap",
            default="1",
            module="MD Manager - Run Plans",
            description="Select plan to remap",
        )
        raw = remap_recorded_index_by_key(self.processor, plans, lambda e: e['name'], str(raw))
        try:
            plan = rp.load_plan(plans[int(raw) - 1]['name'])
            annotate_recorded_key(self.processor, plans[int(raw) - 1]['name'])
        except (ValueError, IndexError):
            self.console.print("[red]Invalid selection[/red]")
            return
        profiles = [e for e in cp.list_profiles() if not e['name'].startswith('_')]
        if not profiles:
            self.console.print("[yellow]No cluster profiles available[/yellow]")
            return
        seen = {}
        for e in profiles:
            seen.setdefault(e['name'], e)
        profiles = list(seen.values())
        self.console.print("\nTarget cluster:")
        for i, e in enumerate(profiles, 1):
            self.console.print(f"  {i}. {e['name']}")
        raw = prompt_with_context(
            self.processor,
            "Pick target",
            default="1",
            module="MD Manager - Run Plans",
            description="Select target cluster for remap",
        )
        raw = remap_recorded_index_by_key(self.processor, profiles, lambda e: e['name'], str(raw))
        try:
            target = cp.load_profile(profiles[int(raw) - 1]['name'])
            annotate_recorded_key(self.processor, profiles[int(raw) - 1]['name'])
        except (ValueError, IndexError):
            self.console.print("[red]Invalid selection[/red]")
            return
        decisions = rp.remap_for_cluster(plan, target)
        self.console.print("\n[bold]Remap preview:[/bold]")
        for d in decisions:
            tag = "[green]✓[/green]" if d['status'] == 'kept' else "[yellow]~[/yellow]"
            self.console.print(
                f"  {tag} {d['step_id']}: {d['current_class']} → "
                f"{d['suggested_class']}  [grey50]({d['reason']})[/grey50]"
            )
        if not confirm_with_context(
            self.processor,
            "\nApply remap?",
            default=True,
            module="MD Manager - Run Plans",
            description="Apply plan remap to target cluster",
        ):
            return
        new_plan = rp.apply_remap(plan, target, decisions)
        rename = prompt_with_context(
            self.processor,
            "Name for remapped plan",
            default=new_plan.name,
            module="MD Manager - Run Plans",
            description="Name for remapped plan",
        ).strip()
        new_plan.name = rename or new_plan.name
        path = rp.save_plan(new_plan, scope="user")
        self.console.print(f"[green]✓ Saved remapped plan: {path}[/green]")

    def _manage_library(self):
        """Manage template and workflow library."""
        while True:
            self.console.print("\n[bold cyan]===== Library Management =====[/bold cyan]")
            self.console.print("\n[bold]Import:[/bold]")
            self.console.print("  1. Import .mdin template")
            self.console.print("  2. Import .json protocol")
            self.console.print("\n[bold]Template Library:[/bold]")
            self.console.print("  3. View template statistics")
            self.console.print("  4. Preview template content")
            self.console.print("  5. Delete user templates")
            self.console.print("\n[bold]Protocol Library:[/bold]")
            self.console.print("  6. View protocol statistics")
            self.console.print("  7. Delete custom protocols")
            self.console.print("\n[bold]Navigation:[/bold]")
            self.console.print("  b. ← Back to main menu")

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3", "4", "5", "6", "7", "b"],
                default="b",
                module="MD Manager - Library Management",
                description="Library management options",
                options_map={
                    "1": "Import .mdin template",
                    "2": "Import .json protocol",
                    "3": "View template statistics",
                    "4": "Preview template content",
                    "5": "Delete user templates",
                    "6": "View protocol statistics",
                    "7": "Delete custom protocols",
                    "b": "← Back to main menu"
                }
            )

            if choice == "1":
                self._import_template()
            elif choice == "2":
                self._import_workflow()
            elif choice == "3":
                self._view_library_statistics()
            elif choice == "4":
                self._preview_library_templates()
            elif choice == "5":
                self._delete_template_from_library()
            elif choice == "6":
                self._view_workflow_statistics()
            elif choice == "7":
                self._delete_workflows_from_library()
            elif choice == "b":
                break

    def _view_workflow_statistics(self):
        """View statistics about the workflow library."""
        self.console.print("\n[bold cyan]Protocol Library Statistics[/bold cyan]\n")

        workflows = self.user_data_manager.list_workflows()

        # Count by source
        builtin_count = sum(1 for w in workflows.values() if w.get('source') == 'builtin')
        # Accept both 'custom' and 'user' for backward compatibility
        custom_count = sum(1 for w in workflows.values() if w.get('source') in ['custom', 'user'])

        from rich.table import Table
        stats_table = Table(title="Protocol Counts")
        stats_table.add_column("Category", style="bright_blue")
        stats_table.add_column("Count", style="green")

        stats_table.add_row("Built-in Protocols", str(builtin_count))
        stats_table.add_row("Custom Protocols", str(custom_count))
        stats_table.add_row("Total Protocols", str(len(workflows)))

        self.console.print(stats_table)

        # Show custom workflows if any
        if custom_count > 0:
            self.console.print("\n[bold]Custom Protocols:[/bold]")
            for workflow_id, metadata in workflows.items():
                # Accept both 'custom' and 'user' for backward compatibility
                if metadata.get('source') in ['custom', 'user']:
                    name = metadata.get('name', workflow_id)
                    step_count = metadata.get('step_count', '?')
                    self.console.print(f"  • {name} ({step_count} steps)")

    def _delete_workflows_from_library(self):
        """Delete custom workflows from the library."""
        self.console.print("\n[bold]Delete Custom Protocols[/bold]")

        # Get only custom workflows (not builtin)
        all_workflows = self.user_data_manager.list_workflows()
        # Accept both 'custom' and 'user' for backward compatibility
        custom_workflows = {k: v for k, v in all_workflows.items() if v.get('source') in ['custom', 'user']}

        if not custom_workflows:
            self.console.print("[yellow]No custom protocols available to delete[/yellow]")
            return

        # Display custom workflows for selection
        self.console.print("\n[bold]Select protocols to delete:[/bold]")
        workflow_choices = {}
        choice_num = 1

        for workflow_id, metadata in custom_workflows.items():
            workflow_choices[str(choice_num)] = workflow_id
            name = metadata.get('name', workflow_id)
            step_count = metadata.get('step_count', '?')
            self.console.print(f"  {choice_num}. {name} ({step_count} steps)")
            choice_num += 1

        self.console.print(f"  0. Cancel")
        self.console.print("[grey50]Tip: Enter comma-separated numbers and/or ranges (e.g., '1,3-5,7') or 'all' to delete all[/grey50]")

        valid_choices = list(workflow_choices.keys()) + ["0", "all"]
        user_input = prompt_with_context(
            self.processor,
            "Select protocol(s) to delete",
            default="0",
            module="MD Manager - Protocol Library",
            description="Select protocols to delete (comma-separated, ranges, 'all', or 0 to cancel)"
        )

        if user_input == "0":
            return

        # Parse comma-separated selections, ranges, or 'all'
        try:
            # Handle 'all' option
            if user_input.lower() == "all":
                selected_workflow_ids = list(workflow_choices.values())
            else:
                selected_choices = self._parse_numeric_selection(user_input, valid_choices)
                if selected_choices is None:
                    return
                selected_workflow_ids = [workflow_choices[c] for c in selected_choices]

            if not selected_workflow_ids:
                return

            # Identify templates used by selected workflows
            self.console.print("\n[grey50]Analyzing template dependencies...[/grey50]")
            all_template_ids = set()
            orphaned_template_ids = []
            shared_template_ids = []

            # Collect all templates used by selected workflows
            for workflow_id in selected_workflow_ids:
                template_ids = self.user_data_manager.get_workflow_template_ids(workflow_id)
                all_template_ids.update(template_ids)

            # Check which templates are orphaned (only used by workflows being deleted)
            for template_id in all_template_ids:
                workflows_using_template = self.user_data_manager.find_template_usage(template_id)
                # Check if ALL workflows using this template are in the deletion list
                if all(wid in selected_workflow_ids for wid in workflows_using_template):
                    orphaned_template_ids.append(template_id)
                else:
                    shared_template_ids.append(template_id)

            # Show what will be deleted and confirm
            self.console.print(f"\n[bold]Protocols to delete:[/bold]")
            for workflow_id in selected_workflow_ids:
                workflow_metadata = custom_workflows[workflow_id]
                self.console.print(f"  - {workflow_metadata.get('name', workflow_id)}")

            if orphaned_template_ids:
                self.console.print(f"\n[bold]Associated templates to delete (orphaned):[/bold]")
                self.console.print(f"  [green]✓[/green] {len(orphaned_template_ids)} template{'s' if len(orphaned_template_ids) != 1 else ''} used only by these protocols")

            if shared_template_ids:
                self.console.print(f"\n[bold]Shared templates (will NOT be deleted):[/bold]")
                self.console.print(f"  [yellow]⚠[/yellow] {len(shared_template_ids)} template{'s' if len(shared_template_ids) != 1 else ''} used by other protocols")

            # Confirm deletion
            deletion_summary = f"{len(selected_workflow_ids)} protocol{'s' if len(selected_workflow_ids) != 1 else ''}"
            if orphaned_template_ids:
                deletion_summary += f" and {len(orphaned_template_ids)} orphaned template{'s' if len(orphaned_template_ids) != 1 else ''}"

            confirm_str = prompt_with_context(
                self.processor,
                f"[red]Delete {deletion_summary}?[/red]",
                choices=["y", "n"],
                default="n",
                module="MD Manager - Protocol Library",
                description="Confirm protocol deletion",
                options_map={"y": "Yes, delete", "n": "No, cancel"}
            )

            if confirm_str.lower() == "y":
                delete_count = 0
                for workflow_id in selected_workflow_ids:
                    workflow_metadata = custom_workflows[workflow_id]
                    try:
                        # Delete workflow file
                        workflow_path = self.user_data_manager.workflow_base_dir / workflow_metadata['workflow_path']
                        if workflow_path.exists():
                            workflow_path.unlink()

                        # Remove from metadata
                        del self.user_data_manager.workflow_metadata[workflow_id]

                        self.console.print(f"[green]✓ Deleted protocol '{workflow_metadata.get('name', workflow_id)}'[/green]")
                        delete_count += 1
                    except Exception as e:
                        self.console.print(f"[red]Error deleting protocol '{workflow_metadata.get('name', workflow_id)}': {e}[/red]")

                # Save workflow metadata changes
                if delete_count > 0:
                    try:
                        self.user_data_manager._save_metadata(
                            self.user_data_manager.workflow_metadata,
                            self.user_data_manager.workflow_metadata_file
                        )
                    except Exception as e:
                        self.console.print(f"[red]Error saving protocol metadata: {e}[/red]")

                # Delete orphaned templates
                template_delete_count = 0
                if orphaned_template_ids:
                    self.console.print(f"\n[grey50]Deleting orphaned templates...[/grey50]")
                    for template_id in orphaned_template_ids:
                        try:
                            # Find and delete template file
                            custom_template_dir = self.user_data_manager.user_template_dir / "custom"
                            if custom_template_dir.exists():
                                # Template files are named: <uuid>_<name>.json
                                for template_file in custom_template_dir.glob(f"{template_id}_*.json"):
                                    if template_file.exists():
                                        template_file.unlink()
                                        template_delete_count += 1
                                        break
                        except Exception as e:
                            self.console.print(f"[yellow]Warning: Could not delete template {template_id}: {e}[/yellow]")

                # Show final summary
                if delete_count > 0 or template_delete_count > 0:
                    summary_parts = []
                    if delete_count > 0:
                        summary_parts.append(f"{delete_count} protocol{'s' if delete_count != 1 else ''}")
                    if template_delete_count > 0:
                        summary_parts.append(f"{template_delete_count} template{'s' if template_delete_count != 1 else ''}")

                    self.console.print(f"[green]✓ Successfully deleted {' and '.join(summary_parts)}[/green]")

                    if shared_template_ids:
                        self.console.print(f"[grey50]  Preserved {len(shared_template_ids)} shared template{'s' if len(shared_template_ids) != 1 else ''}[/grey50]")
            else:
                self.console.print("[grey50]Deletion cancelled[/grey50]")

        except Exception as e:
            self.console.print(f"[red]Error parsing selection: {e}[/red]")

    def _view_library_statistics(self):
        """View statistics about the template library."""
        self.console.print("\n[bold cyan]Template Library Statistics[/bold cyan]\n")
        
        templates = self.user_data_manager.list_templates()
        
        # Categorize templates
        builtin_count = sum(1 for t in templates.values() if t.get('source') == 'builtin')
        # Custom templates include: 'custom', 'from_workflow_step', 'copied_from_builtin', etc.
        custom_count = sum(1 for t in templates.values() if t.get('source') not in ['builtin', 'modified'])
        modified_count = sum(1 for t in templates.values() if t.get('source') == 'modified')
        
        # Count by type
        type_counts = {}
        for metadata in templates.values():
            sim_type = metadata.get('simulation_type', 'unknown')
            type_counts[sim_type] = type_counts.get(sim_type, 0) + 1
            
        # Display statistics
        from rich.table import Table
        
        table = Table(title="Template Library Overview")
        table.add_column("Category", style="bright_blue")
        table.add_column("Count", justify="right")
        
        table.add_row("Built-in Templates", str(builtin_count))
        table.add_row("Custom Templates", str(custom_count))
        table.add_row("Modified Templates", str(modified_count))
        table.add_row("[bold]Total", f"[bold]{len(templates)}")
        
        self.console.print(table)
        
        # Type breakdown
        self.console.print("\n[bold]Templates by Type:[/bold]")
        for sim_type in ["minimization", "heating", "equilibration", "production"]:
            if sim_type in type_counts:
                self.console.print(f"  {sim_type.title()}: {type_counts[sim_type]}")
        if "unknown" in type_counts or "other" in type_counts:
            other_count = type_counts.get("unknown", 0) + type_counts.get("other", 0)
            self.console.print(f"  Other: {other_count}")
            
        # Directory paths
        self.console.print(f"\n[bold]Library Location:[/bold]")
        self.console.print(f"  Base: {self.user_data_manager.template_base_dir}")
        self.console.print(f"  User: {self.user_data_manager.user_template_dir}")
        
    def _browse_for_template_import(self):
        """Browse for .mdin template files to import.

        Thin wrapper over the shared file browser; the recursive `find` command
        is preserved via extra_commands. Returns a Path or None on cancel.
        """
        from pathlib import Path
        from datetime import datetime
        from proprep.utils.file_browser import file_browser

        def _ctime_detail(p):
            try:
                return "created " + datetime.fromtimestamp(os.path.getctime(p)).strftime("%m/%d/%Y")
            except OSError:
                return ""

        extra = {
            "find": ("Search recursively for .mdin files",
                     lambda cur: self._find_mdin_files(Path(cur))),
        }
        return file_browser(
            directory=str(Path.cwd()),
            extensions=[".mdin"],
            console=self.console,
            processor=self.processor,
            label="mdin template",
            entry_detail=_ctime_detail,
            path_factory=Path,
            extra_commands=extra,
            module="MD Manager - Template Browser",
        )


                
    def _find_mdin_files(self, current_dir):
        """Search recursively for .mdin files from current directory."""
        from pathlib import Path
        
        self.console.print(f"\n[bold cyan]Searching for .mdin files in {current_dir}...[/bold cyan]")
        
        # Find all .mdin files recursively
        mdin_files = list(current_dir.rglob("*.mdin"))
        
        if not mdin_files:
            self.console.print("[grey50]No .mdin files found[/grey50]")
            return None
            
        self.console.print(f"[green]Found {len(mdin_files)} .mdin files:[/green]")
        
        # Display found files
        for i, file_path in enumerate(mdin_files, 1):
            relative_path = file_path.relative_to(current_dir)
            creation_time = datetime.fromtimestamp(file_path.stat().st_ctime)
            date_str = creation_time.strftime("%m/%d/%Y")
            self.console.print(f"  {i:2}. {relative_path} (created {date_str})")
            
        # Let user select
        while True:
            choice = prompt_with_context(
                self.processor, f"Select file (1-{len(mdin_files)}) or 'cancel'",
                default="cancel", module="MD Manager - File Search",
                description="Select from recursively-found mdin files",
            ).strip()
            
            if choice.lower() == 'cancel':
                return None
                
            try:
                file_num = int(choice)
                if 1 <= file_num <= len(mdin_files):
                    selected_file = mdin_files[file_num - 1]
                    self.console.print(f"[green]Selected: {selected_file}[/green]")
                    return selected_file
                else:
                    self.console.print(f"[red]Please enter a number between 1 and {len(mdin_files)}[/red]")
            except ValueError:
                self.console.print("[red]Please enter a valid number or 'cancel'[/red]")
                
    def _enhanced_find_for_import(self, current_dir, pattern):
        """Enhanced find specifically for template import."""
        import glob
        
        self.console.print(f"[grey50]Searching for: {pattern}[/grey50]")
        
        matches = []
        try:
            if "*" in pattern or "?" in pattern:
                # Glob pattern search
                for match in current_dir.rglob(pattern):
                    if match.is_file() and match.suffix.lower() == ".mdin":
                        matches.append(match)
            else:
                # Text pattern search
                for match in current_dir.rglob(f"*{pattern}*.mdin"):
                    if match.is_file():
                        matches.append(match)
                        
            if not matches:
                self.console.print("[yellow]No matching .mdin files found[/yellow]")
                return None
                
            if len(matches) == 1:
                # Single match - return directly
                self.console.print(f"[green]Found: {matches[0]}[/green]")
                return matches[0]
            else:
                # Multiple matches - let user choose
                self.console.print(f"\n[bold]Found {len(matches)} matching files:[/bold]")
                for i, match in enumerate(matches[:20], 1):  # Limit to 20 results
                    rel_path = match.relative_to(current_dir)
                    self.console.print(f"  {i}. {rel_path}")
                    
                if len(matches) > 20:
                    self.console.print(f"  ... and {len(matches) - 20} more (showing first 20)")
                    
                self.console.print("  0. Cancel")

                valid_choices = [str(i) for i in range(1, min(21, len(matches)+1))] + ["0"]

                # Build options map
                options_map = {}
                for i, match in enumerate(matches[:20], 1):
                    rel_path = match.relative_to(current_dir)
                    options_map[str(i)] = str(rel_path)
                options_map["0"] = "Cancel"

                choice = prompt_with_context(
                    self.processor,
                    "Select file",
                    choices=valid_choices,
                    default="0",
                    module="MD Manager - File Browser",
                    description="Select file from search results",
                    options_map=options_map
                )
                
                # Replay by basename so changed search results can't mis-pick.
                choice = remap_recorded_index(self.processor, matches, choice)
                if choice == "0":
                    return None
                else:
                    selected_file = matches[int(choice) - 1]
                    annotate_selected_path(self.processor, selected_file)
                    self.console.print(f"[green]Selected: {selected_file}[/green]")
                    return selected_file
                    
        except Exception as e:
            self.console.print(f"[red]Search error: {e}[/red]")
            return None
        
    def _preview_template_content(self, selected_templates):
        """Preview the content of selected templates."""
        if not selected_templates:
            self.console.print("[yellow]No templates selected to preview[/yellow]")
            return
            
        self.console.print("\n[bold]Preview Template Content:[/bold]")
        
        # Display current templates with numbers
        templates = self.user_data_manager.list_templates()
        for i, template_id in enumerate(selected_templates, 1):
            template_name = templates.get(template_id, {}).get('name', template_id)
            self.console.print(f"  {i:2}. {template_name}")
            
        self.console.print(f"\n[grey50]0. Cancel preview[/grey50]")

        # Get user selection
        valid_choices = [str(i) for i in range(1, len(selected_templates)+1)] + ["0"]

        # Build options map
        options_map = {}
        for i, template_id in enumerate(selected_templates, 1):
            template_name = templates.get(template_id, {}).get('name', template_id)
            options_map[str(i)] = template_name
        options_map["0"] = "Cancel preview"

        choice = prompt_with_context(
            self.processor,
            "Select template to preview",
            choices=valid_choices,
            default="0",
            module="MD Manager - Template Library",
            description="Select template to preview content",
            options_map=options_map
        )
        
        if choice == "0":
            return
            
        # Preview the selected template
        template_index = int(choice) - 1
        template_id = selected_templates[template_index]
        template_metadata = templates.get(template_id, {})
        
        # Create a temporary config for preview
        temp_config = SimulationConfig(
            name=template_metadata.get('name', template_id),
            template_id=template_id,
            mdin_path=template_metadata.get('template_path', ''),
            engine="preview"
        )
        
        # Get the controller to use for formatting
        from .amber_controller import AmberController
        controller = AmberController(processor=self.processor)
        
        self._display_template_details(template_metadata, temp_config, controller)
    
    def _select_existing_template(self, controller):
        """Display categorized templates for selection using existing AmberController method."""
        templates = self.user_data_manager.list_templates()
        
        if not templates:
            self.console.print("[red]No templates available.[/red]")
            return None
            
        # Use AmberController's excellent categorized display with numbering
        self.console.print("\n[bold]Available Templates:[/bold]")
        template_choices = controller._display_categorized_templates(templates, show_numbers=True)
        
        # Get user selection (support comma-separated multiple selections)
        self.console.print(f"\n[grey50]0. Cancel selection[/grey50]")
        self.console.print(f"[grey50]Tip: Enter comma-separated numbers for multiple templates (e.g., '1,3,5')[/grey50]")

        valid_choices = list(template_choices.keys()) + ["0"]

        # Build options map
        options_map = {}
        for choice_num, template_id in template_choices.items():
            template_name = templates.get(template_id, {}).get('name', template_id)
            options_map[choice_num] = template_name
        options_map["0"] = "Cancel selection"

        user_input = prompt_with_context(
            self.processor,
            "Select template(s)",
            default="0",
            module="MD Manager - Template Selection",
            description="Select templates (comma-separated for multiple)",
            options_map=options_map
        )
        
        if user_input == "0":
            return None
            
        # Parse comma-separated selections
        selected_ids = []
        try:
            choices = [choice.strip() for choice in user_input.split(',')]
            for choice in choices:
                if choice not in valid_choices:
                    self.console.print(f"[red]Invalid choice: {choice}[/red]")
                    return None
                if choice != "0":
                    template_id = template_choices[choice]
                    selected_ids.append(template_id)
                    
            if not selected_ids:
                return None
                
            # Return list of template IDs (single item if only one selected)
            return selected_ids if len(selected_ids) > 1 else selected_ids[0]
            
        except Exception as e:
            self.console.print(f"[red]Error parsing selection: {e}[/red]")
            return None
        
    def _create_or_modify_template_with_wizard(self, controller):
        """Integrate with existing template creation functionality from AmberController."""
        self.console.print("\n[bold]Create or Modify Template[/bold]")
        self.console.print("[grey50]This will use the interactive template creation/modification process[/grey50]\n")
        
        try:
            # Get initial template count to detect new template creation
            templates_before = set(self.user_data_manager.list_templates().keys())
            
            # Call AmberController's template creation method
            controller._create_single_template()
            
            # Refresh and check if a new template was created
            self.user_data_manager = UserDataManager(console=self.console)  # Refresh the manager
            templates_after = set(self.user_data_manager.list_templates().keys())
            new_templates = templates_after - templates_before
            
            if new_templates:
                # Get the first (and likely only) new template ID
                new_template_id = list(new_templates)[0]
                templates = self.user_data_manager.list_templates()
                template_name = templates.get(new_template_id, {}).get('name', new_template_id)
                
                # Ask if the template should be added to the queue
                add_to_queue_str = prompt_with_context(
                    self.processor,
                    f"\nAdd '{template_name}' to the simulation queue?",
                    choices=["y", "n"],
                    default="y",
                    module="MD Manager - Template Creation",
                    description="Add newly created template to queue",
                    options_map={"y": "Yes, add to queue", "n": "No, just save template"}
                )
                add_to_queue = add_to_queue_str.lower() == "y"
                
                if add_to_queue:
                    self.console.print(f"[green]✓ Added template to queue: {template_name}[/green]")
                    return new_template_id
                else:
                    self.console.print("[grey50]Template created but not added to queue[/grey50]")
                    return None
            else:
                # Even if we can't detect a new template automatically, 
                # refresh the template list and check if any template was created recently
                self.user_data_manager = UserDataManager(console=self.console)  # Refresh the manager
                all_templates = self.user_data_manager.list_templates()
                
                if all_templates:
                    # Sort by creation date to find the most recent
                    sorted_templates = sorted(
                        all_templates.items(),
                        key=lambda x: x[1].get('created_date', '0000-00-00T00:00:00'),
                        reverse=True
                    )
                    
                    most_recent_id, most_recent_meta = sorted_templates[0]
                    template_name = most_recent_meta.get('name', most_recent_id)
                    


                    # Ask the user if this is the template they just created
                    is_new_str = prompt_with_context(
                        self.processor,
                        f"\nWas '{template_name}' the template you just created/modified?",
                        choices=["y", "n"],
                        default="y",
                        module="MD Manager - Template Creation",
                        description="Confirm template just created",
                        options_map={"y": "Yes, that's the one", "n": "No, different template"}
                    )
                    is_new = is_new_str.lower() == "y"
                    
                    if is_new:
                        # Ask if the template should be added to the queue
                        add_to_queue_str = prompt_with_context(
                            self.processor,
                            f"Add '{template_name}' to the simulation queue?",
                            choices=["y", "n"],
                            default="y",
                            module="MD Manager - Template Creation",
                            description="Add template to queue",
                            options_map={"y": "Yes, add to queue", "n": "No, just save template"}
                        )
                        add_to_queue = add_to_queue_str.lower() == "y"
                        
                        if add_to_queue:
                            self.console.print(f"[green]✓ Added template to queue: {template_name}[/green]")
                            return most_recent_id
                        else:
                            self.console.print("[grey50]Template created but not added to queue[/grey50]")
                            return None
                
                self.console.print("[yellow]Template operation was cancelled or failed[/yellow]")
                return None
                
        except Exception as e:
            self.console.print(f"[red]Error with template operation: {e}[/red]")
            return None
            
    def _modify_existing_template(self, controller):
        """Template modification with custom copy creation for builtins."""
        # First select template to modify
        selection = self._select_existing_template(controller)
        if selection is None:
            return None
            
        template_id, metadata = selection
        self.console.print(f"\n[bold]Modifying Template: {metadata.get('name', template_id)}[/bold]")
        
        # Check if this is a builtin template
        is_builtin = template_id.startswith('builtin_')
        
        if is_builtin:
            self.console.print("[yellow]This is a builtin template. Modifications will create a custom copy.[/yellow]")
            
        # Show modification options
        self.console.print("\n[bold]Modification Options:[/bold]")
        self.console.print("  1. Edit metadata only (name, description, priority)")
        self.console.print("  2. Edit template content (launch editor)")
        self.console.print("  3. Edit both metadata and content")
        self.console.print("  4. Cancel")

        choice = prompt_with_context(
            self.processor,
            "Select option",
            choices=["1", "2", "3", "4"],
            default="4",
            module="MD Manager - Template Modification",
            description="Select modification type",
            options_map={
                "1": "Edit metadata only (name, description, priority)",
                "2": "Edit template content (launch editor)",
                "3": "Edit both metadata and content",
                "4": "Cancel"
            }
        )
        
        if choice == "4":
            return None
            
        try:
            # Create a temporary config to work with
            temp_config = SimulationConfig(
                name=f"temp_{template_id}",
                template_id=template_id,
                mdin_path=metadata.get('template_path', ''),
                engine="sander"  # Default, will be changed later
            )
            
            # Make a working copy of metadata
            working_metadata = metadata.copy()
            
            # Perform the modification based on choice
            if choice == "1":
                self._modify_template_metadata(temp_config, working_metadata)
                # Save changes - this handles builtin vs custom properly
                self._save_template_metadata_changes(temp_config, working_metadata)
                
            elif choice == "2":
                self._edit_template_directly(temp_config, working_metadata)
                
            elif choice == "3":
                self._modify_template_metadata(temp_config, working_metadata)
                # Save metadata changes first - this handles builtin vs custom properly
                self._save_template_metadata_changes(temp_config, working_metadata)
                # Then edit content if it's still a custom template
                self._edit_template_directly(temp_config, working_metadata)
                
            # Get the final template ID and metadata after modifications
            # Return the final template ID (may have changed if builtin was converted to custom)
            return temp_config.template_id
            
        except Exception as e:
            self.console.print(f"[red]Error modifying template: {e}[/red]")
            return None
            
    def _configure_individual_engines(self, has_gpus=False):
        """Configure engine and hardware for each queued simulation individually."""
        active_queue = self._get_active_queue()
        if not active_queue:
            self.console.print("[yellow]No simulations in queue to configure[/yellow]")
            return

        # When a cluster profile is loaded, the per-sim editor uses the
        # profile's resource palette instead of the legacy engine picker.
        if self._cluster_profile is not None:
            self._configure_individual_engines_with_profile(active_queue)
            return

        self.console.print("\n[bold]Individual Engine Configuration[/bold]")
        self.console.print("[grey50]Configure engine and hardware for each simulation[/grey50]\n")

        # Get available engines
        available_engines = self._get_available_engines()

        for i, config in enumerate(active_queue):
            self.console.print(f"\n[cyan]Simulation {i+1}: {config.name}[/cyan]")
            
            # Get template name instead of UUID
            template_display = config.template_id
            try:
                template_data = self.user_data_manager.load_custom_template(config.template_id)
                if template_data and 'name' in template_data:
                    template_display = template_data['name']
            except:
                pass  # Keep UUID if template not found
                
            self.console.print(f"Template: {template_display}")
            self.console.print(f"\nCurrent engine: [yellow]{config.engine}[/yellow]")

            # Show engine selection with current engine marked
            self.console.print("\nAvailable engines:")
            for j, engine in enumerate(available_engines, 1):
                if engine == config.engine:
                    self.console.print(f"  {j}. {engine} [grey50](current)[/grey50]")
                else:
                    self.console.print(f"  {j}. {engine}")

            engine_choices = [str(j) for j in range(1, len(available_engines)+1)]

            # Build options map
            options_map = {}
            for j, engine in enumerate(available_engines, 1):
                if engine == config.engine:
                    options_map[str(j)] = f"{engine} (current)"
                else:
                    options_map[str(j)] = engine
            options_map["skip"] = "Keep current engine"

            choice = prompt_with_context(
                self.processor,
                "Select engine",
                choices=engine_choices + ["skip"],
                default="skip",
                module="MD Manager - Engine Configuration",
                description="Select simulation engine",
                options_map=options_map
            )

            if choice != "skip":
                new_engine = available_engines[int(choice)-1]
                if new_engine != config.engine:
                    config.engine = new_engine
                    self.console.print(f"[green]✓ Updated to {new_engine}[/green]")

                    # Configure hardware resources if applicable
                    if "cuda" in new_engine.lower():
                        config.hardware_config = self._configure_gpu_resources()
                    elif "mpi" in new_engine.lower():
                        config.hardware_config = self._configure_mpi_resources()
                else:
                    self.console.print("[grey50]Engine unchanged[/grey50]")
        
        # Force sync to workspace after all changes
        if hasattr(self.simulation_queue, '_sync_to_workspace'):
            self.simulation_queue._sync_to_workspace()
                        
    def _configure_individual_engines_with_profile(self, active_queue):
        """Profile-aware per-step editor: class + wall time + overrides + amber_flags.

        Shows the fully-resolved config for each step before prompting for edits,
        so the user can see exactly what the generated script will use.
        """
        profile = self._cluster_profile
        class_names = list(profile.resource_classes.keys())
        cpu_default = profile.defaults.get("cpu_class", "")
        gpu_default = profile.defaults.get("gpu_class", "")

        self.console.print("\n[bold]Per-step Resource Configuration[/bold]")
        self.console.print(
            f"[grey50]Profile: {profile.name} · "
            f"Available classes: {', '.join(class_names)}[/grey50]\n"
        )

        # Compute the classifier hint per step (CPU vs GPU) so suggested
        # defaults match Amber recommendations.
        npt_indices = []
        npt_counter = 0
        for sim in active_queue:
            params = self._get_simulation_mdin_params(sim)
            try:
                imin = int(params.get('imin', 0))
            except (ValueError, TypeError):
                imin = 0
            try:
                ntp = int(params.get('ntp', 0))
            except (ValueError, TypeError):
                ntp = 0
            if imin == 0 and ntp > 0:
                npt_indices.append(npt_counter)
                npt_counter += 1
            else:
                npt_indices.append(-1)
        npt_total = npt_counter

        for i, sim in enumerate(active_queue):
            target, reason = self._classify_for_engine(
                sim, npt_index=npt_indices[i], npt_total=npt_total)
            suggested_class = cpu_default if target == 'cpu' else gpu_default
            hw = sim.hardware_config or {}
            current_class = hw.get('resource_class', suggested_class)
            current_wall = hw.get('time_limit', '')
            current_overrides = dict(hw.get('overrides', {}))
            current_flags = list(hw.get('amber_flags', []))

            # Show fully resolved config before editing
            self._display_resolved_step_config(
                step_index=i + 1,
                sim=sim,
                current_class=current_class,
                current_wall=current_wall,
                current_overrides=current_overrides,
                current_flags=current_flags,
                classifier_target=target,
                classifier_reason=reason,
                suggested_class=suggested_class,
            )

            self.console.print(
                "[grey50]  k = keep as-is (Enter)  ·  c = change resource class  ·  "
                "w = change wall time  ·  s = skip remaining steps[/grey50]"
            )
            self.console.print(
                "[grey50]  To change partition / memory / GRES / modules, edit the "
                "class itself: Main Menu → clusters → Edit JSON.[/grey50]"
            )
            action = prompt_with_context(
                self.processor,
                "Action",
                choices=["k", "c", "w", "s"],
                default="k",
                module="MD Manager - Step Resources",
                description="Per-step resource action",
                options_map={
                    "k": "Keep as-is, next step",
                    "c": "Change resource class",
                    "w": "Change wall time",
                    "s": "Skip remaining steps",
                },
            )
            if action == "s":
                break
            if action == "k":
                if 'resource_class' not in hw:
                    cls_def = profile.resolve_class(current_class)
                    sim.engine = cls_def.get('binary', sim.engine)
                    sim.hardware_config = {
                        'resource_class': current_class,
                        'cluster_profile': profile.name,
                        'time_limit': current_wall or cls_def.get('default_time', ''),
                    }
                continue
            if action == "c":
                self.console.print("  Available classes:")
                for j, cn in enumerate(class_names, 1):
                    marker = " [grey50](current)[/grey50]" if cn == current_class else ""
                    self.console.print(f"    {j}. {cn}{marker}")
                class_options_map = {str(i + 1): cn for i, cn in enumerate(class_names)}
                raw = prompt_with_context(
                    self.processor,
                    "  Pick class",
                    default=str(class_names.index(current_class) + 1 if current_class in class_names else 1),
                    module="MD Manager - Step Resources",
                    description="Select resource class for step",
                    options_map=class_options_map,
                )
                raw = remap_recorded_index_by_key(self.processor, class_names, lambda c: c, str(raw))
                try:
                    current_class = class_names[int(raw) - 1]
                    annotate_recorded_key(self.processor, current_class)
                    # When class changes, reset wall time to new class's
                    # default unless the user had an explicit non-default
                    # wall time already.
                    if not current_wall:
                        current_wall = profile.resolve_class(current_class).get('default_time', '')
                except (ValueError, IndexError):
                    self.console.print("[red]Invalid selection — keeping current[/red]")
            if action == "w":
                from proprep.md_prep.cluster_profile import parse_slurm_time_to_seconds
                try:
                    cls_for_cap = profile.resolve_class(current_class)
                    default_wall = current_wall or cls_for_cap.get('default_time', '')
                except KeyError:
                    cls_for_cap = {}
                    default_wall = current_wall
                max_time = cls_for_cap.get('max_time', '')
                max_sec = parse_slurm_time_to_seconds(max_time)
                while True:
                    entered = prompt_with_context(
                        self.processor,
                        "  Wall time (HH:MM:SS)",
                        default=default_wall,
                        module="MD Manager - Step Resources",
                        description="Step wall time (HH:MM:SS)",
                    ).strip()
                    entered_sec = parse_slurm_time_to_seconds(entered)
                    if max_sec is not None and entered_sec is not None and entered_sec > max_sec:
                        self.console.print(
                            f"[red]  {entered} exceeds partition max of {max_time} "
                            f"for class '{current_class}' — SLURM will reject this job.[/red]"
                        )
                        continue
                    current_wall = entered
                    break

            # Commit this step's edits
            try:
                cls_def = profile.resolve_class(current_class)
            except KeyError:
                self.console.print(f"[red]Class '{current_class}' no longer exists in profile[/red]")
                continue
            sim.engine = cls_def.get('binary', sim.engine)
            sim.hardware_config = {
                'resource_class': current_class,
                'cluster_profile': profile.name,
                'time_limit': current_wall,
            }

        self.console.print(f"\n[green]✓ Per-step configuration saved[/green]")

    def _display_resolved_step_config(
        self, step_index, sim, current_class, current_wall,
        current_overrides, current_flags,
        classifier_target, classifier_reason, suggested_class,
    ):
        """Show the final resolved config for a step so the user can audit it.

        The first two rows (class, wall time) are step-editable via the action
        prompt. The rest come from the resource class definition and are shown
        here for transparency — to change them, the user must edit the class
        in the 'clusters' submenu.
        """
        from rich.panel import Panel
        profile = self._cluster_profile

        lines = []
        step_name = sim.step_name or sim.name
        lines.append(f"[bold cyan]Step {step_index}:[/bold cyan] {step_name}")
        lines.append(
            f"[grey50]  Classifier: {classifier_target.upper()} — "
            f"{classifier_reason}  ·  suggested: {suggested_class}[/grey50]"
        )
        lines.append("")
        lines.append("[bold]Editable here:[/bold]")
        lines.append(f"  Resource class: [bright_blue]{current_class}[/bright_blue]")
        try:
            cls_def = profile.resolve_class(current_class)
        except KeyError:
            lines.append("  [red](class not found in profile)[/red]")
            self.console.print(Panel("\n".join(lines), expand=False, border_style="cyan"))
            return
        effective_wall = current_wall or cls_def.get('default_time', '')
        lines.append(f"  Wall time:      {effective_wall}")
        lines.append("")
        lines.append("[bold]From the class definition[/bold] [grey50](edit via clusters submenu):[/grey50]")
        lines.append(f"  Binary:         {cls_def.get('binary', '')}")
        lines.append(f"  Partition:      {cls_def.get('partition', '')}")
        if cls_def.get('mode') == 'gpu':
            n = cls_def.get('gpus', 1)
            gres_base = cls_def.get('gres', '')
            gres = gres_base.replace('<n>', str(n)).replace('{n}', str(n))
            lines.append(f"  GRES:           {gres}")
        else:
            ntasks = cls_def.get('ntasks', 1)
            lines.append(f"  Tasks:          {ntasks}")
        lines.append(f"  Memory:         {cls_def.get('memory', '')}")
        launcher = cls_def.get('mpi_launcher', '')
        if launcher:
            lines.append(f"  MPI launcher:   {launcher}")
        env_vars = cls_def.get('env_vars', [])
        if env_vars:
            lines.append(f"  env_vars:       {', '.join(env_vars)}")
        self.console.print()
        self.console.print(Panel("\n".join(lines), expand=False, border_style="cyan"))

    def _configure_bulk_engines(self, has_gpus=False):
        """Set all simulations to the same engine configuration."""
        active_queue = self._get_active_queue()
        if not active_queue:
            self.console.print("[yellow]No simulations in queue to configure[/yellow]")
            return

        self.console.print("\n[bold]Bulk Engine Configuration[/bold]")
        self.console.print("[grey50]Set the same engine for all queued simulations[/grey50]\n")
        
        # Get available engines
        available_engines = self._get_available_engines()
        
        self.console.print("Available engines:")
        for i, engine in enumerate(available_engines, 1):
            self.console.print(f"  {i}. {engine}")

        engine_choices = [str(i) for i in range(1, len(available_engines)+1)]

        # Build options map
        options_map = {}
        for i, engine in enumerate(available_engines, 1):
            options_map[str(i)] = engine
        options_map["cancel"] = "Cancel bulk configuration"

        choice = prompt_with_context(
            self.processor,
            "Select engine for all simulations",
            choices=engine_choices + ["cancel"],
            default="cancel",
            module="MD Manager - Engine Configuration",
            description="Select engine for all queued simulations",
            options_map=options_map
        )
        
        if choice == "cancel":
            return
            
        selected_engine = available_engines[int(choice)-1]
        
        # Configure hardware resources if needed
        hardware_config = {}
        if "cuda" in selected_engine.lower():
            hardware_config = self._configure_gpu_resources()
        elif "mpi" in selected_engine.lower():
            hardware_config = self._configure_mpi_resources()
            
        # Apply to active simulations
        for config in active_queue:
            config.engine = selected_engine
            config.hardware_config = hardware_config

        self.console.print(f"[green]✓ {len(active_queue)} simulations configured with {selected_engine}[/green]")

    def _get_simulation_mdin_params(self, config: SimulationConfig) -> Dict:
        """Get MDIN parameters for a simulation config, from override or template."""
        # Use override content first
        if config.mdin_content_override:
            return self._parse_mdin_params(config.mdin_content_override)

        # Try loading from template
        if hasattr(self, 'user_data_manager') and self.user_data_manager and config.template_id:
            try:
                content, _ = self.user_data_manager.get_template_content(config.template_id)
                if content:
                    params = self._parse_mdin_params(content)
                    # Apply parameter overrides on top
                    if config.parameter_overrides:
                        params.update(config.parameter_overrides)
                    return params
            except Exception:
                pass

        return config.parameter_overrides or {}

    def _classify_for_engine(self, config: SimulationConfig,
                             npt_index: int = 0, npt_total: int = 0) -> tuple:
        """Classify a simulation as CPU or GPU based on MDIN parameters.

        Args:
            config: The simulation configuration to classify.
            npt_index: Zero-based index of this NPT step among all NPT steps.
            npt_total: Total number of NPT steps in the queue.

        Returns:
            (target, reason) where target is 'cpu' or 'gpu' and reason is a display string.
        """
        params = self._get_simulation_mdin_params(config)

        try:
            imin = int(params.get('imin', 0))
        except (ValueError, TypeError):
            imin = 0
        try:
            ntp = int(params.get('ntp', 0))
        except (ValueError, TypeError):
            ntp = 0

        if imin == 1:
            return 'cpu', 'force overflow risk (SPFP)'
        elif ntp > 0:
            # Split NPT across CPU (first step) and GPU (later steps) whenever
            # there are ≥2 NPT steps in the queue. The first NPT step sees the
            # most box shrinkage as density converges, which pmemd.cuda can't
            # handle well (fixed PME grid at job start); subsequent NPT steps
            # run on GPU because the box is by then near-stationary.
            if npt_total > 1 and npt_index > 0:
                return 'gpu', 'box equilibrated, safe for GPU'
            return 'cpu', 'PME grid resizing (early NPT)'
        else:
            return 'gpu', 'safe for GPU'

    def _apply_recommended_engines(self, gpu_info: Dict):
        """Apply recommended engine assignments based on simulation type.

        Minimization steps get pmemd.MPI (CPU) to avoid SPFP force overflow.
        NPT steps: if there are more than 2, only the first gets CPU (box
        equilibrates quickly); otherwise all NPT steps use CPU.
        NVT and production steps get pmemd.cuda (GPU) for speed.
        """
        active_queue = self._get_active_queue()
        if not active_queue:
            self.console.print("[yellow]No simulations in queue to configure[/yellow]")
            return

        from rich.table import Table

        # Pre-pass: identify NPT steps so we can apply the >2 NPT rule
        npt_indices = []
        npt_counter = 0
        for config in active_queue:
            params = self._get_simulation_mdin_params(config)
            try:
                imin = int(params.get('imin', 0))
            except (ValueError, TypeError):
                imin = 0
            try:
                ntp = int(params.get('ntp', 0))
            except (ValueError, TypeError):
                ntp = 0
            if imin == 0 and ntp > 0:
                npt_indices.append(npt_counter)
                npt_counter += 1
            else:
                npt_indices.append(-1)  # not an NPT step

        npt_total = npt_counter

        # Classify all simulations
        classifications = []
        needs_cpu = False
        needs_gpu = False

        for i, config in enumerate(active_queue):
            npt_index = npt_indices[i]
            target, reason = self._classify_for_engine(
                config, npt_index=npt_index, npt_total=npt_total)
            classifications.append((config, target, reason))
            if target == 'cpu':
                needs_cpu = True
            else:
                needs_gpu = True

        # Show classification table
        # If a cluster profile is loaded, drive hardware from the palette
        # instead of prompting for MPI/GPU resources interactively.
        profile = self._cluster_profile
        if profile is not None:
            cpu_class = profile.defaults.get("cpu_class", "")
            gpu_class = profile.defaults.get("gpu_class", "")
            if (needs_cpu and not cpu_class) or (needs_gpu and not gpu_class):
                self.console.print(
                    "[yellow]Profile is missing a default resource class for "
                    "this step type; falling back to manual prompts[/yellow]"
                )
                profile = None

        # Build and show classification table. With a profile, show the
        # resolved resource class + wall time so the user can audit before
        # committing. Group rows by structure_label when the queue spans
        # multiple structures so the table stays readable at scale; global
        # indices are preserved so downstream selection semantics don't shift.
        def new_assignment_table(title=None):
            t = Table(title=title)
            t.add_column("#", style="grey50")
            t.add_column("Step", style="cyan")
            t.add_column("Reason", style="grey50")
            t.add_column("Engine", style="yellow")
            if profile is not None:
                t.add_column("Resource class", style="bright_blue")
                t.add_column("Wall time", style="bright_blue")
            return t

        def row_for(i, config, target, reason):
            display_name = config.step_name or config.name
            engine = "pmemd.MPI" if target == 'cpu' else "pmemd.cuda"
            row = [str(i), display_name, reason, engine]
            if profile is not None:
                class_name = profile.defaults['cpu_class'] if target == 'cpu' else profile.defaults['gpu_class']
                try:
                    wall = profile.resolve_class(class_name).get('default_time', '')
                except KeyError:
                    wall = ''
                row.extend([class_name, wall or '[grey50](unset)[/grey50]'])
            return row

        indexed = list(enumerate(classifications, 1))
        groups: "dict[str, list]" = {}
        for idx, (config, target, reason) in indexed:
            label = getattr(config, 'structure_label', None) or ""
            groups.setdefault(label, []).append((idx, config, target, reason))

        multi_group = len(groups) > 1 and any(label for label in groups)

        if multi_group:
            self.console.print("[bold]Recommended Assignments:[/bold]")
            for label, entries in groups.items():
                display_label = label or "(unlabeled)"
                self.console.print(f"\n  [cyan]{display_label}[/cyan]")
                table = new_assignment_table()
                for idx, config, target, reason in entries:
                    table.add_row(*row_for(idx, config, target, reason))
                self.console.print(table)
        else:
            table = new_assignment_table(title="Recommended Assignments")
            for idx, (config, target, reason) in indexed:
                table.add_row(*row_for(idx, config, target, reason))
            self.console.print(table)

        if profile is not None:
            # Confirm before committing. 'e' hands off to per-step editor.
            self.console.print(
                "\n[grey50]y = apply as shown · e = edit per step · n = cancel[/grey50]"
            )
            choice = prompt_with_context(
                self.processor,
                "Apply these assignments?",
                choices=["y", "n", "e"],
                default="y",
                module="MD Manager - Apply Recommended",
                description="Confirm recommended resource assignments",
                options_map={"y": "Apply as shown", "n": "Cancel", "e": "Edit per step"},
            )
            if choice == "n":
                self.console.print("[yellow]Cancelled — no assignments applied[/yellow]")
                return
            for config, target, reason in classifications:
                class_name = profile.defaults['cpu_class'] if target == 'cpu' else profile.defaults['gpu_class']
                cls_def = profile.resolve_class(class_name)
                config.engine = cls_def.get('binary', 'pmemd.MPI' if target == 'cpu' else 'pmemd.cuda')
                config.hardware_config = {
                    'resource_class': class_name,
                    'cluster_profile': profile.name,
                    'time_limit': cls_def.get('default_time', ''),
                }
            self.console.print(
                f"\n[green]✓ {len(classifications)} simulations configured "
                f"from cluster profile '{profile.name}'[/green]"
            )
            cpu_count = sum(1 for _, t, _ in classifications if t == 'cpu')
            gpu_count = sum(1 for _, t, _ in classifications if t == 'gpu')
            if cpu_count:
                self.console.print(f"  {cpu_count}x {profile.defaults['cpu_class']}")
            if gpu_count:
                self.console.print(f"  {gpu_count}x {profile.defaults['gpu_class']}")
            if choice == "e":
                self._configure_individual_engines(gpu_info['available'] > 0)
            return

        # No profile: original interactive path.
        mpi_config = {}
        if needs_cpu:
            mpi_config = self._configure_mpi_resources()

        gpu_config = {}
        if needs_gpu:
            gpu_config = self._configure_gpu_resources()

        for config, target, reason in classifications:
            if target == 'cpu':
                config.engine = "pmemd.MPI"
                config.hardware_config = mpi_config
            else:
                config.engine = "pmemd.cuda"
                config.hardware_config = gpu_config

        self.console.print(f"\n[green]✓ {len(classifications)} simulations configured[/green]")

        cpu_count = sum(1 for _, t, _ in classifications if t == 'cpu')
        gpu_count = sum(1 for _, t, _ in classifications if t == 'gpu')
        if cpu_count:
            self.console.print(f"  {cpu_count}x pmemd.MPI ({mpi_config.get('mpi_tasks', '?')} processes)")
        if gpu_count:
            self.console.print(f"  {gpu_count}x pmemd.cuda (GPU {gpu_config.get('gpu_ids', '0')})")

    # --- Cluster profile / run plan action handlers (Step 4) -------------

    def _action_load_cluster_profile(self):
        """Interactive profile loader: list profiles, pick one, fill required fields."""
        from proprep.md_prep import cluster_profile as cp

        if self._cluster_profile is not None:
            unload = confirm_with_context(
                self.processor,
                f"A profile is already loaded ('{self._cluster_profile.name}'). "
                "Unload and pick a new one?",
                default=False,
                module="MD Manager - Cluster Profiles",
                description="Unload current profile and pick a new one",
            )
            if not unload:
                return
            self._cluster_profile = None
            self._run_plan = None

        entries = [e for e in cp.list_profiles() if not e['name'].startswith('_')]
        if not entries:
            self.console.print(
                "[yellow]No cluster profiles available. Use option 6 "
                "'Manage cluster profiles' in the MD Manager menu to create "
                "one from the template.[/yellow]"
            )
            return

        # Deduplicate (a name can appear in multiple scopes — we merge at load)
        seen = {}
        for e in entries:
            seen.setdefault(e['name'], e)
        entries = list(seen.values())

        self.console.print("\n[bold]Available cluster profiles:[/bold]")
        for i, e in enumerate(entries, 1):
            self.console.print(f"  {i}. {e['display_name']} [grey50]({e['name']}, {e['source']})[/grey50]")
        self.console.print(
            "[grey50]  Enter a number to preview, or prefix with 'v' (e.g. 'v2') "
            "to view JSON. 'c' to cancel.[/grey50]"
        )

        while True:
            raw = prompt_with_context(
                self.processor,
                "Enter number (or 'c' to cancel)",
                default="1",
                module="MD Manager - Cluster Profiles",
                description="Select cluster profile to load (or 'c' to cancel, 'v#' to view JSON)",
            ).strip().lower()
            if raw == 'c':
                return
            view_raw = False
            if raw.startswith('v'):
                view_raw = True
                raw = raw[1:].strip()
            try:
                idx = int(raw) - 1
                if idx < 0 or idx >= len(entries):
                    raise ValueError
            except ValueError:
                self.console.print("[red]Invalid selection[/red]")
                continue

            name = entries[idx]['name']
            profile = cp.load_profile(name)

            if view_raw:
                self.console.print(
                    f"\n[bold cyan]--- {name}.json ({profile.source}) ---[/bold cyan]"
                )
                self.console.print(json.dumps(profile.to_json(), indent=2))
                continue

            # Human-readable preview
            self._display_cluster_profile_summary(profile)
            use = confirm_with_context(
                self.processor,
                f"Use '{profile.name}' as the active cluster profile?",
                default=True,
                module="MD Manager - Cluster Profiles",
                description=f"Use profile '{profile.name}' as active",
            )
            if use:
                break
            else:
                continue

        # Prompt for missing required fields, then save filled copy to user scope
        missing = profile.missing_required_fields()
        if missing:
            self.console.print(
                f"\n[yellow]Profile '{name}' needs site-specific values:[/yellow]"
            )
            for path in missing:
                value = prompt_with_context(
                    self.processor,
                    f"  {path}",
                    default="",
                    module="MD Manager - Cluster Profiles",
                    description=f"Fill cluster profile field '{path}'",
                )
                # write into cluster dict via dotted path
                keys = path.split(".")
                target = profile.cluster
                # Simple case: single-level field inside cluster dict
                if len(keys) == 1:
                    target[keys[0]] = value
                else:
                    for k in keys[:-1]:
                        target = target.setdefault(k, {})
                    target[keys[-1]] = value
            if confirm_with_context(
                self.processor,
                "Save filled-in profile to ~/.proprep/cluster_profiles/ for next time?",
                default=True,
                module="MD Manager - Cluster Profiles",
                description="Save filled profile to user scope",
            ):
                cp.save_profile(profile, scope="user")
                self.console.print(f"[green]✓ Saved to user scope[/green]")

        errors = cp.validate(profile)
        if errors:
            self.console.print("[red]Profile has validation errors:[/red]")
            for e in errors:
                self.console.print(f"  • {e}")
            return

        self._cluster_profile = profile
        self.console.print(
            f"\n[green]✓ Loaded cluster profile '{profile.name}' "
            f"({len(profile.resource_classes)} resource classes)[/green]"
        )

    def _display_cluster_profile_summary(self, profile):
        """Print a readable preview of a cluster profile's contents."""
        from rich.table import Table
        cluster = profile.cluster or {}

        self.console.print()
        self.console.print(Panel(
            f"[bold]{profile.display_name or profile.name}[/bold] "
            f"[grey50]({profile.name}, {profile.source})[/grey50]\n"
            f"{profile.description or '[grey50](no description)[/grey50]'}",
            title="Profile",
            border_style="cyan",
            expand=False,
            padding=(0, 1),
        ))
        self.console.print(f"  [bold]Scheduler:[/bold] {cluster.get('scheduler', 'slurm')}")
        account = cluster.get('account', '')
        self.console.print(f"  [bold]Account:[/bold]   {account or '[grey50](not set)[/grey50]'}")
        qos = cluster.get('qos', '')
        if qos:
            self.console.print(f"  [bold]QoS:[/bold]       {qos}")
        modules = cluster.get('modules', [])
        if modules:
            self.console.print(f"  [bold]Modules:[/bold]   {', '.join(modules)}")
        conda_env = cluster.get('conda_env', '')
        if conda_env:
            self.console.print(f"  [bold]Conda env:[/bold] {conda_env}")
        pre = cluster.get('pre_commands', [])
        if pre:
            self.console.print(f"  [bold]Pre-commands:[/bold]")
            for cmd in pre:
                self.console.print(f"    {cmd}")
        post = cluster.get('post_commands', [])
        if post:
            self.console.print(f"  [bold]Post-commands:[/bold]")
            for cmd in post:
                self.console.print(f"    {cmd}")
        extra = cluster.get('extra_directives', [])
        if extra:
            self.console.print(f"  [bold]Extra #SBATCH:[/bold]")
            for d in extra:
                self.console.print(f"    {d}")

        # Resource class table
        if profile.resource_classes:
            table = Table(show_header=True, header_style="bold cyan", title="Resource classes")
            table.add_column("Name")
            table.add_column("Mode")
            table.add_column("Partition")
            table.add_column("GRES / CPUs")
            table.add_column("Memory")
            table.add_column("Binary")
            table.add_column("Default wall")
            for name, cdef in profile.resource_classes.items():
                mode = cdef.get('mode', '')
                part = cdef.get('partition', '')
                if mode == 'gpu':
                    gres_cpu = cdef.get('gres', '')
                else:
                    gres_cpu = f"{cdef.get('ntasks', '?')} tasks × {cdef.get('nodes', 1)} node(s)"
                table.add_row(
                    name, mode, part, gres_cpu,
                    str(cdef.get('memory', '')),
                    cdef.get('binary', ''),
                    cdef.get('default_time', ''),
                )
            self.console.print(table)

        defs = profile.defaults or {}
        if defs:
            self.console.print(
                f"  [bold]Defaults:[/bold] cpu → {defs.get('cpu_class', '[grey50]none[/grey50]')}, "
                f"gpu → {defs.get('gpu_class', '[grey50]none[/grey50]')}"
            )

        from proprep.md_prep import cluster_profile as cp
        errors = cp.validate(profile)
        if errors:
            self.console.print("[yellow]Validation warnings:[/yellow]")
            for e in errors:
                self.console.print(f"  • {e}")

    def _action_load_run_plan(self):
        """Interactive run-plan loader, filtered to plans compatible with loaded profile."""
        from proprep.md_prep import run_plan as rp

        all_plans = rp.list_plans()
        if not all_plans:
            self.console.print("[yellow]No saved run plans available[/yellow]")
            return

        # Show all plans, mark those bound to a different cluster
        active_cluster = self._cluster_profile.name
        self.console.print("\n[bold]Available run plans:[/bold]")
        for i, e in enumerate(all_plans, 1):
            marker = "" if e['cluster_name'] == active_cluster else " [grey50 italic](different cluster — will need remap)[/grey50 italic]"
            self.console.print(
                f"  {i}. {e['name']} [grey50]({e['protocol_name']}, cluster={e['cluster_name']})[/grey50]{marker}"
            )

        raw = prompt_with_context(
            self.processor,
            "Enter number (or 'c' to cancel)",
            default="1",
            module="MD Manager - Run Plans",
            description="Select run plan to load",
        )
        if raw.strip().lower() == 'c':
            return
        try:
            idx = int(raw) - 1
            if idx < 0 or idx >= len(all_plans):
                raise ValueError
        except ValueError:
            self.console.print("[red]Invalid selection[/red]")
            return

        plan = rp.load_plan(all_plans[idx]['name'])
        if plan.cluster_name != active_cluster:
            self.console.print(
                f"[yellow]Plan was built for '{plan.cluster_name}', "
                f"but '{active_cluster}' is loaded. Remap not yet implemented "
                f"in this flow — load via the 'plans' submenu to remap.[/yellow]"
            )
            return

        self._run_plan = plan
        self.console.print(f"\n[green]✓ Loaded run plan '{plan.name}'[/green]")

        # Apply plan's assignments to the queue
        queue = self._get_active_queue() or self.simulation_queue.queue
        applied = 0
        for sim in queue:
            step_id = getattr(sim, 'step_id', None) or getattr(sim, 'name', '')
            if step_id in plan.step_resources:
                sr = plan.step_resources[step_id]
                cls = self._cluster_profile.resolve_class(sr.class_name)
                sim.engine = cls.get('binary', sim.engine)
                sim.hardware_config = {
                    'resource_class': sr.class_name,
                    'cluster_profile': self._cluster_profile.name,
                    'time_limit': sr.time_limit,
                    'overrides': dict(sr.overrides),
                }
                applied += 1
        self.console.print(f"  Applied assignments to {applied} simulations")

    def _action_save_run_plan(self):
        """Persist current queue assignments as a run plan in user scope."""
        from proprep.md_prep import run_plan as rp
        from proprep.md_prep.run_plan import RunPlan, StepResource

        queue = self._get_active_queue() or self.simulation_queue.queue
        step_resources = {}
        for sim in queue:
            hw = sim.hardware_config or {}
            cls_name = hw.get('resource_class')
            if not cls_name:
                continue  # skip unassigned or non-profile entries
            step_id = getattr(sim, 'step_id', None) or sim.name
            step_resources[step_id] = StepResource(
                class_name=cls_name,
                time_limit=hw.get('time_limit', ''),
                overrides=dict(hw.get('overrides', {})),
            )
        if not step_resources:
            self.console.print(
                "[yellow]No profile-based assignments to save. "
                "Use option 1 or 2 with a loaded profile first.[/yellow]"
            )
            return

        default_name = f"{self._cluster_profile.name}-plan"
        name = prompt_with_context(
            self.processor,
            "Plan name",
            default=default_name,
            module="MD Manager - Run Plans",
            description="Run plan name to save",
        )
        # Figure out protocol name from active workflow if available
        protocol_name = ""
        protocol_version = ""
        workspace = self.workspace
        if workspace:
            wf = workspace.get('active_workflow') or {}
            protocol_name = wf.get('workflow_name', '')
            protocol_version = wf.get('version', '')
        plan = RunPlan(
            name=name,
            protocol_name=protocol_name,
            protocol_version=protocol_version,
            cluster_name=self._cluster_profile.name,
            step_resources=step_resources,
        )
        path = rp.save_plan(plan, scope="user")
        self._run_plan = plan
        self.console.print(f"[green]✓ Run plan saved: {path}[/green]")

        # Sync to workspace
        if hasattr(self.simulation_queue, '_sync_to_workspace'):
            self.simulation_queue._sync_to_workspace()

    def _configure_slurm_mode(self) -> bool:
        """
        Configure SLURM job script generation for HPC cluster submission.

        Returns:
            True if scripts were generated successfully, False otherwise
        """
        from .slurm_generator import SlurmJobGenerator

        self.console.print("\n[bold cyan]===== SLURM Job Script Generation =====[/bold cyan]\n")
        self.console.print("[grey50]Generate job scripts for HPC cluster submission via SLURM[/grey50]\n")

        # Initialize generator
        generator = SlurmJobGenerator(console=self.console)

        # Interactive configuration
        slurm_config = generator.interactive_configure()

        if not slurm_config:
            self.console.print("[yellow]SLURM configuration cancelled[/yellow]")
            return False

        # Determine output directory
        output_dir = Path.cwd() / "simulations" / "slurm_scripts"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Check if this is a workflow (has dependencies)
        has_workflow = any(
            hasattr(config, 'workflow_id') and config.workflow_id
            for config in self.simulation_queue.queue
        )

        if has_workflow:
            # Generate workflow scripts with dependencies
            self._generate_workflow_slurm_scripts(generator, output_dir)
        else:
            # Generate independent scripts
            self._generate_independent_slurm_scripts(generator, output_dir)

        # Save SLURM configuration
        config_path = output_dir / "slurm_config.json"
        generator.save_config(config_path)
        self.console.print(f"\n[grey50]SLURM configuration saved to: {config_path}[/grey50]")

        # Mark as SLURM mode in workspace
        if self.workspace:
            self.workspace.set('slurm_mode', True)
            self.workspace.set('slurm_output_dir', str(output_dir))

        return True

    def _generate_independent_slurm_scripts(self, generator: 'SlurmJobGenerator', output_dir: Path):
        """Generate independent SLURM scripts (no workflow dependencies)."""
        self.console.print("\n[bold]Generating SLURM job scripts...[/bold]\n")

        generated_scripts = []

        for i, sim_config in enumerate(self.simulation_queue.queue, 1):
            # Use stepN for directory + file naming so SLURM layout matches the
            # batch path. The semantic step name is preserved in progress output.
            step_key = f"step{i}"
            sim_dir = output_dir / step_key
            sim_dir.mkdir(exist_ok=True)

            # Copy topology and coordinate files to sim directory
            import shutil
            prmtop_path = Path(sim_config.prmtop)
            coord_path = Path(sim_config.rst7)

            shutil.copy2(prmtop_path, sim_dir / prmtop_path.name)
            shutil.copy2(coord_path, sim_dir / coord_path.name)

            # Copy MDIN file if exists
            if sim_config.mdin_path and Path(sim_config.mdin_path).exists():
                shutil.copy2(sim_config.mdin_path, sim_dir / "simulation.mdin")
            else:
                # Generate MDIN from template, applying any configured
                # restraints (mirrors what the batch path does).
                template_content = self._resolve_mdin_content(sim_config)
                if template_content:
                    template_content = self._apply_configured_restraints(
                        template_content, sim_config, sim_dir
                    )
                    with open(sim_dir / "simulation.mdin", 'w') as f:
                        f.write(template_content)

            # Update generator job name for this simulation
            original_job_name = generator.config.job_name
            generator.config.job_name = step_key

            # Save location for the sbatch script (inside the sim_dir for
            # independent runs). Passed to generate_script so the emitted
            # `cd` is script-relative rather than an absolute baked path.
            script_path = sim_dir / f"{step_key}.sh"
            script_content = generator.generate_script(
                sim_config=sim_config,
                sim_dir=sim_dir,
                engine=sim_config.engine or 'sander',
                topology_file=sim_dir / prmtop_path.name,
                coordinate_file=sim_dir / coord_path.name,
                output_prefix_override=step_key,
                script_path=script_path,
            )

            # Restore original job name
            generator.config.job_name = original_job_name
            with open(script_path, 'w') as f:
                f.write(script_content)

            script_path.chmod(0o755)
            generated_scripts.append((sim_config.name, script_path))

            self.console.print(f"[green]✓[/green] {sim_config.name}: {script_path.relative_to(output_dir)}")

        # Display submission instructions
        self._display_slurm_submission_instructions(generated_scripts, output_dir, is_workflow=False)

    def _write_slurm_scripts_from_profile(self, submit: bool = False) -> bool:
        """Generate SLURM scripts from the loaded cluster profile + run plan.

        When no run plan is loaded, builds an in-memory plan from each sim's
        hardware_config (which _apply_recommended_engines / _action_load_run_plan
        populated). If submit=True, runs sbatch on submit_workflow.sh.
        """
        from proprep.md_prep.slurm_generator import generate_scripts_from_plan
        from proprep.md_prep.run_plan import RunPlan, StepResource

        profile = self._cluster_profile
        if profile is None:
            self.console.print("[red]No cluster profile loaded[/red]")
            return False

        # Build or reuse the run plan
        plan = self._run_plan
        if plan is None:
            step_resources = {}
            for sim in self.simulation_queue.queue:
                hw = sim.hardware_config or {}
                cls_name = hw.get('resource_class')
                if not cls_name:
                    self.console.print(
                        f"[yellow]Simulation '{sim.name}' has no resource class "
                        f"assigned — go back to Step 4 option 1 or 2[/yellow]"
                    )
                    return False
                step_id = getattr(sim, 'step_id', None) or sim.name
                step_resources[step_id] = StepResource(
                    class_name=cls_name,
                    time_limit=hw.get('time_limit', ''),
                    overrides=dict(hw.get('overrides', {})),
                )
            plan = RunPlan(
                name="ad-hoc",
                protocol_name="",
                cluster_name=profile.name,
                step_resources=step_resources,
            )

        output_dir = Path.cwd() / "simulations" / "slurm_scripts"
        output_dir.mkdir(parents=True, exist_ok=True)

        workflow_sims = self._build_workflow_sims_for_slurm(output_dir)

        try:
            script_paths = generate_scripts_from_plan(
                cluster_profile=profile,
                run_plan=plan,
                workflow_sims=workflow_sims,
                output_dir=output_dir,
                console=self.console,
            )
        except Exception as e:
            self.console.print(f"[red]SLURM generation failed: {e}[/red]")
            return False

        self.console.print(f"\n[bold green]✓ Wrote {len(script_paths)} SLURM scripts to {output_dir}[/bold green]")
        for name, path in script_paths.items():
            self.console.print(f"  [green]✓[/green] {name}: {path.name}")

        submit_script = output_dir / "submit_workflow.sh"
        self.console.print(f"  [green]✓[/green] Master submission: {submit_script.name}")

        # Preview / edit loop. Before anything hits sbatch the user gets a
        # chance to inspect and tweak each generated script.
        proceed = self._preview_generated_scripts(
            script_paths=script_paths,
            submit_script=submit_script,
            require_confirm_to_submit=submit,
        )
        if submit and not proceed:
            self.console.print("[yellow]Submission cancelled — scripts remain on disk for manual review.[/yellow]")
            submit = False

        if submit:
            self.console.print("\n[bold]Submitting via sbatch...[/bold]")
            try:
                result = subprocess.run(
                    ["bash", str(submit_script)],
                    cwd=output_dir,
                    capture_output=True,
                    text=True,
                )
                self.console.print(result.stdout)
                if result.returncode != 0:
                    self.console.print(f"[red]sbatch submission failed (exit {result.returncode})[/red]")
                    self.console.print(result.stderr)
                    return False
            except Exception as e:
                self.console.print(f"[red]Error running submit script: {e}[/red]")
                return False
            self.console.print("[green]✓ Jobs submitted[/green]")

        if self.workspace is not None:
            self.workspace.set('slurm_mode', True)
            self.workspace.set('slurm_output_dir', str(output_dir))
        return True

    def _preview_generated_scripts(
        self,
        script_paths: Dict[str, Path],
        submit_script: Path,
        require_confirm_to_submit: bool,
    ) -> bool:
        """Review loop for generated SLURM scripts (files already on disk).

        The user can view each script, edit in $EDITOR (changes persist to
        disk), accept, or cancel (when require_confirm_to_submit=True,
        cancel aborts sbatch but leaves files in place).
        """
        ordered = list(script_paths.items())
        ordered.append(("__submit_workflow__", submit_script))

        while True:
            self.console.print("\n[bold cyan]Script Review[/bold cyan]")
            self.console.print(
                "[grey50]Files have been written. Edits you make here are saved to the files on disk.[/grey50]"
            )
            for i, (name, path) in enumerate(ordered, 1):
                label = "master submission" if name == "__submit_workflow__" else name
                self.console.print(f"  {i}. {label} [grey50]{path.name}[/grey50]")
            self.console.print(
                "\n[grey50]v <N> = view contents · e <N> = edit in $EDITOR · "
                "d = done · c = cancel[/grey50]"
            )
            raw = prompt_with_context(
                self.processor,
                "Action",
                default="d",
                module="MD Manager - SLURM",
                description="Script review action (v N=view, e N=edit, d=done, c=cancel)",
            ).strip()
            if not raw:
                continue
            low = raw.lower()
            if low in ("d", "done", ""):
                return True
            if low in ("c", "cancel"):
                return False

            # Parse 'v N' or 'e N' (accept legacy 'p N' for view too)
            parts = low.split()
            if len(parts) != 2 or parts[0] not in ("v", "p", "e"):
                self.console.print("[yellow]Unrecognized input. Use 'v 1', 'e 2', 'd', or 'c'.[/yellow]")
                continue
            try:
                idx = int(parts[1]) - 1
                if idx < 0 or idx >= len(ordered):
                    raise ValueError
            except ValueError:
                self.console.print("[red]Invalid script number[/red]")
                continue

            name, path = ordered[idx]
            if parts[0] in ("v", "p"):
                self.console.print()
                self.console.print(Panel(
                    path.read_text(),
                    title=path.name,
                    border_style="cyan",
                    expand=False,
                ))
            else:  # edit
                editor = os.environ.get("EDITOR") or os.environ.get("VISUAL") or "vi"
                try:
                    subprocess.run([editor, str(path)])
                    self.console.print(f"[grey50]Saved edits to {path.name}[/grey50]")
                except FileNotFoundError:
                    self.console.print(
                        f"[red]Editor '{editor}' not found. Edit {path} directly.[/red]"
                    )
                except Exception as e:
                    self.console.print(f"[red]Editor error: {e}[/red]")

    def _build_workflow_sims_for_slurm(self, output_dir: Path) -> List[Dict[str, Any]]:
        """Stage files into per-sim directories; return the workflow_sims list.

        Shared prmtop + initial rst7 are staged once at output_dir; each
        step directory holds only its own mdin (and CpHMD files if any).
        The -c/-p/-ref arguments for each step reference the shared files
        at the parent level, and step N > 1 reads step N-1's output rst7.
        """
        queue = list(self.simulation_queue.queue)
        workflow_sims: List[Dict[str, Any]] = []
        if not queue:
            return workflow_sims

        first = queue[0]
        prmtop_path = Path(first.prmtop)
        coord_path = Path(first.rst7)
        if prmtop_path.exists():
            shutil.copy2(prmtop_path, output_dir / prmtop_path.name)
        if coord_path.exists():
            shutil.copy2(coord_path, output_dir / coord_path.name)

        prev_step_key: Optional[str] = None
        for i, sim_config in enumerate(queue, 1):
            # stepN for dir + file prefix so SLURM layout mirrors batch.
            step_key = f"step{i}"
            sim_dir = output_dir / step_key
            sim_dir.mkdir(exist_ok=True)

            if sim_config.mdin_path and Path(sim_config.mdin_path).exists():
                shutil.copy2(sim_config.mdin_path, sim_dir / "simulation.mdin")
            else:
                template_content = self._resolve_mdin_content(sim_config)
                if template_content:
                    template_content = self._apply_configured_restraints(
                        template_content, sim_config, sim_dir
                    )
                    with open(sim_dir / "simulation.mdin", 'w') as f:
                        f.write(template_content)

            topology_arg = f"../{prmtop_path.name}"
            if prev_step_key is None:
                coord_arg = f"../{coord_path.name}"
            else:
                coord_arg = f"../{prev_step_key}/{prev_step_key}.rst7"
            reference_arg = coord_arg

            cpmd_flags = None
            if self._is_production_step(sim_config):
                workflow = self._get_workflow_for_step(sim_config)
                if workflow and workflow.cpin_file:
                    cpin_path = Path(workflow.cpin_file)
                    if cpin_path.exists():
                        shutil.copy2(cpin_path, sim_dir / cpin_path.name)
                        cpmd_flags = [
                            "-cpin", cpin_path.name,
                            "-cpout", f"{step_key}.cpout",
                            "-cprestrt", f"{step_key}.cprestrt",
                        ]
                    if workflow.cpin_config:
                        mod_prmtop = workflow.cpin_config.get('modified_prmtop')
                        if mod_prmtop and Path(mod_prmtop).exists():
                            shutil.copy2(mod_prmtop, sim_dir / Path(mod_prmtop).name)

            workflow_sims.append({
                'config': sim_config,
                'sim_dir': sim_dir,
                'extra_flags': cpmd_flags,
                'topology_arg': topology_arg,
                'coord_arg': coord_arg,
                'reference_arg': reference_arg,
                'output_prefix': step_key,
                'step_key': step_key,
            })
            prev_step_key = step_key

        return workflow_sims

    def _generate_workflow_slurm_scripts(self, generator: 'SlurmJobGenerator', output_dir: Path):
        """Generate SLURM scripts with workflow dependencies.

        Layout:
            output_dir/
              prmtop                     (one copy, shared)
              initial.rst7               (one copy, shared)
              <step_slug>/
                simulation.mdin          (per-step inputs)
              <step_slug>.sh             (sbatch wrapper, -c refs prior step)
              submit_workflow.sh

        Each step reads its predecessor's output rst7 via a relative path
        from the step's working directory, so step N actually picks up
        where step N-1 left off instead of restarting from the initial
        coordinates every time.
        """
        self.console.print("\n[bold]Generating protocol SLURM job scripts with dependencies...[/bold]\n")

        import shutil
        queue = list(self.simulation_queue.queue)
        if not queue:
            self.console.print("[yellow]No simulations queued[/yellow]")
            return

        # Stage prmtop + initial rst7 once at output_dir. All sims in the
        # workflow share the same topology/initial coords.
        first = queue[0]
        prmtop_path = Path(first.prmtop)
        coord_path = Path(first.rst7)
        if prmtop_path.exists():
            shutil.copy2(prmtop_path, output_dir / prmtop_path.name)
        if coord_path.exists():
            shutil.copy2(coord_path, output_dir / coord_path.name)

        # Prepare per-step directories and the workflow_sims entries.
        workflow_sims = []
        prev_step_key: Optional[str] = None

        for i, sim_config in enumerate(queue, 1):
            # stepN for directory + all output artifact names so the SLURM
            # layout matches the batch path.
            step_key = f"step{i}"
            sim_dir = output_dir / step_key
            sim_dir.mkdir(exist_ok=True)

            # MDIN — still per-step (it's the point of per-step dirs).
            # Apply any user-configured restraints (restraintmask / GROUP /
            # DISANG) so the SLURM mdin matches what the batch path writes.
            if sim_config.mdin_path and Path(sim_config.mdin_path).exists():
                shutil.copy2(sim_config.mdin_path, sim_dir / "simulation.mdin")
            else:
                template_content = self._resolve_mdin_content(sim_config)
                if template_content:
                    template_content = self._apply_configured_restraints(
                        template_content, sim_config, sim_dir
                    )
                    with open(sim_dir / "simulation.mdin", 'w') as f:
                        f.write(template_content)

            # Paths that go into the amber command line, relative to sim_dir.
            topology_arg = f"../{prmtop_path.name}"
            if prev_step_key is None:
                coord_arg = f"../{coord_path.name}"
            else:
                coord_arg = f"../{prev_step_key}/{prev_step_key}.rst7"
            # -ref tracks whatever -c is, matching the batch run_workflow.sh
            # behaviour; fix that separately if users want restraints
            # anchored to the original structure throughout.
            reference_arg = coord_arg

            # CpHMD flags + modified topology (still staged into sim_dir,
            # since those are per-production-step).
            cpmd_flags = None
            if self._is_production_step(sim_config):
                workflow = self._get_workflow_for_step(sim_config)
                if workflow and workflow.cpin_file:
                    cpin_path = Path(workflow.cpin_file)
                    if cpin_path.exists():
                        shutil.copy2(cpin_path, sim_dir / cpin_path.name)
                        cpmd_flags = [
                            "-cpin", cpin_path.name,
                            "-cpout", f"{step_key}.cpout",
                            "-cprestrt", f"{step_key}.cprestrt",
                        ]
                    if workflow.cpin_config:
                        mod_prmtop = workflow.cpin_config.get('modified_prmtop')
                        if mod_prmtop and Path(mod_prmtop).exists():
                            shutil.copy2(mod_prmtop, sim_dir / Path(mod_prmtop).name)

            workflow_sims.append({
                'config': sim_config,
                'sim_dir': sim_dir,
                'extra_flags': cpmd_flags,
                'topology_arg': topology_arg,
                'coord_arg': coord_arg,
                'reference_arg': reference_arg,
                'output_prefix': step_key,
                'step_key': step_key,
            })
            prev_step_key = step_key

        # Generate scripts with dependency chaining
        script_paths = generator.generate_workflow_scripts(workflow_sims, output_dir)

        # Display generated scripts
        for sim_name, script_path in script_paths.items():
            self.console.print(f"[green]✓[/green] {sim_name}: {script_path.relative_to(output_dir)}")

        # Generate master submission script
        submit_script = generator.generate_submission_script(workflow_sims, script_paths, output_dir)
        self.console.print(f"\n[bold green]✓[/bold green] Master submission script: {submit_script.relative_to(output_dir)}")

        # Display workflow submission instructions
        self._display_workflow_slurm_instructions(submit_script, script_paths, output_dir)

    def _display_slurm_submission_instructions(self, scripts: List, output_dir: Path, is_workflow: bool):
        """Display instructions for submitting SLURM jobs."""
        self.console.print(f"\n[bold cyan]═══ SLURM Job Submission Instructions ═══[/bold cyan]\n")

        self.console.print(f"[bold]Scripts Location:[/bold] {output_dir}\n")

        self.console.print("[bold]To submit jobs:[/bold]")
        self.console.print("1. Transfer scripts to your HPC cluster", highlight=False)
        self.console.print("2. Navigate to the scripts directory", highlight=False)
        self.console.print("3. Submit each job with:\n", highlight=False)

        for sim_name, script_path in scripts[:3]:  # Show first 3 as examples
            relative_path = script_path.relative_to(output_dir)
            self.console.print(f"   [cyan]sbatch {relative_path}[/cyan]")

        if len(scripts) > 3:
            self.console.print(f"   [grey50]... and {len(scripts) - 3} more[/grey50]")

        self.console.print("\n[bold]To check job status:[/bold]")
        self.console.print("   [cyan]squeue -u $USER[/cyan]")

        self.console.print("\n[bold]To cancel a job:[/bold]")
        self.console.print("   [cyan]scancel <job_id>[/cyan]")

    def _display_workflow_slurm_instructions(self, submit_script: Path, script_paths: Dict, output_dir: Path):
        """Display instructions for submitting workflow with dependencies."""
        self.console.print(f"\n[bold cyan]═══ Protocol Submission Instructions ═══[/bold cyan]\n")

        self.console.print(f"[bold]Scripts Location:[/bold] {output_dir}\n")

        self.console.print(
            "[grey50]Layout: prmtop + initial rst7 are at the top level; each "
            "step has its own subdirectory with simulation.mdin. Step N>1 "
            "reads the prior step's output rst7 via a relative path "
            "(../<prev_step>/<prev>.rst7).[/grey50]\n"
        )
        self.console.print(
            "[grey50]Note: batch_*/ directories (if present) are a separate "
            "local-run layout written during protocol setup; slurm_scripts/ "
            "is the HPC layout. Pick one.[/grey50]\n"
        )

        self.console.print("[bold]Option 1: Use master submission script (Recommended)[/bold]")
        self.console.print("This script submits all jobs with proper dependencies:\n")
        self.console.print(f"   [cyan]./{submit_script.name}[/cyan]\n")

        self.console.print("[bold]Option 2: Submit individually[/bold]")
        self.console.print("Submit jobs one at a time, noting each job ID for dependencies:\n")

        sim_names = list(script_paths.keys())
        for i, (sim_name, script_path) in enumerate(script_paths.items()):
            relative_path = script_path.relative_to(output_dir)
            if i == 0:
                self.console.print(f"   [cyan]JOBID1=$(sbatch --parsable {relative_path})[/cyan]")
            else:
                prev_idx = i
                self.console.print(f"   [cyan]JOBID{i+1}=$(sbatch --parsable --dependency=afterok:$JOBID{i} {relative_path})[/cyan]")

        self.console.print("\n[bold]To monitor protocol:[/bold]")
        self.console.print("   [cyan]squeue -u $USER[/cyan]  (shows all your jobs)")
        self.console.print("   [cyan]watch squeue -u $USER[/cyan]  (auto-refresh every 2s)")

    def _show_hardware_suggestions(self):
        """Display hardware suggestions for different simulation types."""
        self.console.print("\n[bold cyan]===== Hardware Suggestions =====[/bold cyan]")
        
        # Engine overview
        self.console.print("\n[bold]Available AMBER Engines:[/bold]")
        self.console.print("  • [bold]pmemd.cuda[/bold] - GPU-accelerated, fastest for most simulations")
        self.console.print("  • [bold]pmemd.MPI[/bold] - Parallel CPU version for multi-core systems")
        self.console.print("  • [bold]pmemd[/bold] - Serial optimized CPU version")
        self.console.print("  • [bold]sander[/bold] - Standard CPU version with full feature support")
        
        # Best practices
        self.console.print("\n[bold]Recommended Engine Selection:[/bold]")
        self.console.print("\n  [cyan]For GPU systems:[/cyan]")
        self.console.print("    • Use [green]pmemd.cuda[/green] for heating and production simulations")
        self.console.print("    • Use [yellow]CPU engines[/yellow] for minimization (numerical stability)")
        self.console.print("    • Use [yellow]CPU engines[/yellow] for NPT equilibration (pressure/density stability)")
        
        self.console.print("\n  [cyan]For CPU-only systems:[/cyan]")
        self.console.print("    • Use [green]pmemd.MPI[/green] for production simulations (if available)")
        self.console.print("    • Use [green]pmemd[/green] for single-core optimized performance")
        self.console.print("    • Use [green]sander[/green] when special features are needed")
        
        # Simulation-specific recommendations
        self.console.print("\n[bold]Simulation-Specific Recommendations:[/bold]")
        self.console.print("\n  [cyan]Minimization:[/cyan]")
        self.console.print("    Prefer CPU engines (pmemd, sander) for numerical precision")
        
        self.console.print("\n  [cyan]NPT Equilibration:[/cyan]")
        self.console.print("    Prefer CPU engines for pressure/density stability")
        self.console.print("    GPU calculations may show pressure fluctuations")
        
        self.console.print("\n  [cyan]NVT Heating/Production:[/cyan]")
        self.console.print("    pmemd.cuda provides best performance (if GPU available)")
        self.console.print("    pmemd.MPI for parallel CPU performance")
        
        # Current queue analysis
        if self.simulation_queue:
            self.console.print("\n[bold]Your Current Queue:[/bold]")
            min_count = sum(1 for c in self.simulation_queue if 'minim' in c.name.lower())
            heat_count = sum(1 for c in self.simulation_queue if 'heat' in c.name.lower())
            npt_count = sum(1 for c in self.simulation_queue if 'npt' in c.name.lower() or 'equil' in c.name.lower())
            prod_count = sum(1 for c in self.simulation_queue if 'prod' in c.name.lower())
            
            if min_count > 0:
                self.console.print(f"    • {min_count} minimization(s) - [yellow]Consider CPU engine[/yellow]")
            if npt_count > 0:
                self.console.print(f"    • {npt_count} NPT equilibration(s) - [yellow]Consider CPU engine[/yellow]")
            if heat_count > 0:
                self.console.print(f"    • {heat_count} heating simulation(s) - [green]GPU recommended if available[/green]")
            if prod_count > 0:
                self.console.print(f"    • {prod_count} production run(s) - [green]GPU recommended if available[/green]")
        
        self.console.print("\n[grey50]Press Enter to continue...[/grey50]")
        input()
            
    def _get_available_engines(self):
        """Get list of available AMBER engines."""
        # This would normally detect available engines on the system
        # For now, return a standard list
        return ["sander", "pmemd", "pmemd.MPI", "pmemd.cuda"]
        
    def _configure_gpu_resources(self):
        """Configure GPU resources for CUDA engines."""
        self.console.print("\n[bold]GPU Configuration:[/bold]")

        # Get GPU info from system
        gpu_info = self._get_gpu_info()
        gpu_count = gpu_info.get('available', 1)

        if gpu_count == 0:
            self.console.print("[yellow]No GPUs detected on this machine. Defaulting to GPU 0.[/yellow]")
            self.console.print("[grey50]  This is expected if setting up on a login node to run on a compute node.[/grey50]")
            return {"gpu_ids": "0"}
        elif gpu_count == 1:
            self.console.print("Using GPU 0 (only GPU available)")
            return {"gpu_ids": "0"}
        else:
            # Prompt for GPU IDs (can be comma-separated for multiple GPUs)
            gpu_ids_str = prompt_with_context(
                self.processor,
                f"GPU IDs to use (comma-separated, 0-{gpu_count-1})",
                default="0",
                module="MD Manager - Hardware Configuration",
                description="Enter GPU device IDs for CUDA execution"
            )
            return {"gpu_ids": gpu_ids_str}
            
    def _configure_mpi_resources(self):
        """Configure MPI resources for parallel engines."""
        self.console.print("\n[bold]MPI Configuration:[/bold]")
        
        # Get available cores from system
        cpu_info = self._get_cpu_info()
        max_cores = cpu_info.get('available', 8)

        # Default to reasonable number of cores (16 or half of available, whichever is smaller)
        default_cores = min(max_cores // 2, 16) if max_cores > 1 else 1

        num_processes_str = prompt_with_context(
            self.processor,
            f"Number of MPI processes (1-{max_cores})",
            default=str(default_cores),
            module="MD Manager - Hardware Configuration",
            description="Configure MPI process count"
        )
        num_processes = int(num_processes_str)
        return {"mpi_tasks": num_processes}

    # ========================================================================
    # TRAJECTORY ANALYSIS METHODS
    # ========================================================================

    def _browse_for_trajectory_files(self, start_dir=None):
        """Browse for .nc trajectory files (multi-select for concatenation).

        Thin wrapper over the shared file browser. Returns a list of Paths
        (single picks come back as a one-element list), None on cancel, or the
        result of the recursive `find` command.
        """
        from pathlib import Path
        from datetime import datetime
        from proprep.utils.file_browser import file_browser

        start = Path(start_dir) if start_dir else self.working_directory

        def _nc_detail(p):
            try:
                size_mb = os.path.getsize(p) / (1024 ** 2)
                date_str = datetime.fromtimestamp(os.path.getmtime(p)).strftime("%m/%d/%Y %H:%M")
                return f"{size_mb:.1f} MB, {date_str}"
            except OSError:
                return ""

        extra = {
            "find": ("Recursive search for .nc files",
                     lambda cur: self._find_trajectory_files(Path(cur))),
        }
        return file_browser(
            directory=str(start),
            extensions=[".nc"],
            console=self.console,
            processor=self.processor,
            multi=True,
            label="trajectory file",
            entry_detail=_nc_detail,
            path_factory=Path,
            extra_commands=extra,
            module="MD Manager - Trajectory Browser",
        )



    def _find_trajectory_files(self, current_dir):
        """Recursive search for .nc trajectory files."""
        self.console.print(f"\n[bold cyan]Searching for .nc files...[/bold cyan]")

        nc_files = list(current_dir.rglob("*.nc"))

        if not nc_files:
            self.console.print("[grey50]No .nc files found[/grey50]")
            return None

        self.console.print(f"[green]Found {len(nc_files)} files:[/green]")

        for i, file_path in enumerate(nc_files[:20], 1):  # Limit display to 20
            relative = file_path.relative_to(current_dir)
            size_mb = file_path.stat().st_size / (1024**2)
            self.console.print(f"  {i:2}. {relative} ({size_mb:.1f} MB)")

        if len(nc_files) > 20:
            self.console.print(f"  [grey50]... and {len(nc_files) - 20} more[/grey50]")

        choice = prompt_with_context(
            self.processor, f"Select file (1-{min(20, len(nc_files))}) or 'cancel'",
            default="cancel", module="MD Manager - File Search",
            description="Select from recursively-found trajectory files",
        ).strip()

        if choice.lower() == 'cancel':
            return None

        try:
            file_num = int(choice)
            if 1 <= file_num <= min(20, len(nc_files)):
                selected = nc_files[file_num - 1]
                self.console.print(f"[green]Selected: {selected.name}[/green]")
                return [selected]
        except ValueError:
            pass

        return None

    def _concatenate_trajectories(self, nc_files: list, prmtop: Path):
        """
        Load and concatenate multiple trajectory segments.

        Uses pytraj to load multiple files as single trajectory.
        Mirrors _concatenate_monitors() for mdout files.

        Args:
            nc_files: List of trajectory file paths (in order)
            prmtop: Topology file path

        Returns:
            TrajectoryAnalyzer: Analyzer with concatenated trajectory
        """
        from .trajectory_analyzer import TrajectoryAnalyzer
        import pytraj as pt

        self.console.print(f"\n[bold cyan]Loading {len(nc_files)} Trajectory Segments[/bold cyan]")
        for i, f in enumerate(nc_files, 1):
            size_mb = f.stat().st_size / (1024**2)
            self.console.print(f"  {i}. {f.name} ({size_mb:.1f} MB)")

        self.console.print(f"\n[grey50]Concatenating segments...[/grey50]")

        try:
            # pytraj can load multiple files directly
            traj = pt.load([str(f) for f in nc_files], top=str(prmtop))

            total_frames = len(traj)
            self.console.print(f"[green]✓ Loaded {total_frames:,} total frames[/green]")

            # Create analyzer with pre-loaded trajectory
            analyzer = TrajectoryAnalyzer(str(prmtop), traj_object=traj)

            return analyzer

        except Exception as e:
            self.console.print(f"[red]Error concatenating trajectories: {e}[/red]")
            raise

    def _get_analysis_region_selection(self, analyzer, analysis_type: str) -> dict:
        """
        Interactive region selection for trajectory analysis.

        Similar to MD Restraint Manager's residue selection but for atoms/residues.

        Args:
            analyzer: TrajectoryAnalyzer instance
            analysis_type: Name of analysis (for display)

        Returns:
            dict: {'mask': str, 'description': str}
        """
        self.console.print(f"\n[bold cyan]Select Region for {analysis_type}[/bold cyan]\n")

        # Get number of residues for protein selection
        n_residues = analyzer.system_info.get('n_residues', 0)

        # Build options based on system
        # Detect protein residue range (exclude common solvent)
        # Common solvent residues: WAT, Na+, Cl-, etc.
        # For now, we'll use a heuristic: protein is typically first ~90% of residues
        # User can use custom selection for exact ranges
        if n_residues > 100:
            # Likely has solvent - estimate protein range
            protein_end = int(n_residues * 0.8)  # Conservative estimate
            protein_desc = f"Protein residues 1-{protein_end} (no H, excludes solvent)"
            protein_mask = f':1-{protein_end}&!@H='
        else:
            # Small system, likely no solvent
            protein_desc = f"Protein residues 1-{n_residues} (no H)"
            protein_mask = f':1-{n_residues}&!@H='

        options_map = {
            "1": "Backbone atoms (@C,CA,N,O&!:WAT)",
            "2": "C-alpha atoms (@CA)",
            "3": "All atoms (*)",
            "4": protein_desc,
            "5": "Specific residues (custom)",
            "6": "Custom AMBER mask (advanced)"
        }

        self.console.print("[bold]Common Selections:[/bold]")
        for key, desc in options_map.items():
            self.console.print(f"  {key}. {desc}")

        choice = prompt_with_context(
            self.processor,
            "Select region",
            choices=["1", "2", "3", "4", "5", "6"],
            default="2",
            module="MD Manager - Trajectory Analysis",
            description=f"Region for {analysis_type}",
            options_map=options_map
        )

        if choice == "1":
            return {
                'mask': '@C,CA,N,O&!:WAT',
                'description': 'Backbone atoms'
            }

        elif choice == "2":
            return {
                'mask': '@CA',
                'description': 'C-alpha atoms'
            }

        elif choice == "3":
            return {
                'mask': '*',
                'description': 'All atoms'
            }

        elif choice == "4":
            return {
                'mask': protein_mask,
                'description': protein_desc
            }

        elif choice == "5":
            return self._get_specific_residue_mask(analyzer)

        else:  # choice == "6"
            return self._get_custom_mask()

    def _get_specific_residue_mask(self, analyzer) -> dict:
        """Get residue specification from user."""
        self.console.print("\n[bold]Residue Specification Examples:[/bold]")
        self.console.print("  42            - Single residue")
        self.console.print("  10-50         - Residue range")
        self.console.print("  10-50,75-100  - Multiple ranges")
        self.console.print("  10,15,20,42   - Specific residues")

        residue_spec = prompt_with_context(
            self.processor,
            "Enter residue specification",
            module="MD Manager - Trajectory Analysis",
            description="Specify residues"
        )

        # Build AMBER mask
        mask = f":{residue_spec}"

        # Ask if user wants specific atoms within these residues
        atom_choice = prompt_with_context(
            self.processor,
            "Atom selection within these residues",
            choices=["1", "2", "3"],
            default="2",
            module="MD Manager - Trajectory Analysis",
            description="Select atoms",
            options_map={
                "1": "All atoms in these residues",
                "2": "C-alpha only (recommended)",
                "3": "Custom atom selection"
            }
        )

        if atom_choice == "1":
            pass  # Use mask as-is
        elif atom_choice == "2":
            mask += "&@CA"
        else:
            atom_mask = prompt_with_context(
                self.processor,
                "Enter atom selection (e.g., @CA,C,N,O)",
                module="MD Manager - Trajectory Analysis",
                description="Atom mask"
            )
            mask += f"&{atom_mask}"

        return {
            'mask': mask,
            'description': f'Residues {residue_spec}'
        }

    def _get_custom_mask(self) -> dict:
        """Get custom AMBER mask from user with examples."""
        self.console.print("\n[bold]AMBER Mask Syntax:[/bold]")
        self.console.print("[grey50]Format: :[residues]@[atoms]&[operators][/grey50]\n")

        self.console.print("[bold]Examples:[/bold]")
        self.console.print("  @CA                - All C-alpha atoms")
        self.console.print("  :1-50              - All atoms in residues 1-50")
        self.console.print("  :ALA,GLY           - All alanine and glycine residues")
        self.console.print("  :1-50@CA           - C-alpha in residues 1-50")
        self.console.print("  :1-50&!@H=         - Residues 1-50, exclude hydrogens")
        self.console.print("  @CA,C,N,O          - Multiple atom types")
        self.console.print("  :ALA@CB            - CB atoms of alanines")

        self.console.print("\n[bold]Operators:[/bold]")
        self.console.print("  &  - AND")
        self.console.print("  |  - OR")
        self.console.print("  !  - NOT")

        mask = prompt_with_context(
            self.processor,
            "Enter AMBER mask",
            module="MD Manager - Trajectory Analysis",
            description="Custom AMBER mask"
        )

        return {
            'mask': mask,
            'description': f'Custom mask: {mask}'
        }

    def _analyze_trajectory(self, nc_files: list, prmtop: Path, sim_name: str):
        """
        Main trajectory analysis workflow.

        Args:
            nc_files: List of trajectory file paths
            prmtop: Topology file path
            sim_name: Simulation name for display
        """
        from .trajectory_analyzer import TrajectoryAnalyzer

        try:
            # Load trajectory (single or concatenated)
            if len(nc_files) == 1:
                self.console.print(f"\n[bold cyan]Loading Trajectory[/bold cyan]")
                size_mb = nc_files[0].stat().st_size / (1024**2)
                self.console.print(f"[grey50]File: {nc_files[0].name} ({size_mb:.1f} MB)[/grey50]")

                analyzer = TrajectoryAnalyzer(str(prmtop), str(nc_files[0]))
            else:
                analyzer = self._concatenate_trajectories(nc_files, prmtop)

            self.console.print(f"[green]✓ Loaded successfully[/green]")

            # Display trajectory info
            self._display_trajectory_info(analyzer, sim_name)

            # Analysis loop - keep showing menu until user exits
            while True:
                # Analysis selection menu
                selected_analyses = self._get_analysis_selection_menu()

                if selected_analyses is None:
                    # User selected exit (option 0)
                    self.console.print("[cyan]Returning to MD Manager menu...[/cyan]")
                    return

                if not selected_analyses:
                    self.console.print("[yellow]No analyses selected[/yellow]")
                    continue  # Show menu again

                # Perform selected analyses
                self._execute_trajectory_analyses(analyzer, selected_analyses, sim_name)

                # Export option
                self._offer_trajectory_export(analyzer, sim_name)

                # Loop back to menu for more analyses

        except ImportError as e:
            self.console.print(f"[red]Error: {e}[/red]")
            self.console.print(f"[yellow]Install pytraj: conda install -c conda-forge pytraj ambertools[/yellow]")
        except Exception as e:
            self.console.print(f"[red]Error in trajectory analysis: {e}[/red]")
            import traceback
            traceback.print_exc()

    def _display_trajectory_info(self, analyzer, sim_name: str):
        """Display trajectory information."""
        from rich.table import Table

        self.console.print("\n" + "="*70)
        self.console.print(f"TRAJECTORY ANALYSIS: {sim_name}")
        self.console.print("="*70)

        table = Table(title="Trajectory Information", show_header=False)
        table.add_column("Property", style="bright_blue")
        table.add_column("Value", style="white")

        table.add_row("Frames", f"{analyzer.system_info['n_frames']:,}")
        table.add_row("Atoms", f"{analyzer.system_info['n_atoms']:,}")
        table.add_row("Residues", f"{analyzer.system_info['n_residues']:,}")

        if analyzer.frame_times:
            time_range = f"{analyzer.frame_times[0]:.1f} - {analyzer.frame_times[-1]:.1f} ps"
            duration = analyzer.frame_times[-1] - analyzer.frame_times[0]
            duration_ns = duration / 1000.0
            table.add_row("Time range", f"{time_range} ({duration_ns:.1f} ns)")

        table.add_row("Box type", analyzer.system_info.get('box_type', 'unknown'))

        self.console.print(table)

    def _get_analysis_selection_menu(self) -> list:
        """
        Display categorized menu for selecting trajectory analyses.

        Returns:
            list: Selected analysis IDs, or None to exit
        """
        self.console.print("\n[bold bright_blue]" + "="*70 + "[/bold bright_blue]")
        self.console.print("[bold bright_blue]                TRAJECTORY ANALYSIS MENU[/bold bright_blue]")
        self.console.print("[bold bright_blue]" + "="*70 + "[/bold bright_blue]\n")

        analyses = {
            "1": "RMSD - Root mean square deviation",
            "2": "RMSF - Per-residue fluctuations",
            "3": "Contacts - Contact maps & native contacts",
            "4": "Salt bridges - Electrostatic interactions",
            "5": "SASA - Solvent accessible surface area",
            "6": "DSSP - Secondary structure (helix, sheet, coil)",
            "7": "Ramachandran - Backbone dihedral angles (φ, ψ)",
            "8": "B-factors - Pseudo B-factors from MD",
            "9": "Distance - Atom-atom distances",
            "10": "Angle - Three-atom angles",
            "11": "Dihedral - Four-atom torsion angles",
            "12": "Vector - Orientation tracking",
            "13": "PCA - Principal component analysis",
            "14": "Clustering - Conformational clustering",
            "15": "Autocorrelation - Temporal correlations",
            "16": "Pairwise RMSD - All-vs-all RMSD matrix",
            "17": "Hydrogen bonds - H-bond network",
            "18": "Water RDF - Radial distribution function",
            "19": "Water shells - Hydration layer analysis",
            "20": "Density maps - Spatial density distributions",
            "21": "Radius of gyration - Structural compactness",
            "22": "Contact frequency - Per-residue contact analysis",
            "0": "Exit to MD Manager menu"
        }

        # Display categorized menu — styling matches the main-menu convention:
        # bold-blue "══ section ══" headers and a bold-blue name scan column,
        # with highlight=False so Rich's ReprHighlighter doesn't tint the numbers.
        sections = [
            ("STRUCTURAL ANALYSIS", ["1", "2", "3", "4", "5"]),
            ("SECONDARY STRUCTURE & GEOMETRY", ["6", "7", "8"]),
            ("PAIRWISE MEASUREMENTS", ["9", "10", "11", "12"]),
            ("DYNAMICS & CORRELATION", ["13", "14", "15", "16"]),
            ("SOLVATION ANALYSIS", ["17", "18", "19", "20"]),
            ("GEOMETRIC PROPERTIES", ["21", "22"]),
        ]
        for title, keys in sections:
            self.console.print(f"\n[bold blue]══ {title} ══[/bold blue]", highlight=False)
            for key in keys:
                self.console.print(f"  {key:>2}. [bold blue]{analyses[key]}[/bold blue]", highlight=False)

        self.console.print(f"\n  {'0':>2}. [bold blue]{analyses['0']}[/bold blue]", highlight=False)

        selection = prompt_with_context(
            self.processor,
            "\nSelect analyses (comma-separated, e.g., 1,2,6)",
            module="MD Manager - Trajectory Analysis",
            description="Select analysis types"
        )

        # Parse selection
        selected = []
        if '0' in selection:
            # Exit option
            return None  # Signal to exit
        else:
            selected = [s.strip() for s in selection.split(',') if s.strip() in analyses and s.strip() != '0']

        if selected:
            selected_names = [analyses[s].split('-')[0].strip() for s in selected]
            self.console.print(f"\n[green]Selected:[/green] {', '.join(selected_names)}")

        return selected

    def _execute_trajectory_analyses(self, analyzer, selected_analyses: list, sim_name: str):
        """Execute selected trajectory analyses."""
        total = len(selected_analyses)

        for i, analysis_id in enumerate(selected_analyses, 1):
            self.console.print(f"\n{'='*70}")
            self.console.print(f"Analysis {i}/{total}")
            self.console.print("="*70)

            # Structural Analysis
            if analysis_id == "1":
                self._analyze_rmsd(analyzer)
            elif analysis_id == "2":
                self._analyze_rmsf(analyzer)
            elif analysis_id == "3":
                self._analyze_contacts(analyzer)
            elif analysis_id == "4":
                self._analyze_salt_bridges(analyzer)
            elif analysis_id == "5":
                self._analyze_sasa(analyzer)

            # Secondary Structure & Geometry
            elif analysis_id == "6":
                self._analyze_dssp(analyzer)
            elif analysis_id == "7":
                self._analyze_ramachandran(analyzer)
            elif analysis_id == "8":
                self._analyze_bfactors(analyzer)

            # Pairwise Measurements
            elif analysis_id == "9":
                self._analyze_distance(analyzer)
            elif analysis_id == "10":
                self._analyze_angle(analyzer)
            elif analysis_id == "11":
                self._analyze_dihedral(analyzer)
            elif analysis_id == "12":
                self._analyze_vector(analyzer)

            # Dynamics & Correlation
            elif analysis_id == "13":
                self._analyze_pca(analyzer)
            elif analysis_id == "14":
                self._analyze_clustering(analyzer)
            elif analysis_id == "15":
                self._analyze_autocorrelation(analyzer)
            elif analysis_id == "16":
                self._analyze_pairwise_rmsd(analyzer)

            # Solvation Analysis
            elif analysis_id == "17":
                self._analyze_hbonds(analyzer)
            elif analysis_id == "18":
                self._analyze_water_rdf(analyzer)
            elif analysis_id == "19":
                self._analyze_water_shells(analyzer)
            elif analysis_id == "20":
                self._analyze_density_maps(analyzer)

            # Geometric Properties
            elif analysis_id == "21":
                self._analyze_radius_of_gyration(analyzer)
            elif analysis_id == "22":
                self._analyze_contact_frequency_per_residue(analyzer)

    def _analyze_rmsd(self, analyzer):
        """Perform RMSD analysis."""
        from rich.table import Table

        self.console.print("\n[bold cyan]RMSD Analysis[/bold cyan]")

        # Get region selection
        region = self._get_analysis_region_selection(analyzer, "RMSD")

        # Get reference frame
        self.console.print("\n[bold]Select Reference Frame:[/bold]")
        self.console.print("  1. Frame 0 (default)")
        self.console.print("  2. Average structure")
        self.console.print("  3. Specific frame")

        ref_choice = prompt_with_context(
            self.processor,
            "Select reference frame",
            choices=["1", "2", "3"],
            default="1",
            module="MD Manager - RMSD",
            description="Reference frame selection",
            options_map={
                "1": "Frame 0 (default)",
                "2": "Average structure",
                "3": "Specific frame"
            }
        )

        reference = 0
        if ref_choice == "2":
            reference = -1  # pytraj uses -1 for average
        elif ref_choice == "3":
            frame_str = prompt_with_context(
                self.processor,
                f"Enter frame number (0-{analyzer.system_info['n_frames']-1})",
                module="MD Manager - RMSD",
                description="Specify reference frame"
            )
            reference = int(frame_str)

        # Calculate
        self.console.print(f"\n[grey50]Calculating RMSD for {region['description']}...[/grey50]")
        rmsd = analyzer.calculate_rmsd(
            mask=region['mask'],
            reference=reference,
            label=f"rmsd_{region['description'].replace(' ', '_')}"
        )
        self.console.print(f"[green]✓ Complete[/green]")

        # Display results
        self._display_rmsd_results(analyzer, region, reference)

        # Offer additional RMSD analyses
        additional = prompt_with_context(
            self.processor,
            "Analyze additional regions? (y/n)",
            choices=["y", "n"],
            default="n",
            module="MD Manager - RMSD",
            description="Additional RMSD",
            options_map={"y": "Yes", "n": "No"}
        )

        if additional == "y":
            self._analyze_rmsd(analyzer)  # Recursive for additional regions

    def _display_rmsd_results(self, analyzer, region: dict, reference: int):
        """Display RMSD analysis results."""
        from rich.table import Table
        import numpy as np

        # Get the most recent RMSD data
        rmsd_keys = [k for k in analyzer.data.keys() if k.startswith('rmsd_') and not k.endswith('_mask') and not k.endswith('_reference')]
        if not rmsd_keys:
            self.console.print("[red]Error: No RMSD data found[/red]")
            return

        latest_key = rmsd_keys[-1]
        rmsd_values = analyzer.data.get(latest_key)

        if not rmsd_values:
            self.console.print(f"[red]Error: No data found for key '{latest_key}'[/red]")
            return

        # Statistics
        stats = analyzer.get_statistics(latest_key)

        if not stats:
            self.console.print("[red]Error: Could not calculate statistics for RMSD data[/red]")
            return

        # Display
        self.console.print(f"\n[bold]RMSD Analysis: {region['description']}[/bold]")
        self.console.print(f"[grey50]Mask: {region['mask']}[/grey50]")
        ref_desc = f"Frame {reference}" if reference >= 0 else "Average structure"
        self.console.print(f"[grey50]Reference: {ref_desc}[/grey50]\n")

        # Statistics table
        table = Table(title="RMSD Statistics", show_header=False)
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        table.add_row("Mean", f"{stats['mean']:.2f} ± {stats['std']:.2f} Å")
        table.add_row("Range", f"{stats['min']:.2f} - {stats['max']:.2f} Å")
        table.add_row("Final", f"{rmsd_values[-1]:.2f} Å")

        self.console.print(table)

        # ASCII plot
        self.console.print("\n" + "="*70)
        self.console.print("RMSD vs Time")
        self.console.print("="*70)

        plot = self._create_ascii_plot(
            rmsd_values,
            title=f"RMSD: {region['description']}",
            xlabel="Time (ps)",
            ylabel="RMSD (Å)",
            x_values=analyzer.frame_times if analyzer.frame_times else None
        )
        self.console.print(plot)

    def _analyze_rmsf(self, analyzer):
        """Perform RMSF analysis."""
        from rich.table import Table

        self.console.print("\n[bold cyan]RMSF Analysis[/bold cyan]")

        # Get region selection
        region = self._get_analysis_region_selection(analyzer, "RMSF")

        # Calculate
        self.console.print(f"\n[grey50]Calculating RMSF for {region['description']}...[/grey50]")
        rmsf = analyzer.calculate_rmsf(
            mask=region['mask'],
            label=f"rmsf_{region['description'].replace(' ', '_')}"
        )
        self.console.print(f"[green]✓ Complete[/green]")

        # Display results
        self._display_rmsf_results(analyzer, region)

    def _display_rmsf_results(self, analyzer, region: dict):
        """Display RMSF analysis results."""
        from rich.table import Table
        import numpy as np

        # Get the most recent RMSF data (exclude metadata keys)
        rmsf_keys = [k for k in analyzer.data.keys()
                     if k.startswith('rmsf_')
                     and not k.endswith('_mask')
                     and not k.endswith('_residue_indices')]
        if not rmsf_keys:
            self.console.print("[red]Error: No RMSF data found[/red]")
            return

        latest_key = rmsf_keys[-1]
        rmsf_values = analyzer.data.get(latest_key)

        if not rmsf_values:
            self.console.print(f"[red]Error: No data found for key '{latest_key}'[/red]")
            return

        # Get residue indices if available (for per-residue RMSF)
        residue_indices_key = f'{latest_key}_residue_indices'
        residue_indices = analyzer.data.get(residue_indices_key)

        # Statistics
        stats = analyzer.get_statistics(latest_key)

        # Check if statistics are available
        if not stats:
            self.console.print("[red]Error: Could not calculate statistics for RMSF data[/red]")
            return

        # Display
        self.console.print(f"\n[bold]RMSF Analysis: {region['description']}[/bold]")
        self.console.print(f"[grey50]Mask: {region['mask']}[/grey50]")
        self.console.print(f"[grey50]Calculation: Per-residue (mass-weighted average)[/grey50]\n")

        # Statistics table
        table = Table(title="RMSF Statistics", show_header=False)
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        table.add_row("Mean", f"{stats['mean']:.2f} ± {stats['std']:.2f} Å")
        table.add_row("Range", f"{stats['min']:.2f} - {stats['max']:.2f} Å")
        table.add_row("Median", f"{stats['median']:.2f} Å")

        self.console.print(table)

        # Find high-fluctuation regions (just reporting, no threshold)
        rmsf_array = np.array(rmsf_values)

        # Handle both 1D and 2D arrays (pytraj can return different formats)
        if rmsf_array.ndim > 1:
            # If 2D, take the mean across frames for each residue
            rmsf_array = np.mean(rmsf_array, axis=0)

        # Report top 10% as "high fluctuation" for information
        percentile_90 = np.percentile(rmsf_array, 90)

        high_fluct_indices = np.where(rmsf_array >= percentile_90)[0]
        if len(high_fluct_indices) > 0:
            self.console.print(f"\n[bold]High Fluctuation Residues[/bold] (top 10%, n={len(high_fluct_indices)}):")
            for idx in high_fluct_indices:
                # Use actual residue number if available, otherwise use sequential
                if residue_indices:
                    residue_num = residue_indices[idx]
                else:
                    residue_num = idx + 1
                self.console.print(f"  Residue {residue_num}: {rmsf_array[idx]:.2f} Å")

        # ASCII plot
        self.console.print("\n" + "="*70)
        self.console.print("Per-Residue RMSF")
        self.console.print("="*70)

        # Create x-values (residue numbers)
        if residue_indices:
            x_values = residue_indices
        else:
            x_values = list(range(1, len(rmsf_array) + 1))

        plot = self._create_ascii_plot(
            rmsf_array.tolist() if hasattr(rmsf_array, 'tolist') else rmsf_array,
            title=f"RMSF: {region['description']}",
            xlabel="Residue Number",
            ylabel="RMSF (Å)",
            x_values=x_values
        )
        self.console.print(plot)

    def _analyze_distance(self, analyzer):
        """Perform distance measurement analysis."""
        from rich.table import Table

        self.console.print("\n[bold cyan]Distance Measurement[/bold cyan]")

        # Get two atom selections
        self.console.print("\n[bold]Select First Atom/Group:[/bold]")
        mask1 = prompt_with_context(
            self.processor,
            "Enter first atom mask (e.g., :42@CA)",
            module="MD Manager - Distance",
            description="First atom selection"
        )

        self.console.print("\n[bold]Select Second Atom/Group:[/bold]")
        mask2 = prompt_with_context(
            self.processor,
            "Enter second atom mask (e.g., :89@CA)",
            module="MD Manager - Distance",
            description="Second atom selection"
        )

        # Optional label
        label = prompt_with_context(
            self.processor,
            "Enter description for this distance (optional)",
            default="distance",
            module="MD Manager - Distance",
            description="Label for analysis"
        )

        # Calculate
        self.console.print(f"\n[grey50]Calculating distance between {mask1} and {mask2}...[/grey50]")
        distances = analyzer.calculate_distance(mask1, mask2, label=label)
        self.console.print(f"[green]✓ Complete[/green]")

        # Display results
        self._display_distance_results(analyzer, mask1, mask2, label)

    def _display_distance_results(self, analyzer, mask1: str, mask2: str, label: str):
        """Display distance analysis results."""
        from rich.table import Table

        # Get the distance data
        dist_keys = [k for k in analyzer.data.keys() if label in k and not k.endswith('_mask1') and not k.endswith('_mask2')]
        if not dist_keys:
            self.console.print("[red]Error: No distance data found[/red]")
            return

        latest_key = dist_keys[-1]
        dist_values = analyzer.data.get(latest_key)

        if not dist_values:
            self.console.print(f"[red]Error: No data found for key '{latest_key}'[/red]")
            return

        # Statistics
        stats = analyzer.get_statistics(latest_key)

        if not stats:
            self.console.print("[red]Error: Could not calculate statistics for distance data[/red]")
            return

        # Display
        self.console.print(f"\n[bold]Distance Analysis: {label}[/bold]")
        self.console.print(f"[grey50]From: {mask1}[/grey50]")
        self.console.print(f"[grey50]To: {mask2}[/grey50]\n")

        # Statistics table
        table = Table(title="Distance Statistics", show_header=False)
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        table.add_row("Mean", f"{stats['mean']:.2f} ± {stats['std']:.2f} Å")
        table.add_row("Range", f"{stats['min']:.2f} - {stats['max']:.2f} Å")
        table.add_row("Median", f"{stats['median']:.2f} Å")

        self.console.print(table)

        # ASCII plot - time series
        self.console.print("\n" + "="*70)
        self.console.print("Distance vs Time")
        self.console.print("="*70)

        plot = self._create_ascii_plot(
            dist_values,
            title=f"Distance: {label}",
            xlabel="Time (ps)",
            ylabel="Distance (Å)",
            x_values=analyzer.frame_times if analyzer.frame_times else None
        )
        self.console.print(plot)

        # Distribution histogram
        self.console.print("\n" + "="*70)
        self.console.print("Distance Distribution")
        self.console.print("="*70)

        hist_plot = self._create_ascii_histogram(
            dist_values,
            title="Distance Distribution",
            xlabel="Distance (Å)"
        )
        self.console.print(hist_plot, highlight=False)

    def _analyze_angle(self, analyzer):
        """Perform angle measurement analysis."""
        from rich.table import Table

        self.console.print("\n[bold cyan]Angle Measurement[/bold cyan]")

        # Get three atom selections
        self.console.print("\n[bold]Select Three Atoms:[/bold]")
        mask1 = prompt_with_context(
            self.processor,
            "First atom (e.g., :42@CA)",
            module="MD Manager - Angle",
            description="First atom"
        )

        mask2 = prompt_with_context(
            self.processor,
            "Second atom - vertex (e.g., :43@CA)",
            module="MD Manager - Angle",
            description="Vertex atom"
        )

        mask3 = prompt_with_context(
            self.processor,
            "Third atom (e.g., :44@CA)",
            module="MD Manager - Angle",
            description="Third atom"
        )

        # Optional label
        label = prompt_with_context(
            self.processor,
            "Enter description (optional)",
            default="angle",
            module="MD Manager - Angle",
            description="Label"
        )

        # Calculate
        self.console.print(f"\n[grey50]Calculating angle...[/grey50]")
        angles = analyzer.calculate_angle(mask1, mask2, mask3, label=label)
        self.console.print(f"[green]✓ Complete[/green]")

        # Display results
        self._display_angle_results(analyzer, mask1, mask2, mask3, label)

    def _display_angle_results(self, analyzer, mask1: str, mask2: str, mask3: str, label: str):
        """Display angle analysis results."""
        from rich.table import Table

        # Get the angle data
        angle_keys = [k for k in analyzer.data.keys() if label in k and not k.endswith('_masks')]
        if not angle_keys:
            self.console.print("[red]Error: No angle data found[/red]")
            return

        latest_key = angle_keys[-1]
        angle_values = analyzer.data.get(latest_key)

        if not angle_values:
            self.console.print(f"[red]Error: No data found for key '{latest_key}'[/red]")
            return

        # Statistics
        stats = analyzer.get_statistics(latest_key)

        if not stats:
            self.console.print("[red]Error: Could not calculate statistics for angle data[/red]")
            return

        # Display
        self.console.print(f"\n[bold]Angle Analysis: {label}[/bold]")
        self.console.print(f"[grey50]Atoms: {mask1} - {mask2} - {mask3}[/grey50]\n")

        # Statistics table
        table = Table(title="Angle Statistics", show_header=False)
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        table.add_row("Mean", f"{stats['mean']:.2f} ± {stats['std']:.2f}°")
        table.add_row("Range", f"{stats['min']:.2f} - {stats['max']:.2f}°")
        table.add_row("Median", f"{stats['median']:.2f}°")

        self.console.print(table)

        # ASCII plot
        self.console.print("\n" + "="*70)
        self.console.print("Angle vs Time")
        self.console.print("="*70)

        plot = self._create_ascii_plot(
            angle_values,
            title=f"Angle: {label}",
            xlabel="Time (ps)",
            ylabel="Angle (degrees)",
            x_values=analyzer.frame_times if analyzer.frame_times else None
        )
        self.console.print(plot)

    def _analyze_dihedral(self, analyzer):
        """Perform dihedral angle measurement analysis."""
        from rich.table import Table

        self.console.print("\n[bold cyan]Dihedral Angle Measurement[/bold cyan]")

        # Get four atom selections
        self.console.print("\n[bold]Select Four Atoms:[/bold]")
        mask1 = prompt_with_context(
            self.processor,
            "First atom (e.g., :42@C)",
            module="MD Manager - Dihedral",
            description="First atom"
        )

        mask2 = prompt_with_context(
            self.processor,
            "Second atom (e.g., :43@N)",
            module="MD Manager - Dihedral",
            description="Second atom"
        )

        mask3 = prompt_with_context(
            self.processor,
            "Third atom (e.g., :43@CA)",
            module="MD Manager - Dihedral",
            description="Third atom"
        )

        mask4 = prompt_with_context(
            self.processor,
            "Fourth atom (e.g., :43@C)",
            module="MD Manager - Dihedral",
            description="Fourth atom"
        )

        # Optional label
        label = prompt_with_context(
            self.processor,
            "Enter description (optional, e.g., 'Psi angle Res43')",
            default="dihedral",
            module="MD Manager - Dihedral",
            description="Label"
        )

        # Calculate
        self.console.print(f"\n[grey50]Calculating dihedral...[/grey50]")
        dihedrals = analyzer.calculate_dihedral(mask1, mask2, mask3, mask4, label=label)
        self.console.print(f"[green]✓ Complete[/green]")

        # Display results
        self._display_dihedral_results(analyzer, mask1, mask2, mask3, mask4, label)

    def _display_dihedral_results(self, analyzer, mask1: str, mask2: str,
                                  mask3: str, mask4: str, label: str):
        """Display dihedral analysis results."""
        from rich.table import Table

        # Get the dihedral data
        dih_keys = [k for k in analyzer.data.keys() if label in k and not k.endswith('_masks')]
        if not dih_keys:
            self.console.print("[red]Error: No dihedral data found[/red]")
            return

        latest_key = dih_keys[-1]
        dih_values = analyzer.data.get(latest_key)

        if not dih_values:
            self.console.print(f"[red]Error: No data found for key '{latest_key}'[/red]")
            return

        # Statistics
        stats = analyzer.get_statistics(latest_key)

        if not stats:
            self.console.print("[red]Error: Could not calculate statistics for dihedral data[/red]")
            return

        # Display
        self.console.print(f"\n[bold]Dihedral Analysis: {label}[/bold]")
        self.console.print(f"[grey50]Atoms: {mask1} - {mask2} - {mask3} - {mask4}[/grey50]\n")

        # Statistics table
        table = Table(title="Dihedral Statistics", show_header=False)
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        table.add_row("Mean", f"{stats['mean']:.2f} ± {stats['std']:.2f}°")
        table.add_row("Range", f"{stats['min']:.2f} - {stats['max']:.2f}°")
        table.add_row("Median", f"{stats['median']:.2f}°")

        self.console.print(table)
        self.console.print("[grey50]Note: Dihedral angles range from -180° to +180°[/grey50]")

        # ASCII plot
        self.console.print("\n" + "="*70)
        self.console.print("Dihedral vs Time")
        self.console.print("="*70)

        plot = self._create_ascii_plot(
            dih_values,
            title=f"Dihedral: {label}",
            xlabel="Time (ps)",
            ylabel="Dihedral (degrees)",
            x_values=analyzer.frame_times if analyzer.frame_times else None
        )
        self.console.print(plot)

    def _analyze_hbonds(self, analyzer):
        """Perform hydrogen bond analysis."""
        from rich.table import Table

        self.console.print("\n[bold cyan]Hydrogen Bond Analysis[/bold cyan]")

        # Get parameters
        self.console.print("\n[bold]H-Bond Criteria:[/bold]")

        donor_mask = prompt_with_context(
            self.processor,
            "Donor mask (leave empty for all)",
            default="",
            module="MD Manager - H-Bonds",
            description="Donor selection"
        )
        if not donor_mask:
            donor_mask = None

        acceptor_mask = prompt_with_context(
            self.processor,
            "Acceptor mask (leave empty for all)",
            default="",
            module="MD Manager - H-Bonds",
            description="Acceptor selection"
        )
        if not acceptor_mask:
            acceptor_mask = None

        distance_cutoff = float(prompt_with_context(
            self.processor,
            "Distance cutoff (Angstroms)",
            default="3.0",
            module="MD Manager - H-Bonds",
            description="Max donor-acceptor distance"
        ))

        angle_cutoff = float(prompt_with_context(
            self.processor,
            "Angle cutoff (degrees)",
            default="135",
            module="MD Manager - H-Bonds",
            description="Min donor-H-acceptor angle"
        ))

        # Calculate
        self.console.print(f"\n[grey50]Calculating H-bonds...[/grey50]")
        self.console.print(f"[grey50]Note: H-bond analysis may take some time for large trajectories[/grey50]")

        try:
            hbond_data = analyzer.calculate_hbonds(
                donor_mask=donor_mask,
                acceptor_mask=acceptor_mask,
                distance_cutoff=distance_cutoff,
                angle_cutoff=angle_cutoff
            )
            self.console.print(f"[green]✓ Complete[/green]")

            # Display results
            self._display_hbond_results(analyzer, distance_cutoff, angle_cutoff)
        except Exception as e:
            self.console.print(f"[yellow]H-bond calculation encountered issues: {e}[/yellow]")
            self.console.print(f"[grey50]Showing basic statistics...[/grey50]")
            self._display_hbond_results(analyzer, distance_cutoff, angle_cutoff)

    def _display_hbond_results(self, analyzer, distance_cutoff: float, angle_cutoff: float):
        """Display H-bond analysis results."""
        from rich.table import Table

        # Get H-bond data (exclude metadata keys)
        hbond_keys = [k for k in analyzer.data.keys()
                      if k.startswith('hbonds')
                      and not k.endswith('_distance_cutoff')
                      and not k.endswith('_angle_cutoff')
                      and not k.endswith('_pair_occupancies')]
        if not hbond_keys:
            self.console.print("[yellow]No H-bond data available[/yellow]")
            return

        latest_key = hbond_keys[-1]
        hbond_counts = analyzer.data.get(latest_key)

        if not hbond_counts:
            self.console.print(f"[red]Error: No data found for key '{latest_key}'[/red]")
            return

        # Get H-bond pair occupancies if available
        occupancies_key = f'{latest_key}_pair_occupancies'
        pair_occupancies = analyzer.data.get(occupancies_key, [])

        # Statistics
        stats = analyzer.get_statistics(latest_key)

        if not stats:
            self.console.print("[red]Error: Could not calculate statistics for H-bond data[/red]")
            return

        # Display
        self.console.print(f"\n[bold]Hydrogen Bond Analysis[/bold]")
        self.console.print(f"[grey50]Distance cutoff: {distance_cutoff} Å[/grey50]")
        self.console.print(f"[grey50]Angle cutoff: {angle_cutoff}°[/grey50]")
        if pair_occupancies:
            self.console.print(f"[grey50]Unique donor-acceptor pairs: {len(pair_occupancies)}[/grey50]\n")
        else:
            self.console.print()

        # Statistics table
        table = Table(title="H-Bond Statistics", show_header=False)
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        table.add_row("Mean", f"{stats['mean']:.1f} ± {stats['std']:.1f} bonds")
        table.add_row("Range", f"{stats['min']:.0f} - {stats['max']:.0f} bonds")
        table.add_row("Median", f"{stats['median']:.1f} bonds")

        self.console.print(table)

        # Show most persistent H-bonds
        if pair_occupancies:
            self.console.print(f"\n[bold]Most Persistent H-Bonds[/bold] (top 30 by occupancy):")
            self.console.print(f"[grey50]Occupancy = % of frames where H-bond is present[/grey50]\n")

            # Show top 30 most persistent
            top_n = min(30, len(pair_occupancies))
            for i, hb_data in enumerate(pair_occupancies[:top_n], 1):
                pair_name = hb_data['pair']
                occupancy = hb_data['occupancy']

                # Create occupancy bar (50 chars max)
                bar_length = int(occupancy / 2)  # Scale to 50 chars max
                bar = "█" * bar_length

                self.console.print(f"  {i:2d}. {pair_name:40s} {occupancy:5.1f}% {bar}")

        # ASCII plot
        self.console.print("\n" + "="*70)
        self.console.print("H-Bonds vs Time")
        self.console.print("="*70)

        plot = self._create_ascii_plot(
            hbond_counts,
            title="Hydrogen Bonds over Time",
            xlabel="Time (ps)",
            ylabel="Number of H-Bonds",
            x_values=analyzer.frame_times if analyzer.frame_times else None
        )
        self.console.print(plot)

    def _analyze_radius_of_gyration(self, analyzer):
        """Perform radius of gyration analysis."""
        from rich.table import Table

        self.console.print("\n[bold cyan]Radius of Gyration Analysis[/bold cyan]")

        # Get region selection
        region = self._get_analysis_region_selection(analyzer, "Radius of Gyration")

        # Calculate
        self.console.print(f"\n[grey50]Calculating radius of gyration for {region['description']}...[/grey50]")
        rgyr = analyzer.calculate_radius_of_gyration(
            mask=region['mask'],
            label=f"rgyr_{region['description'].replace(' ', '_')}"
        )
        self.console.print(f"[green]✓ Complete[/green]")

        # Display results
        self._display_rgyr_results(analyzer, region)

    def _display_rgyr_results(self, analyzer, region: dict):
        """Display Rg analysis results."""
        from rich.table import Table

        # Get Rg data (exclude metadata keys)
        rgyr_keys = [k for k in analyzer.data.keys() if k.startswith('rgyr_') and not k.endswith('_mask')]
        if not rgyr_keys:
            self.console.print("[red]Error: No Rg data found[/red]")
            return

        latest_key = rgyr_keys[-1]
        rgyr_values = analyzer.data.get(latest_key)

        if not rgyr_values:
            self.console.print(f"[red]Error: No data found for key '{latest_key}'[/red]")
            return

        # Statistics
        stats = analyzer.get_statistics(latest_key)

        if not stats:
            self.console.print("[red]Error: Could not calculate statistics for Rg data[/red]")
            return

        # Display
        self.console.print(f"\n[bold]Radius of Gyration: {region['description']}[/bold]")
        self.console.print(f"[grey50]Mask: {region['mask']}[/grey50]\n")

        # Statistics table
        table = Table(title="Rg Statistics", show_header=False)
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        table.add_row("Mean", f"{stats['mean']:.2f} ± {stats['std']:.2f} Å")
        table.add_row("Range", f"{stats['min']:.2f} - {stats['max']:.2f} Å")
        table.add_row("Median", f"{stats['median']:.2f} Å")

        self.console.print(table)

        # ASCII plot
        self.console.print("\n" + "="*70)
        self.console.print("Radius of Gyration vs Time")
        self.console.print("="*70)

        plot = self._create_ascii_plot(
            rgyr_values,
            title=f"Rg: {region['description']}",
            xlabel="Time (ps)",
            ylabel="Rg (Å)",
            x_values=analyzer.frame_times if analyzer.frame_times else None
        )
        self.console.print(plot)

    def _analyze_water_rdf(self, analyzer):
        """Perform water radial distribution function analysis."""
        from rich.table import Table

        self.console.print("\n[bold cyan]Water Radial Distribution Function[/bold cyan]")

        # Get solute selection
        self.console.print("\n[bold]Solute Selection:[/bold]")
        solute_mask = prompt_with_context(
            self.processor,
            "Enter solute mask (e.g., :42-45 for residues 42-45)",
            module="MD Manager - RDF",
            description="Solute selection"
        )

        # Solvent selection (default to water oxygens)
        solvent_mask = prompt_with_context(
            self.processor,
            "Solvent mask",
            default=":WAT@O",
            module="MD Manager - RDF",
            description="Solvent selection"
        )

        # RDF parameters
        max_distance = float(prompt_with_context(
            self.processor,
            "Maximum distance (Angstroms)",
            default="10.0",
            module="MD Manager - RDF",
            description="Max distance for RDF"
        ))

        bin_spacing = float(prompt_with_context(
            self.processor,
            "Bin spacing (Angstroms)",
            default="0.1",
            module="MD Manager - RDF",
            description="Histogram bin width"
        ))

        # Calculate
        self.console.print(f"\n[grey50]Calculating RDF...[/grey50]")
        try:
            distances, gr_values = analyzer.calculate_water_radial_distribution(
                solute_mask=solute_mask,
                solvent_mask=solvent_mask,
                max_distance=max_distance,
                bin_spacing=bin_spacing
            )
            self.console.print(f"[green]✓ Complete[/green]")

            # Display results
            self._display_rdf_results(analyzer, solute_mask, solvent_mask, distances, gr_values)
        except Exception as e:
            self.console.print(f"[red]Error calculating RDF: {e}[/red]")

    def _display_rdf_results(self, analyzer, solute_mask: str, solvent_mask: str,
                            distances, gr_values):
        """Display RDF analysis results."""
        from rich.table import Table
        import numpy as np

        # Display
        self.console.print(f"\n[bold]Radial Distribution Function[/bold]")
        self.console.print(f"[grey50]Solute: {solute_mask}[/grey50]")
        self.console.print(f"[grey50]Solvent: {solvent_mask}[/grey50]\n")

        # Find peaks (simple peak finding)
        gr_array = np.array(gr_values)
        dist_array = np.array(distances)

        # First shell peak (within first 5 Å)
        first_shell_mask = dist_array < 5.0
        if np.any(first_shell_mask):
            first_shell_gr = gr_array[first_shell_mask]
            first_shell_dist = dist_array[first_shell_mask]
            first_peak_idx = np.argmax(first_shell_gr)

            table = Table(title="Solvation Shell Information", show_header=False)
            table.add_column("Property", style="bright_blue")
            table.add_column("Value", style="white")

            table.add_row("First shell peak", f"{first_shell_dist[first_peak_idx]:.2f} Å")
            table.add_row("Peak g(r)", f"{first_shell_gr[first_peak_idx]:.2f}")

            self.console.print(table)

        # ASCII plot of g(r)
        self.console.print("\n" + "="*70)
        self.console.print("Radial Distribution Function g(r)")
        self.console.print("="*70)

        plot = self._create_ascii_plot(
            gr_values,
            title="g(r) vs Distance",
            ylabel="g(r)",
            x_values=distances
        )
        self.console.print(plot)

    def _analyze_sasa(self, analyzer):
        """Perform solvent accessible surface area analysis."""
        from rich.table import Table

        self.console.print("\n[bold cyan]Solvent Accessible Surface Area[/bold cyan]")

        # Get region selection
        region = self._get_analysis_region_selection(analyzer, "SASA")

        # Probe radius
        probe_radius = float(prompt_with_context(
            self.processor,
            "Probe radius (Angstroms)",
            default="1.4",
            module="MD Manager - SASA",
            description="Solvent probe radius (1.4 for water)"
        ))

        # Calculate
        self.console.print(f"\n[grey50]Calculating SASA for {region['description']}...[/grey50]")
        try:
            sasa = analyzer.calculate_sasa(
                mask=region['mask'],
                probe_radius=probe_radius,
                label=f"sasa_{region['description'].replace(' ', '_')}"
            )
            self.console.print(f"[green]✓ Complete[/green]")

            # Display results
            self._display_sasa_results(analyzer, region, probe_radius)
        except Exception as e:
            self.console.print(f"[red]Error calculating SASA: {e}[/red]")

    def _display_sasa_results(self, analyzer, region: dict, probe_radius: float):
        """Display SASA analysis results."""
        from rich.table import Table

        # Get SASA data (exclude metadata keys)
        sasa_keys = [k for k in analyzer.data.keys() if k.startswith('sasa_') and not k.endswith('_mask') and not k.endswith('_probe')]
        if not sasa_keys:
            self.console.print("[red]Error: No SASA data found[/red]")
            return

        latest_key = sasa_keys[-1]
        sasa_values = analyzer.data.get(latest_key)

        if not sasa_values:
            self.console.print(f"[red]Error: No data found for key '{latest_key}'[/red]")
            return

        # Statistics
        stats = analyzer.get_statistics(latest_key)

        if not stats:
            self.console.print("[red]Error: Could not calculate statistics for SASA data[/red]")
            return

        # Display
        self.console.print(f"\n[bold]SASA Analysis: {region['description']}[/bold]")
        self.console.print(f"[grey50]Mask: {region['mask']}[/grey50]")
        self.console.print(f"[grey50]Probe radius: {probe_radius} Å[/grey50]\n")

        # Statistics table
        table = Table(title="SASA Statistics", show_header=False)
        table.add_column("Metric", style="bright_blue")
        table.add_column("Value", style="white")

        table.add_row("Mean", f"{stats['mean']:.1f} ± {stats['std']:.1f} Ų")
        table.add_row("Range", f"{stats['min']:.1f} - {stats['max']:.1f} Ų")
        table.add_row("Median", f"{stats['median']:.1f} Ų")

        self.console.print(table)

        # ASCII plot
        self.console.print("\n" + "="*70)
        self.console.print("SASA vs Time")
        self.console.print("="*70)

        plot = self._create_ascii_plot(
            sasa_values,
            title=f"SASA: {region['description']}",
            xlabel="Time (ps)",
            ylabel="SASA (Ų)",
            x_values=analyzer.frame_times if analyzer.frame_times else None
        )
        self.console.print(plot)

    def _analyze_dssp(self, analyzer):
        """Perform DSSP secondary structure analysis."""
        from rich.table import Table

        self.console.print("\n[bold cyan]Secondary Structure Analysis (DSSP)[/bold cyan]")

        # Get region selection
        region = self._get_analysis_region_selection(analyzer, "DSSP")

        self.console.print("\nCalculating secondary structure...")
        self.console.print("[grey50]Note: This requires dssp executable in PATH[/grey50]")
        self.console.print("[grey50]Install via: conda install -c salilab dssp[/grey50]")

        try:
            dssp_data = analyzer.calculate_dssp(mask=region['mask'])
            self.console.print(f"[green]✓ Complete[/green]")

            # Display results
            self._display_dssp_results(analyzer, region)
        except Exception as e:
            self.console.print(f"[red]DSSP calculation failed: {e}[/red]")
            self.console.print(f"[yellow]Make sure dssp is installed and in PATH[/yellow]")

    def _display_dssp_results(self, analyzer, region):
        """Display DSSP secondary structure results."""
        from rich.table import Table
        import numpy as np

        # Get DSSP data (exclude metadata keys)
        dssp_keys = [k for k in analyzer.data.keys()
                     if k == 'dssp_helix_pct' or k.startswith('dssp_') and k.endswith('_helix_pct')]

        if not dssp_keys:
            self.console.print("[red]Error: No DSSP data found[/red]")
            return

        # Extract the base key (remove '_helix_pct' suffix)
        base_key = dssp_keys[0].replace('_helix_pct', '')

        # Get all secondary structure data
        helix_pct = analyzer.data.get(f'{base_key}_helix_pct')
        sheet_pct = analyzer.data.get(f'{base_key}_sheet_pct')
        turn_pct = analyzer.data.get(f'{base_key}_turn_pct')
        coil_pct = analyzer.data.get(f'{base_key}_coil_pct')

        if not all([helix_pct, sheet_pct, turn_pct, coil_pct]):
            self.console.print("[red]Error: Incomplete DSSP data[/red]")
            return

        # Calculate statistics for each SS type
        helix_stats = analyzer.get_statistics(f'{base_key}_helix_pct')
        sheet_stats = analyzer.get_statistics(f'{base_key}_sheet_pct')
        turn_stats = analyzer.get_statistics(f'{base_key}_turn_pct')
        coil_stats = analyzer.get_statistics(f'{base_key}_coil_pct')

        # Display
        self.console.print(f"\n[bold]Secondary Structure Analysis: {region['description']}[/bold]")
        self.console.print(f"[grey50]Mask: {region['mask']}[/grey50]\n")

        # Summary statistics table
        table = Table(title="Secondary Structure Statistics", show_header=True)
        table.add_column("Structure", style="bright_blue")
        table.add_column("Mean %", style="white")
        table.add_column("Range %", style="white")

        table.add_row("α-Helix", f"{helix_stats['mean']:.1f} ± {helix_stats['std']:.1f}",
                     f"{helix_stats['min']:.1f} - {helix_stats['max']:.1f}")
        table.add_row("β-Sheet", f"{sheet_stats['mean']:.1f} ± {sheet_stats['std']:.1f}",
                     f"{sheet_stats['min']:.1f} - {sheet_stats['max']:.1f}")
        table.add_row("Turn", f"{turn_stats['mean']:.1f} ± {turn_stats['std']:.1f}",
                     f"{turn_stats['min']:.1f} - {turn_stats['max']:.1f}")
        table.add_row("Coil", f"{coil_stats['mean']:.1f} ± {coil_stats['std']:.1f}",
                     f"{coil_stats['min']:.1f} - {coil_stats['max']:.1f}")

        self.console.print(table)

        # Secondary structure evolution plot
        self.console.print("\n" + "="*70)
        self.console.print("Secondary Structure Evolution")
        self.console.print("="*70)

        # Create stacked area visualization - show all four types
        # We'll create individual plots for each type
        x_values = analyzer.frame_times if analyzer.frame_times else list(range(len(helix_pct)))

        # α-Helix
        self.console.print("\n[bold cyan]α-Helix Content[/bold cyan]")
        helix_plot = self._create_ascii_plot(
            helix_pct,
            title="α-Helix %",
            xlabel="Time (ps)",
            ylabel="% α-Helix",
            x_values=x_values
        )
        self.console.print(helix_plot, highlight=False)

        # β-Sheet
        self.console.print("\n[bold cyan]β-Sheet Content[/bold cyan]")
        sheet_plot = self._create_ascii_plot(
            sheet_pct,
            title="β-Sheet %",
            xlabel="Time (ps)",
            ylabel="% β-Sheet",
            x_values=x_values
        )
        self.console.print(sheet_plot, highlight=False)

        # Show composition summary
        self.console.print("\n[bold]Average Composition:[/bold]")
        total_structured = helix_stats['mean'] + sheet_stats['mean']
        total_flexible = turn_stats['mean'] + coil_stats['mean']

        self.console.print(f"  Structured regions (helix + sheet): {total_structured:.1f}%")
        self.console.print(f"  Flexible regions (turn + coil):     {total_flexible:.1f}%")

        # Visual bar chart of average composition
        self.console.print("\n[bold]Average Secondary Structure Composition:[/bold]")
        max_pct = max(helix_stats['mean'], sheet_stats['mean'], turn_stats['mean'], coil_stats['mean'])
        scale = 50.0 / max_pct if max_pct > 0 else 1.0

        helix_bar = "█" * int(helix_stats['mean'] * scale)
        sheet_bar = "█" * int(sheet_stats['mean'] * scale)
        turn_bar = "█" * int(turn_stats['mean'] * scale)
        coil_bar = "█" * int(coil_stats['mean'] * scale)

        self.console.print(f"  α-Helix: {helix_bar} {helix_stats['mean']:.1f}%")
        self.console.print(f"  β-Sheet: {sheet_bar} {sheet_stats['mean']:.1f}%")
        self.console.print(f"  Turn:    {turn_bar} {turn_stats['mean']:.1f}%")
        self.console.print(f"  Coil:    {coil_bar} {coil_stats['mean']:.1f}%")

    def _analyze_ramachandran(self, analyzer):
        """Perform Ramachandran/multidihedral analysis."""
        from rich.table import Table

        self.console.print("\n[bold cyan]Ramachandran / Backbone Dihedral Analysis[/bold cyan]")

        # Ask user what type of analysis
        self.console.print("\n[bold]Select Dihedral Type:[/bold]")
        self.console.print("  1. Phi-Psi (Ramachandran plot)")
        self.console.print("  2. Chi1 (sidechain rotamers)")
        self.console.print("  3. Omega (peptide bond planarity)")

        dihedral_choice = prompt_with_context(
            self.processor,
            "Select dihedral type",
            default="1",
            module="MD Manager - Ramachandran",
            description="Select dihedral type to analyze"
        )

        dihedral_map = {
            "1": "phi-psi",
            "2": "chi1",
            "3": "omega"
        }

        dihedral_type = dihedral_map.get(dihedral_choice, "phi-psi")

        # Ask for residue selection
        res_choice = prompt_with_context(
            self.processor,
            "Residue range (e.g., 1-50, or press Enter for all protein)",
            default="",
            module="MD Manager - Ramachandran",
            description="Specify residue range"
        )

        residue_selection = res_choice if res_choice.strip() else None

        self.console.print(f"\nCalculating {dihedral_type} dihedrals...")
        try:
            rama_data = analyzer.calculate_ramachandran(
                residue_selection=residue_selection,
                dihedral_type=dihedral_type
            )
            self.console.print(f"[green]✓ Complete[/green]")

            # Display results
            if dihedral_type == "phi-psi":
                self._display_ramachandran_results(analyzer, rama_data)
            else:
                self._display_dihedral_results_general(analyzer, rama_data, dihedral_type)

        except Exception as e:
            self.console.print(f"[red]Dihedral calculation failed: {e}[/red]")

    def _display_ramachandran_results(self, analyzer, rama_data):
        """Display Ramachandran plot and phi-psi analysis."""
        from rich.table import Table
        import numpy as np

        phi = rama_data['phi']  # shape: (n_residues, n_frames)
        psi = rama_data['psi']
        resrange = rama_data['resrange']

        n_residues, n_frames = phi.shape

        self.console.print(f"\n[bold]Ramachandran Analysis[/bold]")
        self.console.print(f"[grey50]Residues: {resrange}[/grey50]")
        self.console.print(f"[grey50]Frames: {n_frames}[/grey50]\n")

        # Classify residues into Ramachandran regions
        # α-helix: φ ≈ -60°, ψ ≈ -45°
        # β-sheet: φ ≈ -120°, ψ ≈ +120°
        # Left-handed helix: φ ≈ +60°, ψ ≈ +45°

        # Calculate mean phi/psi for each residue
        phi_mean = np.mean(phi, axis=1)
        psi_mean = np.mean(psi, axis=1)
        phi_std = np.std(phi, axis=1)
        psi_std = np.std(psi, axis=1)

        # Classify into regions
        helix_count = 0
        sheet_count = 0
        other_count = 0

        for i in range(n_residues):
            # α-helix region: φ in [-90, -30], ψ in [-75, -15]
            if -90 <= phi_mean[i] <= -30 and -75 <= psi_mean[i] <= -15:
                helix_count += 1
            # β-sheet region: φ in [-180, -90], ψ in [90, 180]
            elif -180 <= phi_mean[i] <= -90 and 90 <= psi_mean[i] <= 180:
                sheet_count += 1
            else:
                other_count += 1

        # Summary statistics
        table = Table(title="Ramachandran Classification", show_header=True)
        table.add_column("Region", style="bright_blue")
        table.add_column("Residues", style="white")
        table.add_column("Percentage", style="white")

        total = n_residues
        table.add_row("α-helix region", str(helix_count), f"{(helix_count/total)*100:.1f}%")
        table.add_row("β-sheet region", str(sheet_count), f"{(sheet_count/total)*100:.1f}%")
        table.add_row("Other", str(other_count), f"{(other_count/total)*100:.1f}%")

        self.console.print(table)

        # ASCII Ramachandran plot (aggregated over all frames)
        self.console.print("\n" + "="*70)
        self.console.print("Ramachandran Plot (φ vs ψ)")
        self.console.print("="*70)
        self.console.print("[grey50]All residues, all frames aggregated[/grey50]\n")

        # Create 2D ASCII scatter plot
        # Bin the data into a grid
        phi_flat = phi.flatten()
        psi_flat = psi.flatten()

        # Create 2D histogram
        phi_bins = np.linspace(-180, 180, 36)  # 10° bins
        psi_bins = np.linspace(-180, 180, 36)

        hist, _, _ = np.histogram2d(phi_flat, psi_flat, bins=[phi_bins, psi_bins])

        # Display as ASCII heatmap
        width = 60
        height = 20

        # Downsample histogram to fit display
        from scipy import ndimage
        hist_display = ndimage.zoom(hist, (height/hist.shape[0], width/hist.shape[1]), order=0)

        # Normalize and convert to characters
        max_count = hist_display.max()
        if max_count > 0:
            hist_display = hist_display / max_count

        chars = ' ░▒▓█'
        self.console.print("ψ (deg)")
        for j in range(height-1, -1, -1):
            line = ""
            for i in range(width):
                density = hist_display[j, i]
                char_idx = min(int(density * (len(chars)-1)), len(chars)-1)
                line += chars[char_idx]

            # Add y-axis labels
            y_val = -180 + (j / (height-1)) * 360
            self.console.print(f"{y_val:>6.0f}° │{line}")

        # X-axis
        self.console.print("       └" + "─" * width)
        self.console.print("        " + "-180°" + " " * 20 + "0°" + " " * 20 + "180°")
        self.console.print(" " * 30 + "φ (deg)\n")

        self.console.print("[grey50]█ = high density  ░ = low density[/grey50]")

        # Show residues with unusual dihedrals (outliers)
        self.console.print("\n[bold]Residues with High Flexibility:[/bold]")
        self.console.print("[grey50](Large φ/ψ fluctuations, std > 30°)[/grey50]\n")

        flexible_residues = []
        for i in range(min(10, n_residues)):  # Show up to 10
            if phi_std[i] > 30 or psi_std[i] > 30:
                res_num = int(resrange.split('-')[0]) + i
                flexible_residues.append((res_num, phi_mean[i], psi_mean[i], phi_std[i], psi_std[i]))

        if flexible_residues:
            flex_table = Table(show_header=True)
            flex_table.add_column("Residue", style="bright_blue")
            flex_table.add_column("φ (deg)", style="white")
            flex_table.add_column("ψ (deg)", style="white")

            for res_num, pm, psm, ps, pss in flexible_residues[:10]:
                flex_table.add_row(
                    str(res_num),
                    f"{pm:.1f} ± {ps:.1f}",
                    f"{psm:.1f} ± {pss:.1f}"
                )

            self.console.print(flex_table)
        else:
            self.console.print("  No highly flexible residues found")

    def _display_dihedral_results_general(self, analyzer, dihedral_data, dihedral_type):
        """Display results for chi or omega angles."""
        from rich.table import Table
        import numpy as np

        if 'chi' in dihedral_data:
            angles = dihedral_data['chi']
        elif 'omega' in dihedral_data:
            angles = dihedral_data['omega']
        else:
            self.console.print("[red]Error: No dihedral data found[/red]")
            return

        resrange = dihedral_data['resrange']
        n_residues, n_frames = angles.shape

        self.console.print(f"\n[bold]{dihedral_type.upper()} Dihedral Analysis[/bold]")
        self.console.print(f"[grey50]Residues: {resrange}[/grey50]")
        self.console.print(f"[grey50]Frames: {n_frames}[/grey50]\n")

        # Calculate statistics
        angle_mean = np.mean(angles, axis=1)
        angle_std = np.std(angles, axis=1)

        # Summary
        self.console.print(f"[bold]Summary Statistics:[/bold]")
        self.console.print(f"  Mean angle across all residues: {np.mean(angle_mean):.1f}°")
        self.console.print(f"  Overall std: {np.mean(angle_std):.1f}°\n")

        # For omega, check planarity
        if dihedral_type == 'omega':
            # Omega should be close to 180° (trans) or 0° (cis)
            # Deviation from 180° indicates non-planarity
            trans_count = np.sum(np.abs(angle_mean - 180) < 30)
            cis_count = np.sum(np.abs(angle_mean) < 30)

            self.console.print(f"[bold]Peptide Bond Geometry:[/bold]")
            self.console.print(f"  Trans peptides (ω ≈ 180°): {trans_count}/{n_residues}")
            self.console.print(f"  Cis peptides (ω ≈ 0°): {cis_count}/{n_residues}")

            # Find non-planar peptide bonds
            non_planar = []
            for i in range(n_residues):
                deviation = min(abs(angle_mean[i] - 180), abs(angle_mean[i]))
                if deviation > 30:
                    res_num = int(resrange.split('-')[0]) + i
                    non_planar.append((res_num, angle_mean[i], deviation))

            if non_planar:
                self.console.print(f"\n[yellow]Non-planar peptide bonds found:[/yellow]")
                for res_num, omega, dev in non_planar[:5]:
                    self.console.print(f"  Residue {res_num}: ω = {omega:.1f}° (deviation: {dev:.1f}°)")

    def _analyze_contacts(self, analyzer):
        """Perform contact map and native contact analysis."""

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]CONTACT MAP ANALYSIS[/bold cyan]")
        self.console.print("="*70)

        # Ask for contact type
        self.console.print("\n[bold]Select Contact Type:[/bold]")
        self.console.print("  1. Intra-protein contacts (C-alpha atoms)")
        self.console.print("  2. Inter-domain contacts (specify two regions)")
        self.console.print("  3. Custom selection (advanced)")

        contact_type = prompt_with_context(
            self.processor,
            "Select contact type",
            choices=["1", "2", "3"],
            default="1",
            module="MD Manager - Contact Maps",
            description="Select contact type to analyze",
            options_map={
                "1": "Intra-protein contacts (C-alpha atoms)",
                "2": "Inter-domain contacts (specify two regions)",
                "3": "Custom selection (advanced)"
            }
        )

        # Get mask selections based on choice
        if contact_type == "1":
            # Default C-alpha contacts
            mask1 = "@CA"
            mask2 = "@CA"
            self.console.print("\n[grey50]Using C-alpha atoms for contact analysis[/grey50]")

        elif contact_type == "2":
            # Inter-domain contacts
            self.console.print("\n[bold]Define Domain 1:[/bold]")
            self.console.print("Example: :1-50 (residues 1-50)")

            domain1 = prompt_with_context(
                self.processor,
                "Domain 1 selection (AMBER mask)",
                module="MD Manager - Contact Maps",
                description="Domain 1 selection"
            )

            self.console.print("\n[bold]Define Domain 2:[/bold]")
            self.console.print("Example: :51-100 (residues 51-100)")

            domain2 = prompt_with_context(
                self.processor,
                "Domain 2 selection (AMBER mask)",
                module="MD Manager - Contact Maps",
                description="Domain 2 selection"
            )

            mask1 = f"{domain1} & @CA"
            mask2 = f"{domain2} & @CA"

        else:
            # Custom masks
            self.console.print("\n[bold]Custom Selection:[/bold]")
            self.console.print("Examples:")
            self.console.print("  @CA          - C-alpha atoms")
            self.console.print("  :1-50@CA     - C-alpha in residues 1-50")
            self.console.print("  @C,CA,N,O    - Backbone atoms")

            mask1 = prompt_with_context(
                self.processor,
                "Selection 1 (AMBER mask)",
                module="MD Manager - Contact Maps",
                description="Selection 1"
            )

            mask2 = prompt_with_context(
                self.processor,
                "Selection 2 (AMBER mask, or press Enter for same as selection 1)",
                module="MD Manager - Contact Maps",
                description="Selection 2"
            )

            if not mask2.strip():
                mask2 = mask1

        # Get distance cutoff
        self.console.print("\n[bold]Distance Cutoff:[/bold]")
        self.console.print("  Default: 4.5 Å (typical for C-alpha contacts)")
        self.console.print("  Heavy atoms: 3.5-4.0 Å")
        self.console.print("  Include sidechain: 6.0-7.0 Å")

        cutoff_input = prompt_with_context(
            self.processor,
            "Distance cutoff (Å)",
            default="4.5",
            module="MD Manager - Contact Maps",
            description="Distance cutoff for contacts"
        )

        try:
            distance_cutoff = float(cutoff_input)
        except ValueError:
            distance_cutoff = 4.5

        # Get reference frame
        self.console.print("\n[bold]Reference Frame for Native Contacts:[/bold]")
        self.console.print("  1. Frame 0 (first frame)")
        self.console.print("  2. Last frame")
        self.console.print("  3. Specific frame number")

        ref_choice = prompt_with_context(
            self.processor,
            "Reference frame",
            choices=["1", "2", "3"],
            default="1",
            module="MD Manager - Contact Maps",
            description="Reference frame selection",
            options_map={
                "1": "Frame 0 (first frame)",
                "2": "Last frame",
                "3": "Specific frame number"
            }
        )

        if ref_choice == "1":
            reference_frame = 0
        elif ref_choice == "2":
            reference_frame = analyzer.traj.n_frames - 1
        else:
            frame_input = prompt_with_context(
                self.processor,
                "Frame number",
                default="0",
                module="MD Manager - Contact Maps",
                description="Reference frame number"
            )
            try:
                reference_frame = int(frame_input)
            except ValueError:
                reference_frame = 0

        # Perform analysis
        self.console.print(f"\n[cyan]Calculating contacts...[/cyan]")
        self.console.print(f"  Selection 1: {mask1}")
        self.console.print(f"  Selection 2: {mask2}")
        self.console.print(f"  Distance cutoff: {distance_cutoff} Å")
        self.console.print(f"  Reference frame: {reference_frame}")

        try:
            contact_data = analyzer.calculate_contacts(
                mask1=mask1,
                mask2=mask2,
                distance_cutoff=distance_cutoff,
                reference_frame=reference_frame
            )

            self.console.print("[green]✓ Analysis complete[/green]")

            # Display results
            self._display_contact_results(analyzer, contact_data)

        except Exception as e:
            self.console.print(f"[red]Error during contact analysis: {e}[/red]")
            logger.error(f"Contact analysis failed: {e}", exc_info=True)

    def _display_contact_results(self, analyzer, contact_data):
        """Display contact map analysis results."""
        from rich.table import Table
        import numpy as np

        q_values = contact_data['q_values']
        persistent_contacts = contact_data['persistent_contacts']
        n_native_contacts = contact_data['n_native_contacts']
        mask1 = contact_data['mask1']
        mask2 = contact_data['mask2']
        contact_frequency = contact_data['contact_frequency']

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]CONTACT ANALYSIS RESULTS[/bold cyan]")
        self.console.print("="*70)

        # Summary statistics
        self.console.print(f"\n[bold]Contact Summary:[/bold]")
        self.console.print(f"  Selection 1: {mask1}")
        self.console.print(f"  Selection 2: {mask2}")
        self.console.print(f"  Native contacts (reference): {n_native_contacts}")
        self.console.print(f"  Persistent contacts (>50% occupancy): {len(persistent_contacts)}")

        # Q-value statistics
        q_stats = {
            'mean': float(np.mean(q_values)),
            'std': float(np.std(q_values)),
            'min': float(np.min(q_values)),
            'max': float(np.max(q_values))
        }

        self.console.print(f"\n[bold]Q-value (Native Contact Retention):[/bold]")
        self.console.print(f"  Mean: {q_stats['mean']:.3f} ± {q_stats['std']:.3f}")
        self.console.print(f"  Range: {q_stats['min']:.3f} - {q_stats['max']:.3f}")

        # ASCII plot of Q-value over time
        self.console.print("\n" + "="*70)
        self.console.print("Q-VALUE VS TIME")
        self.console.print("="*70)

        # Get frame times
        frame_times = [analyzer.traj.time[i] for i in range(len(q_values))]

        q_plot = self._create_ascii_plot(
            q_values,
            title="Q-value (Fraction of Native Contacts)",
            xlabel="Time (ps)",
            ylabel="Q-value",
            x_values=frame_times
        )

        self.console.print(q_plot, highlight=False)

        # Display most persistent contacts
        if persistent_contacts:
            self.console.print("\n" + "="*70)
            self.console.print("[bold]MOST PERSISTENT CONTACTS (>50% occupancy)[/bold]")
            self.console.print("="*70)

            # Create table
            contact_table = Table(show_header=True, header_style="bold bright_blue")
            contact_table.add_column("Residue 1", style="bright_blue")
            contact_table.add_column("Residue 2", style="bright_blue")
            contact_table.add_column("Atoms", style="white")
            contact_table.add_column("Occupancy", style="green")
            contact_table.add_column("Native", style="yellow")

            # Show top 20 contacts
            for contact in persistent_contacts[:20]:
                res1_str = f"{contact['resname_i']}{contact['res_i']}"
                res2_str = f"{contact['resname_j']}{contact['res_j']}"
                atoms_str = f"{contact['atom_i']} - {contact['atom_j']}"
                occupancy_str = f"{contact['occupancy']*100:.1f}%"
                native_str = "✓" if contact['is_native'] else ""

                contact_table.add_row(
                    res1_str,
                    res2_str,
                    atoms_str,
                    occupancy_str,
                    native_str
                )

            self.console.print(contact_table)

            if len(persistent_contacts) > 20:
                self.console.print(f"\n[grey50]... and {len(persistent_contacts)-20} more persistent contacts[/grey50]")

        # Contact frequency heatmap (if not too large)
        if contact_frequency.shape[0] <= 50 and contact_frequency.shape[1] <= 50:
            self.console.print("\n" + "="*70)
            self.console.print("[bold]CONTACT FREQUENCY MAP[/bold]")
            self.console.print("="*70)

            # Create ASCII heatmap
            self._display_contact_heatmap(contact_frequency)
        else:
            self.console.print(f"\n[grey50]Contact matrix too large for display ({contact_frequency.shape[0]}x{contact_frequency.shape[1]})[/grey50]")

    def _display_contact_heatmap(self, contact_frequency):
        """Display contact frequency as ASCII heatmap."""
        import numpy as np

        # Downsample if needed
        max_size = 40
        if contact_frequency.shape[0] > max_size or contact_frequency.shape[1] > max_size:
            from scipy import ndimage
            zoom_factor = min(max_size / contact_frequency.shape[0], max_size / contact_frequency.shape[1])
            contact_frequency = ndimage.zoom(contact_frequency, zoom_factor, order=0)

        height, width = contact_frequency.shape

        # Characters for density levels
        chars = ' ░▒▓█'

        self.console.print(f"\n[grey50]Residue j →[/grey50]")
        for i in range(height):
            line = ""
            for j in range(width):
                freq = contact_frequency[i, j]
                char_idx = min(int(freq * (len(chars)-1)), len(chars)-1)
                line += chars[char_idx]

            if i == 0:
                self.console.print(f"Res i ↓ │{line}")
            else:
                self.console.print(f"       │{line}")

        self.console.print(f"\n[grey50]Legend: {chars[0]}=never  {chars[-1]}=always[/grey50]")

    def _analyze_salt_bridges(self, analyzer):
        """Perform salt bridge analysis."""

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]SALT BRIDGE ANALYSIS[/bold cyan]")
        self.console.print("="*70)

        self.console.print("\n[bold]Salt Bridge Detection:[/bold]")
        self.console.print("Salt bridges are electrostatic interactions between:")
        self.console.print("  • Acidic residues: ASP, GLU")
        self.console.print("  • Basic residues: LYS, ARG, HIS")

        # Get distance cutoff
        self.console.print("\n[bold]Distance Cutoff:[/bold]")
        self.console.print("  Default: 4.0 Å (standard for salt bridges)")
        self.console.print("  Relaxed: 5.0 Å (includes weaker interactions)")
        self.console.print("  Strict: 3.2 Å (only strong interactions)")

        cutoff_input = prompt_with_context(
            self.processor,
            "Distance cutoff (Å)",
            default="4.0",
            module="MD Manager - Salt Bridges",
            description="Distance cutoff for salt bridge detection"
        )

        try:
            distance_cutoff = float(cutoff_input)
        except ValueError:
            distance_cutoff = 4.0

        # Perform analysis
        self.console.print(f"\n[cyan]Analyzing salt bridges...[/cyan]")
        self.console.print(f"  Distance cutoff: {distance_cutoff} Å")

        try:
            salt_bridge_data = analyzer.calculate_salt_bridges(
                distance_cutoff=distance_cutoff
            )

            self.console.print("[green]✓ Analysis complete[/green]")

            # Display results
            self._display_salt_bridge_results(analyzer, salt_bridge_data)

        except Exception as e:
            self.console.print(f"[red]Error during salt bridge analysis: {e}[/red]")
            logger.error(f"Salt bridge analysis failed: {e}", exc_info=True)

    def _display_salt_bridge_results(self, analyzer, salt_bridge_data):
        """Display salt bridge analysis results."""
        from rich.table import Table
        import numpy as np

        salt_bridges = salt_bridge_data['salt_bridges']
        cutoff = salt_bridge_data['cutoff']
        n_bridges = salt_bridge_data['n_bridges']

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]SALT BRIDGE RESULTS[/bold cyan]")
        self.console.print("="*70)

        # Summary
        self.console.print(f"\n[bold]Summary:[/bold]")
        self.console.print(f"  Distance cutoff: {cutoff} Å")
        self.console.print(f"  Total salt bridges detected: {n_bridges}")

        if n_bridges == 0:
            self.console.print("\n[yellow]No salt bridges found in the trajectory[/yellow]")
            return

        # Count persistent salt bridges (>50% occupancy)
        persistent = [sb for sb in salt_bridges if sb['occupancy'] > 0.5]
        transient = [sb for sb in salt_bridges if sb['occupancy'] <= 0.5]

        self.console.print(f"  Persistent (>50% occupancy): {len(persistent)}")
        self.console.print(f"  Transient (≤50% occupancy): {len(transient)}")

        # Display persistent salt bridges
        if persistent:
            self.console.print("\n" + "="*70)
            self.console.print("[bold]PERSISTENT SALT BRIDGES (>50% occupancy)[/bold]")
            self.console.print("="*70)

            sb_table = Table(show_header=True, header_style="bold bright_blue")
            sb_table.add_column("Acidic", style="red")
            sb_table.add_column("Basic", style="blue")
            sb_table.add_column("Occupancy", style="green")
            sb_table.add_column("Distance (Å)", style="white")
            sb_table.add_column("Std (Å)", style="grey50")

            for sb in persistent[:15]:
                acidic_str = f"{sb['acidic_name']}{sb['acidic_res']}"
                basic_str = f"{sb['basic_name']}{sb['basic_res']}"
                occupancy_str = f"{sb['occupancy']*100:.1f}%"
                dist_str = f"{sb['mean_distance']:.2f}"
                std_str = f"{sb['std_distance']:.2f}"

                sb_table.add_row(
                    acidic_str,
                    basic_str,
                    occupancy_str,
                    dist_str,
                    std_str
                )

            self.console.print(sb_table)

            if len(persistent) > 15:
                self.console.print(f"\n[grey50]... and {len(persistent)-15} more persistent salt bridges[/grey50]")

        # Display transient salt bridges (top 10)
        if transient:
            self.console.print("\n" + "="*70)
            self.console.print("[bold]TRANSIENT SALT BRIDGES (top 10 by occupancy)[/bold]")
            self.console.print("="*70)

            trans_table = Table(show_header=True, header_style="bold bright_blue")
            trans_table.add_column("Acidic", style="red")
            trans_table.add_column("Basic", style="blue")
            trans_table.add_column("Occupancy", style="yellow")
            trans_table.add_column("Distance (Å)", style="white")

            for sb in transient[:10]:
                acidic_str = f"{sb['acidic_name']}{sb['acidic_res']}"
                basic_str = f"{sb['basic_name']}{sb['basic_res']}"
                occupancy_str = f"{sb['occupancy']*100:.1f}%"
                dist_str = f"{sb['mean_distance']:.2f}"

                trans_table.add_row(
                    acidic_str,
                    basic_str,
                    occupancy_str,
                    dist_str
                )

            self.console.print(trans_table)

        # Plot most persistent salt bridge distance over time
        if persistent:
            top_bridge = persistent[0]

            self.console.print("\n" + "="*70)
            self.console.print(f"[bold]DISTANCE TRAJECTORY: Most Persistent Salt Bridge[/bold]")
            self.console.print("="*70)

            self.console.print(f"\n[cyan]{top_bridge['acidic_name']}{top_bridge['acidic_res']} ← → "
                             f"{top_bridge['basic_name']}{top_bridge['basic_res']}[/cyan]")
            self.console.print(f"[grey50]Occupancy: {top_bridge['occupancy']*100:.1f}%[/grey50]")

            # Get frame times
            frame_times = [analyzer.traj.time[i] for i in range(len(top_bridge['distances']))]

            # Create ASCII plot
            dist_plot = self._create_ascii_plot(
                top_bridge['distances'],
                title=f"Distance (Å)",
                xlabel="Time (ps)",
                ylabel="Distance (Å)",
                x_values=frame_times
            )

            self.console.print(dist_plot, highlight=False)

            # Add cutoff reference
            self.console.print(f"\n[grey50]Cutoff threshold: {cutoff} Å[/grey50]")

        # Occupancy distribution histogram
        if n_bridges > 0:
            self.console.print("\n" + "="*70)
            self.console.print("[bold]SALT BRIDGE OCCUPANCY DISTRIBUTION[/bold]")
            self.console.print("="*70)

            occupancies = [sb['occupancy'] for sb in salt_bridges]

            # Create histogram
            hist_plot = self._create_ascii_histogram(
                occupancies,
                title="Salt Bridge Occupancy",
                xlabel="Occupancy",
                bins=10
            )

            self.console.print(hist_plot, highlight=False)

    def _analyze_pca(self, analyzer):
        """Perform Principal Component Analysis."""

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]PRINCIPAL COMPONENT ANALYSIS (PCA)[/bold cyan]")
        self.console.print("="*70)

        self.console.print("\n[bold]PCA Overview:[/bold]")
        self.console.print("PCA identifies collective motions by reducing dimensionality.")
        self.console.print("Principal components (PCs) capture the largest variance in motion.")

        # Get region selection
        region = self._get_analysis_region_selection()

        # Get number of components
        self.console.print("\n[bold]Number of Components:[/bold]")
        self.console.print("  Typical: 3-5 components (covers major motions)")
        self.console.print("  Extended: 10+ components (detailed analysis)")

        n_components_input = prompt_with_context(
            self.processor,
            "Number of principal components",
            default="3",
            module="MD Manager - PCA",
            description="Number of PCs to calculate"
        )

        try:
            n_components = int(n_components_input)
            if n_components < 1:
                n_components = 3
        except ValueError:
            n_components = 3

        # Perform analysis
        self.console.print(f"\n[cyan]Calculating PCA...[/cyan]")
        self.console.print(f"  Region: {region}")
        self.console.print(f"  Components: {n_components}")

        try:
            pca_data = analyzer.calculate_pca(
                mask=region,
                n_components=n_components
            )

            self.console.print("[green]✓ Analysis complete[/green]")

            # Display results
            self._display_pca_results(analyzer, pca_data)

        except Exception as e:
            self.console.print(f"[red]Error during PCA: {e}[/red]")
            logger.error(f"PCA failed: {e}", exc_info=True)

    def _display_pca_results(self, analyzer, pca_data):
        """Display PCA results."""
        from rich.table import Table
        import numpy as np

        variance_explained = pca_data['variance_explained']
        cumulative_variance = pca_data['cumulative_variance']
        projections = pca_data['projections']
        mask = pca_data['mask']
        n_components = pca_data['n_components']

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]PCA RESULTS[/bold cyan]")
        self.console.print("="*70)

        # Summary
        self.console.print(f"\n[bold]Analysis Summary:[/bold]")
        self.console.print(f"  Selection: {mask}")
        self.console.print(f"  Components: {n_components}")
        self.console.print(f"  Frames: {len(projections['PC1'])}")

        # Variance explained table
        self.console.print("\n" + "="*70)
        self.console.print("[bold]VARIANCE EXPLAINED BY PRINCIPAL COMPONENTS[/bold]")
        self.console.print("="*70)

        var_table = Table(show_header=True, header_style="bold bright_blue")
        var_table.add_column("Component", style="bright_blue")
        var_table.add_column("Variance (%)", style="green")
        var_table.add_column("Cumulative (%)", style="yellow")
        var_table.add_column("Bar", style="blue")

        for i in range(min(n_components, len(variance_explained))):
            pc_name = f"PC{i+1}"
            var_pct = variance_explained[i]
            cum_pct = cumulative_variance[i]

            # Create bar visualization
            bar_length = int(var_pct / 2)  # Scale to fit
            bar = "█" * bar_length

            var_table.add_row(
                pc_name,
                f"{var_pct:.2f}",
                f"{cum_pct:.2f}",
                bar
            )

        self.console.print(var_table)

        # Key insights
        self.console.print(f"\n[bold]Key Insights:[/bold]")
        self.console.print(f"  • PC1 captures {variance_explained[0]:.1f}% of motion variance")
        if len(variance_explained) >= 3:
            total_3pc = sum(variance_explained[:3])
            self.console.print(f"  • First 3 PCs capture {total_3pc:.1f}% of total variance")

        # Scree plot
        self.console.print("\n" + "="*70)
        self.console.print("[bold]SCREE PLOT (Variance by Component)[/bold]")
        self.console.print("="*70)

        # Create bar chart
        max_var = max(variance_explained[:min(10, len(variance_explained))])
        chart_width = 50

        for i in range(min(10, len(variance_explained))):
            var_pct = variance_explained[i]
            bar_len = int((var_pct / max_var) * chart_width)
            bar = "█" * bar_len

            self.console.print(f"  PC{i+1:>2} │{bar} {var_pct:.2f}%")

        # PC1 vs time
        self.console.print("\n" + "="*70)
        self.console.print("[bold]PC1 PROJECTION VS TIME[/bold]")
        self.console.print("="*70)

        pc1_values = projections['PC1']
        frame_times = [analyzer.traj.time[i] for i in range(len(pc1_values))]

        pc1_plot = self._create_ascii_plot(
            pc1_values,
            title="PC1 Projection",
            xlabel="Time (ps)",
            ylabel="PC1",
            x_values=frame_times
        )

        self.console.print(pc1_plot, highlight=False)

        # PC1 vs PC2 scatter plot (if PC2 exists)
        if 'PC2' in projections:
            self.console.print("\n" + "="*70)
            self.console.print("[bold]PC1 vs PC2 PROJECTION[/bold]")
            self.console.print("="*70)

            self._display_pca_scatter(projections['PC1'], projections['PC2'])

        # PC distribution histograms
        self.console.print("\n" + "="*70)
        self.console.print("[bold]PC1 DISTRIBUTION[/bold]")
        self.console.print("="*70)

        pc1_hist = self._create_ascii_histogram(
            pc1_values,
            title="PC1 Distribution",
            xlabel="PC1 Value",
            bins=20
        )

        self.console.print(pc1_hist)

    def _display_pca_scatter(self, pc1_values, pc2_values):
        """Display PC1 vs PC2 as ASCII scatter plot."""
        import numpy as np

        # Create 2D histogram for density
        pc1_array = np.array(pc1_values)
        pc2_array = np.array(pc2_values)

        # Define plot dimensions
        height = 30
        width = 60

        # Get ranges
        pc1_min, pc1_max = np.min(pc1_array), np.max(pc1_array)
        pc2_min, pc2_max = np.min(pc2_array), np.max(pc2_array)

        # Create 2D histogram
        hist, xedges, yedges = np.histogram2d(
            pc1_array, pc2_array,
            bins=[width, height],
            range=[[pc1_min, pc1_max], [pc2_min, pc2_max]]
        )

        # Characters for density levels
        chars = ' ·○●'

        # Plot
        self.console.print(f"\n[grey50]PC2 ↑[/grey50]")

        for j in range(height-1, -1, -1):
            line = ""
            for i in range(width):
                density = hist[i, j]
                if density == 0:
                    char_idx = 0
                else:
                    # Log scale for better visualization
                    log_density = np.log1p(density)
                    max_log = np.log1p(np.max(hist))
                    char_idx = min(int((log_density / max_log) * (len(chars)-1)) + 1, len(chars)-1)
                line += chars[char_idx]

            # Add axis labels
            if j == height // 2:
                y_val = pc2_min + (j / height) * (pc2_max - pc2_min)
                self.console.print(f"{y_val:>8.1f} │{line}")
            else:
                self.console.print(f"         │{line}")

        # X-axis
        x_axis = "─" * width
        self.console.print(f"         └{x_axis}→ PC1")

        # X-axis labels
        self.console.print(f"         {pc1_min:>8.1f}{' '*(width-16)}{pc1_max:>8.1f}")

        self.console.print(f"\n[grey50]Density: {chars[0]}=low  {chars[-1]}=high[/grey50]")

    def _analyze_clustering(self, analyzer):
        """Perform conformational clustering analysis."""

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]CONFORMATIONAL CLUSTERING[/bold cyan]")
        self.console.print("="*70)

        self.console.print("\n[bold]Clustering Overview:[/bold]")
        self.console.print("Groups similar conformations to identify representative structures.")
        self.console.print("Useful for analyzing conformational diversity and transitions.")

        # Get region selection
        region = self._get_analysis_region_selection()

        # Get clustering algorithm
        self.console.print("\n[bold]Clustering Algorithm:[/bold]")
        self.console.print("  1. K-means (fast, good for well-separated clusters)")
        self.console.print("  2. Hierarchical (comprehensive, better for overlapping clusters)")

        algo_choice = prompt_with_context(
            self.processor,
            "Select algorithm",
            choices=["1", "2"],
            default="1",
            module="MD Manager - Clustering",
            description="Select clustering algorithm",
            options_map={
                "1": "K-means",
                "2": "Hierarchical"
            }
        )

        algorithm = 'kmeans' if algo_choice == "1" else 'hierarchical'

        # Get number of clusters
        self.console.print("\n[bold]Number of Clusters:[/bold]")
        self.console.print("  Typical: 5-10 clusters")
        self.console.print("  Tip: Start with fewer clusters, refine if needed")

        n_clusters_input = prompt_with_context(
            self.processor,
            "Number of clusters",
            default="5",
            module="MD Manager - Clustering",
            description="Number of clusters to identify"
        )

        try:
            n_clusters = int(n_clusters_input)
            if n_clusters < 2:
                n_clusters = 5
        except ValueError:
            n_clusters = 5

        # Perform analysis
        self.console.print(f"\n[cyan]Performing clustering...[/cyan]")
        self.console.print(f"  Region: {region}")
        self.console.print(f"  Algorithm: {algorithm}")
        self.console.print(f"  Clusters: {n_clusters}")

        try:
            cluster_data = analyzer.calculate_clustering(
                mask=region,
                n_clusters=n_clusters,
                algorithm=algorithm
            )

            self.console.print("[green]✓ Clustering complete[/green]")

            # Display results
            self._display_clustering_results(analyzer, cluster_data)

        except Exception as e:
            self.console.print(f"[red]Error during clustering: {e}[/red]")
            logger.error(f"Clustering failed: {e}", exc_info=True)

    def _display_clustering_results(self, analyzer, cluster_data):
        """Display clustering analysis results."""
        from rich.table import Table
        import numpy as np

        assignments = cluster_data['assignments']
        populations = cluster_data['populations']
        representatives = cluster_data['representatives']
        cluster_rmsd = cluster_data['cluster_rmsd']
        n_clusters = cluster_data['n_clusters']
        algorithm = cluster_data['algorithm']
        mask = cluster_data['mask']

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]CLUSTERING RESULTS[/bold cyan]")
        self.console.print("="*70)

        # Summary
        self.console.print(f"\n[bold]Analysis Summary:[/bold]")
        self.console.print(f"  Selection: {mask}")
        self.console.print(f"  Algorithm: {algorithm}")
        self.console.print(f"  Number of clusters: {n_clusters}")
        self.console.print(f"  Total frames: {len(assignments)}")

        # Cluster populations table
        self.console.print("\n" + "="*70)
        self.console.print("[bold]CLUSTER POPULATIONS[/bold]")
        self.console.print("="*70)

        cluster_table = Table(show_header=True, header_style="bold bright_blue")
        cluster_table.add_column("Cluster", style="bright_blue")
        cluster_table.add_column("Frames", style="green")
        cluster_table.add_column("Population (%)", style="yellow")
        cluster_table.add_column("Representative", style="white")
        cluster_table.add_column("Avg RMSD (Å)", style="grey50")
        cluster_table.add_column("Bar", style="blue")

        total_frames = len(assignments)

        for i in range(n_clusters):
            count = populations[i]
            pct = (count / total_frames) * 100
            rep_frame = representatives[i]
            avg_rmsd = cluster_rmsd[i]['mean']

            # Create bar visualization
            bar_length = int(pct / 2)  # Scale to fit
            bar = "█" * bar_length

            cluster_table.add_row(
                str(i),
                str(count),
                f"{pct:.1f}",
                f"Frame {rep_frame}",
                f"{avg_rmsd:.2f}",
                bar
            )

        self.console.print(cluster_table)

        # Find most and least populated clusters
        max_pop_cluster = max(populations.items(), key=lambda x: x[1])
        min_pop_cluster = min(populations.items(), key=lambda x: x[1])

        self.console.print(f"\n[bold]Key Insights:[/bold]")
        self.console.print(f"  • Most populated: Cluster {max_pop_cluster[0]} ({max_pop_cluster[1]} frames, "
                         f"{(max_pop_cluster[1]/total_frames)*100:.1f}%)")
        self.console.print(f"  • Least populated: Cluster {min_pop_cluster[0]} ({min_pop_cluster[1]} frames, "
                         f"{(min_pop_cluster[1]/total_frames)*100:.1f}%)")

        # Cluster assignment timeline
        self.console.print("\n" + "="*70)
        self.console.print("[bold]CLUSTER ASSIGNMENT TIMELINE[/bold]")
        self.console.print("="*70)

        # Sample timeline (subsample if too many frames)
        max_timeline_points = 100
        if len(assignments) > max_timeline_points:
            step = len(assignments) // max_timeline_points
            sampled_assignments = assignments[::step]
            sampled_times = [analyzer.traj.time[i] for i in range(0, len(assignments), step)]
        else:
            sampled_assignments = assignments
            sampled_times = [analyzer.traj.time[i] for i in range(len(assignments))]

        # Create timeline visualization
        timeline_height = n_clusters
        timeline_width = min(len(sampled_assignments), 80)

        self.console.print(f"\n[grey50]Time →[/grey50]")

        # Build timeline grid
        for cluster_id in range(n_clusters-1, -1, -1):
            line = ""
            for i in range(timeline_width):
                idx = int((i / timeline_width) * len(sampled_assignments))
                if sampled_assignments[idx] == cluster_id:
                    line += "█"
                else:
                    line += " "

            self.console.print(f"Cluster {cluster_id} │{line}")

        # Time axis
        x_axis = "─" * timeline_width
        self.console.print(f"           └{x_axis}")

        if len(sampled_times) > 0:
            start_time = sampled_times[0]
            end_time = sampled_times[-1]
            self.console.print(f"           {start_time:>8.0f} ps{' '*(timeline_width-20)}{end_time:>8.0f} ps")

        # Population distribution pie chart (ASCII)
        self.console.print("\n" + "="*70)
        self.console.print("[bold]POPULATION DISTRIBUTION[/bold]")
        self.console.print("="*70)

        # Create horizontal bar chart
        max_count = max(populations.values())
        chart_width = 50

        for i in range(n_clusters):
            count = populations[i]
            pct = (count / total_frames) * 100
            bar_len = int((count / max_count) * chart_width)
            bar = "█" * bar_len

            self.console.print(f"  Cluster {i} │{bar} {pct:.1f}% ({count} frames)")

        # RMSD compactness
        self.console.print("\n" + "="*70)
        self.console.print("[bold]CLUSTER COMPACTNESS (Avg RMSD to Centroid)[/bold]")
        self.console.print("="*70)

        # Bar chart of RMSD
        rmsd_values = [cluster_rmsd[i]['mean'] for i in range(n_clusters)]
        max_rmsd = max(rmsd_values) if rmsd_values else 1.0

        for i in range(n_clusters):
            mean_rmsd = cluster_rmsd[i]['mean']
            std_rmsd = cluster_rmsd[i]['std']

            if max_rmsd > 0:
                bar_len = int((mean_rmsd / max_rmsd) * 40)
            else:
                bar_len = 0

            bar = "█" * bar_len

            self.console.print(f"  Cluster {i} │{bar} {mean_rmsd:.2f} ± {std_rmsd:.2f} Å")

        self.console.print(f"\n[grey50]Lower RMSD = more compact cluster[/grey50]")

    def _analyze_pairwise_rmsd(self, analyzer):
        """Perform pairwise RMSD matrix analysis."""

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]PAIRWISE RMSD MATRIX[/bold cyan]")
        self.console.print("="*70)

        self.console.print("\n[bold]Pairwise RMSD Overview:[/bold]")
        self.console.print("Computes all-vs-all RMSD to map conformational similarity.")
        self.console.print("Useful for identifying transitions and revisited conformations.")

        n_frames = analyzer.traj.n_frames

        # Warn about computational cost for large trajectories
        estimated_comparisons = (n_frames * (n_frames - 1)) // 2

        self.console.print(f"\n[bold]Trajectory Size:[/bold]")
        self.console.print(f"  Frames: {n_frames}")
        self.console.print(f"  Comparisons needed: {estimated_comparisons:,}")

        # Recommend subsampling for large trajectories
        subsample = None
        if n_frames > 500:
            self.console.print(f"\n[yellow]Large trajectory detected![/yellow]")
            self.console.print("  Recommended: Subsample to reduce computation time")

            self.console.print("\n[bold]Subsampling:[/bold]")
            self.console.print("  1. No subsampling (use all frames)")
            self.console.print(f"  2. Every 2nd frame ({n_frames//2} frames)")
            self.console.print(f"  3. Every 5th frame ({n_frames//5} frames)")
            self.console.print(f"  4. Every 10th frame ({n_frames//10} frames)")
            self.console.print("  5. Custom")

            subsample_choice = prompt_with_context(
                self.processor,
                "Subsampling option",
                choices=["1", "2", "3", "4", "5"],
                default="3",
                module="MD Manager - Pairwise RMSD",
                description="Select subsampling strategy",
                options_map={
                    "1": "No subsampling",
                    "2": "Every 2nd frame",
                    "3": "Every 5th frame",
                    "4": "Every 10th frame",
                    "5": "Custom"
                }
            )

            if subsample_choice == "1":
                subsample = None
            elif subsample_choice == "2":
                subsample = 2
            elif subsample_choice == "3":
                subsample = 5
            elif subsample_choice == "4":
                subsample = 10
            else:
                custom_input = prompt_with_context(
                    self.processor,
                    "Subsample every Nth frame",
                    default="5",
                    module="MD Manager - Pairwise RMSD",
                    description="Custom subsampling interval"
                )
                try:
                    subsample = int(custom_input)
                except ValueError:
                    subsample = 5

        # Get region selection
        region = self._get_analysis_region_selection()

        # Perform analysis
        self.console.print(f"\n[cyan]Calculating pairwise RMSD matrix...[/cyan]")
        self.console.print(f"  Region: {region}")
        if subsample:
            self.console.print(f"  Subsampling: every {subsample}th frame")
        self.console.print("[yellow]This may take a while for large trajectories...[/yellow]")

        try:
            rmsd_data = analyzer.calculate_pairwise_rmsd(
                mask=region,
                subsample=subsample
            )

            self.console.print("[green]✓ Analysis complete[/green]")

            # Display results
            self._display_pairwise_rmsd_results(analyzer, rmsd_data)

        except Exception as e:
            self.console.print(f"[red]Error during pairwise RMSD analysis: {e}[/red]")
            logger.error(f"Pairwise RMSD failed: {e}", exc_info=True)

    def _display_pairwise_rmsd_results(self, analyzer, rmsd_data):
        """Display pairwise RMSD matrix results."""
        from rich.table import Table
        import numpy as np

        matrix = rmsd_data['matrix']
        avg_per_frame = rmsd_data['avg_per_frame']
        most_central = rmsd_data['most_central']
        central_rmsd = rmsd_data['central_rmsd']
        most_extreme = rmsd_data['most_extreme']
        extreme_rmsd = rmsd_data['extreme_rmsd']
        overall_mean = rmsd_data['overall_mean']
        overall_std = rmsd_data['overall_std']
        overall_max = rmsd_data['overall_max']
        mask = rmsd_data['mask']
        n_frames = rmsd_data['n_frames']

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]PAIRWISE RMSD RESULTS[/bold cyan]")
        self.console.print("="*70)

        # Summary statistics
        self.console.print(f"\n[bold]Overall Statistics:[/bold]")
        self.console.print(f"  Selection: {mask}")
        self.console.print(f"  Frames analyzed: {n_frames}")
        self.console.print(f"  Mean RMSD: {overall_mean:.2f} ± {overall_std:.2f} Å")
        self.console.print(f"  Maximum RMSD: {overall_max:.2f} Å")

        # Central and extreme frames
        self.console.print(f"\n[bold]Key Frames:[/bold]")
        self.console.print(f"  Most central frame: {most_central} (avg RMSD: {central_rmsd:.2f} Å)")
        self.console.print(f"  Most extreme frame: {most_extreme} (avg RMSD: {extreme_rmsd:.2f} Å)")

        # RMSD matrix heatmap (if not too large)
        if n_frames <= 80:
            self.console.print("\n" + "="*70)
            self.console.print("[bold]RMSD MATRIX HEATMAP[/bold]")
            self.console.print("="*70)

            self._display_rmsd_matrix_heatmap(matrix)
        else:
            self.console.print(f"\n[grey50]Matrix too large for full heatmap display ({n_frames}x{n_frames})[/grey50]")
            self.console.print("[grey50]Showing downsampled heatmap...[/grey50]")

            # Downsample matrix for display
            from scipy import ndimage
            max_display = 60
            zoom_factor = max_display / n_frames
            downsampled = ndimage.zoom(matrix, (zoom_factor, zoom_factor), order=0)

            self.console.print("\n" + "="*70)
            self.console.print("[bold]RMSD MATRIX HEATMAP (downsampled)[/bold]")
            self.console.print("="*70)

            self._display_rmsd_matrix_heatmap(downsampled)

        # Average RMSD per frame plot
        self.console.print("\n" + "="*70)
        self.console.print("[bold]AVERAGE RMSD PER FRAME[/bold]")
        self.console.print("="*70)

        # Create ASCII plot
        frame_indices = list(range(n_frames))
        avg_plot = self._create_ascii_plot(
            avg_per_frame,
            title="Avg RMSD to all other frames",
            xlabel="Frame",
            ylabel="Avg RMSD (Å)",
            x_values=frame_indices
        )

        self.console.print(avg_plot, highlight=False)

        self.console.print("\n[grey50]Lower values = more representative/central conformations[/grey50]")

        # Distribution of pairwise RMSDs
        self.console.print("\n" + "="*70)
        self.console.print("[bold]PAIRWISE RMSD DISTRIBUTION[/bold]")
        self.console.print("="*70)

        # Flatten upper triangle
        upper_triangle = matrix[np.triu_indices(n_frames, k=1)]

        rmsd_hist = self._create_ascii_histogram(
            upper_triangle,
            title="Distribution of Pairwise RMSDs",
            xlabel="RMSD (Å)",
            bins=20
        )

        self.console.print(rmsd_hist)

    def _display_rmsd_matrix_heatmap(self, matrix):
        """Display RMSD matrix as ASCII heatmap."""
        import numpy as np

        n = matrix.shape[0]

        # Characters for RMSD levels
        chars = ' ░▒▓█'

        # Get max RMSD for scaling
        max_rmsd = np.max(matrix)

        self.console.print(f"\n[grey50]Frame j →[/grey50]")

        # Display heatmap
        for i in range(min(n, 50)):  # Limit rows
            line = ""
            for j in range(min(n, 60)):  # Limit columns
                rmsd = matrix[i, j]
                # Scale to character range
                char_idx = min(int((rmsd / max_rmsd) * (len(chars)-1)), len(chars)-1)
                line += chars[char_idx]

            if i == 0:
                self.console.print(f"Frame i ↓ │{line}")
            else:
                self.console.print(f"          │{line}")

        self.console.print(f"\n[grey50]Legend: {chars[0]}=0Å (identical)  {chars[-1]}={max_rmsd:.1f}Å (max)[/grey50]")
        self.console.print("[grey50]Diagonal = self-comparison (0 Å)[/grey50]")

    def _analyze_bfactors(self, analyzer):
        """Perform B-factor (pseudo B-factor from MD) analysis."""

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]B-FACTOR ANALYSIS[/bold cyan]")
        self.console.print("="*70)

        self.console.print("\n[bold]B-factor Overview:[/bold]")
        self.console.print("Calculates pseudo B-factors from atomic fluctuations (RMSF).")
        self.console.print("Formula: B = (8π²/3) × RMSF²")
        self.console.print("Can be compared to crystallographic B-factors.")

        # Get region selection
        region = self._get_analysis_region_selection()

        # Perform analysis
        self.console.print(f"\n[cyan]Calculating B-factors...[/cyan]")
        self.console.print(f"  Region: {region}")

        try:
            bfactor_data = analyzer.calculate_bfactors(
                mask=region,
                by_residue=True
            )

            self.console.print("[green]✓ Analysis complete[/green]")

            # Display results
            self._display_bfactor_results(analyzer, bfactor_data)

        except Exception as e:
            self.console.print(f"[red]Error during B-factor analysis: {e}[/red]")
            logger.error(f"B-factor analysis failed: {e}", exc_info=True)

    def _display_bfactor_results(self, analyzer, bfactor_data):
        """Display B-factor analysis results."""
        from rich.table import Table
        import numpy as np

        bfactors = bfactor_data['bfactors']
        rmsf = bfactor_data['rmsf']
        residue_ids = bfactor_data['residue_ids']
        residue_names = bfactor_data['residue_names']
        mask = bfactor_data['mask']

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]B-FACTOR RESULTS[/bold cyan]")
        self.console.print("="*70)

        # Statistics
        mean_bf = np.mean(bfactors)
        std_bf = np.std(bfactors)
        min_bf = np.min(bfactors)
        max_bf = np.max(bfactors)

        self.console.print(f"\n[bold]Overall Statistics:[/bold]")
        self.console.print(f"  Selection: {mask}")
        self.console.print(f"  Residues: {len(residue_ids)}")
        self.console.print(f"  Mean B-factor: {mean_bf:.2f} ± {std_bf:.2f} Å²")
        self.console.print(f"  Range: {min_bf:.2f} - {max_bf:.2f} Å²")

        # B-factor vs residue plot
        self.console.print("\n" + "="*70)
        self.console.print("[bold]B-FACTORS VS RESIDUE[/bold]")
        self.console.print("="*70)

        bf_plot = self._create_ascii_plot(
            bfactors,
            title="Per-Residue B-factors",
            xlabel="Residue",
            ylabel="B-factor (Å²)",
            x_values=residue_ids
        )

        self.console.print(bf_plot, highlight=False)

        # High B-factor residues
        threshold = mean_bf + std_bf
        high_bf_indices = np.where(bfactors > threshold)[0]

        if len(high_bf_indices) > 0:
            self.console.print("\n" + "="*70)
            self.console.print(f"[bold]HIGH B-FACTOR RESIDUES (> {threshold:.2f} Å²)[/bold]")
            self.console.print("="*70)

            bf_table = Table(show_header=True, header_style="bold bright_blue")
            bf_table.add_column("Residue", style="bright_blue")
            bf_table.add_column("B-factor (Å²)", style="yellow")
            bf_table.add_column("RMSF (Å)", style="green")
            bf_table.add_column("Deviation", style="grey50")

            # Sort by B-factor (highest first)
            sorted_indices = sorted(high_bf_indices, key=lambda i: bfactors[i], reverse=True)

            for idx in sorted_indices[:15]:  # Show top 15
                resid = residue_ids[idx]
                resname = residue_names[idx]
                bf = bfactors[idx]
                rmsf_val = rmsf[idx]
                deviation = (bf - mean_bf) / std_bf

                bf_table.add_row(
                    f"{resname}{resid}",
                    f"{bf:.2f}",
                    f"{rmsf_val:.2f}",
                    f"{deviation:.1f}σ"
                )

            self.console.print(bf_table)

            if len(high_bf_indices) > 15:
                self.console.print(f"\n[grey50]... and {len(high_bf_indices)-15} more high B-factor residues[/grey50]")

        # Low B-factor residues (rigid regions)
        low_threshold = mean_bf - std_bf
        low_bf_indices = np.where(bfactors < low_threshold)[0]

        if len(low_bf_indices) > 0:
            self.console.print("\n" + "="*70)
            self.console.print(f"[bold]LOW B-FACTOR RESIDUES (< {low_threshold:.2f} Å²)[/bold]")
            self.console.print("="*70)
            self.console.print("[grey50]These residues are relatively rigid[/grey50]\n")

            low_bf_table = Table(show_header=True, header_style="bold bright_blue")
            low_bf_table.add_column("Residue", style="bright_blue")
            low_bf_table.add_column("B-factor (Å²)", style="green")
            low_bf_table.add_column("RMSF (Å)", style="white")

            # Sort by B-factor (lowest first)
            sorted_low = sorted(low_bf_indices, key=lambda i: bfactors[i])

            for idx in sorted_low[:10]:  # Show top 10 most rigid
                resid = residue_ids[idx]
                resname = residue_names[idx]
                bf = bfactors[idx]
                rmsf_val = rmsf[idx]

                low_bf_table.add_row(
                    f"{resname}{resid}",
                    f"{bf:.2f}",
                    f"{rmsf_val:.2f}"
                )

            self.console.print(low_bf_table)

        # B-factor distribution
        self.console.print("\n" + "="*70)
        self.console.print("[bold]B-FACTOR DISTRIBUTION[/bold]")
        self.console.print("="*70)

        bf_hist = self._create_ascii_histogram(
            bfactors,
            title="Distribution of B-factors",
            xlabel="B-factor (Å²)",
            bins=20
        )

        self.console.print(bf_hist)

        # Interpretation guide
        self.console.print("\n" + "="*70)
        self.console.print("[bold]INTERPRETATION GUIDE[/bold]")
        self.console.print("="*70)
        self.console.print("\n[bold]B-factor ranges (general guidelines):[/bold]")
        self.console.print("  • < 20 Å²:  Very rigid (core, secondary structure)")
        self.console.print("  • 20-40 Å²: Moderate flexibility")
        self.console.print("  • 40-60 Å²: High flexibility (loops, termini)")
        self.console.print("  • > 60 Å²:  Very flexible (disordered regions)")
        self.console.print("\n[grey50]Note: These are typical ranges; interpretation depends on system[/grey50]")

    def _analyze_water_shells(self, analyzer):
        """Perform water shell (hydration layer) analysis."""

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]WATER SHELL ANALYSIS[/bold cyan]")
        self.console.print("="*70)

        self.console.print("\n[bold]Water Shell Overview:[/bold]")
        self.console.print("Analyzes hydration layers around solute (protein/molecule).")
        self.console.print("Counts water molecules in concentric shells at different distances.")

        # Get solute selection
        self.console.print("\n[bold]Solute Selection:[/bold]")
        self.console.print("  1. Protein (default)")
        self.console.print("  2. Custom selection")

        solute_choice = prompt_with_context(
            self.processor,
            "Select solute",
            choices=["1", "2"],
            default="1",
            module="MD Manager - Water Shells",
            description="Select solute for hydration analysis",
            options_map={
                "1": "Protein",
                "2": "Custom selection"
            }
        )

        if solute_choice == "2":
            solute_mask = prompt_with_context(
                self.processor,
                "Solute mask (AMBER syntax)",
                default=":1-100",
                module="MD Manager - Water Shells",
                description="Custom solute selection"
            )
        else:
            solute_mask = None  # Will use default

        # Get shell parameters
        self.console.print("\n[bold]Shell Parameters:[/bold]")
        self.console.print("  Shell width: Thickness of each hydration layer")
        self.console.print("  Max distance: How far from solute to analyze")

        shell_width_input = prompt_with_context(
            self.processor,
            "Shell width (Å)",
            default="2.0",
            module="MD Manager - Water Shells",
            description="Shell width in Angstroms"
        )

        max_distance_input = prompt_with_context(
            self.processor,
            "Maximum distance (Å)",
            default="10.0",
            module="MD Manager - Water Shells",
            description="Maximum distance to analyze"
        )

        try:
            shell_width = float(shell_width_input)
            max_distance = float(max_distance_input)
        except ValueError:
            shell_width = 2.0
            max_distance = 10.0

        # Perform analysis
        self.console.print(f"\n[cyan]Analyzing water shells...[/cyan]")
        self.console.print(f"  Shell width: {shell_width} Å")
        self.console.print(f"  Max distance: {max_distance} Å")

        try:
            shell_data = analyzer.calculate_water_shells(
                solute_mask=solute_mask,
                shell_width=shell_width,
                max_distance=max_distance
            )

            self.console.print("[green]✓ Analysis complete[/green]")

            # Display results
            self._display_water_shell_results(analyzer, shell_data)

        except Exception as e:
            self.console.print(f"[red]Error during water shell analysis: {e}[/red]")
            logger.error(f"Water shell analysis failed: {e}", exc_info=True)

    def _display_water_shell_results(self, analyzer, shell_data):
        """Display water shell analysis results."""
        from rich.table import Table
        import numpy as np

        populations = shell_data['populations']
        stats = shell_data['stats']
        shell_boundaries = shell_data['shell_boundaries']
        solute_mask = shell_data['solute_mask']
        n_shells = shell_data['n_shells']

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]WATER SHELL RESULTS[/bold cyan]")
        self.console.print("="*70)

        # Summary
        self.console.print(f"\n[bold]Analysis Summary:[/bold]")
        self.console.print(f"  Solute: {solute_mask}")
        self.console.print(f"  Number of shells: {n_shells}")
        self.console.print(f"  Frames analyzed: {populations.shape[0]}")

        # Shell statistics table
        self.console.print("\n" + "="*70)
        self.console.print("[bold]HYDRATION SHELL STATISTICS[/bold]")
        self.console.print("="*70)

        shell_table = Table(show_header=True, header_style="bold bright_blue")
        shell_table.add_column("Shell", style="bright_blue")
        shell_table.add_column("Distance (Å)", style="white")
        shell_table.add_column("Mean Waters", style="green")
        shell_table.add_column("Std", style="grey50")
        shell_table.add_column("Range", style="yellow")

        for i, stat in enumerate(stats):
            shell_table.add_row(
                str(i+1),
                stat['range'],
                f"{stat['mean']:.1f}",
                f"{stat['std']:.1f}",
                f"{stat['min']}-{stat['max']}"
            )

        self.console.print(shell_table)

        # Key insights
        first_shell_mean = stats[0]['mean']
        self.console.print(f"\n[bold]Key Insights:[/bold]")
        self.console.print(f"  • First hydration shell: {first_shell_mean:.1f} waters on average")
        self.console.print(f"  • First shell range: {stats[0]['range']}")

        # Average shell populations (bar chart)
        self.console.print("\n" + "="*70)
        self.console.print("[bold]AVERAGE WATER POPULATION BY SHELL[/bold]")
        self.console.print("="*70)

        max_mean = max(stat['mean'] for stat in stats)
        chart_width = 40

        for i, stat in enumerate(stats):
            mean_val = stat['mean']
            if max_mean > 0:
                bar_len = int((mean_val / max_mean) * chart_width)
            else:
                bar_len = 0

            bar = "█" * bar_len

            self.console.print(f"  Shell {i+1} ({stat['range']}) │{bar} {mean_val:.1f} waters")

        # First shell population over time
        self.console.print("\n" + "="*70)
        self.console.print("[bold]FIRST SHELL POPULATION VS TIME[/bold]")
        self.console.print("="*70)

        first_shell_pop = populations[:, 0]
        frame_times = [analyzer.traj.time[i] for i in range(len(first_shell_pop))]

        shell_plot = self._create_ascii_plot(
            first_shell_pop,
            title="First Hydration Shell (0.0-2.0 Å)",
            xlabel="Time (ps)",
            ylabel="# Waters",
            x_values=frame_times
        )

        self.console.print(shell_plot, highlight=False)

        # Shell population distribution for first shell
        self.console.print("\n" + "="*70)
        self.console.print("[bold]FIRST SHELL POPULATION DISTRIBUTION[/bold]")
        self.console.print("="*70)

        first_shell_hist = self._create_ascii_histogram(
            first_shell_pop,
            title="First Shell Water Count",
            xlabel="# Waters",
            bins=15
        )

        self.console.print(first_shell_hist)

    def _analyze_density_maps(self, analyzer):
        """Perform density map analysis."""

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]DENSITY MAP ANALYSIS[/bold cyan]")
        self.console.print("="*70)

        self.console.print("\n[bold]Density Map Overview:[/bold]")
        self.console.print("Creates 3D spatial density distribution.")
        self.console.print("Shows where selected atoms spend most time.")
        self.console.print("Useful for: water sites, ion binding, ligand positions")

        # Get selection
        self.console.print("\n[bold]Select atoms for density map:[/bold]")
        self.console.print("Examples:")
        self.console.print("  :WAT@O       - Water oxygens")
        self.console.print("  :NA+,CL-     - Ions")
        self.console.print("  :LIG         - Ligand")

        selection = prompt_with_context(
            self.processor,
            "Selection (AMBER mask)",
            default=":WAT@O",
            module="MD Manager - Density Maps",
            description="Atom selection for density map"
        )

        # Get grid spacing
        self.console.print("\n[bold]Grid Spacing:[/bold]")
        self.console.print("  Fine: 0.5 Å (high resolution, more memory)")
        self.console.print("  Medium: 1.0 Å (balanced)")
        self.console.print("  Coarse: 2.0 Å (fast, less detail)")

        spacing_input = prompt_with_context(
            self.processor,
            "Grid spacing (Å)",
            default="1.0",
            module="MD Manager - Density Maps",
            description="Grid spacing for density calculation"
        )

        try:
            grid_spacing = float(spacing_input)
        except ValueError:
            grid_spacing = 1.0

        # Perform analysis
        self.console.print(f"\n[cyan]Calculating density map...[/cyan]")
        self.console.print(f"  Selection: {selection}")
        self.console.print(f"  Grid spacing: {grid_spacing} Å")

        try:
            density_data = analyzer.calculate_density_map(
                selection_mask=selection,
                grid_spacing=grid_spacing
            )

            self.console.print("[green]✓ Analysis complete[/green]")

            # Display results
            self._display_density_results(analyzer, density_data)

        except Exception as e:
            self.console.print(f"[red]Error during density map analysis: {e}[/red]")
            logger.error(f"Density map analysis failed: {e}", exc_info=True)

    def _display_density_results(self, analyzer, density_data):
        """Display density map results."""
        from rich.table import Table
        import numpy as np

        density_slice = density_data['density_slice']
        grid_dims = density_data['grid_dims']
        grid_min = density_data['grid_min']
        grid_max = density_data['grid_max']
        grid_spacing = density_data['grid_spacing']
        peaks = density_data['peaks']
        max_density = density_data['max_density']
        selection = density_data['selection']

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]DENSITY MAP RESULTS[/bold cyan]")
        self.console.print("="*70)

        # Summary
        self.console.print(f"\n[bold]Grid Information:[/bold]")
        self.console.print(f"  Selection: {selection}")
        self.console.print(f"  Grid dimensions: {grid_dims[0]} × {grid_dims[1]} × {grid_dims[2]}")
        self.console.print(f"  Grid spacing: {grid_spacing:.2f} Å")
        self.console.print(f"  Grid bounds:")
        self.console.print(f"    X: {grid_min[0]:.1f} to {grid_max[0]:.1f} Å")
        self.console.print(f"    Y: {grid_min[1]:.1f} to {grid_max[1]:.1f} Å")
        self.console.print(f"    Z: {grid_min[2]:.1f} to {grid_max[2]:.1f} Å")

        # Density statistics
        self.console.print(f"\n[bold]Density Statistics:[/bold]")
        self.console.print(f"  Maximum density: {max_density:.4f}")

        # Peak locations
        if peaks:
            self.console.print("\n" + "="*70)
            self.console.print("[bold]PEAK DENSITY LOCATIONS[/bold]")
            self.console.print("="*70)

            peak_table = Table(show_header=True, header_style="bold bright_blue")
            peak_table.add_column("Rank", style="bright_blue")
            peak_table.add_column("X (Å)", style="white")
            peak_table.add_column("Y (Å)", style="white")
            peak_table.add_column("Z (Å)", style="white")
            peak_table.add_column("Density", style="green")

            for i, peak in enumerate(peaks, 1):
                x, y, z = peak['coords']
                dens = peak['density']

                peak_table.add_row(
                    str(i),
                    f"{x:.1f}",
                    f"{y:.1f}",
                    f"{z:.1f}",
                    f"{dens:.4f}"
                )

            self.console.print(peak_table)

        # 2D slice visualization
        self.console.print("\n" + "="*70)
        self.console.print("[bold]DENSITY SLICE (XY plane at mid-Z)[/bold]")
        self.console.print("="*70)

        # Downsample if needed for display
        max_display = 60
        if density_slice.shape[0] > max_display or density_slice.shape[1] > max_display:
            from scipy import ndimage
            zoom_x = min(1.0, max_display / density_slice.shape[0])
            zoom_y = min(1.0, max_display / density_slice.shape[1])
            density_slice = ndimage.zoom(density_slice, (zoom_x, zoom_y), order=0)

        # Create ASCII heatmap
        chars = ' ░▒▓█'
        height, width = density_slice.shape

        self.console.print(f"\n[grey50]Y ↑[/grey50]")

        for j in range(height-1, -1, -1):
            line = ""
            for i in range(width):
                dens = density_slice[i, j]

                if max_density > 0:
                    normalized = dens / max_density
                    char_idx = min(int(normalized * (len(chars)-1)), len(chars)-1)
                else:
                    char_idx = 0

                line += chars[char_idx]

            if j == height // 2:
                self.console.print(f"   │{line}")
            else:
                self.console.print(f"   │{line}")

        # X axis
        x_axis = "─" * width
        self.console.print(f"   └{x_axis}→ X")

        self.console.print(f"\n[grey50]Legend: {chars[0]}=no density  {chars[-1]}=max density[/grey50]")
        self.console.print(f"[grey50]Slice shown at Z = {(grid_min[2] + grid_max[2])/2:.1f} Å[/grey50]")

    def _analyze_vector(self, analyzer):
        """Perform vector orientation analysis."""

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]VECTOR ORIENTATION ANALYSIS[/bold cyan]")
        self.console.print("="*70)

        self.console.print("\n[bold]Vector Analysis Overview:[/bold]")
        self.console.print("Tracks orientation of molecular vectors over time.")
        self.console.print("Examples: helix tilt, bond rotation, molecular axis")

        # Get vector definition
        self.console.print("\n[bold]Define Vector:[/bold]")
        self.console.print("Vector = (Selection 2) - (Selection 1)")
        self.console.print("Example: Helix axis from N-terminal to C-terminal")

        atom1 = prompt_with_context(
            self.processor,
            "Selection 1 (vector start)",
            default=":1-10@CA",
            module="MD Manager - Vector",
            description="Vector start selection"
        )

        atom2 = prompt_with_context(
            self.processor,
            "Selection 2 (vector end)",
            default=":40-50@CA",
            module="MD Manager - Vector",
            description="Vector end selection"
        )

        # Perform analysis
        self.console.print(f"\n[cyan]Analyzing vector orientation...[/cyan]")

        try:
            vector_data = analyzer.calculate_vector_analysis(
                atom1_mask=atom1,
                atom2_mask=atom2
            )

            self.console.print("[green]✓ Analysis complete[/green]")

            # Display results
            self._display_vector_results(analyzer, vector_data)

        except Exception as e:
            self.console.print(f"[red]Error during vector analysis: {e}[/red]")
            logger.error(f"Vector analysis failed: {e}", exc_info=True)

    def _display_vector_results(self, analyzer, vector_data):
        """Display vector analysis results."""
        import numpy as np

        magnitudes = vector_data['magnitudes']
        polar_angles = vector_data['polar_angles']
        azimuthal_angles = vector_data['azimuthal_angles']

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]VECTOR RESULTS[/bold cyan]")
        self.console.print("="*70)

        # Statistics
        self.console.print(f"\n[bold]Vector Magnitude:[/bold]")
        self.console.print(f"  Mean: {np.mean(magnitudes):.2f} ± {np.std(magnitudes):.2f} Å")
        self.console.print(f"  Range: {np.min(magnitudes):.2f} - {np.max(magnitudes):.2f} Å")

        self.console.print(f"\n[bold]Polar Angle (from Z-axis):[/bold]")
        self.console.print(f"  Mean: {np.mean(polar_angles):.1f} ± {np.std(polar_angles):.1f}°")

        # Plot polar angle over time
        frame_times = [analyzer.traj.time[i] for i in range(len(polar_angles))]
        polar_plot = self._create_ascii_plot(
            polar_angles,
            title="Polar Angle vs Time",
            xlabel="Time (ps)",
            ylabel="Angle (°)",
            x_values=frame_times
        )

        self.console.print("\n" + "="*70)
        self.console.print(polar_plot, highlight=False)

    def _analyze_autocorrelation(self, analyzer):
        """Perform autocorrelation analysis on RMSD."""
        self.console.print("\n[cyan]Note: Autocorrelation requires pre-calculated data.[/cyan]")
        self.console.print("[cyan]Calculating RMSD first...[/cyan]")

        # Calculate RMSD for autocorrelation
        region = self._get_analysis_region_selection(analyzer, "Autocorrelation")
        rmsd = pt.rmsd(analyzer.traj, mask=region)

        # Calculate autocorrelation
        autocorr_data = analyzer.calculate_autocorrelation(rmsd)

        # Display
        autocorr = autocorr_data['autocorr']
        lag_times = list(range(len(autocorr)))

        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]AUTOCORRELATION RESULTS[/bold cyan]")
        self.console.print("="*70)

        autocorr_plot = self._create_ascii_plot(
            autocorr,
            title="Autocorrelation Function",
            xlabel="Lag (frames)",
            ylabel="Correlation",
            x_values=lag_times
        )

        self.console.print(autocorr_plot, highlight=False)

    def _analyze_contact_frequency_per_residue(self, analyzer):
        """Perform per-residue contact frequency analysis."""
        from rich.table import Table

        region = self._get_analysis_region_selection(analyzer, "Contact Frequency")

        self.console.print(f"\n[cyan]Calculating contact frequencies...[/cyan]")

        freq_data = analyzer.calculate_contact_frequency(mask=region)

        # Display top residues
        self.console.print("\n" + "="*70)
        self.console.print("[bold cyan]PER-RESIDUE CONTACT FREQUENCY[/bold cyan]")
        self.console.print("="*70)

        freq_table = Table(show_header=True, header_style="bold bright_blue")
        freq_table.add_column("Rank", style="bright_blue")
        freq_table.add_column("Residue", style="white")
        freq_table.add_column("Avg Contacts", style="green")

        for i, (resid, freq) in enumerate(zip(freq_data['residues'][:15], freq_data['frequencies'][:15]), 1):
            freq_table.add_row(str(i), str(resid), f"{freq:.1f}")

        self.console.print(freq_table)

    def _offer_trajectory_export(self, analyzer, sim_name: str):
        """Offer to export trajectory analysis data."""
        self.console.print("\n" + "="*70)
        self.console.print("EXPORT DATA")
        self.console.print("="*70)

        self.console.print("\n[bold]Export Options:[/bold]")
        self.console.print("  1. CSV format (for Excel, Origin, plotting)")
        self.console.print("  2. JSON format (for Python, R)")
        self.console.print("  3. Both formats")
        self.console.print("  4. Don't export")

        export_choice = prompt_with_context(
            self.processor,
            "Export analysis results?",
            choices=["1", "2", "3", "4"],
            default="4",
            module="MD Manager - Export",
            description="Export trajectory data",
            options_map={
                "1": "CSV format (for Excel, Origin, plotting)",
                "2": "JSON format (for Python, R)",
                "3": "Both formats",
                "4": "Don't export"
            }
        )

        if export_choice == "4":
            return

        # Determine output directory
        output_dir = self.working_directory
        base_name = f"{sim_name}_trajectory_analysis"

        formats_to_export = []
        if export_choice == "1":
            formats_to_export = ['csv']
        elif export_choice == "2":
            formats_to_export = ['json']
        elif export_choice == "3":
            formats_to_export = ['csv', 'json']

        # Export each format
        for fmt in formats_to_export:
            output_file = output_dir / f"{base_name}.{fmt}"
            try:
                exported_path = analyzer.export_data(
                    output_format=fmt,
                    output_file=output_file
                )
                self.console.print(f"[green]✓ Exported to {exported_path}[/green]")
            except Exception as e:
                self.console.print(f"[red]Error exporting {fmt}: {e}[/red]")