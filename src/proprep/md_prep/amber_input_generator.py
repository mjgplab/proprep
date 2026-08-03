"""
AMBER Input Generator Module - Comprehensive Implementation

Generates optimized mdin files for sander or pmemd based on simulation requirements.
Guides users through parameter selection with educational explanations.
Fixed all identified issues and implemented missing methods.
"""

import os
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import json
import copy
import importlib.util
from datetime import datetime
import glob

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
    float_prompt_with_context,
)
from proprep.utils.file_browser import remap_recorded_index, annotate_selected_path

from proprep.utils.paths import get_package_dir
from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.md_prep.amber_controller import AmberController


# Temporarily disabled - replaced by MolecularDynamicsManager
# @register_module
class AmberInputGenerator(ProcessingModule):
    """
    Generate AMBER MD input files using unified template-wizard system.

    This redesigned system provides:
    - Template-centric workflow (always start/end with annotated templates)
    - Choice of configuration method (wizard vs direct edit)
    - Individual simulation phase configuration
    - Complete workflow assembly
    - Custom template and preset management
    """

    NAME = "AMBER Input Generator"
    CATEGORY = "MD Preparation"
    DESCRIPTION = "Create optimized mdin files using annotated templates and comprehensive wizard"
    VERSION = "2.0.0"

    def __init__(self):
        super().__init__()

        # Package directory for template storage
        self.package_dir = get_package_dir() / "md_prep"
        
        # Initialize the unified controller
        self.controller = None
        
    def set_processor(self, processor):
        """Set the processor reference and initialize controller."""
        self.processor = processor
        
        # Initialize controller with processor context
        self.controller = AmberController(processor, self.package_dir)

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
        """Get redesigned menu options for simulation-phase-based workflow."""
        if self.controller:
            return self.controller.get_menu_options()
        else:
            # Fallback menu before controller is initialized
            return {
                "minimization": "Configure Minimization Settings",
                "heating": "Configure Heating/Thermal Equilibration", 
                "equilibration": "Configure Pressure Equilibration",
                "production": "Configure Production MD Settings",
                "assemble_workflow": "Assemble Complete Workflow",
                "load_workflow": "Load Workflow Preset",
                "save_workflow": "Save Current Workflow as Preset",
            }

    def get_workspace_requirements(self) -> List[str]:
        """No specific workspace requirements - can work standalone."""
        return []

    def can_process(self, workspace) -> bool:
        """Can always process."""
        return True

    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection using the unified controller."""
        try:
            if not self.controller:
                self.console.print("[red]Controller not initialized. Please set processor first.[/red]")
                return False
                
            return self.controller.handle_menu_option(option)
            
        except Exception as e:
            self.console.print(f"[red]Error executing option '{option}': {e}[/red]")
            return False
        """Ensure all required directories exist."""
        for directory in [
            self.user_data_dir,
            self.templates_dir,
            self.presets_dir,
            self.builtin_templates_dir,
            self.builtin_presets_dir,
        ]:
            directory.mkdir(parents=True, exist_ok=True)

    def _initialize_builtin_content(self):
        """Initialize builtin templates and presets."""
        self._create_builtin_templates()
        self._create_builtin_presets()

    def _create_builtin_templates(self):
        """Create default builtin templates."""
        templates = {
            "basic_md": {
                "name": "Basic MD Simulation",
                "description": "Standard molecular dynamics simulation",
                "simulation_type": "md",
                "config": {
                    "imin": 0,
                    "ntx": 5,
                    "irest": 1,
                    "nstlim": 2500000,
                    "dt": 0.002,
                    "ntpr": 1000,
                    "ntwx": 1000,
                    "ntwr": 10000,
                    "ntxo": 2,
                    "ioutfm": 1,
                    "ntt": 3,
                    "temp0": 300.0,
                    "gamma_ln": 2.0,
                    "ntp": 1,
                    "pres0": 1.0,
                    "barostat": 2,
                    "ntb": 2,
                    "ntc": 2,
                    "ntf": 2,
                    "cut": 10.0,
                }
            },
            "minimization": {
                "name": "Energy Minimization",
                "description": "Standard energy minimization",
                "simulation_type": "minimization",
                "config": {
                    "imin": 1,
                    "ntx": 1,
                    "maxcyc": 10000,
                    "ncyc": 1000,
                    "ntmin": 1,
                    "dx0": 0.01,
                    "drms": 0.0001,
                    "ntpr": 100,
                    "ntwr": 500,
                    "cut": 10.0,
                }
            },
            "heating": {
                "name": "Heating Protocol",
                "description": "Gradual heating from 0 to 300K",
                "simulation_type": "heating",
                "config": {
                    "imin": 0,
                    "ntx": 1,
                    "irest": 0,
                    "nstlim": 50000,
                    "dt": 0.001,
                    "ntpr": 500,
                    "ntwx": 500,
                    "ntwr": 2500,
                    "ntxo": 2,
                    "ioutfm": 1,
                    "ntt": 3,
                    "tempi": 0.0,
                    "temp0": 300.0,
                    "gamma_ln": 2.0,
                    "ntr": 1,
                    "restraint_wt": 10.0,
                    "restraintmask": "!@H=",
                    "ntb": 1,
                    "ntc": 2,
                    "ntf": 2,
                    "cut": 10.0,
                    "nmropt": 1,
                }
            }
        }
        
        for template_name, template_data in templates.items():
            template_file = self.builtin_templates_dir / f"{template_name}.json"
            if not template_file.exists():
                with open(template_file, 'w') as f:
                    json.dump(template_data, f, indent=2)

    def _create_builtin_presets(self):
        """Create builtin workflow presets."""
        presets = {
            "standard_protein": {
                "name": "Standard Protein MD Workflow",
                "description": "Complete protein MD simulation workflow",
                "steps": {
                    "minimization": {
                        "description": "Energy minimization",
                        "filename": "01_min.in",
                        "config": {
                            "imin": 1,
                            "ntx": 1,
                            "maxcyc": 10000,
                            "ncyc": 1000,
                            "ntmin": 1,
                            "ntpr": 100,
                            "ntwr": 500,
                            "cut": 10.0,
                            "ntr": 1,
                            "restraint_wt": 500.0,
                            "restraintmask": "!@H=",
                        }
                    },
                    "heating": {
                        "description": "Heating from 0 to 300K",
                        "filename": "02_heat.in",
                        "config": {
                            "imin": 0,
                            "ntx": 1,
                            "irest": 0,
                            "nstlim": 50000,
                            "dt": 0.001,
                            "ntpr": 500,
                            "ntwx": 500,
                            "ntwr": 2500,
                            "ntxo": 2,
                            "ioutfm": 1,
                            "ntt": 3,
                            "tempi": 0.0,
                            "temp0": 300.0,
                            "gamma_ln": 2.0,
                            "ntr": 1,
                            "restraint_wt": 10.0,
                            "restraintmask": "!@H=",
                            "ntb": 1,
                            "ntc": 2,
                            "ntf": 2,
                            "cut": 10.0,
                            "nmropt": 1,
                        },
                        "nmr_section": "&wt type='TEMP0', istep1=0, istep2=25000,\n    value1=0.0, value2=300.0 /\n&wt type='END' /"
                    },
                    "equilibration": {
                        "description": "NPT equilibration",
                        "filename": "03_equil.in",
                        "config": {
                            "imin": 0,
                            "ntx": 5,
                            "irest": 1,
                            "nstlim": 100000,
                            "dt": 0.002,
                            "ntpr": 1000,
                            "ntwx": 1000,
                            "ntwr": 5000,
                            "ntxo": 2,
                            "ioutfm": 1,
                            "ntt": 3,
                            "temp0": 300.0,
                            "gamma_ln": 2.0,
                            "ntp": 1,
                            "pres0": 1.0,
                            "barostat": 2,
                            "ntb": 2,
                            "ntr": 1,
                            "restraint_wt": 1.0,
                            "restraintmask": "!@H=",
                            "ntc": 2,
                            "ntf": 2,
                            "cut": 10.0,
                        }
                    },
                    "production": {
                        "description": "Production MD simulation",
                        "filename": "04_prod.in",
                        "config": {
                            "imin": 0,
                            "ntx": 5,
                            "irest": 1,
                            "nstlim": 2500000,
                            "dt": 0.002,
                            "ntpr": 1000,
                            "ntwx": 1000,
                            "ntwr": 10000,
                            "ntxo": 2,
                            "ioutfm": 1,
                            "ntt": 3,
                            "temp0": 300.0,
                            "gamma_ln": 2.0,
                            "ntp": 1,
                            "pres0": 1.0,
                            "barostat": 2,
                            "ntb": 2,
                            "ntc": 2,
                            "ntf": 2,
                            "cut": 10.0,
                        }
                    }
                }
            }
        }
        
        for preset_name, preset_data in presets.items():
            preset_file = self.builtin_presets_dir / f"{preset_name}.json"
            if not preset_file.exists():
                with open(preset_file, 'w') as f:
                    json.dump(preset_data, f, indent=2)

    def _load_all_workflow_presets(self) -> Dict:
        """Load all available workflow presets."""
        presets = {}
        
        # Load builtin presets
        for preset_file in self.builtin_presets_dir.glob("*.json"):
            try:
                with open(preset_file, 'r') as f:
                    preset_data = json.load(f)
                    presets[preset_file.stem] = preset_data
            except Exception as e:
                self.console.print(f"[yellow]Warning: Could not load preset {preset_file}: {e}[/yellow]")
        
        # Load user presets
        for preset_file in self.presets_dir.glob("*.json"):
            if preset_file.parent != self.builtin_presets_dir:
                try:
                    with open(preset_file, 'r') as f:
                        preset_data = json.load(f)
                        presets[preset_file.stem] = preset_data
                except Exception as e:
                    self.console.print(f"[yellow]Warning: Could not load preset {preset_file}: {e}[/yellow]")
        
        return presets

    # ============================================================================
    # MAIN WORKFLOW METHODS
    # ============================================================================

    def _generate_input_workflow(self) -> bool:
        """Main workflow for generating single input file - comprehensive 15-step process."""
        try:
            self.console.print("\n[bold cyan]===== AMBER Input Generator =====[/bold cyan]")
            self.console.print("Let's configure your MD simulation step by step.\n")

            # Step 1: Simulation Type
            if not self._select_simulation_type():
                return False

            # Step 2: Engine Selection
            if not self._select_engine():
                return False

            # Step 3: Performance Options
            if not self._configure_performance():
                return False

            # Step 4: System Setup (periodic vs non-periodic)
            if not self._configure_system_type():
                return False

            # Step 5: Input/Restart Configuration
            if not self._configure_input_restart():
                return False

            # Step 6: Simulation Parameters (context-dependent)
            if not self._configure_simulation_parameters():
                return False

            # Step 7: Time Control (for MD)
            if self.simulation_type != "minimization":
                if not self._configure_time_control():
                    return False

            # Step 8: Temperature Control
            if not self._configure_temperature():
                return False

            # Step 9: Pressure Control (if needed)
            if not self._configure_pressure():
                return False

            # Step 10: Constraints
            if not self._configure_constraints():
                return False

            # Step 11: Restraints (optional)
            if not self._configure_restraints():
                return False

            # Step 12: Non-bonded Interactions
            if not self._configure_nonbonded():
                return False

            # Step 13: Output Control
            if not self._configure_output():
                return False

            # Step 14: Advanced Options
            if not self._configure_advanced_options():
                return False

            # Step 15: Review and Generate
            return self._review_and_generate()

        except KeyboardInterrupt:
            self.console.print("\n[yellow]Configuration cancelled by user[/yellow]")
            return False
        except Exception as e:
            self.console.print(f"\n[red]Error in workflow: {e}[/red]")
            return False

    def _select_simulation_type(self) -> bool:
        """Step 1: Enhanced simulation type selection with all imin options."""
        self.console.print("[bold]Step 1: Simulation Type Selection[/bold]")

        options = {
            "1": ("minimization", "Energy Minimization (imin=1)"),
            "2": ("md", "Molecular Dynamics (imin=0)"),
            "3": ("trajectory_analysis", "Trajectory Analysis (imin=5)"),
            "4": ("md_analysis", "MD Analysis of Trajectory (imin=6)"),
            "5": ("socket_server", "Socket Server Mode (imin=7)"),
            "6": ("custom", "Custom Setup"),
        }

        for key, (_, desc) in options.items():
            self.console.print(f"[cyan]{key}.[/cyan] {desc}")

        sim_type_opts = {k: v[1] for k, v in options.items()}
        choice = prompt_with_context(
            self.processor,
            "Select simulation type",
            choices=list(options.keys()),
            default="2",
            module="AMBER Input Generator",
            description="Select AMBER simulation type",
            options_map=sim_type_opts,
        )
        
        self.simulation_type, description = options[choice]
        
        # Set imin based on choice
        if self.simulation_type == "minimization":
            self.config["imin"] = 1
        elif self.simulation_type == "trajectory_analysis":
            self.config["imin"] = 5
        elif self.simulation_type == "md_analysis":
            self.config["imin"] = 6
        elif self.simulation_type == "socket_server":
            self.config["imin"] = 7
        else:
            self.config["imin"] = 0

        self.console.print(f"[green]Selected: {description}[/green]")
        
        # Provide educational information
        if self.simulation_type == "trajectory_analysis":
            self.console.print("[cyan]ℹ️  Trajectory analysis reads structures and performs operations like energy calculation[/cyan]")
        elif self.simulation_type == "socket_server":
            self.console.print("[cyan]ℹ️  Socket server mode for use with external drivers like i-PI[/cyan]")
        elif self.simulation_type == "md":
            self.console.print("[cyan]ℹ️  MD includes heating, equilibration, production, and other dynamics[/cyan]")
        
        return True

    def _select_engine(self) -> bool:
        """Step 2: Engine Selection with CPU/GPU configuration."""
        self.console.print("\n[bold]Step 2: Engine Selection[/bold]")

        self.console.print("\n[bold cyan]2a.[/bold cyan] Select MD engine:")
        self.console.print("[cyan]1.[/cyan] sander (CPU-based, versatile)")
        self.console.print("[cyan]2.[/cyan] pmemd (optimized, CPU)")
        self.console.print("[cyan]3.[/cyan] pmemd.cuda (GPU-accelerated)")

        choice = prompt_with_context(
            self.processor,
            "Select MD engine",
            choices=["1", "2", "3"],
            default="2",
            module="AMBER Input Generator",
            description="Select MD engine",
            options_map={"1": "sander (CPU)", "2": "pmemd (optimized CPU)", "3": "pmemd.cuda (GPU)"},
        )

        engines = {"1": "sander", "2": "pmemd", "3": "pmemd.cuda"}
        self.engine_choice = engines[choice]
        
        # Configure compute resources based on engine choice
        if choice == "3":  # GPU
            self.config["_use_gpu"] = True
            self.console.print("\n[bold cyan]2b.[/bold cyan] GPU Configuration:")
            gpu_ids = prompt_with_context(
                self.processor,
                "GPU indices (comma-separated, e.g., 0,1)",
                default="0",
                module="AMBER Input Generator",
                description="GPU indices (comma-separated)",
            )
            self.config["_gpu_ids"] = gpu_ids
            self.console.print(f"[green]GPU acceleration enabled on GPU(s): {gpu_ids}[/green]")
            self.console.print("[cyan]ℹ️  Run script will set CUDA_VISIBLE_DEVICES[/cyan]")
        else:  # CPU
            self.config["_use_gpu"] = False
            self.console.print("\n[bold cyan]2b.[/bold cyan] CPU Configuration:")
            ncpus = int_prompt_with_context(
                self.processor,
                "Number of CPU cores",
                default=1,
                module="AMBER Input Generator",
                description="Number of CPU cores",
            )
            self.config["_ncpus"] = ncpus
            
            if ncpus > 1:
                self.config["_use_mpi"] = True
                self.console.print(f"[green]MPI parallelization enabled with {ncpus} cores[/green]")
                self.console.print("[cyan]ℹ️  Run script will use mpirun[/cyan]")
            else:
                self.config["_use_mpi"] = False

        self.console.print(f"\n[green]Selected: {self.engine_choice}[/green]")
        return True

    def _configure_performance(self) -> bool:
        """Step 3: Performance Options."""
        self.console.print("\n[bold]Step 3: Performance Options[/bold]")

        # NRESPA multiple time stepping
        self.console.print("\n[bold cyan]3a.[/bold cyan] Multiple time stepping:")
        use_respa = confirm_with_context(
            self.processor,
            "Use multiple time stepping (NRESPA)?",
            default=False,
            module="AMBER Input Generator",
            description="Use multiple time stepping (NRESPA)",
        )
        if use_respa:
            nrespa = int_prompt_with_context(
                self.processor,
                "NRESPA factor",
                default=2,
                module="AMBER Input Generator",
                description="NRESPA factor",
            )
            self.config["nrespa"] = nrespa
            self.console.print(f"[cyan]ℹ️  Slow forces evaluated every {nrespa} steps[/cyan]")

        # Center of mass motion removal
        if self.simulation_type not in ["minimization", "trajectory_analysis"]:
            self.console.print("\n[bold cyan]3b.[/bold cyan] Center-of-mass motion removal:")
            nscm = int_prompt_with_context(
                self.processor,
                "Center-of-mass motion removal frequency (nscm)",
                default=1000,
                module="AMBER Input Generator",
                description="NSCM: COM motion removal frequency",
            )
            self.config["nscm"] = nscm
            if nscm == 0:
                self.console.print("[yellow]⚠️  COM motion removal disabled[/yellow]")

        return True

    def _configure_system_type(self) -> bool:
        """Step 4: System Setup (periodic vs non-periodic)."""
        self.console.print("\n[bold]Step 4: System Setup[/bold]")

        self.console.print("[cyan]1.[/cyan] Vacuum/Gas phase (ntb=0)")
        self.console.print("[cyan]2.[/cyan] Periodic boundary conditions (ntb=1)")
        self.console.print("[cyan]3.[/cyan] Constant pressure periodic (ntb=2)")

        choice = prompt_with_context(
            self.processor,
            "Select boundary conditions",
            choices=["1", "2", "3"],
            default="2",
            module="AMBER Input Generator",
            description="Select boundary conditions",
            options_map={
                "1": "Vacuum/Gas phase (ntb=0)",
                "2": "Periodic boundary conditions (ntb=1)",
                "3": "Constant pressure periodic (ntb=2)",
            },
        )

        ntb_values = {"1": 0, "2": 1, "3": 2}
        self.config["ntb"] = ntb_values[choice]

        if choice == "1":
            self.console.print("[cyan]ℹ️  No periodic boundaries - suitable for gas phase or small molecules[/cyan]")
        else:
            self.console.print("[cyan]ℹ️  Periodic boundary conditions - typical for condensed phase[/cyan]")

        return True

    def _configure_input_restart(self) -> bool:
        """Step 5: Input/Restart Configuration."""
        self.console.print("\n[bold]Step 5: Input/Restart Configuration[/bold]")

        if self.simulation_type in ["md", "custom"] or self.config["imin"] == 0:
            self.console.print("\n[bold cyan]5a.[/bold cyan] Input file type:")
            self.console.print("[cyan]1.[/cyan] New simulation (ntx=1, irest=0)")
            self.console.print("[cyan]2.[/cyan] Restart from previous run (ntx=5, irest=1)")

            choice = prompt_with_context(
                self.processor,
                "Select input type",
                choices=["1", "2"],
                default="1",
                module="AMBER Input Generator",
                description="Input file type (new vs restart)",
                options_map={"1": "New simulation (ntx=1, irest=0)", "2": "Restart from previous run"},
            )

            if choice == "1":
                self.config["ntx"] = 1
                self.config["irest"] = 0
                self.console.print("[cyan]ℹ️  Starting new simulation - velocities will be generated[/cyan]")
            else:
                self.config["ntx"] = 5
                self.config["irest"] = 1
                self.console.print("[cyan]ℹ️  Restarting - coordinates and velocities read from restart file[/cyan]")
        else:
            # Minimization typically starts fresh
            self.config["ntx"] = 1
            self.config["irest"] = 0
            self.console.print("\n[bold cyan]5a.[/bold cyan] Input configuration:")
            self.console.print("[cyan]ℹ️  Minimization uses ntx=1, irest=0 (coordinates only)[/cyan]")

        return True

    def _configure_simulation_parameters(self) -> bool:
        """Step 6: Simulation-specific parameters."""
        self.console.print(f"\n[bold]Step 6: {self.simulation_type.title()} Settings[/bold]")

        if self.simulation_type == "minimization":
            return self._configure_minimization_complete()
        elif self.simulation_type in ["trajectory_analysis", "md_analysis"]:
            return self._configure_analysis_mode()
        elif self.simulation_type == "socket_server":
            return self._configure_socket_mode()
        else:
            return self._configure_dynamics_complete()

    def _configure_minimization_complete(self) -> bool:
        """Complete minimization configuration."""
        # Maximum cycles
        self.console.print("\n[bold cyan]6a.[/bold cyan] Minimization cycles:")
        maxcyc = int_prompt_with_context(
            self.processor,
            "Maximum minimization cycles (maxcyc)",
            default=10000,
            module="AMBER Input Generator",
            description="Max minimization cycles (maxcyc)",
        )
        self.config["maxcyc"] = maxcyc

        # Method selection
        self.console.print("\n[bold cyan]6b.[/bold cyan] Minimization method:")
        self.console.print("[cyan]0.[/cyan] Full conjugate gradient")
        self.console.print("[cyan]1.[/cyan] Steepest descent → conjugate gradient (recommended)")
        self.console.print("[cyan]2.[/cyan] Steepest descent only")
        self.console.print("[cyan]3.[/cyan] XMIN method")
        self.console.print("[cyan]4.[/cyan] LMOD method")

        ntmin = prompt_with_context(
            self.processor,
            "Select method",
            choices=["0", "1", "2", "3", "4"],
            default="1",
            module="AMBER Input Generator",
            description="Minimization method (ntmin)",
            options_map={
                "0": "Full conjugate gradient",
                "1": "Steepest descent → conjugate gradient",
                "2": "Steepest descent only",
                "3": "XMIN method",
                "4": "LMOD method",
            },
        )
        self.config["ntmin"] = int(ntmin)

        if ntmin == "1":
            self.console.print("\n[bold cyan]6c.[/bold cyan] Steepest descent configuration:")
            ncyc = int_prompt_with_context(
                self.processor,
                "Steepest descent cycles before switching (ncyc)",
                default=1000,
                module="AMBER Input Generator",
                description="Steepest-descent cycles before switching (ncyc)",
            )
            self.config["ncyc"] = ncyc

        # Convergence criteria
        self.console.print("\n[bold cyan]6d.[/bold cyan] Convergence criteria:")
        drms = float_prompt_with_context(
            self.processor,
            "Convergence criterion (drms) [kcal/mol/Å]",
            default=0.0001,
            module="AMBER Input Generator",
            description="Convergence criterion drms [kcal/mol/Å]",
        )
        self.config["drms"] = drms

        # Initial step size
        self.console.print("\n[bold cyan]6e.[/bold cyan] Step size:")
        dx0 = float_prompt_with_context(
            self.processor,
            "Initial step size (dx0)",
            default=0.01,
            module="AMBER Input Generator",
            description="Initial minimization step size (dx0)",
        )
        self.config["dx0"] = dx0

        return True

    def _configure_analysis_mode(self) -> bool:
        """Configure trajectory analysis mode."""
        if self.simulation_type == "trajectory_analysis":
            self.console.print("[cyan]ℹ️  Trajectory analysis mode - specify trajectory file with -y flag[/cyan]")
            
            # For analysis, often want single point energies
            single_point = confirm_with_context(
                self.processor,
                "Single point energy calculation only?",
                default=True,
                module="AMBER Input Generator",
                description="Single-point energy calculation only",
            )
            if single_point:
                self.config["maxcyc"] = 1
            else:
                maxcyc = int_prompt_with_context(
                    self.processor,
                    "Maximum cycles per frame",
                    default=100,
                    module="AMBER Input Generator",
                    description="Max cycles per frame (trajectory analysis)",
                )
                self.config["maxcyc"] = maxcyc

        elif self.simulation_type == "md_analysis":
            self.console.print("[cyan]ℹ️  MD analysis mode - runs MD from each trajectory frame[/cyan]")

            nstlim = int_prompt_with_context(
                self.processor,
                "MD steps per frame (0 for single point)",
                default=0,
                module="AMBER Input Generator",
                description="MD steps per frame for MD-analysis mode",
            )
            self.config["nstlim"] = nstlim

            if nstlim > 0:
                dt = float_prompt_with_context(
                    self.processor,
                    "Time step [ps]",
                    default=0.002,
                    module="AMBER Input Generator",
                    description="MD time step [ps]",
                )
                self.config["dt"] = dt

        return True

    def _configure_socket_mode(self) -> bool:
        """Configure socket server mode."""
        self.console.print("[cyan]ℹ️  Socket server mode for external drivers (e.g., i-PI)[/cyan]")
        self.console.print("[cyan]ℹ️  Use -host and -port command line arguments to specify connection[/cyan]")
        
        # Socket mode doesn't need many parameters
        self.console.print("[green]Socket mode configured - minimal parameters needed[/green]")
        return True

    def _configure_dynamics_complete(self) -> bool:
        """Configure molecular dynamics parameters."""
        self.console.print("\n[bold cyan]6a.[/bold cyan] MD simulation setup:")
        self.console.print("[cyan]ℹ️  MD parameters will be configured in the following steps[/cyan]")
        return True

    def _configure_time_control(self) -> bool:
        """Step 7: Time Control."""
        self.console.print("\n[bold]Step 7: Time Control[/bold]")

        # Time step
        self.console.print("\n[bold cyan]7a.[/bold cyan] Time step configuration:")
        if self.simulation_type == "md":
            default_dt = 0.002
        else:
            default_dt = 0.001

        dt = float_prompt_with_context(
            self.processor,
            "Time step (dt) [ps]",
            default=default_dt,
            module="AMBER Input Generator",
            description="MD time step dt [ps]",
        )
        self.config["dt"] = dt

        if dt > 0.002:
            self.console.print("[yellow]⚠️  Large time step - ensure proper constraints[/yellow]")

        # Number of steps
        self.console.print("\n[bold cyan]7b.[/bold cyan] Simulation length:")
        if self.simulation_type == "md":
            default_steps = 2500000
        else:
            default_steps = 100000

        nstlim = int_prompt_with_context(
            self.processor,
            "Number of MD steps (nstlim)",
            default=default_steps,
            module="AMBER Input Generator",
            description="Number of MD steps (nstlim)",
        )
        self.config["nstlim"] = nstlim

        # Calculate and show runtime
        runtime_ps = nstlim * dt
        runtime_ns = runtime_ps / 1000
        self.console.print(f"[cyan]📊 Total simulation time: {runtime_ps:.1f} ps ({runtime_ns:.3f} ns)[/cyan]")

        # Start time (optional)
        self.console.print("\n[bold cyan]7c.[/bold cyan] Start time (optional):")
        if confirm_with_context(
            self.processor,
            "Set custom start time?",
            default=False,
            module="AMBER Input Generator",
            description="Set custom start time",
        ):
            t = float_prompt_with_context(
                self.processor,
                "Start time (t) [ps]",
                default=0.0,
                module="AMBER Input Generator",
                description="MD start time t [ps]",
            )
            self.config["t"] = t

        return True

    def _configure_temperature(self) -> bool:
        """Step 8: Temperature Control."""
        self.console.print("\n[bold]Step 8: Temperature Control[/bold]")

        if self.simulation_type == "minimization":
            self.console.print("[cyan]ℹ️  No temperature control needed for minimization[/cyan]")
            return True

        self.console.print("\n[bold cyan]8a.[/bold cyan] Thermostat selection:")
        self.console.print("[cyan]0.[/cyan] No temperature control (NVE)")
        self.console.print("[cyan]1.[/cyan] Weak coupling (not recommended)")
        self.console.print("[cyan]2.[/cyan] Andersen thermostat")
        self.console.print("[cyan]3.[/cyan] Langevin dynamics (recommended)")
        self.console.print("[cyan]11.[/cyan] Bussi thermostat")

        ntt = prompt_with_context(
            self.processor,
            "Select thermostat",
            choices=["0", "1", "2", "3", "11"],
            default="3",
            module="AMBER Input Generator",
            description="Thermostat selection (ntt)",
            options_map={
                "0": "No temperature control (NVE)",
                "1": "Weak coupling (not recommended)",
                "2": "Andersen",
                "3": "Langevin (recommended)",
                "11": "Bussi",
            },
        )
        self.config["ntt"] = int(ntt)

        if ntt != "0":
            # Target temperature
            self.console.print("\n[bold cyan]8b.[/bold cyan] Target temperature:")
            temp0 = float_prompt_with_context(
                self.processor,
                "Target temperature (temp0) [K]",
                default=300.0,
                module="AMBER Input Generator",
                description="Target temperature temp0 [K]",
            )
            self.config["temp0"] = temp0

            if temp0 > 350:
                self.console.print("[yellow]⚠️  High temperature - consider reducing time step[/yellow]")

            # Initial temperature for new simulations
            if self.config.get("irest", 0) == 0:
                self.console.print("\n[bold cyan]8c.[/bold cyan] Initial temperature:")
                tempi = float_prompt_with_context(
                    self.processor,
                    "Initial temperature (tempi) [K]",
                    default=temp0,
                    module="AMBER Input Generator",
                    description="Initial temperature tempi [K]",
                )
                self.config["tempi"] = tempi

            # Thermostat-specific parameters
            self.console.print("\n[bold cyan]8d.[/bold cyan] Thermostat parameters:")
            if ntt == "1":
                tautp = float_prompt_with_context(
                    self.processor,
                    "Heat bath time constant (tautp) [ps]",
                    default=1.0,
                    module="AMBER Input Generator",
                    description="Heat bath time constant tautp [ps] (weak-coupling)",
                )
                self.config["tautp"] = tautp
                self.console.print("[yellow]⚠️  Weak coupling can cause problems - consider Langevin[/yellow]")

            elif ntt == "2":
                vrand = int_prompt_with_context(
                    self.processor,
                    "Collision frequency (vrand steps)",
                    default=1000,
                    module="AMBER Input Generator",
                    description="Andersen collision frequency vrand [steps]",
                )
                self.config["vrand"] = vrand

            elif ntt == "3":
                gamma_ln = float_prompt_with_context(
                    self.processor,
                    "Collision frequency (gamma_ln) [ps⁻¹]",
                    default=2.0,
                    module="AMBER Input Generator",
                    description="Langevin collision frequency gamma_ln [ps⁻¹]",
                )
                self.config["gamma_ln"] = gamma_ln
                if gamma_ln < 0.1:
                    self.console.print("[cyan]ℹ️  Low collision frequency - good for enhanced sampling[/cyan]")

            elif ntt == "11":
                tautp = float_prompt_with_context(
                    self.processor,
                    "Thermostat time constant (tautp) [ps]",
                    default=1.0,
                    module="AMBER Input Generator",
                    description="Bussi thermostat time constant tautp [ps]",
                )
                self.config["tautp"] = tautp

            # Random seed
            self.console.print("\n[bold cyan]8e.[/bold cyan] Random number generation:")
            ig = int_prompt_with_context(
                self.processor,
                "Random seed (ig, -1 for time-based)",
                default=-1,
                module="AMBER Input Generator",
                description="Random seed ig (-1 = time-based)",
            )
            self.config["ig"] = ig

        return True

    def _configure_pressure(self) -> bool:
        """Step 9: Pressure Control."""
        self.console.print("\n[bold]Step 9: Pressure Control[/bold]")

        self.console.print("\n[bold cyan]9a.[/bold cyan] Pressure control setup:")
        if self.config.get("ntb", 1) < 2:
            use_pressure = confirm_with_context(
                self.processor,
                "Enable constant pressure (requires ntb=2)?",
                default=False,
                module="AMBER Input Generator",
                description="Enable constant pressure (sets ntb=2)",
            )
            if use_pressure:
                self.config["ntb"] = 2
            else:
                self.console.print("[cyan]ℹ️  Constant volume simulation[/cyan]")
                return True

        if self.config.get("ntb", 1) == 2 or self.simulation_type in ["md", "custom"]:
            self.console.print("\n[bold cyan]9b.[/bold cyan] Pressure scaling method:")
            self.console.print("[cyan]0.[/cyan] No pressure control")
            self.console.print("[cyan]1.[/cyan] Isotropic scaling")
            self.console.print("[cyan]2.[/cyan] Anisotropic scaling")
            self.console.print("[cyan]3.[/cyan] Semi-isotropic scaling")

            ntp = prompt_with_context(
                self.processor,
                "Select pressure control",
                choices=["0", "1", "2", "3"],
                default="1",
                module="AMBER Input Generator",
                description="Pressure scaling method (ntp)",
                options_map={
                    "0": "No pressure control",
                    "1": "Isotropic scaling",
                    "2": "Anisotropic scaling",
                    "3": "Semi-isotropic scaling",
                },
            )
            self.config["ntp"] = int(ntp)

            if ntp != "0":
                # Reference pressure
                self.console.print("\n[bold cyan]9c.[/bold cyan] Pressure parameters:")
                pres0 = float_prompt_with_context(
                    self.processor,
                    "Reference pressure (pres0) [bar]",
                    default=1.0,
                    module="AMBER Input Generator",
                    description="Reference pressure pres0 [bar]",
                )
                self.config["pres0"] = pres0

                # Barostat type
                self.console.print("\n[bold cyan]9d.[/bold cyan] Barostat selection:")
                self.console.print("[cyan]1.[/cyan] Berendsen")
                self.console.print("[cyan]2.[/cyan] Monte Carlo (recommended)")

                barostat = prompt_with_context(
                    self.processor,
                    "Select barostat",
                    choices=["1", "2"],
                    default="2",
                    module="AMBER Input Generator",
                    description="Select barostat",
                    options_map={"1": "Berendsen", "2": "Monte Carlo (recommended)"},
                )
                self.config["barostat"] = int(barostat)

                if barostat == "2":
                    self.console.print("\n[bold cyan]9e.[/bold cyan] Monte Carlo barostat settings:")
                    mcbarint = int_prompt_with_context(
                        self.processor,
                        "MC barostat interval",
                        default=100,
                        module="AMBER Input Generator",
                        description="Monte Carlo barostat interval (mcbarint)",
                    )
                    self.config["mcbarint"] = mcbarint

                # Pressure relaxation time
                self.console.print("\n[bold cyan]9f.[/bold cyan] Coupling parameters:")
                taup = float_prompt_with_context(
                    self.processor,
                    "Pressure relaxation time (taup) [ps]",
                    default=2.0,
                    module="AMBER Input Generator",
                    description="Pressure relaxation time taup [ps]",
                )
                self.config["taup"] = taup

                # Compressibility
                self.console.print("\n[bold cyan]9g.[/bold cyan] Compressibility:")
                comp = float_prompt_with_context(
                    self.processor,
                    "Compressibility [10⁻⁶ bar⁻¹]",
                    default=44.6,
                    module="AMBER Input Generator",
                    description="Compressibility [10⁻⁶ bar⁻¹]",
                )
                self.config["comp"] = comp

        return True

    def _configure_constraints(self) -> bool:
        """Step 10: Bond Constraints."""
        self.console.print("\n[bold]Step 10: Bond Constraints[/bold]")

        if self.simulation_type == "minimization":
            self.console.print("\n[bold cyan]10a.[/bold cyan] SHAKE during minimization:")
            use_shake = confirm_with_context(
                self.processor,
                "Use SHAKE during minimization?",
                default=False,
                module="AMBER Input Generator",
                description="Use SHAKE during minimization",
            )
            if not use_shake:
                self.console.print("[cyan]ℹ️  SHAKE typically not used for minimization[/cyan]")
                return True

        self.console.print("\n[bold cyan]10a.[/bold cyan] SHAKE constraint options:")
        self.console.print("[cyan]1.[/cyan] No constraints")
        self.console.print("[cyan]2.[/cyan] Constrain bonds to hydrogen (recommended)")
        self.console.print("[cyan]3.[/cyan] Constrain all bonds")

        choice = prompt_with_context(
            self.processor,
            "Select constraints",
            choices=["1", "2", "3"],
            default="2",
            module="AMBER Input Generator",
            description="SHAKE constraint option",
            options_map={
                "1": "No constraints",
                "2": "Constrain bonds to hydrogen (recommended)",
                "3": "Constrain all bonds",
            },
        )

        ntc_values = {"1": 1, "2": 2, "3": 3}
        self.config["ntc"] = ntc_values[choice]

        if choice != "1":
            # SHAKE tolerance
            self.console.print("\n[bold cyan]10b.[/bold cyan] SHAKE tolerance:")
            tol = float_prompt_with_context(
                self.processor,
                "SHAKE tolerance",
                default=0.00001,
                module="AMBER Input Generator",
                description="SHAKE tolerance (tol)",
            )
            self.config["tol"] = tol

            # Force evaluation - typically matches ntc
            if choice == "2":
                self.config["ntf"] = 2  # Skip H-bond forces
                self.console.print("[cyan]ℹ️  Hydrogen bond forces will be omitted (ntf=2)[/cyan]")
            elif choice == "3":
                self.config["ntf"] = 3  # Skip all bond forces
                self.console.print("[cyan]ℹ️  All bond forces will be omitted (ntf=3)[/cyan]")

            # Water model handling
            if choice == "2":
                self.console.print("\n[bold cyan]10c.[/bold cyan] Water SHAKE configuration:")
                jfastw = int_prompt_with_context(
                    self.processor,
                    "Fast water SHAKE (jfastw, 0=auto, 4=disabled)",
                    default=0,
                    module="AMBER Input Generator",
                    description="Fast water SHAKE jfastw",
                )
                if jfastw != 0:
                    self.config["jfastw"] = jfastw

        return True

    def _configure_restraints(self) -> bool:
        """Step 11: Positional Restraints."""
        self.console.print("\n[bold]Step 11: Positional Restraints[/bold]")

        self.console.print("\n[bold cyan]11a.[/bold cyan] Restraint method:")
        self.console.print("[cyan]1.[/cyan] No restraints")
        self.console.print("[cyan]2.[/cyan] Simple mask-based restraints (ntr=1)")
        self.console.print("[cyan]3.[/cyan] Advanced GROUP specification restraints")

        restraint_choice = prompt_with_context(
            self.processor,
            "Select restraint method",
            choices=["1", "2", "3"],
            default="1",
            module="AMBER Input Generator",
            description="Restraint method selection",
        )

        if restraint_choice == "1":
            return True
        elif restraint_choice == "2":
            return self._configure_simple_restraints()
        else:
            return self._configure_group_restraints()

    def _configure_simple_restraints(self) -> bool:
        """Configure simple mask-based restraints."""
        self.config["ntr"] = 1

        # Restraint weight
        self.console.print("\n[bold cyan]11b.[/bold cyan] Restraint force constant:")
        restraint_wt = float_prompt_with_context(
            self.processor,
            "Restraint force constant [kcal/mol/Å²]",
            default=10.0,
            module="AMBER Input Generator",
            description="Restraint force constant [kcal/mol/Å²]",
        )
        self.config["restraint_wt"] = restraint_wt

        # Restraint mask
        self.console.print("\n[bold cyan]11c.[/bold cyan] Restraint mask:")
        self.console.print("Common restraint masks:")
        self.console.print("!@H=           - All non-hydrogen atoms")
        self.console.print("@CA,C,N        - Protein backbone atoms")
        self.console.print(":1-10          - Residues 1 to 10")
        self.console.print("!:WAT,Na+,Cl-  - All except water and ions")

        default_mask = "!@H=" if self.simulation_type == "md" else "@CA,C,N"
        restraintmask = prompt_with_context(
            self.processor,
            "Restraint mask",
            default=default_mask,
            module="AMBER Input Generator",
            description="AMBER restraint mask",
        )
        self.config["restraintmask"] = restraintmask

        self.console.print("[cyan]ℹ️  Reference coordinates will be read from refc file[/cyan]")
        return True

    def _configure_group_restraints(self) -> bool:
        """Configure advanced GROUP specification restraints."""
        self.console.print("\n11b. GROUP specification restraints:")
        self.console.print("[cyan]ℹ️  GROUP restraints allow multiple groups with different force constants[/cyan]")
        
        groups = []
        group_num = 1
        
        while True:
            self.console.print(f"\n11b.{group_num}. Group {group_num} configuration:")
            
            # Group title
            title = prompt_with_context(
                self.processor,
                "Group title",
                default=f"Restraint group {group_num}",
                module="AMBER Input Generator",
                description=f"GROUP restraint title (group {group_num})",
            )
            
            # Force constant
            force_constant = float_prompt_with_context(
                self.processor,
                "Force constant [kcal/mol/Å²]",
                default=10.0,
                module="AMBER Input Generator",
                description="GROUP restraint force constant",
            )
            
            # Selection method
            self.console.print("\nSelection method:")
            self.console.print("1. Residue range (RES)", highlight=False)
            self.console.print("2. Atom range (ATOM)", highlight=False)
            self.console.print("3. Advanced filter (FIND)", highlight=False)
            
            method_choice = prompt_with_context(
                self.processor,
                "Select method [1-3]",
                choices=["1", "2", "3"],
                default="1",
                module="AMBER Input Generator",
                description="GROUP restraint selection method",
                options_map={
                    "1": "Residue range",
                    "2": "Atom range",
                    "3": "FIND search",
                },
            )
            
            group_lines = [title, str(force_constant)]
            
            if method_choice == "1":
                start_res = int_prompt_with_context(
                    self.processor,
                    "Start residue number",
                    default=1,
                    module="AMBER Input Generator",
                    description="GROUP restraint start residue",
                )
                end_res = int_prompt_with_context(
                    self.processor,
                    "End residue number",
                    default=10,
                    module="AMBER Input Generator",
                    description="GROUP restraint end residue",
                )
                group_lines.extend([
                    f"RES {start_res} {end_res}",
                    "END"
                ])
            elif method_choice == "2":
                start_atom = int_prompt_with_context(
                    self.processor,
                    "Start atom number",
                    default=1,
                    module="AMBER Input Generator",
                    description="GROUP restraint start atom",
                )
                end_atom = int_prompt_with_context(
                    self.processor,
                    "End atom number",
                    default=100,
                    module="AMBER Input Generator",
                    description="GROUP restraint end atom",
                )
                group_lines.extend([
                    f"ATOM {start_atom} {end_atom}",
                    "END"
                ])
            else:  # FIND method
                self.console.print("\nFilter specification (use * for wildcards):")
                atom_name = prompt_with_context(
                    self.processor,
                    "Atom name",
                    default="CA",
                    module="AMBER Input Generator",
                    description="FIND atom name",
                )
                atom_type = prompt_with_context(
                    self.processor,
                    "Atom type (* for any)",
                    default="*",
                    module="AMBER Input Generator",
                    description="FIND atom type",
                )
                tree_name = prompt_with_context(
                    self.processor,
                    "Tree name (* for any)",
                    default="*",
                    module="AMBER Input Generator",
                    description="FIND tree name",
                )
                residue_name = prompt_with_context(
                    self.processor,
                    "Residue name (* for any)",
                    default="*",
                    module="AMBER Input Generator",
                    description="FIND residue name",
                )
                
                start_res = int_prompt_with_context(
                    self.processor,
                    "Start residue for search",
                    default=1,
                    module="AMBER Input Generator",
                    description="FIND start residue",
                )
                end_res = int_prompt_with_context(
                    self.processor,
                    "End residue for search",
                    default=999,
                    module="AMBER Input Generator",
                    description="FIND end residue",
                )
                
                group_lines.extend([
                    "FIND",
                    f"{atom_name} {atom_type} {tree_name} {residue_name}",
                    "SEARCH",
                    f"RES {start_res} {end_res}",
                    "END"
                ])
            
            groups.append(group_lines)
            
            # Ask for more groups
            if not confirm_with_context(
                self.processor,
                "\nAdd another restraint group?",
                default=False,
                module="AMBER Input Generator",
                description="Add another GROUP restraint",
            ):
                break
            
            group_num += 1
        
        # Store groups for later inclusion in mdin file
        self.config["_restraint_groups"] = groups
        self.console.print(f"\n[green]Configured {len(groups)} restraint group(s)[/green]")
        self.console.print("[cyan]ℹ️  Reference coordinates will be read from refc file[/cyan]")
        
        return True

    def _configure_nonbonded(self) -> bool:
        """Step 12: Non-bonded Interactions."""
        self.console.print("\n[bold]Step 12: Non-bonded Interactions[/bold]")

        self.console.print("\n[bold cyan]12a.[/bold cyan] Cutoff distance:")
        if self.config.get("ntb", 1) > 0:
            # Periodic system - use PME
            cut = float_prompt_with_context(
                self.processor,
                "Nonbonded cutoff (cut) [Å]",
                default=10.0,
                module="AMBER Input Generator",
                description="Nonbonded cutoff cut [Å]",
            )
            self.config["cut"] = cut

            if cut < 8.0:
                self.console.print("[yellow]⚠️  Short cutoff - may affect accuracy[/yellow]")
            elif cut > 12.0:
                self.console.print("[yellow]⚠️  Long cutoff - will be expensive[/yellow]")

            self.console.print("[cyan]ℹ️  PME will be used automatically for electrostatics[/cyan]")

        else:
            # Non-periodic system - usually no cutoff or very large cutoff
            use_cutoff = confirm_with_context(
                self.processor,
                "Use nonbonded cutoff?",
                default=False,
                module="AMBER Input Generator",
                description="Use nonbonded cutoff (non-periodic)",
            )
            if use_cutoff:
                cut = float_prompt_with_context(
                    self.processor,
                    "Nonbonded cutoff [Å]",
                    default=999.0,
                    module="AMBER Input Generator",
                    description="Non-periodic nonbonded cutoff cut [Å]",
                )
                self.config["cut"] = cut
            else:
                self.config["cut"] = 999.0
                self.console.print("[cyan]ℹ️  No cutoff - all interactions calculated[/cyan]")

        return True

    def _configure_output(self) -> bool:
        """Step 13: Output Control."""
        self.console.print("\n[bold]Step 13: Output Control[/bold]")

        # Energy output frequency
        self.console.print("\n[bold cyan]13a.[/bold cyan] Energy output:")
        if self.simulation_type == "minimization":
            ntpr = int_prompt_with_context(
                self.processor,
                "Energy print frequency (ntpr)",
                default=100,
                module="AMBER Input Generator",
                description="Energy print frequency ntpr (minimization)",
            )
        else:
            ntpr = int_prompt_with_context(
                self.processor,
                "Energy print frequency (ntpr)",
                default=1000,
                module="AMBER Input Generator",
                description="Energy print frequency ntpr (MD)",
            )
        self.config["ntpr"] = ntpr

        # Restart file frequency (applies to both minimization and MD)
        self.console.print("\n[bold cyan]13b.[/bold cyan] Restart file output:")
        if self.simulation_type == "minimization":
            ntwr = int_prompt_with_context(
                self.processor,
                "Restart write frequency (ntwr)",
                default=500,
                module="AMBER Input Generator",
                description="Restart write frequency ntwr (minimization)",
            )
        else:
            ntwr = int_prompt_with_context(
                self.processor,
                "Restart write frequency (ntwr)",
                default=10000,
                module="AMBER Input Generator",
                description="Restart write frequency ntwr (MD)",
            )
        self.config["ntwr"] = ntwr

        # Final coordinate format (applies to both)
        self.console.print("\n[bold cyan]13c.[/bold cyan] Final coordinate format:")
        self.console.print("[cyan]1.[/cyan] ASCII")
        self.console.print("[cyan]2.[/cyan] NetCDF (recommended)")
        ntxo = prompt_with_context(
            self.processor,
            "Select format",
            choices=["1", "2"],
            default="2",
            module="AMBER Input Generator",
            description="Final coordinate format (ntxo)",
            options_map={"1": "ASCII", "2": "NetCDF (recommended)"},
        )
        self.config["ntxo"] = int(ntxo)

        # MD-specific output options
        if self.simulation_type != "minimization":
            # Trajectory output
            self.console.print("\n[bold cyan]13d.[/bold cyan] Trajectory output:")
            write_trajectory = confirm_with_context(
                self.processor,
                "Write coordinate trajectory?",
                default=True,
                module="AMBER Input Generator",
                description="Write coordinate trajectory",
            )
            if write_trajectory:
                ntwx = int_prompt_with_context(
                    self.processor,
                    "Trajectory write frequency (ntwx)",
                    default=1000,
                    module="AMBER Input Generator",
                    description="Trajectory write frequency ntwx",
                )
                self.config["ntwx"] = ntwx

                # Trajectory format
                self.console.print("\n[bold cyan]13e.[/bold cyan] Trajectory format:")
                self.console.print("[cyan]0.[/cyan] ASCII (large files)")
                self.console.print("[cyan]1.[/cyan] NetCDF binary (recommended)")

                ioutfm = prompt_with_context(
                    self.processor,
                    "Select format",
                    choices=["0", "1"],
                    default="1",
                    module="AMBER Input Generator",
                    description="Trajectory format (ioutfm)",
                    options_map={"0": "ASCII (large files)", "1": "NetCDF binary (recommended)"},
                )
                self.config["ioutfm"] = int(ioutfm)

                # Coordinate wrapping
                self.console.print("\n[bold cyan]13f.[/bold cyan] Coordinate wrapping:")
                iwrap = confirm_with_context(
                    self.processor,
                    "Wrap coordinates to primary box?",
                    default=False,
                    module="AMBER Input Generator",
                    description="Wrap coordinates to primary box (iwrap)",
                )
                if iwrap:
                    self.config["iwrap"] = 1
                    self.console.print("[cyan]ℹ️  Coordinates will be wrapped for visualization[/cyan]")

            # Additional output options
            self.console.print("\n[bold cyan]13g.[/bold cyan] Velocity trajectory:")
            write_velocities = confirm_with_context(
                self.processor,
                "Write velocity trajectory?",
                default=False,
                module="AMBER Input Generator",
                description="Write velocity trajectory",
            )
            if write_velocities:
                ntwv_choice = prompt_with_context(
                    self.processor,
                    "Velocity output: 1=separate file, -1=combined with coordinates",
                    choices=["1", "-1"],
                    default="1",
                    module="AMBER Input Generator",
                    description="Velocity output mode (ntwv)",
                    options_map={"1": "Separate file", "-1": "Combined with coordinates"},
                )
                self.config["ntwv"] = int(ntwv_choice)
            
            # Force trajectory
            self.console.print("\n[bold cyan]13h.[/bold cyan] Force trajectory:")
            write_forces = confirm_with_context(
                self.processor,
                "Write force trajectory?",
                default=False,
                module="AMBER Input Generator",
                description="Write force trajectory",
            )
            if write_forces:
                ntwf_choice = prompt_with_context(
                    self.processor,
                    "Force output: 1=separate file, -1=combined with coordinates",
                    choices=["1", "-1"],
                    default="1",
                    module="AMBER Input Generator",
                    description="Force output mode (ntwf)",
                    options_map={"1": "Separate file", "-1": "Combined with coordinates"},
                )
                self.config["ntwf"] = int(ntwf_choice)
            
            # Energy file
            self.console.print("\n[bold cyan]13i.[/bold cyan] Compact energy file:")
            write_energies = confirm_with_context(
                self.processor,
                "Write compact energy file (mden)?",
                default=False,
                module="AMBER Input Generator",
                description="Write compact energy file (mden)",
            )
            if write_energies:
                ntwe = int_prompt_with_context(
                    self.processor,
                    "Energy write frequency (ntwe)",
                    default=1000,
                    module="AMBER Input Generator",
                    description="Energy write frequency ntwe",
                )
                self.config["ntwe"] = ntwe

        return True

    def _configure_advanced_options(self) -> bool:
        """Step 14: Advanced Options."""
        self.console.print("\n[bold]Step 14: Advanced Options[/bold]")

        # NMR restraints and varying conditions
        if self.simulation_type in ["md", "custom"]:
            self.console.print("\n[bold cyan]14a.[/bold cyan] Varying conditions:")
            use_nmr = confirm_with_context(
                self.processor,
                "Use varying conditions (&wt blocks)?",
                default=False,
                module="AMBER Input Generator",
                description="Enable &wt varying conditions (nmropt=1)",
            )
            if use_nmr:
                self.config["nmropt"] = 1
                self.console.print("[cyan]ℹ️  nmropt=1 enables &wt varying conditions[/cyan]")
                
                # Configure &wt blocks
                self._configure_wt_blocks()

        # Energy decomposition
        self.console.print("\n[bold cyan]14b.[/bold cyan] Energy decomposition:")
        use_decomp = confirm_with_context(
            self.processor,
            "Enable energy decomposition?",
            default=False,
            module="AMBER Input Generator",
            description="Enable energy decomposition",
        )
        if use_decomp:
            self.console.print("Decomposition options:")
            self.console.print("[cyan]1.[/cyan] Per-residue (1-4 terms with internal)")
            self.console.print("[cyan]2.[/cyan] Per-residue (1-4 terms with non-bonded)")
            self.console.print("[cyan]3.[/cyan] Pairwise per-residue (1-4 with internal)")
            self.console.print("[cyan]4.[/cyan] Pairwise per-residue (1-4 with non-bonded)")

            idecomp = prompt_with_context(
                self.processor,
                "Select decomposition",
                choices=["1", "2", "3", "4"],
                default="1",
                module="AMBER Input Generator",
                description="Energy decomposition mode (idecomp)",
                options_map={
                    "1": "Per-residue (1-4 terms with internal)",
                    "2": "Per-residue (1-4 terms with non-bonded)",
                    "3": "Pairwise per-residue (1-4 with internal)",
                    "4": "Pairwise per-residue (1-4 with non-bonded)",
                },
            )
            self.config["idecomp"] = int(idecomp)

        # Belly dynamics (legacy)
        self.console.print("\n[bold cyan]14c.[/bold cyan] Belly dynamics:")
        if confirm_with_context(
            self.processor,
            "Use belly dynamics (subset of moving atoms)?",
            default=False,
            module="AMBER Input Generator",
            description="Use belly dynamics",
        ):
            self.config["ibelly"] = 1
            bellymask = prompt_with_context(
                self.processor,
                "Belly mask (moving atoms)",
                default="!:WAT",
                module="AMBER Input Generator",
                description="Belly mask (moving-atom AMBER mask)",
            )
            self.config["bellymask"] = bellymask
            self.console.print("[yellow]⚠️  Belly dynamics is legacy - consider restraints instead[/yellow]")

        return True

    def _configure_wt_blocks(self) -> bool:
        """Configure &wt varying conditions blocks."""
        self.console.print("\n[bold cyan]14a.1.[/bold cyan] &wt Block Configuration:")
        self.console.print("[cyan]ℹ️  &wt blocks allow parameters to change during simulation[/cyan]")
        
        # Main configuration options
        self.console.print("\n&wt Configuration options:")
        self.console.print("[cyan]1.[/cyan] Temperature ramping (heating protocols)")
        self.console.print("[cyan]2.[/cyan] Restraint ramping (gradual release)")
        self.console.print("[cyan]3.[/cyan] Temperature + restraint ramping (combined)")
        self.console.print("[cyan]4.[/cyan] Custom &wt blocks (advanced)")
        
        config_choice = prompt_with_context(
            self.processor,
            "Select configuration",
            choices=["1", "2", "3", "4"],
            default="1",
            module="AMBER Input Generator",
            description="&wt block configuration type",
        )
        
        wt_conditions = []
        
        if config_choice == "1":
            wt_conditions = self._configure_temperature_ramping()
        elif config_choice == "2":
            wt_conditions = self._configure_restraint_ramping()
        elif config_choice == "3":
            wt_conditions = self._configure_combined_ramping()
        else:
            wt_conditions = self._configure_custom_wt_blocks()
        
        # Store conditions
        self.config["_wt_conditions"] = wt_conditions
        
        self.console.print(f"\n[green]Configured {len(wt_conditions)} &wt condition(s)[/green]")
        return True

    def _configure_temperature_ramping(self) -> List[str]:
        """Configure temperature ramping protocols."""
        self.console.print("\n[bold cyan]Temperature Ramping Configuration[/bold cyan]")
        
        # Get total simulation steps
        nstlim = self.config.get("nstlim", 50000)
        
        wt_blocks = []
        ramp_num = 1
        last_step = 0
        
        while True:
            self.console.print(f"\n[bold cyan]Temperature Ramp {ramp_num}:[/bold cyan]")
            
            # Start step
            if ramp_num == 1:
                istep1 = int_prompt_with_context(
                    self.processor,
                    "Start step",
                    default=0,
                    module="AMBER Input Generator",
                    description=f"Temperature ramp {ramp_num} start step",
                )
            else:
                istep1 = int_prompt_with_context(
                    self.processor,
                    "Start step",
                    default=last_step,
                    module="AMBER Input Generator",
                    description=f"Temperature ramp {ramp_num} start step",
                )
            
            # End step with intelligent defaults
            remaining_steps = nstlim - istep1
            if ramp_num == 1:
                default_end = min(nstlim // 2, istep1 + remaining_steps // 2)
            else:
                default_end = nstlim
            
            istep2 = int_prompt_with_context(
                self.processor,
                f"End step (max: {nstlim})",
                default=default_end,
                module="AMBER Input Generator",
                description=f"Temperature ramp {ramp_num} end step",
            )
            
            # Validate step range
            if istep2 > nstlim:
                self.console.print(f"[yellow]⚠️  End step cannot exceed simulation length ({nstlim}). Using {nstlim}.[/yellow]")
                istep2 = nstlim
            
            # Temperature values
            if ramp_num == 1:
                temp1 = float_prompt_with_context(
                    self.processor,
                    "Initial temperature [K]",
                    default=0.0,
                    module="AMBER Input Generator",
                    description=f"Temperature ramp {ramp_num} initial temp",
                )
            else:
                temp1 = float_prompt_with_context(
                    self.processor,
                    "Initial temperature [K]",
                    default=300.0,
                    module="AMBER Input Generator",
                    description=f"Temperature ramp {ramp_num} initial temp",
                )

            temp2 = float_prompt_with_context(
                self.processor,
                "Final temperature [K]",
                default=300.0,
                module="AMBER Input Generator",
                description=f"Temperature ramp {ramp_num} final temp",
            )
            
            # Create &wt block
            if istep2 > istep1:
                wt_block = (
                    f"&wt type='TEMP0', istep1={istep1}, istep2={istep2},\n"
                    f"    value1={temp1}, value2={temp2} /"
                )
            else:
                wt_block = (
                    f"&wt type='TEMP0', istep1={istep1}, istep2={nstlim},\n"
                    f"    value1={temp1}, value2={temp2} /"
                )
                istep2 = nstlim
            
            wt_blocks.append(wt_block)
            last_step = istep2
            
            # Show preview
            self.console.print(f"\n[green]Created temperature ramp:[/green]")
            self.console.print(f"[grey50]{wt_block}[/grey50]")
            
            # Check if we've reached the end
            if istep2 >= nstlim:
                self.console.print(f"[cyan]ℹ️  Reached end of simulation at step {nstlim}[/cyan]")
                break
            
            # Ask for more ramps
            remaining_steps = nstlim - istep2
            if remaining_steps <= 0:
                break
                
            if not confirm_with_context(
                self.processor,
                f"\nAdd another temperature ramp? ({remaining_steps} steps remaining)",
                default=False,
                module="AMBER Input Generator",
                description="Add another temperature ramp",
            ):
                break
            
            ramp_num += 1
        
        return wt_blocks

    def _configure_restraint_ramping(self) -> List[str]:
        """Configure restraint ramping protocols."""
        self.console.print("\n[bold cyan]Restraint Ramping Configuration[/bold cyan]")
        
        wt_blocks = []
        ramp_num = 1
        
        while True:
            self.console.print(f"\n[bold cyan]Restraint Ramp {ramp_num}:[/bold cyan]")
            
            # Get parameters for this ramp
            istep1 = int_prompt_with_context(
                self.processor,
                "Start step",
                default=0,
                module="AMBER Input Generator",
                description=f"Restraint ramp {ramp_num} start step",
            )
            istep2 = int_prompt_with_context(
                self.processor,
                "End step (0 = end of simulation)",
                default=0,
                module="AMBER Input Generator",
                description=f"Restraint ramp {ramp_num} end step",
            )
            weight1 = float_prompt_with_context(
                self.processor,
                "Initial restraint weight",
                default=10.0,
                module="AMBER Input Generator",
                description=f"Restraint ramp {ramp_num} initial weight",
            )
            weight2 = float_prompt_with_context(
                self.processor,
                "Final restraint weight",
                default=1.0,
                module="AMBER Input Generator",
                description=f"Restraint ramp {ramp_num} final weight",
            )
            
            # Create &wt block
            if istep2 > 0:
                wt_block = (
                    f"&wt type='REST', istep1={istep1}, istep2={istep2},\n"
                    f"    value1={weight1}, value2={weight2} /"
                )
            else:
                wt_block = (
                    f"&wt type='REST', istep1={istep1}, istep2=0,\n"
                    f"    value1={weight1} /"
                )
            
            wt_blocks.append(wt_block)
            
            # Show preview
            self.console.print(f"\n[green]Created restraint ramp:[/green]")
            self.console.print(f"[grey50]{wt_block}[/grey50]")
            
            # Ask for more ramps
            if not confirm_with_context(
                self.processor,
                "\nAdd another restraint ramp?",
                default=False,
                module="AMBER Input Generator",
                description="Add another restraint ramp",
            ):
                break
            
            ramp_num += 1
        
        return wt_blocks

    def _configure_combined_ramping(self) -> List[str]:
        """Configure combined temperature and restraint ramping."""
        self.console.print("\n[bold cyan]Combined Temperature + Restraint Ramping[/bold cyan]")
        
        wt_blocks = []
        
        # Temperature ramps
        self.console.print("\n[bold]Temperature Configuration:[/bold]")
        temp_blocks = self._configure_temperature_ramping()
        wt_blocks.extend(temp_blocks)
        
        # Restraint ramps
        self.console.print("\n[bold]Restraint Configuration:[/bold]")
        restraint_blocks = self._configure_restraint_ramping()
        wt_blocks.extend(restraint_blocks)
        
        return wt_blocks

    def _configure_custom_wt_blocks(self) -> List[str]:
        """Configure custom &wt blocks."""
        self.console.print("\n[bold cyan]Custom &wt Block Configuration[/bold cyan]")
        
        wt_blocks = []
        block_num = 1
        
        while True:
            self.console.print(f"\n[bold cyan]&wt Block {block_num}:[/bold cyan]")
            
            # Parameter type
            self.console.print("\nParameter to vary:")
            self.console.print("[cyan]1.[/cyan] TEMP0 - Target temperature")
            self.console.print("[cyan]2.[/cyan] REST - All NMR restraint weights")  
            self.console.print("[cyan]3.[/cyan] BOND - Bond energy weights")
            self.console.print("[cyan]4.[/cyan] ANGLE - Angle energy weights")
            self.console.print("[cyan]5.[/cyan] TORSION - Torsion energy weights")
            self.console.print("[cyan]6.[/cyan] VDW - van der Waals weights")
            self.console.print("[cyan]7.[/cyan] ELEC - Electrostatic weights")
            self.console.print("[cyan]8.[/cyan] NB - All non-bonded weights")
            self.console.print("[cyan]9.[/cyan] TAUTP - Temperature coupling")
            self.console.print("[cyan]10.[/cyan] Custom parameter")
            
            param_choice = prompt_with_context(
                self.processor,
                "Select parameter",
                choices=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
                default="1",
                module="AMBER Input Generator",
                description="&wt parameter to vary",
                options_map={
                    "1": "TEMP0 - Target temperature",
                    "2": "REST - NMR restraint weights",
                    "3": "BOND - Bond energy",
                    "4": "ANGLE - Angle energy",
                    "5": "TORSION - Torsion energy",
                    "6": "VDW - van der Waals",
                    "7": "ELEC - Electrostatic",
                    "8": "NB - All non-bonded",
                    "9": "TAUTP - Temperature coupling",
                    "10": "Custom parameter",
                },
            )
            
            param_types = {
                "1": "TEMP0", "2": "REST", "3": "BOND", "4": "ANGLE",
                "5": "TORSION", "6": "VDW", "7": "ELEC", "8": "NB", "9": "TAUTP"
            }
            
            if param_choice == "10":
                param_type = prompt_with_context(
                    self.processor,
                    "Parameter name",
                    module="AMBER Input Generator",
                    description="Custom &wt parameter name",
                ).upper()
            else:
                param_type = param_types[param_choice]

            # Step range
            istep1 = int_prompt_with_context(
                self.processor,
                "Start step (istep1)",
                default=0,
                module="AMBER Input Generator",
                description=f"&wt {param_type} istep1",
            )
            istep2 = int_prompt_with_context(
                self.processor,
                "End step (istep2, 0=until end)",
                default=0,
                module="AMBER Input Generator",
                description=f"&wt {param_type} istep2",
            )

            # Values
            value1 = float_prompt_with_context(
                self.processor,
                "Initial value (value1)",
                default=0.0,
                module="AMBER Input Generator",
                description=f"&wt {param_type} value1",
            )

            if istep2 > 0:
                value2 = float_prompt_with_context(
                    self.processor,
                    "Final value (value2)",
                    default=1.0,
                    module="AMBER Input Generator",
                    description=f"&wt {param_type} value2",
                )

                # Interpolation type
                linear = confirm_with_context(
                    self.processor,
                    "Linear interpolation (vs multiplicative)?",
                    default=True,
                    module="AMBER Input Generator",
                    description=f"&wt {param_type} linear interpolation",
                )
                imult = 0 if linear else 1

                # Step size
                iinc = int_prompt_with_context(
                    self.processor,
                    "Step increment (0=continuous)",
                    default=0,
                    module="AMBER Input Generator",
                    description=f"&wt {param_type} step increment",
                )
                
                wt_block = (
                    f"&wt type='{param_type}', istep1={istep1}, istep2={istep2},\n"
                    f"    value1={value1}, value2={value2}, imult={imult}, iinc={iinc} /"
                )
            else:
                wt_block = (
                    f"&wt type='{param_type}', istep1={istep1}, istep2=0,\n"
                    f"    value1={value1} /"
                )
            
            wt_blocks.append(wt_block)
            
            # Show what was created
            self.console.print(f"\n[green]Created &wt block:[/green]")
            self.console.print(f"[grey50]{wt_block}[/grey50]")
            
            # Ask for more blocks
            if not confirm_with_context(
                self.processor,
                "\nAdd another &wt block?",
                default=False,
                module="AMBER Input Generator",
                description="Add another custom &wt block",
            ):
                break
            
            block_num += 1
        
        return wt_blocks

    def _generate_editable_template(self) -> str:
        """Generate an editable template with comments explaining each section."""
        content = self._build_mdin_content()
        
        # Add explanatory header
        header = """# AMBER Input File Template - Generated by ProPrep
# Edit the parameters below as needed. Changes will be preserved.
# Lines starting with # are comments and will be ignored.
# 
# Common parameters:
#   imin: 0=MD, 1=minimization, 5=trajectory analysis
#   nstlim: Number of MD steps  
#   dt: Time step in ps
#   temp0/tempi: Target/initial temperature in K
#   ntpr: Print energies every N steps
#   ntwx: Write coordinates every N steps
#
# For more details, see the AMBER manual or use 'ambmask -help'

"""
        return header + content

    def _parse_editable_template(self, template_content: str) -> bool:
        """Parse the edited template back into configuration."""
        try:
            # Clear current config but preserve special keys
            special_keys = {k: v for k, v in self.config.items() if k.startswith('_')}
            self.config.clear()
            self.config.update(special_keys)
            
            # Parse the template content
            current_section = None
            
            for line in template_content.split('\n'):
                line = line.strip()
                
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue
                
                # Detect sections
                if line.startswith('&'):
                    current_section = line
                    continue
                elif line == '/':
                    current_section = None
                    continue
                
                # Parse parameters in &cntrl section
                if current_section == '&cntrl':
                    if '=' in line:
                        # Handle comma-separated parameters on one line
                        parts = line.rstrip(',').split(',')
                        for part in parts:
                            if '=' in part:
                                key, value = part.split('=', 1)
                                key = key.strip()
                                value = value.strip().rstrip(',')
                                
                                # Convert to appropriate type
                                try:
                                    # Try integer first
                                    if '.' not in value:
                                        self.config[key] = int(value)
                                    else:
                                        self.config[key] = float(value)
                                except ValueError:
                                    # Keep as string, removing quotes if present
                                    if value.startswith('"') and value.endswith('"'):
                                        value = value[1:-1]
                                    elif value.startswith("'") and value.endswith("'"):
                                        value = value[1:-1]
                                    self.config[key] = value
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error parsing template: {e}[/red]")
            return False

    def _edit_template_interactively(self, template: str) -> str:
        """Allow user to edit the template interactively."""
        try:
            import tempfile
            import subprocess
            import os
            
            # Create temporary file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.in', delete=False) as tmp:
                tmp.write(template)
                tmp_path = tmp.name
            
            # Get editor preference
            editor = os.environ.get('EDITOR', 'nano')
            
            try:
                # Open editor
                self.console.print(f"[cyan]Opening template in {editor}...[/cyan]")
                self.console.print("[grey50]Save and exit the editor when done.[/grey50]")
                
                result = subprocess.run([editor, tmp_path], check=True)
                
                # Read back the edited content
                with open(tmp_path, 'r') as f:
                    edited_content = f.read()
                
                return edited_content
                
            finally:
                # Clean up
                if os.path.exists(tmp_path):
                    os.unlink(tmp_path)
                    
        except Exception as e:
            self.console.print(f"[red]Error opening editor: {e}[/red]")
            self.console.print("[yellow]Falling back to current template[/yellow]")
            return template

    def _review_and_generate(self) -> bool:
        """Step 15: Final Review and Generation using editable template."""
        self.console.print("\n[bold]Step 15: Review and Generate[/bold]")
        
        # Generate editable template
        template = self._generate_editable_template()
        
        # Show the template
        from rich.syntax import Syntax
        self.console.print("\n[bold cyan]Generated AMBER Input Template[/bold cyan]")
        self.console.print("[yellow]You can edit this template directly. Changes will be preserved.[/yellow]")
        self.console.print("─" * 60)
        
        syntax = Syntax(template, "fortran", theme="monokai", line_numbers=True)
        self.console.print(syntax)
        self.console.print("─" * 60)
        
        # Edit template if user wants
        if confirm_with_context(
            self.processor,
            "Edit the template?",
            default=True,
            module="AMBER Input Generator",
            description="Edit the generated template",
        ):
            self.console.print("[yellow]⚠️  Once you edit the template, you're responsible for parameter correctness.[/yellow]")
            self.console.print("[yellow]   ProPrep will save your changes as-is without validation.[/yellow]")
            
            if confirm_with_context(
                self.processor,
                "Proceed with editing?",
                default=True,
                module="AMBER Input Generator",
                description="Proceed with editing template",
            ):
                template = self._edit_template_interactively(template)
                self.console.print("[green]✓ Template edited successfully[/green]")
        
        # Final options
        self.console.print("\n[bold]Final Options:[/bold]")
        self.console.print("1. Save input file", highlight=False)
        self.console.print("2. Save as template", highlight=False)
        self.console.print("3. Generate run script", highlight=False)
        self.console.print("4. All of the above", highlight=False)
        self.console.print("5. Return to main menu", highlight=False)

        choice = prompt_with_context(
            self.processor,
            "Enter choice [1-5]",
            choices=["1", "2", "3", "4", "5"],
            default="1",
            module="AMBER Input Generator",
            description="Final action (save/template/script/all/return)",
        )

        success = True
        
        if choice in ["1", "4"]:
            success &= self._save_final_input_file(template)
            
        if choice in ["2", "4"]:
            success &= self._save_template_file(template)
            
        if choice in ["3", "4"]:
            success &= self._export_run_script()
            
        return success

    def _save_final_input_file(self, template_content: str) -> bool:
        """Save the final input file."""
        default_name = f"{self.simulation_type}.in"
        filename = prompt_with_context(
            self.processor,
            "Output filename",
            default=default_name,
            module="AMBER Input Generator",
            description="Output filename",
        )

        if not filename.endswith(".in"):
            filename += ".in"

        # Check if file exists
        if os.path.exists(filename):
            if not confirm_with_context(
                self.processor,
                f"File {filename} exists. Overwrite?",
                default=False,
                module="AMBER Input Generator",
                description=f"Overwrite existing file {filename}",
            ):
                self.console.print("[yellow]File save cancelled[/yellow]")
                return False

        try:
            with open(filename, 'w') as f:
                f.write(template_content)
            
            self.console.print(f"[green]✓ Input file saved: {filename}[/green]")
            self._show_usage_info(filename)
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error writing file: {e}[/red]")
            return False

    def _save_template_file(self, template_content: str) -> bool:
        """Save as a reusable template."""
        default_name = f"{self.simulation_type}_template.in"
        filename = prompt_with_context(
            self.processor,
            "Template filename",
            default=default_name,
            module="AMBER Input Generator",
            description="Template output filename",
        )

        if not filename.endswith(".in"):
            filename += ".in"

        try:
            with open(filename, 'w') as f:
                f.write(template_content)
            
            self.console.print(f"[green]✓ Template saved: {filename}[/green]")
            self.console.print("[cyan]This template can be reused and modified for future simulations[/cyan]")
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error writing template: {e}[/red]")
            return False

    def _display_comprehensive_summary(self):
        """Display comprehensive configuration summary."""
        table = Table(title="Configuration Summary")
        table.add_column("Category", style="cyan")
        table.add_column("Parameter", style="blue")
        table.add_column("Value", style="green")

        # Organize by categories
        categories = {
            "Simulation": ["imin", "simulation_type", "engine_choice"],
            "Time": ["nstlim", "dt", "t"],
            "Temperature": ["ntt", "temp0", "tempi", "gamma_ln"],
            "Pressure": ["ntp", "ntb", "pres0", "barostat"],
            "Constraints": ["ntc", "ntf", "tol"],
            "Output": ["ntpr", "ntwx", "ntwr", "ioutfm"],
            "Advanced": ["nmropt", "idecomp", "nrespa"]
        }

        for category, params in categories.items():
            for i, param in enumerate(params):
                if param == "simulation_type":
                    value = self.simulation_type or "None"
                elif param == "engine_choice":
                    value = f"{self.engine_choice}" + (" + GPU" if self.config.get("_use_gpu") else "")
                elif param in self.config:
                    value = str(self.config[param])
                else:
                    continue
                    
                cat_display = category if i == 0 else ""
                table.add_row(cat_display, param, value)

        self.console.print(table)

        # Runtime calculation for MD
        if self.simulation_type != "minimization" and "nstlim" in self.config and "dt" in self.config:
            runtime_ps = self.config["nstlim"] * self.config["dt"]
            runtime_ns = runtime_ps / 1000
            self.console.print(f"\n[bold]Total simulation time:[/bold] {runtime_ps:.1f} ps ({runtime_ns:.3f} ns)")

        # Parameter count
        param_count = len([k for k in self.config.keys() if not k.startswith("_")])
        self.console.print(f"[grey50]Total parameters configured: {param_count}[/grey50]")

    def _generate_mdin_file(self) -> bool:
        """Generate the final mdin file."""
        if not self.config:
            self.console.print("[red]No configuration available[/red]")
            return False

        # Get filename with safety check
        default_name = f"{self.simulation_type}.in"
        filename = prompt_with_context(
            self.processor,
            "Output filename",
            default=default_name,
            module="AMBER Input Generator",
            description="Output filename",
        )

        if not filename.endswith(".in"):
            filename += ".in"

        # Check if file exists
        if os.path.exists(filename):
            if not confirm_with_context(
                self.processor,
                f"File {filename} exists. Overwrite?",
                default=False,
                module="AMBER Input Generator",
                description=f"Overwrite existing file {filename}",
            ):
                self.console.print("[yellow]File generation cancelled[/yellow]")
                return False

        # Generate content
        content = self._build_mdin_content()

        try:
            with open(filename, 'w') as f:
                f.write(content)
            
            self.console.print(f"[green]✓ Input file generated: {filename}[/green]")
            
            # Show what was included
            components = ["&cntrl section"]
            
            if self.config.get("_wt_conditions"):
                components.append("&wt varying conditions")
            
            if self.config.get("_file_redirections"):
                components.append("file redirection")
                
            if self.config.get("_restraint_groups"):
                components.append("GROUP specification")
            
            if self.config.get("_nmr_protocol"):
                components.append("NMR heating protocol")
            
            self.console.print(f"[cyan]Components included: {', '.join(components)}[/cyan]")
            
            # Show usage info
            self._show_usage_info(filename)
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error writing file: {e}[/red]")
            return False

    def _build_mdin_content(self) -> str:
        """Build the complete mdin file content."""
        lines = []
        
        # Title
        title = f"{self.simulation_type.title()} simulation - Generated by ProPrep AMBER Input Generator"
        lines.append(title)
        
        # &cntrl section
        lines.append("&cntrl")
        
        # Add parameters in logical order
        ordered_params = [
            # Basic control
            "imin", "ntx", "irest",
            # Time control
            "nstlim", "dt", "t",
            # Minimization
            "maxcyc", "ncyc", "ntmin", "dx0", "drms",
            # Temperature
            "ntt", "temp0", "tempi", "gamma_ln", "tautp", "vrand", "ig",
            # Pressure
            "ntp", "ntb", "pres0", "barostat", "mcbarint", "taup", "comp",
            # Constraints
            "ntc", "ntf", "tol", "jfastw",
            # Restraints
            "ntr", "restraint_wt", "restraintmask",
            # Belly
            "ibelly", "bellymask",
            # Nonbonded
            "cut",
            # Output
            "ntpr", "ntwx", "ntwr", "ntxo", "ioutfm", "iwrap",
            # Advanced
            "nscm", "nrespa", "nmropt", "idecomp",
        ]
        
        for param in ordered_params:
            if param in self.config and not param.startswith("_"):
                value = self.config[param]
                if isinstance(value, str) and " " in value:
                    lines.append(f"  {param} = '{value}',")
                else:
                    lines.append(f"  {param} = {value},")
        
        lines.append("/")
        
        # Add &wt section if present
        wt_section = self._generate_wt_section()
        if wt_section:
            lines.append("")
            lines.append(wt_section)
        
        # Add file redirection if present
        redirection_section = self._generate_file_redirection_section()
        if redirection_section:
            lines.append("")
            lines.append(redirection_section)
        
        # Add GROUP specification if present
        group_section = self._generate_group_specification()
        if group_section:
            lines.append("")
            lines.append(group_section)
        
        return "\n".join(lines) + "\n"

    def _generate_wt_section(self) -> str:
        """Generate &wt section for varying conditions."""
        if not self.config.get("nmropt", 0):
            return ""
        
        # Use configured wt conditions if available
        wt_conditions = self.config.get("_wt_conditions", [])
        if wt_conditions:
            lines = []
            for condition in wt_conditions:
                lines.append(condition)
            lines.append("&wt type='END' /")
            return "\n".join(lines)
        
        # Default heating protocol for heating simulations (fallback)
        if self.simulation_type == "heating" and self.config.get("tempi", 0) == 0:
            temp0 = self.config.get("temp0", 300.0)
            nstlim = self.config.get("nstlim", 50000)
            
            return (
                f"&wt type='TEMP0', istep1=0, istep2={nstlim//2},\n"
                f"    value1=0.0, value2={temp0} /\n"
                "&wt type='END' /"
            )
        
        return ""

    def _generate_file_redirection_section(self) -> str:
        """Generate file redirection section."""
        redirections = self.config.get("_file_redirections", {})
        if not redirections:
            return ""
        
        lines = []
        for redirect_type, filename in redirections.items():
            lines.append(f"{redirect_type} = {filename}")
        
        return "\n".join(lines)

    def _generate_group_specification(self) -> str:
        """Generate GROUP specification."""
        groups = self.config.get("_restraint_groups", [])
        if not groups:
            return ""
        
        lines = []
        for group in groups:
            lines.extend(group)
        lines.append("END")
        
        return "\n".join(lines)

    def _show_usage_info(self, filename: str):
        """Show usage information for the generated file."""
        self.console.print(f"\n[bold]Usage Information:[/bold]")
        
        command = self._get_engine_command(filename)
        self.console.print(f"[cyan]Run command:[/cyan] {command}")
        
        # Required files
        required_files = ["prmtop", "inpcrd"]
        if self.config.get("irest", 0) == 1:
            required_files.append("restrt (restart file)")
        if self.config.get("ntr", 0) == 1:
            required_files.append("refc (reference coordinates)")
        
        self.console.print(f"[cyan]Required files:[/cyan] {', '.join(required_files)}")
        
        # Optional files
        optional_files = []
        if self.config.get("nmropt", 0) > 0:
            optional_files.append("DISANG (restraint file)")
        
        if optional_files:
            self.console.print(f"[cyan]Optional files:[/cyan] {', '.join(optional_files)}")

    def _get_engine_command(self, mdin_file: str) -> str:
        """Get the appropriate engine command with proper MPI/GPU setup."""
        base_files = f"-i {mdin_file} -o mdout -p prmtop -c inpcrd"
        
        if self.config.get("irest", 0) == 1:
            base_files += " -r restrt"
        else:
            base_files += " -r restrt"
        
        if self.config.get("ntwx", 0) > 0:
            base_files += " -x mdcrd"
        
        if self.config.get("ntr", 0) == 1 or self.config.get("_restraint_groups"):
            base_files += " -ref refc"
        
        # Build command with proper parallelization
        if self.config.get("_use_gpu", False):
            # GPU command with CUDA_VISIBLE_DEVICES
            gpu_ids = self.config.get("_gpu_ids", "0")
            return f"export CUDA_VISIBLE_DEVICES={gpu_ids} && {self.engine_choice} {base_files}"
        elif self.config.get("_use_mpi", False):
            # MPI command
            ncpus = self.config.get("_ncpus", 1)
            return f"mpirun -np {ncpus} {self.engine_choice}.MPI {base_files}"
        else:
            # Serial command
            return f"{self.engine_choice} {base_files}"

    # ============================================================================
    # WORKFLOW GENERATION METHODS
    # ============================================================================

    def _generate_workflow_inputs(self) -> bool:
        """Generate complete workflow input files using editable templates."""
        try:
            self.console.print("[bold]Generate Complete AMBER Workflow[/bold]")
            
            # Show available presets
            if not self.workflow_presets:
                self.console.print("[yellow]No workflow presets available[/yellow]")
                return False
            
            self.console.print("\nAvailable workflow presets:")
            preset_names = list(self.workflow_presets.keys())
            
            table = Table()
            table.add_column("Option", style="cyan")
            table.add_column("Preset Name", style="green")
            table.add_column("Description", style="white")
            
            for i, preset_key in enumerate(preset_names, 1):
                preset = self.workflow_presets[preset_key]
                table.add_row(str(i), preset["name"], preset["description"])
            
            self.console.print(table)
            
            # Get user choice
            choices = [str(i) for i in range(1, len(preset_names) + 1)]
            preset_options_map = {
                str(i + 1): self.workflow_presets[k]["name"]
                for i, k in enumerate(preset_names)
            }
            choice = prompt_with_context(
                self.processor,
                f"Select preset [1-{len(preset_names)}]",
                choices=choices,
                default="1",
                module="AMBER Input Generator",
                description="Select workflow preset",
                options_map=preset_options_map,
            )
            
            preset_key = preset_names[int(choice) - 1]
            preset = self.workflow_presets[preset_key]
            
            self.console.print(f"\n[green]Selected: {preset['name']}[/green]")
            self.console.print(f"Description: {preset['description']}")
            
            # Store workflow information
            self.simulation_type = "workflow"
            self.current_workflow = {
                "preset_key": preset_key,
                "preset": preset,
                "step_templates": {},
                "step_contents": {}
            }
            
            # Generate and edit templates for each workflow step
            return self._workflow_template_editor()
            
        except Exception as e:
            self.console.print(f"[red]Error generating workflow: {e}[/red]")
            return False

    def _workflow_template_editor(self) -> bool:
        """Handle the workflow template editing process."""
        preset = self.current_workflow["preset"]
        step_names = list(preset["steps"].keys())
        
        self.console.print(f"\n[bold cyan]Workflow Template Editor[/bold cyan]")
        self.console.print(f"[yellow]You will edit templates for {len(step_names)} workflow steps:[/yellow]")
        
        for i, step_name in enumerate(step_names, 1):
            step_config = preset["steps"][step_name]
            self.console.print(f"  {i}. {step_name}: {step_config['description']}")
        
        # Generate templates for each step
        for step_name in step_names:
            step_config = preset["steps"][step_name]
            self.console.print(f"\n[bold]Step {step_names.index(step_name) + 1}: {step_name}[/bold]")
            self.console.print(f"[grey50]{step_config['description']}[/grey50]")
            
            # Generate initial template
            template = self._generate_workflow_step_template(step_name, step_config)
            self.current_workflow["step_templates"][step_name] = template
            
            # Show template
            from rich.syntax import Syntax
            syntax = Syntax(template, "fortran", theme="monokai", line_numbers=True)
            self.console.print(syntax)
            self.console.print("─" * 60)
            
            # Edit template if user wants
            if confirm_with_context(
                self.processor,
                f"Edit template for {step_name}?",
                default=True,
                module="AMBER Input Generator",
                description=f"Edit template for {step_name}",
            ):
                self.console.print("[yellow]⚠️  Once you edit the template, you're responsible for parameter correctness.[/yellow]")
                self.console.print("[yellow]   ProPrep will save your changes as-is without validation.[/yellow]")
                
                if confirm_with_context(
                self.processor,
                "Proceed with editing?",
                default=True,
                module="AMBER Input Generator",
                description="Proceed with editing template",
            ):
                    edited_template = self._edit_template_interactively(template)
                    self.current_workflow["step_templates"][step_name] = edited_template
                    self.console.print(f"[green]✓ {step_name} template edited successfully[/green]")
            
            # Store final content
            self.current_workflow["step_contents"][step_name] = self.current_workflow["step_templates"][step_name]
        
        # Final workflow options
        return self._workflow_final_options()

    def _generate_workflow_step_template(self, step_name: str, step_config: dict) -> str:
        """Generate an editable template for a workflow step."""
        # Build template content from step configuration
        content = self._generate_mdin_content(
            step_config["config"],
            step_config["description"],
            step_config.get("nmr_section", ""),
            step_config.get("group_restraints")
        )
        
        # Add helpful header comments for workflow context
        header = f"""# AMBER Workflow Step: {step_name}
# {step_config['description']}
# Generated by ProPrep - Edit as needed
#
# This is step {list(self.current_workflow['preset']['steps'].keys()).index(step_name) + 1} of {len(self.current_workflow['preset']['steps'])} in the workflow
# Input files will be chained automatically between steps
#
# Key parameters for this step:
"""
        
        # Add step-specific parameter explanations
        config = step_config["config"]
        if config.get("imin") == 1:
            header += "#   imin=1: Energy minimization\n"
            header += "#   maxcyc: Maximum minimization cycles\n"
            header += "#   ncyc: Steepest descent cycles before switching to conjugate gradient\n"
        elif config.get("ntt"):
            if "tempi" in config and "temp0" in config and config["tempi"] != config["temp0"]:
                header += "#   Heating: tempi → temp0 over nstlim steps\n"
            else:
                header += "#   MD simulation at constant temperature\n"
            header += "#   nstlim: Number of MD steps\n"
            header += "#   dt: Time step in picoseconds\n"
        
        header += "\n"
        
        return header + content

    def _workflow_final_options(self) -> bool:
        """Handle final options for workflow generation."""
        self.console.print("\n[bold]Workflow Final Options:[/bold]")
        self.console.print("1. Save workflow input files (separate .in files)", highlight=False)
        self.console.print("2. Save as templates (reusable)", highlight=False)
        self.console.print("3. Generate run scripts", highlight=False)
        self.console.print("4. All of the above", highlight=False)
        self.console.print("5. Return to main menu", highlight=False)

        choice = prompt_with_context(
            self.processor,
            "Enter choice [1-5]",
            choices=["1", "2", "3", "4", "5"],
            default="1",
            module="AMBER Input Generator",
            description="Final action (save/template/script/all/return)",
        )

        success = True
        
        if choice in ["1", "4"]:
            success &= self._save_workflow_input_files()
            
        if choice in ["2", "4"]:
            success &= self._save_workflow_templates()
            
        if choice in ["3", "4"]:
            success &= self._workflow_run_script_options()
            
        return success

    def _save_workflow_input_files(self) -> bool:
        """Save workflow input files with auto-generated names."""
        try:
            preset = self.current_workflow["preset"]
            step_contents = self.current_workflow["step_contents"]
            saved_files = []
            
            self.console.print("\n[bold]Saving workflow input files...[/bold]")
            
            for step_name, content in step_contents.items():
                # Use preset filename or generate one
                step_config = preset["steps"][step_name]
                default_filename = step_config.get("filename", f"{step_name}.in")
                
                # Check if file exists
                if os.path.exists(default_filename):
                    if not confirm_with_context(
                        self.processor,
                        f"File {default_filename} exists. Overwrite?",
                        default=False,
                        module="AMBER Input Generator",
                        description=f"Overwrite existing file {default_filename}",
                    ):
                        continue
                
                # Save file
                with open(default_filename, 'w') as f:
                    f.write(content)
                
                saved_files.append(default_filename)
                self.console.print(f"  ✓ Saved {default_filename}")
            
            self.console.print(f"\n[green]Successfully saved {len(saved_files)} workflow files![/green]")
            
            # Show saved files summary
            table = Table(title="Saved Workflow Files")
            table.add_column("File", style="cyan")
            table.add_column("Description", style="white")
            
            for step_name in step_contents.keys():
                step_config = preset["steps"][step_name]
                filename = step_config.get("filename", f"{step_name}.in")
                if filename in saved_files:
                    table.add_row(filename, step_config["description"])
            
            self.console.print(table)
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error saving workflow files: {e}[/red]")
            return False

    def _save_workflow_templates(self) -> bool:
        """Save workflow templates for reuse."""
        try:
            step_contents = self.current_workflow["step_contents"]
            saved_templates = []
            
            self.console.print("\n[bold]Saving workflow templates...[/bold]")
            
            for step_name, content in step_contents.items():
                template_filename = f"{step_name}_template.in"
                
                # Check if file exists
                if os.path.exists(template_filename):
                    if not confirm_with_context(
                        self.processor,
                        f"Template {template_filename} exists. Overwrite?",
                        default=False,
                        module="AMBER Input Generator",
                        description=f"Overwrite existing template {template_filename}",
                    ):
                        continue
                
                # Save template
                with open(template_filename, 'w') as f:
                    f.write(content)
                
                saved_templates.append(template_filename)
                self.console.print(f"  ✓ Saved {template_filename}")
            
            self.console.print(f"\n[green]Successfully saved {len(saved_templates)} template files![/green]")
            self.console.print("[cyan]These templates can be loaded and modified for future workflows.[/cyan]")
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error saving templates: {e}[/red]")
            return False

    def _workflow_run_script_options(self) -> bool:
        """Let user choose run script options."""
        self.console.print("\n[bold]Run Script Options:[/bold]")
        self.console.print("1. Master script (runs entire workflow)", highlight=False)
        self.console.print("2. Individual scripts (separate script for each step)", highlight=False)
        self.console.print("3. Both master and individual scripts", highlight=False)
        
        choice = prompt_with_context(
            self.processor,
            "Select run script type [1-3]",
            choices=["1", "2", "3"],
            default="1",
            module="AMBER Input Generator",
            description="Run script type",
        )
        
        success = True
        
        if choice in ["1", "3"]:
            success &= self._generate_master_workflow_script()
            
        if choice in ["2", "3"]:
            success &= self._generate_individual_workflow_scripts()
            
        return success

    def _generate_master_workflow_script(self) -> bool:
        """Generate a master workflow script that runs all steps."""
        try:
            # Get topology and coordinate files from user
            topology_file, coord_file = self._discover_workflow_files()
            if not topology_file or not coord_file:
                return False
            
            preset = self.current_workflow["preset"]
            step_contents = self.current_workflow["step_contents"]
            
            script_name = "run_workflow.sh"
            
            # Check if file exists
            if os.path.exists(script_name):
                if not confirm_with_context(
                    self.processor,
                    f"Script {script_name} exists. Overwrite?",
                    default=False,
                    module="AMBER Input Generator",
                    description=f"Overwrite existing script {script_name}",
                ):
                    return False
            
            # Build script content
            lines = [
                "#!/bin/bash",
                "# AMBER Workflow Master Script",
                f"# Generated by ProPrep - Workflow: {preset['name']}",
                "",
                "set -e  # Exit on any error",
                "",
                "# Check for required files",
                f"if [ ! -f {topology_file} ]; then",
                f"    echo 'Error: Topology file {topology_file} not found'",
                "    exit 1",
                "fi",
                "",
                f"if [ ! -f {coord_file} ]; then",
                f"    echo 'Error: Coordinate file {coord_file} not found'",
                "    exit 1",
                "fi",
                "",
                "# Set AMBER engine (modify as needed)",
                "ENGINE=${ENGINE:-pmemd}",
                "",
                f"echo 'Starting {preset['name']} workflow...'",
                "echo 'Using topology: {topology_file}'",
                f"echo 'Using coordinates: {coord_file}'",
                "echo",
                ""
            ]
            
            # Add steps
            steps = list(preset["steps"].items())
            prev_coord = coord_file
            
            for i, (step_name, step_config) in enumerate(steps):
                step_num = i + 1
                input_file = step_config.get("filename", f"{step_name}.in")
                output_coord = f"{step_name}.rst7"
                output_file = f"{step_name}.out"

                lines.extend([
                    f"# Step {step_num}: {step_config['description']}",
                    f"echo 'Step {step_num}: {step_config['description']}'",
                    f"echo 'Running: $ENGINE -O -i {input_file} -o {output_file} -p {topology_file} -c {prev_coord} -r {output_coord}'",
                    f"$ENGINE -O -i {input_file} -o {output_file} -p {topology_file} -c {prev_coord} -r {output_coord}",
                    f"EXIT_CODE=$?",
                    f"wait  # Ensure all background processes complete",
                    f"if [ $EXIT_CODE -ne 0 ]; then",
                    f"    echo 'Error in step {step_num}: {step_name} (exit code: '$EXIT_CODE')'",
                    "    exit 1",
                    "fi",
                    f"echo 'Step {step_num} completed successfully'",
                    "echo",
                    ""
                ])

                prev_coord = output_coord
            
            lines.extend([
                "echo 'Workflow completed successfully!'",
                f"echo 'Final coordinates: {prev_coord}'",
                ""
            ])
            
            # Write script
            with open(script_name, 'w') as f:
                f.write("\n".join(lines))
            
            # Make executable
            os.chmod(script_name, 0o755)
            
            self.console.print(f"[green]✓ Master workflow script generated: {script_name}[/green]")
            self.console.print(f"[cyan]Usage: ./{script_name}[/cyan]")
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error generating master workflow script: {e}[/red]")
            return False

    def _generate_individual_workflow_scripts(self) -> bool:
        """Generate individual run scripts for each workflow step."""
        try:
            # Get topology and coordinate files from user
            topology_file, coord_file = self._discover_workflow_files()
            if not topology_file or not coord_file:
                return False
            
            preset = self.current_workflow["preset"]
            generated_scripts = []
            
            self.console.print("\n[bold]Generating individual step scripts...[/bold]")
            
            steps = list(preset["steps"].items())
            prev_coord = coord_file
            
            for i, (step_name, step_config) in enumerate(steps):
                step_num = i + 1
                script_name = f"run_{step_name}.sh"
                input_file = step_config.get("filename", f"{step_name}.in")
                output_coord = f"{step_name}.rst7"
                output_file = f"{step_name}.out"
                
                # Check if script exists
                if os.path.exists(script_name):
                    if not confirm_with_context(
                    self.processor,
                    f"Script {script_name} exists. Overwrite?",
                    default=False,
                    module="AMBER Input Generator",
                    description=f"Overwrite existing script {script_name}",
                ):
                        continue
                
                # Build script content
                lines = [
                    "#!/bin/bash",
                    f"# AMBER Step {step_num} Script: {step_name}",
                    f"# {step_config['description']}",
                    f"# Generated by ProPrep",
                    "",
                    "set -e  # Exit on any error",
                    "",
                    "# Check for required files",
                    f"if [ ! -f {topology_file} ]; then",
                    f"    echo 'Error: Topology file {topology_file} not found'",
                    "    exit 1",
                    "fi",
                    "",
                    f"if [ ! -f {prev_coord} ]; then",
                    f"    echo 'Error: Coordinate file {prev_coord} not found'",
                    "    exit 1",
                    "fi",
                    "",
                    f"if [ ! -f {input_file} ]; then",
                    f"    echo 'Error: Input file {input_file} not found'",
                    "    exit 1",
                    "fi",
                    "",
                    "# Set AMBER engine (modify as needed)",
                    "ENGINE=${ENGINE:-pmemd}",
                    "",
                    f"echo 'Running step {step_num}: {step_config['description']}'",
                    f"echo 'Input: {input_file}'",
                    f"echo 'Coordinates: {prev_coord}'",
                    f"echo 'Output: {output_file}'",
                    "echo",
                    "",
                    "# Run AMBER",
                    f"$ENGINE -O -i {input_file} -o {output_file} -p {topology_file} -c {prev_coord} -r {output_coord}",
                    "",
                    f"echo 'Step {step_num} completed successfully'",
                    f"echo 'Output coordinates: {output_coord}'",
                    ""
                ]
                
                # Write script
                with open(script_name, 'w') as f:
                    f.write("\n".join(lines))
                
                # Make executable
                os.chmod(script_name, 0o755)
                
                generated_scripts.append(script_name)
                self.console.print(f"  ✓ Generated {script_name}")
                
                prev_coord = output_coord
            
            if generated_scripts:
                self.console.print(f"\n[green]Successfully generated {len(generated_scripts)} individual scripts![/green]")
                
                # Show summary
                table = Table(title="Generated Individual Scripts")
                table.add_column("Script", style="cyan")
                table.add_column("Step", style="white")
                table.add_column("Description", style="grey50")
                
                for i, (step_name, step_config) in enumerate(steps):
                    script_name = f"run_{step_name}.sh"
                    if script_name in generated_scripts:
                        table.add_row(script_name, f"Step {i+1}", step_config["description"])
                
                self.console.print(table)
                
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error generating individual workflow scripts: {e}[/red]")
            return False

    def _discover_workflow_files(self) -> Tuple[str, str]:
        """Discover and let user select topology and coordinate files for workflow."""
        try:
            # Find topology files
            topology_patterns = ["*.prmtop", "*.parm7", "*.top"]
            topology_files = []
            for pattern in topology_patterns:
                topology_files.extend(glob.glob(pattern))
            
            if not topology_files:
                self.console.print("[red]No topology files found (.prmtop, .parm7, .top)[/red]")
                return None, None
            
            # Find coordinate files
            coord_patterns = ["*.rst7", "*.inpcrd", "*.crd"]
            coord_files = []
            for pattern in coord_patterns:
                coord_files.extend(glob.glob(pattern))
            
            if not coord_files:
                self.console.print("[red]No coordinate files found (.rst7, .inpcrd, .crd)[/red]")
                return None, None
            
            # Let user select topology file
            self.console.print("\n[bold]Select topology file:[/bold]")
            for i, f in enumerate(topology_files, 1):
                self.console.print(f"{i}. {f}")
            
            while True:
                try:
                    choice = int(prompt_with_context(
                        self.processor,
                        f"Select topology file [1-{len(topology_files)}]",
                        module="AMBER Input Generator",
                        description="Select topology (prmtop) file",
                    ))
                    if 1 <= choice <= len(topology_files):
                        topology_file = topology_files[choice - 1]
                        break
                    else:
                        self.console.print("[red]Invalid choice[/red]")
                except ValueError:
                    self.console.print("[red]Please enter a number[/red]")
            
            # Let user select coordinate file
            self.console.print(f"\n[bold]Select initial coordinate file:[/bold]")
            for i, f in enumerate(coord_files, 1):
                self.console.print(f"{i}. {f}")
            
            while True:
                try:
                    choice = int(prompt_with_context(
                        self.processor,
                        f"Select coordinate file [1-{len(coord_files)}]",
                        module="AMBER Input Generator",
                        description="Select coordinate (rst7/inpcrd) file",
                    ))
                    if 1 <= choice <= len(coord_files):
                        coord_file = coord_files[choice - 1]
                        break
                    else:
                        self.console.print("[red]Invalid choice[/red]")
                except ValueError:
                    self.console.print("[red]Please enter a number[/red]")
            
            return topology_file, coord_file
            
        except Exception as e:
            self.console.print(f"[red]Error discovering files: {e}[/red]")
            return None, None

    def _generate_preset_files(self, preset: Dict) -> bool:
        """Generate all input files from a preset."""
        try:
            self.console.print(f"\n[bold]Generating workflow files...[/bold]")
            
            generated_files = []
            
            for step_name, step_config in preset["steps"].items():
                filename = step_config["filename"]
                
                # Check if file exists
                if os.path.exists(filename):
                    if not confirm_with_context(
                self.processor,
                f"File {filename} exists. Overwrite?",
                default=False,
                module="AMBER Input Generator",
                description=f"Overwrite existing file {filename}",
            ):
                        self.console.print("[yellow]Workflow generation cancelled[/yellow]")
                        return False
                
                # Generate file content
                content = self._generate_mdin_content(
                    step_config["config"],
                    step_config["description"],
                    step_config.get("nmr_section", ""),
                    step_config.get("group_restraints")
                )
                
                # Write file
                with open(filename, 'w') as f:
                    f.write(content)
                generated_files.append(filename)
                self.console.print(f"  ✅ Generated {filename}")
            
            self.console.print(f"\n[green]Successfully generated {len(generated_files)} input files![/green]")
            
            # Show generated files
            table = Table(title="Generated Workflow Files")
            table.add_column("File", style="cyan")
            table.add_column("Description", style="white")
            
            for step_name, step_config in preset["steps"].items():
                table.add_row(step_config["filename"], step_config["description"])
            
            self.console.print(table)
            
            # Offer to generate run script
            generate_script = confirm_with_context(
                self.processor,
                "\nGenerate run script?",
                default=True,
                module="AMBER Input Generator",
                description="Generate run script",
            )
            if generate_script:
                self._generate_workflow_run_script(preset)
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error generating preset files: {e}[/red]")
            return False

    def _generate_mdin_content(self, config: Dict, description: str, nmr_section: str = "", group_restraints: List[Dict] = None) -> str:
        """Generate mdin file content from configuration.

        Args:
            config: Dictionary of AMBER parameters for &cntrl namelist
            description: Description line for the mdin file
            nmr_section: Optional NMR restraints section (e.g., &wt namelists)
            group_restraints: Optional list of GROUP restraint specifications
        """
        lines = [description, "&cntrl"]

        # Add parameters in organized order
        # Note: restraintmask and restraint_wt should NOT be present if using GROUP
        ordered_params = [
            "imin", "ntx", "irest", "nstlim", "dt", "maxcyc", "ncyc", "ntmin",
            "ntpr", "ntwx", "ntwr", "ntxo", "ioutfm",
            "ntt", "temp0", "tempi", "gamma_ln", "ntp", "ntb", "pres0", "barostat",
            "ntc", "ntf", "ntr", "restraint_wt", "restraintmask", "cut", "nmropt"
        ]

        for param in ordered_params:
            if param in config:
                value = config[param]
                if isinstance(value, str) and " " in value:
                    lines.append(f"  {param} = '{value}',")
                else:
                    lines.append(f"  {param} = {value},")

        # Add any remaining parameters
        for key, value in config.items():
            if key not in ordered_params:
                if isinstance(value, str) and " " in value:
                    lines.append(f"  {key} = '{value}',")
                else:
                    lines.append(f"  {key} = {value},")

        lines.append("/")

        # Add NMR section if present
        if nmr_section:
            lines.append("")
            lines.append(nmr_section)

        # Add GROUP restraints if present
        # IMPORTANT: GROUP specification must have NO indentation!
        if group_restraints:
            lines.append("")
            for group in group_restraints:
                lines.append(group.get("title", "GROUP restraint"))
                lines.append(str(group.get("force_constant", 10.0)))

                # Add FIND section if specified
                if group.get("find_criteria"):
                    lines.append("FIND")
                    # find_criteria is a list of criteria strings
                    for criterion in group["find_criteria"]:
                        lines.append(criterion)
                    lines.append("SEARCH")

                # Add RES lines with residue ranges
                if group.get("residue_ranges"):
                    # Format: RES start1 end1 start2 end2 ... (up to 7 pairs per line)
                    ranges = group["residue_ranges"]
                    res_line = "RES"
                    for start, end in ranges:
                        res_line += f" {start} {end}"
                    lines.append(res_line)

                lines.append("END")

            # Final END to terminate all GROUP specifications
            lines.append("END")

        return "\n".join(lines) + "\n"

    def _generate_workflow_run_script(self, preset: Dict) -> bool:
        """Generate run script for workflow."""
        try:
            script_name = "run_workflow.sh"
            
            if os.path.exists(script_name):
                if not confirm_with_context(
                    self.processor,
                    f"Script {script_name} exists. Overwrite?",
                    default=False,
                    module="AMBER Input Generator",
                    description=f"Overwrite existing script {script_name}",
                ):
                    return False
            
            script_content = self._build_workflow_script_content(preset)
            
            with open(script_name, 'w') as f:
                f.write(script_content)
            
            # Make executable
            os.chmod(script_name, 0o755)
            
            self.console.print(f"[green]Run script generated: {script_name}[/green]")
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error generating run script: {e}[/red]")
            return False

    def _build_workflow_script_content(self, preset: Dict) -> str:
        """Build the workflow run script content."""
        lines = [
            "#!/bin/bash",
            "# AMBER Workflow Run Script",
            f"# Generated by ProPrep AMBER Input Generator",
            f"# Workflow: {preset['name']}",
            "",
            "# Check for required files",
            "if [ ! -f prmtop ]; then",
            "    echo 'Error: prmtop file not found'",
            "    exit 1",
            "fi",
            "",
            "if [ ! -f inpcrd ]; then",
            "    echo 'Error: inpcrd file not found'",
            "    exit 1",
            "fi",
            "",
            "# Set engine (modify as needed)",
            "ENGINE=pmemd",
            "",
            "echo 'Starting AMBER workflow...'",
            ""
        ]
        
        steps = list(preset["steps"].items())
        
        for i, (step_name, step_config) in enumerate(steps):
            filename = step_config["filename"]
            description = step_config["description"]
            
            lines.extend([
                f"# Step {i+1}: {description}",
                f"echo 'Running {description}...'",
            ])
            
            # Determine input and output files
            if i == 0:
                # First step
                input_coord = "inpcrd"
                restart_flag = ""
            else:
                # Subsequent steps
                prev_step = steps[i-1][1]["filename"].replace(".in", "")
                input_coord = f"{prev_step}.rst"
                restart_flag = "-r"
            
            output_base = filename.replace(".in", "")
            
            # Build command
            cmd_parts = [
                f"$ENGINE -i {filename}",
                f"-o {output_base}.out",
                "-p prmtop",
                f"-c {input_coord}",
                f"-r {output_base}.rst"
            ]
            
            # Add trajectory output if specified
            if step_config["config"].get("ntwx", 0) > 0:
                cmd_parts.append(f"-x {output_base}.nc")
            
            # Add reference coordinates for restraints
            if step_config["config"].get("ntr", 0) == 1:
                if i == 0:
                    cmd_parts.append("-ref inpcrd")
                else:
                    cmd_parts.append(f"-ref {input_coord}")
            
            command = " \\\n    ".join(cmd_parts) + " < /dev/null"

            lines.extend([
                command,
                "EXIT_CODE=$?",
                "wait  # Ensure all background processes complete",
                "if [ $EXIT_CODE -ne 0 ]; then",
                f"    echo 'Error in {description} (exit code: '$EXIT_CODE')'",
                "    exit 1",
                "fi",
                f"echo '{description} completed successfully'",
                ""
            ])
        
        lines.extend([
            "echo 'Workflow completed successfully!'",
            "echo 'Output files:'",
        ])
        
        for step_name, step_config in preset["steps"].items():
            output_base = step_config["filename"].replace(".in", "")
            description = step_config["description"]
            lines.append(f"echo '  {output_base}.out - {description} output'")
            lines.append(f"echo '  {output_base}.rst - {description} restart'")
            if step_config["config"].get("ntwx", 0) > 0:
                lines.append(f"echo '  {output_base}.nc - {description} trajectory'")
        
        return "\n".join(lines) + "\n"

    def _customize_workflow_preset(self, preset: Dict) -> bool:
        """Allow user to customize workflow preset parameters."""
        try:
            self.console.print(f"\n[bold]Customizing: {preset['name']}[/bold]")
            
            # Create a copy of the preset to modify
            custom_preset = copy.deepcopy(preset)
            
            # Show current steps and allow selection for modification
            steps = list(custom_preset["steps"].keys())
            
            while True:
                self.console.print("\nWorkflow steps:")
                table = Table()
                table.add_column("Option", style="cyan")
                table.add_column("Step", style="green")
                table.add_column("Description", style="white")
                
                for i, step_name in enumerate(steps, 1):
                    step_config = custom_preset["steps"][step_name]
                    table.add_row(str(i), step_name.replace('_', ' ').title(), step_config["description"])
                
                table.add_row(str(len(steps) + 1), "Modify global settings", "Change settings that affect all steps")
                table.add_row(str(len(steps) + 2), "Generate with current settings", "Create files with current parameters")
                
                self.console.print(table)
                
                choices = [str(i) for i in range(1, len(steps) + 3)]
                choice = prompt_with_context(
                    self.processor,
                    f"Select step to customize [1-{len(steps) + 2}]",
                    choices=choices,
                    default=str(len(steps) + 2),
                    module="AMBER Input Generator",
                    description="Select step to customize",
                )
                
                if choice == str(len(steps) + 2):  # Generate with current settings
                    break
                elif choice == str(len(steps) + 1):  # Modify global settings
                    self._modify_global_settings(custom_preset)
                else:
                    # Modify specific step
                    step_index = int(choice) - 1
                    step_name = steps[step_index]
                    self._modify_step_settings(custom_preset, step_name)
                
                # Ask if done customizing
                done = confirm_with_context(
                    self.processor,
                    "\nFinished customizing?",
                    default=False,
                    module="AMBER Input Generator",
                    description="Finished customizing steps",
                )
                if done:
                    break
            
            # Generate files with customized preset
            return self._generate_preset_files(custom_preset)
            
        except Exception as e:
            self.console.print(f"[red]Error customizing workflow: {e}[/red]")
            return False

    def _modify_global_settings(self, preset: Dict) -> bool:
        """Modify global settings that affect all steps."""
        try:
            self.console.print("\n[bold]Global Settings Modification[/bold]")
            
            # Common global changes
            options = {
                "1": "Change target temperature",
                "2": "Change time step",
                "3": "Change output frequencies",
                "4": "Change cutoff distance",
                "5": "Return to step selection"
            }
            
            for key, desc in options.items():
                self.console.print(f"{key}. {desc}")
            
            choice = prompt_with_context(
                self.processor,
                "Select global change [1-5]",
                choices=list(options.keys()),
                default="5",
                module="AMBER Input Generator",
                description="Select global workflow-wide change",
            )
            
            if choice == "1":
                new_temp = float_prompt_with_context(
                    self.processor,
                    "New target temperature [K]",
                    default=300.0,
                    module="AMBER Input Generator",
                    description="New workflow-wide target temperature [K]",
                )
                for step_config in preset["steps"].values():
                    if "temp0" in step_config["config"]:
                        step_config["config"]["temp0"] = new_temp
                self.console.print(f"[green]Updated temperature to {new_temp}K in all applicable steps[/green]")
            
            elif choice == "2":
                new_dt = float_prompt_with_context(
                    self.processor,
                    "New time step [ps]",
                    default=0.002,
                    module="AMBER Input Generator",
                    description="New workflow-wide time step [ps]",
                )
                for step_config in preset["steps"].values():
                    if "dt" in step_config["config"]:
                        step_config["config"]["dt"] = new_dt
                self.console.print(f"[green]Updated time step to {new_dt}ps in all applicable steps[/green]")
            
            elif choice == "3":
                new_ntpr = int_prompt_with_context(
                    self.processor,
                    "New energy print frequency",
                    default=1000,
                    module="AMBER Input Generator",
                    description="New workflow-wide energy print frequency",
                )
                new_ntwx = int_prompt_with_context(
                    self.processor,
                    "New trajectory write frequency",
                    default=1000,
                    module="AMBER Input Generator",
                    description="New workflow-wide trajectory write frequency",
                )
                for step_config in preset["steps"].values():
                    if "ntpr" in step_config["config"]:
                        step_config["config"]["ntpr"] = new_ntpr
                    if "ntwx" in step_config["config"]:
                        step_config["config"]["ntwx"] = new_ntwx
                self.console.print(f"[green]Updated output frequencies in all steps[/green]")
            
            elif choice == "4":
                new_cut = float_prompt_with_context(
                    self.processor,
                    "New cutoff distance [Å]",
                    default=10.0,
                    module="AMBER Input Generator",
                    description="New workflow-wide cutoff [Å]",
                )
                for step_config in preset["steps"].values():
                    if "cut" in step_config["config"]:
                        step_config["config"]["cut"] = new_cut
                self.console.print(f"[green]Updated cutoff to {new_cut}Å in all steps[/green]")
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error modifying global settings: {e}[/red]")
            return False

    def _modify_step_settings(self, preset: Dict, step_name: str) -> bool:
        """Modify settings for a specific step."""
        try:
            step_config = preset["steps"][step_name]
            
            self.console.print(f"\n[bold]Customizing: {step_name.replace('_', ' ').title()}[/bold]")
            self.console.print(f"Description: {step_config['description']}")
            
            # Show current parameters
            self.console.print("\nCurrent parameters:")
            param_table = Table()
            param_table.add_column("Parameter", style="cyan")
            param_table.add_column("Value", style="white")
            
            for param, value in step_config["config"].items():
                param_table.add_row(param, str(value))
            
            self.console.print(param_table)
            
            # Ask which parameters to modify
            modify = confirm_with_context(
                self.processor,
                f"\nModify parameters for {step_name}?",
                default=False,
                module="AMBER Input Generator",
                description=f"Modify parameters for step {step_name}",
            )
            if not modify:
                return True
            
            # Step-specific customizations
            if step_name == "minimization":
                self._customize_minimization_step(step_config)
            elif step_name == "heating":
                self._customize_heating_step(step_config)
            elif "equilibration" in step_name:
                self._customize_equilibration_step(step_config, step_name)
            elif step_name == "production":
                self._customize_production_step(step_config)
            else:
                self._customize_generic_step(step_config)
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error modifying step settings: {e}[/red]")
            return False

    def _customize_minimization_step(self, step_config: Dict):
        """Customize minimization step parameters."""
        config = step_config["config"]
        
        # Common minimization parameters
        MOD = "AMBER Input Generator"
        if confirm_with_context(self.processor, "Change maximum cycles?", default=False,
                                module=MOD, description="Change maxcyc"):
            config["maxcyc"] = int_prompt_with_context(self.processor, "Maximum cycles", default=config.get("maxcyc", 10000),
                                                       module=MOD, description="Maxcyc (minimization customize)")

        if confirm_with_context(self.processor, "Change steepest descent cycles?", default=False,
                                module=MOD, description="Change ncyc"):
            config["ncyc"] = int_prompt_with_context(self.processor, "Steepest descent cycles", default=config.get("ncyc", 1000),
                                                     module=MOD, description="Ncyc (minimization customize)")

        if confirm_with_context(self.processor, "Change restraint weight?", default=False,
                                module=MOD, description="Change minimization restraint weight"):
            config["restraint_wt"] = float_prompt_with_context(self.processor, "Restraint weight [kcal/mol/Å²]",
                                                               default=config.get("restraint_wt", 500.0),
                                                               module=MOD, description="Minimization restraint weight")

    def _customize_heating_step(self, step_config: Dict):
        """Customize heating step parameters."""
        config = step_config["config"]
        MOD = "AMBER Input Generator"

        if confirm_with_context(self.processor, "Change simulation length?", default=False,
                                module=MOD, description="Change heating nstlim"):
            config["nstlim"] = int_prompt_with_context(self.processor, "Number of steps", default=config.get("nstlim", 50000),
                                                       module=MOD, description="Heating nstlim")

        if confirm_with_context(self.processor, "Change initial temperature?", default=False,
                                module=MOD, description="Change heating tempi"):
            config["tempi"] = float_prompt_with_context(self.processor, "Initial temperature [K]", default=config.get("tempi", 0.0),
                                                        module=MOD, description="Heating initial temperature tempi")

        if confirm_with_context(self.processor, "Change final temperature?", default=False,
                                module=MOD, description="Change heating temp0"):
            config["temp0"] = float_prompt_with_context(self.processor, "Final temperature [K]", default=config.get("temp0", 300.0),
                                                        module=MOD, description="Heating final temperature temp0")

        if confirm_with_context(self.processor, "Change restraint weight?", default=False,
                                module=MOD, description="Change heating restraint weight"):
            config["restraint_wt"] = float_prompt_with_context(self.processor, "Restraint weight [kcal/mol/Å²]",
                                                               default=config.get("restraint_wt", 10.0),
                                                               module=MOD, description="Heating restraint weight")

    def _customize_equilibration_step(self, step_config: Dict, step_name: str):
        """Customize equilibration step parameters."""
        config = step_config["config"]
        MOD = "AMBER Input Generator"

        if confirm_with_context(self.processor, "Change simulation length?", default=False,
                                module=MOD, description="Change equilibration nstlim"):
            config["nstlim"] = int_prompt_with_context(self.processor, "Number of steps", default=config.get("nstlim", 100000),
                                                       module=MOD, description="Equilibration nstlim")

        if confirm_with_context(self.processor, "Change pressure control?", default=False,
                                module=MOD, description="Change equilibration pressure control"):
            if confirm_with_context(self.processor, "Enable pressure control?", default=config.get("ntp", 0) > 0,
                                    module=MOD, description="Enable pressure control"):
                config["ntp"] = 1
                config["ntb"] = 2
                config["pres0"] = float_prompt_with_context(self.processor, "Target pressure [bar]", default=config.get("pres0", 1.0),
                                                            module=MOD, description="Equilibration target pressure [bar]")
            else:
                config["ntp"] = 0
                config["ntb"] = 1

        if confirm_with_context(self.processor, "Change restraints?", default=False,
                                module=MOD, description="Change equilibration restraints"):
            if confirm_with_context(self.processor, "Apply restraints?", default=config.get("ntr", 0) > 0,
                                    module=MOD, description="Apply equilibration restraints"):
                config["ntr"] = 1
                config["restraint_wt"] = float_prompt_with_context(self.processor, "Restraint weight [kcal/mol/Å²]",
                                                                   default=config.get("restraint_wt", 1.0),
                                                                   module=MOD, description="Equilibration restraint weight")
                config["restraintmask"] = prompt_with_context(self.processor, "Restraint mask",
                                                              default=config.get("restraintmask", "!@H="),
                                                              module=MOD, description="Equilibration restraint mask")
            else:
                config["ntr"] = 0

    def _customize_production_step(self, step_config: Dict):
        """Customize production step parameters."""
        config = step_config["config"]
        MOD = "AMBER Input Generator"

        if confirm_with_context(self.processor, "Change simulation length?", default=False,
                                module=MOD, description="Change production nstlim"):
            config["nstlim"] = int_prompt_with_context(self.processor, "Number of steps", default=config.get("nstlim", 2500000),
                                                       module=MOD, description="Production nstlim")

            # Calculate and show runtime
            dt = config.get("dt", 0.002)
            runtime_ps = config["nstlim"] * dt
            runtime_ns = runtime_ps / 1000
            self.console.print(f"[cyan]Total runtime: {runtime_ps:.1f} ps ({runtime_ns:.3f} ns)[/cyan]")

        if confirm_with_context(self.processor, "Change output frequencies?", default=False,
                                module=MOD, description="Change production output frequencies"):
            config["ntpr"] = int_prompt_with_context(self.processor, "Energy print frequency",
                                                     default=config.get("ntpr", 1000),
                                                     module=MOD, description="Production ntpr")
            config["ntwx"] = int_prompt_with_context(self.processor, "Trajectory write frequency",
                                                     default=config.get("ntwx", 1000),
                                                     module=MOD, description="Production ntwx")
            config["ntwr"] = int_prompt_with_context(self.processor, "Restart write frequency",
                                                     default=config.get("ntwr", 10000),
                                                     module=MOD, description="Production ntwr")

    def _customize_generic_step(self, step_config: Dict):
        """Customize generic step parameters."""
        config = step_config["config"]
        
        # Show all parameters and allow selection
        params = list(config.keys())
        if not params:
            self.console.print("[yellow]No parameters to customize[/yellow]")
            return
        
        self.console.print("Available parameters:")
        for i, param in enumerate(params, 1):
            self.console.print(f"{i}. {param} = {config[param]}")
        
        while True:
            param_choices = [str(i) for i in range(1, len(params) + 1)] + ["0"]
            param_options_map = {str(i + 1): p for i, p in enumerate(params)}
            param_options_map["0"] = "Finish"
            choice = prompt_with_context(
                self.processor,
                f"Select parameter to modify [1-{len(params)}] (0 to finish)",
                choices=param_choices,
                default="0",
                module="AMBER Input Generator",
                description="Select generic step parameter to modify",
                options_map=param_options_map,
            )

            if choice == "0":
                break

            param = params[int(choice) - 1]
            current_value = config[param]

            # Determine input type
            if isinstance(current_value, bool):
                new_value = confirm_with_context(
                    self.processor,
                    f"New value for {param}",
                    default=current_value,
                    module="AMBER Input Generator",
                    description=f"New value for {param}",
                )
            elif isinstance(current_value, int):
                new_value = int_prompt_with_context(
                    self.processor,
                    f"New value for {param}",
                    default=current_value,
                    module="AMBER Input Generator",
                    description=f"New value for {param}",
                )
            elif isinstance(current_value, float):
                new_value = float_prompt_with_context(
                    self.processor,
                    f"New value for {param}",
                    default=current_value,
                    module="AMBER Input Generator",
                    description=f"New value for {param}",
                )
            else:
                new_value = prompt_with_context(
                    self.processor,
                    f"New value for {param}",
                    default=str(current_value),
                    module="AMBER Input Generator",
                    description=f"New value for {param}",
                )
            
            config[param] = new_value
            self.console.print(f"[green]Updated {param} to {new_value}[/green]")

    # ============================================================================
    # TEMPLATE MANAGEMENT METHODS
    # ============================================================================

    def _load_template(self) -> bool:
        """Load configuration from template."""
        try:
            templates = self._list_available_templates()
            
            if not templates:
                self.console.print("[yellow]No templates available[/yellow]")
                return False
            
            # Show available templates
            table = Table(title="Available Templates")
            table.add_column("Option", style="cyan")
            table.add_column("Name", style="green")
            table.add_column("Source", style="blue")
            table.add_column("Description", style="white")
            
            template_keys = list(templates.keys())
            for i, template_key in enumerate(template_keys, 1):
                template_info = templates[template_key]
                template_data = template_info["data"]
                source = template_info["source"]
                name = template_data.get("name", template_key)
                description = template_data.get("description", "No description")
                
                table.add_row(str(i), name, source, description)
            
            self.console.print(table)
            
            # Get user choice
            choices = [str(i) for i in range(1, len(template_keys) + 1)]
            choice = prompt_with_context(
                self.processor,
                f"Select template [1-{len(template_keys)}]",
                choices=choices,
                module="AMBER Input Generator",
                description="Select custom template",
            )
            
            selected_key = template_keys[int(choice) - 1]
            selected_template = templates[selected_key]["data"]
            
            # Load configuration
            self.config = selected_template.get("config", {}).copy()
            self.simulation_type = selected_template.get("simulation_type", "md")
            
            self.console.print(f"[green]Template loaded: {selected_template.get('name')}[/green]")
            
            # Display loaded configuration
            self._display_loaded_template_config(selected_template)
            
            # Offer review and editing options
            return self._review_and_edit_loaded_template()
            
        except Exception as e:
            self.console.print(f"[red]Error loading template: {e}[/red]")
            return False

    def _list_available_templates(self) -> Dict:
        """List all available templates."""
        templates = {}
        
        # Load builtin templates
        for template_file in self.builtin_templates_dir.glob("*.json"):
            try:
                with open(template_file, 'r') as f:
                    template_data = json.load(f)
                    templates[template_file.stem] = {
                        "data": template_data,
                        "source": "builtin"
                    }
            except Exception as e:
                self.console.print(f"[yellow]Warning: Could not load template {template_file}: {e}[/yellow]")
        
        # Load user templates
        for template_file in self.templates_dir.glob("*.json"):
            if template_file.parent != self.builtin_templates_dir:
                try:
                    with open(template_file, 'r') as f:
                        template_data = json.load(f)
                        templates[template_file.stem] = {
                            "data": template_data,
                            "source": "user"
                        }
                except Exception as e:
                    self.console.print(f"[yellow]Warning: Could not load template {template_file}: {e}[/yellow]")
        
        return templates

    def _display_loaded_template_config(self, template_data: Dict) -> None:
        """Display the configuration that was just loaded from template."""
        self.console.print(f"\n[bold]Loaded Configuration Summary[/bold]")
        
        # Basic info
        info_table = Table()
        info_table.add_column("Property", style="cyan")
        info_table.add_column("Value", style="green")
        
        info_table.add_row("Template Name", template_data.get("name", "Unknown"))
        info_table.add_row("Description", template_data.get("description", "N/A"))
        info_table.add_row("Simulation Type", template_data.get("simulation_type", "Unknown"))
        info_table.add_row("Parameters", str(len(template_data.get("config", {}))))
        
        self.console.print(info_table)

    def _review_and_edit_loaded_template(self) -> bool:
        """Review and optionally edit loaded template configuration."""
        try:
            while True:
                self.console.print(f"\n[bold]Template Configuration Options[/bold]")
                self.console.print("1. Use template as-is", highlight=False)
                self.console.print("2. Review configuration", highlight=False)
                self.console.print("3. Edit parameters", highlight=False)
                self.console.print("4. Generate input file", highlight=False)
                self.console.print("5. Save as new template", highlight=False)
                self.console.print("6. Return to main menu", highlight=False)
                
                choice = prompt_with_context(
                    self.processor,
                    "Select option [1-6]",
                    choices=["1", "2", "3", "4", "5", "6"],
                    default="1",
                    module="AMBER Input Generator",
                    description="Template/workflow management action",
                )
                
                if choice == "1":
                    self.console.print("[green]Template configuration accepted[/green]")
                    return True
                elif choice == "2":
                    self._display_summary()
                elif choice == "3":
                    self._modify_current_settings()
                elif choice == "4":
                    return self._generate_mdin_file()
                elif choice == "5":
                    return self._save_template()
                else:
                    return True
                    
        except Exception as e:
            self.console.print(f"[red]Error in template review: {e}[/red]")
            return False

    def _save_template(self) -> bool:
        """Save current configuration as template."""
        try:
            if not self.config:
                self.console.print("[yellow]No configuration to save[/yellow]")
                return False

            # Handle workflow vs single configuration differently
            if self.simulation_type == "workflow":
                return self._save_workflow_template()
            else:
                return self._save_single_template()
                
        except Exception as e:
            self.console.print(f"[red]Error saving template: {e}[/red]")
            return False

    def _save_single_template(self) -> bool:
        """Save single configuration template."""
        try:
            # Get template info
            name = prompt_with_context(
                self.processor,
                "Template name",
                module="AMBER Input Generator",
                description="Custom template name",
            )
            description = prompt_with_context(
                self.processor,
                "Template description",
                default="User-created template",
                module="AMBER Input Generator",
                description="Custom template description",
            )
            
            # Create template data
            template_data = {
                "name": name,
                "description": description,
                "simulation_type": self.simulation_type,
                "config": {k: v for k, v in self.config.items() if not k.startswith("_")},
                "created": datetime.now().isoformat(),
                "version": "1.0"
            }
            
            # Save to file
            safe_name = "".join(c for c in name if c.isalnum() or c in (' ', '-', '_')).rstrip()
            safe_name = safe_name.replace(' ', '_').lower()
            template_file = self.templates_dir / f"{safe_name}.json"
            
            if template_file.exists():
                if not confirm_with_context(
                    self.processor,
                    f"Template {safe_name} exists. Overwrite?",
                    default=False,
                    module="AMBER Input Generator",
                    description=f"Overwrite existing template {safe_name}",
                ):
                    return False
            
            with open(template_file, 'w') as f:
                json.dump(template_data, f, indent=2)
            
            self.console.print(f"[green]Template saved: {template_file}[/green]")
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error saving single template: {e}[/red]")
            return False

    def _save_workflow_template(self) -> bool:
        """Save workflow configuration as preset."""
        try:
            # Get preset info
            name = prompt_with_context(
                self.processor,
                "Workflow preset name",
                module="AMBER Input Generator",
                description="Custom workflow preset name",
            )
            description = prompt_with_context(
                self.processor,
                "Workflow description",
                default="User-created workflow",
                module="AMBER Input Generator",
                description="Custom workflow description",
            )
            
            # Get the final preset from config
            final_preset = self.config.get("_final_preset")
            if not final_preset:
                self.console.print("[yellow]No workflow preset to save[/yellow]")
                return False
            
            # Create preset data
            preset_data = {
                "name": name,
                "description": description,
                "steps": final_preset["steps"],
                "created": datetime.now().isoformat(),
                "version": "1.0"
            }
            
            # Save to file
            safe_name = "".join(c for c in name if c.isalnum() or c in (' ', '-', '_')).rstrip()
            safe_name = safe_name.replace(' ', '_').lower()
            preset_file = self.presets_dir / f"{safe_name}.json"
            
            if preset_file.exists():
                if not confirm_with_context(
                    self.processor,
                    f"Preset {safe_name} exists. Overwrite?",
                    default=False,
                    module="AMBER Input Generator",
                    description=f"Overwrite existing preset {safe_name}",
                ):
                    return False
            
            with open(preset_file, 'w') as f:
                json.dump(preset_data, f, indent=2)
            
            self.console.print(f"[green]Workflow preset saved: {preset_file}[/green]")
            
            # Reload presets
            self.workflow_presets = self._load_all_workflow_presets()
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error saving workflow template: {e}[/red]")
            return False

    # ============================================================================
    # CONFIGURATION MANAGEMENT METHODS
    # ============================================================================

    def _review_configuration(self) -> bool:
        """Review current configuration."""
        try:
            if not self.config:
                self.console.print("[yellow]No configuration set[/yellow]")
                return False

            # Handle workflow vs single configuration differently
            if self.simulation_type == "workflow":
                return self._review_workflow_configuration()
            else:
                return self._review_single_configuration()
                
        except Exception as e:
            self.console.print(f"[red]Error reviewing configuration: {e}[/red]")
            return False

    def _review_workflow_configuration(self) -> bool:
        """Review workflow configuration."""
        try:
            self.console.print("[bold]Workflow Configuration Review[/bold]")
            
            # Basic workflow info
            info_table = Table(title="Workflow Information")
            info_table.add_column("Property", style="cyan")
            info_table.add_column("Value", style="white")
            
            info_table.add_row("Workflow Name", self.config.get("_workflow_name", "Unknown"))
            info_table.add_row("Description", self.config.get("_workflow_description", "N/A"))
            info_table.add_row("Base Preset", self.config.get("_workflow_preset", "Unknown"))
            info_table.add_row("Steps", ", ".join(self.config.get("_workflow_steps", [])))
            
            self.console.print(info_table)
            
            # Show detailed step configuration
            show_details = confirm_with_context(
                self.processor,
                "\nShow detailed step configurations?",
                default=False,
                module="AMBER Input Generator",
                description="Show detailed step configurations",
            )
            if show_details:
                final_preset = self.config.get("_final_preset")
                if final_preset and "steps" in final_preset:
                    for step_name, step_config in final_preset["steps"].items():
                        self.console.print(f"\n[bold]{step_name.replace('_', ' ').title()}:[/bold]")
                        self.console.print(f"  Description: {step_config.get('description', 'N/A')}")
                        self.console.print(f"  Filename: {step_config.get('filename', 'N/A')}")
                        
                        # Show key parameters
                        config = step_config.get("config", {})
                        if config:
                            key_params = ["nstlim", "dt", "temp0", "cut", "ntpr"]
                            for param in key_params:
                                if param in config:
                                    self.console.print(f"  {param}: {config[param]}")
                else:
                    self.console.print("[yellow]No detailed step configuration available[/yellow]")
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error reviewing workflow configuration: {e}[/red]")
            return False

    def _review_single_configuration(self) -> bool:
        """Review single input configuration."""
        try:
            self._display_summary()

            # Show full config
            show_full = confirm_with_context(
                self.processor,
                "Show full configuration?",
                default=False,
                module="AMBER Input Generator",
                description="Show full template configuration",
            )
            if show_full:
                self.console.print("\n[bold]Full Configuration:[/bold]")
                for key, value in sorted(self.config.items()):
                    if not key.startswith("_"):  # Skip internal keys
                        self.console.print(f"  {key} = {value}")

            return True
            
        except Exception as e:
            self.console.print(f"[red]Error reviewing single configuration: {e}[/red]")
            return False

    def _display_summary(self):
        """Display configuration summary."""
        table = Table(title="Configuration Summary")
        table.add_column("Parameter", style="cyan")
        table.add_column("Value", style="green")
        table.add_column("Description", style="white")

        # Key parameters to highlight
        highlights = {
            "Simulation Type": self.simulation_type.title(),
            "Engine": f"{self.engine_choice}" + (" + GPU" if self.config.get("_use_gpu") else ""),
            "System": "Periodic" if self.config.get("ntb", 0) > 0 else "Non-periodic",
        }

        if self.simulation_type != "minimization":
            runtime = self.config.get("nstlim", 0) * self.config.get("dt", 0.002)
            highlights["Runtime"] = f"{runtime:.2f} ps ({runtime/1000:.3f} ns)"

            if self.config.get("temp0"):
                highlights["Temperature"] = f"{self.config['temp0']} K"

            if self.config.get("ntp", 0) > 0:
                highlights["Pressure"] = f"{self.config.get('pres0', 1.0)} bar"

        for param, value in highlights.items():
            table.add_row(param, str(value), "")

        self.console.print(table)

    def _modify_current_settings(self) -> bool:
        """Modify current settings."""
        try:
            if not self.config:
                self.console.print("[yellow]No configuration to modify[/yellow]")
                return False
            
            self.console.print(f"\n[bold cyan]Modify Current Configuration[/bold cyan]")
            self.console.print(f"[bold]Simulation Type:[/bold] {self.simulation_type}")
            self.console.print(f"[bold]Engine:[/bold] {self.engine_choice}")
            
            while True:
                self.console.print(f"\n[bold]Configuration Modification Options[/bold]")
                self.console.print("1. Edit individual parameters", highlight=False)
                self.console.print("2. Add new parameters", highlight=False)
                self.console.print("3. Remove parameters", highlight=False)
                self.console.print("4. Show full configuration", highlight=False)
                self.console.print("5. Restart workflow from beginning", highlight=False)
                self.console.print("6. Return to review", highlight=False)
                
                choice = prompt_with_context(
                    self.processor,
                    "Select option [1-6]",
                    choices=["1", "2", "3", "4", "5", "6"],
                    default="1",
                    module="AMBER Input Generator",
                    description="Template/workflow management action",
                )
                
                if choice == "1":
                    self._edit_existing_parameters()
                elif choice == "2":
                    self._add_new_parameters()
                elif choice == "3":
                    self._remove_parameters()
                elif choice == "4":
                    self._display_current_full_config()
                elif choice == "5":
                    self.config.clear()
                    return self._generate_input_workflow()
                else:
                    return True
                    
        except Exception as e:
            self.console.print(f"[red]Error modifying settings: {e}[/red]")
            return False

    def _edit_existing_parameters(self) -> None:
        """Edit existing parameters."""
        try:
            if not self.config:
                self.console.print("[yellow]No parameters to edit[/yellow]")
                return
            
            # Show current parameters
            table = Table(title="Current Parameters")
            table.add_column("Option", style="cyan")
            table.add_column("Parameter", style="blue")
            table.add_column("Current Value", style="green")
            
            param_list = [k for k in self.config.keys() if not k.startswith("_")]
            for i, param in enumerate(param_list, 1):
                table.add_row(str(i), param, str(self.config[param]))
            
            self.console.print(table)
            
            # Get parameter to edit
            choices = [str(i) for i in range(1, len(param_list) + 1)] + ["0"]
            choice = prompt_with_context(
                self.processor,
                f"Select parameter to edit [1-{len(param_list)}] (0 to cancel)",
                choices=choices,
                module="AMBER Input Generator",
                description="Select parameter to edit",
            )
            
            if choice == "0":
                return
            
            param_name = param_list[int(choice) - 1]
            current_value = self.config[param_name]
            
            self.console.print(f"\n[bold]Editing parameter: {param_name}[/bold]")
            self.console.print(f"Current value: [green]{current_value}[/green]")
            
            # Determine input type based on current value
            try:
                if isinstance(current_value, bool):
                    new_value = confirm_with_context(
                        self.processor,
                        "New value",
                        default=current_value,
                        module="AMBER Input Generator",
                        description=f"New value for {param_name}",
                    )
                elif isinstance(current_value, int):
                    new_value = int_prompt_with_context(
                        self.processor,
                        "New value",
                        default=current_value,
                        module="AMBER Input Generator",
                        description=f"New value for {param_name}",
                    )
                elif isinstance(current_value, float):
                    new_value = float_prompt_with_context(
                        self.processor,
                        "New value",
                        default=current_value,
                        module="AMBER Input Generator",
                        description=f"New value for {param_name}",
                    )
                else:
                    new_value = prompt_with_context(
                        self.processor,
                        "New value",
                        default=str(current_value),
                        module="AMBER Input Generator",
                        description=f"New value for {param_name}",
                    )
                    # Try to convert string numbers back to appropriate type
                    if str(current_value).replace('.', '').replace('-', '').isdigit():
                        try:
                            if '.' in str(current_value):
                                new_value = float(new_value)
                            else:
                                new_value = int(new_value)
                        except ValueError:
                            pass  # Keep as string
            except (ValueError, KeyboardInterrupt):
                self.console.print("[yellow]Edit cancelled[/yellow]")
                return
            
            # Confirm change
            if confirm_with_context(
                self.processor,
                f"Change {param_name} from {current_value} to {new_value}?",
                default=True,
                module="AMBER Input Generator",
                description=f"Confirm change {param_name}",
            ):
                self.config[param_name] = new_value
                self.console.print(f"[green]Parameter {param_name} updated to {new_value}[/green]")
            else:
                self.console.print("[yellow]Change cancelled[/yellow]")
                
        except Exception as e:
            self.console.print(f"[red]Error editing parameters: {e}[/red]")

    def _add_new_parameters(self) -> None:
        """Add new parameters to the configuration."""
        try:
            self.console.print(f"\n[bold]Add New Parameter[/bold]")
            
            param_name = prompt_with_context(
                self.processor,
                "Parameter name",
                module="AMBER Input Generator",
                description="New parameter name",
            ).strip()
            if not param_name:
                self.console.print("[yellow]Invalid parameter name[/yellow]")
                return
            
            if param_name in self.config:
                self.console.print(f"[yellow]Parameter {param_name} already exists. Use edit option to modify it.[/yellow]")
                return
            
            param_value = prompt_with_context(
                self.processor,
                "Parameter value",
                module="AMBER Input Generator",
                description="New parameter value",
            ).strip()
            
            # Try to convert to appropriate type
            converted_value = self._convert_parameter_value(param_value)
            
            self.console.print(f"Adding: {param_name} = {converted_value} (type: {type(converted_value).__name__})")
            
            if confirm_with_context(
                self.processor,
                "Add this parameter?",
                default=True,
                module="AMBER Input Generator",
                description="Confirm add new parameter",
            ):
                self.config[param_name] = converted_value
                self.console.print(f"[green]Parameter {param_name} added[/green]")
            else:
                self.console.print("[yellow]Parameter not added[/yellow]")
                
        except Exception as e:
            self.console.print(f"[red]Error adding parameter: {e}[/red]")

    def _remove_parameters(self) -> None:
        """Remove parameters from the configuration."""
        try:
            param_list = [k for k in self.config.keys() if not k.startswith("_")]
            if not param_list:
                self.console.print("[yellow]No parameters to remove[/yellow]")
                return
            
            # Show current parameters
            table = Table(title="Current Parameters")
            table.add_column("Option", style="cyan")
            table.add_column("Parameter", style="blue")
            table.add_column("Value", style="green")
            
            for i, param in enumerate(param_list, 1):
                table.add_row(str(i), param, str(self.config[param]))
            
            self.console.print(table)
            
            # Get parameter to remove
            choices = [str(i) for i in range(1, len(param_list) + 1)] + ["0"]
            choice = prompt_with_context(
                self.processor,
                f"Select parameter to remove [1-{len(param_list)}] (0 to cancel)",
                choices=choices,
                module="AMBER Input Generator",
                description="Select parameter to remove",
            )
            
            if choice == "0":
                return
            
            param_name = param_list[int(choice) - 1]
            
            if confirm_with_context(
                self.processor,
                f"Remove parameter {param_name}?",
                default=False,
                module="AMBER Input Generator",
                description=f"Confirm remove parameter {param_name}",
            ):
                removed_value = self.config.pop(param_name)
                self.console.print(f"[green]Parameter {param_name} (value: {removed_value}) removed[/green]")
            else:
                self.console.print("[yellow]Removal cancelled[/yellow]")
                
        except Exception as e:
            self.console.print(f"[red]Error removing parameters: {e}[/red]")

    def _display_current_full_config(self) -> None:
        """Display the complete current configuration."""
        try:
            if not self.config:
                self.console.print("[yellow]No configuration loaded[/yellow]")
                return
            
            self.console.print(f"\n[bold cyan]Current Full Configuration[/bold cyan]")
            
            table = Table(title="All Parameters")
            table.add_column("Parameter", style="cyan")
            table.add_column("Value", style="green")
            table.add_column("Type", style="blue")
            
            for param, value in sorted(self.config.items()):
                if not param.startswith("_"):  # Skip internal keys
                    table.add_row(param, str(value), type(value).__name__)
            
            self.console.print(table)
            
            # Show summary info
            param_count = len([k for k in self.config.keys() if not k.startswith("_")])
            self.console.print(f"\n[grey50]Total parameters: {param_count}[/grey50]")
            if self.simulation_type:
                self.console.print(f"[grey50]Simulation type: {self.simulation_type}[/grey50]")
                
        except Exception as e:
            self.console.print(f"[red]Error displaying configuration: {e}[/red]")

    def _convert_parameter_value(self, value_str: str):
        """Convert a string parameter value to the most appropriate type."""
        value_str = value_str.strip()
        
        # Handle boolean values
        if value_str.lower() in ['true', 'yes', 'on', '1']:
            return True
        elif value_str.lower() in ['false', 'no', 'off', '0']:
            return False
        
        # Handle numeric values
        try:
            if '.' in value_str:
                return float(value_str)
            else:
                return int(value_str)
        except ValueError:
            pass
        
        # Handle quoted strings (remove quotes)
        if (value_str.startswith('"') and value_str.endswith('"')) or \
           (value_str.startswith("'") and value_str.endswith("'")):
            return value_str[1:-1]
        
        # Return as string
        return value_str

    def _clear_configuration(self) -> bool:
        """Clear current configuration."""
        try:
            if not self.config:
                self.console.print("[yellow]No configuration to clear[/yellow]")
                return False
            
            if confirm_with_context(
                self.processor,
                "Clear all configuration?",
                default=False,
                module="AMBER Input Generator",
                description="Clear all configuration",
            ):
                self.config.clear()
                self.simulation_type = None
                self.engine_choice = None
                self.console.print("[green]Configuration cleared[/green]")
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error clearing configuration: {e}[/red]")
            return False

    def _export_run_script(self) -> bool:
        """Export run script - discovers files and lets user select."""
        try:
            script_name = prompt_with_context(
                self.processor,
                "Script filename",
                default="run_amber.sh",
                module="AMBER Input Generator",
                description="Run script filename",
            )
            
            if not script_name.endswith(".sh"):
                script_name += ".sh"

            if os.path.exists(script_name):
                if not confirm_with_context(
                    self.processor,
                    f"Script {script_name} exists. Overwrite?",
                    default=False,
                    module="AMBER Input Generator",
                    description=f"Overwrite existing script {script_name}",
                ):
                    return False

            # Discover and select files
            file_selections = self._discover_and_select_files()
            if not file_selections:
                return False
                
            script_content = self._build_file_based_run_script_content(file_selections)

            with open(script_name, 'w') as f:
                f.write(script_content)

            # Make executable
            os.chmod(script_name, 0o755)

            self.console.print(f"[green]✓ Run script generated: {script_name}[/green]")
            self._show_script_summary(file_selections)
            return True

        except Exception as e:
            self.console.print(f"[red]Error exporting script: {e}[/red]")
            return False

    def _discover_and_select_files(self) -> dict:
        """Discover files in directory and let user select what to use."""
        import glob
                
        self.console.print("\n[bold]Discovering files for run script...[/bold]")
        
        file_selections = {}
        
        # 1. Input file (.in)
        input_files = glob.glob("*.in")
        if input_files:
            self.console.print(f"\n[cyan]Found {len(input_files)} AMBER input files:[/cyan]")
            for i, file in enumerate(input_files, 1):
                self.console.print(f"  {i}. {file}")
            
            if len(input_files) == 1:
                file_selections['input'] = input_files[0]
                self.console.print(f"[green]Using: {input_files[0]}[/green]")
            else:
                input_map = {str(i + 1): f for i, f in enumerate(input_files)}
                choice = prompt_with_context(
                    self.processor,
                    "Select input file",
                    choices=[str(i) for i in range(1, len(input_files)+1)],
                    default="1",
                    module="AMBER Input Generator",
                    description="Select .in input file for run script",
                    options_map=input_map,
                )
                choice = remap_recorded_index(self.processor, input_files, str(choice))
                file_selections['input'] = input_files[int(choice)-1]
                annotate_selected_path(self.processor, file_selections['input'])
        else:
            self.console.print("[red]No .in files found in directory[/red]")
            return {}

        # 2. Topology file (.prmtop)
        topo_files = glob.glob("*.prmtop")
        if topo_files:
            self.console.print(f"\n[cyan]Found {len(topo_files)} topology files:[/cyan]")
            for i, file in enumerate(topo_files, 1):
                self.console.print(f"  {i}. {file}")
            
            if len(topo_files) == 1:
                file_selections['prmtop'] = topo_files[0]
                self.console.print(f"[green]Using: {topo_files[0]}[/green]")
            else:
                topo_map = {str(i + 1): f for i, f in enumerate(topo_files)}
                choice = prompt_with_context(
                    self.processor,
                    "Select topology file",
                    choices=[str(i) for i in range(1, len(topo_files)+1)],
                    default="1",
                    module="AMBER Input Generator",
                    description="Select .prmtop topology file for run script",
                    options_map=topo_map,
                )
                choice = remap_recorded_index(self.processor, topo_files, str(choice))
                file_selections['prmtop'] = topo_files[int(choice)-1]
                annotate_selected_path(self.processor, file_selections['prmtop'])
        else:
            self.console.print("[yellow]No .prmtop files found - script will check for generic 'prmtop'[/yellow]")
            file_selections['prmtop'] = 'prmtop'

        # 3. Coordinate file (.rst7, .inpcrd)
        coord_files = glob.glob("*.rst7") + glob.glob("*.inpcrd")
        if coord_files:
            self.console.print(f"\n[cyan]Found {len(coord_files)} coordinate files:[/cyan]")
            for i, file in enumerate(coord_files, 1):
                self.console.print(f"  {i}. {file}")
            
            if len(coord_files) == 1:
                file_selections['inpcrd'] = coord_files[0]
                self.console.print(f"[green]Using: {coord_files[0]}[/green]")
            else:
                coord_map = {str(i + 1): f for i, f in enumerate(coord_files)}
                choice = prompt_with_context(
                    self.processor,
                    "Select coordinate file",
                    choices=[str(i) for i in range(1, len(coord_files)+1)],
                    default="1",
                    module="AMBER Input Generator",
                    description="Select .rst7/.inpcrd coordinate file for run script",
                    options_map=coord_map,
                )
                choice = remap_recorded_index(self.processor, coord_files, str(choice))
                file_selections['inpcrd'] = coord_files[int(choice)-1]
                annotate_selected_path(self.processor, file_selections['inpcrd'])
        else:
            self.console.print("[yellow]No .rst7/.inpcrd files found - script will check for generic 'inpcrd'[/yellow]")
            file_selections['inpcrd'] = 'inpcrd'

        # 4. Reference file (optional, for restraints)
        ref_files = glob.glob("*refc*") + glob.glob("*ref*.rst7")
        if ref_files:
            self.console.print(f"\n[cyan]Found {len(ref_files)} potential reference files:[/cyan]")
            for i, file in enumerate(ref_files, 1):
                self.console.print(f"  {i}. {file}")
            self.console.print(f"  {len(ref_files)+1}. None (no restraints)")
            
            ref_map = {str(i + 1): f for i, f in enumerate(ref_files)}
            ref_map[str(len(ref_files) + 1)] = "None (no restraints)"
            choice = prompt_with_context(
                self.processor,
                "Select reference file (for restraints)",
                choices=[str(i) for i in range(1, len(ref_files)+2)],
                default=str(len(ref_files)+1),
                module="AMBER Input Generator",
                description="Select reference file for restraints",
                options_map=ref_map,
            )
            # Remap a recorded file pick by basename; the 'None' pseudo-option
            # carries no basename and passes through to the <= len guard.
            choice = remap_recorded_index(self.processor, ref_files, str(choice))
            if int(choice) <= len(ref_files):
                file_selections['refc'] = ref_files[int(choice)-1]
                annotate_selected_path(self.processor, file_selections['refc'])

        return file_selections

    def _build_file_based_run_script_content(self, file_selections: dict) -> str:
        """Build run script content based on selected files."""
        lines = [
            "#!/bin/bash",
            "# AMBER Run Script - Generated by ProPrep",
            f"# Simulation type: {self.simulation_type}",
            f"# Input file: {file_selections['input']}",
            "",
            "# File checks"
        ]

        # Add file existence checks
        for file_type, filename in file_selections.items():
            if file_type == 'refc':
                continue  # Handle separately
            lines.extend([
                f"if [ ! -f {filename} ]; then",
                f"    echo 'Error: {filename} not found'",
                f"    exit 1",
                f"fi"
            ])

        # Handle reference file check (optional)
        if 'refc' in file_selections:
            lines.extend([
                "",
                f"if [ ! -f {file_selections['refc']} ]; then",
                f"    echo 'Warning: Reference file {file_selections['refc']} not found'",
                f"    echo 'Restraints may not work properly'",
                f"fi"
            ])

        # Add engine selection and execution
        lines.extend([
            "",
            f"# Engine selection",
            f"ENGINE={self.engine_choice}",
            "",
            f"echo 'Running {self.simulation_type} with $ENGINE...'",
            f"echo 'Input: {file_selections['input']}'",
            f"echo 'Topology: {file_selections['prmtop']}'", 
            f"echo 'Coordinates: {file_selections['inpcrd']}'",
            "",
            f"# Run simulation",
            f"$ENGINE -O -i {file_selections['input']} -o mdout -p {file_selections['prmtop']} -c {file_selections['inpcrd']} -r restrt"
        ])

        # Add reference file to command if needed
        if 'refc' in file_selections:
            lines[-1] += f" -ref {file_selections['refc']}"

        lines.extend([
            "",
            "if [ $? -eq 0 ]; then",
            "    echo 'Simulation completed successfully'",
            "else",
            "    echo 'Simulation failed'",
            "    exit 1",
            "fi"
        ])

        return "\n".join(lines)

    def _show_script_summary(self, file_selections: dict):
        """Show summary of generated script."""
        self.console.print("\n[bold]Script Summary:[/bold]")
        self.console.print(f"[cyan]Engine:[/cyan] {self.engine_choice}")
        for file_type, filename in file_selections.items():
            type_name = {'input': 'Input file', 'prmtop': 'Topology', 'inpcrd': 'Coordinates', 'refc': 'Reference'}.get(file_type, file_type)
            self.console.print(f"[cyan]{type_name}:[/cyan] {filename}")
        self.console.print("[grey50]Make the script executable with: chmod +x run_amber.sh[/grey50]")

    def _build_single_run_script_content(self) -> str:
        """Build run script content for single simulation."""
        lines = [
            "#!/bin/bash",
            "# AMBER Run Script",
            f"# Generated by ProPrep AMBER Input Generator",
            f"# Simulation type: {self.simulation_type}",
            "",
            "# Check for required files",
            "if [ ! -f prmtop ]; then",
            "    echo 'Error: prmtop file not found'",
            "    exit 1",
            "fi",
            "",
            "if [ ! -f inpcrd ]; then",
            "    echo 'Error: inpcrd file not found'",
            "    exit 1",
            "fi",
            "",
        ]

        # Add restart file check if needed
        if self.config.get("irest", 0) == 1:
            lines.extend([
                "if [ ! -f restrt ]; then",
                "    echo 'Error: restrt file not found (required for restart)'",
                "    exit 1",
                "fi",
                ""
            ])

        # Add reference file check if needed
        if self.config.get("ntr", 0) == 1:
            lines.extend([
                "if [ ! -f refc ]; then",
                "    echo 'Error: refc file not found (required for restraints)'",
                "    exit 1",
                "fi",
                ""
            ])

        # Set engine
        lines.extend([
            f"# Set engine",
            f"ENGINE={self.engine_choice}",
            "",
            f"echo 'Running {self.simulation_type} with $ENGINE...'",
            ""
        ])

        # Build command
        mdin_file = f"{self.simulation_type}.in"
        cmd_parts = [
            f"$ENGINE -i {mdin_file}",
            f"-o {self.simulation_type}.out",
            "-p prmtop",
            "-c inpcrd",
            f"-r {self.simulation_type}.rst"
        ]

        if self.config.get("ntwx", 0) > 0:
            cmd_parts.append(f"-x {self.simulation_type}.nc")

        if self.config.get("ntr", 0) == 1:
            cmd_parts.append("-ref refc")

        command = " \\\n    ".join(cmd_parts)

        lines.extend([
            command,
            "",
            "if [ $? -eq 0 ]; then",
            f"    echo '{self.simulation_type.title()} completed successfully!'",
            "else",
            f"    echo 'Error in {self.simulation_type}'",
            "    exit 1",
            "fi"
        ])

        return "\n".join(lines) + "\n"