"""
AMBER Configuration Wizard - Standalone Module

Comprehensive 13-step wizard for AMBER MD parameter configuration.
This module serves as a template configuration assistant that handles
all simulation parameters. Hardware/engine selection and restraints are
handled separately by the MD Manager.
"""

import os
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import json
import copy

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich import box

from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
    float_prompt_with_context,
)

from proprep.md_prep.amber_parameter_database import build_parameter_database, Parameter
from proprep.md_prep.amber_parameter_display import (
    Verbosity,
    format_parameter_hint,
    print_reference_guide,
)


class AmberWizard:
    """
    Comprehensive AMBER parameter configuration wizard.
    
    Runs the full 15-step parameter selection process and returns
    a complete configuration dictionary suitable for template population.
    """
    
    def __init__(self, simulation_type: str = None, initial_config: Dict[str, Any] = None,
                 console: Console = None, verbosity: Verbosity = None, processor=None):
        """
        Initialize the AMBER wizard.

        Args:
            simulation_type: Target simulation type (minimization, md, heating, etc.)
            initial_config: Starting configuration to build upon
            console: Rich console for output
            verbosity: Educational content level (None = prompt at start)
            processor: Optional ProPrep processor for session recording context
        """
        self.simulation_type = simulation_type
        self.config = initial_config.copy() if initial_config else {}
        self.console = console or Console()
        self.verbosity = verbosity
        self.processor = processor

        # Educational parameter database
        self.param_db = build_parameter_database()

        # Store original config for comparison
        self.initial_config = initial_config.copy() if initial_config else {}
        
    @classmethod
    def configure(cls, simulation_type: str = None, initial_config: Dict[str, Any] = None,
                  console: Console = None, verbosity: Verbosity = None, processor=None) -> Dict[str, Any]:
        """
        Factory method to run complete wizard configuration.

        Args:
            simulation_type: Target simulation type (None = prompt in Step 1)
            initial_config: Optional starting configuration
            console: Rich console for output
            verbosity: Educational content level (None = prompt at start)
            processor: Optional ProPrep processor for session recording context

        Returns:
            Complete AMBER configuration dictionary (includes _simulation_type key)
        """
        wizard = cls(simulation_type, initial_config, console, verbosity, processor=processor)
        return wizard.run_complete_workflow()
        
    def run_complete_workflow(self) -> Dict[str, Any]:
        """
        Execute the comprehensive 13-step AMBER configuration wizard.

        Note: Hardware/engine configuration and restraints are handled by the
        MD Manager (Step 4: Hardware and Step 3: Restraints respectively).

        Returns:
            Complete configuration dictionary
        """
        try:
            self.console.print("\n[bold cyan]===== AMBER Configuration Wizard =====[/bold cyan]")

            if self.simulation_type:
                self.console.print(f"[green]Configuring: {self.simulation_type} simulation[/green]")
            else:
                self.console.print("Let's configure your MD simulation step by step.")

            # Display workflow overview
            overview = (
                "[bold]13-Step Template Configuration:[/bold]\n\n"
                "[cyan] 1.[/cyan] Simulation Type\n"
                "[cyan] 2.[/cyan] Performance Options\n"
                "[cyan] 3.[/cyan] System Setup\n"
                "[cyan] 4.[/cyan] Input/Restart\n"
                "[cyan] 5.[/cyan] Simulation Parameters\n"
                "[cyan] 6.[/cyan] Constraints\n"
                "[cyan] 7.[/cyan] Time Control\n"
                "[cyan] 8.[/cyan] Temperature Control\n"
                "[cyan] 9.[/cyan] Pressure Control\n"
                "[cyan]10.[/cyan] Non-bonded Interactions\n"
                "[cyan]11.[/cyan] Output Control\n"
                "[cyan]12.[/cyan] Advanced Options\n"
                "[cyan]13.[/cyan] Review & Summary\n\n"
                "[grey50]Note: Hardware/engine and restraints are configured\n"
                "separately in the MD Manager workflow.[/grey50]"
            )
            self.console.print(Panel(overview, title="Workflow Overview", expand=False, border_style="cyan"))

            # Educational content verbosity selection
            if self.verbosity is None:
                self.verbosity = self._select_verbosity()

            # Offer reference guide before starting
            if confirm_with_context(
                self.processor,
                "View parameter reference guide before starting?",
                default=False,
                module="AMBER Wizard",
                description="Show the AMBER parameter reference guide before the wizard starts?",
            ):
                print_reference_guide(self.param_db, self.console)

            # Step 1: Simulation Type (always show, even if pre-specified)
            if not self._select_simulation_type():
                return self.config

            # Step 2: Performance Options
            if not self._configure_performance():
                return self.config

            # Step 3: System Setup (periodic vs non-periodic)
            if not self._configure_system_type():
                return self.config

            # Step 4: Input/Restart Configuration
            if not self._configure_input_restart():
                return self.config

            # Step 5: Simulation Parameters (context-dependent)
            if not self._configure_simulation_parameters():
                return self.config

            # Step 6: Constraints
            if not self._configure_constraints():
                return self.config

            # Step 7: Time Control
            if not self._configure_time_control():
                return self.config

            # Step 8: Temperature Control
            if not self._configure_temperature():
                return self.config

            # Step 9: Pressure Control (if needed)
            if not self._configure_pressure():
                return self.config

            # Step 10: Non-bonded Interactions
            if not self._configure_nonbonded():
                return self.config

            # Step 11: Output Control
            if not self._configure_output():
                return self.config

            # Step 12: Advanced Options
            if not self._configure_advanced_options():
                return self.config

            # Step 13: Final review
            self._show_configuration_summary()

            # Include simulation type in config for callers that need it
            self.config["_simulation_type"] = self.simulation_type

            return self.config

        except KeyboardInterrupt:
            self.console.print("\n[yellow]Configuration cancelled by user[/yellow]")
            return self.initial_config
        except Exception as e:
            self.console.print(f"\n[red]Error in wizard: {e}[/red]")
            return self.initial_config
            
    def _select_verbosity(self) -> Verbosity:
        """Prompt user for educational content verbosity level."""
        self.console.print("\n[bold cyan]Educational Content Level:[/bold cyan]")
        self.console.print("[cyan]1.[/cyan] Terse   - Parameter name, options, and recommended value")
        self.console.print("[cyan]2.[/cyan] Compact - Terse + what the parameter physically does")
        self.console.print("[cyan]3.[/cyan] Verbose - Compact + why the default, when to change")

        choice = prompt_with_context(
            self.processor,
            "Select verbosity",
            choices=["1", "2", "3"],
            default="2",
            module="AMBER Wizard",
            description="Educational content verbosity level for parameter hints",
            options_map={
                "1": "Terse (name, options, recommended value)",
                "2": "Compact (adds physical meaning)",
                "3": "Verbose (adds rationale and when to change)",
            },
        )
        return {
            "1": Verbosity.TERSE,
            "2": Verbosity.COMPACT,
            "3": Verbosity.VERBOSE,
        }[choice]

    def _show_param_hint(self, param_name: str, recommended: Any = None) -> None:
        """Display educational content for a parameter if it exists in the database."""
        param = self.param_db.get(param_name)
        if param:
            format_parameter_hint(param, self.verbosity, self.console, recommended)

    def _select_simulation_type(self) -> bool:
        """Step 1: Simulation type selection with all imin options."""
        self.console.print("\n[bold]Step 1: Simulation Type Selection[/bold]")

        self._show_param_hint("imin")

        options = {
            "1": ("minimization", "Energy Minimization (imin=1)", 1),
            "2": ("md", "Molecular Dynamics (imin=0)", 0),
            "3": ("trajectory_analysis", "Trajectory Analysis (imin=5)", 5),
            "4": ("md_analysis", "MD Analysis of Trajectory (imin=6)", 6),
            "5": ("socket_server", "Socket Server Mode (imin=7)", 7),
            "6": ("custom", "Custom Setup", None),
        }

        # Determine default based on pre-selection (if any)
        default_choice = "1"
        if self.simulation_type:
            # Map simulation_type to menu choice
            type_to_choice = {
                "minimization": "1",
                "md": "2", "heating": "2", "equilibration": "2", "production": "2",
                "trajectory_analysis": "3",
                "md_analysis": "4",
                "socket_server": "5",
                "custom": "6",
            }
            default_choice = type_to_choice.get(self.simulation_type, "1")

        for key, (_, desc, _) in options.items():
            self.console.print(f"[cyan]{key}.[/cyan] {desc}")

        choice = prompt_with_context(
            self.processor,
            "Select simulation type",
            choices=list(options.keys()),
            default=default_choice,
            module="AMBER Wizard",
            description="Step 1: top-level AMBER simulation mode (imin)",
            options_map={k: desc for k, (_, desc, _) in options.items()},
        )

        self.simulation_type, description, imin_val = options[choice]

        # Set imin based on choice
        if imin_val is not None:
            self.config["imin"] = imin_val
        else:
            imin = int_prompt_with_context(
                self.processor,
                "Enter imin value",
                default=0,
                module="AMBER Wizard",
                description="Custom imin value (AMBER integration mode)",
            )
            self.config["imin"] = imin

        self.console.print(f"[cyan]imin={self.config['imin']} ({self.simulation_type})[/cyan]")

        return True

    def _configure_performance(self) -> bool:
        """Step 2: Performance Options."""
        self.console.print("\n[bold]Step 2: Performance Options[/bold]")

        # NRESPA multiple time stepping
        self.console.print("\n[bold cyan]2a.[/bold cyan] Multiple time stepping:")
        self._show_param_hint("nrespa")
        use_respa = confirm_with_context(
            self.processor,
            "Use multiple time stepping (NRESPA)?",
            default=False,
            module="AMBER Wizard",
            description="Step 2a: enable NRESPA multiple time stepping for slow forces",
        )
        if use_respa:
            nrespa = int_prompt_with_context(
                self.processor,
                "NRESPA factor",
                default=2,
                module="AMBER Wizard",
                description="NRESPA factor: slow forces evaluated every N steps",
            )
            self.config["nrespa"] = nrespa
            self.console.print(f"[cyan]ℹ️  Slow forces evaluated every {nrespa} steps[/cyan]")

        # Center of mass motion removal
        if self.simulation_type not in ["minimization", "trajectory_analysis"]:
            self.console.print("\n[bold cyan]2b.[/bold cyan] Center-of-mass motion removal:")
            self._show_param_hint("nscm")
            nscm = int_prompt_with_context(
                self.processor,
                "Center-of-mass motion removal frequency (nscm)",
                default=1000,
                module="AMBER Wizard",
                description="Step 2b: nscm — COM motion removal frequency (0 disables)",
            )
            self.config["nscm"] = nscm
            if nscm == 0:
                self.console.print("[yellow]⚠️  COM motion removal disabled[/yellow]")
                
        return True
        
    def _configure_system_type(self) -> bool:
        """Step 3: System setup and boundary conditions."""
        self.console.print("\n[bold]Step 3: System Setup[/bold]")
        self.console.print("\n[bold cyan]3a.[/bold cyan] Boundary conditions:")
        self._show_param_hint("ntb", recommended=2)
        self.console.print("[cyan]0.[/cyan] No periodicity (ntb=0)")
        self.console.print("[cyan]1.[/cyan] Constant volume periodic (ntb=1)")
        self.console.print("[cyan]2.[/cyan] Constant pressure periodic (ntb=2)")

        ntb = prompt_with_context(
            self.processor,
            "Select boundary condition",
            choices=["0", "1", "2"],
            default="2",
            module="AMBER Wizard",
            description="Step 3a: periodic boundary conditions (ntb)",
            options_map={
                "0": "No periodicity (ntb=0)",
                "1": "Constant volume periodic (ntb=1)",
                "2": "Constant pressure periodic (ntb=2)",
            },
        )
        self.config["ntb"] = int(ntb)

        if ntb == "0":
            self.console.print("[yellow]⚠️  Non-periodic system - pressure control disabled[/yellow]")
            self.config["ntp"] = 0  # No pressure control

            # Offer implicit solvent for non-periodic systems
            self.console.print("\n[bold cyan]3b.[/bold cyan] Implicit solvent (Generalized Born):")
            use_gb = confirm_with_context(
                self.processor,
                "Use implicit solvent (GB)?",
                default=False,
                module="AMBER Wizard",
                description="Step 3b: enable Generalized Born implicit solvent for non-periodic system",
            )

            if use_gb:
                self._configure_implicit_solvent()
        elif ntb == "1":
            self.console.print("[cyan]ℹ️  Constant volume - consider pressure equilibration first[/cyan]")
        else:  # ntb == "2"
            self.console.print("[green]✓ Constant pressure - good for production MD[/green]")

        # Water cap option (for non-periodic or special setups)
        if ntb == "0":
            self.console.print("\n[bold cyan]3c.[/bold cyan] Water cap (spherical boundary):")
            use_cap = confirm_with_context(
                self.processor,
                "Use water cap?",
                default=False,
                module="AMBER Wizard",
                description="Step 3c: enable spherical water cap boundary",
            )
            if use_cap:
                self._configure_water_cap()

        return True

    def _configure_implicit_solvent(self) -> bool:
        """Configure Generalized Born implicit solvent options."""
        self.console.print("\n[bold cyan]Implicit Solvent Configuration[/bold cyan]")

        self._show_param_hint("igb", recommended=2)
        self.console.print("\nGeneralized Born model selection:")
        self.console.print("[cyan]1.[/cyan] HCT (igb=1) - Hawkins, Cramer, Truhlar")
        self.console.print("[cyan]2.[/cyan] OBC I (igb=2) - Onufriev, Bashford, Case (recommended)")
        self.console.print("[cyan]5.[/cyan] OBC II (igb=5) - Modified OBC")
        self.console.print("[cyan]6.[/cyan] Vacuum (igb=6) - No solvent")
        self.console.print("[cyan]7.[/cyan] GBn (igb=7) - GBneck")
        self.console.print("[cyan]8.[/cyan] GBn2 (igb=8) - GBneck2 (recommended for nucleic acids)")

        igb_choice = prompt_with_context(
            self.processor,
            "Select GB model",
            choices=["1", "2", "5", "6", "7", "8"],
            default="2",
            module="AMBER Wizard",
            description="Generalized Born implicit solvent model (igb)",
            options_map={
                "1": "HCT (igb=1)",
                "2": "OBC I (igb=2, recommended)",
                "5": "OBC II (igb=5)",
                "6": "Vacuum (igb=6)",
                "7": "GBn (igb=7)",
                "8": "GBn2 (igb=8, nucleic acids)",
            },
        )
        igb_map = {"1": 1, "2": 2, "5": 5, "6": 6, "7": 7, "8": 8}
        self.config["igb"] = igb_map[igb_choice]

        if igb_choice != "6":  # Not vacuum
            # Salt concentration
            self._show_param_hint("saltcon")
            saltcon = float_prompt_with_context(
                self.processor,
                "Salt concentration (saltcon) [M]",
                default=0.0,
                module="AMBER Wizard",
                description="GB salt concentration (saltcon) in Molar",
            )
            if saltcon > 0:
                self.config["saltcon"] = saltcon

            # GB interaction cutoff
            self._show_param_hint("rgbmax")
            rgbmax = float_prompt_with_context(
                self.processor,
                "GB interaction cutoff (rgbmax) [Å]",
                default=25.0,
                module="AMBER Wizard",
                description="GB Born radius calculation cutoff (rgbmax, Å)",
            )
            if rgbmax != 25.0:
                self.config["rgbmax"] = rgbmax

            # Surface area term for GBSA
            self._show_param_hint("gbsa")
            self.console.print("\nSurface area (SA) term:")
            self.console.print("[cyan]0.[/cyan] No SA term")
            self.console.print("[cyan]1.[/cyan] LCPO surface area")
            self.console.print("[cyan]2.[/cyan] Recursive surface area")

            gbsa = prompt_with_context(
                self.processor,
                "Include SA term?",
                choices=["0", "1", "2"],
                default="0",
                module="AMBER Wizard",
                description="GB surface area term (gbsa)",
                options_map={
                    "0": "No SA term",
                    "1": "LCPO surface area",
                    "2": "Recursive surface area",
                },
            )
            self.config["gbsa"] = int(gbsa)

            if gbsa != "0":
                self._show_param_hint("surften")
                surften = float_prompt_with_context(
                    self.processor,
                    "Surface tension (surften) [cal/mol·Å²]",
                    default=0.005,
                    module="AMBER Wizard",
                    description="GBSA surface tension (surften, cal/mol·Å²)",
                )
                self.config["surften"] = surften

        self.console.print(f"[green]✓ Configured GB model igb={self.config['igb']}[/green]")
        return True

    def _configure_water_cap(self) -> bool:
        """Configure water cap (spherical boundary) options."""
        self.console.print("\n[bold cyan]Water Cap Configuration[/bold cyan]")

        # Cap type
        self._show_param_hint("ivcap")
        self.console.print("\nCap type:")
        self.console.print("[cyan]1.[/cyan] Water cap with restraining potential")
        self.console.print("[cyan]2.[/cyan] Orthorhombic virtual box")

        ivcap = prompt_with_context(
            self.processor,
            "Select cap type",
            choices=["1", "2"],
            default="1",
            module="AMBER Wizard",
            description="Water cap type (ivcap)",
            options_map={
                "1": "Water cap with restraining potential",
                "2": "Orthorhombic virtual box",
            },
        )
        self.config["ivcap"] = int(ivcap)

        # Cap force constant
        self._show_param_hint("fcap")
        fcap = float_prompt_with_context(
            self.processor,
            "Cap force constant (fcap) [kcal/mol·Å²]",
            default=1.5,
            module="AMBER Wizard",
            description="Water cap restraining force constant (fcap)",
        )
        self.config["fcap"] = fcap

        # Cap radius
        self._show_param_hint("cutcap")
        cutcap = float_prompt_with_context(
            self.processor,
            "Cap radius (cutcap) [Å]",
            default=15.0,
            module="AMBER Wizard",
            description="Water cap radius (cutcap, Å)",
        )
        self.config["cutcap"] = cutcap

        # Cap center
        self.console.print("\nCap center coordinates (leave as 0 for center of mass):")
        xcap = float_prompt_with_context(
            self.processor,
            "X center (xcap)",
            default=0.0,
            module="AMBER Wizard",
            description="Water cap center X coordinate (0 = center of mass)",
        )
        ycap = float_prompt_with_context(
            self.processor,
            "Y center (ycap)",
            default=0.0,
            module="AMBER Wizard",
            description="Water cap center Y coordinate (0 = center of mass)",
        )
        zcap = float_prompt_with_context(
            self.processor,
            "Z center (zcap)",
            default=0.0,
            module="AMBER Wizard",
            description="Water cap center Z coordinate (0 = center of mass)",
        )

        if xcap != 0.0:
            self.config["xcap"] = xcap
        if ycap != 0.0:
            self.config["ycap"] = ycap
        if zcap != 0.0:
            self.config["zcap"] = zcap

        self.console.print(f"[green]✓ Configured water cap with radius {cutcap} Å[/green]")
        return True
        
    def _configure_input_restart(self) -> bool:
        """Step 4: Input and restart configuration."""
        self.console.print("\n[bold]Step 4: Input/Restart Configuration[/bold]")
        
        if self.simulation_type == "minimization":
            self.config["ntx"] = 1
            self.config["irest"] = 0
            self.console.print("[cyan]ℹ️  Minimization: ntx=1, irest=0 (coordinates only)[/cyan]")
            return True
            
        self.console.print("\n[bold cyan]4a.[/bold cyan] Input coordinate reading:")
        self._show_param_hint("ntx", recommended=5)
        self.console.print("[cyan]1.[/cyan] Coordinates only (ntx=1)")
        self.console.print("[cyan]5.[/cyan] Coordinates and velocities (ntx=5)")
        
        ntx = prompt_with_context(
            self.processor,
            "Select input type",
            choices=["1", "5"],
            default="5",
            module="AMBER Wizard",
            description="Step 4a: input coordinate reading mode (ntx)",
            options_map={
                "1": "Coordinates only (ntx=1)",
                "5": "Coordinates and velocities (ntx=5)",
            },
        )
        self.config["ntx"] = int(ntx)
        
        self.console.print("\n[bold cyan]4b.[/bold cyan] Restart flag:")
        self._show_param_hint("irest")
        if ntx == "5":
            irest = confirm_with_context(
                self.processor,
                "Restart from previous simulation (irest=1)?",
                default=True,
                module="AMBER Wizard",
                description="Step 4b: restart flag (irest) — continue from prior simulation",
            )
            self.config["irest"] = 1 if irest else 0
        else:
            self.config["irest"] = 0
            self.console.print("[cyan]ℹ️  New simulation: irest=0[/cyan]")
            
        return True
        
    def _configure_simulation_parameters(self) -> bool:
        """Step 5: Simulation-specific parameters (context-dependent)."""
        if self.simulation_type == "minimization":
            self.console.print("\n[bold]Step 5: Simulation Parameters[/bold]")
            return self._configure_minimization_complete()
        elif self.simulation_type in ["trajectory_analysis", "md_analysis"]:
            self.console.print("\n[bold]Step 5: Simulation Parameters[/bold]")
            return self._configure_analysis_mode()
        elif self.simulation_type == "socket_server":
            self.console.print("\n[bold]Step 5: Simulation Parameters[/bold]")
            return self._configure_socket_mode()
        else:  # md, heating, custom
            self.console.print("\n[bold]Step 5: Simulation Parameters[/bold]")
            self._show_param_hint("imin", recommended=0)
            self.config["imin"] = 0
            self.console.print(f"[cyan]imin=0 (molecular dynamics, {self.simulation_type})[/cyan]")
            return True
            
    def _configure_minimization_complete(self) -> bool:
        """Complete minimization configuration."""
        # Maximum cycles
        self.console.print("\n[bold cyan]5a.[/bold cyan] Minimization cycles:")
        self._show_param_hint("maxcyc")
        maxcyc = int_prompt_with_context(
            self.processor,
            "Maximum minimization cycles (maxcyc)",
            default=10000,
            module="AMBER Wizard",
            description="Step 5a: maximum minimization cycles (maxcyc)",
        )
        self.config["maxcyc"] = maxcyc

        # Method selection
        self.console.print("\n[bold cyan]5b.[/bold cyan] Minimization method:")
        self._show_param_hint("ntmin", recommended=1)
        self.console.print("[cyan]0.[/cyan] Full conjugate gradient")
        self.console.print("[cyan]1.[/cyan] Steepest descent → conjugate gradient (recommended)")
        self.console.print("[cyan]2.[/cyan] Steepest descent only")
        self.console.print("[cyan]3.[/cyan] XMIN method")
        self.console.print("[cyan]4.[/cyan] LMOD method")
        self.console.print("[cyan]5.[/cyan] TNCG (Truncated Newton conjugate gradient)")

        ntmin = prompt_with_context(
            self.processor,
            "Select method",
            choices=["0", "1", "2", "3", "4", "5"],
            default="1",
            module="AMBER Wizard",
            description="Step 5b: minimization method (ntmin)",
            options_map={
                "0": "Full conjugate gradient",
                "1": "Steepest descent → conjugate gradient (recommended)",
                "2": "Steepest descent only",
                "3": "XMIN method",
                "4": "LMOD method",
                "5": "TNCG (Truncated Newton CG)",
            },
        )
        self.config["ntmin"] = int(ntmin)

        if ntmin == "1":
            self.console.print("\n[bold cyan]5c.[/bold cyan] Steepest descent configuration:")
            self._show_param_hint("ncyc")
            ncyc = int_prompt_with_context(
                self.processor,
                "Steepest descent cycles before switching (ncyc)",
                default=1000,
                module="AMBER Wizard",
                description="Step 5c: ncyc — steepest descent cycles before switching to CG",
            )
            self.config["ncyc"] = ncyc

        # Convergence criteria
        self.console.print("\n[bold cyan]5d.[/bold cyan] Convergence criteria:")
        self._show_param_hint("drms")
        drms = float_prompt_with_context(
            self.processor,
            "Convergence criterion (drms) [kcal/mol/Å]",
            default=0.0001,
            module="AMBER Wizard",
            description="Step 5d: drms — minimization convergence criterion",
        )
        self.config["drms"] = drms

        # Initial step size
        self.console.print("\n[bold cyan]5e.[/bold cyan] Step size:")
        self._show_param_hint("dx0")
        dx0 = float_prompt_with_context(
            self.processor,
            "Initial step size (dx0)",
            default=0.01,
            module="AMBER Wizard",
            description="Step 5e: dx0 — initial minimizer step size",
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
                module="AMBER Wizard",
                description="Trajectory analysis: single-point energies vs. re-minimization per frame",
            )
            if single_point:
                self.config["maxcyc"] = 1
            else:
                maxcyc = int_prompt_with_context(
                    self.processor,
                    "Maximum cycles per frame",
                    default=100,
                    module="AMBER Wizard",
                    description="Trajectory analysis: maximum minimization cycles per frame",
                )
                self.config["maxcyc"] = maxcyc

        elif self.simulation_type == "md_analysis":
            self.console.print("[cyan]ℹ️  MD analysis mode - runs MD from each trajectory frame[/cyan]")
            
            nstlim = int_prompt_with_context(
                self.processor,
                "MD steps per frame (0 for single point)",
                default=0,
                module="AMBER Wizard",
                description="MD analysis: nstlim MD steps per trajectory frame",
            )
            self.config["nstlim"] = nstlim

            if nstlim > 0:
                dt = float_prompt_with_context(
                    self.processor,
                    "Time step [ps]",
                    default=0.002,
                    module="AMBER Wizard",
                    description="MD analysis: integration time step (dt, ps)",
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
        
    def _configure_time_control(self) -> bool:
        """Step 7: Time step and simulation length."""
        self.console.print("\n[bold]Step 7: Time Control[/bold]")
        
        if self.simulation_type == "minimization":
            self.console.print("[cyan]ℹ️  No time control needed for minimization[/cyan]")
            return True
            
        self.console.print("\n[bold cyan]7a.[/bold cyan] Time step:")
        self._show_param_hint("dt", recommended=0.002)
        if self.config.get("ntc", 1) >= 2:  # SHAKE constraints
            dt_default = 0.002
            self.console.print("[cyan]ℹ️  SHAKE constraints allow 2 fs time step[/cyan]")
        else:
            dt_default = 0.001
            self.console.print("[yellow]⚠️  No constraints - consider smaller time step[/yellow]")
            
        dt = float_prompt_with_context(
            self.processor,
            "Time step (dt) [ps]",
            default=dt_default,
            module="AMBER Wizard",
            description="Step 7a: integration time step (dt, ps)",
        )
        self.config["dt"] = dt
        
        if dt > 0.002:
            self.console.print("[yellow]⚠️  Large time step - ensure appropriate constraints[/yellow]")
            
        self.console.print("\n[bold cyan]7b.[/bold cyan] Simulation length:")
        self._show_param_hint("nstlim")

        # Suggest based on simulation type
        if self.simulation_type == "heating":
            default_steps = 50000  # 100 ps for heating
        elif "equil" in self.simulation_type.lower():
            default_steps = 100000  # 200 ps for equilibration
        else:
            default_steps = 2500000  # 5 ns for production
            
        nstlim = int_prompt_with_context(
            self.processor,
            "Number of steps (nstlim)",
            default=default_steps,
            module="AMBER Wizard",
            description="Step 7b: total number of MD steps (nstlim)",
        )
        self.config["nstlim"] = nstlim
        
        # Calculate and display total time
        total_time = nstlim * dt
        if total_time < 1:
            time_str = f"{total_time*1000:.1f} ps"
        else:
            time_str = f"{total_time:.2f} ns"
            
        self.console.print(f"[green]Total simulation time: {time_str}[/green]")
        
        # Start time (optional)
        self.console.print("\n[bold cyan]7c.[/bold cyan] Start time (optional):")
        if confirm_with_context(
            self.processor,
            "Set custom start time?",
            default=False,
            module="AMBER Wizard",
            description="Step 7c: override default start time (t)?",
        ):
            t = float_prompt_with_context(
                self.processor,
                "Start time (t) [ps]",
                default=0.0,
                module="AMBER Wizard",
                description="Simulation start time (t, ps)",
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
        self._show_param_hint("ntt", recommended=3)
        self.console.print("[cyan]0.[/cyan] No temperature control (NVE)")
        self.console.print("[cyan]1.[/cyan] Weak coupling (Berendsen, not recommended)")
        self.console.print("[cyan]2.[/cyan] Andersen thermostat")
        self.console.print("[cyan]3.[/cyan] Langevin dynamics (recommended)")
        self.console.print("[cyan]9.[/cyan] Optimized isokinetic Nose-Hoover (RESPA compatible)")
        self.console.print("[cyan]10.[/cyan] Stochastic isokinetic Nose-Hoover")
        self.console.print("[cyan]11.[/cyan] Bussi velocity rescaling")

        ntt = prompt_with_context(
            self.processor,
            "Select thermostat",
            choices=["0", "1", "2", "3", "9", "10", "11"],
            default="3",
            module="AMBER Wizard",
            description="Step 8a: thermostat type (ntt)",
            options_map={
                "0": "No temperature control (NVE)",
                "1": "Weak coupling (Berendsen)",
                "2": "Andersen thermostat",
                "3": "Langevin dynamics (recommended)",
                "9": "Optimized isokinetic Nose-Hoover (RESPA)",
                "10": "Stochastic isokinetic Nose-Hoover",
                "11": "Bussi velocity rescaling",
            },
        )
        self.config["ntt"] = int(ntt)
        
        if ntt != "0":
            # Target temperature
            self.console.print("\n[bold cyan]8b.[/bold cyan] Target temperature:")
            self._show_param_hint("temp0")
            temp0 = float_prompt_with_context(
                self.processor,
                "Target temperature (temp0) [K]",
                default=300.0,
                module="AMBER Wizard",
                description="Step 8b: target temperature (temp0, K)",
            )
            self.config["temp0"] = temp0
            
            if temp0 > 350:
                self.console.print("[yellow]⚠️  High temperature - consider reducing time step[/yellow]")
                
            # Initial temperature for new simulations
            if self.config.get("irest", 0) == 0:
                self.console.print("\n[bold cyan]8c.[/bold cyan] Initial temperature:")
                self._show_param_hint("tempi")
                if self.simulation_type == "heating":
                    tempi = float_prompt_with_context(
                        self.processor,
                        "Initial temperature (tempi) [K]",
                        default=0.0,
                        module="AMBER Wizard",
                        description="Step 8c: initial temperature (tempi, K) — heating protocol",
                    )
                else:
                    tempi = float_prompt_with_context(
                        self.processor,
                        "Initial temperature (tempi) [K]",
                        default=temp0,
                        module="AMBER Wizard",
                        description="Step 8c: initial temperature (tempi, K)",
                    )
                self.config["tempi"] = tempi
                
            # Thermostat-specific parameters
            self.console.print("\n[bold cyan]8d.[/bold cyan] Thermostat parameters:")
            
            if ntt == "1":
                self._show_param_hint("tautp")
                tautp = float_prompt_with_context(
                    self.processor,
                    "Heat bath time constant (tautp) [ps]",
                    default=1.0,
                    module="AMBER Wizard",
                    description="Berendsen heat bath coupling time (tautp, ps)",
                )
                self.config["tautp"] = tautp
                self.console.print("[yellow]⚠️  Weak coupling can cause problems - consider Langevin[/yellow]")
                
            elif ntt == "2":
                self._show_param_hint("vrand")
                vrand = int_prompt_with_context(
                    self.processor,
                    "Collision frequency (vrand steps)",
                    default=1000,
                    module="AMBER Wizard",
                    description="Andersen thermostat velocity randomization frequency (vrand)",
                )
                self.config["vrand"] = vrand
                
            elif ntt == "3":
                self._show_param_hint("gamma_ln", recommended=2.0)
                gamma_ln = float_prompt_with_context(
                    self.processor,
                    "Collision frequency (gamma_ln) [ps⁻¹]",
                    default=2.0,
                    module="AMBER Wizard",
                    description="Langevin collision frequency (gamma_ln, ps⁻¹)",
                )
                self.config["gamma_ln"] = gamma_ln
                if gamma_ln < 0.1:
                    self.console.print("[cyan]ℹ️  Low collision frequency - good for enhanced sampling[/cyan]")
                    
            elif ntt == "9":
                # Optimized isokinetic Nose-Hoover
                self._show_param_hint("nkija")
                nkija = int_prompt_with_context(
                    self.processor,
                    "Number of NH chain iterations (nkija)",
                    default=1,
                    module="AMBER Wizard",
                    description="Optimized isokinetic Nose-Hoover: NH chain iterations (nkija)",
                )
                self.config["nkija"] = nkija

                self.console.print("\nInitial velocity distribution:")
                self.console.print("[cyan]0.[/cyan] Uniform")
                self.console.print("[cyan]1.[/cyan] Gaussian")
                idistr = prompt_with_context(
                    self.processor,
                    "Distribution type (idistr)",
                    choices=["0", "1"],
                    default="0",
                    module="AMBER Wizard",
                    description="Initial velocity distribution (idistr) for ntt=9",
                    options_map={"0": "Uniform", "1": "Gaussian"},
                )
                self.config["idistr"] = int(idistr)

            elif ntt == "10":
                # Stochastic isokinetic Nose-Hoover
                self._show_param_hint("nkija")
                nkija = int_prompt_with_context(
                    self.processor,
                    "Number of NH chain iterations (nkija)",
                    default=1,
                    module="AMBER Wizard",
                    description="Stochastic isokinetic Nose-Hoover: NH chain iterations (nkija)",
                )
                self.config["nkija"] = nkija

                self._show_param_hint("sinrtau")
                sinrtau = float_prompt_with_context(
                    self.processor,
                    "Stochastic time constant (sinrtau) [ps]",
                    default=1.0,
                    module="AMBER Wizard",
                    description="Stochastic isokinetic Nose-Hoover time constant (sinrtau, ps)",
                )
                self.config["sinrtau"] = sinrtau

                self.console.print("\nInitial velocity distribution:")
                self.console.print("[cyan]0.[/cyan] Uniform")
                self.console.print("[cyan]1.[/cyan] Gaussian")
                idistr = prompt_with_context(
                    self.processor,
                    "Distribution type (idistr)",
                    choices=["0", "1"],
                    default="0",
                    module="AMBER Wizard",
                    description="Initial velocity distribution (idistr) for ntt=10",
                    options_map={"0": "Uniform", "1": "Gaussian"},
                )
                self.config["idistr"] = int(idistr)

            elif ntt == "11":
                self._show_param_hint("tautp")
                tautp = float_prompt_with_context(
                    self.processor,
                    "Thermostat time constant (tautp) [ps]",
                    default=1.0,
                    module="AMBER Wizard",
                    description="Bussi velocity rescaling coupling time (tautp, ps)",
                )
                self.config["tautp"] = tautp

            # Middle scheme integrator (alternative to standard Verlet)
            self.console.print("\n[bold cyan]8e.[/bold cyan] Integration scheme:")
            self._show_param_hint("ischeme")
            use_middle = confirm_with_context(
                self.processor,
                "Use 'middle' scheme integrator?",
                default=False,
                module="AMBER Wizard",
                description="Step 8e: enable 'middle' scheme integrator (ischeme)",
            )

            if use_middle:
                self.config["ischeme"] = 1  # Middle scheme

                self._show_param_hint("ithermostat", recommended=2)
                self.console.print("\nMiddle scheme thermostat:")
                self.console.print("[cyan]0.[/cyan] No thermostat in middle scheme")
                self.console.print("[cyan]1.[/cyan] Langevin (uses gamma_ln)")
                self.console.print("[cyan]2.[/cyan] Velocity rescaling (Bussi-like)")

                ithermostat = prompt_with_context(
                    self.processor,
                    "Middle scheme thermostat",
                    choices=["0", "1", "2"],
                    default="2",
                    module="AMBER Wizard",
                    description="Middle scheme thermostat (ithermostat)",
                    options_map={
                        "0": "No thermostat",
                        "1": "Langevin (uses gamma_ln)",
                        "2": "Velocity rescaling (Bussi-like)",
                    },
                )
                self.config["ithermostat"] = int(ithermostat)

                if ithermostat == "2":
                    self._show_param_hint("therm_par")
                    therm_par = float_prompt_with_context(
                        self.processor,
                        "Thermostat coupling parameter (therm_par)",
                        default=3.0,
                        module="AMBER Wizard",
                        description="Middle scheme thermostat coupling parameter (therm_par)",
                    )
                    self.config["therm_par"] = therm_par

                self.console.print("[cyan]ℹ️  Middle scheme can improve sampling and stability[/cyan]")

            # Velocity limiting (stability)
            self.console.print("\n[bold cyan]8f.[/bold cyan] Velocity control:")
            self._show_param_hint("vlimit")
            use_vlimit = confirm_with_context(
                self.processor,
                "Enable velocity limiting for stability?",
                default=True,
                module="AMBER Wizard",
                description="Step 8f: enable velocity limiting (vlimit)?",
            )
            if use_vlimit:
                vlimit = float_prompt_with_context(
                    self.processor,
                    "Velocity limit (vlimit) [Å/ps]",
                    default=20.0,
                    module="AMBER Wizard",
                    description="Velocity limit (vlimit, Å/ps)",
                )
                self.config["vlimit"] = vlimit
                
            # LES particles (if applicable)
            use_les = confirm_with_context(
                self.processor,
                "Using LES particles?",
                default=False,
                module="AMBER Wizard",
                description="Step 8f: is this a Locally Enhanced Sampling (LES) simulation?",
            )
            if use_les:
                temp0les = float_prompt_with_context(
                    self.processor,
                    "LES particle temperature (temp0les) [K]",
                    default=-1,
                    module="AMBER Wizard",
                    description="LES particle target temperature (temp0les, K; <0 = same as temp0)",
                )
                if temp0les >= 0:
                    self.config["temp0les"] = temp0les
                
            # Random seed
            self.console.print("\n[bold cyan]8g.[/bold cyan] Random number generation:")
            self._show_param_hint("ig")
            ig = int_prompt_with_context(
                self.processor,
                "Random seed (ig, -1 for time-based)",
                default=-1,
                module="AMBER Wizard",
                description="Step 8g: random number seed (ig)",
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
                module="AMBER Wizard",
                description="Step 9a: switch to constant pressure (will set ntb=2)?",
            )
            if use_pressure:
                self.config["ntb"] = 2
            else:
                self.console.print("[cyan]ℹ️  Constant volume simulation[/cyan]")
                return True
                
        if self.config.get("ntb", 1) == 2 or self.simulation_type in ["md", "custom"]:
            self.console.print("\n[bold cyan]9b.[/bold cyan] Pressure scaling method:")
            self._show_param_hint("ntp", recommended=1)
            self.console.print("[cyan]0.[/cyan] No pressure control")
            self.console.print("[cyan]1.[/cyan] Isotropic scaling")
            self.console.print("[cyan]2.[/cyan] Anisotropic scaling")
            self.console.print("[cyan]3.[/cyan] Semi-isotropic (x,y coupled, z independent)")
            self.console.print("[cyan]4.[/cyan] Semi-isotropic membrane (for lipid bilayers)")

            ntp = prompt_with_context(
                self.processor,
                "Select pressure control",
                choices=["0", "1", "2", "3", "4"],
                default="1",
                module="AMBER Wizard",
                description="Step 9b: pressure scaling method (ntp)",
                options_map={
                    "0": "No pressure control",
                    "1": "Isotropic scaling",
                    "2": "Anisotropic scaling",
                    "3": "Semi-isotropic (x,y coupled, z independent)",
                    "4": "Semi-isotropic membrane (lipid bilayers)",
                },
            )
            self.config["ntp"] = int(ntp)

            # Membrane-specific options for ntp=4
            if ntp == "4":
                self.console.print("\n[bold cyan]9c.[/bold cyan] Membrane orientation:")
                self._show_param_hint("baroscalingdir", recommended=3)
                self.console.print("[cyan]1.[/cyan] Membrane normal along X")
                self.console.print("[cyan]2.[/cyan] Membrane normal along Y")
                self.console.print("[cyan]3.[/cyan] Membrane normal along Z (standard)")

                baroscalingdir = prompt_with_context(
                    self.processor,
                    "Membrane normal direction",
                    choices=["1", "2", "3"],
                    default="3",
                    module="AMBER Wizard",
                    description="Step 9c: membrane normal direction (baroscalingdir) for ntp=4",
                    options_map={
                        "1": "Membrane normal along X",
                        "2": "Membrane normal along Y",
                        "3": "Membrane normal along Z (standard)",
                    },
                )
                self.config["baroscalingdir"] = int(baroscalingdir)
                self.console.print("[cyan]ℹ️  Semi-isotropic: membrane plane scales together, normal independently[/cyan]")
            
            if ntp != "0":
                # Reference pressure
                self.console.print("\n[bold cyan]9d.[/bold cyan] Pressure parameters:")
                self._show_param_hint("pres0")
                pres0 = float_prompt_with_context(
                    self.processor,
                    "Reference pressure (pres0) [bar]",
                    default=1.0,
                    module="AMBER Wizard",
                    description="Step 9d: target pressure (pres0, bar)",
                )
                self.config["pres0"] = pres0
                
                # Barostat type
                self.console.print("\n[bold cyan]9e.[/bold cyan] Barostat selection:")
                self._show_param_hint("barostat", recommended=2)
                self.console.print("[cyan]1.[/cyan] Berendsen")
                self.console.print("[cyan]2.[/cyan] Monte Carlo (recommended)")
                
                barostat = prompt_with_context(
                    self.processor,
                    "Select barostat",
                    choices=["1", "2"],
                    default="2",
                    module="AMBER Wizard",
                    description="Step 9e: barostat algorithm (barostat)",
                    options_map={
                        "1": "Berendsen",
                        "2": "Monte Carlo (recommended)",
                    },
                )
                self.config["barostat"] = int(barostat)
                
                if barostat == "2":
                    self.console.print("\n[bold cyan]9f.[/bold cyan] Monte Carlo barostat settings:")
                    self._show_param_hint("mcbarint")
                    mcbarint = int_prompt_with_context(
                        self.processor,
                        "MC barostat interval",
                        default=100,
                        module="AMBER Wizard",
                        description="Step 9f: Monte Carlo barostat attempt interval (mcbarint)",
                    )
                    self.config["mcbarint"] = mcbarint

                # Pressure relaxation time
                self.console.print("\n[bold cyan]9g.[/bold cyan] Coupling parameters:")
                self._show_param_hint("taup")
                taup = float_prompt_with_context(
                    self.processor,
                    "Pressure relaxation time (taup) [ps]",
                    default=2.0,
                    module="AMBER Wizard",
                    description="Step 9g: pressure coupling relaxation time (taup, ps)",
                )
                self.config["taup"] = taup

                # Compressibility
                self.console.print("\n[bold cyan]9h.[/bold cyan] Compressibility:")
                self._show_param_hint("comp")
                comp = float_prompt_with_context(
                    self.processor,
                    "Compressibility [10⁻⁶ bar⁻¹]",
                    default=44.6,
                    module="AMBER Wizard",
                    description="Step 9h: isothermal compressibility (comp, 10⁻⁶ bar⁻¹)",
                )
                self.config["comp"] = comp
                    
        return True
        
    def _configure_constraints(self) -> bool:
        """Step 6: Bond and angle constraints."""
        self.console.print("\n[bold]Step 6: Constraints[/bold]")
        
        if self.simulation_type == "minimization":
            self.console.print("[cyan]ℹ️  No constraints needed for minimization[/cyan]")
            return True
            
        self.console.print("\n[bold cyan]6a.[/bold cyan] Bond constraints:")
        self._show_param_hint("ntc", recommended=2)
        self.console.print("[cyan]1.[/cyan] No constraints (all bonds flexible)")
        self.console.print("[cyan]2.[/cyan] SHAKE (bonds with hydrogen)")
        self.console.print("[cyan]3.[/cyan] All bonds constrained")
        
        ntc = prompt_with_context(
            self.processor,
            "Select constraints",
            choices=["1", "2", "3"],
            default="2",
            module="AMBER Wizard",
            description="Step 6a: bond constraints (ntc)",
            options_map={
                "1": "No constraints (all bonds flexible)",
                "2": "SHAKE (bonds with hydrogen)",
                "3": "All bonds constrained",
            },
        )
        self.config["ntc"] = int(ntc)
        
        if ntc != "1":
            # Force evaluation
            self.console.print("\n[bold cyan]6b.[/bold cyan] Force evaluation:")
            self._show_param_hint("ntf")
            self.console.print("[cyan]1.[/cyan] Complete force evaluation")
            self.console.print("[cyan]2.[/cyan] Omit bonds with hydrogen (if ntc=2)")
            self.console.print("[cyan]3.[/cyan] Omit all bond forces (if ntc=3)")
            
            if ntc == 2:
                ntf = prompt_with_context(
                    self.processor,
                    "Force evaluation",
                    choices=["1", "2"],
                    default="2",
                    module="AMBER Wizard",
                    description="Step 6b: force evaluation flag (ntf) with SHAKE on H bonds",
                    options_map={
                        "1": "Complete force evaluation",
                        "2": "Omit bonds with hydrogen",
                    },
                )
            else:  # ntc == 3
                ntf = prompt_with_context(
                    self.processor,
                    "Force evaluation",
                    choices=["1", "3"],
                    default="3",
                    module="AMBER Wizard",
                    description="Step 6b: force evaluation flag (ntf) with all bonds constrained",
                    options_map={
                        "1": "Complete force evaluation",
                        "3": "Omit all bond forces",
                    },
                )
                
            self.config["ntf"] = int(ntf)
            
            # SHAKE tolerance
            if ntc >= 2:
                self._show_param_hint("tol")
                tol = float_prompt_with_context(
                    self.processor,
                    "SHAKE tolerance",
                    default=0.00001,
                    module="AMBER Wizard",
                    description="SHAKE geometric tolerance (tol)",
                )
                self.config["tol"] = tol
                
        # Water constraints
        if self.config.get("ntc", 1) >= 2:
            self.console.print("\n[bold cyan]6c.[/bold cyan] Water geometry:")
            self._show_param_hint("jfastw")
            settle = confirm_with_context(
                self.processor,
                "Use SETTLE for water molecules?",
                default=True,
                module="AMBER Wizard",
                description="Step 6c: use SETTLE for water geometry (jfastw=4)?",
            )
            if settle:
                self.config["jfastw"] = 4  # SETTLE algorithm
                
        return True

    def _configure_nonbonded(self) -> bool:
        """Step 10: Non-bonded interactions."""
        self.console.print("\n[bold]Step 10: Non-bonded Interactions[/bold]")

        self.console.print("\n[bold cyan]10a.[/bold cyan] Cutoff distance:")
        self._show_param_hint("cut", recommended=10.0)
        cut = float_prompt_with_context(
            self.processor,
            "Non-bonded cutoff [Å]",
            default=10.0,
            module="AMBER Wizard",
            description="Step 10a: non-bonded cutoff distance (cut, Å)",
        )
        self.config["cut"] = cut

        if cut > 12.0:
            self.console.print("[yellow]⚠️  Large cutoff - ensure box size is adequate[/yellow]")
        elif cut < 8.0:
            self.console.print("[yellow]⚠️  Small cutoff - may affect electrostatics[/yellow]")

        # Force switching for vdW
        self.console.print("\n[bold cyan]10b.[/bold cyan] Force switching:")
        self._show_param_hint("fswitch")
        use_fswitch = confirm_with_context(
            self.processor,
            "Use force switching for vdW interactions?",
            default=False,
            module="AMBER Wizard",
            description="Step 10b: enable vdW force switching (fswitch)?",
        )
        if use_fswitch:
            fswitch = float_prompt_with_context(
                self.processor,
                "Force switch distance (fswitch) [Å]",
                default=cut - 2.0,
                module="AMBER Wizard",
                description="vdW force switch distance (fswitch, Å)",
            )
            self.config["fswitch"] = fswitch
            self.console.print("[cyan]ℹ️  vdW forces smoothly switched between fswitch and cut[/cyan]")

        # Dielectric constant
        self.console.print("\n[bold cyan]10c.[/bold cyan] Dielectric constant:")
        self._show_param_hint("dielc", recommended=1.0)
        use_dielc = confirm_with_context(
            self.processor,
            "Modify dielectric constant (default=1.0)?",
            default=False,
            module="AMBER Wizard",
            description="Step 10c: change dielectric constant (dielc) from 1.0?",
        )
        if use_dielc:
            dielc = float_prompt_with_context(
                self.processor,
                "Dielectric constant (dielc)",
                default=1.0,
                module="AMBER Wizard",
                description="Dielectric constant (dielc)",
            )
            self.config["dielc"] = dielc

        # Neighbor list update frequency
        self.console.print("\n[bold cyan]10d.[/bold cyan] Neighbor list:")
        self._show_param_hint("nsnb", recommended=0)
        nsnb = int_prompt_with_context(
            self.processor,
            "Neighbor list update frequency (nsnb, 0=auto)",
            default=0,
            module="AMBER Wizard",
            description="Step 10d: neighbor list update frequency (nsnb)",
        )
        if nsnb > 0:
            self.config["nsnb"] = nsnb

        # Long-range electrostatics
        if self.config.get("ntb", 1) > 0:  # Periodic system
            self.console.print("\n[bold cyan]10e.[/bold cyan] Long-range electrostatics:")

            use_pme = confirm_with_context(
                self.processor,
                "Use Particle Mesh Ewald (PME)?",
                default=True,
                module="AMBER Wizard",
                description="Step 10e: use Particle Mesh Ewald for long-range electrostatics?",
            )

            if use_pme:
                self.console.print("[green]✓ PME will be used for long-range electrostatics[/green]")

                # PME parameters (optional advanced settings)
                if confirm_with_context(
                    self.processor,
                    "Configure PME parameters manually?",
                    default=False,
                    module="AMBER Wizard",
                    description="Step 10e: customize PME grid/Ewald settings manually?",
                ):
                    self.console.print("\n[bold cyan]10f.[/bold cyan] PME grid settings:")

                    # Ewald type
                    self._show_param_hint("ew_type", recommended=1)
                    self.console.print("\nEwald method:")
                    self.console.print("[cyan]0.[/cyan] Regular Ewald (slow)")
                    self.console.print("[cyan]1.[/cyan] PME (default, fast)")
                    ew_type = prompt_with_context(
                        self.processor,
                        "Ewald type (ew_type)",
                        choices=["0", "1"],
                        default="1",
                        module="AMBER Wizard",
                        description="Step 10f: Ewald method (ew_type)",
                        options_map={
                            "0": "Regular Ewald (slow)",
                            "1": "PME (fast, default)",
                        },
                    )
                    if ew_type == "0":
                        self.config["ew_type"] = 0

                    # B-spline order
                    self._show_param_hint("order", recommended=4)
                    order = int_prompt_with_context(
                        self.processor,
                        "B-spline interpolation order (order, 4-6)",
                        default=4,
                        module="AMBER Wizard",
                        description="PME B-spline interpolation order (order, 4–6)",
                    )
                    if order != 4:
                        self.config["order"] = order

                    # Grid points
                    self._show_param_hint("nfft1", recommended=0)
                    nfft1 = int_prompt_with_context(
                        self.processor,
                        "Grid points X-direction (nfft1, 0=auto)",
                        default=0,
                        module="AMBER Wizard",
                        description="PME grid points in X-direction (nfft1, 0=auto)",
                    )
                    if nfft1 > 0:
                        self.config["nfft1"] = nfft1
                        nfft2 = int_prompt_with_context(
                            self.processor,
                            "Grid points Y-direction (nfft2)",
                            default=nfft1,
                            module="AMBER Wizard",
                            description="PME grid points in Y-direction (nfft2)",
                        )
                        nfft3 = int_prompt_with_context(
                            self.processor,
                            "Grid points Z-direction (nfft3)",
                            default=nfft1,
                            module="AMBER Wizard",
                            description="PME grid points in Z-direction (nfft3)",
                        )
                        self.config["nfft2"] = nfft2
                        self.config["nfft3"] = nfft3

                    self._show_param_hint("dsum_tol", recommended=0.000001)
                    dsum_tol = float_prompt_with_context(
                        self.processor,
                        "PME direct sum tolerance",
                        default=0.000001,
                        module="AMBER Wizard",
                        description="PME direct sum tolerance (dsum_tol)",
                    )
                    self.config["dsum_tol"] = dsum_tol
            else:
                self.console.print("[yellow]⚠️  No PME - results may be poor for charged systems[/yellow]")

        # Electric field application
        self.console.print("\n[bold cyan]10g.[/bold cyan] External electric field:")
        self._show_param_hint("efx", recommended=0.0)
        use_efield = confirm_with_context(
            self.processor,
            "Apply external electric field?",
            default=False,
            module="AMBER Wizard",
            description="Step 10g: apply external electric field?",
        )
        if use_efield:
            self.console.print("[cyan]ℹ️  Electric field units: kcal/(mol·Å·e)[/cyan]")
            efx = float_prompt_with_context(
                self.processor,
                "Field in X direction (efx)",
                default=0.0,
                module="AMBER Wizard",
                description="External E-field X component (efx, kcal/mol·Å·e)",
            )
            self._show_param_hint("efy", recommended=0.0)
            efy = float_prompt_with_context(
                self.processor,
                "Field in Y direction (efy)",
                default=0.0,
                module="AMBER Wizard",
                description="External E-field Y component (efy, kcal/mol·Å·e)",
            )
            self._show_param_hint("efz", recommended=0.0)
            efz = float_prompt_with_context(
                self.processor,
                "Field in Z direction (efz)",
                default=0.0,
                module="AMBER Wizard",
                description="External E-field Z component (efz, kcal/mol·Å·e)",
            )
            if efx != 0.0:
                self.config["efx"] = efx
            if efy != 0.0:
                self.config["efy"] = efy
            if efz != 0.0:
                self.config["efz"] = efz
            self.console.print("[yellow]⚠️  Electric field adds constant force on charged atoms[/yellow]")

        return True
        
    def _configure_output(self) -> bool:
        """Step 11: Output control."""
        self.console.print("\n[bold]Step 11: Output Control[/bold]")
        
        self.console.print("\n[bold cyan]11a.[/bold cyan] Output frequencies:")

        # Energy output
        if self.simulation_type == "minimization":
            ntpr_default = 100
        else:
            ntpr_default = 1000

        self._show_param_hint("ntpr", recommended=ntpr_default)
        ntpr = int_prompt_with_context(
            self.processor,
            "Energy output frequency (ntpr)",
            default=ntpr_default,
            module="AMBER Wizard",
            description="Step 11a: energy/status output frequency (ntpr)",
        )
        self.config["ntpr"] = ntpr
        
        # Coordinate output (MD only)
        if self.simulation_type != "minimization":
            self._show_param_hint("ntwx", recommended=ntpr)
            ntwx = int_prompt_with_context(
                self.processor,
                "Coordinate output frequency (ntwx)",
                default=ntpr,
                module="AMBER Wizard",
                description="Step 11a: trajectory coordinate output frequency (ntwx)",
            )
            self.config["ntwx"] = ntwx
            
            # Additional output options
            self.console.print("\n[bold cyan]11b.[/bold cyan] Additional trajectory outputs:")
            self._show_param_hint("ntwv", recommended=0)
            write_velocities = confirm_with_context(
                self.processor,
                "Write velocity trajectory?",
                default=False,
                module="AMBER Wizard",
                description="Step 11b: write velocity trajectory (ntwv)?",
            )
            if write_velocities:
                ntwv_choice = prompt_with_context(
                    self.processor,
                    "Velocity output: 1=separate file, -1=combined with coordinates",
                    choices=["1", "-1"],
                    default="1",
                    module="AMBER Wizard",
                    description="Step 11b: velocity output destination (ntwv)",
                    options_map={
                        "1": "Separate velocity file",
                        "-1": "Combined with coordinate trajectory",
                    },
                )
                self.config["ntwv"] = int(ntwv_choice)
            
            # Force trajectory
            self._show_param_hint("ntwf", recommended=0)
            write_forces = confirm_with_context(
                self.processor,
                "Write force trajectory?",
                default=False,
                module="AMBER Wizard",
                description="Step 11b: write force trajectory (ntwf)?",
            )
            if write_forces:
                ntwf_choice = prompt_with_context(
                    self.processor,
                    "Force output: 1=separate file, -1=combined with coordinates",
                    choices=["1", "-1"],
                    default="1",
                    module="AMBER Wizard",
                    description="Step 11b: force output destination (ntwf)",
                    options_map={
                        "1": "Separate force file",
                        "-1": "Combined with coordinate trajectory",
                    },
                )
                self.config["ntwf"] = int(ntwf_choice)
            
            # Energy file
            self._show_param_hint("ntwe", recommended=0)
            write_energies = confirm_with_context(
                self.processor,
                "Write compact energy file (mden)?",
                default=False,
                module="AMBER Wizard",
                description="Step 11b: write compact mden energy file (ntwe)?",
            )
            if write_energies:
                ntwe = int_prompt_with_context(
                    self.processor,
                    "Energy write frequency (ntwe)",
                    default=1000,
                    module="AMBER Wizard",
                    description="mden compact energy output frequency (ntwe)",
                )
                self.config["ntwe"] = ntwe
            
        # Restart file frequency
        self.console.print("\n[bold cyan]11c.[/bold cyan] Restart file frequency:")
        self._show_param_hint("ntwr", recommended=ntpr * 10)
        ntwr = int_prompt_with_context(
            self.processor,
            "Restart file frequency (ntwr)",
            default=ntpr * 10,
            module="AMBER Wizard",
            description="Step 11c: restart file write frequency (ntwr)",
        )
        self.config["ntwr"] = ntwr

        # Additional output control parameters
        self.console.print("\n[bold cyan]11d.[/bold cyan] Additional output parameters:")
        
        # Running averages
        if self.simulation_type != "minimization":
            self._show_param_hint("ntave", recommended=0)
            use_averages = confirm_with_context(
                self.processor,
                "Output running averages?",
                default=False,
                module="AMBER Wizard",
                description="Step 11d: print running averages (ntave)?",
            )
            if use_averages:
                ntave = int_prompt_with_context(
                    self.processor,
                    "Running average frequency (ntave)",
                    default=ntpr,
                    module="AMBER Wizard",
                    description="Running average window size (ntave)",
                )
                self.config["ntave"] = ntave
        
        # Coordinate wrapping
        if self.config.get("ntb", 0) > 0:  # Only for periodic systems
            self._show_param_hint("iwrap", recommended=0)
            iwrap = confirm_with_context(
                self.processor,
                "Wrap coordinates to primary box?",
                default=False,
                module="AMBER Wizard",
                description="Step 11d: wrap coordinates into primary unit cell (iwrap)?",
            )
            if iwrap:
                self.config["iwrap"] = 1
        
        # Velocity output type
        if "ntwv" in self.config and self.config["ntwv"] > 0:
            self._show_param_hint("ionstepvelocities", recommended=0)
            ionstepvelocities = confirm_with_context(
                self.processor,
                "Use on-step velocities (vs half-step)?",
                default=False,
                module="AMBER Wizard",
                description="Step 11d: write on-step (vs leapfrog half-step) velocities (ionstepvelocities)",
            )
            if ionstepvelocities:
                self.config["ionstepvelocities"] = 1
        
        # Trajectory subset (for large systems)
        if self.simulation_type != "minimization":
            self._show_param_hint("ntwprt", recommended=0)
            use_subset = confirm_with_context(
                self.processor,
                "Output only subset of atoms?",
                default=False,
                module="AMBER Wizard",
                description="Step 11d: limit trajectory to first N atoms (ntwprt)?",
            )
            if use_subset:
                natoms_prompt = "Number of atoms to output (ntwprt, 0=all)"
                ntwprt = int_prompt_with_context(
                    self.processor,
                    natoms_prompt,
                    default=0,
                    module="AMBER Wizard",
                    description="Number of atoms to write to trajectory (ntwprt)",
                )
                if ntwprt > 0:
                    self.config["ntwprt"] = ntwprt
            
        self.console.print("\n[bold cyan]11e.[/bold cyan] File formats:")

        # Coordinate format
        self._show_param_hint("ioutfm", recommended=1)
        self.console.print("[cyan]1.[/cyan] ASCII (human readable)")
        self.console.print("[cyan]2.[/cyan] NetCDF (binary, smaller)")

        ioutfm = prompt_with_context(
            self.processor,
            "Coordinate format",
            choices=["1", "2"],
            default="2",
            module="AMBER Wizard",
            description="Step 11e: trajectory/restart file format (ioutfm, ntxo)",
            options_map={
                "1": "ASCII (human readable)",
                "2": "NetCDF (binary, recommended)",
            },
        )
        if ioutfm == "2":
            self.config["ioutfm"] = 1  # NetCDF for trajectory
            self.config["ntxo"] = 2    # NetCDF for restart
            
        return True
        
    def _configure_advanced_options(self) -> bool:
        """Step 12: Advanced Options."""
        self.console.print("\n[bold]Step 12: Advanced Options[/bold]")

        # NMR restraints and varying conditions
        self.console.print("\n[bold cyan]12a.[/bold cyan] Varying conditions:")
        if self.simulation_type in ["md", "custom", "heating"]:
            self._show_param_hint("nmropt", recommended=0)
            use_nmr = confirm_with_context(
                self.processor,
                "Use varying conditions (&wt blocks)?",
                default=False,
                module="AMBER Wizard",
                description="Step 12a: enable &wt varying-conditions blocks (nmropt)?",
            )

            if use_nmr:
                self.config["nmropt"] = 1
                self.console.print("[cyan]ℹ️  nmropt=1 enables &wt varying conditions[/cyan]")

                # Configure &wt blocks
                self._configure_wt_blocks()

                # File redirection commands (available with nmropt=1)
                self.console.print("\n[bold cyan]12b.[/bold cyan] File redirection (nmropt output):")
                self.console.print("[grey50]These commands control output for varying conditions and restraints[/grey50]")

                use_listout = confirm_with_context(
                    self.processor,
                    "Write restraint/condition listing (LISTOUT)?",
                    default=False,
                    module="AMBER Wizard",
                    description="Step 12b: enable LISTOUT restraint/condition listing?",
                )
                if use_listout:
                    listout_file = prompt_with_context(
                        self.processor,
                        "LISTOUT filename",
                        default="restraint_list.out",
                        module="AMBER Wizard",
                        description="LISTOUT output filename",
                    )
                    self.config["_listout"] = listout_file

                use_dumpave = confirm_with_context(
                    self.processor,
                    "Write time-averaged values (DUMPAVE)?",
                    default=False,
                    module="AMBER Wizard",
                    description="Step 12b: enable DUMPAVE time-averaged values?",
                )
                if use_dumpave:
                    dumpave_file = prompt_with_context(
                        self.processor,
                        "DUMPAVE filename",
                        default="restraint_ave.out",
                        module="AMBER Wizard",
                        description="DUMPAVE output filename",
                    )
                    self.config["_dumpave"] = dumpave_file

                use_listin = confirm_with_context(
                    self.processor,
                    "Read restraint info from file (LISTIN)?",
                    default=False,
                    module="AMBER Wizard",
                    description="Step 12b: read restraint info via LISTIN?",
                )
                if use_listin:
                    self.console.print("[cyan]1.[/cyan] POUT - Echo restraints to mdout")
                    self.console.print("[cyan]2.[/cyan] Specify input file")
                    listin_choice = prompt_with_context(
                        self.processor,
                        "LISTIN option",
                        choices=["1", "2"],
                        default="1",
                        module="AMBER Wizard",
                        description="LISTIN source selection",
                        options_map={
                            "1": "POUT — echo restraints to mdout",
                            "2": "Specify input file",
                        },
                    )
                    if listin_choice == "1":
                        self.config["_listin"] = "POUT"
                    else:
                        listin_file = prompt_with_context(
                            self.processor,
                            "LISTIN filename",
                            default="restraint_list.in",
                            module="AMBER Wizard",
                            description="LISTIN input filename",
                        )
                        self.config["_listin"] = listin_file
        else:
            self.console.print("[grey50]Varying conditions not applicable for this simulation type[/grey50]")

        # PMEMD-specific tuning options
        self.console.print("\n[bold cyan]12c.[/bold cyan] PMEMD-specific options:")
        use_pmemd_opts = confirm_with_context(
            self.processor,
            "Configure PMEMD-specific parameters?",
            default=False,
            module="AMBER Wizard",
            description="Step 12c: configure PMEMD-specific tuning parameters?",
        )

        if use_pmemd_opts:
            # Output flush interval
            mdout_flush = int_prompt_with_context(
                self.processor,
                "Output flush interval (mdout_flush_interval)",
                default=100,
                module="AMBER Wizard",
                description="PMEMD mdout flush interval (mdout_flush_interval)",
            )
            if mdout_flush != 100:
                self.config["mdout_flush_interval"] = mdout_flush

            # Separate cutoffs for ES and vdW
            self.console.print("\nSeparate cutoffs (PMEMD only):")
            use_sep_cutoffs = confirm_with_context(
                self.processor,
                "Use separate ES/vdW cutoffs?",
                default=False,
                module="AMBER Wizard",
                description="PMEMD: use distinct electrostatic and vdW cutoffs?",
            )
            if use_sep_cutoffs:
                es_cutoff = float_prompt_with_context(
                    self.processor,
                    "Electrostatic cutoff (es_cutoff) [Å]",
                    default=self.config.get("cut", 10.0),
                    module="AMBER Wizard",
                    description="PMEMD electrostatic cutoff (es_cutoff, Å)",
                )
                vdw_cutoff = float_prompt_with_context(
                    self.processor,
                    "van der Waals cutoff (vdw_cutoff) [Å]",
                    default=self.config.get("cut", 10.0),
                    module="AMBER Wizard",
                    description="PMEMD van der Waals cutoff (vdw_cutoff, Å)",
                )
                self.config["es_cutoff"] = es_cutoff
                self.config["vdw_cutoff"] = vdw_cutoff

            # FFT grid spacing
            fft_grids = float_prompt_with_context(
                self.processor,
                "FFT grids per angstrom (fft_grids_per_ang)",
                default=1.0,
                module="AMBER Wizard",
                description="PMEMD FFT grid density (fft_grids_per_ang)",
            )
            if fft_grids != 1.0:
                self.config["fft_grids_per_ang"] = fft_grids

            self.console.print("[cyan]ℹ️  These options only apply when running with pmemd[/cyan]")
                
        # Energy decomposition
        self.console.print("\n[bold cyan]12d.[/bold cyan] Energy decomposition:")
        self._show_param_hint("idecomp", recommended=0)
        use_decomp = confirm_with_context(
            self.processor,
            "Enable energy decomposition?",
            default=False,
            module="AMBER Wizard",
            description="Step 12d: enable per-residue energy decomposition (idecomp)?",
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
                module="AMBER Wizard",
                description="Step 12d: energy decomposition type (idecomp)",
                options_map={
                    "1": "Per-residue (1-4 with internal)",
                    "2": "Per-residue (1-4 with non-bonded)",
                    "3": "Pairwise per-residue (1-4 with internal)",
                    "4": "Pairwise per-residue (1-4 with non-bonded)",
                },
            )
            self.config["idecomp"] = int(idecomp)
            
        # Belly dynamics (legacy)
        self.console.print("\n[bold cyan]12e.[/bold cyan] Belly dynamics:")
        self._show_param_hint("ibelly", recommended=0)
        if confirm_with_context(
            self.processor,
            "Use belly dynamics (subset of moving atoms)?",
            default=False,
            module="AMBER Wizard",
            description="Step 12e: enable legacy belly dynamics (ibelly)?",
        ):
            self.config["ibelly"] = 1
            bellymask = prompt_with_context(
                self.processor,
                "Belly mask (moving atoms)",
                default="!:WAT",
                module="AMBER Wizard",
                description="AMBER mask selecting atoms allowed to move (bellymask)",
            )
            self.config["bellymask"] = bellymask
            self.console.print("[yellow]⚠️  Belly dynamics is legacy - consider restraints instead[/yellow]")
        
        # SHAKE customization
        if self.config.get("ntc", 1) > 1:  # If SHAKE is enabled
            self.console.print("\n[bold cyan]12f.[/bold cyan] SHAKE customization:")
            
            # Custom water names
            custom_water = confirm_with_context(
                self.processor,
                "Customize water residue/atom names?",
                default=False,
                module="AMBER Wizard",
                description="Step 12f: override default SHAKE water residue/atom names?",
            )
            if custom_water:
                watnam = prompt_with_context(
                    self.processor,
                    "Water residue name (WATNAM)",
                    default="WAT",
                    module="AMBER Wizard",
                    description="SHAKE water residue name (WATNAM)",
                )
                self.config["WATNAM"] = f"'{watnam}'"

                owtnm = prompt_with_context(
                    self.processor,
                    "Oxygen atom name (OWTNM)",
                    default="O",
                    module="AMBER Wizard",
                    description="SHAKE water oxygen atom name (OWTNM)",
                )
                self.config["OWTNM"] = f"'{owtnm}'"

                hwtnm1 = prompt_with_context(
                    self.processor,
                    "First hydrogen name (HWTNM1)",
                    default="H1",
                    module="AMBER Wizard",
                    description="SHAKE water first hydrogen name (HWTNM1)",
                )
                self.config["HWTNM1"] = f"'{hwtnm1}'"

                hwtnm2 = prompt_with_context(
                    self.processor,
                    "Second hydrogen name (HWTNM2)",
                    default="H2",
                    module="AMBER Wizard",
                    description="SHAKE water second hydrogen name (HWTNM2)",
                )
                self.config["HWTNM2"] = f"'{hwtnm2}'"
            
            # No-shake mask
            self._show_param_hint("noshakemask")
            use_noshake = confirm_with_context(
                self.processor,
                "Exclude atoms from SHAKE?",
                default=False,
                module="AMBER Wizard",
                description="Step 12f: exclude atoms from SHAKE (noshakemask)?",
            )
            if use_noshake:
                noshakemask = prompt_with_context(
                    self.processor,
                    "No-shake mask (atoms to exclude)",
                    default=":1-10",
                    module="AMBER Wizard",
                    description="AMBER mask of atoms to exclude from SHAKE (noshakemask)",
                )
                self.config["noshakemask"] = noshakemask
                self.console.print("[yellow]⚠️  noshakemask will set ntf=1 automatically[/yellow]")

        # Reference guide access
        self.console.print("\n[bold cyan]12g.[/bold cyan] Parameter Reference Guide:")
        if confirm_with_context(
            self.processor,
            "View full parameter reference guide?",
            default=False,
            module="AMBER Wizard",
            description="Step 12g: display the full parameter reference guide now?",
        ):
            print_reference_guide(self.param_db, self.console)

        return True
        
    def _show_configuration_summary(self) -> None:
        """Display final configuration summary (Step 13)."""
        self.console.print("\n[bold green]═══ Configuration Complete! ═══[/bold green]")
        self.console.print(f"[cyan]Simulation type: {self.simulation_type}[/cyan]\n")

        # Comprehensive parameter descriptions organized by category
        param_info = {
            # Simulation Control
            "imin": ("Simulation Control", "Minimization flag (0=MD, 1=minimize)"),
            "ntx": ("Simulation Control", "Coordinate reading format"),
            "irest": ("Simulation Control", "Restart flag (0=new, 1=restart)"),
            "nstlim": ("Simulation Control", "Number of MD steps"),
            "t": ("Simulation Control", "Initial time (ps)"),

            # Minimization
            "maxcyc": ("Minimization", "Maximum minimization cycles"),
            "ncyc": ("Minimization", "Steepest descent cycles before CG"),
            "ntmin": ("Minimization", "Minimization method (0-5)"),
            "drms": ("Minimization", "Convergence criterion (kcal/mol/Å)"),
            "dx0": ("Minimization", "Initial step size"),

            # Time Control
            "dt": ("Time Control", "Time step (ps)"),
            "nrespa": ("Time Control", "RESPA multiple timestep factor"),
            "nscm": ("Time Control", "COM motion removal frequency"),

            # Temperature Control
            "ntt": ("Temperature", "Thermostat type (0-3,9-11)"),
            "temp0": ("Temperature", "Target temperature (K)"),
            "tempi": ("Temperature", "Initial temperature (K)"),
            "temp0les": ("Temperature", "LES particle temperature (K)"),
            "gamma_ln": ("Temperature", "Langevin collision frequency (ps⁻¹)"),
            "tautp": ("Temperature", "Berendsen/Bussi coupling time (ps)"),
            "vrand": ("Temperature", "Velocity randomization frequency"),
            "vlimit": ("Temperature", "Velocity limit (default -1=off)"),
            "ig": ("Temperature", "Random seed (-1=based on time)"),
            "nkija": ("Temperature", "Nose-Hoover chain iterations (ntt=9,10)"),
            "idistr": ("Temperature", "Initial velocity distribution (0=uniform)"),
            "sinrtau": ("Temperature", "Stochastic time constant (ntt=10)"),

            # Middle Scheme Integrator
            "ischeme": ("Integrator", "Integration scheme (0=std, 1=middle)"),
            "ithermostat": ("Integrator", "Middle scheme thermostat (0-2)"),
            "therm_par": ("Integrator", "Thermostat coupling parameter"),

            # Pressure Control
            "ntp": ("Pressure", "Pressure control (0-4)"),
            "pres0": ("Pressure", "Target pressure (bar)"),
            "taup": ("Pressure", "Pressure coupling time (ps)"),
            "barostat": ("Pressure", "Barostat type (1=Berendsen, 2=MC)"),
            "mcbarint": ("Pressure", "MC barostat attempt interval"),
            "comp": ("Pressure", "Compressibility (bar⁻¹)"),
            "baroscalingdir": ("Pressure", "Membrane normal direction (1=X,2=Y,3=Z)"),

            # System Setup
            "ntb": ("System Setup", "Periodic boundaries (0=none, 1=const V, 2=const P)"),
            "iwrap": ("System Setup", "Wrap coordinates into box"),
            "ifbox": ("System Setup", "Box type"),

            # Implicit Solvent (GB)
            "igb": ("Implicit Solvent", "GB model (0=off, 1=HCT, 2=OBC, 5-8)"),
            "saltcon": ("Implicit Solvent", "Salt concentration (M)"),
            "rgbmax": ("Implicit Solvent", "GB interaction cutoff (Å)"),
            "gbsa": ("Implicit Solvent", "Surface area term (0=off, 1=LCPO)"),
            "surften": ("Implicit Solvent", "Surface tension (cal/mol·Å²)"),

            # Water Cap
            "ivcap": ("Water Cap", "Cap type (1=cap, 2=virtual box)"),
            "fcap": ("Water Cap", "Cap force constant (kcal/mol·Å²)"),
            "cutcap": ("Water Cap", "Cap radius (Å)"),
            "xcap": ("Water Cap", "Cap center X coordinate"),
            "ycap": ("Water Cap", "Cap center Y coordinate"),
            "zcap": ("Water Cap", "Cap center Z coordinate"),

            # Constraints
            "ntc": ("Constraints", "SHAKE constraint type"),
            "ntf": ("Constraints", "Force evaluation exclusions"),
            "tol": ("Constraints", "SHAKE tolerance"),
            "jfastw": ("Constraints", "Fast water (4=SETTLE)"),
            "noshakemask": ("Constraints", "Atoms excluded from SHAKE"),

            # Non-bonded
            "cut": ("Non-bonded", "Nonbonded cutoff (Å)"),
            "fswitch": ("Non-bonded", "Force switching distance (Å)"),
            "dielc": ("Non-bonded", "Dielectric constant"),
            "nsnb": ("Non-bonded", "Neighbor list update frequency"),

            # PME / Electrostatics
            "ew_type": ("Electrostatics", "Ewald type (0=regular, 1=PME)"),
            "order": ("Electrostatics", "B-spline interpolation order"),
            "nfft1": ("Electrostatics", "PME grid X dimension"),
            "nfft2": ("Electrostatics", "PME grid Y dimension"),
            "nfft3": ("Electrostatics", "PME grid Z dimension"),
            "dsum_tol": ("Electrostatics", "PME direct sum tolerance"),
            "eedmeth": ("Electrostatics", "Electrostatic method"),
            "scnb": ("Electrostatics", "1-4 vdW scaling"),
            "scee": ("Electrostatics", "1-4 electrostatic scaling"),

            # Electric Field
            "efx": ("Electric Field", "Field in X direction (kcal/mol·Å·e)"),
            "efy": ("Electric Field", "Field in Y direction (kcal/mol·Å·e)"),
            "efz": ("Electric Field", "Field in Z direction (kcal/mol·Å·e)"),

            # Output Control
            "ntpr": ("Output", "Energy output frequency"),
            "ntwx": ("Output", "Coordinate output frequency"),
            "ntwr": ("Output", "Restart file frequency"),
            "ntwe": ("Output", "Energy file frequency"),
            "ntwv": ("Output", "Velocity output frequency"),
            "ntwf": ("Output", "Force output frequency"),
            "ntave": ("Output", "Averaging frequency"),
            "ioutfm": ("Output", "Output format (0=ASCII, 1=NetCDF)"),
            "ntxo": ("Output", "Restart format (1=ASCII, 2=NetCDF)"),
            "ntwprt": ("Output", "Atoms to write (0=all)"),
            "ionstepvelocities": ("Output", "On-step velocities (1=yes)"),

            # PMEMD-specific
            "mdout_flush_interval": ("PMEMD", "Output flush interval"),
            "es_cutoff": ("PMEMD", "Electrostatic cutoff (Å)"),
            "vdw_cutoff": ("PMEMD", "van der Waals cutoff (Å)"),
            "fft_grids_per_ang": ("PMEMD", "FFT grids per angstrom"),

            # Advanced
            "nmropt": ("Advanced", "NMR/restraint options"),
            "pencut": ("Advanced", "Cutoff for 1-4 printing"),
            "ibelly": ("Advanced", "Belly dynamics flag"),
            "bellymask": ("Advanced", "Belly atom mask"),
            "idecomp": ("Advanced", "Energy decomposition"),

            # File Redirection (internal markers)
            "_listout": ("File Redirection", "LISTOUT restraint/condition listing"),
            "_dumpave": ("File Redirection", "DUMPAVE time-averaged values"),
            "_listin": ("File Redirection", "LISTIN restraint input (or POUT)"),
        }

        # Group parameters by category
        categories = {}
        for param, value in self.config.items():
            # Skip internal parameters (prefixed with _)
            if param.startswith("_"):
                continue
            if param in param_info:
                category, description = param_info[param]
            else:
                category = "Other"
                description = ""
            if category not in categories:
                categories[category] = []
            categories[category].append((param, value, description))

        # Define category order
        category_order = [
            "Simulation Control", "Minimization", "Time Control", "Temperature",
            "Integrator", "Pressure", "System Setup", "Implicit Solvent",
            "Water Cap", "Constraints", "Non-bonded", "Electrostatics",
            "Electric Field", "Output", "PMEMD", "Advanced",
            "File Redirection", "Other"
        ]

        # Display one table per category for clean grouping
        for category in category_order:
            if category not in categories:
                continue
            params = categories[category]
            if not params:
                continue

            table = Table(title=f"[bold yellow]{category}[/bold yellow]",
                          box=box.SIMPLE_HEAVY, expand=False, show_lines=False)
            table.add_column("Parameter", style="cyan")
            table.add_column("Value", style="green")
            table.add_column("Description", style="white")
            table.add_column("Commonality", width=18)

            for param, value, description in sorted(params):
                # Format value nicely
                if isinstance(value, float):
                    value_str = f"{value:g}"
                elif isinstance(value, list):
                    value_str = ", ".join(str(v) for v in value)
                else:
                    value_str = str(value)

                # Commonality tag from parameter database
                if param in self.param_db:
                    db_param = self.param_db[param]
                    if db_param.common:
                        tag = "[green]commonly changed[/green]"
                    else:
                        tag = "rarely changed"
                else:
                    tag = ""
                table.add_row(param, value_str, description, tag)

            self.console.print(table)
            self.console.print()

        # Show &wt blocks if configured
        if "_wt_blocks" in self.config:
            self.console.print("\n[bold yellow]&wt Varying Conditions:[/bold yellow]")
            for i, block in enumerate(self.config["_wt_blocks"], 1):
                self.console.print(f"  Block {i}: [grey50]{block}[/grey50]")
        
    def _configure_wt_blocks(self) -> bool:
        """Configure &wt varying conditions blocks."""
        self.console.print("\n[bold cyan]&wt Block Configuration:[/bold cyan]")
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
            module="AMBER Wizard",
            description="&wt configuration preset",
            options_map={
                "1": "Temperature ramping",
                "2": "Restraint ramping",
                "3": "Temperature + restraint ramping",
                "4": "Custom &wt blocks (advanced)",
            },
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
                    module="AMBER Wizard",
                    description=f"Temperature ramp {ramp_num}: start step (istep1)",
                )
            else:
                istep1 = int_prompt_with_context(
                    self.processor,
                    "Start step",
                    default=last_step,
                    module="AMBER Wizard",
                    description=f"Temperature ramp {ramp_num}: start step (istep1)",
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
                module="AMBER Wizard",
                description=f"Temperature ramp {ramp_num}: end step (istep2)",
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
                    module="AMBER Wizard",
                    description=f"Temperature ramp {ramp_num}: initial temperature (value1, K)",
                )
            else:
                temp1 = float_prompt_with_context(
                    self.processor,
                    "Initial temperature [K]",
                    default=300.0,
                    module="AMBER Wizard",
                    description=f"Temperature ramp {ramp_num}: initial temperature (value1, K)",
                )

            temp2 = float_prompt_with_context(
                self.processor,
                "Final temperature [K]",
                default=300.0,
                module="AMBER Wizard",
                description=f"Temperature ramp {ramp_num}: final temperature (value2, K)",
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
                module="AMBER Wizard",
                description="Add another temperature ramp segment?",
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
                module="AMBER Wizard",
                description=f"Restraint ramp {ramp_num}: start step (istep1)",
            )
            istep2 = int_prompt_with_context(
                self.processor,
                "End step (0 = end of simulation)",
                default=0,
                module="AMBER Wizard",
                description=f"Restraint ramp {ramp_num}: end step (istep2, 0=end)",
            )
            weight1 = float_prompt_with_context(
                self.processor,
                "Initial restraint weight",
                default=10.0,
                module="AMBER Wizard",
                description=f"Restraint ramp {ramp_num}: initial weight (value1)",
            )
            weight2 = float_prompt_with_context(
                self.processor,
                "Final restraint weight",
                default=1.0,
                module="AMBER Wizard",
                description=f"Restraint ramp {ramp_num}: final weight (value2)",
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
                module="AMBER Wizard",
                description="Add another restraint ramp segment?",
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
            self.console.print("[cyan]9.[/cyan] TAUTP - Temperature coupling time")
            self.console.print("[cyan]10.[/cyan] GAMMA_LN - Langevin collision frequency")
            self.console.print("[cyan]11.[/cyan] DUMPFREQ - DUMPAVE output frequency")
            self.console.print("[cyan]12.[/cyan] NSTEP0 - Shift step counter")
            self.console.print("[cyan]13.[/cyan] DISAVE - Distance restraint averaging")
            self.console.print("[cyan]14.[/cyan] RSTAR - Reference distance (NOESY)")
            self.console.print("[cyan]15.[/cyan] Custom parameter")

            param_choice = prompt_with_context(
                self.processor,
                "Select parameter",
                choices=["1", "2", "3", "4", "5", "6", "7", "8", "9",
                         "10", "11", "12", "13", "14", "15"],
                default="1",
                module="AMBER Wizard",
                description=f"Custom &wt block {block_num}: parameter to vary",
                options_map={
                    "1": "TEMP0 — Target temperature",
                    "2": "REST — NMR restraint weights",
                    "3": "BOND — Bond energy weights",
                    "4": "ANGLE — Angle energy weights",
                    "5": "TORSION — Torsion energy weights",
                    "6": "VDW — van der Waals weights",
                    "7": "ELEC — Electrostatic weights",
                    "8": "NB — All non-bonded weights",
                    "9": "TAUTP — Temperature coupling time",
                    "10": "GAMMA_LN — Langevin collision frequency",
                    "11": "DUMPFREQ — DUMPAVE output frequency",
                    "12": "NSTEP0 — Shift step counter",
                    "13": "DISAVE — Distance restraint averaging",
                    "14": "RSTAR — Reference distance (NOESY)",
                    "15": "Custom parameter",
                },
            )

            param_types = {
                "1": "TEMP0", "2": "REST", "3": "BOND", "4": "ANGLE",
                "5": "TORSION", "6": "VDW", "7": "ELEC", "8": "NB", "9": "TAUTP",
                "10": "GAMMA_LN", "11": "DUMPFREQ", "12": "NSTEP0",
                "13": "DISAVE", "14": "RSTAR"
            }
            
            if param_choice == "15":
                param_type = prompt_with_context(
                    self.processor,
                    "Parameter name",
                    module="AMBER Wizard",
                    description="Custom &wt parameter name",
                ).upper()
            else:
                param_type = param_types[param_choice]
            
            # Step range
            istep1 = int_prompt_with_context(
                self.processor,
                "Start step (istep1)",
                default=0,
                module="AMBER Wizard",
                description=f"Custom &wt block {block_num}: start step (istep1)",
            )
            istep2 = int_prompt_with_context(
                self.processor,
                "End step (istep2, 0=until end)",
                default=0,
                module="AMBER Wizard",
                description=f"Custom &wt block {block_num}: end step (istep2, 0=until end)",
            )

            # Values
            value1 = float_prompt_with_context(
                self.processor,
                "Initial value (value1)",
                default=0.0,
                module="AMBER Wizard",
                description=f"Custom &wt block {block_num}: initial value (value1)",
            )

            if istep2 > 0:
                value2 = float_prompt_with_context(
                    self.processor,
                    "Final value (value2)",
                    default=1.0,
                    module="AMBER Wizard",
                    description=f"Custom &wt block {block_num}: final value (value2)",
                )

                # Interpolation type
                linear = confirm_with_context(
                    self.processor,
                    "Linear interpolation (vs multiplicative)?",
                    default=True,
                    module="AMBER Wizard",
                    description=f"Custom &wt block {block_num}: linear (imult=0) vs multiplicative (imult=1)?",
                )
                imult = 0 if linear else 1

                # Step size
                iinc = int_prompt_with_context(
                    self.processor,
                    "Step increment (0=continuous)",
                    default=0,
                    module="AMBER Wizard",
                    description=f"Custom &wt block {block_num}: step increment (iinc, 0=continuous)",
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
                module="AMBER Wizard",
                description="Add another custom &wt block?",
            ):
                break
            
            block_num += 1
        
        return wt_blocks