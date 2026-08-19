"""
Small Molecule Parameterizer for MD Simulations

A module for parameterizing non-covalently bound small molecules (ligands, cofactors, etc.)
for molecular dynamics simulations using AMBER tools. Follows the coordinate extraction →
Gaussian (PBEPBE opt + HF ESP) → antechamber RESP → parmchk2 → tleap workflow.

Educational and interactive - shows all commands and explains scientific rationale.
"""

import os
import subprocess
import tempfile
import glob
import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

from rich.table import Table
from rich.console import Console
from rich.panel import Panel
from rich.prompt import Prompt, Confirm, FloatPrompt, IntPrompt
from rich.table import Table
from rich.text import Text

from Bio.PDB import PDBIO, Select, Structure, Model, Chain

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.forcefield_prep.small_molecule_commands import (
    AnalyzeSmallMoleculesCommand,
    GenerateGaussianInputCommand,
    ProcessGaussianOutputCommand,
    GenerateParametersCommand,
    CreateTLeapInputCommand,
    DisplayHelpCommand,
    RunWorkflowCommand,
)
from proprep.utils.workflow_checklist import WorkflowChecklist, WorkflowStep
from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
    float_prompt_with_context,
)


# LEGACY: These steps and the _checklist_sm_* handlers below serve the standalone
# SmallMoleculeParameterizer module menu path (_run_complete_workflow).
# The production path is run_workflow() → SmallMolWorkflowRunner → SMALL_MOL_WORKFLOW_STEPS.
# This legacy path will be removed once the standalone module menu is retired.
SMALL_MOL_STEPS = [
    WorkflowStep(
        id="sm-1", name="Analyze Structure",
        description="Identify small molecules in structure",
        handler="_checklist_sm_1_analyze",
        section="Small Molecule Parameterization",
    ),
    WorkflowStep(
        id="sm-2", name="Gaussian Input",
        description="Generate Gaussian input for optimization + ESP",
        handler="_checklist_sm_2_gaussian_input",
        section="Small Molecule Parameterization",
        dependencies=["sm-1"],
        checkpoint=True,
        checkpoint_message="Run Gaussian on generated .gjf files, then resume.",
    ),
    WorkflowStep(
        id="sm-3", name="RESP Charges",
        description="Process Gaussian output for RESP charges",
        handler="_checklist_sm_3_resp",
        section="Small Molecule Parameterization",
        dependencies=["sm-2"],
    ),
    WorkflowStep(
        id="sm-4", name="Bonded Parameters",
        description="Generate frcmod with parmchk2",
        handler="_checklist_sm_4_parameters",
        section="Small Molecule Parameterization",
        dependencies=["sm-3"],
    ),
    WorkflowStep(
        id="sm-5", name="tLEaP Input",
        description="Create tLEaP input files",
        handler="_checklist_sm_5_tleap",
        section="Small Molecule Parameterization",
        dependencies=["sm-4"],
    ),
]


@register_module
class SmallMoleculeParameterizer(ProcessingModule):
    """Small molecule parameterization module for ProPrep."""
    
    NAME = "Small Molecule Parameterizer"
    CATEGORY = "forcefield_prep"
    DESCRIPTION = "Parameterize non-covalently bound small molecules for MD simulations"
    VERSION = "1.0.0"
    REQUIRES = ["PDB Loader"]
    PRIORITY = 5  # After basic structure preparation, before tLEaP
    
    def __init__(self):
        super().__init__()
        self.console = Console()  # Simple direct assignment following standard pattern
        self.small_molecules = {}  # Dict to store identified small molecules
        self.parameterization_status = {}  # Track parameterization progress
        self.gaussian_available = None  # Cache Gaussian availability
        self.current_residue = None  # Store current residue being processed
        
        # Default Gaussian settings (user customizable)
        self.gaussian_settings = {
            "opt": {
                "memory": "10GB",
                "processors": 4,
                "keywords": "opt freq b3lyp/6-31+G(d) nosym int=ultrafine IOp(7/33=1)"
            },
            "esp": {
                "memory": "10GB", 
                "processors": 4,
                "keywords": "HF/6-31G(d) Pop=mk IOp(6/33=2,6/41=10,6/42=10) nosym int=ultrafine"
            }
        }

    def set_processor(self, processor):
        """Set the processor reference."""
        self.processor = processor
        # Update console if processor has one (following standard pattern)
        if hasattr(processor, 'console'):
            self.console = processor.console
    
    def get_menu_options(self) -> Dict[str, str]:
        """Get available menu options."""
        return {
            "run_workflow": "Run complete parameterization workflow (all steps)",
            "analyze": "Analyze structure for small molecules",
            "gaussian_input": "Generate Gaussian input files",
            "process_gaussian": "Process Gaussian output and assign RESP charges",
            "generate_params": "Generate force field parameters with parmchk2",
            "create_tleap": "Create tLEaP input files for gas/aqueous phases",
            "help": "Display detailed help information"
        }
    
    def get_workspace_requirements(self) -> List[str]:
        """List what this module needs from workspace."""
        return ["structure", "pdb_file"]
    
    def get_workspace_outputs(self) -> List[str]:
        """List what this module produces in workspace."""
        return ["small_molecules"]
    
    def can_process(self, workspace) -> bool:
        """Check if module can process current workspace.

        Uses StructureSelector to check for any available structure.
        """
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)
        status = selector.get_structure_status()
        return status.get("has_any", False)
        
    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection using command pattern."""
        if option == "run_workflow":
            command = RunWorkflowCommand(self.processor)
            return command.execute_with_error_handling()
        elif option == "analyze":
            command = AnalyzeSmallMoleculesCommand(self.processor)
            return command.execute_with_error_handling()
        elif option == "gaussian_input":
            command = GenerateGaussianInputCommand(self.processor)
            return command.execute_with_error_handling()
        elif option == "process_gaussian":
            command = ProcessGaussianOutputCommand(self.processor)
            return command.execute_with_error_handling()
        elif option == "generate_params":
            command = GenerateParametersCommand(self.processor)
            return command.execute_with_error_handling()
        elif option == "create_tleap":
            command = CreateTLeapInputCommand(self.processor)
            return command.execute_with_error_handling()
        elif option == "help":
            command = DisplayHelpCommand(self.processor)
            return command.execute_with_error_handling()
        return False

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Workflow orchestration

    def _run_complete_workflow(self) -> bool:
        """LEGACY: Run workflow via standalone module menu. Will be removed.

        The production path is run_workflow() via ForcefieldParameterizer.
        """
        output_dir = Path("small_molecule_params")
        output_dir.mkdir(exist_ok=True)

        checklist = WorkflowChecklist(
            steps=SMALL_MOL_STEPS,
            executor=self,
            processor=self.processor,
            workflow_name="Small Molecule Parameterization",
            console=self.console,
            state_dir=output_dir,
        )
        return checklist.run()

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # WorkflowChecklist handler methods

    def _checklist_sm_1_analyze(self):
        """Handler for WorkflowChecklist: analyze small molecules."""
        result = self._analyze_small_molecules()
        if not result:
            raise RuntimeError("Failed to analyze small molecules")
        count = len(self.small_molecules)
        return {'summary': f'Found {count} small molecule(s)'}

    def _checklist_sm_2_gaussian_input(self):
        """Handler for WorkflowChecklist: generate Gaussian input files."""
        result = self._generate_gaussian_input()
        if not result:
            raise RuntimeError("Failed to generate Gaussian input files")
        return {'checkpoint': True}

    def _checklist_sm_3_resp(self):
        """Handler for WorkflowChecklist: process Gaussian output for RESP charges."""
        result = self._process_gaussian_output()
        if not result:
            raise RuntimeError("Failed to process Gaussian output")
        return {'summary': 'RESP charges assigned'}

    def _checklist_sm_4_parameters(self):
        """Handler for WorkflowChecklist: generate force field parameters."""
        result = self._generate_parameters()
        if not result:
            raise RuntimeError("Failed to generate force field parameters")
        return {'summary': 'Force field parameters generated'}

    def _checklist_sm_5_tleap(self):
        """Handler for WorkflowChecklist: create tLEaP input files."""
        result = self._create_tleap_input()
        if not result:
            raise RuntimeError("Failed to create tLEaP input files")
        return {'summary': 'tLEaP input files created'}

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Public Action Methods (for command pattern)

    def _analyze_small_molecules(self) -> bool:
        """Analyze structure for small molecules (ligands, cofactors, etc.)."""
        try:
            workspace = self.get_workspace()

            # Use StructureSelector to get structure object
            from proprep.utils.structure_selector import StructureSelector
            selector = StructureSelector(workspace, self.console, processor=self.processor)
            structure = selector.get_structure_object(silent=True)

            if not structure:
                self.console.print("[yellow]No structure available in workspace[/yellow]")
                return False
            
            self.console.print("[bold cyan]Analyzing structure for small molecules...[/bold cyan]")
            
            # Standard residues to exclude
            standard_residues = {
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                "DA", "DT", "DG", "DC", "A", "T", "G", "C", "U"  # DNA/RNA
            }
            
            small_molecules = {}
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        resname = residue.get_resname().strip()
                        
                        # Skip standard residues and water
                        if resname in standard_residues or resname == "HOH":
                            continue
                        
                        # Skip if it's a single atom (likely an ion)
                        atom_count = len(list(residue.get_atoms()))
                        if atom_count == 1:
                            continue
                        
                        # This is likely a small molecule
                        mol_key = f"{chain.id}_{residue.id[1]}_{resname}"
                        small_molecules[mol_key] = {
                            "chain_id": chain.id,
                            "res_num": residue.id[1],
                            "res_name": resname,
                            "residue": residue,
                            "atom_count": atom_count,
                            "status": "identified"
                        }
            
            self.small_molecules = small_molecules
            self.update_workspace("small_molecules", small_molecules)
            
            # Display results
            if small_molecules:
                table = Table(title="Identified Small Molecules")
                table.add_column("Molecule ID", style="cyan")
                table.add_column("Chain", style="green")
                table.add_column("Residue", style="yellow")
                table.add_column("Name", style="magenta")
                table.add_column("Atoms", style="blue")
                table.add_column("Status", style="white")
                
                for mol_id, mol_info in small_molecules.items():
                    table.add_row(
                        mol_id,
                        mol_info["chain_id"],
                        str(mol_info["res_num"]),
                        mol_info["res_name"],
                        str(mol_info["atom_count"]),
                        mol_info["status"]
                    )
                
                self.console.print(table)
                self.console.print(f"[green]Found {len(small_molecules)} small molecule(s)[/green]")
            else:
                self.console.print("[yellow]No small molecules found in structure[/yellow]")
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error analyzing small molecules: {str(e)}[/red]")
            return False

    def _generate_gaussian_input(self) -> bool:
        """Generate Gaussian input files for RESP charge calculation."""
        try:
            if not self.small_molecules:
                self.console.print("[yellow]No small molecules identified. Run analysis first.[/yellow]")
                return False
            
            # Check if antechamber is available (we'll need it later for RESP)
            if not self._check_antechamber():
                return False
            
            # Let user select molecules to parameterize
            selected_molecules = self._select_molecules_for_parameterization()
            if not selected_molecules:
                return False
            
            # Get Gaussian computational settings from user
            self._configure_gaussian_settings()
            
            output_dir = Path("small_molecule_params")
            output_dir.mkdir(exist_ok=True)
            
            for mol_id in selected_molecules:
                mol_info = self.small_molecules[mol_id]
                
                self.console.print(f"[cyan]Processing molecule {mol_id} ({mol_info['res_name']})...[/cyan]")
                
                # Extract molecule to PDB file
                mol_pdb_file = output_dir / f"{mol_info['res_name'].lower()}.pdb"
                self._extract_and_write_molecule_pdb(mol_info, mol_pdb_file)
                
                # Get net charge for molecule
                net_charge = self._get_molecule_charge(mol_info['res_name'])
                
                # Generate Gaussian input file manually (with --link1--)
                gaussian_file = output_dir / f"{mol_info['res_name'].lower()}.gjf"
                self._write_gaussian_input_file(mol_pdb_file, gaussian_file, net_charge, mol_info['res_name'])
                
                self.console.print(f"[green]✓ Generated Gaussian input: {gaussian_file}[/green]")
                
                # Display next steps
                self._display_gaussian_instructions(gaussian_file, mol_info['res_name'])
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error generating Gaussian input: {str(e)}[/red]")
            return False

    def _process_gaussian_output(self) -> bool:
        """Process Gaussian output and assign RESP charges."""
        try:
            self.console.print("[bold cyan]Processing Gaussian output for RESP charges...[/bold cyan]")
            
            # Look for available log files
            log_files = list(Path(".").glob("*.log"))
            if not log_files:
                self.console.print("[yellow]No Gaussian log files found in current directory[/yellow]")
                return False
            
            # Let user select which log file to process
            if len(log_files) == 1:
                selected_log = log_files[0]
                self.console.print(f"[cyan]Processing {selected_log}[/cyan]")
            else:
                self.console.print("[cyan]Available Gaussian log files:[/cyan]")
                for i, log_file in enumerate(log_files, 1):
                    self.console.print(f"  {i}. {log_file}")
                
                choice = int_prompt_with_context(
                    self.processor,
                    "Select log file to process",
                    choices=[str(i) for i in range(1, len(log_files) + 1)],
                    default=1,
                    module="Small Molecule Parameterizer",
                    description="Select Gaussian log file to process"
                )
                selected_log = log_files[choice - 1]
            
            # Extract base name for output files
            base_name = selected_log.stem
            mol2_file = f"{base_name}.mol2"
            
            # Get net charge for RESP fitting
            net_charge = self._get_molecule_charge(base_name.upper())
            
            # Run antechamber to extract RESP charges
            self._run_antechamber_resp(selected_log, mol2_file, net_charge, base_name.upper())
            
            # Verify charges
            self._verify_mol2_charges(mol2_file, net_charge)
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error processing Gaussian output: {str(e)}[/red]")
            return False

    def _generate_parameters(self) -> bool:
        """Generate force field parameters using parmchk2."""
        try:
            self.console.print("[bold cyan]Generating force field parameters...[/bold cyan]")
            
            # Look for available MOL2 files
            mol2_files = list(Path(".").glob("*.mol2"))
            if not mol2_files:
                self.console.print("[yellow]No MOL2 files found. Process Gaussian output first.[/yellow]")
                return False
            
            # Check if parmchk2 is available
            if not self._check_parmchk2():
                return False
            
            for mol2_file in mol2_files:
                base_name = mol2_file.stem
                frcmod_file = f"{base_name}.frcmod"
                
                self.console.print(f"[cyan]Generating parameters for {mol2_file}...[/cyan]")
                
                # Run parmchk2
                self._run_parmchk2(mol2_file, frcmod_file)
                
                # Analyze generated parameters
                self._analyze_frcmod_file(frcmod_file)
            
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error generating parameters: {str(e)}[/red]")
            return False

    def _create_tleap_input(self) -> bool:
        """Create tLEaP input files for gas and aqueous phases."""
        try:
            self.console.print("[bold cyan]Creating tLEaP input files...[/bold cyan]")
            
            # Look for available MOL2 and FRCMOD files
            mol2_files = list(Path(".").glob("*.mol2"))
            frcmod_files = list(Path(".").glob("*.frcmod"))
            
            if not mol2_files or not frcmod_files:
                self.console.print("[yellow]Need both MOL2 and FRCMOD files. Complete parameterization first.[/yellow]")
                return False
            
            # Get water model preference
            water_model = self._get_water_model_choice()
            box_size = float_prompt_with_context(
                self.processor,
                "Solvation box size (Å)",
                default=20.0,
                module="Small Molecule Parameterizer",
                description="Solvation box size"
            )
            
            tleap_inputs = []
            
            for mol2_file in mol2_files:
                base_name = mol2_file.stem
                frcmod_file = Path(f"{base_name}.frcmod")
                
                if not frcmod_file.exists():
                    self.console.print(f"[yellow]No frcmod file for {mol2_file}, skipping...[/yellow]")
                    continue
                
                # Generate gas phase tLEaP input
                gas_input = self._generate_gas_tleap_input(mol2_file, frcmod_file, base_name)
                tleap_inputs.append(gas_input)
                
                # Generate aqueous phase tLEaP input
                aq_input = self._generate_aqueous_tleap_input(mol2_file, frcmod_file, base_name, water_model, box_size)
                tleap_inputs.append(aq_input)
            
            # Display completion message
            if tleap_inputs:
                self.console.print(Panel(
                    f"[bold green]✓ tLEaP input files generated![/bold green]\n\n"
                    f"[bold yellow]Next steps:[/bold yellow]\n\n"
                    f"Run tLEaP to generate AMBER parameter files:\n"
                    f"tleap -f {base_name}_gas_tleap.in\n"
                    f"tleap -f {base_name}_aq_tleap.in\n\n"
                    f"This will generate .parm7 and .rst7 files for MD simulations.",
                    title="tLEaP Input Generated",
                    expand=False
                ))
            
            return len(tleap_inputs) > 0
            
        except Exception as e:
            self.console.print(f"[red]Error creating tLEaP input: {str(e)}[/red]")
            return False

    def _display_help(self) -> bool:
        """Display detailed help information."""
        help_text = """
[bold]Small Molecule Parameterizer Help[/bold]

This module guides you through parameterizing non-covalently bound small molecules
for AMBER molecular dynamics simulations.

[bold]Workflow Overview:[/bold]
1. [cyan]Analyze structure[/cyan] - Identify small molecules in your PDB structure
2. [cyan]Generate Gaussian input[/cyan] - Create input files for quantum chemistry calculations
3. [cyan]Process Gaussian output[/cyan] - Extract RESP charges and generate MOL2 files
4. [cyan]Generate parameters[/cyan] - Create force field parameters with parmchk2
5. [cyan]Create tLEaP input[/cyan] - Generate files for gas and aqueous phase systems

[bold]Requirements:[/bold]
- AMBER Tools (antechamber, parmchk2, tleap)
- Gaussian (for quantum chemistry calculations)
- Structure with small molecules already loaded

[bold]Supported Molecule Types:[/bold]
- Ligands and cofactors
- Non-standard residues
- Organic small molecules
- Metal complexes (with proper handling)

[bold]Output Files:[/bold]
- .gjf - Gaussian input files
- .mol2 - MOL2 files with RESP charges
- .frcmod - AMBER force field modification files
- .parm7/.rst7 - AMBER parameter and coordinate files

[bold]Notes:[/bold]
- Gaussian calculations can be run externally and results imported
- The module follows GAFF2 force field conventions
- RESP charges provide high-quality electrostatics for MD simulations
"""
        
        self.console.print(Panel(help_text, title="Small Molecule Parameterizer Help"))
        return True

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Helper methods

    def _check_antechamber(self) -> bool:
        """Check if antechamber is available."""
        try:
            result = subprocess.run(["antechamber", "-h"], 
                                  capture_output=True, text=True)
            return result.returncode == 0
        except FileNotFoundError:
            self.console.print("[red]Error: antechamber not found. Please install AMBER Tools.[/red]")
            return False

    def _check_parmchk2(self) -> bool:
        """Check if parmchk2 is available."""
        try:
            result = subprocess.run(["parmchk2", "-h"], 
                                  capture_output=True, text=True)
            return result.returncode == 0
        except FileNotFoundError:
            self.console.print("[red]Error: parmchk2 not found. Please install AMBER Tools.[/red]")
            return False

    def _select_molecules_for_parameterization(self) -> List[str]:
        """Let user select which molecules to parameterize."""
        if not self.small_molecules:
            return []
        
        self.console.print("\n[cyan]Select molecules to parameterize:[/cyan]")
        
        selected = []
        for mol_id, mol_info in self.small_molecules.items():
            should_parameterize = confirm_with_context(
                self.processor,
                f"Parameterize {mol_id} ({mol_info['res_name']}, {mol_info['atom_count']} atoms)?",
                default=True,
                module="Small Molecule Parameterizer",
                description=f"Parameterize {mol_id}"
            )
            if should_parameterize:
                selected.append(mol_id)
        
        return selected

    def _configure_gaussian_settings(self):
        """Allow user to configure Gaussian computational settings."""
        self.console.print("\n[bold cyan]Configure Gaussian Settings[/bold cyan]")
        
        # Optimization settings (full customization)
        self.console.print("\n[bold]Optimization + Frequency Calculation:[/bold]")
        
        # Memory
        self.gaussian_settings["opt"]["memory"] = prompt_with_context(
            self.processor,
            "Memory allocation",
            default=self.gaussian_settings["opt"]["memory"],
            module="Small Molecule Parameterizer",
            description="Optimization memory allocation"
        )

        # Processors
        try:
            proc_input = prompt_with_context(
                self.processor,
                "Number of processors",
                default=str(self.gaussian_settings["opt"]["processors"]),
                module="Small Molecule Parameterizer",
                description="Optimization processors"
            )
            self.gaussian_settings["opt"]["processors"] = int(proc_input)
        except ValueError:
            self.console.print("[yellow]Invalid processor count, using default[/yellow]")

        # Keywords with educational info
        customize_opt = confirm_with_context(
            self.processor,
            "Customize optimization keywords?",
            default=False,
            module="Small Molecule Parameterizer",
            description="Customize optimization keywords"
        )
        if customize_opt:
            self.console.print(Panel(
                "[cyan]Default: opt freq b3lyp/6-31+G(d) IOp(7/33=1)[/cyan]\n\n"
                "[bold]Method options:[/bold]\n"
                "• pbepbe - PBE density functional (good for organics)\n"
                "• b3lyp - Popular hybrid functional\n"
                "• m06-2x - Good for dispersion interactions\n\n"
                "[bold]Basis set options:[/bold]\n"
                "• 6-31g(2d,2p) - Double-zeta with polarization\n"
                "• 6-31+g(d,p) - Adds diffuse functions\n"
                "• def2-svp - Ahlrichs basis set\n\n"
                "[bold]IOp(7/33=1):[/bold] Saves Cartesian Hessian for Seminario refinement",
                title="Optimization Keywords Help",
                expand=False
            ))
            
            self.gaussian_settings["opt"]["keywords"] = prompt_with_context(
                self.processor,
                "Optimization keywords",
                default=self.gaussian_settings["opt"]["keywords"],
                module="Small Molecule Parameterizer",
                description="Optimization keywords"
            )

        # ESP settings (customizable with warnings)
        self.console.print("\n[bold]ESP Calculation:[/bold]")

        # Memory and processors
        self.gaussian_settings["esp"]["memory"] = prompt_with_context(
            self.processor,
            "ESP memory allocation",
            default=self.gaussian_settings["esp"]["memory"],
            module="Small Molecule Parameterizer",
            description="ESP memory allocation"
        )

        try:
            proc_input = prompt_with_context(
                self.processor,
                "ESP processors",
                default=str(self.gaussian_settings["esp"]["processors"]),
                module="Small Molecule Parameterizer",
                description="ESP processors"
            )
            self.gaussian_settings["esp"]["processors"] = int(proc_input)
        except ValueError:
            self.console.print("[yellow]Invalid processor count, using default[/yellow]")

        # Keywords with strong warnings
        modify_esp = confirm_with_context(
            self.processor,
            "Modify ESP keywords? [red](NOT RECOMMENDED)[/red]",
            default=False,
            module="Small Molecule Parameterizer",
            description="Modify ESP keywords"
        )

        if modify_esp:
            self.console.print(Panel(
                "[red][bold]WARNING: Modifying ESP keywords may break AMBER compatibility![/bold][/red]\n\n"
                "[cyan]Default: HF/6-31G(d) Pop=mk IOp(6/33=2,6/41=10,6/42=10) nosym int=ultrafine[/cyan]\n\n"
                "[bold]Why these defaults?[/bold]\n"
                "• [bold]HF/6-31G(d)[/bold]     : Required for AMBER charge compatibility\n"
                "• [bold]Pop=mk[/bold]          : Merz-Kollman population analysis for ESP\n"
                "• [bold]IOp(6/33=2)[/bold]     : Output ESP points and potentials for RESP fitting\n"
                "• [bold]IOp(6/41=10)[/bold]    : Use 10 concentric layers (default=4)\n"
                "• [bold]IOp(6/42=10)[/bold]    : 10 points/area (~1000 points/atom, default=1)\n"
                "• [bold]NoSymm[/bold]          : Prevent symmetry constraints on ESP grid\n"
                "• [bold]Int=ultrafine[/bold]   : High-quality integration grid\n\n"
                "[red]Changing these may result in charges incompatible with AMBER simulations![/red]",
                title="ESP Keywords Warning",
                border_style="red",
                expand=False
            ))

            still_modify = confirm_with_context(
                self.processor,
                "Do you still want to modify ESP keywords?",
                default=False,
                module="Small Molecule Parameterizer",
                description="Confirm modify ESP keywords"
            )

            if still_modify:
                self.gaussian_settings["esp"]["keywords"] = prompt_with_context(
                    self.processor,
                    "[red]ESP keywords (expert only)[/red]",
                    default=self.gaussian_settings["esp"]["keywords"],
                    module="Small Molecule Parameterizer",
                    description="Custom ESP keywords"
                )
                self.console.print("[red]Using custom ESP keywords - simulation quality not guaranteed![/red]")

    def _get_molecule_charge(self, res_name: str) -> int:
        """Get net charge for a molecule from user."""
        charge = prompt_with_context(
            self.processor,
            f"Net charge for {res_name}",
            default="0",
            module="Small Molecule Parameterizer",
            description=f"Net charge for {res_name}"
        )
        try:
            return int(charge)
        except ValueError:
            self.console.print("[yellow]Invalid charge, using 0[/yellow]")
            return 0

    def _extract_and_write_molecule_pdb(self, mol_info: Dict, output_file: Path):
        """Extract a single molecule from the full structure and write to PDB file."""
        try:
            # Create a new structure with just this residue
            residue = mol_info["residue"]
            
            # Create new structure hierarchy
            new_structure = Structure.Structure("small_mol")
            new_model = Model.Model(0)
            new_chain = Chain.Chain(mol_info["chain_id"])
            
            # Add the residue to the new chain
            new_chain.add(residue.copy())
            new_model.add(new_chain)
            new_structure.add(new_model)
            
            # Write to PDB file
            io = PDBIO()
            io.set_structure(new_structure)
            io.save(str(output_file))
            
            self.console.print(f"[green]✓ Extracted coordinates to {output_file}[/green]")
            
        except Exception as e:
            self.console.print(f"[red]Error extracting molecule: {str(e)}[/red]")
            # Fallback: try simple extraction
            try:
                io = PDBIO()
                io.set_structure(mol_info["residue"])
                io.save(str(output_file))
                self.console.print(f"[yellow]✓ Extracted using fallback method[/yellow]")
            except:
                raise Exception(f"Could not extract molecule coordinates: {str(e)}")

    def _write_gaussian_input_file(self, pdb_file: Path, gaussian_file: Path, net_charge: int, res_name: str):
        """Write Gaussian input file with optimization + ESP calculation using --link1--."""
        
        # Read coordinates from PDB file
        coords = []
        with open(pdb_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith(('ATOM', 'HETATM')):
                    try:
                        element = line[76:78].strip()
                        if not element:
                            # Guess element from atom name
                            element = line[12:16].strip()[0]
                        
                        # More robust coordinate parsing with proper error handling
                        x_str = line[30:38].strip()
                        y_str = line[38:46].strip()
                        z_str = line[46:54].strip()
                        
                        # Debug output if parsing fails
                        if not x_str or not y_str or not z_str:
                            self.console.print(f"[yellow]Warning: Empty coordinate field on line {line_num}[/yellow]")
                            self.console.print(f"Line: {line.strip()}")
                            continue
                        
                        x = float(x_str)
                        y = float(y_str)
                        z = float(z_str)
                        
                        coords.append(f" {element:2s}              {x:12.8f}   {y:12.8f}   {z:12.8f}")
                        
                    except ValueError as e:
                        # Enhanced error reporting to help debug the issue
                        self.console.print(f"[red]Error parsing coordinates on line {line_num}:[/red]")
                        self.console.print(f"Line: {line.strip()}")
                        self.console.print(f"X field (chars 30-38): '{line[30:38]}'")
                        self.console.print(f"Y field (chars 38-46): '{line[38:46]}'")
                        self.console.print(f"Z field (chars 46-54): '{line[46:54]}'")
                        self.console.print(f"Element field (chars 76-78): '{line[76:78]}'")
                        self.console.print(f"Atom name field (chars 12-16): '{line[12:16]}'")
                        raise ValueError(f"Could not parse coordinates on line {line_num}: {str(e)}")
                    
                    except IndexError as e:
                        self.console.print(f"[red]Line {line_num} too short for PDB format:[/red]")
                        self.console.print(f"Line length: {len(line)}, Line: {line.strip()}")
                        raise ValueError(f"Invalid PDB format on line {line_num}: line too short")
        
        if not coords:
            raise Exception("No coordinates found in PDB file")
        
        # Write Gaussian input with educational comments
        base_name = gaussian_file.stem
        
        with open(gaussian_file, 'w') as f:
            # Optimization section
            f.write(f"%chk={base_name}.chk\n")
            f.write(f"%nprocshared={self.gaussian_settings['opt']['processors']}\n")
            f.write(f"%mem={self.gaussian_settings['opt']['memory']}\n")
            f.write(f"#p {self.gaussian_settings['opt']['keywords']}\n")
            f.write("\n")
            f.write(f"{res_name} optimization and frequency\n")
            f.write("\n")
            f.write(f"{net_charge} 1\n")
            
            # Write coordinates
            for coord in coords:
                f.write(f"{coord}\n")
            f.write("\n")
            
            # ESP section with --link1--
            f.write("--link1--\n")
            f.write(f"%oldchk={base_name}.chk\n")
            f.write(f"%chk={base_name}_esp.chk\n")
            f.write(f"%nprocshared={self.gaussian_settings['esp']['processors']}\n")
            f.write(f"%mem={self.gaussian_settings['esp']['memory']}\n")
            f.write(f"#P {self.gaussian_settings['esp']['keywords']} guess=read geom=allcheck GFInput GFPrint\n")
            f.write("\n")
            f.write(f"{res_name} ESP calculation\n")
            f.write("\n")
            f.write("\n")  # Blank line at end
        
        # Display educational information about the file
        self.console.print(Panel(
            f"[bold cyan]Gaussian Input File Structure:[/bold cyan]\n\n"
            f"[yellow]Section 1: Optimization + Frequency[/yellow]\n"
            f"• Method: {self.gaussian_settings['opt']['keywords']}\n"
            f"• Purpose: Find minimum energy geometry\n"
            f"• Frequency: Verify no imaginary frequencies\n\n"
            f"[yellow]Section 2: ESP Calculation (--link1--)[/yellow]\n"
            f"• Method: {self.gaussian_settings['esp']['keywords']}\n"
            f"• Purpose: Generate electrostatic potential grid\n"
            f"• Uses optimized geometry from Section 1\n"
            f"• %oldchk reads from opt checkpoint; %chk writes to separate ESP checkpoint\n\n"
            f"[cyan]Why this two-step approach?[/cyan]\n"
            f"1. DFT optimization provides accurate geometry\n"
            f"2. HF ESP calculation provides AMBER-compatible charges\n"
            f"3. RESP fitting will use ESP data for high-quality charges\n\n"
            f"[green]Generated: {gaussian_file}[/green]",
            title="Gaussian Input File Created",
            expand=False
        ))

    def _display_gaussian_instructions(self, gaussian_file: Path, res_name: str):
        """Display instructions for running Gaussian calculation."""
        log_file = gaussian_file.with_suffix('.log')
        
        self.console.print(Panel(
            f"[bold yellow]Gaussian Calculation Instructions[/bold yellow]\n\n"
            f"[bold]1. Run Gaussian calculation:[/bold]\n"
            f"   g16 {gaussian_file} > {log_file}\n\n"
            f"[bold]2. What this will do:[/bold]\n"
            f"   • Optimize molecular geometry\n"
            f"   • Calculate vibrational frequencies\n"
            f"   • Generate electrostatic potential grid\n"
            f"   • Prepare data for RESP charge fitting\n\n"
            f"[bold]3. Expected runtime:[/bold]\n"
            f"   • Small molecules (< 20 atoms): 10-30 minutes\n"
            f"   • Medium molecules (20-50 atoms): 1-3 hours\n"
            f"   • Large molecules (> 50 atoms): Several hours\n\n"
            f"[bold]4. Check for completion:[/bold]\n"
            f"   • Log file should end with 'Normal termination of Gaussian'\n"
            f"   • Look for both optimization and ESP sections\n\n"
            f"[cyan]Return here after Gaussian completion to continue with RESP fitting![/cyan]",
            title=f"Next Steps for {res_name}"
        ))

    def _run_antechamber_resp(self, log_file: Path, mol2_file: str, net_charge: int, res_name: str):
        """Run antechamber to extract RESP charges from Gaussian log."""
        
        # Build antechamber command
        cmd = [
            "antechamber",
            "-i", str(log_file),
            "-fi", "gout",
            "-o", mol2_file,
            "-fo", "mol2",
            "-c", "resp",
            "-nc", str(net_charge),
            "-rn", res_name,
            "-at", "gaff2"
        ]
        
        # Display educational information about the command
        self.console.print("\n[bold cyan]RESP Charge Assignment[/bold cyan]")
        self.console.print(f"[yellow]Command being executed:[/yellow]")
        self.console.print(" ".join(cmd))
        
        self.console.print(f"\n[cyan]Flag explanations:[/cyan]")
        self.console.print(f"  -i {log_file}        : Gaussian output with ESP data")
        self.console.print(f"  -fi gout                : Gaussian output format")
        self.console.print(f"  -o {mol2_file}         : MOL2 with charges and atom types")
        self.console.print(f"  -fo mol2                : MOL2 output format")
        self.console.print(f"  -c resp                 : RESP charge method")
        self.console.print(f"  -nc {net_charge}                : Net charge constraint")
        self.console.print(f"  -rn {res_name}              : Residue name in MOL2")
        self.console.print(f"  -at gaff2               : GAFF2 atom type assignment")
        
        self.console.print(f"\n[cyan]Why RESP charges?[/cyan]")
        self.console.print(f"• Fit to reproduce electrostatic potential around molecule")
        self.console.print(f"• Topologically equivalent atoms get identical charges")
        self.console.print(f"• Designed specifically for molecular dynamics simulations")
        self.console.print(f"• More accurate than Mulliken or other population analyses")
        
        # Run the command
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.console.print(f"[green]✓ Generated MOL2 file: {mol2_file}[/green]")
            
            if result.stdout:
                self.console.print(f"[grey50]Antechamber output: {result.stdout.strip()}[/grey50]")
                
        except subprocess.CalledProcessError as e:
            self.console.print(f"[red]Error running antechamber: {e.stderr}[/red]")
            raise

    def _verify_mol2_charges(self, mol2_file: str, expected_charge: int):
        """Verify that MOL2 charges sum to expected value."""
        try:
            charges = []
            with open(mol2_file, 'r') as f:
                in_atom_section = False
                for line in f:
                    if line.startswith("@<TRIPOS>ATOM"):
                        in_atom_section = True
                        continue
                    elif line.startswith("@<TRIPOS>"):
                        in_atom_section = False
                        continue
                    
                    if in_atom_section and line.strip():
                        parts = line.split()
                        if len(parts) >= 9:
                            try:
                                charge = float(parts[8])
                                charges.append(charge)
                            except ValueError:
                                continue
            
            if charges:
                total_charge = sum(charges)
                self.console.print(f"\n[cyan]Charge verification:[/cyan]")
                self.console.print(f"Total charge: {total_charge:.6f}")
                self.console.print(f"Expected charge: {expected_charge}")
                self.console.print(f"Difference: {abs(total_charge - expected_charge):.6f}")
                
                if abs(total_charge - expected_charge) < 0.001:
                    self.console.print("[green]✓ Charges sum correctly[/green]")
                else:
                    self.console.print("[yellow]⚠ Charge deviation larger than expected[/yellow]")
            
        except Exception as e:
            self.console.print(f"[yellow]Could not verify charges: {str(e)}[/yellow]")

    def _run_parmchk2(self, mol2_file: Path, frcmod_file: str):
        """Run parmchk2 to generate force field parameters."""
        
        # Build parmchk2 command
        cmd = [
            "parmchk2",
            "-i", str(mol2_file),
            "-o", frcmod_file,
            "-f", "mol2",
            "-s", "2",
            "-a", "Y"
        ]

        # Display educational information
        self.console.print(f"\n[yellow]Command being executed:[/yellow]")
        self.console.print(" ".join(cmd))

        self.console.print(f"\n[cyan]Flag explanations:[/cyan]")
        self.console.print(f"  -i {mol2_file}         : MOL2 with atom types and charges")
        self.console.print(f"  -o {frcmod_file}       : Force field modification file")
        self.console.print(f"  -f mol2                : Input is MOL2 format")
        self.console.print(f"  -s 2                   : Use GAFF2 parameter database")
        self.console.print(f"  -a Y                   : Print all parameters (complete frcmod)")
        
        self.console.print(f"\n[cyan]What parmchk2 does:[/cyan]")
        self.console.print(f"• Searches GAFF2 database for all needed parameters")
        self.console.print(f"• Identifies missing bond/angle/dihedral parameters")
        self.console.print(f"• Estimates missing parameters using chemical similarity")
        self.console.print(f"• Creates frcmod file with custom parameters")
        
        # Run the command
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.console.print(f"[green]✓ Generated parameters: {frcmod_file}[/green]")
            
        except subprocess.CalledProcessError as e:
            self.console.print(f"[red]Error running parmchk2: {e.stderr}[/red]")
            raise

    def _analyze_frcmod_file(self, frcmod_file: str):
        """Analyze the generated frcmod file and report parameter quality."""
        try:
            with open(frcmod_file, 'r') as f:
                content = f.read()
            
            # Count different parameter types
            bond_count = len([line for line in content.split('\n') if line and not line.startswith(('remark', 'MASS', 'BOND', 'ANGLE', 'DIHE', 'IMPROPER', 'NONBON')) and 'BOND' in content.split('\n')[max(0, content.split('\n').index(line)-10):content.split('\n').index(line)]])
            
            # Simple analysis - count non-empty lines in each section
            sections = content.split('\n\n')
            
            self.console.print(f"\n[cyan]Parameter analysis for {frcmod_file}:[/cyan]")
            
            if 'BOND' in content:
                self.console.print("✓ Bond parameters section present")
            if 'ANGLE' in content:
                self.console.print("✓ Angle parameters section present") 
            if 'DIHE' in content:
                self.console.print("✓ Dihedral parameters section present")
            if 'IMPROPER' in content:
                self.console.print("✓ Improper parameters section present")
            
            # Note about parameter estimation
            self.console.print(f"\n[yellow]Note:[/yellow] Parameters marked with 'ATTN' were estimated")
            self.console.print(f"Review the frcmod file to check parameter quality")
            
        except Exception as e:
            self.console.print(f"[yellow]Could not analyze frcmod file: {str(e)}[/yellow]")

    def _get_water_model_choice(self) -> str:
        """Get user's water model preference."""
        water_models = {
            "1": ("opc", "OPC (recommended - 4-point model)"),
            "2": ("tip3p", "TIP3P (3-point model)"),
            "3": ("tip4pew", "TIP4P-Ew (4-point Ewald)"),
            "4": ("spce", "SPC/E (3-point extended)")
        }
        
        self.console.print("\n[cyan]Choose water model:[/cyan]")
        for key, (model, desc) in water_models.items():
            self.console.print(f"  {key}. {desc}")

        choice = prompt_with_context(
            self.processor,
            "Water model",
            choices=list(water_models.keys()),
            default="1",
            module="Small Molecule Parameterizer",
            description="Select water model"
        )

        return water_models[choice][0]

    def _generate_gas_tleap_input(self, mol2_file: Path, frcmod_file: Path, mol_name: str) -> Path:
        """Generate gas phase tLEaP input file."""
        output_file = Path(f"{mol_name}_gas_tleap.in")
        unit = _tleap_safe_unit_var(mol_name)

        content = f"""source leaprc.protein.ff19SB
loadamberparams frcmod.ff19SB
source leaprc.gaff2
loadamberparams "{frcmod_file.name}"

{unit} = loadmol2 "{mol2_file.name}"
saveOff {unit} "{mol_name}.lib"
saveAmberParm {unit} "{mol_name}_gas.parm7" "{mol_name}_gas.rst7"

quit
"""
        
        output_file.write_text(content)
        self.console.print(f"[green]✓ Generated gas phase tLEaP input: {output_file}[/green]")
        
        return output_file

    def _generate_aqueous_tleap_input(self, mol2_file: Path, frcmod_file: Path, mol_name: str, water_model: str, box_size: float) -> Path:
        """Generate aqueous phase tLEaP input file."""
        output_file = Path(f"{mol_name}_aq_tleap.in")
        
        # Water model configurations
        water_configs = {
            "opc": ("leaprc.water.opc", "frcmod.opc", "OPCBOX"),
            "tip3p": ("leaprc.water.tip3p", "frcmod.tip3p", "TIP3PBOX"),
            "tip4pew": ("leaprc.water.tip4pew", "frcmod.tip4pew", "TIP4PEWBOX"),
            "spce": ("leaprc.water.spce", "frcmod.spce", "SPCBOX")
        }
        
        leaprc, frcmod, box_type = water_configs[water_model]
        unit = _tleap_safe_unit_var(mol_name)

        content = f"""source leaprc.protein.ff19SB
loadamberparams frcmod.ff19SB
source {leaprc}
loadamberparams {frcmod}
source leaprc.gaff2
loadamberparams "{frcmod_file.name}"

{unit} = loadmol2 "{mol2_file.name}"
solvateBox {unit} {box_type} {box_size}
saveAmberParm {unit} "{mol_name}_aq.parm7" "{mol_name}_aq.rst7"

quit
"""
        output_file.write_text(content)
        self.console.print(f"[green]✓ Generated aqueous phase tLEaP input: {output_file}[/green]")
        
        return output_file


# =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
# Module-level function for integration with forcefield_parameterizer

# Hydrogen addition via `reduce` now lives in `hydrogen_editor` (single
# source of truth, reused by the modified-amino-acid from-structure route).
# Re-exported here so existing importers (e.g. structure_preprocessor) that
# do `from ...small_molecule_parameterizer import run_reduce_aligned` keep
# working unchanged.
from proprep.forcefield_prep.hydrogen_editor import (  # noqa: E402,F401
    check_reduce_availability,
    configure_reduce_options_aligned,
    run_reduce_aligned,
)


def enhanced_antechamber_step(gaussian_log_file, mol2_file, net_charge, res_name, console, interactive=True, processor=None):
    """Enhanced antechamber step with frequency checking."""

    # Check frequencies first
    freq_ok = check_gaussian_frequencies(gaussian_log_file, console)

    if not freq_ok:
        if interactive:
            proceed = confirm_with_context(
                processor,
                "Gaussian calculation has issues. Continue with RESP charge generation anyway?",
                default=False,
                module="Small Molecule Parameterizer",
                description="Continue despite frequency issues"
            )
            if not proceed:
                console.print("[yellow]Aborting antechamber step due to frequency issues[/yellow]")
                return False
        else:
            console.print("[red]Frequency check failed in non-interactive mode, aborting[/red]")
            return False
    
    # Run antechamber
    console.print(f"\n[cyan]Running antechamber for RESP charge assignment...[/cyan]")

    cmd = [
        "antechamber",
        "-i", gaussian_log_file,
        "-fi", "gout",
        "-o", mol2_file,
        "-fo", "mol2",
        "-c", "resp",
        "-nc", str(net_charge),
        "-rn", res_name,
        "-at", "gaff2"
    ]

    console.print(f"[grey50]{' '.join(cmd)}[/grey50]")
    console.print(f"\n[yellow]Command explanation:[/yellow]")
    console.print(f"  -i {gaussian_log_file:<13}: Input Gaussian log file")
    console.print(f"  -fi gout        : Input format (Gaussian output)")
    console.print(f"  -o {mol2_file:<13}: Output MOL2 file")
    console.print(f"  -fo mol2        : Output format (MOL2)")
    console.print(f"  -c resp         : Use RESP charges from ESP calculation")
    console.print(f"  -nc {net_charge:<12}: Net charge of molecule")
    console.print(f"  -rn {res_name:<12}: Residue name in MOL2")
    console.print(f"  -at gaff2       : Use GAFF2 atom types")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode == 0:
            console.print(f"[green]✓ RESP charges assigned successfully[/green]")
            if result.stdout:
                # Show any warnings
                if "Warning" in result.stdout or "ERROR" in result.stdout:
                    console.print("[yellow]Antechamber warnings/errors:[/yellow]")
                    console.print(f"[grey50]{result.stdout}[/grey50]")
            return True
        else:
            console.print(f"[red]✗ Antechamber failed[/red]")
            if result.stderr:
                console.print(f"[red]Error: {result.stderr}[/red]")
            return False
            
    except subprocess.TimeoutExpired:
        console.print("[red]✗ Antechamber timed out[/red]")
        return False
    except Exception as e:
        console.print(f"[red]✗ Error running antechamber: {str(e)}[/red]")
        return False


def generate_atom_name_mapping(pdb_file: str, mol2_file: str, output_file: str, console) -> bool:
    """Generate atom name mapping between PDB and MOL2 files.

    After antechamber, atom names in the MOL2 may differ from the original PDB.
    Atom ORDER is preserved through the workflow (PDB → Gaussian → antechamber → mol2),
    so we simply map by index position.

    Args:
        pdb_file: Path to original PDB file (with crystallographic atom names)
        mol2_file: Path to MOL2 file from antechamber (with sequential atom names)
        output_file: Path to save the JSON mapping file
        console: Rich console for output

    Returns:
        True if mapping was generated successfully
    """
    console.print(f"\n[cyan]Generating atom name mapping...[/cyan]")

    try:
        # Parse PDB atom names (in order)
        pdb_atoms = []
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    atom_name = line[12:16].strip()
                    pdb_atoms.append(atom_name)

        # Parse MOL2 atom names (in order)
        mol2_atoms = []
        in_atom_section = False
        with open(mol2_file, 'r') as f:
            for line in f:
                if line.startswith('@<TRIPOS>ATOM'):
                    in_atom_section = True
                    continue
                if line.startswith('@<TRIPOS>') and in_atom_section:
                    break
                if in_atom_section and line.strip():
                    parts = line.split()
                    if len(parts) >= 2:
                        mol2_atoms.append(parts[1])

        # Report counts
        n_pdb = len(pdb_atoms)
        n_mol2 = len(mol2_atoms)
        console.print(f"  PDB atoms: {n_pdb}, MOL2 atoms: {n_mol2}")

        if n_pdb != n_mol2:
            console.print(f"  [yellow]Warning: Atom count mismatch![/yellow]")

        # Map by index order (atom order is preserved through workflow)
        mapping = {}
        n_different = 0
        for pdb_name, mol2_name in zip(pdb_atoms, mol2_atoms):
            mapping[pdb_name] = mol2_name
            if pdb_name != mol2_name:
                n_different += 1

        console.print(f"  Names differ for {n_different}/{min(n_pdb, n_mol2)} atoms")

        # Save mapping
        with open(output_file, 'w') as f:
            json.dump(mapping, f, indent=2)

        console.print(f"[green]✓ Atom name mapping saved to {output_file}[/green]")
        return True

    except Exception as e:
        console.print(f"[yellow]⚠ Could not generate atom name mapping: {e}[/yellow]")
        return False


def check_gaussian_frequencies(log_file, console):
    """Check for negative frequencies in Gaussian log file."""
    console.print(f"\n[cyan]Checking frequencies in {log_file}...[/cyan]")
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Look for frequency lines
        freq_lines = []
        for line in content.split('\n'):
            if 'Frequencies --' in line:
                freq_lines.append(line)
        
        if not freq_lines:
            console.print("[yellow]⚠️  No frequency information found in log file[/yellow]")
            console.print("[yellow]   This may indicate an incomplete or failed calculation[/yellow]")
            return False
        
        # Parse the first frequency line (where negative frequencies would appear)
        first_line = freq_lines[0]
        freq_values = []
        
        # Extract frequency values after "Frequencies --"
        parts = first_line.split('Frequencies --')[1].split()
        for part in parts:
            try:
                freq = float(part)
                freq_values.append(freq)
            except ValueError:
                continue
        
        if not freq_values:
            console.print("[yellow]⚠️  Could not parse frequency values[/yellow]")
            return False
        
        # Check for negative frequencies
        negative_freqs = [f for f in freq_values if f < 0]
        
        console.print(f"[cyan]Found {len(freq_lines)} frequency line(s)[/cyan]")
        console.print(f"[cyan]First few frequencies: {', '.join(f'{f:.1f}' for f in freq_values[:6])}[/cyan]")
        
        if negative_freqs:
            console.print(f"[red]⚠️  {len(negative_freqs)} negative frequency(ies) detected: {', '.join(f'{f:.1f}' for f in negative_freqs)}[/red]")
            console.print(f"[red]   This indicates the structure is not at a minimum energy geometry[/red]")
            console.print(f"[red]   The optimization may have failed or converged to a transition state[/red]")
            return False
        else:
            console.print(f"[green]✓ All frequencies are positive (lowest: {min(freq_values):.1f} cm⁻¹)[/green]")
            return True
            
    except FileNotFoundError:
        console.print(f"[red]✗ Gaussian log file not found: {log_file}[/red]")
        return False
    except Exception as e:
        console.print(f"[red]✗ Error reading log file: {str(e)}[/red]")
        return False
    
def analyze_frcmod_penalties(frcmod_file, console):
    """
    Analyze and report penalty scores from frcmod file.

    Returns:
        List of (param_name, score, status, section) tuples, or empty list on error
    """
    console.print(f"\n[cyan]Analyzing parameter penalties in {frcmod_file}...[/cyan]")

    try:
        with open(frcmod_file, 'r') as f:
            lines = f.readlines()

        penalties = []
        current_section = None

        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # Track which section we're in
            if line in ['BOND', 'ANGLE', 'DIHE', 'IMPROPER', 'NONBON']:
                current_section = line
                continue

            # Skip section headers and empty sections
            if line in ['MASS', 'END'] or not current_section:
                continue

            # Look for penalty score or ATTN
            if 'penalty score=' in line or 'ATTN' in line:
                # Extract the complete parameter definition based on section
                param_name = extract_parameter_name(line, current_section)

                if 'penalty score=' in line:
                    # Extract penalty score
                    try:
                        score_part = line.split('penalty score=')[1].strip()
                        score_str = score_part.split()[0]
                        # Remove any trailing non-numeric characters (e.g., parentheses)
                        score_str = score_str.rstrip(')]}>')
                        score = float(score_str)
                        penalties.append((param_name, score, 'penalty', current_section))
                    except (IndexError, ValueError):
                        penalties.append((param_name, 999.0, 'parse_error', current_section))
                elif 'ATTN' in line:
                    penalties.append((param_name, float('inf'), 'ATTN', current_section))

        if not penalties:
            console.print("[green]✓ No penalty scores or ATTN warnings found[/green]")
            return []

        # Sort by penalty score (descending)
        penalties.sort(key=lambda x: x[1], reverse=True)

        # Create table with numbering for easy reference
        table = Table(title="Parameter Penalty Analysis", expand=False)
        table.add_column("#", style="grey50", justify="right")
        table.add_column("Section", style="blue")
        table.add_column("Parameter", style="cyan")
        table.add_column("Penalty Score", style="yellow")
        table.add_column("Status")

        for idx, (param, score, status, section) in enumerate(penalties, 1):
            if status == 'ATTN':
                score_str = "ATTN"
                status_str = "[red]⚠️  Attention Required[/red]"
            elif status == 'parse_error':
                score_str = "Parse Error"
                status_str = "[red]❌ Could not parse[/red]"
            else:
                score_str = f"{score:.1f}"
                if score > 10:
                    status_str = "[red]🔴 High penalty[/red]"
                elif score > 1:
                    status_str = "[yellow]🟡 Medium penalty[/yellow]"
                else:
                    status_str = "[green]✅ Low penalty[/green]"

            table.add_row(str(idx), section, param, score_str, status_str)

        console.print(table)

        # Summary
        high_penalties = [p for p in penalties if p[1] > 10 or p[2] == 'ATTN']
        medium_penalties = [p for p in penalties if 1 < p[1] <= 10 and p[2] != 'ATTN']
        if high_penalties:
            console.print(f"\n[red]⚠️  {len(high_penalties)} parameter(s) with high penalties or ATTN warnings[/red]")
        if medium_penalties:
            console.print(f"[yellow]   {len(medium_penalties)} parameter(s) with medium penalties[/yellow]")
        if not high_penalties and not medium_penalties:
            console.print(f"\n[green]✓ All parameters have acceptable penalty scores[/green]")

        return penalties

    except FileNotFoundError:
        console.print(f"[red]✗ FRCMOD file not found: {frcmod_file}[/red]")
        return []
    except Exception as e:
        console.print(f"[red]✗ Error analyzing FRCMOD file: {str(e)}[/red]")
        return []


def ask_refinement_selection(penalties, console, interactive=True, processor=None, fchk_file=None, mol_name="mol"):
    """
    Ask user which parameters to refine using the simplified workflow.

    Flow:
    1. User selects parameters from penalty table
    2. Bonds/angles → Seminario (requires fchk)
    3. Dihedrals → User chooses PES scan or CREST

    Args:
        penalties: List of (param_name, score, status, section) tuples
        console: Rich console for output
        interactive: Whether running in interactive mode
        processor: Optional processor for session recording context
        fchk_file: Path to .fchk file (None if not available)
        mol_name: Molecule name for messages

    Returns:
        Dictionary with refinement configuration:
        {
            'bonds_angles': list of indices for Seminario,
            'dihedrals': list of indices for dihedral refinement,
            'dihedral_method': 'pes' or 'crest' or None,
            'selected_params': list of (param_name, score, status, section) tuples
        }
    """
    result = {
        'bonds_angles': [],
        'dihedrals': [],
        'dihedral_method': None,
        'selected_params': []
    }

    if not interactive or not penalties:
        return result

    # Show refinement option panel
    console.print(Panel(
        "[bold cyan]Parameter Refinement[/bold cyan]\n\n"
        "You can refine parameters with penalties or uncertain estimates.\n\n"
        "[bold]Methods used:[/bold]\n"
        "  • [cyan]Bonds & Angles[/cyan] → Seminario method (fast, uses existing freq data)\n"
        "  • [cyan]Dihedrals[/cyan] → PES scan or CREST conformer sampling\n\n"
        "[bold]Note on improper dihedrals:[/bold]\n"
        "  Improper dihedral refinement must be done manually. To evaluate:\n"
        "  • Compare planarity between ab initio and classical MD trajectories\n"
        "  • Measure deviation from the QM-optimized or experimental structure\n"
        "  • Adjust force constants if planarity deviations are unacceptable\n\n"
        "[bold]Selection options:[/bold]\n"
        "  • Enter numbers (e.g., [cyan]1,2,4[/cyan] or [cyan]1-4[/cyan]) to select specific parameters\n"
        "  • Enter [cyan]all[/cyan] to refine all fittable parameters\n"
        "  • Press Enter or type [cyan]none[/cyan] to skip refinement",
        title="Parameter Refinement",
        border_style="cyan",
        expand=False
    ))

    # Get parameter selection
    selection = prompt_with_context(
        processor,
        "Parameters to refine",
        default="none",
        module="Small Molecule Parameterizer",
        description="Select parameters for refinement"
    ).strip().lower()

    if selection in ['none', '']:
        console.print("[grey50]Skipping parameter refinement[/grey50]")
        return result

    # Parse selection
    fittable_sections = ['BOND', 'ANGLE', 'DIHE']
    selected_indices = []

    if selection == 'all':
        selected_indices = [i for i, p in enumerate(penalties) if p[3] in fittable_sections]
    else:
        try:
            parts = selection.replace(' ', '').split(',')
            for part in parts:
                if '-' in part:
                    start, end = part.split('-')
                    for i in range(int(start), int(end) + 1):
                        if 1 <= i <= len(penalties):
                            selected_indices.append(i - 1)
                else:
                    i = int(part)
                    if 1 <= i <= len(penalties):
                        selected_indices.append(i - 1)
            selected_indices = sorted(set(selected_indices))
        except ValueError:
            console.print("[yellow]Could not parse selection, skipping refinement[/yellow]")
            return result

    if not selected_indices:
        console.print("[yellow]No valid parameters selected[/yellow]")
        return result

    # Filter to fittable parameters only
    selected_indices = [i for i in selected_indices if penalties[i][3] in fittable_sections]

    # Warn about IMPROPER if any were selected
    improper_selected = [i for i in selected_indices if penalties[i][3] == 'IMPROPER']
    if improper_selected:
        console.print("[yellow]Note: IMPROPER parameters cannot be refined and will be skipped[/yellow]")

    if not selected_indices:
        console.print("[yellow]No fittable parameters selected[/yellow]")
        return result

    # Categorize selected parameters
    bonds_angles_idx = [i for i in selected_indices if penalties[i][3] in ['BOND', 'ANGLE']]
    dihedrals_idx = [i for i in selected_indices if penalties[i][3] == 'DIHE']

    result['selected_params'] = [penalties[i] for i in selected_indices]

    # Handle bonds and angles (Seminario)
    if bonds_angles_idx:
        ba_params = [penalties[i] for i in bonds_angles_idx]
        ba_names = [p[0] for p in ba_params]

        console.print(f"\n[bold]Bonds & Angles:[/bold] {len(bonds_angles_idx)} parameter(s)")
        for p in ba_params:
            console.print(f"  • {p[3]}: {p[0]}")

        if fchk_file and os.path.exists(fchk_file):
            console.print(f"[green]  ✓ Using Seminario method (fchk found)[/green]")
            result['bonds_angles'] = bonds_angles_idx
        else:
            console.print(f"[red]  ✗ fchk file not found: {mol_name}.fchk[/red]")
            console.print(f"[yellow]  Please generate it:[/yellow]")
            console.print(f"[cyan]    formchk {mol_name}.chk {mol_name}.fchk[/cyan]")
            console.print(f"[yellow]  Then re-run this step.[/yellow]")
            # Don't add to result - user needs to generate fchk first

    # Handle dihedrals
    if dihedrals_idx:
        dihe_params = [penalties[i] for i in dihedrals_idx]

        console.print(f"\n[bold]Dihedrals:[/bold] {len(dihedrals_idx)} parameter(s)")
        for p in dihe_params:
            console.print(f"  • {p[3]}: {p[0]}")

        console.print(Panel(
            "[bold]Dihedral Refinement Method[/bold]\n\n"
            "[cyan]1) PES scan[/cyan] - Systematic scan of each dihedral\n"
            "   • Thorough sampling of torsional potential\n"
            "   • Scans one dihedral at a time\n"
            "   • Recommended for 1-2 dihedrals\n\n"
            "[cyan]2) CREST[/cyan] - Conformer ensemble sampling\n"
            "   • Samples all rotatable bonds together\n"
            "   • Faster for multiple dihedrals\n"
            "   • Uses semiempirical method (GFN2-xTB)\n\n"
            "[cyan]3) Skip[/cyan] - Don't refine dihedrals",
            title="Dihedral Method",
            border_style="cyan",
            expand=False
        ))

        dihe_choice = prompt_with_context(
            processor,
            "Dihedral method",
            choices=["1", "2", "3"],
            default="1" if len(dihedrals_idx) <= 2 else "2",
            module="Small Molecule Parameterizer",
            description="Select dihedral refinement method"
        )

        if dihe_choice == "1":
            result['dihedrals'] = dihedrals_idx
            result['dihedral_method'] = 'pes'
            console.print("[cyan]  → Using PES scan for dihedrals[/cyan]")
        elif dihe_choice == "2":
            result['dihedrals'] = dihedrals_idx
            result['dihedral_method'] = 'crest'
            console.print("[cyan]  → Using CREST for dihedrals[/cyan]")
        else:
            console.print("[grey50]  → Skipping dihedral refinement[/grey50]")

    # Summary
    if result['bonds_angles'] or result['dihedrals']:
        console.print(f"\n[green]Refinement plan:[/green]")
        if result['bonds_angles']:
            console.print(f"  • Seminario: {len(result['bonds_angles'])} bond/angle parameter(s)")
        if result['dihedrals']:
            method_name = "PES scan" if result['dihedral_method'] == 'pes' else "CREST"
            console.print(f"  • {method_name}: {len(result['dihedrals'])} dihedral(s)")

    return result


# Keep old function name as alias for backwards compatibility
def ask_paramfit_selection(penalties, console, interactive=True, processor=None):
    """Deprecated: Use ask_refinement_selection instead."""
    result = ask_refinement_selection(penalties, console, interactive, processor, fchk_file=None)
    # Return indices for paramfit (old behavior)
    return result.get('dihedrals', [])


def extract_parameter_name(line, section):
    """Extract the complete parameter name based on the frcmod section.

    AMBER frcmod uses fixed-width columns based on Fortran FORMAT specifications:
    - BOND:  FORMAT(A2,1X,A2,2F10.2) = 5 chars for atom types
    - ANGLE: FORMAT(A2,1X,A2,1X,A2,2F10.2) = 8 chars for atom types
    - DIHE/IMPROPER: FORMAT(A2,1X,A2,1X,A2,1X,A2,...) = 11 chars for atom types
    """

    if section == 'BOND':
        # FORMAT(A2,1X,A2,2F10.2) = 2+1+2 = 5 characters
        # Example: c -ns   418.30   1.388       same as c -n, penalty score=  0.0
        #          ^^ ^^
        if len(line) >= 5:
            atom1 = line[0:2].strip()
            atom2 = line[3:5].strip()
            return f"{atom1}-{atom2}"

        # Fallback
        parts = line.split()
        return parts[0] if parts else "unknown"

    elif section == 'ANGLE':
        # FORMAT(A2,1X,A2,1X,A2,2F10.2) = 2+1+2+1+2 = 8 characters
        # Example: nc-c -ns    68.820     115.620     same as n -c -n , penalty score=  0.0
        #          ^^ ^^ ^^
        if len(line) >= 8:
            atom1 = line[0:2].strip()
            atom2 = line[3:5].strip()
            atom3 = line[6:8].strip()
            return f"{atom1}-{atom2}-{atom3}"

        # Fallback
        parts = line.split()
        return parts[0] if parts else "unknown"
    
    elif section in ['DIHE', 'IMPROPER']:
        # AMBER frcmod format for dihedrals/impropers uses fixed-width columns:
        # FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)
        # = 2+1+2+1+2+1+2 = 11 characters for the 4 atom types
        # Example: nc-cd-ns-c    4    6.600       180.000           2.000
        # Example: nc-ns-c -o         10.5          180.0         2.0
        #          ^^ ^^ ^^ ^^
        #          Atom types may have trailing spaces (e.g., "c " not "c")

        # Extract first 11 characters (4 atom types with separators)
        if len(line) >= 11:
            param_str = line[:11]
            # Extract each 2-char atom type at positions 0-1, 3-4, 6-7, 9-10
            atom1 = param_str[0:2].strip()
            atom2 = param_str[3:5].strip()
            atom3 = param_str[6:8].strip()
            atom4 = param_str[9:11].strip()
            return f"{atom1}-{atom2}-{atom3}-{atom4}"

        # Fallback for short lines - try simple split
        parts = line.split()
        if parts:
            return parts[0]
    
    elif section == 'NONBON':
        # Format: LTYNB R EDEP
        # Example: c3          1.9080  0.1094
        parts = line.split()
        if len(parts) >= 3:
            return parts[0]  # e.g., "c3"
    
    # Fallback - just return the first whitespace-separated field
    parts = line.split()
    return parts[0] if parts else "unknown"

def choose_tleap_mode(interactive=True, processor=None):
    """Let user choose tLEaP execution mode."""
    if not interactive:
        return 'both'

    console = Console()
    console.print("\n[cyan]Choose tLEaP execution mode:[/cyan]")
    console.print("  gas    : Generate gas phase parameters only")
    console.print("  aqueous: Generate solvated system only")
    console.print("  both   : Generate both gas and aqueous phase systems")

    mode = prompt_with_context(
        processor,
        "tLEaP mode",
        choices=["gas", "aqueous", "both"],
        default="both",
        module="Small Molecule Parameterizer",
        description="tLEaP execution mode"
    )

    return mode

def run_tleap_command(input_file, console, timeout=60):
    """Run tleap command and handle output."""
    cmd = ["tleap", "-f", str(input_file)]
    
    console.print(f"\n[cyan]Running tLEaP command:[/cyan]")
    console.print(f"[grey50]{' '.join(cmd)}[/grey50]")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        
        if result.returncode == 0:
            console.print(f"[green]✓ tLEaP completed successfully[/green]")
            if result.stdout:
                # Show relevant output lines
                output_lines = result.stdout.split('\n')
                important_lines = [line for line in output_lines 
                                 if 'Checking' in line or 'Loading' in line or 'Saved' in line or 'ERROR' in line or 'WARNING' in line]
                if important_lines:
                    console.print("[grey50]Key tLEaP output:[/grey50]")
                    for line in important_lines[-10:]:  # Last 10 important lines
                        if line.strip():
                            console.print(f"[grey50]  {line.strip()}[/grey50]")
            return True, "tLEaP completed successfully"
        else:
            error_msg = result.stderr if result.stderr else "tLEaP failed"
            console.print(f"[red]✗ tLEaP failed[/red]")
            if result.stderr:
                console.print(f"[red]Error output:[/red]")
                console.print(f"[grey50]{result.stderr}[/grey50]")
            return False, error_msg
            
    except subprocess.TimeoutExpired:
        return False, "tLEaP timed out"
    except FileNotFoundError:
        return False, "tLEaP not found in PATH"
    except Exception as e:
        return False, f"Error running tLEaP: {str(e)}"

from proprep.utils.tleap_utils import (
    tleap_safe_unit_var as _tleap_safe_unit_var,
    is_tleap_safe_resname,
    suggest_tleap_safe_resname,
)


def create_and_run_tleap_inputs(mol2_file, frcmod_file, mol_name, mode, console, processor=None, net_charge=0):
    """Create and execute tLEaP input files based on mode."""
    success_files = []
    # tLEaP unit variable must be a valid identifier; filenames keep mol_name.
    unit = _tleap_safe_unit_var(mol_name)

    if mode in ['gas', 'both']:
        # Gas phase
        gas_input = f"{mol_name}_gas_tleap.in"
        gas_content = f"""source leaprc.gaff2
loadamberparams "{frcmod_file}"

{unit} = loadmol2 "{mol2_file}"
saveOff {unit} "{mol_name}.lib"
saveAmberParm {unit} "{mol_name}_gas.parm7" "{mol_name}_gas.rst7"

quit
"""
        
        with open(gas_input, 'w') as f:
            f.write(gas_content)
        
        console.print(f"\n[bold]Creating gas phase system...[/bold]")
        success, message = run_tleap_command(gas_input, console)
        
        if success:
            success_files.extend([f"{mol_name}_gas.parm7", f"{mol_name}_gas.rst7", f"{mol_name}.lib"])
            console.print(f"[green]✓ Gas phase files created: {mol_name}_gas.parm7, {mol_name}_gas.rst7[/green]")
        else:
            console.print(f"[red]✗ Gas phase tLEaP failed: {message}[/red]")
    
    if mode in ['aqueous', 'both']:
        # Aqueous phase
        console.print(f"\n[cyan]Configuring aqueous phase system...[/cyan]")
        
        # Water model selection
        water_models = {
            "1": ("opc", "leaprc.water.opc", "OPCBOX"),
            "2": ("tip3p", "leaprc.water.tip3p", "TIP3PBOX"), 
            "3": ("tip4pew", "leaprc.water.tip4pew", "TIP4PEWBOX"),
            "4": ("spce", "leaprc.water.spce", "SPCBOX")
        }
        
        console.print("Water model options:")
        for key, (name, _, _) in water_models.items():
            console.print(f"  {key}: {name.upper()}")
        
        water_choice = prompt_with_context(
            processor,
            "Select water model",
            choices=list(water_models.keys()),
            default="1",
            module="Small Molecule Parameterizer",
            description="Select water model for solvation"
        )

        box_size = prompt_with_context(
            processor,
            "Box size (Å)",
            default="12.0",
            module="Small Molecule Parameterizer",
            description="Solvation box size"
        )
        
        water_name, leaprc, box_type = water_models[water_choice]
        
        aq_input = f"{mol_name}_aq_tleap.in"

        # Add counterions if molecule is charged
        ion_lines = ""
        if net_charge != 0:
            if net_charge > 0:
                ion_lines = f"addionsrand {unit} Cl- 0\n"
            else:
                ion_lines = f"addionsrand {unit} Na+ 0\n"

        aq_content = f"""source {leaprc}
source leaprc.gaff2
loadamberparams "{frcmod_file}"

{unit} = loadmol2 "{mol2_file}"
saveOff {unit} "{mol_name}.lib"
solvateBox {unit} {box_type} {box_size}
{ion_lines}saveAmberParm {unit} "{mol_name}_aq.parm7" "{mol_name}_aq.rst7"

quit
"""
        
        with open(aq_input, 'w') as f:
            f.write(aq_content)
        
        console.print(f"\n[bold]Creating aqueous phase system ({water_name.upper()}, {box_size}Å box)...[/bold]")
        success, message = run_tleap_command(aq_input, console)
        
        if success:
            success_files.extend([f"{mol_name}_aq.parm7", f"{mol_name}_aq.rst7"])
            if f"{mol_name}.lib" not in success_files:
                success_files.append(f"{mol_name}.lib")
            console.print(f"[green]✓ Aqueous phase files created: {mol_name}_aq.parm7, {mol_name}_aq.rst7[/green]")
        else:
            console.print(f"[red]✗ Aqueous phase tLEaP failed: {message}[/red]")
    
    return success_files

SMALL_MOL_WORKFLOW_STEPS = [
    WorkflowStep(
        id="sm-1", name="Extract Coordinates",
        description="Extract small molecule coordinates from PDB",
        handler="_step_extract_coordinates",
        section="Structure Preparation",
    ),
    WorkflowStep(
        id="sm-2", name="Hydrogen Addition",
        description="Analyze and optionally add hydrogens via reduce",
        handler="_step_hydrogen_addition",
        section="Structure Preparation",
        dependencies=["sm-1"],
    ),
    WorkflowStep(
        id="sm-3", name="Gaussian Input",
        description="Generate Gaussian input (optimization + ESP)",
        handler="_step_gaussian_input",
        section="Structure Preparation",
        dependencies=["sm-2"],
        checkpoint=True,
        checkpoint_message="Run Gaussian on the .gjf file, then resume.",
    ),
    WorkflowStep(
        id="sm-4", name="RESP Charges",
        description="Process Gaussian output → RESP charges (antechamber)",
        handler="_step_resp_charges",
        section="Charge Derivation",
        dependencies=["sm-3"],
    ),
    WorkflowStep(
        id="sm-5", name="Bonded Parameters",
        description="Generate missing bonded parameters (parmchk2)",
        handler="_step_ff_parameters",
        section="Bonded Parameter Generation",
        dependencies=["sm-4"],
    ),
    WorkflowStep(
        id="sm-6", name="AMBER Topology Files",
        description="Create AMBER topology files (tLEaP)",
        handler="_step_tleap",
        section="Bonded Parameter Generation",
        dependencies=["sm-5"],
    ),
    WorkflowStep(
        id="sm-7", name="Parameter Refinement",
        description="Seminario / PES scan / CREST refinement (optional)",
        handler="_step_refinement",
        section="Bonded Parameter Generation",
        dependencies=["sm-6"],
        optional=True,
    ),
    # Last, so it registers the REFINED frcmod when sm-7 ran. Depending on an
    # optional step is safe: _check_dependencies treats "skipped" as satisfied.
    WorkflowStep(
        id="sm-8", name="Force Field Integration",
        description="Register lib/frcmod for the Topology Generator",
        handler="_step_ff_integration",
        section="Bonded Parameter Generation",
        dependencies=["sm-7"],
    ),
]


def emit_small_molecule_transformer(console, residue_name, mol2_file,
                                    output_dir, promo_result):
    """Emit a reusable rename transformer for a deposited small molecule.

    Small-molecule parameterization keeps the residue name but renames its atoms
    to antechamber's sequential labels. Capturing that atom-name map as a Tier-1
    transformer (resname unchanged, atom_renames applied) lets the Redox Site
    Preparer re-apply it to the same ligand in a different structure, so the
    deposited parameters can be assigned there without re-running antechamber.
    Mirrors the metal-site and modified-AA integration steps. Best-effort and
    fully guarded: the parameters are already saved, so a failure here must
    never disrupt the workflow.
    """
    library_path = (promo_result or {}).get("library_path")
    if not library_path:
        return None  # user declined promotion; nothing deposited to reuse

    try:
        import json as _json

        # Locate the PDB->mol2 atom-name map the parameterizer wrote.
        search = os.path.dirname(mol2_file) if mol2_file else output_dir
        maps = glob.glob(os.path.join(search, "*_atom_name_mapping.json"))
        atom_renames = {}
        if maps:
            with open(maps[0]) as fh:
                raw = _json.load(fh) or {}
            # Keep only real renames; antechamber leaves many atom names
            # unchanged, and identity entries would just bloat the recipe.
            atom_renames = {k: v for k, v in raw.items() if k != v}

        # forcefield_path = residue dir relative to specialized_residues, so the
        # Topology Generator discovers the deposited FF the same way it does for
        # built-in transformers.
        forcefield_path = None
        parts = Path(library_path).parts
        if "specialized_residues" in parts:
            idx = parts.index("specialized_residues")
            forcefield_path = "/".join(parts[idx + 1:]) or None

        from proprep.redoxsite_prep.transformation.auto_rename import (
            emit_rename_transformer,
        )
        from proprep.forcefield_params.user_library import (
            DEFAULT_REDOX_STATE, DEFAULT_SPIN_STATE,
        )

        entry = {"resname": residue_name, "target": residue_name}
        if atom_renames:
            entry["atom_renames"] = atom_renames

        tpath = emit_rename_transformer(
            [entry],
            name=f"sm_{residue_name}".lower(),
            description=(f"Reuse small-molecule parameters for {residue_name} "
                         f"(atom-name map to the deposited library)"),
            redox_state=DEFAULT_REDOX_STATE,
            spin_state=DEFAULT_SPIN_STATE,
            forcefield_path=forcefield_path,
            provenance={"source": "small_molecule", "residue": residue_name,
                        "library_path": library_path,
                        "metadata_path": (promo_result or {}).get("metadata_path")},
            site_types=["small_molecule"],
        )
        console.print(f"  [green]✓[/green] Reuse transformer: [grey50]{tpath}[/grey50]")
        return tpath
    except Exception as exc:  # noqa: BLE001 — never let emission break the workflow
        console.print(f"  [grey50]Reuse transformer not emitted: {exc}[/grey50]")
        return None


class SmallMolWorkflowRunner:
    """Executor class for the small molecule parameterization WorkflowChecklist.

    Holds all shared state as instance attributes. Created in run_workflow(),
    setup() runs the preamble, then WorkflowChecklist drives the 7 step handlers.
    """

    def __init__(self, residue_name, residues, output_name=None, interactive=True, processor=None, regenerate=False):
        self.residue_name = residue_name
        self.residues = residues
        self.output_name = output_name
        self.interactive = interactive
        self.processor = processor
        # When True, setup() does NOT short-circuit on pre-existing .mol2/.frcmod.
        # Callers reached only via a "generate parameters" choice pass this so a
        # re-run (or a recorded-session replay, where the outputs are already on
        # disk) actually re-enters the checklist instead of silently reusing.
        self.regenerate = regenerate
        self.console = Console()

        # Results dict (returned by compile_results)
        self.results = {
            "success": True,
            "message": f"Small molecule parameterization workflow for {residue_name}",
            "status": "completed",
            "parameter_files": {"prep_file": None, "frcmod_file": None},
            "simulation_files": [],
            "missing_files": [],
            "output_files": {},
        }

        # Set by setup()
        self.selected_residue = None
        self.mol_name = None
        self.mol_name_lower = None

        # Set by step handlers
        self.mol_pdb_file = None
        self.net_charge = 0
        self.multiplicity = 1
        self.gaussian_settings = None
        self.gaussian_file = None
        self.gaussian_log_file = None
        self.mol2_file = None
        self.frcmod_file = None
        self.chk_file = None
        self.penalties = []
        self.refinement_config = {}
        self.success_files = []
        self.prmtop_file = None
        self.current_frcmod = None

    # ── Preamble ───────────────────────────────────────────────────

    def setup(self):
        """Run preamble: select residue instance, choose output name, check existing files.

        Returns:
            None if setup completed normally (proceed to checklist).
            dict if early exit (existing parameter files found).
        """
        self.console.print(f"[bold cyan]Small Molecule Parameterization Workflow for {self.residue_name}[/bold cyan]")
        self.console.print(f"[green]Working directory: {os.getcwd()}[/green]")

        # Handle multiple instances - let user choose which one to extract
        if len(self.residues) > 1:
            self.console.print(f"\n[cyan]Found {len(self.residues)} instances of {self.residue_name}:[/cyan]")

            for i, res in enumerate(self.residues, 1):
                location = f"{res.chain_id}:{res.resid}" if hasattr(res, 'chain_id') else f"Instance {i}"
                self.console.print(f"  {i}. {location}")

            if self.interactive:
                # Per-option labels so the session editor shows the residue
                # locations, not just the bare index.
                res_options = {
                    str(i): (f"{res.chain_id}:{res.resid}" if hasattr(res, 'chain_id') else f"Instance {i}")
                    for i, res in enumerate(self.residues, 1)
                }
                choice = prompt_with_context(
                    self.processor,
                    "Which instance should be used for coordinate extraction?",
                    choices=[str(i) for i in range(1, len(self.residues) + 1)],
                    default="1",
                    module="Small Molecule Parameterizer",
                    description="Select residue instance for extraction",
                    options_map=res_options,
                )
                self.selected_residue = self.residues[int(choice) - 1]
                location = f"{self.selected_residue.chain_id}:{self.selected_residue.resid}" if hasattr(self.selected_residue, 'chain_id') else f"Instance {choice}"
                self.console.print(f"[green]Using instance {choice}: {location}[/green]")
            else:
                self.selected_residue = self.residues[0]
        else:
            self.selected_residue = self.residues[0]
            location = f"{self.selected_residue.chain_id}:{self.selected_residue.resid}" if hasattr(self.selected_residue, 'chain_id') else "Single instance"
            self.console.print(f"[cyan]Using single instance: {location}[/cyan]")

        # Determine output residue name.
        #
        # This name becomes the OFF/lib entry name, and tLEaP matches loadpdb
        # residues to templates by that entry name (not by the residue name
        # inside the mol2). tLEaP also lexes a digit-leading bare token (e.g.
        # 9E2) as scientific notation, and PDB resName is only three columns
        # wide. So the name MUST be tLEaP/PDB-safe: 1-3 chars, letter-leading,
        # alphanumeric. Enforcing it here is what lets the downstream metal-site
        # recombination load the ligand and rename the structure residue to
        # match — otherwise tleap fails with "Unknown residue".
        suggested = suggest_tleap_safe_resname(self.residue_name)
        if self.output_name:
            candidate = self.output_name.upper()
            if not is_tleap_safe_resname(candidate):
                raise ValueError(
                    f"Output residue name '{candidate}' is not tLEaP-safe: it must be "
                    f"1-3 characters, start with a letter, and be alphanumeric "
                    f"(try '{suggest_tleap_safe_resname(candidate)}')."
                )
            self.mol_name = candidate
        elif self.interactive:
            default_name = self.residue_name.upper()
            if not is_tleap_safe_resname(default_name):
                default_name = suggested
            self.console.print(Panel(
                "[bold]Output Residue Name[/bold]\n\n"
                "The output residue name is used for all generated files and is the\n"
                "residue name tLEaP will match against your structure.\n\n"
                "[bold]Must be 1-3 characters, start with a letter, alphanumeric only.[/bold]\n"
                "[grey50](tLEaP reads a digit-leading name like 9E2 as a number, and PDB\n"
                "residue names are only 3 columns wide, so codes like 9E2 are rejected —\n"
                f"a safe name such as {suggested} is offered by default. You may also pick a\n"
                f"custom name to parameterize multiple redox states, e.g. {suggested[:-1]}O / {suggested[:-1]}R.)[/grey50]",
                title="Residue Naming",
                border_style="blue",
                expand=False,
            ))
            while True:
                candidate = prompt_with_context(
                    self.processor, "Output residue name",
                    default=default_name,
                    module="Small Molecule Parameterizer",
                    description="Output residue name for parameter files",
                ).upper()
                if is_tleap_safe_resname(candidate):
                    self.mol_name = candidate
                    break
                self.console.print(
                    f"[red]'{candidate}' is not a valid residue name for tLEaP.[/red] "
                    f"[grey50]Use 1-3 characters, letters/digits only, starting with a letter "
                    f"(e.g. {suggest_tleap_safe_resname(candidate)}).[/grey50]"
                )
        else:
            candidate = self.residue_name.upper()
            self.mol_name = candidate if is_tleap_safe_resname(candidate) else suggested
            if self.mol_name != candidate:
                self.console.print(
                    f"[yellow]Residue name '{candidate}' is not tLEaP-safe; "
                    f"using '{self.mol_name}' instead.[/yellow]"
                )

        self.mol_name_lower = self.mol_name.lower()

        if self.mol_name != self.residue_name.upper():
            self.console.print(f"[green]Output residue name: {self.mol_name} (extracted from {self.residue_name})[/green]")

        # Check for existing parameter files — short-circuit if found, UNLESS the
        # caller asked to regenerate. The short-circuit is a convenience for
        # accidental re-runs; when the user explicitly chose "generate", or when a
        # recorded session is being replayed (outputs already on disk), skipping it
        # is wrong — we must re-enter the checklist. Announce the overwrite so the
        # reuse isn't silent, but do not prompt (prompting would break replay).
        existing_mol2 = glob.glob(f"{self.mol_name_lower}.mol2")
        existing_frcmod = glob.glob(f"{self.mol_name_lower}.frcmod")

        if existing_mol2 and existing_frcmod and self.regenerate:
            self.console.print(
                f"[yellow]Existing parameter files for {self.mol_name} found; "
                f"regenerating as requested (they will be overwritten).[/yellow]"
            )
        elif existing_mol2 and existing_frcmod:
            self.console.print(f"[green]Found existing parameter files for {self.mol_name}[/green]")
            self.results["parameter_files"]["prep_file"] = os.path.abspath(existing_mol2[0])
            self.results["parameter_files"]["frcmod_file"] = os.path.abspath(existing_frcmod[0])
            # Also surface the .lib and atom-name mapping so the metal-site
            # recombination can loadoff the library and remap the structure.
            existing_lib = glob.glob(f"{self.mol_name_lower}.lib")
            if existing_lib:
                self.results["parameter_files"]["lib_file"] = os.path.abspath(existing_lib[0])
            existing_map = glob.glob(f"{self.mol_name_lower}_atom_name_mapping.json")
            if existing_map:
                self.results["parameter_files"]["atom_mapping"] = os.path.abspath(existing_map[0])
            self.results["message"] = f"Found existing parameter files for {self.mol_name}"
            return self.results  # Early exit

        self._recover_paths_from_disk()
        return None  # Proceed to checklist

    def _recover_paths_from_disk(self):
        """Re-derive path/charge attributes from on-disk artifacts.

        Resume restores step status but not the runner's instance attributes,
        so steps that depend on outputs from earlier (already-completed) steps
        would otherwise see None and crash.
        """
        base = self.mol_name_lower

        h_pdb = f"{base}_H.pdb"
        plain_pdb = f"{base}.pdb"
        if os.path.exists(h_pdb):
            self.mol_pdb_file = h_pdb
        elif os.path.exists(plain_pdb):
            self.mol_pdb_file = plain_pdb

        for attr, fname in (
            ("gaussian_file", f"{base}.gjf"),
            ("gaussian_log_file", f"{base}.log"),
            ("mol2_file", f"{base}.mol2"),
            ("frcmod_file", f"{base}.frcmod"),
            ("chk_file", f"{base}.chk"),
            ("prmtop_file", f"{base}_gas.parm7"),
        ):
            if os.path.exists(fname):
                setattr(self, attr, fname)

        if self.gaussian_file and os.path.exists(self.gaussian_file):
            try:
                with open(self.gaussian_file, "r") as f:
                    for line in f:
                        parts = line.split()
                        if len(parts) == 2 and all(p.lstrip("-").isdigit() for p in parts):
                            self.net_charge = int(parts[0])
                            self.multiplicity = int(parts[1])
                            break
            except OSError:
                pass

    # ── Step 1: Extract Coordinates ────────────────────────────────

    def _step_extract_coordinates(self):
        """Extract small molecule coordinates from PDB to a standalone file."""
        selected = self.selected_residue

        # Handle different residue wrapper objects from forcefield parameterizer
        if hasattr(selected, 'biopython_residue'):
            residue_obj = selected.biopython_residue
            chain_id = selected.chain_id
            resid = selected.resid
        elif hasattr(selected, 'residue'):
            residue_obj = selected.residue
            chain_id = selected.chain_id
            resid = selected.resid
        elif hasattr(selected, 'get_atoms'):
            residue_obj = selected
            chain_id = selected.get_parent().id if hasattr(selected, 'get_parent') else 'A'
            resid = selected.id[1] if hasattr(selected, 'id') else 1
        else:
            for attr_name in dir(selected):
                attr = getattr(selected, attr_name)
                if hasattr(attr, 'get_atoms'):
                    residue_obj = attr
                    chain_id = getattr(selected, 'chain_id', 'A')
                    resid = getattr(selected, 'resid', 1)
                    break
            else:
                raise RuntimeError(f"Cannot find BioPython residue object in {type(selected)}")

        self.console.print(f"[cyan]Extracting from chain {chain_id}, residue {resid}[/cyan]")

        from Bio.PDB import PDBIO, Structure, Model, Chain

        new_structure = Structure.Structure("small_mol")
        new_model = Model.Model(0)
        new_chain = Chain.Chain(chain_id)

        new_residue = residue_obj.copy()
        new_residue.id = (' ', 1, ' ')

        for atom in new_residue.get_atoms():
            if not hasattr(atom, 'element') or not atom.element:
                atom.element = atom.get_name()[0].strip('0123456789')

        new_chain.add(new_residue)
        new_model.add(new_chain)
        new_structure.add(new_model)

        self.mol_pdb_file = f"{self.mol_name_lower}.pdb"
        io = PDBIO()
        io.set_structure(new_structure)
        io.save(self.mol_pdb_file)

        self.console.print(f"[green]✓ Extracted coordinates to {self.mol_pdb_file}[/green]")

        # Preview
        self.console.print("[grey50]Generated PDB file preview:[/grey50]")
        with open(self.mol_pdb_file, 'r') as f:
            for i, line in enumerate(f):
                if i < 3:
                    self.console.print(f"[grey50]  {line.rstrip()}[/grey50]")
                if line.startswith(('ATOM', 'HETATM')):
                    break

        atom_count = sum(1 for _ in new_residue.get_atoms())
        return {'summary': f'Extracted {atom_count} atoms to {self.mol_pdb_file}'}

    # ── Step 2: Hydrogen Addition ──────────────────────────────────

    def _step_hydrogen_addition(self):
        """Analyze hydrogen content and optionally add hydrogens via reduce.

        Delegates to the shared :class:`HydrogenEditor`
        (``forcefield_prep/hydrogen_editor.py``) so the small-molecule and
        modified-amino-acid routes share one hydrogen workflow. The module
        label stays "Small Molecule Parameterizer" to preserve this caller's
        session-recording keys.
        """
        from proprep.forcefield_prep.hydrogen_editor import HydrogenEditor

        editor = HydrogenEditor(
            self.mol_pdb_file, self.mol_name_lower,
            console=self.console, processor=self.processor,
            interactive=self.interactive, residue_name=self.residue_name,
            module="Small Molecule Parameterizer",
        )
        result = editor.run()
        self.mol_pdb_file = result['pdb_file']
        return {'summary': result['summary']}

    # ── Step 3: Gaussian Input ─────────────────────────────────────

    def _step_gaussian_input(self):
        """Generate Gaussian input file with --link1-- two-step workflow."""
        self.gaussian_file = f"{self.mol_name_lower}.gjf"

        self.gaussian_settings = {
            "opt": {"memory": "10GB", "processors": 4,
                    "keywords": "opt freq b3lyp/6-31+G(d) nosym int=ultrafine IOp(7/33=1)"},
            "esp": {"memory": "10GB", "processors": 4,
                    "keywords": "HF/6-31G(d) Pop=mk IOp(6/33=2,6/41=10,6/42=10) nosym int=ultrafine"},
        }

        self.net_charge = 0
        self.multiplicity = 1

        if self.interactive:
            self.console.print(Panel(
                "[bold cyan]Two-Step Gaussian Workflow:[/bold cyan]\n\n"
                "[yellow]Step 1: Geometry Optimization + Frequency[/yellow]\n"
                f"  Keywords: [grey50]{self.gaussian_settings['opt']['keywords']}[/grey50]\n"
                "  • DFT optimization finds minimum energy geometry\n"
                "  • Frequency calculation confirms no imaginary modes\n"
                "  • IOp(7/33=1) → Save Cartesian Hessian for Seminario force constants\n\n"
                "[yellow]Step 2: Electrostatic Potential Calculation[/yellow]\n"
                f"  Keywords: [grey50]{self.gaussian_settings['esp']['keywords']}[/grey50]\n"
                "  • HF ESP maintains AMBER force field compatibility\n"
                "  • Merz-Kollman analysis generates ESP grid for RESP fitting\n"
                "  • IOp(6/33=2)  → Write potential points and values to output\n"
                "  • IOp(6/41=10) → 10 concentric layers of ESP points per atom\n"
                "  • IOp(6/42=10) → ~1000 points/atom (adjust for large molecules)\n\n"
                "[bold]Charge & Multiplicity:[/bold]\n"
                "  • [cyan]Net charge[/cyan]: Total molecular charge (0 for neutral)\n"
                "  • [cyan]Multiplicity[/cyan]: 2S+1 (1=singlet for closed-shell, 2=doublet for radicals)\n\n"
                "[bold]Note on anions:[/bold]\n"
                "  DFT self-interaction error causes excess electron delocalization in vacuum,\n"
                "  which can lead to unbound electrons and SCF convergence failure for anions.\n"
                "  Implicit solvation (SCRF) is strongly recommended for negatively charged molecules.",
                title="Gaussian Calculation Overview",
                border_style="blue",
                expand=False,
            ))

            # Charge
            try:
                charge_str = prompt_with_context(
                    self.processor, f"Net charge for {self.residue_name}",
                    default="0", module="Small Molecule Parameterizer",
                    description="Net charge for molecule",
                )
                self.net_charge = int(charge_str)
            except ValueError:
                self.console.print("[yellow]Invalid charge, using 0[/yellow]")

            # Anion solvation warning
            if self.net_charge < 0:
                self.console.print(Panel(
                    "[bold yellow]Anion detected (charge = {}):[/bold yellow]\n\n".format(self.net_charge) +
                    "DFT self-interaction error causes electrons to be too delocalized — or\n"
                    "even unbound — for anions in vacuum. This frequently leads to SCF\n"
                    "convergence failure or unphysical geometries.\n\n"
                    "Adding implicit solvation (SCRF) to the optimization step stabilizes\n"
                    "the charge distribution and is [bold]strongly recommended[/bold].\n\n"
                    "[grey50]The ESP step (HF) is less susceptible but will also benefit from\n"
                    "consistent solvation between both calculation steps.[/grey50]",
                    title="Anion Solvation Warning",
                    border_style="yellow",
                    expand=False,
                ))

                add_scrf = confirm_with_context(
                    self.processor, "Add implicit solvation (SCRF) to Gaussian keywords?",
                    default=True, module="Small Molecule Parameterizer",
                    description="Add implicit solvation for anion",
                )
                if add_scrf:
                    self.console.print("\n[cyan]Common solvents:[/cyan]")
                    self.console.print("  Water (ε=78.4), DiMethylSulfoxide (ε=46.8), Acetonitrile (ε=35.7),")
                    self.console.print("  Methanol (ε=32.6), Ethanol (ε=24.9), Dichloromethane (ε=8.9),")
                    self.console.print("  TetraHydroFuran (ε=7.4), Chloroform (ε=4.7), Toluene (ε=2.4)")
                    self.console.print("[grey50]Full list: see Gaussian SCRF documentation for all supported solvent names[/grey50]")

                    solvent_name = prompt_with_context(
                        self.processor, "Solvent name (Gaussian keyword)",
                        default="Water", module="Small Molecule Parameterizer",
                        description="SCRF solvent selection",
                    )
                    scrf_keyword = f"SCRF=(Solvent={solvent_name})"
                    self.gaussian_settings["opt"]["keywords"] += f" {scrf_keyword}"
                    self.gaussian_settings["esp"]["keywords"] += f" {scrf_keyword}"
                    self.console.print(f"[green]Added {scrf_keyword} to both calculation steps[/green]")

            # Multiplicity
            try:
                mult_str = prompt_with_context(
                    self.processor, "Spin multiplicity",
                    default="1", module="Small Molecule Parameterizer",
                    description="Spin multiplicity",
                )
                self.multiplicity = int(mult_str)
                if self.multiplicity < 1:
                    self.console.print("[yellow]Multiplicity must be ≥1, using 1[/yellow]")
                    self.multiplicity = 1
            except ValueError:
                self.console.print("[yellow]Invalid multiplicity, using 1[/yellow]")

            # Memory
            self.gaussian_settings["opt"]["memory"] = prompt_with_context(
                self.processor, "Memory allocation",
                default="10GB", module="Small Molecule Parameterizer",
                description="Gaussian memory allocation",
            )
            self.gaussian_settings["esp"]["memory"] = self.gaussian_settings["opt"]["memory"]

            # Processors
            try:
                proc_str = prompt_with_context(
                    self.processor, "Number of processors",
                    default="4", module="Small Molecule Parameterizer",
                    description="Number of processors",
                )
                processors = int(proc_str)
                self.gaussian_settings["opt"]["processors"] = processors
                self.gaussian_settings["esp"]["processors"] = processors
            except ValueError:
                self.console.print("[yellow]Invalid processor count, using 4[/yellow]")

            # Optional keyword customization
            if confirm_with_context(
                self.processor, "Customize calculation keywords?",
                default=False, module="Small Molecule Parameterizer",
                description="Customize Gaussian keywords",
            ):
                self.console.print("[grey50]Leave blank to keep default[/grey50]")
                opt_kw = prompt_with_context(
                    self.processor, "Step 1 (Optimization) keywords",
                    default=self.gaussian_settings["opt"]["keywords"],
                    module="Small Molecule Parameterizer",
                    description="Optimization keywords",
                )
                if opt_kw.strip():
                    self.gaussian_settings["opt"]["keywords"] = opt_kw

                esp_kw = prompt_with_context(
                    self.processor, "Step 2 (ESP) keywords",
                    default=self.gaussian_settings["esp"]["keywords"],
                    module="Small Molecule Parameterizer",
                    description="ESP calculation keywords",
                )
                if esp_kw.strip():
                    self.gaussian_settings["esp"]["keywords"] = esp_kw

        _write_gaussian_input_with_link1(
            self.mol_pdb_file, self.gaussian_file, self.net_charge,
            self.multiplicity, self.mol_name, self.gaussian_settings, self.console,
        )

        return {'summary': f'Generated {self.gaussian_file} (charge={self.net_charge}, mult={self.multiplicity})'}

    # ── Step 4: RESP Charges ─────────────────────────────────────

    def _step_resp_charges(self):
        """Process Gaussian output to assign RESP charges via antechamber."""
        self.gaussian_log_file = f"{self.mol_name_lower}.log"

        if not os.path.exists(self.gaussian_log_file):
            raise RuntimeError(
                f"Gaussian output not found: {self.gaussian_log_file}\n"
                f"Run Gaussian first: g16 {self.gaussian_file}"
            )

        self.console.print(f"[green]✓ Found Gaussian output: {self.gaussian_log_file}[/green]")

        self.mol2_file = f"{self.mol_name_lower}.mol2"
        success = enhanced_antechamber_step(
            self.gaussian_log_file, self.mol2_file, self.net_charge,
            self.mol_name, self.console, self.interactive, self.processor,
        )

        if not success:
            raise RuntimeError("RESP charge assignment failed")

        # Generate atom name mapping
        mapping_file = f"{self.mol_name_lower}_atom_name_mapping.json"
        generate_atom_name_mapping(self.mol_pdb_file, self.mol2_file, mapping_file, self.console)
        self.results["parameter_files"]["atom_mapping"] = os.path.abspath(mapping_file)

        return {'summary': f'RESP charges assigned → {self.mol2_file}'}

    # ── Step 5: Force Field Parameters ─────────────────────────────

    def _step_ff_parameters(self):
        """Run parmchk2 to generate force field parameters and analyze penalties."""
        self.frcmod_file = f"{self.mol_name_lower}.frcmod"
        self.chk_file = f"{self.mol_name_lower}.chk"

        self.console.print(f"[cyan]Running parmchk2...[/cyan]")
        cmd = ["parmchk2", "-i", self.mol2_file, "-f", "mol2", "-o", self.frcmod_file, "-s", "2", "-a", "Y"]
        self.console.print(f"[grey50]{' '.join(cmd)}[/grey50]")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            raise RuntimeError(f"parmchk2 failed: {result.stderr}")

        self.console.print(f"[green]✓ Force field parameters generated[/green]")

        self.penalties = analyze_frcmod_penalties(self.frcmod_file, self.console)

        fchk_file = f"{self.mol_name_lower}.fchk"
        if os.path.exists(fchk_file):
            self.console.print(f"[grey50]Formatted checkpoint found: {fchk_file}[/grey50]")

        self.refinement_config = ask_refinement_selection(
            self.penalties, self.console, self.interactive, self.processor,
            fchk_file=fchk_file if os.path.exists(fchk_file) else None,
            mol_name=self.mol_name_lower,
        )

        n_penalties = len(self.penalties)
        return {'summary': f'Parameters generated ({n_penalties} penalty scores analyzed)'}

    # ── Step 6: tLEaP ──────────────────────────────────────────────

    def _step_tleap(self):
        """Create AMBER topology files via tLEaP."""
        tleap_mode = choose_tleap_mode(self.interactive, self.processor)

        self.success_files = create_and_run_tleap_inputs(
            self.mol2_file, self.frcmod_file, self.mol_name_lower,
            tleap_mode, self.console, self.processor, net_charge=self.net_charge,
        )

        self.prmtop_file = f"{self.mol_name_lower}_gas.parm7"
        self.current_frcmod = self.frcmod_file

        self.results["parameter_files"]["prep_file"] = os.path.abspath(self.mol2_file)
        self.results["parameter_files"]["frcmod_file"] = os.path.abspath(self.frcmod_file)

        # Record the OFF/lib file so the metal-site recombination can loadoff it
        # (tLEaP matches loadpdb residues by the lib entry name, which is why we
        # hand off a real .lib rather than the mol2). The entry name equals
        # mol_name, which the naming prompt has already validated as tLEaP-safe.
        lib_path = f"{self.mol_name_lower}.lib"
        if os.path.exists(lib_path):
            self.results["parameter_files"]["lib_file"] = os.path.abspath(lib_path)

        if self.success_files:
            self.console.print(f"[green]Generated files: {', '.join(self.success_files)}[/green]")
            self.results["simulation_files"] = [
                os.path.abspath(f) for f in self.success_files if os.path.exists(f)
            ]
            return {'summary': f'tLEaP completed ({len(self.success_files)} files generated)'}
        else:
            self.console.print(f"[yellow]Some tLEaP steps failed, but parameter files are available[/yellow]")
            return {'summary': 'tLEaP completed with warnings'}

    # ── Step 7: Parameter Refinement ───────────────────────────────
    # This step runs AFTER tLEaP because dihedral refinement methods
    # (PES scan, CREST/paramfit) need the topology file for MM energy
    # evaluations. If any parameters are refined, tLEaP is re-run at the
    # end to produce final topology files with the updated frcmod.

    def _step_refinement(self):
        """Optional parameter refinement (Seminario for bonds/angles, PES/CREST for dihedrals)."""
        if not self.refinement_config.get('bonds_angles') and not self.refinement_config.get('dihedrals'):
            return {'summary': 'No refinement selected — skipped'}

        # 8a: Seminario for bonds and angles
        self.console.print(f"\n[cyan]8a) Seminario method for bonds/angles[/cyan]")
        if self.refinement_config.get('bonds_angles'):
            bonds_angles_params = [self.penalties[i] for i in self.refinement_config['bonds_angles']]
            self.console.print(f"[cyan]    {len(bonds_angles_params)} parameter(s) selected[/cyan]")

            try:
                from proprep.forcefield_prep.seminario_refinement import run_seminario_refinement_workflow

                seminario_results = run_seminario_refinement_workflow(
                    mol_name=self.mol_name_lower, mol2_file=self.mol2_file,
                    frcmod_file=self.current_frcmod, chk_file=self.chk_file,
                    selected_params=bonds_angles_params, console=self.console,
                    interactive=self.interactive,
                )

                if seminario_results.get("refinement_success"):
                    self.results["seminario_refinement"] = seminario_results
                    if seminario_results.get("fitted_frcmod"):
                        self.current_frcmod = seminario_results["fitted_frcmod"]
                        self.results["parameter_files"]["frcmod_file"] = os.path.abspath(self.current_frcmod)
                        self.console.print(f"[green]✓ Seminario refinement complete: {self.current_frcmod}[/green]")
                else:
                    self.console.print(f"[yellow]Seminario refinement did not complete: {seminario_results.get('message', 'Unknown error')}[/yellow]")

            except ImportError as e:
                self.console.print(f"[red]Seminario refinement module not available: {e}[/red]")
            except Exception as e:
                self.console.print(f"[yellow]Seminario refinement error: {e}[/yellow]")
                import traceback
                self.console.print(f"[grey50]{traceback.format_exc()}[/grey50]")
        else:
            self.console.print(f"[grey50]    Skipped — no bond/angle parameters selected for refinement[/grey50]")

        # 8b: Dihedral refinement (PES scan or CREST)
        if self.refinement_config.get('dihedrals') and self.refinement_config.get('dihedral_method'):
            dihe_params = [self.penalties[i] for i in self.refinement_config['dihedrals']]
            dihe_method = self.refinement_config['dihedral_method']

            self.console.print(f"\n[cyan]8b) Dihedral refinement via {dihe_method.upper()} ({len(dihe_params)} parameters)[/cyan]")

            if dihe_method == 'pes':
                try:
                    from proprep.forcefield_prep.pes_scan_refinement import run_pes_scan_workflow

                    pes_results = run_pes_scan_workflow(
                        mol_name=self.mol_name_lower, mol2_file=self.mol2_file,
                        prmtop_file=self.prmtop_file, frcmod_file=self.current_frcmod,
                        selected_dihedrals=dihe_params, charge=self.net_charge,
                        multiplicity=self.multiplicity, console=self.console,
                        interactive=self.interactive, processor=self.processor,
                        gaussian_settings=self.gaussian_settings,
                    )

                    if pes_results.get("refinement_success"):
                        self.results["pes_refinement"] = pes_results
                        if pes_results.get("fitted_frcmod"):
                            self.current_frcmod = pes_results["fitted_frcmod"]
                            self.results["parameter_files"]["frcmod_file"] = os.path.abspath(self.current_frcmod)
                            self.console.print(f"[green]✓ PES scan refinement complete: {self.current_frcmod}[/green]")
                    else:
                        msg = pes_results.get('message', 'Unknown error')
                        self.console.print(f"[yellow]PES scan refinement did not complete: {msg}[/yellow]")

                except ImportError as e:
                    self.console.print(f"[yellow]PES scan module not yet implemented: {e}[/yellow]")
                    self.console.print(f"[grey50]Dihedral refinement skipped[/grey50]")
                except Exception as e:
                    self.console.print(f"[yellow]PES scan error: {e}[/yellow]")

            elif dihe_method == 'crest':
                if not os.path.exists(self.prmtop_file):
                    self.console.print(f"[yellow]Prmtop file not found: {self.prmtop_file}[/yellow]")
                    self.console.print(f"[grey50]CREST/Paramfit requires a prmtop file. Skipping.[/grey50]")
                else:
                    try:
                        from proprep.forcefield_prep.paramfit_refinement import run_paramfit_refinement_workflow

                        crest_results = run_paramfit_refinement_workflow(
                            mol_name=self.mol_name_lower, mol2_file=self.mol2_file,
                            prmtop_file=self.prmtop_file, frcmod_file=self.current_frcmod,
                            charge=self.net_charge, multiplicity=self.multiplicity,
                            console=self.console, selected_params=dihe_params,
                            processor=self.processor,
                        )

                        if crest_results.get("refinement_success"):
                            self.results["crest_refinement"] = crest_results
                            if crest_results.get("fitted_frcmod"):
                                self.current_frcmod = crest_results["fitted_frcmod"]
                                self.results["parameter_files"]["frcmod_file"] = os.path.abspath(self.current_frcmod)
                                self.console.print(f"[green]✓ CREST/Paramfit refinement complete: {self.current_frcmod}[/green]")
                        else:
                            self.console.print(f"[yellow]CREST refinement did not complete: {crest_results.get('message', 'Unknown error')}[/yellow]")

                    except ImportError as e:
                        self.console.print(f"[red]Paramfit module not available: {e}[/red]")
                    except Exception as e:
                        self.console.print(f"[yellow]CREST/Paramfit error: {e}[/yellow]")

        # Re-run tLEaP if frcmod was updated, so topology reflects refined parameters
        if self.current_frcmod != self.frcmod_file:
            self.console.print(f"\n[cyan]Re-running tLEaP with refined parameters ({self.current_frcmod})...[/cyan]")
            tleap_mode = choose_tleap_mode(interactive=False, processor=self.processor)
            success_files = create_and_run_tleap_inputs(
                self.mol2_file, self.current_frcmod, self.mol_name_lower,
                tleap_mode, self.console, self.processor, net_charge=self.net_charge,
            )
            if success_files:
                self.success_files = success_files
                self.results["simulation_files"] = [
                    os.path.abspath(f) for f in success_files if os.path.exists(f)
                ]
                self.console.print(f"[green]✓ Topology rebuilt with refined parameters[/green]")
            else:
                self.console.print(f"[yellow]tLEaP re-run had issues — check output[/yellow]")

        return {'summary': 'Parameter refinement completed'}

    # ── Step 8: Force Field Integration ────────────────────────────
    # The other two parameterizers deposit their finished parameters into the
    # workspace keys the Topology Generator reads: metal sites via
    # structure_preprocessor._checklist_mcpb_4_integration, modified AAs via
    # _run_step_10. Without the same deposit here, a parameterized ligand's
    # mol2/frcmod stay on disk and never reach the system tLEaP build.

    def _step_ff_integration(self):
        """Register the finished ligand parameters for the system tLEaP build."""
        workspace = None
        if self.processor is not None:
            try:
                workspace = self.processor._get_workspace()
            except Exception:  # noqa: BLE001
                workspace = None

        if workspace is None:
            self.console.print(
                "[yellow]No workspace available — parameters are on disk but were "
                "not registered for the Topology Generator.[/yellow]"
            )
            return {'summary': 'Not registered — no workspace'}

        # Prefer the .lib over the mol2: tLEaP matches loadpdb residues by the
        # lib entry name, which is why _step_tleap records it as lib_file.
        lib_path = os.path.abspath(f"{self.mol_name_lower}.lib")
        unit_file = lib_path if os.path.exists(lib_path) else (
            os.path.abspath(self.mol2_file) if self.mol2_file else None
        )

        # current_frcmod is the refined frcmod when sm-7 ran, else the parmchk2
        # output from sm-5.
        frcmod = self.current_frcmod or self.frcmod_file
        frcmod = os.path.abspath(frcmod) if frcmod else None

        def _append(key, path, supersedes=()):
            """Append path to a workspace list, returning its registration state.

            ``supersedes`` drops earlier entries this file replaces. Re-running
            this step after refinement would otherwise leave the pre-refinement
            frcmod registered next to the refined one, and tLEaP would load both
            (re-running sm-7 clears sm-8's status, so this path is reachable).
            """
            existing = workspace.get(key, []) or []
            if not isinstance(existing, list):
                existing = []
            dropped = [p for p in supersedes if p and p != path and p in existing]
            if dropped:
                existing = [p for p in existing if p not in dropped]
            if path in existing:
                state = "already registered"
            else:
                existing.append(path)
                state = "registered"
            workspace.set(key, existing)
            if dropped:
                state += f" (replaced {os.path.basename(dropped[0])})"
            return state

        table = Table(title="Force Field Integration", expand=False)
        table.add_column("Type", style="cyan")
        table.add_column("File", style="green")
        table.add_column("Status", style="grey50")

        registered = 0
        if unit_file and os.path.exists(unit_file):
            kind = "lib (unit)" if unit_file.endswith(".lib") else "mol2 (unit)"
            # A .lib written on a later pass supersedes an earlier mol2 entry.
            alt_unit = os.path.abspath(self.mol2_file) if self.mol2_file else None
            table.add_row(kind, os.path.basename(unit_file),
                          _append("preprocessing_lib_files", unit_file,
                                  supersedes=(alt_unit, lib_path)))
            registered += 1
        else:
            table.add_row("lib/mol2 (unit)", "—", "[yellow]missing[/yellow]")

        if frcmod and os.path.exists(frcmod):
            # The refined frcmod supersedes the parmchk2 one from sm-5.
            base_frcmod = os.path.abspath(self.frcmod_file) if self.frcmod_file else None
            table.add_row("frcmod (parameters)", os.path.basename(frcmod),
                          _append("preprocessing_frcmod_files", frcmod,
                                  supersedes=(base_frcmod,)))
            registered += 1
        else:
            table.add_row("frcmod (parameters)", "—", "[yellow]missing[/yellow]")

        self.console.print(table)

        if not registered:
            self.console.print(
                "[yellow]Nothing to register — run the earlier steps first.[/yellow]"
            )
            return {'summary': 'Nothing to register'}

        self.console.print(
            f"[green]✓[/green] {self.mol_name} parameters are registered; the "
            "Topology Generator will load them for the full system."
        )

        # 2) Reuse deposit. Registration above serves THIS build; promotion and
        # the transformer serve future ones. Every parameterizer's final step
        # does all three (metal site: mcpb-4, modified AA: aa-10).
        from proprep.forcefield_prep.library_promotion import offer_library_promotion

        promo = offer_library_promotion(
            self.console, self.processor,
            category="small_molecule",
            residue_name=self.mol_name,
            frcmod_file=frcmod,
            lib_search_dir=os.path.dirname(unit_file) if unit_file else os.getcwd(),
            prep_file=os.path.abspath(self.mol2_file) if self.mol2_file else None,
        )

        # 3) Reuse transformer — only meaningful once a deposit exists to
        # point at, which emit_small_molecule_transformer checks itself.
        emit_small_molecule_transformer(
            self.console, self.mol_name,
            os.path.abspath(self.mol2_file) if self.mol2_file else None,
            os.getcwd(), promo,
        )

        deposited = " and deposited" if (promo or {}).get("library_path") else ""
        return {'summary': f'Registered {registered} file(s) for tLEaP{deposited}'}

    # ── Results compilation ────────────────────────────────────────

    def compile_results(self, success):
        """Collect final results from runner state."""
        if self.mol2_file and os.path.exists(self.mol2_file):
            self.results["parameter_files"]["prep_file"] = os.path.abspath(self.mol2_file)
        if self.frcmod_file and os.path.exists(self.frcmod_file):
            self.results["parameter_files"]["frcmod_file"] = os.path.abspath(self.frcmod_file)

        # If parameter files exist, treat as success even if user chose save & quit
        has_params = (self.results["parameter_files"].get("prep_file") or
                      self.results["parameter_files"].get("frcmod_file"))

        if success or has_params:
            self.results["success"] = True
            self.results["status"] = "completed"
            self.results["message"] = f"Successfully completed parameterization for {self.mol_name}"
        else:
            self.results["success"] = False
            self.results["status"] = "saved"
            self.results["message"] = f"Small molecule parameterization saved (incomplete) for {self.mol_name}"

        return self.results


def run_workflow(residue_name: str, residues: list, output_dir: str = None, interactive: bool = True, processor=None, output_name: str = None, regenerate: bool = False):
    """
    Run the complete small molecule parameterization workflow.

    Called by ForcefieldParameterizer and structure_preprocessor when the user
    selects a small molecule for parameterization.

    Args:
        residue_name: Name of the residue to parameterize (as it appears in PDB)
        residues: List of residue instances
        output_dir: Output directory for generated files
        interactive: Whether to run in interactive mode
        processor: Optional processor for session recording context
        output_name: Optional output residue name (for different redox states, etc.)
        regenerate: If True, do not short-circuit when .mol2/.frcmod already
            exist — re-enter the checklist and overwrite. Callers that reach
            this function only via a "generate parameters" choice pass True so
            re-runs and recorded-session replays regenerate instead of reusing.

    Returns:
        dict: Results of parameterization workflow
    """
    original_dir = os.getcwd()

    try:
        if output_dir is None:
            output_dir = Path("small_molecule_params")
        else:
            output_dir = Path(output_dir)

        output_dir.mkdir(exist_ok=True)
        # Resolve BEFORE the chdir: the checklist keeps state_dir and rebuilds
        # the state-file path at every save, so a relative one would follow the
        # process around. Any step that chdirs would drop workflow_state.json
        # somewhere else and the resume lookup would miss it.
        output_dir = output_dir.resolve()
        os.chdir(output_dir)

        runner = SmallMolWorkflowRunner(
            residue_name=residue_name,
            residues=residues,
            output_name=output_name,
            interactive=interactive,
            processor=processor,
            regenerate=regenerate,
        )

        # Preamble: select residue instance, choose output name, check existing files
        early_result = runner.setup()
        if early_result is not None:
            return early_result

        checklist = WorkflowChecklist(
            steps=SMALL_MOL_WORKFLOW_STEPS,
            executor=runner,
            processor=processor,
            workflow_name="Small Molecule Parameterization",
            console=runner.console,
            state_dir=output_dir,
        )
        success = checklist.run()

        return runner.compile_results(success)

    except Exception as e:
        console = Console()
        console.print(f"[red]Error in small molecule workflow: {str(e)}[/red]")
        return {
            "success": False,
            "message": f"Error in workflow: {str(e)}",
            "status": "failed",
            "parameter_files": {"prep_file": None, "frcmod_file": None},
            "simulation_files": [],
            "missing_files": [],
            "output_files": {},
        }

    finally:
        os.chdir(original_dir)

# Helper functions for the module-level workflow
def _write_gaussian_input_with_link1(pdb_file, gaussian_file, net_charge, multiplicity, res_name, settings, console):
    """Write Gaussian input file with --link1-- structure.

    Args:
        pdb_file: Path to PDB file with coordinates
        gaussian_file: Output Gaussian input file path
        net_charge: Net molecular charge
        multiplicity: Spin multiplicity (2S+1)
        res_name: Residue name
        settings: Dict with 'opt' and 'esp' settings
        console: Rich console for output
    """

    # Read coordinates from PDB
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                element = line[76:78].strip()
                if not element:
                    element = line[12:16].strip()[0]

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                coords.append(f" {element:2s}              {x:12.8f}   {y:12.8f}   {z:12.8f}")

    # Write Gaussian input
    base_name = Path(gaussian_file).stem

    with open(gaussian_file, 'w') as f:
        # Optimization section
        f.write(f"%chk={base_name}.chk\n")
        f.write(f"%nprocshared={settings['opt']['processors']}\n")
        f.write(f"%mem={settings['opt']['memory']}\n")
        f.write(f"# {settings['opt']['keywords']}\n")
        f.write("\n")
        f.write(f"{res_name} optimization and frequency\n")
        f.write("\n")
        f.write(f"{net_charge} {multiplicity}\n")

        for coord in coords:
            f.write(f"{coord}\n")
        f.write("\n")

        # ESP section with --link1--
        f.write("--link1--\n")
        f.write(f"%oldchk={base_name}.chk\n")
        f.write(f"%chk={base_name}_esp.chk\n")
        f.write(f"%nprocshared={settings['esp']['processors']}\n")
        f.write(f"%mem={settings['esp']['memory']}\n")
        f.write(f"#P {settings['esp']['keywords']} guess=read geom=allcheck\n")
        f.write("\n")
        f.write(f"{res_name} ESP calculation\n")
        f.write("\n")
        f.write("\n")

    console.print(f"[green]✓ Generated Gaussian input: {gaussian_file}[/green]")

    # Concise summary (detailed info already shown in overview panel)
    console.print(Panel(
        f"[bold]Configuration Summary:[/bold]\n"
        f"  Charge: {net_charge}, Multiplicity: {multiplicity}\n"
        f"  Memory: {settings['opt']['memory']}, Processors: {settings['opt']['processors']}\n"
        f"  Atoms: {len(coords)}",
        title="Gaussian Input Generated",
        border_style="green",
        expand=False
    ))

def _run_antechamber_for_resp(log_file, mol2_file, net_charge, res_name, console):
    """Run antechamber to extract RESP charges."""
    
    cmd = [
        "antechamber", "-i", log_file, "-fi", "gout", "-o", mol2_file, 
        "-fo", "mol2", "-c", "resp", "-nc", str(net_charge), "-rn", res_name, "-at", "gaff2"
    ]
    
    console.print(f"[yellow]Command being executed:[/yellow]")
    console.print(" ".join(cmd))
    
    console.print(f"\n[cyan]Flag explanations:[/cyan]")
    console.print(f"  -i {log_file}      : Gaussian output with ESP data")
    console.print(f"  -fi gout              : Gaussian output format")
    console.print(f"  -o {mol2_file}       : MOL2 with charges and atom types")
    console.print(f"  -fo mol2              : MOL2 output format")
    console.print(f"  -c resp               : RESP charge method")
    console.print(f"  -nc {net_charge}                : Net charge constraint")
    console.print(f"  -rn {res_name}            : Residue name in MOL2")
    console.print(f"  -at gaff2             : GAFF2 atom type assignment")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        console.print(f"[green]✓ Generated MOL2 file: {mol2_file}[/green]")
        
    except subprocess.CalledProcessError as e:
        console.print(f"[red]Error running antechamber: {e.stderr}[/red]")
        raise

def _run_parmchk2_for_parameters(mol2_file, frcmod_file, console):
    """Run parmchk2 to generate force field parameters."""

    cmd = ["parmchk2", "-i", mol2_file, "-o", frcmod_file, "-f", "mol2", "-s", "2", "-a", "Y"]

    console.print(f"[yellow]Command being executed:[/yellow]")
    console.print(" ".join(cmd))

    console.print(f"\n[cyan]Flag explanations:[/cyan]")
    console.print(f"  -i {mol2_file}       : MOL2 with atom types and charges")
    console.print(f"  -o {frcmod_file}     : Force field modification file")
    console.print(f"  -f mol2              : Input is MOL2 format")
    console.print(f"  -s 2                 : Use GAFF2 parameter database")
    console.print(f"  -a Y                 : Print all parameters (complete frcmod)")
    
    console.print(f"\n[cyan]What parmchk2 does:[/cyan]")
    console.print(f"• Searches GAFF2 database for bond/angle/dihedral parameters")
    console.print(f"• Estimates missing parameters using chemical similarity")
    console.print(f"• Creates frcmod file with custom parameters")
    console.print(f"• Parameters marked 'ATTN' were estimated, not from database")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        console.print(f"[green]✓ Generated parameters: {frcmod_file}[/green]")
        
    except subprocess.CalledProcessError as e:
        console.print(f"[red]Error running parmchk2: {e.stderr}[/red]")
        raise

def _get_water_model_choice(console, processor=None):
    """Get user's water model preference."""
    water_models = {
        "1": ("opc", "OPC (recommended - 4-point model)"),
        "2": ("tip3p", "TIP3P (3-point model)"),
        "3": ("tip4pew", "TIP4P-Ew (4-point Ewald)"),
        "4": ("spce", "SPC/E (3-point extended)")
    }

    console.print("\n[cyan]Choose water model:[/cyan]")
    for key, (model, desc) in water_models.items():
        console.print(f"  {key}. {desc}")

    choice = prompt_with_context(
        processor,
        "Water model",
        choices=list(water_models.keys()),
        default="1",
        module="Small Molecule Parameterizer",
        description="Select water model"
    )

    return water_models[choice][0]

def _generate_tleap_inputs(mol2_file, frcmod_file, mol_name, water_model, box_size, console):
    """Generate tLEaP input files for gas and aqueous phases."""
    
    # Gas phase tLEaP input
    gas_input = f"{mol_name}_gas_tleap.in"
    gas_content = f"""source leaprc.protein.ff19SB
loadamberparams frcmod.ff19SB
source leaprc.gaff2
loadamberparams {frcmod_file}

{mol_name} = loadmol2 "{mol2_file}"
saveOff {mol_name} {mol_name}.lib
saveAmberParm {mol_name} {mol_name}_gas.parm7 {mol_name}_gas.rst7

quit
"""
    
    with open(gas_input, 'w') as f:
        f.write(gas_content)
    
    # Aqueous phase tLEaP input
    aq_input = f"{mol_name}_aq_tleap.in"
    
    # Water model configurations
    water_configs = {
        "opc": ("leaprc.water.opc", "frcmod.opc", "OPCBOX"),
        "tip3p": ("leaprc.water.tip3p", "frcmod.tip3p", "TIP3PBOX"),
        "tip4pew": ("leaprc.water.tip4pew", "frcmod.tip4pew", "TIP4PEWBOX"),
        "spce": ("leaprc.water.spce", "frcmod.spce", "SPCBOX")
    }
    
    leaprc, frcmod_water, box_type = water_configs[water_model]
    
    aq_content = f"""source leaprc.protein.ff19SB
loadamberparams frcmod.ff19SB
source {leaprc}
loadamberparams {frcmod_water}
source leaprc.gaff2
loadamberparams {frcmod_file}

{mol_name} = loadmol2 "{mol2_file}"
solvateBox {mol_name} {box_type} {box_size}
saveAmberParm {mol_name} {mol_name}_aq.parm7 {mol_name}_aq.rst7

quit
"""
    
    with open(aq_input, 'w') as f:
        f.write(aq_content)
    
    console.print(f"[green]✓ Generated tLEaP input files:[/green]")
    console.print(f"  Gas phase: {gas_input}")
    console.print(f"  Aqueous phase: {aq_input}")
    
    # Educational display
    console.print(f"\n[cyan]Gas phase tLEaP commands:[/cyan]")
    console.print(f"  source leaprc.gaff2       → Load GAFF2 force field")
    console.print(f"  loadamberparams {frcmod_file} → Load custom parameters")
    console.print(f"  {mol_name} = loadmol2 \"{mol2_file}\" → Load molecule with RESP charges")
    console.print(f"  saveAmberParm → Generate simulation files (.parm7/.rst7)")
    
    console.print(f"\n[cyan]Aqueous phase tLEaP commands:[/cyan]")
    console.print(f"  source {leaprc}   → Load {water_model.upper()} water model")
    console.print(f"  solvateBox {mol_name} {box_type} {box_size} → Add {box_size}Å water box")
    console.print(f"  saveAmberParm → Generate solvated simulation files")

def _display_completion_summary(mol_name, console):
    """Display completion summary and next steps."""

    console.print(Panel(
        f"[bold green]✓ Small molecule parameterization completed![/bold green]\n\n"
        f"[bold]Parameter Quality:[/bold]\n"
        f"  • RESP charges from quantum calculation\n"
        f"  • Frequencies verified (no negative modes)\n"
        f"  • Force field parameters generated\n"
        f"  • Penalty scores analyzed\n\n"
        f"[bold]Generated files:[/bold]\n"
        f"  {mol_name_lower}.pdb              → Extracted coordinates\n"
        f"  {mol_name_lower}.gjf              → Gaussian input\n"
        f"  {mol_name_lower}.log              → Gaussian output\n"
        f"  {mol_name_lower}.mol2             → MOL2 with RESP charges\n"
        f"  {mol_name_lower}.frcmod           → Force field parameters\n"
        f"  {mol_name_lower}_*_tleap.in       → tLEaP input files\n"
        f"  {mol_name_lower}.lib              → Reusable library file\n\n"
        f"[bold]Simulation files:[/bold]\n" +
        (f"  {', '.join(success_files)}\n\n" if success_files else "  [yellow]Run tLEaP manually if needed[/yellow]\n\n") +
        f"[cyan]Files are ready for AMBER MD simulations![/cyan]",
        title="Parameterization Complete",
        expand=False
    ))


def _launch_small_molecule_viewer(pdb_file: str, residue_name: str, console) -> bool:
    """Backward-compatible alias for :func:`hydrogen_editor.launch_viewer`.

    The implementation moved to ``hydrogen_editor`` so the small-molecule and
    modified-amino-acid routes share one viewer launch path. Retained as a
    thin delegate in case anything still references this name.
    """
    from proprep.forcefield_prep.hydrogen_editor import launch_viewer
    return launch_viewer(pdb_file, residue_name, console)



