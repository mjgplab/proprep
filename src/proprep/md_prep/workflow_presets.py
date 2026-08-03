"""
AMBER Workflow Preset System

Reads and manages workflow presets based on existing mdin template files.
Each preset contains a sequence of simulation steps with real parameters from templates.
"""

import json
import os
import re
from typing import Dict, List, Any, Optional
from pathlib import Path
from dataclasses import dataclass, asdict
from datetime import datetime

from rich.console import Console
from rich.table import Table
from rich.panel import Panel

from proprep.utils.paths import get_package_dir
from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
)

@dataclass
class WorkflowStep:
    """Definition of a single workflow step with actual mdin content."""
    id: str
    name: str
    type: str  # minimization, heating, equilibration, production
    description: str
    expected_runtime_gpu: str
    expected_runtime_cpu: str
    
    # File paths and content
    mdin_file: str
    mdin_content: str = ""
    
    # File handling
    input_coord: str = "inpcrd"
    output_coord: str = None
    dependencies: List[str] = None
    
    def __post_init__(self):
        if self.dependencies is None:
            self.dependencies = []
        if self.output_coord is None:
            self.output_coord = f"{self.id}.rst"

@dataclass 
class WorkflowPreset:
    """Complete workflow preset definition."""
    workflow_name: str
    description: str
    version: str
    category: str  # protein, membrane, small_molecule
    
    steps: List[WorkflowStep]
    
    # Global settings
    global_settings: Dict[str, Any]
    estimated_time: Dict[str, str]
    
    # Metadata
    created: str = None
    mdin_source_dir: str = None
    
    def __post_init__(self):
        if self.created is None:
            self.created = datetime.now().isoformat()

class WorkflowPresetSystem:
    """System for managing workflow presets based on real mdin files."""
    
    def __init__(self, console: Console = None, mdin_dir: Path = None, processor=None):
        self.console = console or Console()
        self.processor = processor
        self.presets: Dict[str, WorkflowPreset] = {}
        
        # Set mdin directory path
        if mdin_dir is None:
            self.mdin_dir = get_package_dir() / "mdins"
        else:
            self.mdin_dir = Path(mdin_dir)
            
        self._load_builtin_presets()
    
    def _parse_mdin_file(self, filepath: Path) -> str:
        """Read and return content of mdin file."""
        try:
            with open(filepath, 'r') as f:
                return f.read()
        except Exception as e:
            self.console.print(f"[red]Error reading {filepath}: {e}[/red]")
            return ""
    
    def _extract_simulation_info(self, mdin_content: str, filename: str) -> Dict[str, Any]:
        """Extract key simulation parameters from mdin content."""
        info = {
            "is_minimization": False,
            "is_heating": False,
            "has_restraints": False,
            "restraint_weight": 0.0,
            "nstlim": 0,
            "dt": 0.001,
            "simulation_time_ps": 0
        }
        
        # Check if minimization
        if re.search(r'imin\s*=\s*1', mdin_content):
            info["is_minimization"] = True
            
        # Check for heating (nmropt=1 with TEMP0 ramping)
        if re.search(r'nmropt\s*=\s*1', mdin_content) and 'TEMP0' in mdin_content:
            info["is_heating"] = True
            
        # Extract restraint weight
        restraint_match = re.search(r'restraint_wt\s*=\s*([\d.]+)', mdin_content)
        if restraint_match:
            info["has_restraints"] = True
            info["restraint_weight"] = float(restraint_match.group(1))
            
        # Extract nstlim and dt
        nstlim_match = re.search(r'nstlim\s*=\s*(\d+)', mdin_content)
        if nstlim_match:
            info["nstlim"] = int(nstlim_match.group(1))
            
        dt_match = re.search(r'dt\s*=\s*([\d.]+)', mdin_content)
        if dt_match:
            info["dt"] = float(dt_match.group(1))
            
        # Calculate simulation time
        if not info["is_minimization"]:
            info["simulation_time_ps"] = info["nstlim"] * info["dt"]
            
        return info
    
    def _load_builtin_presets(self):
        """Load built-in workflow presets from actual mdin files."""
        
        if not self.mdin_dir.exists():
            self.console.print(f"[yellow]Warning: mdin directory not found: {self.mdin_dir}[/yellow]")
            return
            
        # Define the protein equilibration workflow sequence
        protein_workflow = [
            {
                "filename": "000_Min.mdin",
                "id": "000_Min",
                "name": "Initial Solvent Minimization",
                "type": "minimization", 
                "description": "Energy minimize solvent with solute restrained",
                "runtime_gpu": "30 sec",
                "runtime_cpu": "2 min",
                "input_coord": "inpcrd",
                "dependencies": []
            },
            {
                "filename": "001_NPT.mdin",
                "id": "001_NPT", 
                "name": "Initial NPT Equilibration",
                "type": "equilibration",
                "description": "Initial density equilibration (5 ps)",
                "runtime_gpu": "2 min",
                "runtime_cpu": "15 min", 
                "input_coord": "000_Min.rst",
                "dependencies": ["000_Min"]
            },
            {
                "filename": "002_NPT.mdin",
                "id": "002_NPT",
                "name": "Continued NPT Equilibration", 
                "type": "equilibration",
                "description": "Continued density equilibration (5 ps)",
                "runtime_gpu": "2 min",
                "runtime_cpu": "15 min",
                "input_coord": "001_NPT.rst",
                "dependencies": ["001_NPT"]
            },
            {
                "filename": "003_NPT.mdin",
                "id": "003_NPT",
                "name": "Extended NPT Equilibration",
                "type": "equilibration",
                "description": "Extended density equilibration (190 ps)", 
                "runtime_gpu": "15 min",
                "runtime_cpu": "2 hours",
                "input_coord": "002_NPT.rst",
                "dependencies": ["002_NPT"]
            },
            {
                "filename": "01_Min.mdin",
                "id": "01_Min",
                "name": "Solvent Minimization",
                "type": "minimization",
                "description": "Energy minimize solvent with solute restrained",
                "runtime_gpu": "30 sec",
                "runtime_cpu": "2 min",
                "input_coord": "003_NPT.rst",
                "dependencies": ["003_NPT"]
            },
            {
                "filename": "02_Heat.mdin", 
                "id": "02_Heat",
                "name": "Solvent Heating",
                "type": "heating",
                "description": "Heat from 0→300K over 600ps + 1ns hold",
                "runtime_gpu": "45 min",
                "runtime_cpu": "6 hours",
                "input_coord": "01_Min.rst",
                "dependencies": ["01_Min"]
            },
            {
                "filename": "03_Equil_NPT_solvent.mdin",
                "id": "03_Equil_NPT",
                "name": "Solvent NPT Equilibration", 
                "type": "equilibration",
                "description": "NPT equilibration at 300K (5 ns)",
                "runtime_gpu": "3 hours",
                "runtime_cpu": "24 hours",
                "input_coord": "02_Heat.rst",
                "dependencies": ["02_Heat"]
            },
            {
                "filename": "04_Min_solute.mdin",
                "id": "04_Min_solute",
                "name": "Solute Minimization",
                "type": "minimization", 
                "description": "Minimize solute (no restraints)",
                "runtime_gpu": "30 sec",
                "runtime_cpu": "2 min",
                "input_coord": "03_Equil_NPT.rst",
                "dependencies": ["03_Equil_NPT"]
            },
            {
                "filename": "05_Heat_solute_Rwt25.mdin",
                "id": "05_Heat_Rwt25",
                "name": "Solute Heating (25 kcal/mol·Å²)",
                "type": "heating",
                "description": "Heat solute with 25 kcal/mol·Å² restraints",
                "runtime_gpu": "45 min", 
                "runtime_cpu": "6 hours",
                "input_coord": "04_Min_solute.rst",
                "dependencies": ["04_Min_solute"]
            },
            {
                "filename": "06_Equil_solute_Rwt25.mdin",
                "id": "06_Equil_Rwt25",
                "name": "Solute Equilibration (25 kcal/mol·Å²)",
                "type": "equilibration",
                "description": "NPT equilibration with 25 kcal/mol·Å² restraints (500 ps)",
                "runtime_gpu": "30 min",
                "runtime_cpu": "4 hours",
                "input_coord": "05_Heat_Rwt25.rst",
                "dependencies": ["05_Heat_Rwt25"]
            },
            {
                "filename": "07_Equil_solute_Rwt10.mdin", 
                "id": "07_Equil_Rwt10",
                "name": "Solute Equilibration (10 kcal/mol·Å²)",
                "type": "equilibration", 
                "description": "NPT equilibration with 10 kcal/mol·Å² restraints (500 ps)",
                "runtime_gpu": "30 min",
                "runtime_cpu": "4 hours",
                "input_coord": "06_Equil_Rwt25.rst",
                "dependencies": ["06_Equil_Rwt25"]
            },
            {
                "filename": "08_Equil_solute_Rwt5.mdin",
                "id": "08_Equil_Rwt5", 
                "name": "Solute Equilibration (5 kcal/mol·Å²)",
                "type": "equilibration",
                "description": "NPT equilibration with 5 kcal/mol·Å² restraints (500 ps)",
                "runtime_gpu": "30 min",
                "runtime_cpu": "4 hours",
                "input_coord": "07_Equil_Rwt10.rst",
                "dependencies": ["07_Equil_Rwt10"]
            },
            {
                "filename": "09_Equil_solute_Rwt2.mdin",
                "id": "09_Equil_Rwt2",
                "name": "Solute Equilibration (2 kcal/mol·Å²)",
                "type": "equilibration",
                "description": "Final NPT equilibration with 2 kcal/mol·Å² restraints (500 ps)",
                "runtime_gpu": "30 min",
                "runtime_cpu": "4 hours", 
                "input_coord": "08_Equil_Rwt5.rst",
                "dependencies": ["08_Equil_Rwt5"]
            }
        ]
        
        # Create workflow steps by reading actual mdin files
        steps = []
        total_simulation_time = 0.0
        
        for step_def in protein_workflow:
            mdin_path = self.mdin_dir / step_def["filename"]
            mdin_content = self._parse_mdin_file(mdin_path)
            
            if mdin_content:
                # Extract simulation info
                sim_info = self._extract_simulation_info(mdin_content, step_def["filename"])
                
                # Update description with actual simulation time if available
                description = step_def["description"]
                if not sim_info["is_minimization"] and sim_info["simulation_time_ps"] > 0:
                    time_ps = sim_info["simulation_time_ps"]
                    if time_ps >= 1000:
                        time_str = f"{time_ps/1000:.1f} ns"
                    else:
                        time_str = f"{time_ps:.0f} ps"
                    description += f" [{time_str}]"
                    total_simulation_time += time_ps / 1000  # Convert to ns
                
                step = WorkflowStep(
                    id=step_def["id"],
                    name=step_def["name"],
                    type=step_def["type"],
                    description=description,
                    expected_runtime_gpu=step_def["runtime_gpu"],
                    expected_runtime_cpu=step_def["runtime_cpu"],
                    mdin_file=step_def["filename"],
                    mdin_content=mdin_content,
                    input_coord=step_def["input_coord"],
                    output_coord=f"{step_def['id']}.rst",
                    dependencies=step_def["dependencies"]
                )
                steps.append(step)
            else:
                self.console.print(f"[yellow]Warning: Could not read {mdin_path}[/yellow]")
                
        if steps:
            protein_preset = WorkflowPreset(
                workflow_name="Standard Protein Equilibration",
                description="Complete solvent + solute equilibration protocol for proteins in water",
                version="1.0",
                category="protein",
                steps=steps,
                global_settings={
                    "system_type": "protein_in_water",
                    "temperature": 300,
                    "pressure": 1.0,
                    "default_restraint_mask": "!@H= & !:NA,CL,MG,WAT",
                    "cutoff": 12.0,
                    "time_step": 0.001
                },
                estimated_time={
                    "gpu": "6-8 hours", 
                    "cpu": "48-72 hours",
                    "total_simulation_time": f"{total_simulation_time:.1f} ns"
                },
                mdin_source_dir=str(self.mdin_dir)
            )
            
            self.presets["protein_equilibration"] = protein_preset
    
    def get_available_presets(self) -> Dict[str, WorkflowPreset]:
        """Get all available workflow presets."""
        return self.presets.copy()
    
    def display_preset_menu(self) -> Optional[str]:
        """Display preset selection menu and return selected preset key."""
        self.console.print("\n[bold cyan]===== Available Workflow Presets =====[/bold cyan]")
        
        if not self.presets:
            self.console.print("[yellow]No presets available[/yellow]")
            return None
            
        # Built-in presets
        self.console.print("\n[bold]Built-in Presets:[/bold]")
        preset_choices = {}
        choice_num = 1
        
        for key, preset in self.presets.items():
            if not key.startswith("custom_"):
                total_time = preset.estimated_time.get("total_simulation_time", "unknown")
                
                self.console.print(f"  {choice_num}. {preset.workflow_name} ({len(preset.steps)} steps)")
                self.console.print(f"     • {preset.description}")
                self.console.print(f"     • Duration: ~{total_time} total simulation time")
                self.console.print(f"     • Restraint weights: 50 → 25 → 10 → 5 → 2 kcal/mol·Å²")
                
                preset_choices[str(choice_num)] = key
                choice_num += 1
        
        # Custom presets
        custom_presets = {k: v for k, v in self.presets.items() if k.startswith("custom_")}
        if custom_presets:
            self.console.print("\n[bold]Custom Presets:[/bold]")
            for key, preset in custom_presets.items():
                self.console.print(f"  {choice_num}. {preset.workflow_name} (user created)")
                preset_choices[str(choice_num)] = key
                choice_num += 1
                
        if not preset_choices:
            return None
        
        # Add return option
        self.console.print(f"  {choice_num}. ← Return to AMBER Input Generator")
        preset_choices[str(choice_num)] = None  # Return option
            
        choices = list(preset_choices.keys())
        preset_opts = {
            k: (p.name if p else "Return to AMBER Input Generator")
            for k, p in preset_choices.items()
        }
        choice = prompt_with_context(
            self.processor,
            "Select preset",
            choices=choices,
            default="1",
            module="MD Manager - Workflow Presets",
            description="Select workflow preset",
            options_map=preset_opts,
        )

        return preset_choices.get(choice)
    
    def display_preset_overview(self, preset: WorkflowPreset) -> bool:
        """Display detailed overview of a workflow preset."""
        self.console.print(f"\n[bold cyan]===== {preset.workflow_name} =====[/bold cyan]")
        
        # Create workflow overview table
        table = Table(title="Workflow Steps")
        table.add_column("Step", style="cyan", width=13)
        table.add_column("Type", style="green", width=15)
        table.add_column("Time", style="yellow", width=10)
        table.add_column("Restraints", style="magenta", width=12)  
        table.add_column("Weight", style="red", width=10)
        
        for step in preset.steps:
            # Extract time from description or parse mdin
            time_match = re.search(r'\\[(.*?)\\]', step.description)
            if time_match:
                time_str = time_match.group(1)
            elif step.type == "minimization":
                time_str = "-"
            else:
                # Try to extract from mdin content
                sim_info = self._extract_simulation_info(step.mdin_content, step.mdin_file)
                if sim_info["simulation_time_ps"] > 0:
                    time_ps = sim_info["simulation_time_ps"]
                    time_str = f"{time_ps/1000:.1f} ns" if time_ps >= 1000 else f"{time_ps:.0f} ps"
                else:
                    time_str = "?"
            
            # Extract restraint info from mdin content
            sim_info = self._extract_simulation_info(step.mdin_content, step.mdin_file)
            if sim_info["has_restraints"]:
                weight = sim_info["restraint_weight"]
                restraint_type = "Solvent" if weight >= 25 else "Solute"
                weight_str = f"{weight:.1f}"
            else:
                restraint_type = "None"
                weight_str = "0.0"
                
            table.add_row(
                step.id,
                step.type.title(),
                time_str,
                restraint_type,
                weight_str
            )
            
        self.console.print(table)
        
        # Show estimated runtime
        total_time = preset.estimated_time.get("total_simulation_time", "~8.6 ns")
        self.console.print(f"\nTotal simulation time: {total_time}")
        self.console.print(f"Estimated runtime: {preset.estimated_time['gpu']} (GPU) / {preset.estimated_time['cpu']} (CPU)")
        
        return True
    
    def display_step_menu(self, preset: WorkflowPreset) -> Optional[int]:
        """Display step selection menu for editing."""
        self.console.print(f"\n[bold cyan]===== Edit Individual Workflow Steps =====[/bold cyan]")
        
        self.console.print("Select step to modify:")
        step_choices = {}
        
        for i, step in enumerate(preset.steps, 1):
            # Get time info for display
            sim_info = self._extract_simulation_info(step.mdin_content, step.mdin_file)
            
            if step.type == "minimization":
                time_info = ""
            elif sim_info["simulation_time_ps"] > 0:
                time_ps = sim_info["simulation_time_ps"]
                time_str = f"{time_ps/1000:.1f} ns" if time_ps >= 1000 else f"{time_ps:.0f} ps"
                time_info = f" ({time_str})"
            else:
                time_info = ""
                
            # Get restraint info
            restraint_info = ""
            if sim_info["has_restraints"]:
                restraint_info = f" [{sim_info['restraint_weight']:.0f} kcal/mol·Å²]"
                
            self.console.print(f"  {i}. {step.id:<12} - {step.name}{time_info}{restraint_info}")
            step_choices[str(i)] = i - 1  # Convert to 0-based index
            
        self.console.print(f"  {len(preset.steps) + 1}. Add new step")
        self.console.print(f"  {len(preset.steps) + 2}. Remove step")
        self.console.print(f"  {len(preset.steps) + 3}. Reorder steps")
        
        # Add special choices
        step_choices[str(len(preset.steps) + 1)] = -1  # Add new
        step_choices[str(len(preset.steps) + 2)] = -2  # Remove
        step_choices[str(len(preset.steps) + 3)] = -3  # Reorder
        
        choices = list(step_choices.keys())
        choice = prompt_with_context(
            self.processor,
            "Select option",
            choices=choices,
            default="1",
            module="MD Manager - Workflow Presets",
            description="Edit preset step (select, add, remove, or reorder)",
        )

        return step_choices.get(choice)
    
    def save_preset_to_file(self, preset: WorkflowPreset, filepath: Path):
        """Save preset to JSON file."""
        preset_dict = asdict(preset) 
        with open(filepath, 'w') as f:
            json.dump(preset_dict, f, indent=2)
    
    def load_preset_from_file(self, filepath: Path) -> WorkflowPreset:
        """Load preset from JSON file."""
        with open(filepath, 'r') as f:
            preset_dict = json.load(f)
            
        # Reconstruct WorkflowStep objects
        steps = []
        for step_dict in preset_dict['steps']:
            steps.append(WorkflowStep(**step_dict))
            
        preset_dict['steps'] = steps
        return WorkflowPreset(**preset_dict)
    
    def generate_workflow_sequence_json(self, preset: WorkflowPreset, output_dir: Path = None) -> Path:
        """Generate JSON workflow sequence file for the Workflow Manager."""
        if output_dir is None:
            output_dir = Path.cwd()
            
        sequence_file = output_dir / f"{preset.workflow_name.lower().replace(' ', '_')}_workflow.json"
        
        sequence_data = {
            "workflow_name": preset.workflow_name,
            "description": preset.description,
            "version": preset.version,
            "created": datetime.now().isoformat(),
            "total_steps": len(preset.steps),
            "estimated_time": preset.estimated_time,
            "global_settings": preset.global_settings,
            "steps": [],
            "file_paths": {
                "workflow_dir": str(output_dir),
                "mdin_files": [],
                "execution_scripts": [
                    "run_step.sh",
                    "run_workflow.sh"
                ]
            },
            "workspace_integration": {
                "processor": "amber_input_generator",
                "generated_timestamp": datetime.now().isoformat(),
                "source": "workflow_preset_system"
            }
        }
        
        for step in preset.steps:
            step_data = {
                "id": step.id,
                "name": step.name,
                "type": step.type,
                "mdin_file": step.mdin_file,
                "dependencies": step.dependencies,
                "input_coord": step.input_coord,
                "output_coord": step.output_coord,
                "expected_runtime_gpu": step.expected_runtime_gpu,
                "expected_runtime_cpu": step.expected_runtime_cpu,
                "description": step.description
            }
            sequence_data["steps"].append(step_data)
            sequence_data["file_paths"]["mdin_files"].append(step.mdin_file)
            
        with open(sequence_file, 'w') as f:
            json.dump(sequence_data, f, indent=2)
            
        return sequence_file