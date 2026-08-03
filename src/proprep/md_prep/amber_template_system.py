"""
AMBER Annotated Template System

Provides rich, annotated templates for AMBER simulations with parameter explanations
and seamless integration with the AmberWizard for configuration.
"""

import json
import os
from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path
from dataclasses import dataclass, asdict
from enum import Enum

from rich.console import Console
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text

from proprep.utils.paths import get_package_dir


class SimulationType(Enum):
    """AMBER simulation types."""
    MINIMIZATION = "minimization"
    HEATING = "heating" 
    EQUILIBRATION = "equilibration"
    PRODUCTION = "production"
    CUSTOM = "custom"


@dataclass
class ParameterAnnotation:
    """Annotation for a single AMBER parameter."""
    name: str
    value: Any
    description: str
    typical_values: List[str]
    warnings: List[str] = None
    dependencies: List[str] = None
    category: str = "general"
    
    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []
        if self.dependencies is None:
            self.dependencies = []


@dataclass 
class AmberTemplate:
    """
    Rich AMBER simulation template with annotations and metadata.
    """
    name: str
    description: str
    simulation_type: SimulationType
    parameters: Dict[str, ParameterAnnotation]
    nmr_section: Optional[str] = None
    comments: List[str] = None
    metadata: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.comments is None:
            self.comments = []
        if self.metadata is None:
            self.metadata = {}
            
    def get_config_dict(self) -> Dict[str, Any]:
        """Extract just the parameter values for AMBER input file generation."""
        return {name: param.value for name, param in self.parameters.items()}
        
    def update_from_wizard_config(self, wizard_config: Dict[str, Any]) -> None:
        """Update template parameters from wizard configuration results."""
        for param_name, value in wizard_config.items():
            if param_name in self.parameters:
                # Update existing parameter
                self.parameters[param_name].value = value
            else:
                # Add new parameter with basic annotation
                self.parameters[param_name] = ParameterAnnotation(
                    name=param_name,
                    value=value,
                    description=f"Parameter {param_name} configured by wizard",
                    typical_values=[str(value)],
                    category="wizard_configured"
                )
                
    def generate_mdin_content(self) -> str:
        """Generate complete AMBER mdin file content from template."""
        lines = []
        
        # Title and description
        lines.append(f"{self.name}")
        lines.append("&cntrl")
        
        # Group parameters by category for organized output
        categories = {}
        for param in self.parameters.values():
            cat = param.category
            if cat not in categories:
                categories[cat] = []
            categories[cat].append(param)
            
        # Write parameters by category with comments
        for category, params in categories.items():
            if params:
                lines.append(f"! {category.replace('_', ' ').title()} Parameters")
                
                for param in params:
                    # Add parameter with inline comment
                    comment = param.description[:50] + "..." if len(param.description) > 50 else param.description
                    lines.append(f"  {param.name}={param.value}, ! {comment}")
                    
                lines.append("")  # Blank line between categories
                
        lines.append("&end")
        
        # Add NMR section if present
        if self.nmr_section:
            lines.append("")
            lines.append("! Varying conditions (&wt blocks)")
            lines.extend(self.nmr_section.split('\n'))
            
        return '\n'.join(lines)


class AmberTemplateSystem:
    """
    Manages AMBER simulation templates with rich annotations and wizard integration.
    """
    
    def __init__(self, package_dir: Path = None, console: Console = None):
        """
        Initialize the template system.
        
        Args:
            package_dir: Base package directory for template storage
            console: Rich console for output
        """
        self.package_dir = package_dir or (get_package_dir() / "md_prep")
        self.console = console or Console()
        
        # Template storage paths
        self.templates_dir = self.package_dir / "templates" / "annotated"
        self.user_templates_dir = self.templates_dir / "user"
        self.builtin_templates_dir = self.templates_dir / "builtin"
        
        # Ensure directories exist
        self._ensure_directories()
        
        # Initialize with default templates
        self._create_default_templates()
        
    def _ensure_directories(self):
        """Ensure all template directories exist."""
        for directory in [self.templates_dir, self.user_templates_dir, self.builtin_templates_dir]:
            directory.mkdir(parents=True, exist_ok=True)
            
    def _create_default_templates(self):
        """Create comprehensive default templates for each simulation type."""
        
        # Minimization template
        minimization = AmberTemplate(
            name="Energy Minimization",
            description="Standard energy minimization with restraints",
            simulation_type=SimulationType.MINIMIZATION,
            parameters={
                "imin": ParameterAnnotation(
                    name="imin",
                    value=1,
                    description="Enable minimization (imin=1)",
                    typical_values=["1"],
                    category="simulation_control"
                ),
                "ntx": ParameterAnnotation(
                    name="ntx",
                    value=1,
                    description="Read coordinates only (no velocities)",
                    typical_values=["1"],
                    category="input_control"
                ),
                "maxcyc": ParameterAnnotation(
                    name="maxcyc",
                    value=10000,
                    description="Maximum number of minimization cycles",
                    typical_values=["5000", "10000", "50000"],
                    warnings=["Very large values may indicate convergence problems"],
                    category="minimization"
                ),
                "ncyc": ParameterAnnotation(
                    name="ncyc",
                    value=1000,
                    description="Steepest descent cycles before switching to conjugate gradient",
                    typical_values=["500", "1000", "2000"],
                    dependencies=["maxcyc"],
                    category="minimization"
                ),
                "ntmin": ParameterAnnotation(
                    name="ntmin",
                    value=2,
                    description="Minimization algorithm (1=steepest descent only, 2=steepest+conjugate)",
                    typical_values=["1", "2"],
                    category="minimization"
                ),
                "dx0": ParameterAnnotation(
                    name="dx0", 
                    value=0.01,
                    description="Initial step size for minimization",
                    typical_values=["0.01", "0.1", "1.0"],
                    warnings=["Large values may cause instability"],
                    category="minimization"
                ),
                "drms": ParameterAnnotation(
                    name="drms",
                    value=0.0001,
                    description="RMS gradient convergence criterion",
                    typical_values=["0.0001", "0.001", "0.01"],
                    category="minimization"
                ),
                "cut": ParameterAnnotation(
                    name="cut",
                    value=10.0,
                    description="Nonbonded cutoff distance in Angstroms",
                    typical_values=["8.0", "10.0", "12.0"],
                    warnings=["Values < 8.0 may affect electrostatics", "Values > 12.0 require large boxes"],
                    category="nonbonded"
                ),
                "ntpr": ParameterAnnotation(
                    name="ntpr",
                    value=100,
                    description="Print energy every ntpr steps",
                    typical_values=["50", "100", "500"],
                    category="output_control"
                ),
                "ntr": ParameterAnnotation(
                    name="ntr",
                    value=1,
                    description="Enable positional restraints (1=yes, 0=no)",
                    typical_values=["0", "1"],
                    dependencies=["restraint_wt", "restraintmask"],
                    category="restraints"
                ),
                "restraint_wt": ParameterAnnotation(
                    name="restraint_wt", 
                    value=500.0,
                    description="Restraint force constant in kcal/mol/Å²",
                    typical_values=["1.0", "10.0", "100.0", "500.0"],
                    dependencies=["ntr"],
                    category="restraints"
                ),
                "restraintmask": ParameterAnnotation(
                    name="restraintmask",
                    value="!@H=",
                    description="Atom selection for restraints (!@H= means all heavy atoms)",
                    typical_values=["!@H=", "@CA,C,N,O", "!:WAT", ":1-100@CA"],
                    dependencies=["ntr"],
                    category="restraints"
                )
            },
            comments=[
                "Standard minimization protocol with restraints on heavy atoms",
                "Adjust restraint_wt and restraintmask as needed for your system",
                "Consider reducing restraints or removing them entirely for final minimization"
            ]
        )
        
        # Heating template
        heating = AmberTemplate(
            name="Heating Protocol", 
            description="Gradual heating from 0 to 300K with restraints",
            simulation_type=SimulationType.HEATING,
            parameters={
                "imin": ParameterAnnotation(
                    name="imin",
                    value=0,
                    description="Molecular dynamics simulation (imin=0)",
                    typical_values=["0"],
                    category="simulation_control"
                ),
                "ntx": ParameterAnnotation(
                    name="ntx", 
                    value=1,
                    description="Read coordinates only (new simulation)",
                    typical_values=["1", "5"],
                    dependencies=["irest"],
                    category="input_control"
                ),
                "irest": ParameterAnnotation(
                    name="irest",
                    value=0,
                    description="Do not restart simulation (new start)",
                    typical_values=["0", "1"],
                    dependencies=["ntx"],
                    category="input_control"
                ),
                "nstlim": ParameterAnnotation(
                    name="nstlim",
                    value=50000,
                    description="Number of MD steps (50k steps = 50 ps with dt=0.001)",
                    typical_values=["25000", "50000", "100000"],
                    dependencies=["dt"],
                    category="time_control"
                ),
                "dt": ParameterAnnotation(
                    name="dt",
                    value=0.001,
                    description="Time step in ps (1 fs for heating)",
                    typical_values=["0.001", "0.002"],
                    warnings=["Use 0.001 ps for heating phases", "0.002 ps requires SHAKE constraints"],
                    category="time_control"
                ),
                "ntt": ParameterAnnotation(
                    name="ntt",
                    value=3,
                    description="Langevin thermostat (recommended)",
                    typical_values=["1", "2", "3", "11"],
                    dependencies=["gamma_ln", "temp0", "tempi"],
                    category="temperature"
                ),
                "tempi": ParameterAnnotation(
                    name="tempi",
                    value=0.0,
                    description="Initial temperature in K (start from 0K)",
                    typical_values=["0.0", "100.0", "200.0"],
                    dependencies=["nmropt"],
                    category="temperature"
                ),
                "temp0": ParameterAnnotation(
                    name="temp0",
                    value=300.0, 
                    description="Target temperature in K",
                    typical_values=["300.0", "310.0", "350.0"],
                    warnings=["High temperatures (>350K) may require smaller time steps"],
                    category="temperature"
                ),
                "gamma_ln": ParameterAnnotation(
                    name="gamma_ln",
                    value=2.0,
                    description="Langevin collision frequency in ps⁻¹",
                    typical_values=["1.0", "2.0", "5.0"],
                    dependencies=["ntt"],
                    category="temperature"
                ),
                "ntb": ParameterAnnotation(
                    name="ntb",
                    value=1,
                    description="Constant volume periodic boundary (NVT)",
                    typical_values=["1", "2"],
                    warnings=["Use ntb=1 for heating, ntb=2 for production"],
                    category="system_setup"
                ),
                "ntc": ParameterAnnotation(
                    name="ntc",
                    value=2,
                    description="SHAKE constraint algorithm (bonds with hydrogen)",
                    typical_values=["1", "2", "3"],
                    dependencies=["ntf", "tol"],
                    category="constraints"
                ),
                "ntf": ParameterAnnotation(
                    name="ntf", 
                    value=2,
                    description="Omit bond forces with hydrogen (when ntc=2)",
                    typical_values=["1", "2", "3"], 
                    dependencies=["ntc"],
                    category="constraints"
                ),
                "cut": ParameterAnnotation(
                    name="cut",
                    value=10.0,
                    description="Nonbonded cutoff distance in Angstroms",
                    typical_values=["8.0", "10.0", "12.0"],
                    warnings=["Values < 8.0 may affect electrostatics"],
                    category="nonbonded"
                ),
                "nmropt": ParameterAnnotation(
                    name="nmropt",
                    value=1,
                    description="Enable varying conditions (&wt blocks for temperature ramping)",
                    typical_values=["0", "1"],
                    dependencies=["tempi", "temp0"],
                    category="advanced"
                ),
                "ntr": ParameterAnnotation(
                    name="ntr",
                    value=1,
                    description="Enable positional restraints",
                    typical_values=["0", "1"],
                    dependencies=["restraint_wt", "restraintmask"],
                    category="restraints"
                ),
                "restraint_wt": ParameterAnnotation(
                    name="restraint_wt",
                    value=10.0,
                    description="Restraint force constant in kcal/mol/Å²",
                    typical_values=["1.0", "10.0", "25.0"],
                    category="restraints"
                ),
                "restraintmask": ParameterAnnotation(
                    name="restraintmask",
                    value="!@H=",
                    description="Restrain all heavy atoms during heating",
                    typical_values=["!@H=", "@CA,C,N,O", "!:WAT"],
                    category="restraints"
                ),
                "ntpr": ParameterAnnotation(
                    name="ntpr",
                    value=500,
                    description="Print energy every 500 steps",
                    typical_values=["250", "500", "1000"],
                    category="output_control"
                ),
                "ntwx": ParameterAnnotation(
                    name="ntwx",
                    value=500,
                    description="Write trajectory every 500 steps",
                    typical_values=["250", "500", "1000"],
                    category="output_control"
                ),
                "ntwr": ParameterAnnotation(
                    name="ntwr",
                    value=2500,
                    description="Write restart file every 2500 steps",
                    typical_values=["1000", "2500", "5000"],
                    category="output_control"
                ),
                "ntxo": ParameterAnnotation(
                    name="ntxo",
                    value=2,
                    description="NetCDF format for restart file",
                    typical_values=["1", "2"],
                    category="output_control"
                ),
                "ioutfm": ParameterAnnotation(
                    name="ioutfm",
                    value=1,
                    description="NetCDF format for trajectory file",
                    typical_values=["0", "1"], 
                    category="output_control"
                )
            },
            nmr_section="&wt type='TEMP0', istep1=0, istep2=25000,\n    value1=0.0, value2=300.0 /\n&wt type='END' /",
            comments=[
                "Gradual heating protocol from 0K to 300K over 50 ps",
                "Uses Langevin thermostat with moderate restraints",
                "Temperature ramping controlled by &wt blocks",
                "Adjust heating rate by modifying nstlim and &wt blocks"
            ]
        )
        
        # Save default templates
        self._save_template(minimization, self.builtin_templates_dir)
        self._save_template(heating, self.builtin_templates_dir)
        
        # Create equilibration and production templates similarly...
        # (I'll add these in the next step to keep this manageable)
        
    def _save_template(self, template: AmberTemplate, directory: Path):
        """Save a template to disk."""
        filename = f"{template.simulation_type.value}.json"
        filepath = directory / filename
        
        # Convert to serializable format
        template_dict = {
            "name": template.name,
            "description": template.description,
            "simulation_type": template.simulation_type.value,
            "parameters": {name: asdict(param) for name, param in template.parameters.items()},
            "nmr_section": template.nmr_section,
            "comments": template.comments,
            "metadata": template.metadata
        }
        
        with open(filepath, 'w') as f:
            json.dump(template_dict, f, indent=2)
            
    def load_template(self, simulation_type: SimulationType, template_name: str = None) -> AmberTemplate:
        """
        Load a template by simulation type and optional name.
        
        Args:
            simulation_type: Type of simulation
            template_name: Specific template name (defaults to builtin)
            
        Returns:
            Loaded template
        """
        # Try user templates first, then builtin
        for directory in [self.user_templates_dir, self.builtin_templates_dir]:
            filename = f"{simulation_type.value}.json"
            filepath = directory / filename
            
            if filepath.exists():
                with open(filepath, 'r') as f:
                    data = json.load(f)
                    
                # Reconstruct template
                parameters = {}
                for name, param_data in data["parameters"].items():
                    parameters[name] = ParameterAnnotation(**param_data)
                    
                return AmberTemplate(
                    name=data["name"],
                    description=data["description"],
                    simulation_type=SimulationType(data["simulation_type"]),
                    parameters=parameters,
                    nmr_section=data.get("nmr_section"),
                    comments=data.get("comments", []),
                    metadata=data.get("metadata", {})
                )
                
        raise FileNotFoundError(f"Template for {simulation_type.value} not found")
        
    def display_template_overview(self, template: AmberTemplate) -> None:
        """Display rich overview of template parameters."""
        # Title panel
        title = f"[bold cyan]{template.name}[/bold cyan]"
        subtitle = f"[grey50]{template.description}[/grey50]"
        
        self.console.print(Panel(f"{title}\n{subtitle}", 
                                title="AMBER Template Overview",
                                border_style="cyan"))
        
        # Parameter table by category
        categories = {}
        for param in template.parameters.values():
            cat = param.category.replace("_", " ").title()
            if cat not in categories:
                categories[cat] = []
            categories[cat].append(param)
            
        for category, params in categories.items():
            table = Table(title=category, show_header=True, header_style="bold blue")
            table.add_column("Parameter", style="cyan", width=12)
            table.add_column("Value", style="green", width=10)
            table.add_column("Description", style="white", width=40)
            table.add_column("Typical", style="yellow", width=15)
            
            for param in params:
                typical = ", ".join(param.typical_values[:3])  # Show first 3 typical values
                table.add_row(param.name, str(param.value), param.description[:40] + "...", typical)
                
            self.console.print(table)
            self.console.print()  # Blank line
            
        # Show comments if present
        if template.comments:
            comment_text = "\n".join([f"• {comment}" for comment in template.comments])
            self.console.print(Panel(comment_text, title="Notes", border_style="yellow"))
            
    def display_template_content(self, template: AmberTemplate) -> None:
        """Display the generated AMBER input file content."""
        content = template.generate_mdin_content()
        syntax = Syntax(content, "fortran", theme="monokai", line_numbers=True)
        
        self.console.print(Panel(syntax, 
                                title=f"Generated Input File: {template.name}",
                                border_style="green"))