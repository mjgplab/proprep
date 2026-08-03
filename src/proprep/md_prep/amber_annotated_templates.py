"""
AMBER Annotated Templates - Simple Clean System

Just like the original template system, but with inline comments added to the generated mdin files.
Exactly what the user wanted - no complex tables, just helpful comments after each parameter.
"""

import json
import logging
import os
from typing import Dict, Any, List, Optional
from pathlib import Path
from enum import Enum

from rich.console import Console
from rich.panel import Panel
from rich.syntax import Syntax

from proprep.utils.paths import get_package_dir

logger = logging.getLogger(__name__)


class SimulationType(Enum):
    """AMBER simulation types."""
    MINIMIZATION = "minimization"
    HEATING = "heating" 
    EQUILIBRATION = "equilibration"
    PRODUCTION = "production"


class AmberAnnotatedTemplate:
    """
    Simple AMBER template with inline comment generation.
    Same structure as original, just adds helpful comments to mdin output.
    """
    
    def __init__(self, name: str, description: str, simulation_type: SimulationType, 
                 config: Dict[str, Any], nmr_section: str = ""):
        self.name = name
        self.description = description
        self.simulation_type = simulation_type
        self.config = config
        self.nmr_section = nmr_section
        self.metadata = {}  # Initialize metadata for wizard tracking
        
    def generate_mdin_content(self) -> str:
        """Generate mdin file content with inline comments - like user's example."""
        lines = [
            self.description,
            "&cntrl"
        ]
        
        # Add parameters in organized order with inline comments
        ordered_params = [
            "imin", "ntx", "irest", "nstlim", "dt", "maxcyc", "ncyc", "ntmin",
            "ntpr", "ntwx", "ntwr", "ntxo", "ioutfm",
            "ntt", "temp0", "tempi", "gamma_ln", "ntp", "ntb", "pres0", "barostat",
            "ntc", "ntf", "ntr", "restraint_wt", "restraintmask", "cut", "nmropt"
        ]
        
        for param in ordered_params:
            if param in self.config:
                value = self.config[param]
                comment = self._get_param_comment(param)
                
                if isinstance(value, str) and " " in value:
                    lines.append(f"  {param}='{value}', ! {comment}")
                else:
                    lines.append(f"  {param}={value}, ! {comment}")
        
        # Add any remaining parameters
        for key, value in self.config.items():
            if key not in ordered_params:
                comment = self._get_param_comment(key)
                if isinstance(value, str) and " " in value:
                    lines.append(f"  {key}='{value}', ! {comment}")
                else:
                    lines.append(f"  {key}={value}, ! {comment}")
        
        lines.append("/")
        
        # Add NMR section if present
        if self.nmr_section:
            lines.append("")
            lines.extend(self.nmr_section.split('\n'))
        
        return "\n".join(lines) + "\n"
        
    def _get_param_comment(self, param: str) -> str:
        """Get inline comment for parameter - like user's example."""
        comments = {
            "imin": "0=MD, 1=minimization",
            "ntx": "1=coordinates only, 5=coordinates+velocities", 
            "irest": "0=new simulation, 1=restart",
            "nstlim": "Number of MD steps",
            "dt": "Time step in ps",
            "maxcyc": "Maximum minimization cycles",
            "ncyc": "Steepest descent cycles before conjugate gradient",
            "ntmin": "0=conjugate gradient, 1=steepest descent then conjugate, 2=steepest descent only",
            "ntpr": "Print energy every N steps",
            "ntwx": "Write trajectory every N steps",
            "ntwr": "Write restart file every N steps", 
            "ntxo": "1=formatted, 2=NetCDF restart",
            "ioutfm": "0=formatted, 1=NetCDF trajectory",
            "ntt": "0=NVE, 1=weak coupling, 3=Langevin",
            "temp0": "Target temperature (K)",
            "tempi": "Initial temperature (K)",
            "gamma_ln": "Langevin collision frequency (ps-1)",
            "ntp": "0=no pressure control, 1=isotropic scaling",
            "ntb": "1=constant volume, 2=constant pressure",
            "pres0": "Reference pressure (bar)",
            "barostat": "1=Berendsen, 2=Monte Carlo",
            "ntc": "1=no SHAKE, 2=bonds with H, 3=all bonds",
            "ntf": "1=complete forces, 2=omit bonds with H",
            "ntr": "0=no restraints, 1=positional restraints",
            "restraint_wt": "Restraint force constant (kcal/mol/A^2)",
            "restraintmask": "Atom selection for restraints",
            "cut": "Nonbonded cutoff distance (A)",
            "nmropt": "0=no varying conditions, 1=&wt blocks enabled",
            "dx0": "Initial step size for minimization (A)",
            "drms": "RMS gradient convergence criterion (kcal/mol/A)"
        }
        return comments.get(param, "")
        
    def get_config_dict(self) -> Dict[str, Any]:
        """Get the current configuration dictionary for wizard initialization."""
        return self.config.copy()
        
    def update_from_wizard_config(self, wizard_config: Dict[str, Any]) -> None:
        """Update template config from wizard results."""
        self.config.update(wizard_config)


class AmberAnnotatedTemplateSystem:
    """
    Simple annotated template system - just adds inline comments.
    Same JSON structure as original, just better mdin output.
    """
    
    def __init__(self, package_dir: Path = None, console: Console = None):
        self.package_dir = package_dir or (get_package_dir() / "md_prep")
        self.console = console or Console()

        # Built-in templates ship inside the package (read-only). User templates
        # must live in a user-writable location: the package install dir is
        # frequently read-only (e.g. a shared conda env), so creating or saving
        # under it raises PermissionError. Mirror UserDataManager's ~/.proprep
        # layout for user data.
        self.user_data_dir = self.package_dir / "user_data"
        self.templates_dir = self.user_data_dir / "templates"
        self.builtin_templates_dir = self.templates_dir / "builtin"
        self.user_templates_dir = Path.home() / ".proprep" / "md_annotated_templates" / "user"

        # Ensure directories exist
        self._ensure_directories()

        # Create default templates if needed
        self._create_default_templates()

    def _ensure_directories(self):
        """Ensure the user-writable template directory exists.

        Only the user templates dir is created — it lives under ~/.proprep and is
        always writable. The builtin dir ships inside the (possibly read-only)
        package, so we never mkdir it here.
        """
        self.user_templates_dir.mkdir(parents=True, exist_ok=True)
            
    def _create_default_templates(self):
        """Create default templates - same structure as original."""
        
        templates = {
            "minimization": {
                "name": "Energy Minimization",
                "description": "Energy Minimization Stage in Explicit Solvent",
                "simulation_type": "minimization",
                "config": {
                    "imin": 1,
                    "ntx": 1,
                    "maxcyc": 10000,
                    "ncyc": 1000,
                    "ntmin": 2,
                    "dx0": 0.01,
                    "drms": 0.0001,
                    "cut": 10.0,
                    "ntpr": 100,
                    "ntwr": 500,
                    "ntr": 1,
                    "restraint_wt": 500.0,
                    "restraintmask": "!@H="
                }
            },
            "heating": {
                "name": "Heating Protocol",
                "description": "Heating Stage from 0 to 300K",
                "simulation_type": "heating",
                "config": {
                    "imin": 0,
                    "ntx": 1,
                    "irest": 0,
                    "nstlim": 50000,
                    "dt": 0.001,
                    "ntt": 3,
                    "tempi": 0.0,
                    "temp0": 300.0,
                    "gamma_ln": 2.0,
                    "ntb": 1,
                    "ntc": 2,
                    "ntf": 2,
                    "cut": 10.0,
                    "nmropt": 1,
                    "ntr": 1,
                    "restraint_wt": 10.0,
                    "restraintmask": "!@H=",
                    "ntpr": 500,
                    "ntwx": 500,
                    "ntwr": 2500,
                    "ntxo": 2,
                    "ioutfm": 1
                },
                "nmr_section": "&wt type='TEMP0', istep1=0, istep2=25000,\n    value1=0.0, value2=300.0 /\n&wt type='END' /"
            },
            "equilibration": {
                "name": "NPT Equilibration",
                "description": "Pressure Equilibration Stage at 300K",
                "simulation_type": "equilibration",
                "config": {
                    "imin": 0,
                    "ntx": 5,
                    "irest": 1,
                    "nstlim": 100000,
                    "dt": 0.002,
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
                    "ntr": 1,
                    "restraint_wt": 1.0,
                    "restraintmask": "!@H=",
                    "ntpr": 1000,
                    "ntwx": 1000,
                    "ntwr": 5000,
                    "ntxo": 2,
                    "ioutfm": 1
                }
            },
            "production": {
                "name": "Production MD",
                "description": "Production Molecular Dynamics",
                "simulation_type": "production",
                "config": {
                    "imin": 0,
                    "ntx": 5,
                    "irest": 1,
                    "nstlim": 2500000,
                    "dt": 0.002,
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
                    "ntpr": 1000,
                    "ntwx": 1000,
                    "ntwr": 10000,
                    "ntxo": 2,
                    "ioutfm": 1
                }
            }
        }
        
        # Save templates in original JSON format. Builtin templates ship with the
        # package, so this only fills gaps. On a read-only install the package dir
        # can't be written — that's fine (the shipped files are already present),
        # so tolerate PermissionError/OSError instead of crashing.
        for template_name, template_data in templates.items():
            template_file = self.builtin_templates_dir / f"{template_name}.json"
            if not template_file.exists():
                try:
                    self.builtin_templates_dir.mkdir(parents=True, exist_ok=True)
                    with open(template_file, 'w') as f:
                        json.dump(template_data, f, indent=2)
                except OSError as e:
                    logger.debug(
                        "Skipping builtin template regeneration for %s "
                        "(package dir not writable): %s", template_name, e)
                    
    def load_template(self, simulation_type: SimulationType) -> AmberAnnotatedTemplate:
        """Load template from JSON file - same format as original."""
        # Try user templates first, then builtin
        for directory in [self.user_templates_dir, self.builtin_templates_dir]:
            filepath = directory / f"{simulation_type.value}.json"
            if filepath.exists():
                with open(filepath, 'r') as f:
                    data = json.load(f)
                    
                return AmberAnnotatedTemplate(
                    name=data["name"],
                    description=data["description"],
                    simulation_type=SimulationType(data["simulation_type"]),
                    config=data["config"],
                    nmr_section=data.get("nmr_section", "")
                )
                
        raise FileNotFoundError(f"Template for {simulation_type.value} not found")
        
    def save_template(self, template: AmberAnnotatedTemplate, directory: Path = None):
        """Save template to JSON file."""
        if directory is None:
            directory = self.user_templates_dir
            
        template_data = {
            "name": template.name,
            "description": template.description,
            "simulation_type": template.simulation_type.value,
            "config": template.config,
            "nmr_section": template.nmr_section
        }
        
        filepath = directory / f"{template.simulation_type.value}.json"
        with open(filepath, 'w') as f:
            json.dump(template_data, f, indent=2)
        
    def display_template_overview(self, template: AmberAnnotatedTemplate) -> None:
        """Show simple template overview - no tables, just clean info."""
        # This method is deprecated - use display_template_content instead
        pass
        
    def display_template_content(self, template: AmberAnnotatedTemplate) -> None:
        """Display the clean template content with syntax highlighting."""
        content = template.generate_mdin_content()
        
        # Use tLEaP-style display: title + dashed lines + syntax + dashed lines
        self.console.print(f"\n[bold]Current {template.name}:[/bold]")
        self.console.print("[grey50]" + "─" * 60 + "[/grey50]")
        
        syntax = Syntax(content, "fortran", theme="monokai", line_numbers=True)
        self.console.print(syntax)
        
        self.console.print("[grey50]" + "─" * 60 + "[/grey50]")