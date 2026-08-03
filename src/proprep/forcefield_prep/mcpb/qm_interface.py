"""
QM Interface for MCPB Step 2
Manages QM calculation workflow for bonded parameter generation.

© 2024 ProPrep Developer. All rights reserved.

This module handles:
- QM input file generation from PDB coordinates
- QM software support (Gaussian; ORCA under consideration)
- Guided and external execution modes
- QM output validation

PROPRIETARY SOFTWARE: This file contains proprietary code belonging to the ProPrep project.
Unauthorized copying, distribution, or modification is strictly prohibited.
"""

from enum import Enum
from pathlib import Path
from typing import Optional, Dict, Any, Tuple
from rich.console import Console
from rich.panel import Panel


class QMSoftware(Enum):
    """Supported QM software packages."""
    GAUSSIAN = "gaussian"
    GAMESS = "gamess"
    ORCA = "orca"


class QMCalculationMode(Enum):
    """QM calculation execution modes."""
    EXTERNAL = "external"  # User provides completed output files
    GUIDED = "guided"      # ProPrep generates inputs, user executes manually


class QMInterface:
    """
    Manages QM calculation workflow for force constant generation.

    Supports two modes:
    1. GUIDED: ProPrep writes input files, user runs calculations externally
    2. EXTERNAL: User provides completed QM output files

    QM calculations required:
    - Geometry optimization (optional, if structure not pre-optimized)
    - Frequency calculation to obtain Hessian matrix
    """

    def __init__(self, software: QMSoftware, mode: QMCalculationMode,
                 console: Optional[Console] = None):
        """
        Initialize QM interface.

        Args:
            software: QM software to use (Gaussian, GAMESS, ORCA)
            mode: Execution mode (GUIDED or EXTERNAL)
            console: Rich console for output
        """
        self.software = software
        self.mode = mode
        self.console = console or Console()

    def generate_input_files(self, pdb_file: Path, output_dir: Path,
                            charge: int, multiplicity: int,
                            memory_gb: int = 4, n_processors: int = 4,
                            functional: str = "B3LYP", basis_set: str = "6-31G*",
                            basis_groups: list = None, ecp_specs: list = None,
                            job_type: str = "Opt Freq",
                            additional_keywords: str = "",
                            title_card: str = "MCPB Step 2 - Small model",
                            output_name: str = None,
                            route_line: str = None,
                            include_seminario_iop: bool = True,
                            metal_radii: Optional[Dict[str, float]] = None,
                            frozen_atoms: Optional[set] = None) -> Dict[str, Path]:
        """
        Generate QM input files from PDB coordinates.

        Creates input files for QM software with appropriate settings for
        MCPB parameter generation.

        Args:
            pdb_file: PDB file with small model coordinates
            output_dir: Directory for QM input files
            charge: Total charge of system
            multiplicity: Spin multiplicity (2S+1)
            memory_gb: Memory allocation in GB
            n_processors: Number of processors
            functional: DFT functional (e.g., B3LYP)
            basis_set: Basis set name (e.g., 6-31G*, GenECP)
            basis_groups: List of (atoms_list, basis_name) tuples for Gen/GenECP
            ecp_specs: List of (atoms_list, ecp_name) tuples for GenECP
            job_type: Job specification (e.g., "Opt Freq", "Freq", "SP")
            additional_keywords: Additional Gaussian keywords
            title_card: Title card text
            output_name: Custom output file base name (default: {pdb_stem}_freq)
            route_line: Full route section override (if provided, ignores functional/basis_set/job_type)
            include_seminario_iop: Include IOp(7/33=1) for Seminario method (default True)
            metal_radii: Dict of {element: radius} for MK ESP with ReadRadii (e.g., {'Zn': 1.395})
            frozen_atoms: Set of 0-based atom indices to freeze during geometry optimization

        Returns:
            Dict mapping file types to paths:
                - "input_file": Main input file
                - "checkpoint_file": Checkpoint file path (Gaussian)
                - "log_file": Expected output log file path
                - "fchk_file": Expected formatted checkpoint (Gaussian)
        """
        if self.software == QMSoftware.GAUSSIAN:
            return self._generate_gaussian_input(
                pdb_file, output_dir, charge, multiplicity,
                memory_gb, n_processors, functional, basis_set,
                basis_groups, ecp_specs, job_type, additional_keywords, title_card,
                output_name, route_line, include_seminario_iop, metal_radii, frozen_atoms
            )
        elif self.software == QMSoftware.GAMESS:
            return self._generate_gamess_input(pdb_file, output_dir, charge, multiplicity)
        elif self.software == QMSoftware.ORCA:
            return self._generate_orca_input(pdb_file, output_dir, charge, multiplicity)
        else:
            raise NotImplementedError(f"{self.software.value} not yet supported")

    def _parse_pdb_coordinates(self, pdb_file: Path) -> list:
        """
        Parse coordinates from PDB file.

        Args:
            pdb_file: Path to PDB file

        Returns:
            List of (element, x, y, z) tuples
        """
        # Known metal elements for fallback detection
        from proprep.forcefield_prep.structure_preprocessor import METAL_ELEMENTS

        coords = []
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # PDB format: columns 13-16 = atom name, 17-20 = resname,
                    # 31-38 = x, 39-46 = y, 47-54 = z, 77-78 = element
                    element = line[76:78].strip()
                    atom_name = line[12:16].strip()
                    resname = line[17:20].strip()

                    if not element:
                        # Fallback: extract from atom name
                        element = ''.join([c for c in atom_name if c.isalpha()])[:2]

                    # BioPython can mis-parse metal elements (e.g., ZN → N)
                    # Only apply metal fallback for single-atom residues where
                    # the residue name matches a metal (e.g., ZN residue with element "N")
                    resname_title = resname[0].upper() + resname[1:].lower() if len(resname) > 1 else resname.upper()
                    if resname_title in METAL_ELEMENTS and resname.upper() == atom_name.strip().upper():
                        element = resname_title

                    # Capitalize properly (e.g., 'ZN' -> 'Zn', 'H' -> 'H')
                    if len(element) == 2:
                        element = element[0].upper() + element[1].lower()
                    else:
                        element = element.upper()

                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append((element, x, y, z))
        return coords

    def _generate_gaussian_input(self, pdb_file: Path, output_dir: Path,
                                 charge: int, multiplicity: int,
                                 memory_gb: int = 4, n_processors: int = 4,
                                 functional: str = "B3LYP", basis_set: str = "6-31G*",
                                 basis_groups: list = None, ecp_specs: list = None,
                                 job_type: str = "Opt Freq",
                                 additional_keywords: str = "",
                                 title_card: str = "MCPB Step 2 - Small model",
                                 output_name: str = None,
                                 route_line: str = None,
                                 include_seminario_iop: bool = True,
                                 metal_radii: Optional[Dict[str, float]] = None,
                                 frozen_atoms: Optional[set] = None) -> Dict[str, Path]:
        """
        Generate Gaussian input template with customizable settings.

        Args:
            output_name: Custom base name for output files (default: {pdb_stem}_freq)
            route_line: Full route section override (if provided, ignores functional/basis_set/job_type)
            include_seminario_iop: Include IOp(7/33=1) for Seminario method (default True)
            metal_radii: Dict of {element: radius} for MK ESP with ReadRadii
            frozen_atoms: Set of 0-based atom indices to freeze (written with -1 flag)

        Supports Gen/GenECP for mixed basis sets.
        USER MUST REVIEW AND EDIT before running!
        """
        # Determine output file naming
        base_name = output_name if output_name else f"{pdb_file.stem}_freq"
        gjf_file = output_dir / f"{base_name}.gjf"

        # Parse PDB coordinates
        coords = self._parse_pdb_coordinates(pdb_file)

        # Write Gaussian input template
        with open(gjf_file, 'w') as f:
            # Link 0 commands
            f.write(f"%Chk={base_name}.chk\n")
            f.write(f"%Mem={memory_gb}GB\n")
            f.write(f"%NProcShared={n_processors}\n")

            # Route section
            if route_line:
                # Full route override provided - use as-is
                final_route = route_line if route_line.startswith("#") else f"# {route_line}"
            else:
                # Build route from components
                final_route = f"# {functional}/{basis_set} {job_type}"
                # Standard quality keywords
                final_route += " Geom=PrintInputOrient Integral=(Grid=UltraFine)"
                # Seminario method requires IOp(7/33=1) for Cartesian force constants
                if include_seminario_iop:
                    final_route += " IOp(7/33=1)"
                # Add user's additional keywords
                if additional_keywords:
                    final_route += f" {additional_keywords}"

            f.write(f"{final_route}\n")
            f.write("\n")

            # Title card
            f.write(f"{title_card}\n")
            f.write("\n")

            # Charge and multiplicity
            f.write(f"{charge} {multiplicity}\n")

            # Write coordinates (with freeze flag: 0=free, -1=frozen)
            for i, (element, x, y, z) in enumerate(coords):
                freeze_flag = -1 if frozen_atoms and i in frozen_atoms else 0
                f.write(f"{element:2s}  {freeze_flag:2d}  {x:12.8f}  {y:12.8f}  {z:12.8f}\n")

            f.write("\n")

            # Handle Gen/GenECP basis set specifications
            if basis_set.lower() in ["gen", "genecp"] and basis_groups:
                # Write basis set specifications
                for atoms, basis_name in basis_groups:
                    # Write atom list
                    atom_line = " ".join(atoms) + " 0"
                    f.write(f"{atom_line}\n")
                    # Write basis set name
                    f.write(f"{basis_name}\n")
                    f.write("****\n")

                f.write("\n")

                # Write ECP specifications if GenECP
                if basis_set.lower() == "genecp" and ecp_specs:
                    for ecp_atoms, ecp_basis in ecp_specs:
                        # Write atom list
                        ecp_atom_line = " ".join(ecp_atoms) + " 0"
                        f.write(f"{ecp_atom_line}\n")
                        # Write ECP basis set name
                        f.write(f"{ecp_basis}\n")

                    f.write("\n")

            # Write metal radii for MK ESP with ReadRadii
            # Gaussian requires radii written TWICE with blank line separators
            # Format: radii_block, blank, radii_block, blank, blank
            if metal_radii:
                for element, radius in metal_radii.items():
                    f.write(f"{element} {radius:.3f}\n")
                f.write("\n")
                for element, radius in metal_radii.items():
                    f.write(f"{element} {radius:.3f}\n")
                f.write("\n")

            # Final blank line
            f.write("\n")

        self.console.print(f"[green]✅ Generated Gaussian input: {gjf_file.name}[/green]")

        # Show execution instructions
        self._show_gaussian_execution_guide(gjf_file)

        return {
            "input_file": gjf_file,
            "checkpoint_file": output_dir / f"{base_name}.chk",
            "log_file": output_dir / f"{base_name}.log",
            "fchk_file": output_dir / f"{base_name}.fchk"
        }

    def _generate_gamess_input(self, pdb_file: Path, output_dir: Path,
                               charge: int, multiplicity: int) -> Dict[str, Path]:
        """Generate GAMESS-US input template for Hessian calculation."""
        inp_file = output_dir / f"{pdb_file.stem}_freq.inp"

        # Parse PDB coordinates
        coords = self._parse_pdb_coordinates(pdb_file)

        # Convert to GAMESS format
        with open(inp_file, 'w') as f:
            f.write(f"! MCPB Step 2 - Small model\n")
            f.write(f" $CONTRL SCFTYP=RHF RUNTYP=HESSIAN\n")
            f.write(f"  ICHARG={charge} MULT={multiplicity} $END\n")
            f.write(f" $SYSTEM MWORDS=500 $END\n")
            f.write(f" $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n")
            f.write(f" $STATPT OPTTOL=0.0001 NSTEP=200 $END\n")
            f.write(f" $DATA\n")
            f.write(f"Small model\n")
            f.write(f"C1\n")

            # Write atomic coordinates
            for element, x, y, z in coords:
                f.write(f"{element:2s}  {self._get_atomic_number(element):2.1f}  "
                       f"{x:12.8f}  {y:12.8f}  {z:12.8f}\n")

            f.write(f" $END\n")

        self.console.print(f"[green]✅ Generated GAMESS input: {inp_file.name}[/green]")

        # Show execution instructions
        self._show_gamess_execution_guide(inp_file)

        return {
            "input_file": inp_file,
            "log_file": output_dir / f"{pdb_file.stem}_freq.log",
            "dat_file": output_dir / f"{pdb_file.stem}_freq.dat"
        }

    def _generate_orca_input(self, pdb_file: Path, output_dir: Path,
                            charge: int, multiplicity: int) -> Dict[str, Path]:
        """Generate ORCA input template for optimization + frequency calculation."""
        inp_file = output_dir / f"{pdb_file.stem}_freq.inp"

        # Parse PDB coordinates
        coords = self._parse_pdb_coordinates(pdb_file)

        with open(inp_file, 'w') as f:
            f.write(f"# MCPB Step 2 - Small model\n")
            f.write(f"! B3LYP 6-31G* TightSCF Grid5 GridX5\n")
            # Always optimization + frequency for bonded parameters
            f.write(f"! Opt Freq\n")
            f.write(f"\n")
            f.write(f"%pal nprocs 4 end\n")
            f.write(f"%maxcore 1000\n")
            f.write(f"\n")
            f.write(f"* xyz {charge} {multiplicity}\n")

            for element, x, y, z in coords:
                f.write(f"{element:2s}  {x:12.8f}  {y:12.8f}  {z:12.8f}\n")

            f.write(f"*\n")

        self.console.print(f"[green]✅ Generated ORCA input: {inp_file.name}[/green]")

        # Show execution instructions
        self._show_orca_execution_guide(inp_file)

        return {
            "input_file": inp_file,
            "log_file": output_dir / f"{pdb_file.stem}_freq.out",
            "hess_file": output_dir / f"{pdb_file.stem}_freq.hess"
        }

    def _show_gaussian_execution_guide(self, input_file: Path):
        """Display execution instructions for Gaussian."""
        output_dir = input_file.parent
        guide_text = (
            f"Open a separate terminal session to inspect and run:\n\n"
            f"1. Navigate to directory:\n"
            f"   [cyan]cd {output_dir}[/cyan]\n\n"
            f"2. Review and edit (if needed):\n"
            f"   {input_file.absolute()}\n\n"
            f"3. Run the calculation from within that directory:\n"
            f"   [cyan]g16 {input_file.name}[/cyan]\n\n"
            f"4. Convert checkpoint file:\n"
            f"   [cyan]formchk {input_file.stem}.chk {input_file.stem}.fchk[/cyan]\n"
        )
        self.console.print(Panel(guide_text, title="Next Steps", border_style="yellow", expand=False))

    def _show_gamess_execution_guide(self, input_file: Path):
        """Display execution instructions for GAMESS."""
        guide_text = (
            f"[bold yellow]Input Template Location:[/bold yellow]\n"
            f"  {input_file.absolute()}\n\n"
            f"[bold yellow]IMPORTANT: Review and edit this file before running![/bold yellow]\n"
            f"  • Adjust basis set for metals\n"
            f"  • Modify memory settings for your system\n"
            f"  • Review calculation settings\n\n"
            f"[bold]After editing, run:[/bold]\n"
            f"  [cyan]rungms {input_file.name} > {input_file.stem}.log[/cyan]\n"
        )
        self.console.print(Panel(guide_text, title="Template Generated", border_style="yellow"))

    def _show_orca_execution_guide(self, input_file: Path):
        """Display execution instructions for ORCA."""
        guide_text = (
            f"[bold yellow]Input Template Location:[/bold yellow]\n"
            f"  {input_file.absolute()}\n\n"
            f"[bold yellow]IMPORTANT: Review and edit this file before running![/bold yellow]\n"
            f"  • Adjust basis set for metals (consider ECP for heavy metals)\n"
            f"  • Modify processors/memory for your system\n"
            f"  • Add solvent model if needed\n"
            f"  • Review method/basis set appropriateness\n\n"
            f"[bold]After editing, run:[/bold]\n"
            f"  [cyan]orca {input_file.name} > {input_file.stem}.out[/cyan]\n"
        )
        self.console.print(Panel(guide_text, title="Template Generated", border_style="yellow"))

    def _get_atomic_number(self, element: str) -> float:
        """Get atomic number for element symbol."""
        atomic_numbers = {
            'H': 1.0, 'C': 6.0, 'N': 7.0, 'O': 8.0, 'S': 16.0,
            'P': 15.0, 'F': 9.0, 'Cl': 17.0, 'Br': 35.0,
            'Fe': 26.0, 'Cu': 29.0, 'Zn': 30.0, 'Mg': 12.0,
            'Ca': 20.0, 'Mn': 25.0, 'Co': 27.0, 'Ni': 28.0
        }
        return atomic_numbers.get(element, 0.0)

    def check_calculation_complete(self, log_file: Path) -> Tuple[bool, str]:
        """
        Check if QM calculation completed successfully.

        Args:
            log_file: Path to QM output log file

        Returns:
            (success, message) tuple
        """
        if not log_file.exists():
            return False, f"Log file not found: {log_file}"

        try:
            with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()

            if self.software == QMSoftware.GAUSSIAN:
                if "Normal termination" in content:
                    if "Frequencies --" in content:
                        return True, "Gaussian calculation completed successfully"
                    else:
                        return False, "Calculation completed but no frequencies found"
                elif "Error termination" in content:
                    return False, "Gaussian calculation terminated with error"
                else:
                    return False, "Gaussian calculation did not complete"

            elif self.software == QMSoftware.GAMESS:
                if "EXECUTION OF GAMESS TERMINATED NORMALLY" in content:
                    if "HESSIAN MATRIX" in content:
                        return True, "GAMESS calculation completed successfully"
                    else:
                        return False, "Calculation completed but no Hessian found"
                else:
                    return False, "GAMESS calculation did not complete normally"

            elif self.software == QMSoftware.ORCA:
                if "ORCA TERMINATED NORMALLY" in content:
                    if "VIBRATIONAL FREQUENCIES" in content:
                        return True, "ORCA calculation completed successfully"
                    else:
                        return False, "Calculation completed but no frequencies found"
                else:
                    return False, "ORCA calculation did not complete normally"

            return False, "Could not determine calculation status"

        except Exception as e:
            return False, f"Error checking log file: {str(e)}"

    def generate_fchk(self, chk_file: Path) -> Optional[Path]:
        """
        Generate formatted checkpoint file from binary checkpoint.
        Required for Hessian extraction from Gaussian calculations.

        Args:
            chk_file: Path to binary checkpoint file

        Returns:
            Path to fchk file if successful, None otherwise
        """
        if self.software != QMSoftware.GAUSSIAN:
            return None

        if not chk_file.exists():
            self.console.print(f"[red]❌ Checkpoint file not found: {chk_file}[/red]")
            return None

        fchk_file = chk_file.with_suffix('.fchk')

        import subprocess
        cmd = ["formchk", str(chk_file), str(fchk_file)]

        try:
            self.console.print(f"[cyan]Running formchk to generate {fchk_file.name}...[/cyan]")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

            if result.returncode == 0 and fchk_file.exists():
                self.console.print(f"[green]✅ Generated fchk file: {fchk_file.name}[/green]")
                return fchk_file
            else:
                self.console.print(f"[red]❌ formchk failed:[/red]")
                if result.stderr:
                    self.console.print(f"[red]{result.stderr}[/red]")
                return None

        except FileNotFoundError:
            self.console.print("[red]❌ formchk command not found in PATH[/red]")
            self.console.print("[yellow]Please ensure Gaussian is properly installed[/yellow]")
            return None

        except subprocess.TimeoutExpired:
            self.console.print("[red]❌ formchk timed out[/red]")
            return None

        except Exception as e:
            self.console.print(f"[red]❌ Error running formchk: {str(e)}[/red]")
            return None

    def get_required_files(self) -> Dict[str, str]:
        """
        Get list of required output files for this QM software.

        Returns:
            Dict mapping file type to description
        """
        if self.software == QMSoftware.GAUSSIAN:
            return {
                "log_file": "Gaussian output log (.log)",
                "chk_file": "Binary checkpoint file (.chk)",
                "fchk_file": "Formatted checkpoint with Hessian (.fchk)"
            }
        elif self.software == QMSoftware.GAMESS:
            return {
                "log_file": "GAMESS output log with Hessian (.log)"
            }
        elif self.software == QMSoftware.ORCA:
            return {
                "log_file": "ORCA output file (.out)",
                "hess_file": "ORCA Hessian file (.hess)"
            }
        return {}
