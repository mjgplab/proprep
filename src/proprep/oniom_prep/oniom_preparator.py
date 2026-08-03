"""
ONIOM Preparator - Main Orchestrator for ONIOM Setup

High-level interface for creating complete ONIOM QM/MM input files from RedoxSite objects.
Coordinates layer assignment, link atom placement, validation, and file writing.

© 2024 ProPrep Developer. All rights reserved.
"""

import logging
from typing import Dict, List, Optional, Tuple
from datetime import datetime
from pathlib import Path
from rich.console import Console

from proprep.structure_prep.comprehensive_redox_detector import RedoxSite
from proprep.forcefield_prep.mcpb.ff_parameter_reader import FFParameterReader

from .data_structures import ONIOMSetup
from .layer_classifier import LayerClassifier
from .link_atom_placer import LinkAtomPlacer
from .oniom_writer import ONIOMWriter
from .oniom_validator import ONIOMValidator, ValidationResult


logger = logging.getLogger(__name__)


class ONIOMPreparator:
    """
    Main class for preparing ONIOM QM/MM calculations from RedoxSite objects.

    Workflow:
    1. Initialize with RedoxSite and configuration
    2. Assign atoms to ONIOM layers (HIGH/MEDIUM/LOW)
    3. Place link atoms at QM/MM boundaries
    4. Load atom types and charges (from MCPB or defaults)
    5. Load force field parameters
    6. Validate setup for errors
    7. Write Gaussian ONIOM input file

    This class serves as the high-level interface that coordinates all
    sub-modules while showing the user what is happening at each step.
    """

    def __init__(
        self,
        redox_site: RedoxSite,
        logger_instance=None,
        console=None
    ):
        """
        Initialize ONIOM preparator.

        Args:
            redox_site: RedoxSite containing structure and bonding
            logger_instance: Optional logger
            console: Optional Rich Console for output
        """
        self.redox_site = redox_site
        self.logger = logger_instance or logger
        self.console = console or Console()

        # ONIOM setup object (populated during preparation)
        self.oniom_setup: Optional[ONIOMSetup] = None

        # Sub-modules
        self.layer_classifier: Optional[LayerClassifier] = None
        self.link_atom_placer: Optional[LinkAtomPlacer] = None
        self.ff_reader: Optional[FFParameterReader] = None
        self.validator: Optional[ONIOMValidator] = None
        self.writer: Optional[ONIOMWriter] = None

    def prepare_oniom_setup(
        self,
        # Layer selection (parmed residue indices)
        high_residues: Optional[List[int]] = None,
        medium_residues: Optional[List[int]] = None,
        n_layers: int = 2,

        # Freezing options
        freeze_distance_cutoff: Optional[float] = None,
        frozen_residues: Optional[List[int]] = None,

        # Link atom options
        custom_scale_factors: Optional[Dict[str, float]] = None,

        # QM settings
        qm_functional: str = "B3LYP",
        qm_basis_set: str = "6-31G*",
        qm_charge: int = 0,
        qm_multiplicity: int = 1,

        # MM settings
        mm_forcefield: str = "AMBER",
        force_field_name: str = "ff19SB",
        standard_dat_files: Optional[List[str]] = None,
        standard_frcmod_files: Optional[List[str]] = None,
        custom_frcmod_files: Optional[List[str]] = None,
        amberhome: Optional[str] = None,
        prmtop_path: Optional[str] = None,

        # Medium layer settings (3-layer only)
        medium_method: str = "HF/3-21G",
        medium_charge: int = 0,
        medium_multiplicity: int = 1,

        # Job settings
        job_type: str = "Opt",
        n_processors: int = 4,
        memory_gb: int = 8,
        additional_keywords: str = "",

        # Atom types and charges (keyed by parmed atom index)
        atom_types: Optional[Dict[int, str]] = None,
        charges: Optional[Dict[int, float]] = None,

        # Pre-computed layer assignments (keyed by parmed atom index)
        layer_assignments: Optional[Dict[int, 'LayerAssignment']] = None,

        # Boundary atoms (parmed atom indices)
        boundary_atom_indices: Optional[List[int]] = None,

        # parmed structure
        parm = None,
    ) -> ONIOMSetup:
        """
        Prepare complete ONIOM setup from configuration.

        This is the main workflow function that coordinates all steps,
        showing the user what is happening and why at each stage.

        Args:
            high_residues: List of (chain, resid, insertion_code) for HIGH layer
            medium_residues: List for MEDIUM layer (3-layer only)
            n_layers: 2 or 3
            freeze_distance_cutoff: Distance (Å) beyond which atoms frozen
            frozen_residues: Explicitly frozen residues
            custom_scale_factors: Link atom scale factors (e.g., {"C-C": 0.723})
            qm_functional: QM functional (e.g., "B3LYP", "M06-2X")
            qm_basis_set: QM basis set (e.g., "6-31G*", "6-311G(d,p)")
            qm_charge: Charge of QM region
            qm_multiplicity: Spin multiplicity of QM region
            mm_forcefield: MM forcefield keyword (usually "AMBER")
            force_field_name: AMBER forcefield name (e.g., "ff19SB") or custom leaprc names
            custom_frcmod_files: Optional list of custom .frcmod files from forcefield_params directory
            amberhome: Path to AMBER installation
            medium_method: Method for medium layer (3-layer only)
            medium_charge: Charge of medium layer
            medium_multiplicity: Multiplicity of medium layer
            job_type: Gaussian job type ("Opt", "SP", "Freq", etc.)
            n_processors: Number of processors
            memory_gb: Memory allocation (GB)
            additional_keywords: Additional Gaussian keywords
            atom_types: Atom type assignments (from MCPB)
            charges: Atomic charges (from MCPB)
            layer_assignments: Pre-computed layer assignments (bypasses layer classification)

        Returns:
            Complete ONIOMSetup object

        Raises:
            ValueError: If configuration is invalid
        """
        self.logger.info("=" * 70)
        self.logger.info("Starting ONIOM QM/MM Preparation")
        self.logger.info("=" * 70)

        # ========================================
        # STEP 1: Create ONIOMSetup object
        # ========================================
        self.console.print("\n[bold cyan]Initializing ONIOM Configuration[/bold cyan]")
        self.console.print(f"  QM Method: {qm_functional}/{qm_basis_set}")
        self.console.print(f"  QM Charge: {qm_charge}, Multiplicity: {qm_multiplicity}")
        self.console.print(f"  MM Forcefield: {mm_forcefield}")
        self.console.print(f"  Number of layers: {n_layers}")
        self.console.print(f"  Job type: {job_type}")

        self.oniom_setup = ONIOMSetup(
            redox_site=self.redox_site,
            parm=parm,
            qm_functional=qm_functional,
            qm_basis_set=qm_basis_set,
            qm_charge=qm_charge,
            qm_multiplicity=qm_multiplicity,
            mm_forcefield=mm_forcefield,
            medium_method=medium_method,
            medium_charge=medium_charge,
            medium_multiplicity=medium_multiplicity,
            job_type=job_type,
            n_processors=n_processors,
            memory_gb=memory_gb,
            additional_keywords=additional_keywords,
            n_layers=n_layers,
            setup_timestamp=datetime.now().isoformat()
        )

        self.logger.info(f"Created ONIOM configuration: {n_layers}-layer, QM={qm_functional}/{qm_basis_set}")

        # ========================================
        # STEP 2: Assign layers
        # ========================================
        self.console.print("\n[bold cyan]Step 1/6: Assigning ONIOM Layers[/bold cyan]")

        if layer_assignments is not None:
            # Use pre-computed layer assignments (from MM configuration)
            self.console.print("  Using pre-computed layer assignments from MM configuration")
            self.console.print("  (HIGH region already expanded to include C=O and N-H from flanking residues)", highlight=False)
            self.oniom_setup.layer_assignments = layer_assignments

            # Apply distance-based freezing to pre-computed assignments.
            # An atom stays active if it is within freeze_distance_cutoff of
            # ANY HIGH atom (nearest-atom shell), not from a single centroid.
            if freeze_distance_cutoff is not None and parm is not None:
                from .data_structures import ONIOMLayer, FreezeFlag

                high_coords = [
                    (parm.atoms[idx].xx, parm.atoms[idx].xy, parm.atoms[idx].xz)
                    for idx, a in layer_assignments.items()
                    if a.layer == ONIOMLayer.HIGH
                ]
                if high_coords:
                    cutoff_sq = freeze_distance_cutoff * freeze_distance_cutoff
                    n_frozen = 0
                    for idx, assignment in layer_assignments.items():
                        if assignment.layer in (ONIOMLayer.HIGH, ONIOMLayer.MEDIUM):
                            continue
                        atom = parm.atoms[idx]
                        ax, ay, az = atom.xx, atom.xy, atom.xz
                        # Frozen unless within cutoff of at least one HIGH atom
                        within_shell = False
                        for hx, hy, hz in high_coords:
                            dx = ax - hx
                            dy = ay - hy
                            dz = az - hz
                            if dx * dx + dy * dy + dz * dz <= cutoff_sq:
                                within_shell = True
                                break
                        if not within_shell:
                            assignment.freeze = FreezeFlag.FROZEN
                            n_frozen += 1
                    self.console.print(
                        f"  Froze {n_frozen} atoms beyond {freeze_distance_cutoff}Å from any HIGH atom",
                        highlight=False
                    )
        else:
            # Compute layer assignments from scratch
            self.console.print(f"  HIGH layer residues: {len(high_residues)}")
            if medium_residues:
                self.console.print(f"  MEDIUM layer residues: {len(medium_residues)}")
            self.console.print("  LOW layer: all remaining atoms")
            self.console.print("\n  Applying alpha carbon bridging to prevent peptide bond cutting...")
            self.console.print("  Merging nearby fragments to create contiguous QM regions...")

            self.layer_classifier = LayerClassifier(
                redox_site=self.redox_site,
                logger_instance=self.logger,
                console=self.console
            )

            layer_assignments = self.layer_classifier.assign_layers(
                high_residues=high_residues,
                medium_residues=medium_residues,
                freeze_distance_cutoff=freeze_distance_cutoff,
                frozen_residues=frozen_residues,
                n_layers=n_layers
            )

            self.oniom_setup.layer_assignments = layer_assignments

        # Count atoms per layer
        n_high = sum(1 for a in self.oniom_setup.layer_assignments.values() if a.layer.value == "H")
        n_medium = sum(1 for a in self.oniom_setup.layer_assignments.values() if a.layer.value == "M")
        n_low = sum(1 for a in self.oniom_setup.layer_assignments.values() if a.layer.value == "L")

        self.console.print(f"\n  [green]✓[/green] Layer assignment complete:")
        self.console.print(f"    HIGH: {n_high} atoms")
        if n_medium > 0:
            self.console.print(f"    MEDIUM: {n_medium} atoms")
        self.console.print(f"    LOW: {n_low} atoms")

        # ========================================
        # STEP 3: Place link atoms
        # ========================================
        self.console.print("\n[bold cyan]Step 2/6: Placing Link Atoms at QM/MM Boundaries[/bold cyan]")
        self.console.print("  Link atoms are hydrogen atoms placed along bonds crossing the QM/MM boundary.")
        self.console.print("  They maintain proper valence at the boundary while allowing separate treatment")
        self.console.print("  of QM and MM regions.")

        self.link_atom_placer = LinkAtomPlacer(
            redox_site=self.redox_site,
            layer_assignments=self.oniom_setup.layer_assignments,
            boundary_atom_indices=boundary_atom_indices,
            parm=parm,
            logger_instance=self.logger,
            console=self.console
        )

        link_atoms = self.link_atom_placer.place_link_atoms(
            custom_scale_factors=custom_scale_factors
        )

        self.oniom_setup.link_atoms = link_atoms

        if link_atoms:
            self.console.print(f"\n  [green]✓[/green] Placed {len(link_atoms)} link atoms")
            # Classify by boundary type
            boundary_counts = {}
            for la in link_atoms:
                boundary_counts[la.boundary_type] = boundary_counts.get(la.boundary_type, 0) + 1
            for btype, count in sorted(boundary_counts.items()):
                self.console.print(f"    {btype}: {count}")
        else:
            self.console.print("  [yellow]⚠[/yellow] No link atoms placed (no QM/MM boundary bonds found)")

        # ========================================
        # STEP 3.5: Build connectivity table (NEW!)
        # ========================================
        if parm is not None:
            self.console.print("\n[bold cyan]Building Connectivity Table[/bold cyan]")
            self.console.print("  Reading bonds from prmtop topology.")

            try:
                from proprep.oniom_prep.connectivity_builder import ConnectivityBuilder

                conn_builder = ConnectivityBuilder(
                    oniom_setup=self.oniom_setup,
                    parm=parm,
                    logger_instance=self.logger
                )

                connectivity_table = conn_builder.build_connectivity()
                self.oniom_setup.connectivity = connectivity_table
                self.oniom_setup.use_explicit_connectivity = True

                n_bonds = len(conn_builder.bonds)
                n_link = conn_builder.stats.get('link_bonds', 0)
                self.console.print(f"\n  [green]✓[/green] Built connectivity table: {n_bonds} bonds")
                if n_link:
                    self.console.print(f"    {n_link} link atom bonds added")

            except Exception as e:
                self.console.print(f"\n  [yellow]⚠[/yellow] Connectivity building failed: {e}")
                self.console.print("    Gaussian will infer connectivity from inter-atomic distances")
                self.logger.warning(f"Connectivity building failed: {e}", exc_info=True)
                self.oniom_setup.use_explicit_connectivity = False
        else:
            self.console.print("\n[grey50]Skipping connectivity building (no prmtop available)[/grey50]")
            self.oniom_setup.use_explicit_connectivity = False

        # ========================================
        # STEP 4: Load atom types and charges
        # ========================================
        self.console.print("\n[bold cyan]Step 3/6: Loading Atom Types and Charges[/bold cyan]")
        self.console.print("  Atom types define the force field parameters for each atom.")
        self.console.print("  Charges are atomic partial charges from RESP or other methods.")

        if atom_types is not None:
            self.oniom_setup.atom_types = atom_types
            n_unique_types = len(set(atom_types.values()))
            self.console.print(f"\n  [green]✓[/green] Loaded {len(atom_types)} atom type assignments ({n_unique_types} unique types)")
        else:
            self.console.print("\n  [yellow]⚠[/yellow] No atom types provided, Gaussian will use defaults")

        if charges is not None:
            self.oniom_setup.charges = charges
            total_charge = sum(charges.values())
            self.console.print(f"  [green]✓[/green] Loaded {len(charges)} atomic charges (total: {total_charge:.2f})")
        else:
            self.console.print("  [yellow]⚠[/yellow] No charges provided, using 0.0 as default")

        # ========================================
        # STEP 5: Load force field parameters
        # ========================================
        self.console.print("\n[bold cyan]Step 4/6: Loading Force Field Parameters[/bold cyan]")
        self.console.print("  Parameters include bond lengths, angles, VDW radii, etc.")
        self.console.print("  These are written to the ONIOM input file (Gaussian only has AMBER94 built-in).")

        try:
            if prmtop_path:
                # Read all parameters directly from the prmtop (preferred path)
                from proprep.forcefield_prep.metal_site_parameterizer import PrmtopParameterProvider
                provider = PrmtopParameterProvider(prmtop_path, console=self.console)

                # Build an adapter that exposes the interface the ONIOM writer expects:
                # .get_nonbonded_parameter(), .get_bond_parameter(), .get_angle_parameter(),
                # .dihedral_params (dict with string keys and DihedralParam objects)
                provider.dihedral_params = provider.as_dihedral_params_dict()
                self.ff_reader = provider

                stats = provider.get_statistics()
                self.console.print(f"\n  [green]✓[/green] Parameters extracted from prmtop:")
                self.console.print(f"    {stats['n_bond']} bond types")
                self.console.print(f"    {stats['n_angle']} angle types")
                self.console.print(f"    {stats['n_dihedral']} dihedral types")
                self.console.print(f"    {stats['n_nonbonded']} VDW types")
            else:
                # Legacy path: load from AMBERHOME dat/frcmod files
                from proprep.forcefield_prep.mcpb.ff_parameter_reader import FFParameterReader
                import os

                self.ff_reader = object.__new__(FFParameterReader)
                self.ff_reader.force_field = "custom"
                self.ff_reader.logger = self.logger
                self.ff_reader.amberhome = amberhome or os.getenv('AMBERHOME', '')

                self.ff_reader.mass_params = {}
                self.ff_reader.nonbonded_params = {}
                self.ff_reader.bond_params = {}
                self.ff_reader.angle_params = {}
                self.ff_reader.dihedral_params = {}

                if standard_dat_files:
                    self.console.print(f"\n  [cyan]Loading base force field parameters:[/cyan]")
                    for dat_path in standard_dat_files:
                        dat_name = Path(dat_path).name
                        self.console.print(f"    • {dat_name}")
                        self.ff_reader._parse_parameter_file(dat_path, source=dat_name, overwrite=True)

                if standard_frcmod_files:
                    self.console.print(f"\n  [cyan]Loading standard force field parameters:[/cyan]")
                    for param_path in standard_frcmod_files:
                        param_name = Path(param_path).name
                        self.console.print(f"    • {param_name}")
                        self.ff_reader._parse_parameter_file(param_path, source=param_name, overwrite=True)

                if custom_frcmod_files:
                    self.console.print(f"\n  [cyan]Loading custom force field parameters:[/cyan]")
                    for frcmod_path in custom_frcmod_files:
                        frcmod_name = Path(frcmod_path).name
                        self.console.print(f"    • {frcmod_name}")
                        self.ff_reader._parse_frcmod_file(frcmod_path, source=frcmod_name, overwrite=True)

                stats = self.ff_reader.get_statistics()
                self.console.print(f"\n  [green]✓[/green] Loaded force field parameters:")
                self.console.print(f"    {stats['n_bonds']} bond types")
                self.console.print(f"    {stats['n_angles']} angle types")
                self.console.print(f"    {stats['n_dihedrals']} dihedral types")

        except Exception as e:
            self.console.print(f"\n  [yellow]⚠[/yellow] Could not load force field parameters: {e}")
            self.console.print("    ONIOM input will be generated without FF parameters")
            self.console.print("    Gaussian will use its built-in AMBER94 parameters")
            self.logger.error(f"Force field parameter loading failed: {e}", exc_info=True)
            import traceback
            traceback.print_exc()
            self.ff_reader = None

        # ========================================
        # STEP 6: Validate setup
        # ========================================
        self.console.print("\n[bold cyan]Step 5/6: Validating ONIOM Setup[/bold cyan]")
        self.console.print("  Checking for common errors and potential issues...")

        self.validator = ONIOMValidator(
            oniom_setup=self.oniom_setup,
            logger_instance=self.logger,
            console=self.console
        )

        validation_result = self.validator.validate_setup()
        self.oniom_setup.validation_passed = validation_result.passed

        if validation_result.passed:
            self.console.print(f"\n  [green]✓[/green] Validation passed")
            self.console.print(f"    {len(validation_result.warnings)} warnings, {len(validation_result.info)} info messages")
            if validation_result.warnings:
                self.console.print("    Review warnings before running calculation")
        else:
            self.console.print(f"\n  [red]✗[/red] Validation failed with {len(validation_result.errors)} errors")
            self.console.print("    Please review and fix errors before using this setup")

        # ========================================
        # STEP 7: Summary
        # ========================================
        self.console.print("\n[bold cyan]Step 6/6: Preparation Complete[/bold cyan]")
        self.console.print("\n  Setup is ready. Use write_input_file() to generate Gaussian input.")
        self.console.print("  Use print_validation_report() to see detailed validation results.")
        self.console.print("=" * 70)

        return self.oniom_setup

    def write_input_file(
        self,
        output_path: str,
        title: Optional[str] = None,
        include_comments: bool = True,
        include_ff_parameters: bool = True
    ) -> bool:
        """
        Write Gaussian ONIOM input file.

        Args:
            output_path: Path to output .com file
            title: Optional custom title
            include_comments: Include informative comments in file
            include_ff_parameters: Include force field parameter section

        Returns:
            True if successful, False otherwise

        Raises:
            RuntimeError: If setup not prepared yet
        """
        if self.oniom_setup is None:
            raise RuntimeError(
                "ONIOM setup not prepared. Call prepare_oniom_setup() first."
            )

        self.console.print(f"\n[bold cyan]Writing ONIOM Input File[/bold cyan]")
        self.console.print(f"  Output: {output_path}")

        # Create writer
        self.writer = ONIOMWriter(
            oniom_setup=self.oniom_setup,
            ff_parameter_reader=self.ff_reader if include_ff_parameters else None,
            logger_instance=self.logger,
            console=self.console
        )

        # Write file
        success = self.writer.write_input_file(
            output_path=output_path,
            title=title,
            include_comments=include_comments,
            include_ff_parameters=include_ff_parameters
        )

        if success:
            output_file = Path(output_path)
            file_size = output_file.stat().st_size / 1024  # KB
            self.console.print(f"\n  [green]✓[/green] Successfully wrote ONIOM input file ({file_size:.1f} KB)")
            self.console.print(f"    File: {output_path}")
            self.console.print("\n  This file is ready to submit to Gaussian for ONIOM calculation.")
        else:
            self.console.print(f"\n  [red]✗[/red] Failed to write ONIOM input file")

        return success

    def write_diagnostic_file(
        self,
        output_path: str,
        title: Optional[str] = None,
    ) -> bool:
        """
        Write a diagnostic Gaussian ONIOM input file.

        Uses ONIOM=InputFiles so Gaussian writes separate input files for
        each sub-calculation (real-low, model-high, model-low). Also uses
        the Test keyword so Gaussian does not run the SCF.

        This lets you verify layer assignments, atom counts, and charges
        for each ONIOM sub-calculation before running an expensive job.

        Returns:
            True if successful, False otherwise
        """
        if self.oniom_setup is None:
            raise RuntimeError(
                "ONIOM setup not prepared. Call prepare_oniom_setup() first."
            )

        self.console.print(f"\n[bold cyan]Writing ONIOM Diagnostic File[/bold cyan]")
        self.console.print(f"  Output: {output_path}")

        # Create writer if not already created
        if self.writer is None:
            self.writer = ONIOMWriter(
                oniom_setup=self.oniom_setup,
                ff_parameter_reader=self.ff_reader,
                logger_instance=self.logger,
                console=self.console
            )

        success = self.writer.write_diagnostic_file(
            output_path=output_path,
            title=title,
        )

        if success:
            output_file = Path(output_path)
            file_size = output_file.stat().st_size / 1024
            self.console.print(f"\n  [green]✓[/green] Wrote diagnostic file ({file_size:.1f} KB)")
            self.console.print(f"    File: {output_path}")
            self.console.print(
                "\n  Run this file with Gaussian. It will:"
                "\n    1. Write separate .com files for each ONIOM sub-calculation"
                "\n    2. Stop before the expensive SCF calculation (Test keyword)"
                "\n  Check the generated sub-files to verify layer assignments and charges."
            )
        else:
            self.console.print(f"\n  [red]✗[/red] Failed to write diagnostic file")

        return success

    def write_summary_report(self, output_path: str) -> bool:
        """
        Write human-readable summary report.

        This report explains the ONIOM setup in detail for documentation
        and record-keeping purposes.

        Args:
            output_path: Path to output .txt file

        Returns:
            True if successful, False otherwise
        """
        if self.oniom_setup is None or self.writer is None:
            raise RuntimeError(
                "ONIOM setup not prepared. Call prepare_oniom_setup() and "
                "write_input_file() first."
            )

        try:
            report = self.writer.generate_summary_report()

            output_file = Path(output_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)

            with open(output_file, 'w') as f:
                f.write(report)

            self.console.print(f"\n[green]✓[/green] Wrote summary report: {output_path}")
            return True

        except Exception as e:
            self.console.print(f"[red]✗[/red] Failed to write summary report: {e}")
            self.logger.error(f"Failed to write summary report: {e}")
            return False

    def print_validation_report(self):
        """
        Print detailed validation report to console.

        This shows all errors, warnings, and informational messages
        from the validation process.
        """
        if self.validator is None:
            self.console.print("[yellow]Validation not performed yet[/yellow]")
            return

        self.validator.print_validation_report()

    def get_validation_result(self) -> Optional[ValidationResult]:
        """
        Get validation result for programmatic access.

        Returns:
            ValidationResult or None if not validated yet
        """
        if self.validator is None:
            return None

        return self.validator.result
