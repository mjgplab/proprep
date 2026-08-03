"""
ONIOM Validator - Quality Control for ONIOM Setup

Validates ONIOMSetup configurations for common errors and potential issues.
Uses parmed integer atom indices as keys throughout.
"""

import logging
from typing import List, Optional
from collections import defaultdict
from rich.console import Console

from .data_structures import (
    ONIOMSetup,
    ONIOMLayer,
    FreezeFlag,
    calculate_layer_statistics,
)
from .structure_utils import calculate_distance


logger = logging.getLogger(__name__)


class ValidationResult:
    """Results from ONIOM setup validation."""

    def __init__(self):
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.info: List[str] = []

    @property
    def passed(self) -> bool:
        return len(self.errors) == 0

    @property
    def has_warnings(self) -> bool:
        return len(self.warnings) > 0

    def add_error(self, message: str):
        self.errors.append(message)
        logger.error(f"Validation error: {message}")

    def add_warning(self, message: str):
        self.warnings.append(message)
        logger.warning(f"Validation warning: {message}")

    def add_info(self, message: str):
        self.info.append(message)
        logger.info(f"Validation info: {message}")


class ONIOMValidator:
    """Validates ONIOMSetup configurations."""

    def __init__(self, oniom_setup: ONIOMSetup, logger_instance=None, console=None):
        self.oniom_setup = oniom_setup
        self.logger = logger_instance or logger
        self.console = console or Console()

    def print_validation_report(self):
        """Print a formatted validation report to the console."""
        result = self.validate_setup()

        if result.errors:
            self.console.print("\n[bold red]Errors:[/bold red]")
            for msg in result.errors:
                self.console.print(f"  [red]✗ {msg}[/red]")

        if result.warnings:
            self.console.print("\n[bold yellow]Warnings:[/bold yellow]")
            for msg in result.warnings:
                self.console.print(f"  [yellow]⚠ {msg}[/yellow]")

        if result.info:
            self.console.print("\n[bold cyan]Info:[/bold cyan]")
            for msg in result.info:
                self.console.print(f"  {msg}")

        if result.passed:
            self.console.print("\n[bold green]✓ Validation passed[/bold green]")
        else:
            self.console.print("\n[bold red]✗ Validation failed[/bold red]")

    def validate_setup(self) -> ValidationResult:
        """Run all validation checks."""
        result = ValidationResult()

        self._validate_layer_assignments(result)
        self._validate_geometry(result)
        self._validate_link_atoms(result)
        self._validate_atom_types_and_charges(result)
        self._validate_charge_balance(result)
        self._validate_qm_settings(result)

        return result

    def _validate_layer_assignments(self, result: ValidationResult):
        """Validate layer assignments."""
        assignments = self.oniom_setup.layer_assignments

        if not assignments:
            result.add_error("No layer assignments found")
            return

        stats = calculate_layer_statistics(self.oniom_setup)

        result.add_info(f"Total atoms: {stats.n_high + stats.n_medium + stats.n_low}")
        result.add_info(f"HIGH: {stats.n_high}, MEDIUM: {stats.n_medium}, LOW: {stats.n_low}")
        result.add_info(f"Link atoms: {stats.n_link}")
        result.add_info(f"Frozen: {stats.n_frozen}, Active: {stats.n_active}")

        if stats.n_high == 0:
            result.add_error("No atoms in HIGH layer (QM region is empty)")

        if stats.n_low == 0:
            result.add_warning("No atoms in LOW layer (entire system is QM)")

        if self.oniom_setup.n_layers == 3 and stats.n_medium == 0:
            result.add_warning("3-layer ONIOM requested but MEDIUM layer is empty")

        if stats.n_high > 500:
            result.add_warning(
                f"Large QM region ({stats.n_high} atoms) may be computationally expensive"
            )

    def _validate_geometry(self, result: ValidationResult):
        """Check for short interatomic distances that would crash Gaussian."""
        parm = self.oniom_setup.parm
        if parm is None:
            return

        CLASH_THRESHOLD = 0.5  # Angstroms — Gaussian rejects atoms this close

        # Check all bonded pairs for unreasonably short distances
        clashes = []
        for bond in parm.bonds:
            idx1 = bond.atom1.idx
            idx2 = bond.atom2.idx
            if idx1 not in self.oniom_setup.layer_assignments:
                continue
            if idx2 not in self.oniom_setup.layer_assignments:
                continue

            a1 = bond.atom1
            a2 = bond.atom2
            dist = calculate_distance(
                (a1.xx, a1.xy, a1.xz), (a2.xx, a2.xy, a2.xz))

            if dist < CLASH_THRESHOLD:
                clashes.append((idx1, idx2, dist))

        if clashes:
            for idx1, idx2, dist in clashes:
                a1 = parm.atoms[idx1]
                a2 = parm.atoms[idx2]
                result.add_warning(
                    f"Short distance ({dist:.3f} A): "
                    f"{a1.residue.name}{a1.residue.idx + 1}:{a1.name} — "
                    f"{a2.residue.name}{a2.residue.idx + 1}:{a2.name}. "
                    f"Gaussian may reject this. Consider minimizing the structure first."
                )

    def _validate_link_atoms(self, result: ValidationResult):
        """Validate link atom placements."""
        link_atoms = self.oniom_setup.link_atoms
        assignments = self.oniom_setup.layer_assignments

        if not link_atoms:
            if len(assignments) > 0:
                n_high = sum(1 for a in assignments.values() if a.layer == ONIOMLayer.HIGH)
                n_low = sum(1 for a in assignments.values() if a.layer == ONIOMLayer.LOW)
                if n_high > 0 and n_low > 0:
                    result.add_warning("No link atoms but both HIGH and LOW layers present")
            return

        for i, la in enumerate(link_atoms):
            # Verify QM parent is in HIGH/MEDIUM layer
            qm_assignment = assignments.get(la.qm_parent_idx)
            if qm_assignment is None:
                result.add_error(f"Link atom {i+1}: QM parent idx {la.qm_parent_idx} not in assignments")
            elif qm_assignment.layer not in (ONIOMLayer.HIGH, ONIOMLayer.MEDIUM):
                result.add_error(
                    f"Link atom {i+1}: QM parent is in {qm_assignment.layer.value} "
                    f"layer (expected HIGH or MEDIUM)"
                )

            # Verify MM parent is in LOW layer
            mm_assignment = assignments.get(la.mm_parent_idx)
            if mm_assignment is None:
                result.add_error(f"Link atom {i+1}: MM parent idx {la.mm_parent_idx} not in assignments")
            elif mm_assignment.layer != ONIOMLayer.LOW:
                result.add_error(
                    f"Link atom {i+1}: MM parent is in {mm_assignment.layer.value} "
                    f"layer (expected LOW)"
                )

            # Check link atom distance from QM parent
            dist_to_qm = calculate_distance(la.coords, la.qm_parent_coords)
            if dist_to_qm < 0.3:
                result.add_warning(f"Link atom {i+1}: very close to QM parent ({dist_to_qm:.2f} A)")
            elif dist_to_qm > 2.0:
                result.add_warning(f"Link atom {i+1}: far from QM parent ({dist_to_qm:.2f} A)")

        # Check for duplicate link atoms
        for i in range(len(link_atoms)):
            for j in range(i + 1, len(link_atoms)):
                dist = calculate_distance(link_atoms[i].coords, link_atoms[j].coords)
                if dist < 0.01:
                    result.add_warning(f"Duplicate link atoms at indices {i+1} and {j+1}")

    def _validate_atom_types_and_charges(self, result: ValidationResult):
        """Validate atom types and charges are present."""
        atom_types = self.oniom_setup.atom_types
        charges = self.oniom_setup.charges
        assignments = self.oniom_setup.layer_assignments

        if not atom_types:
            result.add_warning("No atom types assigned (Gaussian will use defaults)")
            return

        # Check all assigned atoms have types and charges
        missing_types = 0
        missing_charges = 0
        for atom_idx in assignments:
            if atom_idx not in atom_types:
                missing_types += 1
            if atom_idx not in charges:
                missing_charges += 1

        if missing_types > 0:
            result.add_warning(f"{missing_types} atoms missing atom type assignments")
        if missing_charges > 0:
            result.add_warning(f"{missing_charges} atoms missing charge assignments")

        # Check for unusual charges
        for atom_idx, charge in charges.items():
            if abs(charge) > 3.0:
                assignment = assignments.get(atom_idx)
                if assignment:
                    result.add_warning(
                        f"Unusually large charge ({charge:.4f}) on "
                        f"{assignment.atom_name} in {assignment.residue_name}{assignment.residue_number}"
                    )

    def _validate_charge_balance(self, result: ValidationResult):
        """Validate total charge is close to an integer."""
        charges = self.oniom_setup.charges
        if not charges:
            return

        total = sum(charges.values())
        nearest_int = round(total)
        deviation = abs(total - nearest_int)

        result.add_info(f"Total system charge: {total:.6f} (nearest integer: {nearest_int})")

        if deviation > 0.01:
            result.add_warning(
                f"Total charge ({total:.4f}) deviates from nearest integer "
                f"({nearest_int}) by {deviation:.4f}"
            )

    def _validate_qm_settings(self, result: ValidationResult):
        """Validate QM calculation settings."""
        setup = self.oniom_setup

        if not setup.qm_functional:
            result.add_error("No QM functional specified")
        if not setup.qm_basis_set:
            result.add_error("No QM basis set specified")

        if setup.qm_multiplicity < 1:
            result.add_error(f"Invalid QM multiplicity: {setup.qm_multiplicity}")

        if setup.n_layers == 3:
            if not setup.medium_method:
                result.add_warning("3-layer ONIOM but no MEDIUM layer method specified")
