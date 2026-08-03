#!/usr/bin/env python3
"""
Site Requirements Collector

Interactive collector for defining site matching requirements
(centers, atoms/residues, bonds)
"""

import logging
from typing import Dict, List, Optional

from ..data_models import (
    SiteRequirements,
    CenterRequirements,
    AtomRequirements,
    ResidueRequirement,
    BondRequirements,
    BondGroup,
    CenterTypeEnum
)

logger = logging.getLogger(__name__)


class SiteRequirementsCollector:
    """Collect site requirements interactively"""

    def __init__(self):
        pass

    def collect(self, existing: Optional[SiteRequirements] = None) -> SiteRequirements:
        """
        Collect site requirements from user

        Args:
            existing: Existing requirements to edit (from template)

        Returns:
            Complete SiteRequirements object
        """
        print("\n" + "="*70)
        print("STEP 1: SITE REQUIREMENTS")
        print("="*70)
        print("\nDefine what makes a redox site match this transformer.")
        print("This includes: centers, required residues, and bond patterns.\n")

        # Collect center requirements
        centers = self._collect_center_requirements(
            existing.centers if existing else None
        )

        # Ask if user wants to specify atom requirements
        print("\nDo you want to specify required residues/atoms?")
        print("  1. Yes - define specific residue requirements")
        print("  2. No - skip (less restrictive)")
        choice = input("Choice [1/2]: ").strip()

        atoms = None
        if choice == '1':
            atoms = self._collect_atom_requirements(
                existing.atoms if existing else None
            )

        # Ask if user wants to specify bond requirements
        print("\nDo you want to specify required bonds?")
        print("  1. Yes - define bond patterns (advanced)")
        print("  2. No - skip (recommended for first iteration)")
        choice = input("Choice [1/2]: ").strip()

        bonds = None
        if choice == '1':
            bonds = self._collect_bond_requirements(
                existing.bonds if existing else None
            )

        return SiteRequirements(
            centers=centers,
            atoms=atoms,
            bonds=bonds
        )

    def _collect_center_requirements(
        self,
        existing: Optional[CenterRequirements] = None
    ) -> CenterRequirements:
        """Collect center requirements"""
        print("\n--- Center Requirements ---")

        # Number of centers
        default_count = existing.required_count if existing else 1
        print(f"\nHow many centers are required? [{default_count}]")
        count_str = input("> ").strip()
        count = int(count_str) if count_str else default_count

        # Center types
        print("\nWhat types of centers are allowed?")
        print("  1. Metal ions (Fe, Cu, Zn, etc. - isolated single atoms)")
        print("  2. Organometallic cofactors (heme, Fe-S clusters - metal in multi-atom residue)")
        print("  3. Organic cofactors (non-standard residues without metal)")
        print("  4. Redox-active amino acids (TYR, TRP, CYS, etc.)")
        print("  5. Multiple types (comma-separated: 1,2)")

        default_choice = "1"
        if existing and existing.center_types:
            type_nums = []
            for ct in existing.center_types:
                if ct == CenterTypeEnum.METAL_ION:
                    type_nums.append("1")
                elif ct == CenterTypeEnum.ORGANOMETALLIC_COFACTOR:
                    type_nums.append("2")
                elif ct == CenterTypeEnum.ORGANIC_COFACTOR:
                    type_nums.append("3")
                elif ct == CenterTypeEnum.REDOX_AMINO_ACID:
                    type_nums.append("4")
            default_choice = ",".join(type_nums)

        choice = input(f"Choice [{default_choice}]: ").strip() or default_choice

        center_types = []
        for c in choice.split(','):
            c = c.strip()
            if c == '1':
                center_types.append(CenterTypeEnum.METAL_ION)
            elif c == '2':
                center_types.append(CenterTypeEnum.ORGANOMETALLIC_COFACTOR)
            elif c == '3':
                center_types.append(CenterTypeEnum.ORGANIC_COFACTOR)
            elif c == '4':
                center_types.append(CenterTypeEnum.REDOX_AMINO_ACID)

        # Element specification (for metal ions)
        elements = None
        if CenterTypeEnum.METAL_ION in center_types:
            print("\nSpecify required element(s) for metal centers?")
            print("  Examples: FE, CU, ZN")
            print("  Leave blank for any metal")

            default_elem = ""
            if existing and existing.elements:
                default_elem = ",".join(existing.elements)

            elem_input = input(f"Element(s) [{default_elem or 'any'}]: ").strip()
            if elem_input:
                elements = [e.strip().upper() for e in elem_input.split(',')]
            elif default_elem:
                elements = [e.strip().upper() for e in default_elem.split(',')]

        return CenterRequirements(
            required_count=count,
            center_types=center_types,
            elements=elements
        )

    def _collect_atom_requirements(
        self,
        existing: Optional[AtomRequirements] = None
    ) -> AtomRequirements:
        """Collect atom/residue requirements"""
        print("\n--- Residue Requirements ---")
        print("Specify which residues must be present in the site.")
        print("Format: RESNAME min_count max_count")
        print("Example: HIS 2 2 (exactly 2 histidines)")
        print("Example: CYS 1 None (at least 1 cysteine)")
        print("Enter blank line when done.\n")

        required_residues = {}
        if existing and existing.required_residues:
            required_residues = dict(existing.required_residues)
            print("Existing requirements:")
            for resname, req in required_residues.items():
                max_str = str(req.max_count) if req.max_count else "None"
                print(f"  {resname}: {req.min_count} - {max_str}")
            print()

        print("Add/modify requirements (blank line to finish):")
        while True:
            line = input("> ").strip()
            if not line:
                break

            parts = line.split()
            if len(parts) < 2:
                print("Invalid format. Use: RESNAME min_count [max_count]")
                continue

            resname = parts[0].upper()
            min_count = int(parts[1])
            max_count = int(parts[2]) if len(parts) > 2 and parts[2].lower() != 'none' else None

            required_residues[resname] = ResidueRequirement(
                min_count=min_count,
                max_count=max_count
            )
            max_str = str(max_count) if max_count else "any"
            print(f"  Added: {resname} ({min_count} - {max_str})")

        # Alternative groups
        alternative_groups = None
        if existing and existing.alternative_groups:
            alternative_groups = existing.alternative_groups

        print("\nDefine alternative residue groups?")
        print("(Groups where only ONE from the group is required)")
        print("Example: HEM,HEC (either HEM or HEC, not both)")
        choice = input("Add alternative groups? [y/N]: ").strip().lower()

        if choice == 'y':
            alternative_groups = []
            print("Enter groups (comma-separated, blank line to finish):")
            while True:
                line = input("> ").strip()
                if not line:
                    break
                group = [r.strip().upper() for r in line.split(',')]
                alternative_groups.append(group)
                print(f"  Added group: {group}")

        return AtomRequirements(
            required_residues=required_residues,
            alternative_groups=alternative_groups
        )

    def _collect_bond_requirements(
        self,
        existing: Optional[BondRequirements] = None
    ) -> BondRequirements:
        """Collect bond requirements (simplified for now)"""
        print("\n--- Bond Requirements ---")
        print("Bond requirements are advanced and optional.")
        print("For most transformers, residue requirements are sufficient.\n")

        print("Define bond requirements?")
        print("  1. Yes - define bond patterns (complex)")
        print("  2. No - skip for now")
        choice = input("Choice [2]: ").strip() or "2"

        if choice == '2':
            return BondRequirements(
                required_bond_groups=[],
                require_one_group=False
            )

        # Simplified bond collection
        print("\nNote: Full bond requirement specification is complex.")
        print("Consider defining this directly in the generated code.\n")

        # Return minimal structure
        return BondRequirements(
            required_bond_groups=[],
            require_one_group=False
        )
