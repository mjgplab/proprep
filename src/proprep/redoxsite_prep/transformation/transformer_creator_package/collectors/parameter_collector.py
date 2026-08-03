#!/usr/bin/env python3
"""
Parameter Collector

Interactive collector for defining parameters and parameter mappings
"""

import logging
from typing import Dict, List, Optional
from itertools import product

from ..data_models import (
    Parameters,
    ParameterDefinition,
    ParameterType
)

logger = logging.getLogger(__name__)


class ParameterCollector:
    """Collect parameter definitions and mappings interactively"""

    def __init__(self):
        pass

    def collect(self, existing: Optional[Parameters] = None) -> Optional[Parameters]:
        """
        Collect parameters and mappings

        Args:
            existing: Existing parameters (from template)

        Returns:
            Complete Parameters object, or None if no parameters needed
        """
        print("\n" + "="*70)
        print("STEP 4: PARAMETERS AND MAPPINGS (Optional)")
        print("="*70)
        print("\nParameters allow users to control transformer behavior.")
        print("Common use: redox_state and spin_state for redox-active cofactors.\n")

        print("Does this transformer need user-configurable parameters?")
        choice = input("Add parameters? [y/N]: ").strip().lower()

        if choice != 'y':
            return None

        # Collect parameter definitions
        definitions = self._collect_parameter_definitions(
            existing.definitions if existing else {}
        )

        if not definitions:
            return None

        # Generate parameter combinations
        param_combinations = self._generate_combinations(definitions)
        print(f"\nParameter combinations: {param_combinations}")

        # Collect mappings
        mappings = self._collect_parameter_mappings(
            param_combinations,
            existing.mappings if existing else {}
        )

        return Parameters(
            definitions=definitions,
            mappings=mappings
        )

    def _collect_parameter_definitions(
        self,
        existing: Dict[str, ParameterDefinition]
    ) -> Dict[str, ParameterDefinition]:
        """Collect parameter definitions"""
        print("\n--- Parameter Definitions ---")

        definitions = dict(existing)

        if existing:
            print("\nExisting parameters:")
            for name, defn in existing.items():
                print(f"  {name}: {defn.type} - {defn.description}")
            print()

        print("Add parameters (blank name to finish):")

        while True:
            name = input("\n  Parameter name (e.g., redox_state, spin_state): ").strip()
            if not name:
                break

            # Check if editing existing
            if name in definitions:
                print(f"  Editing existing parameter: {name}")

            # Description
            description = input("  Description: ").strip()
            if not description:
                description = name.replace('_', ' ')

            # Type
            print("  Type:")
            print("    1. Choice (user selects from options)")
            print("    2. Fixed (always the same value)")
            type_choice = input("  Choice [1]: ").strip() or "1"
            param_type = ParameterType.CHOICE if type_choice == "1" else ParameterType.FIXED

            # Type-specific fields
            if param_type == ParameterType.CHOICE:
                print("  Options (comma-separated):")
                print("    Example: reduced,oxidized")
                options_str = input("  > ").strip()
                options = [opt.strip() for opt in options_str.split(',')]

                default_str = input(f"  Default [{options[0]}]: ").strip()
                default = default_str if default_str else options[0]

                definitions[name] = ParameterDefinition(
                    description=description,
                    type=param_type,
                    options=options,
                    default=default
                )

            else:  # FIXED
                value = input("  Value: ").strip()
                note = input("  Note (optional): ").strip() or None

                definitions[name] = ParameterDefinition(
                    description=description,
                    type=param_type,
                    value=value,
                    note=note
                )

            print(f"  Added parameter: {name}")

        return definitions

    def _generate_combinations(
        self,
        definitions: Dict[str, ParameterDefinition]
    ) -> List[str]:
        """Generate all parameter combinations"""
        # Get choice parameters only
        choice_params = []
        for name, defn in definitions.items():
            if defn.type == ParameterType.CHOICE and defn.options:
                choice_params.append((name, defn.options))

        if not choice_params:
            return []

        # Generate combinations
        combinations = []
        param_names = [name for name, _ in choice_params]
        param_options = [options for _, options in choice_params]

        for combo in product(*param_options):
            combo_str = "_".join(combo)
            combinations.append(combo_str)

        return combinations

    def _collect_parameter_mappings(
        self,
        param_combinations: List[str],
        existing: Dict[str, Dict[str, str]]
    ) -> Dict[str, Dict[str, str]]:
        """Collect parameter mappings"""
        print("\n--- Parameter Mappings ---")
        print("Define how parameters map to residue names.\n")

        print("Example:")
        print("  Combination 'reduced_low_spin' might map to:")
        print("    heme_name: HCR")
        print("    proximal_ligand_name: PHR")
        print("    distal_ligand_name: DHR\n")

        mappings = dict(existing)

        # Get mapping keys (what to map)
        print("What should be mapped? (comma-separated)")
        print("  Common: heme_name, proximal_ligand_name, distal_ligand_name")
        keys_str = input("  Keys: ").strip()
        if not keys_str:
            return mappings

        mapping_keys = [k.strip() for k in keys_str.split(',')]

        # Collect mappings for each combination
        for combo in param_combinations:
            print(f"\n  Mappings for '{combo}':")

            if combo in mappings:
                print(f"    Existing: {mappings[combo]}")
                print("    Options:")
                print("      1. Keep existing")
                print("      2. Edit")
                choice = input("    Choice [1]: ").strip() or "1"
                if choice == "1":
                    continue

            combo_mappings = {}
            for key in mapping_keys:
                default = ""
                if combo in existing and key in existing[combo]:
                    default = existing[combo][key]

                value = input(f"    {key} [{default}]: ").strip()
                if value:
                    combo_mappings[key] = value.upper()
                elif default:
                    combo_mappings[key] = default

            if combo_mappings:
                mappings[combo] = combo_mappings

        return mappings
