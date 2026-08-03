#!/usr/bin/env python3
"""
Component Matcher Collector

Interactive collector for defining component matching rules using multi-stage approach
"""

import logging
from typing import Dict, List, Optional, Any

from ..data_models import (
    ComponentMatching,
    ComponentMatchingRule,
    IntrinsicFilterStage,
    PositionalStage,
    ArithmeticStage,
    StructuralInfoStage,
    ExcludeMatchedStage,
    CustomStage,
    MatchingStage
)

logger = logging.getLogger(__name__)


class ComponentMatcherCollector:
    """Collect component matching rules interactively"""

    def __init__(self):
        self.matched_components = []  # Track what's been matched

    def collect(
        self,
        component_names: List[str],
        existing: Optional[ComponentMatching] = None
    ) -> ComponentMatching:
        """
        Collect component matching rules

        Args:
            component_names: List of component names that need matching rules
            existing: Existing matching rules (from template)

        Returns:
            Complete ComponentMatching object
        """
        print("\n" + "="*70)
        print("STEP 2: COMPONENT MATCHING RULES")
        print("="*70)
        print("\nDefine HOW to automatically identify components in detected sites.")
        print("Each component gets a multi-stage matching rule.\n")

        print("Stage types available:")
        print("  1. Filter by properties (residue name, coordinating atom, distance)")
        print("  2. Sort and select by position (e.g., closest, first by resid)")
        print("  3. Calculate from another component (e.g., c_ring_cys_id + 1)")
        print("  4. Exclude already matched")
        print("  5. Custom logic\n")

        rules = {}
        existing_rules = existing.rules if existing else {}

        for component in component_names:
            print(f"\n--- Define matching rule for: {component} ---")

            # Check if exists
            if component in existing_rules:
                print(f"\nExisting rule found for '{component}':")
                self._display_rule(existing_rules[component])
                print("\nOptions:")
                print("  1. Keep existing rule")
                print("  2. Edit existing rule")
                print("  3. Start fresh")
                choice = input("Choice [1]: ").strip() or "1"

                if choice == "1":
                    rules[component] = existing_rules[component]
                    self.matched_components.append(component)
                    continue
                elif choice == "2":
                    rule = self._edit_rule(component, existing_rules[component])
                else:
                    rule = self._create_new_rule(component)
            else:
                rule = self._create_new_rule(component)

            rules[component] = rule
            self.matched_components.append(component)

        return ComponentMatching(rules=rules)

    def _create_new_rule(self, component_name: str) -> ComponentMatchingRule:
        """Create a new matching rule for a component"""
        print(f"\nCreate matching rule for: {component_name}")

        # Get description
        description = input(f"Description [{component_name}]: ").strip()
        if not description:
            description = component_name

        # Collect stages
        stages = []
        print("\nAdd stages to narrow down candidates (blank to finish):")

        stage_num = 1
        while True:
            print(f"\n  Stage {stage_num}:")
            stage = self._collect_stage(component_name)
            if stage is None:
                break

            stages.append(stage)
            self._display_stage(stage, stage_num)
            stage_num += 1

        if not stages:
            print("Warning: No stages defined. This component won't be matched.")

        return ComponentMatchingRule(
            description=description,
            stages=stages
        )

    def _edit_rule(
        self,
        component_name: str,
        existing_rule: ComponentMatchingRule
    ) -> ComponentMatchingRule:
        """Edit an existing matching rule"""
        print(f"\nEdit matching rule for: {component_name}")

        # Edit description
        print(f"\nCurrent description: {existing_rule.description}")
        new_desc = input("New description (blank to keep): ").strip()
        description = new_desc if new_desc else existing_rule.description

        # Edit stages
        stages = list(existing_rule.stages)
        while True:
            print("\n  Current stages:")
            for i, stage in enumerate(stages, 1):
                self._display_stage(stage, i)

            print("\nOptions:")
            print("  1. Add stage")
            print("  2. Remove stage")
            print("  3. Done")
            choice = input("Choice [3]: ").strip() or "3"

            if choice == "1":
                new_stage = self._collect_stage(component_name)
                if new_stage:
                    stages.append(new_stage)
            elif choice == "2":
                if stages:
                    idx_str = input(f"Remove stage number (1-{len(stages)}): ").strip()
                    try:
                        idx = int(idx_str) - 1
                        if 0 <= idx < len(stages):
                            removed = stages.pop(idx)
                            print(f"Removed stage {idx+1}")
                    except ValueError:
                        print("Invalid number")
            elif choice == "3":
                break

        return ComponentMatchingRule(
            description=description,
            stages=stages
        )

    def _collect_stage(self, component_name: str) -> Optional[MatchingStage]:
        """Collect a single matching stage"""
        print("\n  Stage type:")
        print("    1. Filter by properties")
        print("    2. Sort and select by position")
        print("    3. Calculate from another component")
        print("    4. Exclude already matched")
        print("    5. Custom logic")
        print("    (blank to finish)")

        choice = input("  Type: ").strip()

        if not choice:
            return None
        elif choice == "1":
            return self._collect_intrinsic_filter()
        elif choice == "2":
            return self._collect_positional()
        elif choice == "3":
            return self._collect_arithmetic()
        elif choice == "4":
            return ExcludeMatchedStage()
        elif choice == "5":
            return self._collect_custom()
        else:
            print("Invalid choice")
            return None

    def _collect_intrinsic_filter(self) -> IntrinsicFilterStage:
        """Collect intrinsic filter stage"""
        print("\n  Filter by properties:")
        stage_data = {}

        # Residue name
        print("    Filter by residue name? (e.g., HIS, or HEM,HEC for multiple)")
        resname = input("    Residue name(s) [leave blank to skip]: ").strip().upper()
        if resname:
            if ',' in resname:
                stage_data['resnames'] = [r.strip() for r in resname.split(',')]
            else:
                stage_data['resname'] = resname

        # Coordinating atom
        print("    Filter by coordinating atom? (e.g., NE2, or OD1,OD2 for multi-dentate)")
        coord_atom = input("    Coordinating atom(s) [leave blank to skip]: ").strip().upper()
        if coord_atom:
            if ',' in coord_atom:
                stage_data['coord_atoms'] = [a.strip() for a in coord_atom.split(',')]
            else:
                stage_data['coord_atom'] = coord_atom

        # Distance range
        print("    Filter by distance range?")
        dist_choice = input("    Add distance filter? [y/N]: ").strip().lower()
        if dist_choice == 'y':
            min_dist = float(input("    Min distance (Å) [0]: ").strip() or "0")
            max_dist = float(input("    Max distance (Å) [3.0]: ").strip() or "3.0")
            stage_data['min_distance'] = min_dist
            stage_data['max_distance'] = max_dist

        return IntrinsicFilterStage(**stage_data)

    def _collect_positional(self) -> PositionalStage:
        """Collect positional stage"""
        print("\n  Sort and select by position:")

        # Sort by
        print("    Sort by:")
        print("      1. Residue ID")
        print("      2. Distance")
        print("      3. Chain")
        sort_choice = input("    Choice [1]: ").strip() or "1"

        sort_by_map = {'1': 'resid', '2': 'distance', '3': 'chain'}
        sort_by = sort_by_map.get(sort_choice, 'resid')

        # Reference for distance sorting
        reference_component = None
        reference_atom = None
        if sort_by == 'distance':
            print("    Measure distance from which component?")
            if self.matched_components:
                print(f"    Matched so far: {', '.join(self.matched_components)}")
            reference_component = input("    Component name: ").strip()
            reference_atom = input("    Atom name (optional): ").strip() or None

        # Position
        print("    Select position (0=first, -1=last):")
        position = int(input("    Position [0]: ").strip() or "0")

        # Reverse
        reverse = input("    Reverse sort order? [y/N]: ").strip().lower() == 'y'

        return PositionalStage(
            sort_by=sort_by,
            position=position,
            reverse=reverse,
            reference_component=reference_component,
            reference_atom=reference_atom
        )

    def _collect_arithmetic(self) -> ArithmeticStage:
        """Collect arithmetic stage"""
        print("\n  Calculate from another component:")

        # Base component
        print("    Base component:")
        if self.matched_components:
            print(f"    Available: {', '.join(self.matched_components)}")
        base_component = input("    Component name: ").strip()

        # Operation
        print("    Operation:")
        print("      1. Add")
        print("      2. Subtract")
        op_choice = input("    Choice [1]: ").strip() or "1"
        operation = "add" if op_choice == "1" else "subtract"

        # Value
        value = int(input("    Value: ").strip())

        return ArithmeticStage(
            base_component=base_component,
            operation=operation,
            value=value
        )

    def _collect_custom(self) -> CustomStage:
        """Collect custom stage"""
        print("\n  Custom logic:")
        description = input("    Description: ").strip()
        print("    Python code (single line):")
        code = input("    > ").strip()

        return CustomStage(
            code=code,
            description=description
        )

    def _display_rule(self, rule: ComponentMatchingRule):
        """Display a matching rule"""
        print(f"  Description: {rule.description}")
        print(f"  Stages ({len(rule.stages)}):")
        for i, stage in enumerate(rule.stages, 1):
            self._display_stage(stage, i)

    def _display_stage(self, stage: MatchingStage, stage_num: int):
        """Display a single stage"""
        stage_dict = stage.model_dump()
        stage_type = stage_dict.get('type', 'unknown')

        print(f"    {stage_num}. {stage_type}:", end="")

        if stage_type == "intrinsic_filter":
            parts = []
            if 'resname' in stage_dict and stage_dict['resname']:
                parts.append(f"resname={stage_dict['resname']}")
            if 'resnames' in stage_dict and stage_dict['resnames']:
                parts.append(f"resnames={stage_dict['resnames']}")
            if 'coord_atom' in stage_dict and stage_dict['coord_atom']:
                parts.append(f"coord_atom={stage_dict['coord_atom']}")
            if 'min_distance' in stage_dict and stage_dict['min_distance'] is not None:
                parts.append(f"distance={stage_dict['min_distance']}-{stage_dict.get('max_distance', 'inf')}")
            print(" " + ", ".join(parts) if parts else " (no filters)")

        elif stage_type == "positional":
            print(f" sort_by={stage_dict['sort_by']}, position={stage_dict['position']}, reverse={stage_dict.get('reverse', False)}")

        elif stage_type == "arithmetic":
            print(f" {stage_dict['base_component']} {stage_dict['operation']} {stage_dict['value']}")

        elif stage_type == "exclude_matched":
            print(" (exclude already matched components)")

        elif stage_type == "custom":
            print(f" {stage_dict.get('description', 'custom logic')}")

        else:
            print(f" {stage_dict}")
