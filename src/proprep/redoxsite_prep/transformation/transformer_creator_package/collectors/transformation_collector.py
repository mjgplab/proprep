#!/usr/bin/env python3
"""
Transformation Collector

Interactive collector for defining the transformation sequence
"""

import logging
from typing import Dict, List, Optional

from ..data_models import (
    TransformationSequence,
    Transformation,
    TransformationSelector,
    TransformationAction
)

logger = logging.getLogger(__name__)


class TransformationCollector:
    """Collect transformation sequence interactively"""

    def __init__(self):
        self.component_names = []

    def collect(
        self,
        component_names: List[str],
        existing: Optional[TransformationSequence] = None
    ) -> TransformationSequence:
        """
        Collect transformation sequence

        Args:
            component_names: Available component names for selectors
            existing: Existing transformations (from template)

        Returns:
            Complete TransformationSequence object
        """
        self.component_names = component_names

        print("\n" + "="*70)
        print("STEP 3: TRANSFORMATION SEQUENCE")
        print("="*70)
        print("\nDefine the transformations to apply to matched components.")
        print("Each transformation has a SELECTOR (what atoms) and ACTION (what to do).\n")

        transformations = []
        if existing and existing.transformations:
            print("Existing transformations found:")
            for i, trans in enumerate(existing.transformations, 1):
                print(f"  {i}. {trans.id}: {trans.description}")

            print("\nOptions:")
            print("  1. Keep all existing transformations")
            print("  2. Edit existing transformations")
            print("  3. Start fresh")
            choice = input("Choice [2]: ").strip() or "2"

            if choice == "1":
                return existing
            elif choice == "2":
                transformations = self._edit_transformations(existing.transformations)
            # else: start fresh with empty list

        # Add new transformations
        print("\nAdd transformations (blank ID to finish):")
        trans_num = len(transformations) + 1

        while True:
            trans = self._collect_transformation(trans_num)
            if trans is None:
                break

            transformations.append(trans)
            print(f"\n  Added transformation {trans_num}: {trans.id}")
            self._display_transformation(trans)
            trans_num += 1

        return TransformationSequence(transformations=transformations)

    def _edit_transformations(
        self,
        existing_trans: List[Transformation]
    ) -> List[Transformation]:
        """Edit existing transformations"""
        transformations = list(existing_trans)

        while True:
            print("\nCurrent transformations:")
            for i, trans in enumerate(transformations, 1):
                print(f"  {i}. {trans.id}: {trans.description}")

            print("\nOptions:")
            print("  1. Add transformation")
            print("  2. Remove transformation")
            print("  3. Done")
            choice = input("Choice [3]: ").strip() or "3"

            if choice == "1":
                trans = self._collect_transformation(len(transformations) + 1)
                if trans:
                    transformations.append(trans)
            elif choice == "2":
                if transformations:
                    idx_str = input(f"Remove number (1-{len(transformations)}): ").strip()
                    try:
                        idx = int(idx_str) - 1
                        if 0 <= idx < len(transformations):
                            removed = transformations.pop(idx)
                            print(f"Removed: {removed.id}")
                    except ValueError:
                        print("Invalid number")
            elif choice == "3":
                break

        return transformations

    def _collect_transformation(self, trans_num: int) -> Optional[Transformation]:
        """Collect a single transformation"""
        print(f"\n--- Transformation {trans_num} ---")

        # ID
        trans_id = input("  ID (e.g., rename_hem, move_cys_sidechain): ").strip()
        if not trans_id:
            return None

        # Description
        description = input("  Description: ").strip()
        if not description:
            description = trans_id.replace('_', ' ')

        # Selector
        print("\n  Selector (which atoms to transform):")
        selector = self._collect_selector()

        # Action
        print("\n  Action (what to do):")
        action = self._collect_action()

        return Transformation(
            id=trans_id,
            description=description,
            selector=selector,
            action=action
        )

    def _collect_selector(self) -> TransformationSelector:
        """Collect transformation selector"""
        selector_data = {}

        print("    Available components:", ", ".join(self.component_names))

        # Chain ID
        chain = input("    Chain ID (or component field like 'center_chain'): ").strip()
        if chain:
            selector_data['chain_id'] = chain

        # Residue ID
        resid = input("    Residue ID (or component field like 'center_id'): ").strip()
        if resid:
            # Try to parse as int, otherwise keep as string (component reference)
            try:
                selector_data['residue_id'] = int(resid)
            except ValueError:
                selector_data['residue_id'] = resid

        # Residue name
        resname = input("    Residue name (e.g., HIS, CYS): ").strip().upper()
        if resname:
            selector_data['residue_name'] = resname

        # Atom names
        atoms = input("    Atom names (comma-separated, or blank for all): ").strip().upper()
        if atoms:
            selector_data['atom_names'] = [a.strip() for a in atoms.split(',')]

        return TransformationSelector(**selector_data)

    def _collect_action(self) -> TransformationAction:
        """Collect transformation action"""
        action_data = {}

        print("    Available actions:")
        print("      1. Change residue name")
        print("      2. Change residue ID")
        print("      3. Change chain ID")
        print("      4. Rename atoms")
        print("      5. Convert to HETATM")
        print("      6. Multiple actions")

        choice = input("    Choice [6]: ").strip() or "6"

        if choice == "1" or choice == "6":
            new_resname = input("    New residue name: ").strip().upper()
            if new_resname:
                action_data['change_residue_name'] = new_resname

        if choice == "2" or choice == "6":
            new_resid = input("    New residue ID (or component field): ").strip()
            if new_resid:
                try:
                    action_data['change_residue_id'] = int(new_resid)
                except ValueError:
                    action_data['change_residue_id'] = new_resid

        if choice == "3" or choice == "6":
            new_chain = input("    New chain ID (or component field): ").strip()
            if new_chain:
                action_data['change_chain_id'] = new_chain

        if choice == "4" or choice == "6":
            print("    Rename atoms (format: OLD=NEW, comma-separated):")
            print("    Example: CB=CBB2,SG=SGB2")
            atom_renames = input("    Atom renames: ").strip()
            if atom_renames:
                rename_dict = {}
                for pair in atom_renames.split(','):
                    if '=' in pair:
                        old, new = pair.split('=', 1)
                        rename_dict[old.strip().upper()] = new.strip().upper()
                if rename_dict:
                    action_data['rename_atoms'] = rename_dict

        if choice == "5" or choice == "6":
            hetatm = input("    Convert to HETATM? [y/N]: ").strip().lower()
            if hetatm == 'y':
                action_data['convert_to_hetatm'] = True

        return TransformationAction(**action_data)

    def _display_transformation(self, trans: Transformation):
        """Display a transformation"""
        print(f"    Selector: {trans.selector.model_dump(exclude_none=True)}")
        print(f"    Action: {trans.action.model_dump(exclude_none=True)}")
