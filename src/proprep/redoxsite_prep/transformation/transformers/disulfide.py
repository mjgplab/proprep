#!/usr/bin/env python3
"""
Disulfide Bond Transformer for RedoxSite System

This transformer handles disulfide-bonded cysteine residues by renaming CYS to CYX.
This is required for proper handling by TLeaP during MD simulation preparation.

Background:
- CYS: Cysteine with free thiol (-SH)
- CYX: Cysteine involved in disulfide bond (no hydrogen on sulfur)
- TLeaP requires CYX naming to properly recognize disulfide bonds

Use cases:
- All disulfide bonds detected by comprehensive_redox_detector
- Automatically assigned to disulfide-bonded CYS pairs

Author: Claude Code Implementation
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    RedoxSiteTransformerBase,
    TransformerEvaluation,
    TransformerEvaluationDetail,
    register_redox_transformer
)

# Import RedoxSite components
from proprep.structure_prep.comprehensive_redox_detector import CenterType

logger = logging.getLogger(__name__)

@register_redox_transformer
class DisulfideTransformer(RedoxSiteTransformerBase):
    """Transformer for disulfide-bonded cysteine residues (CYS → CYX renaming)"""

    # Transformer identification
    TRANSFORMER_NAME = "disulfide"
    DESCRIPTION = "Rename disulfide-bonded CYS residues to CYX for TLeaP compatibility"
    SUPPORTED_SITE_TYPES = ["disulfide"]
    FORCEFIELD_PATH = None  # No custom forcefield needed - standard AMBER naming

    @classmethod
    def evaluate_redox_site(cls, redox_site) -> TransformerEvaluation:
        """
        Evaluate if this site is a valid disulfide bond requiring CYS→CYX renaming.

        Requirements:
        - Exactly 2 centers (the two CYS residues)
        - Both centers are CYS residues
        - Both centers have is_disulfide_bonded property set to True
        """
        details = []
        requirements_met = 0
        total_requirements = 3

        # Requirement 1: Exactly 2 centers
        center_count = len(redox_site.centers) if hasattr(redox_site, 'centers') else 0
        center_count_check = (center_count == 2)
        details.append(
            TransformerEvaluationDetail(
                description="Must have exactly 2 centers (disulfide pair)",
                passed=center_count_check,
                value_found=center_count,
                value_expected=2,
                error_message=None if center_count_check else f"Found {center_count} centers, expected 2"
            )
        )
        if center_count_check:
            requirements_met += 1

        # Requirement 2: Both centers are CYS residues
        cys_count = 0
        if hasattr(redox_site, 'centers'):
            for center in redox_site.centers:
                if hasattr(center, 'resname') and center.resname in ['CYS', 'CYX']:
                    cys_count += 1

        cys_check = (cys_count == 2)
        details.append(
            TransformerEvaluationDetail(
                description="Both centers must be CYS residues",
                passed=cys_check,
                value_found=cys_count,
                value_expected=2,
                error_message=None if cys_check else f"Found {cys_count} CYS residues, expected 2"
            )
        )
        if cys_check:
            requirements_met += 1

        # Requirement 3: Both centers marked as disulfide-bonded
        disulfide_bonded_count = 0
        if hasattr(redox_site, 'centers'):
            for center in redox_site.centers:
                if hasattr(center, 'properties') and center.properties.get('is_disulfide_bonded', False):
                    disulfide_bonded_count += 1

        disulfide_check = (disulfide_bonded_count == 2)
        details.append(
            TransformerEvaluationDetail(
                description="Both centers must have is_disulfide_bonded property",
                passed=disulfide_check,
                value_found=disulfide_bonded_count,
                value_expected=2,
                error_message=None if disulfide_check else f"Found {disulfide_bonded_count} disulfide-bonded centers, expected 2"
            )
        )
        if disulfide_check:
            requirements_met += 1

        # Calculate confidence
        confidence = requirements_met / total_requirements if total_requirements > 0 else 0.0
        is_valid = (requirements_met == total_requirements)

        return TransformerEvaluation(
            is_valid=is_valid,
            confidence=confidence,
            description=cls.DESCRIPTION,
            requirements_met=requirements_met,
            total_requirements=total_requirements,
            details=details
        )

    @classmethod
    def get_site_requirements(cls) -> Dict[str, Any]:
        """
        Site requirements for disulfide bonds.
        """
        return {
            "centers": {
                "required_count": 2,
                "center_types": [CenterType.REDOX_AMINO_ACID],
                "residue_names": ["CYS", "CYX"]
            },
            "atoms": {
                "required_residues": {
                    "CYS": {"min_count": 0, "max_count": 2},
                    "CYX": {"min_count": 0, "max_count": 2}
                },
                "alternative_groups": [
                    ["CYS", "CYX"]  # Either CYS or CYX, or mix
                ]
            },
            "bonds": {
                "required_bond_groups": [
                    {
                        "description": "Disulfide bond",
                        "min_count": 1,
                        "bond_types": {
                            "disulfide": [
                                (("CYS", "SG"), ("CYS", "SG")),
                                (("CYX", "SG"), ("CYX", "SG")),
                                (("CYS", "SG"), ("CYX", "SG"))
                            ]
                        }
                    }
                ],
                "require_one_group": False  # Optional bond validation
            }
        }

    @classmethod
    def match_components(cls, redox_site) -> Tuple[Dict[str, Any], List[str]]:
        """
        Match the two CYS residues in the disulfide bond.

        Returns:
            matched_components: Dict with cys1 and cys2 chain/resid info
            missing_components: List of missing components (if any)
        """
        matched_components = {}
        missing_components = []

        # Get the two CYS centers
        if not hasattr(redox_site, 'centers') or len(redox_site.centers) < 2:
            missing_components.extend(["cys1_chain", "cys1_id", "cys2_chain", "cys2_id"])
            return matched_components, missing_components

        cys_centers = [c for c in redox_site.centers if hasattr(c, 'resname') and c.resname in ['CYS', 'CYX']]

        if len(cys_centers) >= 2:
            # Assign first and second CYS
            matched_components["cys1_chain"] = cys_centers[0].chain
            matched_components["cys1_id"] = cys_centers[0].resid
            matched_components["cys2_chain"] = cys_centers[1].chain
            matched_components["cys2_id"] = cys_centers[1].resid
        else:
            missing_components.extend(["cys1_chain", "cys1_id", "cys2_chain", "cys2_id"])

        return matched_components, missing_components

    @classmethod
    def get_transformation_sequence(cls, components: Dict[str, Any], parameters: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Generate transformation sequence to rename CYS to CYX.

        This is a simple transformation that renames both CYS residues to CYX
        so that TLeaP can properly recognize them as disulfide-bonded cysteines.

        Args:
            components: Matched site components (cys1 and cys2 info)
            parameters: User-specified parameters (none needed for this transformer)

        Returns:
            List of transformation operations
        """
        transformations = []

        # Only perform transformation if both CYS residues are present
        if not all(key in components for key in ["cys1_chain", "cys1_id", "cys2_chain", "cys2_id"]):
            logger.warning("Missing CYS component information - cannot generate transformations")
            return transformations

        # Transformation 1: Rename first CYS to CYX
        transformations.append({
            "id": "rename_cys1_to_cyx",
            "description": "Rename first CYS in disulfide bond to CYX",
            "selector": {
                "chain_id": components["cys1_chain"],
                "residue_name": "CYS",
                "residue_id": components["cys1_id"]
            },
            "action": {
                "change_residue_name": "CYX"
            }
        })

        # Transformation 2: Rename second CYS to CYX
        transformations.append({
            "id": "rename_cys2_to_cyx",
            "description": "Rename second CYS in disulfide bond to CYX",
            "selector": {
                "chain_id": components["cys2_chain"],
                "residue_name": "CYS",
                "residue_id": components["cys2_id"]
            },
            "action": {
                "change_residue_name": "CYX"
            }
        })

        return transformations

    @classmethod
    def get_required_residue_count(cls) -> int:
        """
        Disulfide bonds don't need new residue IDs - they use existing CYS positions.

        Returns:
            1 (minimum space needed to preserve the site)
        """
        return 1

    @classmethod
    def get_residue_space_plan(cls, components: Dict[str, Any]) -> Dict[str, int]:
        """
        No residue space planning needed - just renaming in place.

        Args:
            components: Matched site components (ignored)

        Returns:
            Empty dict - no residue allocation needed
        """
        return {}

    @classmethod
    def get_parameter_definitions(cls) -> Dict[str, Any]:
        """
        No parameters required - transformer is fully automatic.

        Returns:
            Empty dict - no user parameters needed
        """
        return {}

    @classmethod
    def validate_parameters(cls, parameters: Dict[str, Any]) -> Tuple[bool, str]:
        """
        Always valid - no parameters to validate.

        Args:
            parameters: Parameter values (ignored)

        Returns:
            Tuple of (True, "No parameters required")
        """
        return True, "No parameters required"

    @classmethod
    def get_parameter_mappings(cls, parameters: Dict[str, Any]) -> Dict[str, str]:
        """
        No parameter mappings needed.

        Args:
            parameters: User parameters (ignored)

        Returns:
            Empty dict - no mappings needed
        """
        return {}
