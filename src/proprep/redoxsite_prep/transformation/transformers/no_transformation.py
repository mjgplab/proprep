#!/usr/bin/env python3
"""
No-Transformation Transformer for RedoxSite System

This transformer is a pass-through for redox sites that are already forcefield-compliant
and require no modifications to residue names, atom names, or other PDB metadata.

Use cases:
- Standard metal binding sites (e.g., Ca2+ sites with standard naming)
- Cofactor sites already properly named
- Any detected redox site that needs tracking but no transformation

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

logger = logging.getLogger(__name__)

@register_redox_transformer
class NoTransformationTransformer(RedoxSiteTransformerBase):
    """Pass-through transformer for sites that are already forcefield-compliant"""

    # Transformer identification
    TRANSFORMER_NAME = "no_transformation"
    DESCRIPTION = "Pass-through for sites that are already forcefield-compliant (no modifications needed)"
    SUPPORTED_SITE_TYPES = ["*"]  # Compatible with all site types
    FORCEFIELD_PATH = None  # No forcefield files needed

    @classmethod
    def evaluate_redox_site(cls, redox_site) -> TransformerEvaluation:
        """
        Always returns valid - this transformer accepts any site.

        The philosophy: if a user explicitly assigns this transformer to a site,
        they're declaring that the site is already correct and needs no transformation.
        """
        # Count atoms in the site for the requirements count
        atom_count = len(redox_site.atoms) if hasattr(redox_site, 'atoms') else 0

        details = [
            TransformerEvaluationDetail(
                description="Site validation - no transformation needed",
                passed=True,
                value_found=atom_count,
                value_expected="any",
                error_message=None
            )
        ]

        return TransformerEvaluation(
            is_valid=True,  # Always valid
            confidence=1.0,  # Always 100% confident
            description=cls.DESCRIPTION,
            requirements_met=1,
            total_requirements=1,
            details=details
        )

    @classmethod
    def get_site_requirements(cls) -> Dict[str, Any]:
        """
        No requirements - accepts any site structure.
        """
        return {
            "centers": [],
            "residues": [],
            "bonds": []
        }

    @classmethod
    def match_components(cls, redox_site) -> Tuple[Dict[str, Any], List[str]]:
        """
        Returns empty components since no transformation is needed.

        Since this transformer doesn't modify anything, there's no need
        to identify specific components. The site is accepted as-is.

        Returns:
            Empty dict and empty missing list
        """
        matched_components = {}
        missing_components = []

        return matched_components, missing_components

    @classmethod
    def get_transformation_sequence(cls, components: Dict[str, Any], parameters: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Returns empty transformation sequence - no modifications needed.

        Args:
            components: Matched site components (ignored)
            parameters: User-specified parameters (ignored)

        Returns:
            Empty list - no transformations to apply
        """
        return []

    @classmethod
    def get_required_residue_count(cls) -> int:
        """
        Requires 1 residue ID to preserve the existing site residue.

        The site already exists as 1 residue and needs to occupy
        1 residue ID in the transformed structure as well.

        Returns:
            1 (one residue space for the existing site)
        """
        return 1

    @classmethod
    def get_residue_space_plan(cls, components: Dict[str, Any]) -> Dict[str, int]:
        """
        No residue space planning needed - no modifications.

        Args:
            components: Matched site components (ignored)

        Returns:
            Empty dict - no residue space allocation needed
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
