#!/usr/bin/env python3
"""Shared base class for organic-cofactor transformers (FMN/FAD/NADP/NAD/H4B
family). These cofactors are single-residue substitutions: the transformer
renames the residue to match the user's chosen redox state, and that's it.
Hydrogen toggles (HN1/HN5/etc.) are deferred to the Hydrogen Addition step
(structure preprocessor 0c) using the new residue's .lib as the template.

Subclasses override class attributes:
    TRANSFORMER_NAME, DESCRIPTION, FORCEFIELD_PATH (required)
    CANONICAL_INPUT_RESNAMES (set of resnames detector reports for this family,
        e.g. {"FMN"} for FlavinFMNTransformer; {"FAD"} for FlavinFADTransformer)
    REDOX_STATE_OPTIONS (list of (redox_state_key, output_resname, description))
    ATOM_NAME_RENAMES (per-state dict of PDB-input-atom-name → lib-atom-name
        renames; identity renames omitted)
"""

import logging
from typing import Any, Dict, List, Tuple

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    RedoxSiteTransformerBase,
    TransformerEvaluation,
    TransformerEvaluationDetail,
)
from proprep.structure_prep.comprehensive_redox_detector import CenterType

logger = logging.getLogger(__name__)


class OrganicCofactorTransformerBase(RedoxSiteTransformerBase):
    """Base for single-residue organic-cofactor transformers."""

    SUPPORTED_SITE_TYPES = ["organic_cofactor"]

    # ── Subclass overrides ────────────────────────────────────────────────
    CANONICAL_INPUT_RESNAMES: set = set()  # e.g. {"FMN"} or {"FAD"}
    # List of (redox_state_key, output_resname, description, is_default).
    # Order in the list defines the order shown to the user.
    REDOX_STATE_OPTIONS: List[Tuple[str, str, str, bool]] = []
    # Per-output-resname atom-name renames applied during transformation.
    # Empty dict means atom names are already canonical.
    ATOM_NAME_RENAMES: Dict[str, Dict[str, str]] = {}

    # ── Standard transformer interface ────────────────────────────────────

    @classmethod
    def evaluate_redox_site(cls, redox_site) -> TransformerEvaluation:
        details = []
        met = 0
        total = 2

        # R1: exactly 1 center
        n = len(redox_site.centers) if hasattr(redox_site, 'centers') else 0
        passed1 = n == 1
        details.append(TransformerEvaluationDetail(
            description="Must have exactly 1 center (single organic cofactor residue)",
            passed=passed1, value_found=n, value_expected=1,
            error_message=None if passed1 else f"Found {n} centers, expected 1",
        ))
        if passed1:
            met += 1

        # R2: center is ORGANIC_COFACTOR with canonical resname
        passed2 = False
        if n >= 1:
            c = redox_site.centers[0]
            ct = getattr(c, 'center_type', None)
            ct_value = ct.value if hasattr(ct, 'value') else str(ct).lower()
            resname_ok = getattr(c, 'resname', None) in cls.CANONICAL_INPUT_RESNAMES
            type_ok = ct_value == 'organic_cofactor'
            passed2 = resname_ok and type_ok
        details.append(TransformerEvaluationDetail(
            description=f"Center must be ORGANIC_COFACTOR with resname in {sorted(cls.CANONICAL_INPUT_RESNAMES)}",
            passed=passed2,
            value_found=getattr(redox_site.centers[0], 'resname', '?') if n >= 1 else 'no center',
            value_expected=sorted(cls.CANONICAL_INPUT_RESNAMES),
        ))
        if passed2:
            met += 1

        return TransformerEvaluation(
            is_valid=(met == total),
            confidence=met / total if total else 0.0,
            description=cls.DESCRIPTION,
            requirements_met=met,
            total_requirements=total,
            details=details,
        )

    @classmethod
    def get_site_requirements(cls) -> Dict[str, Any]:
        return {
            "centers": {
                "required_count": 1,
                "center_types": [CenterType.ORGANIC_COFACTOR],
                "residue_names": sorted(cls.CANONICAL_INPUT_RESNAMES),
            },
            "atoms": {
                "required_residues": {
                    rn: {"min_count": 1, "max_count": 1}
                    for rn in cls.CANONICAL_INPUT_RESNAMES
                },
            },
        }

    @classmethod
    def match_components(cls, redox_site) -> Tuple[Dict[str, Any], List[str]]:
        if not hasattr(redox_site, 'centers') or not redox_site.centers:
            return {}, ["center_id", "center_chain"]
        c = redox_site.centers[0]
        return {
            "center_id": c.resid,
            "center_chain": c.chain,
            "center_resname": c.resname,
        }, []

    @classmethod
    def get_required_residue_count(cls) -> int:
        return 1

    @classmethod
    def get_residue_space_plan(cls, components: Dict[str, Any]) -> Dict[str, int]:
        return {}

    @classmethod
    def get_parameter_definitions(cls) -> Dict[str, Any]:
        opts = [rs for rs, _, _, _ in cls.REDOX_STATE_OPTIONS]
        default = next(
            (rs for rs, _, _, is_default in cls.REDOX_STATE_OPTIONS if is_default),
            opts[0] if opts else None,
        )
        return {
            "redox_state": {
                "description": "Redox/protonation state of the cofactor",
                "type": "choice",
                "options": opts,
                "default": default,
            },
            "spin_state": {
                "description": "Spin state (set by the redox state; not user-configurable)",
                "type": "fixed",
                "value": "default",
            },
        }

    @classmethod
    def get_option_description(cls, param_name: str, option_value: Any) -> str:
        if param_name != "redox_state":
            return ""
        for rs, _resname, descr, _is_default in cls.REDOX_STATE_OPTIONS:
            if rs == option_value:
                return descr
        return ""

    @classmethod
    def validate_parameters(cls, parameters: Dict[str, Any]) -> Tuple[bool, str]:
        if "redox_state" not in parameters:
            return False, "Missing required parameter: redox_state"
        valid_states = [rs for rs, _, _, _ in cls.REDOX_STATE_OPTIONS]
        if parameters["redox_state"] not in valid_states:
            return False, (
                f"Invalid redox_state: {parameters['redox_state']}. "
                f"Options: {valid_states}"
            )
        return True, "Parameters valid"

    @classmethod
    def get_parameter_mappings(cls, parameters: Dict[str, Any]) -> Dict[str, str]:
        redox_state = parameters.get("redox_state")
        for rs, output_resname, _, _ in cls.REDOX_STATE_OPTIONS:
            if rs == redox_state:
                return {"output_residue_name": output_resname}
        # Fallback: first option
        if cls.REDOX_STATE_OPTIONS:
            return {"output_residue_name": cls.REDOX_STATE_OPTIONS[0][1]}
        return {}

    @classmethod
    def get_transformation_sequence(
        cls, components: Dict[str, Any], parameters: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        mappings = cls.get_parameter_mappings(parameters)
        new_resname = mappings.get("output_residue_name")
        if not new_resname:
            return []

        # Identity rename (e.g., FMN → FMN when redox_state=oxidized) is still
        # worth emitting so the transformation manager records the site as
        # transformed and so any atom-name renames in ATOM_NAME_RENAMES still
        # apply. Cost is zero (a no-op residue rename has no side effects).
        old_resname = components.get("center_resname")
        center_chain = components.get("center_chain")
        center_id = components.get("center_id")

        transformations: List[Dict[str, Any]] = []

        # Residue rename
        transformations.append({
            "id": f"rename_{old_resname}_to_{new_resname}",
            "description": f"Rename {old_resname} → {new_resname} (redox state: {parameters.get('redox_state')})",
            "selector": {
                "chain_id": center_chain,
                "residue_name": old_resname,
                "residue_id": center_id,
            },
            "action": {
                "change_residue_name": new_resname,
            },
        })

        # Atom-name renames (if any defined for this output state)
        atom_renames = cls.ATOM_NAME_RENAMES.get(new_resname, {})
        if atom_renames:
            transformations.append({
                "id": f"rename_atoms_in_{new_resname}",
                "description": f"Normalize PDB atom names to library convention for {new_resname}",
                "selector": {
                    "chain_id": center_chain,
                    # Selector now matches the new residue name (after the
                    # rename above has been applied by the manager). All atoms
                    # in the residue are eligible; the rename_atoms map is
                    # keyed by old atom name.
                    "residue_name": new_resname,
                    "residue_id": center_id,
                },
                "action": {
                    "rename_atoms": atom_renames,
                },
            })

        return transformations
