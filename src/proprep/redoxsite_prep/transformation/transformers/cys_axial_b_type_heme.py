#!/usr/bin/env python3
"""
Cys-Axial B-Type Heme Transformer for RedoxSite System

Handles b-type hemes with a single proximal Cys thiolate ligand — the
P450 / NOS oxygenase-domain class of hemes. Distinct from bis-histidine
b-type hemes (no axial His, only the proximal Cys), and from c-type hemes
(no covalent thioether links to the protein, only the Fe-S coordination).

Splits the HEM residue into the macrocycle (renamed to HTO ferric / HTR
ferrous) plus two propionate side chains (PRN, same constph residue
the bis-his transformers use). The proximal Cys is renamed to CTO/CTR
to swap in the NOS-fit Cys-thiolate charges. The Fe-Sγ coordination bond
and the two propionate-attachment bonds (HEM.C2A → PRN.CA, HEM.C3D →
PRN.CA) are emitted explicitly per the upstream parameterization handoff.

Parameterized for the NOS oxygenase-domain substrate-bound geometry (rat
nNOS Cys415, PDB 2G6H). Force constants are calibrated to that
conformational basin; for 6-coordinate Fe-O₂/NO/CO states, additional
parameterization would be needed.
"""

import logging
from typing import Any, Dict, List, Tuple

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    RedoxSiteTransformerBase,
    TransformerEvaluation,
    register_redox_transformer,
)
from proprep.structure_prep.comprehensive_redox_detector import CenterType

logger = logging.getLogger(__name__)


@register_redox_transformer
class CysAxialBTypeHemeTransformer(RedoxSiteTransformerBase):
    """Transformer for b-type hemes with a single proximal Cys thiolate ligand."""

    TRANSFORMER_NAME = "heme_cys_axial_b_type"
    DESCRIPTION = (
        "Cys-axial b-type heme (NOS oxygenase / P450-class; "
        "5-coordinate with proximal Cys thiolate)"
    )
    SUPPORTED_SITE_TYPES = ["heme_cys_axial_b_type"]
    FORCEFIELD_PATH = "heme/cys_axial_b_type"

    # ────────────────────────────────────────────────────────────────────
    # Site evaluation — uses the standard requirements-driven validators.
    # ────────────────────────────────────────────────────────────────────

    @classmethod
    def evaluate_redox_site(cls, redox_site) -> TransformerEvaluation:
        all_details = []
        met = 0
        total = 0

        for validator in (
            cls.validate_centers_from_requirements,
            cls.validate_atoms_from_requirements,
            cls.validate_bonds_from_requirements,
        ):
            m, t, details = validator(redox_site)
            all_details.extend(details)
            met += m
            total += t

        confidence = met / total if total else 0.0
        # Exact criteria — see bis_his_c_type_heme for the full note.
        is_valid = met == total

        return TransformerEvaluation(
            is_valid=is_valid,
            confidence=confidence,
            description=cls.DESCRIPTION,
            requirements_met=met,
            total_requirements=total,
            details=all_details,
        )

    # ────────────────────────────────────────────────────────────────────
    # Component matching
    # ────────────────────────────────────────────────────────────────────

    @classmethod
    def match_components(cls, redox_site) -> Tuple[Dict[str, Any], List[str]]:
        """Locate the HEM center and its single proximal CYS ligand.

        Unlike bis-his (2 His → ambiguity needs user resolution), Cys-axial
        is unambiguous: at most one Cys can sit within Fe-Sγ coordination
        distance on either face. We pick the unique Cys whose SG is bonded
        to FE in the redox_site's bond list.
        """
        matched: Dict[str, Any] = {}
        missing: List[str] = []

        # Step 1: heme center
        heme_centers = []
        for c in redox_site.centers:
            if not hasattr(c, "center_type"):
                continue
            ct = c.center_type
            type_value = ct.value if hasattr(ct, "value") else str(ct).lower()
            if type_value in ("organic_cofactor", "organometallic_cofactor"):
                if getattr(c, "resname", None) == "HEM":
                    heme_centers.append(c)

        if not heme_centers:
            logger.warning(f"No HEM center in site {redox_site.site_id}")
            missing.extend(["center_id", "center_chain"])
            return matched, missing

        heme = heme_centers[0]
        matched["center_id"] = heme.resid
        matched["center_chain"] = heme.chain

        # Step 2: proximal Cys via Fe-SG coordination bond
        cys_ligands = []
        for bond in redox_site.bonds:
            if bond.chemical_type != "coordinate":
                continue
            a1 = redox_site.coord_to_pdb.get(bond.atom1_coords, {})
            a2 = redox_site.coord_to_pdb.get(bond.atom2_coords, {})
            cys_atom = None
            if (a1.get("resname") == "CYS" and a1.get("atom_name") == "SG"
                    and a2.get("element") == "FE"):
                cys_atom = a1
            elif (a2.get("resname") == "CYS" and a2.get("atom_name") == "SG"
                    and a1.get("element") == "FE"):
                cys_atom = a2
            if cys_atom is None:
                continue
            already = any(
                c["resid"] == cys_atom["resid"] and c["chain"] == cys_atom["chain"]
                for c in cys_ligands
            )
            if not already:
                cys_ligands.append({
                    "chain": cys_atom["chain"],
                    "resid": cys_atom["resid"],
                })

        if len(cys_ligands) != 1:
            logger.warning(
                f"Expected 1 proximal Cys, found {len(cys_ligands)} in site "
                f"{redox_site.site_id}"
            )
            missing.extend(["proximal_cys_id", "proximal_cys_chain"])
        else:
            matched["proximal_cys_chain"] = cys_ligands[0]["chain"]
            matched["proximal_cys_id"] = cys_ligands[0]["resid"]

        # Step 3: propionate residue IDs (same convention as bis-his)
        matched["prop_a_id"] = matched["center_id"] + 1
        matched["prop_a_chain"] = matched["center_chain"]
        matched["prop_d_id"] = matched["center_id"] + 2
        matched["prop_d_chain"] = matched["center_chain"]

        return matched, missing

    # ────────────────────────────────────────────────────────────────────
    # Residue-space allocation
    # ────────────────────────────────────────────────────────────────────

    @classmethod
    def update_components_with_id_mapping(
        cls, components: Dict[str, Any],
        id_mapping: Dict[Tuple[str, int], int],
    ) -> Dict[str, Any]:
        updated = super().update_components_with_id_mapping(components, id_mapping)
        if "center_id" in updated:
            new_id = updated["center_id"]
            updated["prop_a_id"] = new_id + 1
            updated["prop_d_id"] = new_id + 2
            if "center_chain" in updated:
                updated["prop_a_chain"] = updated["center_chain"]
                updated["prop_d_chain"] = updated["center_chain"]
        return updated

    @classmethod
    def get_required_residue_count(cls) -> int:
        """Cys-axial heme uses 3 new residue IDs: heme + 2 propionates.
        The proximal Cys keeps its original residue ID (only resname changes)."""
        return 3

    @classmethod
    def get_residue_space_plan(cls, components: Dict[str, Any]) -> Dict[str, int]:
        return {
            "center": 0,
            "propionate_a": 1,
            "propionate_d": 2,
        }

    # ────────────────────────────────────────────────────────────────────
    # Parameters
    # ────────────────────────────────────────────────────────────────────

    @classmethod
    def get_parameter_definitions(cls) -> Dict[str, Any]:
        defs = {
            "redox_state": {
                "description": "Iron oxidation state",
                "type": "choice",
                "options": ["oxidized", "reduced"],
                "default": "oxidized",
            },
            "spin_state": {
                "description": (
                    "Electronic spin state (Cys-axial NOS heme is high-spin in "
                    "both ferric and ferrous; only one parameter set per redox state)"
                ),
                "type": "fixed",
                "value": "high_spin",
                "note": (
                    "Ferric is S=5/2 sextet, ferrous is S=2 quintet — both "
                    "high-spin in the substrate-bound active site"
                ),
            },
        }
        # Generic pH-treatment fork (constant_pH vs fixed_pH) + per-propionate
        # protomer choices, derived from this cofactor's metadata. Yields {} until
        # a fixed-pH forcefield set is declared (see Phase 6 metadata patch and
        # RedoxSiteTransformerBase.protonation_parameter_definitions), so pre-patch
        # behavior is unchanged. Added after redox/spin so the manager configures
        # those first (ph_treatment + protonation options depend on them).
        defs.update(cls.protonation_parameter_definitions())
        return defs

    @classmethod
    def validate_parameters(cls, parameters: Dict[str, Any]) -> Tuple[bool, str]:
        # Validate against get_valid_options (not the static options list) so that
        # gated params behave correctly: a protonation_<role> choice is only
        # required/valid under fixed_pH (its valid options are empty otherwise),
        # and its options come from the chosen set's protonation_model — not the
        # intentionally-empty static `options`. Params with no valid options in
        # the current configuration are not applicable and are skipped.
        param_defs = cls.get_parameter_definitions()
        for name, defn in param_defs.items():
            others = {k: v for k, v in parameters.items() if k != name}
            valid = cls.get_valid_options(name, others)
            if not valid:
                continue  # gated off / not applicable in this configuration
            if name not in parameters:
                return False, f"Missing required parameter: {name}"
            if defn["type"] == "choice" and parameters[name] not in valid:
                return False, f"Invalid {name}: {parameters[name]}. Options: {valid}"
        return True, "Parameters valid"

    @classmethod
    def get_parameter_mappings(cls, parameters: Dict[str, Any]) -> Dict[str, str]:
        redox = parameters.get("redox_state", "oxidized")
        mappings = {
            "oxidized": {"heme_name": "HTO", "proximal_cys_name": "CTO"},
            "reduced":  {"heme_name": "HTR", "proximal_cys_name": "CTR"},
        }
        return mappings.get(redox, mappings["oxidized"])

    # ────────────────────────────────────────────────────────────────────
    # Transformation sequence
    # ────────────────────────────────────────────────────────────────────

    @classmethod
    def get_transformation_sequence(
        cls, components: Dict[str, Any], parameters: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Build the 4-step transformation sequence:

          1. Extract propionate A from HEM → new PRN residue at C2A
          2. Extract propionate D from HEM → new PRN residue at C3D
          3. Rename HEM → HTO/HTR (state-specific)
          4. Rename proximal CYS → CTO/CTR (state-specific)

        The propionate-extraction atom set matches the bis-his transformer
        (CAA/CBA/CGA/O1A/O2A and CAD/CBD/CGD/O1D/O2D). The HEM's bare
        β-carbons C2A and C3D stay in the macrocycle as sp2 CB atoms (the
        NOS parameterization re-types them so they're ready for cross-
        residue propionate bonds via PRN.CA).
        """
        # All output residue names come from metadata, keyed by role (no hardcoded
        # residue codes). center → HTO/HTR, proximal_cys → CTO/CTR, and the two
        # propionate roles → PRN (constant-pH) or PRD/PRP per ring (fixed-pH),
        # per the chosen forcefield set's protonation_model + the user's
        # ph_treatment / protonation_<role> choices. Structural details below
        # (atom selection, atom renames, residue-IDs, chain-IDs, bonds) are
        # unchanged — only the residue-name *values* are resolved.
        resolved = cls.resolve_output_residue_names(parameters)
        prop_a_name = resolved["propionate_a"]
        prop_d_name = resolved["propionate_d"]
        heme_name = resolved["center"]
        proximal_cys_name = resolved["proximal_cys"]
        transformations: List[Dict[str, Any]] = []

        # Step 1: propionate A (pyrrole A side chain → propionate residue)
        transformations.append({
            "id": "extract_propionate_a",
            "description": f"Create {prop_a_name} residue from propionate A group (HEM C2A side chain)",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEM",
                "residue_id": components["center_id"],
                "atom_names": ["CAA", "CBA", "CGA", "O1A", "O2A"],
            },
            "action": {
                "change_residue_name": prop_a_name,
                "change_residue_id": components["prop_a_id"],
                "change_chain_id": components["center_chain"],
                "rename_atoms": {
                    "CAA": "CA", "CBA": "CB", "CGA": "CG",
                    "O1A": "O1", "O2A": "O2",
                },
            },
        })

        # Step 2: propionate D (pyrrole D side chain → propionate residue)
        transformations.append({
            "id": "extract_propionate_d",
            "description": f"Create {prop_d_name} residue from propionate D group (HEM C3D side chain)",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEM",
                "residue_id": components["center_id"],
                "atom_names": ["CAD", "CBD", "CGD", "O1D", "O2D"],
            },
            "action": {
                "change_residue_name": prop_d_name,
                "change_residue_id": components["prop_d_id"],
                "change_chain_id": components["center_chain"],
                "rename_atoms": {
                    "CAD": "CA", "CBD": "CB", "CGD": "CG",
                    "O1D": "O1", "O2D": "O2",
                },
            },
        })

        # Step 3: HEM → HTO/HTR
        transformations.append({
            "id": "rename_heme_redox_specific",
            "description": f"Rename HEM to {heme_name}",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEM",
                "residue_id": components["center_id"],
            },
            "action": {"change_residue_name": heme_name},
        })

        # Step 4: proximal CYS → CTO/CTR
        if "proximal_cys_id" in components and "proximal_cys_chain" in components:
            transformations.append({
                "id": "rename_proximal_cys",
                "description": (
                    f"Rename proximal CYS at "
                    f"{components['proximal_cys_chain']}:"
                    f"{components['proximal_cys_id']} to {proximal_cys_name}"
                ),
                "selector": {
                    "chain_id": components["proximal_cys_chain"],
                    "residue_name": "CYS",
                    "residue_id": components["proximal_cys_id"],
                },
                "action": {"change_residue_name": proximal_cys_name},
            })

        return transformations

    # ────────────────────────────────────────────────────────────────────
    # Site-requirements declaration
    # ────────────────────────────────────────────────────────────────────

    @classmethod
    def get_site_requirements(cls) -> Dict[str, Any]:
        """Cys-axial b-type heme:
          * 1 HEM center
          * Exactly 1 CYS within Fe-Sγ coordination distance
          * 1 Fe-SG coordinate bond
          * 2 propionate carve-out bonds (C2A-CAA, C3D-CAD); intra-HEM in
            the input PDB but inter-residue after the transformer extracts
            CAA/CAD into the PRN side residues."""
        return {
            "centers": {
                "required_count": 1,
                "center_types": [CenterType.ORGANOMETALLIC_COFACTOR],
                "residue_names": ["HEM"],
            },
            "atoms": {
                "required_residues": {
                    "HEM": {"min_count": 1, "max_count": 1},
                    "CYS": {"min_count": 1, "max_count": 1},
                },
            },
            "bonds": {
                # Bond requirements reflect bonds that will be INTER-residue
                # AFTER the transformer's atom migration — the bonds the user
                # must define during detection so the tleap input gets explicit
                # cross-residue `bond` directives for them.
                #
                # Inter-residue bonds in the final topology:
                #   HEM: C2A-CAA, C3D-CAD  -> HTO/HTR - PRN_A/PRN_D  (propionate)
                #   HEM-CYS: FE-SG         -> HTO/HTR - CTO/CTR      (axial coordination)
                #
                # The Cys is NOT absorbed into the heme residue here (unlike the
                # bis-his c-type), so CB-SG stays intra-CTO/CTR and is declared
                # by the CTO/CTR lib template — no requirement needed for it.
                "required_bond_groups": [
                    {
                        "description": (
                            "HEM bonds (inter-residue post-transformation): "
                            "Fe-Sγ axial coordination + 2 propionate carve-outs"
                        ),
                        "min_count": 3,
                        "bond_types": {
                            "coordinate": [
                                (("HEM", "FE"), ("CYS", "SG")),   # -> HTO/HTR - CTO/CTR
                            ],
                            "covalent": [
                                (("HEM", "C2A"), ("HEM", "CAA")),  # -> HTO/HTR - PRN_A
                                (("HEM", "C3D"), ("HEM", "CAD")),  # -> HTO/HTR - PRN_D
                            ],
                        },
                    },
                ],
            },
        }
