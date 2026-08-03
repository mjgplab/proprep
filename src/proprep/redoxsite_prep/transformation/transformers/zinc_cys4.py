#!/usr/bin/env python3
"""
Zn(Cys)4 Tetrathiolate Transformer for RedoxSite System

Handles tetrahedral Zn(II) sites coordinated by 4 cysteinate ligands
(structural zinc fingers, zinc-binding motifs, etc.). Single state, no
redox or spin choice. All 4 Cys ligands are renamed to a single canonical
CYZ residue (averaged-symmetry parameters); the Zn residue is renamed to
ZNC. The bundled parameters were RESP-fit with cross-Cys side-chain
equivalence so the 4 ligands are chemically interchangeable by
construction.
"""

import logging
from typing import Any, Dict, List, Tuple

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    RedoxSiteTransformerBase,
    TransformerEvaluation,
    TransformerEvaluationDetail,
    register_redox_transformer,
)

logger = logging.getLogger(__name__)


# Target canonical residue names emitted by this transformer.
ZN_TARGET_RESNAME = "ZNC"
CYS_TARGET_RESNAME = "CYZ"


@register_redox_transformer
class ZincCys4Transformer(RedoxSiteTransformerBase):
    """Transformer for tetrahedral Zn(II)(Cys-)4 structural sites."""

    TRANSFORMER_NAME = "zinc_cys4"
    DESCRIPTION = "Tetrahedral Zn(II) coordinated by 4 cysteinate ligands (structural zinc-thiolate site)"
    SUPPORTED_SITE_TYPES = ["zinc_cys4"]
    FORCEFIELD_PATH = "zinc/cys4"

    @classmethod
    def evaluate_redox_site(cls, redox_site) -> TransformerEvaluation:
        """Validate that the site contains 1 Zn + 4 ligating CYS via coordinate Zn-SG bonds."""
        details: List[TransformerEvaluationDetail] = []
        requirements_met = 0
        total_requirements = 3

        zn_atoms = [a for a in redox_site.atoms
                    if getattr(a, "element", None) and a.element.upper() == "ZN"]
        cys_sg_atoms = [a for a in redox_site.atoms
                        if a.resname == "CYS" and a.atom_name == "SG"]
        cys_residues = {(a.chain, a.resid) for a in redox_site.atoms if a.resname == "CYS"}

        # 1. Exactly 1 Zn atom
        ok = len(zn_atoms) == 1
        details.append(TransformerEvaluationDetail(
            description="1 Zn atom in site",
            passed=ok, value_found=len(zn_atoms), value_expected=1,
        ))
        if ok:
            requirements_met += 1

        # 2. 4 ligating CYS residues with SG
        ok = len(cys_residues) == 4 and len(cys_sg_atoms) == 4
        details.append(TransformerEvaluationDetail(
            description="4 CYS residues with SG",
            passed=ok, value_found=len(cys_residues), value_expected=4,
        ))
        if ok:
            requirements_met += 1

        # 3. At least 4 coordinate Zn-SG bonds
        coord_zn_sg = 0
        for bond in redox_site.bonds:
            if getattr(bond, "chemical_type", None) != "coordinate":
                continue
            a1 = redox_site.coord_to_pdb.get(bond.atom1_coords, {}) or {}
            a2 = redox_site.coord_to_pdb.get(bond.atom2_coords, {}) or {}
            zn_side = sg_side = None
            if (a1.get("element") or "").upper() == "ZN" and a2.get("resname") == "CYS" and a2.get("atom_name") == "SG":
                zn_side, sg_side = a1, a2
            elif (a2.get("element") or "").upper() == "ZN" and a1.get("resname") == "CYS" and a1.get("atom_name") == "SG":
                zn_side, sg_side = a2, a1
            if zn_side and sg_side:
                coord_zn_sg += 1

        ok = coord_zn_sg >= 4
        details.append(TransformerEvaluationDetail(
            description="4 coordinate Zn-SG bonds",
            passed=ok, value_found=coord_zn_sg, value_expected=4,
        ))
        if ok:
            requirements_met += 1

        confidence = requirements_met / total_requirements if total_requirements else 0.0
        return TransformerEvaluation(
            is_valid=requirements_met == total_requirements,
            confidence=confidence,
            description=cls.DESCRIPTION,
            requirements_met=requirements_met,
            total_requirements=total_requirements,
            details=details,
        )

    @classmethod
    def match_components(cls, redox_site) -> Tuple[Dict[str, Any], List[str]]:
        """Collect the Zn residue and the 4 ligating Cys (chain, resid) tuples.

        All 4 Cys are chemically equivalent under our averaged-symmetry
        parameters, so order does not matter — we just collect them.
        """
        matched: Dict[str, Any] = {}
        missing: List[str] = []

        zn_atoms = [a for a in redox_site.atoms
                    if getattr(a, "element", None) and a.element.upper() == "ZN"]
        if len(zn_atoms) != 1:
            missing.append("zn_atom")
            return matched, missing
        zn = zn_atoms[0]
        matched["zn_chain"] = zn.chain
        matched["zn_id"] = zn.resid
        matched["zn_source_resname"] = zn.resname

        cys_residues = sorted({(a.chain, a.resid) for a in redox_site.atoms if a.resname == "CYS"})
        if len(cys_residues) != 4:
            missing.append("4_cys_residues")
            return matched, missing
        for n, (chain, resid) in enumerate(cys_residues, start=1):
            matched[f"cys_{n}_chain"] = chain
            matched[f"cys_{n}_id"] = resid

        return matched, missing

    @classmethod
    def get_required_residue_count(cls) -> int:
        """Zn residue keeps its ID; 4 ligating Cys keep theirs. 1 slot for the cluster center."""
        return 1

    @classmethod
    def get_residue_space_plan(cls, components: Dict[str, Any]) -> Dict[str, int]:
        return {"center": 0}

    @classmethod
    def update_components_with_id_mapping(
        cls, components: Dict[str, Any],
        id_mapping: Dict[Tuple[str, int], int],
    ) -> Dict[str, Any]:
        """Remap Zn(Cys)4-specific component IDs after the ID-mapping step.

        The base class implementation only knows about the standard
        ``center_id`` / ``site_id`` / ``main_residue_id`` naming. This
        transformer uses ``zn_id`` for the Zn ion and ``cys_<n>_id`` for the
        four ligating cysteines, so the base implementation silently no-ops
        and the components dict keeps the pre-mapping IDs. Downstream, the
        rename selector then targets a residue that doesn't exist post-
        mapping, the framework's coordinate-match finds nothing, and the
        Zn rename silently no-ops — leaving the residue as ``ZN`` in the
        final PDB. tleap then matches it to ions12.lib's built-in ``ZN``
        unit (atom type ``Zn2+``), and every Zn-S bond/angle parameter
        lookup fails with ``Could not find ... SZ - Zn2+``.

        Override remaps both the Zn ion and the 4 Cys ligands if they're in
        the id_mapping (typically the Zn is, the Cys ligands aren't —
        but the lookup is safe either way).
        """
        updated = super().update_components_with_id_mapping(components, id_mapping)

        zn_chain = updated.get("zn_chain")
        zn_id = updated.get("zn_id")
        if zn_chain is not None and zn_id is not None:
            key = (zn_chain, zn_id)
            if key in id_mapping:
                updated["zn_id"] = id_mapping[key]

        for n in (1, 2, 3, 4):
            cys_chain = updated.get(f"cys_{n}_chain")
            cys_id = updated.get(f"cys_{n}_id")
            if cys_chain is not None and cys_id is not None:
                key = (cys_chain, cys_id)
                if key in id_mapping:
                    updated[f"cys_{n}_id"] = id_mapping[key]

        return updated

    @classmethod
    def get_parameter_definitions(cls) -> Dict[str, Any]:
        """Tetrahedral Zn(Cys)4 has no redox or spin variation (closed-shell d10
        singlet, charge -2 always). Still declare redox_state and spin_state as
        fixed parameters so they thread through the workflow's parameters dict
        — downstream code (FF picker, atom-types resolver) reads
        parameters['redox_state'] / parameters['spin_state'] and falls back to
        'unknown' if they're missing, which then breaks the
        `<cofactor_path>/<redox_state>/<spin_state>/` lookup in the loader.
        """
        return {
            "redox_state": {
                "description": "Redox state (Zn(Cys)4 has only one)",
                "type": "fixed",
                "value": "single_state",
                "note": (
                    "Zn(II) d10 closed shell; not redox-active. Single "
                    "parameter set covers all tetrahedral Zn(Cys)4 sites."
                ),
            },
            "spin_state": {
                "description": "Electronic spin state (always singlet)",
                "type": "fixed",
                "value": "default",
                "note": "Closed-shell d10 singlet (S=0, multiplicity=1).",
            },
        }

    @classmethod
    def validate_parameters(cls, parameters: Dict[str, Any]) -> Tuple[bool, str]:
        return True, "No parameters required"

    @classmethod
    def get_parameter_mappings(cls, parameters: Dict[str, Any]) -> Dict[str, str]:
        return {
            "zn_residue_name": ZN_TARGET_RESNAME,
            "cys_residue_name": CYS_TARGET_RESNAME,
        }

    @classmethod
    def get_transformation_sequence(
        cls, components: Dict[str, Any], parameters: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        mappings = cls.get_parameter_mappings(parameters)
        transformations: List[Dict[str, Any]] = []

        # 1. Rename the Zn residue to ZNC
        zn_selector: Dict[str, Any] = {
            "chain_id": components["zn_chain"],
            "residue_id": components["zn_id"],
        }
        source_resname = components.get("zn_source_resname")
        if source_resname:
            zn_selector["residue_name"] = source_resname
        transformations.append({
            "id": "rename_zn_residue",
            "description": f"Rename {source_resname or 'Zn'} residue to {mappings['zn_residue_name']}",
            "selector": zn_selector,
            "action": {"change_residue_name": mappings["zn_residue_name"]},
        })

        # 2-5. Rename each ligating Cys to CYZ (all 4 get the same name).
        for n in (1, 2, 3, 4):
            transformations.append({
                "id": f"rename_ligating_cys_{n}",
                "description": (
                    f"Rename CYS at {components[f'cys_{n}_chain']}:"
                    f"{components[f'cys_{n}_id']} (ligating Zn) to {mappings['cys_residue_name']}"
                ),
                "selector": {
                    "chain_id": components[f"cys_{n}_chain"],
                    "residue_name": "CYS",
                    "residue_id": components[f"cys_{n}_id"],
                },
                "action": {"change_residue_name": mappings["cys_residue_name"]},
            })

        return transformations

    @classmethod
    def get_site_requirements(cls) -> Dict[str, Any]:
        return {
            "centers": {
                "required_count": 1,
                "elements": ["Zn"],
            },
            "atoms": {
                "required_residues": {
                    "CYS": {"min_count": 4, "max_count": 4},
                },
            },
            "bonds": {
                "required_bond_groups": [
                    {
                        "description": "4 coordinate Zn-SG bonds (ligating Cys)",
                        "min_count": 4,
                        "bond_types": {
                            "coordinate": [
                                (("ANY", "ZN"), ("CYS", "SG")),
                            ],
                        },
                    },
                ],
            },
        }
