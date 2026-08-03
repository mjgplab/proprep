#!/usr/bin/env python3
"""
Bis-Histidine B-Type Heme Transformer for RedoxSite System

This transformer handles bis-histidine ligated B-type hemes (non-covalently bound).
Unlike C-type hemes, B-type hemes have no thioether bonds to cysteine residues.

The two axial histidine ligands are designated 'front' and 'rear' based on the
parameterization model (PDB 1CYO). There is no general structural rule to
distinguish them in an arbitrary protein, so the transformer declares the
front/rear assignment as an ``_ambiguous`` block — the user resolves it
interactively via :py:meth:`RedoxTransformationManager._resolve_ambiguities`.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    RedoxSiteTransformerBase,
    TransformerEvaluation,
    TransformerEvaluationDetail,
    register_redox_transformer
)

from proprep.structure_prep.comprehensive_redox_detector import CenterType

logger = logging.getLogger(__name__)

@register_redox_transformer
class BisHisBTypeHemeTransformer(RedoxSiteTransformerBase):
    """Transformer for bis-histidine ligated B-type hemes (no thioether bonds)"""

    TRANSFORMER_NAME = "heme_bis_his_b_type"
    DESCRIPTION = "Bis-histidine ligated B-type heme (non-covalently bound)"
    SUPPORTED_SITE_TYPES = ["heme_bis_his_b_type"]
    FORCEFIELD_PATH = "heme/bis_his_b_type"

    @classmethod
    def evaluate_redox_site(cls, redox_site) -> TransformerEvaluation:
        """Evaluate site compatibility for bis-his b-type heme"""
        all_details = []
        total_requirements_met = 0
        total_requirements = 0

        # Validate centers
        center_met, center_total, center_details = cls.validate_centers_from_requirements(redox_site)
        all_details.extend(center_details)
        total_requirements_met += center_met
        total_requirements += center_total

        # Validate atoms/residues
        atom_met, atom_total, atom_details = cls.validate_atoms_from_requirements(redox_site)
        all_details.extend(atom_details)
        total_requirements_met += atom_met
        total_requirements += atom_total

        # Validate bonds
        bond_met, bond_total, bond_details = cls.validate_bonds_from_requirements(redox_site)
        all_details.extend(bond_details)
        total_requirements_met += bond_met
        total_requirements += bond_total

        confidence = total_requirements_met / total_requirements if total_requirements > 0 else 0.0
        # Exact criteria — see bis_his_c_type_heme for the full note.
        is_valid = total_requirements_met == total_requirements

        return TransformerEvaluation(
            is_valid=is_valid,
            confidence=confidence,
            description=cls.DESCRIPTION,
            requirements_met=total_requirements_met,
            total_requirements=total_requirements,
            details=all_details
        )

    @classmethod
    def match_components(cls, redox_site) -> Tuple[Dict[str, Any], List[str]]:
        """
        Match components for a bis-his b-type heme site.

        Simpler than c-type: no cysteine identification needed.
        Finds heme center and two coordinated histidines, then assigns
        front/rear by residue ID (lower = front, higher = rear).
        """
        matched_components = {}
        missing_components = []

        # Step 1: Find the heme center
        heme_centers = []
        for c in redox_site.centers:
            if hasattr(c, 'center_type'):
                center_type = c.center_type
                if hasattr(center_type, 'value'):
                    type_value = center_type.value
                else:
                    type_value = str(center_type).lower()

                if type_value in ("organic_cofactor", "organometallic_cofactor"):
                    if hasattr(c, 'resname') and c.resname in ['HEM', 'HEC']:
                        heme_centers.append(c)

        if not heme_centers:
            logger.warning(f"No heme centers found in site {redox_site.site_id}")
            missing_components.extend(["center_id", "center_chain"])
            return matched_components, missing_components

        heme_center = heme_centers[0]
        matched_components["center_id"] = heme_center.resid
        matched_components["center_chain"] = heme_center.chain

        # Step 2: Find coordinated histidine ligands via coordination bonds
        his_ligands = []

        for bond in redox_site.bonds:
            if bond.chemical_type == "coordinate":
                atom1_info = redox_site.coord_to_pdb.get(bond.atom1_coords, {})
                atom2_info = redox_site.coord_to_pdb.get(bond.atom2_coords, {})

                his_atom = None
                if (atom1_info and atom1_info.get('resname') == 'HIS' and
                    atom2_info and atom2_info.get('element') == 'FE'):
                    his_atom = atom1_info
                elif (atom2_info and atom2_info.get('resname') == 'HIS' and
                      atom1_info and atom1_info.get('element') == 'FE'):
                    his_atom = atom2_info

                if his_atom:
                    # Avoid duplicates (same His can have multiple bonds listed)
                    already_found = any(
                        h['resid'] == his_atom['resid'] and h['chain'] == his_atom['chain']
                        for h in his_ligands
                    )
                    if not already_found:
                        his_ligands.append({
                            'chain': his_atom['chain'],
                            'resid': his_atom['resid'],
                            'atom_name': his_atom['atom_name']
                        })

        # Step 3: Defer front/rear assignment to the user via the resolver.
        # No general structural rule distinguishes the two axial His ligands of
        # a bis-His b-type heme, so we list both as candidates and let the
        # transformation manager prompt for the assignment.
        if len(his_ligands) >= 2:
            his_ligands.sort(key=lambda h: (h['chain'], h['resid']))
            candidates = [
                {
                    'chain': h['chain'],
                    'resid': h['resid'],
                    'display': f"HIS {h['chain']}:{h['resid']}",
                }
                for h in his_ligands[:2]
            ]
            matched_components.setdefault("_ambiguous", []).append({
                "label": "Axial histidine assignment",
                "description": (
                    "Bis-His b-type heme is symmetric: which His is 'front' "
                    "(matches FHO/FHR in the parameter set) and which is 'rear' "
                    "(matches RHO/RHR)?"
                ),
                "candidates": candidates,
                "roles": {"front_ligand": 1, "rear_ligand": 1},
            })
        else:
            logger.warning(
                f"Expected 2 coordinated histidines, found {len(his_ligands)} "
                f"in site {redox_site.site_id}"
            )
            missing_components.extend([
                "front_ligand_id", "front_ligand_chain",
                "rear_ligand_id", "rear_ligand_chain"
            ])

        # Step 4: Calculate propionate residue IDs
        if "center_id" in matched_components:
            matched_components["prop_a_id"] = matched_components["center_id"] + 1
            matched_components["prop_a_chain"] = matched_components["center_chain"]
            matched_components["prop_d_id"] = matched_components["center_id"] + 2
            matched_components["prop_d_chain"] = matched_components["center_chain"]

        return matched_components, missing_components

    @classmethod
    def update_components_with_id_mapping(cls, components: Dict[str, Any],
                                        id_mapping: Dict[Tuple[str, int], int]) -> Dict[str, Any]:
        """Recalculate propionate IDs when center ID is remapped"""
        updated_components = super().update_components_with_id_mapping(components, id_mapping)

        if "center_id" in updated_components:
            new_center_id = updated_components["center_id"]
            updated_components["prop_a_id"] = new_center_id + 1
            updated_components["prop_d_id"] = new_center_id + 2

            if "center_chain" in updated_components:
                updated_components["prop_a_chain"] = updated_components["center_chain"]
                updated_components["prop_d_chain"] = updated_components["center_chain"]

        return updated_components

    @classmethod
    def get_required_residue_count(cls) -> int:
        """B-type heme needs 3 residue IDs: heme + 2 propionates"""
        return 3

    @classmethod
    def get_residue_space_plan(cls, components: Dict[str, Any]) -> Dict[str, int]:
        """Define how the 3 residue IDs will be allocated"""
        return {
            "center": 0,
            "propionate_a": 1,
            "propionate_d": 2
        }

    @classmethod
    def get_parameter_definitions(cls) -> Dict[str, Any]:
        """Parameters for b-type hemes"""
        return {
            "redox_state": {
                "description": "Iron oxidation state",
                "type": "choice",
                "options": ["reduced", "oxidized"],
                "default": "reduced"
            },
            "spin_state": {
                "description": "Electronic spin state (bis-histidine coordination enforces low-spin)",
                "type": "fixed",
                "value": "low_spin",
                "note": "Bis-histidine ligated hemes are always low-spin due to strong field ligands"
            }
        }

    @classmethod
    def validate_parameters(cls, parameters: Dict[str, Any]) -> Tuple[bool, str]:
        """Validate parameter values"""
        param_defs = cls.get_parameter_definitions()

        for param_name, param_def in param_defs.items():
            if param_name not in parameters:
                return False, f"Missing required parameter: {param_name}"

            value = parameters[param_name]
            if param_def["type"] == "choice":
                if value not in param_def["options"]:
                    return False, f"Invalid {param_name}: {value}. Options: {param_def['options']}"

        return True, "Parameters valid"

    @classmethod
    def get_parameter_mappings(cls, parameters: Dict[str, Any]) -> Dict[str, str]:
        """Map user parameters to forcefield-specific residue names"""
        redox_state = parameters.get("redox_state", "reduced")
        spin_state = parameters.get("spin_state", "low_spin")

        mappings = {
            "reduced_low_spin": {
                "heme_name": "HBR",
                "front_ligand_name": "FHR",
                "rear_ligand_name": "RHR",
            },
            "oxidized_low_spin": {
                "heme_name": "HBO",
                "front_ligand_name": "FHO",
                "rear_ligand_name": "RHO",
            }
        }

        param_key = f"{redox_state}_{spin_state}"
        return mappings.get(param_key, mappings["reduced_low_spin"])

    @classmethod
    def get_transformation_sequence(cls, components: Dict[str, Any], parameters: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Generate transformation sequence for bis-his b-type heme.

        Simpler than c-type: no cysteine sidechain migration or CYO conversion.
        Steps:
          1. Extract propionate A from HEM into PRN residue
          2. Extract propionate D from HEM into PRN residue
          3. Rename HEM to redox-specific heme name (HBO/HBR)
          4. Rename front histidine to redox-specific name (FHO/FHR)
          5. Rename rear histidine to redox-specific name (RHO/RHR)
        """
        mappings = cls.get_parameter_mappings(parameters)

        transformations = []

        # Step 1: Extract propionate A
        transformations.append({
            "id": "extract_propionate_a",
            "description": "Create PRN residue from propionate A group",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEM",
                "residue_id": components["center_id"],
                "atom_names": ["CAA", "CBA", "CGA", "O1A", "O2A"]
            },
            "action": {
                "change_residue_name": "PRN",
                "change_residue_id": components["prop_a_id"],
                "change_chain_id": components["center_chain"],
                "rename_atoms": {"CAA": "CA", "CBA": "CB", "CGA": "CG", "O1A": "O1", "O2A": "O2"}
            }
        })

        # Step 2: Extract propionate D
        transformations.append({
            "id": "extract_propionate_d",
            "description": "Create PRN residue from propionate D group",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEM",
                "residue_id": components["center_id"],
                "atom_names": ["CAD", "CBD", "CGD", "O1D", "O2D"]
            },
            "action": {
                "change_residue_name": "PRN",
                "change_residue_id": components["prop_d_id"],
                "change_chain_id": components["center_chain"],
                "rename_atoms": {"CAD": "CA", "CBD": "CB", "CGD": "CG", "O1D": "O1", "O2D": "O2"}
            }
        })

        # Step 3: Apply redox-specific heme name
        transformations.append({
            "id": "apply_redox_specific_heme_name",
            "description": f"Rename HEM to {mappings['heme_name']}",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEM",
                "residue_id": components["center_id"]
            },
            "action": {
                "change_residue_name": mappings["heme_name"]
            }
        })

        # Step 4: Rename front histidine
        transformations.append({
            "id": "apply_redox_specific_front_his",
            "description": f"Rename front His to {mappings['front_ligand_name']}",
            "selector": {
                "chain_id": components["front_ligand_chain"],
                "residue_id": components["front_ligand_id"]
            },
            "action": {
                "change_residue_name": mappings["front_ligand_name"]
            }
        })

        # Step 5: Rename rear histidine
        transformations.append({
            "id": "apply_redox_specific_rear_his",
            "description": f"Rename rear His to {mappings['rear_ligand_name']}",
            "selector": {
                "chain_id": components["rear_ligand_chain"],
                "residue_id": components["rear_ligand_id"]
            },
            "action": {
                "change_residue_name": mappings["rear_ligand_name"]
            }
        })

        return transformations

    @classmethod
    def get_site_requirements(cls) -> Dict[str, Any]:
        """
        Site requirements for b-type bis-his heme.

        Key difference from c-type: NO cysteine requirements.
        B-type hemes are held in place by the protein's hydrophobic pocket
        and axial histidine coordination, not thioether bonds.
        """
        return {
            "centers": {
                "required_count": 1,
                "center_types": [CenterType.ORGANOMETALLIC_COFACTOR],
                "residue_names": ["HEM"]
            },
            "atoms": {
                "required_residues": {
                    "HEM": {"min_count": 1, "max_count": 1},
                    "HIS": {"min_count": 2, "max_count": 2},
                    # A b-type heme is held by the pocket and the two axial
                    # His; it has no thioether linkage. Stating the absence
                    # explicitly is what separates this transformer from the
                    # bis-His c-type one when a c-type heme is deposited (or
                    # mis-annotated) as HEM rather than HEC — without it,
                    # every check here passes on a c-type site.
                    "CYS": {"min_count": 0, "max_count": 0}
                }
            },
            "bonds": {
                "required_bond_groups": [
                    {
                        "description": "Fe-His coordination bonds",
                        "min_count": 2,
                        "bond_types": {
                            "coordinate": [
                                (("HEM", "FE"), ("HIS", "NE2")),
                                (("HEM", "FE"), ("HIS", "NE2"))
                            ]
                        }
                    }
                ],
                "require_one_group": True
            }
        }
