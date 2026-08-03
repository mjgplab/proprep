#!/usr/bin/env python3
"""
Bis-Histidine C-Type Heme Transformer for RedoxSite System

This transformer handles bis-histidine ligated C-type hemes with thioether-linked cysteines.
Part of the general RedoxSite transformation framework.

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

# Import RedoxSite components - must use same CenterType as RedoxSite objects
from proprep.structure_prep.comprehensive_redox_detector import CenterType

logger = logging.getLogger(__name__)

@register_redox_transformer
class BisHisCTypeHemeTransformer(RedoxSiteTransformerBase):
    """Transformer for bis-histidine ligated C-type hemes with thioether bonds to cysteines"""
    
    # Transformer identification
    TRANSFORMER_NAME = "heme_bis_his_c_type"
    DESCRIPTION = "Bis-histidine ligated C-type heme with thioether-linked cysteines"
    SUPPORTED_SITE_TYPES = ["heme_bis_his_c_type"]
    FORCEFIELD_PATH = "heme/bis_his_c_type"
    
    @classmethod
    def evaluate_redox_site(cls, redox_site) -> TransformerEvaluation:
        """Comprehensive evaluation for bis-his c-type heme compatibility using declarative validation"""
        all_details = []
        total_requirements_met = 0
        total_requirements = 0
        
        # Validate centers using generic validation
        center_met, center_total, center_details = cls.validate_centers_from_requirements(redox_site)
        all_details.extend(center_details)
        total_requirements_met += center_met
        total_requirements += center_total
        
        # Validate atoms/residues using generic validation
        atom_met, atom_total, atom_details = cls.validate_atoms_from_requirements(redox_site)
        all_details.extend(atom_details)
        total_requirements_met += atom_met
        total_requirements += atom_total
        
        # Validate bonds using generic validation
        bond_met, bond_total, bond_details = cls.validate_bonds_from_requirements(redox_site)
        all_details.extend(bond_details)
        total_requirements_met += bond_met
        total_requirements += bond_total
        
        # Calculate confidence based on percentage of requirements met
        confidence = total_requirements_met / total_requirements if total_requirements > 0 else 0.0
        
        # Every declared requirement must be met. The axial-ligand composition
        # IS the definition of this site type — a bis-His c-type heme has two
        # His and two thioether Cys, and a site missing either is a different
        # cofactor needing different parameters, not a near miss. A partial-
        # credit threshold let sibling heme transformers pass on each other's
        # sites, because the ligand identity checks are only a couple of votes
        # among many and an absent residue scores a free pass on its max_count
        # check. confidence keeps the ratio for display and ranking.
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
        Match components using RedoxSite bond and atom information
        """
        import logging
        logger = logging.getLogger(__name__)

        matched_components = {}
        missing_components = []

        # Step 1: Find the center (heme organic cofactor)
        # Use value-based comparison to handle enum serialization issues
        heme_centers = []
        for c in redox_site.centers:
            # Debug logging to understand what we're getting
            has_center_type = hasattr(c, 'center_type')
            logger.debug(f"Center {c.resname if hasattr(c, 'resname') else 'unknown'} - has_center_type: {has_center_type}")

            if has_center_type:
                center_type = c.center_type
                has_value = hasattr(center_type, 'value')
                logger.debug(f"  center_type: {center_type}, has_value: {has_value}")

                # Handle both enum objects and string values
                if has_value:
                    type_value = center_type.value
                else:
                    # If it's already a string (shouldn't happen but be defensive)
                    type_value = str(center_type).lower()

                logger.debug(f"  type_value: {type_value}")

                if type_value in ("organic_cofactor", "organometallic_cofactor"):
                    # Verify it's a heme by checking resname
                    if hasattr(c, 'resname') and c.resname in ['HEM', 'HEC']:
                        heme_centers.append(c)
                        logger.debug(f"  Found heme center: {c.resname}")

        if not heme_centers:
            logger.warning(f"No heme centers found in site {redox_site.site_id}. Total centers: {len(redox_site.centers)}")
            for i, c in enumerate(redox_site.centers):
                logger.warning(f"  Center {i}: resname={getattr(c, 'resname', 'N/A')}, "
                             f"center_type={getattr(c, 'center_type', 'N/A')}")
            missing_components.extend(["center_id", "center_chain"])
            return matched_components, missing_components

        heme_center = heme_centers[0]
        matched_components["center_id"] = heme_center.resid
        matched_components["center_chain"] = heme_center.chain
        
        # Step 2: Identify cysteine attachments using distance-based assignment
        # Find CAB and CAC atoms in the heme
        cab_atom = None
        cac_atom = None
        
        for atom in redox_site.atoms:
            if atom.resname in ['HEM', 'HEC']:
                if atom.atom_name == 'CAB':
                    cab_atom = atom
                elif atom.atom_name == 'CAC':
                    cac_atom = atom
        
        
        # Find all CYS SG atoms
        cys_sg_atoms = []
        for atom in redox_site.atoms:
            if atom.resname == 'CYS' and atom.atom_name == 'SG':
                cys_sg_atoms.append(atom)
        
        
        # Step 3: Assign b-ring and c-ring cysteines based on closest distances
        b_ring_cys = None
        c_ring_cys = None
        
        if cab_atom and cys_sg_atoms:
            # Find CYS SG closest to CAB (b-ring)
            min_distance = float('inf')
            closest_cys = None
            
            for cys_atom in cys_sg_atoms:
                dx = cys_atom.coords[0] - cab_atom.coords[0]
                dy = cys_atom.coords[1] - cab_atom.coords[1]
                dz = cys_atom.coords[2] - cab_atom.coords[2]
                distance = (dx*dx + dy*dy + dz*dz) ** 0.5
                
                
                if distance < min_distance:
                    min_distance = distance
                    closest_cys = cys_atom
            
            if closest_cys:
                b_ring_cys = {
                    'cys_chain': closest_cys.chain,
                    'cys_resid': closest_cys.resid,
                    'distance': min_distance
                }
        
        if cac_atom and cys_sg_atoms:
            # Find CYS SG closest to CAC (c-ring)
            min_distance = float('inf')
            closest_cys = None
            
            for cys_atom in cys_sg_atoms:
                dx = cys_atom.coords[0] - cac_atom.coords[0]
                dy = cys_atom.coords[1] - cac_atom.coords[1]
                dz = cys_atom.coords[2] - cac_atom.coords[2]
                distance = (dx*dx + dy*dy + dz*dz) ** 0.5
                
                
                if distance < min_distance:
                    min_distance = distance
                    closest_cys = cys_atom
            
            if closest_cys:
                c_ring_cys = {
                    'cys_chain': closest_cys.chain,
                    'cys_resid': closest_cys.resid,
                    'distance': min_distance
                }
        
        if b_ring_cys:
            matched_components["b_ring_cys_id"] = b_ring_cys['cys_resid']
            matched_components["b_ring_cys_chain"] = b_ring_cys['cys_chain']
        else:
            missing_components.extend(["b_ring_cys_id", "b_ring_cys_chain"])
        
        if c_ring_cys:
            matched_components["c_ring_cys_id"] = c_ring_cys['cys_resid'] 
            matched_components["c_ring_cys_chain"] = c_ring_cys['cys_chain']
        else:
            missing_components.extend(["c_ring_cys_id", "c_ring_cys_chain"])
        
        # Step 4: Identify histidine ligands.
        #
        # We used to look for Fe-NE2 coordinate bonds in redox_site.bonds, but
        # those bonds are no longer required in the user-defined bond set —
        # they become intra-residue after the transformer migrates the
        # imidazole into the heme residue, and the conste heme template
        # already declares them. So we walk site.atoms instead: any HIS
        # residue in the site is an axial ligand by construction (the user
        # added it during template refinement specifically because it
        # coordinates Fe).
        his_ligand_keys = set()
        for atom in redox_site.atoms:
            if atom.resname == 'HIS':
                his_ligand_keys.add((atom.chain, atom.resid))

        his_ligands = [
            {'chain': chain, 'resid': resid, 'atom_name': None}
            for chain, resid in sorted(his_ligand_keys)
        ]
        
        # Step 5: Apply site-specific assignment rules
        if len(his_ligands) >= 2 and c_ring_cys:
            # Rule: proximal His = c-ring Cys + 1
            proximal_candidate_id = c_ring_cys['cys_resid'] + 1
            proximal_candidate_chain = c_ring_cys['cys_chain']
            
            proximal_his = None
            distal_his = None
            
            for his in his_ligands:
                if (his['resid'] == proximal_candidate_id and 
                    his['chain'] == proximal_candidate_chain):
                    proximal_his = his
                else:
                    distal_his = his  # By elimination
            
            if proximal_his:
                matched_components["proximal_ligand_id"] = proximal_his['resid']
                matched_components["proximal_ligand_chain"] = proximal_his['chain']
            else:
                missing_components.extend(["proximal_ligand_id", "proximal_ligand_chain"])
            
            if distal_his:
                matched_components["distal_ligand_id"] = distal_his['resid']
                matched_components["distal_ligand_chain"] = distal_his['chain']
            else:
                missing_components.extend(["distal_ligand_id", "distal_ligand_chain"])
        else:
            missing_components.extend([
                "proximal_ligand_id", "proximal_ligand_chain",
                "distal_ligand_id", "distal_ligand_chain"
            ])
        
        # Step 6: Calculate propionate residue IDs (by convention)
        if "center_id" in matched_components:
            matched_components["prop_a_id"] = matched_components["center_id"] + 1
            matched_components["prop_a_chain"] = matched_components["center_chain"]
            matched_components["prop_d_id"] = matched_components["center_id"] + 2  
            matched_components["prop_d_chain"] = matched_components["center_chain"]
        
        return matched_components, missing_components
    
    @classmethod
    def update_components_with_id_mapping(cls, components: Dict[str, Any], 
                                        id_mapping: Dict[Tuple[str, int], int]) -> Dict[str, Any]:
        """
        Override to handle propionate ID dependencies on center ID mapping.
        
        The propionate IDs are calculated as center_id + 1 and center_id + 2,
        so when the center gets ID mapped, we need to recalculate the propionates
        based on the new center ID, not look for explicit propionate mappings.
        """
        # First, apply the default ID mapping updates
        updated_components = super().update_components_with_id_mapping(components, id_mapping)
        
        # Then, handle the special propionate ID dependencies
        if "center_id" in updated_components:
            # Recalculate propionate IDs based on the potentially mapped center ID
            new_center_id = updated_components["center_id"]
            updated_components["prop_a_id"] = new_center_id + 1
            updated_components["prop_d_id"] = new_center_id + 2
            
            # Propionates use the same chain as the center
            if "center_chain" in updated_components:
                updated_components["prop_a_chain"] = updated_components["center_chain"]
                updated_components["prop_d_chain"] = updated_components["center_chain"]
        
        return updated_components
    
    @classmethod
    def get_required_residue_count(cls) -> int:
        """C-type heme needs 3 residue IDs: heme + 2 propionates"""
        return 3
    
    @classmethod
    def get_residue_space_plan(cls, components: Dict[str, Any]) -> Dict[str, int]:
        """Define how the 3 residue IDs will be allocated"""
        return {
            "center": 0,           # Original heme position
            "propionate_a": 1,     # center_id + 1  
            "propionate_d": 2      # center_id + 2
        }
    
    @classmethod
    def get_parameter_definitions(cls) -> Dict[str, Any]:
        """Parameters specific to c-type hemes"""
        defs = {
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
        # Generic pH-treatment fork (constant_pH PRN vs fixed_pH PRD/PRP) + per-ring
        # protomer choices, derived from this cofactor's metadata. Mirrors the NOS
        # cys_axial transformer; yields {} for any cofactor that declares only one
        # treatment, so behavior is unchanged where no fixed_pH set exists.
        defs.update(cls.protonation_parameter_definitions())
        return defs

    @classmethod
    def validate_parameters(cls, parameters: Dict[str, Any]) -> Tuple[bool, str]:
        """Validate parameter values against get_valid_options so gated params
        (protonation_<role>, valid only under fixed_pH) behave correctly — same
        contract as the cys_axial transformer."""
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
        """Map user parameters to forcefield-specific names"""
        redox_state = parameters.get("redox_state", "reduced")
        spin_state = parameters.get("spin_state", "low_spin")
        
        # Parameter mappings. The conste-style HEH/HCO/HCR lib absorbs both
        # axial His side chains into the heme residue, leaving identical
        # backbone-only HIO stubs at the proximal and distal positions —
        # so both proximal_ligand_name and distal_ligand_name resolve to HIO.
        mappings = {
            "reduced_low_spin": {
                "heme_name": "HCR",
                "proximal_ligand_name": "HIO",
                "distal_ligand_name": "HIO",
            },
            "oxidized_low_spin": {
                "heme_name": "HCO",
                "proximal_ligand_name": "HIO",
                "distal_ligand_name": "HIO",
            }
        }
        
        param_key = f"{redox_state}_{spin_state}"
        return mappings.get(param_key, mappings["reduced_low_spin"])
    
    @classmethod
    def get_transformation_sequence(cls, components: Dict[str, Any], parameters: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Generate the specific transformation sequence for bis-his c-type heme"""
        
        # Get parameter mappings for residue names (heme + His/Cys stubs — these
        # are pH-independent: HCO/HCR, HIO, CYO).
        mappings = cls.get_parameter_mappings(parameters)

        # Propionate residue names are pH-treatment dependent: resolved from
        # metadata by role — PRN (constant_pH) or PRD/PRP per ring (fixed_pH),
        # per the chosen forcefield set's protonation_model + the user's
        # ph_treatment / protonation_<role> choices. Falls back to PRN if a
        # cofactor declares no protonation_model (pre-split behavior).
        resolved = cls.resolve_output_residue_names(parameters)
        prop_a_name = resolved.get("propionate_a", "PRN")
        prop_d_name = resolved.get("propionate_d", "PRN")

        transformations = []
        
        # Step 1: Rename HEM to HEC if needed
        transformations.append({
            "id": "rename_hem_to_hec",
            "description": "Rename HEM to HEC (if needed)",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEM",
                "residue_id": components["center_id"]
            },
            "action": {
                "change_residue_name": "HEC"
            }
        })
        
        # Step 2: Move b-ring cysteine sidechain to heme
        # (conste convention: 4 atoms per Cys absorbed into HEH — CB, HB2, HB3, SG)
        transformations.append({
            "id": "transform_b_ring_cys_sidechain",
            "description": "Move b-ring Cys sidechain atoms (CB/HB2/HB3/SG) to heme",
            "selector": {
                "chain_id": components["b_ring_cys_chain"],
                "residue_name": "CYS",
                "residue_id": components["b_ring_cys_id"],
                "atom_names": ["CB", "HB2", "HB3", "SG"]
            },
            "action": {
                "change_residue_name": "HEC",
                "change_residue_id": components["center_id"],
                "change_chain_id": components["center_chain"],
                "change_insertion_code": "",
                "rename_atoms": {"CB": "CBB2", "HB2": "HB2B", "HB3": "HB3B", "SG": "SGB2"},
                "convert_to_hetatm": True
            }
        })

        # Step 3: Move c-ring cysteine sidechain to heme
        # (conste convention: 4 atoms per Cys absorbed into HEH — CB, HB2, HB3, SG)
        transformations.append({
            "id": "transform_c_ring_cys_sidechain",
            "description": "Move c-ring Cys sidechain atoms (CB/HB2/HB3/SG) to heme",
            "selector": {
                "chain_id": components["c_ring_cys_chain"],
                "residue_name": "CYS",
                "residue_id": components["c_ring_cys_id"],
                "atom_names": ["CB", "HB2", "HB3", "SG"]
            },
            "action": {
                "change_residue_name": "HEC",
                "change_residue_id": components["center_id"],
                "change_chain_id": components["center_chain"],
                "change_insertion_code": "",
                "rename_atoms": {"CB": "CBC1", "HB2": "HB2C", "HB3": "HB3C", "SG": "SGC1"},
                "convert_to_hetatm": True
            }
        })

        # Step 4: Convert b-ring cysteine backbone to CYO
        # (conste convention: 6-atom CYO stub — N, H, CA, HA, C, O)
        transformations.append({
            "id": "convert_b_ring_cys_to_cyo",
            "description": "Change thioether-linked Cys residues to CYO for backbone atoms (b-ring)",
            "selector": {
                "chain_id": components["b_ring_cys_chain"],
                "residue_name": "CYS",
                "residue_id": components["b_ring_cys_id"],
                "atom_names": ["N", "H", "CA", "HA", "C", "O"]
            },
            "action": {
                "change_residue_name": "CYO"
            }
        })

        # Step 5: Convert c-ring cysteine backbone to CYO
        # (conste convention: 6-atom CYO stub — N, H, CA, HA, C, O)
        transformations.append({
            "id": "convert_c_ring_cys_to_cyo",
            "description": "Change thioether-linked Cys residues to CYO for backbone atoms (c-ring)",
            "selector": {
                "chain_id": components["c_ring_cys_chain"],
                "residue_name": "CYS",
                "residue_id": components["c_ring_cys_id"],
                "atom_names": ["N", "H", "CA", "HA", "C", "O"]
            },
            "action": {
                "change_residue_name": "CYO"
            }
        })
        
        # Step 6: Extract propionate A
        transformations.append({
            "id": "extract_propionate_a",
            "description": f"Create {prop_a_name} residue from propionate A group",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEC",
                "residue_id": components["center_id"],
                "atom_names": ["CAA", "CBA", "CGA", "O1A", "O2A"]
            },
            "action": {
                "change_residue_name": prop_a_name,
                "change_residue_id": components["prop_a_id"],
                "change_chain_id": components["center_chain"],
                "rename_atoms": {"CAA": "CA", "CBA": "CB", "CGA": "CG", "O1A": "O1", "O2A": "O2"}
            }
        })

        # Step 7: Extract propionate D
        transformations.append({
            "id": "extract_propionate_d",
            "description": f"Create {prop_d_name} residue from propionate D group",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEC",
                "residue_id": components["center_id"],
                "atom_names": ["CAD", "CBD", "CGD", "O1D", "O2D"]
            },
            "action": {
                "change_residue_name": prop_d_name,
                "change_residue_id": components["prop_d_id"],
                "change_chain_id": components["center_chain"],
                "rename_atoms": {"CAD": "CA", "CBD": "CB", "CGD": "CG", "O1D": "O1", "O2D": "O2"}
            }
        })
        
        # NOTE: the HEC -> redox-specific heme name (HCO/HCR) rename is
        # DEFERRED to the very end. All atom migrations into the heme residue
        # (Cys sidechains above; His sidechains below) MUST happen while the
        # heme residue is still named HEC, otherwise they'd land in a phantom
        # HEC residue with the same residue id as the renamed HCO/HCR.
        #
        # His sidechain atoms migrated into the heme residue (conste convention):
        # CB/HB2/HB3 + the entire imidazole ring (CG/ND1/HD1/CE1/HE1/NE2/CD2/HD2),
        # 11 atoms total per His. Inside HEC they're given the conste atom names
        # CB1/HB21/HB31/CG1/ND11/HD11/CE11/HE11/NE21/CD21/HD21 for the proximal
        # His and CB2/HB22/HB32/CG2/ND12/HD12/CE12/HE12/NE22/CD22/HD22 for the
        # distal His.
        his_sidechain_atoms = [
            "CB", "HB2", "HB3", "CG", "ND1", "HD1",
            "CE1", "HE1", "NE2", "CD2", "HD2",
        ]
        his_backbone_atoms = ["N", "H", "CA", "HA", "C", "O"]

        def his_sidechain_rename_map(idx: int) -> Dict[str, str]:
            i = str(idx)
            return {
                "CB":  f"CB{i}",
                "HB2": f"HB2{i}",
                "HB3": f"HB3{i}",
                "CG":  f"CG{i}",
                "ND1": f"ND1{i}",
                "HD1": f"HD1{i}",
                "CE1": f"CE1{i}",
                "HE1": f"HE1{i}",
                "NE2": f"NE2{i}",
                "CD2": f"CD2{i}",
                "HD2": f"HD2{i}",
            }

        # Step 9a: Move proximal His sidechain to heme
        transformations.append({
            "id": "transform_proximal_his_sidechain",
            "description": "Move proximal His sidechain (CB + imidazole) to heme as CB1/CG1/ND11/HD11/CE11/HE11/NE21/CD21/HD21",
            "selector": {
                "chain_id": components["proximal_ligand_chain"],
                "residue_name": "HIS",
                "residue_id": components["proximal_ligand_id"],
                "atom_names": his_sidechain_atoms,
            },
            "action": {
                "change_residue_name": "HEC",
                "change_residue_id": components["center_id"],
                "change_chain_id": components["center_chain"],
                "change_insertion_code": "",
                "rename_atoms": his_sidechain_rename_map(1),
                "convert_to_hetatm": True,
            }
        })

        # Step 9b: Rename proximal His backbone stub to HIO
        transformations.append({
            "id": "apply_redox_specific_proximal_his",
            "description": "Rename proximal ligand based on redox state",
            "selector": {
                "chain_id": components["proximal_ligand_chain"],
                "residue_name": "HIS",
                "residue_id": components["proximal_ligand_id"],
                "atom_names": his_backbone_atoms,
            },
            "action": {
                "change_residue_name": mappings["proximal_ligand_name"]
            }
        })

        # Step 10a: Move distal His sidechain to heme
        transformations.append({
            "id": "transform_distal_his_sidechain",
            "description": "Move distal His sidechain (CB + imidazole) to heme as CB2/CG2/ND12/HD12/CE12/HE12/NE22/CD22/HD22",
            "selector": {
                "chain_id": components["distal_ligand_chain"],
                "residue_name": "HIS",
                "residue_id": components["distal_ligand_id"],
                "atom_names": his_sidechain_atoms,
            },
            "action": {
                "change_residue_name": "HEC",
                "change_residue_id": components["center_id"],
                "change_chain_id": components["center_chain"],
                "change_insertion_code": "",
                "rename_atoms": his_sidechain_rename_map(2),
                "convert_to_hetatm": True,
            }
        })

        # Step 10b: Rename distal His backbone stub to HIO
        transformations.append({
            "id": "apply_redox_specific_distal_his",
            "description": "Rename distal ligand based on redox state",
            "selector": {
                "chain_id": components["distal_ligand_chain"],
                "residue_name": "HIS",
                "residue_id": components["distal_ligand_id"],
                "atom_names": his_backbone_atoms,
            },
            "action": {
                "change_residue_name": mappings["distal_ligand_name"]
            }
        })

        # Step 11: Final HEC -> HCO/HCR rename. Deferred to the end so all
        # atom migrations into the heme residue (Cys and His sidechains above)
        # land under the HEC name; this single step then promotes the merged
        # residue to its state-specific final name.
        transformations.append({
            "id": "apply_redox_specific_heme_name",
            "description": "Rename HEC to redox-specific heme name",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEC",
                "residue_id": components["center_id"]
            },
            "action": {
                "change_residue_name": mappings["heme_name"]
            }
        })

        return transformations
    
    @classmethod
    def get_site_requirements(cls) -> Dict[str, Any]:
        """Return site requirements"""
        return {
            "centers": {
                "required_count": 1,
                "center_types": [CenterType.ORGANOMETALLIC_COFACTOR],  # Heme contains Fe
                "residue_names": ["HEM", "HEC"]
            },
            "atoms": {
                "required_residues": {
                    "HEM": {"min_count": 0, "max_count": 1},
                    "HEC": {"min_count": 0, "max_count": 1}, 
                    "HIS": {"min_count": 2, "max_count": 2},
                    "CYS": {"min_count": 2, "max_count": 2}
                },
                "alternative_groups": [
                    ["HEM", "HEC"]  # Either HEM or HEC, not both
                ]
            },
            "bonds": {
                # Bond requirements reflect the bonds that will be INTER-residue
                # AFTER the transformer's atom migration (and therefore the bonds
                # the user must define during detection so the tleap input gets
                # explicit cross-residue `bond` directives for them).
                #
                # Each entry below corresponds to one inter-residue bond in the
                # final topology:
                #   HEM/HEC: C2A-CAA, C3D-CAD  -> HCO-PRN_A, HCO-PRN_D  (propionate)
                #   CYS:    CA-CB ×2           -> CYO-HCO              (Cys CA-CB-S sidechain to heme)
                #   HIS:    CA-CB ×2           -> HIO-HCO              (His CA-imidazole sidechain to heme)
                #
                # The Fe-NE2 coordination bonds and CAB/CAC-SG thioether bonds
                # are DELIBERATELY NOT required here: they're inter-residue in
                # the input PDB but the transformer migrates the NE2 / SG atoms
                # into the heme residue, making them intra-residue in the final
                # topology where the conste heme template already declares them.
                # Requiring them would force the user to define bonds that the
                # tleap generator then duplicates against the lib template,
                # triggering tleap's "1-4: cannot add bond" fatal error.
                "required_bond_groups": [
                    {
                        "description": "HEM variant bonds (inter-residue post-transformation)",
                        "min_count": 6,
                        "bond_types": {
                            "covalent": [
                                (("HEM", "C2A"), ("HEM", "CAA")),   # -> HCO-PRN_A
                                (("HEM", "C3D"), ("HEM", "CAD")),   # -> HCO-PRN_D
                                (("CYS", "CA"), ("CYS", "CB")),     # -> CYO-HCO (b-ring Cys)
                                (("CYS", "CA"), ("CYS", "CB")),     # -> CYO-HCO (c-ring Cys)
                                (("HIS", "CA"), ("HIS", "CB")),     # -> HIO-HCO (proximal His)
                                (("HIS", "CA"), ("HIS", "CB"))      # -> HIO-HCO (distal His)
                            ]
                        }
                    },
                    {
                        "description": "HEC variant bonds (inter-residue post-transformation)",
                        "min_count": 6,
                        "bond_types": {
                            "covalent": [
                                (("HEC", "C2A"), ("HEC", "CAA")),
                                (("HEC", "C3D"), ("HEC", "CAD")),
                                (("CYS", "CA"), ("CYS", "CB")),
                                (("CYS", "CA"), ("CYS", "CB")),
                                (("HIS", "CA"), ("HIS", "CB")),
                                (("HIS", "CA"), ("HIS", "CB"))
                            ]
                        }
                    }
                ],
                "require_one_group": True  # At least one group must meet its min_count
            }
        }