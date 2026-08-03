#!/usr/bin/env python3
"""
4Fe-4S Iron-Sulfur Cluster Transformer for RedoxSite System

Handles cubane [4Fe-4S] clusters with 4 cysteinate ligands. Supports three
oxidation states (HiPIP-ox / HiPIP-red == Fd-ox / Fd-red) with multiple
broken-symmetry parameter sets per state. The 4 cysteinates are matched to
the 4 cluster Fe atoms by a deterministic distance rule: for each FE<n> in
the cluster residue, the closest CYS SG becomes the n-th modified
cysteinate residue (e.g. FO1: FE1 -> C01, FE2 -> C02, FE3 -> C03, FE4 -> C04).
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    RedoxSiteTransformerBase,
    TransformerEvaluation,
    TransformerEvaluationDetail,
    register_redox_transformer,
)

logger = logging.getLogger(__name__)


# Source PDB residue names recognized as 4Fe-4S clusters.
SUPPORTED_SOURCE_RESNAMES = {"SF4", "F4S", "FS4", "CFM", "ICS"}

# Cluster atoms expected inside the source residue.
EXPECTED_FE_NAMES = ["FE1", "FE2", "FE3", "FE4"]
EXPECTED_S_NAMES = ["S1", "S2", "S3", "S4"]


# Per-set Cys ligand residue names. Keyed by spin_state code (which equals
# the cluster residue name in metadata.json). Source: README in
# forcefield_params/specialized_residues/iron_sulfur/4fe4s/.
LIGAND_CYS_RESIDUES: Dict[str, List[str]] = {
    # Oxidized HiPIP (FO1-FO6); averaged set (FO0) shares one ligand resname across all 4 positions
    "FO1": ["C01", "C02", "C03", "C04"],
    "FO2": ["C05", "C06", "C07", "C08"],
    "FO3": ["C09", "C10", "C11", "C12"],
    "FO4": ["C13", "C14", "C15", "C16"],
    "FO5": ["C17", "C18", "C19", "C20"],
    "FO6": ["C21", "C22", "C23", "C24"],
    "FO0": ["CXO", "CXO", "CXO", "CXO"],
    # Reduced HiPIP / Oxidized Fd (FR1-FR3); averaged set (FR0)
    "FR1": ["C25", "C26", "C27", "C28"],
    "FR2": ["C29", "C30", "C31", "C32"],
    "FR3": ["C33", "C34", "C35", "C36"],
    "FR0": ["CXR", "CXR", "CXR", "CXR"],
    # Reduced Fd (FD1-FD6); averaged set (FD0)
    "FD1": ["C37", "C38", "C39", "C40"],
    "FD2": ["C41", "C42", "C43", "C44"],
    "FD3": ["C45", "C46", "C47", "C48"],
    "FD4": ["C49", "C50", "C51", "C52"],
    "FD5": ["C53", "C54", "C55", "C56"],
    "FD6": ["C57", "C58", "C59", "C60"],
    "FD0": ["CXD", "CXD", "CXD", "CXD"],
}

# Allowed (redox_state, spin_state) combinations. Code "F<state>0" denotes the
# averaged set (single atom-type-collapsed parameterization across BS guesses).
VALID_STATES: Dict[str, List[str]] = {
    "oxidized":      ["FO1", "FO2", "FO3", "FO4", "FO5", "FO6", "FO0"],
    "reduced_hipip": ["FR1", "FR2", "FR3", "FR0"],
    "reduced_fd":    ["FD1", "FD2", "FD3", "FD4", "FD5", "FD6", "FD0"],
}


def _distance_sq(a, b) -> float:
    dx = a.coords[0] - b.coords[0]
    dy = a.coords[1] - b.coords[1]
    dz = a.coords[2] - b.coords[2]
    return dx * dx + dy * dy + dz * dz


@register_redox_transformer
class IronSulfur4Fe4STransformer(RedoxSiteTransformerBase):
    """Transformer for cubane [4Fe-4S] clusters with 4 cysteinate ligands."""

    TRANSFORMER_NAME = "iron_sulfur_4fe4s"
    DESCRIPTION = "4Fe-4S cubane cluster with 4 cysteinate ligands (HiPIP / ferredoxin)"
    SUPPORTED_SITE_TYPES = ["iron_sulfur_4fe4s", "Fe4S4 cluster"]
    FORCEFIELD_PATH = "iron_sulfur/4fe4s"

    @classmethod
    def evaluate_redox_site(cls, redox_site) -> TransformerEvaluation:
        """Validate that the site contains a 4Fe-4S cluster + 4 ligating CYS."""
        details: List[TransformerEvaluationDetail] = []
        requirements_met = 0
        total_requirements = 4

        fe_atoms = [a for a in redox_site.atoms
                    if getattr(a, "element", None) and a.element.upper() == "FE"]
        s_inorganic = [a for a in redox_site.atoms
                       if getattr(a, "element", None) and a.element.upper() == "S"
                       and a.resname != "CYS"]
        cys_residues = {(a.chain, a.resid) for a in redox_site.atoms if a.resname == "CYS"}
        cys_sg_atoms = [a for a in redox_site.atoms
                        if a.resname == "CYS" and a.atom_name == "SG"]

        # 1. Exactly 4 Fe atoms in the site
        ok = len(fe_atoms) == 4
        details.append(TransformerEvaluationDetail(
            description="4 Fe atoms in cluster",
            passed=ok, value_found=len(fe_atoms), value_expected=4,
        ))
        if ok:
            requirements_met += 1

        # 2. 4 inorganic (non-Cys) sulfur atoms (the bridging sulfides)
        ok = len(s_inorganic) == 4
        details.append(TransformerEvaluationDetail(
            description="4 bridging sulfide atoms",
            passed=ok, value_found=len(s_inorganic), value_expected=4,
        ))
        if ok:
            requirements_met += 1

        # 3. 4 ligating CYS residues with SG
        ok = len(cys_residues) == 4 and len(cys_sg_atoms) == 4
        details.append(TransformerEvaluationDetail(
            description="4 CYS residues with SG",
            passed=ok, value_found=len(cys_residues), value_expected=4,
        ))
        if ok:
            requirements_met += 1

        # 4. At least 4 coordinate Fe-SG bonds
        coord_fe_sg = 0
        for bond in redox_site.bonds:
            if getattr(bond, "chemical_type", None) != "coordinate":
                continue
            a1 = redox_site.coord_to_pdb.get(bond.atom1_coords, {}) or {}
            a2 = redox_site.coord_to_pdb.get(bond.atom2_coords, {}) or {}
            fe_side = sg_side = None
            if (a1.get("element") or "").upper() == "FE" and a2.get("resname") == "CYS" and a2.get("atom_name") == "SG":
                fe_side, sg_side = a1, a2
            elif (a2.get("element") or "").upper() == "FE" and a1.get("resname") == "CYS" and a1.get("atom_name") == "SG":
                fe_side, sg_side = a2, a1
            if fe_side and sg_side:
                coord_fe_sg += 1

        ok = coord_fe_sg >= 4
        details.append(TransformerEvaluationDetail(
            description="4 coordinate Fe-SG bonds",
            passed=ok, value_found=coord_fe_sg, value_expected=4,
        ))
        if ok:
            requirements_met += 1

        confidence = requirements_met / total_requirements if total_requirements else 0.0
        is_valid = requirements_met == total_requirements
        return TransformerEvaluation(
            is_valid=is_valid,
            confidence=confidence,
            description=cls.DESCRIPTION,
            requirements_met=requirements_met,
            total_requirements=total_requirements,
            details=details,
        )

    @classmethod
    def match_components(cls, redox_site) -> Tuple[Dict[str, Any], List[str]]:
        """Match cluster + 4 ligating Cys deterministically.

        For each cluster Fe (named FE1-FE4), find the closest CYS SG. That
        Cys becomes the n-th ligand position. Each Cys is assigned to at
        most one Fe.
        """
        matched: Dict[str, Any] = {}
        missing: List[str] = []

        # Find Fe atoms in the site, grouped by atom name.
        fe_by_name = {}
        for atom in redox_site.atoms:
            if getattr(atom, "element", None) and atom.element.upper() == "FE":
                fe_by_name.setdefault(atom.atom_name, atom)

        # Identify the cluster residue (residue containing the Fe atoms).
        # The expected case is a single residue (SF4/F4S) with all 4 Fe.
        fe_residues = {(a.chain, a.resid, a.resname) for a in fe_by_name.values()}
        if len(fe_residues) == 1:
            chain, resid, resname = next(iter(fe_residues))
            matched["center_chain"] = chain
            matched["center_id"] = resid
            matched["cluster_source_resname"] = resname
        else:
            logger.warning(
                "Fe4S4 site has %d distinct residues holding Fe atoms; "
                "expected 1 merged cluster residue. Residues: %s",
                len(fe_residues), fe_residues,
            )
            missing.extend(["center_id", "center_chain"])
            return matched, missing

        # Verify the cluster has the canonical FE1-FE4 atom names.
        if not all(name in fe_by_name for name in EXPECTED_FE_NAMES):
            logger.warning(
                "Cluster missing expected Fe atom names. Found: %s, expected: %s",
                sorted(fe_by_name.keys()), EXPECTED_FE_NAMES,
            )
            missing.append("expected_fe_atom_names")
            return matched, missing

        # Collect candidate Cys SG atoms.
        cys_sg = [a for a in redox_site.atoms
                  if a.resname == "CYS" and a.atom_name == "SG"]
        if len(cys_sg) < 4:
            logger.warning("Site has %d CYS SG atoms; expected 4.", len(cys_sg))
            missing.extend([f"cys_{n}_id" for n in (1, 2, 3, 4)])
            missing.extend([f"cys_{n}_chain" for n in (1, 2, 3, 4)])
            return matched, missing

        # For each Fe, find the closest unassigned CYS SG.
        used: set = set()
        for n in (1, 2, 3, 4):
            fe = fe_by_name[f"FE{n}"]
            best = None
            best_d2 = None
            for sg in cys_sg:
                key = (sg.chain, sg.resid)
                if key in used:
                    continue
                d2 = _distance_sq(fe, sg)
                if best is None or d2 < best_d2:
                    best, best_d2 = sg, d2
            if best is None:
                missing.append(f"cys_{n}_id")
                missing.append(f"cys_{n}_chain")
                continue
            matched[f"cys_{n}_id"] = best.resid
            matched[f"cys_{n}_chain"] = best.chain
            used.add((best.chain, best.resid))

        return matched, missing

    @classmethod
    def get_required_residue_count(cls) -> int:
        """The cluster keeps its existing ID; ligating Cys keep theirs. 1 is enough."""
        return 1

    @classmethod
    def get_residue_space_plan(cls, components: Dict[str, Any]) -> Dict[str, int]:
        return {"center": 0}

    @classmethod
    def get_parameter_definitions(cls) -> Dict[str, Any]:
        return {
            "redox_state": {
                "description": "Cluster oxidation state",
                "type": "choice",
                "options": list(VALID_STATES.keys()),
                "default": "reduced_hipip",
            },
            "spin_state": {
                "description": (
                    "Broken-symmetry parameter set. Valid codes per redox state: "
                    "oxidized=FO1-FO6, reduced_hipip=FR1-FR3, reduced_fd=FD1-FD6."
                ),
                "type": "choice",
                "options": [code for codes in VALID_STATES.values() for code in codes],
                "default": "FR1",
            },
        }

    @classmethod
    def validate_parameters(cls, parameters: Dict[str, Any]) -> Tuple[bool, str]:
        redox = parameters.get("redox_state")
        spin = parameters.get("spin_state")
        if redox is None or spin is None:
            return False, "Both 'redox_state' and 'spin_state' are required"
        if redox not in VALID_STATES:
            return False, f"Invalid redox_state '{redox}'. Options: {list(VALID_STATES)}"
        if spin not in VALID_STATES[redox]:
            return False, (
                f"spin_state '{spin}' is not valid for redox_state '{redox}'. "
                f"Valid codes: {VALID_STATES[redox]}"
            )
        if spin not in LIGAND_CYS_RESIDUES:
            return False, f"No ligand-Cys mapping defined for spin_state '{spin}'"
        return True, "Parameters valid"

    @classmethod
    def get_valid_options(cls, param_name: str,
                          current_parameters: Dict[str, Any]) -> List[Any]:
        """Filter spin_state options by the redox_state already chosen."""
        if param_name == "spin_state":
            redox = current_parameters.get("redox_state")
            if redox in VALID_STATES:
                return list(VALID_STATES[redox])
            # If redox_state hasn't been chosen yet, return the full union so
            # the manager has something to render. The manager should always
            # configure redox_state before spin_state given dict insertion
            # order in get_parameter_definitions.
            return [code for codes in VALID_STATES.values() for code in codes]
        return super().get_valid_options(param_name, current_parameters)

    @classmethod
    def get_option_description(cls, param_name: str, option_value: Any) -> str:
        """Surface descriptions from the Fe4S4 metadata.json."""
        try:
            from proprep.forcefield_params import load_forcefield_metadata
            md = load_forcefield_metadata(cls.FORCEFIELD_PATH)
        except Exception:
            return ""

        if param_name == "redox_state":
            return md.get("redox_states", {}).get(option_value, {}).get("description", "")

        if param_name == "spin_state":
            for redox_data in md.get("redox_states", {}).values():
                spin_data = redox_data.get("spin_states", {}).get(option_value)
                if spin_data:
                    return spin_data.get("description", "")
            return ""

        return ""

    @classmethod
    def get_parameter_mappings(cls, parameters: Dict[str, Any]) -> Dict[str, str]:
        spin = parameters.get("spin_state", "FR1")
        cys_names = LIGAND_CYS_RESIDUES.get(spin, ["C01", "C02", "C03", "C04"])
        return {
            "cluster_residue_name": spin,
            "cys_1_residue_name": cys_names[0],
            "cys_2_residue_name": cys_names[1],
            "cys_3_residue_name": cys_names[2],
            "cys_4_residue_name": cys_names[3],
        }

    @classmethod
    def get_transformation_sequence(
        cls, components: Dict[str, Any], parameters: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        mappings = cls.get_parameter_mappings(parameters)
        source_resname = components.get("cluster_source_resname")

        transformations: List[Dict[str, Any]] = []

        # 1. Rename the cluster residue to the chosen FOn/FRn/FDn code.
        cluster_selector: Dict[str, Any] = {
            "chain_id": components["center_chain"],
            "residue_id": components["center_id"],
        }
        if source_resname:
            cluster_selector["residue_name"] = source_resname
        transformations.append({
            "id": "rename_cluster_residue",
            "description": (
                f"Rename {source_resname or 'cluster'} residue to "
                f"{mappings['cluster_residue_name']}"
            ),
            "selector": cluster_selector,
            "action": {"change_residue_name": mappings["cluster_residue_name"]},
        })

        # 2-5. Rename each ligating Cys to its position-specific Cn code.
        for n in (1, 2, 3, 4):
            cys_id_key = f"cys_{n}_id"
            cys_chain_key = f"cys_{n}_chain"
            new_name = mappings[f"cys_{n}_residue_name"]
            transformations.append({
                "id": f"rename_ligating_cys_{n}",
                "description": (
                    f"Rename CYS at {components[cys_chain_key]}:"
                    f"{components[cys_id_key]} (ligating FE{n}) to {new_name}"
                ),
                "selector": {
                    "chain_id": components[cys_chain_key],
                    "residue_name": "CYS",
                    "residue_id": components[cys_id_key],
                },
                "action": {"change_residue_name": new_name},
            })

        return transformations

    @classmethod
    def get_site_requirements(cls) -> Dict[str, Any]:
        return {
            "centers": {
                "required_count": 0,
                "elements": ["Fe"],
            },
            "atoms": {
                "required_residues": {
                    "CYS": {"min_count": 4, "max_count": 4},
                },
            },
            "bonds": {
                "required_bond_groups": [
                    {
                        "description": "4 coordinate Fe-SG bonds (ligating Cys)",
                        "min_count": 4,
                        "bond_types": {
                            "coordinate": [
                                (("ANY", "FE1"), ("CYS", "SG")),
                                (("ANY", "FE2"), ("CYS", "SG")),
                                (("ANY", "FE3"), ("CYS", "SG")),
                                (("ANY", "FE4"), ("CYS", "SG")),
                            ],
                        },
                    },
                ],
            },
        }
