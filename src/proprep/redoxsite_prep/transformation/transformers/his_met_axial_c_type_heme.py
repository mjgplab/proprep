#!/usr/bin/env python3
"""
His/Met-Axial C-Type Heme Transformer for RedoxSite System

Handles the His/Met-axial c-type heme of cytochrome c2 (and the wider
cytochrome-c/c2 family that uses a His(NE2)/Met(SD) axial pair). Like the
bis-his c-type heme, the macrocycle is covalently anchored to the protein
through two thioether linkages (the CXXCH motif), and every redox-active
atom is carried in a single heme residue (HMO ferric / HMR ferrous): the
porphyrin plus the side chains of both thioether Cys, the axial His
imidazole, and the axial Met thioether are absorbed into the heme unit.
The four source protein residues are reduced to backbone-only stubs
(CYO ×2, HIO, MEO) and the two propionates are carved out as PRN, all
re-attached with explicit cross-residue bonds.

Distinct from bis-his c-type (this set swaps one axial His for an axial
Met, so it migrates a Met side chain and leaves an MEO stub) and from the
Cys-axial b-type (no thioether links, single proximal Cys).

The Fe-NE2(His) and Fe-SD(Met) axial coordination bonds and the two
thioether S-C bonds are INTERNAL to HMO/HMR (declared in the lib after
migration) — they are deliberately not required or re-emitted here.
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
class HisMetAxialCTypeHemeTransformer(RedoxSiteTransformerBase):
    """Transformer for His/Met-axial c-type hemes (cytochrome c2 family)."""

    TRANSFORMER_NAME = "heme_his_met_axial_c_type"
    DESCRIPTION = (
        "His/Met-axial c-type heme (cytochrome c2; thioether-linked CXXCH "
        "motif, six-coordinate with axial His-NE2 and Met-SD)"
    )
    SUPPORTED_SITE_TYPES = ["heme_his_met_axial_c_type"]
    FORCEFIELD_PATH = "heme/his_met_axial_c_type"

    # ────────────────────────────────────────────────────────────────────
    # Site evaluation — standard requirements-driven validators.
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
        # Exact criteria: one axial His AND one axial Met, both required. Under
        # the previous 80% threshold this transformer matched bis-His sites —
        # the missing Met failed its min_count check but scored a free pass on
        # its max_count check, and the surplus His likewise split 1-1, leaving
        # the total above threshold. See bis_his_c_type_heme for the full note.
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
        """Locate the heme center, its two thioether Cys, the axial His, and
        the axial Met.

        Ligand identification follows the established heme-transformer
        conventions:
          * The two thioether Cys are disambiguated geometrically by which
            heme vinyl-alpha carbon their SG sits closest to — the Cys near
            CAB (ring B) maps to the heme's CBC1/SGC1 slot, the Cys near CAC
            (ring C) maps to CBC2/SGC2. (Fe-S coordination is no longer in
            the required bond set, so we use coordinates, not bonds.)
          * The axial His is identified by the CXXCH sequence rule: it is the
            residue immediately following the ring-C thioether Cys
            (c-ring Cys resid + 1). Falls back to the unique HIS in the site.
          * The axial Met is the unique MET present in the curated site (it
            was added during template refinement precisely because its SD
            coordinates the Fe — same by-presence logic the bis-his
            transformer uses for its axial His residues).
        """
        matched: Dict[str, Any] = {}
        missing: List[str] = []

        # Step 1: heme center (organic/organometallic cofactor, HEM or HEC)
        heme_centers = []
        for c in redox_site.centers:
            if not hasattr(c, "center_type"):
                continue
            ct = c.center_type
            type_value = ct.value if hasattr(ct, "value") else str(ct).lower()
            if type_value in ("organic_cofactor", "organometallic_cofactor"):
                if getattr(c, "resname", None) in ("HEM", "HEC"):
                    heme_centers.append(c)

        if not heme_centers:
            logger.warning(f"No heme center in site {redox_site.site_id}")
            missing.extend(["center_id", "center_chain"])
            return matched, missing

        heme = heme_centers[0]
        matched["center_id"] = heme.resid
        matched["center_chain"] = heme.chain

        # Step 2: heme vinyl-alpha carbons CAB (ring B) and CAC (ring C)
        cab_atom = None
        cac_atom = None
        for atom in redox_site.atoms:
            if atom.resname in ("HEM", "HEC"):
                if atom.atom_name == "CAB":
                    cab_atom = atom
                elif atom.atom_name == "CAC":
                    cac_atom = atom

        # Step 3: thioether Cys, assigned by SG proximity to CAB / CAC
        cys_sg_atoms = [
            atom for atom in redox_site.atoms
            if atom.resname == "CYS" and atom.atom_name == "SG"
        ]

        def _closest_cys(ref_atom):
            if ref_atom is None or not cys_sg_atoms:
                return None
            best = None
            best_d = float("inf")
            for sg in cys_sg_atoms:
                dx = sg.coords[0] - ref_atom.coords[0]
                dy = sg.coords[1] - ref_atom.coords[1]
                dz = sg.coords[2] - ref_atom.coords[2]
                d = (dx * dx + dy * dy + dz * dz) ** 0.5
                if d < best_d:
                    best_d = d
                    best = sg
            return best

        bring = _closest_cys(cab_atom)   # ring B Cys -> CBC1/SGC1
        cring = _closest_cys(cac_atom)   # ring C Cys -> CBC2/SGC2

        if bring is not None and cring is not None and (
            bring.chain == cring.chain and bring.resid == cring.resid
        ):
            logger.warning(
                f"Both thioether Cys resolved to the same residue in site "
                f"{redox_site.site_id}; cannot disambiguate ring B vs ring C."
            )
            bring = cring = None

        if bring is not None:
            matched["bring_cys_id"] = bring.resid
            matched["bring_cys_chain"] = bring.chain
        else:
            missing.extend(["bring_cys_id", "bring_cys_chain"])

        if cring is not None:
            matched["cring_cys_id"] = cring.resid
            matched["cring_cys_chain"] = cring.chain
        else:
            missing.extend(["cring_cys_id", "cring_cys_chain"])

        # Step 4: axial His by CXXCH sequence rule (c-ring Cys + 1),
        #          falling back to the unique HIS in the site.
        his_keys = sorted({
            (atom.chain, atom.resid)
            for atom in redox_site.atoms if atom.resname == "HIS"
        })
        his_choice = None
        if cring is not None:
            candidate = (cring.chain, cring.resid + 1)
            if candidate in his_keys:
                his_choice = candidate
        if his_choice is None and len(his_keys) == 1:
            his_choice = his_keys[0]
            if cring is not None:
                logger.info(
                    "Axial His did not match the CXXCH sequence rule "
                    "(c-ring Cys + 1); using the unique HIS in the site."
                )
        if his_choice is not None:
            matched["his_chain"], matched["his_id"] = his_choice
        else:
            missing.extend(["his_id", "his_chain"])

        # Step 5: axial Met — the unique MET in the curated site.
        met_keys = sorted({
            (atom.chain, atom.resid)
            for atom in redox_site.atoms if atom.resname == "MET"
        })
        if len(met_keys) == 1:
            matched["met_chain"], matched["met_id"] = met_keys[0]
        else:
            logger.warning(
                f"Expected exactly 1 axial MET, found {len(met_keys)} in site "
                f"{redox_site.site_id}"
            )
            missing.extend(["met_id", "met_chain"])

        # Step 6: propionate residue IDs (same convention as the other hemes)
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
        """3 new residue IDs: heme + 2 propionates. The thioether Cys, axial
        His, and axial Met keep their original residue IDs (only their
        resnames change to the CYO/HIO/MEO backbone stubs)."""
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
                    "Electronic spin state (the strong His/Met axial field "
                    "enforces low-spin in both redox states)"
                ),
                "type": "fixed",
                "value": "low_spin",
                "note": (
                    "Ferric is S=1/2 (low-spin doublet), ferrous is S=0 "
                    "(low-spin closed-shell singlet)"
                ),
            },
        }
        # Generic pH-treatment fork (constant_pH PRN vs fixed_pH PRD/PRP) + per-ring
        # protomer choices, from this cofactor's metadata. Mirrors cys_axial.
        defs.update(cls.protonation_parameter_definitions())
        return defs

    @classmethod
    def validate_parameters(cls, parameters: Dict[str, Any]) -> Tuple[bool, str]:
        # Validate against get_valid_options so gated protonation_<role> params
        # (valid only under fixed_pH) behave correctly — same as cys_axial.
        param_defs = cls.get_parameter_definitions()
        for name, defn in param_defs.items():
            others = {k: v for k, v in parameters.items() if k != name}
            valid = cls.get_valid_options(name, others)
            if not valid:
                continue
            if name not in parameters:
                return False, f"Missing required parameter: {name}"
            if defn["type"] == "choice" and parameters[name] not in valid:
                return False, f"Invalid {name}: {parameters[name]}. Options: {valid}"
        return True, "Parameters valid"

    @classmethod
    def get_parameter_mappings(cls, parameters: Dict[str, Any]) -> Dict[str, str]:
        redox = parameters.get("redox_state", "oxidized")
        mappings = {
            "oxidized": {
                "heme_name": "HMO",
                "thioether_cys_name": "CYO",
                "axial_his_name": "HIO",
                "axial_met_name": "MEO",
            },
            "reduced": {
                "heme_name": "HMR",
                "thioether_cys_name": "CYO",
                "axial_his_name": "HIO",
                "axial_met_name": "MEO",
            },
        }
        return mappings.get(redox, mappings["oxidized"])

    # ────────────────────────────────────────────────────────────────────
    # Transformation sequence
    # ────────────────────────────────────────────────────────────────────

    # Migrated side-chain atom maps (PDB name -> heme-unit name), matching
    # the bundled HMO/HMR lib roster exactly. HMR is atom-for-atom identical
    # to HMO, so the same maps serve both redox states.
    _BRING_CYS_MAP = {"CB": "CBC1", "SG": "SGC1", "HB2": "HB21", "HB3": "HB31"}
    _CRING_CYS_MAP = {"CB": "CBC2", "SG": "SGC2", "HB2": "HB22", "HB3": "HB32"}
    # Axial His is the HID tautomer (proton on ND1; NE2 free to coordinate Fe).
    # Only CB and its hydrogens + CG are renamed; the imidazole ring atoms keep
    # their standard names (they already exist under those names in the lib).
    _HIS_SIDECHAIN = [
        "CB", "HB2", "HB3", "CG", "ND1", "HD1",
        "CE1", "HE1", "NE2", "CD2", "HD2",
    ]
    _HIS_MAP = {
        "CB": "CBH", "HB2": "HB2H", "HB3": "HB3H", "CG": "CGH",
        "ND1": "ND1", "HD1": "HD1", "CE1": "CE1", "HE1": "HE1",
        "NE2": "NE2", "CD2": "CD2", "HD2": "HD2",
    }
    _MET_SIDECHAIN = [
        "CB", "HB2", "HB3", "CG", "HG2", "HG3", "SD", "CE", "HE1", "HE2", "HE3",
    ]
    _MET_MAP = {
        "CB": "CBM", "HB2": "HB2M", "HB3": "HB3M",
        "CG": "CGM", "HG2": "HG2M", "HG3": "HG3M",
        "SD": "SD", "CE": "CE",
        "HE1": "HE1M", "HE2": "HE2M", "HE3": "HE3M",
    }
    _BACKBONE = ["N", "H", "CA", "HA", "C", "O"]

    @classmethod
    def get_transformation_sequence(
        cls, components: Dict[str, Any], parameters: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Build the transformation sequence.

        All side-chain migrations land in the heme residue while it is still
        named HEC; the HEC -> HMO/HMR rename is deferred to the very end so
        the migrated atoms don't land in a phantom residue (same discipline
        as the bis-his c-type transformer).

          1.  HEM -> HEC (if needed)
          2.  ring-B Cys side chain -> heme (CBC1/SGC1/HB21/HB31)
          3.  ring-C Cys side chain -> heme (CBC2/SGC2/HB22/HB32)
          4.  ring-B Cys backbone -> CYO
          5.  ring-C Cys backbone -> CYO
          6.  axial His side chain -> heme (CBH/CGH + imidazole)
          7.  axial His backbone -> HIO
          8.  axial Met side chain -> heme (CBM/CGM/SD/CE + H)
          9.  axial Met backbone -> MEO
          10. propionate A -> PRN
          11. propionate D -> PRN
          12. HEC -> HMO/HMR
        """
        m = cls.get_parameter_mappings(parameters)
        # Propionate names are pH-treatment dependent — resolved from metadata by
        # role (PRN constant_pH, or PRD/PRP per ring fixed_pH). Falls back to PRN.
        resolved = cls.resolve_output_residue_names(parameters)
        prop_a_name = resolved.get("propionate_a", "PRN")
        prop_d_name = resolved.get("propionate_d", "PRN")
        seq: List[Dict[str, Any]] = []

        # 1. HEM -> HEC (no-op if already HEC)
        seq.append({
            "id": "rename_hem_to_hec",
            "description": "Rename HEM to HEC (if needed)",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEM",
                "residue_id": components["center_id"],
            },
            "action": {"change_residue_name": "HEC"},
        })

        # Helper to build a side-chain migration step
        def migrate_sidechain(step_id, desc, chain, resid, src_resname,
                              atom_names, rename_map):
            seq.append({
                "id": step_id,
                "description": desc,
                "selector": {
                    "chain_id": chain,
                    "residue_name": src_resname,
                    "residue_id": resid,
                    "atom_names": atom_names,
                },
                "action": {
                    "change_residue_name": "HEC",
                    "change_residue_id": components["center_id"],
                    "change_chain_id": components["center_chain"],
                    "change_insertion_code": "",
                    "rename_atoms": rename_map,
                    "convert_to_hetatm": True,
                },
            })

        def convert_backbone(step_id, desc, chain, resid, src_resname, new_name):
            seq.append({
                "id": step_id,
                "description": desc,
                "selector": {
                    "chain_id": chain,
                    "residue_name": src_resname,
                    "residue_id": resid,
                    "atom_names": cls._BACKBONE,
                },
                "action": {"change_residue_name": new_name},
            })

        # 2-3. thioether Cys side chains -> heme
        migrate_sidechain(
            "transform_bring_cys_sidechain",
            "Move ring-B thioether Cys side chain (CB/SG/HB2/HB3) to heme as "
            "CBC1/SGC1/HB21/HB31",
            components["bring_cys_chain"], components["bring_cys_id"], "CYS",
            list(cls._BRING_CYS_MAP), cls._BRING_CYS_MAP,
        )
        migrate_sidechain(
            "transform_cring_cys_sidechain",
            "Move ring-C thioether Cys side chain (CB/SG/HB2/HB3) to heme as "
            "CBC2/SGC2/HB22/HB32",
            components["cring_cys_chain"], components["cring_cys_id"], "CYS",
            list(cls._CRING_CYS_MAP), cls._CRING_CYS_MAP,
        )

        # 4-5. thioether Cys backbones -> CYO stubs
        convert_backbone(
            "convert_bring_cys_to_cyo",
            "Reduce ring-B thioether Cys to CYO backbone stub",
            components["bring_cys_chain"], components["bring_cys_id"], "CYS",
            m["thioether_cys_name"],
        )
        convert_backbone(
            "convert_cring_cys_to_cyo",
            "Reduce ring-C thioether Cys to CYO backbone stub",
            components["cring_cys_chain"], components["cring_cys_id"], "CYS",
            m["thioether_cys_name"],
        )

        # 6-7. axial His side chain -> heme, backbone -> HIO
        migrate_sidechain(
            "transform_axial_his_sidechain",
            "Move axial His side chain (CB + imidazole) to heme as "
            "CBH/CGH + ND1/HD1/CE1/HE1/NE2/CD2/HD2",
            components["his_chain"], components["his_id"], "HIS",
            cls._HIS_SIDECHAIN, cls._HIS_MAP,
        )
        convert_backbone(
            "convert_axial_his_to_hio",
            "Reduce axial His to HIO backbone stub",
            components["his_chain"], components["his_id"], "HIS",
            m["axial_his_name"],
        )

        # 8-9. axial Met side chain -> heme, backbone -> MEO
        migrate_sidechain(
            "transform_axial_met_sidechain",
            "Move axial Met side chain (CB/CG/SD/CE + H) to heme as "
            "CBM/CGM/SD/CE + HB2M/HB3M/HG2M/HG3M/HE1M/HE2M/HE3M",
            components["met_chain"], components["met_id"], "MET",
            cls._MET_SIDECHAIN, cls._MET_MAP,
        )
        convert_backbone(
            "convert_axial_met_to_meo",
            "Reduce axial Met to MEO backbone stub",
            components["met_chain"], components["met_id"], "MET",
            m["axial_met_name"],
        )

        # 10-11. propionates -> PRN (constant_pH) or PRD/PRP per ring (fixed_pH)
        seq.append({
            "id": "extract_propionate_a",
            "description": f"Create {prop_a_name} residue from propionate A group",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEC",
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
        seq.append({
            "id": "extract_propionate_d",
            "description": f"Create {prop_d_name} residue from propionate D group",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEC",
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

        # 12. final HEC -> HMO/HMR
        seq.append({
            "id": "apply_redox_specific_heme_name",
            "description": f"Rename HEC to {m['heme_name']}",
            "selector": {
                "chain_id": components["center_chain"],
                "residue_name": "HEC",
                "residue_id": components["center_id"],
            },
            "action": {"change_residue_name": m["heme_name"]},
        })

        return seq

    # ────────────────────────────────────────────────────────────────────
    # Site-requirements declaration
    # ────────────────────────────────────────────────────────────────────

    @classmethod
    def get_site_requirements(cls) -> Dict[str, Any]:
        """His/Met-axial c-type heme:
          * 1 heme center (HEM or HEC)
          * 2 thioether CYS, 1 axial HIS, 1 axial MET
          * The bonds that become INTER-residue after the transformer migrates
            the side chains into the heme: each ligand's CA-CB (-> stub-heme
            splice) plus the two propionate carve-outs (C2A-CAA, C3D-CAD).

        The Fe-NE2(His), Fe-SD(Met) coordination bonds and the two thioether
        S-C bonds are deliberately NOT required: the migration pulls NE2/SD/SG
        into the heme residue, where the lib template already declares them.
        Requiring them would make the tleap generator duplicate bonds the
        template owns, triggering tleap's "cannot add bond" fatal error.
        """
        return {
            "centers": {
                "required_count": 1,
                "center_types": [CenterType.ORGANOMETALLIC_COFACTOR],
                "residue_names": ["HEM", "HEC"],
            },
            "atoms": {
                "required_residues": {
                    "HEM": {"min_count": 0, "max_count": 1},
                    "HEC": {"min_count": 0, "max_count": 1},
                    "CYS": {"min_count": 2, "max_count": 2},
                    "HIS": {"min_count": 1, "max_count": 1},
                    "MET": {"min_count": 1, "max_count": 1},
                },
                "alternative_groups": [
                    ["HEM", "HEC"],  # either HEM or HEC, not both
                ],
            },
            "bonds": {
                "required_bond_groups": [
                    {
                        "description": (
                            "HEM-variant bonds (inter-residue post-transformation): "
                            "2 thioether-Cys + axial-His + axial-Met side-chain "
                            "splices + 2 propionate carve-outs"
                        ),
                        "min_count": 6,
                        "bond_types": {
                            "covalent": [
                                (("HEM", "C2A"), ("HEM", "CAA")),   # -> heme-PRN_A
                                (("HEM", "C3D"), ("HEM", "CAD")),   # -> heme-PRN_D
                                (("CYS", "CA"), ("CYS", "CB")),     # -> CYO-heme (ring B)
                                (("CYS", "CA"), ("CYS", "CB")),     # -> CYO-heme (ring C)
                                (("HIS", "CA"), ("HIS", "CB")),     # -> HIO-heme
                                (("MET", "CA"), ("MET", "CB")),     # -> MEO-heme
                            ],
                        },
                    },
                    {
                        "description": (
                            "HEC-variant bonds (inter-residue post-transformation)"
                        ),
                        "min_count": 6,
                        "bond_types": {
                            "covalent": [
                                (("HEC", "C2A"), ("HEC", "CAA")),
                                (("HEC", "C3D"), ("HEC", "CAD")),
                                (("CYS", "CA"), ("CYS", "CB")),
                                (("CYS", "CA"), ("CYS", "CB")),
                                (("HIS", "CA"), ("HIS", "CB")),
                                (("MET", "CA"), ("MET", "CB")),
                            ],
                        },
                    },
                ],
                "require_one_group": True,
            },
        }
