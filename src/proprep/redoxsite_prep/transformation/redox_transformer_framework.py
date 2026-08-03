#!/usr/bin/env python3
"""
RedoxSite Transformation Framework

Enhanced transformation system using RedoxSite objects with coordinate-based tracking.
Handles all types of redox sites: metal centers, organic cofactors, and redox-active amino acids.

Author: Claude Code Implementation
"""

import copy
import logging
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from rich.console import Console
from rich.panel import Panel
from rich.prompt import Confirm, Prompt
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn

from proprep.utils.prompts import prompt_with_context, confirm_with_context

logger = logging.getLogger(__name__)

# ===== EVALUATION RESULT CLASSES =====

@dataclass
class TransformerEvaluationDetail:
    """Individual requirement check result"""
    description: str
    passed: bool
    value_found: Any = None
    value_expected: Any = None
    error_message: Optional[str] = None

@dataclass 
class TransformerEvaluation:
    """Complete evaluation result from a transformer"""
    is_valid: bool
    confidence: float  # 0.0 to 1.0
    description: str
    requirements_met: int
    total_requirements: int
    details: List[TransformerEvaluationDetail]
    error_message: Optional[str] = None

# ===== BASE TRANSFORMER CLASS =====

class RedoxSiteTransformerBase:
    """
    Base class that all site-specific transformers must implement.

    A transformer is a reproducible recipe that says: "whenever you see this type of site
    (defined by requirements), identify the components (via matching logic), and apply these
    exact transformations" - making the same forcefield-compatible modifications every time,
    regardless of which protein structure contains the site.
    """

    # Each transformer must define these class attributes
    TRANSFORMER_NAME: str = "base"
    DESCRIPTION: str = "Base transformer"
    SUPPORTED_SITE_TYPES: List[str] = []
    FORCEFIELD_PATH: Optional[str] = None
    
    @classmethod
    def evaluate_redox_site(cls, redox_site) -> TransformerEvaluation:
        """
        Comprehensive evaluation of RedoxSite compatibility
        Must be implemented by each transformer
        
        Args:
            redox_site: RedoxSite object to evaluate
            
        Returns:
            TransformerEvaluation with compatibility assessment
        """
        raise NotImplementedError("Each transformer must implement evaluate_redox_site")
    
    @classmethod
    def get_site_requirements(cls) -> Dict[str, Any]:
        """
        Define what the site must contain
        
        Returns:
            Dictionary defining site requirements (centers, residues, bonds)
        """
        raise NotImplementedError("Each transformer must define site requirements")
    
    @classmethod
    def match_components(cls, redox_site) -> Tuple[Dict[str, Any], List[str]]:
        """
        Site-specific logic to identify component roles from RedoxSite structure

        Args:
            redox_site: RedoxSite object

        Returns:
            Tuple of (matched_components, missing_components)

        Ambiguous role assignments
        --------------------------
        When a transformer can identify a set of candidate residues but cannot
        choose which candidate fills which role from geometry alone (e.g.
        proximal vs distal His on a symmetric bis-His heme), it may declare an
        ``"_ambiguous"`` key in ``matched_components``. The transformation
        manager will prompt the user to resolve each ambiguity and then strip
        the ``"_ambiguous"`` key before downstream code sees the dict.

        ``_ambiguous`` is a list of blocks; each block has shape::

            {
                "label": str,                # short title shown to user
                "description": str,          # explanation (one sentence)
                "candidates": [
                    {"chain": "A", "resid": 23, "display": "HIS A23"},
                    ...
                ],
                "roles": {<role_name>: <count>, ...},
            }

        For each role the user picks ``count`` candidates from the remaining
        pool. Resolved values are written back as ``<role>_id`` / ``<role>_chain``
        when ``count == 1`` and ``<role>_id_<n>`` / ``<role>_chain_<n>``
        (1-indexed) when ``count > 1``. The total slots requested
        (``sum(roles.values())``) must equal ``len(candidates)``; otherwise the
        block is skipped with a warning.
        """
        raise NotImplementedError("Each transformer must implement component matching")
    
    @classmethod
    def get_transformation_sequence(cls, components: Dict[str, Any], parameters: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Generate ordered transformation sequence for PDB modification
        
        Args:
            components: Matched site components (from match_components)
            parameters: User-specified parameters (redox state, etc.)
            
        Returns:
            List of transformation dictionaries in execution order
        """
        raise NotImplementedError("Each transformer must define transformation sequence")
    
    @classmethod
    def get_required_residue_count(cls) -> int:
        """
        How many residue IDs this transformer needs (including original)
        
        Returns:
            Total residue IDs required for this transformation
        """
        raise NotImplementedError("Each transformer must specify residue space requirements")
    
    @classmethod 
    def get_residue_space_plan(cls, components: Dict[str, Any]) -> Dict[str, int]:
        """
        Detailed breakdown of residue space allocation
        
        Args:
            components: Matched components
            
        Returns:
            Dict mapping component roles to relative residue offsets
        """
        raise NotImplementedError("Each transformer must define residue space plan")
    
    @classmethod
    def get_parameter_definitions(cls) -> Dict[str, Any]:
        """
        Define required parameters and their validation
        
        Returns:
            Dictionary defining parameter structure and validation rules
        """
        raise NotImplementedError("Each transformer must define parameters")
    
    @classmethod
    def validate_parameters(cls, parameters: Dict[str, Any]) -> Tuple[bool, str]:
        """
        Validate parameter values

        Args:
            parameters: Parameter values to validate

        Returns:
            Tuple of (is_valid, error_message)
        """
        raise NotImplementedError("Each transformer must implement parameter validation")

    @classmethod
    def get_valid_options(cls, param_name: str,
                          current_parameters: Dict[str, Any]) -> List[Any]:
        """
        Return the valid option values for ``param_name`` given the parameters
        already chosen for this site.

        Default implementation: returns the static ``options`` list from
        :py:meth:`get_parameter_definitions` for choice-typed parameters,
        ``[value]`` for fixed-typed parameters, and ``[]`` for unknown params
        or types. Override in transformers whose parameter options depend on
        other parameter values (e.g. spin_state options gated by redox_state).

        Args:
            param_name: Name of the parameter being configured
            current_parameters: Parameter values already assigned for this site
                (does not include ``param_name``)

        Returns:
            List of valid option values for this parameter on this site.
        """
        # Generic pH-treatment / protonation handling (active only for cofactors
        # whose metadata declares them; see protonation_parameter_definitions).
        if param_name == cls.PH_TREATMENT_PARAM:
            redox = current_parameters.get("redox_state")
            spin = current_parameters.get("spin_state")
            if redox and spin:
                treatments = cls.available_ph_treatments(redox, spin)
                if treatments:
                    return treatments
        elif param_name.startswith(cls.PROTONATION_PARAM_PREFIX):
            # A per-site protomer choice is valid ONLY under fixed_pH; otherwise
            # return [] so the manager (whose static options are also empty for
            # these params) skips it. Under fixed_pH, options are that site's
            # variant keys from the selected set's protonation_model.
            if current_parameters.get(cls.PH_TREATMENT_PARAM) != "fixed_pH":
                return []
            role = param_name[len(cls.PROTONATION_PARAM_PREFIX):]
            for site in cls._protonation_sites(current_parameters):
                if site.get("role") == role and "variants" in site:
                    return list(site["variants"].keys())
            return []

        defs = cls.get_parameter_definitions()
        param_def = defs.get(param_name)
        if not param_def:
            return []
        if param_def.get("type") == "choice":
            return list(param_def.get("options", []))
        if param_def.get("type") == "fixed":
            return [param_def.get("value")]
        return []

    @classmethod
    def get_option_description(cls, param_name: str, option_value: Any) -> str:
        """
        Return a human-readable description for one option of a parameter.

        The transformation manager surfaces this in the parameter-configuration
        prompt so users can see what each option means (e.g. "BS partition
        {Fe1,Fe2}|{Fe3,Fe4}, spin variant a") before picking. Default returns
        ``""`` (no description shown). Override to provide context.

        Args:
            param_name: Name of the parameter
            option_value: One of the option values for that parameter

        Returns:
            Short description string (single line preferred).
        """
        # Generic descriptions for the metadata-driven pH-treatment fork and the
        # per-ring protomer choices, so every fixed-pH-capable cofactor explains
        # this (consequential) choice in the parameter prompt without each
        # transformer re-stating it. Subclasses that override this for their own
        # params should call super() to keep these.
        if param_name == cls.PH_TREATMENT_PARAM:
            return {
                "fixed_pH": (
                    "Static propionates (PRP protonated / PRD deprotonated) bundled in the "
                    "cofactor library; keeps the modern ff14SB/ff19SB protein backbone. Use "
                    "for ordinary fixed-protonation MD - you pick each ring's protonation."
                ),
                "constant_pH": (
                    "Titratable PRN propionates from AMBER's constant-(pH,E) libraries "
                    "(requires leaprc.constph + leaprc.conste, ff10 backbone). Use for "
                    "constant-pH / constant-(pH,E) MD - protonation is sampled dynamically, "
                    "not chosen here."
                ),
            }.get(option_value, "")
        if param_name.startswith(cls.PROTONATION_PARAM_PREFIX):
            return {
                "deprotonated": "Deprotonated carboxylate (COO-, net -1) - standard above the propionate pKa (~4.8).",
                "protonated": "Protonated carboxylic acid (COOH, net 0, syn conformer).",
            }.get(option_value, "")
        return ""
    
    @classmethod
    def get_parameter_mappings(cls, parameters: Dict[str, Any]) -> Dict[str, str]:
        """
        Map user parameters to forcefield-specific names
        
        Args:
            parameters: User parameter values
            
        Returns:
            Dictionary mapping component types to residue names
        """
        raise NotImplementedError("Each transformer must implement parameter mappings")

    # ────────────────────────────────────────────────────────────────────
    # Protonation / pH-treatment infrastructure (generic, metadata-driven)
    #
    # These helpers let any cofactor whose force-field sets declare a
    # ph_treatment / protonation_model (see forcefield_params.loader) expose a
    # "constant_pH vs fixed_pH" fork plus per-site protomer choices, and resolve
    # every output residue name from metadata by role — without the transformer
    # hardcoding residue codes. Cofactors with no such metadata are unaffected:
    # every method degrades to "no treatments / no sites / center+ligand only".
    #
    # PARAM NAMING CONVENTION:
    #   "ph_treatment"            → constant_pH | fixed_pH (the fork)
    #   "protonation_<role>"      → one choice per fixed-pH protonation site
    # ────────────────────────────────────────────────────────────────────

    PROTONATION_PARAM_PREFIX = "protonation_"
    PH_TREATMENT_PARAM = "ph_treatment"

    @classmethod
    def _forcefield_sets(cls, redox_state: str, spin_state: str) -> List[Dict[str, Any]]:
        """All discovered forcefield sets for a state ([] if none / on error)."""
        if not cls.FORCEFIELD_PATH:
            return []
        try:
            from proprep.forcefield_params import discover_forcefield_files
            return discover_forcefield_files(cls.FORCEFIELD_PATH, redox_state, spin_state)
        except Exception as e:
            logger.debug("discover_forcefield_files failed for %s/%s/%s: %s",
                         cls.FORCEFIELD_PATH, redox_state, spin_state, e)
            return []

    @classmethod
    def _all_redox_spin_pairs(cls) -> List[Tuple[str, str]]:
        """Every (redox_state, spin_state) pair declared in this cofactor's metadata."""
        if not cls.FORCEFIELD_PATH:
            return []
        try:
            from proprep.forcefield_params import load_forcefield_metadata
            meta = load_forcefield_metadata(cls.FORCEFIELD_PATH)
        except Exception:
            return []
        pairs = []
        for rstate, rdata in (meta.get("redox_states", {}) or {}).items():
            for sstate in (rdata.get("spin_states", {}) or {}).keys():
                pairs.append((rstate, sstate))
        return pairs

    @classmethod
    def available_ph_treatments(cls, redox_state: str, spin_state: str) -> List[str]:
        """Distinct, order-preserved ph_treatment values among a state's sets
        (excluding sets that don't declare one)."""
        seen, treatments = set(), []
        for s in cls._forcefield_sets(redox_state, spin_state):
            t = s.get("ph_treatment")
            if t and t not in seen:
                seen.add(t)
                treatments.append(t)
        return treatments

    @classmethod
    def select_forcefield_set_name(cls, parameters: Dict[str, Any]) -> Optional[str]:
        """Pick a representative set whose names this state+treatment imply.

        Residue names (center/ligand/protomer) don't depend on the charge-method
        axis (RESP vs CM5) — only on redox/spin + ph_treatment — so for NAME
        resolution any set matching the chosen treatment is equivalent. An
        explicit ``forcefield_set`` parameter wins if present; otherwise prefer
        the is_default set, else the first. The specific lib (charge method) is
        chosen later in the Topology Generator.
        """
        redox = parameters.get("redox_state")
        spin = parameters.get("spin_state")
        if not redox or not spin:
            return None
        sets = cls._forcefield_sets(redox, spin)
        if not sets:
            return None

        explicit = parameters.get("forcefield_set")
        if explicit:
            for s in sets:
                if s["name"] == explicit:
                    return explicit

        treatment = parameters.get(cls.PH_TREATMENT_PARAM)
        if treatment:
            sets = [s for s in sets if s.get("ph_treatment") == treatment] or sets

        for s in sets:
            if s.get("is_default"):
                return s["name"]
        return sets[0]["name"]

    @classmethod
    def _protonation_sites(cls, parameters: Dict[str, Any]) -> List[Dict[str, Any]]:
        """protonation_model sites for the set the parameters select ([] if none)."""
        set_name = cls.select_forcefield_set_name(parameters)
        if not set_name:
            return []
        try:
            from proprep.forcefield_params import get_protonation_model
            model = get_protonation_model(
                cls.FORCEFIELD_PATH, parameters["redox_state"],
                parameters["spin_state"], set_name)
        except Exception:
            return []
        return list((model or {}).get("sites", []) or [])

    @classmethod
    def protonation_parameter_definitions(cls) -> Dict[str, Any]:
        """Build the generic ph_treatment + per-site protonation parameter defs
        for a transformer to merge into ``get_parameter_definitions``.

        - ``ph_treatment`` is emitted only when ≥2 distinct treatments exist
          across the cofactor's states (otherwise there's no fork — every
          current cofactor returns {}).
        - one ``protonation_<role>`` choice per fixed-pH site role found in any
          fixed-pH set. Static ``options`` is intentionally EMPTY so the
          parameter is *gated* by get_valid_options: it shows only when
          ph_treatment == fixed_pH (the manager skips a param whose static
          options are empty and whose get_valid_options returns []).
        """
        defs: Dict[str, Any] = {}

        # ph_treatment fork — union across all states
        all_treatments, seen = [], set()
        roles_defaults: Dict[str, Tuple[str, str]] = {}  # role -> (label, default)
        for redox, spin in cls._all_redox_spin_pairs():
            for s in cls._forcefield_sets(redox, spin):
                t = s.get("ph_treatment")
                if t and t not in seen:
                    seen.add(t)
                    all_treatments.append(t)
                model = s.get("protonation_model") or {}
                if s.get("ph_treatment") == "fixed_pH":
                    for site in model.get("sites", []) or []:
                        role = site.get("role")
                        if role and role not in roles_defaults and "variants" in site:
                            roles_defaults[role] = (
                                site.get("label", role),
                                site.get("default") or next(iter(site["variants"]), None),
                            )

        if len(all_treatments) < 2:
            return defs  # no fork → no generic params

        default_treatment = "constant_pH" if "constant_pH" in all_treatments else all_treatments[0]
        defs[cls.PH_TREATMENT_PARAM] = {
            "description": "Propionate/titratable-site pH treatment",
            "type": "choice",
            "options": list(all_treatments),
            "default": default_treatment,
        }

        for role, (label, default) in roles_defaults.items():
            defs[f"{cls.PROTONATION_PARAM_PREFIX}{role}"] = {
                "description": f"Protonation state of {label}",
                "type": "choice",
                "options": [],  # empty ⇒ gated by get_valid_options (fixed_pH only)
                "default": default,
            }
        return defs

    @classmethod
    def resolve_output_residue_names(cls, parameters: Dict[str, Any]) -> Dict[str, str]:
        """Resolve {role: residue_name} for all output residues from metadata.

        Single source of truth used by ``get_transformation_sequence`` so the
        transformer hardcodes no residue codes. Reads the chosen set's
        protonation_model + the state's center/ligand names. ``protonation_<role>``
        params become the per-site protomer choices.
        """
        set_name = cls.select_forcefield_set_name(parameters)
        if not set_name:
            return {}
        protonation_choices = {
            k[len(cls.PROTONATION_PARAM_PREFIX):]: v
            for k, v in parameters.items()
            if k.startswith(cls.PROTONATION_PARAM_PREFIX) and v is not None
        }
        from proprep.forcefield_params import resolve_residue_names
        return resolve_residue_names(
            cls.FORCEFIELD_PATH, parameters["redox_state"],
            parameters["spin_state"], set_name, protonation_choices)

    @classmethod
    def update_components_with_id_mapping(cls, components: Dict[str, Any],
                                        id_mapping: Dict[Tuple[str, int], int]) -> Dict[str, Any]:
        """
        Update component IDs based on ID mapping results
        
        Args:
            components: Original matched components
            id_mapping: Mapping from (chain, original_id) to new_id
            
        Returns:
            Updated components with new residue IDs
        """
        # Default implementation - transformers can override for custom logic
        updated_components = components.copy()
        
        # Update main site ID if mapped (works for any center type)
        main_id_keys = ["center_id", "site_id", "main_residue_id"]
        main_chain_keys = ["center_chain", "site_chain", "main_residue_chain"]
        
        for id_key, chain_key in zip(main_id_keys, main_chain_keys):
            if id_key in components and chain_key in components:
                original_key = (components[chain_key], components[id_key])
                if original_key in id_mapping:
                    new_id = id_mapping[original_key]
                    updated_components[id_key] = new_id
                    
                    # Update dependent IDs based on space plan
                    space_plan = cls.get_residue_space_plan(components)
                    for role, offset in space_plan.items():
                        dependent_key = f"{role}_id"
                        if dependent_key in updated_components and offset > 0:
                            updated_components[dependent_key] = new_id + offset
                break
        
        return updated_components
    
    @classmethod
    def handle_errors(cls, error_type: str, context: Dict[str, Any]) -> Tuple[bool, str]:
        """
        Site-specific error handling and fallback strategies
        
        Args:
            error_type: Type of error encountered
            context: Error context information
            
        Returns:
            Tuple of (can_continue, error_message)
        """
        return False, f"Unhandled error: {error_type}"
    
    @classmethod
    def validate_centers_from_requirements(cls, redox_site) -> Tuple[int, int, List[TransformerEvaluationDetail]]:
        """
        Generic center validation using get_site_requirements()
        
        Args:
            redox_site: RedoxSite object to validate
            
        Returns:
            Tuple of (requirements_met, total_requirements, validation_details)
        """
        details = []
        requirements_met = 0
        total_requirements = 0
        
        try:
            requirements = cls.get_site_requirements()
            centers_section = requirements.get("centers", {})
            
            # Validate required count
            required_count = centers_section.get("required_count", 0)
            if required_count > 0:
                total_requirements += 1
                centers_found = len(redox_site.centers)
                passed = centers_found == required_count
                
                details.append(TransformerEvaluationDetail(
                    description=f"Required center count",
                    passed=passed,
                    value_found=centers_found,
                    value_expected=required_count
                ))
                
                if passed:
                    requirements_met += 1
            
            # Validate center types
            required_center_types = centers_section.get("center_types", [])
            if required_center_types:
                total_requirements += 1
                matching_types = [c for c in redox_site.centers if c.center_type in required_center_types]
                passed = len(matching_types) >= required_count
                
                details.append(TransformerEvaluationDetail(
                    description=f"Required center types: {required_center_types}",
                    passed=passed,
                    value_found=len(matching_types),
                    value_expected=required_count
                ))
                
                if passed:
                    requirements_met += 1
            
            # Validate elements
            required_elements = centers_section.get("elements", [])
            if required_elements:
                total_requirements += 1
                matching_elements = [c for c in redox_site.centers if c.element in required_elements]
                passed = len(matching_elements) >= required_count

                details.append(TransformerEvaluationDetail(
                    description=f"Required elements: {required_elements}",
                    passed=passed,
                    value_found=len(matching_elements),
                    value_expected=required_count
                ))

                if passed:
                    requirements_met += 1

            # Validate residue names
            required_residue_names = centers_section.get("residue_names", [])
            if required_residue_names:
                total_requirements += 1
                matching_residues = [c for c in redox_site.centers if c.resname in required_residue_names]
                passed = len(matching_residues) >= required_count

                details.append(TransformerEvaluationDetail(
                    description=f"Required residue names: {required_residue_names}",
                    passed=passed,
                    value_found=len(matching_residues),
                    value_expected=required_count
                ))

                if passed:
                    requirements_met += 1

        except Exception as e:
            details.append(TransformerEvaluationDetail(
                description="Center validation error",
                passed=False,
                error_message=str(e)
            ))
        
        return requirements_met, total_requirements, details
    
    @classmethod
    def validate_atoms_from_requirements(cls, redox_site) -> Tuple[int, int, List[TransformerEvaluationDetail]]:
        """
        Generic atom/residue validation using get_site_requirements()
        
        Args:
            redox_site: RedoxSite object to validate
            
        Returns:
            Tuple of (requirements_met, total_requirements, validation_details)
        """
        details = []
        requirements_met = 0
        total_requirements = 0
        
        try:
            requirements = cls.get_site_requirements()
            atoms_section = requirements.get("atoms", {})
            required_residues = atoms_section.get("required_residues", {})
            alternative_groups = atoms_section.get("alternative_groups", [])
            
            # Validate individual residue requirements
            for resname, req in required_residues.items():
                min_count = req.get("min_count", 0)
                max_count = req.get("max_count", float('inf'))
                
                # Count residues of this type
                residues = [atom for atom in redox_site.atoms if atom.resname == resname]
                unique_residues = set((atom.chain, atom.resid) for atom in residues)
                found_count = len(unique_residues)
                
                # Check min count requirement
                if min_count > 0:
                    total_requirements += 1
                    passed_min = found_count >= min_count
                    
                    details.append(TransformerEvaluationDetail(
                        description=f"{resname} residues (min: {min_count})",
                        passed=passed_min,
                        value_found=found_count,
                        value_expected=min_count
                    ))
                    
                    if passed_min:
                        requirements_met += 1
                
                # Check max count requirement  
                if max_count < float('inf'):
                    total_requirements += 1
                    passed_max = found_count <= max_count
                    
                    details.append(TransformerEvaluationDetail(
                        description=f"{resname} residues (max: {max_count})",
                        passed=passed_max,
                        value_found=found_count,
                        value_expected=max_count
                    ))
                    
                    if passed_max:
                        requirements_met += 1
            
            # Validate alternative groups (e.g., either HEM or HEC, not both)
            for group in alternative_groups:
                total_requirements += 1
                group_counts = []
                for resname in group:
                    residues = [atom for atom in redox_site.atoms if atom.resname == resname]
                    unique_residues = set((atom.chain, atom.resid) for atom in residues)
                    group_counts.append(len(unique_residues))
                
                # Exactly one residue type from the group should be present
                present_types = sum(1 for count in group_counts if count > 0)
                passed = present_types == 1
                
                details.append(TransformerEvaluationDetail(
                    description=f"Alternative group {group} (exactly one type)",
                    passed=passed,
                    value_found=present_types,
                    value_expected=1
                ))
                
                if passed:
                    requirements_met += 1
                    
        except Exception as e:
            details.append(TransformerEvaluationDetail(
                description="Atom validation error",
                passed=False,
                error_message=str(e)
            ))
        
        return requirements_met, total_requirements, details
    
    @classmethod
    def validate_bonds_from_requirements(cls, redox_site) -> Tuple[int, int, List[TransformerEvaluationDetail]]:
        """
        Generic bond validation using get_site_requirements()
        Supports both old format (required_bond_types) and new format (required_bond_groups)
        
        Args:
            redox_site: RedoxSite object to validate
            
        Returns:
            Tuple of (bonds_found, bonds_required, validation_details)
        """
        details = []
        total_bonds_found = 0
        total_bonds_required = 0
        
        try:
            requirements = cls.get_site_requirements()
            bonds_section = requirements.get("bonds", {})
            
            # Check if using new bond groups format
            if "required_bond_groups" in bonds_section:
                return cls._validate_bond_groups(redox_site, bonds_section)
            
            # Fall back to old format for backward compatibility
            required_bond_types = bonds_section.get("required_bond_types", {})
            
            for bond_type, bond_info in required_bond_types.items():
                min_count = bond_info.get("min_count", 0)
                atom_pairs = bond_info.get("atom_pairs", [])
                description = bond_info.get("description", f"{bond_type} bonds")
                
                # Count matching bonds in the RedoxSite. Credit is clamped to
                # the requirement for the same reason as in _validate_bond_groups:
                # surplus bonds must not offset a failed requirement elsewhere
                # in the caller's met/total tally.
                bonds_found = cls._count_matching_bonds(redox_site, bond_type, atom_pairs)
                total_bonds_found += min(bonds_found, min_count)
                total_bonds_required += min_count
                
                # Create evaluation detail
                passed = bonds_found >= min_count
                details.append(TransformerEvaluationDetail(
                    description=description,
                    passed=passed,
                    value_found=bonds_found,
                    value_expected=min_count
                ))
                
        except Exception as e:
            details.append(TransformerEvaluationDetail(
                description="Bond validation error",
                passed=False,
                error_message=str(e)
            ))
        
        return total_bonds_found, total_bonds_required, details
    
    @classmethod
    def _validate_bond_groups(cls, redox_site, bonds_section) -> Tuple[int, int, List[TransformerEvaluationDetail]]:
        """
        Validate bond requirements using the new bond groups format
        
        Args:
            redox_site: RedoxSite object to validate
            bonds_section: Bonds section from get_site_requirements()
            
        Returns:
            Tuple of (bonds_found, bonds_required, validation_details)
        """
        details = []
        bond_groups = bonds_section.get("required_bond_groups", [])
        require_one_group = bonds_section.get("require_one_group", False)
        
        group_results = []
        
        # Evaluate each bond group
        for group in bond_groups:
            group_description = group.get("description", "Bond group")
            group_min_count = group.get("min_count", 0)
            bond_types = group.get("bond_types", {})
            
            group_bonds_found = 0
            group_details = []
            
            # Count bonds for each type within the group
            for bond_type, atom_pairs in bond_types.items():
                bonds_found = cls._count_matching_bonds(redox_site, bond_type, atom_pairs)
                group_bonds_found += bonds_found
                
                group_details.append(TransformerEvaluationDetail(
                    description=f"{group_description} - {bond_type} bonds",
                    passed=bonds_found > 0,
                    value_found=bonds_found,
                    value_expected=f"part of {group_min_count} total"
                ))
            
            # Check if this group meets its requirements.
            #
            # The satisfied count is clamped to the requirement. Callers add
            # these numbers into the same running met/total tally they use for
            # pass/fail requirement checks, so an over-satisfied bond group
            # would otherwise contribute surplus credit that pays for a FAILED
            # composition check elsewhere — e.g. a site with 2 Cys scoring 6/7
            # on residue checks plus 4-found/3-required on bonds lands at
            # 10 == 10 and reads as fully valid. Clamping makes "met == total"
            # mean what it says: every check passed and every bond satisfied.
            group_passed = group_bonds_found >= group_min_count
            group_credit = min(group_bonds_found, group_min_count)
            group_results.append((group_passed, group_credit, group_min_count,
                                  group_bonds_found))
            
            # Add group summary
            details.append(TransformerEvaluationDetail(
                description=f"{group_description} (total)",
                passed=group_passed,
                value_found=group_bonds_found,
                value_expected=group_min_count
            ))
            
            # Add individual bond type details
            details.extend(group_details)
        
        # Calculate overall results
        if require_one_group:
            # At least one group must pass. Report that group's own numbers —
            # selecting among the groups that PASSED, not simply the one with
            # the largest count: a failing group can carry more bonds than a
            # passing one (a bigger min_count it fell short of), and picking it
            # would sink an evaluation that a satisfied group had already met.
            passing = [r for r in group_results if r[0]]
            if passing:
                best_group = max(passing, key=lambda r: r[3])
                return best_group[1], best_group[2], details
            else:
                # No group passed
                total_required = max(result[2] for result in group_results) if group_results else 0
                return 0, total_required, details
        else:
            # All groups must pass
            total_found = sum(result[1] for result in group_results)
            total_required = sum(result[2] for result in group_results)
            return total_found, total_required, details
    
    @classmethod
    def _count_matching_bonds(cls, redox_site, bond_type: str, atom_pairs: List[Tuple]) -> int:
        """
        Count bonds in RedoxSite that match the specified atom pairs
        
        Args:
            redox_site: RedoxSite object
            bond_type: Type of bond (coordinate, covalent, etc.)
            atom_pairs: List of ((resname1, atom1), (resname2, atom2)) tuples
            
        Returns:
            Number of matching bonds found
        """
        count = 0
        
        for bond in redox_site.bonds:
            # Only check bonds of the specified type
            if bond.chemical_type != bond_type:
                continue
                
            # Get atom info for both ends of the bond
            atom1_info = redox_site.get_current_pdb_info(bond.atom1_coords)
            atom2_info = redox_site.get_current_pdb_info(bond.atom2_coords)
            
            if not atom1_info or not atom2_info:
                continue
                
            # Check if this bond matches any of the required atom pairs
            for pair in atom_pairs:
                (resname1, atomname1), (resname2, atomname2) = pair
                
                # Check both directions (order agnostic)
                if cls._bond_matches_pair(atom1_info, atom2_info, resname1, atomname1, resname2, atomname2):
                    count += 1
                    break  # Don't double-count the same bond
                    
        return count
    
    @classmethod
    def _bond_matches_pair(cls, atom1_info: Dict, atom2_info: Dict, 
                          resname1: str, atomname1: str, resname2: str, atomname2: str) -> bool:
        """
        Check if a bond matches a specific residue/atom pair (order agnostic)
        
        Args:
            atom1_info, atom2_info: Atom information dictionaries
            resname1, atomname1: First residue/atom pair
            resname2, atomname2: Second residue/atom pair
            
        Returns:
            True if the bond matches the specified pair
        """
        # Forward direction
        if (atom1_info.get('resname') == resname1 and atom1_info.get('atom_name') == atomname1 and
            atom2_info.get('resname') == resname2 and atom2_info.get('atom_name') == atomname2):
            return True
            
        # Reverse direction (order agnostic)
        if (atom1_info.get('resname') == resname2 and atom1_info.get('atom_name') == atomname2 and
            atom2_info.get('resname') == resname1 and atom2_info.get('atom_name') == atomname1):
            return True
            
        return False

# ===== TRANSFORMER REGISTRY =====

class RedoxSiteTransformerRegistry:
    """Registry for all available RedoxSite transformers"""
    
    def __init__(self):
        self._transformers: Dict[str, type] = {}
    
    def register(self, transformer_class: type):
        """Register a transformer class"""
        if not issubclass(transformer_class, RedoxSiteTransformerBase):
            raise ValueError(f"Transformer must inherit from RedoxSiteTransformerBase")
        
        name = transformer_class.TRANSFORMER_NAME
        if name in self._transformers:
            logger.warning(f"Overriding existing transformer: {name}")
        
        self._transformers[name] = transformer_class
        logger.debug(f"Registered RedoxSite transformer: {name}")
    
    def get_transformer(self, name: str) -> Optional[type]:
        """Get transformer class by name"""
        return self._transformers.get(name)
    
    def get_all_transformers(self) -> Dict[str, type]:
        """Get all registered transformers"""
        return self._transformers.copy()
    
    def list_transformer_names(self) -> List[str]:
        """Get list of registered transformer names"""
        return list(self._transformers.keys())

# Global registry instance
redox_transformer_registry = RedoxSiteTransformerRegistry()

def register_redox_transformer(transformer_class: type):
    """Decorator to register a RedoxSite transformer"""
    redox_transformer_registry.register(transformer_class)
    return transformer_class

# ===== TRANSFORMER SELECTOR =====

class RedoxSiteTransformerSelector:
    """Central module that evaluates RedoxSites against all available transformers"""

    def __init__(self, console: Console = None, processor=None):
        self.console = console or Console()
        self.processor = processor
        self.registry = redox_transformer_registry
    
    def get_compatible_transformers(self, redox_site) -> List[Dict[str, Any]]:
        """
        Get all compatible transformers for a RedoxSite (without user interaction)

        Args:
            redox_site: RedoxSite object to evaluate

        Returns:
            List of dicts with transformer info: [{"name": str, "confidence": float, "evaluation": TransformerEvaluation}, ...]
        """
        compatible = []

        for transformer_name, transformer_class in self.registry.get_all_transformers().items():
            try:
                evaluation = transformer_class.evaluate_redox_site(redox_site)
                if evaluation.is_valid:
                    compatible.append({
                        "name": transformer_name,
                        "confidence": evaluation.confidence,
                        "evaluation": evaluation
                    })
            except Exception as e:
                logger.warning(f"Transformer {transformer_name} evaluation failed: {e}")

        return compatible

    def select_transformer_for_site(self, redox_site) -> Optional[str]:
        """
        Evaluate RedoxSite against all transformers and let user select

        Args:
            redox_site: RedoxSite object to evaluate

        Returns:
            Selected transformer name or None if no valid transformers
        """
        # Step 1: Get candidates from all transformers
        candidates = []

        for transformer_name, transformer_class in self.registry.get_all_transformers().items():
            try:
                # Each transformer handles its own evaluation
                evaluation = transformer_class.evaluate_redox_site(redox_site)
                if evaluation.is_valid:
                    candidates.append((transformer_name, transformer_class, evaluation))
            except Exception as e:
                self.console.print(f"[yellow]Warning: {transformer_name} evaluation failed: {e}[/yellow]")
                logger.warning(f"Transformer {transformer_name} evaluation failed: {e}")

        # Step 2: Present candidates to user
        if not candidates:
            self.console.print(f"[red]No transformers compatible with site {redox_site.site_id}[/red]")
            return None
        elif len(candidates) == 1:
            name, cls, eval_result = candidates[0]
            self.console.print(f"[green]Auto-selected {name} (confidence: {eval_result.confidence:.2f})[/green]")
            return name
        else:
            return self._interactive_selection(redox_site, candidates)
    
    def _interactive_selection(self, redox_site, candidates) -> str:
        """Present multiple candidates for user selection"""
        self.console.print(f"\n[bold]Multiple transformers available for site {redox_site.site_id}:[/bold]")
        
        table = Table(title="Compatible Transformers")
        table.add_column("#", style="cyan")
        table.add_column("Transformer", style="green") 
        table.add_column("Confidence", style="yellow")
        table.add_column("Description", style="blue")
        table.add_column("Requirements Met", style="magenta")
        
        for i, (name, cls, evaluation) in enumerate(candidates, 1):
            requirements_met = f"{evaluation.requirements_met}/{evaluation.total_requirements}"
            table.add_row(
                str(i), 
                name, 
                f"{evaluation.confidence:.2f}",
                evaluation.description,
                requirements_met
            )
        
        self.console.print(table)
        
        # Show detailed evaluation for each candidate
        for i, (name, cls, evaluation) in enumerate(candidates, 1):
            if evaluation.details:
                self.console.print(f"\n[bold]{i}. {name} Details:[/bold]")
                for detail in evaluation.details:
                    status = "[green]✓[/green]" if detail.passed else "[red]✗[/red]"
                    self.console.print(f"  {status} {detail.description}")
        
        # Build options_map for the choices
        options_map = {}
        for i, (name, cls, evaluation) in enumerate(candidates, 1):
            options_map[str(i)] = name

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Select transformer",
            choices=[str(i) for i in range(1, len(candidates) + 1)],
            default="1",
            module="Redox Site Transformation",
            description="Select transformer for site",
            options_map=options_map
        )

        selected_name = candidates[int(choice) - 1][0]
        self.console.print(f"[green]Selected: {selected_name}[/green]")
        return selected_name

# ===== SPACE ANALYZER =====

class RedoxSiteSpaceAnalyzer:
    """Analyzes RedoxSite transformation space requirements"""
    
    def __init__(self, console: Console = None):
        self.console = console or Console()
        self.registry = redox_transformer_registry
    
    def analyze_transformation_space(self, redox_sites: List, 
                                   selected_transformers: Dict[str, str]) -> Dict[str, Any]:
        """
        Analyze space requirements for all RedoxSite transformations
        
        Args:
            redox_sites: List of RedoxSite objects
            selected_transformers: site_id -> transformer_name mapping
            
        Returns:
            Analysis results including conflicts and space requirements
        """
        space_requirements = {}
        conflicts = []
        
        for site in redox_sites:
            transformer_name = selected_transformers.get(site.site_id)
            if not transformer_name:
                continue
                
            transformer_class = self.registry.get_transformer(transformer_name)
            if not transformer_class:
                continue
            
            # Get components for this site
            try:
                components, missing = transformer_class.match_components(site)
                if missing:
                    logger.warning(f"Missing components for site {site.site_id}: {missing}")
                    continue
            except Exception as e:
                logger.error(f"Component matching failed for site {site.site_id}: {e}")
                continue
            
            # Calculate space needed
            required_count = transformer_class.get_required_residue_count()
            space_plan = transformer_class.get_residue_space_plan(components)
            
            # Find the anchor residue (main center of the site)
            anchor_residue = self._find_anchor_residue(site, components)
            if not anchor_residue:
                continue
            
            anchor_id = anchor_residue['resid']
            anchor_chain = anchor_residue['chain']
            
            # Calculate required ID range
            max_offset = max(space_plan.values()) if space_plan else 0
            required_range = (anchor_id, anchor_id + max_offset)
            
            space_requirements[site.site_id] = {
                'anchor_chain': anchor_chain,
                'anchor_id': anchor_id,
                'required_count': required_count,
                'required_range': required_range,
                'space_plan': space_plan,
                'transformer': transformer_name,
                'components': components
            }
        
        # Detect conflicts
        conflicts = self._detect_space_conflicts(space_requirements)
        
        return {
            'space_requirements': space_requirements,
            'conflicts': conflicts,
            'needs_id_mapping': len(conflicts) > 0
        }
    
    def _find_anchor_residue(self, redox_site, components: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Find the anchor residue for space calculation - works for any center type"""
        # Look for various center component patterns
        center_patterns = [
            ("center_id", "center_chain"),
            ("site_id", "site_chain"), 
            ("main_residue_id", "main_residue_chain"),
            ("metal_center_id", "metal_center_chain"),  # Legacy metal support
            ("cofactor_id", "cofactor_chain"),
            ("residue_id", "residue_chain")
        ]
        
        for id_key, chain_key in center_patterns:
            if id_key in components and chain_key in components:
                return {
                    'resid': components[id_key],
                    'chain': components[chain_key]
                }
        
        # Fall back to first center of any type
        if redox_site.centers:
            center = redox_site.centers[0]
            return {
                'resid': center.resid,
                'chain': center.chain
            }
        
        return None
    
    def _detect_space_conflicts(self, space_requirements: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Detect overlapping residue ID requirements"""
        conflicts = []
        
        # Group by chain
        by_chain = {}
        for site_id, req in space_requirements.items():
            chain = req['anchor_chain']
            if chain not in by_chain:
                by_chain[chain] = []
            by_chain[chain].append((site_id, req))
        
        # Check for overlaps within each chain
        for chain, chain_reqs in by_chain.items():
            for i, (site1_id, req1) in enumerate(chain_reqs):
                for site2_id, req2 in chain_reqs[i+1:]:
                    range1 = req1['required_range']
                    range2 = req2['required_range']
                    
                    # Check if ranges overlap
                    if self._ranges_overlap(range1, range2):
                        conflicts.append({
                            'type': 'site_overlap',
                            'site1': site1_id,
                            'site2': site2_id,
                            'chain': chain,
                            'range1': range1,
                            'range2': range2
                        })
        
        return conflicts
    
    def _ranges_overlap(self, range1: Tuple[int, int], range2: Tuple[int, int]) -> bool:
        """Check if two ID ranges overlap"""
        start1, end1 = range1
        start2, end2 = range2
        return not (end1 < start2 or end2 < start1)

# ===== TRANSFORMATION EXECUTOR =====

class TransformationExecutor:
    """
    Executes transformation sequences using RedoxSite coordinate mapping

    Args:
        console: Rich Console for output (optional)
        temp_dir: Directory for temporary files (optional)
        verbose: Enable verbose output (default: True)
        save_intermediate_structures: Save PDB after each transformation step for debugging (default: False)
                                      Warning: Generates many files - only enable when debugging specific issues
    """

    def __init__(self, console: Console = None, temp_dir = None, verbose: bool = True, save_intermediate_structures: bool = False):
        self.console = console or Console()
        self.temp_dir = temp_dir
        self.verbose = verbose
        self.save_intermediate_structures = save_intermediate_structures
    
    def apply_transformation_sequence(self, redox_site, 
                                    transformations: List[Dict[str, Any]], 
                                    pdb_lines: List[str]) -> Tuple[List[str], Any]:
        """
        Apply ordered transformation sequence to PDB structure
        
        Args:
            redox_site: RedoxSite with coordinate mapping
            transformations: Ordered list from transformer
            pdb_lines: Current PDB content
            
        Returns:
            (modified_pdb_lines, updated_redox_site)
        """
        # Set current site ID for intermediate structure saving
        self.current_site_id = getattr(redox_site, 'site_id', 'unknown')
        
        current_lines = pdb_lines.copy()
        
        for i, transform in enumerate(transformations):
            logger.debug(f"Applying transformation {i+1}/{len(transformations)}: {transform['description']}")
            
            # Apply transformation using coordinate-based selection
            # (RedoxSite metadata is updated within _apply_single_transformation)
            current_lines, lines_modified = self._apply_single_transformation(
                transform, redox_site, current_lines
            )
            
            # Save intermediate structure after each transformation step for debugging
            self._save_intermediate_structure(current_lines, redox_site, i+1, transform['description'])
            
            # Display progress
            if self.verbose:
                self.console.print(f"⚙️  Applying {transform['description']} ({i+1}/{len(transformations)}, {(i+1)/len(transformations)*100:.1f}%)")
                if lines_modified > 0:
                    self.console.print(f"   → Modified {lines_modified} lines (total: {sum(t.get('lines_modified', 0) for t in transformations[:i+1])})")
                else:
                    self.console.print("   → No changes made")
        
        return current_lines, redox_site
    
    def _apply_single_transformation(self, transform: Dict[str, Any], 
                                   redox_site, pdb_lines: List[str]) -> Tuple[List[str], int]:
        """Apply a single transformation using RedoxSite coordinate mapping"""
        
        # Get target coordinates from RedoxSite based on selector
        target_coords = self._resolve_transformation_targets(transform["selector"], redox_site)
        
        
        if not target_coords:
            logger.debug(f"No target coordinates found for transformation: {transform['description']}")
            return pdb_lines, 0
        
        # Apply action to all matching PDB lines and track coordinate updates
        modified_lines = []
        lines_modified = 0
        first_few_checked = 0
        coord_updates = {}  # Track what coordinates got new metadata
        
        for line in pdb_lines:
            if line.startswith(('ATOM', 'HETATM')):
                line_coords = self._extract_coordinates_from_pdb_line(line)
                
                # DEBUG: Print first few coordinate comparisons
                # Debug output removed for cleaner interface
                
                # Use coordinate-based matching with tolerance (like original MetalSite system)
                if self._coordinate_matches_any_target(line_coords, target_coords):
                    # Apply the transformation action
                    modified_line = self._apply_action_to_pdb_line(line, transform["action"])
                    modified_lines.append(modified_line)
                    lines_modified += 1
                    
                    # Extract new metadata from the modified line
                    new_metadata = self._extract_metadata_from_pdb_line(modified_line)
                    coord_updates[line_coords] = new_metadata
                else:
                    # DEBUG: Show coordinate mismatch for first few lines
                    if first_few_checked < 3 and line_coords:
                        # Show closest target coord
                        if target_coords:
                            closest_target = min(target_coords, key=lambda t: self._calculate_distance(line_coords, t))
                            distance = self._calculate_distance(line_coords, closest_target)
                        first_few_checked += 1
                    modified_lines.append(line)
            else:
                modified_lines.append(line)
        
        # Update RedoxSite metadata immediately after transformation
        if coord_updates:
            if self.verbose:
                self.console.print(f"[cyan]🔄 Updating RedoxSite metadata for {len(coord_updates)} coordinates[/cyan]")
            
            # DEBUG: Show what residue groups exist before update
            before_groups = list(redox_site.residue_groups.keys())
            # self.console.print(f"[grey50]🔍 DEBUG: Residue groups before update: {before_groups}[/grey50]")
            
            # DEBUG: Check if coordinates are found in residue groups BEFORE update
            for coords, metadata in coord_updates.items():
                found_in_group = None
                for group_key, coord_list in redox_site.residue_groups.items():
                    if coords in coord_list:
                        found_in_group = group_key
                        break
                if found_in_group:
                    pass  # Debug output suppressed
                else:
                    pass  # Debug output suppressed
            
            # Fix coordinate precision mismatch by mapping to closest stored coordinates
            corrected_coord_updates = {}
            for coords, metadata in coord_updates.items():
                # Find the closest stored coordinate in RedoxSite
                closest_stored_coord = self._find_closest_stored_coordinate(coords, redox_site)
                if closest_stored_coord:
                    corrected_coord_updates[closest_stored_coord] = metadata
                    if closest_stored_coord != coords:
                        pass  # Silently correct coordinate precision
                else:
                    self.console.print(f"[red]❌ No matching stored coordinate found for {coords}[/red]")
            
            if corrected_coord_updates:
                redox_site.update_atom_metadata(corrected_coord_updates)
                # Silently update RedoxSite coordinates
            else:
                self.console.print(f"[red]❌ No coordinates could be updated due to precision mismatch[/red]")
            
            # DEBUG: Show what residue groups exist after update
            after_groups = list(redox_site.residue_groups.keys())
            # self.console.print(f"[grey50]🔍 DEBUG: Residue groups after update: {after_groups}[/grey50]")
            
            # DEBUG: Check if coordinates are found in residue groups AFTER update
            for coords, metadata in coord_updates.items():
                found_in_group = None
                for group_key, coord_list in redox_site.residue_groups.items():
                    if coords in coord_list:
                        found_in_group = group_key
                        break
                if found_in_group:
                    pass  # Debug output suppressed
                else:
                    pass  # Debug output suppressed
                    
                # self.console.print(f"[grey50]🔍 DEBUG: Updated {coords} → chain={metadata.get('chain')}, resname={metadata.get('resname')}, resid={metadata.get('resid')}[/grey50]")
        else:
            self.console.print(f"[yellow]⚠️ No coordinate updates to apply[/yellow]")
        
        # Store modification count for progress tracking
        transform['lines_modified'] = lines_modified
        
        return modified_lines, lines_modified
    
    def _save_intermediate_structure(self, pdb_lines: List[str], redox_site, step_num: int, description: str):
        """Save intermediate structure after each transformation step for debugging"""
        # Skip if intermediate structure saving is disabled
        if not self.save_intermediate_structures:
            return

        try:
            # Debug: Print temp_dir status
            if not hasattr(self, 'temp_dir') or self.temp_dir is None:
                self.console.print(f"[red]ERROR: temp_dir not set for intermediate structure saving![/red]")
                return
            
            # Create step-specific filename
            safe_description = "".join(c for c in description if c.isalnum() or c in (' ', '-', '_')).rstrip()
            safe_description = "_".join(safe_description.split())[:50]  # Limit length
            filename = f"step_{step_num:02d}_{safe_description}.pdb"
            
            # Get site-specific directory from current transformation
            if hasattr(self, 'current_site_id') and self.current_site_id:
                site_dir = self.temp_dir / "intermediate_structures" / self.current_site_id
            else:
                site_dir = self.temp_dir / "intermediate_structures" / "unknown_site"
            
            site_dir.mkdir(parents=True, exist_ok=True)
            output_path = site_dir / filename
            
            # Write structure to file (fix formatting - lines already have newlines)
            with open(output_path, 'w') as f:
                for line in pdb_lines:
                    if line.endswith('\n'):
                        f.write(line)
                    else:
                        f.write(line + '\n')
            
            # Silently save intermediate structure for debugging
            
            # Also save a brief summary
            summary_path = site_dir / f"step_{step_num:02d}_summary.txt"
            with open(summary_path, 'w') as f:
                f.write(f"Transformation Step {step_num}\n")
                f.write(f"Description: {description}\n")
                f.write(f"Total lines: {len(pdb_lines)}\n")
                f.write(f"ATOM/HETATM lines: {sum(1 for line in pdb_lines if line.startswith(('ATOM', 'HETATM')))}\n")
                f.write(f"Site ID: {getattr(redox_site, 'site_id', 'unknown')}\n")
                f.write(f"Structure saved to: {filename}\n")
            
        except Exception as e:
            # Don't fail transformation if saving fails
            self.console.print(f"[red]ERROR saving intermediate structure for step {step_num}: {e}[/red]")
            logger.warning(f"Failed to save intermediate structure for step {step_num}: {e}")

    def cleanup_intermediate_structures(self):
        """Clean up all intermediate structure files to free disk space"""
        if not hasattr(self, 'temp_dir') or self.temp_dir is None:
            return

        try:
            import shutil
            intermediate_dir = self.temp_dir / "intermediate_structures"

            if intermediate_dir.exists():
                # Count files before deletion for reporting
                file_count = sum(1 for _ in intermediate_dir.rglob('*') if _.is_file())

                # Remove the entire intermediate structures directory
                shutil.rmtree(intermediate_dir)

                if self.verbose:
                    self.console.print(f"[grey50]🧹 Cleaned up {file_count} intermediate structure files[/grey50]")
                logger.info(f"Cleaned up {file_count} intermediate structure files from {intermediate_dir}")
        except Exception as e:
            # Don't fail if cleanup fails, just log it
            logger.warning(f"Failed to cleanup intermediate structures: {e}")
            if self.verbose:
                self.console.print(f"[yellow]⚠️  Failed to cleanup intermediate structures: {e}[/yellow]")

    def _resolve_transformation_targets(self, selector: Dict, redox_site) -> Set[Tuple[float, float, float]]:
        """Convert selector criteria to coordinate sets using RedoxSite data"""
        target_coords = set()

        # Handle different selector types
        if 'chain_id' in selector and 'residue_id' in selector:
            # Select by residue
            chain_id = selector['chain_id']
            residue_id = selector['residue_id']
            insertion_code = selector.get('insertion_code', '')
            
            # DEBUG: Print what we're looking for
            # Debug output suppressed
            
            # DEBUG: Print available residues in RedoxSite
            available_residues = list(redox_site.residue_groups.keys())
            # self.console.print(f"[grey50]🔍 DEBUG: Available residues in RedoxSite: {available_residues}[/grey50]")
            
            # DEBUG: Show detailed atom list for the target residue
            # self.console.print(f"[grey50]🔍 DEBUG: Atoms currently in target residue ({chain_id}, {residue_id}):[/grey50]")
            target_atoms = redox_site.get_atoms_by_residue(chain_id, residue_id, insertion_code)
            if not target_atoms and insertion_code == '':
                target_atoms = redox_site.get_atoms_by_residue(chain_id, residue_id, ' ')
            # Detailed atom listing suppressed for cleaner output
            
            # Try to get atoms with exact insertion code first
            atoms = redox_site.get_atoms_by_residue(chain_id, residue_id, insertion_code)
            
            # If not found and insertion code is empty, try with space (common PDB format difference)
            if not atoms and insertion_code == '':
                atoms = redox_site.get_atoms_by_residue(chain_id, residue_id, ' ')
                if atoms:
                    pass  # Debug output suppressed
            
            # If still not found and insertion code is space, try with empty
            if not atoms and insertion_code == ' ':
                atoms = redox_site.get_atoms_by_residue(chain_id, residue_id, '')
                if atoms:
                    pass  # Debug output suppressed
            
            
            # Filter by residue name if specified
            if 'residue_name' in selector:
                expected_resname = selector['residue_name']
                # Use the current residue name from coord_to_pdb mapping, not atom.resname
                # This handles atoms that have been moved between residues
                filtered_atoms = []
                for atom in atoms:
                    current_metadata = redox_site.coord_to_pdb.get(atom.coords, {})
                    current_resname = current_metadata.get('resname', atom.resname)
                    if current_resname == expected_resname:
                        filtered_atoms.append(atom)
                atoms = filtered_atoms
            
            # Filter by atom names if specified
            if 'atom_names' in selector:
                allowed_names = set(selector['atom_names'])
                atoms = [atom for atom in atoms if atom.atom_name in allowed_names]
            
            # Collect coordinates
            atom_coords = [atom.coords for atom in atoms]
            target_coords.update(atom_coords)
            
            # DEBUG: Print coordinates we're looking for
            # self.console.print(f"[grey50]🔍 DEBUG: Target coordinates: {list(atom_coords)[:3]}{'...' if len(atom_coords) > 3 else ''}[/grey50]")
        
        
        # DEBUG: If this is step 8 (apply_redox_specific_heme_name), show ALL atoms in the target residue
        if 'residue_name' in selector and selector.get('residue_name') == 'HEC' and 'residue_id' in selector:
            target_resid = selector['residue_id']
            target_chain = selector['chain_id']
            count = 0
            for coord, atom_info in redox_site.coord_to_pdb.items():
                if (atom_info.get('chain') == target_chain and 
                    atom_info.get('resid') == target_resid):
                    count += 1
        
        return target_coords
    
    def _coordinate_matches_any_target(self, line_coords: Tuple[float, float, float], 
                                     target_coords: Set[Tuple[float, float, float]], 
                                     tolerance: float = 0.001) -> bool:
        """Check if line coordinates match any target coordinates within tolerance"""
        for target in target_coords:
            if self._calculate_distance(line_coords, target) <= tolerance:
                return True
        return False
    
    def _calculate_distance(self, coords1: Tuple[float, float, float], 
                          coords2: Tuple[float, float, float]) -> float:
        """Calculate Euclidean distance between two coordinate tuples"""
        import numpy as np
        return np.sqrt(np.sum((np.array(coords1) - np.array(coords2)) ** 2))
    
    def _find_closest_stored_coordinate(self, target_coords: Tuple[float, float, float], 
                                      redox_site, tolerance: float = 0.001) -> Optional[Tuple[float, float, float]]:
        """Find the closest stored coordinate in RedoxSite within tolerance"""
        best_coord = None
        best_distance = float('inf')
        
        # Search through all stored coordinates in RedoxSite
        all_stored_coords = set()
        for atom in redox_site.atoms:
            all_stored_coords.add(atom.coords)
        for center in redox_site.centers:
            all_stored_coords.add(center.coords)
        
        for stored_coord in all_stored_coords:
            distance = self._calculate_distance(target_coords, stored_coord)
            if distance < best_distance and distance <= tolerance:
                best_distance = distance
                best_coord = stored_coord
        
        return best_coord
    
    def _extract_coordinates_from_pdb_line(self, line: str) -> Tuple[float, float, float]:
        """Extract coordinates from PDB line as raw floats (no rounding)"""
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip()) 
            z = float(line[46:54].strip())
            return (x, y, z)  # Keep original precision for distance-based matching
        except (ValueError, IndexError):
            return (0.0, 0.0, 0.0)
    
    def _extract_metadata_from_pdb_line(self, line: str) -> Dict[str, Any]:
        """Extract atom metadata from PDB line for RedoxSite updates"""
        try:
            # Extract insertion code properly - preserve original values
            # Now that get_atoms_by_residue treats '' and ' ' as equivalent for no-insertion cases,
            # we preserve the original insertion codes rather than normalizing them.
            insertion_code = line[26] if len(line) > 26 else ' '
            metadata = {
                'chain': line[21],
                'resname': line[17:20].strip(),
                'resid': int(line[22:26].strip()),
                'insertion_code': insertion_code,  # Keep space as space for consistency
                'atom_name': line[12:16].strip()
            }
            return metadata
        except (ValueError, IndexError) as e:
            self.console.print(f"[red]ERROR: Failed to extract metadata from PDB line: {e}[/red]")
            self.console.print(f"[red]Line: {repr(line)}[/red]")
            return {}
    
    def _apply_action_to_pdb_line(self, line: str, action: Dict[str, Any]) -> str:
        """Apply transformation action to a PDB line"""
        modified_line = line
        
        # Apply residue name change
        if 'change_residue_name' in action:
            new_resname = action['change_residue_name']
            old_resname = line[17:20].strip()
            atom_name = line[12:16].strip()
            resid = line[22:26].strip()
            chain = line[21:22].strip()
            
            # DEBUG: Print residue name change details
            
            modified_line = modified_line[:17] + f"{new_resname:>3}" + modified_line[20:]
        
        # Apply residue ID change
        if 'change_residue_id' in action:
            new_resid = action['change_residue_id']
            modified_line = modified_line[:22] + f"{new_resid:>4}" + modified_line[26:]
        
        # Apply chain ID change
        if 'change_chain_id' in action:
            new_chain = action['change_chain_id']
            modified_line = modified_line[:21] + new_chain + modified_line[22:]
        
        # Apply insertion code change
        if 'change_insertion_code' in action:
            new_insertion = action['change_insertion_code']
            old_insertion = line[26:27]
            atom_name = line[12:16].strip()
            modified_line = modified_line[:26] + new_insertion + modified_line[27:]
        
        # Apply atom name changes
        if 'rename_atoms' in action:
            current_atom_name = line[12:16].strip()
            atom_mapping = action['rename_atoms']
            if current_atom_name in atom_mapping:
                new_atom_name = atom_mapping[current_atom_name]
                modified_line = modified_line[:12] + f"{new_atom_name:>4}" + modified_line[16:]
        
        # Convert ATOM to HETATM if requested
        if action.get('convert_to_hetatm', False) and line.startswith('ATOM'):
            modified_line = 'HETATM' + modified_line[6:]
        
        return modified_line
    
