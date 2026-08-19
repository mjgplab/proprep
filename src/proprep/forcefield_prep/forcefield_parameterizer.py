"""
Force Field Parameterizer for MD Simulations

A module for handling forcefield parameterization of non-standard residues for MD simulations.
Acts as an interface between the PDB processor and specialized parameterization modules.
"""

import glob
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from collections import defaultdict, OrderedDict
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from rich.console import Console
from rich.panel import Panel
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from rich.table import Table

from proprep.application.menu_commands import prompt_with_context
import proprep.structure_prep.chem_comp_dict_fetcher as ccd
from proprep.structure_prep.comprehensive_redox_detector import METALS, CenterType
from proprep.utils.module_registry import ProcessingModule, register_module
from .forcefield_worker import NonStandardResidue
from .forcefield_commands import (
    AnalyzeNonstandardResiduesCommand,
    ParameterizeAACommand,
    ReclassifyResidueCommand,
    ManageMappingsCommand,
    ConfigureClassificationSettingsCommand,
    DisplayHelpCommand,
)

from . import modified_amino_acid_parameterizer
from . import metal_site_parameterizer

logger = logging.getLogger(__name__)

def _modaa_output_dirname(residue_name):
    """Working-directory name for one modified-amino-acid parameterization.

    Matches the sibling parameterizers' convention — ``small_molecule_params_E4Z``,
    ``metal_site_params_MN_A_202`` — so all three read alike in a run directory.
    (The modified-AA route used to nest under ``parameterized_residues/<name>/``,
    which both broke that convention and collided by name with the
    ``parameterized_residues`` WORKSPACE key, a dict of parameter data that has
    nothing to do with this directory.)
    """
    return f"modified_aa_params_{residue_name}"


def _imported_library_atom_names(result: dict) -> set:
    """Non-hydrogen atom names in the library that was just deposited."""
    from proprep.forcefield_prep.library_promotion import library_atom_names

    state_dir = result.get("state_dir") or result.get("library_path")
    if not state_dir:
        return set()
    names = set()
    root = Path(state_dir)
    for lib in sorted(root.glob("*.lib")) + sorted(root.glob("*.off")):
        names |= library_atom_names(lib)
    return names


def _rank_sites_for_library(redox_sites, residue_name, lib_atom_names):
    """Detected sites ranked by how well they match an imported library.

    Returns ``[(site, matched_resname, overlap_fraction), ...]``, best first.

    Matching on the library's RESIDUE NAME alone is wrong in the very case a
    transformer exists for. Parameters for an oxidized FAD might be named FAO
    while the structure says FAD -- renaming FAD -> FAO is the transformer's
    whole purpose -- so looking for a site containing FAO finds nothing and
    would send the user off to define a site that cannot exist.

    Atom names survive the rename. An FAO library derived from FAD still
    carries FAD's atom names, so composition identifies the target residue
    whatever it is called. Verified on a real import: the deposited FAD
    library's 53 non-hydrogen names matched the structure's FAD residue exactly.

    An outright name match short-circuits, being the strongest signal when the
    names do agree.
    """
    wanted = str(residue_name).strip().upper() if residue_name else None
    ranked = []

    for site in redox_sites:
        by_residue = {}
        for atom in getattr(site, "atoms", None) or []:
            key = (getattr(atom, "resname", "") or "").strip().upper()
            if key:
                by_residue.setdefault(key, set()).add(
                    (getattr(atom, "atom_name", "") or "").strip().upper())

        best_name, best_score = None, 0.0
        for resname, atom_names in by_residue.items():
            if wanted and resname == wanted:
                best_name, best_score = resname, 1.0
                break
            if lib_atom_names and atom_names:
                score = len(atom_names & lib_atom_names) / max(len(atom_names), 1)
                if score > best_score:
                    best_name, best_score = resname, score
        if best_name and best_score > 0:
            ranked.append((site, best_name, best_score))

    ranked.sort(key=lambda item: item[2], reverse=True)
    return ranked


@register_module
class ForcefieldParameterizer(ProcessingModule):
    """Module for parameterizing non-standard residues for MD simulations"""

    NAME = "Force Field Parameterizer"
    DESCRIPTION = "Parameterize organic molecules, modified amino acids, and metal sites"
    VERSION = "1.0.0"
    CATEGORY = "preparation"
    REQUIRES = ["PDB Loader"]
    PRIORITY = 5

    def __init__(self):
        """Initialize the forcefield parameterizer module"""
        super().__init__()
        self.console = Console()
        self.ccd_parser = ccd.CCDParser(use_cache=True)
        self.non_standard_residues = []

        self.classification_settings = {
            "max_small_molecule_atoms": 200,
            "min_small_molecule_atoms": 2,
        }

        self.standard_aa = {
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        }

        self.modified_aa_map = {
            "2AS": "ASP",
            "3AH": "HIS",
            "5HP": "GLU",
            "5OW": "LYS",
            "ACL": "ARG",
            "AGM": "ARG",
            "AIB": "ALA",
            "ALM": "ALA",
            "ALO": "THR",
            "ALY": "LYS",
            "ARM": "ARG",
            "ASA": "ASP",
            "ASB": "ASP",
            "ASK": "ASP",
            "ASL": "ASP",
            "ASQ": "ASP",
            "ASH": "ASP",  # Protonated aspartate
            "AYA": "ALA",
            "BCS": "CYS",
            "BHD": "ASP",
            "BMT": "THR",
            "BNN": "ALA",
            "BUC": "CYS",
            "BUG": "LEU",
            "C5C": "CYS",
            "C6C": "CYS",
            "CAS": "CYS",
            "CCS": "CYS",
            "CEA": "CYS",
            "CGU": "GLU",
            "CHG": "ALA",
            "CLE": "LEU",
            "CME": "CYS",
            "CSD": "ALA",
            "CSO": "CYS",
            "CSP": "CYS",
            "CSS": "CYS",
            "CSW": "CYS",
            "CSX": "CYS",
            "CXM": "MET",
            "CY1": "CYS",
            "CY3": "CYS",
            "CYG": "CYS",
            "CYM": "CYS",
            "CYQ": "CYS",
            "DAH": "PHE",
            "DAL": "ALA",
            "DAR": "ARG",
            "DAS": "ASP",
            "DCY": "CYS",
            "DGL": "GLU",
            "DGN": "GLN",
            "DHA": "ALA",
            "DHI": "HIS",
            "DIL": "ILE",
            "DIV": "VAL",
            "DLE": "LEU",
            "DLY": "LYS",
            "DNP": "ALA",
            "DPN": "PHE",
            "DPR": "PRO",
            "DSN": "SER",
            "DSP": "ASP",
            "DTH": "THR",
            "DTR": "TRP",
            "DTY": "TYR",
            "DVA": "VAL",
            "EFC": "CYS",
            "FLA": "ALA",
            "FME": "MET",
            "GGL": "GLU",
            "GL3": "GLY",
            "GLZ": "GLY",
            "GMA": "GLU",
            "GLH": "GLU",  # Protonated glutamate
            "GSC": "GLY",
            "HAC": "ALA",
            "HAR": "ARG",
            "HIC": "HIS",
            "HID": "HIS",  # Histidine with hydrogen on delta nitrogen
            "HIE": "HIS",  # Histidine with hydrogen on epsilon nitrogen
            "HIP": "HIS",  # Doubly protonated histidine
            "HMR": "ARG",
            "HPQ": "PHE",
            "HTR": "TRP",
            "HYP": "PRO",
            "IAS": "ASP",
            "IIL": "ILE",
            "IYR": "TYR",
            "KCX": "LYS",
            "LLP": "LYS",
            "LLY": "LYS",
            "LTR": "TRP",
            "LYM": "LYS",
            "LYZ": "LYS",
            "MAA": "ALA",
            "MEN": "ASN",
            "MHS": "HIS",
            "MIS": "SER",
            "MK8": "LEU",
            "MLE": "LEU",
            "MPQ": "GLY",
            "MSA": "GLY",
            "MSE": "MET",
            "MVA": "VAL",
            "NEM": "HIS",
            "NEP": "HIS",
            "NLE": "LEU",
            "NLN": "LEU",
            "NLP": "LEU",
            "NMC": "GLY",
            "OAS": "SER",
            "OCS": "CYS",
            "OMT": "MET",
            "PAQ": "TYR",
            "PCA": "GLU",
            "PEC": "CYS",
            "PHI": "PHE",
            "PHL": "PHE",
            "PR3": "CYS",
            "PRR": "ALA",
            "PTR": "TYR",
            "PYX": "CYS",
            "SAC": "SER",
            "SAR": "GLY",
            "SCH": "CYS",
            "SCS": "CYS",
            "SCY": "CYS",
            "SEL": "SER",
            "SEP": "SER",
            "SET": "SER",
            "SHC": "CYS",
            "SHR": "LYS",
            "SMC": "CYS",
            "SOC": "CYS",
            "STY": "TYR",
            "SVA": "SER",
            "TIH": "ALA",
            "TPL": "TRP",
            "TPO": "THR",
            "TPQ": "ALA",
            "TRG": "LYS",
            "TRO": "TRP",
            "TYB": "TYR",
            "TYI": "TYR",
            "TYQ": "TYR",
            "TYS": "TYR",
            "TYY": "TYR",
        }

        # Metal elements - use comprehensive set from comprehensive_redox_detector
        # METALS is imported at module level from proprep.structure_prep.comprehensive_redox_detector

        # User-defined residue classifications (persistent across sessions)
        self.user_residue_classifications = {}
        
        self._load_classification_settings()
        self._load_user_residue_classifications()

        # Load user-defined mappings
        self._load_modified_aa_map()

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Helper methods to centralize interactions with the workspace
    def get_workspace(self):
        """
        Helper method to get the appropriate workspace object

        Returns:
            The current workspace object
        """
        return self.processor.workspace

    def get_from_workspace(self, key, default=None):
        """
        Helper method to get values from the processor's workspace

        Args:
            key: Key to retrieve
            default: Default value if key not found

        Returns:
            The value for the given key or default
        """
        return self.processor.workspace.get(key, default)

    # Priority order for forcefield parameterization (tLEaP convention)
    _PRIORITY_KEYS = [
        "transformed_pdb_file",
        "protonation_pdb_file",
        "structure_with_prot_resnames",  # Legacy key
        "repaired_pdb_file",
        "filtered_pdb_file",
        "hstripped_pdb_file",
        "rcsb_pdb_file",
        "local_pdb_file",
        "alphafold_pdb_file",
        "alphafill_pdb_file",
        "alphafold_homolog_pdb_file",
    ]

    # Mapping from workspace key to display name
    _KEY_TO_NAME = {
        "transformed_pdb_file": "transformed",
        "protonation_pdb_file": "protonated",
        "structure_with_prot_resnames": "protonated",
        "repaired_pdb_file": "repaired",
        "filtered_pdb_file": "filtered",
        "hstripped_pdb_file": "H-stripped",
        "rcsb_pdb_file": "RCSB",
        "local_pdb_file": "local",
        "alphafold_pdb_file": "AlphaFold",
        "alphafill_pdb_file": "AlphaFill",
        "alphafold_homolog_pdb_file": "AlphaFold homolog",
    }

    def _get_structure_with_priority(self):
        """
        Get structure from workspace using priority system.

        Priority (following tLEaP convention):
        1. transformed_pdb_file - After redox transformation
        2. protonation_pdb_file - After protonation state assignment
        3. repaired_pdb_file - After structure repair (MODELLER)
        4. filtered_pdb_file - After PDB filtering
        5. Structure Loader files - Direct from loading (rcsb, local, alphafold, alphafill)

        Checks file paths first (preferred), then falls back to BioPython Structure objects.
        Uses StructureSelector for consistent selection behavior.

        Returns:
            Tuple of (structure, structure_type_name, is_file_path)
        """
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(
            self.processor.workspace, self.console, processor=self.processor
        )

        # Try file paths first (preferred for forcefield generation)
        result = selector.get_structure(
            priority_override=self._PRIORITY_KEYS,
            return_key=True,
            silent=True,
        )

        if result is not None:
            structure_file, workspace_key = result
            name = self._KEY_TO_NAME.get(workspace_key, "selected")
            self.console.print(f"[green]Using {name} structure: {structure_file}[/green]")
            return structure_file, name, True

        # Fallback to BioPython Structure objects
        result = selector.get_structure_object(
            priority_override=self._PRIORITY_KEYS,
            return_key=True,
            silent=True,
        )

        if result is not None:
            structure, structure_key = result
            # Map structure key to display name
            file_key = structure_key.replace("_structure", "_pdb_file")
            name = self._KEY_TO_NAME.get(file_key, "selected")
            self.console.print(f"[yellow]Using {name} structure object (BioPython)[/yellow]")
            return structure, name, False

        return None, None, False

    def _get_structure_file_for_preprocessing(self) -> str:
        """
        Get a PDB file path for structure preprocessing.

        Uses same priority as _get_structure_with_priority() but only returns
        file paths (not BioPython objects). Uses StructureSelector for consistent
        selection behavior.

        Returns:
            PDB file path or None
        """
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(
            self.processor.workspace, self.console, processor=self.processor
        )

        result = selector.get_structure(
            priority_override=self._PRIORITY_KEYS,
            return_key=True,
            silent=True,
        )

        if result is not None:
            structure_file, workspace_key = result
            name = self._KEY_TO_NAME.get(workspace_key, "selected")
            self.console.print(
                f"[green]Using {name} structure for preprocessing: {structure_file}[/green]"
            )
            return structure_file

        self.console.print("[yellow]No PDB file found in workspace[/yellow]")
        return None

    def _get_structure_object(self):
        """
        Get a BioPython Structure object from workspace using priority system.

        This is a convenience method for code that needs direct structure access.
        Priority order matches _get_structure_with_priority().
        Uses StructureSelector for consistent selection behavior.

        Returns:
            BioPython Structure object or None
        """
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(
            self.processor.workspace, self.console, processor=self.processor
        )

        result = selector.get_structure_object(
            priority_override=self._PRIORITY_KEYS,
            silent=True,
        )

        return result

    def _extract_nonstandard_from_redox_sites(self):
        """
        Extract non-standard residues from RedoxSite objects.

        RedoxSites contain redox-active centers which are often non-standard residues
        (hemes, metal ions, cofactors, modified amino acids, etc.)

        Returns:
            List of NonStandardResidue objects extracted from RedoxSites
        """
        redox_sites = self.get_from_workspace("detected_redox_sites", [])

        if not redox_sites:
            return []


        # Track unique non-standard residues by (chain, resid, resname)
        unique_residues = {}

        for site in redox_sites:
            # Extract from centers (these are the key redox-active components)
            for center in site.centers:
                # Skip standard amino acids
                if center.resname in self.standard_aa:
                    continue

                key = (center.chain, center.resid, center.resname)

                if key not in unique_residues:
                    # Category is decided later by content (_classify_unit); the
                    # detector's center_type is only a hint, never the authority.
                    ns_res = NonStandardResidue(
                        name=center.resname,
                        chain_id=center.chain,
                        resid=center.resid,
                        category="unknown",
                    )

                    # Link to source RedoxSite for MCPB parameterization
                    # The full RedoxSite contains metal + coordinating residues
                    ns_res.redox_site_id = site.site_id
                    ns_res.source_redox_site = site  # Full RedoxSite object
                    ns_res.redox_center_data = center
                    ns_res.is_redox_component = True

                    # Store elements from RedoxSite atoms for this residue
                    # This allows metal detection even without BioPython residue
                    residue_elements = set()
                    for atom in site.atoms:
                        if (atom.chain == center.chain and
                            atom.resid == center.resid and
                            atom.resname == center.resname):
                            if atom.element:
                                residue_elements.add(atom.element.upper())
                    ns_res.redox_site_elements = residue_elements

                    unique_residues[key] = ns_res

        return list(unique_residues.values())
    def get_from_workspace_obj(self, workspace, key, default=None):
        """
        Helper method to get values from a workspace object

        Args:
            workspace: Workspace object
            key: Key to retrieve
            default: Default value if key not found

        Returns:
            The value for the given key or default
        """
        return workspace.get(key, default)

    def update_workspace(self, key, value):
        """
        Helper method to update the processor's workspace

        Args:
            key: Key to update
            value: New value
        """
        old_value = self.processor.workspace.get(key)
        self.processor.workspace.set(key, value)
        debug_enabled = self.processor.workspace.get("debug", False)

        # Debug output if enabled
        if debug_enabled and self.processor.console:
            self._debug_value(key, value, old_value, "updated")

    def update_workspace_obj(self, workspace, key, value):
        """
        Helper method to update a workspace object

        Args:
            workspace: Workspace object
            key: Key to update
            value: New value

        Returns:
            Updated workspace object
        """
        debug_enabled = workspace.get("debug", False)
        old_value = workspace.get(key)

        workspace.set(key, value)

        # Debug output if enabled
        if debug_enabled and self.processor.console:
            self._debug_value(key, value, old_value, "updated")

        return workspace

    def _debug_value(self, key, value, old_value, action):
        """
        Helper method to print debug information about a workspace value

        Args:
            key: Key being updated
            value: New value
            old_value: Previous value
            action: Action being performed (e.g., "updated", "added")
        """
        value_type = type(value).__name__
        value_info = f"{value_type}"
        if hasattr(value, "__len__") and not isinstance(value, (str, bytes)):
            value_info = f"{value_type}({len(value)})"

        if old_value is None:
            self.processor.console.print(
                f"[grey50]DEBUG: {action.capitalize()} value for '{key}': {value_info}[/grey50]"
            )
        else:
            old_type = type(old_value).__name__
            old_info = f"{old_type}"
            if hasattr(old_value, "__len__") and not isinstance(
                old_value, (str, bytes)
            ):
                old_info = f"{old_type}({len(old_value)})"

            self.processor.console.print(
                f"[grey50]DEBUG: {action.capitalize()} '{key}' from {old_info} to {value_info}[/grey50]"
            )

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        menu = OrderedDict()
        menu["configure"] = "Analyze and configure residue classifications"
        # Import is structure-independent (registers loose .frcmod/.lib files);
        # place it ahead of Parameterize so it reads as an available entry
        # point before the structure-gated workflow steps.
        menu["import_library"] = "Import existing parameters into your library"
        menu["parameterize"] = "Parameterize residues (new or resume)"
        menu["status"] = "View parameterization status"
        menu["help"] = "Display detailed help information"
        return menu

    def get_enhanced_menu_options(self, workspace):
        """
        Get menu options with enhanced status information.

        Args:
            workspace: Current workspace

        Returns:
            List of MenuOption objects with status
        """
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus

        options = []

        # Check workspace state
        non_standard_residues = workspace.get("non_standard_residues", [])
        has_analysis = len(non_standard_residues) > 0
        pending_parameterizations = workspace.get("pending_parameterizations", {})
        has_pending = len(pending_parameterizations) > 0

        # Option 1: Analyze and configure - needs a loaded structure
        if has_analysis:
            ff_status, ff_dep = OptionStatus.COMPLETED, ""
        elif self.can_process(workspace):
            ff_status, ff_dep = OptionStatus.AVAILABLE, ""
        else:
            ff_status = OptionStatus.BLOCKED
            ff_dep = self.availability_note(workspace) or "Load a structure first"
        options.append(MenuOption(
            key="1",
            description="Analyze and configure residue classifications",
            status=ff_status,
            dependency_text=ff_dep,
        ))

        # Option 2: Import existing parameters - always available (structure-
        # independent; registers loose .frcmod/.lib files into the user library).
        # Placed ahead of Parameterize as an available entry point.
        options.append(MenuOption(
            key="2",
            description="Import existing parameters into your library",
            status=OptionStatus.AVAILABLE
        ))

        # Option 3: Parameterize - requires analysis OR pending; ● once done (nothing pending)
        if workspace.get("parameterized_residues") and not has_pending:
            status = OptionStatus.COMPLETED
            dep_text = ""
        elif has_analysis or has_pending:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to analyze residues first] ○"

        options.append(MenuOption(
            key="3",
            description="Parameterize residues (new or resume)",
            status=status,
            dependency_text=dep_text
        ))

        # Option 4: View status - requires analysis OR pending parameterizations
        if has_analysis or has_pending:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to analyze residues first] ○"

        options.append(MenuOption(
            key="4",
            description="View parameterization status",
            status=status,
            dependency_text=dep_text
        ))

        # Option 5: Help - always available
        options.append(MenuOption(
            key="5",
            description="Display detailed help information",
            status=OptionStatus.AVAILABLE
        ))

        return options

    def get_menu_suggestion(self, workspace):
        """
        Get a suggestion for the next recommended action.

        Args:
            workspace: Current workspace

        Returns:
            Suggestion text or None
        """
        non_standard_residues = workspace.get("non_standard_residues", [])
        has_analysis = len(non_standard_residues) > 0
        pending_parameterizations = workspace.get("pending_parameterizations", {})
        pending_count = len(pending_parameterizations)

        # Check for pending work to resume (highest priority)
        # Option numbers below must track get_menu_options: 1 analyze,
        # 2 import, 3 parameterize, 4 status, 5 help. They were left at the
        # pre-"import" numbering and pointed one option short of each action.
        if pending_count > 0:
            return f"⚠ Found {pending_count} pending parameterization{'s' if pending_count != 1 else ''} to resume. Use option 3 to continue, or view status with option 4"

        if not has_analysis:
            if not self.can_process(workspace):
                return f"{self.availability_note(workspace) or 'A structure is required'}. Load one via the Structure Loader."
            return "Start by analyzing residue classifications (option 1) to identify non-standard residues"
        else:
            # Count PARAMETERIZATION UNITS, not residues. Every member of a
            # metal-site unit is stamped "metal_site", the coordinating ligand
            # included, so counting residues reported 4UHX's three metal sites
            # as four (its MoCo site contributes both MOS and its MTE ligand).
            def _units(category):
                return len({self._unit_key(res) for res in non_standard_residues
                            if getattr(res, 'category', None) == category})

            metal_sites = _units("metal_site")
            modified_aas = _units("modified_amino_acid")
            small_molecules = _units("small_molecule")

            total = metal_sites + modified_aas + small_molecules

            if total == 0:
                return "✓ No non-standard residues requiring parameterization. View help (option 5) or press [m] to return to the main menu"

            # Build breakdown
            parts = []
            if metal_sites > 0:
                parts.append(f"{metal_sites} metal site{'s' if metal_sites != 1 else ''}")
            if modified_aas > 0:
                parts.append(f"{modified_aas} modified amino acid{'s' if modified_aas != 1 else ''}")
            if small_molecules > 0:
                parts.append(f"{small_molecules} small molecule{'s' if small_molecules != 1 else ''}")

            breakdown = ", ".join(parts)
            return f"Found {breakdown}. Parameterize with option 3, view status with option 4, or get help with option 5"

    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection using command pattern."""
        if option == "configure":
            return self.configure_classifications_submenu()
        elif option == "parameterize":
            return self.parameterize_with_auto_resume()
        elif option == "status":
            self.display_parameterization_status()
            return True
        elif option == "help":
            command = DisplayHelpCommand(self.processor)
            return command.execute_with_error_handling()
        elif option == "import_library":
            return self._import_into_user_library()
        return False

    def _import_into_user_library(self) -> bool:
        """Standalone import of previously-developed parameters into the library."""
        try:
            from proprep.forcefield_prep.library_promotion import run_import_wizard
            result = run_import_wizard(self.console, self.processor)
        except Exception as e:  # noqa: BLE001 — keep the menu alive on any error
            logger.debug("Import wizard error: %s", e)
            self.console.print(f"[red]Could not import parameters: {e}[/red]")
            return True

        if result:
            self._offer_transformer_for_import(result)
        return True

    @staticmethod
    def _import_forcefield_seed(result: dict) -> dict:
        """Turn a promotion result into a transformer force-field link.

        ``path`` is relative to ``specialized_residues``, which is how a
        transformer names a deposited library; ``state_dir`` ends
        ``.../<redox>/<spin>``.
        """
        seed = {}
        lib_path = result.get("library_path")
        if lib_path:
            parts = Path(lib_path).parts
            if "specialized_residues" in parts:
                idx = parts.index("specialized_residues")
                rel = "/".join(parts[idx + 1:])
                if rel:
                    seed["path"] = rel
        state_dir = result.get("state_dir")
        if state_dir:
            state_parts = Path(state_dir).parts
            if len(state_parts) >= 2:
                seed["redox_state"] = state_parts[-2]
                seed["spin_state"] = state_parts[-1]
        # For a small molecule or modified AA the entry name IS the residue the
        # library parameterizes, so a detected site holding that residue can be
        # identified. For a METAL SITE it is a site identifier ("4hux_fe2s2"),
        # naming no residue -- every detected metal site stays a plausible
        # target there, so no residue name is recorded.
        residue_name = result.get("residue_name")
        if not residue_name and seed.get("path"):
            family, _, entry = seed["path"].partition("/")
            if family in ("small_molecules", "modified_aa") and entry:
                residue_name = entry.split("/")[0]
        if residue_name:
            seed["residue_name"] = residue_name
        return seed

    def _offer_transformer_for_import(self, result: dict) -> None:
        """Offer to build the transformer the imported parameters need.

        Deposited parameters are inert on their own. A transformer is what
        BINDS them to a site: the Topology Generator resolves parameters
        through ``transformer.FORCEFIELD_PATH``, so without one the library is
        unreachable no matter how it is named.

        Renaming is a separate job the same object happens to do. MCPB output
        needs both (CYS -> CM1, plus the library); a cofactor already named as
        the library names it needs only the binding, and saves as a
        pass-through with no edits. Describing the transformer as "the thing
        that renames" made that second case look like it needed nothing, when
        it needs a transformer just as much.

        The editor works on a DETECTED redox site, and the import wizard runs
        without a structure loaded, so this can only be offered when sites are
        present. When they are not, say what is needed rather than launching
        into a failure.
        """
        seed = self._import_forcefield_seed(result)

        workspace = self.processor._get_workspace() if self.processor else None
        redox_sites = (workspace.get("detected_redox_sites") if workspace else None) or []

        self.console.print(
            "\n[bold]Imported parameters need a transformer[/bold]")
        self.console.print(
            "[grey50]A transformer is what binds a deposited library to a site: "
            "the Topology Generator finds parameters through it. If the residues "
            "also need renaming to the library's names it does that too, and if "
            "they are already named correctly it is saved as a pass-through with "
            "no edits — either way the binding is what makes the library "
            "reachable.[/grey50]")

        # Rank the detected sites by how well they match what was imported
        # and SAY so, but never refuse on it. Whether a structure residue
        # should be treated as this library's residue is user knowledge: the
        # FAD in a structure may well be the FAO these parameters describe, and
        # that rename is exactly what the transformer is for.
        if redox_sites:
            ranked = _rank_sites_for_library(
                redox_sites, seed.get("residue_name"),
                _imported_library_atom_names(result))
            strong = [r for r in ranked if r[2] >= 0.8]
            if strong:
                for site, matched, score in strong:
                    self.console.print(
                        f"[grey50]  {getattr(site, 'site_id', '?')} contains "
                        f"{matched} - {score:.0%} of its atoms match the "
                        f"imported library[/grey50]")
            else:
                self.console.print(
                    "[yellow]  No detected site resembles these parameters."
                    "[/yellow]")
                self.console.print(
                    "[grey50]  If the residue they describe is not a site yet: "
                    "run the Redox Site Detector to define one, then Redox Site "
                    "Preparer -> 'Create transformer'. With nothing to rename, "
                    "type 'save' straight away and accept the pass-through."
                    + (f"\n  Link it to: {seed['path']}" if seed.get("path") else "")
                    + "[/grey50]")

        if not redox_sites:
            self.console.print(
                "[yellow]No detected redox sites in this session, and the "
                "transformer editor works on one.[/yellow]")
            self.console.print(
                "[grey50]  To add it later: load the structure, run the Redox "
                "Site Detector, then Redox Site Preparer -> "
                "'Create transformer (interactive PDB editor)'. With nothing to "
                "rename, type 'save' straight away and accept the pass-through."
                + (f"\n  Link it to: {seed['path']}" if seed.get("path") else "")
                + "[/grey50]")
            return

        if not confirm_with_context(
            self.processor,
            "Create the transformer for these parameters now?",
            default=True,
            module="Force Field Parameterizer",
            description="Create a transformer for imported parameters",
        ):
            self.console.print(
                "[grey50]Skipped. Redox Site Preparer → 'Create transformer' "
                "when you are ready"
                + (f"; link it to {seed['path']}" if seed.get("path") else "")
                + ".[/grey50]")
            return

        try:
            module = self.processor.get_module_instance("Redox Site Preparer")
            if module is None:
                raise RuntimeError("Redox Site Preparer module unavailable")
            module.create_custom_transformer(forcefield_default=seed)
        except Exception as e:  # noqa: BLE001 — the import already succeeded
            logger.debug("Transformer creation after import failed: %s", e)
            self.console.print(
                f"[yellow]Could not open the transformer creator ({e}). The "
                f"parameters are imported; create the transformer from the "
                f"Redox Site Preparer menu.[/yellow]")

    def configure_classifications_submenu(self) -> bool:
        """Sub-menu for analyzing and configuring residue classifications."""
        while True:
            self.console.print("\n[bold cyan]Analyze and Configure Residue Classifications[/bold cyan]")
            self.console.print("1. Analyze non-standard residues in structure", highlight=False)
            self.console.print("2. Change classification of a specific residue", highlight=False)
            self.console.print("3. Configure classification settings (auto-classify rules)", highlight=False)
            self.console.print("4. Update modified amino acid mapping table", highlight=False)
            self.console.print("5. Return to Force Field Parameterizer menu", highlight=False)

            choice = prompt_with_context(
                self.processor,
                "\nSelect option",
                choices=["1", "2", "3", "4", "5"],
                default="1",
                module="Force Field Parameterizer",
                description="Configure residue classifications",
                options_map={
                    "1": "Analyze non-standard residues",
                    "2": "Change classification",
                    "3": "Configure settings",
                    "4": "Update mapping table",
                    "5": "Return to menu"
                }
            )

            if choice == "1":
                # Analyze non-standard residues
                command = AnalyzeNonstandardResiduesCommand(self.processor)
                command.execute_with_error_handling()
            elif choice == "2":
                # Change classification
                command = ReclassifyResidueCommand(self.processor)
                command.execute_with_error_handling()
            elif choice == "3":
                # Configure settings
                command = ConfigureClassificationSettingsCommand(self.processor)
                command.execute_with_error_handling()
            elif choice == "4":
                # Update mapping table
                command = ManageMappingsCommand(self.processor)
                command.execute_with_error_handling()
            elif choice == "5":
                # Return to main menu
                return True

        return True

    def parameterize_with_auto_resume(self) -> bool:
        """Parameterize residues with automatic detection of pending workflows."""
        # Check for pending parameterizations
        pending_parameterizations = self.get_from_workspace("pending_parameterizations", {})

        if pending_parameterizations:
            self.console.print(f"\n[yellow]Found {len(pending_parameterizations)} pending parameterization(s)[/yellow]")

            # Show pending workflows
            self.console.print("\n[bold]Pending Parameterizations:[/bold]")
            table = Table()
            table.add_column("No.", style="cyan")
            table.add_column("Residue", style="magenta")
            table.add_column("Type", style="yellow")
            table.add_column("Status", style="green")

            pending_list = list(pending_parameterizations.items())
            for i, (residue_name, data) in enumerate(pending_list, 1):
                param_type = data.get("type", "unknown")
                missing_files = data.get("missing_files", [])
                status = f"Needs {len(missing_files)} calculation(s)"

                table.add_row(str(i), residue_name, param_type.replace("_", " ").title(), status)

            self.console.print(table)

            # Ask user what to do
            self.console.print("\n[bold]Options:[/bold]")
            self.console.print("  [cyan]1[/cyan] Resume a pending workflow")
            self.console.print("  [cyan]2[/cyan] Start new parameterization")
            self.console.print("  [cyan]3[/cyan] Return to main menu")

            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3"],
                default="1",
                module="Force Field Parameterizer",
                description="Resume or start parameterization",
                options_map={
                    "1": "Resume pending workflow",
                    "2": "Start new parameterization",
                    "3": "Return to main menu"
                }
            )

            if choice == "1":
                # Resume workflow
                return self._resume_selected_workflow(pending_list)
            elif choice == "2":
                # Start new parameterization
                return self.parameterize_residue()
            else:
                # Return to main menu
                return True
        else:
            # No pending workflows, go straight to parameterization
            return self.parameterize_residue()

    def _resume_selected_workflow(self, pending_list: list) -> bool:
        """Helper to resume a selected workflow from the pending list."""
        # Get user selection
        resume_options_map = {
            str(i + 1): f"{name} ({data.get('type', 'unknown')})"
            for i, (name, data) in enumerate(pending_list)
        }
        resume_options_map["q"] = "Cancel"
        choice = prompt_with_context(
            self.processor,
            "Select parameterization to resume (or 'q' to cancel)",
            choices=[str(i) for i in range(1, len(pending_list) + 1)] + ["q"],
            default="1",
            module="Force Field Parameterizer",
            description="Select pending parameterization to resume",
            options_map=resume_options_map,
        )

        if choice == "q":
            self.console.print("[yellow]Resume cancelled[/yellow]")
            return True

        # Resume selected workflow
        selected_idx = int(choice) - 1
        residue_name, data = pending_list[selected_idx]

        # Halo every instance of this residue type so the user sees what
        # the resumed parameterization will (eventually) cover. The
        # pending list is keyed by residue name only — no specific
        # (chain, resid) — so the halo covers all matching residues.
        self._halo_residues_by_name(residue_name)

        self.console.print(f"\n[cyan]Resuming parameterization for {residue_name}...[/cyan]")

        # Route to appropriate resume function based on type
        param_type = data.get("type", "unknown")
        if param_type == "modified_amino_acid":
            self._resume_modified_amino_acid_workflow(residue_name, data)
        elif param_type == "small_molecule":
            self._resume_small_molecule_workflow(residue_name, data)
        elif param_type == "metal_site":
            self._resume_metal_site_workflow(residue_name, data)
        else:
            self.console.print(f"[red]Unknown parameterization type: {param_type}[/red]")

        return True

    def reclassify_residue(self):
        """Interactive function to change the classification of a non-standard residue"""
        if not self.non_standard_residues:
            # Try to get from workspace
            self.non_standard_residues = self.get_from_workspace(
                "non_standard_residues", []
            )

            if not self.non_standard_residues:
                self.console.print(
                    "[yellow]No non-standard residues analyzed yet. Running analysis...[/yellow]"
                )
                self.analyze_nonstandard_residues()

                if not self.non_standard_residues:
                    self.console.print(
                        "[red]No non-standard residues found to reclassify[/red]"
                    )
                    return

        # Display all non-standard residues
        self.console.print("\n[bold]Select a residue to reclassify:[/bold]")
        table = Table(title="Available Non-Standard Residues")
        table.add_column("No.", style="cyan", justify="right")
        table.add_column("Name", style="magenta")
        table.add_column("Location", style="green")
        table.add_column("Current Category", style="yellow")
        table.add_column("Parent/Notes", style="blue")

        # Track residues we've already displayed to avoid duplicates
        displayed_residues = set()

        # Add rows for each residue
        displayed_list = []
        for res in self.non_standard_residues:
            # Create a unique identifier for this residue
            res_key = (res.name, res.chain_id, res.resid)

            # Skip if we've already displayed this residue
            if res_key in displayed_residues:
                continue

            displayed_residues.add(res_key)
            displayed_list.append(res)

            category_display = res.category.replace("_", " ").capitalize()
            location = res.get_location_str() if hasattr(res, 'get_location_str') else f"{res.chain_id}:{res.resid}"

            # Add parent info if available
            notes = ""
            if res.category == "modified_amino_acid" and hasattr(res, 'parent_residue') and res.parent_residue:
                notes = f"Parent: {res.parent_residue}"
            elif res.category == "small_molecule":
                atom_count = getattr(res, 'atom_count', 'unknown')
                notes = f"Small molecule ({atom_count} atoms)"
            elif res.category == "metal_site":
                notes = "Metal coordination site"
            elif hasattr(res, 'notes') and res.notes:
                notes = res.notes

            table.add_row(
                str(len(displayed_list)), res.name, location, category_display, notes
            )

        self.console.print(table)

        # Get selection
        reclass_options_map = {
            str(i + 1): f"{r.name} ({r.category})"
            for i, r in enumerate(displayed_list)
        }
        reclass_options_map["q"] = "Cancel"
        choice = prompt_with_context(
            self.processor,
            "Enter residue number to reclassify (or 'q' to cancel)",
            choices=[str(i) for i in range(1, len(displayed_list) + 1)] + ["q"],
            default="q",
            module="Force Field Parameterizer",
            description="Select residue to reclassify",
            options_map=reclass_options_map,
        )

        if choice == "q":
            self.console.print("[yellow]Reclassification cancelled[/yellow]")
            return

        # Get selected residue
        selected_idx = int(choice) - 1
        selected_residue = displayed_list[selected_idx]

        # Halo the picked residue so the user has visual context for the
        # category they're about to choose. Sits on top of the per-
        # classification halos already drawn after the analyze step.
        self._halo_residue_instance(selected_residue)

        # Display current classification
        self.console.print(
            f"\nSelected: {selected_residue.name} ({location})"
        )
        self.console.print(
            f"Current classification: {selected_residue.category.replace('_', ' ').capitalize()}"
        )
        if (
            selected_residue.category == "modified_amino_acid"
            and hasattr(selected_residue, 'parent_residue') and selected_residue.parent_residue
        ):
            self.console.print(
                f"Current parent residue: {selected_residue.parent_residue}"
            )

        # Show available categories with descriptions
        self.console.print("\n[cyan]Available classifications:[/cyan]")
        
        categories = [
            "modified_amino_acid",
            "small_molecule", 
            "metal_site",
            "unknown"
        ]

        category_descriptions = {
            "modified_amino_acid": "Modified amino acid (derived from standard amino acid)",
            "small_molecule": "Small molecule (ligand, cofactor, or organic compound)", 
            "metal_site": "Metal site (metal ion or metal-containing compound)",
            "unknown": "Unknown (requires manual classification)"
        }

        for i, category in enumerate(categories, 1):
            category_display = category.replace("_", " ").capitalize()
            description = category_descriptions[category]
            current_marker = " ← Current" if selected_residue.category == category else ""
            self.console.print(f"  {i}. {category_display} - {description}{current_marker}")

        # Get new classification
        cat_options_map = {str(i + 1): c.replace("_", " ").capitalize() for i, c in enumerate(categories)}
        cat_options_map["q"] = "Cancel"
        cat_choice = prompt_with_context(
            self.processor,
            "Select new classification",
            choices=[str(i) for i in range(1, len(categories) + 1)] + ["q"],
            default="q",
            module="Force Field Parameterizer",
            description="New classification for residue",
            options_map=cat_options_map,
        )

        if cat_choice == "q":
            self.console.print("[yellow]Reclassification cancelled[/yellow]")
            return

        # Update residue classification
        new_category = categories[int(cat_choice) - 1]
        old_category = selected_residue.category
        selected_residue.category = new_category

        self.console.print(
            f"[green]Successfully reclassified {selected_residue.name} from '{old_category.replace('_', ' ').capitalize()}' to '{new_category.replace('_', ' ').capitalize()}'[/green]"
        )

        # Handle specific classification requirements
        if new_category == "modified_amino_acid":
            # Check if we have this in our mapping table (use uppercase for comparison)
            if selected_residue.name.upper() in self.modified_aa_map:
                mapped_parent = self.modified_aa_map[selected_residue.name.upper()]
                self.console.print(
                    f"[green]Found in mapping table: {selected_residue.name} is derived from {mapped_parent}[/green]"
                )

                use_mapped = confirm_with_context(
                    self.processor,
                    f"Use mapped parent residue ({mapped_parent})?",
                    default=True,
                    module="Force Field Parameterizer",
                    description=f"Use mapped parent residue {mapped_parent}",
                )

                if use_mapped:
                    selected_residue.parent_residue = mapped_parent
                    selected_residue.notes = f"Modified {mapped_parent} - requires parameterization"
                    self.console.print(
                        f"[green]Set parent residue to {mapped_parent}[/green]"
                    )
                else:
                    # Manual parent selection
                    standard_aas = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
                                "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", 
                                "THR", "TRP", "TYR", "VAL"]
                    
                    self.console.print(f"\n[cyan]Select parent amino acid for {selected_residue.name}:[/cyan]")
                    for i, aa in enumerate(standard_aas, 1):
                        self.console.print(f"  {i:2d}. {aa}")
                    
                    parent_aa_map = {str(i + 1): aa for i, aa in enumerate(standard_aas)}
                    parent_choice = prompt_with_context(
                        self.processor,
                        "Choose parent amino acid",
                        choices=[str(i) for i in range(1, len(standard_aas) + 1)],
                        module="Force Field Parameterizer",
                        description="Select parent amino acid for modified AA",
                        options_map=parent_aa_map,
                    )
                    
                    parent_aa = standard_aas[int(parent_choice) - 1]
                    selected_residue.parent_residue = parent_aa
                    selected_residue.notes = f"Modified {parent_aa} - requires parameterization"
                    self.console.print(f"[green]Set parent residue to {parent_aa}[/green]")
            else:
                # Not in mapping table, ask user to select parent
                standard_aas = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
                            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", 
                            "THR", "TRP", "TYR", "VAL"]
                
                self.console.print(f"\n[cyan]Select parent amino acid for {selected_residue.name}:[/cyan]")
                for i, aa in enumerate(standard_aas, 1):
                    self.console.print(f"  {i:2d}. {aa}")
                
                parent_aa_map2 = {str(i + 1): aa for i, aa in enumerate(standard_aas)}
                parent_choice = prompt_with_context(
                    self.processor,
                    "Choose parent amino acid",
                    choices=[str(i) for i in range(1, len(standard_aas) + 1)],
                    module="Force Field Parameterizer",
                    description="Select parent amino acid for modified AA",
                    options_map=parent_aa_map2,
                )
                
                parent_aa = standard_aas[int(parent_choice) - 1]
                selected_residue.parent_residue = parent_aa
                selected_residue.notes = f"Modified {parent_aa} - requires parameterization"
                self.console.print(f"[green]Set parent residue to {parent_aa}[/green]")
                
                # Offer to add to mapping table
                add_to_map = confirm_with_context(
                    self.processor,
                    f"Add {selected_residue.name} → {parent_aa} to mapping table for future use?",
                    default=True,
                    module="Force Field Parameterizer",
                    description=f"Add {selected_residue.name} to modified AA mapping table",
                )
                if add_to_map:
                    self.modified_aa_map[selected_residue.name.upper()] = parent_aa
                    self._save_modified_aa_map()
                    self.console.print(f"[green]✓ Added {selected_residue.name} → {parent_aa} to mapping table[/green]")

        elif new_category == "small_molecule":
            # Count atoms if possible - look for atom_count attribute or count from structure
            atom_count = getattr(selected_residue, 'atom_count', 'unknown')
            if atom_count == 'unknown':
                # Try to get from the residue object if available
                try:
                    if hasattr(selected_residue, 'residue') and hasattr(selected_residue.residue, 'get_atoms'):
                        atom_count = len(list(selected_residue.residue.get_atoms()))
                    elif hasattr(selected_residue, 'biopython_residue') and hasattr(selected_residue.biopython_residue, 'get_atoms'):
                        atom_count = len(list(selected_residue.biopython_residue.get_atoms()))
                except:
                    atom_count = 'unknown'
            
            selected_residue.notes = f"Small molecule ({atom_count} atoms) - can be parameterized with GAFF2/RESP"
            self.console.print(f"[green]Classified as small molecule with {atom_count} atoms[/green]")
            self.console.print("[cyan]This residue can now be parameterized using the Small Molecule Parameterizer[/cyan]")

        elif new_category == "metal_site":
            selected_residue.notes = "Metal site - requires metal parameterization"
            self.console.print(f"[yellow]Classified as metal site - metal parameterization will be needed[/yellow]")

        elif new_category == "unknown":
            selected_residue.notes = "Classification uncertain - manual review needed"
            self.console.print(f"[yellow]Classified as unknown - manual classification may be needed[/yellow]")

        # Store user classification for future reference
        if not hasattr(self, 'user_residue_classifications'):
            self.user_residue_classifications = {}

        self.user_residue_classifications[selected_residue.name.upper()] = {
            "category": new_category,
            "notes": selected_residue.notes,
            "user_defined": True
        }

        # Show updated information
        self.console.print(f"\n[bold]Updated classification for {selected_residue.name}:[/bold]")
        self.console.print(f"  Category: {new_category.replace('_', ' ').capitalize()}")
        if hasattr(selected_residue, 'notes') and selected_residue.notes:
            self.console.print(f"  Notes: {selected_residue.notes}")
        if new_category == "modified_amino_acid" and hasattr(selected_residue, 'parent_residue') and selected_residue.parent_residue:
            self.console.print(f"  Parent residue: {selected_residue.parent_residue}")

        # Update workspace
        self.update_workspace("non_standard_residues", self.non_standard_residues)
        self.update_workspace("user_residue_classifications", self.user_residue_classifications)

        self.console.print("[green]✓ Classification updated and saved to workspace[/green]")

        return True

    def update_aa_mapping(self):
        """Function to add or update an entry in the modified amino acid mapping table"""
        self.console.print("\n[bold]Update Modified Amino Acid Mapping Table[/bold]")

        # Display current mapping table
        self._display_aa_mapping_table()

        # Get the modified amino acid code
        mod_aa = prompt_with_context(
            self.processor,
            "Enter the 3-letter code for the modified amino acid",
            default="",
            module="Force Field Parameterizer",
            description="Modified amino acid 3-letter code",
        ).upper()

        if not mod_aa:
            self.console.print("[yellow]Operation cancelled[/yellow]")
            return

        # Check if it's already in the mapping
        if mod_aa in self.modified_aa_map:
            current_parent = self.modified_aa_map[mod_aa]
            self.console.print(
                f"[yellow]Note: {mod_aa} is already mapped to {current_parent}[/yellow]"
            )

        # List of standard amino acids for selection
        standard_aa_list = list(self.standard_aa)
        standard_aa_list.sort()

        self.console.print("\n[bold]Select parent residue:[/bold]")
        for i, aa in enumerate(standard_aa_list, 1):
            self.console.print(f"  {i}. {aa}")

        parent_choice_map = {str(i + 1): aa for i, aa in enumerate(standard_aa_list)}
        parent_choice_map["c"] = "Custom 3-letter code"
        aa_choice = prompt_with_context(
            self.processor,
            "Select parent amino acid (or enter 'c' for custom 3-letter code)",
            choices=[str(i) for i in range(1, len(standard_aa_list) + 1)] + ["c"],
            default="c",
            module="Force Field Parameterizer",
            description="Select parent amino acid",
            options_map=parent_choice_map,
        )

        if aa_choice == "c":
            # Custom entry
            parent_aa = prompt_with_context(
                self.processor,
                "Enter 3-letter code for parent amino acid",
                default="UNK",
                module="Force Field Parameterizer",
                description="Custom parent amino acid 3-letter code",
            ).upper()
        else:
            # Selection from list
            parent_aa = standard_aa_list[int(aa_choice) - 1]

        # Update the mapping
        self.modified_aa_map[mod_aa] = parent_aa
        self._save_modified_aa_map()

        self.console.print(
            f"[green]Successfully added/updated mapping: {mod_aa} -> {parent_aa}[/green]"
        )

        # Display the updated table
        self._display_aa_mapping_table()

    def _display_aa_mapping_table(self):
        """Display the current modified amino acid mapping table"""
        if not self.modified_aa_map:
            self.console.print(
                "[yellow]The modified amino acid mapping table is empty[/yellow]"
            )
            return

        table = Table(title="Modified Amino Acid Mapping")
        table.add_column("Modified AA", style="cyan")
        table.add_column("Parent AA", style="green")

        # Sort by modified AA code for consistent display
        sorted_keys = sorted(self.modified_aa_map.keys())

        for mod_aa in sorted_keys:
            parent_aa = self.modified_aa_map[mod_aa]
            table.add_row(mod_aa, parent_aa)

        self.console.print(table)
        self.console.print(f"Total entries: {len(self.modified_aa_map)}")

    def _save_modified_aa_map(self):
        """Save the modified amino acid mapping table to a config file"""
        try:
            # Create config directory if it doesn't exist
            config_dir = os.path.join(os.path.expanduser("~"), ".proprep")
            os.makedirs(config_dir, exist_ok=True)

            # Save the mapping to a file
            config_file = os.path.join(config_dir, "modified_aa_map.json")
            with open(config_file, "w") as f:
                json.dump(self.modified_aa_map, f, indent=2)

            self.console.print(f"[green]Saved mapping table to {config_file}[/green]")
        except Exception as e:
            self.console.print(f"[red]Error saving mapping table: {str(e)}[/red]")

    def _load_modified_aa_map(self):
        """Load the modified amino acid mapping table from a config file"""
        try:
            config_file = os.path.join(
                os.path.expanduser("~"), ".proprep", "modified_aa_map.json"
            )
            if os.path.exists(config_file):
                with open(config_file, "r") as f:
                    user_map = json.load(f)

                # Update the built-in mapping with user entries
                # Note: This preserves the built-in mappings while adding user ones
                self.modified_aa_map.update(user_map)
                # logger.info(
                #     f"Loaded user-defined mapping table with {len(user_map)} entries"
                # )
        except Exception as e:
            self.console.print(
                f"[yellow]Error loading mapping table: {str(e)}[/yellow]"
            )

    def _prompt_parameterization(self, residue):
        """Helper method to prompt for parameterization of a residue"""
        if residue.category == "modified_amino_acid" and residue.parent_residue:
            parameterize_now = confirm_with_context(
                self.processor,
                f"Parameterize {residue.name} now?",
                default=True,
                module="Force Field Parameterizer",
                description=f"Parameterize {residue.name} now",
            )

            if parameterize_now:
                self.parameterize_modified_amino_acid(residue.name, [residue])

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements - needs at least one structure loaded"""
        return [
            "rcsb_pdb_file | local_pdb_file | alphafold_pdb_file | alphafill_pdb_file | alphafold_homolog_pdb_file"
        ]

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return [
            "parameterized_residues",
            "non_standard_residues",
            "user_residue_classifications",
            "pending_parameterizations",
            "global_atom_registry_data",
        ]

    def can_process(self, workspace: Dict[str, Any]) -> bool:
        """Check if the module can process the current workspace"""
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)
        status = selector.get_structure_status()
        return status.get("has_any", False)

    def display_help(self):
        """
        Display detailed help information for the Forcefield Parameterization module.
        """
        help_text = """
[bold cyan]Forcefield Parameterization Module[/bold cyan]
[grey50]═══════════════════════════════════[/grey50]

This module identifies and parameterizes non-standard residues for molecular dynamics
simulations. The workflow is:

[bold]1.[/bold] Use [cyan]"Analyze non-standard residues"[/cyan] to scan your structure
[bold]2.[/bold] Review automatic classifications:
   • Modified amino acids (residues derived from standard amino acids)
   • Small molecules (ligands, cofactors, organic compounds)
   • Metal sites (metal ions and coordination environment)
   • Unknown components (for manual classification)
[bold]3.[/bold] Reclassify residues if needed using [cyan]"Change classification"[/cyan]
[bold]4.[/bold] Use [cyan]"Parameterize a specific residue"[/cyan] to generate parameters
[bold]5.[/bold] Manage modified amino acid mappings via [cyan]"Add/update mapping table"[/cyan]

[bold cyan]━━━ Small Molecule Parameterization ━━━[/bold cyan]

[bold]Workflow (8 steps):[/bold]
  [cyan]1.[/cyan] Extract coordinates from PDB structure
  [cyan]2.[/cyan] Add hydrogen atoms (optional, via reduce)
  [cyan]3.[/cyan] Generate Gaussian input file (two-step workflow)
  [cyan]4.[/cyan] Gaussian calculation checkpoint (interactive or deferred)
  [cyan]5.[/cyan] Process Gaussian output → RESP charges (antechamber)
  [cyan]6.[/cyan] Generate force field parameters (parmchk2)
  [cyan]7.[/cyan] Create AMBER topology files (tLEaP)
  [cyan]8.[/cyan] Parameter refinement (optional)

[bold]Gaussian Two-Step Workflow:[/bold]
  • [yellow]Step 1:[/yellow] B3LYP/6-31+G(d) optimization + frequency
    - Finds minimum energy geometry
    - IOp(7/33=1) saves Hessian for Seminario refinement
  • [yellow]Step 2:[/yellow] HF/6-31G(d) ESP calculation (via --link1--)
    - Merz-Kollman ESP for RESP fitting
    - HF level maintains AMBER charge compatibility

[bold]Parameter Refinement Options:[/bold]
  • [cyan]Seminario method[/cyan] (bonds/angles): Derives force constants from QM Hessian
    - Uses existing frequency calculation data
    - Fast (seconds), no additional QM needed
  • [cyan]PES scan[/cyan] (dihedrals): Systematic torsional scans
    - Default: 24 points at 15° increments (360° scan)
    - Uses relaxed scans (opt=modredundant)
    - Fits parameters with AMBER paramfit
  • [cyan]CREST[/cyan] (dihedrals): Conformer sampling with paramfit
    - GFN2-xTB metadynamics for conformer generation
    - Samples all rotatable bonds simultaneously
    - Faster than PES for molecules with many dihedrals

[bold]Output Files:[/bold]
  • molecule.mol2 - Coordinates, atom types, RESP charges
  • molecule.frcmod - Force field parameters
  • molecule.lib - AMBER library file for tLEaP
  • molecule_seminario.frcmod - Seminario-refined parameters (if used)

[bold cyan]━━━ Modified Amino Acid Parameterization ━━━[/bold cyan]

[bold]Workflow (10 steps):[/bold]
  [cyan]1.[/cyan] Generate ACE-XXX-NME capped tripeptide (tLEaP)
  [cyan]2.[/cyan] Set backbone conformations with cpptraj
      - α-helix: φ=-60°, ψ=-45°
      - β-sheet: φ=-135°, ψ=135°
  [cyan]3.[/cyan] Create Gaussian input files
  [cyan]4.[/cyan] Run Gaussian calculations (external)
  [cyan]5.[/cyan] Analyze calculations and extract structures
  [cyan]6.[/cyan] Generate ESP from optimized structures
  [cyan]7.[/cyan] Create AC file from lowest-energy structure
  [cyan]8.[/cyan] Run residuegen for multi-conformation RESP fitting
  [cyan]9.[/cyan] Generate bonded parameters (parmchk2)
  [cyan]10.[/cyan] Create AMBER library file

[bold]Conformational Sampling Options:[/bold]
  • [cyan]Simple workflow[/cyan]: 2 conformations (α-helix + β-sheet)
  • [cyan]PES scan workflow[/cyan]: Additional backbone dihedral scans
    - Broader conformational coverage for charge fitting
    - Recommended for flexible modifications

[bold]Protonation State Changes:[/bold]
  When the modified residue has different protonation than parent:
  1. Run PES scans in [yellow]neutral[/yellow] state (avoid vacuum artifacts)
  2. Extract structures, modify to target protonation
  3. Re-optimize in [yellow]charged[/yellow] state with constraints:
     - Dihedral freeze: Lock backbone φ/ψ angles
     - Atom freeze: Freeze all except modified group

[bold]Gaussian Level of Theory:[/bold]
  • Uses HF/6-31G* throughout (single-step, not two-step)
  • Maintains compatibility with standard AMBER backbone charges
  • Multi-conformation ESP → charges valid across secondary structures

[bold]Output Files:[/bold]
  • residue.prep - AMBER prep format with RESP charges
  • residue.frcmod - Force field parameters
  • residue.lib - AMBER OFF library file

[bold cyan]━━━ Metal Site Parameterization ━━━[/bold cyan]

  • Step 1: Structure preparation and validation ✅ Available
  • Step 2: Metal center parameter building 🚧 Coming soon
  • Step 3: Force field integration 🚧 Coming soon
  • Step 4: Simulation setup 🚧 Coming soon

[bold]Features:[/bold]
  • Interactive structure validation and fixing
  • Automatic component identification
  • Integration with small molecule parameterizer for ligands
  • Educational interface with command explanations

[grey50]Note: Metal site parameterization requires a repaired structure (no missing
atoms/residues) due to the sensitivity of metal coordination geometry.[/grey50]
        """

        self.console.print(Panel(help_text, title="Forcefield Parameterization Help", expand=False))


    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    def configure_classification_settings(self) -> bool:
        """Display classification strategy and allow atom count configuration."""
        self.console.print(Panel(
            "[bold cyan]Classification Strategy[/bold cyan]\n\n"
            "Shows how ProPrep automatically classifies non-standard residues.\n"
            "Only atom count range for small molecules is configurable.",
            title="Classification Settings",
            expand=False
        ))

        # Display current settings
        self._display_current_settings()

        # Explain the classification process
        self.console.print("\n[bold]═══ Classification Process ═══[/bold]\n")

        min_atoms = self.classification_settings.get("min_small_molecule_atoms", 2)
        max_atoms = self.classification_settings.get("max_small_molecule_atoms", 200)

        process = f"""[cyan]Priority 1: User Manual Classification[/cyan]
   • Residues you've manually classified are always respected
   • Set via: 'Change classification of specific residue'

[cyan]Priority 2: RedoxSite Classification[/cyan]
   • Residues identified in metal/redox sites by RedoxDetector
   • Categories: metal_site, small_molecule, modified_amino_acid

[cyan]Priority 3: Metal Ion Detection[/cyan]
   • Uses comprehensive METALS dictionary from periodic table
   • Source: comprehensive_redox_detector.METALS (always enabled)

[cyan]Priority 4: Modified Amino Acid Detection[/cyan]
   • Checked in order:
     a. Local mapping table (modified_aa_map)
     b. CCD parent residue field (mon_nstd_parent_comp_id)
     c. CCD peptide linking classification
     d. Protein backbone atoms (N, CA, C, O)

[cyan]Priority 5: Small Molecule Detection[/cyan]
   • Atom count range: [green]{min_atoms}-{max_atoms} atoms[/green] (configurable)
   • CCD NON-POLYMER check (always enabled when CCD available)
   • Applied ONLY if not classified above

[cyan]Priority 6: Unknown[/cyan]
   • Residues that don't match any criteria
   • Require manual classification
"""

        self.console.print(process)

        # Ask if user wants to modify atom count settings
        if confirm_with_context(
            self.processor,
            "\nModify small molecule atom count range?",
            default=False,
            module="Force Field Parameterizer",
            description="Modify small molecule atom count classification range",
        ):
            self.classification_settings["min_small_molecule_atoms"] = int(prompt_with_context(
                self.processor,
                "Minimum atoms for small molecule",
                default=str(min_atoms),
                module="Force Field Parameterizer",
                description="Minimum atoms for small-molecule classification",
            ))

            self.classification_settings["max_small_molecule_atoms"] = int(prompt_with_context(
                self.processor,
                "Maximum atoms for small molecule",
                default=str(max_atoms),
                module="Force Field Parameterizer",
                description="Maximum atoms for small-molecule classification",
            ))

            # Save settings
            self._save_classification_settings()

            # Display updated settings
            self.console.print("\n[green]✓ Settings updated![/green]")
            self._display_current_settings()

        return True

    def manage_mappings(self) -> bool:
        """Manage amino acid and residue classification mappings."""
        while True:
            self.console.print(Panel(
                "[bold cyan]Manage Classification Mappings[/bold cyan]\n\n"
                "View and edit mappings for:\n"
                "• Modified amino acids (e.g., MSE → MET)\n"
                "• User-defined residue classifications",
                title="Mapping Management"
            ))
            
            options = {
                "view_aa": "View amino acid mappings",
                "edit_aa": "Edit amino acid mappings", 
                "view_user": "View user residue classifications",
                "edit_user": "Edit user residue classifications",
                "clear_user": "Clear user classifications",
                "back": "Back to main menu"
            }
            
            self.console.print("\n[bold]Available Actions:[/bold]")
            for key, desc in options.items():
                self.console.print(f"  {key}: {desc}")
            
            choice = prompt_with_context(None,
                "Choose action",
                choices=list(options.keys()),
                default="back"
            )
            
            if choice == "view_aa":
                self._display_aa_mappings()
            elif choice == "edit_aa":
                self._edit_aa_mappings()
            elif choice == "view_user":
                self._display_user_classifications()
            elif choice == "edit_user":
                self._edit_user_classifications()
            elif choice == "clear_user":
                self._clear_user_classifications()
            elif choice == "back":
                break
        
        return True

    def _display_current_settings(self):
        """Display current classification settings."""
        table = Table(title="Current Classification Settings")
        table.add_column("Setting", style="cyan")
        table.add_column("Value", style="green")
        table.add_column("Description", style="yellow")

        table.add_row(
            "Small Molecule Atom Range",
            f"{self.classification_settings['min_small_molecule_atoms']}-{self.classification_settings['max_small_molecule_atoms']} atoms",
            "Atom count range for small molecules"
        )

        table.add_row(
            "Metal Detection",
            "✓ (always enabled)",
            "Uses comprehensive METALS dictionary"
        )

        table.add_row(
            "CCD Integration",
            "✓ (always enabled)",
            "Uses Chemical Component Dictionary when available"
        )

        table.add_row(
            "Backbone Detection",
            "✓ (always enabled)",
            "Detects modified AAs via N, CA, C, O atoms"
        )

        self.console.print(table)

    def _display_aa_mappings(self):
        """Display amino acid mappings."""
        if not self.modified_aa_map:
            self.console.print("[yellow]No amino acid mappings found[/yellow]")
            return
        
        table = Table(title="Modified Amino Acid Mappings")
        table.add_column("Modified AA", style="cyan")
        table.add_column("Parent AA", style="green")
        
        for modified, parent in sorted(self.modified_aa_map.items()):
            table.add_row(modified, parent)
        
        self.console.print(table)

    def _edit_aa_mappings(self):
        """Edit amino acid mappings."""
        self.console.print("\n[bold]Edit Amino Acid Mappings[/bold]")
        
        modified_aa = prompt_with_context(
            self.processor,
            "Modified amino acid code (3-letter)",
            module="Force Field Parameterizer",
            description="Modified amino acid 3-letter code to classify",
        ).upper().strip()
        if not modified_aa or len(modified_aa) != 3:
            self.console.print("[red]Please enter a valid 3-letter code[/red]")
            return
        
        standard_aas = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        
        self.console.print(f"\nSelect parent amino acid for {modified_aa}:")
        for i, aa in enumerate(standard_aas, 1):
            self.console.print(f"  {i:2d}. {aa}")
        
        choice = prompt_with_context(None,
            "Choose parent amino acid",
            choices=[str(i) for i in range(1, len(standard_aas) + 1)]
        )
        
        parent_aa = standard_aas[int(choice) - 1]
        
        self.modified_aa_map[modified_aa] = parent_aa
        self._save_modified_aa_map()
        
        self.console.print(f"[green]✓ {modified_aa} → {parent_aa} mapping added[/green]")

    def _display_user_classifications(self):
        """Display user-defined residue classifications."""
        if not self.user_residue_classifications:
            self.console.print("[yellow]No user-defined classifications found[/yellow]")
            return
        
        table = Table(title="User-Defined Residue Classifications")
        table.add_column("Residue", style="cyan")
        table.add_column("Classification", style="green")
        table.add_column("Notes", style="yellow")
        
        for residue, info in self.user_residue_classifications.items():
            table.add_row(
                residue,
                info.get("category", "unknown"),
                info.get("notes", "")
            )
        
        self.console.print(table)

    def _edit_user_classifications(self):
        """Edit user-defined residue classifications."""
        residue = prompt_with_context(
            self.processor,
            "Residue name to classify",
            module="Force Field Parameterizer",
            description="Residue name to classify",
        ).upper().strip()
        
        if not residue:
            return
        
        categories = ["small_molecule", "metal_site", "modified_amino_acid", "unknown", "ignore"]
        
        self.console.print(f"\nClassify {residue} as:")
        for i, cat in enumerate(categories, 1):
            desc = {
                "small_molecule": "Small molecule (will be sent to Small Molecule Parameterizer)",
                "metal_site": "Metal site (metal coordination)",
                "modified_amino_acid": "Modified amino acid (will be sent to AA parameterizer)",
                "unknown": "Unknown (no automatic parameterization)",
                "ignore": "Ignore (skip in analysis)"
            }
            self.console.print(f"  {i}. {cat.replace('_', ' ').title()} - {desc[cat]}")
        
        choice = prompt_with_context(None,
            "Choose classification",
            choices=[str(i) for i in range(1, len(categories) + 1)]
        )
        
        selected_category = categories[int(choice) - 1]
        notes = prompt_with_context(
            self.processor,
            "Notes (optional)",
            default="",
            module="Force Field Parameterizer",
            description="Classification notes for residue",
        )
        
        self.user_residue_classifications[residue] = {
            "category": selected_category,
            "notes": notes
        }
        
        self._save_user_residue_classifications()
        self.console.print(f"[green]✓ {residue} classified as {selected_category}[/green]")

    def _clear_user_classifications(self):
        """Clear all user-defined classifications."""
        if not self.user_residue_classifications:
            self.console.print("[yellow]No user classifications to clear[/yellow]")
            return
        
        confirm = confirm_with_context(
            self.processor,
            f"Clear all {len(self.user_residue_classifications)} user classifications?",
            default=False,
            module="Force Field Parameterizer",
            description="Clear all user residue classifications",
        )
        
        if confirm:
            self.user_residue_classifications.clear()
            self._save_user_residue_classifications()
            self.console.print("[green]✓ All user classifications cleared[/green]")

    def _load_classification_settings(self):
        """Load classification settings from config file."""
        try:
            config_file = os.path.join(
                os.path.expanduser("~"), ".proprep", "classification_settings.json"
            )
            if os.path.exists(config_file):
                with open(config_file, "r") as f:
                    saved_settings = json.load(f)
                    self.classification_settings.update(saved_settings)
        except Exception as e:
            # Use defaults if loading fails
            pass

    def _save_classification_settings(self):
        """Save classification settings to config file."""
        try:
            config_dir = os.path.join(os.path.expanduser("~"), ".proprep")
            os.makedirs(config_dir, exist_ok=True)
            
            config_file = os.path.join(config_dir, "classification_settings.json")
            with open(config_file, "w") as f:
                json.dump(self.classification_settings, f, indent=2)
        except Exception as e:
            self.console.print(f"[yellow]Warning: Could not save settings: {str(e)}[/yellow]")

    def _load_user_residue_classifications(self):
        """Load user residue classifications from config file."""
        try:
            config_file = os.path.join(
                os.path.expanduser("~"), ".proprep", "user_residue_classifications.json"
            )
            if os.path.exists(config_file):
                with open(config_file, "r") as f:
                    self.user_residue_classifications = json.load(f)
        except Exception as e:
            self.user_residue_classifications = {}

    def _save_user_residue_classifications(self):
        """Save user residue classifications to config file."""
        try:
            config_dir = os.path.join(os.path.expanduser("~"), ".proprep")
            os.makedirs(config_dir, exist_ok=True)
            
            config_file = os.path.join(config_dir, "user_residue_classifications.json")
            with open(config_file, "w") as f:
                json.dump(self.user_residue_classifications, f, indent=2)
        except Exception as e:
            self.console.print(f"[yellow]Warning: Could not save classifications: {str(e)}[/yellow]")

    def _is_simple_metal_ion(self, res_code, atom_count):
        """Check if residue is a simple metal ion (residue name = metal name, single atom)."""
        return res_code.upper() in self.metals and atom_count == 1

    def _get_atom_count_from_structure(self, res_code, chain_id, res_id):
        """Get atom count from the actual structure."""
        try:
            structure = self._get_structure_object()

            if structure:
                for model in structure:
                    for chain in model:
                        if chain.id == chain_id:
                            for residue in chain:
                                if residue.id[1] == res_id and residue.get_resname().strip() == res_code:
                                    return len(list(residue.get_atoms()))
            
            # Fallback - conservative default
            return 10
            
        except Exception:
            return 10  # Safe default

    def _is_small_molecule(self, res_code, ccd_data, atom_count):
        """
        Determine if a residue should be classified as a small molecule based on current settings.
        NO HARD-CODED RESIDUE NAMES - only use user classifications and CCD data.
        """
        # Check if small molecule detection is enabled
        if not self.classification_settings["small_molecule_detection"]:
            return False
        
        # Check user-defined classification first (highest priority)
        if res_code.upper() in self.user_residue_classifications:
            user_class = self.user_residue_classifications[res_code.upper()]
            return user_class.get("category") == "small_molecule"
        
        # Check atom count range
        min_atoms = self.classification_settings["min_small_molecule_atoms"]
        max_atoms = self.classification_settings["max_small_molecule_atoms"]
        
        if not (min_atoms <= atom_count <= max_atoms):
            return False
        
        # Rule out water (always)
        if res_code.upper() in ["HOH", "WAT", "H2O"]:
            return False
        
        # Check CCD NON-POLYMER classification if enabled and available
        if self.classification_settings["use_ccd_non_polymer"] and ccd_data:
            try:
                if "type" in ccd_data and ccd_data["type"]:
                    ccd_type = ccd_data["type"].upper()
                    if "NON-POLYMER" in ccd_type:
                        return True
            except:
                pass
        
        # If no CCD data available, default to unknown (will need user classification)
        return False

    def _check_metal_sites_availability(self):
        """
        Check if metal sites have been detected and provide informational message.
        
        This is purely informational - the analysis can proceed without metal sites.
        
        Returns:
            bool: Always returns True (analysis can always proceed)
        """
        
        metal_sites = self.get_from_workspace("metal_sites", [])
        
        if not metal_sites:
            self.console.print("\n[blue]ℹ️  Metal Site Detection Info:[/blue]")
            self.console.print("[blue]No metal sites found in workspace. If your structure contains metal sites,[/blue]")
            self.console.print("[blue]run 'Metal Site Preparation' → 'Detect metal coordination sites' first[/blue]")
            self.console.print("[blue]for more accurate metal site classification in this analysis.[/blue]")
            self.console.print("[blue]For non-metalloprotein structures, this message can be ignored.[/blue]")
        else:
            self.console.print(f"\n[green]✅ Found {len(metal_sites)} metal site(s) in workspace - will be used for classification[/green]")
        
        return True

    # ========================================================================
    # Viewer integration helpers (dispatcher-level)
    # ========================================================================
    #
    # Per-classification hex colours used by ``_halo_classified_residues``.
    # Picked from the ColorBrewer Paired-12 set the coordinator already uses
    # so the halos sit in the same visual family as Chunk A/B's overlays.
    _CATEGORY_HALO = {
        "metal_site":          ("#e31a1c", "ffparam_metal_site"),       # dark red
        "modified_amino_acid": ("#33a02c", "ffparam_modified_aa"),      # dark green
        "small_molecule":      ("#ff7f00", "ffparam_small_mol"),        # dark orange
        "unknown":             ("#6a3d9a", "ffparam_unknown"),          # dark purple
    }

    def _halo_priority_structure(self):
        """Common helper: ensure the priority structure is loaded in the viewer.

        Returns the file path so the caller knows there's a structure to
        attach halos to (None if the workspace has no priority PDB).
        """
        try:
            from proprep.utils.structure_selector import get_priority_pdb_file
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        except Exception:
            return None
        path = get_priority_pdb_file(self.processor.workspace, silent=True)
        if not path:
            return None
        _viewer.show_structure(path)
        return path

    def _halo_classified_residues(self, residues):
        """Halo every non-standard residue, coloured by its classification.

        One halo rep per category (4 categories: metal_site / modified_amino_acid
        / small_molecule / unknown), each under its own stable label so a later
        re-analyze cleanly replaces the prior set. Prints a one-line legend
        to the terminal so the user knows what each colour represents.
        """
        if not residues or not self._halo_priority_structure():
            return
        from proprep.structure_prep.viewer_coordinator import viewer as _viewer

        by_category: Dict[str, List] = {}
        for res in residues:
            cat = getattr(res, "category", None) or "unknown"
            by_category.setdefault(cat, []).append(res)

        legend_parts = []
        for cat, items in by_category.items():
            color, label = self._CATEGORY_HALO.get(cat, ("#ffff00", f"ffparam_{cat}"))
            clauses = []
            for r in items:
                ch = getattr(r, "chain_id", "") or ""
                rid = getattr(r, "resid", None)
                if rid is None:
                    continue
                clauses.append(f"(:{ch} and {rid})")
            if not clauses:
                continue
            _viewer.highlight(
                " or ".join(clauses),
                style="halo",
                color=color,
                label=label,
            )
            legend_parts.append(
                f"[{color}]●[/{color}] {cat.replace('_', ' ')} ({len(items)})"
            )

        if legend_parts:
            self.console.print(
                "\n[grey50]Viewer halos:[/grey50] " + " · ".join(legend_parts)
            )

    def _halo_residue_instance(self, residue, *, label="ffparam_current",
                                color="#ffff00"):
        """Halo a single residue instance the user just picked from a table.

        Stable label ``ffparam_current`` so re-firing replaces the previous
        pick — every selection moves the highlight to the freshly chosen
        residue, leaving the per-category halos in place underneath.
        """
        if not residue or not self._halo_priority_structure():
            return
        ch = getattr(residue, "chain_id", "") or ""
        rid = getattr(residue, "resid", None)
        if rid is None:
            return
        from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        _viewer.highlight(
            f":{ch} and {rid}",
            style="halo",
            color=color,
            label=label,
        )

    def _halo_residues_by_name(self, residue_name, *, label="ffparam_current",
                                color="#ffff00"):
        """Halo every instance of a residue type — used by the resume picker.

        The resume table only knows the residue name (e.g. "HEM"), not a
        specific (chain, resid). Halo every instance via NGL's bracket
        syntax so the user sees all the residues that the resumed
        parameterization will (eventually) cover.
        """
        if not residue_name or not self._halo_priority_structure():
            return
        from proprep.structure_prep.viewer_coordinator import viewer as _viewer
        _viewer.highlight(
            f"[{residue_name}]",
            style="halo",
            color=color,
            label=label,
        )

    @staticmethod
    def _coord_key(coords):
        """Rounded coordinate tuple used to identify a RedoxSite atom.

        A site's atom names and residue numbers change under transformation but
        its coordinates never do, so coords are the permanent identifier (see
        :class:`RedoxSiteAtom`). 3 decimals matches PDB precision and the
        rounding used by ``sync_redox_sites_from_pdb``.
        """
        return tuple(round(float(c), 3) for c in coords)

    @staticmethod
    def _residue_key(res):
        """(chain, resid, resname) identity shared by members and site atoms."""
        return (res.chain_id, int(res.resid), (res.name or "").strip())

    def _site_residue_by_coords(self, site):
        """Map every site atom's coord key -> (chain, resid, resname)."""
        index = {}
        for atom in getattr(site, "atoms", []) or []:
            try:
                index[self._coord_key(atom.coords)] = (
                    atom.chain, int(atom.resid), (atom.resname or "").strip())
            except (TypeError, ValueError):
                continue
        return index

    def _bond_endpoint_residue(self, coords, residue_info, coord_index):
        """Resolve one bond endpoint to (chain, resid, resname), or None.

        Resolves through the site's atoms by coordinate first; ``residue_info``
        is only a fallback because it is not reliably a dict — a site
        deserialized by ``structure_preprocessor`` can carry ``""`` there.
        """
        hit = coord_index.get(self._coord_key(coords)) if coords else None
        if hit:
            return hit
        if isinstance(residue_info, dict) and residue_info.get("resname"):
            try:
                return (residue_info.get("chain"), int(residue_info["resid"]),
                        str(residue_info["resname"]).strip())
            except (TypeError, ValueError, KeyError):
                return None
        return None

    def _covalent_partner_members(self, site, members):
        """Site residues covalently bonded to a unit's members but not in it.

        :meth:`_extract_nonstandard_from_redox_sites` seeds a unit only from a
        site's redox CENTERS, and skips centers that are standard amino acids.
        A covalent adduct (a Cys plus the inhibitor bonded to its SG, say)
        therefore reaches classification as the ligand alone, and reads as a
        lone small molecule instead of the modified amino acid it is. Pull back
        in every residue the site joins to a member by an inter-residue COVALENT
        bond, so the unit is the whole conjugate.

        Only covalent bonds are followed: ``classify_bond_types`` marks a bond
        "covalent" exactly when neither atom is a metal, so coordinate and
        metal-metal contacts (a metal site's coordination sphere, which MCPB
        re-detects from the site itself) are left alone. Bonds are followed
        transitively, but never outside the site.

        Returns new :class:`NonStandardResidue` objects in discovery order.
        These belong to the unit only — they are deliberately NOT added to
        ``self.non_standard_residues``, since a standard Cys pulled in as a
        covalent partner is not itself a non-standard residue.
        """
        coord_index = self._site_residue_by_coords(site)

        edges = []
        for bond in getattr(site, "bonds", []) or []:
            if getattr(bond, "chemical_type", "") != "covalent":
                continue
            if getattr(bond, "bond_type", "") != "interresidue":
                continue
            a = self._bond_endpoint_residue(
                getattr(bond, "atom1_coords", None),
                getattr(bond, "atom1_residue_info", None), coord_index)
            b = self._bond_endpoint_residue(
                getattr(bond, "atom2_coords", None),
                getattr(bond, "atom2_residue_info", None), coord_index)
            if a and b and a != b:
                edges.append((a, b))
        if not edges:
            return []

        seeds = {self._residue_key(m) for m in members}
        reached = set(seeds)
        frontier = list(seeds)
        discovered = []  # keeps output order deterministic
        while frontier:
            nxt = []
            for a, b in edges:
                for x, y in ((a, b), (b, a)):
                    if x in frontier and y not in reached:
                        reached.add(y)
                        nxt.append(y)
                        discovered.append(y)
            frontier = nxt

        # Elements per residue, so _member_has_metal works without a structure.
        elements = defaultdict(set)
        for atom in getattr(site, "atoms", []) or []:
            try:
                key = (atom.chain, int(atom.resid), (atom.resname or "").strip())
            except (TypeError, ValueError):
                continue
            if atom.element:
                elements[key].add(atom.element.upper())

        partners = []
        for chain, resid, resname in discovered:
            ns_res = NonStandardResidue(
                name=resname, chain_id=chain, resid=resid, category="unknown")
            ns_res.redox_site_id = site.site_id
            ns_res.source_redox_site = site
            ns_res.is_redox_component = True
            ns_res.is_covalent_partner = True
            ns_res.redox_site_elements = elements.get((chain, resid, resname), set())
            partners.append(ns_res)
        return partners

    def _group_residues_into_display_units(self):
        """Group non-standard residues into parameterization units, classified.

        A unit is a RedoxSite (its member residues, plus any residue the site
        covalently bonds them to) or a single standalone residue. Each unit is
        then classified by content via :meth:`_classify_unit`, so every unit
        carries exactly one category.

        Returns a list of unit dicts in stable input order, each of the form
        ``{"members", "site_id", "category", "procedure", "metals", "ligands",
        "parent_residue", "unsupported_kind", "basis", ...}``.
        """
        structure = self._get_structure_object()

        # Bucket residues by the redox site they belong to (preserving order),
        # keeping loose residues (no site) in the same positional stream.
        site_members = OrderedDict()   # site_id -> [ns_res, ...]
        order = []                     # sequence of ("site", site_id) / ("loose", ns_res)
        loose_by_key = {}              # (chain, resid, resname) -> loose ns_res
        for res in self.non_standard_residues:
            site = getattr(res, "source_redox_site", None)
            site_id = getattr(res, "redox_site_id", None) or getattr(site, "site_id", None)
            if site is not None and site_id:
                if site_id not in site_members:
                    site_members[site_id] = []
                    order.append(("site", site_id))
                site_members[site_id].append(res)
            else:
                order.append(("loose", res))
                loose_by_key[self._residue_key(res)] = res

        # Complete each metal-free site with the residues it covalently bonds
        # its members to. A metal site is already whole: MCPB re-detects its own
        # coordination sphere from the site, and extraction seeds it with every
        # metal/ligand CENTER. Only a metal-free site can be missing a member,
        # because extraction skips a covalent partner that is a standard AA.
        absorbed = set()
        for site_id, members in site_members.items():
            site = next((getattr(m, "source_redox_site", None) for m in members
                         if getattr(m, "source_redox_site", None) is not None), None)
            if site is None or any(self._member_has_metal(m) for m in members):
                continue
            for partner in self._covalent_partner_members(site, members):
                key = self._residue_key(partner)
                existing = loose_by_key.get(key)
                if existing is not None:
                    # A NON-standard partner the structure scan already picked
                    # up on its own: it belongs to the conjugate, so take the
                    # scanned object (it carries ccd_data / atom counts) and
                    # drop its separate unit rather than parameterize it twice.
                    members.append(existing)
                    absorbed.add(key)
                else:
                    members.append(partner)

        units = []
        emitted_sites = set()
        for kind, ref in order:
            if kind == "loose":
                if self._residue_key(ref) in absorbed:
                    continue
                members, site_id = [ref], None
            else:
                if ref in emitted_sites:
                    continue
                emitted_sites.add(ref)
                members, site_id = site_members[ref], ref
            unit = {"members": members, "site_id": site_id}
            unit.update(self._classify_unit(members, structure))
            units.append(unit)
        return units

    def _format_unit(self, unit):
        """Return (name_desc, category_display, status) for a unit row.

        One place so the analysis table and the selection menu label a unit
        identically, keyed only on its content-derived category.
        """
        members = unit["members"]
        members_str = ", ".join(f"{m.name}({m.chain_id}:{m.resid})" for m in members)
        category = unit["category"]

        if category == "metal_site":
            metals = unit["metals"] or members
            metals_str = "+".join(f"{m.name}({m.chain_id}:{m.resid})" for m in metals)
            name_desc = f"Metal site: {metals_str}"
            ligands = unit.get("ligands") or []
            if ligands:
                name_desc += "\n  ligand: " + ", ".join(
                    f"{l.name} ({l.chain_id}:{l.resid})" for l in ligands)
            # Nuclearity counts metal ATOMS, not metal-bearing residues: an
            # Fe2S2 cluster is one FES residue holding two Fe, and counting
            # residues labelled it mononuclear.
            n = sum(self._count_metal_atoms(m) for m in metals) or len(metals)
            nuclearity = {1: "mononuclear", 2: "binuclear", 3: "trinuclear",
                          4: "tetranuclear"}.get(n, f"{n}-nuclear")
            return (name_desc, "Metal Site",
                    f"[green]MCPB ({nuclearity}{', +ligand' if ligands else ''})[/green]")

        if category == "modified_amino_acid":
            name_desc = members_str if len(members) > 1 else \
                f"{members[0].name} ({members[0].chain_id}:{members[0].resid})"
            proc = "conjugate, from structure" if unit.get("procedure") == "from_structure" else "de-novo"
            return (name_desc, "Modified Amino Acid",
                    f"[green]Modified AA Parameterizer ({proc})[/green]")

        if category == "small_molecule":
            name_desc = f"{members[0].name} ({members[0].chain_id}:{members[0].resid})" \
                if len(members) == 1 else members_str
            return (name_desc, "Small Molecule", "[green]Small Molecule Parameterizer available[/green]")

        if category == "unsupported":
            return (members_str, f"Unsupported ({unit.get('unsupported_kind','polymer')})",
                    "[yellow]Dedicated parameterizer not yet implemented[/yellow]")

        # unknown
        name_desc = f"{members[0].name} ({members[0].chain_id}:{members[0].resid})" \
            if len(members) == 1 else members_str
        return (name_desc, "Unknown", "[yellow]Needs classification first[/yellow]")

    def analyze_nonstandard_residues(self):
        """
        Analyze the structure for non-standard residues.

        New workflow:
        1. Extract non-standard residues from RedoxSites (if available)
        2. Scan structure for additional non-standard residues
        3. Merge and deduplicate
        4. Classify and enrich with CCD data
        """
        self.console.print("\n[bold]═══ Analyzing Non-Standard Residues ═══[/bold]\n")

        # Display classification methodology
        self._display_classification_methodology()

        # Get structure using priority system
        structure_data, structure_type, is_file_path = self._get_structure_with_priority()

        if not structure_data:
            self.console.print("[red]No structure available in workspace[/red]")
            return

        self.console.print(f"[grey50]Using structure: {structure_type}[/grey50]\n")

        # Parse structure if it's a file path
        structure = self._ensure_structure_object(structure_data, is_file_path)
        if not structure:
            return

        # Step 1: Extract from RedoxSites
        self.console.print("[bold cyan]Step 1: Extracting from RedoxSites[/bold cyan]")
        redox_nonstandard = self._extract_nonstandard_from_redox_sites()
        if redox_nonstandard:
            self.console.print(f"  [green]✓ Found {len(redox_nonstandard)} residue(s) from RedoxSite analysis[/green]")
            self.console.print("  [grey50]  (categories are decided per unit in Step 3, by content)[/grey50]")
            for res in redox_nonstandard:
                self.console.print(f"    • {res.name} ({res.chain_id}:{res.resid})")
        else:
            self.console.print("  [grey50]No RedoxSites detected in workspace[/grey50]")

        # Step 2: Scan structure for additional non-standard residues
        self.console.print("\n[bold cyan]Step 2: Scanning Structure[/bold cyan]")
        self.console.print("[grey50]  Iterating through all residues, skipping standard amino acids and water.[/grey50]")
        # Pass already-found residues to skip them in scan
        already_found_keys = {(res.chain_id, res.resid, res.name) for res in redox_nonstandard}
        structure_nonstandard = self._scan_structure_for_nonstandard(structure, skip_keys=already_found_keys)
        if structure_nonstandard:
            self.console.print(f"  [green]✓ Found {len(structure_nonstandard)} additional non-standard residue(s)[/green]")
            for res in structure_nonstandard:
                self.console.print(f"    • {res.name} ({res.chain_id}:{res.resid}) - {res.atom_count} atoms")
        else:
            self.console.print("  [grey50]No additional non-standard residues found[/grey50]")

        # Step 3: Classify residues
        self.console.print("\n[bold cyan]Step 3: Classification[/bold cyan]")
        self.console.print("[grey50]  Applying classification tests in priority order (see methodology above).[/grey50]")

        # Merge and deduplicate
        self.non_standard_residues = self._merge_residue_lists(
            redox_nonstandard,
            structure_nonstandard
        )

        # Enrich with CCD data
        self._enrich_with_ccd_data()

        # Store in workspace
        self.update_workspace("non_standard_residues", self.non_standard_residues)

        # Halo every detected non-standard residue, coloured by classification
        # category, so the user sees the analysis result in the structure pane
        # (not just in the terminal table).
        self._halo_classified_residues(self.non_standard_residues)

        # Display summary
        self.console.print()
        self.display_analysis_summary()

    def _display_classification_methodology(self):
        """Display the analysis process and classification methodology"""
        from rich.panel import Panel

        methodology = (
            "[bold]Analysis Process:[/bold]\n"
            "  1. Build parameterization units: each RedoxSite's residues are grouped\n"
            "     into ONE unit (so a covalent AA+ligand adduct is analyzed together);\n"
            "     every other non-standard residue is its own unit.\n"
            "  2. Classify each unit by its STRUCTURAL CONTENT (atoms / backbone / CCD),\n"
            "     not by any label carried over from redox detection.\n\n"
            "[bold]Classification order (first match wins):[/bold]\n\n"
            "[cyan]1. User classification[/cyan] [grey50](highest priority)[/grey50]\n"
            "   A manual category you've set for a residue name overrides every test below.\n\n"
            "[cyan]2. Metal content[/cyan] → [yellow]metal_site[/yellow] (MCPB)\n"
            "   Test: any residue in the unit contains a metal ATOM (per-atom element,\n"
            "   or a CCD / periodic-table metal-ion match). One rule covers lone ions\n"
            "   AND organometallic cofactors (heme, Fe-S).\n\n"
            "[cyan]3. Peptide backbone[/cyan] → [yellow]modified_amino_acid[/yellow]\n"
            "   Test: a residue carries a full N, CA, C, O backbone.\n"
            "     • Unit has >1 residue (a RedoxSite bundled the AA with a covalent\n"
            "       partner) → from-structure (conjugate) route.\n"
            "     • Single residue → de-novo route; the parent AA is then looked up\n"
            "       (local PTM map, else CCD parent) to build the capped model.\n\n"
            "[cyan]4. Recognized biopolymer[/cyan] → [yellow]unsupported[/yellow]\n"
            "   Test: CCD classifies it as nucleic-acid or carbohydrate. Flagged rather\n"
            "   than forced through GAFF — there is no dedicated parameterizer yet.\n\n"
            "[cyan]5. Small organic molecule[/cyan] → [yellow]small_molecule[/yellow] (GAFF2)\n"
            "   Test: atom count within the configurable window AND CCD is non-polymer /\n"
            "   not an explicit peptide·DNA·RNA·saccharide.\n\n"
            "[cyan]6. Otherwise[/cyan] → [yellow]unknown[/yellow] (requires manual classification)"
        )

        self.console.print(Panel(methodology, title="Analysis Methodology",
                                  border_style="blue", expand=False))
        self.console.print()

    def _ensure_structure_object(self, structure_data, is_file_path):
        """Convert file path to BioPython Structure object if needed"""
        if is_file_path:
            from Bio.PDB import PDBParser
            parser = PDBParser(QUIET=True)
            try:
                return parser.get_structure("protein", structure_data)
            except Exception as e:
                self.console.print(f"[red]Error parsing structure file: {e}[/red]")
                return None
        else:
            return structure_data

    def _scan_structure_for_nonstandard(self, structure, skip_keys=None):
        """Scan BioPython structure for non-standard residues.

        Args:
            structure: BioPython structure object
            skip_keys: Optional set of (chain_id, resid, resname) tuples to skip
                       (e.g., residues already found from RedoxSite extraction)
        """
        nonstandard = {}
        skip_keys = skip_keys or set()

        for model in structure:
            for chain in model:
                chain_id = chain.id
                for residue in chain:
                    resname = residue.get_resname().strip()
                    resid = residue.id[1]

                    # Skip standard amino acids and water
                    if resname in self.standard_aa or resname in ["HOH", "WAT"]:
                        continue

                    key = (chain_id, resid, resname)

                    # Skip if already found (e.g., from RedoxSite extraction)
                    if key in skip_keys:
                        continue

                    if key not in nonstandard:
                        # Count atoms
                        atom_count = len(list(residue.get_atoms()))

                        ns_res = NonStandardResidue(
                            name=resname,
                            chain_id=chain_id,
                            resid=resid,
                            category="unknown",  # Will be classified later
                        )
                        ns_res.atom_count = atom_count
                        ns_res.biopython_residue = residue  # Store for backbone detection
                        nonstandard[key] = ns_res

            break  # Only process first model

        return list(nonstandard.values())

    def _merge_residue_lists(self, redox_list, structure_list):
        """Merge and deduplicate residue lists, prioritizing redox-derived data"""
        merged = {}

        # Add redox-derived residues first (higher priority)
        for res in redox_list:
            key = (res.chain_id, res.resid, res.name)
            merged[key] = res

        # Add structure-derived residues if not already present
        for res in structure_list:
            key = (res.chain_id, res.resid, res.name)
            if key not in merged:
                merged[key] = res

        return list(merged.values())

    def _enrich_with_ccd_data(self):
        """Fetch CCD data, then stamp each residue with its UNIT's category.

        A residue's category has to come from the unit it will actually be
        parameterized in, not from the residue in isolation: a covalent
        adduct's ligand on its own looks like a small molecule, and only the
        unit (the ligand plus the amino acid it is bonded to) shows the
        modified amino acid. Classifying per unit here keeps the per-residue
        readers (the classification editor, workspace storage, halos) agreeing
        with the unit table rather than contradicting it.

        CCD data is fetched first because :meth:`_classify_unit` reads it.
        """
        for res in self.non_standard_residues:
            if getattr(res, "ccd_data", None) is None:
                res.ccd_data = self.ccd_parser.get_residue_data(res.name)

        self.console.print()
        for unit in self._group_residues_into_display_units():
            category, basis = unit["category"], unit.get("basis", "")
            for res in unit["members"]:
                # A covalent partner belongs to the unit only, never to
                # self.non_standard_residues — there is nothing to stamp.
                if getattr(res, "is_covalent_partner", False):
                    continue
                res.category = category
                res.classification_basis = basis
                if category == "modified_amino_acid" and unit.get("parent_residue"):
                    res.parent_residue = unit["parent_residue"]
            label = ", ".join(
                f"{m.name} ({m.chain_id}:{m.resid})" for m in unit["members"])
            self.console.print(
                f"  [green]✓[/green] {label}: "
                f"[yellow]{category.replace('_', ' ')}[/yellow] "
                f"[grey50]— {basis}[/grey50]"
            )
    def _is_metal_ion(self, resname, ccd_data):
        """Check if residue is a metal ion using comprehensive METALS dictionary"""
        return resname.upper() in METALS

    def _contains_metal_atoms(self, res) -> Tuple[bool, List[str]]:
        """
        Check if a residue contains any metal atoms.

        This is used to detect metal-containing small molecules (like heme)
        that need MCPB-style parameterization instead of pure GAFF2.

        Args:
            res: NonStandardResidue object (may have biopython_residue or redox_site_elements)

        Returns:
            Tuple of (contains_metal: bool, metal_elements: List[str])
        """
        metal_elements_found = []

        # Method 1: Check redox_site_elements (from RedoxSite extraction)
        if hasattr(res, 'redox_site_elements') and res.redox_site_elements:
            for element in res.redox_site_elements:
                if element.upper() in METALS:
                    if element.upper() not in metal_elements_found:
                        metal_elements_found.append(element.upper())

        # Method 2: Check BioPython residue object (from structure scan)
        if not metal_elements_found and hasattr(res, 'biopython_residue') and res.biopython_residue is not None:
            for atom in res.biopython_residue.get_atoms():
                # Get element from atom - BioPython stores it in atom.element
                element = getattr(atom, 'element', None)
                if element is None:
                    # Fallback: infer from atom name (first 1-2 chars)
                    atom_name = atom.get_name().strip()
                    element = atom_name[0:2].strip() if len(atom_name) >= 2 else atom_name[0:1]

                element_upper = element.upper().strip()

                if element_upper in METALS:
                    if element_upper not in metal_elements_found:
                        metal_elements_found.append(element_upper)

        return (len(metal_elements_found) > 0, metal_elements_found)

    def _count_metal_atoms(self, res) -> int:
        """Number of metal ATOMS in a residue.

        Nuclearity is a property of the metal centre, not of how many residues
        carry one: an Fe2S2 cluster is a single FES residue holding two Fe, and
        counting residues calls it mononuclear. _contains_metal_atoms cannot
        answer this — it returns the distinct metal ELEMENTS, so FES yields
        ['FE'] whether it holds one iron or four.

        Sources in order: the RedoxSite's atom list (exact, and available
        without a structure), then the BioPython residue, then a floor of 1 for
        a residue known to hold a metal but whose atoms we cannot enumerate.
        """
        site = getattr(res, "source_redox_site", None)
        if site is not None:
            try:
                want = (res.chain_id, int(res.resid), (res.name or "").strip())
            except (TypeError, ValueError):
                want = None
            if want is not None:
                count = 0
                for atom in getattr(site, "atoms", []) or []:
                    try:
                        key = (atom.chain, int(atom.resid), (atom.resname or "").strip())
                    except (TypeError, ValueError):
                        continue
                    if key == want and atom.element and atom.element.upper() in METALS:
                        count += 1
                if count:
                    return count

        bp = getattr(res, "biopython_residue", None)
        if bp is not None:
            count = 0
            for atom in bp.get_atoms():
                element = getattr(atom, "element", None)
                if element is None:
                    name = atom.get_name().strip()
                    element = name[0:2].strip() if len(name) >= 2 else name[0:1]
                if element and element.upper().strip() in METALS:
                    count += 1
            if count:
                return count

        return 1 if self._member_has_metal(res) else 0

    def _unit_key(self, res) -> str:
        """Key identifying the parameterization unit a residue belongs to.

        Every member of a metal-site unit is stamped ``category="metal_site"``,
        the coordinating ligand included, so counting residues by category
        overcounts sites — 4UHX reported "4 metal sites" for three (its MoCo
        site contributes both MOS and its MTE ligand). Residues extracted from
        a RedoxSite share its id; a standalone residue is its own unit.
        """
        site_id = getattr(res, "redox_site_id", None)
        if site_id:
            return f"site:{site_id}"
        return f"res:{self._residue_key(res)}"

    def _is_small_molecule_from_ccd(self, ccd_data, atom_count):
        """Check if residue is a small molecule based on settings and CCD data"""
        # Check atom count range first
        min_atoms = self.classification_settings.get("min_small_molecule_atoms", 2)
        max_atoms = self.classification_settings.get("max_small_molecule_atoms", 200)

        if not (min_atoms <= atom_count <= max_atoms):
            return False

        # If CCD data available, check type
        if ccd_data:
            ccd_type = ccd_data.get("type", "").upper()
            if "NON-POLYMER" in ccd_type or "HETAIN" in ccd_type:
                return True
            # If CCD explicitly says it's something else (e.g., PEPTIDE), trust it
            if any(keyword in ccd_type for keyword in ["PEPTIDE", "DNA", "RNA", "SACCHARIDE"]):
                return False

        # No CCD or inconclusive - use atom count alone
        return True

    def _has_peptide_backbone(self, residue):
        """
        Check if a residue has the canonical peptide backbone atoms (N, CA, C, O).
        This indicates it's likely a modified amino acid rather than a small molecule.

        Args:
            residue: BioPython residue object

        Returns:
            bool: True if all backbone atoms are present
        """
        try:
            backbone_atoms = {"N", "CA", "C", "O"}
            residue_atoms = {atom.get_name().strip().upper() for atom in residue.get_atoms()}
            return backbone_atoms.issubset(residue_atoms)
        except Exception:
            return False

    # ── Content-based unit classification ──────────────────────────────────
    # A "unit" is a RedoxSite (one or more members) or a single standalone
    # residue. Grouping is decided elsewhere (the RedoxSite bundles; a lone
    # residue is its own unit). _classify_unit decides the CATEGORY by
    # structural content — never by the detector's center_type, which is only a
    # hint. Exactly one of: metal_site, modified_amino_acid, small_molecule
    # (the three parameterizers), plus unknown ("auto couldn't decide, you
    # pick") and unsupported (recognized nucleic acid / carbohydrate, no
    # dedicated parameterizer yet).

    def _member_biopython_residue(self, res, structure=None):
        """BioPython residue for a NonStandardResidue: its cached
        ``biopython_residue`` (structure-scanned) or a lookup by (chain, resid)
        (RedoxSite-extracted members carry no cached residue)."""
        bp = getattr(res, "biopython_residue", None)
        if bp is not None:
            return bp
        structure = structure or self._get_structure_object()
        if not structure:
            return None
        try:
            model = next(iter(structure))
        except StopIteration:
            return None
        for chain in model:
            if chain.id != res.chain_id:
                continue
            for r in chain:
                if r.id[1] == int(res.resid) and r.resname.strip() == res.name.strip():
                    return r
        return None

    def _member_has_backbone(self, res, structure=None):
        """True if the residue has a full N,CA,C,O peptide backbone."""
        bp = self._member_biopython_residue(res, structure)
        return self._has_peptide_backbone(bp) if bp is not None else False

    def _member_has_metal(self, res):
        """True if the residue is a metal ion or contains any metal atom."""
        if self._is_metal_ion(res.name, getattr(res, "ccd_data", None)):
            return True
        has_metal, _ = self._contains_metal_atoms(res)
        return has_metal

    def _member_atom_count(self, res, structure=None):
        """Heavy+H atom count, from the cached count or a structure lookup."""
        n = getattr(res, "atom_count", 0)
        if n:
            return n
        bp = self._member_biopython_residue(res, structure)
        return len(list(bp.get_atoms())) if bp is not None else 0

    def _lookup_parent_aa(self, res):
        """Parent standard AA for a modified residue: static mapping table, then
        a live CCD ``mon_nstd_parent_comp_id`` query, else None (prompt later)."""
        if res.name.upper() in self.modified_aa_map:
            return self.modified_aa_map[res.name.upper()]
        ccd = getattr(res, "ccd_data", None) or self.ccd_parser.get_residue_data(res.name)
        if ccd:
            parent = (ccd.get("mon_nstd_parent_comp_id") or "").strip()
            if parent in self.standard_aa:
                return parent
        return None

    def _ccd_polymer_kind(self, res):
        """'nucleic acid' / 'carbohydrate' if CCD marks this an unsupported
        biopolymer residue, else None."""
        ccd = getattr(res, "ccd_data", None) or self.ccd_parser.get_residue_data(res.name)
        if not ccd:
            return None
        ccd_type = (ccd.get("type") or "").upper()
        if "RNA" in ccd_type or "DNA" in ccd_type:
            return "nucleic acid"
        if "SACCHARIDE" in ccd_type:
            return "carbohydrate"
        return None

    def _classify_unit(self, members, structure=None):
        """Classify a parameterization unit by structural content.

        ``members`` is a RedoxSite's residues (>=1) or a single standalone
        residue. Returns a dict: category, procedure (modAA only),
        metals/ligands (metal_site routing), parent_residue (modAA de-novo),
        unsupported_kind, basis.
        """
        structure = structure or self._get_structure_object()

        # 0. User manual override wins (first overridden member).
        for m in members:
            override = self.user_residue_classifications.get(m.name.upper())
            if override:
                cat = override.get("category", "unknown")
                return {
                    "category": cat,
                    "procedure": "de_novo" if cat == "modified_amino_acid" else None,
                    "metals": [x for x in members if self._member_has_metal(x)],
                    "ligands": [], "parent_residue": override.get("parent_residue"),
                    "unsupported_kind": None,
                    "basis": "User classification (manual override)",
                }

        # 1. Metal content → metal_site (MCPB). Covers lone ions and
        #    organometallic cofactors (heme, Fe-S) alike.
        metals = [m for m in members if self._member_has_metal(m)]
        if metals:
            ligands = [m for m in members if m not in metals]
            return {
                "category": "metal_site", "procedure": None,
                "metals": metals, "ligands": ligands, "parent_residue": None,
                "unsupported_kind": None,
                "basis": "Contains metal atom(s) → MCPB",
            }

        # 2. Peptide backbone → modified amino acid. A multi-member unit means
        #    the RedoxSite bundled the AA with a covalent partner → conjugate →
        #    from-structure; a lone modified residue → de-novo (parent lookup).
        aa_members = [m for m in members if self._member_has_backbone(m, structure)]
        if aa_members:
            if len(members) > 1:
                procedure, parent = "from_structure", None
                basis = "Modified AA conjugate (RedoxSite bundles AA + ligand)"
            else:
                procedure = "de_novo"
                parent = self._lookup_parent_aa(aa_members[0])
                basis = "Modified amino acid (peptide backbone)"
            return {
                "category": "modified_amino_acid", "procedure": procedure,
                "metals": [], "ligands": [], "parent_residue": parent,
                "aa_member": aa_members[0], "unsupported_kind": None, "basis": basis,
            }

        # 3. Organic. Recognized biopolymers (nucleic acid / carbohydrate) have
        #    no dedicated parameterizer yet — flag rather than mis-GAFF them.
        for m in members:
            kind = self._ccd_polymer_kind(m)
            if kind:
                return {
                    "category": "unsupported", "procedure": None,
                    "metals": [], "ligands": [], "parent_residue": None,
                    "unsupported_kind": kind,
                    "basis": f"CCD type: {kind} (no dedicated parameterizer yet)",
                }

        # 4. Small organic molecule, else unknown (auto couldn't decide).
        if all(self._is_small_molecule_from_ccd(
                getattr(m, "ccd_data", None), self._member_atom_count(m, structure))
                for m in members):
            return {
                "category": "small_molecule", "procedure": None,
                "metals": [], "ligands": [], "parent_residue": None,
                "unsupported_kind": None, "basis": "Organic small molecule",
            }

        return {
            "category": "unknown", "procedure": None,
            "metals": [], "ligands": [], "parent_residue": None,
            "unsupported_kind": None,
            "basis": "Auto-classification could not decide",
        }

    def _offer_workspace_save_and_exit(self):
        """Offer user options to save workspace and exit or background ProPrep while Gaussian runs"""
        from proprep.utils.prompts import prompt_with_context
        import sys

        self.console.print("\n[bold cyan]═══ Gaussian Calculation Required ═══[/bold cyan]")
        self.console.print("\n[yellow]ProPrep has paused and needs external Gaussian calculations.[/yellow]")
        self.console.print("\nWhat would you like to do?\n")
        self.console.print("  [cyan]1.[/cyan] Save workspace and exit ProPrep")
        self.console.print("     → Run Gaussian externally, restart ProPrep later, load workspace, resume")
        self.console.print("\n  [cyan]2.[/cyan] Keep ProPrep running (minimize terminal)")
        self.console.print("     → Run Gaussian in another terminal, return here when done")
        self.console.print("\n  [cyan]3.[/cyan] Continue (return to menu)")
        self.console.print("     → Handle the Gaussian calculations yourself\n")

        choice = prompt_with_context(
            self.processor,
            "Select option",
            choices=["1", "2", "3"],
            default="1",
            module="Force Field Parameterizer",
            description="Gaussian calculation workflow options",
            options_map={
                "1": "Save workspace and exit",
                "2": "Keep ProPrep running",
                "3": "Continue to menu"
            }
        )

        if choice == "1":
            # Save workspace and exit
            filename = prompt_with_context(
                self.processor,
                "\nEnter filename to save workspace",
                default="workspace.json",
                module="Force Field Parameterizer",
                description="Workspace filename for saving"
            )

            # Save workspace
            if hasattr(self.processor, 'save_workspace'):
                success = self.processor.save_workspace(filename)
                if success:
                    self.console.print(f"\n[green]✓ Workspace saved to {filename}[/green]")
                    self.console.print("\n[bold]To resume later:[/bold]")
                    self.console.print("  1. Run your Gaussian calculations")
                    self.console.print(f"  2. Restart ProPrep")
                    self.console.print(f"  3. Load workspace: Main Menu → Workspace Options → Load Workspace")
                    self.console.print(f"  4. Resume: Force Field Parameterizer → Resume pending parameterization\n")

                    # Ask if they want to exit now
                    exit_now = confirm_with_context(
                        self.processor,
                        "Exit ProPrep now?",
                        default=True,
                        module="Force Field Parameterizer",
                        description="Exit ProPrep now",
                    )
                    if exit_now:
                        self.console.print("[cyan]Goodbye! 👋[/cyan]")
                        sys.exit(0)
                else:
                    self.console.print("[red]Failed to save workspace[/red]")

        elif choice == "2":
            # Keep running - just inform user
            self.console.print("\n[green]✓ ProPrep will remain running[/green]")
            self.console.print("\n[bold]Instructions:[/bold]")
            self.console.print("  1. Open a new terminal")
            self.console.print("  2. Run your Gaussian calculations")
            self.console.print("  3. Return to this terminal when done")
            self.console.print("  4. Use option 7 (Resume pending parameterization) to continue\n")
            input("Press Enter to return to menu...")

        elif choice == "3":
            # Just return
            self.console.print("\n[yellow]Returning to menu. Remember to complete Gaussian calculations.[/yellow]")
            input("Press Enter to continue...")

    # ===== LEGACY METAL SITE HANDLING - DEPRECATED =====
    # The methods below were part of the old MetalSite-based workflow.
    # They are kept for reference but are no longer called.
    # All functionality now uses RedoxSite objects from detected_redox_sites.

    def display_analysis_summary(self):
        """Display analysis summary with enhanced formatting for mixed residue types"""
        if not self.non_standard_residues:
            self.console.print("[yellow]No non-standard residues found[/yellow]")
            return

        # One unit = one RedoxSite (its members) or one standalone residue,
        # each classified into a single category by content.
        units = self._group_residues_into_display_units()

        counts = defaultdict(int)
        for u in units:
            counts[u["category"]] += 1

        self.console.print(f"\n[bold]Found {len(units)} parameterization unit(s):[/bold]")
        self.console.print(f"  Modified amino acids: {counts['modified_amino_acid']}")
        self.console.print(f"  Small molecules: {counts['small_molecule']}")
        self.console.print(f"  Metal sites: {counts['metal_site']}")
        if counts["unsupported"]:
            self.console.print(f"  Unsupported (nucleic acid / carbohydrate): {counts['unsupported']}")
        self.console.print(f"  Unknown: {counts['unknown']}")

        # Create enhanced table
        table = Table(title="Available Non-Standard Residues")
        table.add_column("No.", style="cyan", justify="right")
        table.add_column("Residue", style="magenta")
        table.add_column("Category", style="yellow")
        table.add_column("Classification Basis", style="grey50")
        table.add_column("Status", style="blue")

        for i, unit in enumerate(units, 1):
            name_desc, category_display, status = self._format_unit(unit)
            basis = unit.get("basis", "")
            if len(basis) > 35:
                basis = basis[:32] + "..."
            table.add_row(str(i), name_desc, category_display, basis, status)

        self.console.print(table)

        # Add explanatory information
        self.console.print("\n[bold]Parameterization Options:[/bold]")

        if counts["modified_amino_acid"] > 0:
            self.console.print(
                "[green]• Modified amino acids[/green] → modified AA workflow "
                "(de-novo, or from-structure for a covalent conjugate)"
            )
        if counts["small_molecule"] > 0:
            self.console.print(
                "[green]• Small molecules[/green] → small molecule (GAFF) workflow"
            )
        if counts["metal_site"] > 0:
            self.console.print(
                "[green]• Metal sites[/green] → MCPB metal site workflow "
                "(a poly-nuclear site and its coordinating ligands are one unit)"
            )
        if counts["unsupported"] > 0:
            self.console.print(
                "[yellow]• Unsupported[/yellow] → nucleic acid / carbohydrate; "
                "no dedicated parameterizer yet"
            )
        if counts["unknown"] > 0:
            self.console.print(
                "[yellow]• Unknown[/yellow] → classify first "
                "(choose one of the three parameterizers)"
            )
    def parameterize_residue(self):
        """Interactive function to parameterize a specific non-standard residue - central hub for all parameterization strategies"""
        # Try to get residues from workspace if not already loaded
        if not self.non_standard_residues:
            self.non_standard_residues = self.get_from_workspace(
                "non_standard_residues", []
            )

            if not self.non_standard_residues:
                self.console.print(
                    "[yellow]No non-standard residues analyzed yet. Running analysis...[/yellow]"
                )
                self.analyze_nonstandard_residues()

                if not self.non_standard_residues:
                    self.console.print(
                        "[red]No non-standard residues found to parameterize[/red]"
                    )
                    return

        # Create selection table - no more grouping for metal sites
        self.console.print("\n[bold]Select a residue to parameterize:[/bold]")
        table = Table(title="Available Non-Standard Residues")
        table.add_column("No.", style="cyan", justify="right")
        table.add_column("Name/Description", style="magenta")
        table.add_column("Category", style="yellow")
        table.add_column("Instances", style="green", justify="right")
        table.add_column("Status", style="blue")

        # One selectable entry per parameterization unit (a RedoxSite's members
        # collapse to one row; a standalone residue is its own row).
        units = self._group_residues_into_display_units()
        options = [("unit", unit) for unit in units]

        for i, (_, unit) in enumerate(options, 1):
            name_desc, category_display, status = self._format_unit(unit)
            instances = str(len(unit["members"]))
            table.add_row(str(i), name_desc, category_display, instances, status)

        self.console.print(table)

        # The selection names the SITES to parameterize. It does not group them:
        # residues are grouped into sites by the Redox Site Detector, which is
        # the authority on what a unit is. Each selected site goes to the
        # parameterizer its own category implies, so a mixed selection is
        # ordinary rather than a conflict.
        self.console.print(
            "\n[cyan]Select the sites to parameterize — one (e.g. '3') or "
            "several (e.g. '1,2,3'). Each goes to its own parameterizer.[/cyan]")
        choice = prompt_with_context(
            self.processor,
            "Enter site number(s) to parameterize (or 'q' to cancel)",
            default="q",
            module="Force Field Parameterizer",
            description="Select site number(s) to parameterize",
        )

        if choice.lower() == "q":
            self.console.print("[yellow]Parameterization cancelled[/yellow]")
            return

        # '+' is still accepted as a separator: it was the combine syntax, and
        # old session logs replay through here.
        indices = []
        for token in choice.replace("+", ",").split(","):
            token = token.strip()
            if not token:
                continue
            try:
                idx = int(token) - 1
            except ValueError:
                self.console.print(
                    f"[red]Invalid input: '{token}'. Enter site numbers like "
                    f"'3' or '1,2,3'[/red]")
                return
            if idx < 0 or idx >= len(options):
                self.console.print(f"[red]Invalid selection: {token} is out of range[/red]")
                return
            if idx not in indices:
                indices.append(idx)

        if not indices:
            self.console.print("[yellow]Nothing selected[/yellow]")
            return

        self._parameterize_selected_units([options[i][1] for i in indices], units)

    def _parameterize_selected_units(self, selected_units, all_units):
        """Send each selected site to the parameterizer its category implies.

        Metal sites are handled as a group because MCPB has to number their
        M*/Y* atom types against each other; everything else is independent, so
        it is routed one unit at a time. Nothing is combined here — a unit is
        already whole.

        Args:
            selected_units: The units the user chose.
            all_units: Every displayed unit, so unselected metal sites can be
                named (the topology build needs every metal site to have
                parameters, whether derived here or reused).
        """
        metal_units = [u for u in selected_units
                       if u["category"] in ("metal_site", "metal_ion")]
        other_units = [u for u in selected_units if u not in metal_units]

        if metal_units:
            self._parameterize_metal_sites(metal_units, all_units)

        for unit in other_units:
            label = ", ".join(f"{m.name} ({m.chain_id}:{m.resid})"
                              for m in unit["members"])
            self.console.print(f"\n[bold cyan]Parameterizing: {label}[/bold cyan]")
            try:
                self._route_unit(unit)
            except Exception as exc:
                # One site failing must not cancel the sites selected after it.
                self.console.print(
                    f"[red]Parameterization of {label} failed: {exc}[/red]")

    def _route_unit(self, unit):
        """Dispatch one parameterization unit to the correct parameterizer.

        HANDOFF INVARIANT: once a unit reaches a parameterizer, behavior is
        exactly as it is today — this method only chooses the door. A unit is
        already whole (the Redox Site Detector grouped it), so nothing is
        combined or split here. A metal site goes through _parameterize_metal_sites
        so its site id is recorded and MCPB parameterizes it and no others.
        """
        category = unit["category"]
        members = unit["members"]

        if category == "metal_site":
            # Through the same door as a multi-site selection, so the site id is
            # recorded and MCPB parameterizes this site and no others.
            self._parameterize_metal_sites([unit])

        elif category == "modified_amino_acid":
            if unit.get("procedure") == "from_structure":
                members_name = "+".join(m.name for m in members)
                capped = self._build_capped_conjugate_pdb(
                    members, members_name,
                    redox_site=getattr(members[0], "source_redox_site", None),
                )
                if capped:
                    conformer_pdbs, mod_resname = capped
                    self.console.print(
                        f"[cyan]Launching modified AA parameterization for covalent "
                        f"adduct {members_name} (residue {mod_resname})...[/cyan]"
                    )
                    self.parameterize_modified_amino_acid(
                        mod_resname, members, combined_pdb=conformer_pdbs
                    )
                else:
                    self.console.print(
                        "[yellow]Could not build a capped model compound "
                        "(chain terminus / missing neighbour?).[/yellow]"
                    )
            else:
                # De-novo single modified residue (Route A) — unchanged handoff.
                aa = unit.get("aa_member", members[0])
                if unit.get("parent_residue"):
                    aa.parent_residue = unit["parent_residue"]
                self.console.print(
                    f"[cyan]Launching modified AA parameterization for {aa.name}...[/cyan]")
                self.parameterize_modified_amino_acid(aa.name, [aa])

        elif category == "small_molecule":
            for member in members:
                self.console.print(
                    f"[cyan]Launching small molecule parameterization for {member.name}...[/cyan]")
                self.parameterize_small_molecule(member.name, [member])

        elif category == "unsupported":
            kind = unit.get("unsupported_kind", "biopolymer")
            names = ", ".join(m.name for m in members)
            self.console.print(
                f"[yellow]{names}: recognized as {kind}, but a dedicated "
                f"parameterizer is not yet implemented.[/yellow]")

        else:  # unknown
            names = ", ".join(m.name for m in members)
            self.console.print(
                f"[yellow]{names} needs classification first.[/yellow]\n"
                "[yellow]Use 'Change classification' to assign a parameterizer, "
                "then try again.[/yellow]")

    def _offer_library_promotion(self, *, category, residue_name, frcmod_file,
                                 lib_search_dir, prep_file=None, lib_file=None,
                                 atom_types=None):
        """Offer to save a finished parameterization into the user library.

        Best-effort and fully guarded: a failure here must never break the
        parameterization workflow, which has already succeeded and stored its
        results in the workspace. Returns the deposit result dict (with
        ``library_path``/``metadata_path``) when a promotion occurred, else
        None, so callers can locate the deposited library (e.g. to emit a reuse
        transformer that points at it).

        Thin wrapper over library_promotion.offer_library_promotion, which the
        parameterizer checklists' final steps call directly.
        """
        from proprep.forcefield_prep.library_promotion import offer_library_promotion
        return offer_library_promotion(
            self.console, self.processor,
            category=category, residue_name=residue_name,
            frcmod_file=frcmod_file, lib_search_dir=lib_search_dir,
            prep_file=prep_file, lib_file=lib_file, atom_types=atom_types,
        )

    def _save_combined_residues_to_pdb(self, residues: List[NonStandardResidue], combined_name: str) -> Optional[str]:
        """Extract combined residues from structure and save to PDB file"""
        try:
            from Bio.PDB import PDBIO, Select

            # Get structure
            structure = self._get_structure_object()
            if not structure:
                self.console.print("[red]No structure available in workspace[/red]")
                return None

            # Create custom selector for the specific residues
            class CombinedResidueSelector(Select):
                def __init__(self, target_residues):
                    self.targets = {(r.chain_id, r.resid) for r in target_residues}

                def accept_residue(self, residue):
                    chain_id = residue.get_parent().id
                    resid = residue.id[1]
                    return (chain_id, resid) in self.targets

            # Create output directory
            output_dir = Path("combined_residues")
            output_dir.mkdir(exist_ok=True)

            # Save PDB
            pdb_filename = f"{combined_name.replace('+', '_')}.pdb"
            pdb_path = output_dir / pdb_filename

            io = PDBIO()
            io.set_structure(structure)
            io.save(str(pdb_path), CombinedResidueSelector(residues))

            self.console.print(f"[green]✓ Saved combined residues to: {pdb_path}[/green]")
            return str(pdb_path)

        except Exception as e:
            self.console.print(f"[red]Error saving combined residues: {e}[/red]")
            return None

    # Heavy-atom cap layout the sep-bond detector expects post-hydrogenation:
    #   ACE (first): C, O, CH3  (+ reduce adds HH31/HH32/HH33 -> 6 atoms)
    #   NME (last):  N, CH3     (+ reduce adds H, HH31/HH32/HH33 -> 6 atoms)
    _BACKBONE = ("N", "CA", "C", "O")

    def _name_conjugate_residue(self, src_resname, aa_resid, residues):
        """Pick the 3-character residue name for a covalent AA↔ligand adduct.

        The adduct is a new chemical entity, so it must not inherit the parent's
        name: a library unit called CYS would shadow standard cysteine, and the
        working directory / file prefixes / AC file would all claim to be plain
        Cys. The default follows the same MCPB convention the integration step
        uses (first letter + last letter + counter, so CYS → CS1), checked
        against the names tLEaP already knows and against the residues in this
        conjugate, then offered for the user to override with something
        meaningful.

        Returns an uppercase name of at most 3 characters.
        """
        from proprep.forcefield_prep.mcpb.integration_utils import (
            generate_unique_residue_names)
        from proprep.forcefield_prep.structure_preprocessor import STANDARD_RESIDUES
        from proprep.forcefield_params.loader import get_registered_residue_names

        # Names that must not be reused: everything tLEaP already resolves, plus
        # every residue making up this conjugate (the ligand keeps its own code
        # in the model file, so the adduct cannot take it either).
        reserved = set(STANDARD_RESIDUES) | {"HOH", "WAT", "HEM"}
        reserved |= {(r.name or "").strip().upper() for r in residues}

        # Names already claimed by residues in the forcefield library (bundled +
        # ~/.proprep). Unlike `reserved`, these are a SOFT constraint: the same
        # adduct re-parameterized legitimately wants its old name back, so a
        # collision here warns and asks rather than hard-rejecting. But we still
        # want the DEFAULT to skip past them, so a second Cys adduct from another
        # protein does not silently default to CS1 again and clobber the first at
        # deposit time. {NAME: owning-cofactor-path}; empty if the lib is unreadable.
        library_names = get_registered_residue_names()

        # NB: generate_unique_residue_names ADDS what it generates to the set it
        # is handed (so a multi-residue call avoids self-collisions). Give it a
        # copy, or `reserved` would come back containing the very default we are
        # about to offer — and the validation below would reject it forever.
        # Include the library names so the counter advances past them (CS1 taken
        # by an existing adduct -> the default becomes CS2).
        # Seed the fallback from a letter, never a digit: an empty/degenerate
        # src_resname must not yield "1" (a leading-digit tLEaP hazard).
        seed = "".join(c for c in src_resname if c.isalpha())[:2] or "XX"
        default = generate_unique_residue_names(
            [(aa_resid, src_resname)],
            existing_names=set(reserved) | set(library_names),
        ).get((aa_resid, src_resname)) or f"{seed}1"

        self.console.print(
            f"\n[grey50]The covalent adduct of {src_resname} and "
            f"{', '.join(r.name for r in residues if (r.name or '').strip().upper() != src_resname)} "
            f"is a new residue — it cannot be called {src_resname}, which tLEaP already "
            f"resolves to the standard amino acid.[/grey50]"
        )
        while True:
            name = prompt_with_context(
                self.processor,
                "Residue name for the conjugate (max 3 characters)",
                default=default,
                module="Force Field Parameterizer",
                description="Residue name for the covalent adduct",
            ).strip().upper()
            if not name or len(name) > 3 or not name.isalnum():
                self.console.print(
                    "[red]Use 1-3 alphanumeric characters.[/red]")
                continue
            if name[0].isdigit():
                # A leading digit forces tLEaP unit-name quoting and reads as a
                # count, not a residue — reject rather than emit a fragile name.
                self.console.print(
                    "[red]A residue name must start with a letter "
                    "(a leading digit breaks tLEaP unit references).[/red]")
                continue
            if name in reserved:
                self.console.print(
                    f"[red]{name} is already taken (standard residue or part of this "
                    f"conjugate) — tLEaP would resolve it to the wrong unit.[/red]")
                continue
            if name in library_names:
                # A library entry already owns this name. Reusing it is only right
                # when this IS that same residue (a re-parameterization); for a
                # different adduct it would shadow the existing unit at deposit.
                # Warn with the owner and require an explicit confirmation.
                self.console.print(
                    f"[yellow]{name} is already parameterized in your forcefield "
                    f"library ([grey50]{library_names[name]}[/grey50]). Reuse it only "
                    f"if this is the same residue — otherwise the new parameters will "
                    f"overwrite the existing ones.[/yellow]")
                if not confirm_with_context(
                    self.processor,
                    f"Reuse the existing library name {name}?",
                    default=False,
                    module="Force Field Parameterizer",
                    description="Confirm reuse of an existing library residue name",
                ):
                    continue
            return name

    def _build_capped_conjugate_pdb(self, residues, combined_name, redox_site=None):
        """Build an ACE/NME-capped model compound for a covalent AA↔ligand adduct.

        Unlike the de-novo route (which builds ACE-<parent>-NME from the ff14SB
        library and discards the real structure), this keeps the
        crystallographic coordinates of the whole conjugate — the modifying
        ligand comes along already attached — and caps the peptide backbone
        using the real i-1/i+1 neighbour residues, trimmed to ACE/NME so the
        junction peptide geometry stays crystallographic.

        The covalent AA↔ligand bond recorded on the RedoxSite is honoured: a
        CONECT is written for it and the H it formally replaces (e.g. a Cys
        HG) is dropped. Heavy atoms only are emitted here; hydrogens are then
        added and curated through the shared :class:`HydrogenEditor`.

        Args:
            residues: NonStandardResidue objects making up the conjugate.
            combined_name: display/base name for the combined unit.
            redox_site: RedoxSite carrying the covalent bond (falls back to
                ``residues[0].source_redox_site``).

        Returns:
            Path to the capped, hydrogen-curated PDB, or None on failure.
        """
        import numpy as np

        structure = self._get_structure_object()
        if not structure:
            self.console.print("[red]No structure available in workspace[/red]")
            return None

        model = next(iter(structure))  # first model

        # Per-chain residue maps (by sequence number) + the conjugate residues.
        targets = {(r.chain_id, int(r.resid)) for r in residues}
        chain_maps = {}
        conjugate = []  # BioPython residues belonging to the conjugate
        for chain in model:
            cmap = {}
            for res in chain:
                # standard/hetero residues only keyed by seq num (icode ignored)
                if res.id[2] not in (" ", ""):
                    self.console.print(
                        f"[yellow]⚠ Insertion code on {res.resname} {chain.id}{res.id[1]}"
                        f"{res.id[2]} — neighbour lookup ignores icodes[/yellow]"
                    )
                cmap[res.id[1]] = res
                if (chain.id, res.id[1]) in targets:
                    conjugate.append((chain.id, res))
            chain_maps[chain.id] = cmap

        if not conjugate:
            self.console.print("[red]Could not locate the conjugate residues in the structure[/red]")
            return None

        # The amino-acid partner = the conjugate residue carrying a peptide backbone.
        aa_chain = aa_res = None
        for chain_id, res in conjugate:
            if all(name in res for name in self._BACKBONE):
                aa_chain, aa_res = chain_id, res
                break
        if aa_res is None:
            self.console.print(
                "[red]No conjugate residue has a full N,CA,C,O backbone — cannot cap as a "
                "modified amino acid. (Is this really an AA↔ligand adduct?)[/red]"
            )
            return None

        aa_resid = aa_res.id[1]
        src_resname = aa_res.resname.strip()[:3]
        cmap = chain_maps[aa_chain]
        prev_res = cmap.get(aa_resid - 1)
        next_res = cmap.get(aa_resid + 1)

        # Caps need real neighbour geometry: prev C,O,CA and next N,CA.
        if prev_res is None or not all(n in prev_res for n in ("C", "O", "CA")):
            self.console.print(
                f"[red]Missing/incomplete i-1 neighbour for {src_resname} {aa_chain}{aa_resid} "
                f"— cannot build the ACE cap from real geometry (chain terminus?).[/red]"
            )
            return None
        if next_res is None or not all(n in next_res for n in ("N", "CA")):
            self.console.print(
                f"[red]Missing/incomplete i+1 neighbour for {src_resname} {aa_chain}{aa_resid} "
                f"— cannot build the NME cap from real geometry (chain terminus?).[/red]"
            )
            return None

        # The adduct is its own chemical entity, NOT the parent amino acid: a
        # Cys with an inhibitor on its SG is not CYS. Name it now, once capping
        # is known to be viable, so the working directory, every file prefix, the
        # AC file and the deposited library all carry one identity. Only the
        # returned identity changes — the model PDB below deliberately keeps the
        # SOURCE residue names so reduce can protonate components by name.
        resname_out = self._name_conjugate_residue(src_resname, aa_resid, residues)

        def _rec(name, resname, resid, coord, element):
            return {"name": name, "resname": resname, "resid": resid,
                    "x": float(coord[0]), "y": float(coord[1]), "z": float(coord[2]),
                    "element": element}

        # --- ACE cap (from i-1): carbonyl C/O retained, CA becomes the methyl C ---
        ace = [
            _rec("C", "ACE", aa_resid - 1, prev_res["C"].get_coord(), "C"),
            _rec("O", "ACE", aa_resid - 1, prev_res["O"].get_coord(), "O"),
            _rec("CH3", "ACE", aa_resid - 1, prev_res["CA"].get_coord(), "C"),
        ]

        # --- conjugate body: every heavy atom of every conjugate residue.
        #     Each atom KEEPS its source residue's own name/number here. That is
        #     what lets `reduce` protonate the model: it looks components up in
        #     its residue + HETATM-connection dictionaries by residue name (e.g.
        #     RBF -> riboflavin), so a well-known modifier is fully hydrogenated
        #     automatically; anything reduce doesn't recognise falls through to
        #     manual curation in the viewer. The name is transient — downstream
        #     QM/antechamber bin atoms positionally and stamp the single modified
        #     residue name, so the finished library is still one residue.
        #     Atom names stay unique across the whole conjugate (antechamber later
        #     treats it as one residue), and the thiol H the covalent bond
        #     replaces (Cys HG) is dropped.
        #
        #     Atoms are split into the amino-acid body and the ligand body so the
        #     file can be emitted ACE - AA - NME - ligand (see the assembly
        #     below): reduce only builds the NME cap hydrogens when the
        #     ACE-AA-NME peptide is a contiguous, uninterrupted run — a ligand
        #     residue wedged between the AA and NME makes reduce treat NME as
        #     disconnected and leave it bare. The ligand therefore trails after
        #     NME. This ordering is invisible downstream: capped_pdb_to_gaussian
        #     re-bins by residue name and re-emits ACE - middle - NME, so the
        #     .gjf/.ac/.prep that tleap consumes still have ACE first / NME last
        #     for head/tail detection. ---
        HG_NAMES = {"HG", "HG1"}  # thiol H replaced by the covalent bond

        # --- alternate conformations of the conjugate --------------------
        # A disordered conjugate (e.g. a flavin modelled in two altlocs) is
        # emitted as one capped model PER altloc letter, so each crystallographic
        # conformer gets its own QM optimisation + ESP and the modified-AA
        # checklist fits RESP charges across all of them jointly. Atoms with no
        # altloc are shared verbatim by every conformer; only the disordered
        # atoms swap coordinates, so the atom roster and order stay identical
        # across conformers (a hard requirement for multi-conformer RESP). The
        # crystal model is heavy-atom only, so hydrogens never carry an altloc
        # and are ignored during detection.
        altloc_letters = set()
        for _cid, res in conjugate:
            for atom in res:
                if not atom.is_disordered():
                    continue
                for child in atom.disordered_get_list():
                    al = (child.get_altloc() or "").strip()
                    el = (child.element or "").strip()
                    nm = child.get_name().strip()
                    is_h = el in ("H", "D") or (not el and nm[:1] == "H")
                    if al and not is_h:
                        altloc_letters.add(al)
        conformer_letters = sorted(altloc_letters)
        multi_conformer = len(conformer_letters) > 1

        def _build_bodies(altloc, warn=False):
            """Heavy-atom AA/ligand bodies for one conformer. altloc=None (or a
            letter absent on a given site) falls back to the Biopython-selected
            child, so ordered atoms are shared by every conformer."""
            aa_body, lig_body = [], []
            seen_names = set()
            for _chain_id, res in conjugate:
                res_name = res.resname.strip()[:3]
                res_id = res.id[1]
                dest = aa_body if res is aa_res else lig_body
                for atom in res:
                    if atom.is_disordered():
                        if altloc is not None and atom.disordered_has_id(altloc):
                            src = atom.disordered_get(altloc)
                        else:
                            src = atom.disordered_get_list()[0]
                    else:
                        src = atom
                    el = (src.element or "").strip()
                    nm = src.get_name().strip()
                    is_h = el in ("H", "D") or (not el and nm[:1] == "H")
                    if is_h:
                        continue  # crystal is H-less; drop any stray H (incl. Cys HG)
                    if res is aa_res and nm in HG_NAMES:
                        continue
                    if nm in seen_names:
                        if warn:
                            self.console.print(
                                f"[yellow]⚠ Duplicate atom name '{nm}' across conjugate residues — "
                                f"keeping first; rename to disambiguate if this is wrong.[/yellow]"
                            )
                        continue
                    seen_names.add(nm)
                    dest.append(_rec(nm, res_name, res_id,
                                     src.get_coord(), el or nm[:1]))
            return aa_body, lig_body

        # --- NME cap (from i+1): amide N retained, CA becomes the methyl C ---
        nme = [
            _rec("N", "NME", aa_resid + 1, next_res["N"].get_coord(), "N"),
            _rec("CH3", "NME", aa_resid + 1, next_res["CA"].get_coord(), "C"),
        ]
        assert len(ace) == 3 and len(nme) == 2, "cap heavy-atom counts drifted"

        def _assemble(aa_body, lig_body):
            # ACE - AA - NME - ligand: keeps the peptide backbone contiguous for
            # reduce; downstream canonicalises by resname so NME still ends up last.
            atoms = ace + aa_body + nme + lig_body
            for i, a in enumerate(atoms, start=1):
                a["serial"] = i
            return atoms

        # The reference conformer (first altloc letter, or the selected atoms
        # when the conjugate is ordered) defines the atom roster + serials; the
        # covalent-bond CONECT is computed once and shared, because every
        # conformer emits the same atoms in the same order.
        ref_letter = conformer_letters[0] if multi_conformer else None
        ref_aa_body, ref_lig_body = _build_bodies(ref_letter, warn=True)
        atoms = _assemble(ref_aa_body, ref_lig_body)

        # --- honour the covalent bond: CONECT between its two endpoints ---
        conect_pairs = []
        bonds = []
        rsite = redox_site or getattr(residues[0], "source_redox_site", None)
        if rsite is not None:
            bonds = [b for b in getattr(rsite, "bonds", [])
                     if getattr(b, "chemical_type", "") == "covalent"]

        def _serial_at(coord, tol=0.6):
            best, best_d = None, tol
            for a in atoms:
                d = ((a["x"] - coord[0]) ** 2 + (a["y"] - coord[1]) ** 2
                     + (a["z"] - coord[2]) ** 2) ** 0.5
                if d < best_d:
                    best, best_d = a["serial"], d
            return best

        for b in bonds:
            s1 = _serial_at(b.atom1_coords)
            s2 = _serial_at(b.atom2_coords)
            if s1 and s2 and s1 != s2:
                conect_pairs.append((s1, s2))
        if bonds and not conect_pairs:
            self.console.print(
                "[yellow]⚠ Covalent bond defined on the RedoxSite but its atoms did not match "
                "the extracted coordinates — no CONECT written; verify in the viewer.[/yellow]"
            )

        # --- write one capped heavy-atom PDB per conformer ---------------
        output_dir = Path("combined_residues")
        output_dir.mkdir(exist_ok=True)
        base = combined_name.replace('+', '_')

        def _atom_name_field(name):
            # 4-wide name field; short non-metal names get a leading space (col 13)
            return name.ljust(4) if len(name) >= 4 else (" " + name).ljust(4)

        def _write_capped_pdb(pdb_atoms, pdb_path):
            with open(pdb_path, "w") as f:
                f.write(f"REMARK  Capped conjugate model compound for {combined_name}\n")
                f.write("REMARK  Order: ACE(head) - modified residue - NME(tail); heavy atoms only\n")
                for a in pdb_atoms:
                    f.write(
                        # Cols: 13-16 atom name, 17 altLoc (blank), 18-20 resName,
                        # 22 chain. The blank between {name} and {resname} is the
                        # altLoc column — without it resName shifts left and reduce
                        # misreads it (ACE->"CE") and refuses to protonate.
                        "ATOM  {serial:5d} {name} {resname:>3s} A{resid:4d}    "
                        "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n".format(
                            serial=a["serial"], name=_atom_name_field(a["name"]),
                            resname=a["resname"], resid=a["resid"],
                            x=a["x"], y=a["y"], z=a["z"], el=a["element"][:2],
                        )
                    )
                for s1, s2 in conect_pairs:
                    f.write(f"CONECT{s1:5d}{s2:5d}\n")
                    f.write(f"CONECT{s2:5d}{s1:5d}\n")
                f.write("END\n")

        # Returns a list of (conformer_label, abspath), reference conformer first.
        conformer_pdbs = []
        n_heavy = len(ref_aa_body) + len(ref_lig_body)
        if multi_conformer:
            for letter in conformer_letters:
                aa_body, lig_body = _build_bodies(letter)
                atoms_c = _assemble(aa_body, lig_body)
                pdb_path = output_dir / f"{base}_capped_alt{letter}.pdb"
                _write_capped_pdb(atoms_c, pdb_path)
                conformer_pdbs.append((f"xtal{letter}", os.path.abspath(str(pdb_path))))
            self.console.print(
                f"[green]✓ Built {len(conformer_letters)} capped conjugate models "
                f"for alternate conformations {', '.join(conformer_letters)}:[/green]\n"
                f"  [grey50]each ACE(3) + {n_heavy} conjugate heavy atoms + NME(2); "
                f"{len(conect_pairs)} covalent CONECT bond(s); modified residue = "
                f"{resname_out}. All conformers share one atom roster; the checklist fits "
                f"RESP charges across them jointly.[/grey50]"
            )
        else:
            pdb_path = output_dir / f"{base}_capped.pdb"
            _write_capped_pdb(atoms, pdb_path)
            conformer_pdbs.append(("xtal", os.path.abspath(str(pdb_path))))
            self.console.print(
                f"[green]✓ Built capped conjugate model:[/green] {pdb_path}\n"
                f"  [grey50]ACE(3) + {n_heavy} conjugate heavy atoms + NME(2); "
                f"{len(conect_pairs)} covalent CONECT bond(s); modified residue = "
                f"{resname_out}[/grey50]"
            )

        # Hydrogens are NOT added here: the heavy-atom capped model(s) are handed
        # to the modified-AA checklist, whose Step 1 (from-structure branch) runs
        # the shared HydrogenEditor on the reference conformer and reconciles the
        # rest. This keeps every interactive step inside the checklist UI,
        # consistent with the metal / small-molecule / Route-A flows.
        return conformer_pdbs, resname_out

    def _parameterize_modaa_from_structure(self, residue_name, residues, conformer_pdbs):
        """Route B: parameterize a modified AA from capped model compound(s)
        extracted from the real structure. Bypasses the de-novo parent build.

        conformer_pdbs is a list of (conformer_label, path) tuples — one capped
        model per crystallographic alternate conformation. The reference (first)
        conformer's path is passed as ``starting_pdb`` for backward compatibility;
        the full set drives per-conformer QM and joint multi-conformer RESP.
        """
        output_dir = os.path.join(os.getcwd(), _modaa_output_dirname(residue_name))
        os.makedirs(output_dir, exist_ok=True)
        self.console.print(
            f"[green]Launching Modified Amino Acid Parameterizer "
            f"(from-structure route) for {residue_name}...[/green]"
        )
        conformer_pdbs = [(label, os.path.abspath(path)) for label, path in conformer_pdbs]
        try:
            workflow_result = modified_amino_acid_parameterizer.run_workflow(
                amino_acid=residue_name,
                output_dir=output_dir,
                interactive=True,
                processor=self.processor,
                starting_pdb=conformer_pdbs[0][1],
                conformer_pdbs=conformer_pdbs,
                conformer_mode="from_structure",
                # Original RedoxSite members, so the workflow's final FF-integration
                # step can name/rename residues + emit a reuse transformer.
                source_residues=[
                    {"name": r.name, "chain_id": r.chain_id, "resid": r.resid}
                    for r in residues
                ],
            )
        except Exception as e:
            import traceback

            traceback.print_exc()
            self.console.print(f"[red]Error during parameterization: {str(e)}[/red]")
            return

        if not workflow_result or not workflow_result.get("success"):
            # Paused-for-QM and genuine failures both land here; the workflow
            # step has already printed the specifics (e.g. the gjf to run).
            msg = (workflow_result or {}).get(
                "message", "Workflow did not complete (see messages above)."
            )
            self.console.print(f"[grey50]{msg}[/grey50]")
            return

        prep_file = workflow_result["parameter_files"].get("prep_file")
        frcmod_file = workflow_result["parameter_files"].get("frcmod_file")
        if prep_file and frcmod_file:
            if not self.get_workspace().has("parameterized_residues"):
                self.update_workspace("parameterized_residues", {})
            param_residues = self.get_from_workspace("parameterized_residues", {})
            param_residues[residue_name] = {
                "prep_file": prep_file,
                "frcmod_file": frcmod_file,
                "parent_residue": None,
                "from_structure": True,
                "success": True,
            }
            self.update_workspace("parameterized_residues", param_residues)
            # FF integration (deposit + prepared-PDB rename + workspace keys +
            # reuse transformer) now runs as the workflow's final checklist step
            # (_run_step_10_from_structure), mirroring the metal-site parameterizer.
        else:
            self.console.print(
                "[yellow]Workflow completed but no parameter files were generated[/yellow]"
            )

    def parameterize_modified_amino_acid(
        self, residue_name: str, residues: List[NonStandardResidue], combined_pdb: Optional[str] = None
    ):
        """Parameterize a modified amino acid residue

        Args:
            residue_name: Name of the residue or combined name (e.g., "L1R+ORE")
            residues: List of NonStandardResidue objects
            combined_pdb: Optional path to PDB file containing combined residues
        """
        if not residues:
            self.console.print("[red]No residues to parameterize[/red]")
            return

        # Route B: capped model compound(s) extracted from the real structure
        # were supplied. Parameterize directly instead of building ACE-<parent>-NME.
        # combined_pdb is either a single path (string, one conformer) or a list
        # of (conformer_label, path) tuples — one capped model per crystallographic
        # alternate conformation.
        if combined_pdb:
            if isinstance(combined_pdb, str):
                conformer_pdbs = [("xtal", combined_pdb)]
            else:
                conformer_pdbs = list(combined_pdb)
            if len(conformer_pdbs) > 1:
                self.console.print(
                    f"[cyan]Capped model compounds ({len(conformer_pdbs)} alternate "
                    f"conformations):[/cyan] {', '.join(p for _, p in conformer_pdbs)}")
            else:
                self.console.print(f"[cyan]Capped model compound: {conformer_pdbs[0][1]}[/cyan]")
            self.console.print(f"[cyan]Using it as the QM starting structure (from-structure route)[/cyan]")
            self._parameterize_modaa_from_structure(residue_name, residues, conformer_pdbs)
            return

        # Get parent residue from first instance
        parent_residue = residues[0].parent_residue
        ccd_data = residues[0].ccd_data

        self.console.print(
            f"\n[bold]Parameterizing Modified Amino Acid: {residue_name}[/bold]"
        )

        if ccd_data and "name" in ccd_data:
            self.console.print(f"Full Name: {ccd_data['name']}")

            # Show warning about CCD entries potentially being different compounds
            self.console.print(
                "\n[yellow]⚠ Note:[/yellow] The CCD may have an entry for this residue code that differs "
                "from your intended compound."
            )
            self.console.print(
                "[grey50]  For example, 'TYO' in the CCD is a ring-opened oxidative degradation product,[/grey50]"
            )
            self.console.print(
                "[grey50]  not a tyrosyl radical. Verify the Full Name above matches your intended molecule.[/grey50]"
            )

        # If no parent residue, prompt user to specify it
        if not parent_residue:
            self.console.print(
                f"[yellow]No parent residue found for {residue_name}.[/yellow]"
            )
            self.console.print("[yellow]Please specify the parent standard amino acid.[/yellow]")

            # List of standard amino acids for selection
            standard_aa_list = list(self.standard_aa)
            standard_aa_list.sort()

            self.console.print("\n[bold]Select parent residue:[/bold]")
            for i, aa in enumerate(standard_aa_list, 1):
                self.console.print(f"  {i}. {aa}")

            aa_choice = prompt_with_context(None,
                "Select parent amino acid number (or 'q' to cancel)",
                choices=[str(i) for i in range(1, len(standard_aa_list) + 1)] + ["q"],
            )

            if aa_choice == "q":
                self.console.print("[yellow]Parameterization cancelled[/yellow]")
                return

            parent_residue = standard_aa_list[int(aa_choice) - 1]

            # Update the residue objects with the parent
            for res in residues:
                res.parent_residue = parent_residue

            # Optionally add to mapping table
            add_to_map = confirm_with_context(
                self.processor,
                f"Add {residue_name} → {parent_residue} to mapping table for future use?",
                default=True,
                module="Force Field Parameterizer",
                description=f"Add {residue_name} to modified AA mapping table",
            )
            if add_to_map:
                self.modified_aa_map[residue_name.upper()] = parent_residue
                self._save_modified_aa_map()
                self.console.print(f"[green]✓ Added to mapping table[/green]")

        self.console.print(f"Parent Standard Residue: {parent_residue}")

        # Allow user to accept or change the parent residue
        confirm_parent = confirm_with_context(
            self.processor,
            f"Use {parent_residue} as parent residue?",
            default=True,
            module="Force Field Parameterizer",
            description=f"Use {parent_residue} as parent residue",
        )

        if not confirm_parent:
            # List of standard amino acids for selection
            standard_aa_list = list(self.standard_aa)
            standard_aa_list.sort()

            self.console.print("\n[bold]Select parent residue:[/bold]")
            for i, aa in enumerate(standard_aa_list, 1):
                self.console.print(f"  {i}. {aa}")

            aa_choice = prompt_with_context(None,
                "Select parent amino acid (or enter 'c' for custom 3-letter code)",
                choices=[str(i) for i in range(1, len(standard_aa_list) + 1)] + ["c"],
                default="c",
            )

            if aa_choice == "c":
                # Custom entry
                parent_residue = prompt_with_context(
                    self.processor,
                    "Enter 3-letter code for parent amino acid",
                    default=parent_residue,
                    module="Force Field Parameterizer",
                    description="Custom parent amino acid 3-letter code",
                ).upper()
            else:
                # Selection from list
                parent_residue = standard_aa_list[int(aa_choice) - 1]

        # Important: Clearly explain the workflow to the user
        self.console.print()
        self.console.print("[yellow]IMPORTANT WORKFLOW INFORMATION[/yellow]")
        self.console.print(
            f"The parameterization workflow will use the [bold]parent amino acid ({parent_residue})[/bold] structure as a starting point."
        )
        self.console.print(
            f"During the workflow, you will need to modify this structure to create the [bold]modified amino acid ({residue_name})[/bold]."
        )
        self.console.print(
            f"The output files will be saved in a directory named "
            f"[bold]{_modaa_output_dirname(residue_name)}[/bold]."
        )

        # Create output directory for this residue if it doesn't exist
        output_dir = os.path.join(os.getcwd(), _modaa_output_dirname(residue_name))
        os.makedirs(output_dir, exist_ok=True)

        # Change to the output directory - important for the workflow
        original_dir = os.getcwd()
        os.chdir(output_dir)

        self.console.print(f"[green]Output directory: {output_dir}[/green]")
        self.console.print(
            "[green]Launching Modified Amino Acid Parameterizer workflow...[/green]"
        )

        try:
            # Check if parameter files already exist
            existing_prep = glob.glob(f"{residue_name.lower()}.prep")
            existing_frcmod = glob.glob(f"{residue_name.lower()}.frcmod")

            if existing_prep and existing_frcmod:
                self.console.print(
                    f"[yellow]Parameter files already exist for {residue_name}[/yellow]"
                )
                overwrite = prompt_with_context(
                    self.processor,
                    "Overwrite existing parameter files?",
                    choices=["y", "n"],
                    default="n",
                    module="Force Field Parameterizer",
                    description="Overwrite existing parameter files",
                    options_map={"y": "Yes, overwrite", "n": "No, keep existing"},
                )

                if overwrite.lower() != "y":
                    self.console.print(
                        f"Using existing parameter files for {residue_name}"
                    )

                    # Store the result in the workspace
                    if not self.get_workspace().has("parameterized_residues"):
                        self.update_workspace("parameterized_residues", {})

                    param_residues = self.get_from_workspace(
                        "parameterized_residues", {}
                    )
                    param_residues[residue_name] = {
                        "prep_file": os.path.abspath(existing_prep[0]),
                        "frcmod_file": os.path.abspath(existing_frcmod[0]),
                        "parent_residue": parent_residue,
                        "reused_existing": True,
                        "success": True,
                    }
                    self.update_workspace("parameterized_residues", param_residues)

                    # Return to original directory
                    os.chdir(original_dir)
                    return

            # Run the full workflow from modified_amino_acid_parameterizer
            # Pass the PARENT amino acid name to the workflow
            workflow_result = modified_amino_acid_parameterizer.run_workflow(
                amino_acid=parent_residue,  # IMPORTANT: Pass parent residue here, not the modified name
                output_dir=output_dir,
                interactive=True,
            )

            # Process workflow results
            if workflow_result["success"]:
                # Get parameter files
                prep_file = workflow_result["parameter_files"].get("prep_file")
                frcmod_file = workflow_result["parameter_files"].get("frcmod_file")

                # Check workflow status
                if (
                    workflow_result.get("status") == "paused"
                    or workflow_result.get("status") == "pending_calculations"
                ):
                    # Workflow needs Gaussian calculations
                    missing_files = workflow_result.get("missing_files", [])
                    self.console.print(
                        f"[yellow]Parameterization requires Gaussian calculations to complete[/yellow]"
                    )
                    self.console.print(f"The following files need to be created:")
                    for file in missing_files:
                        self.console.print(f"  - {file}")

                    # Store in workspace as a pending workflow
                    if not self.get_workspace().has("pending_parameterizations"):
                        self.update_workspace("pending_parameterizations", {})

                    pending_param = self.get_from_workspace(
                        "pending_parameterizations", {}
                    )
                    pending_param[residue_name] = {
                        "output_dir": output_dir,
                        "missing_files": missing_files,
                        "status": workflow_result.get("status"),
                        "parent_residue": parent_residue,
                    }
                    self.update_workspace("pending_parameterizations", pending_param)

                    self.console.print(
                        "[yellow]The parameterization process has been paused.[/yellow]"
                    )
                    self.console.print(
                        "You will need to perform the required Gaussian calculations"
                    )
                    self.console.print(
                        "and then restart the parameterization workflow."
                    )

                    # Offer to save workspace and exit/background
                    self._offer_workspace_save_and_exit()

                elif prep_file and frcmod_file:
                    # Workflow completed successfully with parameter files
                    # Store in workspace (output already shown in modified_amino_acid_parameterizer)
                    if not self.get_workspace().has("parameterized_residues"):
                        self.update_workspace("parameterized_residues", {})

                    param_residues = self.get_from_workspace(
                        "parameterized_residues", {}
                    )
                    param_residues[residue_name] = {
                        "prep_file": prep_file,
                        "frcmod_file": frcmod_file,
                        "parent_residue": parent_residue,
                        "success": True,
                    }
                    self.update_workspace("parameterized_residues", param_residues)

                    # Library promotion happens in the workflow's final step
                    # (aa-10 "Force Field Integration") for both routes. The
                    # from-structure route already deposited there via
                    # promote_state, so promoting again here double-deposited.
                else:
                    # Workflow completed but didn't create parameter files
                    self.console.print(
                        "[yellow]Workflow completed but no parameter files were generated[/yellow]"
                    )
            else:
                # Workflow failed
                error = workflow_result.get(
                    "error", workflow_result.get("message", "Unknown error")
                )
                self.console.print(
                    f"[red]Error in parameterization workflow: {error}[/red]"
                )

            # Return to original directory
            os.chdir(original_dir)

        except ImportError:
            self.console.print(
                "[yellow]Modified Amino Acid Parameterizer module not found[/yellow]"
            )
            self.console.print(
                f"[italic]This would launch the parameterization for {residue_name} ({parent_residue})[/italic]"
            )
        except Exception as e:
            import traceback

            traceback.print_exc()
            self.console.print(f"[red]Error during parameterization: {str(e)}[/red]")

            # Make sure we return to original directory
            if "original_dir" in locals():
                os.chdir(original_dir)

    def parameterize_small_molecule(self, residue_name: str, residues: List[NonStandardResidue], combined_pdb: Optional[str] = None):
        """Parameterize a small molecule residue - parallels parameterize_modified_amino_acid method.

        Args:
            residue_name: Name of the residue or combined name (e.g., "L1R+ORE")
            residues: List of NonStandardResidue objects
            combined_pdb: Optional path to PDB file containing combined residues
        """
        if not residues:
            self.console.print("[red]No residues to parameterize[/red]")
            return

        self.console.print(f"\n[bold]Parameterizing Small Molecule: {residue_name}[/bold]")

        # Display info about combined residues if applicable
        if combined_pdb:
            self.console.print(f"[cyan]Combined residue structure saved to: {combined_pdb}[/cyan]")
            self.console.print(f"[cyan]This PDB will be used as the starting structure for parameterization[/cyan]")
        
        # Display residue information
        locations = [f"{res.chain_id}:{res.resid}" for res in residues[:3]]
        if len(residues) > 3:
            locations.append(f"... and {len(residues) - 3} more")
        self.console.print(f"Found {len(residues)} instance(s) at: {', '.join(locations)}")

        # Get the actual BioPython residue objects from the structure
        structure = self._get_structure_object()
        if not structure:
            self.console.print("[red]No structure available in workspace[/red]")
            return
        
        # Convert NonStandardResidue objects to actual BioPython residue objects
        biopython_residues = []
        for ns_residue in residues:
            # Find the actual BioPython residue in the structure
            for model in structure:
                for chain in model:
                    if chain.id == ns_residue.chain_id:
                        for residue in chain:
                            if (residue.id[1] == ns_residue.resid and 
                                residue.get_resname().strip() == ns_residue.name):
                                # Create a wrapper object that has both the BioPython residue and metadata
                                class ResidueWrapper:
                                    def __init__(self, biopython_residue, chain_id, resid):
                                        self.biopython_residue = biopython_residue
                                        self.chain_id = chain_id
                                        self.resid = resid
                                    
                                    def get_atoms(self):
                                        return self.biopython_residue.get_atoms()
                                    
                                    def get_resname(self):
                                        return self.biopython_residue.get_resname()
                                
                                wrapper = ResidueWrapper(residue, ns_residue.chain_id, ns_residue.resid)
                                biopython_residues.append(wrapper)
                                break
        
        if not biopython_residues:
            self.console.print("[red]Could not find BioPython residue objects in structure[/red]")
            return

        # Create output directory
        original_dir = os.getcwd()
        output_dir = f"small_molecule_params_{residue_name}"
        
        try:
            # Import the small molecule parameterizer
            from . import small_molecule_parameterizer
            
            self.console.print(f"[cyan]Starting small molecule parameterization workflow...[/cyan]")
            
            # Run the workflow - now passes BioPython residue objects.
            # parameterize_small_molecule is a generate-only action (no
            # "I already have parameters" branch), so regenerate rather than
            # silently reuse pre-existing .mol2/.frcmod on a re-run or replay.
            workflow_result = small_molecule_parameterizer.run_workflow(
                residue_name=residue_name,
                residues=biopython_residues,  # Pass BioPython residue objects instead
                output_dir=output_dir,
                interactive=True,
                processor=self.processor,
                regenerate=True,
            )
            
            # Process workflow results (matching the modified amino acid parameterizer pattern)
            if workflow_result["success"]:
                # Get parameter files
                mol2_file = workflow_result["parameter_files"].get("prep_file")  # MOL2 file stored as prep_file
                frcmod_file = workflow_result["parameter_files"].get("frcmod_file")

                # Check workflow status
                if (
                    workflow_result.get("status") == "paused"
                    or workflow_result.get("status") == "pending_calculations"
                ):
                    # Workflow needs Gaussian calculations
                    missing_files = workflow_result.get("missing_files", [])

                    # Store in workspace as a pending workflow
                    if not self.get_workspace().has("pending_parameterizations"):
                        self.update_workspace("pending_parameterizations", {})

                    pending_param = self.get_from_workspace(
                        "pending_parameterizations", {}
                    )
                    pending_param[residue_name] = {
                        "output_dir": output_dir,
                        "missing_files": missing_files,
                        "status": workflow_result.get("status"),
                        "type": "small_molecule",  # Mark as small molecule type
                    }
                    self.update_workspace("pending_parameterizations", pending_param)

                    # Offer to save workspace and exit/background
                    self._offer_workspace_save_and_exit()

                elif mol2_file and frcmod_file:
                    # Workflow completed successfully - store in workspace
                    # (Summary already displayed by small_molecule_parameterizer)
                    if not self.get_workspace().has("parameterized_residues"):
                        self.update_workspace("parameterized_residues", {})

                    param_residues = self.get_from_workspace(
                        "parameterized_residues", {}
                    )
                    param_residues[residue_name] = {
                        "mol2_file": mol2_file,
                        "frcmod_file": frcmod_file,
                        "type": "small_molecule",
                        "success": True,
                    }
                    self.update_workspace("parameterized_residues", param_residues)

                    # Registration, library promotion and transformer emission
                    # all happen in the workflow's final step (sm-8 "Force Field
                    # Integration"), so every entry point into run_workflow gets
                    # them — not just this one. Doing them here too would
                    # double-prompt for the deposit.
            else:
                # Workflow failed or requires action (message already displayed)
                pass

            # Return to original directory
            os.chdir(original_dir)

        except ImportError:
            self.console.print(
                "[red]Small Molecule Parameterizer module not found[/red]"
            )
            self.console.print(
                f"[italic]This would launch the small molecule parameterization for {residue_name}[/italic]"
            )
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.console.print(f"[red]Error during small molecule parameterization: {str(e)}[/red]")

            # Make sure we return to original directory
            if "original_dir" in locals():
                os.chdir(original_dir)

    def validate_metal_site_requirements(self) -> Tuple[bool, str]:
        """
        Validate that requirements are met for metal site parameterization.
        
        Returns:
            Tuple[bool, str]: (requirements_met, message)
        """
        
        # Check for structure (prefer repaired > filtered > loader structures)
        repaired_structure = self.get_from_workspace("repaired_structure")
        filtered_structure = self.get_from_workspace("filtered_structure")

        if repaired_structure:
            return True, "Repaired structure available - ready for metal site parameterization"
        elif filtered_structure:
            # Accept filtered structure but warn about completeness
            return True, "Using filtered structure - please ensure it is complete (no missing atoms/residues/hydrogens)"

        # Check for Structure Loader keys as fallback
        loader_keys = ["rcsb_structure", "local_structure", "alphafold_structure",
                       "alphafill_structure", "alphafold_homolog_structure"]
        for key in loader_keys:
            if self.get_from_workspace(key):
                return True, "Using loaded structure - strongly recommend running Structure Completeness module first"

        return False, "No structure found in workspace. Please load a structure first."
        
    def get_metal_site_status_summary(self) -> Dict[str, Any]:
        """
        Get status summary of metal site parameterizations.
        
        Returns:
            Dict containing status information for display
        """
        
        parameterized = self.get_from_workspace("parameterized_residues", {})
        
        metal_sites = {
            name: info for name, info in parameterized.items() 
            if info.get("type") == "metal_site"
        }
        
        summary = {
            "total_metal_sites": len(metal_sites),
            "completed_step1": 0,
            "ready_for_step2": 0,
            "fully_parameterized": 0,
            "details": {}
        }
        
        for name, info in metal_sites.items():
            status_detail = {
                "status": info.get("status", "unknown"),
                "step1_completed": info.get("step1_completed", False),
                "workflow_directory": info.get("workflow_directory"),
                "component_summary": info.get("metal_site_results", {}).get("component_summary", {})
            }
            
            summary["details"][name] = status_detail
            
            if info.get("step1_completed"):
                summary["completed_step1"] += 1
                
            if info.get("status") == "partial":
                summary["ready_for_step2"] += 1
            elif info.get("status") == "completed":
                summary["fully_parameterized"] += 1
        
        return summary

    def display_metal_site_status(self):
        """
        Display current metal site parameterization status.
        
        This is a new method that can be called from the main menu or help system.
        """
        
        status = self.get_metal_site_status_summary()
        
        if status["total_metal_sites"] == 0:
            self.console.print("[yellow]No metal sites have been parameterized yet[/yellow]")
            return
        
        self.console.print(f"\n[bold cyan]Metal Site Parameterization Status[/bold cyan]")
        
        # Overall summary
        table = Table(title="Summary")
        table.add_column("Metric", style="cyan")
        table.add_column("Count", justify="right", style="green")
        
        table.add_row("Total Metal Sites", str(status["total_metal_sites"]))
        table.add_row("Step 1 Completed", str(status["completed_step1"]))
        table.add_row("Ready for Step 2", str(status["ready_for_step2"]))
        table.add_row("Fully Parameterized", str(status["fully_parameterized"]))
        
        self.console.print(table)
        
        # Detailed status
        if status["details"]:
            detail_table = Table(title="Detailed Status")
            detail_table.add_column("Metal Site", style="magenta")
            detail_table.add_column("Status", style="yellow")
            detail_table.add_column("Components", style="blue")
            detail_table.add_column("Directory", style="grey50")
            
            for name, details in status["details"].items():
                status_str = details["status"].title()
                if details["step1_completed"]:
                    status_str += " ✅"
                
                # Format component summary
                comp_summary = details.get("component_summary", {})
                components = []
                if comp_summary.get("ligands", 0) > 0:
                    components.append(f"L:{comp_summary['ligands']}")
                if comp_summary.get("metals", 0) > 0:
                    components.append(f"M:{comp_summary['metals']}")
                if comp_summary.get("waters", 0) > 0:
                    components.append(f"W:{comp_summary['waters']}")
                
                comp_str = ", ".join(components) if components else "None"
                
                # Format directory path
                work_dir = details.get("workflow_directory", "")
                if work_dir:
                    work_dir = str(Path(work_dir).name)  # Just show directory name
                
                detail_table.add_row(name, status_str, comp_str, work_dir)
            
            self.console.print(detail_table)

    def parameterize_metal_site(self, residue_name, selected_residues):
        """
        Parameterize metal site residues using the metal site parameterizer.
        
        Args:
            residue_name: Name of the metal site (e.g., "FE", "HEM_FE")
            selected_residues: List of NonStandardResidue objects
        """
        try:
            # Validate requirements with more flexible structure handling
            requirements_met, message = self.validate_metal_site_requirements()
            
            if not requirements_met:
                self.console.print(f"[red]❌ Cannot proceed with metal site parameterization:[/red]")
                self.console.print(f"[red]{message}[/red]")
                return
            
            # Show structure completeness warning if not using repaired structure
            repaired_structure = self.get_from_workspace("repaired_structure")
            if not repaired_structure:
                self.console.print(f"[yellow]⚠️  {message}[/yellow]")
                self.console.print("[yellow]Metal site parameterization is sensitive to structure completeness.[/yellow]")
                
                if not confirm_with_context(
                    self.processor,
                    "Continue with current structure?",
                    default=True,
                    module="Force Field Parameterizer",
                    description="Continue despite structure issues",
                ):
                    self.console.print("[cyan]Consider running Structure Completeness module first[/cyan]")
                    return

            # Import the metal site parameterizer
            from proprep.forcefield_prep import metal_site_parameterizer
            
            self.console.print(f"[cyan]Starting metal site parameterization workflow for {residue_name}...[/cyan]")
                            
            # Convert NonStandardResidue objects to BioPython residues for the workflow
            biopython_residues = []
            structure = self._get_structure_object()
            for nsr in selected_residues:
                # Get the actual BioPython residue from the structure
                if structure:
                    for model in structure:
                        for chain in model:
                            if chain.id == nsr.chain_id:
                                for residue in chain:
                                    if residue.id[1] == nsr.resid:
                                        biopython_residues.append(residue)
                                        break
            
            if not biopython_residues:
                self.console.print("[red]❌ Could not find corresponding residues in structure[/red]")
                return
            
            # Set up output directory
            output_dir = Path("metal_site_parameterization")
            
            # Set workspace manager for the metal site parameterizer
            # This is a bit tricky since we need to pass the workspace context
            
            # Pass workspace manager to the workflow
            workflow_result = metal_site_parameterizer.run_workflow(
                residue_name=residue_name,
                residues=biopython_residues,
                output_dir=str(output_dir),
                interactive=True,
                workspace_manager=self  # Pass self as workspace manager
            )
                                        
            # Process workflow results (similar to other parameterizers)
            if workflow_result["success"]:
                self.console.print("[green]✅ Metal site parameterization completed successfully[/green]")
                
                # Get result files
                step1_results = workflow_result.get("metal_site_results", {})
                assembled_pdb = workflow_result["parameter_files"].get("prep_file")  # Step 1 PDB stored as prep_file
                
                # Check workflow status
                status = workflow_result.get("status", "completed")
                
                if status == "partial" or status == "completed":
                    if assembled_pdb:
                        self.console.print(f"[green]Structure preparation completed:[/green]")
                        self.console.print(f"  Prepared structure: {Path(assembled_pdb).name}")
                        
                        # Show component summary
                        component_summary = step1_results.get("component_summary", {})
                        if component_summary:
                            self.console.print(f"  Components processed:")
                            if component_summary.get("ligands", 0) > 0:
                                self.console.print(f"    • Ligands: {component_summary['ligands']}")
                            if component_summary.get("metals", 0) > 0:
                                self.console.print(f"    • Metal ions: {component_summary['metals']}")
                            if component_summary.get("waters", 0) > 0:
                                self.console.print(f"    • Waters: {component_summary['waters']}")
                        
                        # Store in workspace
                        if not self.get_workspace().has("parameterized_residues"):
                            self.update_workspace("parameterized_residues", {})
                        
                        param_residues = self.get_from_workspace("parameterized_residues", {})
                        param_residues[residue_name] = {
                            "prep_file": assembled_pdb,
                            "type": "metal_site",
                            "success": True,
                            "status": status,
                            "step1_completed": True,
                            "metal_site_results": step1_results,
                            "workflow_directory": step1_results.get("work_directory")
                        }
                        self.update_workspace("parameterized_residues", param_residues)
                        
                        # Show next steps
                        if status == "partial":
                            self.console.print("\n[yellow]Next Steps:[/yellow]")
                            self.console.print("  • Step 1 (Structure Preparation) completed ✅")
                            self.console.print("  • Step 2 (Metal Center Parameter Building) - Coming soon")
                            self.console.print("  • Step 3 (Force Field Integration) - Coming soon")
                            self.console.print("  • Step 4 (Simulation Setup) - Coming soon")
                            
                            self.console.print(f"\n[cyan]Current Status:[/cyan]")
                            self.console.print(f"  The structure has been prepared and is ready for metal center")
                            self.console.print(f"  parameter building when Steps 2-4 are implemented.")
                        else:
                            self.console.print("\n[green]All available workflow steps completed![/green]")
                    
                    else:
                        self.console.print("[yellow]Workflow completed but no structure file was generated[/yellow]")
                
                elif status == "failed":
                    error = workflow_result.get("message", "Unknown error")
                    self.console.print(f"[red]❌ Metal site parameterization failed: {error}[/red]")
                    
                    failed_step = workflow_result.get("failed_step")
                    if failed_step:
                        self.console.print(f"[red]Failed at Step {failed_step}[/red]")
                    
            else:
                # Workflow failed
                error = workflow_result.get("message", "Unknown error")
                self.console.print(f"[red]❌ Error in metal site parameterization workflow: {error}[/red]")
                
        except ImportError:
            self.console.print("[yellow]Metal Site Parameterizer module not found[/yellow]")
            self.console.print(f"[italic]This would launch the metal site parameterization for {residue_name}[/italic]")
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.console.print(f"[red]Error during metal site parameterization: {str(e)}[/red]")

    def _parameterize_metal_sites(self, metal_units, all_units=None):
        """Parameterize exactly the selected metal sites, in one MCPB pass.

        MCPB numbers M*/Y* atom types across sites, so the selected sites go
        through a single checklist run rather than one run each. Only the
        selected sites are parameterized: the site ids are recorded in the
        workspace and mcpb-1 filters its redox-site list to them.

        Args:
            metal_units: The selected metal-site units.
            all_units: Every displayed unit, used to name the metal sites that
                were NOT selected.
        """
        def _site_label(unit):
            metals = unit.get("metals") or unit["members"]
            label = ", ".join(f"{m.name} ({m.chain_id}:{m.resid})" for m in metals)
            ligands = unit.get("ligands") or []
            if ligands:
                label += " + ligand " + ", ".join(
                    f"{l.name} ({l.chain_id}:{l.resid})" for l in ligands)
            return label

        n = len(metal_units)
        self.console.print(
            f"\n[bold cyan]Metal Site Parameterization: "
            f"{n} site{'s' if n != 1 else ''} selected[/bold cyan]")
        for unit in metal_units:
            self.console.print(f"  [green]will run:[/green] {_site_label(unit)}")

        # Name the metal sites left out. They are not an error — reusing one
        # site's parameters on an equivalent site is the point of the reuse
        # transformer mcpb-4 emits — but the topology build needs every metal
        # site to have parameters from somewhere.
        skipped = [u for u in (all_units or [])
                   if u["category"] in ("metal_site", "metal_ion")
                   and u not in metal_units]
        if skipped:
            for unit in skipped:
                self.console.print(f"  [grey50]not run: [/grey50] {_site_label(unit)}")
            self.console.print(
                "\n[yellow]Every metal site needs parameters before the topology "
                "build will succeed.[/yellow]")
            self.console.print(
                "[grey50]  For a site equivalent to one being run, apply the reuse "
                "transformer emitted at the end of this run instead of "
                "re-deriving it.[/grey50]")

        # Record the selection so mcpb-1 parameterizes these sites and no others.
        selected_ids = [u.get("site_id") for u in metal_units if u.get("site_id")]
        self.update_workspace("mcpb_selected_site_ids", selected_ids)

        # The seed residue names the output directory and anchors the run.
        seed_unit = metal_units[0]
        metal_site_residue = (seed_unit.get("metals") or seed_unit["members"])[0]
        self._parameterize_metal_site(metal_site_residue, announce=False)

    def _parameterize_metal_site(self, metal_site_residue, announce: bool = True):
        """
        Handle parameterization of an individual metal site.

        This function routes to the metal_site_parameterizer workflow which:
        1. Uses the source RedoxSite (if available) or finds matching one from workspace
        2. Runs MCPB atom typing and fingerprint generation
        3. (Future) Runs RESP charge fitting and parameter building

        Args:
            metal_site_residue: The seed residue; names the output directory.
            announce: False when _parameterize_metal_sites has already listed
                every selected site, so the seed is not re-announced as though
                it were the only one.
        """

        if announce:
            self.console.print(f"\n[bold cyan]Parameterizing Metal Site: {metal_site_residue.name}[/bold cyan]")
            self.console.print(f"[grey50]Location: {metal_site_residue.chain_id}:{metal_site_residue.resid}[/grey50]")

        # Check if we have a direct RedoxSite reference (preferred)
        source_redox_site = getattr(metal_site_residue, 'source_redox_site', None)

        if source_redox_site:
            # Only when this is the whole story. On a multi-site selection the
            # caller has already listed every site, and describing the seed's
            # RedoxSite here made it look like the only one being run.
            if announce:
                self.console.print(f"[green]✅ Using RedoxSite: {source_redox_site.site_id}[/green]")

                # Show what's included in this RedoxSite
                if hasattr(source_redox_site, 'centers') and source_redox_site.centers:
                    self.console.print("[cyan]This RedoxSite includes:[/cyan]")
                    for center in source_redox_site.centers:
                        marker = "→" if (center.chain == metal_site_residue.chain_id and
                                         center.resid == metal_site_residue.resid) else " "
                        self.console.print(f"  {marker} {center.resname} ({center.chain}:{center.resid})")
        else:
            # Fallback: try to find matching RedoxSite from workspace
            detected_redox_sites = self.get_from_workspace("detected_redox_sites")
            if detected_redox_sites:
                # Find the RedoxSite that contains this residue
                for site in detected_redox_sites:
                    for center in site.centers:
                        if (center.chain == metal_site_residue.chain_id and
                            center.resid == metal_site_residue.resid and
                            center.resname == metal_site_residue.name):
                            source_redox_site = site
                            break
                    if source_redox_site:
                        break

                if source_redox_site:
                    self.console.print(f"[green]✅ Found matching RedoxSite: {source_redox_site.site_id}[/green]")
                else:
                    self.console.print(f"[yellow]⚠ No matching RedoxSite found for {metal_site_residue.name}[/yellow]")
                    self.console.print("[grey50]Redox sites will be detected during preprocessing (Step 0f)[/grey50]")
            else:
                self.console.print("[grey50]No redox sites detected yet -- will be detected during preprocessing (Step 0f)[/grey50]")

        # Confirm before proceeding
        if not confirm_with_context(
            self.processor,
            "\nProceed with metal site parameterization?",
            default=True,
            module="Force Field Parameterizer",
            description="Proceed with metal site parameterization",
        ):
            self.console.print("[yellow]Metal site parameterization cancelled[/yellow]")
            return

        # Create output directory with unique identifier
        site_name = metal_site_residue.name.replace(" ", "_").replace(":", "_").replace("(", "").replace(")", "")
        chain_resid = f"{metal_site_residue.chain_id}_{metal_site_residue.resid}"
        output_dir = f"metal_site_params_{site_name}_{chain_resid}"
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        self.console.print(f"[cyan]Output directory: {output_dir}[/cyan]")

        # === STRUCTURE PREPROCESSING ===
        # Run the preprocessing workflow before metal site parameterization
        # This handles: triage, H addition, exclusion, recombination, redox detection, atom typing
        self.console.print("\n[bold cyan]═══ Structure Preprocessing ═══[/bold cyan]")

        try:
            from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor
        except ImportError:
            self.console.print("[red]Structure Preprocessor not available[/red]")
            return

        # Get the PDB file to preprocess
        pdb_file = self._get_structure_file_for_preprocessing()
        if not pdb_file:
            self.console.print("[red]No structure file available for preprocessing[/red]")
            return

        # Run preprocessing with checklist interface
        # This provides visual progress tracking, state persistence, and
        # the ability to resume after Gaussian calculations
        preprocessor = StructurePreprocessor(self.processor, self.console)
        preprocessing_result = preprocessor.run_checklist(
            pdb_file=pdb_file,
            output_dir=output_path,
            interactive=True
        )

        # The checklist handles the entire workflow (Steps 0a-0g preprocessing + Steps 1-4 parameterization)
        # It provides visual progress tracking, state persistence for Gaussian checkpoints,
        # and the ability to resume workflows later
        if preprocessing_result.success:
            self.console.print("\n[green]Workflow completed successfully![/green]")
            if preprocessing_result.prepared_pdb:
                self.console.print(f"[grey50]Final structure: {preprocessing_result.prepared_pdb}[/grey50]")
        else:
            self.console.print(f"[yellow]Workflow incomplete or cancelled[/yellow]")
            if preprocessing_result.error_message:
                self.console.print(f"[grey50]{preprocessing_result.error_message}[/grey50]")

    def _store_completed_parameterization(self, residue_name: str, workflow_result: dict, param_type: str):
        """Store completed parameterization results in workspace"""
        if not self.get_workspace().has("parameterized_residues"):
            self.update_workspace("parameterized_residues", {})
        
        param_residues = self.get_from_workspace("parameterized_residues", {})
        param_residues[residue_name] = {
            "type": param_type,
            "success": True,
            "status": "completed",
            "output_dir": workflow_result.get("output_dir"),
            "output_files": workflow_result.get("output_files", {}),
            "completion_time": datetime.now().isoformat()
        }
        
        # Add metal site specific data
        if param_type == "metal_site":
            param_residues[residue_name].update({
                "completed_steps": workflow_result.get("completed_steps", []),
                "step_results": workflow_result.get("step_results", {}),
                "metal_site_results": workflow_result.get("metal_site_results", {})
            })
        
        self.update_workspace("parameterized_residues", param_residues)

    def _store_partial_parameterization(self, residue_name: str, workflow_result: dict, param_type: str):
        """Store partial parameterization results in workspace"""
        if not self.get_workspace().has("parameterized_residues"):
            self.update_workspace("parameterized_residues", {})
        
        param_residues = self.get_from_workspace("parameterized_residues", {})
        param_residues[residue_name] = {
            "type": param_type,
            "success": True,
            "status": "partial",
            "completed_steps": workflow_result.get("completed_steps", []),
            "output_dir": workflow_result.get("output_dir"),
            "step_results": workflow_result.get("step_results", {}),
            "last_updated": datetime.now().isoformat()
        }
        
        # Add metal site specific data
        if param_type == "metal_site":
            param_residues[residue_name].update({
                "workflow_metadata": workflow_result.get("workflow_metadata", {}),
                "metal_site_results": workflow_result.get("metal_site_results", {})
            })
        
        self.update_workspace("parameterized_residues", param_residues)
        
    def _display_metal_site_details(self, metal_site):
        """Display detailed information about the metal site"""
        
        # Extract metal information
        metal_element = self._get_site_attribute(metal_site, "metal_element", "?")
        metal_chain = self._get_site_attribute(metal_site, "metal_chain", "?")
        metal_resid = self._get_site_attribute(metal_site, "metal_resid", "?")
        metal_resname = self._get_site_attribute(metal_site, "metal_resname", "?")
        coordination_number = self._get_site_attribute(metal_site, "coordination_number", "?")
        
        # Create detailed information panel
        details = []
        details.append(f"[cyan]Metal Center:[/cyan] {metal_element} in {metal_resname} {metal_chain}:{metal_resid}")
        details.append(f"[cyan]Coordination Number:[/cyan] {coordination_number}")
        
        # Display ligands
        ligands = self._get_site_attribute(metal_site, "ligands", [])
        if ligands:
            details.append(f"[cyan]Coordinating Ligands:[/cyan]")
            for i, ligand in enumerate(ligands, 1):
                if isinstance(ligand, dict):
                    chain = ligand.get("chain", "?")
                    resname = ligand.get("resname", "?")
                    resid = ligand.get("resid", "?")
                    atom_name = ligand.get("atom_name", "?")
                    distance = ligand.get("distance", "?")
                    details.append(f"  {i}. {resname} {chain}:{resid} {atom_name} ({distance:.2f} Å)")
        
        # Display additional properties
        is_heme = self._get_site_attribute(metal_site, "is_heme", False)
        if is_heme:
            details.append(f"[yellow]Special Properties:[/yellow] Heme-containing site")
        
        self.console.print(Panel("\n".join(details), title="Metal Site Details", border_style="blue"))

    def _display_metal_site_parameterization_info(self):
        """Display information about metal site parameterization workflow"""
        
        info_text = """
    [bold cyan]Metal Site Parameterization Workflow[/bold cyan]

    The metal site parameterizer handles comprehensive parameterization of metal-containing systems:

    [bold]Step 1:[/bold] Structure preparation and validation
    • Validates metal coordination geometry
    • Ensures complete structure (no missing atoms)
    • Identifies coordinating ligands and water molecules

    [bold]Step 2:[/bold] Metal center parameter building
    • Uses AmberTools Metal Center Parameter Building methodology
    • Generates metal ion parameters with proper charges
    • Handles different oxidation states and coordination environments

    [bold]Step 3:[/bold] Ligand parameterization
    • Parameterizes coordinating ligands using appropriate methods
    • Integrates with existing small molecule parameterizer for non-standard ligands
    • Handles protein residues (His, Cys, Met, etc.) coordination

    [bold]Step 4:[/bold] Force field integration
    • Generates complete parameter sets for AMBER
    • Creates bond definitions for TLEaP input
    • Ensures compatibility with standard AMBER force fields

    [bold cyan]Note:[/bold cyan] Metal site parameterization requires a repaired structure with no missing atoms due to the sensitivity of metal coordination geometry.
    """
        
        self.console.print(Panel(info_text, title="Metal Site Parameterization Information", expand=False))

    def display_parameterization_status(self):
        """Display status of all parameterization workflows in workspace"""
        self.console.print("\n[bold]Parameterization Status Summary[/bold]")
        
        # Check for completed parameterizations
        parameterized_residues = self.get_from_workspace("parameterized_residues", {})
        pending_parameterizations = self.get_from_workspace("pending_parameterizations", {})
        
        if not parameterized_residues and not pending_parameterizations:
            self.console.print("[yellow]No parameterizations found in workspace[/yellow]")
            return
        
        # Create status table
        table = Table(title="Parameterization Status")
        table.add_column("Residue", style="cyan")
        table.add_column("Type", style="magenta")
        table.add_column("Status", style="yellow")
        table.add_column("Details", style="green")
        
        # Add completed parameterizations
        for residue_name, data in parameterized_residues.items():
            param_type = data.get("type", "unknown")
            status = data.get("status", "completed") if data.get("success", False) else "failed"
            
            # Create details based on type and status
            details = ""
            if param_type == "metal_site":
                if status == "completed":
                    details = "✅ Complete workflow"
                elif status == "partial":
                    completed_steps = data.get("completed_steps", [])
                    details = f"⚠️ Steps {', '.join(map(str, completed_steps))} done"
                else:
                    details = "❌ Failed"
            elif param_type == "modified_amino_acid":
                if data.get("prep_file") and data.get("frcmod_file"):
                    details = "✅ PREP + FRCMOD files"
                else:
                    details = "⚠️ Incomplete"
            elif param_type == "small_molecule":
                if data.get("mol2_file") and data.get("frcmod_file"):
                    details = "✅ MOL2 + FRCMOD files"
                else:
                    details = "⚠️ Incomplete"
            
            table.add_row(residue_name, param_type.replace("_", " ").title(), status.title(), details)
        
        # Add pending parameterizations
        for residue_name, data in pending_parameterizations.items():
            param_type = data.get("type", "unknown")
            status = "Pending"
            missing_files = data.get("missing_files", [])
            details = f"⏳ Needs {len(missing_files)} calculation(s)"
            
            table.add_row(residue_name, param_type.replace("_", " ").title(), status, details)
        
        self.console.print(table)
        
        # Show actionable items
        if pending_parameterizations:
            self.console.print("\n[yellow]Action Required:[/yellow]")
            self.console.print("• Complete pending Gaussian calculations")
            self.console.print("• Use 'Resume workflow' option to continue parameterizations")

    def resume_pending_parameterization(self):
        """Resume a pending parameterization workflow"""
        pending_parameterizations = self.get_from_workspace("pending_parameterizations", {})
        
        if not pending_parameterizations:
            self.console.print("[yellow]No pending parameterizations found[/yellow]")
            return
        
        # Show pending workflows
        self.console.print("\n[bold]Pending Parameterizations:[/bold]")
        table = Table()
        table.add_column("No.", style="cyan")
        table.add_column("Residue", style="magenta")
        table.add_column("Type", style="yellow")
        table.add_column("Status", style="green")
        
        pending_list = list(pending_parameterizations.items())
        for i, (residue_name, data) in enumerate(pending_list, 1):
            param_type = data.get("type", "unknown")
            missing_files = data.get("missing_files", [])
            status = f"Needs {len(missing_files)} calculation(s)"
            
            table.add_row(str(i), residue_name, param_type.replace("_", " ").title(), status)
        
        self.console.print(table)
        
        # Get user selection
        choice = prompt_with_context(None,
            "Select parameterization to resume (or 'q' to cancel)",
            choices=[str(i) for i in range(1, len(pending_list) + 1)] + ["q"],
            default="q"
        )
        
        if choice == "q":
            self.console.print("[yellow]Resume cancelled[/yellow]")
            return
        
        # Resume selected workflow
        selected_idx = int(choice) - 1
        residue_name, data = pending_list[selected_idx]
        
        self.console.print(f"\n[cyan]Resuming parameterization for {residue_name}...[/cyan]")
        
        # Route to appropriate resume function based on type
        param_type = data.get("type", "unknown")
        if param_type == "modified_amino_acid":
            self._resume_modified_amino_acid_workflow(residue_name, data)
        elif param_type == "small_molecule":
            self._resume_small_molecule_workflow(residue_name, data)
        elif param_type == "metal_site":
            self._resume_metal_site_workflow(residue_name, data)
        else:
            self.console.print(f"[red]Unknown parameterization type: {param_type}[/red]")

    def _get_site_attribute(self, site, attribute, default=None):
        """Read an attribute off a site that may be an object or a dict.

        Defined on PDBProcessor, and called here seven times without existing
        on this class -- every call an AttributeError waiting for the metal-site
        summary to render. Delegates when the processor has it, so the two
        cannot drift, and falls back to the same behaviour when it does not.
        """
        processor_impl = getattr(self.processor, "_get_site_attribute", None)
        if callable(processor_impl):
            return processor_impl(site, attribute, default)

        if isinstance(site, dict):
            return site.get(attribute, default)
        return getattr(site, attribute, default)

    def _resume_modified_amino_acid_workflow(self, residue_name: str,
                                             workflow_data: dict):
        """Resume a paused modified-amino-acid parameterization.

        This method did not exist while being dispatched to, so choosing
        "Resume pending workflow" for a modAA entry raised AttributeError.
        ``resume_paused_workflow`` is the real entry point: it checks whether
        the awaited QM logs have appeared and re-enters the workflow if so.
        """
        output_dir = workflow_data.get("output_dir")
        if not output_dir:
            self.console.print("[red]No output directory recorded for this workflow[/red]")
            return

        try:
            from proprep.forcefield_prep.modified_amino_acid_parameterizer import (
                resume_paused_workflow,
            )
        except ImportError as e:
            self.console.print(
                f"[red]Modified amino acid parameterizer not available: {e}[/red]")
            return

        self.console.print(
            f"[cyan]Resuming modified amino acid workflow in {output_dir}...[/cyan]")
        try:
            result = resume_paused_workflow(residue_name, output_dir)
        except Exception as e:  # noqa: BLE001 - keep the menu alive
            logger.debug("modAA resume failed: %s", e)
            self.console.print(f"[red]Error resuming workflow: {e}[/red]")
            return

        self._record_resume_result(residue_name, result, "modified_amino_acid")

    def _resume_small_molecule_workflow(self, residue_name: str,
                                        workflow_data: dict):
        """Resume a paused small-molecule parameterization.

        This method did not exist while being dispatched to, so choosing
        "Resume pending workflow" for a small molecule -- the case that pauses
        for Gaussian and is therefore the most likely to be resumed -- raised
        AttributeError.

        The parameterizer is checklist-driven and keeps its own state in the
        output directory, so re-entering ``run_workflow`` there resumes rather
        than restarts. ``regenerate`` is left False so completed steps are
        reused; that is what resuming means.
        """
        output_dir = workflow_data.get("output_dir")
        if not output_dir:
            self.console.print("[red]No output directory recorded for this workflow[/red]")
            return

        try:
            from proprep.forcefield_prep.small_molecule_parameterizer import (
                run_workflow,
            )
        except ImportError as e:
            self.console.print(
                f"[red]Small molecule parameterizer not available: {e}[/red]")
            return

        self.console.print(
            f"[cyan]Resuming small molecule workflow in {output_dir}...[/cyan]")
        try:
            result = run_workflow(
                residue_name=residue_name,
                residues=workflow_data.get("residues", []) or [],
                output_dir=output_dir,
                interactive=True,
                processor=self.processor,
            )
        except Exception as e:  # noqa: BLE001 - keep the menu alive
            logger.debug("small molecule resume failed: %s", e)
            self.console.print(f"[red]Error resuming workflow: {e}[/red]")
            return

        self._record_resume_result(residue_name, result, "small_molecule")

    def _record_resume_result(self, residue_name: str, result, param_type: str):
        """Report a resume outcome and clear the entry once it completes."""
        if not isinstance(result, dict):
            self.console.print("[yellow]Workflow returned no result[/yellow]")
            return

        if result.get("success"):
            self.console.print("[green]Workflow resumed successfully[/green]")
            self._update_parameterization_results(residue_name, result, param_type)
            if result.get("status") == "completed":
                pending = self.get_from_workspace("pending_parameterizations", {})
                if residue_name in pending:
                    del pending[residue_name]
                    self.update_workspace("pending_parameterizations", pending)
        else:
            message = result.get("message", "Resume failed")
            self.console.print(f"[yellow]{message}[/yellow]")
            for missing in result.get("missing_files", []) or []:
                self.console.print(f"  [red]missing:[/red] {missing}")

    def _resume_metal_site_workflow(self, residue_name: str, workflow_data: dict):
        """Point at where a metal site actually resumes.

        This used to import ``proprep.ff_prep.metal_site_parameterizer``, a
        path that predates the rename to ``forcefield_prep`` and no longer
        exists -- along with MetalSiteParameterizationWorkflow and
        get_workflow_status, which exist nowhere. The ImportError was caught
        and reported as "parameterizer not available", so the branch has been
        dead rather than working.

        A metal site resumes through the parameterization checklist, which
        keeps its own workflow_state.json in the site directory and offers to
        resume on entry. Say so instead of re-implementing it.
        """
        output_dir = workflow_data.get("output_dir")
        if not output_dir:
            self.console.print("[red]No output directory recorded for this workflow[/red]")
            return

        self.console.print(
            "[cyan]Metal site parameterization resumes through its checklist."
            "[/cyan]")
        try:
            state_file = Path(output_dir) / "workflow_state.json"
            if state_file.exists():
                self.console.print(f"[grey50]  Saved state: {state_file}[/grey50]")
            else:
                self.console.print(
                    f"[grey50]  No saved state found in {output_dir}[/grey50]")
        except Exception as e:  # noqa: BLE001 - a path problem must not crash the menu
            logger.debug("Could not check metal-site state file: %s", e)

        self.console.print(
            "[grey50]  Choose 'Parameterize residues (new or resume)' from this "
            "menu; the checklist offers to resume from that state.[/grey50]")

    def _update_parameterization_results(self, residue_name: str, results: dict, param_type: str):
        """Update workspace with parameterization results"""
        if not self.get_workspace().has("parameterized_residues"):
            self.update_workspace("parameterized_residues", {})
        
        param_residues = self.get_from_workspace("parameterized_residues", {})
        
        # Update or create entry
        if residue_name not in param_residues:
            param_residues[residue_name] = {"type": param_type}
        
        # Merge results
        param_residues[residue_name].update({
            "success": results.get("success", False),
            "status": results.get("status", "unknown"),
            "last_updated": datetime.now().isoformat(),
            "output_files": results.get("output_files", [])
        })
        
        # Add type-specific data
        if param_type == "metal_site":
            param_residues[residue_name].update({
                "completed_steps": results.get("completed_steps", []),
                "step_results": results.get("step_results", {})
            })
        
        self.update_workspace("parameterized_residues", param_residues)

    def validate_workspace_for_parameterization(self):
        """Validate that workspace has required data for parameterization"""
        issues = []

        # Check for structure (workflow-processed or Structure Loader)
        structure_available = False
        struct_keys = [
            "repaired_structure", "filtered_structure",
            "rcsb_structure", "local_structure", "alphafold_structure",
            "alphafill_structure", "alphafold_homolog_structure"
        ]
        for struct_key in struct_keys:
            if self.get_from_workspace(struct_key):
                structure_available = True
                break

        if not structure_available:
            issues.append("No structure available in workspace")
        
        # Check for filter selections
        filter_selections = self.get_from_workspace("filter_selections")
        if not filter_selections:
            issues.append("No filter selections available - run PDB Filter first")
        
        # Check for metal sites if parameterizing metal sites
        metal_sites = self.get_from_workspace("metal_sites", [])
        
        return {
            "valid": len(issues) == 0,
            "issues": issues,
            "has_metal_sites": len(metal_sites) > 0,
            "structure_type": self._get_best_available_structure_type()
        }

    def _get_best_available_structure_type(self):
        """Get the best available structure type in workspace"""
        for struct_key, desc in [
            ("repaired_structure", "repaired"),
            ("filtered_structure", "filtered"),
            # Structure Loader keys as fallback
            ("rcsb_structure", "RCSB"),
            ("local_structure", "local"),
            ("alphafold_structure", "AlphaFold"),
            ("alphafill_structure", "AlphaFill"),
            ("alphafold_homolog_structure", "AlphaFold homolog"),
        ]:
            if self.get_from_workspace(struct_key):
                return desc
        return None

    def process(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Process the workspace"""
        # Convert to a local workspace reference
        workspace_obj = workspace

        # Check if we have filter selections to work with
        if not self.get_from_workspace_obj(workspace_obj, "filter_selections"):
            return workspace_obj

        # Update non_standard_residues in workspace if not already present
        if not self.get_from_workspace_obj(workspace_obj, "non_standard_residues"):
            # For automatic processing without user input
            if self.processor:
                # Store the processor's workspace for context
                # The processor's workspace is always available
                self.analyze_nonstandard_residues()

                # Update workspace with any changes made during analysis
                workspace_obj = self.update_workspace_obj(
                    workspace_obj, "non_standard_residues", self.non_standard_residues
                )

        return workspace_obj

    def get_workspace_display(self, workspace_key, value, console):
        """
        Custom display method for Force Field Parameterizer data

        Args:
            workspace_key: The key in the workspace dictionary
            value: The value stored in the workspace
            console: The rich console object for output

        Returns:
            bool: True if the module handled the display, False otherwise
        """
        if workspace_key == "non_standard_residues" and value:
            # Display non-standard residues with a custom format
            console.print(Panel("[bold]Non-Standard Residues[/bold]", style="green"))

            # Group by type for cleaner display
            by_category = defaultdict(list)
            for res in value:
                by_category[res.category].append(res)

            # Display each category
            for category, residues in by_category.items():
                category_display = category.replace("_", " ").capitalize()
                console.print(f"[bold]{category_display} ({len(residues)}):[/bold]")

                # Group by residue name
                by_name = defaultdict(list)
                for res in residues:
                    by_name[res.name].append(res)

                # Display each residue type
                for name, instances in by_name.items():
                    locations = [f"{res.chain_id}:{res.resid}" for res in instances[:3]]
                    location_str = ", ".join(locations)
                    if len(instances) > 3:
                        location_str += f", ... ({len(instances)-3} more)"

                    parent_info = ""
                    if (
                        category == "modified_amino_acid"
                        and instances[0].parent_residue
                    ):
                        parent_info = f" (derived from {instances[0].parent_residue})"

                    console.print(f"  {name}{parent_info}: {location_str}")

            return True

        elif workspace_key == "parameterized_residues" and value:
            # Display parameterized residues
            console.print(Panel("[bold]Parameterized Residues[/bold]", style="green"))

            # Create table for overview
            table = Table()
            table.add_column("Residue", style="magenta")
            table.add_column("Type", style="cyan")
            table.add_column("Status", style="green")
            table.add_column("Notes", style="blue")

            for res_name, res_data in value.items():
                status = (
                    "[green]Success[/green]"
                    if res_data.get("success", False)
                    else "[red]Failed[/red]"
                )
                res_type = (
                    "Modified amino acid"
                    if res_data.get("parent_residue")
                    else "Unknown"
                )

                notes = ""
                if "parent_residue" in res_data:
                    notes = f"Derived from {res_data['parent_residue']}"
                if res_data.get("reused_existing", False):
                    notes += " (used existing parameter files)"

                table.add_row(res_name, res_type, status, notes)

            console.print(table)

            # Display tleap commands
            console.print("\nTo use these parameterized residues in tleap:")
            for res_name, res_data in value.items():
                if "prep_file" in res_data:
                    console.print(
                        f'  loadAmberPrep "{os.path.basename(res_data["prep_file"])}"'
                    )
                if "frcmod_file" in res_data:
                    console.print(
                        f'  loadAmberParams "{os.path.basename(res_data["frcmod_file"])}"'
                    )

            return True

        return False

    # Enhanced startup validation method
    def initialize_module(self):
        """Initialize the module and validate workspace"""
        self.console.print("[cyan]Initializing Force Field Parameterizer...[/cyan]")
        
        # Validate workspace
        validation = self.validate_workspace_for_parameterization()
        
        if not validation["valid"]:
            self.console.print("[yellow]⚠️  Workspace validation issues found:[/yellow]")
            for issue in validation["issues"]:
                self.console.print(f"  • {issue}")
            self.console.print("\n[grey50]Some functionality may be limited until these issues are resolved.[/grey50]")
        else:
            self.console.print("[green]✅ Workspace validation passed[/green]")
            
            # Show helpful context
            if validation["has_metal_sites"]:
                metal_sites = self.get_from_workspace("metal_sites", [])
                self.console.print(f"[cyan]Found {len(metal_sites)} metal site(s) for parameterization[/cyan]")
            
            structure_type = validation["structure_type"]
            if structure_type:
                self.console.print(f"[cyan]Using {structure_type} structure[/cyan]")
                
                if structure_type == "original":
                    self.console.print("[yellow]💡 Consider running Structure Completeness for better parameterization results[/yellow]")

    # Enhanced display method for better user experience
    def display_enhanced_summary(self):
        """Display enhanced summary with actionable information"""
        if not self.non_standard_residues:
            self.console.print("[yellow]No non-standard residues analyzed yet.[/yellow]")
            self.console.print("[cyan]💡 Use option 1 to analyze your structure[/cyan]")
            return
        
        # Standard summary
        self.display_analysis_summary()
        
        # Add actionable guidance
        self.console.print("\n[bold cyan]Quick Actions:[/bold cyan]")
        
        # Count PARAMETERIZATION UNITS, not residues: a metal site stamps the
        # category onto its ligand members too, so counting residues overcounts
        # sites. Option numbers track get_menu_options (1 analyze, 2 import,
        # 3 parameterize/resume, 4 status, 5 help) — they had drifted, and the
        # resume line pointed at an option 6 that does not exist.
        def _units(category):
            return len({self._unit_key(res) for res in self.non_standard_residues
                        if res.category == category})

        metal_sites = _units("metal_site")
        modified_aas = _units("modified_amino_acid")
        small_molecules = _units("small_molecule")
        unknown = _units("unknown")

        if metal_sites > 0:
            self.console.print(f"[green]• Use option 3 to parameterize any of {metal_sites} metal site(s)[/green]")

        if modified_aas > 0:
            self.console.print(f"[green]• Use option 3 to parameterize any of {modified_aas} modified amino acid(s)[/green]")

        if small_molecules > 0:
            self.console.print(f"[green]• Use option 3 to parameterize any of {small_molecules} small molecule(s)[/green]")

        if unknown > 0:
            self.console.print(f"[yellow]• Use option 1 to classify {unknown} unknown residue(s)[/yellow]")

        # Check for pending work
        pending = self.get_from_workspace("pending_parameterizations", {})
        if pending:
            self.console.print(f"[yellow]• Use option 3 to resume {len(pending)} pending parameterization(s)[/yellow]")

    # Enhanced error handling wrapper
    def safe_execute_with_context(self, operation_name: str, operation_func, *args, **kwargs):
        """Execute an operation with enhanced error handling and user context"""
        try:
            self.console.print(f"[cyan]Starting {operation_name}...[/cyan]")
            
            # Pre-operation validation
            validation = self.validate_workspace_for_parameterization()
            critical_issues = [issue for issue in validation["issues"] if "structure" in issue.lower()]
            
            if critical_issues:
                self.console.print(f"[red]Cannot proceed with {operation_name}:[/red]")
                for issue in critical_issues:
                    self.console.print(f"  • {issue}")
                return False
            
            # Execute operation
            result = operation_func(*args, **kwargs)
            
            if result:
                self.console.print(f"[green]✅ {operation_name} completed successfully[/green]")
            else:
                self.console.print(f"[yellow]⚠️  {operation_name} completed with warnings[/yellow]")
                
            return result
            
        except KeyboardInterrupt:
            self.console.print(f"\n[yellow]⚠️  {operation_name} interrupted by user[/yellow]")
            return False
            
        except Exception as e:
            self.console.print(f"[red]❌ Error during {operation_name}: {str(e)}[/red]")
            
            # Show helpful troubleshooting
            if "import" in str(e).lower():
                self.console.print("[yellow]💡 This may be a missing dependency issue[/yellow]")
            elif "file" in str(e).lower() or "path" in str(e).lower():
                self.console.print("[yellow]💡 This may be a file path or permissions issue[/yellow]")
            elif "workspace" in str(e).lower():
                self.console.print("[yellow]💡 Try reloading your structure or running earlier preparation steps[/yellow]")
                
            return False
        
def main():
    """Run the module as a standalone script"""
    console = Console()
    console.print("[bold]Force Field Parameterizer[/bold]")
    console.print(
        "[italic]This module is designed to run within the PDB Processor framework[/italic]"
    )
    console.print("For standalone usage, please use the PDB Processor main interface.")

if __name__ == "__main__":
    main()
