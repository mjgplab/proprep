"""
Structure Preprocessor for Metal Site Parameterization

This module orchestrates the preprocessing workflow for metal site parameterization:

Step 0a: Structure Filtering (optional)
         - Launches PDB Filter to select chains, remove unwanted residues

Step 0b: Structure Triage + Force Field Selection
         - Categorizes residues: A (protein), B (non-standard), C (water), D (metal)
         - Selects force fields for each component type
         - Stores selections in workspace for use by later steps

Step 0c: Hydrogen Addition
         - Protein (Category A) → protonation_state_analyzer + tLEaP
         - Non-standard (Category B) → reduce

Step 0d: Atom Exclusion + H Capping
         - User specifies atoms/residues to exclude
         - Cut bonds are capped with H atoms at 1.09Å

Step 0e: Structure Recombination
         - Merges all H-added pieces into prepared_structure.pdb

Step 0f: Redox Site Sync
         - Updates existing RedoxSite objects to match prepared structure

Step 0g: Pure Atom Typing
         - Uses FF selections from Step 0b (no prompts)
         - Assigns original_type to all atoms

The key design principle is LAUNCHING existing modules
(protonation_state_analyzer, small_molecule_parameterizer, redox_detector, etc.)
and collecting their results, not reimplementing their logic.

By the end of preprocessing:
- We have a final structure with all H atoms
- We have RedoxSite objects synced to the final structure
- Every atom has its original_type assigned
- All FF selections are stored in workspace

Author: ProPrep Development Team
"""

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Set, TYPE_CHECKING

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
)

from proprep.forcefield_prep.ff_types import (
    AtomTypeAssignment,
    AtomSource,
    TerminalType
)

if TYPE_CHECKING:
    from proprep.structure_prep.comprehensive_redox_detector import RedoxSite

from proprep.utils.workflow_checklist import WorkflowStep, WorkflowChecklist

# Two except-handlers already called logger.debug without one being defined,
# which would have raised NameError instead of reporting the thing they caught.
logger = logging.getLogger(__name__)

# ``!entry.<unit>.unit.atoms table`` only. The naive ``'.unit.atoms' in line``
# also matches ``.unit.atomspertinfo``, whose rows repeat every atom name.
_LIB_ATOMS_TABLE_RE = re.compile(r"!entry\.[^.]+\.unit\.atoms\s+table", re.IGNORECASE)


def _is_hydrogen_name(atom_name: str) -> bool:
    """Whether a PDB/library atom name denotes a hydrogen.

    PDB names may lead with the digit of a branch index (``1HB``), so the first
    ALPHABETIC character decides.
    """
    for char in (atom_name or "").strip():
        if char.isalpha():
            return char.upper() == "H"
    return False




def _unwrap_serialized(value):
    """Recursively unwrap workspace serialization wrappers ({__type__, value})."""
    if isinstance(value, dict):
        if "__type__" in value and "value" in value:
            type_name = value["__type__"]
            inner = value["value"]
            if type_name == "tuple":
                return tuple(_unwrap_serialized(v) for v in inner)
            if type_name == "str":
                return str(inner)
            if type_name == "Path":
                return str(inner)
            # For enums like CenterType, extract the _value_ field
            if isinstance(inner, dict) and "_value_" in inner:
                return inner["_value_"]
            # For other wrapped objects, recursively unwrap the inner dict
            return _unwrap_serialized(inner)
        return {k: _unwrap_serialized(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_unwrap_serialized(v) for v in value]
    return value


def _ensure_redox_site_objects(sites: list) -> list:
    """Convert dicts back to RedoxSite objects if needed (e.g., after JSON resume)."""
    if not sites:
        return sites
    # Already proper objects?
    if hasattr(sites[0], 'centers'):
        return sites
    # Need reconstruction
    from proprep.structure_prep.comprehensive_redox_detector import (
        RedoxSite, RedoxCenter, RedoxSiteAtom, RedoxSiteBond, CenterType
    )
    # Unwrap all serialization wrappers first
    sites = _unwrap_serialized(sites)

    result = []
    for sd in sites:
        site = RedoxSite(sd["site_id"], sd.get("structure_id", "imported"))
        site.site_type = sd.get("site_type", "")
        for cd in sd.get("centers", []):
            coords = cd.get("coords", cd.get("coordinates", (0, 0, 0)))
            if isinstance(coords, (list, tuple)):
                coords = tuple(round(float(x), 3) for x in coords)
            center_type_val = cd.get("center_type", "metal_ion")
            center = RedoxCenter(
                chain=cd["chain"], resname=cd["resname"], resid=cd["resid"],
                atom_name=cd.get("atom_name"),
                insertion_code=cd.get("insertion_code", ""),
                altloc=cd.get("altloc", ""),
                coords=coords,
                center_type=CenterType(center_type_val),
                element=cd.get("element"),
            )
            site.add_center(center)
        for ad in sd.get("atoms", []):
            coords = ad.get("coords", ad.get("coordinates", (0, 0, 0)))
            if isinstance(coords, (list, tuple)):
                coords = tuple(round(float(x), 3) for x in coords)
            atom = RedoxSiteAtom(
                chain=ad["chain"], resname=ad["resname"], resid=ad["resid"],
                atom_name=ad["atom_name"],
                coords=coords,
                element=ad["element"],
            )
            site.add_atom(atom)
        for bd in sd.get("bonds", []):
            a1c = bd.get("atom1_coords", bd.get("atom1_coordinates", (0, 0, 0)))
            a2c = bd.get("atom2_coords", bd.get("atom2_coordinates", (0, 0, 0)))
            if isinstance(a1c, (list, tuple)):
                a1c = tuple(round(float(x), 3) for x in a1c)
            if isinstance(a2c, (list, tuple)):
                a2c = tuple(round(float(x), 3) for x in a2c)
            bond = RedoxSiteBond(
                atom1_coords=a1c,
                atom2_coords=a2c,
                bond_type=bd["bond_type"],
                chemical_type=bd.get("chemical_type", "unknown"),
                distance=float(bd["distance"]),
                atom1_element=bd.get("atom1_element", ""),
                atom2_element=bd.get("atom2_element", ""),
                atom1_residue_info=bd.get("atom1", bd.get("atom1_residue_info", "")),
                atom2_residue_info=bd.get("atom2", bd.get("atom2_residue_info", "")),
                treatment=bd.get("treatment", "bonded"),
            )
            site.bonds.append(bond)
        result.append(site)
    return result


def _redox_site_has_metal(site) -> bool:
    """True if a redox site has a metal center (an isolated metal ion or an
    organometallic cofactor with an embedded metal).

    MCPB parameterization only applies to metal sites; purely organic cofactors
    (e.g. an inhibitor/ligand) are handled by the Small Molecule Parameterizer
    and must be skipped by every MCPB step, or the workflow stalls trying to
    build a QM model around a metal that isn't there. Accepts either a RedoxSite
    object or its serialized-dict form.
    """
    metal_values = {"metal_ion", "organometallic_cofactor"}
    centers = site.get("centers", []) if isinstance(site, dict) else getattr(site, "centers", [])
    for c in centers:
        ct = c.get("center_type") if isinstance(c, dict) else getattr(c, "center_type", None)
        ct_value = ct.value if hasattr(ct, "value") else str(ct)
        if ct_value in metal_values:
            return True
    return False


# =============================================================================
# Constants
# =============================================================================

# Standard protein residues
STANDARD_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # Protonation variants
    'HIE', 'HID', 'HIP', 'CYX', 'CYM', 'ASH', 'GLH', 'LYN',
    # Terminal variants
    'ACE', 'NME', 'NHE',
}

# Water residues
WATER_RESIDUES = {'HOH', 'WAT', 'TIP', 'TIP3', 'TP3', 'SPC', 'T3P', 'SOL'}

# Metal elements for detecting organometallic residues
METAL_ELEMENTS = {
    'Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
    'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo',
    'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 'Ba', 'La', 'Ce',
    'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb',
    'Bi', 'Po', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu',
}


@dataclass
class MetalInfo:
    """
    Tracks an atom for removal before tLEaP and reinsertion after.

    Used for isolated metal ions (category D), embedded metals in
    organometallic residues (category C), and every atom of a pure inorganic
    metal cluster (category F).

    cluster_id groups the atoms of one pure-cluster residue (e.g. an Fe2S2's
    two Fe and two bridging S) so reinsertion recreates them as a SINGLE
    residue rather than one residue per atom. It is None for single metals.
    """
    atom_name: str
    element: str
    coords: Tuple[float, float, float]
    original_chain: str
    original_resid: int
    original_resname: str
    is_isolated: bool  # True = own residue (category D or F), False = embedded in C
    cluster_id: Optional[str] = None  # set for category-F cluster atoms; groups them

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'atom_name': self.atom_name,
            'element': self.element,
            'coords': list(self.coords),
            'original_chain': self.original_chain,
            'original_resid': self.original_resid,
            'original_resname': self.original_resname,
            'is_isolated': self.is_isolated,
            'cluster_id': self.cluster_id,
        }

    @classmethod
    def from_dict(cls, d: dict) -> 'MetalInfo':
        """Create from dictionary.

        Tolerates checklist-state serialization. Coordinates come from
        BioPython as numpy floats; older state files stored each one as
        ``{"__type__": "str", "value": "-46.078"}`` because the serializer had
        no case for a numpy scalar, so ``coords`` arrived as a list of dicts of
        strings and PDBIO refused to write the atom. Unwrap and coerce, so a
        state file written before that fix still resumes.
        """
        from proprep.utils.workspace import unwrap_serialized

        d = unwrap_serialized(d)

        def _xyz(raw):
            out = []
            for component in (raw or ()):
                try:
                    out.append(float(component))
                except (TypeError, ValueError):
                    out.append(0.0)
            return tuple(out)

        return cls(
            atom_name=d['atom_name'],
            element=d['element'],
            coords=_xyz(d.get('coords')),
            original_chain=d['original_chain'],
            original_resid=int(d['original_resid']),
            original_resname=d['original_resname'],
            is_isolated=bool(d['is_isolated']),
            cluster_id=d.get('cluster_id'),
        )


# =============================================================================
# Workflow Step Definitions
# =============================================================================

# Define all preprocessing and parameterization steps declaratively.
# This makes it easy to add/remove/rearrange steps.
#
# Triage Categories:
#   A: Protein residues
#   B: Organic small molecules (no metal)
#   C: Water
#   D: Isolated metal ions (single atom, coordinated to protein)
#   E: Organometallic small molecules (contains embedded metal)

PREPROCESSING_STEPS = [
    # -------------------------------------------------------------------------
    # STRUCTURE PREPARATION
    # -------------------------------------------------------------------------
    WorkflowStep(
        id="prep-1",
        name="Structure Filtering",
        description="Select chains, remove unwanted residues (optional)",
        handler="_checklist_prep_1_filtering",
        section="Structure Preparation",
        dependencies=[],
        optional=True,
    ),
    WorkflowStep(
        id="prep-1b",
        name="Structure Completeness",
        description="Resolve alternate conformations, missing atoms/residues (optional)",
        handler="_checklist_prep_1b_completeness",
        section="Structure Preparation",
        dependencies=[],
        optional=True,
    ),
    WorkflowStep(
        id="prep-2",
        name="Structure Triage",
        description="Categorize: protein/organic/water/isolated metal/organometallic",
        handler="_checklist_prep_2_triage",
        section="Structure Preparation",
        dependencies=[],
    ),

    # -------------------------------------------------------------------------
    # COMPONENT PARAMETERIZATION
    # -------------------------------------------------------------------------
    WorkflowStep(
        id="param-1",
        name="Protein: FF + Protonation",
        description="Select force field, analyze protonation states (HIS, ASP, etc.)",
        handler="_checklist_param_1_protein",
        section="Component Parameterization",
        dependencies=["prep-2"],
    ),
    WorkflowStep(
        id="param-2",
        name="Water: Model Selection",
        description="Select water model for coordinating waters",
        handler="_checklist_param_2_water",
        section="Component Parameterization",
        dependencies=["prep-2"],
    ),
    WorkflowStep(
        id="param-3",
        name="Organic Small Molecules",
        description="Get lib/frcmod or generate parameters",
        handler="_checklist_param_3_organic",
        section="Component Parameterization",
        dependencies=["prep-2"],
    ),
    WorkflowStep(
        id="param-4",
        name="Organometallic Small Molecules",
        description="Edit residue if needed, get lib/frcmod or generate",
        handler="_checklist_param_4_organometallic",
        section="Component Parameterization",
        dependencies=["prep-2"],
    ),
    # Metal clusters listed before isolated metals: a cluster is the more
    # complete metal object, the lone ion the degenerate case. (Display order
    # follows this list; the param-5/param-6 ids are internal only.)
    WorkflowStep(
        id="param-6",
        name="Metal Clusters",
        description="Pure inorganic clusters (Fe-S, etc.): withhold whole cluster for MCPB",
        handler="_checklist_param_6_metal_clusters",
        section="Component Parameterization",
        dependencies=["prep-2"],
    ),
    WorkflowStep(
        id="param-5",
        name="Isolated Metal Ions",
        description="Acknowledge metals for removal/reinsertion",
        handler="_checklist_param_5_isolated_metals",
        section="Component Parameterization",
        dependencies=["prep-2"],
    ),

    # -------------------------------------------------------------------------
    # STRUCTURE ASSEMBLY
    # -------------------------------------------------------------------------
    WorkflowStep(
        id="assembly-1",
        name="Structure Recombination",
        description="tLEaP processing + metal reinsertion + type collection",
        handler="_checklist_assembly_1_recombination",
        section="Structure Assembly",
        dependencies=["param-1", "param-2", "param-3", "param-4", "param-5"],
    ),
    WorkflowStep(
        id="assembly-2",
        name="Redox Site Sync",
        description="Update site definitions for new coordinates",
        handler="_checklist_assembly_2_redox_sync",
        section="Structure Assembly",
        dependencies=["assembly-1"],
    ),

    # -------------------------------------------------------------------------
    # METAL SITE PARAMETERIZATION
    # -------------------------------------------------------------------------
    WorkflowStep(
        id="mcpb-1",
        name="MCPB Atom Typing",
        description="M*/Y* renaming, build small/large models",
        handler="_checklist_mcpb_1_typing",
        section="Metal Site Parameterization",
        dependencies=["assembly-2"],
    ),
    WorkflowStep(
        id="mcpb-2",
        name="Bonded Parameters",
        description="Gaussian optimization + Seminario method",
        handler="_checklist_mcpb_2_bonded",
        section="Metal Site Parameterization",
        dependencies=["mcpb-1"],
        checkpoint=True,
        checkpoint_message="Gaussian geometry optimization required. Run Gaussian, then resume.",
    ),
    WorkflowStep(
        id="mcpb-3",
        name="RESP Charges",
        description="Gaussian ESP calculation for charges",
        handler="_checklist_mcpb_3_resp",
        section="Metal Site Parameterization",
        dependencies=["mcpb-2"],
        checkpoint=True,
        checkpoint_message="Gaussian ESP calculation required. Run Gaussian, then resume.",
    ),
    WorkflowStep(
        id="mcpb-4",
        name="Force Field Integration",
        description="Combine parameters into tLEaP-ready files",
        handler="_checklist_mcpb_4_integration",
        section="Metal Site Parameterization",
        dependencies=["mcpb-3"],
    ),
]


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class PreprocessingResult:
    """
    Result of structure preprocessing.

    This dataclass contains everything Step 1 needs:
    - prepared_pdb: The final combined structure with H atoms
    - redox_sites: Fresh RedoxSite objects detected on final structure
    - type_assignments: Every atom's original_type (Step 1 adds renamed_type)
    """
    success: bool
    prepared_pdb: Optional[Path] = None

    # Fresh RedoxSites detected on final structure (Step 0f)
    redox_sites: List['RedoxSite'] = field(default_factory=list)

    # The key output: original_type for ALL atoms (Step 0g)
    type_assignments: Dict[Tuple[float, float, float], AtomTypeAssignment] = field(
        default_factory=dict
    )

    # Supporting data
    ff_data: Any = None  # ForceFieldData
    residue_map: Dict[str, str] = field(default_factory=dict)  # PDB_resname → FF_resname
    atom_map: Dict[Tuple[str, str], str] = field(default_factory=dict)  # (PDB_resname, PDB_atom) → FF_atom
    small_mol_results: Dict[str, Dict] = field(default_factory=dict)  # resname → {mol2, frcmod, types, charges}
    excluded_atoms: List[Tuple[float, float, float]] = field(default_factory=list)
    h_caps: List[Dict] = field(default_factory=list)
    triage: Dict[str, str] = field(default_factory=dict)  # residue_id → category
    error_message: str = ""


# =============================================================================
# Constants
# =============================================================================

STANDARD_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
    'TYR', 'VAL', 'HIE', 'HID', 'HIP', 'CYX', 'ACE', 'NME'
}

WATER_RESIDUES = {'HOH', 'WAT', 'TIP', 'TIP3', 'SPC'}

# Note: METAL_ELEMENTS is defined at the top of the file (comprehensive list)


# =============================================================================
# Main Class
# =============================================================================

class StructurePreprocessor:
    """
    Prepares structures for metal site parameterization.

    This module orchestrates the preprocessing workflow:
    0a. Structure filtering (optional, via PDB Filter module)
    0b. Structure triage + force field selection
    0c. Hydrogen addition (protein via tLEaP, non-standard via reduce)
    0d. Atom exclusion + H capping
    0e. Structure recombination → prepared_structure.pdb
    0f. Redox site sync (updates existing sites to match prepared structure)
    0g. Atom typing (uses FF selections from Step 0b)

    The key design principle is LAUNCHING existing modules
    (protonation_state_analyzer, small_molecule_parameterizer, redox_detector, etc.)
    and collecting their results, not reimplementing their logic.

    By the end of preprocessing:
    - We have a final structure with all H atoms
    - We have RedoxSite objects synced to the final structure
    - Every atom has its original_type assigned
    - All FF selections are stored in workspace
    """

    def __init__(self, processor, console: Console = None):
        """
        Initialize the preprocessor.

        Args:
            processor: The PDBProcessor instance (for workspace access and module launching)
            console: Rich Console for output
        """
        self.processor = processor
        self.console = console or Console()
        self.workspace = processor._get_workspace() if processor else None

        # State accumulated through preprocessing
        self.ff_data = None
        self.type_assignments: Dict[Tuple[float, float, float], AtomTypeAssignment] = {}
        self.residue_map: Dict[str, str] = {}  # PDB_resname → FF_resname
        self.atom_map: Dict[Tuple[str, str], str] = {}  # (PDB_res, PDB_atom) → FF_atom
        self.triage_results: Dict[str, str] = {}  # residue_id → category
        self.excluded_atoms: List[Tuple[float, float, float]] = []
        self.h_caps: List[Dict] = []  # {kept_atom_coords, h_coords, h_name}
        self.small_mol_results: Dict[str, Dict] = {}  # resname → {mol2, frcmod, types, charges}
        self.redox_sites: List['RedoxSite'] = []  # Fresh from detector on final structure

        # Workflow state (set by checklist or run_preprocessing)
        self._pdb_file: Optional[str] = None
        self._output_dir: Optional[Path] = None
        self._interactive: bool = True
        self._h_results: Dict[str, Path] = {}  # Category/res_key → PDB with H
        self._final_pdb: Optional[Path] = None

    # =========================================================================
    # Checklist-Based Workflow Entry Point
    # =========================================================================

    def run_checklist(
        self,
        pdb_file: str,
        output_dir: Path,
        interactive: bool = True
    ) -> PreprocessingResult:
        """
        Run the workflow with interactive checklist interface.

        This is the recommended entry point. It provides:
        - Visual progress tracking with checkmarks
        - Ability to re-run steps
        - State persistence for resumption (e.g., after Gaussian calculations)
        - Clear overview of all steps

        Args:
            pdb_file: Path to input PDB file
            output_dir: Directory for output files
            interactive: Whether to prompt user for input

        Returns:
            PreprocessingResult with prepared structure and all metadata
        """
        # Store workflow state for step handlers to access
        self._pdb_file = pdb_file
        self._output_dir = Path(output_dir).resolve()  # Absolute path to avoid issues after chdir
        self._output_dir.mkdir(parents=True, exist_ok=True)
        self._interactive = interactive

        # Store in workspace for modules to find
        if self.workspace:
            self.workspace.set("structure_pdb_file", pdb_file)
            self.workspace.set("preprocessing_output_dir", str(self._output_dir))

        # Create and run the checklist
        checklist = WorkflowChecklist(
            steps=PREPROCESSING_STEPS,
            executor=self,
            processor=self.processor,
            workflow_name="Metal Site Parameterization",
            console=self.console,
            state_dir=self._output_dir,  # Save state in output directory for resume
        )

        success = checklist.run(pdb_file=pdb_file)

        # Build and return result
        if success and self._final_pdb:
            return PreprocessingResult(
                success=True,
                prepared_pdb=self._final_pdb,
                redox_sites=self.redox_sites,
                type_assignments=self.type_assignments,
                ff_data=self.ff_data,
                residue_map=self.residue_map,
                atom_map=self.atom_map,
                small_mol_results=self.small_mol_results,
                excluded_atoms=self.excluded_atoms,
                h_caps=self.h_caps,
                triage=self.triage_results
            )
        else:
            return PreprocessingResult(
                success=False,
                error_message="Workflow incomplete or cancelled"
            )

    # =========================================================================
    # Checklist Step Handlers
    # =========================================================================
    # These wrapper methods are called by the WorkflowChecklist.
    # Handler naming: _checklist_{section}_{number}_{description}

    # -------------------------------------------------------------------------
    # STRUCTURE PREPARATION
    # -------------------------------------------------------------------------

    def _checklist_prep_1_filtering(self) -> dict:
        """Checklist handler for prep-1: Structure Filtering."""
        result = self._step_0a_structure_filtering(self._pdb_file, self._interactive)
        self._pdb_file = result  # Update for subsequent steps
        return {"summary": f"Using {Path(result).name}"}

    def _checklist_prep_1b_completeness(self) -> dict:
        """Checklist handler for prep-1b: Structure Completeness."""
        result = self._step_0a2_structure_completeness(self._pdb_file, self._interactive)
        self._pdb_file = result  # Update for subsequent steps
        return {"summary": f"Using {Path(result).name}"}

    def _ensure_triage_results(self) -> Dict[str, str]:
        """Triage categories, restoring them after a resume.

        ``triage_results`` is populated by step 3 and lives on the instance. A
        resumed run builds a new preprocessor with step 3 already marked
        complete, so it never runs and the attribute stays empty -- and every
        step that reads it concludes the structure has no such residues.

        Restores from the workspace, falling back to re-running triage, which
        is deterministic and needs only the structure.
        """
        if self.triage_results:
            return self.triage_results

        if self.workspace:
            saved = self.workspace.get("preprocessing_triage", None)
            if isinstance(saved, dict) and saved:
                self.triage_results = dict(saved)
                self.console.print(
                    f"[grey50]Restored triage for {len(self.triage_results)} "
                    f"residue(s) from workspace[/grey50]")
                return self.triage_results

        if self._pdb_file:
            self.console.print(
                "[grey50]Triage results unavailable; re-running "
                "categorization[/grey50]")
            try:
                self.triage_results = self._run_triage_only(self._pdb_file)
            except Exception as exc:  # noqa: BLE001
                logger.debug("Could not re-run triage: %s", exc)

        return self.triage_results

    def _checklist_prep_2_triage(self) -> dict:
        """Checklist handler for prep-2: Structure Triage (categorization only)."""
        # Run triage without FF selection (FF selection moved to param steps)
        self.triage_results = self._run_triage_only(self._pdb_file)

        # Count categories
        cat_counts = {}
        for cat in self.triage_results.values():
            cat_counts[cat] = cat_counts.get(cat, 0) + 1

        summary_parts = []
        if 'A' in cat_counts:
            summary_parts.append(f"{cat_counts['A']} protein")
        if 'B' in cat_counts:
            summary_parts.append(f"{cat_counts['B']} organic")
        if 'C' in cat_counts:
            summary_parts.append(f"{cat_counts['C']} organometallic")
        if 'F' in cat_counts:
            summary_parts.append(f"{cat_counts['F']} metal cluster")
        if 'D' in cat_counts:
            summary_parts.append(f"{cat_counts['D']} isolated metal")
        if 'E' in cat_counts:
            summary_parts.append(f"{cat_counts['E']} water")

        return {"summary": ", ".join(summary_parts)}

    # -------------------------------------------------------------------------
    # COMPONENT PARAMETERIZATION
    # -------------------------------------------------------------------------

    def _checklist_param_1_protein(self) -> dict:
        """Checklist handler for param-1: Protein FF Selection + Protonation Analysis."""
        protein_residues = [k for k, v in self._ensure_triage_results().items() if v == 'A']
        if not protein_residues:
            return {"summary": "No protein residues - skipped"}

        self.console.print(Panel(
            f"[bold]Protein Force Field Selection[/bold]\n"
            f"{len(protein_residues)} protein residues found",
            border_style="cyan",
            expand=False
        ))

        # Select protein FF
        protein_ff = self._select_protein_forcefield()
        if not protein_ff:
            return {"summary": "No force field selected"}

        # Store selection (will be used by tLEaP in assembly step)
        if self.workspace:
            self.workspace.set("preprocessing_protein_ff", protein_ff)

        ff_name = protein_ff.split('.')[-1] if '.' in protein_ff else protein_ff

        # Launch protonation state analyzer
        protonation_done = False
        if self.processor:
            prot_analyzer = self.processor.get_module_instance("Protonation State Analyzer")
            if prot_analyzer:
                self.console.print("\n[bold]═══ Protonation State Analysis ═══[/bold]")
                self.console.print("[grey50]Determine protonation states for HIS, ASP, GLU, etc.[/grey50]\n")

                # Analyze protonation states
                if hasattr(prot_analyzer, 'analyze_protonation_states'):
                    prot_analyzer.analyze_protonation_states()
                    protonation_done = True
                elif hasattr(prot_analyzer, 'process'):
                    prot_analyzer.process(self.workspace)
                    protonation_done = True

                # Set residue names based on protonation states
                if protonation_done and hasattr(prot_analyzer, 'set_residue_names'):
                    self.console.print("\n[grey50]Setting residue names based on protonation states...[/grey50]")
                    prot_analyzer.set_residue_names()
            else:
                self.console.print("[yellow]Protonation State Analyzer not available[/yellow]")

        summary = f"Selected {ff_name}"
        if protonation_done:
            summary += " + protonation analyzed"
        return {"summary": summary}

    def _checklist_param_2_water(self) -> dict:
        """Checklist handler for param-2: Water Model Selection (no tLEaP yet).

        Always prompts for water model selection because:
        1. The system will likely be solvated later
        2. Water model choice affects ion parameters
        3. Must be compatible with protein forcefield
        """
        water_residues = [k for k, v in self._ensure_triage_results().items() if v == 'E']

        if water_residues:
            self.console.print(Panel(
                f"[bold]Water Model Selection[/bold]\n"
                f"{len(water_residues)} water molecules in structure",
                border_style="cyan",
                expand=False
            ))
        else:
            self.console.print(Panel(
                "[bold]Water Model Selection[/bold]\n"
                "No water in structure, but model needed for solvation",
                border_style="cyan",
                expand=False
            ))

        # Always select water model (needed for solvation and ion parameters)
        water_model = self._select_water_model()
        if not water_model:
            return {"summary": "No water model selected"}

        # Store selection
        if self.workspace:
            self.workspace.set("preprocessing_water_model", water_model)

        model_name = water_model.split('.')[-1] if '.' in water_model else water_model
        return {"summary": f"Selected {model_name}"}

    def _checklist_param_3_organic(self) -> dict:
        """Checklist handler for param-3: Organic Small Molecules."""
        organic_keys = [k for k, v in self._ensure_triage_results().items() if v == 'B']
        if not organic_keys:
            self.console.print("[grey50]No organic small molecules in structure - skipping[/grey50]")
            return {"summary": "No organic small molecules - skipped"}

        self.console.print(Panel(
            f"[bold]Organic Small Molecules[/bold]\n"
            f"{len(organic_keys)} residues to parameterize",
            border_style="cyan",
            expand=False
        ))

        organic_ff = {}
        for res_key in organic_keys:
            result = self._process_organic_residue(res_key)
            if result:
                resname = res_key.split(':')[2]
                organic_ff[resname] = result

        # Store selections
        if self.workspace:
            self.workspace.set("preprocessing_organic_ff", organic_ff)

        return {"summary": f"Configured {len(organic_ff)}/{len(organic_keys)} residues"}

    def _checklist_param_4_organometallic(self) -> dict:
        """Checklist handler for param-4: Organometallic Small Molecules."""
        orgmet_keys = [k for k, v in self._ensure_triage_results().items() if v == 'C']
        if not orgmet_keys:
            self.console.print("[grey50]No organometallic residues in structure - skipping[/grey50]")
            return {"summary": "No organometallic residues - skipped"}

        self.console.print(Panel(
            f"[bold]Organometallic Small Molecules[/bold]\n"
            f"{len(orgmet_keys)} residues containing metals",
            border_style="cyan",
            expand=False
        ))

        orgmet_ff = {}
        for res_key in orgmet_keys:
            result = self._process_organometallic_residue(res_key)
            if result:
                resname = res_key.split(':')[2]
                orgmet_ff[resname] = result

        # Store selections
        if self.workspace:
            self.workspace.set("preprocessing_organometallic_ff", orgmet_ff)

        return {"summary": f"Configured {len(orgmet_ff)}/{len(orgmet_keys)} residues"}

    def _checklist_param_5_isolated_metals(self) -> dict:
        """Checklist handler for param-5: Isolated Metal Ions."""
        metal_keys = [k for k, v in self._ensure_triage_results().items() if v == 'D']
        if not metal_keys:
            self.console.print("[grey50]No isolated metal ions in structure - skipping[/grey50]")
            return {"summary": "No isolated metal ions - skipped"}

        self.console.print(Panel(
            f"[bold]Isolated Metal Ions[/bold]\n"
            f"{len(metal_keys)} metals will be removed before tLEaP, reinserted after",
            border_style="cyan",
            expand=False
        ))

        isolated_metals = {}
        for res_key in metal_keys:
            metal_info = self._extract_metal_info_from_key(res_key, is_isolated=True)
            if metal_info:
                isolated_metals[res_key] = {'metal_info': metal_info.to_dict()}
                self.console.print(f"  [magenta]{res_key}[/magenta]: {metal_info.element} at "
                                   f"({metal_info.coords[0]:.2f}, {metal_info.coords[1]:.2f}, {metal_info.coords[2]:.2f})")

        # Store for later use
        if self.workspace:
            self.workspace.set("preprocessing_isolated_metals", isolated_metals)

        return {"summary": f"{len(isolated_metals)} metals tracked for reinsertion"}

    def _checklist_param_6_metal_clusters(self) -> dict:
        """Checklist handler for param-6: Pure inorganic metal clusters (Fe-S, etc.).

        A pure cluster has no organic fragment, so unlike an organometallic
        cofactor it is NOT split. Every atom of the cluster residue (metals plus
        bridging sulfides) is withheld from the standard-FF tLEaP pass — none of
        them have standard templates — and reinserted afterward as one residue.
        MCPB owns the cluster's internal and coordinating parameters via the
        detected redox site, so no per-fragment FF is collected here.
        """
        cluster_keys = [k for k, v in self._ensure_triage_results().items() if v == 'F']
        if not cluster_keys:
            self.console.print("[grey50]No pure metal clusters in structure - skipping[/grey50]")
            return {"summary": "No metal clusters - skipped"}

        self.console.print(Panel(
            f"[bold]Metal Clusters[/bold]\n"
            f"{len(cluster_keys)} pure inorganic cluster(s) — whole residue withheld "
            f"for MCPB, reinserted after tLEaP",
            border_style="cyan",
            expand=False
        ))

        # Offer hydrogens BEFORE the atoms are collected: everything downstream
        # (withholding, reinsertion, redox re-detection, typing, the QM models)
        # reads self._pdb_file, so an H added there needs no special casing.
        self._offer_cluster_hydrogens(cluster_keys)

        metal_clusters = {}
        for res_key in cluster_keys:
            cluster_atoms = self._extract_cluster_atoms_from_key(res_key)
            if cluster_atoms:
                metal_clusters[res_key] = {
                    'atoms': [a.to_dict() for a in cluster_atoms],
                }
                metals = [a for a in cluster_atoms if a.element.title() in METAL_ELEMENTS]
                self.console.print(
                    f"  [magenta]{res_key}[/magenta]: {len(cluster_atoms)} atoms "
                    f"({len(metals)} metal, {len(cluster_atoms) - len(metals)} bridging)")

        if self.workspace:
            self.workspace.set("preprocessing_metal_clusters", metal_clusters)

        return {"summary": f"{len(metal_clusters)} cluster(s) tracked for reinsertion"}

    # -------------------------------------------------------------------------
    # STRUCTURE ASSEMBLY
    # -------------------------------------------------------------------------

    def _checklist_assembly_1_recombination(self) -> dict:
        """
        Checklist handler for assembly-1: Structure Recombination.

        This is the main assembly step that:
        1. Removes metals that need removal (isolated + organometallic with organic-only params)
        2. Runs Topology Generator with all collected FF info
        3. Reinserts metals into the tLEaP output
        4. Collects atom types from prmtop
        """
        self.console.print(Panel(
            "[bold]Structure Recombination[/bold]\n"
            "Running tLEaP with all components, then reinserting metals",
            border_style="cyan",
            expand=False
        ))

        # Determine input PDB in order of preference:
        # 1. protonation_pdb_file - has MD residue names (HID, HIE, ASH, etc.)
        # 2. filtered_pdb_file - from PDB filter step
        # 3. self._pdb_file - original input
        input_pdb = self._pdb_file
        if self.workspace:
            protonation_pdb = self.workspace.get("protonation_pdb_file")
            filtered_pdb = self.workspace.get("filtered_pdb_file")

            if protonation_pdb and Path(protonation_pdb).exists():
                input_pdb = protonation_pdb
                self.console.print(f"[grey50]Using protonation-updated structure: {Path(protonation_pdb).name}[/grey50]")
            elif filtered_pdb and Path(filtered_pdb).exists():
                input_pdb = filtered_pdb
                self.console.print(f"[grey50]Using filtered structure: {Path(filtered_pdb).name}[/grey50]")

        # 1. Collect all metals to remove
        metals_to_remove = self._collect_metals_to_remove()
        if metals_to_remove:
            self.console.print(f"[grey50]Will remove {len(metals_to_remove)} metal(s) before tLEaP[/grey50]")

        # 2. Create metal-free structure
        if metals_to_remove:
            metal_free_pdb = self._remove_metals_from_structure(metals_to_remove, input_pdb)
            if self.workspace:
                self.workspace.set("preprocessing_metal_free_pdb", str(metal_free_pdb))
        else:
            metal_free_pdb = Path(input_pdb)

        # 3. Apply atom name mappings (PDB names → lib names)
        remapped_pdb = self._apply_atom_name_mappings(metal_free_pdb)
        if remapped_pdb:
            metal_free_pdb = remapped_pdb

        # 4. Build residue sequence map BEFORE tLEaP (for later sync)
        residue_sequence_map = self._build_residue_sequence_map(metal_free_pdb)
        if self.workspace and residue_sequence_map:
            self.workspace.set("preprocessing_residue_sequence_map", residue_sequence_map)

        # 5. Configure tLEaP with all FF info
        self._configure_tleap_for_assembly(metal_free_pdb)

        # 6. Run Topology Generator
        tleap_success = self._run_tleap_assembly()

        if not tleap_success:
            # Mark the step FAILED, not completed: the steps after this one read
            # the prepared structure it did not produce, and would otherwise
            # fail on a missing file instead of on the real cause.
            return {"summary": "tLEaP failed", "success": False}

        # 7. Get tLEaP output and convert to PDB
        parm7 = self.workspace.get("parm7_file") if self.workspace else None
        rst7 = self.workspace.get("rst7_file") if self.workspace else None

        if not parm7 or not rst7:
            return {"summary": "tLEaP did not produce output files"}

        tleap_pdb = self._output_dir / "tleap_output.pdb"
        try:
            self._convert_amber_to_pdb(parm7, rst7, tleap_pdb)
        except RuntimeError as exc:
            # Everything after this reads the file that was not written.
            self.console.print(f"[red]{exc}[/red]")
            return {"summary": f"Could not convert tLEaP output: {exc}",
                    "success": False}

        # 8. Insert metals back
        if metals_to_remove:
            final_pdb = self._insert_metals(tleap_pdb, metals_to_remove)
        else:
            final_pdb = tleap_pdb

        # Store final structure
        self._final_pdb = final_pdb
        if self.workspace:
            self.workspace.set("prepared_pdb", str(final_pdb))

        # 9. Collect atom types and charges from prmtop
        # The prmtop IS the force field - no need to load leaprcs separately
        atom_data = self._extract_atom_data_from_prmtop(parm7)
        if self.workspace:
            self.workspace.set("preprocessing_atom_data", atom_data)

        self.console.print(f"[green]✓ Created {final_pdb.name}[/green]")
        if metals_to_remove:
            self.console.print(f"[green]✓ Reinserted {len(metals_to_remove)} metal(s)[/green]")
        self.console.print(f"[green]✓ Collected {len(atom_data)} atoms (types + charges)[/green]")

        return {"summary": f"Created {final_pdb.name} with {len(atom_data)} typed atoms"}

    def _checklist_assembly_2_redox_sync(self) -> dict:
        """Checklist handler for assembly-2: Redox Site Sync."""
        # Restore _final_pdb from workspace if resuming from saved state
        if self._final_pdb is None and self.workspace:
            prepared_pdb = self.workspace.get("prepared_pdb")
            if prepared_pdb:
                self._final_pdb = Path(prepared_pdb)

        self.redox_sites = self._step_0f_redox_site_sync(self._final_pdb, self._interactive)
        return {"summary": f"{len(self.redox_sites)} redox site(s) synchronized"}

    # -------------------------------------------------------------------------
    # METAL SITE PARAMETERIZATION
    # -------------------------------------------------------------------------

    def _scan_prior_mcpb_types(self) -> Dict[str, Dict[str, List[int]]]:
        """M*/Y* type positions already used, per prior site directory.

        The fingerprint files an earlier run wrote are the record of which
        names it consumed (``202-MN-MN  1  MN -> M1``). Reading them makes the
        numbering recoverable in a fresh ProPrep session, where the workspace
        is empty — without asking anyone to count Y types across sites, which
        is error-prone precisely because the labels stop matching the count at
        ``Y9``/``YA``.

        Scans sibling ``metal_site_params_*`` directories, so it sees work done
        for this structure under a different site's output directory.
        """
        from proprep.forcefield_prep.metal_site_parameterizer import mcpb_type_index

        used: Dict[str, Dict[str, List[int]]] = {}
        roots = []
        if self._output_dir:
            roots.append(Path(self._output_dir))
            parent = Path(self._output_dir).parent
            roots.extend(sorted(parent.glob("metal_site_params_*")))

        seen_dirs = set()
        for root in roots:
            for fp in sorted(Path(root).glob("site_*/models/standard.fingerprint")):
                key = str(fp.parent.parent)
                if key in seen_dirs:
                    continue
                seen_dirs.add(key)
                metals, ligands = [], []
                try:
                    for line in fp.read_text().splitlines():
                        if "->" not in line:
                            continue
                        name = line.split("->")[-1].strip()
                        pos = mcpb_type_index(name)
                        if pos is None:
                            continue
                        (metals if name.upper().startswith("M") else ligands).append(pos)
                except OSError:
                    continue
                if metals or ligands:
                    used[key] = {"metal": sorted(set(metals)),
                                 "ligand": sorted(set(ligands))}
        return used

    def _offer_prior_type_reuse(self, site_output_dir: Path, site_id: str):
        """Offer to reuse the M*/Y* names a previous run gave THIS site.

        Returns the (metal_start, ligand_start) that reproduces the earlier
        naming, or None to allocate fresh names after the high-water mark.
        """
        from proprep.forcefield_prep.metal_site_parameterizer import (
            mcpb_type_index, MCPB_METAL_TYPE_NAMES, MCPB_LIGAND_TYPE_NAMES,
        )

        fp = Path(site_output_dir) / "models" / "standard.fingerprint"
        if not fp.exists():
            return None

        metals, ligands = [], []
        try:
            for line in fp.read_text().splitlines():
                if "->" not in line:
                    continue
                name = line.split("->")[-1].strip()
                pos = mcpb_type_index(name)
                if pos is None:
                    continue
                (metals if name.upper().startswith("M") else ligands).append(pos)
        except OSError:
            return None

        if not metals and not ligands:
            return None

        def _names(positions, table):
            picked = sorted(set(positions))
            return ", ".join(table[p] for p in picked if p < len(table)) or "-"

        self.console.print(
            f"\n[yellow]{site_id} was parameterized before, in "
            f"{Path(site_output_dir).name}:[/yellow]")
        self.console.print(
            f"  [grey50]metals {_names(metals, MCPB_METAL_TYPE_NAMES)} · "
            f"ligating atoms {_names(ligands, MCPB_LIGAND_TYPE_NAMES)}[/grey50]")

        if not self._interactive:
            return None

        if confirm_with_context(
            self.processor,
            "Reuse those atom-type names (replaces the earlier entry)?",
            default=True,
            module="MCPB Atom Typing",
            description="Reuse a re-run site's previous M*/Y* names",
        ):
            return (min(metals) if metals else 0,
                    min(ligands) if ligands else 0)
        return None

    def _seed_mcpb_type_offsets(self) -> Tuple[int, int]:
        """Starting M*/Y* positions for this run.

        Precedence: the workspace high-water mark (sites parameterized earlier
        in this ProPrep session), then the fingerprints of prior runs on disk (a
        fresh session), then zero. What was found is shown and can be
        overridden, so the number is confirmed rather than computed by hand.
        """
        from proprep.forcefield_prep.metal_site_parameterizer import (
            MCPB_METAL_TYPE_NAMES, MCPB_LIGAND_TYPE_NAMES,
        )

        metal_next = ligand_next = 0
        source = None

        stored = self.workspace.get("mcpb_type_offsets") if self.workspace else None
        if isinstance(stored, dict) and (stored.get("metal") or stored.get("ligand")):
            metal_next = int(stored.get("metal") or 0)
            ligand_next = int(stored.get("ligand") or 0)
            source = "this session"

        prior = self._scan_prior_mcpb_types()
        if prior:
            scan_metal = max((max(v["metal"]) + 1 for v in prior.values() if v["metal"]),
                             default=0)
            scan_ligand = max((max(v["ligand"]) + 1 for v in prior.values() if v["ligand"]),
                              default=0)
            if scan_metal > metal_next or scan_ligand > ligand_next:
                metal_next = max(metal_next, scan_metal)
                ligand_next = max(ligand_next, scan_ligand)
                source = "earlier parameterization on disk"

        if not source or (metal_next == 0 and ligand_next == 0):
            return 0, 0

        def _describe(positions, names):
            if not positions:
                return "-"
            picked = [names[p] for p in positions if p < len(names)]
            return f"{picked[0]}-{picked[-1]}" if len(picked) > 1 else picked[0]

        self.console.print("\n[bold]Metal atom-type numbering[/bold]")
        self.console.print(f"[grey50]Found {source}:[/grey50]")
        for site_dir, v in sorted(prior.items()):
            self.console.print(
                f"  [grey50]{Path(site_dir).parent.name}/{Path(site_dir).name}   "
                f"{_describe(v['metal'], MCPB_METAL_TYPE_NAMES)}   "
                f"{_describe(v['ligand'], MCPB_LIGAND_TYPE_NAMES)}[/grey50]")

        next_metal = (MCPB_METAL_TYPE_NAMES[metal_next]
                      if metal_next < len(MCPB_METAL_TYPE_NAMES) else "exhausted")
        next_ligand = (MCPB_LIGAND_TYPE_NAMES[ligand_next]
                       if ligand_next < len(MCPB_LIGAND_TYPE_NAMES) else "exhausted")
        self.console.print(
            f"  Continuing from: metals at [bold]{next_metal}[/bold], "
            f"ligating atoms at [bold]{next_ligand}[/bold]")

        if self._interactive and confirm_with_context(
            self.processor,
            "Override these atom-type offsets?",
            default=False,
            module="MCPB Atom Typing",
            description="Override the M*/Y* starting offsets",
        ):
            metal_next = int_prompt_with_context(
                self.processor,
                "Metal types already used (M positions to skip)",
                default=metal_next,
                module="MCPB Atom Typing",
                description="Metal atom-type offset",
            )
            ligand_next = int_prompt_with_context(
                self.processor,
                "Ligating-atom types already used (Y positions to skip)",
                default=ligand_next,
                module="MCPB Atom Typing",
                description="Ligand atom-type offset",
            )

        return max(0, metal_next), max(0, ligand_next)

    def _checklist_mcpb_1_typing(self) -> dict:
        """
        Checklist handler for mcpb-1: MCPB Atom Typing.

        Gets all data from workspace:
        - redox_sites (from assembly-2)
        - preprocessing_atom_data (from assembly-1) - types AND charges from prmtop

        The prmtop IS the force field, so we don't need to load leaprcs.
        Metal charges are collected from user in Step 1 (metals aren't in prmtop).

        Then calls MetalSiteWorkflowManager._run_step1() which:
        1. Prompts user for metal charge/spin
        2. Applies M*/Y* renaming
        3. Builds small/large QM models
        4. Generates MCPB fingerprint files
        """
        from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager

        # Get redox sites from instance variable or workspace (for resume case)
        redox_sites = self.redox_sites
        if not redox_sites and self.workspace:
            redox_sites = _ensure_redox_site_objects(self.workspace.get("detected_redox_sites", []))
            if redox_sites:
                self.redox_sites = redox_sites  # Cache for future use
                self.console.print(f"[grey50]Restored {len(redox_sites)} redox site(s) from workspace[/grey50]")

        if not redox_sites:
            self.console.print("[red]No redox sites available for MCPB Step 1[/red]")
            return {"summary": "Failed - no redox sites"}

        # Build type_assignments from workspace (types + charges from prmtop)
        if not self.type_assignments:
            self._build_mcpb_data_from_workspace()

        # Process each redox site
        results = []

        # Helper to get site_id whether site is a dict or RedoxSite object
        def get_site_id(site, idx):
            if isinstance(site, dict):
                return site.get('site_id', f'site_{idx + 1}')
            return site.site_id

        # MCPB is only for metal sites (an isolated metal ion or an
        # organometallic cofactor with an embedded metal). Purely organic
        # cofactors (e.g. an inhibitor/ligand) are parameterized by the Small
        # Molecule Parameterizer and must NOT be pushed through MCPB atom typing
        # — otherwise the workflow stalls trying to build a QM model around a
        # site that has no metal.
        metal_site_count = sum(1 for s in redox_sites if _redox_site_has_metal(s))

        # Only the selected sites are parameterized. The selection is recorded
        # by the Force Field Parameterizer; an empty/absent list means every
        # site, which is what a resumed run and any other entry point get.
        selected_ids = set()
        if self.workspace:
            selected_ids = {s for s in (self.workspace.get("mcpb_selected_site_ids") or []) if s}

        # Running M*/Y* atom-type offsets so metal sites get globally unique
        # type names (site 1 -> M1/Y1..Yn, site 2 -> M2/Y(n+1)..). Restarting at
        # zero would collide once every site's frcmod/mol2 loads into one tLEaP
        # session. Because sites can now be parameterized in separate runs, the
        # starting point is seeded rather than assumed to be zero.
        metal_type_offset, ligand_type_offset = self._seed_mcpb_type_offsets()

        for idx, site in enumerate(redox_sites):
            site_id = get_site_id(site, idx)

            # Skip rather than filter the list: the site_N output directory is
            # named from this index, so a site keeps the same directory whether
            # or not its neighbours were selected this run.
            if selected_ids and site_id not in selected_ids:
                self.console.print(
                    f"\n[grey50]Skipping {site_id}: not selected for this run[/grey50]")
                continue

            # Skip non-metal sites (parameterized separately by the small-molecule path)
            if not _redox_site_has_metal(site):
                self.console.print(f"\n[grey50]Skipping {site_id}: organic cofactor (no metal) — parameterized by the Small Molecule Parameterizer, not MCPB[/grey50]")
                continue

            self.console.print(f"\n[bold cyan]Processing site {idx + 1}/{len(redox_sites)}: {site_id}[/bold cyan]")

            # Create MetalSiteWorkflowManager
            workflow = MetalSiteWorkflowManager(console=self.console, processor=self.processor)

            # Provide the RedoxSite directly
            workflow.provided_redox_site = site

            # Create PreprocessingResult with our data
            # ff_data is None - we get types/charges directly from prmtop
            workflow.preprocessing_result = PreprocessingResult(
                success=True,
                prepared_pdb=self._final_pdb,
                ff_data=None,  # Not needed - prmtop has all data
                type_assignments=self.type_assignments
            )

            # Create site-specific output directory
            site_output_dir = self._output_dir / f"site_{idx + 1}"
            site_output_dir.mkdir(parents=True, exist_ok=True)

            # Re-running a site already parameterized (to correct its charges,
            # say) should replace its entry rather than allocate a second set of
            # types on top of the high-water mark and strand the first.
            site_metal_start, site_ligand_start = metal_type_offset, ligand_type_offset
            reused = self._offer_prior_type_reuse(site_output_dir, site_id)
            if reused is not None:
                site_metal_start, site_ligand_start = reused

            # Run Step 1 — thread the running M*/Y* offsets so this site's
            # types don't collide with prior sites'.
            result = workflow._run_step1(
                residue_name=site_id,
                residues=[],  # Not used when provided_redox_site is set
                output_dir=site_output_dir,
                interactive=self._interactive,
                metal_type_start=site_metal_start,
                ligand_type_start=site_ligand_start,
            )

            # Advance the offsets by however many types this site consumed so
            # the next metal site continues the numbering instead of restarting.
            # max(), not assignment: a site that REUSED its earlier names ends
            # below the high-water mark, and taking its end verbatim would hand
            # the next site names that are already spoken for.
            metal_type_offset = max(
                metal_type_offset,
                getattr(workflow, "type_offset_metal_end", metal_type_offset))
            ligand_type_offset = max(
                ligand_type_offset,
                getattr(workflow, "type_offset_ligand_end", ligand_type_offset))

            results.append(result)

            if result.get("success"):
                stats = result.get("atom_summary", {})
                self.console.print(f"[green]✓ Site {idx + 1}: {stats.get('total_atoms', 0)} atoms typed, "
                                   f"{stats.get('renamed_atoms', 0)} renamed to M*/Y*[/green]")
            else:
                self.console.print(f"[red]✗ Site {idx + 1} failed: {result.get('message', 'Unknown error')}[/red]")

        # Carry the high-water mark, so a metal site parameterized later in this
        # same ProPrep session continues the numbering instead of restarting.
        # A fresh session recovers it from the fingerprints on disk instead.
        if self.workspace:
            self.workspace.set("mcpb_type_offsets",
                               {"metal": metal_type_offset,
                                "ligand": ligand_type_offset})

        # Summarize
        successful = sum(1 for r in results if r.get("success"))
        total_renamed = sum(r.get("atom_summary", {}).get("renamed_atoms", 0) for r in results)

        if not results:
            self.console.print("[yellow]No metal sites required MCPB atom typing.[/yellow]")

        # Count against what was actually attempted: with a selection in play,
        # "2/2 sites" is the honest denominator, not the structure's total.
        attempted = len(results) or metal_site_count
        selection_note = ""
        if selected_ids and attempted < metal_site_count:
            selection_note = (f" ({metal_site_count - attempted} site(s) not "
                              f"selected for this run)")

        return {
            "summary": (f"{successful}/{attempted} metal site(s) processed, "
                        f"{total_renamed} atoms renamed to M*/Y*{selection_note}"),
            "results": results
        }

    def _build_mcpb_data_from_workspace(self) -> None:
        """
        Build type_assignments from workspace for MCPB Step 1.

        The prmtop IS the force field applied to the system. We get atom types
        and charges directly from it, so there's no need to load leaprcs.

        Note: Metals aren't in prmtop (removed before tLEaP, reinserted after).
        Metal charges/spins are collected from the user in Step 1.

        Reads:
        - preprocessing_atom_data: {(chain, resid, resname, atom_name): {'type': str, 'charge': float}}
        - redox_sites (self.redox_sites)

        Sets:
        - self.type_assignments: Dict[coords, AtomTypeAssignment]
        """
        self.console.print("[grey50]Building MCPB data from workspace...[/grey50]")

        # Get atom data (types + charges) from prmtop extraction
        prmtop_data = self.workspace.get("preprocessing_atom_data", {}) if self.workspace else {}
        self.console.print(f"[grey50]  Found {len(prmtop_data)} atoms from prmtop (types + charges)[/grey50]")

        # Build type_assignments by matching prmtop data to redox site atoms
        self.type_assignments = {}

        # Build lookup from (resid, resname, atom_name) -> {'type': ..., 'charge': ...}
        # (tLEaP strips chain IDs, so we ignore chain)
        data_lookup = {}
        for (chain, resid, resname, atom_name), atom_info in prmtop_data.items():
            key = (resid, resname, atom_name)
            data_lookup[key] = atom_info

        # Import METALS for metal identification
        from proprep.structure_prep.comprehensive_redox_detector import METALS

        # Get redox sites from instance variable or workspace (for resume case)
        redox_sites = self.redox_sites
        if not redox_sites and self.workspace:
            redox_sites = _ensure_redox_site_objects(self.workspace.get("detected_redox_sites", []))

        # Match each atom in redox sites to prmtop data
        for site in redox_sites:
            # Build metal_coords for this site from centers
            # - METAL_ION: center coords ARE the metal
            # - ORGANOMETALLIC_COFACTOR: find metal atoms in same residue
            metal_coords = set()
            for center in site.centers:
                if hasattr(center, 'center_type'):
                    if center.center_type.value == 'metal_ion':
                        metal_coords.add(center.coords)
                    elif center.center_type.value == 'organometallic_cofactor':
                        # Find metal atoms in this residue
                        for atom in site.atoms:
                            if (atom.chain == center.chain and
                                atom.resid == center.resid and
                                atom.resname == center.resname and
                                atom.element.upper() in METALS):
                                metal_coords.add(atom.coords)

            # Build ligand_coords from coordinate bonds
            # If a bond has chemical_type == 'coordinate', one end is metal, other is ligand
            ligand_coords = set()
            for bond in getattr(site, 'bonds', []):
                if getattr(bond, 'chemical_type', '') == 'coordinate':
                    if bond.atom1_coords in metal_coords:
                        ligand_coords.add(bond.atom2_coords)
                    elif bond.atom2_coords in metal_coords:
                        ligand_coords.add(bond.atom1_coords)

            # Now create assignments for each atom in this site
            for atom in site.atoms:
                # Try to find data in prmtop
                lookup_key = (atom.resid, atom.resname, atom.atom_name)
                atom_info = data_lookup.get(lookup_key)

                if not atom_info:
                    # Try without exact resid match (prmtop may renumber)
                    for (resid, resname, atom_name), info in data_lookup.items():
                        if resname == atom.resname and atom_name == atom.atom_name:
                            atom_info = info
                            break

                # Extract type and charge from prmtop data
                # Metals won't be in prmtop - their charge is set later by user
                if atom_info:
                    original_type = atom_info['type']
                    charge = atom_info['charge']
                else:
                    # Fallback for metals or missing atoms: a withheld cluster
                    # residue is not in the prmtop at all.
                    #
                    # The element symbol is a fine placeholder for a metal or a
                    # bridging sulfide -- 'MO'/'FE'/'S' are not Amber types, so
                    # they read as provisional and those atoms are renamed to
                    # M*/Y* anyway. HYDROGEN is the exception: 'H' IS a valid
                    # Amber type, the amide/amine one (r* 0.6000), so a hydroxo
                    # proton silently acquired the wrong nonbonded terms where
                    # the hydroxyl convention 'HO' is 0.0000. Type it from the
                    # atom it is bonded to.
                    if (atom.element or '').strip().upper() == 'H':
                        from proprep.forcefield_prep.mcpb.atom_typer import (
                            hydrogen_type_from_neighbors,
                        )
                        original_type = hydrogen_type_from_neighbors(
                            atom.coords,
                            [(a.coords, a.element) for a in site.atoms
                             if (a.chain, a.resid) == (atom.chain, atom.resid)],
                        )
                    else:
                        original_type = atom.element
                    charge = None  # Will be set by user for metals

                # Determine is_center and is_metal_ligand from our computed sets
                is_center = atom.coords in metal_coords
                is_metal_ligand = atom.coords in ligand_coords

                assignment = AtomTypeAssignment(
                    coords=atom.coords,
                    chain=atom.chain,
                    resname=atom.resname,
                    resid=atom.resid,
                    atom_name=atom.atom_name,
                    element=atom.element,
                    original_type=original_type,
                    renamed_type=original_type,
                    charge=charge,
                    is_center=is_center,
                    is_metal_ligand=is_metal_ligand,
                )
                self.type_assignments[atom.coords] = assignment

        # Also include ALL atoms from prepared structure (for gap residues like PHE 91)
        # This ensures gap-filling residues have charges from prmtop
        prepared_pdb = self.workspace.get("prepared_pdb") if self.workspace else None
        if prepared_pdb and Path(prepared_pdb).exists():
            from Bio.PDB import PDBParser
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("prepared", prepared_pdb)

            all_atoms_added = 0
            for model in structure:
                for chain in model:
                    for residue in chain:
                        resid = residue.get_id()[1]
                        resname = residue.get_resname().strip()
                        for atom in residue:
                            coords = tuple(round(x, 3) for x in atom.coord)

                            # Skip if already in type_assignments (RedoxSite atoms)
                            if coords in self.type_assignments:
                                continue

                            # Look up in prmtop data
                            lookup_key = (resid, resname, atom.name.strip())
                            atom_info = data_lookup.get(lookup_key)

                            if atom_info:
                                assignment = AtomTypeAssignment(
                                    coords=coords,
                                    chain=chain.id,
                                    resname=resname,
                                    resid=resid,
                                    atom_name=atom.name.strip(),
                                    element=atom.element.strip() if atom.element else '',
                                    original_type=atom_info['type'],
                                    renamed_type=atom_info['type'],
                                    charge=atom_info['charge'],
                                    is_center=False,
                                    is_metal_ligand=False,
                                )
                                self.type_assignments[coords] = assignment
                                all_atoms_added += 1

            if all_atoms_added > 0:
                self.console.print(f"[grey50]  Added {all_atoms_added} additional atoms from prepared structure[/grey50]")

        # Count metals (atoms without prmtop data - need user input for charge)
        metals_needing_charge = sum(1 for a in self.type_assignments.values() if a.charge is None)
        if metals_needing_charge > 0:
            self.console.print(f"[grey50]  {metals_needing_charge} metal atom(s) need charge from user[/grey50]")

        self.console.print(f"[grey50]  Built {len(self.type_assignments)} type assignments[/grey50]")

    def _checklist_mcpb_2_bonded(self) -> dict:
        """
        Checklist handler for mcpb-2: Bonded Parameters (Seminario Method).

        For each redox site:
        1. Check if Gaussian output files exist (small_freq.log, small_freq.fchk)
        2. If they exist, run Seminario method to extract force constants
        3. If not, show checkpoint message
        """
        from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager

        # Get redox sites from instance variable or workspace (for resume case)
        redox_sites = self.redox_sites
        if not redox_sites and self.workspace:
            redox_sites = _ensure_redox_site_objects(self.workspace.get("detected_redox_sites", []))
            if redox_sites:
                self.redox_sites = redox_sites  # Cache for future use
                self.console.print(f"[grey50]Restored {len(redox_sites)} redox site(s) from workspace[/grey50]")

        if not redox_sites:
            self.console.print("[red]No redox sites available for MCPB Step 2[/red]")
            return {"summary": "Failed - no redox sites"}

        results = []
        all_complete = True
        stale_sites = []

        # Helper to get site_id whether site is a dict or RedoxSite object
        def get_site_id(site, idx):
            if isinstance(site, dict):
                return site.get('site_id', f'site_{idx + 1}')
            return site.site_id

        for idx, site in enumerate(redox_sites):
            # Only metal sites go through MCPB (organic cofactors are handled by
            # the Small Molecule Parameterizer). Skip preserves idx so the ZN
            # site keeps the same site_N directory as in MCPB Step 1.
            if not _redox_site_has_metal(site):
                continue

            site_output_dir = self._output_dir / f"site_{idx + 1}"
            step1_dir = site_output_dir / "models"

            # Check for Gaussian output files
            small_fchk = step1_dir / "small_freq.fchk"
            small_log = step1_dir / "small_freq.log"

            if not small_fchk.exists() or not small_log.exists():
                self.console.print(f"[yellow]Site {idx + 1}: Gaussian output not found[/yellow]")
                self.console.print(f"  Expected: {step1_dir}/small_freq.fchk and small_freq.log")
                all_complete = False
                continue

            site_id = get_site_id(site, idx)
            self.console.print(f"\n[bold cyan]Processing site {idx + 1}/{len(redox_sites)}: {site_id}[/bold cyan]")

            # Step 12 rebuilds the small model too, and Gaussian is run on it
            # manually; a leftover .log gives force constants for the old one.
            stale = self._stale_gaussian_output(step1_dir, stem="small_freq")
            if stale:
                self.console.print(
                    f"[red]✗ Site {idx + 1}: the Gaussian output is for a "
                    f"different model than the input beside it[/red]")
                self.console.print(f"  {stale}")
                self.console.print(
                    f"  The Hessian in that log is for the earlier model, so "
                    f"the Seminario force constants would not describe this "
                    f"site.")
                self.console.print(
                    f"  [yellow]Re-run Gaussian on {step1_dir}/small_freq.gjf, "
                    f"then run this step again.[/yellow]")
                results.append({"success": False, "site_id": site_id,
                                "message": f"stale Gaussian output ({stale})"})
                stale_sites.append(idx + 1)
                all_complete = False
                continue

            self.console.print(f"[green]✓ Found Gaussian output: small_freq.fchk, small_freq.log[/green]")

            # Create MetalSiteWorkflowManager
            workflow = MetalSiteWorkflowManager(console=self.console, processor=self.processor)
            workflow.provided_redox_site = site

            # Load Step 1 results if available
            workflow_state_file = site_output_dir / "workflow_state.json"
            if workflow_state_file.exists():
                import json
                with open(workflow_state_file) as f:
                    state = json.load(f)
                    workflow.step_results = state.get("step_results", {})

            # Step 2a: Generate pre-frcmod with all parameters (NON/YES markers)
            step2a_result = workflow._run_step2a(
                residue_name=site_id,
                residues=[],
                output_dir=site_output_dir,
                interactive=self._interactive
            )

            if not step2a_result.get("success"):
                self.console.print(f"[red]✗ Site {idx + 1} Step 2a failed: {step2a_result.get('message', 'Unknown error')}[/red]")
                results.append(step2a_result)
                continue

            # Step 2b: Run Seminario method and merge with pre-frcmod
            result = workflow._run_step2b(
                residue_name=site_id,
                residues=[],
                output_dir=site_output_dir,
                interactive=self._interactive
            )

            results.append(result)

            if result.get("success"):
                stats = result.get("statistics", {})
                bond_count = stats.get("n_bonds", 0)
                angle_count = stats.get("n_angles", 0)
                self.console.print(f"[green]✓ Site {idx + 1}: {bond_count} bonds, {angle_count} angles parameterized[/green]")
            else:
                self.console.print(f"[red]✗ Site {idx + 1} failed: {result.get('message', 'Unknown error')}[/red]")

        if not all_complete:
            if stale_sites:
                listed = ", ".join(str(n) for n in stale_sites)
                self.console.print(
                    f"\n[yellow]Site(s) {listed}: the Gaussian output describes "
                    f"a different model than the input beside it. Re-run "
                    f"Gaussian for those sites, then resume.[/yellow]")
                return {"summary": f"Stale Gaussian output for site(s) {listed}",
                        "checkpoint": True}
            self.console.print("\n[yellow]Some sites missing Gaussian output. Run Gaussian, then resume.[/yellow]")
            return {"summary": "Waiting for Gaussian output", "checkpoint": True}

        successful = sum(1 for r in results if r.get("success"))
        return {"summary": f"{successful}/{len(redox_sites)} sites parameterized"}

    @staticmethod
    def _parse_gaussian_geometry(lines) -> List[Tuple[str, float, float, float]]:
        """``(element, x, y, z)`` per atom from Gaussian Cartesian input lines.

        Both the .gjf geometry block and the log's ``Symbolic Z-matrix`` echo
        use the same layout, with an optional frozen-atom flag between the
        element and the coordinates::

            N   -1  -36.96000000  -14.48000000  -45.92400000
            N       -36.96        -14.48        -45.924
        """
        geometry = []
        for line in lines:
            parts = line.split()
            if len(parts) == 5:
                symbol, coords = parts[0], parts[2:5]
            elif len(parts) == 4:
                symbol, coords = parts[0], parts[1:4]
            else:
                continue
            try:
                x, y, z = (float(c) for c in coords)
            except ValueError:
                continue
            geometry.append((symbol.upper(), x, y, z))
        return geometry

    @classmethod
    def _gjf_geometry(cls, gjf: Path):
        """The model a Gaussian input describes, or None if unreadable.

        Layout: route, blank, title, blank, "<charge> <mult>", the geometry,
        then a blank line. The ReadRadii entries after that blank line look
        exactly like atom lines ("Fe 1.383"), so the block has to end there.
        """
        try:
            blanks = 0
            in_geometry = False
            body = []
            with open(gjf, errors="ignore") as fh:
                for line in fh:
                    if not line.strip():
                        if in_geometry:
                            break
                        blanks += 1
                        continue
                    if blanks < 2:
                        continue
                    if not in_geometry:
                        # The charge/multiplicity line opens the geometry.
                        if len(line.split()) == 2:
                            in_geometry = True
                        continue
                    body.append(line)
            return cls._parse_gaussian_geometry(body) or None
        except OSError:
            return None

    @classmethod
    def _log_input_geometry(cls, log: Path):
        """The model Gaussian was given, from its ``Symbolic Z-matrix`` echo.

        This is the input verbatim. The ``Input orientation`` table further
        down is NOT interchangeable with it -- Gaussian re-centers and
        reorients the molecule there, so comparing it against the .gjf reports
        differences that are pure rigid-body motion.
        """
        try:
            body = []
            found = False
            with open(log, errors="ignore") as fh:
                for line in fh:
                    if not found:
                        if "Symbolic Z-matrix:" in line:
                            found = True
                        continue
                    if not line.strip():
                        break
                    if "Charge" in line and "Multiplicity" in line:
                        continue
                    body.append(line)
            return cls._parse_gaussian_geometry(body) or None
        except OSError:
            return None

    @classmethod
    def _stale_gaussian_output(cls, models_dir: Path, stem: str = "large_resp",
                               tolerance: float = 1e-3):
        """Why ``<stem>.log`` is superseded by ``<stem>.gjf``, or None.

        Compares the model in the input against the one the log says Gaussian
        was given. Content, not timestamps: mtimes are not evidence about what
        a file contains, and they are rewritten by copying a log back from a
        cluster or checking files out again. Re-running step 12 also rewrites
        an input that is byte-identical to the one that ran, which a timestamp
        test reports as stale when nothing has changed.

        Applies to both QM models, which go stale for the same reason -- step
        12 rebuilds them together, and Gaussian is run manually on each:

        - ``small_freq`` feeds the Seminario force constants (mcpb-2)
        - ``large_resp`` feeds the RESP charges (mcpb-3)

        The failure this prevents is quiet. A superseded log describes the old
        model, so the Hessian or ESP taken from it belongs to a different
        molecule, and the output still looks like force constants or charges.
        Re-running Gaussian is the only repair, so the step refuses.
        """
        models_dir = Path(models_dir)
        log = models_dir / f"{stem}.log"
        gjf = models_dir / f"{stem}.gjf"
        if not log.exists() or not gjf.exists():
            return None

        want = cls._gjf_geometry(gjf)
        got = cls._log_input_geometry(log)
        if want is None or got is None:
            return None

        if len(want) != len(got):
            return (f"large_resp.gjf has {len(want)} atoms but large_resp.log "
                    f"was run on {len(got)}")

        want_elements = [a[0] for a in want]
        got_elements = [a[0] for a in got]
        if want_elements != got_elements:
            changed = sum(1 for a, b in zip(want_elements, got_elements) if a != b)
            return (f"{stem}.gjf and {stem}.log agree on {len(want)} "
                    f"atoms but {changed} differ in element")

        worst = 0.0
        for (_e, *a), (_f, *b) in zip(want, got):
            worst = max(worst, max(abs(p - q) for p, q in zip(a, b)))
        if worst > tolerance:
            return (f"{stem}.gjf and {stem}.log describe the same atoms "
                    f"at different coordinates (up to {worst:.3f} A apart)")

        return None

    @staticmethod
    def _esp_charge_multiplicity(models_dir: Path):
        """(charge, multiplicity) the ESP was actually computed with, or (None, 1).

        RESP fits point charges to a specific electrostatic potential, so its
        total-charge constraint has to be the charge that potential was computed
        under. Any other value is not a worse fit — it is a fit to a different
        molecule, and RESP will spread the difference over the atoms rather than
        fail.

        Read from the Gaussian log first (what the calculation actually used),
        then the .gjf (what it was asked to use). Both live beside the model, so
        neither can be confused with another site's.

        That preference only holds while the log corresponds to the .gjf beside
        it; call ``_stale_gaussian_output`` first, or an obsolete log will hand
        back the charge of a model that has since been regenerated.
        """
        log = Path(models_dir) / "large_resp.log"
        if log.exists():
            try:
                pattern = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")
                with open(log, errors="ignore") as fh:
                    for line in fh:
                        m = pattern.search(line)
                        if m:
                            return int(m.group(1)), int(m.group(2))
            except OSError:
                pass

        gjf = Path(models_dir) / "large_resp.gjf"
        if gjf.exists():
            try:
                # Route, blank, title, blank, then "<charge> <multiplicity>".
                blanks = 0
                with open(gjf, errors="ignore") as fh:
                    for line in fh:
                        if not line.strip():
                            blanks += 1
                            continue
                        if blanks >= 2:
                            parts = line.split()
                            if len(parts) == 2:
                                try:
                                    return int(parts[0]), int(parts[1])
                                except ValueError:
                                    return None, 1
                            return None, 1
            except OSError:
                pass

        return None, 1

    def _checklist_mcpb_3_resp(self) -> dict:
        """
        Checklist handler for mcpb-3: RESP Charges.

        For each redox site:
        1. Check if Gaussian ESP output exists (large_resp.log in step1)
        2. If it exists, run steps 3B, 3C, 3D to fit RESP charges
        3. If not, show checkpoint message

        The integrated workflow generates large_resp.gjf in step1, not step3,
        so we synthesize step_3a results to point to step1 directory.
        """
        from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager

        # Get redox sites from instance variable or workspace (for resume case)
        redox_sites = self.redox_sites
        if not redox_sites and self.workspace:
            redox_sites = _ensure_redox_site_objects(self.workspace.get("detected_redox_sites", []))
            if redox_sites:
                self.redox_sites = redox_sites
                self.console.print(f"[grey50]Restored {len(redox_sites)} redox site(s) from workspace[/grey50]")

        if not redox_sites:
            self.console.print("[red]No redox sites available for MCPB Step 3[/red]")
            return {"summary": "Failed - no redox sites"}

        results = []
        all_complete = True
        stale_sites = []

        # Helper to get site_id whether site is a dict or RedoxSite object
        def get_site_id(site, idx):
            if isinstance(site, dict):
                return site.get('site_id', f'site_{idx + 1}')
            return site.site_id

        for idx, site in enumerate(redox_sites):
            # Only metal sites go through MCPB (see MCPB Step 1/2). Skip preserves
            # idx so the ZN site keeps the same site_N directory.
            if not _redox_site_has_metal(site):
                continue

            site_output_dir = self._output_dir / f"site_{idx + 1}"
            step1_dir = site_output_dir / "models"

            # Check for Gaussian ESP output file (in models, not charge_fit)
            large_resp_log = step1_dir / "large_resp.log"
            large_pdb = step1_dir / "large.pdb"

            if not large_resp_log.exists():
                self.console.print(f"[yellow]Site {idx + 1}: Gaussian ESP output not found[/yellow]")
                self.console.print(f"  Expected: {step1_dir}/large_resp.log")
                all_complete = False
                continue

            site_id = get_site_id(site, idx)
            self.console.print(f"\n[bold cyan]Processing site {idx + 1}/{len(redox_sites)}: {site_id}[/bold cyan]")

            # Present is not the same as current: step 1 may have been re-run
            # since Gaussian last was, leaving a log for a superseded model.
            stale = self._stale_gaussian_output(step1_dir)
            if stale:
                self.console.print(
                    f"[red]✗ Site {idx + 1}: the Gaussian output is for a "
                    f"different model than the input beside it[/red]")
                self.console.print(f"  {stale}")
                self.console.print(
                    f"  The ESP in that log was computed for the earlier model, "
                    f"so fitting to it would produce charges for a molecule "
                    f"this site no longer describes.")
                self.console.print(
                    f"  [yellow]Re-run Gaussian on {step1_dir}/large_resp.gjf, "
                    f"then run this step again.[/yellow]")
                results.append({"success": False, "site_id": site_id,
                                "message": f"stale Gaussian output ({stale})"})
                stale_sites.append(idx + 1)
                all_complete = False
                continue

            self.console.print(f"[green]✓ Found Gaussian ESP output: large_resp.log[/green]")

            # Create MetalSiteWorkflowManager
            workflow = MetalSiteWorkflowManager(console=self.console, processor=self.processor)
            workflow.provided_redox_site = site

            # Load workflow state if available
            workflow_state_file = site_output_dir / "workflow_state.json"
            if workflow_state_file.exists():
                import json
                with open(workflow_state_file) as f:
                    state = json.load(f)
                    workflow.step_results = state.get("step_results", {})

            # Charge and multiplicity for the RESP fit MUST be the ones the ESP
            # was computed with, so they are read from THIS site's Gaussian
            # artifacts rather than from step_results.
            #
            # step_results is now partitioned per site, so step_1 is this
            # site's. Reading the Gaussian artifacts stays the primary source
            # anyway: they record what the calculation ACTUALLY used, which
            # step_1 cannot know if the .gjf was edited by hand or the run was
            # repeated at a different charge.
            #
            # Kept also because it is what caught the original defect, when
            # mcpb_step_results was one dict every site shared: a -1 Fe2S2
            # model fitted against site 2's -3 constraint put +5.6 on a metal.
            large_charge, large_mult = self._esp_charge_multiplicity(step1_dir)

            step1_results = workflow.step_results.get("step_1", {})
            qm_params = step1_results.get("qm_parameters", {}).get("large_model", {})

            if large_charge is None:
                # No Gaussian artifact to read: fall back to step_1, then to the
                # large.pdb REMARK (the SUGGESTED charge, which the user may
                # have overridden at the prompt).
                large_charge = qm_params.get("charge", None)
                large_mult = qm_params.get("multiplicity", 1)
                if large_charge is None:
                    large_charge = 0
                    try:
                        with open(large_pdb) as f:
                            for line in f:
                                if line.startswith("REMARK") and "Total charge:" in line:
                                    charge_str = line.split("Total charge:")[1].strip()
                                    large_charge = round(float(charge_str))
                                    break
                    except Exception:
                        pass
                self.console.print(
                    f"[yellow]  Could not read the charge from this site's Gaussian "
                    f"input/output; using {large_charge:+d} from "
                    f"{'step 1' if qm_params else 'the large.pdb REMARK'}. Verify it "
                    f"matches the ESP calculation.[/yellow]")
            else:
                stored = qm_params.get("charge")
                if stored is not None and stored != large_charge:
                    # Exactly the symptom above — worth naming rather than
                    # silently preferring the right one.
                    self.console.print(
                        f"[yellow]  Step-1 records charge {stored:+d} for this site "
                        f"but its ESP was computed at {large_charge:+d}; using "
                        f"{large_charge:+d} to match the ESP.[/yellow]")
                self.console.print(
                    f"[grey50]  ESP charge {large_charge:+d}, multiplicity "
                    f"{large_mult} (from this site's Gaussian files)[/grey50]")

            # Synthesize step_3a results from step1 data
            # (In integrated workflow, ESP input is generated in step1, not step3)
            # Step 3B expects: step_3a["output_dir"], step_3a["large_pdb"]
            # And looks for large_resp.log in output_dir
            workflow.step_results["step_3a"] = {
                "output_dir": str(step1_dir),  # ESP files are in step1
                "large_pdb": str(large_pdb),
                "charge": large_charge,
                "multiplicity": large_mult,
                "status": "input_generated"
            }

            # Prompt for cross-residue charge equivalence before generating resp.in.
            # Gives the user explicit control over which same-resname ligands
            # (e.g., 4 Cys on a Zn site, 2 axial His on a heme) should share
            # side-chain charges. Stored on the workflow instance for
            # _run_step3b to read.
            workflow.cross_residue_eq_groups = self._prompt_cross_residue_equivalence(
                site, idx + 1
            )

            # Run Step 3B: RESP Input Generation
            result_3b = workflow._run_step3b(
                residue_name=site_id,
                output_dir=site_output_dir,
                interactive=self._interactive
            )

            if not result_3b.get("success"):
                self.console.print(f"[red]✗ Step 3B failed: {result_3b.get('message', 'Unknown error')}[/red]")
                results.append(result_3b)
                continue

            # Run Step 3C: RESP Execution
            result_3c = workflow._run_step3c(
                residue_name=site_id,
                output_dir=site_output_dir
            )

            if not result_3c.get("success"):
                self.console.print(f"[red]✗ Step 3C failed: {result_3c.get('message', 'Unknown error')}[/red]")
                results.append(result_3c)
                continue

            # Run Step 3D: Mol2 File Generation
            result_3d = workflow._run_step3d(
                residue_name=site_id,
                output_dir=site_output_dir
            )

            if not result_3d.get("success"):
                self.console.print(f"[red]✗ Step 3D failed: {result_3d.get('message', 'Unknown error')}[/red]")
                results.append(result_3d)
                continue

            results.append(result_3d)
            self.console.print(f"[green]✓ Site {idx + 1}: RESP charges fitted and mol2 files generated[/green]")

        if not all_complete:
            if stale_sites:
                listed = ", ".join(str(n) for n in stale_sites)
                self.console.print(
                    f"\n[yellow]Site(s) {listed}: the Gaussian output describes "
                    f"a different model than the input beside it. Re-run "
                    f"Gaussian for those sites, then resume.[/yellow]")
                return {"summary": f"Stale Gaussian output for site(s) {listed}",
                        "checkpoint": True}
            self.console.print("\n[yellow]Some sites missing Gaussian ESP output. Run Gaussian, then resume.[/yellow]")
            return {"summary": "Waiting for Gaussian ESP output", "checkpoint": True}

        successful = sum(1 for r in results if r.get("success"))
        return {"summary": f"{successful}/{len(redox_sites)} sites RESP fitted"}

    def _prompt_cross_residue_equivalence(self, site, site_number: int):
        """
        Ask the user which residues should share side-chain charges in RESP stage 2.

        Returns a list of (chain, resid) groups suitable for
        RESPInputGenerator._expand_cross_residue_equivalences. Returns an
        empty list when the user skips or no candidate groupings exist.

        Same-resname is required within each group (atom names must match
        for the equivalence to be well-defined). Single-residue groups
        are silently dropped.
        """
        from collections import OrderedDict
        from rich.panel import Panel
        from rich.table import Table

        from proprep.utils.group_syntax import parse_group_assignments
        from proprep.utils.prompts import prompt_with_context

        # De-duplicate residues from the site's atom list, preserving order.
        residue_keys: "OrderedDict[Tuple[str, int], str]" = OrderedDict()
        for atom in getattr(site, "atoms", []) or []:
            key = (atom.chain, atom.resid)
            if key not in residue_keys:
                residue_keys[key] = atom.resname

        # Drop restrained ligands (e.g. nonbonded waters): they are held by an
        # MD restraint and kept as their standard water model (net-0, no RESP-fit
        # charges), so they are NOT cross-residue equivalence candidates.
        from proprep.forcefield_prep.metal_site_parameterizer import _collect_restrained_ligands
        _, restrained_resids = _collect_restrained_ligands(site)
        if restrained_resids:
            for key in list(residue_keys.keys()):
                if key[1] in restrained_resids:
                    del residue_keys[key]

        ordered_residues = list(residue_keys.items())
        if len(ordered_residues) < 2:
            return []

        # Build a default suggestion: for each resname class with 2+ instances,
        # group those instances together. Multiple classes → semicolon-separated
        # groups in default string.
        by_resname: Dict[str, List[int]] = {}
        for i, (_, resname) in enumerate(ordered_residues, start=1):
            by_resname.setdefault(resname, []).append(i)

        suggestion_groups = [indices for indices in by_resname.values() if len(indices) >= 2]
        if not suggestion_groups:
            # Nothing to suggest; user can still type manually but most won't.
            self.console.print(
                "[grey50]No same-resname residue groups detected; skipping "
                "cross-residue equivalence prompt for this site.[/grey50]"
            )
            return []

        def _format_group(indices: List[int]) -> str:
            indices = sorted(indices)
            # Compress consecutive runs into ranges
            ranges: List[str] = []
            start = indices[0]
            prev = indices[0]
            for x in indices[1:]:
                if x == prev + 1:
                    prev = x
                    continue
                ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
                start = prev = x
            ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
            return ",".join(ranges)

        suggestion = " ".join(
            f"{gid}:{_format_group(indices)}"
            for gid, indices in enumerate(suggestion_groups, start=1)
        )

        # Show residue list
        table = Table(title=f"RESP Charge Equivalence — Site {site_number}", expand=False)
        table.add_column("#", style="cyan")
        table.add_column("Resname", style="green")
        table.add_column("Location", style="yellow")
        for i, ((chain, resid), resname) in enumerate(ordered_residues, start=1):
            table.add_row(str(i), resname, f"{chain}:{resid}")
        self.console.print(table)

        self.console.print(Panel.fit(
            "Group residues that should share side-chain charges across the group.\n"
            "Within a group, atoms with the same name will be equivalenced across\n"
            "residues. Residues in a group must share the same residue name.\n\n"
            "Only for genuinely symmetry-equivalent ligands (e.g. 4 Cys on a Zn,\n"
            "2 axial His on a heme). Different residues or inequivalent metals\n"
            "must NOT be forced to share charges.\n\n"
            "Format: [bold]<group_id>:<residue_list>[/bold] (same as site-grouping syntax)\n"
            "  [grey50]\"1:1-4\"[/grey50]         one group with all 4 ligands\n"
            "  [grey50]\"1:1,2 2:3,4\"[/grey50]    two groups of 2\n"
            "  [grey50]Enter[/grey50]             skip — DEFAULT (no cross-residue equivalence)",
            title="Instructions",
            border_style="grey50",
        ))

        # Equivalence is OFF by default: forcing charge equality is only correct
        # for symmetry-equivalent ligands. Show the same-resname grouping as an
        # opt-in suggestion the user can type in, NOT as the Enter default.
        if suggestion:
            self.console.print(
                f"[grey50]Suggestion, ONLY if these are symmetry-equivalent: "
                f"[bold]{suggestion}[/bold] — type it in to apply, or press Enter to skip.[/grey50]"
            )

        max_idx = len(ordered_residues)
        while True:
            response = prompt_with_context(
                self.processor,
                "Equivalence groups [default: skip]",
                default="",
                module="MCPB RESP",
                description="Cross-residue charge equivalence selection (opt-in)",
            ).strip()

            if response == "" or response.lower() in {"skip", "none", "no"}:
                return []

            parsed = parse_group_assignments(response, max_idx)
            if parsed is None:
                self.console.print(
                    f"[red]Invalid syntax. Expected '[group_id]:[residue_list]' with "
                    f"indices in 1..{max_idx}. Try again or press Enter to skip.[/red]"
                )
                continue

            # Validate: same resname per group, no overlaps, each group >= 2 residues.
            seen_residues = set()
            valid_groups: List[List[Tuple[str, int]]] = []
            error: Optional[str] = None

            for group_id, indices_0based in sorted(parsed.items()):
                if len(indices_0based) < 2:
                    # Silently drop singleton groups
                    continue

                resnames = {ordered_residues[i][1] for i in indices_0based}
                if len(resnames) != 1:
                    error = (
                        f"Group {group_id} contains residues with different resnames "
                        f"({sorted(resnames)}); residues in one group must share the "
                        f"same resname so atom-name correspondence works."
                    )
                    break

                duplicates = [i + 1 for i in indices_0based if i in seen_residues]
                if duplicates:
                    error = (
                        f"Residue(s) {duplicates} appear in multiple groups; each "
                        f"residue can belong to at most one equivalence group."
                    )
                    break

                seen_residues.update(indices_0based)
                valid_groups.append([ordered_residues[i][0] for i in indices_0based])

            if error:
                self.console.print(f"[red]{error}[/red]")
                continue

            if not valid_groups:
                # User entered something but it all collapsed (all singletons); treat as skip.
                return []

            # Echo selection for confirmation
            summary_lines = []
            for group_id, member_keys in enumerate(valid_groups, start=1):
                names = [
                    f"{ordered_residues[next(i for i, (k, _) in enumerate(ordered_residues) if k == mk)][1]} "
                    f"{mk[0]}:{mk[1]}"
                    for mk in member_keys
                ]
                summary_lines.append(f"  Group {group_id}: " + ", ".join(names))
            self.console.print(
                "[green]Cross-residue equivalence groups recorded:[/green]\n"
                + "\n".join(summary_lines)
            )
            return valid_groups

    def _checklist_mcpb_4_integration(self) -> dict:
        """Checklist handler for mcpb-4: Force Field Integration.

        Automates the MCPB integration workflow:
        1. Collects mol2/frcmod deliverables from step 1-3, per metal site
        2. Generates residue names unique across ALL sites (MCPB convention)
        3. Per site: prompts for site type / redox state / spin state and
           deposits one FF parameter library entry under
           ~/.proprep/forcefield_params/, plus one reuse transformer
        4. Renames residues in the prepared PDB
        5. Stores every site's files and types for tLEaP consumption

        Naming is structure-wide (site 2's Cys must not collide with site 1's)
        but the deposit is per-site: a library entry describes one site type,
        and a merged entry can be reused only on a structure that happens to
        carry every site it was built from.
        """
        from rich.panel import Panel
        from rich.table import Table
        from pathlib import Path

        from proprep.forcefield_prep.mcpb.integration_utils import (
            parse_fingerprint,
            generate_unique_residue_names,
            rename_pdb_residues,
            create_ff_library,
        )

        # ================================================================
        # A. Collect MCPB output files, one record per metal site
        # ================================================================
        # A library entry describes ONE site type, so the deliverables are kept
        # partitioned by site_* directory rather than flattened. Only residue
        # NAMING needs a cross-site view (section D), and it takes its union
        # from these records.
        def _site_index(path: Path) -> int:
            """Numeric suffix of a ``site_N`` directory.

            Sorting the glob lexically puts ``site_10`` between ``site_1`` and
            ``site_2``, which would misalign the dir->RedoxSite mapping below.
            """
            suffix = path.name.split("_", 1)[-1]
            return int(suffix) if suffix.isdigit() else 0

        site_dirs = sorted(self._output_dir.glob("site_*"), key=_site_index)
        if not site_dirs:
            return {"summary": "No site directories found"}

        site_records = []
        for site_dir in site_dirs:
            mol2_files = sorted((site_dir / "models").glob("*.mol2"))
            frcmod_files = sorted((site_dir / "bonded_params").glob("*_bonded.frcmod"))
            if not mol2_files and not frcmod_files:
                continue
            fingerprint = site_dir / "models" / "standard.fingerprint"
            assignments = site_dir / "models" / "atom_type_assignments.json"
            site_records.append({
                "dir": site_dir,
                "index": _site_index(site_dir),
                "mol2_files": mol2_files,
                "frcmod_files": frcmod_files,
                "fingerprint": str(fingerprint) if fingerprint.exists() else None,
                "assignments": str(assignments) if assignments.exists() else None,
            })

        if not site_records:
            return {"summary": "No MCPB parameter files found"}

        # Label each record with the site it came from. mcpb-1 names the
        # directory ``site_{idx+1}`` over the FULL redox-site list and skips
        # non-metal sites without creating a directory, so the numeric suffix
        # indexes redox_sites — enumeration order would drift past any skipped
        # organic cofactor.
        labelling_sites = self.redox_sites
        if not labelling_sites and self.workspace:
            labelling_sites = self.workspace.get("detected_redox_sites", [])
        for record in site_records:
            site_obj = None
            pos = record["index"] - 1
            if labelling_sites and 0 <= pos < len(labelling_sites):
                site_obj = labelling_sites[pos]
            site_id = None
            if isinstance(site_obj, dict):
                site_id = site_obj.get("site_id")
            elif site_obj is not None:
                site_id = getattr(site_obj, "site_id", None)
            record["site_obj"] = site_obj
            record["label"] = site_id or record["dir"].name

        # ================================================================
        # B. Display deliverables table
        # ================================================================
        table = Table(title="MCPB Deliverables", expand=False)
        table.add_column("Site", style="bold cyan")
        table.add_column("Type", style="cyan")
        table.add_column("File", style="green")
        table.add_column("Location", style="grey50")

        for record in site_records:
            for mol2 in record["mol2_files"]:
                table.add_row(record["label"], "mol2 (RESP charges)",
                              mol2.name, str(mol2.parent))
            for frcmod in record["frcmod_files"]:
                table.add_row(record["label"], "frcmod (bonded params)",
                              frcmod.name, str(frcmod.parent))

        self.console.print(table)

        # ================================================================
        # C. Site identity — deferred to section E
        # ================================================================
        # type/redox/spin are per-site questions (two sites in one protein can
        # legitimately differ in both), and they are easier to answer once the
        # residue names are settled. Asked in the per-site deposit loop.

        # ================================================================
        # D. Parse EVERY site's fingerprint and generate unique residue names
        # ================================================================
        # A multi-site protein has one fingerprint per site_* directory. Each
        # site's residues and M*/Y* atom-type entries are parsed onto its own
        # record — that slice is what gets deposited. The union of the residue
        # keys is kept alongside for the one job that is genuinely
        # structure-wide: generating residue names unique across all sites.
        parsed_records = [r for r in site_records if r["fingerprint"]]
        if not parsed_records:
            self.console.print("[red]No fingerprint files found in site directories[/red]")
            return {"summary": "Failed — no fingerprint files"}

        combined_residue_keys = []
        seen_keys = set()
        for record in parsed_records:
            fp_data = parse_fingerprint(record["fingerprint"], record["assignments"])
            record["residue_keys"] = list(fp_data["residues"].keys())
            record["atom_type_entries"] = list(fp_data["atom_type_entries"])
            for key in record["residue_keys"]:
                if key not in seen_keys:
                    seen_keys.add(key)
                    combined_residue_keys.append(key)

        if not combined_residue_keys:
            self.console.print("[red]No residues found in fingerprint(s)[/red]")
            return {"summary": "Failed — no residues in fingerprint"}

        # Restrained ligands (e.g. a metal-coordinated water held by an MD
        # distance restraint) stay in the QM model for correct electronics but
        # get no mol2 and no library unit — they load as plain TIP3P WAT. So
        # they must NOT be renamed to an MCPB name: a phantom WT1/WT2 in the
        # PDB would have no tLEaP unit to build from. Drop them before naming.
        from proprep.forcefield_prep.metal_site_parameterizer import _collect_restrained_ligands
        restraint_sites = self.redox_sites
        if not restraint_sites and self.workspace:
            restraint_sites = self.workspace.get("detected_redox_sites", [])
        restrained_resids = set()
        for _site in (restraint_sites or []):
            _, _resids = _collect_restrained_ligands(_site)
            # Normalize to int: fingerprint keys use int resids, but bond
            # residue_info may carry them as strings.
            for _r in _resids:
                try:
                    restrained_resids.add(int(_r))
                except (TypeError, ValueError):
                    pass
        if restrained_resids:
            dropped = [k for k in combined_residue_keys if k[0] in restrained_resids]
            combined_residue_keys = [k for k in combined_residue_keys if k[0] not in restrained_resids]
            for record in parsed_records:
                record["residue_keys"] = [k for k in record["residue_keys"]
                                          if k[0] not in restrained_resids]
            for (resid, resname) in dropped:
                self.console.print(
                    f"  [grey50]Keeping {resname} {resid} as-is "
                    f"(restrained ligand — loads as-is, not renamed)[/grey50]"
                )

        # generate_unique_residue_names threads an internal existing_names set, so
        # calling it once on the union yields names unique across ALL sites.
        residue_name_map = generate_unique_residue_names(combined_residue_keys)

        # First site's fingerprint still passed for logging/back-compat; the
        # merged atom_type_entries below are what create_ff_library actually uses.
        # Display proposed names, attributed to the site each residue belongs
        # to. A residue coordinating two metals (a bridging Cys) is listed
        # under both, and is deposited into both sites' libraries — each entry
        # has to stand alone to be reusable.
        sites_for_key = {}
        for record in parsed_records:
            for key in record["residue_keys"]:
                sites_for_key.setdefault(key, []).append(record["label"])

        name_table = Table(title="Proposed Residue Names", expand=False)
        name_table.add_column("Site", style="cyan")
        name_table.add_column("Original", style="yellow")
        name_table.add_column("New Name", style="green bold")
        name_table.add_column("ResID", style="grey50")

        for (resid, resname), new_name in residue_name_map.items():
            name_table.add_row(
                ", ".join(sites_for_key.get((resid, resname), ["-"])),
                resname, new_name, str(resid),
            )

        self.console.print(name_table)

        if not confirm_with_context(
            self.processor,
            "Accept these residue names?",
            default=True,
            module="MCPB Integration",
            description="Confirm unique residue name mapping",
        ):
            self.console.print("[yellow]Integration cancelled.[/yellow]")
            return {"summary": "Cancelled by user"}

        # ================================================================
        # E. Create one FF parameter library entry per site
        # ================================================================
        # Each site is deposited on its own so it can be reused on its own. A
        # single merged entry keyed by one site_type cannot be applied to a
        # structure that has only one of these sites without dragging in the
        # other's residues, atom types and frcmod.
        if len(parsed_records) > 1:
            self.console.print(
                f"\n[bold]{len(parsed_records)} metal sites — each is named and "
                f"deposited to the library separately.[/bold]"
            )

        claimed_identities = {}   # (type, redox, spin) -> site label that took it
        last_answers = {}         # carried forward as the next site's defaults
        deposits = []
        failed_sites = []

        for record in parsed_records:
            label = record["label"]
            site_keys = record["residue_keys"]
            if not site_keys:
                self.console.print(
                    f"  [yellow]Skipping {label}: no residues left to deposit "
                    f"(all restrained ligands).[/yellow]"
                )
                continue

            site_residue_map = {k: residue_name_map[k] for k in site_keys
                                if k in residue_name_map}
            if not site_residue_map:
                self.console.print(
                    f"  [yellow]Skipping {label}: no named residues.[/yellow]"
                )
                continue

            self.console.print()
            resnames_here = ", ".join(sorted({rn for (_r, rn) in site_residue_map}))
            self.console.print(
                f"[bold cyan]Site {label}[/bold cyan] "
                f"[grey50]({record['dir'].name}: {resnames_here})[/grey50]"
            )

            # Previous site's answers become this one's defaults: sites in one
            # protein are often the same type and state, and a 6-site run
            # should not mean 18 cold prompts.
            while True:
                site_type = prompt_with_context(
                    self.processor,
                    f"Name this site type for {label} (snake_case, e.g., zinc_his3_cys)",
                    default=last_answers.get("site_type"),
                    module="MCPB Integration",
                    description=f"Unique identifier for the metal site type in {record['dir'].name}",
                )
                redox_state = prompt_with_context(
                    self.processor,
                    f"Redox state name for {label} (e.g., oxidized, reduced)",
                    default=last_answers.get("redox_state", "default"),
                    module="MCPB Integration",
                    description="Redox/oxidation state of this metal center",
                )
                spin_state = prompt_with_context(
                    self.processor,
                    f"Spin state name for {label} (e.g., high_spin, low_spin)",
                    default=last_answers.get("spin_state", "default"),
                    module="MCPB Integration",
                    description="Electronic spin state of this metal center",
                )

                identity = (site_type, redox_state, spin_state)
                if identity not in claimed_identities:
                    break
                # Two sites cannot share a library key even when they are
                # chemically equivalent: residue names are unique across the
                # whole structure, so site 1 holds CY1 where site 2 holds CY5.
                # promote_state overwrites a repeated key, which would silently
                # drop the earlier site's residue names from the entry.
                self.console.print(
                    f"  [yellow]{site_type}/{redox_state}/{spin_state} was already "
                    f"used for site {claimed_identities[identity]}. Each site needs "
                    f"its own library key — the residue names differ between "
                    f"sites, so depositing both here would lose one.[/yellow]"
                )

            claimed_identities[identity] = label
            last_answers = {"site_type": site_type, "redox_state": redox_state,
                            "spin_state": spin_state}

            description = f"MCPB-parameterized {site_type} metal site"

            # The MCPB bonded frcmod covers only the metal shell (M*/Y* types).
            # Any organic ligand in the site (e.g. E4Z) needs its own GAFF
            # frcmod deposited too, or its atoms type at reuse with no
            # vdW/torsion params. Scoped to this site's own ligands.
            site_resnames = {resname for (_resid, resname) in site_residue_map}
            ligand_frcmods = self._collect_ligand_frcmods(site_resnames)
            if ligand_frcmods:
                self.console.print(
                    f"  [grey50]Including {len(ligand_frcmods)} ligand GAFF frcmod(s): "
                    f"{', '.join(sorted(ligand_frcmods))}[/grey50]"
                )

            self.console.print("  [bold]Creating FF parameter library...[/bold]")
            # Deposits are per-site now, so each can fail on its own — a site
            # whose bonded frcmod is missing raises here. Report it and keep
            # going, or one incomplete site would cost every other site its
            # already-computed parameters.
            try:
                lib_result = create_ff_library(
                    site_type=site_type,
                    description=description,
                    mol2_files=[str(m) for m in record["mol2_files"]],
                    frcmod_files=[str(f) for f in record["frcmod_files"]],
                    fingerprint_path=record["fingerprint"],
                    assignments_path=record["assignments"],
                    residue_name_map=site_residue_map,
                    redox_state=redox_state,
                    spin_state=spin_state,
                    # This site's own M*/Y* types only. Per-site type offsets
                    # already keep them disjoint, so an entry reused on its own
                    # carries exactly the types its own parameters reference.
                    atom_type_entries=record["atom_type_entries"],
                    extra_frcmod_files=list(ligand_frcmods.values()),
                )
            except Exception as exc:
                self.console.print(
                    f"  [red]✗ Could not deposit {label}: {exc}[/red]"
                )
                failed_sites.append(label)
                continue

            record["site_type"] = site_type
            record["redox_state"] = redox_state
            record["spin_state"] = spin_state
            record["lib_result"] = lib_result
            record["residue_name_map"] = site_residue_map
            deposits.append(record)

            self.console.print(f"  [green]✓[/green] Library: [grey50]{lib_result['library_path']}[/grey50]")
            self.console.print(f"  [green]✓[/green] Metadata: [grey50]{lib_result['metadata_path']}[/grey50]")

        if not deposits:
            self.console.print("[red]No sites were deposited to the library[/red]")
            if failed_sites:
                return {"summary": f"Failed — no site deposited "
                                   f"({', '.join(failed_sites)})"}
            return {"summary": "Failed — nothing to deposit"}

        if failed_sites:
            self.console.print(
                f"\n[yellow]Deposited {len(deposits)} of "
                f"{len(deposits) + len(failed_sites)} sites. Not deposited: "
                f"{', '.join(failed_sites)}.[/yellow]"
            )

        # ================================================================
        # F. Rename residues in the prepared PDB
        # ================================================================
        prepared_pdb = None
        if self.workspace:
            prepared_pdb = self.workspace.get("prepared_pdb")
        if not prepared_pdb:
            prepared_pdb = str(self._final_pdb) if self._final_pdb else None

        if prepared_pdb and Path(prepared_pdb).exists():
            # Build PDB rename map: (chain, resid, old_resname) -> new_name
            # We need chain info — get it from the fingerprint residues and
            # the redox site atoms
            pdb_rename_map = {}
            redox_sites = self.redox_sites
            if not redox_sites and self.workspace:
                redox_sites = self.workspace.get("detected_redox_sites", [])

            for (resid, resname), new_name in residue_name_map.items():
                # Try to find chain from redox site atoms
                chain = self._find_chain_for_residue(resid, resname, redox_sites)
                pdb_rename_map[(chain, resid, resname)] = new_name

            # Write renamed PDB (overwrite in place)
            self.console.print(f"\n[bold]Renaming residues in PDB...[/bold]")
            rename_pdb_residues(prepared_pdb, pdb_rename_map, prepared_pdb)
            self.console.print(f"  [green]✓[/green] Updated: [grey50]{prepared_pdb}[/grey50]")

            renamed_names = ", ".join(
                f"{resname}{resid}→{new_name}"
                for (resid, resname), new_name in residue_name_map.items()
            )
            self.console.print(f"  [grey50]Renames: {renamed_names}[/grey50]")

            # ------------------------------------------------------------
            # Auto-emit a reusable rename transformer.
            # ------------------------------------------------------------
            # The renames we just applied ARE the recipe for reusing these
            # parameters on another instance of this site. Serialize them into
            # a data-only transformer (a rename table + per-residue metal
            # coordination signatures read off the site's bond graph) so the
            # Redox Site Preparer can re-apply them later without re-running
            # MCPB. All matching logic lives in AutoRenameTransformerBase; this
            # only produces data. Best-effort: the parameters are already
            # saved, so a failure here must not abort the integration. Computed
            # BEFORE the RedoxSite rename mutation below so signatures see the
            # original residue names/topology.
            try:
                from proprep.redoxsite_prep.transformation.auto_rename import (
                    emit_rename_transformer,
                    connectivity_signature,
                )

                def _site_containing(chain, resid):
                    for s in (redox_sites or []):
                        if any(getattr(a, "chain", None) == chain
                               and getattr(a, "resid", None) == resid
                               for a in getattr(s, "atoms", [])):
                            return s
                    return None

                # WL connectivity labels per site (cached), so each renamed
                # residue carries the fingerprint of its coordination environment
                # — this is what lets the reused transformer tell same-name
                # residues apart (e.g. the two Mn by their differing ligand sets).
                _label_cache = {}

                def _labels_for(host):
                    sid = id(host)
                    if sid not in _label_cache:
                        _label_cache[sid] = connectivity_signature(host)
                    return _label_cache[sid]

                # Per-resname antechamber atom-name maps ({pdb: mol2}). A small-
                # molecule ligand's deposited library carries antechamber's names,
                # but a reuse structure has the original PDB names; baking these
                # into the transformer lets it rename the reuse structure's atoms
                # to match the library. Source order: the live workspace, then the
                # on-disk mapping the parameterizer always writes (the workspace
                # copy is absent when the organic step was cached/skipped this run).
                _atom_maps_by_resname = self._collect_ligand_atom_maps(
                    {old for (_, _, old) in pdb_rename_map})

                # One transformer per deposited site, built from that site's
                # renames only. A table spanning several sites cannot match
                # anything: evaluate_redox_site (auto_rename.py) counts the
                # table's resnames against a SINGLE redox site and requires
                # met == total, so a merged table fails on every site it was
                # built from — each supplies only its own share of the residues.
                emitted = 0
                for record in deposits:
                    site_keys = set(record["residue_name_map"])
                    site_renames = {
                        (chain, resid, old): new
                        for (chain, resid, old), new in pdb_rename_map.items()
                        if (resid, old) in site_keys
                    }
                    if not site_renames:
                        continue

                    try:
                        rename_table = []
                        for (chain, resid, old_resname), new_name in site_renames.items():
                            entry = {"resname": old_resname, "target": new_name}
                            host = _site_containing(chain, resid)
                            if host is not None:
                                wl = _labels_for(host).get((chain, resid))
                                if wl is not None:
                                    entry["signature"] = wl
                            atom_renames = _atom_maps_by_resname.get(old_resname)
                            if atom_renames:
                                entry["atom_renames"] = atom_renames
                            rename_table.append(entry)

                        site_type = record["site_type"]
                        redox_state = record["redox_state"]
                        spin_state = record["spin_state"]
                        lib_result = record["lib_result"]

                        # Derive the relative cofactor path (under
                        # specialized_residues) so the Topology Generator can
                        # discover the deposited FF the same way built-in
                        # transformers do. library_path is the absolute
                        # ~/.proprep/.../specialized_residues/<cofactor_path> dir.
                        forcefield_path = None
                        lib_path = lib_result.get("library_path")
                        if lib_path:
                            parts = Path(lib_path).parts
                            if "specialized_residues" in parts:
                                idx = parts.index("specialized_residues")
                                forcefield_path = "/".join(parts[idx + 1:]) or None

                        resnames = sorted({old for (_, _, old) in site_renames})
                        tname = (f"mcpb_{site_type}_{'_'.join(resnames)}"
                                 f"_{redox_state}_{spin_state}")
                        tdesc = (
                            f"Reuse MCPB {site_type} parameters "
                            f"(residues: {', '.join(resnames)}; "
                            f"redox={redox_state}, spin={spin_state})"
                        )
                        t_path = emit_rename_transformer(
                            rename_table,
                            name=tname,
                            description=tdesc,
                            redox_state=str(redox_state),
                            spin_state=str(spin_state),
                            forcefield_path=forcefield_path,
                            provenance={
                                "source": "mcpb_metal_site",
                                "site_type": site_type,
                                "site_label": record["label"],
                                "redox_state": str(redox_state),
                                "spin_state": str(spin_state),
                                "library_path": lib_path,
                                "metadata_path": lib_result.get("metadata_path"),
                            },
                            site_types=[site_type],
                        )
                        record["transformer_path"] = str(t_path)
                        emitted += 1
                        self.console.print(
                            f"  [green]✓[/green] Reuse transformer for "
                            f"{record['label']}: [grey50]{t_path}[/grey50]"
                        )
                    except Exception as exc:
                        # The parameters are already deposited, so one site's
                        # transformer failing must not cost the others theirs.
                        self.console.print(
                            f"  [yellow]Note: could not auto-create reuse "
                            f"transformer for {record['label']} ({exc}); "
                            f"parameters are still saved.[/yellow]"
                        )

                if emitted:
                    self.console.print(
                        f"  [grey50]Apply {'them' if emitted > 1 else 'it'} in the "
                        f"Redox Site Preparer to reuse these params on another "
                        f"instance of {'these sites' if emitted > 1 else 'this site'}."
                        f"[/grey50]"
                    )
            except Exception as exc:
                self.console.print(
                    f"  [yellow]Note: could not auto-create reuse transformers "
                    f"({exc}); parameters are still saved.[/yellow]"
                )

            # Update RedoxSite objects to match the renamed PDB.
            # Without this, coord_to_pdb and atoms still carry pre-MCPB
            # names (e.g. HID) while the PDB now has MCPB names (e.g. HD1).
            if redox_sites:
                resid_to_new_name = {resid: new_name for (resid, _), new_name in residue_name_map.items()}
                for site in redox_sites:
                    # Update atoms
                    if hasattr(site, 'atoms'):
                        for atom in site.atoms:
                            if atom.resid in resid_to_new_name:
                                atom.resname = resid_to_new_name[atom.resid]
                    # Update coord_to_pdb
                    if hasattr(site, 'coord_to_pdb'):
                        for coords, pdb_info in site.coord_to_pdb.items():
                            if pdb_info.get('resid') in resid_to_new_name:
                                pdb_info['resname'] = resid_to_new_name[pdb_info['resid']]
                    # Update centers
                    if hasattr(site, 'centers'):
                        for center in site.centers:
                            if center.resid in resid_to_new_name:
                                center.resname = resid_to_new_name[center.resid]
                    # Update residue_groups keys
                    if hasattr(site, 'residue_groups'):
                        new_groups = {}
                        for key, coords_list in site.residue_groups.items():
                            chain, resid, *rest = key
                            if resid in resid_to_new_name:
                                new_key = (chain, resid, *rest) if len(rest) > 0 else (chain, resid)
                                # Replace old resname in key if present
                                new_groups[key] = coords_list
                            else:
                                new_groups[key] = coords_list
                        site.residue_groups = new_groups
                    # Update bond residue_info dicts
                    if hasattr(site, 'bonds'):
                        for bond in site.bonds:
                            if isinstance(bond.atom1_residue_info, dict):
                                if bond.atom1_residue_info.get('resid') in resid_to_new_name:
                                    bond.atom1_residue_info['resname'] = resid_to_new_name[bond.atom1_residue_info['resid']]
                            if isinstance(bond.atom2_residue_info, dict):
                                if bond.atom2_residue_info.get('resid') in resid_to_new_name:
                                    bond.atom2_residue_info['resname'] = resid_to_new_name[bond.atom2_residue_info['resid']]

                # Save updated sites back to workspace
                if self.workspace:
                    self.workspace.set("redox_sites", redox_sites)
                    self.workspace.set("detected_redox_sites", redox_sites)
                self.console.print(f"  [green]✓[/green] Updated RedoxSite residue names to match PDB")
        else:
            self.console.print("[yellow]No prepared PDB found — skipping PDB rename[/yellow]")

        # ================================================================
        # G. Store workspace data for tLEaP
        # ================================================================
        # The library entries are per-site, but the tLEaP run is one session
        # over the whole structure, so every site's files and types are
        # registered together here. Deduped: two sites sharing a bridging
        # residue deposit the same .lib, and a repeated loadoff/addAtomTypes
        # is at best noise in the generated input.
        all_lib_files = []
        all_frcmod_paths = []
        all_atom_types = []
        for record in deposits:
            result = record["lib_result"]
            all_lib_files.extend(result["renamed_mol2_files"])
            all_frcmod_paths.extend(result["frcmod_files"])
            all_atom_types.extend(result.get("atom_type_entries", []))

        def _dedupe(seq):
            seen = set()
            out = []
            for item in seq:
                if item not in seen:
                    seen.add(item)
                    out.append(item)
            return out

        all_lib_files = _dedupe(all_lib_files)
        all_frcmod_paths = _dedupe(all_frcmod_paths)
        all_atom_types = _dedupe(all_atom_types)

        if self.workspace:
            # Append renamed mol2 paths to preprocessing_lib_files
            # (tLEaP generator handles .mol2 -> loadmol2 automatically)
            lib_files = self.workspace.get("preprocessing_lib_files", [])
            if not isinstance(lib_files, list):
                lib_files = []
            lib_files.extend(all_lib_files)
            self.workspace.set("preprocessing_lib_files", lib_files)

            # Append frcmod paths to preprocessing_frcmod_files
            frcmod_files = self.workspace.get("preprocessing_frcmod_files", [])
            if not isinstance(frcmod_files, list):
                frcmod_files = []
            frcmod_files.extend(all_frcmod_paths)
            self.workspace.set("preprocessing_frcmod_files", frcmod_files)

            # Store custom atom types for tLEaP addAtomTypes block
            if all_atom_types:
                existing = self.workspace.get("preprocessing_atom_types", [])
                if not isinstance(existing, list):
                    existing = []
                existing.extend(all_atom_types)
                self.workspace.set("preprocessing_atom_types", existing)

        # ================================================================
        # H. Show advice panel
        # ================================================================
        deposit_table = Table(title="Deposited Library Entries", expand=False)
        deposit_table.add_column("Site", style="cyan")
        deposit_table.add_column("Type / Redox / Spin", style="green bold")
        deposit_table.add_column("Path", style="grey50")
        for record in deposits:
            deposit_table.add_row(
                record["label"],
                f"{record['site_type']} / {record['redox_state']} / {record['spin_state']}",
                str(record["lib_result"]["library_path"]),
            )
        self.console.print()
        self.console.print(deposit_table)

        n_transformers = sum(1 for r in deposits if r.get("transformer_path"))
        plural = "s" if len(deposits) > 1 else ""
        self.console.print(Panel(
            f"[bold green]Force Field Integration Complete[/bold green]\n\n"
            f"[bold]Library entries:[/bold] {len(deposits)} "
            f"(one per metal site)\n"
            f"[bold]Renamed mol2 files:[/bold] {len(all_lib_files)}\n"
            f"[bold]Frcmod files:[/bold] {len(all_frcmod_paths)}\n\n"
            f"[bold]Next steps:[/bold]\n"
            f"  1. Go to the [cyan]Topology Generator[/cyan] to build prmtop/inpcrd\n"
            f"     (mol2 and frcmod files are already registered for tLEaP)\n"
            f"  2. {n_transformers} [cyan]reuse transformer{plural}[/cyan] "
            f"auto-created from these renames —\n"
            f"     apply in the [cyan]Redox Site Preparer[/cyan] to reuse a site's\n"
            f"     parameters on another instance of it (no re-parameterization)",
            title="Force Field Integration",
            border_style="green",
            expand=False,
        ))

        site_names = ", ".join(f"'{r['site_type']}'" for r in deposits)
        summary = (f"Integrated {len(all_lib_files)} mol2, "
                   f"{len(all_frcmod_paths)} frcmod into {len(deposits)} library "
                   f"entr{'ies' if len(deposits) > 1 else 'y'}: {site_names}")
        if failed_sites:
            summary += f" (not deposited: {', '.join(failed_sites)})"
        return {"summary": summary}

    def _find_chain_for_residue(self, resid: int, resname: str,
                                redox_sites) -> str:
        """Find the chain ID for a residue from RedoxSite atoms.

        Falls back to 'A' if the residue cannot be found.
        """
        if not redox_sites:
            return "A"

        for site in redox_sites:
            atoms = []
            if hasattr(site, 'atoms'):
                atoms = site.atoms
            elif isinstance(site, dict):
                atoms = site.get('atoms', [])

            for atom in atoms:
                atom_resid = getattr(atom, 'resid', None)
                atom_resname = getattr(atom, 'resname', None)
                atom_chain = getattr(atom, 'chain', None)
                if isinstance(atom, dict):
                    atom_resid = atom.get('resid')
                    atom_resname = atom.get('resname')
                    atom_chain = atom.get('chain')

                if atom_resid == resid and atom_resname == resname and atom_chain:
                    return atom_chain

        return "A"

    # =========================================================================
    # Component Processing Methods (called by checklist handlers)
    # =========================================================================

    def _run_triage_only(self, pdb_file: str) -> Dict[str, str]:
        """
        Run structure triage (categorization only, no FF selection).

        Categorizes residues into:
        - A: Standard protein residues
        - B: Organic small molecules (no metal)
        - C: Water molecules
        - D: Isolated metal ions (single atom)
        - E: Organometallic small molecules (contains embedded metal)

        Returns:
            Dict mapping residue_key (chain:resid:resname) to category
        """
        from Bio.PDB import PDBParser

        self.console.print(Panel(
            "[bold]Structure Triage[/bold]\n"
            "Categorizing residues by type",
            border_style="cyan",
            expand=False
        ))

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("temp", pdb_file)

        triage = {}
        cat_counts = {'A': 0, 'B': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0}

        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    resid = residue.id[1]
                    res_key = f"{chain.id}:{resid}:{resname}"

                    # Categorize
                    # A=protein, B=organic, C=organometallic, D=isolated metal,
                    # E=water, F=pure inorganic metal cluster (Fe-S etc.)
                    if resname in STANDARD_RESIDUES:
                        triage[res_key] = 'A'
                        cat_counts['A'] += 1
                    elif resname in WATER_RESIDUES:
                        triage[res_key] = 'E'
                        cat_counts['E'] += 1
                    elif self._is_isolated_metal(residue):
                        triage[res_key] = 'D'
                        cat_counts['D'] += 1
                    elif self._is_pure_metal_cluster(residue):
                        # Multi-atom, metal-bearing, NO carbon → no organic
                        # fragment to hand to the small-molecule parameterizer.
                        # The whole cluster is owned by MCPB; check this BEFORE
                        # the organometallic branch (which assumes an organic
                        # scaffold like heme/MoCo pterin).
                        triage[res_key] = 'F'
                        cat_counts['F'] += 1
                    elif self._has_embedded_metal(residue):
                        triage[res_key] = 'C'
                        cat_counts['C'] += 1
                    else:
                        triage[res_key] = 'B'
                        cat_counts['B'] += 1

        # Collect residues by category for the identity column
        cat_residues = {'A': [], 'B': [], 'C': [], 'D': [], 'E': [], 'F': []}
        for res_key, cat in triage.items():
            if cat in cat_residues:
                cat_residues[cat].append(res_key)

        def format_residue_list(res_keys):
            """Format residue list for display."""
            if not res_keys:
                return ""
            return ", ".join(res_keys)

        # Display triage table with identity column
        # Order: A (protein), B (organic), C (organometallic),
        #        F (metal cluster), D (isolated metal), E (water)
        table = Table(title="Structure Triage Results", show_lines=True)
        table.add_column("Category", style="cyan", width=10)
        table.add_column("Count", style="white", width=8)
        table.add_column("Description", style="grey50", width=28)
        table.add_column("Residues", style="green")

        table.add_row("A", str(cat_counts['A']), "Standard protein residues",
                      f"{cat_counts['A']} residues" if cat_counts['A'] > 0 else "")
        table.add_row("B", str(cat_counts['B']), "Organic small molecules",
                      format_residue_list(cat_residues['B']))
        table.add_row("C", str(cat_counts['C']), "Organometallic small molecules",
                      format_residue_list(cat_residues['C']))
        table.add_row("F", str(cat_counts['F']), "Metal clusters (Fe-S, pure inorganic)",
                      format_residue_list(cat_residues['F']))
        table.add_row("D", str(cat_counts['D']), "Isolated metal ions",
                      format_residue_list(cat_residues['D']))
        table.add_row("E", str(cat_counts['E']), "Water molecules",
                      f"{cat_counts['E']} molecules" if cat_counts['E'] > 0 else "")

        self.console.print(table)

        # Store in workspace
        if self.workspace:
            self.workspace.set("preprocessing_triage", triage)

        return triage

    def _is_isolated_metal(self, residue) -> bool:
        """Check if residue is an isolated metal ion (single atom).

        Checks element field first, then falls back to residue name and atom name
        since PDB files often have empty element fields for metal ions.
        """
        atoms = list(residue.get_atoms())
        if len(atoms) == 1:
            # Check element field (normalize to title case for METAL_ELEMENTS)
            element = atoms[0].element.strip().title() if atoms[0].element else ""
            if element:
                # Authoritative when present — the name-based guesses below
                # exist for the blank-field case and must not override it.
                return element in METAL_ELEMENTS
            # Fallback: check residue name (many metal ions have resname = element symbol)
            resname = residue.resname.strip().title()
            if resname in METAL_ELEMENTS:
                return True
            # Fallback: check atom name, honouring the PDB column convention
            if self._element_from_atom_name(atoms[0]) in METAL_ELEMENTS:
                return True
        return False

    @staticmethod
    def _element_from_atom_name(atom) -> str:
        """Element implied by an atom NAME, honouring the PDB column convention.

        A two-letter element starts in column 13 (``FE1 `` is iron); a
        one-letter element is right-justified into column 14, leaving room for
        a remoteness indicator (`` PA `` is phosphorus, `` CA `` an alpha
        carbon). Reading the stripped name instead turns FAD's phosphates PA
        and PB into protactinium and lead, and a protein's CA into calcium.

        Only for atoms whose element field is missing; BioPython keeps the raw
        four-character field in ``fullname``.
        """
        from proprep.utils.pdb_format import element_from_name_field

        full = getattr(atom, "fullname", None) or ""
        if len(full) >= 4:
            return element_from_name_field(full)
        return atom.name.strip().title()

    def _has_embedded_metal(self, residue) -> bool:
        """Check if residue contains an embedded metal (multi-atom with metal).

        For multi-atom residues, checks if any atom is a metal element.
        Falls back to atom name only when the element field is empty.
        """
        atoms = list(residue.get_atoms())
        if len(atoms) <= 1:
            return False  # Single atom residues handled by _is_isolated_metal
        for atom in atoms:
            element = atom.element.strip().title() if atom.element else ""
            if element:
                # The element field is authoritative when present. Consulting
                # the atom name anyway classified FAD as organometallic: its
                # two phosphate atoms are named PA and PB, which title-case to
                # Pa (protactinium) and Pb (lead).
                if element in METAL_ELEMENTS:
                    return True
                continue
            # No element field (common in older/hand-edited PDBs, and the
            # reason this fallback exists — e.g. "FE" in a heme).
            if self._element_from_atom_name(atom) in METAL_ELEMENTS:
                return True
        return False

    def _element_of(self, atom, resname: str = "") -> str:
        """Best-effort element symbol (Title case) for an atom.

        BioPython often leaves the element field blank or mis-parses metals, so
        fall back to the residue name (metal ions frequently share it) and then
        the atom name.
        """
        element = atom.element.strip().title() if atom.element else ""
        if element:
            return element
        if resname:
            guess = resname.strip().title()
            if guess in METAL_ELEMENTS:
                return guess
        return atom.name.strip().title()

    def _is_pure_metal_cluster(self, residue) -> bool:
        """Check if a residue is a pure inorganic metal cluster (Fe-S, etc.).

        A pure cluster is a multi-atom residue that contains at least one metal
        and NO carbon: there is no organic scaffold, so nothing can be split off
        for the small-molecule parameterizer. The whole cluster (metals plus
        bridging atoms such as the sulfides of an Fe2S2/Fe4S4) is owned by MCPB.

        This is the discriminator against an organometallic cofactor (category
        C) like heme or the MoCo pterin, whose carbon-bearing organic part IS
        parameterized separately with the metal removed. The presence of carbon
        is the marker of that organic scaffold.
        """
        atoms = list(residue.get_atoms())
        if len(atoms) <= 1:
            return False  # single-atom metals are category D

        resname = residue.resname
        has_metal = False
        for atom in atoms:
            element = self._element_of(atom, resname)
            if element in METAL_ELEMENTS:
                has_metal = True
            elif element == 'C':
                # An organic scaffold → organometallic (category C), not a
                # pure cluster. (Metals like Ca/Cd/Co/Cu are caught above, so
                # a bare 'C' here is carbon.)
                return False
        return has_metal

    def _get_metal_atoms_in_residue(self, residue) -> List[Tuple[str, str]]:
        """Get list of (atom_name, element) for metal atoms in residue."""
        metals = []
        for atom in residue.get_atoms():
            element = atom.element.strip().title() if atom.element else ""
            if element in METAL_ELEMENTS:
                metals.append((atom.name.strip(), element))
        return metals

    def _get_metal_atoms_in_residue_by_key(self, res_key: str, structure) -> List[Tuple[str, str]]:
        """Get metal atoms in residue by res_key."""
        parts = res_key.split(':')
        chain_id, resid = parts[0], int(parts[1])
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        if residue.id[1] == resid:
                            return self._get_metal_atoms_in_residue(residue)
        return []

    # =========================================================================
    # New Component Processing Methods (for restructured workflow)
    # =========================================================================

    def _process_organic_residue(self, res_key: str) -> Optional[dict]:
        """
        Process an organic small molecule (category B).

        Prompts user for: have params or need to generate.

        Returns:
            Dict with lib/frcmod info, or None if failed
        """
        parts = res_key.split(':')
        resname = parts[2] if len(parts) > 2 else "UNK"

        self.console.print(f"\n[bold cyan]Organic residue: {resname} ({res_key})[/bold cyan]")

        self.console.print("  [1] I have parameters (lib/frcmod files)")
        self.console.print("  [2] I need to generate parameters")

        choice = prompt_with_context(
            self.processor,
            f"Select option for {resname}",
            choices=["1", "2"],
            module="Structure Preprocessor",
            description=f"Parameter source for {resname}"
        )

        if choice == "1":
            # Get lib/frcmod files
            lib_path, frcmod_path, ff_resname = self._prompt_custom_ff_files(resname)
            if lib_path:
                # Handle atom name mapping (PDB names may differ from lib names)
                atom_mapping = self._handle_atom_name_mapping(res_key, lib_path, resname)

                return {
                    'res_key': res_key,
                    'lib_file': str(lib_path),
                    'frcmod_file': str(frcmod_path) if frcmod_path else None,
                    'ff_resname': ff_resname or resname,
                    'atom_name_mapping': atom_mapping,
                }
            return None
        else:
            # Launch Small Molecule Parameterizer
            result = self._launch_small_molecule_parameterizer(resname, res_key, self.triage_results)
            if result:
                # Prefer the generated .lib (loadoff registers a template that
                # matches the structure residue by entry name). Only fall back to
                # the mol2 if no lib was produced.
                return {
                    'res_key': res_key,
                    'source': 'generated',
                    'mol2_file': result.get('mol2_file'),
                    'lib_file': result.get('lib_file') or result.get('mol2_file'),
                    'frcmod_file': result.get('frcmod_file'),
                    'atom_name_mapping': result.get('atom_name_mapping'),
                }
            return None

    def _process_organometallic_residue(self, res_key: str) -> Optional[dict]:
        """
        Process an organometallic small molecule (category C).

        Prompts user for parameter source and handles metal removal tracking.

        Returns:
            Dict with lib/frcmod info and metal tracking, or None if failed
        """
        from Bio.PDB import PDBParser

        parts = res_key.split(':')
        chain_id, resid = parts[0], int(parts[1])
        resname = parts[2] if len(parts) > 2 else "UNK"

        # Get metal info
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("temp", self._pdb_file)
        metals = self._get_metal_atoms_in_residue_by_key(res_key, structure)
        metal_str = ", ".join([f"{m[0]}({m[1]})" for m in metals])

        self.console.print(f"\n[bold yellow]Organometallic residue: {resname} ({res_key})[/bold yellow]")
        self.console.print(f"  Contains metal(s): {metal_str}")

        self.console.print("\n  [1] I have parameters INCLUDING the metal")
        self.console.print("  [2] I have parameters for ORGANIC PART only")
        self.console.print("  [3] I need to generate parameters for the organic part")

        choice = prompt_with_context(
            self.processor,
            f"Select option for {resname}",
            choices=["1", "2", "3"],
            module="Structure Preprocessor",
            description=f"Parameter source for {resname}"
        )

        result = {
            'res_key': res_key,
            'metal_removal_needed': False,
            'metal_info': None,
        }

        if choice == "1":
            # Full params including metal - no removal needed
            lib_path, frcmod_path, ff_resname = self._prompt_custom_ff_files(resname)
            if lib_path:
                # Handle atom name mapping (PDB names may differ from lib names)
                atom_mapping = self._handle_atom_name_mapping(res_key, lib_path, resname)

                result['lib_file'] = str(lib_path)
                result['frcmod_file'] = str(frcmod_path) if frcmod_path else None
                result['ff_resname'] = ff_resname or resname
                result['metal_removal_needed'] = False
                result['atom_name_mapping'] = atom_mapping
                return result
            return None

        elif choice == "2":
            # Organic-only params - need to remove metal
            if confirm_with_context(
                self.processor,
                "Do you need to edit the residue first (atom exclusion + H capping)?",
                default=False,
                module="Structure Preprocessor",
                description="Edit residue (atom exclusion + H capping) before parameterization",
            ):
                # TODO: Implement per-residue atom exclusion
                self.console.print("[yellow]Per-residue editing not yet implemented[/yellow]")

            lib_path, frcmod_path, ff_resname = self._prompt_custom_ff_files(resname)
            if lib_path:
                # Handle atom name mapping (PDB names may differ from lib names)
                atom_mapping = self._handle_atom_name_mapping(res_key, lib_path, resname)

                # Extract metal info for later reinsertion
                metal_info = self._extract_metal_info_from_key(res_key, is_isolated=False)
                result['lib_file'] = str(lib_path)
                result['frcmod_file'] = str(frcmod_path) if frcmod_path else None
                result['ff_resname'] = ff_resname or resname
                result['metal_removal_needed'] = True
                result['metal_info'] = metal_info.to_dict() if metal_info else None
                result['atom_name_mapping'] = atom_mapping
                return result
            return None

        else:
            # Generate params for organic part
            if confirm_with_context(
                self.processor,
                "Do you need to edit the residue first (atom exclusion + H capping)?",
                default=False,
                module="Structure Preprocessor",
                description="Edit residue (atom exclusion + H capping) before parameterization",
            ):
                # TODO: Implement per-residue atom exclusion
                self.console.print("[yellow]Per-residue editing not yet implemented[/yellow]")

            # Extract metal info for later reinsertion
            metal_info = self._extract_metal_info_from_key(res_key, is_isolated=False)

            # Launch parameterizer on organic part (metal will be removed)
            param_result = self._launch_small_molecule_parameterizer(resname, res_key, self.triage_results)
            if param_result:
                result['source'] = 'generated'
                result['mol2_file'] = param_result.get('mol2_file')
                result['lib_file'] = param_result.get('lib_file') or param_result.get('mol2_file')
                result['frcmod_file'] = param_result.get('frcmod_file')
                result['atom_name_mapping'] = param_result.get('atom_name_mapping')
                result['metal_removal_needed'] = True
                result['metal_info'] = metal_info.to_dict() if metal_info else None
                return result
            return None

    def _extract_metal_info_from_key(self, res_key: str, is_isolated: bool) -> Optional[MetalInfo]:
        """Extract MetalInfo for a residue by its key."""
        from Bio.PDB import PDBParser

        parts = res_key.split(':')
        chain_id, resid = parts[0], int(parts[1])
        resname = parts[2] if len(parts) > 2 else "UNK"

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("temp", self._pdb_file)

        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        if residue.id[1] == resid:
                            # Find metal atom(s)
                            for atom in residue.get_atoms():
                                element = atom.element.strip().title() if atom.element else ""
                                # BioPython can mis-parse element for metals (e.g., ZN → N)
                                # Fall back to residue name or atom name
                                if element not in METAL_ELEMENTS:
                                    element = resname.strip().title()
                                if element not in METAL_ELEMENTS:
                                    element = atom.name.strip().title()
                                if element in METAL_ELEMENTS:
                                    return MetalInfo(
                                        atom_name=atom.name.strip(),
                                        element=element,
                                        coords=tuple(atom.coord),
                                        original_chain=chain_id,
                                        original_resid=resid,
                                        original_resname=resname,
                                        is_isolated=is_isolated,
                                    )
        return None

    def _offer_cluster_hydrogens(self, cluster_keys: List[str]) -> None:
        """Offer to add hydrogens to each withheld inorganic cluster.

        Nothing else in the pipeline can put a hydrogen on a cluster: hydrogen
        addition covers protein (category A) and organic residues (category B),
        and ``reduce`` has no chemistry for a Mo-S-O or Fe-S core anyway. So a
        cofactor whose resting state carries a hydroxo — Mo(=O)(=S)(OH) in a
        molybdenum cofactor — reaches the QM model as a bare oxo, with the wrong
        electron count and the wrong charge.

        Editing the generated .gjf by hand is not an alternative. The Gaussian
        input and the model PDB are matched by index: the PDB generates the
        fingerprint, the RESP input and the final mol2, while the Gaussian
        output supplies the Hessian and the ESP. Adding an atom to one of them
        shifts Seminario's atom indices (silently — it validates the Hessian
        only against its own coordinates) and leaves the deposited residue
        template without the hydrogen its charges were fitted with.

        Offered for every cluster rather than guessing which ones want a
        hydrogen; an Fe-S cluster simply declines. The default is no.
        """
        if not self._interactive or not cluster_keys:
            return

        from proprep.forcefield_prep.hydrogen_editor import HydrogenEditor

        for res_key in cluster_keys:
            parts = res_key.split(':')
            if len(parts) < 3:
                continue
            chain_id, resid, resname = parts[0], int(parts[1]), parts[2]

            existing = self._extract_cluster_atoms_from_key(res_key)
            if not existing:
                continue
            composition = ", ".join(sorted({a.element for a in existing}))
            names = ", ".join(a.atom_name for a in existing)

            self.console.print(
                f"\n[bold]{resname} {chain_id}:{resid}[/bold] "
                f"[grey50]({len(existing)} atoms — {names}; elements {composition})[/grey50]")
            self.console.print(
                "[grey50]  A cluster is not hydrogenated anywhere else in the "
                "pipeline. Add one here if this cofactor's resting state carries "
                "a hydroxo/protonated ligand (e.g. Mo-OH); an Fe-S cluster does "
                "not need one.[/grey50]")

            # Put the cluster on screen while it is the subject of the prompt —
            # "which oxygen" is a question about geometry.
            selection = f":{chain_id} and {resid}"
            self._focus_viewer_on_cluster(selection)

            if not confirm_with_context(
                self.processor,
                f"Add hydrogen(s) to {resname} {chain_id}:{resid}?",
                default=False,
                module="Structure Preprocessor",
                description="Add hydrogens to an inorganic cluster",
            ):
                continue

            work_dir = Path(self._output_dir) if self._output_dir else Path.cwd()
            work_dir.mkdir(parents=True, exist_ok=True)
            residue_pdb = work_dir / f"cluster_{resname}_{chain_id}_{resid}.pdb"
            try:
                self._extract_single_residue_to_pdb(
                    self._pdb_file, chain_id, resid, resname, residue_pdb)
            except Exception as exc:  # noqa: BLE001
                self.console.print(
                    f"[yellow]Could not extract {res_key} for editing ({exc}); "
                    f"skipping.[/yellow]")
                continue

            editor = HydrogenEditor(
                str(residue_pdb), f"{resname}_{chain_id}_{resid}",
                console=self.console, processor=self.processor,
                interactive=True, residue_name=resname,
                module="Structure Preprocessor",
            )
            try:
                if editor.add_interactive():
                    added = self._merge_cluster_hydrogens(
                        str(residue_pdb), chain_id, resid, resname)
                    self.console.print(
                        f"  [green]✓ {added} hydrogen(s) merged into the "
                        f"structure[/green]")
                    # The file changed on disk; re-serve it so the viewer shows
                    # the hydrogen that was just added rather than the state
                    # before it.
                    self._focus_viewer_on_cluster(selection, refresh=True)
                    self.console.print(
                        "[grey50]  Their formal charge is asked for with the rest "
                        "of the cluster's core atoms during MCPB atom typing (a "
                        "hydroxo O is -1 where an oxo O is -2).[/grey50]")
            except Exception as exc:  # noqa: BLE001
                self.console.print(
                    f"[yellow]Hydrogen editing failed for {res_key} ({exc}); "
                    f"the cluster is unchanged.[/yellow]")

    def _focus_viewer_on_cluster(self, selection: str, *, refresh: bool = False) -> None:
        """Show the cluster residue in the viewer, optionally re-reading the file.

        ``show_structure`` points the viewer at the structure being edited (a
        no-op if it is already the one displayed). ``refresh`` is what makes an
        added hydrogen visible: the file was rewritten in place, and re-issuing
        the same path alone would not re-read it.

        Best-effort throughout — the viewer is an aid, and a headless or closed
        one must not interrupt the prompt flow.
        """
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer

            if self._pdb_file:
                _viewer.show_structure(str(self._pdb_file))
            if refresh:
                _viewer.refresh_structure()

            _viewer.unhighlight("cluster_h_focus")
            _viewer.highlight(selection, style="ball+stick", color="element",
                              label="cluster_h_focus")
            _viewer.focus_on(selection)
        except Exception as exc:  # noqa: BLE001
            logger.debug("Viewer focus for %s unavailable: %s", selection, exc)

    def _merge_cluster_hydrogens(self, residue_pdb: str, chain_id: str,
                                 resid: int, resname: str) -> int:
        """Copy hydrogens from an edited residue PDB back into the structure.

        Inserted directly after that residue's existing atoms so the residue
        stays contiguous, then every serial in the file is renumbered. Returns
        the number of hydrogens merged.
        """
        from proprep.utils.pdb_format import atom_name_field

        def _is_atom(line):
            return line.startswith(("ATOM", "HETATM"))

        with open(residue_pdb) as fh:
            edited = [ln for ln in fh if _is_atom(ln)]
        new_h = [ln for ln in edited
                 if (ln[76:78].strip().upper() == "H"
                     or (not ln[76:78].strip() and ln[12:16].strip().startswith("H")))]
        if not new_h:
            return 0

        with open(self._pdb_file) as fh:
            lines = fh.readlines()

        # Rewrite each H onto the target residue's identity, so it belongs to
        # the cluster rather than to whatever the extracted file called it.
        rebuilt = []
        for ln in new_h:
            name = ln[12:16].strip()
            rebuilt.append(
                f"HETATM{0:5d} {atom_name_field(name, 'H')} "
                f"{resname:>3.3s} {chain_id}{resid:>4d}    "
                f"{ln[30:54]}"
                f"  1.00  0.00          {'H':>2s}\n"
            )

        last_idx = -1
        for i, ln in enumerate(lines):
            if _is_atom(ln) and ln[21] == chain_id:
                try:
                    if int(ln[22:26]) == resid:
                        last_idx = i
                except ValueError:
                    continue
        if last_idx < 0:
            return 0

        lines[last_idx + 1:last_idx + 1] = rebuilt

        serial = 0
        for i, ln in enumerate(lines):
            if _is_atom(ln):
                serial += 1
                lines[i] = f"{ln[:6]}{serial:5d}{ln[11:]}"

        with open(self._pdb_file, "w") as fh:
            fh.writelines(lines)
        return len(rebuilt)

    def _extract_cluster_atoms_from_key(self, res_key: str) -> List[MetalInfo]:
        """Extract every atom of a pure-cluster residue as MetalInfo.

        All atoms share cluster_id = res_key so reinsertion rebuilds them as one
        residue. Unlike _extract_metal_info_from_key (which returns only the
        first metal), this returns the metals AND the bridging atoms, because
        the whole cluster is withheld from the standard-FF tLEaP pass.
        """
        from Bio.PDB import PDBParser

        parts = res_key.split(':')
        chain_id, resid = parts[0], int(parts[1])
        resname = parts[2] if len(parts) > 2 else "UNK"

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("temp", self._pdb_file)

        cluster_atoms: List[MetalInfo] = []
        for model in structure:
            for chain in model:
                if chain.id != chain_id:
                    continue
                for residue in chain:
                    if residue.id[1] != resid:
                        continue
                    for atom in residue.get_atoms():
                        element = self._element_of(atom, resname)
                        cluster_atoms.append(MetalInfo(
                            atom_name=atom.name.strip(),
                            element=element,
                            coords=tuple(atom.coord),
                            original_chain=chain_id,
                            original_resid=resid,
                            original_resname=resname,
                            is_isolated=True,       # cluster is its own residue
                            cluster_id=res_key,     # groups the atoms on reinsertion
                        ))
        return cluster_atoms

    def _collect_metals_to_remove(self) -> List[MetalInfo]:
        """Collect all metals that need to be removed before tLEaP."""
        metals = []

        # Isolated metals (category D) - all need removal
        isolated = self.workspace.get("preprocessing_isolated_metals", {}) if self.workspace else {}
        for res_key, info in isolated.items():
            if 'metal_info' in info:
                metals.append(MetalInfo.from_dict(info['metal_info']))

        # Organometallic (category C) where metal_removal_needed
        orgmet = self.workspace.get("preprocessing_organometallic_ff", {}) if self.workspace else {}
        for resname, info in orgmet.items():
            if info.get('metal_removal_needed') and info.get('metal_info'):
                metals.append(MetalInfo.from_dict(info['metal_info']))

        # Pure metal clusters (category F) - every atom of the cluster residue
        clusters = self.workspace.get("preprocessing_metal_clusters", {}) if self.workspace else {}
        for res_key, info in clusters.items():
            for atom_dict in info.get('atoms', []):
                metals.append(MetalInfo.from_dict(atom_dict))

        return metals

    def _collect_ligand_atom_maps(self, resnames: Set[str]) -> Dict[str, Dict[str, str]]:
        """Return ``{resname: {pdb_atom: mol2_atom}}`` for antechamber-renamed
        ligands, non-identity entries only.

        The parameterizer renames a ligand's atoms (PDB -> antechamber) and
        records the map; a reuse structure needs the inverse baked into the
        transformer so its atoms match the antechamber-named library. Source
        order: the live workspace first, then the on-disk mapping the small-
        molecule parameterizer always writes
        (``small_molecule_params_<RES>/<res>_atom_name_mapping.json``) — the
        workspace copy is missing whenever the organic step was cached or skipped
        this run, which is exactly when re-emitting a transformer must still
        recover the map.
        """
        import os
        import glob as _glob
        import json as _json

        maps: Dict[str, Dict[str, str]] = {}

        def _keep_diff(m):
            return {pdb: mol2 for pdb, mol2 in (m or {}).items() if pdb != mol2}

        # 1) Workspace (organic + organometallic FF results).
        organic = self.workspace.get("preprocessing_organic_ff", {}) if self.workspace else {}
        orgmet = self.workspace.get("preprocessing_organometallic_ff", {}) if self.workspace else {}
        for resname, info in {**organic, **orgmet}.items():
            diff = _keep_diff(info.get("atom_name_mapping"))
            if resname in resnames and diff:
                maps[resname] = diff

        # 2) On-disk fallback for any still-missing resname.
        base_dirs = []
        for cand in (self.workspace.get("prepared_pdb") if self.workspace else None,
                     str(self._final_pdb) if getattr(self, "_final_pdb", None) else None):
            if cand:
                base_dirs.append(str(Path(cand).parent))
        base_dirs.append(os.getcwd())

        for resname in resnames:
            if resname in maps:
                continue
            for base in base_dirs:
                pattern = os.path.join(
                    base, f"small_molecule_params_{resname}",
                    f"{resname.lower()}_atom_name_mapping.json")
                hits = _glob.glob(pattern)
                if not hits:
                    continue
                try:
                    with open(hits[0]) as f:
                        diff = _keep_diff(_json.load(f))
                except Exception as e:
                    logger.debug("Could not read atom map %s: %s", hits[0], e)
                    continue
                if diff:
                    maps[resname] = diff
                    break

        return maps

    def _collect_ligand_frcmods(self, resnames: Set[str]) -> Dict[str, str]:
        """Return ``{resname: frcmod_path}`` for the site's organic ligands.

        A metal site's organic ligand (e.g. E4Z) is parameterized by the
        small-molecule parameterizer, which writes a parmchk2 GAFF frcmod
        (``small_molecule_params_<RES>/<res>.frcmod``) holding the vdW/torsion
        terms base GAFF lacks. MCPB's bonded frcmod does NOT contain these, so
        the ligand frcmod must be deposited alongside it and loaded at reuse —
        otherwise the ligand's atoms type but resolve no parameters. Sourced
        like the atom-name maps: live workspace first, then the on-disk file the
        parameterizer always writes (present even when the organic step cached).
        """
        import os
        import glob as _glob

        frcmods: Dict[str, str] = {}

        # 1) Workspace (organic + organometallic FF results carry frcmod_file).
        organic = self.workspace.get("preprocessing_organic_ff", {}) if self.workspace else {}
        orgmet = self.workspace.get("preprocessing_organometallic_ff", {}) if self.workspace else {}
        for resname, info in {**organic, **orgmet}.items():
            fm = info.get("frcmod_file")
            if resname in resnames and fm and Path(fm).is_file():
                frcmods[resname] = str(fm)

        # 2) On-disk fallback for any still-missing resname.
        base_dirs = []
        for cand in (self.workspace.get("prepared_pdb") if self.workspace else None,
                     str(self._final_pdb) if getattr(self, "_final_pdb", None) else None):
            if cand:
                base_dirs.append(str(Path(cand).parent))
        base_dirs.append(os.getcwd())

        for resname in resnames:
            if resname in frcmods:
                continue
            for base in base_dirs:
                pattern = os.path.join(
                    base, f"small_molecule_params_{resname}",
                    f"{resname.lower()}.frcmod")
                hits = _glob.glob(pattern)
                if hits:
                    frcmods[resname] = hits[0]
                    break

        return frcmods

    def _remove_metals_from_structure(self, metals: List[MetalInfo], pdb_file: str = None) -> Path:
        """Remove specified metal atoms from structure."""
        from Bio.PDB import PDBParser, PDBIO, Select

        class ExcludeMetals(Select):
            def __init__(self, metals_to_remove):
                self.metals = {
                    (m.original_chain, m.original_resid, m.atom_name)
                    for m in metals_to_remove
                }

            def accept_atom(self, atom):
                res = atom.get_parent()
                chain = res.get_parent()
                key = (chain.id, res.id[1], atom.name.strip())
                return key not in self.metals

        input_file = pdb_file if pdb_file else self._pdb_file
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("struct", input_file)

        output = self._output_dir / "metal_free.pdb"
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output), ExcludeMetals(metals))

        self.console.print(f"[grey50]Created metal-free structure: {output.name}[/grey50]")
        return output

    def _apply_atom_name_mappings(self, pdb_file: Path) -> Optional[Path]:
        """
        Apply atom name mappings to PDB file.

        Renames atoms in non-standard residues so they match the lib file atom names.
        Also fixes residue name case to match the lib file (tLEaP is case-sensitive).

        Args:
            pdb_file: Input PDB file

        Returns:
            Path to remapped PDB file, or None if no mappings to apply
        """
        from Bio.PDB import PDBParser, PDBIO

        # Collect all mappings and lib files from organic and organometallic FF info
        all_mappings = {}  # (chain, resid) -> {'atom_mapping': {...}, 'lib_file': path, 'lib_resname': str}

        organic_ff = self.workspace.get("preprocessing_organic_ff", {}) if self.workspace else {}
        orgmet_ff = self.workspace.get("preprocessing_organometallic_ff", {}) if self.workspace else {}

        for resname, info in {**organic_ff, **orgmet_ff}.items():
            res_key = info.get('res_key', '')
            lib_file = info.get('lib_file')
            if res_key and lib_file:
                parts = res_key.split(':')
                if len(parts) >= 2:
                    chain_id = parts[0]
                    resid = int(parts[1])
                    # Get residue name from lib file (for correct case)
                    lib_resname = self._get_lib_residue_name(lib_file)
                    all_mappings[(chain_id, resid)] = {
                        'atom_mapping': info.get('atom_name_mapping', {}),
                        'lib_file': lib_file,
                        'lib_resname': lib_resname,
                    }

        if not all_mappings:
            return None

        self.console.print(f"[grey50]Applying atom/residue name mappings for {len(all_mappings)} residue(s)[/grey50]")

        # Parse and modify structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("struct", str(pdb_file))

        atoms_renamed = 0
        residues_renamed = 0
        for model in structure:
            for chain in model:
                for residue in chain:
                    key = (chain.id, residue.id[1])
                    if key in all_mappings:
                        info = all_mappings[key]
                        atom_mapping = info['atom_mapping']
                        lib_resname = info['lib_resname']

                        # Fix residue name case to match lib file
                        if lib_resname and residue.resname != lib_resname:
                            old_resname = residue.resname
                            residue.resname = lib_resname
                            residues_renamed += 1
                            self.console.print(f"[grey50]  Residue {old_resname} → {lib_resname}[/grey50]")

                        # Apply atom name mapping. A ligand carried over with a
                        # non-blank altLoc (e.g. a resolved alternate conformation
                        # left at occupancy 0.5) is parsed by BioPython as a
                        # DisorderedAtom; setting .name on that wrapper is silently
                        # dropped on write, so the rename must reach each underlying
                        # child Atom.
                        if atom_mapping:
                            for atom in residue:
                                children = (
                                    atom.disordered_get_list()
                                    if atom.is_disordered()
                                    else [atom]
                                )
                                for child in children:
                                    if child.name in atom_mapping:
                                        new_name = atom_mapping[child.name]
                                        child.name = new_name
                                        child.id = new_name
                                        child.fullname = f" {new_name:<3}"  # PDB format
                                        atoms_renamed += 1

        if atoms_renamed == 0 and residues_renamed == 0:
            return None

        # Save remapped structure
        output = self._output_dir / "remapped.pdb"
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output))

        self.console.print(f"[grey50]Renamed {atoms_renamed} atoms, {residues_renamed} residues → {output.name}[/grey50]")
        return output

    def _get_lib_residue_name(self, lib_file: str) -> Optional[str]:
        """Extract residue name from lib file (preserves case)."""
        try:
            with open(lib_file, 'r') as f:
                for line in f:
                    # Look for: !entry.RESNAME.unit.atoms
                    if '!entry.' in line and '.unit.' in line:
                        # Extract RESNAME from !entry.RESNAME.unit.XXX
                        parts = line.split('.')
                        if len(parts) >= 2:
                            return parts[1]  # Returns the residue name with original case
            return None
        except Exception:
            return None

    def _configure_tleap_for_assembly(self, pdb_file: Path) -> None:
        """Configure tLEaP with all collected FF information."""
        if not self.workspace:
            return

        # Store structure for tLEaP
        self.workspace.set("preprocessing_protein_input", str(pdb_file))

        # For preprocessing, we always use implicit solvent (no bulk water added).
        # Explicit solvation is done later, when the production topology is built.
        # Mark this as PROVISIONAL: it is not a user solvation choice, only a
        # placeholder so the metal-free preprocessing tleap builds without a box.
        # The Topology Generator treats a provisional solvation as unconfigured
        # and still runs the solvation prompt (otherwise this shared key would
        # silently satisfy its `if not solvation_params` guard and the production
        # topology would be built with no periodic box).
        # Coordinating waters (if any) are kept as part of the structure
        self.workspace.set("solvation_parameters",
                           {'solvent_model': 'implicit', 'provisional': True})
        self.console.print("[grey50]Using implicit solvent (no bulk solvation box)[/grey50]")

        # Set up standard forcefields for tLEaP module
        protein_ff = self.workspace.get("preprocessing_protein_ff", "")
        water_model = self.workspace.get("preprocessing_water_model", "")

        tleap_selection = {}
        if protein_ff:
            name = protein_ff.split('.')[-1] if '.' in protein_ff else protein_ff
            tleap_selection["protein"] = {"name": name, "leaprc": protein_ff}
        if water_model:
            name = water_model.split('.')[-1] if '.' in water_model else water_model
            tleap_selection["water"] = {"name": name, "leaprc": water_model, "box": "none"}

        self.workspace.set("selected_standard_forcefields", tleap_selection)

        # Collect lib/frcmod files for small molecules
        lib_files = []
        frcmod_files = []

        # For each residue prefer its .lib (loadoff) over its .mol2: tLEaP matches
        # loadpdb residues to templates by the lib's OFF entry name, and a bare
        # `loadmol2 <file>` with no assignment registers no template at all.
        # Only fall back to the mol2 when no lib is available.
        def _add_component(info):
            lib = info.get('lib_file')
            mol2 = info.get('mol2_file')
            if lib:
                lib_files.append(lib)
            elif mol2:
                lib_files.append(mol2)
            if info.get('frcmod_file'):
                frcmod_files.append(info['frcmod_file'])

        # Organic small molecules
        organic_ff = self.workspace.get("preprocessing_organic_ff", {})
        for resname, info in organic_ff.items():
            _add_component(info)

        # Organometallic (those with full params or generated params)
        orgmet_ff = self.workspace.get("preprocessing_organometallic_ff", {})
        for resname, info in orgmet_ff.items():
            _add_component(info)

        # Store for tLEaP
        if lib_files:
            self.workspace.set("preprocessing_lib_files", lib_files)
        if frcmod_files:
            self.workspace.set("preprocessing_frcmod_files", frcmod_files)

        self.console.print(f"[grey50]Configured tLEaP with {len(lib_files)} lib files, {len(frcmod_files)} frcmod files[/grey50]")

    def _run_tleap_assembly(self) -> bool:
        """Run Topology Generator for structure assembly."""
        if not self.processor:
            self.console.print("[red]No processor available for tLEaP[/red]")
            return False

        tleap_gen = self.processor.get_module_instance("Topology Generator")
        if not tleap_gen:
            self.console.print("[red]Topology Generator not available[/red]")
            return False

        self.console.print("\n[bold]Running Topology Generator[/bold]")

        # Run tLEaP from the output directory so all files land there
        import os
        orig_dir = os.getcwd()
        if self._output_dir:
            os.chdir(self._output_dir)

        # Signal the Topology Generator that THIS whole run is preprocessing's own
        # metal-free build: keep the provisional implicit solvation as-is (no box,
        # no prompt) for the entire assembly. The solvation flow is reached via
        # handle_menu_option("generate_single_state"), NOT only process(), so the
        # flag must cover both — the user's solvation prompt belongs to the later
        # production topology run, where the provisional flag is cleared instead.
        self.workspace.set("_preprocessing_tleap_active", True)
        try:
            # Initial processing
            if hasattr(tleap_gen, 'process'):
                tleap_gen.process(self.workspace)

            # Generate single-state tLEaP input
            if hasattr(tleap_gen, 'handle_menu_option'):
                self.console.print("[cyan]Generating tLEaP input...[/cyan]")
                tleap_gen.handle_menu_option("generate_single_state")

                self.console.print("[cyan]Running tLEaP...[/cyan]")
                tleap_gen.handle_menu_option("generate_topology")

            # Check for output
            parm7 = self.workspace.get("parm7_file") if self.workspace else None
            rst7 = self.workspace.get("rst7_file") if self.workspace else None

            if parm7 and rst7:
                self.console.print(f"[green]✓ tLEaP produced {Path(parm7).name}[/green]")
                return True
            else:
                self.console.print("[yellow]tLEaP did not produce expected output files[/yellow]")
                return False

        except Exception as e:
            # The message alone can be a bare KeyError key — "tLEaP error:
            # 'site_id'" says nothing about where it came from, and the step
            # that follows then fails on the missing output instead. Keep the
            # traceback so the origin is recoverable.
            self.console.print(f"[red]tLEaP error: {e}[/red]")
            logger.exception("Topology Generator raised during structure recombination")
            return False
        finally:
            self.workspace.set("_preprocessing_tleap_active", False)
            os.chdir(orig_dir)

    def _insert_metals(self, tleap_pdb: Path, metals: List[MetalInfo]) -> Path:
        """Insert metals back into the tLEaP output structure."""
        from Bio.PDB import PDBParser, PDBIO
        from Bio.PDB.Atom import Atom
        from Bio.PDB.Residue import Residue

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("struct", tleap_pdb)
        model = structure[0]

        # Get residue sequence map for finding new residue numbers
        seq_map = self.workspace.get("preprocessing_residue_sequence_map", {}) if self.workspace else {}

        # Track metal reinsertion mapping: (orig_chain, orig_resid) -> (new_chain, new_resid)
        metal_reinsertion_map = {}

        # Find max residue number for isolated metals
        max_resid = 0
        for chain in model:
            for residue in chain:
                max_resid = max(max_resid, residue.id[1])

        # Pure metal clusters (category F) reinsert as ONE residue holding all
        # their atoms, not one residue per atom. Group by cluster_id, build the
        # residue, and drop these atoms from the per-atom loop below.
        cluster_groups: Dict[str, List[MetalInfo]] = {}
        singles: List[MetalInfo] = []
        for metal in metals:
            if metal.cluster_id:
                cluster_groups.setdefault(metal.cluster_id, []).append(metal)
            else:
                singles.append(metal)

        for cluster_id, atoms in cluster_groups.items():
            max_resid += 1
            ref = atoms[0]
            if ref.original_chain in [c.id for c in model]:
                chain = model[ref.original_chain]
            else:
                chain = list(model.get_chains())[0]

            new_res = Residue((' ', max_resid, ' '), ref.original_resname, '')
            for a in atoms:
                new_res.add(Atom(
                    a.atom_name, a.coords, 0.0, 1.0, ' ',
                    a.atom_name, max_resid, a.element.upper()
                ))
            chain.add(new_res)
            metal_reinsertion_map[(ref.original_chain, ref.original_resid)] = (chain.id, max_resid)
            self.console.print(
                f"  [green]Inserted metal cluster: {ref.original_resname} "
                f"({len(atoms)} atoms) as residue {max_resid}[/green]")

        for metal in singles:
            if metal.is_isolated:
                # Create new residue for isolated metal
                max_resid += 1

                # Find or use first chain
                if metal.original_chain in [c.id for c in model]:
                    chain = model[metal.original_chain]
                else:
                    chain = list(model.get_chains())[0]

                # Create new residue and atom
                new_res = Residue((' ', max_resid, ' '), metal.original_resname, '')
                new_atom = Atom(
                    metal.atom_name,
                    metal.coords,
                    0.0,  # bfactor
                    1.0,  # occupancy
                    ' ',  # altloc
                    metal.atom_name,
                    max_resid,
                    metal.element.upper()  # BioPython requires uppercase element
                )
                new_res.add(new_atom)
                chain.add(new_res)

                # Track reinsertion mapping
                metal_reinsertion_map[(metal.original_chain, metal.original_resid)] = (chain.id, max_resid)

                self.console.print(f"  [green]Inserted isolated metal: {metal.original_resname} as residue {max_resid}[/green]")

            else:
                # Embedded metal - find the residue and add atom to it
                # Use sequence map to find new resid
                orig_key = f"({metal.original_chain}, {metal.original_resid})"
                new_resid = None

                # Try to find in sequence map
                for key, pos in seq_map.items():
                    if str(key) == orig_key or key == (metal.original_chain, metal.original_resid):
                        # Find residue with this sequence position
                        new_resid = pos
                        break

                if new_resid is None:
                    # Fallback: try to find by resname
                    for chain in model:
                        for residue in chain:
                            if residue.get_resname().strip() == metal.original_resname:
                                new_resid = residue.id[1]
                                break

                if new_resid:
                    # Find residue and add metal atom
                    for chain in model:
                        for residue in chain:
                            if residue.id[1] == new_resid:
                                new_atom = Atom(
                                    metal.atom_name,
                                    metal.coords,
                                    0.0,
                                    1.0,
                                    ' ',
                                    metal.atom_name,
                                    new_resid,
                                    metal.element.upper()  # BioPython requires uppercase
                                )
                                residue.add(new_atom)
                                # Track reinsertion mapping
                                metal_reinsertion_map[(metal.original_chain, metal.original_resid)] = (chain.id, new_resid)
                                self.console.print(f"  [green]Inserted embedded metal: {metal.element} into {metal.original_resname}[/green]")
                                break
                else:
                    self.console.print(f"  [yellow]Could not find residue for {metal.original_resname}[/yellow]")

        # Store metal reinsertion map in workspace for redox site sync
        if self.workspace and metal_reinsertion_map:
            self.workspace.set("preprocessing_metal_reinsertion_map", metal_reinsertion_map)

        output = self._output_dir / "prepared_structure.pdb"
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output))

        return output

    def _extract_atom_data_from_prmtop(self, parm7_file: str) -> Dict[tuple, Dict[str, Any]]:
        """Extract atom types and charges from prmtop file.

        The prmtop IS the force field applied to the system. It contains
        all atom types and charges assigned by tLEaP, so we don't need
        to re-load leaprcs separately.

        Note: Metal atoms won't be in the prmtop (they were removed before
        tLEaP and reinserted after). Metal charges must be obtained from
        the user in MCPB Step 1.

        Returns:
            Dict mapping (chain, resid, resname, atom_name) -> {'type': str, 'charge': float}
        """
        import parmed

        atom_data = {}
        try:
            # Load prmtop (topology only, no coordinates needed)
            parm = parmed.load_file(parm7_file)

            for atom in parm.atoms:
                # Use residue/atom info as key since coordinates aren't in prmtop
                res = atom.residue
                # tLEaP doesn't preserve chain IDs, so use empty string
                key = ('', res.idx + 1, res.name, atom.name)
                atom_data[key] = {
                    'type': atom.type,
                    'charge': atom.charge
                }

        except Exception as e:
            self.console.print(f"[yellow]Could not extract data from prmtop: {e}[/yellow]")

        return atom_data

    # =========================================================================
    # Legacy Component Processing Methods
    # =========================================================================

    def _process_protein_component(self, protein_residues: List[str]) -> Optional[Path]:
        """
        Process protein component: FF selection + protonation + tLEaP.

        This combines FF selection, H addition, and atom typing for protein.
        tLEaP does all three in one step.

        Returns:
            Path to protein_with_H.pdb, or None if failed
        """
        self.console.print(Panel(
            "[bold]Protein Parameterization[/bold]\n"
            "Select force field, analyze protonation states, add hydrogens",
            border_style="cyan",
            expand=False
        ))

        # 1. Select protein force field
        protein_ff = self._select_protein_forcefield()
        if not protein_ff:
            return None

        # Store selection
        if self.workspace:
            self.workspace.set("preprocessing_protein_ff", protein_ff)

        # 2. Extract protein to temporary PDB
        protein_pdb = self._output_dir / "protein_only.pdb"
        self._extract_residues_to_pdb(self._pdb_file, protein_residues, protein_pdb)
        self.console.print(f"[grey50]Extracted {len(protein_residues)} protein residues[/grey50]")

        # 3. Store residue mapping BEFORE tLEaP (which may renumber)
        residue_sequence_map = self._build_residue_sequence_map(protein_pdb)
        if self.workspace and residue_sequence_map:
            self.workspace.set("preprocessing_residue_sequence_map", residue_sequence_map)

        # 4. Set up tLEaP FF selection for the tLEaP module
        self._set_tleap_ff_selection()

        # 5. Store protein PDB in workspace
        if self.workspace:
            self.workspace.set("preprocessing_protein_input", str(protein_pdb))

        # 6. Launch protonation state analyzer
        if self.processor:
            prot_analyzer = self.processor.get_module_instance("Protonation State Analyzer")
            if prot_analyzer:
                self.console.print("\n[bold]Protonation State Analysis[/bold]")
                if hasattr(prot_analyzer, 'analyze_protonation_states'):
                    prot_analyzer.analyze_protonation_states()
                if hasattr(prot_analyzer, 'set_residue_names'):
                    prot_analyzer.set_residue_names()

        # 7. Launch tLEaP
        if self.processor:
            tleap_gen = self.processor.get_module_instance("Topology Generator")
            if tleap_gen:
                self.console.print("\n[bold]tLEaP Processing[/bold]")
                # Preprocessing's own metal-free tleap: keep provisional implicit
                # solvation across the WHOLE run (process + generate_single_state,
                # which is where the solvation flow actually runs). See
                # _run_tleap_assembly for the rationale.
                self.workspace.set("_preprocessing_tleap_active", True)
                try:
                    if hasattr(tleap_gen, 'process'):
                        tleap_gen.process(self.workspace)
                    if hasattr(tleap_gen, 'handle_menu_option'):
                        tleap_gen.handle_menu_option("generate_single_state")
                        tleap_gen.handle_menu_option("generate_topology")
                finally:
                    self.workspace.set("_preprocessing_tleap_active", False)

        # 8. Convert parm7/rst7 to PDB
        parm7 = self.workspace.get("parm7_file") if self.workspace else None
        rst7 = self.workspace.get("rst7_file") if self.workspace else None

        # Clear preprocessing key
        if self.workspace:
            self.workspace.set("preprocessing_protein_input", None)

        if parm7 and rst7:
            protein_with_h = self._output_dir / "protein_with_H.pdb"
            self._convert_amber_to_pdb(parm7, rst7, protein_with_h)
            self.console.print(f"[green]✓ Created {protein_with_h.name}[/green]")
            return protein_with_h
        else:
            self.console.print("[yellow]tLEaP output not found[/yellow]")
            return protein_pdb

    def _process_water_component(self, water_residues: List[str]) -> Optional[Path]:
        """
        Process water component: FF selection + tLEaP.

        For coordinating waters that need explicit parameterization.

        Returns:
            Path to water_with_H.pdb, or None if skipped/failed
        """
        self.console.print(Panel(
            "[bold]Water Parameterization[/bold]\n"
            "Select water model for coordinating waters",
            border_style="cyan",
            expand=False
        ))

        # Select water model
        water_model = self._select_water_model()
        if not water_model:
            self.console.print("[grey50]Skipping water parameterization[/grey50]")
            return None

        # Store selection
        if self.workspace:
            self.workspace.set("preprocessing_water_model", water_model)

        # TODO: Process waters through tLEaP
        # For now, waters are handled with protein or kept as-is
        self.console.print("[grey50]Water processing integrated with protein tLEaP[/grey50]")

        return None

    def _process_nonstandard_residue(self, res_key: str) -> Optional[Path]:
        """
        Process a single non-standard residue.

        Two options:
        1. Have parameters → provide lib/frcmod, run tLEaP for H atoms
        2. Need parameters → launch Small Molecule Parameterizer

        Returns:
            Path to residue_with_H.pdb, or None if failed
        """
        parts = res_key.split(':')
        chain = parts[0]
        resid = int(parts[1])
        resname = parts[2] if len(parts) > 2 else "UNK"

        self.console.print(f"\n[bold cyan]Non-standard residue: {resname} ({chain}:{resid})[/bold cyan]")

        # Two simple options
        self.console.print("  [1] I have parameters (lib/frcmod files)")
        self.console.print("  [2] I need to generate parameters")

        choice = prompt_with_context(
            self.processor,
            f"Select option for {resname}",
            choices=["1", "2"],
            module="Structure Preprocessor",
            description=f"Parameter source for {resname}"
        )

        if choice == "1":
            # Have parameters - get lib/frcmod and run tLEaP
            return self._process_nonstandard_with_params(res_key, resname)
        else:
            # Generate with Small Molecule Parameterizer
            result = self._launch_small_molecule_parameterizer(resname, res_key, self.triage_results)
            if result:
                # Store parameterization result
                if self.workspace:
                    ns_ff = self.workspace.get("preprocessing_nonstandard_ff", {})
                    ns_ff[resname] = {
                        'source': 'generated',
                        'mol2_file': result.get('mol2_file'),
                        'frcmod_file': result.get('frcmod_file'),
                    }
                    self.workspace.set("preprocessing_nonstandard_ff", ns_ff)
                # TODO: Return path to H-added structure from parameterizer
                return None
            return None

    def _process_nonstandard_with_params(self, res_key: str, resname: str) -> Optional[Path]:
        """
        Process non-standard residue with user-provided lib/frcmod files.

        Prompts for lib and frcmod files, then runs tLEaP to add hydrogens.

        Returns:
            Path to residue_with_H.pdb, or None if failed
        """
        # Prompt for lib/frcmod files
        lib_path, frcmod_path, ff_resname = self._prompt_custom_ff_files(resname)

        if not lib_path:
            self.console.print(f"  [yellow]No lib file provided for {resname}[/yellow]")
            return None

        # Store selection
        if self.workspace:
            ns_ff = self.workspace.get("preprocessing_nonstandard_ff", {})
            ns_ff[resname] = {
                'source': 'provided',
                'ff_resname': ff_resname or resname,
                'lib_file': str(lib_path),
                'frcmod_file': str(frcmod_path) if frcmod_path else None,
            }
            self.workspace.set("preprocessing_nonstandard_ff", ns_ff)

        # Extract residue to PDB
        parts = res_key.split(':')
        chain, resid = parts[0], int(parts[1])
        input_pdb = self._output_dir / f"{resname}_{chain}_{resid}_input.pdb"
        self._extract_single_residue_to_pdb(self._pdb_file, chain, resid, resname, input_pdb)
        self.console.print(f"  [grey50]Extracted to {input_pdb.name}[/grey50]")
        self.console.print(f"  [grey50]Parameters: lib={lib_path}, frcmod={frcmod_path}[/grey50]")

        # Run tLEaP with lib/frcmod to add hydrogens
        output_pdb = self._output_dir / f"{resname}_{chain}_{resid}_H.pdb"
        success = self._run_tleap_for_residue(
            input_pdb=input_pdb,
            output_pdb=output_pdb,
            leaprc='leaprc.gaff2',  # Base FF for atom types
            lib_file=str(lib_path),
            frcmod_file=str(frcmod_path) if frcmod_path else None,
            ff_resname=ff_resname,
        )

        if success:
            return output_pdb
        else:
            self.console.print(f"  [yellow]tLEaP failed, returning unprocessed residue[/yellow]")
            return input_pdb

    def _process_metal_ion(self, res_key: str) -> bool:
        """
        Process a metal ion: lookup in MetalIonDatabase.

        Prompts for charge/spin if needed, then looks up VDW parameters
        and atom type from the database.

        Returns:
            True if successfully typed
        """
        from Bio.PDB import PDBParser

        parts = res_key.split(':')
        chain = parts[0]
        resid = int(parts[1])
        resname = parts[2] if len(parts) > 2 else "UNK"

        self.console.print(f"\n  [bold]Metal ion: {resname} ({chain}:{resid})[/bold]")

        # Get atom name from the PDB
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("temp", self._pdb_file)
        atom_name = resname  # Default

        for model in structure:
            for pdb_chain in model:
                if pdb_chain.id == chain:
                    for residue in pdb_chain:
                        if residue.id[1] == resid:
                            atoms = list(residue.get_atoms())
                            if atoms:
                                atom_name = atoms[0].name.strip()
                            break

        # Import MetalIonDatabase
        try:
            from proprep.forcefield_prep.mcpb.metal_ion_database import (
                MetalIonDatabase, MetalConfig, water_model_from_leaprc,
            )

            # Get water model from workspace (if selected). Matched on the
            # leaprc suffix: a substring test sent opc3 to the OPC set.
            water_model = water_model_from_leaprc(
                self.workspace.get("preprocessing_water_model", "") if self.workspace else "",
                logger=logger,
            )

            db = MetalIonDatabase(water_model=water_model)

            # Check if metal is recognized
            if db.is_metal(resname, atom_name):
                self.console.print(f"  [green]✓ Recognized metal: {resname}[/green]")

                # Prompt for charge
                self.console.print("  Common charges for this metal:")
                self.console.print("    Enter formal charge (e.g., 2 for 2+, -1 for anion)")

                charge_str = prompt_with_context(
                    self.processor,
                    f"Formal charge for {resname}",
                    default="2",
                    module="Structure Preprocessor",
                    description=f"Formal charge for {resname}"
                )
                charge = int(charge_str)

                # Prompt for spin
                self.console.print("  Number of unpaired electrons (for high-spin vs low-spin)")
                spin_str = prompt_with_context(
                    self.processor,
                    f"Unpaired electrons for {resname}",
                    default="0",
                    module="Structure Preprocessor",
                    description=f"Spin state for {resname}"
                )
                spin = int(spin_str)

                # Get configuration from database
                config = db.get_metal_config(resname, atom_name, charge, spin)

                if config:
                    self.console.print(f"  [green]✓ Found parameters:[/green]")
                    self.console.print(f"    Element: {config.element}")
                    self.console.print(f"    Atom type: {config.atom_type}")
                    self.console.print(f"    Charge: {config.charge}+")
                    self.console.print(f"    VDW radius: {config.vdw_radius:.3f} Å")
                    self.console.print(f"    VDW epsilon: {config.vdw_epsilon:.6f} kcal/mol")

                    # Store in workspace
                    if self.workspace:
                        metal_types = self.workspace.get("preprocessing_metal_types", {})
                        metal_types[res_key] = {
                            'element': config.element,
                            'atom_name': config.atom_name,
                            'atom_type': config.atom_type,
                            'charge': config.charge,
                            'spin': config.spin,
                            'mass': config.mass,
                            'vdw_radius': config.vdw_radius,
                            'vdw_epsilon': config.vdw_epsilon,
                            'water_model': config.water_model,
                        }
                        self.workspace.set("preprocessing_metal_types", metal_types)

                    return True
                else:
                    self.console.print(f"  [yellow]⚠ No VDW parameters for {resname} with charge {charge}+[/yellow]")
            else:
                self.console.print(f"  [yellow]⚠ Metal {resname} not in database[/yellow]")

        except ImportError as e:
            self.console.print(f"  [yellow]MetalIonDatabase not available: {e}[/yellow]")

        # Fallback: store basic info
        if self.workspace:
            metal_types = self.workspace.get("preprocessing_metal_types", {})
            metal_types[res_key] = {
                'element': resname,
                'atom_name': atom_name,
                'atom_type': resname,  # Use resname as type
                'charge': None,
                'spin': None,
            }
            self.workspace.set("preprocessing_metal_types", metal_types)

        return True

    def _collect_atom_types_from_components(self) -> Dict[Tuple, 'AtomTypeAssignment']:
        """
        Collect atom types from all parameterized components.

        Aggregates types that were assigned during param-1 through param-4.
        Does NOT re-assign types - just collects what was already assigned.

        Returns:
            Dict mapping coordinate tuples to AtomTypeAssignment
        """
        from Bio.PDB import PDBParser

        type_assignments = {}

        self.console.print(Panel(
            "[bold]Collecting Atom Types[/bold]\n"
            "Aggregating types from all parameterized components",
            border_style="cyan",
            expand=False
        ))

        # We need to read the final structure and assign types based on what
        # was determined in the param steps

        if not self._final_pdb or not self._final_pdb.exists():
            self.console.print("[yellow]No final structure available[/yellow]")
            return type_assignments

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("final", str(self._final_pdb))

        # Counters for summary
        protein_count = 0
        water_count = 0
        nonstandard_count = 0
        metal_count = 0
        untyped_count = 0

        # Get stored type information from workspace
        metal_types = self.workspace.get("preprocessing_metal_types", {}) if self.workspace else {}
        nonstandard_ff = self.workspace.get("preprocessing_nonstandard_ff", {}) if self.workspace else {}

        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    resid = residue.id[1]
                    res_key = f"{chain.id}:{resid}:{resname}"

                    for atom in residue.get_atoms():
                        coords = tuple(atom.coord)
                        atom_name = atom.name.strip()
                        element = atom.element.strip() if atom.element else ""

                        # Determine atom type based on residue category
                        atom_type = None
                        source = AtomSource.UNKNOWN

                        # Check if it's a metal
                        if res_key in metal_types:
                            metal_info = metal_types[res_key]
                            atom_type = metal_info.get('atom_type', resname)
                            source = AtomSource.METAL_DATABASE
                            metal_count += 1

                        # Check if it's a standard protein residue
                        elif resname in STANDARD_RESIDUES:
                            # Type comes from protein FF (assigned by tLEaP)
                            # For now, use a placeholder - real types come from prmtop
                            atom_type = f"{resname}_{atom_name}"  # Placeholder
                            source = AtomSource.FORCEFIELD_LIB
                            protein_count += 1

                        # Check if it's water
                        elif resname in WATER_RESIDUES:
                            # Type from water model
                            if atom_name in ('O', 'OW'):
                                atom_type = 'OW'
                            elif atom_name.startswith('H'):
                                atom_type = 'HW'
                            else:
                                atom_type = atom_name
                            source = AtomSource.FORCEFIELD_LIB
                            water_count += 1

                        # Check if it's a parameterized non-standard residue
                        elif resname in nonstandard_ff:
                            ns_info = nonstandard_ff[resname]
                            # Type would come from mol2/lib file
                            atom_type = f"{resname}_{atom_name}"  # Placeholder
                            if ns_info.get('source') == 'generated':
                                source = AtomSource.ANTECHAMBER
                            else:
                                source = AtomSource.FORCEFIELD_LIB
                            nonstandard_count += 1

                        else:
                            # Unknown/untyped
                            atom_type = None
                            source = AtomSource.UNKNOWN
                            untyped_count += 1

                        # Create assignment
                        assignment = AtomTypeAssignment(
                            coords=coords,
                            original_type=atom_type,
                            renamed_type=None,  # Set later in MCPB steps
                            source=source,
                            chain=chain.id,
                            resname=resname,
                            resid=resid,
                            atom_name=atom_name,
                            element=element,
                        )
                        type_assignments[coords] = assignment

        # Display summary
        table = Table(title="Atom Type Collection Summary", show_lines=True)
        table.add_column("Source", style="cyan")
        table.add_column("Count", style="white")

        table.add_row("Protein (FF)", str(protein_count))
        table.add_row("Water", str(water_count))
        table.add_row("Non-standard", str(nonstandard_count))
        table.add_row("Metals", str(metal_count))
        if untyped_count > 0:
            table.add_row("[yellow]Untyped[/yellow]", f"[yellow]{untyped_count}[/yellow]")
        table.add_row("[bold]Total[/bold]", f"[bold]{len(type_assignments)}[/bold]")

        self.console.print(table)

        if untyped_count > 0:
            self.console.print(f"[yellow]⚠ {untyped_count} atoms without assigned types[/yellow]")
            self.console.print("[grey50]These may be H-caps or residues that weren't parameterized[/grey50]")

        return type_assignments

    def _extract_single_residue_to_pdb(
        self,
        pdb_file: str,
        chain_id: str,
        resid: int,
        resname: str,
        output_path: Path
    ) -> None:
        """Extract a single residue to a PDB file."""
        from Bio.PDB import PDBParser, PDBIO, Select

        class ResidueSelect(Select):
            def __init__(self, target_chain, target_resid):
                self.target_chain = target_chain
                self.target_resid = target_resid

            def accept_residue(self, residue):
                return (residue.get_parent().id == self.target_chain and
                        residue.id[1] == self.target_resid)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("temp", pdb_file)

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_path), ResidueSelect(chain_id, resid))

    # =========================================================================
    # Legacy Entry Point (non-checklist)
    # =========================================================================

    def run_preprocessing(
        self,
        pdb_file: str,
        output_dir: Path,
        interactive: bool = True
    ) -> PreprocessingResult:
        """
        Run the complete preprocessing workflow.

        Note: Does NOT take RedoxSite objects as input - they are detected
        on the final structure in Step 0f.

        Args:
            pdb_file: Path to input PDB file
            output_dir: Directory for output files
            interactive: Whether to prompt user for input

        Returns:
            PreprocessingResult containing:
            - prepared_pdb: Path to final prepared structure
            - redox_sites: Fresh RedoxSite objects (detected on final structure)
            - type_assignments: Dict of coords → AtomTypeAssignment (with original_type)
            - ff_data: ForceFieldData
            - small_mol_results: Dict of small molecule param results
        """
        self.console.print(Panel(
            "[bold cyan]Structure Preprocessing[/bold cyan]\n"
            "Preparing structure for metal site parameterization",
            border_style="cyan",
            expand=False
        ))

        # Display preprocessing pipeline overview
        self._display_pipeline_overview()

        try:
            # Step 0a: Structure filtering (optional, via PDB Filter)
            self.console.print("\n[bold]Step 0a: Structure Filtering[/bold]")
            pdb_file = self._step_0a_structure_filtering(pdb_file, interactive)

            # Step 0b: Structure triage (reads Analysis results from workspace)
            self.console.print("\n[bold]Step 0b: Structure Triage[/bold]")
            self.triage_results = self._step_0b_structure_triage(pdb_file, interactive)

            # Step 0c: Hydrogen addition
            self.console.print("\n[bold]Step 0c: Hydrogen Addition[/bold]")
            h_results = self._step_0c_hydrogen_addition(
                pdb_file, self.triage_results, output_dir, interactive
            )

            # Step 0d: Atom exclusion + H capping (manual specification)
            self.console.print("\n[bold]Step 0d: Atom Exclusion + H Capping[/bold]")
            self.excluded_atoms, self.h_caps = self._step_0d_atom_exclusion(
                pdb_file, h_results, output_dir, interactive
            )

            # Step 0e: Structure recombination
            self.console.print("\n[bold]Step 0e: Structure Recombination[/bold]")
            final_pdb = self._step_0e_recombination(
                pdb_file, h_results, self.excluded_atoms, self.h_caps,
                self.triage_results, output_dir
            )

            # Step 0f: Launch Redox Site Detector on final structure
            self.console.print("\n[bold]Step 0f: Redox Site Detection[/bold]")
            self.redox_sites = self._step_0f_redox_site_sync(final_pdb, interactive)

            # Step 0g: FF selection + complete atom typing (on final structure)
            self.console.print("\n[bold]Step 0g: FF Selection + Atom Typing[/bold]")
            self.ff_data, self.type_assignments = self._step_0g_complete_atom_typing(
                final_pdb, self.redox_sites, output_dir, interactive
            )

            # Success!
            self.console.print(Panel(
                f"[green]Preprocessing complete![/green]\n"
                f"Final structure: {final_pdb}\n"
                f"Redox sites: {len(self.redox_sites)}\n"
                f"Atoms typed: {len(self.type_assignments)}",
                border_style="green",
                expand=False
            ))

            return PreprocessingResult(
                success=True,
                prepared_pdb=final_pdb,
                redox_sites=self.redox_sites,
                type_assignments=self.type_assignments,
                ff_data=self.ff_data,
                residue_map=self.residue_map,
                atom_map=self.atom_map,
                small_mol_results=self.small_mol_results,
                excluded_atoms=self.excluded_atoms,
                h_caps=self.h_caps,
                triage=self.triage_results
            )

        except Exception as e:
            self.console.print(f"[red]Preprocessing failed: {e}[/red]")
            return PreprocessingResult(
                success=False,
                error_message=str(e)
            )

    # =========================================================================
    # Pipeline Overview
    # =========================================================================

    def _display_pipeline_overview(self) -> None:
        """Display the preprocessing pipeline steps as a Rich table."""
        table = Table(
            title="Preprocessing Pipeline Overview",
            show_lines=False,
            title_style="bold cyan",
        )
        table.add_column("Step", style="cyan", no_wrap=True, width=6)
        table.add_column("Description", style="white")

        steps = [
            ("0a", "Structure Filtering (optional)"),
            ("0b", "Structure Triage + Force Field Selection"),
            ("0c", "Hydrogen Addition"),
            ("0d", "Atom Exclusion + H Capping"),
            ("0e", "Structure Recombination"),
            ("0f", "Redox Site Sync"),
            ("0g", "Atom Typing"),
        ]
        for step_id, description in steps:
            table.add_row(step_id, description)

        self.console.print(table)

    # =========================================================================
    # Step 0a: Structure Filtering
    # =========================================================================

    def _step_0a_structure_filtering(
        self,
        pdb_file: str,
        interactive: bool
    ) -> str:
        """
        Step 0a: Structure filtering via the PDB Filter module.

        Allows the user to select specific chains and remove unwanted
        components (e.g., extra chains, small molecules, waters) before
        the rest of preprocessing.  Redox site detection is skipped here
        because Step 0f will run fresh detection on the final structure.

        Goes straight into the filter: running this step IS the request to
        filter. The step is optional at the checklist level, where ``<num>s``
        marks it skipped, so a confirmation here only asked the same question
        a second time.

        Args:
            pdb_file: Path to input PDB file
            interactive: Whether to prompt the user

        Returns:
            Path to the (possibly filtered) PDB file for downstream steps
        """
        if not interactive:
            self.console.print("[grey50]Non-interactive mode: skipping structure filtering[/grey50]")
            return pdb_file

        if not self.processor:
            self.console.print("[yellow]Processor not available, skipping filtering[/yellow]")
            return pdb_file

        # Launch PDB Filter with redox detection disabled
        pdb_filter = self.processor.get_module_instance("PDB Filter")
        if not pdb_filter:
            self.console.print(
                "[yellow]PDB Filter module not available, skipping filtering[/yellow]"
            )
            return pdb_file

        self.console.print("\n[bold]═══ PDB Filter ═══[/bold]")
        self.console.print(
            "[grey50]Select chains and components to keep. "
            "Redox detection is handled later in Step 0f.[/grey50]\n"
        )

        pdb_filter.filter_pdb_structure(
            interactive=True,
            skip_redox_detection=True,
        )

        # Check if filtering produced a filtered PDB
        if self.workspace:
            filtered_pdb = self.workspace.get("filtered_pdb_file")
            if filtered_pdb and Path(filtered_pdb).exists():
                self.console.print(
                    f"[green]Using filtered structure: "
                    f"{Path(filtered_pdb).name}[/green]"
                )
                return filtered_pdb

        self.console.print("[grey50]No filtered structure produced, using original[/grey50]")
        return pdb_file

    def _step_0a2_structure_completeness(
        self,
        pdb_file: str,
        interactive: bool,
    ) -> str:
        """
        Step 0a2: Optional structure completeness via the Structure
        Completeness module.

        Drives the shipped ``Structure Completeness`` module end to end --
        the same two menu actions the main menu exposes: ``analyze`` then
        (after a bail prompt) ``process_structure`` (apply repairs / altloc
        selections / caps). We reuse it wholesale rather than reimplementing
        because it already owns the residue/chain mapping that keeps numbering
        consistent when MODELLER fills gaps.

        This resolves alternate conformations (e.g. a bound inhibitor modelled
        in two conformers) before parameterization, so antechamber/MCPB don't
        choke on duplicate atoms downstream. The module writes
        ``repaired_pdb_file``; we adopt it as the working structure. Fresh
        redox detection later (Step 0f) re-runs on these coordinates, so any
        renumbering is naturally absorbed.

        Args:
            pdb_file: Path to the current working PDB (filtered if Step 0a ran)
            interactive: Whether to prompt the user

        Returns:
            Path to the (possibly repaired) PDB for downstream steps
        """
        if not interactive:
            self.console.print(
                "[grey50]Non-interactive mode: skipping structure completeness[/grey50]"
            )
            return pdb_file

        if not self.processor or not self.workspace:
            self.console.print(
                "[yellow]Processor/workspace not available, skipping completeness[/yellow]"
            )
            return pdb_file

        # Registered module NAME is "Structure Fixer" (the completeness module).
        # Import it first so its @register_module decorator has run regardless
        # of app startup import order, then fetch the shared instance.
        try:
            import proprep.structure_prep.structure_completeness  # noqa: F401
        except Exception:
            pass
        module = self.processor.get_module_instance("Structure Fixer")
        if not module:
            self.console.print(
                "[yellow]Structure Fixer module not available, skipping[/yellow]"
            )
            return pdb_file

        # Point the module at the CURRENT working structure. filtered_structure
        # is its priority-1, non-interactive source, so setting it here avoids a
        # mid-workflow "select a structure" prompt and guarantees it analyzes
        # exactly this file (whether or not Step 0a filtering ran).
        try:
            from Bio.PDB import PDBParser
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("completeness_input", pdb_file)
            self.workspace.set("filtered_structure", structure)
        except Exception as exc:
            self.console.print(
                f"[yellow]Could not load structure for completeness: {exc}[/yellow]"
            )
            return pdb_file

        self.console.print("\n[bold]═══ Structure Completeness ═══[/bold]")

        # 1) Analyze -- same as the main menu's "Analyze structure completeness".
        module.handle_menu_option("analyze")

        if not getattr(module, "results", None):
            self.console.print(
                "[grey50]No completeness analysis produced; proceeding.[/grey50]"
            )
            return pdb_file

        # 2) Bail point between analysis and repair. This opens the repair
        # session, where the per-segment and altloc choices are still to be
        # made — "apply" described a decision that has not been taken yet.
        apply_repairs = confirm_with_context(
            self.processor,
            "Repair missing segments and resolve alternate conformations now?",
            default=True,
            module="Structure Preprocessor",
            description="Enter the structure completeness repair session",
        )
        if not apply_repairs:
            self.console.print(
                "[grey50]Analysis kept; structure left unrepaired.[/grey50]")
            return pdb_file

        # 3) Apply -- same as the main menu's "Apply repairs/mutations/caps".
        module.handle_menu_option("process_structure")

        repaired = self.workspace.get("repaired_pdb_file")
        if repaired and Path(repaired).exists():
            self.console.print(
                f"[green]Using repaired structure: {Path(repaired).name}[/green]"
            )
            self.workspace.set("structure_pdb_file", repaired)
            return repaired

        self.console.print("[grey50]No repaired structure produced, using current[/grey50]")
        return pdb_file

    # =========================================================================
    # Step 0b: Structure Triage
    # =========================================================================

    def _step_0b_structure_triage(
        self,
        pdb_file: str,
        interactive: bool
    ) -> Dict[str, str]:
        """
        Step 0b: Structure triage + Force Field selection.

        Part 1 - Triage:
        Reads classification from workspace (set by forcefield_parameterizer Analysis step).
        If not available, LAUNCHES the Analysis step (don't re-implement).

        Categories:
        - A: Standard protein residues → protonation_state_analyzer + tLEaP
        - B: Non-standard residues (ALL) → reduce for H's
        - C: Waters → keep (with option to selectively remove)
        - D: Metal ions → keep, handled in atom typing

        Part 2 - Force Field Selection:
        After triage, selects force fields for each component type.
        Stores selections in workspace for use by Step 0c (tLEaP) and Step 0g (typing).

        Returns:
            Dict mapping residue_id (chain:resid:resname) → category
        """
        # === Part 1: Triage ===
        self.console.print("[bold cyan]Part 1: Structure Classification[/bold cyan]")

        # Try to read Analysis results from workspace
        # The forcefield_parameterizer stores results as "non_standard_residues"
        non_standard_residues = None
        if self.workspace:
            non_standard_residues = self.workspace.get("non_standard_residues")

        if non_standard_residues:
            self.console.print("[green]✓ Using classification from forcefield Analysis step[/green]")
            triage = self._convert_analysis_to_triage(non_standard_residues, pdb_file)
        else:
            # Fallback: LAUNCH the Analysis step (don't re-implement)
            self.console.print("[yellow]No Analysis results found, launching Analysis...[/yellow]")
            triage = self._launch_analysis_step(pdb_file, interactive)

        if not triage:
            return triage

        # === Part 2: Force Field Selection ===
        self.console.print("\n[bold cyan]Part 2: Force Field Selection[/bold cyan]")
        self._step_0b_ff_selection(triage, interactive)

        return triage

    def _convert_analysis_to_triage(
        self, non_standard_residues: List, pdb_file: str
    ) -> Dict[str, str]:
        """
        Convert forcefield_parameterizer Analysis results to triage format.

        The Analysis provides NonStandardResidue objects with category info.
        We scan the PDB to identify which residues actually exist (important
        when a filtered structure is used) and classify them.
        """
        from Bio.PDB import PDBParser

        triage = {}

        # First, scan the structure to get ALL residues that actually exist
        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure("structure", pdb_file)
        except Exception as e:
            self.console.print(f"[red]Failed to parse PDB: {e}[/red]")
            return {}

        # Build set of residue keys present in the (possibly filtered) structure
        residues_in_structure = set()
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.resname.strip()
                    resid = residue.id[1]
                    chain_id = chain.id
                    res_key = f"{chain_id}:{resid}:{resname}"
                    residues_in_structure.add(res_key)

                    # Classify standard residues, waters, metals directly
                    # A=protein, B=organic, C=organometallic, D=isolated metal, E=water
                    if resname in STANDARD_RESIDUES:
                        triage[res_key] = 'A'
                    elif resname in WATER_RESIDUES:
                        triage[res_key] = 'E'
                    elif self._is_metal_residue(residue):
                        triage[res_key] = 'D'
                    # Non-standard will be handled from non_standard_residues list
            break  # Only first model

        # Now process non_standard_residues from Analysis, but ONLY if they
        # exist in the current structure (important for filtered structures)
        for res in non_standard_residues:
            res_key = f"{res.chain_id}:{res.resid}:{res.name}"

            # Skip residues that were filtered out
            if res_key not in residues_in_structure:
                continue

            category = getattr(res, 'category', 'unknown')

            if category == 'metal_site':
                triage[res_key] = 'D'  # Metals
            else:
                # ALL other non-standard (small_molecule, modified_amino_acid, unknown)
                # go to category B for reduce H addition
                triage[res_key] = 'B'

        self._display_triage_summary(triage)
        return triage

    def _launch_analysis_step(self, pdb_file: str, interactive: bool) -> Dict[str, str]:
        """Launch forcefield_parameterizer Analysis step to get classification."""
        # Store PDB in workspace
        if self.workspace:
            self.workspace.set("structure_pdb_file", pdb_file)

        # Try to launch Analysis via forcefield_parameterizer
        if self.processor:
            ff_param = self.processor.get_module_instance("Force Field Parameterizer")
            if ff_param and hasattr(ff_param, 'analyze_nonstandard_residues'):
                self.console.print("\n[bold]═══ Force Field Parameterizer Analysis ═══[/bold]")
                self.console.print("[grey50]Running analysis to classify residues...[/grey50]")
                ff_param.analyze_nonstandard_residues()

                # Check if results are now in workspace (stored as "non_standard_residues")
                if self.workspace:
                    non_standard_residues = self.workspace.get("non_standard_residues")
                    if non_standard_residues:
                        self.console.print(f"[green]✓ Analysis found {len(non_standard_residues)} non-standard residues[/green]")
                        return self._convert_analysis_to_triage(non_standard_residues, pdb_file)

        self.console.print("[red]Analysis step failed - could not retrieve results[/red]")
        return {}

    def _is_metal_residue(self, residue) -> bool:
        """Check if a residue is a metal ion."""
        atoms = list(residue.get_atoms())
        if len(atoms) == 1:
            element = atoms[0].element.strip().upper() if atoms[0].element else ''
            return element in METAL_ELEMENTS
        return False

    def _handle_water_selection(
        self,
        triage: Dict[str, str],
        waters: List[str],
        pdb_file: str
    ) -> Dict[str, str]:
        """Allow user to selectively keep/remove waters."""
        self.console.print(f"\n[cyan]Found {len(waters)} water molecules[/cyan]")

        keep_all = confirm_with_context(
            self.processor,
            "Keep all waters?",
            default=True,
            module="Structure Preprocessor",
            description="Keep all waters",
        )
        if keep_all:
            return triage

        # Let user specify waters to remove
        self.console.print("Enter water residue IDs to REMOVE (format: CHAIN:RESID:HOH)")
        self.console.print("Type 'done' when finished.")

        while True:
            res_spec = prompt_with_context(
                self.processor,
                "Remove water (or 'done')",
                default="done",
                module="Structure Preprocessor",
                description="Water residue ID to remove",
            )
            if res_spec.lower() == 'done':
                break
            if res_spec in waters:
                triage[res_spec] = 'REMOVE'
                self.console.print(f"[yellow]Marked {res_spec} for removal[/yellow]")
            else:
                self.console.print(f"[red]Water not found: {res_spec}[/red]")

        return triage

    def _display_triage_summary(self, triage: Dict[str, str]) -> None:
        """Display a detailed summary table of triage results."""
        # Group residues by category
        # A=protein, B=organic, C=organometallic, D=isolated metal, E=water
        category_residues = {'A': [], 'B': [], 'C': [], 'D': [], 'E': [], 'REMOVE': []}
        for res_key, cat in triage.items():
            if cat in category_residues:
                # res_key format: "chain:resid:resname"
                parts = res_key.split(':')
                if len(parts) == 3:
                    chain, resid, resname = parts
                    category_residues[cat].append(f"{resname}:{chain}{resid}")
                else:
                    category_residues[cat].append(res_key)

        # Category descriptions
        descriptions = {
            'A': 'Standard protein',
            'B': 'Organic small molecules',
            'C': 'Organometallic',
            'D': 'Isolated metal ions',
            'E': 'Waters',
            'REMOVE': 'Marked for removal'
        }

        table = Table(title="Structure Triage Summary", show_lines=True)
        table.add_column("Cat", style="cyan", width=4)
        table.add_column("Description", width=22)
        table.add_column("Residues", style="grey50")
        table.add_column("Count", justify="right", width=5)

        for cat in ['A', 'B', 'C', 'D', 'E', 'REMOVE']:
            residues = category_residues[cat]
            count = len(residues)

            if count == 0:
                continue

            # For category A (protein), just show count - too many to list
            if cat == 'A':
                residue_str = f"[grey50]{count} standard amino acids[/grey50]"
            elif count <= 10:
                # Show all residues if 10 or fewer
                residue_str = ", ".join(sorted(set(residues)))
            else:
                # Show unique residue names with counts
                from collections import Counter
                res_counts = Counter(r.split(':')[0] for r in residues)
                residue_str = ", ".join(f"{name}({c})" for name, c in res_counts.most_common())

            table.add_row(cat, descriptions[cat], residue_str, str(count))

        self.console.print(table)

    # =========================================================================
    # Step 0b Part 2: Force Field Selection (after triage)
    # =========================================================================

    def _step_0b_ff_selection(
        self,
        triage: Dict[str, str],
        interactive: bool
    ) -> None:
        """
        Force field selection based on triage results.

        Called after triage to select force fields for each component type:
        - Protein (Category A): Select protein FF (e.g., ff19SB)
        - Water (Category C): Select water model (e.g., OPC)
        - Non-standard (Category B): Per-residue selection

        Stores selections in workspace for use by Step 0c (tLEaP) and Step 0g (typing).

        NOTE: This always prompts interactively. ProPrep is transparent - no hidden defaults.
        """
        self.console.print(Panel(
            "[cyan]Force Field Selection[/cyan]\n"
            "Select force fields for each component type",
            border_style="cyan",
            expand=False
        ))

        # Check what categories are present
        # A=protein, B=organic, C=organometallic, D=isolated metal, E=water
        has_protein = any(v == 'A' for v in triage.values())
        has_water = any(v == 'E' for v in triage.values())
        has_nonstandard = any(v == 'B' for v in triage.values())

        # 1. Protein force field
        if has_protein:
            self._select_protein_forcefield()

        # 2. Water model
        if has_water:
            self._select_water_model()

        # 3. Non-standard residues (per-residue selection)
        if has_nonstandard:
            nonstandard_residues = self._get_unique_nonstandard_residues(triage)
            self._select_nonstandard_forcefields(nonstandard_residues, triage)

        # 4. Metal ions - inform user (typed from MetalIonDatabase in Step 0g)
        has_metals = any(v == 'D' for v in triage.values())
        if has_metals:
            self.console.print(
                "\n[grey50]Metal ions will be typed from MetalIonDatabase in Step 0g.[/grey50]"
            )

        # Display summary
        self._display_ff_selection_summary()

    def _select_protein_forcefield(self) -> str:
        """Select force field for standard protein residues.

        Renders the curated catalog through the shared force-field menu so this
        picker matches the Topology Generator's exactly (headerless table,
        green names, blue divider, recommendation markers).

        Returns:
            The selected leaprc path (e.g., 'leaprc.protein.ff14SB')
        """
        from proprep.forcefield_params.forcefield_menu import select_single_forcefield

        selected = select_single_forcefield(
            self.console, self.processor, "protein",
            module="Structure Preprocessor",
            prompt_label="Select protein force field",
        )

        selected_leaprc = selected['leaprc']
        selected_name = selected['name']

        # Store in workspace
        if self.workspace:
            self.workspace.set("preprocessing_protein_ff", selected_leaprc)

        self.console.print(f"[green]✓ Selected: {selected_name}[/green]")

        # Ask about modified amino acids (optional)
        self._select_modified_aa_options(selected_name)

        return selected_leaprc

    def _select_modified_aa_options(self, protein_ff_name: str) -> None:
        """Select modified amino acid parameters (optional, multi-select).

        Rendered through the shared force-field menu so it matches the other
        pickers. Options whose name matches the chosen protein FF are marked
        recommended (yellow marker) in place of the old ad-hoc green highlight.
        """
        from proprep.forcefield_params.forcefield_menu import render_forcefield_category
        from proprep.forcefield_params.forcefield_catalog import FORCEFIELD_OPTIONS

        cat_info = FORCEFIELD_OPTIONS.get('modified_aa', {})
        mod_aa_opts = cat_info.get('options', [])
        if not mod_aa_opts:
            return

        # Shallow-copy so FF-matching options can be marked recommended without
        # mutating the shared catalog.
        mod_aa_opts = [dict(opt) for opt in mod_aa_opts]
        for opt in mod_aa_opts:
            if protein_ff_name and protein_ff_name.lower() in opt['name'].lower():
                opt['recommended'] = True
                opt['recommendation_reason'] = f"matches {protein_ff_name}"

        # "None" is row 1 (allow_none); real options follow.
        none_offset = render_forcefield_category(
            self.console,
            f"{cat_info.get('title', 'MODIFIED AMINO ACIDS')}  (optional)",
            cat_info.get('description'),
            mod_aa_opts,
            allow_none=True,
            none_text="None -- no modified AA parameters",
        )
        self.console.print(
            "[grey50]Enter 1 for none, or comma-separated choices (e.g., 2,3), "
            "or press Enter for none[/grey50]"
        )

        # options_map keeps the session recorder legible.
        mod_aa_map = {"1": "None"}
        for i, opt in enumerate(mod_aa_opts, start=1 + none_offset):
            mod_aa_map[str(i)] = opt['name']

        choice = prompt_with_context(
            self.processor,
            "Select modified AA parameters",
            default="1",
            module="Structure Preprocessor",
            description="Select modified AA FF parameters",
            options_map=mod_aa_map,
        )

        # Row 1 / empty / legacy 'n' all mean none.
        if not choice.strip() or choice.strip().lower() in ('1', 'n', 'none'):
            self.console.print("[grey50]No modified AA parameters selected[/grey50]")
            return

        # Parse multi-select (indices are 1-based including the None row).
        selected_leaprcs = []
        for idx_str in choice.split(','):
            try:
                idx = int(idx_str.strip()) - 1 - none_offset
                if 0 <= idx < len(mod_aa_opts):
                    selected_leaprcs.append(mod_aa_opts[idx]['leaprc'])
                    self.console.print(f"[green]✓ Selected: {mod_aa_opts[idx]['name']}[/green]")
            except ValueError:
                pass

        if selected_leaprcs and self.workspace:
            self.workspace.set("preprocessing_modified_aa", selected_leaprcs)

    def _select_water_model(self) -> str:
        """Select water model.

        Renders the curated catalog through the shared force-field menu so this
        picker matches the Topology Generator's exactly.

        Returns:
            Selected leaprc/frcmod string, or None if selection failed.
        """
        from proprep.forcefield_params.forcefield_menu import select_single_forcefield

        selected = select_single_forcefield(
            self.console, self.processor, "water",
            module="Structure Preprocessor",
            prompt_label="Select water model",
            description="Must be compatible with protein force field!",
        )

        # Handle both 'leaprc' and 'frcmod' keys (some water models use frcmod)
        selected_leaprc = selected.get('leaprc') or selected.get('frcmod')
        selected_name = selected['name']

        if self.workspace:
            self.workspace.set("preprocessing_water_model", selected_leaprc)

        self.console.print(f"[green]✓ Selected: {selected_name}[/green]")
        return selected_leaprc

    def _get_unique_nonstandard_residues(self, triage: Dict[str, str]) -> List[Dict]:
        """Get unique non-standard residue names with info."""
        unique = {}  # resname → {chain, resid, res_key}

        for res_key, category in triage.items():
            if category != 'B':
                continue

            parts = res_key.split(':')
            if len(parts) != 3:
                continue

            chain, resid, resname = parts

            if resname not in unique:
                unique[resname] = {
                    'resname': resname,
                    'chain': chain,
                    'resid': resid,
                    'res_key': res_key,
                    'count': 1
                }
            else:
                unique[resname]['count'] += 1

        return list(unique.values())

    def _select_nonstandard_forcefields(
        self,
        nonstandard_residues: List[Dict],
        triage: Dict[str, str]
    ) -> None:
        """Per-residue force field selection for non-standard residues.

        For each non-standard residue, user chooses:
        1. Existing: Use parameters already in AMBER force field
        2. Custom: Provide their own lib/frcmod files
        3. Generate: Launch Small Molecule Parameterizer NOW to generate parameters
        4. Skip: Handle separately (not typed in Step 0g)

        When "Generate" is selected, Small Molecule Parameterizer is launched
        immediately, not deferred to a later step.
        """
        from proprep.forcefield_prep.forcefield_data import ForceFieldData

        self.console.print("\n[bold]Non-Standard Residue Force Fields[/bold]")
        self.console.print(
            "[grey50]For each non-standard residue, select how to obtain force field parameters.[/grey50]\n"
        )

        # Try to load a temporary ForceFieldData to check for existing residues
        try:
            ff_check = ForceFieldData(console=self.console, processor=self.processor)
            # Load common extended force fields to check for heme, etc.
            try:
                ff_check.load_leaprc('leaprc.phosaa10')  # Phosphorylated AA
            except:
                pass
            try:
                ff_check.load_leaprc('leaprc.gaff2')  # GAFF2
            except:
                pass
        except:
            ff_check = None

        selections = {}

        for res_info in nonstandard_residues:
            resname = res_info['resname']
            count = res_info['count']
            chain = res_info['chain']
            resid = res_info['resid']
            res_key = res_info['res_key']

            self.console.print(f"\n[cyan]{resname}[/cyan] ({count} instance(s), e.g., {chain}:{resid})")

            # Check if residue exists in any loaded FF
            exists_in_ff = ff_check.has_residue(resname) if ff_check else False

            # Display options
            options = []
            if exists_in_ff:
                options.append(("existing", f"Use existing FF definition ({resname})", "(Found in force field)"))
            options.append(("custom", "Custom: Provide frcmod + lib files", "(You have parameter files)"))
            options.append(("generate", "Generate: Launch Small Molecule Parameterizer", "(antechamber + GAFF2)"))
            options.append(("skip", "Skip: Handle separately later", "(Not included in typing)"))

            for i, (key, label, note) in enumerate(options, 1):
                if key == "existing":
                    self.console.print(f"  [{i}] {label} [green]{note}[/green]")
                else:
                    self.console.print(f"  [{i}] {label} [grey50]{note}[/grey50]")

            # No default - user must explicitly choose
            choice = prompt_with_context(
                self.processor,
                f"Selection for {resname}",
                choices=[str(i) for i in range(1, len(options) + 1)]
            )

            selected_option = options[int(choice) - 1][0]

            if selected_option == "existing":
                selections[resname] = {
                    'source': 'existing',
                    'ff_resname': resname,
                    'lib_file': None,
                    'frcmod_file': None
                }
                self.console.print(f"  [green]→ Using existing FF definition[/green]")

            elif selected_option == "custom":
                lib_path, frcmod_path, mapped_name = self._prompt_custom_ff_files(resname)
                selections[resname] = {
                    'source': 'custom',
                    'ff_resname': mapped_name or resname,
                    'lib_file': lib_path,
                    'frcmod_file': frcmod_path
                }

            elif selected_option == "generate":
                # Launch Small Molecule Parameterizer NOW
                self.console.print(f"\n[bold cyan]═══ Launching Small Molecule Parameterizer for {resname} ═══[/bold cyan]")
                result = self._launch_small_molecule_parameterizer(resname, res_key, triage)

                if result:
                    selections[resname] = {
                        'source': 'generated',
                        'ff_resname': resname,
                        'lib_file': result.get('lib_file'),
                        'frcmod_file': result.get('frcmod_file'),
                        'mol2_file': result.get('mol2_file')
                    }
                    self.console.print(f"  [green]→ Parameters generated successfully[/green]")
                else:
                    selections[resname] = {
                        'source': 'generate_failed',
                        'ff_resname': None,
                        'lib_file': None,
                        'frcmod_file': None
                    }
                    self.console.print(f"  [red]→ Parameter generation failed[/red]")

            else:  # skip
                selections[resname] = {
                    'source': 'skip',
                    'ff_resname': None,
                    'lib_file': None,
                    'frcmod_file': None
                }
                self.console.print(f"  [grey50]→ Skipped[/grey50]")

        # Store selections
        if self.workspace:
            self.workspace.set("preprocessing_nonstandard_ff", selections)

    def _launch_small_molecule_parameterizer(
        self,
        resname: str,
        res_key: str,
        triage: Dict[str, str]
    ) -> Optional[Dict]:
        """
        Launch Small Molecule Parameterizer for a residue.

        Extracts the residue from the PDB, runs the parameterizer workflow,
        and returns the generated files (mol2, lib, frcmod).

        Args:
            resname: Residue name (e.g., "MNS")
            res_key: Residue key like "A:401" (chain:resid)
            triage: Triage dictionary with residue categories

        Returns:
            Dict with 'mol2_file', 'lib_file', 'frcmod_file' paths, or None if failed
        """
        import os
        from pathlib import Path
        from Bio.PDB import PDBParser

        # Get the current PDB file
        pdb_file = None
        if self.workspace:
            pdb_file = self.workspace.get("structure_pdb_file")
            if not pdb_file:
                pdb_file = self.workspace.get("filtered_pdb_file")

        if not pdb_file:
            self.console.print("[red]No PDB file available in workspace[/red]")
            return None

        pdb_path = Path(pdb_file)
        if not pdb_path.exists():
            self.console.print(f"[red]PDB file not found: {pdb_file}[/red]")
            return None

        # Parse chain and resid from res_key (e.g., "A:862:MNS" or "A:401")
        try:
            parts = res_key.split(':')
            if len(parts) >= 2:
                chain_id = parts[0]
                target_resid = int(parts[1])
            else:
                raise ValueError(f"Expected at least chain:resid, got: {res_key}")
        except ValueError as e:
            self.console.print(f"[red]Invalid residue key format: {res_key} ({e})[/red]")
            return None

        # Parse the PDB and extract the target residue
        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure("temp", str(pdb_path))
        except Exception as e:
            self.console.print(f"[red]Failed to parse PDB: {e}[/red]")
            return None

        # Find the target residue and create wrapper objects
        # (following the pattern from forcefield_parameterizer.py)
        biopython_residues = []

        class ResidueWrapper:
            """Wrapper for BioPython residue with chain/resid metadata."""
            def __init__(self, biopython_residue, chain_id, resid):
                self.biopython_residue = biopython_residue
                self.chain_id = chain_id
                self.resid = resid

            def get_atoms(self):
                return self.biopython_residue.get_atoms()

            def get_resname(self):
                return self.biopython_residue.get_resname()

        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        res_resid = residue.id[1]
                        res_resname = residue.get_resname().strip()
                        if res_resid == target_resid and res_resname == resname:
                            wrapper = ResidueWrapper(residue, chain_id, res_resid)
                            biopython_residues.append(wrapper)

        if not biopython_residues:
            self.console.print(f"[red]Could not find residue {resname} at {res_key} in structure[/red]")
            return None

        self.console.print(f"  Found {len(biopython_residues)} instance(s) of {resname}")

        # Set up output directory
        original_dir = os.getcwd()
        output_dir = f"small_molecule_params_{resname}"

        try:
            # Import and run the small molecule parameterizer workflow
            from proprep.forcefield_prep import small_molecule_parameterizer

            self.console.print(f"[cyan]Starting small molecule parameterization workflow...[/cyan]")

            workflow_result = small_molecule_parameterizer.run_workflow(
                residue_name=resname,
                residues=biopython_residues,
                output_dir=output_dir,
                interactive=True,
                processor=self.processor,
                # This launcher is only reached from the "[2] I need to generate
                # parameters" branch, so never silently reuse pre-existing output
                # (e.g. on a recorded-session replay) — regenerate.
                regenerate=True,
            )

            # Process workflow results
            if workflow_result.get("success"):
                # Get parameter files from the result
                mol2_file = workflow_result.get("parameter_files", {}).get("prep_file")
                frcmod_file = workflow_result.get("parameter_files", {}).get("frcmod_file")

                # Check if workflow paused (needs Gaussian)
                if workflow_result.get("status") == "paused":
                    self.console.print(
                        f"\n[yellow]Parameterization paused - Gaussian calculation required[/yellow]"
                    )
                    self.console.print(
                        f"[grey50]Run Gaussian, then re-run preprocessing to continue[/grey50]"
                    )
                    # Store pending state
                    if self.workspace:
                        pending = self.workspace.get("pending_parameterizations", {})
                        pending[resname] = {
                            "output_dir": output_dir,
                            "status": "pending_gaussian",
                            "res_key": res_key
                        }
                        self.workspace.set("pending_parameterizations", pending)
                    return None

                # Return successful results
                if mol2_file or frcmod_file:
                    self.console.print(f"[green]✓ Parameterization successful for {resname}[/green]")

                    # Load atom name mapping if available
                    atom_mapping = None
                    mapping_path = workflow_result.get("parameter_files", {}).get("atom_mapping")
                    if mapping_path and os.path.exists(mapping_path):
                        import json
                        with open(mapping_path) as f:
                            atom_mapping = json.load(f)

                    # Prefer the real .lib (OFF) file: the recombination loadoff's
                    # it, and tLEaP matches loadpdb residues by the lib entry name.
                    lib_file = workflow_result.get("parameter_files", {}).get("lib_file")

                    return {
                        'mol2_file': mol2_file,
                        'lib_file': lib_file,
                        'frcmod_file': frcmod_file,
                        'atom_name_mapping': atom_mapping,
                    }
                else:
                    self.console.print(f"[yellow]Workflow completed but no parameter files found[/yellow]")
                    return None
            elif workflow_result.get("status") == "saved":
                self.console.print(f"[yellow]Parameterization saved (incomplete) for {resname}[/yellow]")
                return None
            else:
                error_msg = workflow_result.get("message", "Unknown error")
                self.console.print(f"[red]Parameterization failed: {error_msg}[/red]")
                return None

        except ImportError as e:
            self.console.print(f"[red]Failed to import small_molecule_parameterizer: {e}[/red]")
            return None
        except Exception as e:
            self.console.print(f"[red]Error running Small Molecule Parameterizer: {e}[/red]")
            import traceback
            self.console.print(f"[grey50]{traceback.format_exc()}[/grey50]")
            return None
        finally:
            # Return to original directory
            os.chdir(original_dir)

    def _prompt_custom_ff_files(self, resname: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """
        Prompt user for custom force field files using interactive file browser.

        Uses ProPrep's standard file browser pattern with directory navigation.
        """
        import os

        self.console.print(f"\n  [cyan]Custom FF files for {resname}:[/cyan]")

        # FF residue name (may differ from PDB)
        ff_resname = prompt_with_context(
            self.processor,
            f"  FF residue name in parameter files",
            default=resname,
            module="Structure Preprocessor",
            description=f"Force-field residue name for {resname} in parameter files",
        )

        # Get working directory for file browser
        working_dir = os.getcwd()
        if self.workspace:
            working_dir = self.workspace.get("working_directory", working_dir)

        # Browse for lib/off file (required)
        self.console.print("\n  [bold]Select lib/off file:[/bold]")
        lib_path = self._browse_ff_files(
            directory=working_dir,
            extensions=['.lib', '.off'],
            file_type="lib/off"
        )

        if not lib_path:
            self.console.print("  [yellow]No lib file selected[/yellow]")
            return (None, None, ff_resname if ff_resname != resname else None)

        # Browse for frcmod file (optional)
        self.console.print("\n  [bold]Select frcmod file (optional):[/bold]")
        frcmod_path = self._browse_ff_files(
            directory=os.path.dirname(lib_path) if lib_path else working_dir,
            extensions=['.frcmod'],
            file_type="frcmod",
            optional=True
        )

        self.console.print(f"\n  [green]→ Custom files configured[/green]")
        if lib_path:
            self.console.print(f"    lib: {os.path.basename(lib_path)}")
        if frcmod_path:
            self.console.print(f"    frcmod: {os.path.basename(frcmod_path)}")

        return (
            lib_path,
            frcmod_path,
            ff_resname if ff_resname != resname else None
        )

    def _handle_atom_name_mapping(
        self, res_key: str, lib_path: str, resname: str
    ) -> Optional[Dict[str, str]]:
        """
        Handle atom name mapping between PDB and lib file.

        Checks for existing mapping file, or prompts user to create one manually.

        Args:
            res_key: Residue key (chain:resid:resname)
            lib_path: Path to the lib file
            resname: Residue name

        Returns:
            Dict mapping PDB atom names to lib atom names, or None if not needed
        """
        import os
        import json
        from Bio.PDB import PDBParser

        # First, check for existing mapping file
        lib_dir = os.path.dirname(lib_path)
        mapping_file = os.path.join(lib_dir, f"{resname.lower()}_atom_name_mapping.json")

        if os.path.exists(mapping_file):
            self.console.print(f"  [green]Found atom name mapping: {os.path.basename(mapping_file)}[/green]")
            with open(mapping_file, 'r') as f:
                return json.load(f)

        # Get atoms from PDB structure
        pdb_atoms = []
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("temp", self._pdb_file)
        parts = res_key.split(':')
        chain_id, resid = parts[0], int(parts[1])

        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        if residue.id[1] == resid:
                            for atom in residue:
                                pdb_atoms.append((atom.name, atom.element))
                            break

        # Get atoms from lib file
        lib_atoms = self._parse_lib_atoms(lib_path)

        if not lib_atoms:
            self.console.print(f"  [yellow]Could not parse lib file atoms[/yellow]")
            return None

        # Check if names match (no mapping needed)
        pdb_names = set(a[0] for a in pdb_atoms)
        lib_names = set(a[0] for a in lib_atoms)

        if pdb_names == lib_names:
            self.console.print(f"  [green]Atom names match - no mapping needed[/green]")
            return {}

        # A crystal structure usually has no hydrogens while a library always
        # does, so comparing raw name sets reports a mismatch for a library
        # that fits perfectly -- and then asks the user to hand-map 53 atoms
        # against 84. tLEaP adds the missing hydrogens from the template, so
        # only the heavy atoms have to correspond.
        pdb_heavy = {n for n in pdb_names if not _is_hydrogen_name(n)}
        lib_heavy = {n for n in lib_names if not _is_hydrogen_name(n)}

        if pdb_heavy and pdb_heavy == lib_heavy and not (
                pdb_names - pdb_heavy):
            missing = len(lib_names) - len(lib_heavy)
            self.console.print(
                f"  [green]Heavy-atom names match ({len(pdb_heavy)}/{len(pdb_heavy)}) "
                f"- no mapping needed[/green]")
            if missing:
                self.console.print(
                    f"  [grey50]The PDB has no hydrogens; tLEaP will add the "
                    f"{missing} in the library template.[/grey50]")
            return {}

        # Names don't match - prompt for manual mapping
        self.console.print(f"\n  [yellow]Atom name mismatch detected![/yellow]")
        self.console.print(f"  PDB has {len(pdb_atoms)} atoms, lib has {len(lib_atoms)} atoms")
        if pdb_heavy != lib_heavy:
            only_pdb = sorted(pdb_heavy - lib_heavy)
            only_lib = sorted(lib_heavy - pdb_heavy)
            if only_pdb:
                self.console.print(
                    f"  [grey50]heavy atoms only in the PDB: "
                    f"{', '.join(only_pdb[:8])}"
                    f"{' ...' if len(only_pdb) > 8 else ''}[/grey50]")
            if only_lib:
                self.console.print(
                    f"  [grey50]heavy atoms only in the lib: "
                    f"{', '.join(only_lib[:8])}"
                    f"{' ...' if len(only_lib) > 8 else ''}[/grey50]")

        mapping = self._prompt_manual_atom_mapping(pdb_atoms, lib_atoms, resname, lib_dir)
        return mapping

    def _parse_lib_atoms(self, lib_path: str) -> List[Tuple[str, str, str]]:
        """
        Parse atom names from a lib file.

        Returns:
            List of (atom_name, atom_type, element) tuples
        """
        atoms = []
        in_atoms_section = False

        try:
            with open(lib_path, 'r') as f:
                for line in f:
                    # ``.unit.atoms table`` and not ``.unit.atomspertinfo``:
                    # the substring test matched both, and atomspertinfo has one
                    # row per atom carrying the same names, so every library was
                    # read at double length. An 84-atom FAD reported 168, its
                    # listing repeating from index 85.
                    if _LIB_ATOMS_TABLE_RE.search(line):
                        in_atoms_section = True
                        continue
                    if in_atoms_section:
                        if line.startswith('!'):
                            break
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            # Format: "ATOMNAME" "type" ...
                            atom_name = parts[0].strip('"')
                            atom_type = parts[1].strip('"') if len(parts) > 1 else ""
                            # Element is typically first 1-2 chars of atom name
                            element = ''.join(c for c in atom_name if c.isalpha())[:2]
                            atoms.append((atom_name, atom_type, element))
        except Exception as e:
            pass

        return atoms

    def _prompt_manual_atom_mapping(
        self,
        pdb_atoms: List[Tuple[str, str]],
        lib_atoms: List[Tuple[str, str, str]],
        resname: str,
        output_dir: str
    ) -> Dict[str, str]:
        """
        Display atoms from PDB and lib side by side, prompt for manual mapping.

        Args:
            pdb_atoms: List of (name, element) from PDB
            lib_atoms: List of (name, type, element) from lib
            resname: Residue name for output file
            output_dir: Directory to save mapping file

        Returns:
            Dict mapping PDB atom names to lib atom names
        """
        import json
        import os
        from rich.table import Table
        from rich.columns import Columns

        self.console.print("\n  [bold]Manual Atom Name Mapping Required[/bold]")
        self.console.print("  [grey50]Create mapping between PDB atom names and lib file atom names[/grey50]\n")

        # Create two tables side by side
        pdb_table = Table(title="PDB Atoms", show_header=True, header_style="bold cyan")
        pdb_table.add_column("Idx", width=4)
        pdb_table.add_column("Name", width=8)
        pdb_table.add_column("Element", width=8)

        for i, (name, element) in enumerate(pdb_atoms, 1):
            pdb_table.add_row(str(i), name, element or "?")

        lib_table = Table(title="LIB Atoms", show_header=True, header_style="bold green")
        lib_table.add_column("Idx", width=4)
        lib_table.add_column("Name", width=8)
        lib_table.add_column("Type", width=8)

        for i, (name, atype, element) in enumerate(lib_atoms, 1):
            lib_table.add_row(str(i), name, atype)

        # Display tables
        self.console.print(Columns([pdb_table, lib_table], equal=True, expand=True))

        self.console.print("\n  Enter mappings as space-separated pairs: PDB_idx:LIB_idx")
        self.console.print("  Example: 5:9 6:10 7:5  (PDB atom 5 → LIB atom 9, etc.)")
        self.console.print("  Enter 'all' if all atoms have identical names (no mapping needed)")

        mapping_input = prompt_with_context(
            self.processor,
            "  Mappings",
            module="Structure Preprocessor",
            description=f"Atom name mapping for {resname}"
        )

        if mapping_input.lower() == 'all':
            self.console.print("  [green]Using identical names - no mapping needed[/green]")
            return {}

        # Parse the mapping input
        mapping = {}
        for pair in mapping_input.split():
            try:
                pdb_idx, lib_idx = pair.split(':')
                pdb_idx = int(pdb_idx) - 1
                lib_idx = int(lib_idx) - 1

                if 0 <= pdb_idx < len(pdb_atoms) and 0 <= lib_idx < len(lib_atoms):
                    pdb_name = pdb_atoms[pdb_idx][0]
                    lib_name = lib_atoms[lib_idx][0]
                    mapping[pdb_name] = lib_name
                    self.console.print(f"    {pdb_name} → {lib_name}")
            except (ValueError, IndexError):
                self.console.print(f"    [yellow]Invalid pair: {pair}[/yellow]")

        # Save the mapping
        if mapping:
            mapping_file = os.path.join(output_dir, f"{resname.lower()}_atom_name_mapping.json")
            with open(mapping_file, 'w') as f:
                json.dump(mapping, f, indent=2)
            self.console.print(f"  [green]✓ Mapping saved to {os.path.basename(mapping_file)}[/green]")

        return mapping

    def _browse_ff_files(
        self,
        directory: str,
        extensions: List[str],
        file_type: str,
        optional: bool = False
    ) -> Optional[str]:
        """
        Interactive file browser for force field files.

        Thin wrapper over the shared file browser: unified bare-N / q UX and
        filename-based session replay. When ``optional`` is set, the user may
        ``skip``; skip and cancel both return None (callers treat them alike).

        Args:
            directory: Starting directory
            extensions: List of valid file extensions (e.g., ['.lib', '.off'])
            file_type: Description for display (e.g., "lib/off")
            optional: If True, allow skipping without selecting a file

        Returns:
            Selected file path or None if cancelled/skipped
        """
        from proprep.utils.file_browser import file_browser, default_size_detail, SKIP

        result = file_browser(
            directory=directory,
            extensions=extensions,
            console=self.console,
            processor=self.processor,
            label=f"{file_type} file",
            entry_detail=default_size_detail,
            optional=optional,
            module="Structure Preprocessor",
        )
        if result is None or result is SKIP:
            return None
        return result



    def _display_ff_selection_summary(self) -> None:
        """Display a summary of force field selections."""
        if not self.workspace:
            return

        self.console.print("\n[bold]Force Field Selection Summary[/bold]")

        table = Table(show_header=True, header_style="bold", show_lines=True)
        table.add_column("Component", width=20)
        table.add_column("Selection", width=40)

        # Protein FF
        protein_ff = self.workspace.get("preprocessing_protein_ff")
        if protein_ff:
            table.add_row("Protein", protein_ff)

        # Water model
        water_model = self.workspace.get("preprocessing_water_model")
        if water_model:
            table.add_row("Water", water_model)

        # Non-standard
        nonstandard_ff = self.workspace.get("preprocessing_nonstandard_ff")
        if nonstandard_ff:
            for resname, selection in nonstandard_ff.items():
                source = selection.get('source', 'unknown')
                if source == 'existing':
                    ff_resname = selection.get('ff_resname', resname)
                    table.add_row(f"{resname}", f"Existing FF ({ff_resname})")
                elif source == 'custom':
                    lib_file = selection.get('lib_file', 'N/A')
                    table.add_row(f"{resname}", f"Custom: {lib_file}")
                elif source == 'generated':
                    mol2 = selection.get('mol2_file')
                    if mol2:
                        from pathlib import Path
                        table.add_row(f"{resname}", f"[green]Generated: {Path(mol2).name}[/green]")
                    else:
                        table.add_row(f"{resname}", "[green]Generated (antechamber)[/green]")
                elif source == 'generate_failed':
                    table.add_row(f"{resname}", "[red]Generation failed[/red]")
                elif source == 'skip':
                    table.add_row(f"{resname}", "[grey50]Skipped[/grey50]")
                else:
                    table.add_row(f"{resname}", f"[grey50]{source}[/grey50]")

        # Metals
        table.add_row("Metal ions", "[grey50]MetalIonDatabase[/grey50]")

        self.console.print(table)

    def _set_tleap_ff_selection(self) -> None:
        """
        Convert Step 0b FF selections to tLEaP's expected workspace format.

        Topology Generator reads `selected_standard_forcefields` from workspace.
        It expects dicts with 'leaprc' key (and 'box' for water models).
        """
        if not self.workspace:
            return

        protein_ff = self.workspace.get("preprocessing_protein_ff")
        water_model = self.workspace.get("preprocessing_water_model")

        if not protein_ff:
            return  # Let tLEaP prompt if no selection made

        # Build tLEaP-compatible selection dict
        # tLEaP expects: {"protein": {"leaprc": "leaprc.protein.ff19SB"}, ...}
        tleap_selection = {}

        if protein_ff:
            # Extract display name from leaprc (e.g., "leaprc.protein.ff19SB" -> "ff19SB")
            name = protein_ff.split('.')[-1] if '.' in protein_ff else protein_ff
            tleap_selection["protein"] = {"name": name, "leaprc": protein_ff}

        if water_model:
            # Extract display name and determine box type
            name = water_model.split('.')[-1] if '.' in water_model else water_model

            # Map water leaprc to box type
            water_box_map = {
                'leaprc.water.opc': 'OPCBOX',
                'leaprc.water.opc3': 'OPC3BOX',
                'leaprc.water.opc3pol': 'OPC3BOX',
                'leaprc.water.tip4pew': 'TIP4PEWBOX',
                'leaprc.water.tip3p': 'TIP3PBOX',
                'leaprc.water.spce': 'SPCBOX',
                'leaprc.water.spceb': 'SPCBOX',
                'leaprc.water.tip5p': 'TIP5PBOX',
            }
            box_type = water_box_map.get(water_model, 'OPCBOX')

            tleap_selection["water"] = {
                "name": name,
                "leaprc": water_model,
                "box": box_type
            }

        # Store for tLEaP to find
        self.workspace.set("selected_standard_forcefields", tleap_selection)
        self.console.print(f"[grey50]Set tLEaP force field selection from Step 0b[/grey50]")

    # =========================================================================
    # Step 0c: Hydrogen Addition
    # =========================================================================

    def _step_0c_hydrogen_addition(
        self,
        pdb_file: str,
        triage: Dict[str, str],
        output_dir: Path,
        interactive: bool
    ) -> Dict[str, Path]:
        """
        Step 0c: Add hydrogens to structure components.

        Processing paths:
        - Protein (Category A): protonation_state_analyzer + tleap_input_generator
        - ALL non-standard (Category B): reduce program
        - Waters (Category C): keep as-is
        - Metals (Category D): keep as-is

        Note: We use reduce for ALL non-standard residues regardless of whether
        they'll eventually be in FF or need small_molecule_parameterizer.
        That determination happens in Step 0g.

        Returns:
            Dict mapping category/res_key → Path to PDB with H's
        """
        results = {}

        # Separate residues by category
        protein_residues = [k for k, v in triage.items() if v == 'A']
        nonstandard_residues = [k for k, v in triage.items() if v == 'B']

        # === Category A: Protein ===
        if protein_residues:
            self.console.print("\n[bold cyan]Processing Protein (Category A)[/bold cyan]")
            results['protein'] = self._process_protein_hydrogens(
                pdb_file, protein_residues, output_dir, interactive
            )

        # === Category B: ALL non-standard residues → reduce ===
        if nonstandard_residues:
            self.console.print("\n[bold cyan]Processing Non-Standard Residues (Category B)[/bold cyan]")
            self.console.print("[grey50]Using reduce to add hydrogens[/grey50]")
            for res_key in nonstandard_residues:
                chain, resid, resname = res_key.split(':')
                results[res_key] = self._process_reduce_hydrogens(
                    pdb_file, chain, int(resid), resname, output_dir
                )

        # Waters and metals are kept as-is (no H addition needed)

        return results

    def _process_protein_hydrogens(
        self,
        pdb_file: str,
        protein_residues: List[str],
        output_dir: Path,
        interactive: bool
    ) -> Path:
        """
        Process protein through protonation_state_analyzer + tLEaP.

        LAUNCHES existing modules interactively.
        """
        self.console.print("[grey50]Launching protonation state analyzer...[/grey50]")

        # Extract protein to temporary PDB
        protein_pdb = output_dir / "protein_only.pdb"
        self._extract_residues_to_pdb(pdb_file, protein_residues, protein_pdb)
        self.console.print(f"[grey50]Extracted {len(protein_residues)} protein residues to {protein_pdb.name}[/grey50]")

        # Store residue mapping BEFORE tLEaP (which may renumber)
        # Maps (chain, resid) → sequence_position for later sync
        residue_sequence_map = self._build_residue_sequence_map(protein_pdb)
        if self.workspace and residue_sequence_map:
            self.workspace.set("preprocessing_residue_sequence_map", residue_sequence_map)
            self.console.print(f"[grey50]Stored residue sequence mapping ({len(residue_sequence_map)} residues)[/grey50]")

        # Store protein PDB in workspace using preprocessing key
        # This has highest priority in StructureSelector so modules will use it
        if self.workspace:
            self.workspace.set("preprocessing_protein_input", str(protein_pdb))

        # === Launch protonation_state_analyzer ===
        if self.processor:
            prot_analyzer = self.processor.get_module_instance("Protonation State Analyzer")
            if prot_analyzer:
                self.console.print("\n[bold]═══ Protonation State Analyzer ═══[/bold]")
                self.console.print("[grey50]Determine protonation states for HIS, ASP, GLU, etc.[/grey50]\n")

                # Step 1: Analyze protonation states
                if hasattr(prot_analyzer, 'analyze_protonation_states'):
                    prot_analyzer.analyze_protonation_states()
                elif hasattr(prot_analyzer, 'process'):
                    prot_analyzer.process(self.workspace)
                else:
                    self.console.print("[yellow]Could not run protonation analysis[/yellow]")

                # Step 2: Set residue names (HIS→HIE/HID/HIP, etc.)
                if hasattr(prot_analyzer, 'set_residue_names'):
                    self.console.print("\n[grey50]Setting residue names based on protonation states...[/grey50]")
                    prot_analyzer.set_residue_names()
                else:
                    self.console.print("[yellow]Could not set residue names[/yellow]")
            else:
                self.console.print("[yellow]Protonation State Analyzer not available[/yellow]")

        # === Launch tLEaP input generator ===
        if self.processor:
            tleap_gen = self.processor.get_module_instance("Topology Generator")
            if tleap_gen:
                self.console.print("\n[bold]═══ Topology Generator ═══[/bold]")
                self.console.print("[grey50]Generate topology and coordinate files[/grey50]\n")

                # Set tLEaP-compatible FF selection from Step 0b choices
                self._set_tleap_ff_selection()

                # Keep provisional implicit solvation across the WHOLE run
                # (process + generate_single_state, where the solvation flow
                # actually runs). See _run_tleap_assembly for the rationale.
                self.workspace.set("_preprocessing_tleap_active", True)
                try:
                    # First do initial processing (bond definitions, parameters)
                    if hasattr(tleap_gen, 'process'):
                        tleap_gen.process(self.workspace)

                    # Then generate single-state tLEaP input and run it
                    if hasattr(tleap_gen, 'handle_menu_option'):
                        # Generate tLEaP input file
                        self.console.print("[cyan]Generating tLEaP input file...[/cyan]")
                        tleap_gen.handle_menu_option("generate_single_state")

                        # Generate topology files (runs tleap)
                        self.console.print("[cyan]Running tLEaP to generate topology...[/cyan]")
                        tleap_gen.handle_menu_option("generate_topology")
                    else:
                        self.console.print("[yellow]Could not run tLEaP generator[/yellow]")
                finally:
                    self.workspace.set("_preprocessing_tleap_active", False)
            else:
                self.console.print("[yellow]Topology Generator not available[/yellow]")

        # === Convert parm7/rst7 to PDB ===
        parm7 = self.workspace.get("parm7_file") if self.workspace else None
        rst7 = self.workspace.get("rst7_file") if self.workspace else None

        # Clear preprocessing key now that we're done with protein preprocessing
        if self.workspace:
            self.workspace.set("preprocessing_protein_input", None)

        if parm7 and rst7:
            protein_with_h = output_dir / "protein_with_H.pdb"
            self._convert_amber_to_pdb(parm7, rst7, protein_with_h)
            return protein_with_h
        else:
            self.console.print("[yellow]tLEaP output files not found, returning original[/yellow]")
            return protein_pdb

    def _process_reduce_hydrogens(
        self,
        pdb_file: str,
        chain: str,
        resid: int,
        resname: str,
        output_dir: Path
    ) -> Path:
        """
        Add hydrogens to a residue using reduce.

        Uses existing reduce functions from small_molecule_parameterizer.
        """
        # Extract residue to temporary PDB
        residue_pdb = output_dir / f"{resname}_{chain}_{resid}.pdb"
        self._extract_single_residue(pdb_file, chain, resid, residue_pdb)

        try:
            from proprep.forcefield_prep.small_molecule_parameterizer import (
                check_reduce_availability,
                configure_reduce_options_aligned,
                run_reduce_aligned
            )

            # Check reduce availability
            available, info = check_reduce_availability()
            if not available:
                self.console.print(f"[red]reduce not available: {info}[/red]")
                return residue_pdb

            # Configure and run reduce
            options = configure_reduce_options_aligned(
                interactive=True,
                console=self.console,
                processor=self.processor
            )

            output_pdb = output_dir / f"{resname}_{chain}_{resid}_H.pdb"
            success, msg = run_reduce_aligned(
                str(residue_pdb), str(output_pdb), options, self.console
            )

            if success:
                self.console.print(f"[green]Added H atoms to {resname}[/green]")
                return output_pdb
            else:
                self.console.print(f"[red]reduce failed for {resname}: {msg}[/red]")
                return residue_pdb

        except ImportError:
            self.console.print("[yellow]reduce functions not available[/yellow]")
            return residue_pdb

    # =========================================================================
    # Step 0d: Atom Exclusion + H Capping
    # =========================================================================

    def _step_0d_atom_exclusion(
        self,
        pdb_file: str,
        h_addition_results: Dict[str, Path],
        output_dir: Path,
        interactive: bool
    ) -> Tuple[List[Tuple], List[Dict]]:
        """
        Step 0d: Manual atom exclusion and H capping.

        User manually specifies:
        - Residues to exclude entirely
        - Specific atoms to exclude
        - Bonds to cut (specify kept atom and removed atom)

        Returns:
            Tuple of:
            - excluded_atoms: List of coordinate tuples to remove
            - h_caps: List of {kept_coords, h_coords, h_name, h_element}
        """
        excluded_atoms = []
        h_caps = []

        if not interactive:
            return excluded_atoms, h_caps

        self.console.print(Panel(
            "[cyan]Atom Exclusion[/cyan]\n"
            "Specify atoms/residues to exclude from the final structure.\n"
            "Common use: excluding propionate groups, peripheral chains, etc.",
            border_style="cyan",
            expand=False
        ))

        # Display structure summary
        self._display_structure_summary(pdb_file)

        # Option 1: Exclude entire residues
        exclude_residues = confirm_with_context(
            self.processor,
            "Exclude entire residues?",
            default=False,
            module="Structure Preprocessor",
            description="Exclude entire residues from structure",
        )
        if exclude_residues:
            residue_exclusions = self._prompt_residue_exclusions(pdb_file)
            excluded_atoms.extend(residue_exclusions)

        # Option 2: Cut bonds (specify atoms)
        cut_bonds = confirm_with_context(
            self.processor,
            "Cut specific bonds (and cap with H)?",
            default=False,
            module="Structure Preprocessor",
            description="Cut specific bonds and cap with H",
        )
        if cut_bonds:
            bond_cuts, new_caps = self._prompt_bond_cuts_manual(pdb_file, h_addition_results)
            excluded_atoms.extend(bond_cuts)
            h_caps.extend(new_caps)

        # Summary
        if excluded_atoms or h_caps:
            self.console.print(f"\n[green]Excluding {len(excluded_atoms)} atoms[/green]")
            self.console.print(f"[green]Adding {len(h_caps)} H caps[/green]")

        return excluded_atoms, h_caps

    def _display_structure_summary(self, pdb_file: str) -> None:
        """Display a summary of the structure for exclusion selection."""
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('struct', pdb_file)

        table = Table(title="Structure Summary")
        table.add_column("Chain")
        table.add_column("Residue ID")
        table.add_column("Name")
        table.add_column("Atoms", justify="right")

        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    resid = residue.get_id()[1]
                    n_atoms = len(list(residue.get_atoms()))

                    # Only show non-standard residues (most likely candidates for exclusion)
                    if resname not in STANDARD_RESIDUES and resname not in WATER_RESIDUES:
                        table.add_row(chain.id, str(resid), resname, str(n_atoms))

        self.console.print(table)

    def _prompt_residue_exclusions(self, pdb_file: str) -> List[Tuple]:
        """Prompt user to specify residues to exclude entirely."""
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('struct', pdb_file)

        excluded_coords = []

        self.console.print("\nEnter residues to exclude (format: CHAIN:RESID:RESNAME)")
        self.console.print("Example: A:501:PRN or A:502:HEM")
        self.console.print("Type 'done' when finished.\n")

        while True:
            res_spec = prompt_with_context(
                self.processor,
                "Residue to exclude (or 'done')",
                default="done",
                module="Structure Preprocessor",
                description="Residue to exclude (CHAIN:RESID:RESNAME or 'done')",
            )
            if res_spec.lower() == 'done':
                break

            try:
                chain_id, resid_str, resname = res_spec.split(':')
                resid = int(resid_str)

                # Find and add all atoms from this residue
                found = False
                for model in structure:
                    for chain in model:
                        if chain.id == chain_id:
                            for residue in chain:
                                if residue.get_id()[1] == resid:
                                    for atom in residue:
                                        coords = tuple(round(x, 3) for x in atom.coord)
                                        excluded_coords.append(coords)
                                    found = True
                                    self.console.print(f"[green]Marked {resname} for exclusion[/green]")

                if not found:
                    self.console.print(f"[red]Residue not found: {res_spec}[/red]")

            except ValueError:
                self.console.print("[red]Invalid format. Use CHAIN:RESID:RESNAME[/red]")

        return excluded_coords

    def _prompt_bond_cuts_manual(
        self,
        pdb_file: str,
        h_addition_results: Dict[str, Path]
    ) -> Tuple[List[Tuple], List[Dict]]:
        """Prompt user to manually specify bonds to cut."""
        import numpy as np

        excluded_coords = []
        h_caps = []

        self.console.print("\nSpecify bonds to cut by naming the KEPT atom and REMOVED atom.")
        self.console.print("Format: CHAIN:RESID:ATOMNAME for each atom")
        self.console.print("Type 'done' when finished.\n")

        # Build atom lookup from all structures
        atom_lookup = self._build_atom_lookup(pdb_file, h_addition_results)

        while True:
            kept_spec = prompt_with_context(
                self.processor,
                "KEPT atom (or 'done')",
                default="done",
                module="Structure Preprocessor",
                description="Bond-cut KEPT atom (CHAIN:RESID:ATOMNAME or 'done')",
            )
            if kept_spec.lower() == 'done':
                break

            removed_spec = prompt_with_context(
                self.processor,
                "REMOVED atom",
                module="Structure Preprocessor",
                description="Bond-cut REMOVED atom (CHAIN:RESID:ATOMNAME)",
            )

            kept_coords = atom_lookup.get(kept_spec)
            removed_coords = atom_lookup.get(removed_spec)

            if kept_coords and removed_coords:
                # Add removed atom to exclusion list
                excluded_coords.append(removed_coords)

                # Calculate H cap position
                h_cap = self._calculate_h_cap(kept_coords, removed_coords)
                h_caps.append(h_cap)

                self.console.print(f"[green]Bond cut: {kept_spec} --- {removed_spec}[/green]")
            else:
                if not kept_coords:
                    self.console.print(f"[red]KEPT atom not found: {kept_spec}[/red]")
                if not removed_coords:
                    self.console.print(f"[red]REMOVED atom not found: {removed_spec}[/red]")

        return excluded_coords, h_caps

    def _build_atom_lookup(
        self,
        pdb_file: str,
        h_addition_results: Dict[str, Path]
    ) -> Dict[str, Tuple[float, float, float]]:
        """Build lookup table from atom spec to coordinates."""
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)

        lookup = {}

        # Add atoms from original PDB
        structure = parser.get_structure('orig', pdb_file)
        for model in structure:
            for chain in model:
                for residue in chain:
                    resid = residue.get_id()[1]
                    for atom in residue:
                        spec = f"{chain.id}:{resid}:{atom.name.strip()}"
                        lookup[spec] = tuple(round(x, 3) for x in atom.coord)

        # Add atoms from hydrogen-added files
        for res_key, h_path in h_addition_results.items():
            if h_path.exists():
                try:
                    h_struct = parser.get_structure('h', str(h_path))
                    for model in h_struct:
                        for chain in model:
                            for residue in chain:
                                resid = residue.get_id()[1]
                                for atom in residue:
                                    spec = f"{chain.id}:{resid}:{atom.name.strip()}"
                                    lookup[spec] = tuple(round(x, 3) for x in atom.coord)
                except Exception:
                    pass

        return lookup

    def _calculate_h_cap(
        self,
        kept_coords: Tuple[float, float, float],
        removed_coords: Tuple[float, float, float]
    ) -> Dict:
        """Calculate H cap position along bond vector."""
        import numpy as np

        kept = np.array(kept_coords)
        removed = np.array(removed_coords)

        bond_vec = removed - kept
        bond_length = np.linalg.norm(bond_vec)
        unit_vec = bond_vec / bond_length

        H_BOND_LENGTH = 1.09  # Standard C-H bond length
        h_pos = kept + H_BOND_LENGTH * unit_vec

        return {
            'kept_coords': kept_coords,
            'h_coords': tuple(h_pos),
            'h_name': 'H',
            'h_element': 'H'
        }

    # =========================================================================
    # Step 0e: Structure Recombination
    # =========================================================================

    def _step_0e_recombination(
        self,
        original_pdb: str,
        h_addition_results: Dict[str, Path],
        excluded_atoms: List[Tuple],
        h_caps: List[Dict],
        triage: Dict[str, str],
        output_dir: Path
    ) -> Path:
        """
        Step 0e: Merge all pieces into combined structure.

        Process:
        1. Start with protein PDB (with H's)
        2. Add each ligand piece (with H's)
        3. Remove excluded atoms
        4. Add H cap atoms
        5. Add waters and ions from original structure
        6. Write combined PDB with sequential numbering

        Args:
            original_pdb: Path to original input PDB (for waters/ions)
            h_addition_results: Dict mapping category/res_key → Path to H-added PDB
            excluded_atoms: List of coordinate tuples to exclude
            h_caps: List of {kept_coords, h_coords, h_name, h_element} for H caps
            triage: Dict mapping residue_id → category (A/B/C/D/E)
            output_dir: Output directory
        """
        from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue, Atom
        import numpy as np

        # Create new structure
        combined = Structure.Structure('combined')
        model = Model.Model(0)
        combined.add(model)

        atom_serial = 1
        excluded_set = set(excluded_atoms)
        components_added = []

        # Track which residues we've added H caps to (by kept_coords)
        h_cap_residues = {}  # kept_coords → h_cap dict

        # 1. Add protein atoms (if we have processed protein)
        if 'protein' in h_addition_results:
            protein_path = h_addition_results['protein']
            if protein_path.exists():
                start_serial = atom_serial
                atom_serial = self._add_atoms_from_pdb(
                    model, protein_path, excluded_set, atom_serial
                )
                atoms_added = atom_serial - start_serial
                if atoms_added > 0:
                    components_added.append(f"Protein ({atoms_added} atoms)")

        # 2. Add ligand pieces
        ligand_count = 0
        ligand_atoms = 0
        for res_key, pdb_path in h_addition_results.items():
            if res_key == 'protein':
                continue
            if pdb_path.exists():
                start_serial = atom_serial
                atom_serial = self._add_atoms_from_pdb(
                    model, pdb_path, excluded_set, atom_serial
                )
                atoms_added = atom_serial - start_serial
                if atoms_added > 0:
                    ligand_count += 1
                    ligand_atoms += atoms_added

        if ligand_count > 0:
            components_added.append(f"{ligand_count} ligand(s) ({ligand_atoms} atoms)")

        # 3. Add H cap atoms at bond cut positions
        for h_cap in h_caps:
            kept_coords = h_cap.get('kept_coords')
            h_coords = h_cap.get('h_coords')
            h_name = h_cap.get('h_name', 'H')
            h_element = h_cap.get('h_element', 'H')

            if kept_coords and h_coords:
                # Find the residue containing the kept atom
                residue_found = False
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            atom_coords = tuple(round(x, 3) for x in atom.coord)
                            if atom_coords == kept_coords:
                                # Add H cap to this residue
                                new_h = Atom.Atom(
                                    h_name,
                                    np.array(h_coords),
                                    1.0,  # bfactor
                                    1.0,  # occupancy
                                    ' ',  # altloc
                                    f' {h_name} ',  # fullname
                                    atom_serial,
                                    h_element
                                )
                                residue.add(new_h)
                                atom_serial += 1
                                residue_found = True
                                self.console.print(
                                    f"[grey50]Added H cap to {residue.get_resname()} "
                                    f"at {tuple(round(x, 2) for x in h_coords)}[/grey50]"
                                )
                                break
                        if residue_found:
                            break
                    if residue_found:
                        break

        # 4. Add waters and ions from original structure
        # Waters are category 'E', metals are category 'D'
        water_ion_residues = {
            k for k, v in triage.items()
            if v in ('E', 'D')
        }

        if water_ion_residues and original_pdb:
            parser = PDBParser(QUIET=True)
            orig_structure = parser.get_structure('original', str(original_pdb))

            water_count = 0
            ion_count = 0

            for orig_model in orig_structure:
                for orig_chain in orig_model:
                    for orig_residue in orig_chain:
                        res_id = orig_residue.get_id()[1]
                        res_name = orig_residue.get_resname().strip()
                        res_key = f"{orig_chain.id}:{res_id}:{res_name}"

                        if res_key in water_ion_residues:
                            # Get or create chain
                            chain_id = orig_chain.id
                            if chain_id not in [c.id for c in model]:
                                model.add(Chain.Chain(chain_id))
                            target_chain = model[chain_id]

                            # Create new residue
                            new_res_id = orig_residue.id
                            if new_res_id not in [r.id for r in target_chain]:
                                new_res = Residue.Residue(new_res_id, res_name, '')
                                target_chain.add(new_res)
                                target_residue = new_res
                            else:
                                target_residue = target_chain[new_res_id]

                            # Add atoms
                            for orig_atom in orig_residue:
                                coords = tuple(round(x, 3) for x in orig_atom.coord)
                                if coords not in excluded_set:
                                    new_atom = Atom.Atom(
                                        orig_atom.name,
                                        orig_atom.coord,
                                        orig_atom.bfactor,
                                        orig_atom.occupancy,
                                        orig_atom.altloc,
                                        orig_atom.fullname,
                                        atom_serial,
                                        orig_atom.element
                                    )
                                    target_residue.add(new_atom)
                                    atom_serial += 1

                            # Track counts
                            if triage.get(res_key) == 'E':
                                water_count += 1
                            elif triage.get(res_key) == 'D':
                                ion_count += 1

            if water_count > 0:
                components_added.append(f"{water_count} water(s)")
            if ion_count > 0:
                components_added.append(f"{ion_count} ion(s)/metal(s)")

        # Display recombination summary
        if components_added:
            self.console.print(f"[cyan]Combined: {', '.join(components_added)}[/cyan]")
        if excluded_atoms:
            self.console.print(f"[grey50]Excluded {len(excluded_atoms)} atom(s)[/grey50]")
        if h_caps:
            self.console.print(f"[grey50]Added {len(h_caps)} H cap(s)[/grey50]")

        # Write combined PDB
        output_pdb = output_dir / "prepared_structure.pdb"
        io = PDBIO()
        io.set_structure(combined)
        io.save(str(output_pdb))

        self.console.print(f"[green]✓ Combined structure written to {output_pdb.name}[/green]")

        return output_pdb

    def _add_atoms_from_pdb(
        self,
        model,
        pdb_path: Path,
        excluded_set: set,
        start_serial: int
    ) -> int:
        """Add atoms from a PDB file to the model, skipping excluded atoms."""
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('temp', str(pdb_path))

        serial = start_serial

        for src_model in structure:
            for src_chain in src_model:
                # Get or create chain in target
                chain_id = src_chain.id
                if chain_id not in [c.id for c in model]:
                    from Bio.PDB import Chain
                    model.add(Chain.Chain(chain_id))
                target_chain = model[chain_id]

                for src_residue in src_chain:
                    # Get or create residue in target
                    res_id = src_residue.id
                    if res_id not in [r.id for r in target_chain]:
                        from Bio.PDB import Residue
                        new_res = Residue.Residue(res_id, src_residue.resname, '')
                        target_chain.add(new_res)
                    target_residue = target_chain[res_id]

                    for src_atom in src_residue:
                        coords = tuple(round(x, 3) for x in src_atom.coord)
                        if coords not in excluded_set:
                            # Clone atom with new serial
                            from Bio.PDB import Atom
                            new_atom = Atom.Atom(
                                src_atom.name,
                                src_atom.coord,
                                src_atom.bfactor,
                                src_atom.occupancy,
                                src_atom.altloc,
                                src_atom.fullname,
                                serial,
                                src_atom.element
                            )
                            target_residue.add(new_atom)
                            serial += 1

        return serial

    # =========================================================================
    # Step 0f: Launch Redox Site Detector
    # =========================================================================

    def _step_0f_redox_site_sync(
        self,
        final_pdb: Path,
        interactive: bool
    ) -> List['RedoxSite']:
        """
        Step 0f: Synchronize RedoxSite objects with prepared structure.

        Primary path: Update existing RedoxSite objects to match the prepared
        structure. For each residue in the original site:
        - Include all atoms from prepared structure (including new H atoms)
        - Remove atoms that were excluded in Step 0d
        - Update bonds to only include those where both atoms exist

        Fallback: If no existing sites, launch interactive Redox Site Detector.
        """
        self.console.print(Panel(
            "[cyan]Redox Site Synchronization[/cyan]\n"
            "Updating redox sites to match prepared structure",
            border_style="cyan",
            expand=False
        ))

        # Store final PDB in workspace
        if self.workspace:
            self.workspace.set("structure_pdb_file", str(final_pdb))

        # Check for existing redox sites.
        #
        # Sync is destructive: it mutates site objects in place (atoms, centers,
        # coords) and writes them back to detected_redox_sites. To make re-runs
        # non-destructive, snapshot a PRISTINE copy the first time and always
        # sync FROM that pristine copy (on a fresh deep copy so the pristine
        # snapshot itself is never mutated). This way a re-run — or a partial
        # failure — can never corrupt the original site definitions.
        import copy as _copy
        raw_sites = self.workspace.get("detected_redox_sites", []) if self.workspace else []
        if self.workspace:
            pristine = self.workspace.get("redox_sites_pristine")
            if pristine is None and raw_sites:
                self.workspace.set("redox_sites_pristine", _copy.deepcopy(raw_sites))
            elif pristine is not None:
                raw_sites = pristine
        existing_sites = _ensure_redox_site_objects(_copy.deepcopy(raw_sites))

        # Validate that existing sites have actual atoms/centers
        # If sites are empty (e.g., from corrupted saved state), treat as no sites
        if existing_sites:
            valid_sites = []
            for site in existing_sites:
                has_atoms = hasattr(site, 'atoms') and len(site.atoms) > 0
                has_centers = hasattr(site, 'centers') and len(site.centers) > 0
                if has_atoms or has_centers:
                    valid_sites.append(site)
                else:
                    self.console.print(f"[yellow]⚠ Site {getattr(site, 'site_id', 'unknown')} has no atoms/centers, skipping[/yellow]")

            if not valid_sites:
                self.console.print("[yellow]⚠ Existing redox sites are empty (possibly from corrupted saved state)[/yellow]")
                existing_sites = []
            else:
                existing_sites = valid_sites

        if existing_sites:
            # Primary path: sync existing sites with prepared structure
            self.console.print(f"[cyan]Found {len(existing_sites)} existing redox site(s) to synchronize[/cyan]")
            synced_sites = self._sync_redox_sites_with_structure(existing_sites, final_pdb)

            if synced_sites:
                # Store synced sites
                if self.workspace:
                    self.workspace.set("redox_sites", synced_sites)
                    self.workspace.set("detected_redox_sites", synced_sites)

                # Convert bonds to tLEaP commands now (coord_to_pdb is valid)
                # and store in a dedicated workspace key for the tLEaP generator.
                # This mirrors how preprocessing_lib_files / preprocessing_frcmod_files work.
                bond_commands = self._convert_site_bonds_to_tleap_commands(synced_sites)
                if bond_commands and self.workspace:
                    self.workspace.set("preprocessing_bond_commands", bond_commands)

                self.console.print(f"\n[green]✓ Synchronized {len(synced_sites)} redox site(s)[/green]")
                for site in synced_sites:
                    n_atoms = len(site.atoms) if hasattr(site, 'atoms') else 0
                    n_bonds = len(site.bonds) if hasattr(site, 'bonds') else 0
                    site_id = site.site_id if hasattr(site, 'site_id') else 'unknown'
                    self.console.print(f"  - {site_id}: {n_atoms} atoms, {n_bonds} bonds")

                return synced_sites

        # Fallback: no existing sites, launch interactive detection
        self.console.print("[yellow]No existing redox sites found in workspace[/yellow]")
        self.console.print("[grey50]Tip: Run Redox Site Detector before Force Field Parameterizer[/grey50]")

        if interactive and self.processor:
            run_detection = confirm_with_context(
                self.processor,
                "Run interactive redox site detection now?",
                default=True,
                module="Structure Preprocessor",
                description="Run interactive redox site detection",
            )
            if run_detection:
                return self._run_interactive_redox_detection(final_pdb)

        return []

    def _sync_redox_sites_with_structure(
        self,
        redox_sites: List['RedoxSite'],
        final_pdb: Path
    ) -> List['RedoxSite']:
        """
        Synchronize RedoxSite objects with the prepared structure.

        Uses sequence position mapping to handle tLEaP residue renumbering:
        1. Look up original residue's sequence position (from pre-tLEaP mapping)
        2. Find residue at that position in prepared structure (may have new resid)
        3. Include all atoms from that residue (including new H atoms)
        4. Filter bonds to only those where both atoms still exist
        5. Rebuild residue_groups
        """
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("prepared", str(final_pdb))

        # Get the pre-tLEaP sequence mapping from workspace
        # Maps (chain, resid) → global_sequence_position
        original_to_seq_pos = self.workspace.get("preprocessing_residue_sequence_map", {}) if self.workspace else {}

        # Build position-to-residue map for prepared structure
        # Maps global_sequence_position → (chain, resid)
        # Note: tLEaP removes chain IDs and numbers residues globally
        position_to_new_resid = self._build_position_to_residue_map(final_pdb)

        self.console.print(f"[grey50]Sequence map: {len(original_to_seq_pos)} residues, Position map: {len(position_to_new_resid)} residues[/grey50]")

        # Build lookup: (chain, resid) -> list of atoms in prepared structure
        residue_atoms = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    resid = residue.get_id()[1]
                    resname = residue.get_resname().strip()
                    key = (chain.id, resid)
                    if key not in residue_atoms:
                        residue_atoms[key] = []
                    for atom in residue:
                        residue_atoms[key].append({
                            'chain': chain.id,
                            'resid': resid,
                            'resname': resname,
                            'atom_name': atom.name.strip(),
                            'element': atom.element.strip() if atom.element else self._guess_element(atom.name),
                            'coords': tuple(round(x, 3) for x in atom.coord),
                            'altloc': atom.get_altloc() or '',
                            'occupancy': atom.get_occupancy(),
                            'bfactor': atom.get_bfactor(),
                        })

        synced_sites = []

        for site in redox_sites:
            # Get residue keys from original site
            original_residue_keys = set()
            for atom in site.atoms:
                original_residue_keys.add((atom.chain, atom.resid))

            # Also include residues from centers
            if hasattr(site, 'centers'):
                for center in site.centers:
                    original_residue_keys.add((center.chain, center.resid))

            self.console.print(f"\n[grey50]Syncing {site.site_id}: {len(original_residue_keys)} residues[/grey50]")

            # Map original residue keys to new keys using sequence position
            # Note: We use actual tLEaP chain IDs (not restored original) so downstream
            # lookups (terminal classifier, model builder, etc.) work correctly
            new_residue_keys = set()
            # Remember which prepared-structure residue each original residue
            # mapped to, so centers can be updated even when the residue was
            # renamed during preparation (e.g. 9E2 -> x9e).
            orig_to_new_key = {}

            for orig_key in original_residue_keys:
                # First try direct match (for non-protein residues like metals)
                if orig_key in residue_atoms:
                    new_residue_keys.add(orig_key)
                    orig_to_new_key[orig_key] = orig_key
                    continue

                # Check metal reinsertion map (for metals that were removed and reinserted)
                metal_reinsertion_map = self.workspace.get("preprocessing_metal_reinsertion_map", {}) if self.workspace else {}
                if orig_key in metal_reinsertion_map:
                    new_chain, new_resid = metal_reinsertion_map[orig_key]
                    new_key = (new_chain, new_resid)
                    if new_key in residue_atoms:
                        new_residue_keys.add(new_key)
                        orig_to_new_key[orig_key] = new_key
                        self.console.print(f"  [grey50]{orig_key[0]}:{orig_key[1]} → {new_key[0] or '(no chain)'}:{new_key[1]} (reinserted metal)[/grey50]")
                        continue

                # Use sequence position mapping for protein residues
                seq_pos = original_to_seq_pos.get(orig_key)
                if seq_pos is not None:
                    # Find new resid at this global sequence position
                    new_key = position_to_new_resid.get(seq_pos)

                    if new_key and new_key in residue_atoms:
                        new_residue_keys.add(new_key)
                        orig_to_new_key[orig_key] = new_key
                        self.console.print(f"  [grey50]{orig_key[0]}:{orig_key[1]} → {new_key[0] or '(no chain)'}:{new_key[1]} (pos {seq_pos})[/grey50]")
                    else:
                        self.console.print(f"  [yellow]⚠ Residue {orig_key[0]}:{orig_key[1]} (pos {seq_pos}) not found in prepared structure[/yellow]")
                else:
                    self.console.print(f"  [yellow]⚠ Residue {orig_key[0]}:{orig_key[1]} not in sequence map (may be non-protein)[/yellow]")

            # Rebuild atoms list from prepared structure
            new_atoms = []
            new_coord_to_pdb = {}

            # Sort residue keys by resid for deterministic ordering
            for res_key in sorted(new_residue_keys, key=lambda k: (k[1] if isinstance(k[1], int) else int(k[1]) if str(k[1]).isdigit() else 0)):
                if res_key in residue_atoms:
                    # Use actual chain ID from tLEaP output (don't restore original)
                    # This ensures downstream lookups (terminal classifier, etc.) work correctly
                    actual_chain = res_key[0]

                    for atom_info in residue_atoms[res_key]:
                        from proprep.structure_prep.comprehensive_redox_detector import RedoxSiteAtom

                        new_atom = RedoxSiteAtom(
                            chain=actual_chain,  # Use actual tLEaP chain ID
                            resname=atom_info['resname'],
                            resid=atom_info['resid'],
                            atom_name=atom_info['atom_name'],
                            coords=atom_info['coords'],
                            element=atom_info['element'],
                            altloc=atom_info['altloc'],
                            occupancy=atom_info['occupancy'],
                            bfactor=atom_info['bfactor'],
                        )
                        new_atoms.append(new_atom)

                        new_coord_to_pdb[atom_info['coords']] = {
                            'chain': actual_chain,  # Use actual tLEaP chain ID
                            'resname': atom_info['resname'],
                            'resid': atom_info['resid'],
                            'atom_name': atom_info['atom_name'],
                            'element': atom_info['element'],
                            'altloc': atom_info['altloc'],
                        }

            # Track changes
            old_atom_count = len(site.atoms)
            new_atom_count = len(new_atoms)

            # Update site atoms
            site.atoms = new_atoms
            site.coord_to_pdb = new_coord_to_pdb

            # Clear bonds - coordinates have changed (tLEaP translates/rotates),
            # so old coordinate-based bonds are no longer valid
            old_bond_count = len(site.bonds) if hasattr(site, 'bonds') else 0
            site.bonds = []
            new_bond_count = 0

            # Update centers to match the prepared structure.
            # Match by residue KEY (via orig_to_new_key), not by the center's
            # resname: a residue may have been renamed during preparation
            # (e.g. 9E2 -> x9e), which would otherwise leave the center stale
            # and its coords un-refreshed.
            # NOTE on the local names below: `residue_atom_list`/`center_atoms`
            # are deliberately distinct from `residue_atoms` (the shared
            # (chain, resid) -> atoms lookup dict used for EVERY site) — reusing
            # that name previously clobbered it into a list and broke all later
            # sites' lookups.
            if hasattr(site, 'centers'):
                import numpy as np
                for center in site.centers:
                    found = False

                    # Restrict to the atoms of the residue this center maps to.
                    new_key = orig_to_new_key.get((center.chain, center.resid))
                    if new_key is not None:
                        center_atoms = [a for a in new_atoms if (a.chain, a.resid) == new_key]
                    else:
                        # Fallback: match by resname (unrenamed residues)
                        center_atoms = [a for a in new_atoms if a.resname == center.resname]

                    if center.atom_name is not None:
                        # Specific-atom center - match by atom name within the residue
                        for atom in center_atoms:
                            if atom.atom_name == center.atom_name:
                                center.coords = atom.coords
                                center.chain = atom.chain
                                center.resid = atom.resid
                                center.resname = atom.resname  # adopt any renamed residue
                                site.coord_to_pdb[center.coords] = {
                                    'chain': center.chain,
                                    'resname': center.resname,
                                    'resid': center.resid,
                                    'atom_name': center.atom_name,
                                    'element': center.element,
                                    'center_type': center.center_type.value if hasattr(center.center_type, 'value') else str(center.center_type)
                                }
                                found = True
                                break
                    else:
                        # Whole-residue center (atom_name is None) - centroid of the residue
                        if center_atoms:
                            coords_array = np.array([a.coords for a in center_atoms])
                            centroid = tuple(np.mean(coords_array, axis=0))
                            center.coords = centroid
                            center.chain = center_atoms[0].chain
                            center.resid = center_atoms[0].resid
                            center.resname = center_atoms[0].resname  # adopt any renamed residue
                            site.coord_to_pdb[center.coords] = {
                                'chain': center.chain,
                                'resname': center.resname,
                                'resid': center.resid,
                                'atom_name': None,  # Whole residue center
                                'element': None,
                                'center_type': center.center_type.value if hasattr(center.center_type, 'value') else str(center.center_type)
                            }
                            found = True

                    if not found:
                        self.console.print(f"  [yellow]⚠ Center {center.resname}:{center.atom_name or '(whole residue)'} not found in synced atoms[/yellow]")

            # Rebuild residue_groups
            site.residue_groups = {}
            for atom in site.atoms:
                res_key = (atom.chain, atom.resid, getattr(atom, 'insertion_code', ''))
                if res_key not in site.residue_groups:
                    site.residue_groups[res_key] = []
                site.residue_groups[res_key].append(atom.coords)

            # Report atom changes
            atom_diff = new_atom_count - old_atom_count
            atom_change = f"+{atom_diff}" if atom_diff > 0 else str(atom_diff)
            self.console.print(f"  [grey50]Atoms: {old_atom_count} → {new_atom_count} ({atom_change})[/grey50]")

            # Prompt user to redefine bonds (coordinates changed, old bonds invalid)
            if len(site.residue_groups) > 1:
                if old_bond_count > 0:
                    self.console.print(f"  [yellow]Bonds cleared (coordinates changed by tLEaP) - please redefine[/yellow]")
                self._define_site_bonds_interactively(site, structure_file=final_pdb)
            else:
                self.console.print(f"  [grey50]Single residue - no cross-residue bonds needed[/grey50]")

            new_bond_count = len(site.bonds)
            if new_bond_count > 0:
                self.console.print(f"  [green]Defined {new_bond_count} bond(s)[/green]")

            synced_sites.append(site)

        return synced_sites

    def _run_interactive_redox_detection(self, final_pdb: Path) -> List['RedoxSite']:
        """Fallback: run interactive redox site detection."""
        if not self.processor:
            return []

        detector = self.processor.get_module_instance("Redox Site Detector")
        if not detector:
            self.console.print("[red]Redox Site Detector not available[/red]")
            return []

        self.console.print("\n[bold]═══ Redox Site Detector ═══[/bold]")
        self.console.print("[grey50]Detecting metal centers and coordination spheres[/grey50]\n")

        # Use handle_menu_option to run interactive detection
        if hasattr(detector, 'handle_menu_option'):
            detector.handle_menu_option("detect")

        # Retrieve detected sites from workspace
        redox_sites = _ensure_redox_site_objects(self.workspace.get("detected_redox_sites", [])) if self.workspace else []

        if redox_sites:
            self.console.print(f"\n[green]Detected {len(redox_sites)} redox site(s)[/green]")
            for site in redox_sites:
                n_atoms = len(site.atoms) if hasattr(site, 'atoms') else 0
                n_bonds = len(site.bonds) if hasattr(site, 'bonds') else 0
                site_id = site.site_id if hasattr(site, 'site_id') else 'unknown'
                self.console.print(f"  - {site_id}: {n_atoms} atoms, {n_bonds} bonds")
        else:
            self.console.print("[yellow]No redox sites detected[/yellow]")

        return redox_sites

    # =========================================================================
    # Interactive Bond Definition (for Redox Site Sync)
    # =========================================================================

    def _launch_bond_definition_viewer(self, residue_list, structure_file) -> bool:
        """Open the 3D viewer at the bond-definition prompt, colouring each
        numbered site residue to match its table row.

        This is the "who bonds to whom" aid: seeing the metals and their
        coordination sphere in 3D, with table row [N] the same palette colour
        as residue [N], makes the bond pairs obvious.

        The structure tLEaP produced has blank chain IDs and globally
        renumbered residues (see ``_sync_redox_sites_with_structure``), so
        residues are selected by bare residue number, chain-qualified only when
        a chain ID actually survived. User-initiated view => ``force=True`` so
        it launches even in CLI mode. Any viewer failure is swallowed — the
        viewer is a convenience and must never block bond definition.
        """
        try:
            if not structure_file:
                structure_file = str(self._final_pdb) if getattr(self, "_final_pdb", None) else None
            if not structure_file or not Path(structure_file).exists():
                self.console.print("[grey50]No prepared structure available to view[/grey50]")
                return False

            from proprep.structure_prep.viewer_coordinator import viewer as _viewer

            selections = []
            for res_key in residue_list:
                chain, resname, resid, icode = res_key
                if chain and str(chain).strip():
                    selections.append(f"(:{chain} and {resid})")
                else:
                    selections.append(f"{resid}")

            _viewer.show_structure(str(structure_file), force=True)
            _viewer.clear_annotations()
            for i, (res_key, sel) in enumerate(zip(residue_list, selections), 1):
                chain, resname, resid, icode = res_key
                _viewer.highlight(
                    sel,
                    color=f"palette:{i}",
                    style="ball+stick",
                    label=f"[{i}] {resname}{resid}",
                    focused=False,
                    force=True,
                )
            # Zoom to the whole site so the coordination geometry is visible.
            _viewer.focus_on(" or ".join(selections))
            self.console.print(
                "[grey50]3D viewer: each numbered residue is coloured to match the table above.[/grey50]"
            )
            return True
        except Exception as e:  # noqa: BLE001 — viewer must never break the prompt
            self.console.print(f"[grey50]3D viewer unavailable ({e})[/grey50]")
            return False

    def _define_site_bonds_interactively(self, site: 'RedoxSite', structure_file=None) -> None:
        """
        Interactive bond definition for a synced RedoxSite.

        Follows the same UX pattern as comprehensive_redox_detector:
        1. Display numbered table of residues
        2. User enters residue pairs (e.g., "1-2, 1-3, 1-4")
        3. Group by source: source=1, targets=[2,3,4]
        4. Show ALL target atoms from ALL target residues at once
        5. User selects source atom, then target atoms from combined list
        6. Create bonds using site.add_bond_with_classification()

        Supports both inter-residue (1-2) and intra-residue (1-1) bonds.
        """
        import numpy as np
        from rich.table import Table

        # Group atoms by residue
        residue_groups = {}
        for atom in site.atoms:
            res_key = (atom.chain, atom.resname, atom.resid, getattr(atom, 'insertion_code', ''))
            if res_key not in residue_groups:
                residue_groups[res_key] = []
            residue_groups[res_key].append(atom)

        residue_list = list(residue_groups.keys())

        if len(residue_list) < 1:
            self.console.print("[grey50]No residues to define bonds[/grey50]")
            return

        # Display residue table
        table = Table(title="Site Residues")
        table.add_column("#", style="cyan", width=6)
        table.add_column("Residue", style="green")
        table.add_column("Atoms", style="yellow")

        for i, res_key in enumerate(residue_list, 1):
            chain, resname, resid, icode = res_key
            atoms = residue_groups[res_key]
            atom_elements = sorted(set(a.element for a in atoms))
            table.add_row(
                f"[{i}]",
                f"{resname} {chain}:{resid}",
                f"{len(atoms)} atoms ({', '.join(atom_elements)})"
            )

        self.console.print(table)

        # Auto-launch the 3D viewer so the user can SEE which residues are close
        # enough to bond, colour-matched to the table rows above.
        self._launch_bond_definition_viewer(residue_list, structure_file)

        # Combined prompt: enter bond pairs directly, 'v' to re-view, or 'none'.
        # Example is built from THIS site's residue IDs so it matches what the
        # user reads off the table / 3D viewer (no index conversion needed).
        example = ""
        if len(residue_list) >= 2:
            a, b = residue_list[0][2], residue_list[1][2]
            example = f" (e.g. {a}-{b}"
            if len(residue_list) >= 3:
                example += f" {a}-{residue_list[2][2]}"
            example += ")"
        self.console.print("\n[bold]Define bonds between residues[/bold]")
        self.console.print(
            f"[grey50]Enter residue-ID pairs — the numbers shown in the table and 3D viewer{example}. "
            f"Type 'v' to re-open the viewer, or 'none' to skip.[/grey50]"
        )
        while True:
            pairs_input = prompt_with_context(
                self.processor,
                "Bond pairs",
                default="none",
                module="Structure Preprocessor - Bond Definition",
                description="Residue pairs for bond definitions ('v' to view, 'none' to skip)",
                options_map={
                    "none": "Skip bond definition",
                    "v": "Re-open the 3D viewer",
                    "view": "Re-open the 3D viewer",
                },
            ).strip()

            if pairs_input.lower() in ("v", "view"):
                self._launch_bond_definition_viewer(residue_list, structure_file)
                continue
            pairs_input = pairs_input.lower()
            break

        if not pairs_input or pairs_input == "none":
            self.console.print("[grey50]No bonds defined[/grey50]")
            return

        # Parse pairs and group by source residue
        source_to_targets = self._parse_bond_residue_pairs(pairs_input, residue_list)

        if not source_to_targets:
            self.console.print("[yellow]No valid pairs entered[/yellow]")
            return

        # Process each source residue with all its targets
        for source_idx, target_indices in source_to_targets.items():
            self._define_bonds_from_source_residue(
                site, residue_groups, residue_list, source_idx, target_indices
            )

        self.console.print(f"\n[green]Defined {len(site.bonds)} bond(s)[/green]")

    def _parse_bond_residue_pairs(self, input_str: str, residue_list: list) -> Dict[int, List[int]]:
        """
        Parse bond pairs like '185-45 185-182' and group by source residue.

        Each token is interpreted as a RESIDUE ID first (the numbers shown in
        the Site Residues table and the 3D viewer), falling back to a 1-based
        table index only when the token is not a residue ID in this site. This
        lets the user type the residue IDs they read off the viewer instead of
        converting them to row numbers. Post-tLEaP residues are globally
        renumbered and unique, so residue IDs resolve unambiguously; the index
        fallback preserves the old '1-2' behavior.

        Returns:
            Dict mapping source_idx -> list of target_idx (all 0-indexed)
        """
        max_residue = len(residue_list)

        # Residue ID -> 0-based table indices (a single index in the normal,
        # globally-renumbered case; a list guards against a shared resid).
        resid_to_indices: Dict[int, List[int]] = {}
        for idx, res_key in enumerate(residue_list):
            resid_to_indices.setdefault(res_key[2], []).append(idx)

        def _token_to_index(tok: int) -> Optional[int]:
            hits = resid_to_indices.get(tok)
            if hits is not None:
                if len(hits) == 1:
                    return hits[0]
                # Shared residue ID: fall back to the table index if valid.
                if 1 <= tok <= max_residue:
                    return tok - 1
                self.console.print(
                    f"[yellow]Residue ID {tok} is ambiguous (appears {len(hits)}x); "
                    f"use the table index [1-{max_residue}] instead[/yellow]"
                )
                return None
            # Not a residue ID in this site -> treat as a 1-based table index.
            if 1 <= tok <= max_residue:
                return tok - 1
            return None

        source_to_targets: Dict[int, List[int]] = {}

        # Split by comma or space
        pair_strings = [p.strip() for p in input_str.replace(',', ' ').split() if p.strip()]

        for pair_str in pair_strings:
            if '-' not in pair_str:
                self.console.print(f"[yellow]Invalid format: {pair_str} (use e.g. 185-45)[/yellow]")
                continue

            parts = pair_str.split('-')
            if len(parts) != 2:
                continue

            try:
                first = int(parts[0].strip())
                second = int(parts[1].strip())
            except ValueError:
                self.console.print(f"[yellow]Invalid numbers in: {pair_str}[/yellow]")
                continue

            source_idx = _token_to_index(first)
            target_idx = _token_to_index(second)
            if source_idx is None or target_idx is None:
                bad = first if source_idx is None else second
                self.console.print(
                    f"[yellow]{bad} is not a residue ID or table index (1-{max_residue}) "
                    f"in this site[/yellow]"
                )
                continue

            # Group by source
            if source_idx not in source_to_targets:
                source_to_targets[source_idx] = []
            if target_idx not in source_to_targets[source_idx]:
                source_to_targets[source_idx].append(target_idx)

        return source_to_targets

    def _define_bonds_from_source_residue(
        self,
        site: 'RedoxSite',
        residue_groups: dict,
        residue_list: list,
        source_idx: int,
        target_indices: List[int]
    ) -> None:
        """
        Define bonds from a source residue to all its target residues.

        Shows ALL target atoms from ALL target residues in one table,
        sorted by distance from the selected source atom.
        """
        import numpy as np
        from rich.table import Table

        source_key = residue_list[source_idx]
        source_atoms = residue_groups[source_key]
        chain_src, resname_src, resid_src, _ = source_key
        source_label = f"{resname_src} {chain_src}:{resid_src}"

        # Collect ALL target atoms from ALL target residues
        target_atoms_info = []
        target_labels = []
        for target_idx in target_indices:
            target_key = residue_list[target_idx]
            target_residue_atoms = residue_groups[target_key]
            chain_tgt, resname_tgt, resid_tgt, _ = target_key
            target_label = f"{resname_tgt} {chain_tgt}:{resid_tgt}"
            target_labels.append(target_label)

            for atom in target_residue_atoms:
                target_atoms_info.append({
                    'atom': atom,
                    'residue_label': target_label,
                    'target_idx': target_idx
                })

        self.console.print(f"\n[bold underline]DEFINING BONDS FROM {source_label}[/bold underline]")
        self.console.print(f"Target residues: [cyan]{', '.join(target_labels)}[/cyan]")

        # Bond definition loop
        while True:
            # Show source atoms with bond status
            source_table = Table(title=f"Atoms in {source_label}")
            source_table.add_column("#", style="cyan", width=6)
            source_table.add_column("Atom", style="green", width=12)
            source_table.add_column("Status", style="yellow")

            for i, atom in enumerate(source_atoms, 1):
                # Check if this atom already has bonds
                existing_bonds = [b for b in site.bonds if
                                  (b.atom1_coords == atom.coords or b.atom2_coords == atom.coords)]
                if existing_bonds:
                    bond_partners = []
                    for bond in existing_bonds:
                        # Determine which atom is the partner
                        if bond.atom1_coords == atom.coords:
                            partner_info = bond.atom2_residue_info
                        else:
                            partner_info = bond.atom1_residue_info
                        if partner_info:
                            bond_partners.append(
                                f"{partner_info.get('chain', '?')}:{partner_info.get('resname', '?')}"
                                f"{partner_info.get('resid', '?')}:{partner_info.get('atom_name', '?')}"
                            )
                    status = f"bonded to {', '.join(bond_partners)}" if bond_partners else "bonded"
                else:
                    status = "available"

                source_table.add_row(f"[{i}]", atom.atom_name, status)

            self.console.print(source_table)

            # Select source atom
            choice = prompt_with_context(
                self.processor,
                f"Select atom from {source_label} (1-{len(source_atoms)}, atom name, or 'done')",
                default="done",
                module="Structure Preprocessor - Bond Definition"
            ).strip()

            if choice.lower() == 'done':
                break

            # Parse source selection
            source_atom = None
            try:
                idx = int(choice) - 1
                if 0 <= idx < len(source_atoms):
                    source_atom = source_atoms[idx]
            except ValueError:
                # Try matching by atom name
                matches = [a for a in source_atoms if a.atom_name.upper() == choice.upper()]
                if len(matches) == 1:
                    source_atom = matches[0]
                elif len(matches) > 1:
                    self.console.print(f"[yellow]Multiple atoms named '{choice}' - use number[/yellow]")
                    continue

            if source_atom is None:
                self.console.print("[yellow]Invalid selection[/yellow]")
                continue

            # Calculate distances to ALL target atoms
            for info in target_atoms_info:
                target_atom = info['atom']
                dist = np.linalg.norm(np.array(source_atom.coords) - np.array(target_atom.coords))
                info['distance'] = dist

            # Group by residue, sort within each residue by distance
            residue_atom_groups = {}
            for info in target_atoms_info:
                res_label = info['residue_label']
                if res_label not in residue_atom_groups:
                    residue_atom_groups[res_label] = []
                residue_atom_groups[res_label].append(info)

            # Sort atoms within each residue by distance
            for res_label in residue_atom_groups:
                residue_atom_groups[res_label].sort(key=lambda x: x['distance'])

            # Sort residues by their closest atom's distance
            residue_order = sorted(
                residue_atom_groups.keys(),
                key=lambda r: residue_atom_groups[r][0]['distance']
            )

            # Build display list maintaining residue grouping
            display_ordered_atoms = []
            for res_label in residue_order:
                for info in residue_atom_groups[res_label]:
                    display_ordered_atoms.append(info)

            # Display ALL target atoms in one table, grouped by residue
            target_table = Table(title=f"Target atoms (distance from {source_atom.atom_name} {source_atom.element})")
            target_table.add_column("#", style="cyan", width=6)
            target_table.add_column("Residue", style="green", width=12)
            target_table.add_column("Atom", style="yellow", width=8)
            target_table.add_column("Element", style="magenta", width=8)
            target_table.add_column("Distance", style="blue", width=10)

            i = 1
            for residue_idx, res_label in enumerate(residue_order):
                # Add section separator between residues (except first)
                if residue_idx > 0:
                    target_table.add_section()

                for info in residue_atom_groups[res_label]:
                    target_table.add_row(
                        f"[{i}]",
                        res_label,
                        info['atom'].atom_name,
                        info['atom'].element,
                        f"{info['distance']:.2f}Å"
                    )
                    i += 1

            self.console.print(target_table)

            # Use display_ordered_atoms for selection lookup
            targets_with_dist = display_ordered_atoms

            # Select target atoms (can select multiple)
            target_choice = prompt_with_context(
                self.processor,
                f"Select target atom(s) (1-{len(targets_with_dist)}, comma-separated, or 'cancel')",
                default="1",
                module="Structure Preprocessor - Bond Definition"
            ).strip()

            if target_choice.lower() == 'cancel':
                continue

            # Parse target selections
            selected_indices = []
            for part in target_choice.replace(',', ' ').split():
                part = part.strip()
                if not part:
                    continue
                try:
                    idx = int(part) - 1
                    if 0 <= idx < len(targets_with_dist):
                        selected_indices.append(idx)
                    else:
                        self.console.print(f"[yellow]Invalid: {part}[/yellow]")
                except ValueError:
                    # Try atom name
                    matches = [i for i, info in enumerate(targets_with_dist)
                               if info['atom'].atom_name.upper() == part.upper()]
                    if matches:
                        selected_indices.extend(matches)
                    else:
                        self.console.print(f"[yellow]Unknown atom: {part}[/yellow]")

            # Create bonds
            from proprep.structure_prep.comprehensive_redox_detector import METALS
            for sel_idx in selected_indices:
                info = targets_with_dist[sel_idx]
                target_atom = info['atom']
                distance = info['distance']

                # For a metal-water contact, offer to model the water as a
                # RESTRAINED nonbonded ligand (kept TIP3P, held near the metal by
                # an MD distance restraint) instead of a bonded MCPB ligand.
                # Default = bonded (unchanged behavior). Restraining avoids the
                # 1-4 electrostatic distortion of two bonded waters on one metal.
                treatment = "bonded"
                src_metal = (source_atom.element or '').upper() in METALS
                tgt_metal = (target_atom.element or '').upper() in METALS
                src_water = (source_atom.resname or '').upper() in ('HOH', 'WAT')
                tgt_water = (target_atom.resname or '').upper() in ('HOH', 'WAT')
                if (src_metal and tgt_water) or (tgt_metal and src_water):
                    # Name the actual water model chosen in param-2 (leaprc.water.<model>)
                    _wm = (self.workspace.get("preprocessing_water_model", "")
                           if self.workspace else "") or ""
                    _wm_name = _wm.split('.')[-1].upper() if _wm else "the selected water model"
                    if confirm_with_context(
                        self.processor,
                        f"Metal-water contact: model the water as RESTRAINED "
                        f"(keep it as {_wm_name}, held near the metal by an MD distance "
                        f"restraint) instead of a bonded MCPB ligand?",
                        default=False,
                        module="Structure Preprocessor - Bond Definition",
                        description=(
                            "Restrained water stays in the QM models for correct "
                            "electronics but emits no bonded term; it avoids the "
                            "1-4 electrostatic distortion of two bonded waters on "
                            "one metal."
                        ),
                    ):
                        treatment = "restrained"

                site.add_bond_with_classification(
                    source_atom.coords,
                    target_atom.coords,
                    distance,
                    treatment=treatment
                )
                mode_note = "" if treatment == "bonded" else " [yellow](restrained)[/yellow]"
                self.console.print(
                    f"  [green]Created bond: {source_atom.atom_name} ↔ "
                    f"{target_atom.atom_name} in {info['residue_label']} ({distance:.2f}Å)[/green]{mode_note}"
                )

    def _convert_site_bonds_to_tleap_commands(self, sites) -> List[str]:
        """
        Convert RedoxSite bonds to tLEaP bond commands using current coord_to_pdb.

        Called immediately after bond definition when coord_to_pdb is guaranteed
        valid. The resulting commands are stored in workspace so the tLEaP
        generator can use them directly, avoiding fragile coord lookups later.
        """
        commands = []
        for site in sites:
            if not hasattr(site, 'bonds') or not site.bonds:
                continue
            for bond in site.bonds:
                # Restrained contacts are held by an MD distance restraint, not a
                # covalent tLEaP bond — skip so no `bond` command is emitted.
                if getattr(bond, 'treatment', 'bonded') == 'restrained':
                    continue
                atom1_info = site.coord_to_pdb.get(bond.atom1_coords)
                atom2_info = site.coord_to_pdb.get(bond.atom2_coords)
                if not atom1_info or not atom2_info:
                    continue
                resid1 = atom1_info['resid']
                atom1 = atom1_info['atom_name']
                resid2 = atom2_info['resid']
                atom2 = atom2_info['atom_name']
                commands.append(f"bond mol.{resid1}.{atom1} mol.{resid2}.{atom2}")
        if commands:
            self.console.print(f"[green]✓ Generated {len(commands)} tLEaP bond command(s)[/green]")
        return commands

    def _reconstruct_redox_sites_from_dicts(self, site_dicts: List[Dict]) -> List:
        """
        Reconstruct RedoxSite objects from serialized dicts.

        When resuming from saved state, RedoxSite objects are serialized as dicts.
        This reconstructs them back to proper RedoxSite objects.
        """
        from proprep.structure_prep.comprehensive_redox_detector import (
            RedoxSite, RedoxSiteAtom, RedoxSiteBond, RedoxCenter, CenterType
        )

        sites = []
        for site_data in site_dicts:
            # Create RedoxSite
            site = RedoxSite(
                site_data.get("site_id", "site_1"),
                site_data.get("structure_id", "structure")
            )
            site.site_type = site_data.get("site_type", "")

            # Reconstruct centers
            for center_data in site_data.get("centers", []):
                center_type_val = center_data.get("center_type", "organic_cofactor")
                try:
                    center_type = CenterType(center_type_val)
                except ValueError:
                    center_type = CenterType.ORGANIC_COFACTOR

                center = RedoxCenter(
                    chain=center_data.get("chain", ""),
                    resname=center_data.get("resname", ""),
                    resid=center_data.get("resid", 0),
                    atom_name=center_data.get("atom_name"),
                    insertion_code=center_data.get("insertion_code", ""),
                    altloc=center_data.get("altloc", ""),
                    coords=tuple(center_data.get("coords", (0, 0, 0))),
                    center_type=center_type,
                    element=center_data.get("element")
                )
                site.add_center(center)

            # Reconstruct atoms
            for atom_data in site_data.get("atoms", []):
                atom = RedoxSiteAtom(
                    chain=atom_data.get("chain", ""),
                    resname=atom_data.get("resname", ""),
                    resid=atom_data.get("resid", 0),
                    atom_name=atom_data.get("atom_name", ""),
                    coords=tuple(atom_data.get("coords", (0, 0, 0))),
                    element=atom_data.get("element", ""),
                    altloc=atom_data.get("altloc", ""),
                    insertion_code=atom_data.get("insertion_code", ""),
                    occupancy=atom_data.get("occupancy", 1.0),
                    bfactor=atom_data.get("bfactor", 0.0)
                )
                site.add_atom(atom)

            # Reconstruct bonds
            for bond_data in site_data.get("bonds", []):
                bond = RedoxSiteBond(
                    atom1_coords=tuple(bond_data.get("atom1_coords", (0, 0, 0))),
                    atom2_coords=tuple(bond_data.get("atom2_coords", (0, 0, 0))),
                    bond_type=bond_data.get("bond_type", "unknown"),
                    chemical_type=bond_data.get("chemical_type", "unknown"),
                    distance=bond_data.get("distance", 0.0),
                    atom1_element=bond_data.get("atom1_element", ""),
                    atom2_element=bond_data.get("atom2_element", ""),
                    atom1_residue_info=bond_data.get("atom1_residue_info", ""),
                    atom2_residue_info=bond_data.get("atom2_residue_info", ""),
                    treatment=bond_data.get("treatment", "bonded")
                )
                site.bonds.append(bond)

            sites.append(site)

        self.console.print(f"[grey50]Reconstructed {len(sites)} RedoxSite(s) from saved state[/grey50]")
        return sites

    def _guess_element(self, atom_name: str) -> str:
        """Guess element from atom name."""
        name = atom_name.strip().upper()
        if name.startswith('H'):
            return 'H'
        elif name.startswith('C'):
            return 'C'
        elif name.startswith('N'):
            return 'N'
        elif name.startswith('O'):
            return 'O'
        elif name.startswith('S'):
            return 'S'
        elif name in ('FE', 'ZN', 'CU', 'MG', 'CA', 'MN', 'CO', 'NI'):
            return name
        else:
            return name[0] if name else 'X'

    def _build_residue_sequence_map(self, pdb_file: Path) -> Dict[Tuple[str, int], int]:
        """
        Build mapping from (chain, resid) to GLOBAL sequence position.

        This is used to track residue identity across tLEaP processing.
        tLEaP removes chain IDs and numbers residues sequentially across
        all chains, so we use a global position counter.

        Returns:
            Dict mapping (chain_id, resid) -> global_sequence_position (0-indexed)
        """
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure("temp", str(pdb_file))
        except Exception:
            return {}

        residue_map = {}
        global_pos = 0  # Global position across all chains
        for model in structure:
            for chain in model:
                for residue in chain:
                    resid = residue.get_id()[1]
                    residue_map[(chain.id, resid)] = global_pos
                    global_pos += 1

        return residue_map

    def _build_position_to_residue_map(self, pdb_file: Path) -> Dict[int, Tuple[str, int]]:
        """
        Build mapping from global sequence position to (chain, resid).

        Used to find residues in tLEaP output by their sequence position.
        tLEaP removes chain IDs and numbers sequentially, so we use global position.

        Returns:
            Dict mapping global_seq_pos -> (chain_id, resid)
        """
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure("temp", str(pdb_file))
        except Exception:
            return {}

        position_map = {}
        global_pos = 0  # Global position across all chains
        for model in structure:
            for chain in model:
                for residue in chain:
                    resid = residue.get_id()[1]
                    position_map[global_pos] = (chain.id, resid)
                    global_pos += 1

        return position_map

    # =========================================================================
    # Step 0g: FF Selection + Complete Atom Typing
    # =========================================================================

    def _step_0g_complete_atom_typing(
        self,
        final_pdb: Path,
        redox_sites: List['RedoxSite'],
        output_dir: Path,
        interactive: bool
    ) -> Tuple[Any, Dict[Tuple[float, float, float], AtomTypeAssignment]]:
        """
        Step 0g: FF selection and complete atom typing on final structure.

        Flow for each residue:
        1. Metals → MetalIonDatabase + user charge/spin prompts
        2. Standard protein residues → FF lib lookup
        3. Non-standard residues:
           a. Try to map to FF residue (user interactive mapping)
           b. If mapping succeeds → use FF types
           c. If mapping fails → LAUNCH small_molecule_parameterizer → use antechamber types

        Returns:
            Tuple of:
            - ff_data: Loaded ForceFieldData
            - type_assignments: {coords → AtomTypeAssignment} with original_type set
        """
        from Bio.PDB import PDBParser

        self.console.print(Panel(
            "[cyan]Atom Type Assignment[/cyan]\n"
            "Assigning force field atom types to all atoms",
            border_style="cyan",
            expand=False
        ))

        # 1. Interactive FF selection
        ff_data = self._select_force_field(interactive)

        # 2. Load MetalIonDatabase
        try:
            from proprep.forcefield_prep.mcpb.metal_ion_database import MetalIonDatabase
            metal_db = MetalIonDatabase()
        except ImportError:
            metal_db = None
            self.console.print("[yellow]MetalIonDatabase not available[/yellow]")

        # 3. Parse final structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('final', str(final_pdb))

        # 4. First pass: identify non-standard residues that need mapping
        unmapped_residues = self._identify_unmapped_residues(structure, ff_data)

        # 5. For each unmapped residue, try interactive mapping or launch small_mol_param
        for resname in unmapped_residues:
            if self._try_interactive_ff_mapping(resname, ff_data, interactive):
                continue  # Mapping succeeded, residue_map/atom_map updated
            else:
                # Can't map to FF → launch small_molecule_parameterizer
                self._parameterize_unmapped_residue(
                    resname, final_pdb, output_dir, interactive
                )

        # 6. Second pass: assign types to all atoms
        type_assignments = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    resid = residue.get_id()[1]

                    for atom in residue:
                        coords = tuple(round(x, 3) for x in atom.coord)
                        atom_name = atom.name.strip()
                        element = atom.element.strip() if atom.element else self._guess_element(atom_name)

                        assignment = AtomTypeAssignment(
                            coords=coords,
                            chain=chain.id,
                            resname=resname,
                            resid=resid,
                            atom_name=atom_name,
                            element=element
                        )

                        # === TIER 1: Metals ===
                        if metal_db and element.upper() in METAL_ELEMENTS:
                            assignment = self._type_metal_atom(
                                assignment, metal_db, redox_sites, interactive
                            )

                        # === TIER 2: In FF (directly or via mapping) ===
                        elif ff_data and (ff_data.has_residue(resname) or resname in self.residue_map):
                            assignment = self._type_organic_atom(
                                assignment, ff_data, self.residue_map, self.atom_map
                            )

                        # === TIER 3: From small_mol_parameterizer (launched above) ===
                        elif resname in self.small_mol_results:
                            assignment = self._type_from_small_mol(
                                assignment, self.small_mol_results[resname]
                            )

                        # === Fallback: Guess from element ===
                        else:
                            assignment.original_type = element.upper()
                            assignment.source = AtomSource.MANUAL
                            assignment.source_detail = "guessed_from_element"

                        type_assignments[coords] = assignment

        # Summary
        typed_count = sum(1 for a in type_assignments.values() if a.original_type)
        self.console.print(f"\n[green]Assigned types to {typed_count}/{len(type_assignments)} atoms[/green]")

        return ff_data, type_assignments

    def _identify_unmapped_residues(self, structure, ff_data) -> set:
        """Find non-standard residues not in FF."""
        unmapped = set()
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    if resname not in STANDARD_RESIDUES and resname not in WATER_RESIDUES:
                        # Check if it's in FF or already mapped
                        in_ff = ff_data and ff_data.has_residue(resname)
                        already_mapped = resname in self.residue_map
                        if not in_ff and not already_mapped:
                            unmapped.add(resname)
        return unmapped

    def _try_interactive_ff_mapping(self, resname: str, ff_data, interactive: bool) -> bool:
        """Try to map a residue to an FF residue. Returns True if successful."""
        if not interactive:
            return False

        self.console.print(f"\n[cyan]Residue '{resname}' not found in force field[/cyan]")
        can_map = confirm_with_context(
            self.processor,
            f"Can you map '{resname}' to an existing FF residue?",
            default=False,
            module="Structure Preprocessor",
            description=f"Map '{resname}' to existing FF residue",
        )

        if can_map:
            ff_resname = prompt_with_context(
                self.processor,
                "Enter FF residue name to map to",
                module="Structure Preprocessor",
                description=f"FF residue to map '{resname}' to",
            )
            if ff_data and ff_data.has_residue(ff_resname):
                self.residue_map[resname] = ff_resname
                self.console.print(f"[green]Mapped {resname} → {ff_resname}[/green]")
                # TODO: atom-level mapping if needed
                return True
            else:
                self.console.print(f"[red]FF residue '{ff_resname}' not found[/red]")

        return False

    def _parameterize_unmapped_residue(
        self,
        resname: str,
        final_pdb: Path,
        output_dir: Path,
        interactive: bool
    ) -> None:
        """Launch small_molecule_parameterizer for a residue that can't be FF-mapped."""
        self.console.print(f"\n[bold]═══ Small Molecule Parameterizer: {resname} ═══[/bold]")
        self.console.print("[grey50]Generating parameters with antechamber[/grey50]")

        # Extract residue to temporary PDB (first instance)
        residue_pdb = output_dir / f"{resname}_for_param.pdb"
        self._extract_residue_by_name(final_pdb, resname, residue_pdb)

        # Store in workspace
        if self.workspace:
            self.workspace.set("ligand_pdb_file", str(residue_pdb))
            self.workspace.set("ligand_resname", resname)

        # Launch small molecule parameterizer
        if self.processor:
            small_mol = self.processor.get_module_instance("Small Molecule Parameterizer")
            if small_mol:
                # Call the process method
                if hasattr(small_mol, 'process'):
                    small_mol.process(self.workspace)
                else:
                    self.console.print("[yellow]Could not run small molecule parameterizer[/yellow]")
                    return

                # Collect results
                mol2_file = self.workspace.get("ligand_mol2_file") if self.workspace else None
                frcmod_file = self.workspace.get("ligand_frcmod_file") if self.workspace else None

                if mol2_file:
                    self.small_mol_results[resname] = {
                        'mol2': mol2_file,
                        'frcmod': frcmod_file,
                        'types': self._parse_mol2_types(mol2_file),
                        'charges': self._parse_mol2_charges(mol2_file)
                    }
                    self.console.print(f"[green]Parameters generated for {resname}[/green]")
                else:
                    self.console.print(f"[red]Failed to generate parameters for {resname}[/red]")

    def _extract_residue_by_name(
        self,
        pdb_file: Path,
        resname: str,
        output_path: Path
    ) -> None:
        """Extract first instance of a residue by name to a new PDB file."""
        from Bio.PDB import PDBParser, PDBIO, Select

        class FirstResidueSelect(Select):
            def __init__(self, target_resname):
                self.target = target_resname
                self.found = False

            def accept_residue(self, residue):
                if self.found:
                    return False
                if residue.get_resname().strip() == self.target:
                    self.found = True
                    return True
                return False

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('struct', str(pdb_file))

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_path), FirstResidueSelect(resname))

    def _select_force_field(self, interactive: bool) -> 'ForceFieldData':
        """
        Load force field data using selections from Step 0b.

        Priority order:
        1. Preprocessing selections from Step 0b (preprocessing_protein_ff, etc.)
        2. tLEaP selection (selected_standard_forcefields) - backward compatibility
        3. Interactive selection as fallback

        Step 0g is now "pure atom typing" - it should use existing selections,
        not prompt for new ones.

        Args:
            interactive: If True, allow interactive fallback if no prior selection

        Returns:
            ForceFieldData instance with loaded force field data
        """
        from proprep.forcefield_prep.forcefield_data import ForceFieldData

        ff_data = ForceFieldData(console=self.console, processor=self.processor)
        loaded_count = 0

        # === Priority 1: Preprocessing selections from Step 0b ===
        if self.workspace:
            protein_ff = self.workspace.get("preprocessing_protein_ff")
            water_model = self.workspace.get("preprocessing_water_model")
            nonstandard_ff = self.workspace.get("preprocessing_nonstandard_ff")

            if protein_ff or water_model or nonstandard_ff:
                self.console.print("[cyan]Using force field selections from Step 0b:[/cyan]")

                try:
                    # Load protein FF
                    if protein_ff:
                        ff_data.load_leaprc(protein_ff)
                        self.console.print(f"  [grey50]{protein_ff}[/grey50]")
                        loaded_count += 1

                    # Load water model
                    if water_model:
                        ff_data.load_leaprc(water_model)
                        self.console.print(f"  [grey50]{water_model}[/grey50]")
                        loaded_count += 1

                    # Load GAFF2 for generated residues
                    has_generate = nonstandard_ff and any(
                        s.get('source') == 'generate' for s in nonstandard_ff.values()
                    )
                    if has_generate:
                        try:
                            ff_data.load_leaprc('leaprc.gaff2')
                            self.console.print(f"  [grey50]leaprc.gaff2[/grey50]")
                            loaded_count += 1
                        except:
                            pass

                    # Load custom lib/frcmod files for non-standard residues
                    if nonstandard_ff:
                        for resname, selection in nonstandard_ff.items():
                            if selection.get('source') == 'custom':
                                lib_file = selection.get('lib_file')
                                frcmod_file = selection.get('frcmod_file')

                                if lib_file:
                                    ff_data.load_lib_file(lib_file)
                                    self.console.print(f"  [grey50]{lib_file}[/grey50]")
                                    loaded_count += 1

                                if frcmod_file:
                                    ff_data.load_frcmod_file(frcmod_file)
                                    self.console.print(f"  [grey50]{frcmod_file}[/grey50]")

                    if loaded_count > 0:
                        self.console.print(f"[green]✓ Loaded {loaded_count} force field(s)[/green]")
                        return ff_data

                except (FileNotFoundError, RuntimeError) as e:
                    self.console.print(f"[yellow]Warning: Error loading Step 0b selections: {e}[/yellow]")

        # === Priority 2: tLEaP selection (backward compatibility) ===
        prior_selection = self.workspace.get("selected_standard_forcefields") if self.workspace else None

        if prior_selection:
            leaprc_names = self._extract_leaprc_names(prior_selection)
            if leaprc_names:
                self.console.print(f"[cyan]Using force field selection from tLEaP:[/cyan]")
                for name in leaprc_names:
                    self.console.print(f"  [grey50]{name}[/grey50]")

                try:
                    for leaprc in leaprc_names:
                        ff_data.load_leaprc(leaprc)
                    self.console.print(f"[green]✓ Loaded {len(leaprc_names)} force field(s)[/green]")
                    return ff_data
                except (FileNotFoundError, RuntimeError) as e:
                    self.console.print(f"[yellow]Warning: Could not load tLEaP selection: {e}[/yellow]")

        # === Priority 3: Interactive fallback ===
        if interactive:
            self.console.print("[yellow]No prior FF selection found. Running interactive selection...[/yellow]")
            try:
                ff_data.select_and_load(self.console)
            except FileNotFoundError as e:
                self.console.print(f"[red]Error: {e}[/red]")
                return ff_data
            except RuntimeError as e:
                self.console.print(f"[red]Error: {e}[/red]")
                self.console.print(
                    "[yellow]Please ensure AMBERHOME environment variable is set.[/yellow]"
                )
                return ff_data
        else:
            # Non-interactive: load default force fields
            try:
                ff_data.load_leaprc('leaprc.protein.ff14SB')
                ff_data.load_leaprc('leaprc.water.tip3p')
                ff_data.load_leaprc('leaprc.gaff2')
                self.console.print(
                    f"[green]✓ Loaded default force fields: ff14SB, TIP3P, GAFF2[/green]"
                )
            except (FileNotFoundError, RuntimeError) as e:
                self.console.print(f"[yellow]Warning: Could not load defaults: {e}[/yellow]")

        return ff_data

    def _extract_leaprc_names(self, selection: dict) -> List[str]:
        """Extract leaprc names from tLEaP force field selection."""
        leaprc_names = []

        for category, value in selection.items():
            if value is None:
                continue
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, dict) and 'name' in item:
                        leaprc_names.append(item['name'])
            elif isinstance(value, dict) and 'name' in value:
                leaprc_names.append(value['name'])

        return leaprc_names

    def _type_metal_atom(
        self,
        assignment: AtomTypeAssignment,
        metal_db,
        redox_sites: List['RedoxSite'],
        interactive: bool
    ) -> AtomTypeAssignment:
        """Tier 1: Type a metal atom using MetalIonDatabase."""
        element = assignment.element.upper()

        # Check if this metal is a center in any RedoxSite
        is_center = False
        for site in redox_sites:
            if hasattr(site, 'centers'):
                for center in site.centers:
                    center_coords = center.coords if hasattr(center, 'coords') else center
                    if center_coords == assignment.coords:
                        is_center = True
                        break
            if is_center:
                break

        metal_info = metal_db.get_metal_info(element) if metal_db else {}

        if interactive and is_center:
            self.console.print(f"\n[bold cyan]Metal Configuration: {element}[/bold cyan]")

            available_charges = metal_info.get('common_charges', [2])
            charge = int_prompt_with_context(
                self.processor,
                f"Formal charge for {element}",
                default=available_charges[0] if available_charges else 2,
                module="Structure Preprocessor",
                description=f"Formal charge for metal {element}",
            )

            spin = int_prompt_with_context(
                self.processor,
                f"Spin multiplicity (unpaired electrons) for {element}",
                default=metal_info.get('default_spin', 0),
                module="Structure Preprocessor",
                description=f"Spin multiplicity for metal {element}",
            )

            assignment.charge = float(charge)
            assignment.spin = spin
        else:
            assignment.charge = float(metal_info.get('default_charge', 2))
            assignment.spin = metal_info.get('default_spin', 0)

        assignment.original_type = metal_info.get('atom_type', element)
        assignment.mass = metal_info.get('mass')
        assignment.vdw_radius = metal_info.get('vdw_radius')
        assignment.vdw_epsilon = metal_info.get('vdw_epsilon')
        assignment.is_center = is_center
        assignment.source = AtomSource.METAL_DATABASE
        assignment.source_detail = "MetalIonDatabase"

        return assignment

    def _type_organic_atom(
        self,
        assignment: AtomTypeAssignment,
        ff_data,
        residue_map: Dict[str, str],
        atom_map: Dict[Tuple[str, str], str]
    ) -> AtomTypeAssignment:
        """Tier 2: Type an organic atom using FF lib files."""
        resname = assignment.resname
        atom_name = assignment.atom_name

        ff_resname = residue_map.get(resname, resname)
        assignment.ff_resname = ff_resname

        ff_atom_name = atom_map.get((resname, atom_name), atom_name)

        if ff_data:
            atom_type = ff_data.get_atom_type(ff_resname, ff_atom_name)

            if atom_type:
                assignment.original_type = atom_type
                assignment.charge = ff_data.get_atom_charge(ff_resname, ff_atom_name)
                assignment.source = AtomSource.FORCE_FIELD
                assignment.source_detail = ff_data.get_source_file(ff_resname)
            else:
                assignment.original_type = ""
                assignment.source = AtomSource.FORCE_FIELD
                assignment.source_detail = f"NOT_FOUND:{ff_resname}:{ff_atom_name}"
        else:
            assignment.original_type = ""
            assignment.source = AtomSource.FORCE_FIELD
            assignment.source_detail = "no_ff_data"

        return assignment

    def _type_from_small_mol(
        self,
        assignment: AtomTypeAssignment,
        small_mol_result: Dict
    ) -> AtomTypeAssignment:
        """Tier 3: Get type from small_molecule_parameterizer results."""
        atom_name = assignment.atom_name

        types = small_mol_result.get('types', {})
        charges = small_mol_result.get('charges', {})

        if atom_name in types:
            assignment.original_type = types[atom_name]
            assignment.charge = charges.get(atom_name)
            assignment.source = AtomSource.ANTECHAMBER
            assignment.source_detail = small_mol_result.get('mol2', 'antechamber')
        else:
            assignment.original_type = ""
            assignment.source = AtomSource.ANTECHAMBER
            assignment.source_detail = f"NOT_FOUND:{atom_name}"

        return assignment

    def _guess_element(self, atom_name: str) -> str:
        """Guess element from atom name."""
        # Common patterns
        name = atom_name.strip().upper()
        if name.startswith('C'):
            return 'C'
        elif name.startswith('N'):
            return 'N'
        elif name.startswith('O'):
            return 'O'
        elif name.startswith('S'):
            return 'S'
        elif name.startswith('H'):
            return 'H'
        elif name.startswith('FE'):
            return 'FE'
        elif name.startswith('ZN'):
            return 'ZN'
        elif name.startswith('MG'):
            return 'MG'
        elif name.startswith('CA'):
            return 'CA'
        else:
            return name[0] if name else 'X'

    # =========================================================================
    # Utility Methods
    # =========================================================================

    def _extract_residues_to_pdb(
        self,
        pdb_file: str,
        residue_keys: List[str],
        output_path: Path
    ) -> None:
        """Extract specified residues to a new PDB file."""
        from Bio.PDB import PDBParser, PDBIO, Select

        class ResidueSelect(Select):
            def __init__(self, keys):
                self.keys = set(keys)

            def accept_residue(self, residue):
                chain = residue.get_parent()
                resname = residue.get_resname().strip()
                resid = residue.get_id()[1]
                key = f"{chain.id}:{resid}:{resname}"
                return key in self.keys

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('struct', pdb_file)

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_path), ResidueSelect(residue_keys))

    def _extract_single_residue(
        self,
        pdb_file: str,
        chain_id: str,
        resid: int,
        output_path: Path
    ) -> None:
        """Extract a single residue to a new PDB file."""
        from Bio.PDB import PDBParser, PDBIO, Select

        class SingleResidueSelect(Select):
            def __init__(self, chain, resid):
                self.chain = chain
                self.resid = resid

            def accept_residue(self, residue):
                return (residue.get_parent().id == self.chain and
                        residue.get_id()[1] == self.resid)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('struct', pdb_file)

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_path), SingleResidueSelect(chain_id, resid))

    def _run_tleap_for_residue(
        self,
        input_pdb: Path,
        output_pdb: Path,
        leaprc: str,
        lib_file: Optional[str] = None,
        frcmod_file: Optional[str] = None,
        ff_resname: Optional[str] = None,
    ) -> bool:
        """
        Run tLEaP on a single residue to add hydrogens.

        Args:
            input_pdb: Path to input PDB
            output_pdb: Path for output PDB with H atoms
            leaprc: Leaprc file to use (e.g., 'leaprc.gaff2')
            lib_file: Optional path to custom lib/off file
            frcmod_file: Optional path to custom frcmod file
            ff_resname: Force field residue name (if different from PDB)

        Returns:
            True if successful
        """
        import subprocess
        import tempfile

        output_dir = output_pdb.parent
        prefix = output_pdb.stem

        # Build tLEaP script
        lines = [f"source {leaprc}"]

        if lib_file:
            lines.append(f'loadoff "{lib_file}"')
        if frcmod_file:
            lines.append(f'loadamberparams "{frcmod_file}"')

        lines.append(f'mol = loadpdb "{input_pdb}"')
        lines.append("check mol")
        lines.append(f'saveamberparm mol "{output_dir}/{prefix}.parm7" "{output_dir}/{prefix}.rst7"')
        lines.append("quit")

        tleap_input = "\n".join(lines)

        # Write tLEaP script
        tleap_file = output_dir / f"{prefix}_tleap.in"
        tleap_file.write_text(tleap_input)

        self.console.print(f"[grey50]Running tLEaP: {tleap_file}[/grey50]")

        # Run tLEaP
        try:
            result = subprocess.run(
                ['tleap', '-f', str(tleap_file)],
                cwd=str(output_dir),
                capture_output=True,
                text=True,
                timeout=60
            )

            # Check for output files
            parm7 = output_dir / f"{prefix}.parm7"
            rst7 = output_dir / f"{prefix}.rst7"

            if parm7.exists() and rst7.exists():
                self._convert_amber_to_pdb(str(parm7), str(rst7), output_pdb)
                if output_pdb.exists():
                    self.console.print(f"[green]✓ Added hydrogens via tLEaP[/green]")
                    return True
                else:
                    self.console.print("[yellow]tLEaP ran but PDB conversion failed[/yellow]")
            else:
                self.console.print(f"[yellow]tLEaP did not produce output files[/yellow]")
                if result.stderr:
                    for line in result.stderr.split('\n')[:5]:
                        if line.strip():
                            self.console.print(f"[grey50]{line}[/grey50]")

        except subprocess.TimeoutExpired:
            self.console.print("[red]tLEaP timed out[/red]")
        except FileNotFoundError:
            self.console.print("[red]tLEaP not found - is AmberTools installed?[/red]")
        except Exception as e:
            self.console.print(f"[red]tLEaP error: {e}[/red]")

        return False

    def _convert_amber_to_pdb(
        self,
        parm7: str,
        rst7: str,
        output_pdb: Path
    ) -> None:
        """Convert AMBER parm7/rst7 to PDB, preserving the topology's names.

        ``-aatm`` writes atom names as the topology holds them, which is how
        the libraries named them. Without it ambpdb translates to PDB v3
        conventions on the way out -- O1P/O2P become OP1/OP2 and so on.

        That matters because this file is an INTERNAL artifact: the next tLEaP
        pass reloads it against the very libraries whose names were just
        rewritten. An externally supplied FAD library using the v2 phosphate
        names built cleanly on the first pass and then failed on the second
        with "Atom .R<FAD 1311>.A<OP2 85> does not have a type", because by
        then the structure and the library disagreed. The topology was correct
        throughout; only the PDB in between was translated.

        On the reported system the flag changes 2 atoms out of 20479 -- exactly
        the two that broke -- with byte-identical column layout and the element
        column preserved.

        Raises:
            RuntimeError: if the conversion fails. The caller writes this file
                and then reads it back, so continuing would fail later on a
                missing file rather than here on the real cause.
        """
        import subprocess

        # No cpptraj fallback: both ship with AmberTools, so if ambpdb is
        # absent cpptraj is too, and it would write a third naming convention
        # for no benefit.
        try:
            result = subprocess.run(
                ['ambpdb', '-aatm', '-p', parm7, '-c', rst7],
                capture_output=True, text=True, check=True,
            )
        except FileNotFoundError as exc:
            raise RuntimeError(
                "ambpdb not found. AmberTools must be installed and on PATH."
            ) from exc
        except subprocess.CalledProcessError as exc:
            detail = (exc.stderr or "").strip().splitlines()
            raise RuntimeError(
                "ambpdb failed to convert the tLEaP output"
                + (f": {detail[-1]}" if detail else "")
            ) from exc

        if not result.stdout.strip():
            raise RuntimeError("ambpdb produced no output")

        output_pdb.write_text(result.stdout)
        self.console.print(f"[green]Converted to PDB: {output_pdb}[/green]")

    def _parse_mol2_types(self, mol2_file: str) -> Dict[str, str]:
        """Parse atom types from mol2 file."""
        types = {}
        try:
            with open(mol2_file, 'r') as f:
                in_atom_section = False
                for line in f:
                    if '@<TRIPOS>ATOM' in line:
                        in_atom_section = True
                        continue
                    elif '@<TRIPOS>' in line:
                        in_atom_section = False
                        continue

                    if in_atom_section and line.strip():
                        parts = line.split()
                        if len(parts) >= 6:
                            atom_name = parts[1]
                            atom_type = parts[5]
                            types[atom_name] = atom_type
        except Exception:
            pass
        return types

    def _parse_mol2_charges(self, mol2_file: str) -> Dict[str, float]:
        """Parse atom charges from mol2 file."""
        charges = {}
        try:
            with open(mol2_file, 'r') as f:
                in_atom_section = False
                for line in f:
                    if '@<TRIPOS>ATOM' in line:
                        in_atom_section = True
                        continue
                    elif '@<TRIPOS>' in line:
                        in_atom_section = False
                        continue

                    if in_atom_section and line.strip():
                        parts = line.split()
                        if len(parts) >= 9:
                            atom_name = parts[1]
                            charge = float(parts[8])
                            charges[atom_name] = charge
        except Exception:
            pass
        return charges
