"""
Metal Site Parameterizer - Main Orchestrator
© 2024 ProPrep Developer. All rights reserved.

Provides a user-friendly interface to the Metal Center Parameter Building of AmberTools.
Orchestrates the complete workflow across multiple steps while maintaining educational
approach and pause/resume capability.

This module acts as the main entry point for metal site parameterization, similar to
how modified_amino_acid_parameterizer.py works for amino acid parameterization.

PROPRIETARY SOFTWARE: This file contains proprietary code belonging to the ProPrep project.
Unauthorized copying, distribution, or modification is strictly prohibited.
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from proprep.utils.prompts import prompt_with_context
from datetime import datetime

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn

from Bio.PDB import Structure

# Import unified force field data (replaces ForceFieldData)
from .forcefield_data import ForceFieldData
from .nonstandard_residue_handler import NonStandardResidueHandler
from .terminal_classifier import TerminalClassifier, TerminalType
from .mcpb import FingerprintGenerator  # Keep fingerprint generator
from .model_builder import SmallModelBuilder, LargeModelBuilder, ModelResidue, ModelBuilder
from .pdb_writer import PDBWriter
from proprep.structure_prep.comprehensive_redox_detector import RedoxSite, METALS

# Import new three-tier parameterization components
from .mcpb.metal_ion_database import (
    MetalIonDatabase, MetalConfig, water_model_from_leaprc,
)
from .mcpb.ligand_grouping import LigandGroupingInterface, LigandGroup
from .mcpb.antechamber_runner import AntechamberRunner, LigandParameters
from .ff_types import MassParameter, BondParameter, AngleParameter, NonbondedParameter


# MCPB systematic atom-type names: M1-M9 then MA-MZ for metals, Y1-Y9 then
# YA-YZ for their ligating atoms. The Y prefix matches MCPB conventions and
# GlobalAtomTypeRegistry. Module level because the numbering has to be
# continued across sites and across runs, so the preprocessor reads back names
# that a previous run wrote and has to resolve them to the same positions.
MCPB_METAL_TYPE_NAMES = [f"M{i}" for i in range(1, 10)] + [f"M{chr(65 + i)}" for i in range(26)]
MCPB_LIGAND_TYPE_NAMES = [f"Y{i}" for i in range(1, 10)] + [f"Y{chr(65 + i)}" for i in range(26)]


def mcpb_type_index(name: str) -> Optional[int]:
    """Position of an M*/Y* type name in its sequence, or None.

    Note the sequence is not the digit: ``Y9`` is position 8 and ``YA`` is 9,
    so the count of ligating atoms used stops matching the label past nine.
    """
    if not name:
        return None
    name = name.strip().upper()
    for names in (MCPB_METAL_TYPE_NAMES, MCPB_LIGAND_TYPE_NAMES):
        if name in names:
            return names.index(name)
    return None


def _collect_restrained_ligands(site):
    """Return (coords, resids) for the non-metal endpoints of *restrained*
    coordinate bonds in a RedoxSite.

    A restrained ligand (e.g. a metal-coordinated water) stays in the QM models
    for correct electronics but is NOT realized as an MCPB bonded ligand: it
    keeps its standard atom types (no Y* renaming), emits no metal-ligand bond,
    and is treated as RESP scaffolding (net-0, no custom mol2). This helper
    identifies those atoms so the typing/RESP gates can exclude them.

    coords: set of rounded (x,y,z) tuples of the restrained ligand atoms.
    resids: set of the restrained ligand residue numbers (waters have unique
            resids, so matching by number is robust to chain relabeling).
    """
    coords, resids = set(), set()
    for bond in getattr(site, 'bonds', None) or []:
        if getattr(bond, 'treatment', 'bonded') != 'restrained':
            continue
        if (bond.atom1_element or '').upper() in METALS:
            lig_coords, lig_info = bond.atom2_coords, bond.atom2_residue_info
        elif (bond.atom2_element or '').upper() in METALS:
            lig_coords, lig_info = bond.atom1_coords, bond.atom1_residue_info
        else:
            continue
        coords.add(tuple(round(float(x), 3) for x in lig_coords))
        if isinstance(lig_info, dict) and lig_info.get('resid') is not None:
            resids.add(lig_info['resid'])
    return coords, resids


class PrmtopParameterProvider:
    """
    Extract force field parameters from a prmtop file using ParmEd.

    Provides the same interface as ForceFieldData but reads parameters
    from an already-built prmtop file instead of parsing leaprc files.
    This is useful when the preprocessing workflow built the prmtop and
    we need to look up standard parameters for organic bonds/angles.

    Note: The prmtop contains original atom types (CT, N, NB, etc.),
    not the M*/Y* renamed types used in MCPB. The caller must use
    type_assignments to map renamed -> original types for lookups.
    """

    def __init__(self, prmtop_path: str, console: Optional[Console] = None):
        """
        Initialize provider from a prmtop file.

        Args:
            prmtop_path: Path to the prmtop (.parm7) file
            console: Optional Rich console for output
        """
        self.prmtop_path = prmtop_path
        self.console = console or Console()

        # Parameter dictionaries (populated by _load_parameters)
        self.mass_parameters: Dict[str, MassParameter] = {}
        self.bond_parameters: Dict[str, BondParameter] = {}
        self.angle_parameters: Dict[str, AngleParameter] = {}
        self.nonbonded_parameters: Dict[str, NonbondedParameter] = {}

        self._load_parameters()

    def _load_parameters(self):
        """Load all parameters from the prmtop file using ParmEd."""
        import parmed

        try:
            parm = parmed.load_file(self.prmtop_path)
        except Exception as e:
            self.console.print(f"[yellow]Could not load prmtop for parameter lookup: {e}[/yellow]")
            return

        # ================================================================
        # Extract mass parameters (from unique atom types)
        # ================================================================
        seen_types = set()
        for atom in parm.atoms:
            atype = atom.type
            if atype not in seen_types:
                seen_types.add(atype)
                self.mass_parameters[atype] = MassParameter(
                    atom_type=atype,
                    mass=atom.mass,
                    polarizability=0.0,
                    comment=f"From prmtop ({atom.residue.name}:{atom.name})",
                    source="prmtop"
                )

        # ================================================================
        # Extract bond parameters
        # ================================================================
        seen_bonds = set()
        for bond in parm.bonds:
            if bond.type is None:
                continue
            type1, type2 = bond.atom1.type, bond.atom2.type
            # Create canonical key (sorted)
            bond_key = tuple(sorted([type1, type2]))
            if bond_key not in seen_bonds:
                seen_bonds.add(bond_key)
                key_str = f"{bond_key[0]}-{bond_key[1]}"
                self.bond_parameters[key_str] = BondParameter(
                    type1=bond_key[0],
                    type2=bond_key[1],
                    force_constant=bond.type.k,
                    eq_length=bond.type.req,
                    source="prmtop"
                )

        # ================================================================
        # Extract angle parameters
        # ================================================================
        seen_angles = set()
        for angle in parm.angles:
            if angle.type is None:
                continue
            type1, type2, type3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
            # Create canonical key (central atom fixed, endpoints may reverse)
            if type1 <= type3:
                angle_key = (type1, type2, type3)
            else:
                angle_key = (type3, type2, type1)
            if angle_key not in seen_angles:
                seen_angles.add(angle_key)
                key_str = f"{angle_key[0]}-{angle_key[1]}-{angle_key[2]}"
                self.angle_parameters[key_str] = AngleParameter(
                    type1=angle_key[0],
                    type2=angle_key[1],
                    type3=angle_key[2],
                    force_constant=angle.type.k,
                    eq_angle=angle.type.theteq,
                    source="prmtop"
                )

        # ================================================================
        # Extract dihedral parameters (proper)
        # ================================================================
        # Group by canonical type quadruple; collect multi-term entries
        self.dihedral_parameters: Dict[Tuple[str, ...], list] = {}
        for dih in parm.dihedrals:
            if dih.improper:
                continue
            if dih.type is None:
                continue
            t1, t2, t3, t4 = dih.atom1.type, dih.atom2.type, dih.atom3.type, dih.atom4.type
            # Canonical key: central pair (t2,t3) stays, but we can reverse the whole thing
            key = (t1, t2, t3, t4)
            rev = (t4, t3, t2, t1)
            # Use the lexicographically smaller key
            canon = min(key, rev)
            term = (dih.type.phi_k, dih.type.phase, abs(dih.type.per))
            if canon not in self.dihedral_parameters:
                self.dihedral_parameters[canon] = []
            if term not in self.dihedral_parameters[canon]:
                self.dihedral_parameters[canon].append(term)

        # ================================================================
        # Extract improper dihedral parameters
        # ================================================================
        self.improper_parameters: Dict[Tuple[str, ...], list] = {}
        for dih in parm.dihedrals:
            if not dih.improper:
                continue
            if dih.type is None:
                continue
            # AMBER improper format: 3rd atom is central
            t1, t2, t3, t4 = dih.atom1.type, dih.atom2.type, dih.atom3.type, dih.atom4.type
            key = (t1, t2, t3, t4)
            term = (dih.type.phi_k, dih.type.phase, abs(dih.type.per))
            if key not in self.improper_parameters:
                self.improper_parameters[key] = []
            if term not in self.improper_parameters[key]:
                self.improper_parameters[key].append(term)

        # ================================================================
        # Extract VDW (nonbonded) parameters
        # ================================================================
        # Gaussian's AMBER block expects Rmin/2 (the AMBER convention).
        # ParmEd's Atom.rmin already returns Rmin/2, so use it directly.
        for atom in parm.atoms:
            atype = atom.type
            if atype not in self.nonbonded_parameters:
                if hasattr(atom, 'rmin') and atom.rmin is not None:
                    radius = atom.rmin
                elif hasattr(atom, 'sigma') and atom.sigma is not None:
                    # sigma = Rmin / 2^(1/6), so Rmin/2 = sigma * 2^(1/6) / 2
                    radius = atom.sigma * (2 ** (1/6)) / 2.0
                else:
                    radius = 0.0

                epsilon = atom.epsilon if hasattr(atom, 'epsilon') and atom.epsilon else 0.0

                self.nonbonded_parameters[atype] = NonbondedParameter(
                    atom_type=atype,
                    radius=radius,
                    well_depth=epsilon,
                    comment=f"From prmtop",
                    source="prmtop"
                )

    def get_mass_parameter(self, atom_type: str) -> Optional[MassParameter]:
        """Get mass parameter for an atom type."""
        return self.mass_parameters.get(atom_type.strip())

    def get_bond_parameter(self, type1: str, type2: str) -> Optional[BondParameter]:
        """
        Get bond parameter for an atom type pair.

        Handles both orderings (A-B and B-A).
        """
        type1 = type1.strip()
        type2 = type2.strip()

        # Try canonical key (sorted)
        types = sorted([type1, type2])
        key = f"{types[0]}-{types[1]}"
        return self.bond_parameters.get(key)

    def get_angle_parameter(self, type1: str, type2: str, type3: str) -> Optional[AngleParameter]:
        """
        Get angle parameter for an atom type triple.

        Handles both orderings (A-B-C and C-B-A).
        """
        type1 = type1.strip()
        type2 = type2.strip()
        type3 = type3.strip()

        # Try both orderings
        for t1, t3 in [(type1, type3), (type3, type1)]:
            key = f"{t1}-{type2}-{t3}"
            if key in self.angle_parameters:
                return self.angle_parameters[key]

        return None

    def get_nonbonded_parameter(self, atom_type: str) -> Optional[NonbondedParameter]:
        """Get nonbonded (VDW) parameter for an atom type."""
        return self.nonbonded_parameters.get(atom_type.strip())

    def get_dihedral_parameters(self, t1: str, t2: str, t3: str, t4: str) -> Optional[list]:
        """Get proper dihedral parameters for a type quadruple. Returns list of (phi_k, phase, per) terms."""
        key = (t1.strip(), t2.strip(), t3.strip(), t4.strip())
        rev = (key[3], key[2], key[1], key[0])
        canon = min(key, rev)
        return self.dihedral_parameters.get(canon)

    def get_improper_parameters(self, t1: str, t2: str, t3: str, t4: str) -> Optional[list]:
        """Get improper dihedral parameters. 3rd atom is central in AMBER convention."""
        key = (t1.strip(), t2.strip(), t3.strip(), t4.strip())
        return self.improper_parameters.get(key)

    # ------------------------------------------------------------------
    # FFParameterReader-compatible interface for the ONIOM writer
    # ------------------------------------------------------------------

    def as_dihedral_params_dict(self) -> Dict[str, list]:
        """Build a dihedral_params dict with string keys and objects with
        .pn, .phase, .pk, .idivf attributes, matching the format the
        ONIOM writer expects from FFParameterReader.

        Internal storage uses tuple keys and (phi_k, phase, per) tuples.
        This converts to string keys "T1-T2-T3-T4" and lightweight objects.
        """
        from proprep.forcefield_prep.mcpb.ff_parameter_reader import DihedralParam

        result: Dict[str, list] = {}
        for type_tuple, terms in self.dihedral_parameters.items():
            key = "-".join(type_tuple)
            result[key] = [
                DihedralParam(
                    atom_type1=type_tuple[0],
                    atom_type2=type_tuple[1],
                    atom_type3=type_tuple[2],
                    atom_type4=type_tuple[3],
                    idivf=1,
                    pk=phi_k,
                    phase=phase,
                    pn=per,
                    source="prmtop",
                )
                for phi_k, phase, per in terms
            ]
        return result

    def get_statistics(self) -> Dict:
        """Get statistics about loaded parameters."""
        return {
            'force_field': 'prmtop',
            'n_mass': len(self.mass_parameters),
            'n_bond': len(self.bond_parameters),
            'n_angle': len(self.angle_parameters),
            'n_dihedral': len(self.dihedral_parameters),
            'n_improper': len(self.improper_parameters),
            'n_nonbonded': len(self.nonbonded_parameters),
        }


DEFAULT_SITE_KEY = "__default__"


def _site_results_key(site) -> str:
    """Workspace slot for a site's MCPB step results."""
    site_id = getattr(site, "site_id", None)
    if not site_id and isinstance(site, dict):
        site_id = site.get("site_id")
    return str(site_id) if site_id else DEFAULT_SITE_KEY


def _is_legacy_flat_results(stored: dict) -> bool:
    """True for the pre-partition shape: step keys at the top level."""
    return any(str(k).startswith("step_") for k in stored)


def _inferred_element_type(global_assignment, residue_atoms) -> str:
    """Placeholder type for an atom with no library entry.

    The element symbol, except for hydrogen, which is typed from the heavy atom
    it is bonded to. 'MO'/'FE'/'S' are not Amber types so they read as
    placeholders and are renamed to M*/Y*; 'H' IS one -- the amide/amine
    hydrogen -- so leaving it there silently gives a hydroxo or thiol proton the
    wrong nonbonded terms rather than failing.
    """
    element = (global_assignment.element or "").strip()
    if element.upper() != "H":
        return element

    from .mcpb.atom_typer import hydrogen_type_from_neighbors

    key = (global_assignment.chain, global_assignment.resid)
    neighbors = residue_atoms.get(key, [])
    return hydrogen_type_from_neighbors(global_assignment.coords, neighbors)


class MetalSiteWorkflowManager:
    """
    Metal site parameterization step implementations.

    Provides step methods called by the structure_preprocessor's WorkflowChecklist.
    Each step method receives its arguments directly from the checklist handler.

    Step methods:
        _run_step1: MCPB atom typing & fingerprint generation
        _run_step2a: Pre-frcmod generation (before QM)
        _run_step2b: Bonded parameter generation - Seminario method (after QM)
        _run_step3a: ESP calculation setup
        _run_step3b: RESP input generation
        _run_step3c: RESP execution
        _run_step3d: Mol2 file generation
        _run_step4: Force field integration
    """

    def __init__(self, console: Console = None, processor=None):
        self.console = console or Console()
        self.processor = processor
        self.original_working_dir = Path.cwd()
        self.logger = logging.getLogger(self.__class__.__name__)
        # step_results is PER SITE. It used to be restored here from a single
        # workspace key that every site shared, so step_1 belonged to whichever
        # site ran last. That cross-wired site 1 to site 2's standard.fingerprint
        # (every atom typed 'XX', ~150 tleap errors), fitted site 1's RESP
        # against site 2's charge constraint, and still shows up as "Step-1
        # records charge -3 for this site but its ESP was computed at -1".
        #
        # The site is not known at construction -- callers assign
        # provided_redox_site immediately afterwards -- so the restore happens
        # in that setter. Until then this is the standalone workflow's bucket.
        self.step_results: Dict[str, Any] = {}
        self._provided_redox_site = None
        self._site_key = DEFAULT_SITE_KEY
        if processor:
            self._restore_step_results()

    # ------------------------------------------------------------------ #
    # per-site step_results
    # ------------------------------------------------------------------ #

    @property
    def provided_redox_site(self):
        """The site this manager is working on, or None."""
        return self._provided_redox_site

    @provided_redox_site.setter
    def provided_redox_site(self, site):
        """Assigning the site selects which site's step_results to use.

        This is the only point where the manager learns its identity, and
        every caller assigns it right after construction.
        """
        self._provided_redox_site = site
        site_key = _site_results_key(site)
        if site_key != self._site_key:
            self._site_key = site_key
            self._restore_step_results()

    def _restore_step_results(self):
        """Load THIS site's step_results from the workspace."""
        if not self.processor:
            return
        try:
            workspace = self.processor._get_workspace()
            saved = workspace.get("mcpb_step_results", {})
        except Exception:
            return
        if not isinstance(saved, dict):
            return

        bucket = saved.get(self._site_key)
        if bucket is None and _is_legacy_flat_results(saved):
            # Written before results were partitioned: one site's results with
            # no label. Usable only for the standalone workflow, where there is
            # just one site; a labelled site must not adopt them.
            bucket = saved if self._site_key == DEFAULT_SITE_KEY else None

        self.step_results = dict(bucket) if isinstance(bucket, dict) else {}
        if self.step_results:
            self.logger.debug(
                f"Restored step_results for {self._site_key}: "
                f"{list(self.step_results.keys())}")

    def _save_step_results(self):
        """Persist step_results to workspace for resume support."""
        if self.processor:
            try:
                workspace = self.processor._get_workspace()
                # Filter out non-serializable objects (like ff_data)
                serializable = {}
                for key, value in self.step_results.items():
                    if isinstance(value, dict):
                        filtered = {}
                        for k, v in value.items():
                            # Skip large non-serializable objects
                            if k in ('ff_data', 'type_assignments', 'fitted_charges_raw'):
                                continue
                            # Convert Path objects to strings
                            if isinstance(v, Path):
                                filtered[k] = str(v)
                            elif isinstance(v, dict):
                                filtered[k] = {
                                    sk: str(sv) if isinstance(sv, Path) else sv
                                    for sk, sv in v.items()
                                    if not callable(sv)
                                }
                            elif not callable(v):
                                filtered[k] = v
                        serializable[key] = filtered
                    else:
                        serializable[key] = value
                # Merge into this site's slot, leaving other sites alone.
                stored = workspace.get("mcpb_step_results", {})
                if not isinstance(stored, dict) or _is_legacy_flat_results(stored):
                    stored = {}
                stored = dict(stored)
                stored[self._site_key] = serializable
                workspace.set("mcpb_step_results", stored)
            except Exception as e:
                self.logger.debug(f"Could not save step_results: {e}")

    def _get_redox_site_from_workspace(self) -> Optional[RedoxSite]:
        """
        Get RedoxSite object from workspace.

        RedoxSite is stored by pdb_filter under 'detected_redox_sites' key.

        Returns:
            RedoxSite object or None if not found
        """
        if not self.processor:
            self.console.print("[red]❌ Processor not available[/red]")
            return None

        # Get workspace using ProPrep convention
        workspace = self.processor._get_workspace()

        # Get detected redox sites from workspace (stored by pdb_filter)
        detected_redox_sites = workspace.get("detected_redox_sites")

        if detected_redox_sites and len(detected_redox_sites) > 0:
            self.console.print(f"[green]✅ Found {len(detected_redox_sites)} RedoxSite(s) in workspace[/green]")
            if len(detected_redox_sites) > 1:
                self.console.print("[yellow]Multiple RedoxSites found. Using first site.[/yellow]")
                # TODO: Add interactive selection when multiple sites exist
            return detected_redox_sites[0]
        else:
            self.console.print("[red]❌ No RedoxSite found in workspace[/red]")
            self.console.print("[yellow]Please run PDB Filter with redox site detection enabled first[/yellow]")
            return None

    def _get_structure_pdb_file(self) -> Optional[str]:
        """
        Get PDB file path from workspace using StructureSelector.

        Uses the centralized StructureSelector to find the best available
        structure based on priority (repaired > filtered > RCSB > local, etc.)

        Returns:
            Absolute path to PDB file or None
        """
        if not self.processor:
            self.console.print("[red]❌ No processor available for workspace access[/red]")
            return None

        try:
            from proprep.utils.structure_selector import StructureSelector

            # Get workspace using ProPrep convention
            workspace = self.processor._get_workspace()

            # Use StructureSelector for standardized structure access
            selector = StructureSelector(workspace, self.console, self.processor)
            pdb_file = selector.get_structure(silent=False)

            if pdb_file:
                # Ensure absolute path (relative to original working directory if needed)
                p = Path(pdb_file)
                if not p.is_absolute():
                    pdb_file = str((self.original_working_dir / p).resolve())
                return pdb_file

            self.console.print("[red]❌ No valid PDB file found in workspace[/red]")
            return None

        except ImportError:
            self.console.print("[yellow]⚠️  StructureSelector not available, using fallback[/yellow]")
            return self._get_structure_pdb_file_fallback()

    def _get_structure_pdb_file_fallback(self) -> Optional[str]:
        """Fallback method for getting PDB file if StructureSelector unavailable."""
        workspace = self.processor._get_workspace()

        # Check common workspace keys
        pdb_keys = [
            "repaired_pdb_file", "filtered_pdb_file", "rcsb_pdb_file",
            "local_pdb_file", "alphafold_pdb_file", "pdb_file"
        ]

        for key in pdb_keys:
            pdb_path = workspace.get(key)
            if pdb_path and Path(pdb_path).exists():
                self.console.print(f"[green]✅ Using {key}: {Path(pdb_path).name}[/green]")
                p = Path(pdb_path)
                if not p.is_absolute():
                    return str((self.original_working_dir / p).resolve())
                return str(p)

        self.console.print("[red]❌ No valid PDB file found in workspace[/red]")
        return None

    def _populate_serial_numbers(self, redox_site: RedoxSite, pdb_file: str):
        """
        Populate PDB serial numbers for RedoxSite atoms from PDB file.

        This is critical for MCPB fingerprint generation, which uses PDB serial
        numbers as atom IDs for matching with QM model files.

        Args:
            redox_site: RedoxSite to populate serial numbers for
            pdb_file: Path to PDB file

        Raises:
            Warning if atoms don't have serial numbers
        """
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', pdb_file)

        # Coordinate key used to match RedoxSite atoms to PDB atoms. Both sides
        # MUST be normalized identically: coerce to float and round to 3 dp.
        # RedoxSite atom coords can come back from a persisted workflow_state as
        # full-precision *strings* ('30.838055'), which would never compare equal
        # to the PDB's rounded-float keys — leaving every atom without a serial,
        # which silently falls the fingerprint back to synthetic 10000+ IDs and
        # breaks MCPB bond/angle extraction. round(float(...), 3) on both sides
        # tolerates string coords and 6dp-vs-3dp precision alike.
        def _coord_key(coords):
            return tuple(round(float(v), 3) for v in coords)

        # Build coordinate → serial number map from PDB
        coord_to_serial = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue.get_atoms():
                        coord_to_serial[_coord_key(atom.coord)] = atom.serial_number

        # Populate serial numbers in RedoxSite atoms
        atoms_with_serial = 0
        atoms_without_serial = 0

        for atom in redox_site.atoms:
            try:
                key = _coord_key(atom.coords)
            except (TypeError, ValueError):
                key = None
            if key is not None and key in coord_to_serial:
                atom.properties['serial_number'] = coord_to_serial[key]
                atoms_with_serial += 1
            else:
                atoms_without_serial += 1
                self.logger.warning(
                    f"Atom {atom.atom_name} in {atom.resname}:{atom.resid} "
                    f"at coords {atom.coords} not found in PDB file"
                )

        self.console.print(
            f"[green]✓ Populated serial numbers: {atoms_with_serial} atoms[/green]"
        )
        if atoms_without_serial > 0:
            self.console.print(
                f"[yellow]⚠ {atoms_without_serial} atoms missing serial numbers[/yellow]"
            )

        # Debug: Verify serial numbers are actually in properties
        test_atom = redox_site.atoms[0]
        if 'serial_number' in test_atom.properties:
            self.console.print(f"[grey50]Debug: First atom has serial_number = {test_atom.properties['serial_number']}[/grey50]")
        else:
            self.console.print(f"[red]Debug: First atom MISSING serial_number in properties![/red]")

    def _save_atom_type_assignments(self, type_assignments: Dict, output_file: Path):
        """Save atom type assignments to JSON file."""
        # Convert assignments to serializable format
        serialized = {}

        for coords, assignment in type_assignments.items():
            key = f"{coords[0]:.3f},{coords[1]:.3f},{coords[2]:.3f}"

            # Handle both old MCPB format and new dict format
            if isinstance(assignment, dict):
                serialized[key] = assignment
            else:
                # Old MCPB AtomTypeAssignment object format
                serialized[key] = {
                    "chain": assignment.chain,
                    "resname": assignment.resname,
                    "resid": assignment.resid,
                    "atom_name": assignment.atom_name,
                    "element": getattr(assignment, 'element', ''),
                    "original_type": assignment.original_type,
                    "renamed_type": assignment.renamed_type,
                    "is_center": assignment.is_center,
                    "is_metal_ligand": assignment.is_metal_ligand,
                    "terminal_type": assignment.terminal_type.value if hasattr(assignment, 'terminal_type') and hasattr(assignment.terminal_type, 'value') else '',
                    "library_source": getattr(assignment, 'library_source', ''),
                    "charge": getattr(assignment, 'charge', 0.0)
                }

        with open(output_file, 'w') as f:
            json.dump(serialized, f, indent=2)

    def _save_validation_report(self, report: Dict, output_file: Path):
        """Save validation report to JSON file."""
        with open(output_file, 'w') as f:
            json.dump(report, f, indent=2)

        self.console.print(f"[grey50]Saved validation report: {output_file.name}[/grey50]")

    def set_processor(self, processor):
        """Set processor for workspace access."""
        self.processor = processor

    def _run_step1(self, residue_name: str, residues: List, output_dir: Path,
                   interactive: bool, global_types: Dict = None,
                   ff_data=None, metal_type_start: int = 0,
                   ligand_type_start: int = 0) -> Dict[str, Any]:
        """
        Run Step 1: MCPB Atom Typing and Fingerprint Generation.

        Supports three modes for type assignments:
        1. Global types from GlobalAtomTypeRegistry (multi-site coordination)
        2. Preprocessing types from prmtop (preprocessing workflow)
        3. Interactive force field selection (legacy/standalone)

        Process:
        1. Get RedoxSite from workspace (or provided_redox_site)
        2. Get atom types (from global_types, preprocessing, or interactive)
        3. Build small and large QM models
        4. Generate PDB files
        5. Generate MCPB fingerprint files
        6. Generate Gaussian input files

        Args:
            residue_name: Name of metal site
            residues: List of residue objects (not used in new workflow)
            output_dir: Output directory
            interactive: Enable user confirmations
            global_types: Pre-assigned types from GlobalAtomTypeRegistry (optional)
            ff_data: ForceFieldData with charges (optional, used with global_types)

        Returns:
            Dict with success status and output files
        """
        # Determine banner based on type source
        if global_types:
            self.console.print("[bold cyan]▸ MCPB atom typing (using global types)[/bold cyan]")
        else:
            self.console.print("[bold cyan]▸ MCPB atom typing & fingerprint generation[/bold cyan]")

        try:
            # ================================================================
            # 1. Get RedoxSite and PDB file
            # ================================================================
            # Prefer directly-passed RedoxSite over workspace lookup
            redox_site = getattr(self, 'provided_redox_site', None)

            if redox_site:
                self.console.print(f"[green]✓ Using directly-provided RedoxSite: {redox_site.site_id}[/green]")
            else:
                # Fall back to workspace lookup
                redox_site = self._get_redox_site_from_workspace()

            if not redox_site:
                return {
                    "success": False,
                    "message": "No RedoxSite available. Run PDB Filter with redox detection first."
                }

            pdb_file = self._get_structure_pdb_file()

            if not pdb_file:
                return {
                    "success": False,
                    "message": "No PDB file available in workspace"
                }

            self.console.print(f"[cyan]RedoxSite has {len(redox_site.atoms)} atoms[/cyan]")
            self.console.print(f"[cyan]Using PDB: {Path(pdb_file).name}[/cyan]")

            # Populate PDB serial numbers for MCPB fingerprint generation
            self.console.print("\n[bold]Populating PDB serial numbers...[/bold]")
            self._populate_serial_numbers(redox_site, pdb_file)

            # Create models output directory
            step1_dir = output_dir / "models"
            step1_dir.mkdir(parents=True, exist_ok=True)

            # ================================================================
            # 2. Determine type assignment source (3 modes)
            # ================================================================
            # Priority: global_types > preprocessing > interactive

            have_global_types = global_types is not None and len(global_types) > 0

            preprocessing_result = getattr(self, 'preprocessing_result', None)
            have_preprocessing = (
                preprocessing_result is not None and
                getattr(preprocessing_result, 'success', False) and
                getattr(preprocessing_result, 'type_assignments', None) is not None
            )

            # Determine ff_data source
            if have_global_types:
                # Global types path - ff_data passed as parameter
                self.console.print("\n[green]✓ Using global atom types from GlobalAtomTypeRegistry[/green]")
                self.console.print(f"[grey50]  Pre-assigned atoms: {len(global_types)}[/grey50]")
                # ff_data is passed as parameter (may be None)
            elif have_preprocessing:
                # Preprocessing path - types + charges from prmtop
                self.console.print("\n[green]✓ Using atom data from preprocessing (prmtop)[/green]")
                self.console.print(f"[grey50]  Pre-assigned atoms: {len(preprocessing_result.type_assignments)}[/grey50]")
                ff_data = getattr(preprocessing_result, 'ff_data', None)
            else:
                # Legacy/interactive path: need to load force field data
                self.console.print("\n[bold]Force Field Selection[/bold]")
                ff_data = ForceFieldData(console=self.console, processor=self.processor)
                ff_data.select_and_load(self.console)

            # ================================================================
            # 3. Identify residues in RedoxSite
            # ================================================================
            # Key by (chain, resid) so a binuclear site's two MN (and two WAT)
            # show as distinct residues with their IDs, not a collapsed
            # "MN, WAT" set that hides which/how many.
            residues_in_site = {}
            for atom in redox_site.atoms:
                key = (getattr(atom, 'chain', ''), getattr(atom, 'resid', None))
                residues_in_site[key] = atom.resname

            def _res_sort_key(item):
                (chain, resid), _ = item
                return (str(chain), resid if resid is not None else -1)

            residue_labels = []
            for (chain, resid), resname in sorted(residues_in_site.items(), key=_res_sort_key):
                if chain and str(chain).strip():
                    residue_labels.append(f"{resname} {chain}:{resid}")
                elif resid not in (None, ''):
                    residue_labels.append(f"{resname} {resid}")
                else:
                    residue_labels.append(resname)

            self.console.print(f"\n[cyan]Residues in RedoxSite: {', '.join(residue_labels)}[/cyan]")

            # ================================================================
            # 4. Handle missing residues (skip if preprocessing handled this)
            # ================================================================
            if have_preprocessing and getattr(preprocessing_result, 'type_assignments', None):
                self.console.print("[green]✓ Using atom type assignments from preprocessing[/green]")
                self.console.print(f"[grey50]  Pre-assigned types: {len(preprocessing_result.type_assignments)}[/grey50]")
            else:
                # Need to check for and handle missing residues
                missing_residues = set()
                for res_name in residues_in_site:
                    if not ff_data.has_residue(res_name):
                        missing_residues.add(res_name)

                # Initialize MetalIonDatabase to identify metal-containing residues
                metal_db = MetalIonDatabase(water_model=self._water_model(), logger=self.logger)

                # Identify residues that will be handled by three-tier approach
                # These include: metal-containing residues and inorganic ligands
                redox_site_handled_residues = set()
                for atom in redox_site.atoms:
                    # Check if this atom is a metal (Tier 1) by resname/atomname OR element
                    is_metal_residue = metal_db.is_metal(atom.resname, atom.atom_name)
                    is_metal_elem = metal_db.is_metal_element(atom.element)

                    if is_metal_residue or is_metal_elem:
                        redox_site_handled_residues.add(atom.resname)
                    # Check if atom is in a residue not in force field (will be Tier 3)
                    elif atom.resname not in residues_in_site or not ff_data.has_residue(atom.resname):
                        # This is likely an inorganic ligand that will be handled by antechamber
                        redox_site_handled_residues.add(atom.resname)

                # Filter out residues that three-tier approach will handle
                missing_for_handler = missing_residues - redox_site_handled_residues

                if redox_site_handled_residues & missing_residues:
                    self.console.print(
                        f"[cyan]ℹ️  {len(redox_site_handled_residues & missing_residues)} residue(s) will be handled by three-tier approach: "
                        f"{', '.join(sorted(redox_site_handled_residues & missing_residues))}[/cyan]"
                    )

                if missing_for_handler:
                    self.console.print(f"[yellow]⚠️  {len(missing_for_handler)} residue(s) not in force field: {', '.join(sorted(missing_for_handler))}[/yellow]")

                    # Collect PDB atom names for each missing residue (for matching)
                    pdb_residue_atoms = {}
                    for res_name in missing_for_handler:
                        atoms_in_residue = set()
                        for atom in redox_site.atoms:
                            if atom.resname == res_name:
                                atoms_in_residue.add(atom.atom_name)
                        if atoms_in_residue:
                            pdb_residue_atoms[res_name] = atoms_in_residue

                    # Handle missing residues interactively
                    handler = NonStandardResidueHandler(
                        ff_data,
                        self.console,
                        redox_site=redox_site,
                        pdb_file=pdb_file,
                        processor=self.processor,
                    )
                    if not handler.handle_missing_residues(missing_for_handler, pdb_residue_atoms, interactive):
                        return {
                            "success": False,
                            "message": f"Missing force field data for: {', '.join(sorted(missing_for_handler))}"
                        }

            # ================================================================
            # 5. Assign atom types to RedoxSite atoms
            # ================================================================
            self.console.print("\n[bold]Assigning atom types to structure...[/bold]")

            # Priority: global_types > preprocessing > interactive
            if have_global_types:
                # Convert GlobalTypeAssignment to AtomTypeAssignment
                type_assignments = self._convert_global_types(global_types, ff_data)
                self.console.print(f"[green]✓ Converted {len(type_assignments)} global type assignments[/green]")

                # Show assigned types summary
                metal_types = [a.renamed_type for a in type_assignments.values() if a.is_center]
                ligand_types = [a.renamed_type for a in type_assignments.values() if a.is_metal_ligand]
                if metal_types:
                    self.console.print(f"[cyan]Metal types: {', '.join(set(metal_types))}[/cyan]")
                if ligand_types:
                    self.console.print(f"[cyan]Ligand types: {', '.join(set(ligand_types))}[/cyan]")

            elif have_preprocessing and preprocessing_result.type_assignments:
                # Convert preprocessing type_assignments (coord-keyed) to step1 format (tuple-keyed)
                type_assignments = self._convert_preprocessing_types(
                    preprocessing_result.type_assignments, redox_site
                )
                self.console.print(f"[green]✓ Using {len(type_assignments)} pre-assigned types from preprocessing[/green]")

                # Collect metal charges from user (metals weren't in prmtop)
                type_assignments = self._collect_metal_charges(
                    type_assignments, redox_site, interactive
                )
            else:
                # Legacy/interactive path: classify terminal residues for FF lookup
                self.console.print("\n[bold]Classifying terminal residues...[/bold]")
                terminal_classifier = TerminalClassifier(pdb_file, redox_site, self.console)

                type_assignments = self._assign_atom_types(redox_site, ff_data, terminal_classifier, output_dir)
                self.console.print(f"[green]✅ Assigned types for {len(type_assignments)} atoms[/green]")

            # ================================================================
            # 7. Apply systematic renaming for metals and ligands
            # ================================================================
            self.console.print("\n[bold]Applying systematic renaming...[/bold]")

            type_assignments = self._apply_systematic_renaming(
                redox_site, type_assignments,
                metal_start=metal_type_start, ligand_start=ligand_type_start,
            )

            # Count renamed atoms
            renamed_count = sum(1 for assignment in type_assignments.values() if assignment.get('renamed'))

            self.console.print(f"[green]✅ Renamed {renamed_count} metal/ligand atoms[/green]")

            # Show summary. Count actual metal ATOMS (is_center), not centers: an
            # organometallic cofactor / metal cluster is ONE center but can embed
            # several metals (a 2Fe-2S has two Fe), so len(centers) misreports it
            # as "1 metal". Mirror the ligand count just below.
            metal_count = sum(1 for assignment in type_assignments.values() if assignment.get('is_center', False))
            ligand_count = sum(1 for assignment in type_assignments.values() if assignment.get('is_metal_ligand', False))

            self.console.print(f"\n[cyan]Atom Type Summary:[/cyan]")
            self.console.print(f"  • Total atoms: {len(type_assignments)}")
            self.console.print(f"  • Metal atoms: {metal_count}")
            self.console.print(f"  • Ligand atoms: {ligand_count}")
            self.console.print(f"  • Renamed atoms: {renamed_count}")

            # ================================================================
            # 8. Generate fingerprint files
            # ================================================================
            # ================================================================
            # 9. Build QM Models (Small & Large)
            # ================================================================
            self.console.print("\n[bold]Building QM Models...[/bold]")

            # Build small model
            small_model, small_residues = self._build_small_model_interactive(
                redox_site, pdb_file, step1_dir, interactive
            )

            # Build large model
            large_model, large_residues = self._build_large_model_interactive(
                redox_site, pdb_file, step1_dir, interactive, small_residues
            )

            # Generate PDB files for QM calculations
            self.console.print("\n[bold]Generating PDB files for QM calculations...[/bold]")

            # Use unified ForceFieldData (replaces ForceFieldLibrary)
            pdb_writer = PDBWriter(pdb_file, ff_data=ff_data, console=self.console)

            # Create type_assignments for small model (may include gap residues)
            small_type_assignments = type_assignments.copy()

            # Add gap residue atoms from preprocessing (they have charges from prmtop)
            if have_preprocessing and preprocessing_result.type_assignments:
                self._add_gap_residue_atoms_to_assignments(
                    small_model,
                    redox_site,
                    preprocessing_result.type_assignments,
                    small_type_assignments,
                    pdb_file
                )

            small_pdb = step1_dir / "small.pdb"
            pdb_writer.write_pdb(
                small_model,
                small_type_assignments,
                str(small_pdb),
                "Small model for bonded parameters",
                redox_site  # For CONECT records
            )

            # Create type_assignments for large model (includes gap residues)
            large_type_assignments = type_assignments.copy()

            # Add gap residue atoms from preprocessing (they have charges from prmtop)
            if have_preprocessing and preprocessing_result.type_assignments:
                self._add_gap_residue_atoms_to_assignments(
                    large_model,
                    redox_site,
                    preprocessing_result.type_assignments,
                    large_type_assignments,
                    pdb_file
                )

            large_pdb = step1_dir / "large.pdb"
            pdb_writer.write_pdb(
                large_model,
                large_type_assignments,
                str(large_pdb),
                "Large model for RESP charges",
                redox_site  # For CONECT records
            )

            # ================================================================
            # 10. Generate fingerprint files (must match PDB files)
            # ================================================================
            self.console.print("\n[bold]Generating MCPB fingerprint files...[/bold]")

            # A withheld cluster's own bonds (Fe-S inside FES, Mo-S/O inside
            # MOS, the O-H of a hydroxo) are in neither source the Seminario
            # step draws on, so nothing ever derived a force constant for them.
            # Perceive them from the geometry and record them on the site; the
            # LINK writer below is what carries them forward.
            from proprep.structure_prep.comprehensive_redox_detector import (
                perceive_cluster_internal_bonds,
            )
            n_internal = perceive_cluster_internal_bonds(redox_site)
            if n_internal:
                self.console.print(
                    f"[grey50]Perceived {n_internal} bond(s) inside the metal "
                    f"cluster residue(s); these get force constants too[/grey50]")

            # Create fingerprint generator with suppressed logging
            import logging
            fp_logger = logging.getLogger('proprep.forcefield_prep.mcpb.fingerprint_generator')
            fp_logger.setLevel(logging.WARNING)

            fp_gen = FingerprintGenerator()

            # Standard fingerprint (for atom typing - MCPB Steps 3-4)
            # CRITICAL: Uses ORIGINAL residue names (GLY, CYS, etc.), NOT capped versions (ACE, NME)
            # The small model PDB file has caps for QM, but standard fingerprint has original names
            # Atom typing works via PDB serial number lookup between the two files
            standard_fp = step1_dir / "standard.fingerprint"
            self.console.print(f"[grey50]Writing standard fingerprint ({len(redox_site.atoms)} atoms with original residue names)...[/grey50]")

            fp_gen.write_fingerprint(
                site=redox_site,  # Original site without capping groups
                type_assignments=type_assignments,  # Original type assignments
                output_file=str(standard_fp),
                include_links=True
            )

            # Large model fingerprint (for RESP fitting in Step 2)
            large_fp = step1_dir / "large.fingerprint"
            self.console.print(f"[grey50]Writing fingerprint for large model ({len(large_model.model_residues)} residues)...[/grey50]")

            # Add capping group atoms for large model (gap residues already added above)
            self._add_capping_atoms_to_assignments(large_model, pdb_writer, large_type_assignments)

            large_site = self._create_site_from_model(large_model, redox_site, pdb_writer)

            fp_gen.write_large_fingerprint(
                site=large_site,
                output_file=str(large_fp)
            )

            self.console.print(f"[green]✅ Fingerprint files generated:[/green]")
            self.console.print(f"  • Standard (original residues): {standard_fp.name}")
            self.console.print(f"  • Large (with caps): {large_fp.name}")

            # Validate fingerprint
            fp_validation = fp_gen.validate_fingerprint(str(standard_fp))

            if fp_validation['valid']:
                self.console.print(f"[green]✅ Fingerprint validation passed[/green]")
                self.console.print(f"  • Atom lines: {fp_validation['atom_lines']}")
                self.console.print(f"  • LINK lines: {fp_validation['link_lines']}")
                self.console.print(f"  • Renamed atoms: {fp_validation['renamed_atoms']}")
            else:
                self.console.print(f"[red]❌ Fingerprint validation failed[/red]")
                for error in fp_validation['errors']:
                    self.console.print(f"  [red]• {error}[/red]")

                return {
                    "success": False,
                    "message": "Fingerprint validation failed",
                    "validation_errors": fp_validation['errors']
                }

            # ================================================================
            # 11. Generate Gaussian input files for QM calculations
            # ================================================================
            self.console.print("\n[bold cyan]═══ Gaussian Input File Generation ═══[/bold cyan]")
            self.console.print("[grey50]These files are needed for QM calculations in later steps.[/grey50]\n")

            from .mcpb.qm_interface import QMInterface, QMSoftware, QMCalculationMode
            from proprep.utils.prompts import prompt_with_context, int_prompt_with_context, confirm_with_context

            # Read suggested values from PDB files
            small_suggested_charge = 0
            try:
                with open(small_pdb) as f:
                    for line in f:
                        if line.startswith("REMARK") and "Total charge:" in line:
                            charge_str = line.split("Total charge:")[1].strip()
                            small_suggested_charge = round(float(charge_str))
                            break
            except Exception:
                small_suggested_charge = 0

            large_suggested_charge = 0
            try:
                with open(large_pdb) as f:
                    for line in f:
                        if line.startswith("REMARK") and "Total charge:" in line:
                            charge_str = line.split("Total charge:")[1].strip()
                            large_suggested_charge = round(float(charge_str))
                            break
            except Exception:
                large_suggested_charge = small_suggested_charge

            # Calculate suggested multiplicity from metal spins
            total_spin = 0
            for assignment in type_assignments.values():
                if hasattr(assignment, 'spin') and assignment.spin:
                    total_spin += assignment.spin
                elif isinstance(assignment, dict) and assignment.get('spin'):
                    total_spin += assignment['spin']

            total_S = total_spin / 2.0
            suggested_multiplicity = int(2 * abs(total_S) + 1)
            if suggested_multiplicity < 1:
                suggested_multiplicity = 1

            # ----------------------------------------------------------------
            # Small Model: Optimization + Frequency (for Seminario method)
            # ----------------------------------------------------------------
            from rich.panel import Panel
            self.console.print(Panel(
                "[bold cyan]Small Model: Optimization + Frequency Calculation[/bold cyan]\n\n"
                "This calculation optimizes the geometry and computes vibrational frequencies\n"
                "for the metal coordination sphere. The Hessian matrix is used by the Seminario\n"
                "method to derive bond and angle force constants.\n\n"
                "[bold]Default route keywords:[/bold]\n"
                "  • [bold]Opt[/bold]                      Geometry optimization to minimum energy\n"
                "  • [bold]Freq[/bold]                     Frequency calculation (provides Hessian)\n"
                "  • [bold]Geom=PrintInputOrient[/bold]    Print coordinates in standard orientation\n"
                "  • [bold]Integral=(Grid=UltraFine)[/bold] High-quality integration grid\n"
                "  • [bold]IOp(7/33=1)[/bold]              Save Cartesian force constants (REQUIRED)\n\n"
                "[grey50]IOp(7/33=1) is essential - it saves the Hessian in Cartesian coordinates\n"
                "which the Seminario method needs to compute bond/angle force constants.[/grey50]",
                title="Small Model Configuration",
                border_style="cyan",
                expand=False
            ))
            self.console.print()

            small_memory_gb = int_prompt_with_context(
                self.processor,
                "Memory allocation (GB)",
                default=4,
                module="Metal Site Parameterizer",
                description="Small model memory"
            )

            small_nproc = int_prompt_with_context(
                self.processor,
                "Number of processors",
                default=4,
                module="Metal Site Parameterizer",
                description="Small model processors"
            )

            self.console.print(f"\n[grey50]Suggested charge from PDB: {small_suggested_charge}[/grey50]")
            self._print_withheld_cluster_charge_note(redox_site, type_assignments)
            small_charge = int_prompt_with_context(
                self.processor,
                "Total charge of small model",
                default=small_suggested_charge,
                module="Metal Site Parameterizer",
                description="Small model charge"
            )

            self.console.print(f"[grey50]Suggested multiplicity (2S+1): {suggested_multiplicity}[/grey50]")
            small_mult = int_prompt_with_context(
                self.processor,
                "Spin multiplicity",
                default=suggested_multiplicity,
                module="Metal Site Parameterizer",
                description="Small model multiplicity"
            )

            # Anion solvation check
            small_scrf = ""
            if small_charge < 0:
                self.console.print(Panel(
                    f"[bold yellow]Anion detected (charge = {small_charge}):[/bold yellow]\n\n"
                    "DFT self-interaction error causes electrons to be too delocalized — or\n"
                    "even unbound — for anions in vacuum. This frequently leads to SCF\n"
                    "convergence failure or unphysical geometries.\n\n"
                    "Adding implicit solvation (SCRF) stabilizes the charge distribution\n"
                    "and is [bold]strongly recommended[/bold] for negatively charged systems.",
                    title="Anion Solvation Warning",
                    border_style="yellow",
                    expand=False,
                ))

                if confirm_with_context(
                    self.processor,
                    "Add implicit solvation (SCRF)?",
                    default=True,
                    module="Metal Site Parameterizer",
                    description="Add implicit solvation for anion (small model)",
                ):
                    self.console.print("\n[cyan]Common solvents:[/cyan]")
                    self.console.print("  Water (ε=78.4), DiMethylSulfoxide (ε=46.8), Acetonitrile (ε=35.7),")
                    self.console.print("  Methanol (ε=32.6), Ethanol (ε=24.9), Dichloromethane (ε=8.9),")
                    self.console.print("  TetraHydroFuran (ε=7.4), Chloroform (ε=4.7), Toluene (ε=2.4)")
                    self.console.print("[grey50]Full list: see Gaussian SCRF documentation[/grey50]")

                    solvent_name = prompt_with_context(
                        self.processor, "Solvent name (Gaussian keyword)",
                        default="Water",
                        module="Metal Site Parameterizer",
                        description="SCRF solvent selection (small model)",
                    )
                    small_scrf = f"SCRF=(Solvent={solvent_name})"
                    self.console.print(f"[green]Will add {small_scrf} to route[/green]")

            # Method/route customization - panel above already shows defaults
            if confirm_with_context(
                self.processor,
                "Use recommended method and route (B3LYP/6-31G*)?",
                default=True,
                module="Metal Site Parameterizer",
                description="Use recommended small model method"
            ):
                small_functional = "B3LYP"
                small_basis = "6-31G*"
                small_route_line = None
                small_additional = small_scrf
            else:
                small_functional = prompt_with_context(
                    self.processor,
                    "DFT functional",
                    default="B3LYP",
                    module="Metal Site Parameterizer",
                    description="Small model functional"
                )

                small_basis = prompt_with_context(
                    self.processor,
                    "Basis set",
                    default="6-31G*",
                    module="Metal Site Parameterizer",
                    description="Small model basis set"
                )

                self.console.print("\n[yellow]Custom route line[/yellow]")
                self.console.print("[grey50]Include all keywords. IOp(7/33=1) REQUIRED for Seminario method.[/grey50]")
                scrf_part = f" {small_scrf}" if small_scrf else ""
                default_route = f"# {small_functional}/{small_basis} Opt Freq Geom=PrintInputOrient Integral=(Grid=UltraFine) IOp(7/33=1){scrf_part}"
                small_route_line = prompt_with_context(
                    self.processor,
                    "Full route line",
                    default=default_route,
                    module="Metal Site Parameterizer",
                    description="Small model route line"
                )
                small_additional = ""  # SCRF already in custom route line

            # Frozen atoms option for geometry optimization
            self.console.print("\n[bold]Atom Freezing Options[/bold]")
            self.console.print("[grey50]Choose which atoms to freeze during geometry optimization.[/grey50]")
            self.console.print("[grey50]Frozen atoms retain their original coordinates; only free atoms move.[/grey50]\n")
            self.console.print("  [cyan]1[/cyan] Full optimization (no frozen atoms)")
            self.console.print("  [cyan]2[/cyan] Freeze heavy atoms, optimize hydrogens only")
            self.console.print("  [cyan]3[/cyan] Freeze backbone (CA, N, C, O), optimize side chains + H")
            self.console.print("  [cyan]4[/cyan] Custom selection (toggle individual atoms)")
            freeze_choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3", "4"],
                default="1",
                module="Metal Site Parameterizer",
                description="Atom freezing option"
            )

            # Parse PDB to get atom info
            frozen_atoms = None
            atom_info_list = []  # List of (element, x, y, z, atom_name, resname, resid)
            with open(small_pdb) as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        element = line[76:78].strip()
                        if not element:
                            atom_name = line[12:16].strip()
                            element = ''.join([c for c in atom_name if c.isalpha()])[:2]
                            if len(element) == 2:
                                element = element[0].upper() + element[1:].lower()
                            else:
                                element = element.upper()
                        else:
                            atom_name = line[12:16].strip()
                        resname = line[17:20].strip()
                        resid = line[22:26].strip()
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        atom_info_list.append((element, x, y, z, atom_name, resname, resid))

            if freeze_choice != "1":
                frozen_atoms = set()
                if freeze_choice == "2":
                    # Freeze all non-hydrogen atoms
                    for i, (element, _, _, _, _, _, _) in enumerate(atom_info_list):
                        if element.upper() != "H":
                            frozen_atoms.add(i)
                    self.console.print(f"[grey50]  Freezing {len(frozen_atoms)} heavy atoms[/grey50]")
                elif freeze_choice == "3":
                    # Freeze backbone atoms
                    backbone_names = {"CA", "N", "C", "O"}
                    for i, (_, _, _, _, atom_name, _, _) in enumerate(atom_info_list):
                        if atom_name in backbone_names:
                            frozen_atoms.add(i)
                    self.console.print(f"[grey50]  Freezing {len(frozen_atoms)} backbone atoms[/grey50]")
                elif freeze_choice == "4":
                    # Custom selection - show table and let user toggle
                    self.console.print("\n[bold]Custom Atom Selection[/bold]")
                    self.console.print("[grey50]Enter atom numbers to toggle freeze status (space-separated).[/grey50]")
                    self.console.print("[grey50]All atoms start as FREE. Example: '1 5 12' freezes atoms 1, 5, 12.[/grey50]\n")

                    # Display atom table
                    from rich.table import Table
                    table = Table(title="Atoms in Small Model", show_lines=False)
                    table.add_column("#", style="cyan", width=4)
                    table.add_column("Elem", width=4)
                    table.add_column("Name", width=6)
                    table.add_column("Res", width=5)
                    table.add_column("ResID", width=5)
                    table.add_column("Status", width=8)

                    for i, (elem, _, _, _, aname, rname, rid) in enumerate(atom_info_list):
                        status = "[red]FROZEN[/red]" if i in frozen_atoms else "[green]free[/green]"
                        table.add_row(str(i+1), elem, aname, rname, rid, status)

                    self.console.print(table)

                    while True:
                        toggle_input = prompt_with_context(
                            self.processor,
                            "Atoms to toggle (or 'done')",
                            default="done",
                            module="Metal Site Parameterizer",
                            description="Toggle frozen atoms"
                        )
                        if toggle_input.lower() == "done":
                            break
                        try:
                            indices = [int(x) - 1 for x in toggle_input.split()]
                            for idx in indices:
                                if 0 <= idx < len(atom_info_list):
                                    if idx in frozen_atoms:
                                        frozen_atoms.remove(idx)
                                    else:
                                        frozen_atoms.add(idx)
                            # Redisplay table
                            table = Table(title="Atoms in Small Model", show_lines=False)
                            table.add_column("#", style="cyan", width=4)
                            table.add_column("Elem", width=4)
                            table.add_column("Name", width=6)
                            table.add_column("Res", width=5)
                            table.add_column("ResID", width=5)
                            table.add_column("Status", width=8)
                            for i, (elem, _, _, _, aname, rname, rid) in enumerate(atom_info_list):
                                status = "[red]FROZEN[/red]" if i in frozen_atoms else "[green]free[/green]"
                                table.add_row(str(i+1), elem, aname, rname, rid, status)
                            self.console.print(table)
                        except ValueError:
                            self.console.print("[yellow]Invalid input. Enter space-separated numbers or 'done'.[/yellow]")

                    self.console.print(f"[grey50]  Freezing {len(frozen_atoms)} atoms[/grey50]")

                # If no atoms selected, set to None
                if not frozen_atoms:
                    frozen_atoms = None

            # Generate small model input
            qm_interface = QMInterface(QMSoftware.GAUSSIAN, QMCalculationMode.GUIDED, self.console)

            small_gjf_files = qm_interface.generate_input_files(
                pdb_file=small_pdb,
                output_dir=step1_dir,
                charge=small_charge,
                multiplicity=small_mult,
                memory_gb=small_memory_gb,
                n_processors=small_nproc,
                functional=small_functional,
                basis_set=small_basis,
                job_type="Opt Freq",
                additional_keywords=small_additional,
                title_card="MCPB Small Model - Optimization and Frequency",
                output_name="small_freq",
                route_line=small_route_line,
                include_seminario_iop=True,
                frozen_atoms=frozen_atoms
            )
            small_gjf = small_gjf_files["input_file"]

            # ----------------------------------------------------------------
            # Large Model: ESP Calculation (for RESP charges)
            # ----------------------------------------------------------------
            self.console.print("\n[bold]Large Model: ESP Calculation[/bold]")
            self.console.print("[grey50]Used for RESP charge fitting[/grey50]\n")

            if confirm_with_context(
                self.processor,
                f"Use same resources as small model ({small_memory_gb}GB, {small_nproc} proc)?",
                default=True,
                module="Metal Site Parameterizer",
                description="Use same resources"
            ):
                large_memory_gb = small_memory_gb
                large_nproc = small_nproc
            else:
                large_memory_gb = int_prompt_with_context(
                    self.processor,
                    "Memory allocation (GB)",
                    default=small_memory_gb,
                    module="Metal Site Parameterizer",
                    description="Large model memory"
                )
                large_nproc = int_prompt_with_context(
                    self.processor,
                    "Number of processors",
                    default=small_nproc,
                    module="Metal Site Parameterizer",
                    description="Large model processors"
                )

            self.console.print(f"\n[grey50]Suggested charge from PDB: {large_suggested_charge}[/grey50]")
            self._print_withheld_cluster_charge_note(redox_site, type_assignments)
            large_charge = int_prompt_with_context(
                self.processor,
                "Total charge of large model",
                default=large_suggested_charge,
                module="Metal Site Parameterizer",
                description="Large model charge"
            )

            large_mult = int_prompt_with_context(
                self.processor,
                "Spin multiplicity",
                default=small_mult,
                module="Metal Site Parameterizer",
                description="Large model multiplicity"
            )

            # Implicit solvation for the large model.
            #
            # Asked outright rather than inherited. The two models are separate
            # calculations, and copying the small model's SCRF silently meant
            # the large model's solvation was never a visible decision — the
            # run only announced it after the fact.
            #
            # Also asked regardless of charge. The whole block used to sit
            # under `large_charge < 0`, so a neutral or cationic large model
            # dropped the small model's SCRF without comment, leaving the two
            # calculations at different levels of theory for no stated reason.
            large_scrf = ""
            if large_charge < 0 and not small_scrf:
                self.console.print(Panel(
                    f"[bold yellow]Anion detected (charge = {large_charge}):[/bold yellow]\n\n"
                    "Implicit solvation improves SCF convergence and the charge\n"
                    "distribution for anionic systems, which is what the ESP is\n"
                    "fitted to here.",
                    title="Anion Solvation Warning",
                    border_style="yellow",
                    expand=False,
                ))

            if small_scrf:
                self.console.print(
                    f"[grey50]Small model uses {small_scrf}.[/grey50]")

            # Default: match the small model when it used solvation, else
            # recommend it only for an anion.
            if confirm_with_context(
                self.processor,
                "Add implicit solvation (SCRF) to the large model?",
                default=self._default_large_model_solvation(small_scrf, large_charge),
                module="Metal Site Parameterizer",
                description="Add implicit solvation (large model)",
            ):
                use_same = False
                if small_scrf:
                    use_same = confirm_with_context(
                        self.processor,
                        f"Use the same solvent as the small model ({small_scrf})?",
                        default=True,
                        module="Metal Site Parameterizer",
                        description="Match the small model's solvent (large model)",
                    )

                if use_same:
                    large_scrf = small_scrf
                else:
                    self.console.print("\n[cyan]Common solvents:[/cyan]")
                    self.console.print("  Water (ε=78.4), DiMethylSulfoxide (ε=46.8), Acetonitrile (ε=35.7),")
                    self.console.print("  Methanol (ε=32.6), Ethanol (ε=24.9), Dichloromethane (ε=8.9),")
                    self.console.print("  TetraHydroFuran (ε=7.4), Chloroform (ε=4.7), Toluene (ε=2.4)")
                    self.console.print("[grey50]Full list: see Gaussian SCRF documentation[/grey50]")

                    solvent_name = prompt_with_context(
                        self.processor, "Solvent name (Gaussian keyword)",
                        default="Water",
                        module="Metal Site Parameterizer",
                        description="SCRF solvent selection (large model)",
                    )
                    large_scrf = f"SCRF=(Solvent={solvent_name})"

                self.console.print(f"[green]Will add {large_scrf} to route[/green]")
            elif small_scrf:
                # Declining after the small model used solvation is a real
                # divergence between the two calculations; say so once.
                self.console.print(
                    "[yellow]Large model will run without solvation while the "
                    "small model used it — the two models will be at different "
                    "levels of theory.[/yellow]")

            # Default to the small model's level of theory. The MCPB.py tutorial
            # (ambermd.org/tutorials/advanced/tutorial20) states it plainly:
            # "we used the B3LYP/6-31G* level of theory to perform the
            # calculations for both the small and large models". The previous
            # HF/6-31G* default was the generic RESP convention for organics,
            # not what MCPB.py does for a metal site — and it silently differed
            # from the functional the user had just chosen for the small model.
            self.console.print(
                f"[grey50]Small model uses {small_functional}/{small_basis} — "
                f"MCPB.py uses one level of theory for both models.[/grey50]"
            )
            if confirm_with_context(
                self.processor,
                f"Use the same method as the small model "
                f"({small_functional}/{small_basis})?",
                default=True,
                module="Metal Site Parameterizer",
                description="Use recommended large model method"
            ):
                large_functional = small_functional
                large_basis = small_basis
            else:
                large_functional = prompt_with_context(
                    self.processor,
                    "Method for ESP",
                    default=small_functional,
                    module="Metal Site Parameterizer",
                    description="Large model method"
                )

                large_basis = prompt_with_context(
                    self.processor,
                    "Basis set",
                    default=small_basis,
                    module="Metal Site Parameterizer",
                    description="Large model basis set"
                )

            metal_radii = self._collect_metal_radii(type_assignments)
            # Frozen atoms option for geometry optimization (same as small model)
            self.console.print("\n[bold]Atom Freezing Options[/bold]")
            self.console.print("[grey50]Choose which atoms to freeze during geometry optimization.[/grey50]")
            self.console.print("[grey50]Frozen atoms retain their original coordinates; only free atoms move.[/grey50]")
            self.console.print("[grey50]Note: The large model has more atoms, so optimization will take longer.[/grey50]\n")
            self.console.print("  [cyan]1[/cyan] Full optimization (no frozen atoms)")
            self.console.print("  [cyan]2[/cyan] Freeze heavy atoms, optimize hydrogens only")
            self.console.print("  [cyan]3[/cyan] Freeze backbone (CA, N, C, O), optimize side chains + H")
            self.console.print("  [cyan]4[/cyan] Custom selection (toggle individual atoms)")
            large_freeze_choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=["1", "2", "3", "4"],
                default="1",
                module="Metal Site Parameterizer",
                description="Large model atom freezing option"
            )

            # Parse large model PDB to get atom info
            large_frozen_atoms = None
            large_atom_info_list = []  # List of (element, x, y, z, atom_name, resname, resid)
            with open(large_pdb) as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        element = line[76:78].strip()
                        if not element:
                            atom_name = line[12:16].strip()
                            element = ''.join([c for c in atom_name if c.isalpha()])[:2]
                            if len(element) == 2:
                                element = element[0].upper() + element[1:].lower()
                            else:
                                element = element.upper()
                        else:
                            atom_name = line[12:16].strip()
                        resname = line[17:20].strip()
                        resid = line[22:26].strip()
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        large_atom_info_list.append((element, x, y, z, atom_name, resname, resid))

            if large_freeze_choice != "1":
                large_frozen_atoms = set()
                if large_freeze_choice == "2":
                    # Freeze all non-hydrogen atoms
                    for i, (element, _, _, _, _, _, _) in enumerate(large_atom_info_list):
                        if element.upper() != "H":
                            large_frozen_atoms.add(i)
                    self.console.print(f"[grey50]  Freezing {len(large_frozen_atoms)} heavy atoms[/grey50]")
                elif large_freeze_choice == "3":
                    # Freeze backbone atoms
                    backbone_names = {"CA", "N", "C", "O"}
                    for i, (_, _, _, _, atom_name, _, _) in enumerate(large_atom_info_list):
                        if atom_name in backbone_names:
                            large_frozen_atoms.add(i)
                    self.console.print(f"[grey50]  Freezing {len(large_frozen_atoms)} backbone atoms[/grey50]")
                elif large_freeze_choice == "4":
                    # Custom selection - show table and let user toggle
                    self.console.print("\n[bold]Custom Atom Selection[/bold]")
                    self.console.print("[grey50]Enter atom numbers to toggle freeze status (space-separated).[/grey50]")
                    self.console.print("[grey50]All atoms start as FREE. Example: '1 5 12' freezes atoms 1, 5, 12.[/grey50]\n")

                    # Display atom table
                    from rich.table import Table
                    table = Table(title="Atoms in Large Model", show_lines=False)
                    table.add_column("#", style="cyan", width=4)
                    table.add_column("Elem", width=4)
                    table.add_column("Name", width=6)
                    table.add_column("Res", width=5)
                    table.add_column("ResID", width=5)
                    table.add_column("Status", width=8)

                    for i, (elem, _, _, _, aname, rname, rid) in enumerate(large_atom_info_list):
                        status = "[red]FROZEN[/red]" if i in large_frozen_atoms else "[green]free[/green]"
                        table.add_row(str(i+1), elem, aname, rname, rid, status)

                    self.console.print(table)

                    while True:
                        toggle_input = prompt_with_context(
                            self.processor,
                            "Atoms to toggle (or 'done')",
                            default="done",
                            module="Metal Site Parameterizer",
                            description="Toggle frozen atoms (large model)"
                        )
                        if toggle_input.lower() == "done":
                            break
                        try:
                            indices = [int(x) - 1 for x in toggle_input.split()]
                            for idx in indices:
                                if 0 <= idx < len(large_atom_info_list):
                                    if idx in large_frozen_atoms:
                                        large_frozen_atoms.remove(idx)
                                    else:
                                        large_frozen_atoms.add(idx)
                            # Redisplay table
                            table = Table(title="Atoms in Large Model", show_lines=False)
                            table.add_column("#", style="cyan", width=4)
                            table.add_column("Elem", width=4)
                            table.add_column("Name", width=6)
                            table.add_column("Res", width=5)
                            table.add_column("ResID", width=5)
                            table.add_column("Status", width=8)
                            for i, (elem, _, _, _, aname, rname, rid) in enumerate(large_atom_info_list):
                                status = "[red]FROZEN[/red]" if i in large_frozen_atoms else "[green]free[/green]"
                                table.add_row(str(i+1), elem, aname, rname, rid, status)
                            self.console.print(table)
                        except ValueError:
                            self.console.print("[yellow]Invalid input. Enter space-separated numbers or 'done'.[/yellow]")

                    self.console.print(f"[grey50]  Freezing {len(large_frozen_atoms)} atoms[/grey50]")

                # If no atoms selected, set to None
                if not large_frozen_atoms:
                    large_frozen_atoms = None

            # Job type: geometry optimization. The Merz-Kollman ESP (Pop=MK) is
            # then evaluated at the optimized geometry in the same job, so a
            # separate SP keyword is not needed — and "Opt SP" is a contradictory
            # route (Opt and SP are mutually exclusive job types).
            large_job_type = "Opt"

            # Route line customization for ESP
            # Use ReadRadii when we have metal radii to specify
            pop_keyword = "Pop(MK,ReadRadii)" if metal_radii else "Pop=MK"
            scrf_part = f" {large_scrf}" if large_scrf else ""
            default_esp_route = f"# {large_functional}/{large_basis} Opt {pop_keyword} IOp(6/33=2,6/41=10,6/42=17){scrf_part}"

            from rich.panel import Panel
            self.console.print(Panel(
                f"[bold cyan]ESP Calculation Keywords[/bold cyan]\n\n"
                f"Route: [grey50]{default_esp_route}[/grey50]\n\n"
                f"[bold]Keyword explanations:[/bold]\n"
                f"  • [bold]Opt[/bold]             Geometry optimization; the ESP is evaluated at the optimized geometry\n"
                f"  • [bold]{pop_keyword}[/bold]   Merz-Kollman ESP with custom metal radii\n"
                f"  • [bold]IOp(6/33=2)[/bold]     Write ESP points and potentials to output file\n"
                f"  • [bold]IOp(6/41=10)[/bold]    10 concentric layers of ESP points per atom\n"
                f"  • [bold]IOp(6/42=17)[/bold]    ~1700 points per unit area (high density grid)\n\n"
                f"[grey50]The ESP grid is used by RESP to fit atomic partial charges.[/grey50]",
                title="Large Model ESP Route",
                border_style="cyan",
                expand=False
            ))

            # If user already customized method, offer route customization too
            if large_functional != "HF" or large_basis != "6-31G*":
                if not confirm_with_context(
                    self.processor,
                    "Use standard ESP route keywords?",
                    default=True,
                    module="Metal Site Parameterizer",
                    description="Use standard ESP route"
                ):
                    self.console.print("\n[yellow]Custom ESP route line[/yellow]")
                    self.console.print("[grey50]Include all keywords. Must include ESP output options for RESP.[/grey50]")
                    large_route_line = prompt_with_context(
                        self.processor,
                        "Full route line",
                        default=default_esp_route,
                        module="Metal Site Parameterizer",
                        description="Large model route line"
                    )
                    large_esp_keywords = ""
                else:
                    large_route_line = None
                    large_esp_keywords = f"{pop_keyword} IOp(6/33=2,6/41=10,6/42=17){scrf_part}"
            else:
                large_route_line = None
                large_esp_keywords = f"{pop_keyword} IOp(6/33=2,6/41=10,6/42=17){scrf_part}"

            large_gjf_files = qm_interface.generate_input_files(
                pdb_file=large_pdb,
                output_dir=step1_dir,
                charge=large_charge,
                multiplicity=large_mult,
                memory_gb=large_memory_gb,
                n_processors=large_nproc,
                functional=large_functional,
                basis_set=large_basis,
                job_type=large_job_type,
                additional_keywords=large_esp_keywords,
                title_card="MCPB Large Model - Opt + ESP for RESP charges",
                output_name="large_resp",
                route_line=large_route_line,
                include_seminario_iop=False,
                metal_radii=metal_radii if metal_radii else None,
                frozen_atoms=large_frozen_atoms
            )
            large_gjf = large_gjf_files["input_file"]

            # Summary
            self.console.print("\n[bold yellow]Next Steps for QM Calculations:[/bold yellow]")
            self.console.print(f"  1. Run Gaussian on small model:")
            self.console.print(f"     g16 {small_gjf.name}")
            self.console.print(f"     formchk {small_gjf.stem}.chk {small_gjf.stem}.fchk")
            self.console.print(f"  2. Run Gaussian on large model:")
            self.console.print(f"     g16 {large_gjf.name}")
            self.console.print("  3. Return to ProPrep and resume the checklist")

            # ================================================================
            # 12. Save results
            # ================================================================
            # Save atom type assignments to JSON
            assignments_file = step1_dir / "atom_type_assignments.json"
            self._save_atom_type_assignments(type_assignments, assignments_file)
            self.console.print(f"[grey50]Saved atom assignments: {assignments_file.name}[/grey50]")

            self.console.print(f"\n[green]✅ Atom typing & fingerprint generation complete[/green]")

            # Store and persist results
            step1_result = {
                "success": True,
                "step_number": 1,
                "step_description": "MCPB Atom Typing & Fingerprint Generation",
                "next_step": 2,
                "output_files": {
                    "standard_fingerprint": str(standard_fp),
                    "large_fingerprint": str(large_fp),
                    "atom_assignments": str(assignments_file),
                    "small_pdb": str(small_pdb),
                    "large_pdb": str(large_pdb),
                    "small_gjf": str(small_gjf),
                    "large_gjf": str(large_gjf)
                },
                "force_field_info": {
                    "force_field": self._extract_force_field_from_leaprc(ff_data.loaded_leaprcs) if ff_data else "prmtop",
                    "loaded_leaprcs": ff_data.loaded_leaprcs if ff_data else [],
                    "ff_data": ff_data  # May be None if using prmtop directly
                },
                # Store QM parameters for step 3 (RESP)
                "qm_parameters": {
                    "small_model": {
                        "charge": small_charge,
                        "multiplicity": small_mult,
                        "functional": small_functional,
                        "basis_set": small_basis
                    },
                    "large_model": {
                        "charge": large_charge,
                        "multiplicity": large_mult,
                        "functional": large_functional,
                        "basis_set": large_basis
                    }
                },
                # Store type_assignments for Step 3B RESP backbone charges
                "type_assignments": type_assignments,
                "atom_summary": {
                    "total_atoms": len(type_assignments),
                    "metal_atoms": metal_count,
                    "ligand_atoms": ligand_count,
                    "renamed_atoms": renamed_count
                },
                "statistics": {
                    "total_atoms": len(type_assignments),
                    "metal_atoms": metal_count,
                    "ligand_atoms": ligand_count,
                    "renamed_atoms": renamed_count,
                    "atom_lines": fp_validation['atom_lines'],
                    "link_lines": fp_validation['link_lines']
                },
                "fingerprint_validation": fp_validation,
                "redox_site_info": {
                    "total_atoms": len(redox_site.atoms),
                    "metal_centers": len(redox_site.centers),
                    "bonds": len(redox_site.bonds)
                }
            }
            self.step_results["step_1"] = step1_result
            self._save_step_results()
            return step1_result

        except Exception as e:
            import traceback
            traceback.print_exc()

            return {
                "success": False,
                "message": f"Step 1 error: {str(e)}",
                "traceback": traceback.format_exc()
            }

    def _assign_atom_types(self, redox_site: RedoxSite, ff_data: ForceFieldData,
                          terminal_classifier: TerminalClassifier, output_dir: Path) -> Dict:
        """
        Assign atom types using THREE-TIER approach:
        1. Metals: Auto-configured from MetalIonDatabase
        2. Organics: Standard force field lookup
        3. Inorganic Ligands: Antechamber parameterization with user grouping

        Always tracks element information separately from atom type.

        Args:
            redox_site: RedoxSite object with atoms
            ff_data: Force field data collector with atom type definitions
            terminal_classifier: Terminal classifier for getting terminal-aware residue names
            output_dir: Output directory for Step 1 files

        Returns:
            Dictionary mapping atom coordinates to type assignment data
        """
        type_assignments = {}

        # Build metal_coords from centers (handles both METAL_ION and ORGANOMETALLIC_COFACTOR)
        # - METAL_ION: center.coords IS the metal
        # - ORGANOMETALLIC_COFACTOR: find metal atoms in the same residue
        metal_coords = set()
        for center in redox_site.centers:
            if hasattr(center, 'center_type'):
                if center.center_type.value == 'metal_ion':
                    # Isolated metal - center coords ARE the metal
                    metal_coords.add(center.coords)
                elif center.center_type.value == 'organometallic_cofactor':
                    # Embedded metal - find metal atom(s) in this residue
                    for atom in redox_site.atoms:
                        if (atom.chain == center.chain and
                            atom.resid == center.resid and
                            atom.resname == center.resname and
                            atom.element.upper() in METALS):
                            metal_coords.add(atom.coords)

        # Also keep center_coords for is_center check (for whole-residue centers)
        center_coords = metal_coords.copy()

        # Build set of ligand atom coordinates (atoms bonded to metals)
        ligand_coords = set()
        for bond in redox_site.bonds:
            # Restrained contacts (e.g. a nonbonded water held by an MD restraint)
            # are not MCPB bonded ligands, so they must NOT be Y-renamed.
            if bond.chemical_type == 'coordinate' and getattr(bond, 'treatment', 'bonded') == 'bonded':
                # Check which end is the metal
                if bond.atom1_coords in metal_coords:
                    ligand_coords.add(bond.atom2_coords)
                elif bond.atom2_coords in metal_coords:
                    ligand_coords.add(bond.atom1_coords)

        # Pure-cluster internal atoms (bridging sulfides, Mo core O/S) coordinate
        # the metal within the residue and get no 'coordinate' bond above; type
        # them as ligands too so they receive unique Y* types (see
        # _cluster_internal_ligand_coords).
        ligand_coords |= self._cluster_internal_ligand_coords(redox_site, metal_coords)

        # Identify which atoms belong to each tier
        metal_atoms = []       # Tier 1: Metals from MetalIonDatabase
        organic_atoms = []     # Tier 2: Standard force field
        inorganic_ligands = []  # Tier 3: Antechamber

        # Initialize MetalIonDatabase
        metal_db = MetalIonDatabase(water_model=self._water_model(), logger=self.logger)

        # Classify atoms into tiers
        for atom in redox_site.atoms:
            # Tier 1: Check if metal by resname/atomname OR by element
            # This handles both standard metal residues (FE, ZN) and complex residues (SF4)
            is_metal_residue = metal_db.is_metal(atom.resname, atom.atom_name)
            is_metal_elem = metal_db.is_metal_element(atom.element)

            if is_metal_residue or is_metal_elem:
                metal_atoms.append(atom)
            else:
                # Check if in force field
                res_name = atom.resname
                atom_name = atom.atom_name
                terminal_type, _, _ = terminal_classifier.classify_residue(atom.chain, atom.resid, res_name)
                ff_res_name = terminal_classifier.get_force_field_residue_name(res_name, terminal_type)
                atom_def = ff_data.get_atom_definition(ff_res_name, atom_name)

                if atom_def:
                    # Tier 2: Standard force field
                    organic_atoms.append((atom, terminal_type, ff_res_name, atom_def))
                else:
                    # Check if residue exists at all - might need mapping (e.g., HEM → HEH)
                    if not ff_data.has_residue(ff_res_name):
                        # Residue not found - check if we already have an alias
                        alias = ff_data.residue_aliases.get(ff_res_name)
                        if alias:
                            # Try with alias
                            atom_def = ff_data.get_atom_definition(alias, atom_name)
                            if atom_def:
                                organic_atoms.append((atom, terminal_type, ff_res_name, atom_def))
                                continue

                        # Still not found - add to list for mapping prompt later
                        inorganic_ligands.append(atom)
                    else:
                        # Residue exists but atom not found (weird atom name?)
                        inorganic_ligands.append(atom)

        # ================================================================
        # TIER 1: Process Metals (Auto-configuration)
        # ================================================================
        self.console.print(f"\n[cyan]Tier 1: Processing {len(metal_atoms)} metal atoms...[/cyan]")

        for atom in metal_atoms:
            # Always ask user for metal charge and spin (cannot reliably infer)
            self.console.print(f"\n[yellow]Metal: {atom.resname} {atom.atom_name} (Element: {atom.element})[/yellow]")

            while True:
                try:
                    charge_input = prompt_with_context(
                        self.processor, "Enter formal charge (e.g., +2, +3, -1)",
                        module="Metal Site Parameterizer",
                        description="Metal formal charge",
                    ).strip()
                    charge = int(charge_input)
                    break
                except ValueError:
                    self.console.print(f"[red]Invalid charge: {charge_input}. Please enter an integer.[/red]")

            while True:
                try:
                    spin_input = prompt_with_context(
                        self.processor, "Enter number of unpaired electrons (e.g., 0, 1, 5, or -1 for beta)",
                        module="Metal Site Parameterizer",
                        description="Metal unpaired electrons",
                    ).strip()
                    spin = int(spin_input)
                    break
                except ValueError:
                    self.console.print(f"[red]Invalid spin: {spin_input}. Please enter an integer.[/red]")

            # Try direct lookup first (for standard metal residues like FE, ZN)
            metal_config = metal_db.get_metal_config(atom.resname, atom.atom_name, charge, spin)

            if not metal_config:
                # Try element-based lookup (for complex residues like SF4)
                metal_config = metal_db.get_metal_config_by_element(
                    atom.element, atom.atom_name, atom.resname, charge, spin
                )

            if not metal_config:
                # Fallback: manual entry (metal not in database)
                self.console.print(f"[yellow]⚠ Metal {atom.resname}:{atom.atom_name} not in database[/yellow]")
                metal_config = self._manual_metal_entry(atom, charge, spin)

            coords = atom.coords
            is_center = coords in center_coords

            type_assignments[coords] = {
                'chain': atom.chain,
                'resname': atom.resname,
                'ff_resname': atom.resname,  # Metals don't have terminal variants
                'resid': atom.resid,
                'atom_name': atom.atom_name,
                'element': metal_config.element,  # ALWAYS track element
                'original_type': metal_config.atom_type,
                'renamed_type': metal_config.atom_type,
                'charge': float(metal_config.charge),  # Store user-provided charge
                'spin': metal_config.spin,  # Store user-provided spin
                'mass': metal_config.mass,  # Store mass for MASS section
                'terminal_type': 'middle',
                'is_center': is_center,
                'is_metal_ligand': False,  # Metals aren't ligands
                'renamed': False,
                'source': 'metal_database',
                'vdw_params': {
                    'radius': metal_config.vdw_radius,
                    'epsilon': metal_config.vdw_epsilon,
                    'source': metal_config.vdw_source
                }
            }

        # ================================================================
        # TIER 2: Process Organic Atoms (Standard Force Field)
        # ================================================================
        self.console.print(f"[cyan]Tier 2: Processing {len(organic_atoms)} organic atoms...[/cyan]")

        for atom, terminal_type, ff_res_name, atom_def in organic_atoms:
            coords = atom.coords
            is_center = coords in center_coords
            is_ligand = coords in ligand_coords

            type_assignments[coords] = {
                'chain': atom.chain,
                'resname': atom.resname,
                'ff_resname': ff_res_name,
                'resid': atom.resid,
                'atom_name': atom.atom_name,
                'element': atom.element,  # From RedoxSiteAtom
                'original_type': atom_def.atom_type,
                'renamed_type': atom_def.atom_type,
                'charge': atom_def.atom_charge,
                'terminal_type': terminal_type.value,
                'is_center': is_center,
                'is_metal_ligand': is_ligand,
                'renamed': False,
                'source': 'force_field'
            }

        # ================================================================
        # CHECK FOR UNKNOWN RESIDUES - Prompt for mapping before Tier 3
        # ================================================================
        if inorganic_ligands:
            # Find unique residues that are NOT in the force field
            unknown_residues = set()
            for atom in inorganic_ligands:
                if not ff_data.has_residue(atom.resname):
                    unknown_residues.add(atom.resname)

            # Prompt for residue mapping if there are unknown residues
            if unknown_residues:
                self.console.print(f"\n[yellow]Found {len(unknown_residues)} residue(s) not in force field: {list(unknown_residues)}[/yellow]")
                self.console.print("[grey50]These may need to be mapped to existing residues (e.g., HEM → HEH)[/grey50]")

                mappings = ff_data.prompt_residue_mapping(list(unknown_residues))

                if mappings:
                    # Re-classify atoms from mapped residues
                    remaining_inorganic = []
                    for atom in inorganic_ligands:
                        if atom.resname in mappings:
                            # Try to get atom type from mapped residue
                            mapped_res = mappings[atom.resname]
                            atom_def = ff_data.get_atom_definition(atom.resname, atom.atom_name)
                            if atom_def:
                                terminal_type, _, _ = terminal_classifier.classify_residue(atom.chain, atom.resid, atom.resname)
                                organic_atoms.append((atom, terminal_type, atom.resname, atom_def))
                                self.console.print(f"  [green]✓ {atom.resname}-{atom.atom_name} → mapped to {mapped_res}[/green]")
                            else:
                                # Atom name not found even after mapping
                                remaining_inorganic.append(atom)
                                self.console.print(f"  [yellow]⚠ {atom.resname}-{atom.atom_name} not found in {mapped_res}[/yellow]")
                        else:
                            remaining_inorganic.append(atom)

                    inorganic_ligands = remaining_inorganic
                    self.console.print(f"\n[cyan]After mapping: {len(organic_atoms)} organic, {len(inorganic_ligands)} remaining inorganic[/cyan]")

        # ================================================================
        # TIER 3: Process Inorganic Ligands (Antechamber)
        # ================================================================
        if inorganic_ligands:
            self.console.print(f"\n[cyan]Tier 3: Processing {len(inorganic_ligands)} inorganic ligand atoms...[/cyan]")

            # Use LigandGroupingInterface for user-controlled grouping
            grouping_interface = LigandGroupingInterface(console=self.console, logger=self.logger, processor=self.processor)

            # Display table and get user grouping (ONLY for Tier 3 inorganic ligands)
            idx_to_atom = grouping_interface.display_atoms(inorganic_ligands)

            if idx_to_atom:
                # Prompt user for grouping
                ligand_groups = grouping_interface.prompt_user_grouping(idx_to_atom, redox_site)

                # Run antechamber for each group
                antechamber = AntechamberRunner(gaff_version='gaff2', logger=self.logger, processor=self.processor)
                antechamber_output_dir = Path(output_dir) / "antechamber_ligands"

                for group in ligand_groups:
                    self.console.print(f"\n[cyan]Processing {group.group_id}...[/cyan]")

                    # Handle based on parameterization method
                    if group.parameterization_method == 'antechamber':
                        # Special case: single atoms cannot be parameterized with antechamber
                        # (AM1-BCC requires molecular structure for QM calculation)
                        if len(group.atoms) == 1:
                            atom = group.atoms[0]
                            self.console.print(f"[yellow]⚠ Single atom detected. Antechamber requires ≥2 atoms.[/yellow]")
                            self.console.print(f"[cyan]Querying available {atom.element} atom types...[/cyan]\n")

                            # Query force field for this element (includes GAFF2 by default)
                            available_dict = ff_data.get_atom_types_for_element(atom.element, include_gaff2=True)

                            from rich.table import Table

                            # Display force field results
                            ff_results = available_dict.get('force_field', [])
                            if ff_results:
                                ff_table = Table(title=f"{atom.element} Atom Types from Selected Force Field")
                                ff_table.add_column("Residue", style="cyan", no_wrap=True)
                                ff_table.add_column("Atom Name", style="yellow", no_wrap=True)
                                ff_table.add_column("Atom Type", style="green", no_wrap=True)
                                ff_table.add_column("Charge", style="magenta", justify="right")
                                ff_table.add_column("Description", style="grey50")

                                for res, aname, atype, charge, desc in ff_results[:15]:  # Show max 15
                                    ff_table.add_row(res, aname, atype, f"{charge:.4f}", desc)

                                self.console.print(ff_table)
                                if len(ff_results) > 15:
                                    self.console.print(f"[grey50]... and {len(ff_results)-15} more[/grey50]")
                                self.console.print()

                            # Display GAFF2 results
                            gaff2_results = available_dict.get('gaff2', [])
                            if gaff2_results:
                                gaff2_table = Table(title=f"{atom.element} Atom Types from GAFF2 (General Ligands)")
                                gaff2_table.add_column("Residue", style="cyan", no_wrap=True)
                                gaff2_table.add_column("Atom Name", style="yellow", no_wrap=True)
                                gaff2_table.add_column("Atom Type", style="green", no_wrap=True)
                                gaff2_table.add_column("Charge", style="magenta", justify="right")
                                gaff2_table.add_column("Description", style="grey50")

                                for res, aname, atype, charge, desc in gaff2_results[:15]:  # Show max 15
                                    gaff2_table.add_row(res, aname, atype, f"{charge:.4f}", desc)

                                self.console.print(gaff2_table)
                                if len(gaff2_results) > 15:
                                    self.console.print(f"[grey50]... and {len(gaff2_results)-15} more[/grey50]")
                                self.console.print()

                            if not ff_results and not gaff2_results:
                                self.console.print(f"[yellow]No {atom.element} atom types found.[/yellow]\n")

                            # Fall back to manual entry with atom type suggestions
                            params = antechamber.manual_entry_fallback(group.atoms, group.group_id)
                        else:
                            # Run antechamber for multi-atom groups
                            self.console.print(f"[grey50]Running antechamber with charge={group.formal_charge:+d}...[/grey50]")
                            params = antechamber.parameterize_ligand_group(
                                atoms=group.atoms,
                                group_id=group.group_id,
                                formal_charge=group.formal_charge,
                                output_dir=antechamber_output_dir
                            )

                    elif group.parameterization_method == 'existing':
                        # Parse existing mol2/frcmod
                        self.console.print(f"[grey50]Parsing existing files: {Path(group.existing_mol2).name}...[/grey50]")
                        params = antechamber.parse_existing_files(
                            mol2_file=group.existing_mol2,
                            frcmod_file=group.existing_frcmod,
                            atoms=group.atoms,
                            group_id=group.group_id
                        )

                    elif group.parameterization_method == 'manual':
                        # Manual entry
                        self.console.print(f"[grey50]Manual entry for {group.group_id}...[/grey50]")
                        params = antechamber.manual_entry_fallback(group.atoms, group.group_id)

                    else:
                        self.console.print(f"[red]Unknown parameterization method: {group.parameterization_method}[/red]")
                        continue

                    # Store results if successful
                    if params.success:
                        for coords, gaff_type in params.atom_types.items():
                            is_center = coords in center_coords
                            is_ligand = coords in ligand_coords

                            # Find original atom for resname/resid info
                            orig_atom = next((a for a in group.atoms if a.coords == coords), None)

                            if orig_atom:
                                type_assignments[coords] = {
                                    'chain': orig_atom.chain,
                                    'resname': orig_atom.resname,
                                    'ff_resname': orig_atom.resname,
                                    'resid': orig_atom.resid,
                                    'atom_name': orig_atom.atom_name,
                                    'element': params.elements[coords],
                                    'original_type': gaff_type,
                                    'renamed_type': gaff_type,
                                    'charge': params.charges[coords],
                                    'terminal_type': 'middle',
                                    'is_center': is_center,
                                    'is_metal_ligand': is_ligand,
                                    'renamed': False,
                                    'source': f'{group.parameterization_method}_{group.group_id}',
                                    'ligand_group': group.group_id
                                }
                        self.console.print(f"[green]✓ {group.group_id} parameterized successfully[/green]")
                    else:
                        # Failed - ask if user wants to retry with different method
                        self.console.print(f"[red]❌ Failed to parameterize {group.group_id}: {params.error_message}[/red]")
                        self.console.print("[yellow]Falling back to manual entry...[/yellow]")

                        manual_params = antechamber.manual_entry_fallback(group.atoms, group.group_id)

                        for coords, atom_type in manual_params.atom_types.items():
                            orig_atom = next((a for a in group.atoms if a.coords == coords), None)
                            if orig_atom:
                                is_center = coords in center_coords
                                is_ligand = coords in ligand_coords

                                type_assignments[coords] = {
                                    'chain': orig_atom.chain,
                                    'resname': orig_atom.resname,
                                    'ff_resname': orig_atom.resname,
                                    'resid': orig_atom.resid,
                                    'atom_name': orig_atom.atom_name,
                                    'element': manual_params.elements[coords],
                                    'original_type': atom_type,
                                    'renamed_type': atom_type,
                                    'charge': manual_params.charges[coords],
                                    'terminal_type': 'middle',
                                    'is_center': is_center,
                                    'is_metal_ligand': is_ligand,
                                    'renamed': False,
                                    'source': 'manual_entry_fallback'
                                }

        return type_assignments

    def _convert_global_types(
        self,
        global_types: Dict,
        ff_data=None
    ) -> Dict:
        """
        Convert GlobalTypeAssignment objects to Step 1 dict format.

        GlobalTypeAssignment comes from GlobalAtomTypeRegistry for multi-site
        coordination. This method converts to the dict format expected by
        downstream Step 1 processing (renaming, fingerprint generation, etc.).

        Args:
            global_types: Dict[Tuple[float,float,float], GlobalTypeAssignment]
            ff_data: Optional ForceFieldData for charge lookup

        Returns:
            Dict[Tuple[float,float,float], Dict] in Step 1 format
        """
        type_assignments = {}

        # A withheld cluster atom has no library type, so original_type below
        # falls back to the raw ELEMENT. For a metal or a bridging sulfide that
        # is harmless -- 'MO'/'FE'/'S' are not Amber types, and those atoms are
        # renamed to M*/Y* anyway. For a HYDROGEN it is not: 'H' IS a valid
        # Amber type, the amide/amine one, so a hydroxo proton silently
        # acquired its 0.6 A vdW radius instead of the hydroxyl 0.0.
        #
        # Index the residues so an inferred hydrogen can be typed from the atom
        # it is bonded to.
        residue_atoms: Dict[Tuple, List[Tuple]] = {}
        for coords, ga in global_types.items():
            residue_atoms.setdefault((ga.chain, ga.resid), []).append(
                (coords, ga.element))

        # Restrained coordinating ligands (e.g. nonbonded waters) must not be
        # Y-renamed even though the global registry flagged them as ligands.
        restrained_coords, _ = _collect_restrained_ligands(
            getattr(self, 'provided_redox_site', None)
        )

        for coords, global_assignment in global_types.items():
            # Look up charge from ff_data if available
            charge = None
            if ff_data is not None:
                resname = global_assignment.resname
                atom_name = global_assignment.atom_name

                # Try direct lookup
                atom_def = ff_data.get_atom_definition(resname, atom_name)
                if atom_def is not None:
                    charge = atom_def.charge
                else:
                    # Check residue aliases (e.g., HEM → HEH)
                    alias = ff_data.residue_aliases.get(resname)
                    if alias:
                        atom_def = ff_data.get_atom_definition(alias, atom_name)
                        if atom_def is not None:
                            charge = atom_def.charge

            # Build Step 1 format dict
            type_assignments[coords] = {
                'chain': global_assignment.chain,
                'resname': global_assignment.resname,
                'ff_resname': global_assignment.resname,
                'resid': global_assignment.resid,
                'atom_name': global_assignment.atom_name,
                'element': global_assignment.element,
                'original_type': (
                    global_assignment.original_type
                    or _inferred_element_type(global_assignment, residue_atoms)
                    or ""),
                'renamed_type': (
                    global_assignment.global_renamed_type
                    or _inferred_element_type(global_assignment, residue_atoms)
                    or ""),
                'charge': charge,
                'terminal_type': 'internal',  # GlobalTypes are typically internal
                'is_center': global_assignment.is_metal,
                'is_metal_ligand': global_assignment.is_ligand and (
                    tuple(round(float(x), 3) for x in coords) not in restrained_coords
                ),
                'renamed': global_assignment.global_renamed_type is not None,
                'source': 'global_registry',
                'mass': None,
                'spin': None,
                'vdw_params': None
            }

        # Show charge statistics
        atoms_with_charge = sum(1 for a in type_assignments.values() if a.get('charge') is not None)
        atoms_without_charge = len(type_assignments) - atoms_with_charge
        if atoms_with_charge > 0 or atoms_without_charge > 0:
            self.console.print(f"[cyan]Charges: {atoms_with_charge} atoms with FF charges, {atoms_without_charge} without[/cyan]")

        return type_assignments

    def _cluster_internal_ligand_coords(self, redox_site: 'RedoxSite',
                                        metal_coords: set) -> set:
        """Coords of pure-cluster internal atoms that should be typed as ligands.

        A pure inorganic metal cluster (Fe2S2, Fe4S4, a Mo-S-O core) has its
        bridging/terminal atoms bonded to the metal WITHIN the residue. Those
        intra-residue bonds are not emitted as 'coordinate' bonds, so the atoms
        would keep generic element types (S, O) — which can't carry per-atom
        bonded parameters and are inconsistent with the Y-typed external
        ligands. A pure cluster residue (multi-atom, contains a metal, contains
        no carbon → no organic scaffold) has no non-core atoms, so every one of
        its non-metal HEAVY atoms is metal-coordinating core. Return their
        coords so the caller can mark them as ligands (→ unique Y* types).

        Metals themselves are excluded (they stay is_center / M*), and so are
        hydrogens: a hydroxo H added to a Mo-O core is bonded to the oxygen,
        not to the metal, so typing it as a metal ligand would be wrong. The
        "no non-core atoms" assumption held only while a cluster could not
        carry a hydrogen at all.
        """
        # Group the site's atoms by residue.
        by_res: Dict[Tuple, list] = {}
        for atom in redox_site.atoms:
            by_res.setdefault((atom.chain, atom.resid, atom.resname), []).append(atom)

        coords = set()
        for atoms in by_res.values():
            if len(atoms) <= 1:
                continue  # a lone metal ion is not a cluster
            has_metal = any((a.element or '').strip().upper() in METALS for a in atoms)
            has_carbon = any((a.element or '').strip().upper() == 'C' for a in atoms)
            if has_metal and not has_carbon:
                for a in atoms:
                    if a.coords in metal_coords:      # metals stay M*
                        continue
                    if (a.element or '').strip().upper() == 'H':
                        continue                      # bonded to a core atom, not the metal
                    coords.add(a.coords)
        return coords

    def _convert_preprocessing_types(
        self,
        preprocessing_types: Dict,
        redox_site: 'RedoxSite'
    ) -> Dict:
        """
        Convert AtomTypeAssignment objects from preprocessing to Step 1 dict format.

        The preprocessing stage (Step 0g) produces AtomTypeAssignment dataclass
        objects keyed by coordinate tuples. Step 1 expects dicts with certain
        keys for further processing (terminal classification, MCPB renaming).

        Args:
            preprocessing_types: Dict[Tuple[float,float,float], AtomTypeAssignment]
            redox_site: RedoxSite object for filtering (only include atoms in site)

        Returns:
            Dict[Tuple[float,float,float], Dict] in Step 1 format
        """
        type_assignments = {}

        # Build coord sets from RedoxSite for quick lookup
        site_coords = {atom.coords for atom in redox_site.atoms}

        # Build metal_coords from centers (handles both METAL_ION and ORGANOMETALLIC_COFACTOR)
        # - METAL_ION: center.coords IS the metal
        # - ORGANOMETALLIC_COFACTOR: find metal atoms in the same residue
        metal_coords = set()
        for center in redox_site.centers:
            if hasattr(center, 'center_type'):
                if center.center_type.value == 'metal_ion':
                    # Isolated metal - center coords ARE the metal
                    metal_coords.add(center.coords)
                elif center.center_type.value == 'organometallic_cofactor':
                    # Embedded metal - find metal atom(s) in this residue
                    for atom in redox_site.atoms:
                        if (atom.chain == center.chain and
                            atom.resid == center.resid and
                            atom.resname == center.resname and
                            atom.element.upper() in METALS):
                            metal_coords.add(atom.coords)

        # Build ligand coords (atoms bonded to metals via coordinate bonds).
        # Restrained contacts are excluded — they stay nonbonded (no Y-renaming).
        ligand_coords = set()
        for bond in redox_site.bonds:
            if bond.chemical_type == 'coordinate' and getattr(bond, 'treatment', 'bonded') == 'bonded':
                if bond.atom1_coords in metal_coords:
                    ligand_coords.add(bond.atom2_coords)
                elif bond.atom2_coords in metal_coords:
                    ligand_coords.add(bond.atom1_coords)

        # Pure-cluster internal atoms (Fe-S bridging sulfides, Mo-S/O core) are
        # metal-coordinating but bonded within the residue, so they never appear
        # as 'coordinate' bonds above. Treat them as ligands too, for unique Y*
        # types instead of generic element types.
        ligand_coords |= self._cluster_internal_ligand_coords(redox_site, metal_coords)

        # Convert each AtomTypeAssignment in the site to dict format
        for coords, assignment in preprocessing_types.items():
            # Only include atoms that are in this RedoxSite
            if coords not in site_coords:
                continue

            # Determine metal/ligand status
            is_center = coords in metal_coords
            is_ligand = coords in ligand_coords

            # Convert enum to string for terminal_type
            terminal_str = (
                assignment.terminal_type.value
                if hasattr(assignment.terminal_type, 'value')
                else str(assignment.terminal_type)
            )

            # Convert source enum to string
            source_str = (
                assignment.source.value
                if hasattr(assignment.source, 'value')
                else str(assignment.source)
            )

            # Build Step 1 format dict
            type_assignments[coords] = {
                'chain': assignment.chain,
                'resname': assignment.resname,
                'ff_resname': assignment.ff_resname or assignment.resname,
                'resid': assignment.resid,
                'atom_name': assignment.atom_name,
                'element': assignment.element,
                'original_type': assignment.original_type,
                'renamed_type': assignment.renamed_type or assignment.original_type,
                'charge': assignment.charge,
                'terminal_type': terminal_str,
                'is_center': is_center,
                'is_metal_ligand': is_ligand,
                'renamed': assignment.renamed,
                'source': source_str,
                # Include metal-specific fields if present
                'mass': assignment.mass,
                'spin': assignment.spin,
                'vdw_params': {
                    'radius': assignment.vdw_radius,
                    'epsilon': assignment.vdw_epsilon,
                    'source': assignment.vdw_source
                } if assignment.vdw_radius else None
            }

        self.console.print(f"[grey50]Converted {len(type_assignments)} types for RedoxSite atoms[/grey50]")
        return type_assignments

    def _site_models_dir(self) -> Path:
        """THIS site's ``models/`` directory.

        ``step_results`` is restored in ``__init__`` from ``mcpb_step_results``,
        a single workspace key every site shares, so ``step_1`` belongs to
        whichever site wrote it last. The checklist sets ``step_3a`` per site,
        so prefer that and fall back to step_1 only for the standalone workflow.

        Reading step_1 unconditionally is what cross-wired site 1 to site 2's
        standard.fingerprint: their PDB serial ranges do not overlap, every
        atom-type lookup missed, and the mol2 files -- and the libraries built
        from them -- were written with the ``XX`` placeholder type. tleap then
        had no parameter for anything in that site.
        """
        out_dir = (self.step_results.get("step_3a") or {}).get("output_dir")
        if out_dir:
            return Path(out_dir)

        output_files = (self.step_results.get("step_1") or {}).get("output_files") or {}
        for key in ("standard_fingerprint", "large_pdb"):
            recorded = output_files.get(key)
            if recorded:
                return Path(recorded).parent

        raise ValueError("Cannot locate this site's models directory")

    def _water_model(self) -> str:
        """The LJ-table water model for the system being prepared.

        Metal Lennard-Jones parameters are fitted against a specific water
        model, so the frcmod has to use the set matching the model the system
        is solvated in. This was hardcoded to tip3p, which quietly gave every
        non-TIP3P system the wrong metal LJ terms.
        """
        leaprc = ""
        workspace = getattr(self, 'workspace', None)
        if workspace is None and getattr(self, 'processor', None) is not None:
            try:
                workspace = self.processor._get_workspace()
            except Exception:  # noqa: BLE001 — no workspace is not an error here
                workspace = None
        if workspace is not None:
            try:
                leaprc = workspace.get("preprocessing_water_model", "") or ""
            except Exception:  # noqa: BLE001
                leaprc = ""

        return water_model_from_leaprc(leaprc, logger=self.logger)

    def _collect_metal_radii(self, type_assignments: Dict) -> Dict[str, float]:
        """``{element: vdW radius}`` for the site's metals, for MK ReadRadii.

        ReadRadii tells Gaussian to read custom van der Waals radii for the
        Merz-Kollman ESP fit; without them a metal gets no radius of its own and
        its ESP sampling falls back to Gaussian's defaults.

        Sourced from the METAL ATOMS in ``type_assignments``, not from
        ``redox_site.centers``. A center's ``coords`` is the metal atom only for
        a lone ``metal_ion``; for an organometallic cofactor or a pure cluster
        (FES, MOS) the center describes the RESIDUE — its ``element`` is None
        and its ``coords`` are not any atom's — so the old center-based lookup
        matched nothing and every cluster site silently emitted a bare
        ``Pop=MK``. ``type_assignments`` already marks each metal atom
        ``is_center``, carrying its element and the charge the user entered.
        """
        metal_radii: Dict[str, float] = {}
        radius_charges: Dict[str, int] = {}
        metal_db = MetalIonDatabase(water_model=self._water_model(), logger=self.logger)

        for _coords, assignment in type_assignments.items():
            if isinstance(assignment, dict):
                if not assignment.get('is_center'):
                    continue
                element = assignment.get('element')
                charge_val = assignment.get('charge')
            else:
                if not getattr(assignment, 'is_center', False):
                    continue
                element = getattr(assignment, 'element', None)
                charge_val = getattr(assignment, 'charge', None)

            if not element or charge_val is None:
                continue

            # Normalize element symbol (e.g., 'ZN' -> 'Zn')
            element_norm = (element[0].upper() + element[1:].lower()
                            if len(element) > 1 else element.upper())
            radius = metal_db.get_vdw_radius(element_norm, int(charge_val))
            if not radius:
                continue

            # Gaussian's ReadRadii block is per ELEMENT, so two metals of the
            # same element in different oxidation states cannot each carry
            # their own radius. Keep the first and say so.
            if element_norm in metal_radii:
                if radius_charges.get(element_norm) != int(charge_val):
                    self.console.print(
                        f"[yellow]  {element_norm} appears at charge "
                        f"{radius_charges[element_norm]:+d} and {int(charge_val):+d}; "
                        f"ReadRadii is per element, so "
                        f"{metal_radii[element_norm]:.3f} Å is used for both.[/yellow]")
                continue

            metal_radii[element_norm] = radius
            radius_charges[element_norm] = int(charge_val)
            self.console.print(
                f"[grey50]  Metal radius: {element_norm} = {radius:.3f} Å[/grey50]")

        if not metal_radii:
            self.console.print(
                "[yellow]  No metal radii resolved — the ESP will use Gaussian's "
                "defaults (Pop=MK without ReadRadii).[/yellow]")
        return metal_radii

    @staticmethod
    def _default_large_model_solvation(small_scrf: str, large_charge: int) -> bool:
        """Whether to pre-answer yes to solvating the large model.

        Yes when the small model was solvated — running the two models at
        different levels of theory should be a deliberate act, not the result
        of pressing Enter. Otherwise yes only for an anion, where implicit
        solvation improves SCF convergence and the charge distribution the ESP
        is fitted to.

        The user is still asked either way; this only sets the default.
        """
        return bool(small_scrf) or large_charge < 0

    def _collect_metal_charges(
        self,
        type_assignments: Dict,
        redox_site: 'RedoxSite',
        interactive: bool
    ) -> Dict:
        """
        Collect formal charges (and metal spins) for atoms with no prmtop charge.

        Two groups need this, both because they were removed before tLEaP:

        - Metals, which are reinserted afterwards. They need a formal charge and
          a spin state; the charge also keys the van der Waals radius lookup.
        - A withheld cluster's non-metal core atoms — the bridging sulfides of
          an Fe-S cluster, the S/O core of a Mo cofactor. A pure inorganic
          cluster is withheld from the force-field pass as a whole residue, so
          these carry no partial charge either.

        Both are marked the same way: ``charge is None`` on a site atom, with
        ``is_center`` separating the metals from the rest. Collecting only the
        metals left the core atoms at None, which the PDB writer counts as 0.0,
        so the suggested QM charge came out high by their formal charge — an
        Fe2S2 site proposed +2 where ``[Fe2S2(SCys)4]2-`` is -2. Correcting that
        by inflating a metal's charge is worse than it looks: the metal charge is
        also the vdW-radius key and is stored in the deposited library.

        Args:
            type_assignments: Dict of coord -> assignment dict
            redox_site: RedoxSite object
            interactive: Whether to prompt interactively

        Returns:
            Updated type_assignments with charges/spins filled in
        """
        # Split the charge-less site atoms: centers are metals (charge + spin +
        # vdW), everything else is withheld-cluster core (formal charge only).
        metals_needing_charge = []
        core_atoms_needing_charge = []
        for coords, assignment in type_assignments.items():
            if assignment.get('charge') is not None:
                continue
            if assignment.get('is_center'):
                metals_needing_charge.append((coords, assignment))
            else:
                core_atoms_needing_charge.append((coords, assignment))

        if not metals_needing_charge and not core_atoms_needing_charge:
            return type_assignments

        n_metals = len(metals_needing_charge)
        if n_metals:
            self.console.print("\n[bold]Metal Charge Collection[/bold]")
            self.console.print("[grey50]Metals aren't in prmtop - need user input for charge/spin[/grey50]")
            if n_metals > 1:
                self.console.print(
                    f"[grey50]{n_metals} metals in this site are set one at a time; "
                    f"the residue ID below identifies which one.[/grey50]"
                )

        metal_db = MetalIonDatabase(water_model=self._water_model(), logger=self.logger)

        for metal_idx, (coords, assignment) in enumerate(metals_needing_charge, 1):
            element = assignment.get('element', '')
            resname = assignment.get('resname', '')
            atom_name = assignment.get('atom_name', '')
            chain = assignment.get('chain', '')
            resid = assignment.get('resid', '')

            # Disambiguate WHICH metal this is by residue ID (two Mn in a
            # binuclear site otherwise print identically as "MN MN (MN)").
            if chain and str(chain).strip():
                location = f"{chain}:{resid}"
            elif resid not in (None, ''):
                location = f"residue {resid}"
            else:
                location = atom_name  # last-resort identifier
            counter = f" {metal_idx}/{n_metals}" if n_metals > 1 else ""
            self.console.print(
                f"\n[cyan]Metal{counter}: {resname} {location} (atom {atom_name}, {element})[/cyan]"
            )

            if interactive:
                # Prompt for charge — include the residue so the input line is
                # self-identifying even when scrolled back.
                # Name the ATOM, not just the residue. A cluster residue holds
                # several metals, so "formal charge for FES residue 1311" is
                # asked twice, identically, for FE1 and then FE2 — and the two
                # can legitimately differ (a mixed-valence Fe(III)/Fe(II) pair).
                # The core-atom prompts below already name their atom.
                while True:
                    charge_input = prompt_with_context(
                        self.processor,
                        f"Enter formal charge for {resname} {location} atom "
                        f"{atom_name} ({element}) (e.g., +2, +3, -1)",
                        module="Metal Site Parameterizer",
                        description="Metal formal charge",
                    ).strip()
                    try:
                        charge = int(charge_input.replace('+', ''))
                        break
                    except ValueError:
                        self.console.print("[red]Invalid charge format[/red]")

                # Prompt for spin
                while True:
                    spin_input = prompt_with_context(
                        self.processor,
                        f"Enter unpaired electrons for {resname} {location} atom "
                        f"{atom_name} ({element}) (e.g., 0, 1, 5)",
                        module="Metal Site Parameterizer",
                        description="Metal unpaired electrons",
                    ).strip()
                    try:
                        spin = int(spin_input)
                        break
                    except ValueError:
                        self.console.print("[red]Invalid spin format[/red]")
            else:
                # Non-interactive: use defaults
                charge = 2  # Common default
                spin = 0
                self.console.print(f"  [grey50]Using defaults: charge={charge}, spin={spin}[/grey50]")

            # Get metal config from database for VDW parameters.
            # First try the (resname, atom_name) key — works for standalone ions
            # whose resname IS the element (ZN, FE, MN...). For a cluster the
            # resname is the cluster code (FES, SF4) not a metal name, so that
            # lookup misses; fall back to the ELEMENT, which we already know from
            # the atom. Only if the element itself is unknown do we ask the user.
            metal_config = metal_db.get_metal_config(resname, atom_name, charge, spin)
            if metal_config is None and element:
                metal_config = metal_db.get_metal_config_by_element(
                    element, atom_name, resname, charge, spin
                )
            if metal_config is None:
                # Manual fallback (element genuinely not in the database)
                metal_config = self._manual_metal_entry(
                    type('obj', (object,), {'resname': resname, 'atom_name': atom_name, 'element': element})(),
                    charge, spin
                )

            # Update the assignment with all metal properties from database
            assignment['charge'] = float(metal_config.charge)
            assignment['spin'] = metal_config.spin
            assignment['mass'] = metal_config.mass
            if metal_config.vdw_radius:
                assignment['vdw_params'] = {
                    'radius': metal_config.vdw_radius,
                    'epsilon': metal_config.vdw_epsilon,
                    'source': metal_config.vdw_source or 'MetalIonDatabase'
                }
                # These become the metal's Lennard-Jones terms in the frcmod
                # NONB section. When the database had no parameters for this
                # (element, charge) it returns a generic placeholder and logs
                # "VERIFY MANUALLY" — easy to miss next to the ESP radius,
                # which is a DIFFERENT quantity from a different table and may
                # well have resolved. Say plainly which one is a placeholder.
                if 'fallback' in str(metal_config.vdw_source or '').lower():
                    self.console.print(
                        f"  [yellow]⚠ {element} {charge:+d} has no tabulated "
                        f"Lennard-Jones parameters; the frcmod will carry a "
                        f"generic placeholder (r={metal_config.vdw_radius}, "
                        f"eps={metal_config.vdw_epsilon}). This is the force-field "
                        f"vdW term, not the ESP radius — supply a value before "
                        f"using these parameters.[/yellow]")
            else:
                # No radius for this (element, charge). Usually the charge is not
                # a real oxidation state for the element — which is what happens
                # when someone compensates on the metal for a core atom that
                # contributed nothing to the total. Say so here rather than
                # silently leaving the site without a metal radius: this charge
                # is also stored in the deposited library.
                self.console.print(
                    f"  [yellow]⚠ No van der Waals radius for {element or '?'} "
                    f"at charge {charge:+d}. Check that this is the intended "
                    f"oxidation state — the formal charge of a cluster's "
                    f"bridging/core atoms is asked separately below.[/yellow]"
                )

            self.console.print(
                f"  [green]✓ {resname} {location} {atom_name}: "
                f"charge={charge:+d}, spin={spin}[/green]")

        # ------------------------------------------------------------------
        # Withheld-cluster core atoms (bridging sulfides, Mo-S/O core)
        # ------------------------------------------------------------------
        # These have no prmtop charge for the same reason the metals don't: the
        # cluster is withheld from the force-field pass as a whole residue. They
        # need only a formal charge — no spin, no vdW, and nothing downstream
        # reads them as partial charges (step 3B's RESP restraints take their
        # charges from preprocessing_atom_data, not from here).
        if core_atoms_needing_charge:
            self._collect_cluster_core_charges(core_atoms_needing_charge, interactive)

        return type_assignments

    def _collect_cluster_core_charges(self, core_atoms: List, interactive: bool) -> None:
        """Ask for the formal charge of each withheld-cluster core atom.

        Mutates each assignment in place. Left at None, these atoms are summed
        as 0.0 by the PDB writer and the suggested QM charge comes out short by
        their formal charge.

        Args:
            core_atoms: List of (coords, assignment) with charge still None.
            interactive: Whether to prompt. Non-interactive runs leave the
                charges at 0.0 and say so — there is no defensible default
                across a sulfide, an oxo and a hydroxo.
        """
        n_core = len(core_atoms)
        self.console.print("\n[bold]Cluster Core Atom Charges[/bold]")
        self.console.print(
            f"[grey50]{n_core} atom(s) belong to a withheld cluster residue, so they "
            f"carry no prmtop charge. Their formal charge completes the QM total "
            f"(a bridging sulfide is S2-, a terminal oxo O2-).[/grey50]"
        )

        if not interactive:
            for _coords, assignment in core_atoms:
                assignment['charge'] = 0.0
            self.console.print(
                "[yellow]  Non-interactive: left at 0. The suggested QM charge "
                "will be short by these atoms' formal charge.[/yellow]"
            )
            return

        for _coords, assignment in core_atoms:
            resname = assignment.get('resname', '')
            atom_name = assignment.get('atom_name', '')
            element = assignment.get('element', '')
            chain = assignment.get('chain', '')
            resid = assignment.get('resid', '')

            if chain and str(chain).strip():
                location = f"{chain}:{resid}"
            elif resid not in (None, ''):
                location = f"residue {resid}"
            else:
                location = atom_name

            while True:
                raw = prompt_with_context(
                    self.processor,
                    f"Formal charge for {resname} {location} atom {atom_name} "
                    f"({element}) (e.g., -2, -1, 0)",
                    module="Metal Site Parameterizer",
                    description="Withheld-cluster core atom formal charge",
                ).strip()
                try:
                    core_charge = int(raw.replace('+', ''))
                    break
                except ValueError:
                    self.console.print("[red]Invalid charge format[/red]")

            assignment['charge'] = float(core_charge)
            self.console.print(
                f"  [green]✓ {resname} {location} {atom_name}: "
                f"charge={core_charge:+d}[/green]"
            )

    def _manual_metal_entry(self, atom, charge: int, spin: int) -> MetalConfig:
        """
        Manual entry fallback for metals not in database.

        Args:
            atom: RedoxSiteAtom object
            charge: Formal charge (already obtained from user)
            spin: Spin state (already obtained from user)

        Returns:
            MetalConfig with manually entered data
        """
        self.console.print(f"\n[yellow]Manual entry required for metal: {atom.resname} {atom.atom_name}[/yellow]")

        # The element is usually already known from the atom (e.g. Fe of an
        # Fe-S cluster) — don't re-ask for it. Only prompt when it is genuinely
        # missing/unparseable.
        element = (getattr(atom, 'element', '') or '').strip().capitalize()
        if not element:
            element = prompt_with_context(
                self.processor, "Element symbol (e.g., Fe, Cu)",
                module="Metal Site Parameterizer",
                description="Metal element symbol",
            ).strip().capitalize()
        else:
            self.console.print(f"  [grey50]Element: {element} (from structure)[/grey50]")
        mass = float(input(f"  Atomic mass (amu): ").strip())

        # Use default VDW params (charge and spin already provided)
        from .mcpb.metal_ion_database import MetalConfig
        return MetalConfig(
            resname=atom.resname,
            atom_name=atom.atom_name,
            element=element,
            atomic_number=0,  # Unknown
            mass=mass,
            atom_type=atom.resname,
            charge=charge,  # Store user-provided charge
            spin=spin,  # Store user-provided spin
            vdw_radius=1.5,
            vdw_epsilon=0.01,
            vdw_source=f"Manual entry for {element}{charge:+d} (VERIFY)",
            water_model=self._water_model()
        )

    def _print_withheld_cluster_charge_note(self, redox_site, type_assignments) -> None:
        """Warn that the suggested QM charge omits atoms with no charge at all.

        A pure metal cluster (Fe-S, etc.) is removed from the prmtop as a whole
        residue, so neither its metals nor its non-metal core atoms (the bridging
        sulfides, a Mo-S-O core) arrive with a charge. Both are now asked for in
        _collect_metal_charges, so in the normal path nothing is left at None and
        this prints nothing.

        It remains as a safety net for an atom that reaches the model without
        passing through that collection — a gap residue merged in later, or a
        non-interactive run that left the core atoms at 0. Anything still None
        contributes 0.0 to the suggested total, so the suggestion runs high.
        Coordinates key type_assignments and are stable across tLEaP
        renumbering, so a None charge on a site atom is a reliable marker.
        """
        names = []
        for atom in getattr(redox_site, 'atoms', []) or []:
            coords = getattr(atom, 'coords', None)
            if coords is None:
                continue
            assignment = type_assignments.get(coords)
            if assignment is None:
                assignment = type_assignments.get(tuple(round(c, 3) for c in coords))
            if assignment is None:
                continue
            if hasattr(assignment, 'charge'):
                charge = assignment.charge
            elif isinstance(assignment, dict):
                charge = assignment.get('charge')
            else:
                charge = None
            if charge is None:
                names.append(f"{getattr(atom, 'resname', '?')} {getattr(atom, 'atom_name', '?')}")

        if not names:
            return

        shown = ", ".join(names[:6]) + (" ..." if len(names) > 6 else "")
        self.console.print(
            f"[yellow]⚠ This suggestion omits {len(names)} atom(s) ({shown}) that "
            f"still carry no charge and were counted as 0.[/yellow]"
        )
        self.console.print(
            "[yellow]  Their formal charge is not in the sum, so correct the value "
            "to the real net charge of your target oxidation/spin state.[/yellow]"
        )

    def _apply_systematic_renaming(self, redox_site: RedoxSite, type_assignments: Dict,
                                   metal_start: int = 0, ligand_start: int = 0) -> Dict:
        """
        Apply systematic renaming to metal and ligand atoms.

        Renames metal atoms to M1-M9, MA-MZ and ligands to Y1-Y9, YA-YZ.
        Uses Y prefix (not L) to match MCPB conventions and GlobalAtomTypeRegistry.

        For a multi-site protein, the M*/Y* names MUST be globally unique because
        every site's frcmod/mol2 is loaded into the SAME tLEaP session. Passing
        metal_start/ligand_start continues the numbering from where the previous
        site stopped (site 1 -> M1/Y1..Yn, site 2 -> M2/Y(n+1)..). The ending
        offsets are stashed on the instance (type_offset_metal_end /
        type_offset_ligand_end) so the caller can thread them to the next site.

        Args:
            redox_site: RedoxSite object
            type_assignments: Current type assignments
            metal_start: First metal index to use (0-based into metal_names)
            ligand_start: First ligand index to use (0-based into ligand_names)

        Returns:
            Updated type assignments with renamed atoms
        """
        metal_names = MCPB_METAL_TYPE_NAMES
        ligand_names = MCPB_LIGAND_TYPE_NAMES

        metal_index = metal_start
        ligand_index = ligand_start

        for coords, assignment in type_assignments.items():
            if assignment['is_center'] and metal_index < len(metal_names):
                # Rename metal
                assignment['renamed_type'] = metal_names[metal_index]
                assignment['renamed'] = True
                metal_index += 1

            elif assignment['is_metal_ligand'] and ligand_index < len(ligand_names):
                # Rename ligand
                assignment['renamed_type'] = ligand_names[ligand_index]
                assignment['renamed'] = True
                ligand_index += 1

        # Stash ending offsets so a multi-site caller can continue the numbering
        # for the next site (avoids M1/Y1 colliding across sites).
        self.type_offset_metal_end = metal_index
        self.type_offset_ligand_end = ligand_index

        return type_assignments

    def _run_step2a(self, residue_name: str, residues: List, output_dir: Path,
                    interactive: bool) -> Dict[str, Any]:
        """
        Run Step 2a: Generate Pre-frcmod File (Before QM Calculations).

        Following MCPB workflow, this creates a preliminary frcmod file with:
        - NON markers: Parameters that need QM calculations (metal-ligand bonds/angles)
        - YES markers: Parameters available from force field (organic bonds/angles)

        This provides:
        1. Educational checkpoint - shows user what needs QM
        2. Resume capability - can restart from here if QM fails
        3. Template for final frcmod generation

        Process:
        1. Load Step 1 outputs (atom assignments, fingerprint, XYZ)
        2. Parse fingerprint to extract bonds and angles
        3. Initialize force field parameter reader
        4. Generate pre-frcmod with NON/YES markers
        5. Display summary of missing parameters

        Args:
            residue_name: Name of metal site residue
            residues: List of residue objects (not used)
            output_dir: Output directory containing Step 1 results
            interactive: Enable user prompts

        Returns:
            Dict with success status and pre-frcmod file path
        """
        self.console.print("[bold cyan]▸ Pre-frcmod generation (before QM)[/bold cyan]")

        try:
            # ================================================================
            # 1. Verify model outputs exist
            # ================================================================
            step1_dir = output_dir / "models"
            step2_dir = output_dir / "bonded_params"
            step2_dir.mkdir(parents=True, exist_ok=True)

            small_pdb = step1_dir / "small.pdb"
            standard_fp = step1_dir / "standard.fingerprint"
            assignments_file = step1_dir / "atom_type_assignments.json"

            if not all([small_pdb.exists(), standard_fp.exists(), assignments_file.exists()]):
                return {
                    "success": False,
                    "message": "Model outputs not found. Please complete model building first."
                }

            self.console.print(f"[green]✅ Found model outputs[/green]")

            # ================================================================
            # 2. Load atom type assignments
            # ================================================================
            import json
            with open(assignments_file) as f:
                type_assignments_raw = json.load(f)

            # Convert string keys back to tuples
            type_assignments = {}
            for key_str, assignment in type_assignments_raw.items():
                coords = tuple(float(x) for x in key_str.split(','))
                type_assignments[coords] = assignment

            self.console.print(f"[grey50]Loaded {len(type_assignments)} atom type assignments[/grey50]")

            # ================================================================
            # 3. Parse fingerprint file to extract bonds and angles
            # ================================================================
            self.console.print("\n[bold]Parsing fingerprint file...[/bold]")

            bonds, angles, serial_to_coords, boundary_type_assignments = \
                self._parse_fingerprint_for_topology(standard_fp, small_pdb)

            # Merge the standard-typed cap atoms (ACE/NME boundary) into the type
            # assignments so the pre-frcmod generator can type the boundary
            # angle/dihedral terms that cross into the caps. These records use
            # renamed=False, so they add no MASS/NONB entries.
            for coords, record in boundary_type_assignments.items():
                type_assignments.setdefault(coords, record)
            if boundary_type_assignments:
                self.console.print(
                    f"[grey50]Typed {len(boundary_type_assignments)} cap-boundary "
                    f"atom(s) for boundary parameter enumeration[/grey50]")

            self.console.print(f"[grey50]Found {len(bonds)} bonds and {len(angles)} angles[/grey50]")

            # ================================================================
            # 4. Get parameter provider from prmtop
            # ================================================================
            workspace = self.processor._get_workspace() if self.processor else None
            prmtop_path = workspace.get("parm7_file") if workspace else None

            if not prmtop_path or not os.path.exists(prmtop_path):
                return {
                    "success": False,
                    "message": "prmtop file not available. Please complete preprocessing first."
                }

            param_provider = PrmtopParameterProvider(prmtop_path, self.console)
            self.console.print(f"\n[bold]Parameters from: {os.path.basename(prmtop_path)}[/bold]")
            self.console.print(f"[grey50]Loaded: {len(param_provider.bond_parameters)} bonds, "
                             f"{len(param_provider.angle_parameters)} angles[/grey50]")

            # ================================================================
            # 5. Generate pre-frcmod file
            # ================================================================
            from .mcpb.pre_frcmod_generator import PreFrcmodGenerator

            pre_frcmod_gen = PreFrcmodGenerator(param_provider, self.console)

            pre_frcmod_file = step2_dir / f"{residue_name}_pre.frcmod"

            stats = pre_frcmod_gen.generate_pre_frcmod(
                type_assignments,
                bonds,
                angles,
                pre_frcmod_file,
                serial_to_coords
            )

            # ================================================================
            # 7. Save bond topology for reuse in Step 2B and Step 3D
            # ================================================================
            topology_file = step2_dir / "bond_topology.json"
            self._save_bond_topology(bonds, angles, serial_to_coords, topology_file)

            self.console.print(f"[grey50]Saved bond topology: {topology_file.name}[/grey50]")

            # ================================================================
            # 8. Display next steps
            # ================================================================
            self.console.print(f"\n[bold green]✅ Pre-frcmod generation complete[/bold green]")
            self.console.print(f"[cyan]Pre-frcmod file created: {pre_frcmod_file.name}[/cyan]")
            self.console.print(
                f"[grey50]{len(stats['non_bonds'])} bond(s) and {len(stats['non_angles'])} "
                f"angle(s) require QM (Seminario) values[/grey50]")

            return {
                "success": True,
                "step_number": "2a",
                "step_description": "Pre-frcmod Generation (Before QM)",
                "next_step": "2b",
                "output_files": {
                    "pre_frcmod": str(pre_frcmod_file),
                    "bond_topology": str(topology_file)
                },
                "statistics": {
                    "non_bonds": len(stats['non_bonds']),
                    "yes_bonds": len(stats['yes_bonds']),
                    "non_angles": len(stats['non_angles']),
                    "yes_angles": len(stats['yes_angles']),
                    "warnings": len(stats.get('warnings', []))
                },
                "warnings": stats.get('warnings', [])
            }

        except Exception as e:
            import traceback
            traceback.print_exc()

            return {
                "success": False,
                "message": f"Step 2a error: {str(e)}",
                "traceback": traceback.format_exc()
            }

    def _run_step2b(self, residue_name: str, residues: List, output_dir: Path,
                   interactive: bool) -> Dict[str, Any]:
        """
        Run Step 2b: Bonded Parameter Generation via Seminario Method (After QM).

        This step follows Step 2a (pre-frcmod generation) and completes the MCPB workflow
        by computing force constants from QM calculations and merging with the pre-frcmod template.

        Process:
        1. Check for Step 2a outputs (pre-frcmod file)
        2. Check for Step 1 outputs (XYZ file, fingerprint, atom assignments)
        3. Guide user through QM calculation setup
        4. Parse QM outputs (fchk/log files)
        5. Apply Seminario method to generate force constants
        6. Merge computed parameters with pre-frcmod template (replaces NON markers)
        7. Validate parameters

        Args:
            residue_name: Name of metal site residue
            residues: List of residue objects (not used)
            output_dir: Output directory containing Step 1 and 2a results
            interactive: Enable user prompts

        Returns:
            Dict with success status and output files
        """
        self.console.print("[bold cyan]▸ Bonded parameter generation (Seminario method)[/bold cyan]")

        try:
            # ================================================================
            # 1. Verify model outputs exist
            # ================================================================
            step1_dir = output_dir / "models"
            step2_dir = output_dir / "bonded_params"
            step2_dir.mkdir(parents=True, exist_ok=True)

            # Check for required model files
            small_pdb = step1_dir / "small.pdb"
            standard_fp = step1_dir / "standard.fingerprint"
            assignments_file = step1_dir / "atom_type_assignments.json"

            if not all([small_pdb.exists(), standard_fp.exists(), assignments_file.exists()]):
                return {
                    "success": False,
                    "message": "Model outputs not found. Please complete model building first."
                }

            self.console.print(f"[green]✅ Found Step 1 outputs[/green]")

            # Load atom type assignments
            import json
            with open(assignments_file) as f:
                type_assignments_raw = json.load(f)

            # Convert string keys back to tuples
            type_assignments = {}
            for key_str, assignment in type_assignments_raw.items():
                coords = tuple(float(x) for x in key_str.split(','))
                type_assignments[coords] = assignment

            # ================================================================
            # 2. Check for QM Output Files
            # ================================================================
            from .mcpb.qm_interface import QMInterface, QMSoftware, QMCalculationMode
            from proprep.utils.prompts import prompt_with_context, confirm_with_context, int_prompt_with_context
            from proprep.utils.file_browser import remap_recorded_index, annotate_selected_path

            self.console.print("\n[bold]QM Calculation Check[/bold]")
            self.console.print("Step 2 requires Gaussian output from optimization + frequency calculation.\n")

            # Check for expected output files from Step 1's generated input
            # Step 1 generates small_freq.gjf for Seminario method
            small_fchk = step1_dir / "small_freq.fchk"
            small_log = step1_dir / "small_freq.log"

            log_file = None
            fchk_file = None
            qm_mode = None

            if small_fchk.exists() and small_log.exists():
                self.console.print(f"[green]✓ Found QM output files from Step 1:[/green]")
                self.console.print(f"  • {small_log.name}")
                self.console.print(f"  • {small_fchk.name}")

                if confirm_with_context(
                    self.processor,
                    "Use these files for Seminario analysis?",
                    default=True,
                    module="Metal Site Parameterizer",
                    description="Use existing QM files"
                ):
                    log_file = small_log  # Keep as Path
                    fchk_file = small_fchk  # Keep as Path
                    qm_mode = QMCalculationMode.EXTERNAL
            elif small_log.exists() and not small_fchk.exists():
                self.console.print(f"[yellow]⚠ Found {small_log.name} but missing {small_fchk.name}[/yellow]")
                self.console.print("[grey50]Run: formchk small_freq.chk small_freq.fchk[/grey50]\n")
            elif not small_log.exists():
                small_gjf = step1_dir / "small_freq.gjf"
                if small_gjf.exists():
                    self.console.print(f"[yellow]⚠ Gaussian input exists but output not found[/yellow]")
                    self.console.print(f"  Input: {small_gjf.name}")
                    self.console.print(f"  Expected: {small_log.name}, {small_fchk.name}")
                    self.console.print("\n[grey50]Run Gaussian calculation first, then return to this step.[/grey50]\n")

            # If no files found or user declined, offer options
            if qm_mode is None:
                self.console.print("[bold]QM Output Options:[/bold]")
                self.console.print("  [cyan]1[/cyan] - Provide path to existing .log and .fchk files")
                self.console.print("  [cyan]2[/cyan] - Generate new Gaussian input (if Step 1 input was modified)\n")

                mode_choice = prompt_with_context(
                    self.processor,
                    "Select option",
                    default="1",
                    module="Metal Site Parameterizer",
                    description="QM output option"
                )

                if mode_choice == "1":
                    qm_mode = QMCalculationMode.EXTERNAL
                else:
                    qm_mode = QMCalculationMode.GUIDED

            # Currently only Gaussian is supported
            qm_software = QMSoftware.GAUSSIAN
            qm_interface = QMInterface(qm_software, qm_mode, self.console)

            # ================================================================
            # 3. Handle QM Calculation Workflow
            # ================================================================
            if qm_mode == QMCalculationMode.GUIDED:
                # Generate input template
                self.console.print("\n[bold]Generating Gaussian input template...[/bold]\n")

                # ================================================================
                # Computational Resources
                # ================================================================
                self.console.print("[bold]Computational Resources:[/bold]")
                memory_gb = int_prompt_with_context(
                    self.processor, "  Memory allocation (GB)", default=4,
                    module="MCPB Step 2", description="QM memory allocation"
                )
                n_processors = int_prompt_with_context(
                    self.processor, "  Number of processors", default=4,
                    module="MCPB Step 2", description="QM processor count"
                )

                # ================================================================
                # QM System Properties
                # ================================================================
                # Read total charge from small.pdb REMARK record
                suggested_charge = 0
                try:
                    with open(small_pdb) as f:
                        for line in f:
                            if line.startswith("REMARK") and "Total charge:" in line:
                                charge_str = line.split("Total charge:")[1].strip()
                                suggested_charge = round(float(charge_str))
                                break
                except Exception:
                    suggested_charge = 0

                # Calculate suggested multiplicity from metal spins
                # Total S = sum(spin_i / 2), Multiplicity = 2|S| + 1
                total_spin = 0
                for assignment in type_assignments.values():
                    if 'spin' in assignment and assignment.get('element'):
                        # Only count metals (they have 'spin' field)
                        total_spin += assignment['spin']

                # Calculate multiplicity: 2S + 1 where S = total_spin / 2
                total_S = total_spin / 2.0
                suggested_multiplicity = int(2 * abs(total_S) + 1)

                self.console.print("\n[bold]QM System Properties:[/bold]")
                charge = int_prompt_with_context(
                    self.processor, "  Total charge of small model", default=suggested_charge,
                    module="MCPB Step 2", description="QM total charge"
                )
                multiplicity = int_prompt_with_context(
                    self.processor, "  Spin multiplicity (2S+1)", default=suggested_multiplicity,
                    module="MCPB Step 2", description="QM spin multiplicity"
                )

                # ================================================================
                # QM Method
                # ================================================================
                self.console.print("\n[bold]QM Method:[/bold]")
                functional = prompt_with_context(
                    self.processor, "  Functional", default="B3LYP",
                    module="MCPB Step 2", description="QM functional"
                )
                basis_set = prompt_with_context(
                    self.processor, "  Basis set", default="6-31G*",
                    module="MCPB Step 2", description="QM basis set"
                )

                # Handle GenECP basis set specification
                basis_groups = []
                ecp_specs = []
                if basis_set.lower() in ["genecp", "gen"]:
                    # Detect unique elements in small model
                    elements_in_model = set()
                    with open(small_pdb) as f:
                        for line in f:
                            if line.startswith("ATOM") or line.startswith("HETATM"):
                                element = line[76:78].strip()
                                if not element:
                                    atom_name = line[12:16].strip()
                                    element = ''.join([c for c in atom_name if c.isalpha()])[:2]
                                    if len(element) == 2:
                                        element = element[0].upper() + element[1].lower()
                                    else:
                                        element = element.upper()
                                elements_in_model.add(element)

                    elements_list = sorted(elements_in_model)
                    self.console.print(f"\n[cyan]Setting up mixed basis set ({basis_set.upper()})...[/cyan]")
                    self.console.print(f"  Detected elements in structure: {', '.join(elements_list)}")

                    # Collect basis set groups
                    group_num = 1
                    used_atoms = set()

                    while True:
                        self.console.print(f"\n  [bold]Group {group_num}[/bold]")

                        # Suggest remaining atoms
                        remaining = [e for e in elements_list if e not in used_atoms]
                        if not remaining:
                            break

                        atoms_input = prompt_with_context(
                            self.processor, "    Atoms (space-separated)",
                            default=" ".join(remaining) if group_num == 1 else "",
                            module="MCPB Step 2", description="GenECP atoms for basis group"
                        )

                        if not atoms_input.strip():
                            break

                        atoms = atoms_input.strip().split()
                        used_atoms.update(atoms)

                        # Suggest basis set based on atom type
                        default_basis = "LANL2DZ" if any(a in ["Fe", "Cu", "Zn", "Mn", "Co", "Ni", "Mg", "Ca"] for a in atoms) else "6-31G(d)"
                        basis = prompt_with_context(
                            self.processor, "    Basis set for this group", default=default_basis,
                            module="MCPB Step 2", description="GenECP basis set for group"
                        )

                        basis_groups.append((atoms, basis))

                        if not confirm_with_context(
                            self.processor, "  Add another basis set group?", default=False,
                            module="MCPB Step 2", description="Add another GenECP basis group"
                        ):
                            break

                        group_num += 1

                    # Handle ECPs if using GenECP
                    if basis_set.lower() == "genecp":
                        # Suggest metals for ECP
                        metals = [a for atoms, _ in basis_groups for a in atoms
                                 if a in ["Fe", "Cu", "Zn", "Mn", "Co", "Ni", "Mg", "Ca", "Pt", "Pd", "Ag", "Au"]]

                        self.console.print(f"\n  [bold]Effective Core Potentials (ECPs):[/bold]")
                        ecp_atoms_input = prompt_with_context(
                            self.processor, "    Atoms needing ECP (space-separated)",
                            default=" ".join(metals) if metals else "",
                            module="MCPB Step 2", description="ECP atoms selection"
                        )

                        if ecp_atoms_input.strip():
                            ecp_atoms = ecp_atoms_input.strip().split()
                            ecp_basis = prompt_with_context(
                                self.processor, "    ECP basis for these atoms", default="LANL2DZ",
                                module="MCPB Step 2", description="ECP basis set"
                            )
                            ecp_specs.append((ecp_atoms, ecp_basis))

                # ================================================================
                # Job Type
                # ================================================================
                self.console.print("\n[bold]Job Type:[/bold]")
                self.console.print("  [cyan]1[/cyan] Opt Freq - Optimize geometry + frequencies (recommended)")
                self.console.print("  [cyan]2[/cyan] Freq - Frequencies only (if structure pre-optimized)")
                job_type_choice = prompt_with_context(
                    self.processor, "  Select", choices=["1", "2"], default="1",
                    module="MCPB Step 2", description="QM job type selection",
                    options_map={"1": "Opt Freq", "2": "Freq only"}
                )
                job_type = "Opt Freq" if job_type_choice == "1" else "Freq"

                # ================================================================
                # Additional Keywords and Title
                # ================================================================
                self.console.print("\n[grey50]Note: Geom=PrintInputOrient, Integral=(Grid=UltraFine), and IOp(7/33=1) are always included.[/grey50]")
                additional_keywords = prompt_with_context(
                    self.processor, "Additional keywords (optional)", default="",
                    module="MCPB Step 2", description="Additional Gaussian keywords"
                )
                title_card = prompt_with_context(
                    self.processor, "Title card", default="MCPB Step 2 - Small model",
                    module="MCPB Step 2", description="Gaussian title card"
                )

                # Generate input files with all parameters
                qm_files = qm_interface.generate_input_files(
                    pdb_file=small_pdb,
                    output_dir=step2_dir,
                    charge=charge,
                    multiplicity=multiplicity,
                    memory_gb=memory_gb,
                    n_processors=n_processors,
                    functional=functional,
                    basis_set=basis_set,
                    basis_groups=basis_groups,
                    ecp_specs=ecp_specs,
                    job_type=job_type,
                    additional_keywords=additional_keywords,
                    title_card=title_card
                )

                input_file = qm_files["input_file"]
                expected_log = qm_files["log_file"]

                if not confirm_with_context(
                    self.processor, "\nHave you completed these steps?", default=False,
                    module="MCPB Step 2", description="Confirm QM calculation complete"
                ):
                    return {
                        "success": False,
                        "status": "awaiting_qm",
                        "message": "Waiting for QM calculation completion",
                        "input_file": str(input_file),
                        "expected_output": str(expected_log)
                    }

            # ================================================================
            # Select QM output files (both GUIDED and EXTERNAL modes)
            # Skip if files were already selected from Step 1 output
            # ================================================================
            if qm_mode in [QMCalculationMode.GUIDED, QMCalculationMode.EXTERNAL] and log_file is None:
                self.console.print("\n[bold]Select QM output files:[/bold]")

                # List available log files
                log_files = list(step2_dir.glob("*.log")) + list(step2_dir.glob("*.out"))

                if log_files:
                    self.console.print("\n[cyan]Available log files:[/cyan]")
                    for i, f in enumerate(log_files, 1):
                        self.console.print(f"  [{i}] {f.name}")

                    log_choice = prompt_with_context(
                        self.processor, "\nSelect log file number or enter custom path",
                        default="1" if log_files else "",
                        module="MCPB Step 2", description="Select QM log file"
                    )

                    log_choice = remap_recorded_index(self.processor, log_files, str(log_choice))
                    if log_choice.isdigit() and 1 <= int(log_choice) <= len(log_files):
                        log_file = log_files[int(log_choice) - 1]
                        annotate_selected_path(self.processor, log_file)
                    else:
                        log_file = Path(log_choice)
                else:
                    log_file = Path(prompt_with_context(
                        self.processor, "Path to log file (.log or .out)", default="",
                        module="MCPB Step 2", description="QM log file path"
                    ))

                if not log_file.exists():
                    return {"success": False, "message": f"Log file not found: {log_file}"}

            # ================================================================
            # 4. Verify QM calculation completed
            # ================================================================
            self.console.print("\n[bold]Verifying QM calculation...[/bold]")
            success, message = qm_interface.check_calculation_complete(log_file)

            if not success:
                return {
                    "success": False,
                    "message": f"QM calculation check failed: {message}"
                }

            self.console.print(f"[green]✅ {message}[/green]")

            # ================================================================
            # 5. Handle fchk file (Gaussian only)
            # Skip if fchk_file was already selected from Step 1 output
            # ================================================================
            if qm_software == QMSoftware.GAUSSIAN and fchk_file is None:
                # List available fchk files in step2 directory
                fchk_files = list(step2_dir.glob("*.fchk")) + list(step2_dir.glob("*.fch"))

                if fchk_files:
                    self.console.print("\n[cyan]Available fchk files:[/cyan]")
                    for i, f in enumerate(fchk_files, 1):
                        self.console.print(f"  [{i}] {f.name}")

                    fchk_choice = prompt_with_context(
                        self.processor, "\nSelect fchk file number, enter custom path, or type 'convert' to generate from .chk",
                        default="1" if fchk_files else "convert",
                        module="MCPB Step 2", description="Select fchk file"
                    )
                    fchk_choice = remap_recorded_index(self.processor, fchk_files, str(fchk_choice))

                    if fchk_choice.lower() == "convert":
                        # User wants to convert chk to fchk
                        chk_files = list(step2_dir.glob("*.chk"))

                        if chk_files:
                            self.console.print("\n[cyan]Available chk files:[/cyan]")
                            for i, f in enumerate(chk_files, 1):
                                self.console.print(f"  [{i}] {f.name}")

                            chk_choice = prompt_with_context(
                                self.processor, "Select chk file number or enter custom path",
                                default="1" if chk_files else "",
                                module="MCPB Step 2", description="Select chk file for conversion"
                            )

                            chk_choice = remap_recorded_index(self.processor, chk_files, str(chk_choice))
                            if chk_choice.isdigit() and 1 <= int(chk_choice) <= len(chk_files):
                                chk_file = chk_files[int(chk_choice) - 1]
                                annotate_selected_path(self.processor, chk_file)
                            else:
                                chk_file = Path(chk_choice)
                        else:
                            chk_file = Path(prompt_with_context(
                                self.processor, "Path to chk file", default="",
                                module="MCPB Step 2", description="Chk file path"
                            ))

                        if not chk_file.exists():
                            return {"success": False, "message": f"chk file not found: {chk_file}"}

                        # Try to run formchk
                        fchk_file = qm_interface.generate_fchk(chk_file)
                        if not fchk_file:
                            return {
                                "success": False,
                                "message": "Failed to generate fchk file. Please run formchk manually: "
                                           f"formchk {chk_file.name} {chk_file.stem}.fchk"
                            }
                    elif fchk_choice.isdigit() and 1 <= int(fchk_choice) <= len(fchk_files):
                        fchk_file = fchk_files[int(fchk_choice) - 1]
                        annotate_selected_path(self.processor, fchk_file)
                    else:
                        fchk_file = Path(fchk_choice)
                else:
                    # No fchk files found, ask user what to do
                    self.console.print("\n[yellow]No fchk files found in step2 directory[/yellow]")

                    action = prompt_with_context(
                        self.processor, "Enter path to existing fchk file or type 'convert' to generate from .chk",
                        default="convert",
                        module="MCPB Step 2", description="fchk file action"
                    )

                    if action.lower() == "convert":
                        chk_files = list(step2_dir.glob("*.chk"))

                        if chk_files:
                            self.console.print("\n[cyan]Available chk files:[/cyan]")
                            for i, f in enumerate(chk_files, 1):
                                self.console.print(f"  [{i}] {f.name}")

                            chk_choice = prompt_with_context(
                                self.processor, "Select chk file number or enter custom path",
                                default="1" if chk_files else "",
                                module="MCPB Step 2", description="Select chk file for conversion"
                            )

                            chk_choice = remap_recorded_index(self.processor, chk_files, str(chk_choice))
                            if chk_choice.isdigit() and 1 <= int(chk_choice) <= len(chk_files):
                                chk_file = chk_files[int(chk_choice) - 1]
                                annotate_selected_path(self.processor, chk_file)
                            else:
                                chk_file = Path(chk_choice)
                        else:
                            chk_file = Path(prompt_with_context(
                                self.processor, "Path to chk file", default="",
                                module="MCPB Step 2", description="Chk file path"
                            ))

                        if not chk_file.exists():
                            return {"success": False, "message": f"chk file not found: {chk_file}"}

                        fchk_file = qm_interface.generate_fchk(chk_file)
                        if not fchk_file:
                            return {
                                "success": False,
                                "message": "Failed to generate fchk file. Please run formchk manually: "
                                           f"formchk {chk_file.name} {chk_file.stem}.fchk"
                            }
                    else:
                        fchk_file = Path(action)

                if not fchk_file.exists():
                    return {"success": False, "message": f"fchk file not found: {fchk_file}"}

            # ================================================================
            # 6. Extract Hessian matrix and coordinates
            # ================================================================
            self.console.print("\n[bold]Extracting Hessian matrix from QM output...[/bold]")

            from .mcpb.hessian_parser import HessianParser

            try:
                if qm_software == QMSoftware.GAUSSIAN:
                    hessian, coords = HessianParser.parse_gaussian_fchk(fchk_file)
                elif qm_software == QMSoftware.GAMESS:
                    hessian, coords = HessianParser.parse_gamess_log(log_file)
                else:
                    return {"success": False, "message": f"{qm_software.value} not yet supported for Hessian extraction"}

                self.console.print(f"[green]✅ Extracted Hessian: {hessian.shape}[/green]")
                self.console.print(f"[green]✅ Extracted coordinates: {len(coords)//3} atoms[/green]")

                # Validate Hessian
                is_valid, validation_msg = HessianParser.validate_hessian(hessian, coords)
                if not is_valid:
                    return {"success": False, "message": f"Hessian validation failed: {validation_msg}"}

                self.console.print(f"[green]✅ Hessian validation passed[/green]")

                # Show Hessian statistics
                stats = HessianParser.get_hessian_statistics(hessian)
                self.console.print(f"[grey50]Hessian eigenvalue range: {stats['min_eigenvalue']:.2e} to {stats['max_eigenvalue']:.2e}[/grey50]")

            except Exception as e:
                return {"success": False, "message": f"Failed to extract Hessian: {str(e)}"}

            # ================================================================
            # 7. Apply Seminario method
            # ================================================================
            self.console.print("\n[bold]Computing force constants via Seminario method...[/bold]")

            from .mcpb.seminario import SeminarioMethod

            scale_factor_str = prompt_with_context(
                self.processor, "Frequency scaling factor", default="1.0",
                module="MCPB Step 2", description="Seminario frequency scaling factor"
            )
            scale_factor = float(scale_factor_str)

            seminario = SeminarioMethod(hessian, coords, scale_factor, self.console)

            # Extract bonds and angles from fingerprint file
            bond_params, angle_params, bonds = self._extract_parameters_from_fingerprint(
                seminario, standard_fp, type_assignments, small_pdb, coords
            )

            self.console.print(f"[green]✅ Computed {len(bond_params)} bond parameters[/green]")
            self.console.print(f"[green]✅ Computed {len(angle_params)} angle parameters[/green]")

            # Show parameter summary
            self._show_parameter_summary(bond_params, angle_params)

            # ================================================================
            # 8. Generate frcmod file
            # ================================================================
            self.console.print("\n[bold]Generating frcmod file...[/bold]")

            from .mcpb.frcmod_builder import FrcmodBuilder

            frcmod_builder = FrcmodBuilder()

            # Add bond parameters
            for bond in bond_params:
                frcmod_builder.add_bond_parameter(bond)

            # Add angle parameters
            for angle in angle_params:
                frcmod_builder.add_angle_parameter(angle)

            # ================================================================
            # Add force field library parameters for organic bonds/angles
            # ================================================================
            self.console.print("\n[bold]Adding force field library parameters...[/bold]")

            # Get parameter provider from prmtop
            workspace = self.processor._get_workspace() if self.processor else None
            prmtop_path = workspace.get("parm7_file") if workspace else None

            if prmtop_path and os.path.exists(prmtop_path):
                self.console.print(f"[grey50]Using parameters from prmtop: {os.path.basename(prmtop_path)}[/grey50]")
                param_provider = PrmtopParameterProvider(prmtop_path, self.console)
            else:
                self.console.print("[yellow]⚠️  No prmtop available for library parameters[/yellow]")
                param_provider = None

            lib_bond_count, lib_angle_count = self._add_library_parameters(
                frcmod_builder, bond_params, angle_params, bonds_list=bonds,
                ff_data=param_provider
            )

            if lib_bond_count > 0 or lib_angle_count > 0:
                self.console.print(
                    f"[green]✅ Added {lib_bond_count} library bond parameters, "
                    f"{lib_angle_count} library angle parameters[/green]"
                )
            else:
                self.console.print("[grey50]No additional library parameters needed[/grey50]")

            # ================================================================
            # Add MASS and NONB entries for renamed atoms
            # ================================================================
            self.console.print("\n[bold]Adding MASS and NONB parameters...[/bold]")

            mass_count = self._add_mass_and_nonb_parameters(
                frcmod_builder, type_assignments, param_provider
            )

            if mass_count > 0:
                self.console.print(
                    f"[green]✅ Added {mass_count} MASS/NONB entries[/green]"
                )

            # ================================================================
            # Add inherited dihedral/improper parameters for Y* types
            # ================================================================
            self.console.print("\n[bold]Adding inherited dihedral parameters for renamed atom types...[/bold]")
            dihedral_count = self._add_inherited_dihedral_parameters(
                frcmod_builder, type_assignments, param_provider
            )
            if dihedral_count > 0:
                self.console.print(
                    f"[green]✅ Added {dihedral_count} inherited dihedral/improper entries[/green]"
                )
            else:
                self.console.print("[grey50]No additional dihedral parameters needed[/grey50]")

            # ================================================================
            # Check for pre-frcmod from Step 2a and merge
            # ================================================================
            pre_frcmod_file = step2_dir / f"{residue_name}_pre.frcmod"
            frcmod_file = step2_dir / f"{residue_name}_bonded.frcmod"

            if pre_frcmod_file.exists():
                self.console.print(f"\n[cyan]Merging with pre-frcmod template: {pre_frcmod_file.name}[/cyan]")
                frcmod_builder.merge_with_pre_frcmod(pre_frcmod_file, frcmod_file, method="Seminario")
                self.console.print(f"[green]✅ Generated frcmod file (merged with pre-frcmod): {frcmod_file.name}[/green]")
            else:
                self.console.print(f"\n[yellow]⚠️  No pre-frcmod found. Writing direct frcmod file.[/yellow]")
                self.console.print(f"[grey50]Tip: Run Step 2a first to generate pre-frcmod with NON/YES markers[/grey50]")
                frcmod_builder.write_frcmod(frcmod_file, method="Seminario")
                self.console.print(f"[green]✅ Generated frcmod file: {frcmod_file.name}[/green]")

            # ================================================================
            # 9. Validate parameters
            # ================================================================
            is_valid, warnings = frcmod_builder.validate_parameters()

            if warnings:
                self.console.print(f"\n[yellow]⚠️  Parameter validation warnings:[/yellow]")
                for warning in warnings:
                    self.console.print(f"  [yellow]• {warning}[/yellow]")

            # ================================================================
            # 10. Save detailed results
            # ================================================================
            results_json = step2_dir / "step2_results.json"
            self._save_step2_results(
                results_json, bond_params, angle_params, stats, scale_factor
            )

            self.console.print(f"\n[green]✅ Bonded parameter generation complete[/green]")

            return {
                "success": True,
                "step_number": 2,
                "step_description": "Bonded Parameter Generation (Seminario Method)",
                "output_files": {
                    "frcmod_file": str(frcmod_file),
                    "results_json": str(results_json),
                    "log_file": str(log_file),
                    "fchk_file": str(fchk_file) if qm_software == QMSoftware.GAUSSIAN else None
                },
                "statistics": {
                    "n_bonds": len(bond_params),
                    "n_angles": len(angle_params),
                    "scale_factor": scale_factor,
                    "hessian_stats": stats,
                    "validation_warnings": len(warnings)
                },
                "next_step": 3
            }

        except Exception as e:
            import traceback
            traceback.print_exc()

            return {
                "success": False,
                "message": f"Step 2 error: {str(e)}",
                "traceback": traceback.format_exc()
            }

    # Standard ff14SB atom types for the atoms of the ACE/NME QM caps. The caps
    # saturate the valences cut when the QM small model is truncated; in the FULL
    # topology each cap atom stands in for a real backbone atom of the neighbouring
    # residue (NME:N -> next residue's amide N, ACE:C -> previous residue's carbonyl
    # C). A bonded term that crosses that boundary therefore EXISTS in the assembled
    # topology and must carry a parameter in the site frcmod — otherwise, once a
    # ligating atom is retyped (e.g. backbone O -> Y6), tleap has no term for the
    # boundary angle (Y6-C-N) or dihedral (N-C-Y6-M1) and aborts. Typing the cap
    # atom with its standard type lets the pre-frcmod generator resolve those terms
    # to their ordinary ff14SB values.
    _CAP_ATOM_TYPES = {
        ('ACE', 'C'): 'C', ('ACE', 'O'): 'O', ('ACE', 'CH3'): 'CT',
        ('ACE', 'HH31'): 'HC', ('ACE', 'HH32'): 'HC', ('ACE', 'HH33'): 'HC',
        ('NME', 'N'): 'N', ('NME', 'H'): 'H', ('NME', 'CH3'): 'CT',
        ('NME', 'HH31'): 'H1', ('NME', 'HH32'): 'H1', ('NME', 'HH33'): 'H1',
    }
    _CAP_RESNAMES = frozenset(('ACE', 'NME'))

    def _perceive_cap_boundary_bonds(self, pdb_file, fp_serials):
        """
        Perceive bonds crossing a QM-cap (ACE/NME) boundary and type the cap atoms.

        The metal-free preprocessing prmtop never contains the ACE/NME caps (they
        are added only when the QM small model is built), so the peptide bond from
        a fingerprint atom to a cap atom cannot be recovered from the prmtop. It is
        recovered here by covalent-radius distance perception over the small-model
        coordinates, where the cap sits at a real bonded distance (~1.34 A for C-N).

        Only cap atoms that are NOT themselves fingerprint atoms are treated as
        boundary caps (a cap oxygen that directly ligates a metal is a typed
        fingerprint atom and is left alone).

        Args:
            pdb_file: small-model PDB (has cap atoms + coordinates + serials)
            fp_serials: set of PDB serials that belong to the typed fingerprint

        Returns:
            (boundary_bonds, cap_atoms) where
            - boundary_bonds: list of (cap_serial, fp_serial) tuples
            - cap_atoms: dict cap_serial -> {'atom_name','element','renamed_type'}
              for the cap atoms that actually formed a boundary bond
        """
        import math

        # Covalent radii (A); sum + tolerance gives the bond cutoff. Tolerance
        # matches the mol2 intra-residue perception so the two agree.
        radii = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05}
        tol = 0.45

        cap_serials = {}   # serial -> (atom_name, element, resname, coords)
        fp_atoms = {}      # serial -> (element, coords)
        with open(pdb_file) as f:
            for line in f:
                if not (line.startswith('ATOM') or line.startswith('HETATM')):
                    continue
                try:
                    serial = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    resname = line[17:20].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                except (ValueError, IndexError):
                    continue
                element = line[76:78].strip() or atom_name[0]
                element = element[0].upper() + element[1:].lower()
                coords = (x, y, z)
                if serial in fp_serials:
                    fp_atoms[serial] = (element, coords)
                elif resname in self._CAP_RESNAMES:
                    cap_serials[serial] = (atom_name, element, resname, coords)

        boundary_bonds = []
        cap_atoms = {}
        for cap_serial, (atom_name, element, resname, ccoords) in cap_serials.items():
            r_cap = radii.get(element, 0.77)
            for fp_serial, (fp_element, fcoords) in fp_atoms.items():
                r_fp = radii.get(fp_element, 0.77)
                dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(ccoords, fcoords)))
                if dist <= r_cap + r_fp + tol:
                    boundary_bonds.append((cap_serial, fp_serial))
                    cap_type = self._CAP_ATOM_TYPES.get(
                        (resname, atom_name), element.upper()
                    )
                    cap_atoms[cap_serial] = {
                        'atom_name': atom_name,
                        'element': element,
                        'renamed_type': cap_type,
                    }
        return boundary_bonds, cap_atoms

    def _parse_fingerprint_for_topology(
        self,
        fingerprint_file: Path,
        pdb_file: Path
    ) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int, int]], Dict[int, Tuple[float, float, float]], Dict[Tuple[float, float, float], Dict]]:
        """
        Parse fingerprint file to extract bonds and angles topology.

        Simplified version for Step 2a that just extracts connectivity.

        Args:
            fingerprint_file: Path to standard.fingerprint
            pdb_file: PDB file containing atom coordinates and serial numbers

        Returns:
            (bonds, angles, serial_to_coords, boundary_type_assignments) tuple where:
            - bonds: List of (serial1, serial2) tuples (PDB serial numbers)
            - angles: List of (serial1, serial2, serial3) tuples (PDB serial numbers)
            - serial_to_coords: Dict mapping PDB serial number to (x,y,z) coordinates
            - boundary_type_assignments: Dict coords -> type-assignment record for the
              ACE/NME cap atoms that form a boundary bond, so the pre-frcmod
              generator can type the boundary angle/dihedral terms
        """
        # ================================================================
        # 1. Parse fingerprint file to get connectivity
        # ================================================================
        atom_info = {}  # AtomID -> {atom_name, renamed_type}
        link_records = []  # List of (AtomID1, AtomID2) pairs

        with open(fingerprint_file) as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith('LINK'):
                    # Parse LINK record
                    parts = line.split()
                    if len(parts) >= 3:
                        atom1_info = parts[1].split('-')
                        atom2_info = parts[2].split('-')

                        try:
                            atom1_id = int(atom1_info[0])
                            atom2_id = int(atom2_info[0])
                            link_records.append((atom1_id, atom2_id))
                        except (ValueError, IndexError):
                            pass

                elif '-' in line and '->' in line:
                    # Parse atom line
                    parts = line.split()
                    if len(parts) >= 5:
                        try:
                            atom_id = int(parts[1])
                            renamed_type = parts[4].strip()
                            atom_info[atom_id] = {'renamed_type': renamed_type}
                        except (ValueError, IndexError):
                            pass

        # ================================================================
        # 2. Read PDB file to get serial number -> (index, coordinates) mapping
        # ================================================================
        serial_to_index = {}  # PDB serial -> 0-based index
        idx_to_coords = {}    # 0-based index -> (x, y, z)

        index = 0
        with open(pdb_file) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Parse PDB ATOM record
                    # Format: ATOM  serial name resname chain resid    x       y       z
                    try:
                        serial = int(line[6:11].strip())
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        serial_to_index[serial] = index
                        idx_to_coords[index] = (x, y, z)
                        index += 1
                    except (ValueError, IndexError):
                        pass

        # ================================================================
        # 3. Detect all covalent bonds from PDB coordinates
        # ================================================================
        # Build reverse mapping for bond detection
        index_to_serial = {idx: serial for serial, idx in serial_to_index.items()}

        # Create atom_info dict compatible with _detect_distance_bonds
        # Map PDB serial -> atom info from fingerprint
        fp_atom_info = {}

        with open(fingerprint_file) as f:
            for line in f:
                if '-' in line and '->' in line and not line.startswith('LINK'):
                    parts = line.split()
                    if len(parts) >= 5:
                        res_parts = parts[0].split('-')
                        if len(res_parts) >= 3:
                            try:
                                atom_id = int(parts[1])
                                atom_name = res_parts[2]
                                renamed_type = parts[4].strip()
                                fp_atom_info[atom_id] = {
                                    'atom_name': atom_name,
                                    'renamed_type': renamed_type
                                }
                            except (ValueError, IndexError):
                                pass

        # ----------------------------------------------------------------
        # Detect organic bonds: prefer prmtop (accurate) over distance
        # detection (can produce false positives from QM-optimized geometry)
        # ----------------------------------------------------------------
        import numpy as np

        workspace = self.processor._get_workspace() if self.processor else None
        prmtop_path = workspace.get("parm7_file") if workspace else None

        organic_bonds = []
        if prmtop_path and os.path.exists(prmtop_path):
            organic_bonds = self._extract_bonds_from_prmtop(
                prmtop_path, fp_atom_info, serial_to_index
            )
        else:
            # Fallback: distance detection (uses QM-optimized coords, may be distorted)
            self.console.print("[yellow]⚠️  No prmtop available — using distance-based bond detection[/yellow]")
            A_TO_B = 1.0 / 0.529177249
            n_atoms = len(idx_to_coords)
            coords_array = np.zeros(n_atoms * 3)
            for idx, (x, y, z) in idx_to_coords.items():
                coords_array[3*idx] = x * A_TO_B
                coords_array[3*idx + 1] = y * A_TO_B
                coords_array[3*idx + 2] = z * A_TO_B
            organic_bonds = self._detect_distance_bonds(
                coords_array, index_to_serial, fp_atom_info, tolerance=0.40
            )

        # Convert LINK records to same format
        link_bonds = []
        for serial1, serial2 in link_records:
            if serial1 in serial_to_index and serial2 in serial_to_index:
                idx1 = serial_to_index[serial1]
                idx2 = serial_to_index[serial2]
                type1 = fp_atom_info.get(serial1, {}).get('renamed_type', 'XX')
                type2 = fp_atom_info.get(serial2, {}).get('renamed_type', 'XX')
                link_bonds.append((idx1, idx2, type1, type2))

        # ----------------------------------------------------------------
        # Recover the peptide bonds that cross a QM-cap boundary (ACE/NME).
        # The prmtop path above can never contain these bonds because the caps
        # are not part of the metal-free preprocessing prmtop, so we perceive
        # them by distance over the small-model coordinates. Each such bond joins
        # a typed fingerprint atom to a cap atom that we type with its standard
        # ff14SB type below; without it the boundary angle/dihedral to a retyped
        # ligating atom (e.g. Y6-C-N, N-C-Y6-M1) never enumerates.
        # ----------------------------------------------------------------
        boundary_bonds, boundary_cap_atoms = self._perceive_cap_boundary_bonds(
            pdb_file, set(fp_atom_info.keys())
        )
        for cap_serial, fp_serial in boundary_bonds:
            if cap_serial in serial_to_index and fp_serial in serial_to_index:
                cap_idx = serial_to_index[cap_serial]
                fp_idx = serial_to_index[fp_serial]
                cap_type = boundary_cap_atoms[cap_serial]['renamed_type']
                fp_type = fp_atom_info.get(fp_serial, {}).get('renamed_type', 'XX')
                organic_bonds.append((cap_idx, fp_idx, cap_type, fp_type))
        if boundary_cap_atoms:
            self.console.print(
                f"[grey50]Recovered {len(boundary_bonds)} cap-boundary bond(s) "
                f"(ACE/NME) so boundary angle/dihedral terms enumerate[/grey50]")

        # Exclude atoms of restrained ligands (e.g. a nonbonded water held by an
        # MD restraint) from the bonded-term topology. The water stays in the QM
        # model for correct electronics, but it must emit NO bonded terms: its
        # OW-HW bond and HW-OW-HW angle are already defined by leaprc.water, and
        # re-emitting a QM value in the site frcmod would override that generic,
        # type-based term for EVERY water in the system. Every water bond has the
        # coordinating O as an endpoint, so excluding bonds that touch the
        # restrained O coords removes the O-H bonds (and thus the H-O-H angle).
        restrained_coords, _ = _collect_restrained_ligands(
            getattr(self, 'provided_redox_site', None)
        )
        restrained_indices = set()
        if restrained_coords:
            for idx, xyz in idx_to_coords.items():
                if tuple(round(float(c), 3) for c in xyz) in restrained_coords:
                    restrained_indices.add(idx)

        # Combine LINK bonds and organic bonds (remove duplicates)
        bond_set = set()
        n_restrained_skipped = 0
        for idx1, idx2, *_ in link_bonds + organic_bonds:
            if idx1 in restrained_indices or idx2 in restrained_indices:
                n_restrained_skipped += 1
                continue
            bond_set.add((min(idx1, idx2), max(idx1, idx2)))
        if n_restrained_skipped:
            self.console.print(
                f"[grey50]Excluded {n_restrained_skipped} bond(s) of restrained "
                f"(nonbonded) ligand(s) from the frcmod — provided by the water model[/grey50]")

        bonds = [(idx1, idx2) for idx1, idx2 in sorted(bond_set)]

        # ================================================================
        # 4. Generate angles from bonds
        # ================================================================
        # Build adjacency list
        adjacency = {}
        for idx1, idx2 in bonds:
            if idx1 not in adjacency:
                adjacency[idx1] = []
            if idx2 not in adjacency:
                adjacency[idx2] = []
            adjacency[idx1].append(idx2)
            adjacency[idx2].append(idx1)

        # Find all angles (3 consecutive bonded atoms)
        angles = []
        for center_idx in adjacency:
            neighbors = adjacency[center_idx]
            # For each pair of neighbors, create an angle
            for i in range(len(neighbors)):
                for j in range(i+1, len(neighbors)):
                    idx1 = neighbors[i]
                    idx3 = neighbors[j]
                    angles.append((idx1, center_idx, idx3))

        # ================================================================
        # 5. Convert bonds and angles from 0-based indices to PDB serial numbers
        # ================================================================
        # bonds and angles currently use 0-based indices, but need PDB serial numbers
        # for Step 2B (mol2_writer) and Step 3
        bonds_serial = []
        for idx1, idx2 in bonds:
            serial1 = index_to_serial.get(idx1)
            serial2 = index_to_serial.get(idx2)
            if serial1 is not None and serial2 is not None:
                bonds_serial.append((serial1, serial2))

        angles_serial = []
        for idx1, idx2, idx3 in angles:
            serial1 = index_to_serial.get(idx1)
            serial2 = index_to_serial.get(idx2)
            serial3 = index_to_serial.get(idx3)
            if serial1 is not None and serial2 is not None and serial3 is not None:
                angles_serial.append((serial1, serial2, serial3))

        # idx_to_coords needs to be converted to serial_to_coords for consistency
        serial_to_coords = {}
        for idx, coords in idx_to_coords.items():
            serial = index_to_serial.get(idx)
            if serial is not None:
                serial_to_coords[serial] = coords

        # Build type-assignment records for the boundary cap atoms, keyed by the
        # SAME coordinate tuples the pre-frcmod generator looks them up by
        # (serial_to_coords values). renamed=False keeps them out of the MASS/NONB
        # sections (they are standard ff14SB types tleap already knows); their
        # original_type == renamed_type so a boundary angle like Y6-C-N resolves to
        # the ordinary O-C-N force-field value.
        boundary_type_assignments = {}
        for cap_serial, info in boundary_cap_atoms.items():
            coords = serial_to_coords.get(cap_serial)
            if coords is None:
                continue
            cap_type = info['renamed_type']
            boundary_type_assignments[coords] = {
                'renamed_type': cap_type,
                'original_type': cap_type,
                'element': info['element'],
                'renamed': False,
                'is_center': False,
                'is_metal_ligand': False,
            }

        return bonds_serial, angles_serial, serial_to_coords, boundary_type_assignments

    def _extract_parameters_from_fingerprint(
        self,
        seminario: 'SeminarioMethod',
        fingerprint_file: Path,
        type_assignments: Dict,
        pdb_file: Path,
        coords
    ) -> Tuple[List, List, List]:
        """
        Extract bonds and angles from fingerprint file and compute parameters.

        Parses the MCPB fingerprint to get connectivity from LINK records,
        then uses the Seminario method to compute force constants.

        The fingerprint file is the authoritative source for connectivity.

        Args:
            seminario: SeminarioMethod instance
            fingerprint_file: Path to standard.fingerprint
            type_assignments: Dict of atom type assignments
            pdb_file: PDB file containing atom coordinates and serial numbers
            coords: Coordinates array from Hessian (3N,) in Bohr

        Returns:
            (bond_parameters, angle_parameters, bonds) tuple
        """
        from .mcpb.seminario import BondParameter, AngleParameter
        import numpy as np

        # ================================================================
        # 1. Parse fingerprint file to get atoms and connectivity
        # ================================================================

        # Fingerprint format:
        # ResID-ResName-AtomName  AtomID  OrigType -> RenamedType
        # LINK AtomID1-AtomName1 AtomID2-AtomName2

        atom_info = {}  # AtomID -> {resid, resname, atom_name, orig_type, renamed_type}
        link_records = []  # List of (AtomID1, AtomID2) pairs

        with open(fingerprint_file) as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith('LINK'):
                    # Parse LINK record: LINK AtomID1-AtomName1 AtomID2-AtomName2
                    parts = line.split()
                    if len(parts) >= 3:
                        # Extract atom IDs from "ID-Name" format
                        atom1_info = parts[1].split('-')
                        atom2_info = parts[2].split('-')

                        try:
                            atom1_id = int(atom1_info[0])
                            atom2_id = int(atom2_info[0])
                            link_records.append((atom1_id, atom2_id))
                        except (ValueError, IndexError):
                            self.console.print(f"[yellow]⚠️  Failed to parse LINK: {line}[/yellow]")

                elif '-' in line and '->' in line:
                    # Parse atom line: ResID-ResName-AtomName  AtomID  OrigType -> RenamedType
                    parts = line.split()
                    if len(parts) >= 5:
                        # Parse residue info
                        res_parts = parts[0].split('-')
                        if len(res_parts) >= 3:
                            resid = res_parts[0]
                            resname = res_parts[1]
                            atom_name = res_parts[2]

                            try:
                                atom_id = int(parts[1])
                                orig_type = parts[2]
                                renamed_type = parts[4]  # After '->'

                                atom_info[atom_id] = {
                                    'resid': resid,
                                    'resname': resname,
                                    'atom_name': atom_name,
                                    'orig_type': orig_type,
                                    'renamed_type': renamed_type.strip()
                                }
                            except (ValueError, IndexError):
                                self.console.print(f"[yellow]⚠️  Failed to parse atom: {line}[/yellow]")

        self.console.print(f"[grey50]Parsed {len(atom_info)} atoms and {len(link_records)} bonds from fingerprint[/grey50]")

        # ================================================================
        # 2. Read PDB file to map PDB serial numbers to array indices
        # ================================================================

        # Parse PDB file to build serial number -> 0-based index mapping
        # This replaces the old XYZ+mapping file approach
        pdb_id_to_xyz_index = {}  # Map PDB serial number -> 0-based index
        xyz_index_to_pdb_id = {}  # Map 0-based index -> PDB serial number (reverse)

        index = 0
        with open(pdb_file) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Parse PDB ATOM record
                    try:
                        serial = int(line[6:11].strip())
                        pdb_id_to_xyz_index[serial] = index
                        xyz_index_to_pdb_id[index] = serial
                        index += 1
                    except (ValueError, IndexError):
                        pass

        self.console.print(f"[grey50]Loaded {len(pdb_id_to_xyz_index)} atoms from PDB file[/grey50]")

        def atomid_to_xyz_index(atom_id):
            """Convert fingerprint AtomID (PDB serial number) to XYZ index."""
            if atom_id not in pdb_id_to_xyz_index:
                raise KeyError(
                    f"Atom ID {atom_id} not found in mapping file.\n"
                    f"Available IDs: {sorted(pdb_id_to_xyz_index.keys())}"
                )
            return pdb_id_to_xyz_index[atom_id]

        def xyz_index_to_atomtype(xyz_idx):
            """Convert XYZ index to atom type from fingerprint."""
            pdb_id = xyz_index_to_pdb_id.get(xyz_idx)
            if pdb_id is None:
                return None
            atom = atom_info.get(pdb_id)
            if atom is None:
                return None
            return atom['renamed_type']

        # ================================================================
        # 3. Build bond list from LINK records
        # ================================================================

        bonds = []  # List of (idx1, idx2, type1, type2)
        metal_atoms = set()  # Set of XYZ indices that are metals

        for atom1_id, atom2_id in link_records:
            # Convert to XYZ indices
            idx1 = atomid_to_xyz_index(atom1_id)
            idx2 = atomid_to_xyz_index(atom2_id)

            # Get atom types
            if atom1_id in atom_info and atom2_id in atom_info:
                type1 = atom_info[atom1_id]['renamed_type']
                type2 = atom_info[atom2_id]['renamed_type']

                bonds.append((idx1, idx2, type1, type2))

                # Track metal atoms (M1-M9, MA-MZ)
                if type1.startswith('M') and len(type1) == 2:
                    metal_atoms.add(idx1)
                if type2.startswith('M') and len(type2) == 2:
                    metal_atoms.add(idx2)

        self.console.print(f"[grey50]Identified {len(bonds)} bonds from LINK records[/grey50]")

        # ================================================================
        # 3b. Extract organic bonds from prmtop
        # ================================================================
        # The prmtop has all bonds except metal-ligand (metals removed before tLEaP).
        # Combined with LINK records, this gives complete connectivity.
        # This replaces distance-based detection which can produce false positives.

        from rich.table import Table

        # Get prmtop path from workspace
        workspace = self.processor._get_workspace() if self.processor else None
        prmtop_path = workspace.get("parm7_file") if workspace else None

        organic_bonds = []
        if prmtop_path and os.path.exists(prmtop_path):
            organic_bonds = self._extract_bonds_from_prmtop(
                prmtop_path, atom_info, pdb_id_to_xyz_index
            )
            # Filter to only keep bonds involving renamed atoms (M*, L*, Y*, Z*, etc.)
            organic_bonds = self._filter_metal_relevant_bonds(organic_bonds)
        else:
            self.console.print("[yellow]⚠️  No prmtop available - falling back to distance detection[/yellow]")
            # Fallback to distance detection if no prmtop
            all_distance_bonds = self._detect_distance_bonds(
                coords, xyz_index_to_pdb_id, atom_info, tolerance=0.40
            )
            organic_bonds = self._filter_metal_relevant_bonds(all_distance_bonds)

        # Merge LINK bonds with organic bonds for display
        # User selects which bonds to compute force constants for
        selected_bonds = self._merge_and_display_bonds(
            bonds, organic_bonds, atom_info
        )

        # Combine LINK bonds + organic bonds for angle detection
        # LINK = metal-ligand bonds, organic = bonds from prmtop
        all_bonds_for_angles = list(bonds)  # Start with LINK bonds
        existing_keys = {tuple(sorted([b[0], b[1]])) for b in bonds}
        for b in organic_bonds:
            key = tuple(sorted([b[0], b[1]]))
            if key not in existing_keys:
                all_bonds_for_angles.append(b)
                existing_keys.add(key)

        # Update metal atoms set from selected bonds
        metal_atoms = set()
        for idx1, idx2, type1, type2 in selected_bonds:
            if type1.startswith('M') and len(type1) == 2:
                metal_atoms.add(idx1)
            if type2.startswith('M') and len(type2) == 2:
                metal_atoms.add(idx2)

        self.console.print(f"[grey50]Final selection: {len(selected_bonds)} bonds for force constants, {len(metal_atoms)} metals[/grey50]")
        self.console.print(f"[grey50]Using {len(all_bonds_for_angles)} bonds for angle detection[/grey50]")

        # ================================================================
        # 4. Build angle list from ALL bonds involving renamed atoms
        # ================================================================

        # Build adjacency list using ALL bonds (not just selected)
        # This ensures angles like CC-Y2-M1 are detected (requires CC-Y2 bond)
        adjacency = {}  # atom_idx -> list of connected atom indices
        for idx1, idx2, _, _ in all_bonds_for_angles:
            if idx1 not in adjacency:
                adjacency[idx1] = []
            if idx2 not in adjacency:
                adjacency[idx2] = []
            adjacency[idx1].append(idx2)
            adjacency[idx2].append(idx1)

        angles = []  # List of (idx1, idx2, idx3, type1, type2, type3)

        # For each atom that has at least 2 connections, find all angle combinations
        for central_idx, neighbors in adjacency.items():
            if len(neighbors) < 2:
                continue

            # Generate all pairs of neighbors to form angles
            for i in range(len(neighbors)):
                for j in range(i + 1, len(neighbors)):
                    idx1 = neighbors[i]
                    idx2 = central_idx
                    idx3 = neighbors[j]

                    # Get atom types using XYZ index to type mapping
                    type1 = xyz_index_to_atomtype(idx1)
                    type2 = xyz_index_to_atomtype(idx2)
                    type3 = xyz_index_to_atomtype(idx3)

                    if type1 and type2 and type3:
                        # Only include angles involving metals
                        if idx1 in metal_atoms or idx2 in metal_atoms or idx3 in metal_atoms:
                            angles.append((idx1, idx2, idx3, type1, type2, type3))

        self.console.print(f"[grey50]Identified {len(angles)} angles involving metal centers[/grey50]")

        # ================================================================
        # 5. Compute force constants using Seminario method
        # ================================================================

        # Compute bond force constants for SELECTED bonds only (metal-ligand bonds)
        bond_params = []
        for atom1_idx, atom2_idx, atom1_type, atom2_type in selected_bonds:
            try:
                param = seminario.compute_bond_force_constant(
                    atom1_idx, atom2_idx, atom1_type, atom2_type
                )
                bond_params.append(param)
            except Exception as e:
                self.console.print(
                    f"[yellow]⚠️  Failed to compute bond {atom1_type}-{atom2_type}: {str(e)}[/yellow]"
                )

        # Compute angle force constants for ALL metal-involving angles
        angle_params = []
        for atom1_idx, atom2_idx, atom3_idx, atom1_type, atom2_type, atom3_type in angles:
            try:
                param = seminario.compute_angle_force_constant(
                    atom1_idx, atom2_idx, atom3_idx,
                    atom1_type, atom2_type, atom3_type
                )
                angle_params.append(param)
            except ValueError as e:
                # Skip degenerate or linear angles (expected in some cases)
                self.console.print(
                    f"[yellow]⚠️  Skipping angle {atom1_type}-{atom2_type}-{atom3_type}: {str(e)}[/yellow]"
                )
            except Exception as e:
                self.console.print(
                    f"[yellow]⚠️  Failed to compute angle {atom1_type}-{atom2_type}-{atom3_type}: {str(e)}[/yellow]"
                )

        return bond_params, angle_params, selected_bonds

    def _get_atom_type_by_fingerprint_id(self, atom_id: int, atom_info: Dict) -> Optional[str]:
        """
        Get atom type by fingerprint AtomID.

        Args:
            atom_id: Fingerprint AtomID (1-indexed)
            atom_info: Dict mapping AtomID to atom information

        Returns:
            Renamed atom type or None
        """
        if atom_id in atom_info:
            return atom_info[atom_id]['renamed_type']
        return None

    def _detect_distance_bonds(self, coords, xyz_index_to_pdb_id: Dict,
                               atom_info: Dict, tolerance: float = 0.40) -> List[Tuple]:
        """
        Detect bonds by distance using covalent radii (MCPB-style).

        Args:
            coords: (3N,) array of coordinates in Bohr
            xyz_index_to_pdb_id: Mapping from XYZ index to PDB serial number
            atom_info: Atom information from fingerprint
            tolerance: Extra distance added to covalent radii sum (Å)

        Returns:
            List of (xyz_idx1, xyz_idx2, type1, type2) tuples
        """
        import numpy as np

        # Covalent radii (Å)
        COVALENT_RADII = {
            'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05,
            'P': 1.07, 'F': 0.57, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39,
            'Fe': 1.32, 'Cu': 1.32, 'Zn': 1.22, 'Ni': 1.24, 'Co': 1.26,
            'Mn': 1.39, 'Mg': 1.41, 'Ca': 1.76, 'Na': 1.66, 'K': 2.03
        }
        B_TO_A = 0.529177249  # Bohr to Angstrom

        n_atoms = len(coords) // 3
        distance_bonds = []

        for i in range(n_atoms):
            pdb_id_i = xyz_index_to_pdb_id.get(i)
            if not pdb_id_i or pdb_id_i not in atom_info:
                continue

            atom_i = atom_info[pdb_id_i]
            element_i = atom_i['atom_name'][0]  # First character of atom name
            radius_i = COVALENT_RADII.get(element_i, 1.5)
            coords_i = coords[3*i:3*i+3]

            for j in range(i+1, n_atoms):
                pdb_id_j = xyz_index_to_pdb_id.get(j)
                if not pdb_id_j or pdb_id_j not in atom_info:
                    continue

                atom_j = atom_info[pdb_id_j]
                element_j = atom_j['atom_name'][0]
                radius_j = COVALENT_RADII.get(element_j, 1.5)
                coords_j = coords[3*j:3*j+3]

                # Calculate distance in Angstrom
                dist_bohr = np.linalg.norm(coords_j - coords_i)
                dist_angstrom = dist_bohr * B_TO_A

                # Check if bonded
                cutoff = radius_i + radius_j + tolerance
                if dist_angstrom <= cutoff:
                    type_i = atom_i['renamed_type']
                    type_j = atom_j['renamed_type']

                    # Skip if either atom is a metal (M* type)
                    # Metal-ligand bonds should only come from LINK records, not distance detection
                    # This matches MCPB.py behavior in get_mc_blist()
                    is_metal_i = type_i.startswith('M') and len(type_i) == 2
                    is_metal_j = type_j.startswith('M') and len(type_j) == 2
                    if is_metal_i or is_metal_j:
                        continue

                    distance_bonds.append((i, j, type_i, type_j))

        return distance_bonds

    def _extract_bonds_from_prmtop(self, prmtop_path: str, atom_info: Dict,
                                    pdb_id_to_xyz_index: Dict) -> List[Tuple]:
        """
        Extract organic bonds from prmtop file for atoms in the small model.

        The prmtop has all bonds except metal-ligand (metals were removed before tLEaP).
        This is more reliable than distance detection which can produce false positives.

        Matching strategy: PDB serial numbers map directly to prmtop indices.
        PDB serial N corresponds to prmtop atom index N-1 (0-based).

        Args:
            prmtop_path: Path to prmtop file
            atom_info: Dict mapping PDB serial -> {resid, resname, atom_name, renamed_type}
            pdb_id_to_xyz_index: Dict mapping PDB serial -> XYZ index

        Returns:
            List of (xyz_idx1, xyz_idx2, type1, type2) tuples for organic bonds
        """
        import parmed

        try:
            parm = parmed.load_file(prmtop_path)
        except Exception as e:
            self.console.print(f"[yellow]⚠️  Could not load prmtop for bond extraction: {e}[/yellow]")
            return []

        # Build set of prmtop indices for atoms in our small model
        # PDB serial N -> prmtop index N-1
        small_model_prmtop_indices = set()
        prmtop_idx_to_pdb_id = {}  # prmtop_idx -> pdb_serial
        for pdb_id in atom_info.keys():
            prmtop_idx = pdb_id - 1  # PDB serials are 1-based, prmtop is 0-based
            if 0 <= prmtop_idx < len(parm.atoms):
                small_model_prmtop_indices.add(prmtop_idx)
                prmtop_idx_to_pdb_id[prmtop_idx] = pdb_id

        prmtop_bonds = []

        n_hh_skipped = 0
        for bond in parm.bonds:
            idx1 = bond.atom1.idx
            idx2 = bond.atom2.idx

            # Skip hydrogen-hydrogen bonds. AMBER rigid-water templates (TIP3P,
            # SPC/E, etc.) carry an explicit H1-H2 bond used only as a SHAKE
            # constraint; it is not a real chemical bond. Reference MCPB.py never
            # sees it because it perceives bonds by covalent radius (H..H ~1.5 A is
            # beyond the ~1.0 A cutoff). If left in, it seeds a bogus HW-HW bond
            # term and spurious O-H-H angles (e.g. YA-HW-HW) in the pre-frcmod.
            if bond.atom1.atomic_number == 1 and bond.atom2.atomic_number == 1:
                n_hh_skipped += 1
                continue

            # Check if both atoms are in our small model
            if idx1 in small_model_prmtop_indices and idx2 in small_model_prmtop_indices:
                pdb_id1 = prmtop_idx_to_pdb_id[idx1]
                pdb_id2 = prmtop_idx_to_pdb_id[idx2]

                xyz_idx1 = pdb_id_to_xyz_index.get(pdb_id1)
                xyz_idx2 = pdb_id_to_xyz_index.get(pdb_id2)

                if xyz_idx1 is not None and xyz_idx2 is not None:
                    type1 = atom_info[pdb_id1]['renamed_type']
                    type2 = atom_info[pdb_id2]['renamed_type']
                    prmtop_bonds.append((xyz_idx1, xyz_idx2, type1, type2))

        hh_note = f" (skipped {n_hh_skipped} H-H rigid-water bond(s))" if n_hh_skipped else ""
        self.console.print(f"[grey50]Extracted {len(prmtop_bonds)} organic bonds from prmtop{hh_note}[/grey50]")
        return prmtop_bonds

    def _filter_metal_relevant_bonds(self, bonds: List[Tuple]) -> List[Tuple]:
        """
        Filter bonds to only those involving systematically renamed atoms.

        Keeps bonds where at least one atom has a renamed type (M*, L*, Y*, Z*, A*, B*, X*).
        These are the only bonds relevant for metal-involving angles.

        Args:
            bonds: List of (idx1, idx2, type1, type2) tuples

        Returns:
            Filtered list of bonds
        """
        # Patterns for systematically renamed atom types
        RENAMED_PREFIXES = {'M', 'L', 'Y', 'Z', 'A', 'B', 'X'}

        filtered = []
        for idx1, idx2, type1, type2 in bonds:
            # Check if either atom type is renamed (starts with M, L, Y, Z, A, B, X)
            type1_renamed = len(type1) == 2 and type1[0] in RENAMED_PREFIXES
            type2_renamed = len(type2) == 2 and type2[0] in RENAMED_PREFIXES

            if type1_renamed or type2_renamed:
                filtered.append((idx1, idx2, type1, type2))

        return filtered

    def _merge_and_display_bonds(self, link_bonds: List[Tuple],
                                 organic_bonds: List[Tuple],
                                 atom_info: Dict) -> List[Tuple]:
        """
        Merge bonds from LINK records and prmtop, show interactive UI.

        Args:
            link_bonds: Metal-ligand bonds from LINK records [(idx1, idx2, type1, type2), ...]
            organic_bonds: Organic bonds from prmtop (involving renamed atoms)
            atom_info: Atom information from fingerprint

        Returns:
            List of selected bonds
        """
        from rich.table import Table

        # Build dictionary of bonds with source tracking
        # All bonds are pre-selected by default - user can deselect if needed
        bond_dict = {}  # key: (min_idx, max_idx), value: {idx1, idx2, type1, type2, source, selected}

        # Add LINK bonds (metal-ligand coordinate bonds)
        for idx1, idx2, type1, type2 in link_bonds:
            key = tuple(sorted([idx1, idx2]))
            bond_dict[key] = {
                'idx1': idx1, 'idx2': idx2,
                'type1': type1, 'type2': type2,
                'source': 'LINK',
                'selected': True
            }

        # Add organic bonds from prmtop (all pre-selected)
        for idx1, idx2, type1, type2 in organic_bonds:
            key = tuple(sorted([idx1, idx2]))
            if key in bond_dict:
                # Already in LINK bonds - mark as both
                bond_dict[key]['source'] = 'Both'
            else:
                # New bond from prmtop - pre-selected
                bond_dict[key] = {
                    'idx1': idx1, 'idx2': idx2,
                    'type1': type1, 'type2': type2,
                    'source': 'Prmtop',
                    'selected': True
                }

        # Convert to list and sort
        all_bonds = sorted(bond_dict.values(), key=lambda b: (b['idx1'], b['idx2']))

        # Count by source
        n_link = sum(1 for b in all_bonds if b['source'] == 'LINK')
        n_prmtop = sum(1 for b in all_bonds if b['source'] == 'Prmtop')
        n_total = len(all_bonds)

        if n_total == 0:
            self.console.print("[yellow]⚠️ No bonds found for Seminario method[/yellow]")
            return []

        # Show table
        self.console.print("\n[bold]Bonds for Seminario Method:[/bold]\n")

        table = Table(title=f"Found {n_total} bonds for force constant calculation")
        table.add_column("✓", justify="center", style="green", width=3)
        table.add_column("Bond", style="cyan")
        table.add_column("Source", justify="center")
        table.add_column("Notes")

        for bond in all_bonds:
            check = "✓" if bond['selected'] else ""
            bond_str = f"{bond['type1']}-{bond['type2']}"
            source_str = bond['source']

            notes = []
            if bond['source'] == 'LINK':
                notes.append("Metal-ligand")
            elif bond['source'] == 'Prmtop':
                notes.append("Organic (renamed atom)")
            elif bond['source'] == 'Both':
                notes.append("Metal-ligand + prmtop")
            notes_str = " ".join(notes)

            table.add_row(check, bond_str, source_str, notes_str)

        self.console.print(table)
        self.console.print(f"\n[grey50]LINK bonds (metal-ligand): {n_link}[/grey50]")
        self.console.print(f"[grey50]Prmtop bonds (organic with renamed atoms): {n_prmtop}[/grey50]")
        self.console.print(f"[cyan]Total selected: {n_total} bonds[/cyan]\n")

        # Offer options - default is to accept all (already selected)
        from proprep.utils.prompts import prompt_with_context
        choice = prompt_with_context(
            self.processor,
            "Options:\n"
            "  [1] Accept all bonds (recommended)\n"
            "  [2] Metal-ligand (LINK) bonds only\n"
            "  [3] Manually toggle individual bonds\n"
            "Choice",
            choices=["1", "2", "3"],
            default="1",
            module="MCPB Step 2", description="Bond selection action",
            options_map={"1": "Accept all", "2": "LINK only", "3": "Manual toggle"}
        )

        if choice == "2":
            # Only keep LINK bonds
            for bond in all_bonds:
                bond['selected'] = bond['source'] in ('LINK', 'Both')
            n_selected = sum(1 for b in all_bonds if b['selected'])
            self.console.print(f"[green]✅ Selected {n_selected} LINK bonds only[/green]")

        elif choice == "3":
            # Manual toggle
            self._manual_toggle_bonds(all_bonds)

        # Return selected bonds as list of tuples
        selected_bonds = [(b['idx1'], b['idx2'], b['type1'], b['type2'])
                         for b in all_bonds if b['selected']]

        return selected_bonds

    def _manual_toggle_bonds(self, bonds: List[Dict]):
        """Manual bond toggle interface."""
        from proprep.utils.prompts import prompt_with_context
        while True:
            self.console.print("\n[cyan]Current selection:[/cyan]")
            for i, bond in enumerate(bonds, 1):
                status = "[green]✓[/green]" if bond['selected'] else "[ ]"
                bond_str = f"{bond['type1']}-{bond['type2']}"
                self.console.print(f"  [{i:2d}] {status} {bond_str} ({bond['source']})")

            choice = prompt_with_context(
                self.processor,
                "\nEnter bond number to toggle, 'all' to select all, 'none' to deselect all, or 'done' to finish",
                default="done",
                module="MCPB Step 2", description="Manual bond toggle"
            )

            if choice.lower() == 'done':
                break
            elif choice.lower() == 'all':
                for bond in bonds:
                    bond['selected'] = True
            elif choice.lower() == 'none':
                for bond in bonds:
                    bond['selected'] = False
            elif choice.isdigit():
                idx = int(choice) - 1
                if 0 <= idx < len(bonds):
                    bonds[idx]['selected'] = not bonds[idx]['selected']
                else:
                    self.console.print("[red]Invalid bond number[/red]")
            else:
                self.console.print("[red]Invalid input[/red]")

    def _find_atom_type_by_coords(self, target_coords: Tuple[float, float, float],
                                  type_assignments: Dict, tolerance: float = 0.001) -> Optional[str]:
        """
        Find atom type by coordinates with tolerance.

        Args:
            target_coords: Target coordinates (x, y, z)
            type_assignments: Dict mapping coordinates to atom data
            tolerance: Distance tolerance in Angstrom

        Returns:
            Atom type string or None
        """
        import numpy as np

        target = np.array(target_coords)

        for coords, assignment in type_assignments.items():
            coords_array = np.array(coords)
            dist = np.linalg.norm(target - coords_array)

            if dist < tolerance:
                return assignment.get('renamed_type') or assignment.get('original_type')

        return None

    def _show_parameter_summary(self, bond_params: List, angle_params: List):
        """
        Display summary table of computed parameters.

        Args:
            bond_params: List of BondParameter objects
            angle_params: List of AngleParameter objects
        """
        from rich.table import Table

        # Bond parameters table
        if bond_params:
            bond_table = Table(title="Bond Parameters", show_header=True, header_style="bold cyan")
            bond_table.add_column("Bond", style="yellow")
            bond_table.add_column("k (kcal/mol/Å²)", justify="right")
            bond_table.add_column("r₀ (Å)", justify="right")
            bond_table.add_column("Std Dev", justify="right", style="grey50")

            for bond in bond_params:
                bond_table.add_row(
                    f"{bond.atom1_type}-{bond.atom2_type}",
                    f"{bond.force_constant:.1f}",
                    f"{bond.eq_length:.4f}",
                    f"{bond.std_dev:.1f}" if bond.std_dev > 0.01 else "-"
                )

            self.console.print(bond_table)

        # Angle parameters table
        if angle_params:
            angle_table = Table(title="Angle Parameters", show_header=True, header_style="bold cyan")
            angle_table.add_column("Angle", style="yellow")
            angle_table.add_column("k (kcal/mol/rad²)", justify="right")
            angle_table.add_column("θ₀ (deg)", justify="right")
            angle_table.add_column("Std Dev", justify="right", style="grey50")

            for angle in angle_params:
                angle_table.add_row(
                    f"{angle.atom1_type}-{angle.atom2_type}-{angle.atom3_type}",
                    f"{angle.force_constant:.2f}",
                    f"{angle.eq_angle:.2f}",
                    f"{angle.std_dev:.2f}" if angle.std_dev > 0.01 else "-"
                )

            self.console.print(angle_table)

    def _save_step2_results(self, output_file: Path, bond_params: List,
                           angle_params: List, hessian_stats: Dict,
                           scale_factor: float):
        """
        Save detailed Step 2 results to JSON file.

        Args:
            output_file: Output JSON file path
            bond_params: List of BondParameter objects
            angle_params: List of AngleParameter objects
            hessian_stats: Hessian statistics dict
            scale_factor: Frequency scaling factor used
        """
        import json

        results = {
            "step": 2,
            "description": "Bonded Parameter Generation (Seminario Method)",
            "scale_factor": scale_factor,
            "hessian_statistics": hessian_stats,
            "bond_parameters": [
                {
                    "atom1_idx": b.atom1_idx,
                    "atom2_idx": b.atom2_idx,
                    "atom1_type": b.atom1_type,
                    "atom2_type": b.atom2_type,
                    "force_constant": b.force_constant,
                    "eq_length": b.eq_length,
                    "std_dev": b.std_dev
                }
                for b in bond_params
            ],
            "angle_parameters": [
                {
                    "atom1_idx": a.atom1_idx,
                    "atom2_idx": a.atom2_idx,
                    "atom3_idx": a.atom3_idx,
                    "atom1_type": a.atom1_type,
                    "atom2_type": a.atom2_type,
                    "atom3_type": a.atom3_type,
                    "force_constant": a.force_constant,
                    "eq_angle": a.eq_angle,
                    "std_dev": a.std_dev
                }
                for a in angle_params
            ]
        }

        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)

        self.console.print(f"[grey50]Saved detailed results: {output_file.name}[/grey50]")

    def _extract_force_field_from_leaprc(self, loaded_leaprcs: List[str]) -> str:
        """
        Extract force field name from selected leaprc files.

        Parses leaprc filenames to identify the base force field (e.g., ff14SB, ff19SB).
        For special leaprcs like conste/constph that source older FFs, reads the file
        to find the sourced force field.

        Args:
            loaded_leaprcs: List of selected leaprc file names

        Returns:
            Force field name (e.g., 'ff14SB', 'ff19SB', 'ff10')
        """
        import re
        import os

        # Priority order: protein force fields are most informative
        protein_leaprcs = [lrc for lrc in loaded_leaprcs if 'protein' in lrc.lower()]
        if protein_leaprcs:
            leaprc_to_parse = protein_leaprcs[0]
        elif loaded_leaprcs:
            leaprc_to_parse = loaded_leaprcs[0]
        else:
            return 'ff14SB'  # Sensible default

        # Extract force field patterns from filename
        # Matches: ff14SB, ff19SB, ff15ipq, fb15, etc.
        patterns = [
            r'ff\d+[A-Za-z]*',  # ff14SB, ff19SB, ff15ipq, ff15ipq-vac
            r'fb\d+',           # fb15
            r'ff\d+\.\w+',      # ff03.r1
        ]

        for pattern in patterns:
            match = re.search(pattern, leaprc_to_parse)
            if match:
                ff_name = match.group(0)
                # Clean up special cases
                if ff_name == 'ff15ipq-vac':
                    return 'ff15ipq-vac'
                return ff_name

        # Special handling for leaprcs that don't have FF in filename (conste, constph, etc.)
        # Read the file to find what it sources
        special_leaprcs = ['conste', 'constph', 'gaff', 'gaff2']
        if any(s in leaprc_to_parse.lower() for s in special_leaprcs):
            amberhome = os.environ.get('AMBERHOME', '')
            if amberhome:
                leaprc_path = Path(amberhome) / 'dat' / 'leap' / 'cmd' / leaprc_to_parse
                if leaprc_path.exists():
                    try:
                        with open(leaprc_path, 'r') as f:
                            content = f.read()
                        # Look for "source" directive
                        source_match = re.search(r'source\s+(\S+)', content)
                        if source_match:
                            sourced_file = source_match.group(1)
                            # Extract FF from sourced file (e.g., oldff/leaprc.ff10 -> ff10)
                            for pattern in patterns:
                                match = re.search(pattern, sourced_file)
                                if match:
                                    return match.group(0)
                    except Exception:
                        pass

            # GAFF doesn't need protein FF, use ff14SB as compatible default
            if 'gaff' in leaprc_to_parse.lower():
                return 'ff14SB'

        # Fallback: try to extract anything that looks like a force field
        self.logger.warning(
            f"Could not extract force field from leaprc: {leaprc_to_parse}, "
            f"using default ff14SB"
        )
        return 'ff14SB'

    def _get_force_field_from_step1(self) -> str:
        """
        Get force field name from step 1 results.

        Returns:
            Force field name (e.g., 'ff14SB', 'ff19SB')
        """
        if "step_1" in self.step_results:
            step1_results = self.step_results["step_1"]
            if "force_field_info" in step1_results:
                return step1_results["force_field_info"]["force_field"]

        # Fallback: prompt user
        from proprep.utils.prompts import prompt_with_context
        self.console.print(
            "\n[yellow]⚠️  Force field not found in step 1 results[/yellow]"
        )
        force_field = prompt_with_context(
            self.processor, "Which force field are you using?",
            choices=["ff14SB", "ff19SB", "ff99SB", "ff10"],
            default="ff14SB",
            module="MCPB Step 2", description="Force field selection fallback"
        )
        return force_field

    def _add_library_parameters(self, frcmod_builder, seminario_bond_params: List,
                                seminario_angle_params: List, bonds_list: List,
                                ff_data: 'ForceFieldData') -> Tuple[int, int]:
        """
        Add force field library parameters for organic bonds/angles.

        Identifies bonds and angles that don't involve metals and adds
        their parameters from the force field library to the frcmod file.

        Args:
            frcmod_builder: FrcmodBuilder instance
            seminario_bond_params: Bond parameters from Seminario method
            seminario_angle_params: Angle parameters from Seminario method
            bonds_list: List of all bonds (from RedoxSite + distance detection)
            ff_data: ForceFieldData with parameter lookups

        Returns:
            Tuple of (n_bonds_added, n_angles_added)
        """
        if ff_data is None:
            self.console.print(
                "[yellow]⚠️  No ForceFieldData available for library parameter lookup[/yellow]"
            )
            return 0, 0

        # Build set of bond/angle pairs already parameterized by Seminario
        seminario_bonds = set()
        for bond in seminario_bond_params:
            type1, type2 = bond.atom1_type, bond.atom2_type
            seminario_bonds.add(tuple(sorted([type1, type2])))

        seminario_angles = set()
        for angle in seminario_angle_params:
            type1, type2, type3 = angle.atom1_type, angle.atom2_type, angle.atom3_type
            # For angles, order matters but we normalize by keeping central atom in middle
            seminario_angles.add((type1, type2, type3))
            seminario_angles.add((type3, type2, type1))

        # ================================================================
        # Add library bond parameters for organic bonds
        # ================================================================
        lib_bond_count = 0

        for idx1, idx2, type1, type2 in bonds_list:
            # Skip if already parameterized by Seminario
            bond_key = tuple(sorted([type1, type2]))
            if bond_key in seminario_bonds:
                continue

            # Skip if involves metal
            if self._is_metal_type(type1) or self._is_metal_type(type2):
                continue

            # Look up library parameter
            lib_param = ff_data.get_bond_parameter(type1, type2)
            if lib_param:
                # Convert library parameter to Seminario BondParameter format
                from .mcpb.seminario import BondParameter

                bond_param = BondParameter(
                    atom1_idx=idx1,
                    atom2_idx=idx2,
                    atom1_type=type1,
                    atom2_type=type2,
                    force_constant=lib_param.force_constant,
                    eq_length=lib_param.eq_length,
                    std_dev=0.0  # Library parameters don't have std dev
                )
                frcmod_builder.add_bond_parameter(bond_param)
                lib_bond_count += 1

                self.logger.debug(
                    f"Added library bond: {type1}-{type2} "
                    f"(k={lib_param.force_constant:.1f}, r0={lib_param.eq_length:.4f})"
                )

        # ================================================================
        # Add library angle parameters for organic angles
        # ================================================================
        lib_angle_count = 0

        # Generate all possible angles from bonds
        # An angle is formed by three bonded atoms: A-B-C where A-B and B-C are bonded
        bond_dict = {}  # atom_idx -> [(neighbor_idx, type1, type2), ...]

        for idx1, idx2, type1, type2 in bonds_list:
            if idx1 not in bond_dict:
                bond_dict[idx1] = []
            if idx2 not in bond_dict:
                bond_dict[idx2] = []

            bond_dict[idx1].append((idx2, type1, type2))
            bond_dict[idx2].append((idx1, type2, type1))

        # Find all angles
        for central_idx, neighbors in bond_dict.items():
            if len(neighbors) < 2:
                continue

            # Get central atom type (consistent across all neighbor entries)
            type2 = neighbors[0][1]

            # Generate all angle combinations
            for i, (idx1, _, neighbor_type1) in enumerate(neighbors):
                for idx3, _, neighbor_type3 in neighbors[i+1:]:

                    # Skip if already parameterized by Seminario
                    angle_key1 = (neighbor_type1, type2, neighbor_type3)
                    angle_key2 = (neighbor_type3, type2, neighbor_type1)
                    if angle_key1 in seminario_angles or angle_key2 in seminario_angles:
                        continue

                    # Skip if involves metal
                    if (self._is_metal_type(neighbor_type1) or self._is_metal_type(type2) or
                        self._is_metal_type(neighbor_type3)):
                        continue

                    # Look up library parameter
                    lib_param = ff_data.get_angle_parameter(neighbor_type1, type2, neighbor_type3)
                    if lib_param:
                        # Convert library parameter to Seminario AngleParameter format
                        from .mcpb.seminario import AngleParameter

                        angle_param = AngleParameter(
                            atom1_idx=idx1,
                            atom2_idx=central_idx,
                            atom3_idx=idx3,
                            atom1_type=neighbor_type1,
                            atom2_type=type2,
                            atom3_type=neighbor_type3,
                            force_constant=lib_param.force_constant,
                            eq_angle=lib_param.eq_angle,
                            std_dev=0.0  # Library parameters don't have std dev
                        )
                        frcmod_builder.add_angle_parameter(angle_param)
                        lib_angle_count += 1

                        self.logger.debug(
                            f"Added library angle: {neighbor_type1}-{type2}-{neighbor_type3} "
                            f"(k={lib_param.force_constant:.2f}, θ0={lib_param.eq_angle:.2f})"
                        )

        return lib_bond_count, lib_angle_count

    def _add_mass_and_nonb_parameters(self, frcmod_builder, type_assignments: Dict,
                                       ff_data: 'ForceFieldData') -> int:
        """
        Add MASS and NONB entries for systematically renamed atoms.

        Reads mass and VDW parameters from force field files (no hardcoding).

        Args:
            frcmod_builder: FrcmodBuilder instance
            type_assignments: Atom type assignments from step 1
            ff_data: ForceFieldData with parameter lookups

        Returns:
            Number of MASS entries added
        """
        if ff_data is None:
            self.console.print(
                "[yellow]⚠️  No ForceFieldData available for MASS/NONB parameters[/yellow]"
            )
            return 0

        # Collect renamed atoms and look up their parameters
        added_count = 0
        renamed_atoms = {}

        for coords, assignment in type_assignments.items():
            if assignment.get('renamed', False):
                renamed_type = assignment['renamed_type']
                original_type = assignment['original_type'].strip()

                # Skip if already processed
                if renamed_type in renamed_atoms:
                    continue

                # Check if this atom has pre-stored VDW parameters (metals from MetalIonDatabase)
                vdw_params = assignment.get('vdw_params')
                if vdw_params:
                    # Use pre-stored parameters from MetalIonDatabase
                    mass = assignment.get('mass') or 0.0  # Handle None values
                    radius = vdw_params.get('radius') or 0.0
                    epsilon = vdw_params.get('epsilon') or 0.0
                    source = vdw_params.get('source', f'Metal {original_type}')

                    # Add MASS entry for renamed type
                    frcmod_builder.add_mass(
                        atom_type=renamed_type,
                        mass=mass,
                        comment=source
                    )

                    # Add NONB entry for renamed type
                    frcmod_builder.add_nonbonded_parameter(
                        atom_type=renamed_type,
                        radius=radius,
                        well_depth=epsilon,
                        comment=source
                    )

                    renamed_atoms[renamed_type] = {
                        'original': original_type,
                        'mass': mass,
                        'radius': radius,
                        'epsilon': epsilon
                    }
                    added_count += 1

                    self.logger.debug(
                        f"Added MASS/NONB (pre-stored): {renamed_type} <- {original_type} "
                        f"(mass={mass:.2f}, R={radius:.4f}, eps={epsilon:.6f})"
                    )

                else:
                    # Look up parameters from force field (for non-metals)
                    mass_param = ff_data.get_mass_parameter(original_type)
                    if not mass_param:
                        self.logger.warning(
                            f"Could not find MASS parameter for original type '{original_type}', "
                            f"skipping {renamed_type}"
                        )
                        continue

                    # Look up nonbonded parameter from original type
                    nonb_param = ff_data.get_nonbonded_parameter(original_type)
                    if not nonb_param:
                        self.logger.warning(
                            f"Could not find NONB parameter for original type '{original_type}', "
                            f"skipping {renamed_type}"
                        )
                        continue

                    # Add MASS entry for renamed type
                    frcmod_builder.add_mass(
                        atom_type=renamed_type,
                        mass=mass_param.mass,
                        comment=mass_param.comment if mass_param.comment else original_type
                    )

                    # Add NONB entry for renamed type
                    frcmod_builder.add_nonbonded_parameter(
                        atom_type=renamed_type,
                        radius=nonb_param.radius,
                        well_depth=nonb_param.well_depth,
                        comment=nonb_param.comment if nonb_param.comment else original_type
                    )

                    renamed_atoms[renamed_type] = {
                        'original': original_type,
                        'mass': mass_param.mass,
                        'radius': nonb_param.radius,
                        'epsilon': nonb_param.well_depth
                    }
                    added_count += 1

                    self.logger.debug(
                        f"Added MASS/NONB: {renamed_type} <- {original_type} "
                        f"(mass={mass_param.mass:.2f}, R={nonb_param.radius:.4f}, "
                        f"eps={nonb_param.well_depth:.6f})"
                    )

        return added_count

    def _add_inherited_dihedral_parameters(self, frcmod_builder, type_assignments: Dict,
                                          param_provider: PrmtopParameterProvider) -> int:
        """
        Add inherited dihedral (proper + improper) parameters for renamed atom types.

        When MCPB renames coordinating atom types (e.g., NA→Y1), the standard FF
        torsions involving those types no longer match. This method copies the
        parent-type torsion parameters from the prmtop and writes Y*-substituted
        versions to the frcmod.

        Args:
            frcmod_builder: FrcmodBuilder to add parameters to
            type_assignments: Dict mapping coords to type assignment info
            param_provider: PrmtopParameterProvider with extracted dihedrals

        Returns:
            Number of dihedral entries added
        """
        from itertools import product

        # Build original→renamed mapping: {"NA": {"Y1", "Y2"}, "nh": {"Y4"}, ...}
        original_to_renamed = {}
        for coords, assignment in type_assignments.items():
            if assignment.get('renamed'):
                orig = assignment['original_type']
                renamed = assignment['renamed_type']
                if orig not in original_to_renamed:
                    original_to_renamed[orig] = set()
                original_to_renamed[orig].add(renamed)

        if not original_to_renamed:
            return 0

        # Track existing dihedrals to avoid duplicates
        existing_dihe = set()
        for d in frcmod_builder.dihedrals:
            t = d['types']
            existing_dihe.add((t[0], t[1], t[2], t[3]))
            existing_dihe.add((t[3], t[2], t[1], t[0]))

        existing_impr = set()
        for d in frcmod_builder.impropers:
            t = d['types']
            existing_impr.add((t[0], t[1], t[2], t[3]))

        added = 0

        # --- Proper dihedrals ---
        for type_quad, terms in param_provider.dihedral_parameters.items():
            t1, t2, t3, t4 = type_quad

            # Check if any type in the quadruple is a parent of a Y* type
            has_parent = any(t in original_to_renamed for t in (t1, t2, t3, t4))
            if not has_parent:
                continue

            # Build substitution options for each position
            sub_options = []
            for t in (t1, t2, t3, t4):
                if t in original_to_renamed:
                    # Include the ORIGINAL type alongside its renamed variants so
                    # PARTIAL substitutions are generated. When a parent type
                    # occurs in more than one position of a term (e.g. both
                    # carboxylate oxygens O2 in the improper 2C-O2-CO-O2) but only
                    # some positions actually ligate the metal, tleap needs the
                    # MIXED term (2C-Y2-CO-O2 / 2C-O2-CO-Y5), not just the fully
                    # substituted 2C-Y2-CO-Y5 that never occurs. The all-original
                    # combo is skipped below, so keeping the original here is safe.
                    sub_options.append([t] + list(original_to_renamed[t]))
                else:
                    sub_options.append([t])

            # Generate all Y*-substituted combinations
            for combo in product(*sub_options):
                # Skip if identical to original (no Y* substitution happened)
                if combo == (t1, t2, t3, t4):
                    continue

                # Skip if already exists
                if combo in existing_dihe or tuple(reversed(combo)) in existing_dihe:
                    continue

                # Write all terms for this dihedral (multi-term support)
                for i, (phi_k, phase, per) in enumerate(terms):
                    # Negative periodicity signals more terms follow
                    pn = -per if i < len(terms) - 1 else per
                    frcmod_builder.add_dihedral(
                        atom_types=combo,
                        idivf=1,
                        pk=phi_k,
                        phase=phase,
                        pn=pn,
                        comment="prmtop"
                    )
                    added += 1

                existing_dihe.add(combo)
                existing_dihe.add(tuple(reversed(combo)))

        # --- Improper dihedrals ---
        for type_quad, terms in param_provider.improper_parameters.items():
            t1, t2, t3, t4 = type_quad

            has_parent = any(t in original_to_renamed for t in (t1, t2, t3, t4))
            if not has_parent:
                continue

            sub_options = []
            for t in (t1, t2, t3, t4):
                if t in original_to_renamed:
                    # Include the ORIGINAL type alongside its renamed variants so
                    # PARTIAL substitutions are generated. When a parent type
                    # occurs in more than one position of a term (e.g. both
                    # carboxylate oxygens O2 in the improper 2C-O2-CO-O2) but only
                    # some positions actually ligate the metal, tleap needs the
                    # MIXED term (2C-Y2-CO-O2 / 2C-O2-CO-Y5), not just the fully
                    # substituted 2C-Y2-CO-Y5 that never occurs. The all-original
                    # combo is skipped below, so keeping the original here is safe.
                    sub_options.append([t] + list(original_to_renamed[t]))
                else:
                    sub_options.append([t])

            for combo in product(*sub_options):
                if combo == (t1, t2, t3, t4):
                    continue
                if combo in existing_impr:
                    continue

                for phi_k, phase, per in terms:
                    frcmod_builder.add_improper(
                        atom_types=combo,
                        pk=phi_k,
                        phase=phase,
                        pn=per,
                        comment="prmtop"
                    )
                    added += 1

                existing_impr.add(combo)

        return added

    def _is_metal_type(self, atom_type: str) -> bool:
        """Check if atom type represents a metal (starts with M and has digit)."""
        return (len(atom_type) == 2 and
                atom_type[0] == 'M' and
                atom_type[1].isdigit())

    def _run_step3(self, residue_name: str, residues: List, output_dir: Path,
                   interactive: bool) -> Dict[str, Any]:
        """
        Run Step 3: RESP Charge Fitting.

        This is a dispatcher that routes to substeps 3A, 3B, 3C, 3D.
        User should call these substeps individually.
        """
        self.console.print("[yellow]Step 3: RESP Charge Fitting[/yellow]")
        self.console.print("[grey50]Please run substeps individually:[/grey50]")
        self.console.print("[grey50]  • Step 3A: ESP Calculation Setup[/grey50]")
        self.console.print("[grey50]  • Step 3B: RESP Input Generation[/grey50]")
        self.console.print("[grey50]  • Step 3C: RESP Execution[/grey50]")
        self.console.print("[grey50]  • Step 3D: Mol2 File Generation[/grey50]")

        return {
            "success": False,
            "message": "Please run Step 3 substeps (3A, 3B, 3C, 3D) individually",
            "status": "use_substeps"
        }

    def _run_step3a(self, residue_name: str, output_dir: Path, interactive: bool = True) -> Dict[str, Any]:
        """
        MCPB Step 3A: Generate Gaussian input for ESP calculation.

        Workflow:
        1. Check Step 1 completion (NOT Step 2B - we use unoptimized large.pdb)
        2. Load large.pdb from Step 1
        3. Generate large_resp.gjf for ESP calculation
        4. Display instructions for user
        5. Save Step 3A state
        """
        self.console.print("\n[bold cyan]═══ MCPB Step 3A: ESP Calculation Setup ═══[/bold cyan]\n")

        # Check prerequisites - ONLY Step 1 required
        if "step_1" not in self.step_results:
            self.console.print("[red]Error: Step 1 must be completed first[/red]")
            return {"success": False, "message": "Step 1 not completed"}

        # Get output directory from models
        step1_dir = self._site_models_dir()

        # Create charge_fit directory
        step3_dir = step1_dir.parent / "charge_fit"
        step3_dir.mkdir(parents=True, exist_ok=True)

        # Locate large.pdb from Step 1 (NOT optimized from Step 2B)
        large_pdb = step1_dir / "large.pdb"
        if not large_pdb.exists():
            large_pdb = Path(self.step_results["step_1"]["output_files"]["large_pdb"])

        if not large_pdb.exists():
            self.console.print(f"[red]Error: large.pdb not found: {large_pdb}[/red]")
            return {"success": False, "message": f"large.pdb not found: {large_pdb}"}

        self.console.print(f"[grey50]Using large model: {large_pdb.name}[/grey50]")

        # Generate ESP input (large_resp.gjf)
        try:
            from proprep.utils.prompts import int_prompt_with_context

            self.console.print("\n[bold]Molecular Charge and Multiplicity[/bold]")
            charge = int_prompt_with_context(
                self.processor, "Total charge", default=0,
                module="MCPB Step 3A", description="ESP total charge"
            )
            multiplicity = int_prompt_with_context(
                self.processor, "Multiplicity", default=1,
                module="MCPB Step 3A", description="ESP multiplicity"
            )

            from .mcpb.qm_interface import GaussianInputGenerator
            qm_gen = GaussianInputGenerator(logger=self.logger)

            esp_input = qm_gen.generate_esp_input(
                pdb_file=str(large_pdb),
                output_file="large_resp.gjf",
                output_dir=str(step3_dir),
                charge=charge,
                multiplicity=multiplicity,
                method="B3LYP",
                basis_set="6-31G*"
            )

            self.console.print(f"\n[green]✓ ESP input: {Path(esp_input).name}[/green]")

            # Display instructions
            from rich.panel import Panel
            from rich.markdown import Markdown

            instructions = f"""
## Running ESP Calculation

Run Gaussian on the ESP input file:

```bash
g16 < large_resp.gjf > large_resp.log
```

Input file: `{Path(esp_input).name}`
Output file: `large_resp.log`

After calculation completes:
- Check large_resp.log for normal termination
- Return to ProPrep
- Run Step 3B to generate RESP inputs
"""
            panel = Panel(
                Markdown(instructions),
                title="[bold cyan]ESP Calculation Instructions[/bold cyan]",
                border_style="cyan",
                expand=False
            )
            self.console.print(panel)

            # Save state
            self.step_results["step_3a"] = {
                "esp_input": esp_input,
                "large_pdb": str(large_pdb),
                "output_dir": str(step3_dir),
                "charge": charge,
                "multiplicity": multiplicity,
                "status": "input_generated"
            }
            self._save_step_results()
            self._save_workflow_state(output_dir)

            self.console.print("\n[green]✓ Step 3A complete: ESP input generated[/green]")
            self.console.print("[yellow]⚠️  Run Gaussian calculation before proceeding to Step 3B[/yellow]")

            return {
                "success": True,
                "step_number": "3a",
                "step_description": "ESP Calculation Setup",
                "next_step": "3b",
                "output_files": {"esp_input": esp_input},
                "message": "Run Gaussian calculation, then proceed to Step 3B"
            }

        except Exception as e:
            import traceback
            self.logger.error(f"Step 3A error: {e}")
            traceback.print_exc()
            return {"success": False, "message": f"Step 3A error: {str(e)}"}

    def _run_step3b(self, residue_name: str, output_dir: Path, interactive: bool = True) -> Dict[str, Any]:
        """
        MCPB Step 3B: Generate RESP input files with restraints.

        Workflow:
        1. Check Step 3A completion and large_resp.log exists
        2. Extract ESP data from Gaussian log → large_resp.esp
        3. Configure charge restraints (interactive)
        4. Generate resp1.in and resp2.in files
        5. Save Step 3B state
        """
        self.console.print("\n[bold cyan]▸ RESP input generation[/bold cyan]")

        # Check prerequisites
        if "step_3a" not in self.step_results:
            self.console.print("[red]Error: Step 3A must be completed first[/red]")
            return {"success": False, "message": "Step 3A not completed"}

        # Get paths
        output_dir = Path(self.step_results["step_3a"]["output_dir"])
        large_pdb = Path(self.step_results["step_3a"]["large_pdb"])

        # Verify large_resp.log exists
        esp_log = output_dir / "large_resp.log"
        if not esp_log.exists():
            self.console.print(f"[red]Error: ESP calculation output not found: {esp_log}[/red]")
            self.console.print("[yellow]Please run Gaussian calculation first (Step 3A)[/yellow]")
            return {"success": False, "message": f"ESP log not found: {esp_log}"}

        try:
            # Extract ESP data
            from .mcpb.esp_extractor import ESPDataExtractor
            esp_extractor = ESPDataExtractor(logger=self.logger)

            self.console.print(f"[grey50]Extracting ESP data from {esp_log.name}...[/grey50]")
            esp_file = esp_extractor.extract_esp_data(
                gaussian_log=str(esp_log),
                output_file="large_resp.esp",
                output_dir=str(output_dir)
            )
            self.console.print(f"[green]✓ ESP file: {Path(esp_file).name}[/green]")

            # Configure restraints (interactive)
            from proprep.utils.prompts import prompt_with_context
            from rich.panel import Panel

            self.console.print("\n[bold]RESP Restraint Configuration[/bold]\n")

            chgmod_help = Panel(
                """[cyan]Backbone Charge Restraints (chgmod):[/cyan]

[bold]0[/bold]: No restraints (all atoms free)
[bold]1[/bold]: Restrain CA, N, C, O (default - recommended)
[bold]2[/bold]: Restrain CA, H, HA, N, C, O
[bold]3[/bold]: Restrain CA, H, HA, N, C, O, CB

Restraints fix backbone charges to force field values to maintain
transferability and prevent charge overfitting.""",
                border_style="blue",
                expand=False
            )
            self.console.print(chgmod_help)

            chgmod_str = prompt_with_context(
                self.processor, "Select backbone restraint level",
                choices=["0", "1", "2", "3"],
                default="1",
                module="MCPB Step 3B", description="RESP backbone restraint level",
                options_map={"0": "No restraints", "1": "CA,N,C,O", "2": "CA,H,HA,N,C,O", "3": "CA,H,HA,N,C,O,CB"}
            )
            chgmod = int(chgmod_str)

            # Get total charge
            total_charge = self.step_results["step_3a"].get("charge", 0)
            self.console.print(f"\n[grey50]Using molecular charge from Step 3A: {total_charge}[/grey50]")

            # Generate RESP inputs
            from .mcpb.resp_input_generator import RESPInputGenerator
            resp_gen = RESPInputGenerator(logger=self.logger)

            # Get ff_data and type_assignments for backbone charge restraints
            # Try step_results first (in-memory), fall back to workspace (resume case)
            ff_data = None
            type_assignments = {}
            if "step_1" in self.step_results:
                ff_data = self.step_results["step_1"].get("force_field_info", {}).get("ff_data")
                # Copied: the prmtop merge below adds keys, and this dict is
                # the live one in step_results.
                type_assignments = dict(
                    self.step_results["step_1"].get("type_assignments", {}))
            workspace = self.processor._get_workspace() if self.processor else None
            if workspace:
                # Merge in preprocessing_atom_data (prmtop charges) ALWAYS, not
                # only when type_assignments is empty. step_1 stores the
                # site-only assignments, so a gap residue bridging the model —
                # the ARG between two coordinating CYS — is absent from them,
                # and its group constraint would fall back to 0 and neutralize
                # a +1 residue that the ESP was computed with.
                #
                # atom_data keys are ('', resid, resname, atom_name) -> {'type', 'charge'}
                # type_assignments needs coord -> {'resname', 'atom_name', 'charge'}
                # Since we don't have coords, use a dummy key — the RESP generator
                # only uses (resname, atom_name) -> charge from this data
                atom_data = workspace.get("preprocessing_atom_data", {})
                if atom_data:
                    for key, data in atom_data.items():
                        if isinstance(key, (tuple, list)) and len(key) >= 4:
                            _, resid, resname, atom_name = key[:4]
                            dummy_key = (resid, resname, atom_name)  # unique enough
                            if dummy_key in type_assignments:
                                continue  # site assignments win
                            type_assignments[dummy_key] = {
                                'resname': resname,
                                'atom_name': atom_name,
                                'charge': data.get('charge'),
                                'type': data.get('type'),
                            }

            # Get RedoxSite residue IDs — only these get free RESP fitting
            # Everything else (ACE, NME, gap fillers) is scaffolding constrained to 0
            redox_site = getattr(self, 'provided_redox_site', None)
            if redox_site is None:
                redox_site = self._get_redox_site_from_workspace()
            metal_site_resids = set()
            if redox_site:
                # Restrained ligands (nonbonded waters) are RESP scaffolding, not
                # parameterized residues: excluding them here makes the RESP input
                # generator constrain them to net-0 (like a cap) and the Mol2Writer
                # skip them, so they stay standard TIP3P instead of a custom
                # Y-typed residue. They remain in large.pdb for correct ESP.
                _, restrained_resids = _collect_restrained_ligands(redox_site)
                # Collect unique resids, adding chain variants to handle
                # tLEaP stripping chain IDs (PDB may have ' ' or '' for chain)
                resids_only = set()
                for atom in redox_site.atoms:
                    if atom.resid in restrained_resids:
                        continue
                    resids_only.add(atom.resid)
                    metal_site_resids.add((atom.chain, atom.resid))
                for center in redox_site.centers:
                    resids_only.add(center.resid)
                    metal_site_resids.add((center.chain, center.resid))
                # Add blank-chain variants for matching tLEaP output
                for resid in resids_only:
                    metal_site_resids.add((' ', resid))
                    metal_site_resids.add(('', resid))

            self.console.print(f"\n[grey50]Generating RESP input files...[/grey50]")
            if metal_site_resids:
                self.console.print(f"[grey50]  Metal-site residues (free fitting): {len(metal_site_resids)}[/grey50]")
                self.console.print(f"[grey50]  Scaffolding residues (constrained to 0): all others[/grey50]")

            # cross_residue_eq_groups is set by the caller (structure_preprocessor's
            # mcpb-3 prompt) when the user has chosen residues to equivalence
            # across ligand positions. Empty list / unset means no cross-residue
            # equivalence is applied.
            cross_groups = getattr(self, "cross_residue_eq_groups", None) or []

            respin1, respin2 = resp_gen.generate_resp_inputs(
                pdb_file=str(large_pdb),
                esp_file=esp_file,
                output_dir=str(output_dir),
                chgmod=chgmod,
                custom_restraints=None,
                total_charge=total_charge,
                ff_collector=ff_data,  # May be None
                type_assignments=type_assignments,  # Backup source for charges
                metal_site_resids=metal_site_resids,  # Only these get free fitting
                cross_residue_equivalence_groups=cross_groups,
            )

            self.console.print(f"[green]✓ RESP input files: resp1.in, resp2.in[/green]")

            # Save state
            self.step_results["step_3b"] = {
                "esp_file": esp_file,
                "respin1": respin1,
                "respin2": respin2,
                "chgmod": chgmod,
                "total_charge": total_charge,
                "metal_site_resids": [list(r) for r in metal_site_resids],  # serialize tuples
                "status": "inputs_ready"
            }
            self._save_step_results()
            self._save_workflow_state(output_dir)

            self.console.print("\n[green]✓ RESP inputs generated[/green]")

            return {
                "success": True,
                "step_number": "3b",
                "step_description": "RESP Input Generation",
                "next_step": "3c",
                "output_files": {"respin1": respin1, "respin2": respin2, "esp_file": esp_file}
            }

        except Exception as e:
            import traceback
            self.logger.error(f"Step 3B error: {e}")
            traceback.print_exc()
            return {"success": False, "message": f"Step 3B error: {str(e)}"}

    def _run_step3c(self, residue_name: str, output_dir: Path) -> Dict[str, Any]:
        """
        MCPB Step 3C: Execute RESP charge fitting.

        Workflow:
        1. Check Step 3B completion
        2. Run Stage 1 RESP
        3. Run Stage 2 RESP (using Stage 1 output)
        4. Parse final charges from resp2.chg
        5. Display charge summary
        6. Save Step 3C state
        """
        self.console.print("\n[bold cyan]▸ RESP execution[/bold cyan]")

        # Check prerequisites
        if "step_3b" not in self.step_results:
            self.console.print("[red]Error: Step 3B must be completed first[/red]")
            return {"success": False, "message": "Step 3B not completed"}

        try:
            from .mcpb.resp_runner import RESPRunner
            from rich.progress import Progress, SpinnerColumn, TextColumn

            # Initialize RESP runner
            resp_runner = RESPRunner(logger=self.logger)

            # Get paths
            output_dir = Path(self.step_results["step_3a"]["output_dir"])
            esp_file = self.step_results["step_3b"]["esp_file"]
            respin1 = self.step_results["step_3b"]["respin1"]
            respin2 = self.step_results["step_3b"]["respin2"]

            # Run RESP fitting
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                console=self.console
            ) as progress:
                task = progress.add_task("Running RESP fitting...", total=None)

                resp_outputs = resp_runner.run_resp_fitting(
                    esp_file=esp_file,
                    respin1=respin1,
                    respin2=respin2,
                    output_dir=str(output_dir)
                )

                progress.update(task, description="✓ RESP fitting complete")

            # Parse final charges
            charges = resp_runner.parse_final_charges(resp_outputs["stage2_charges"])

            self.console.print(f"\n[green]✓ Fitted {len(charges)} charges[/green]")

            # Display charge summary
            from rich.table import Table

            total = sum(charges)

            table = Table(title="RESP Fitted Charges Summary")
            table.add_column("Statistic", style="cyan")
            table.add_column("Value", justify="right", style="green")

            table.add_row("Total atoms", str(len(charges)))
            table.add_row("Total charge", f"{total:+.4f}")
            table.add_row("Min charge", f"{min(charges):+.4f}")
            table.add_row("Max charge", f"{max(charges):+.4f}")
            table.add_row("Mean |charge|", f"{sum(abs(c) for c in charges)/len(charges):.4f}")

            self.console.print("\n")
            self.console.print(table)

            # Save state
            self.step_results["step_3c"] = {
                "resp_outputs": resp_outputs,
                "fitted_charges": charges,
                "status": "fitting_complete"
            }
            self._save_step_results()
            self._save_workflow_state(output_dir)

            self.console.print("\n[green]✓ RESP fitting finished[/green]")

            return {
                "success": True,
                "step_number": "3c",
                "step_description": "RESP Execution",
                "next_step": "3d",
                "output_files": resp_outputs,
                "statistics": {
                    "n_charges": len(charges),
                    "total_charge": total,
                    "min_charge": min(charges),
                    "max_charge": max(charges)
                }
            }

        except Exception as e:
            import traceback
            self.logger.error(f"Step 3C error: {e}")
            traceback.print_exc()
            return {"success": False, "message": f"Step 3C error: {str(e)}"}

    def _run_step3d(self, residue_name: str, output_dir: Path) -> Dict[str, Any]:
        """
        MCPB Step 3D: Generate mol2 files with atom types and RESP charges.

        Workflow:
        1. Check Step 3C completion
        2. Load bond topology from Step 2A
        3. For each residue:
           a. Load fingerprint (atom types)
           b. Load coordinates from large.pdb
           c. Combine with RESP charges
           d. Write mol2 file using saved bond topology
        4. Display mol2 file locations
        5. Save Step 3D state
        """
        self.console.print("\n[bold cyan]▸ Mol2 file generation[/bold cyan]")

        # Check prerequisites
        if "step_3c" not in self.step_results:
            self.console.print("[red]Error: Step 3C must be completed first[/red]")
            return {"success": False, "message": "Step 3C not completed"}

        # Check for bond topology from bonded parameter generation
        step1_dir = self._site_models_dir()
        step2_dir = step1_dir.parent / "bonded_params"
        bond_topology_file = step2_dir / "bond_topology.json"

        if not bond_topology_file.exists():
            self.console.print("[red]Error: bond_topology.json not found from Step 2A[/red]")
            self.console.print("[yellow]Please run Step 2A first to generate bond topology[/yellow]")
            return {"success": False, "message": "bond_topology.json not found from Step 2A"}

        try:
            from .mcpb.mol2_writer import Mol2Writer

            mol2_writer = Mol2Writer(logger=self.logger)

            # Get paths
            output_dir = Path(self.step_results["step_3a"]["output_dir"])
            large_pdb = Path(self.step_results["step_3a"]["large_pdb"])
            # Must be THIS site's fingerprint: the shared step_1 points at
            # whichever site ran last, whose PDB serials do not overlap.
            standard_fp = step1_dir / "standard.fingerprint"
            if not standard_fp.exists():
                standard_fp = Path(
                    self.step_results["step_1"]["output_files"]["standard_fingerprint"])
            charges = self.step_results["step_3c"]["fitted_charges"]

            # Get metal-site residue IDs for mol2 filtering
            metal_site_resids = None
            step3b_data = self.step_results.get("step_3b", {})
            raw_resids = step3b_data.get("metal_site_resids", [])
            if raw_resids:
                metal_site_resids = {tuple(r) for r in raw_resids}

            self.console.print(f"[bold]Generating mol2 files...[/bold]")

            # Generate mol2 files only for metal-site residues
            mol2_files = mol2_writer.write_mol2_files(
                output_dir=str(output_dir),
                pdb_file=str(large_pdb),
                fingerprint_file=str(standard_fp),
                resp_charges=charges,
                bond_topology_file=str(bond_topology_file),
                metal_site_resids=metal_site_resids
            )

            # Display summary
            from rich.table import Table
            from rich.panel import Panel

            table = Table(title="Generated Mol2 Files")
            table.add_column("Residue", style="cyan")
            table.add_column("File", style="green")
            table.add_column("Net Charge", justify="right", style="yellow")

            total_charge = 0.0
            for residue, mol2_path in mol2_files.items():
                # Sum charges from mol2 atom section
                net_charge = 0.0
                try:
                    in_atoms = False
                    with open(mol2_path) as mf:
                        for line in mf:
                            if '@<TRIPOS>ATOM' in line:
                                in_atoms = True
                                continue
                            if '@<TRIPOS>' in line and in_atoms:
                                break
                            if in_atoms and line.strip():
                                parts = line.split()
                                if len(parts) >= 9:
                                    net_charge += float(parts[8])
                except Exception:
                    pass
                table.add_row(residue, Path(mol2_path).name, f"{net_charge:+.4f}")
                total_charge += net_charge

            table.add_section()
            table.add_row("[bold]Total[/bold]", "", f"[bold]{total_charge:+.4f}[/bold]")

            self.console.print("\n")
            self.console.print(table)

            panel = Panel(
                f"""[green]✓ Mol2 files ready for AMBER parameter generation[/green]

Files location: [cyan]{output_dir}[/cyan]

These files contain:
• Atom names and coordinates
• MCPB atom types (from fingerprints)
• RESP-fitted charges
• Bond connectivity (from Step 2A topology)

[bold]Next step:[/bold] run the [cyan]Force Field Integration[/cyan] checklist step.
It assembles these mol2 files with the bonded .frcmod already produced in
Step 2A into a reusable FF library and registers them for tLEaP — no manual
antechamber/parmchk2 needed. After that, build the system in the Topology
Generator.
""",
                border_style="green",
                title="[bold]Mol2 Generation Complete[/bold]",
                expand=False
            )
            self.console.print(panel)

            # Save state
            self.step_results["step_3d"] = {
                "mol2_files": mol2_files,
                "status": "complete"
            }
            self._save_step_results()
            self._save_workflow_state(output_dir)

            self.console.print("\n[green]✓ Step 3D complete: Mol2 files generated[/green]")
            self.console.print("[bold cyan]→ RESP charge fitting workflow complete![/bold cyan]")

            return {
                "success": True,
                "step_number": "3d",
                "step_description": "Mol2 File Generation",
                "output_files": mol2_files,
                "statistics": {
                    "n_mol2_files": len(mol2_files),
                    "residues": list(mol2_files.keys())
                }
            }

        except Exception as e:
            import traceback
            self.logger.error(f"Step 3D error: {e}")
            traceback.print_exc()
            return {"success": False, "message": f"Step 3D error: {str(e)}"}
    
    def _run_step4(self, residue_name: str, residues: List, output_dir: Path,
                   interactive: bool) -> Dict[str, Any]:
        """Run Step 4: Simulation Setup (Future Implementation)."""
        
        self.console.print("[yellow]Step 4: Simulation Setup[/yellow]")
        self.console.print("[grey50]This step will be implemented in the future[/grey50]")
        
        return {
            "success": False,
            "message": "Step 4 not yet implemented",
            "status": "not_implemented"
        }
    
    def _get_complete_structure_from_workspace(self) -> Optional[Structure.Structure]:
        """Get the best available structure from workspace with completeness warnings.

        Uses StructureSelector for consistent structure access with preference for
        repaired > filtered > original structures.
        """
        if not self.processor:
            self.console.print("[red]❌ Processor not available[/red]")
            return None

        # Get workspace using ProPrep convention
        workspace = self.processor._get_workspace()

        # Use StructureSelector for consistent structure access
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console, processor=self.processor)

        # Try repaired structure first
        repaired_structure = selector.get_structure_by_key(
            "repaired_structure", require_exists=False
        )
        if repaired_structure:
            self.console.print("[green]✅ Using repaired structure from workspace[/green]")
            return repaired_structure

        # Try filtered structure
        filtered_structure = selector.get_structure_by_key(
            "filtered_structure", require_exists=False
        )
        if filtered_structure:
            self.console.print("[yellow]⚠️  Using filtered structure from workspace[/yellow]")
            self.console.print("[yellow]Please ensure this structure is complete:[/yellow]")
            self.console.print("[yellow]  • No missing atoms or residues[/yellow]")
            self.console.print("[yellow]  • Hydrogens added[/yellow]")
            self.console.print("[yellow]  • Alternate locations resolved[/yellow]")

            from proprep.utils.prompts import confirm_with_context
            if not confirm_with_context(
                self.processor, "Is the filtered structure complete?", default=False,
                module="Metal Site Parameterizer", description="Confirm filtered structure completeness"
            ):
                self.console.print("[yellow]Consider running Structure Completeness module first[/yellow]")
                if not confirm_with_context(
                    self.processor, "Continue anyway?", default=False,
                    module="Metal Site Parameterizer", description="Continue with incomplete structure"
                ):
                    return None

            return filtered_structure

        # Try original/any structure
        original_structure = selector.get_structure_object(silent=True)
        if original_structure:
            self.console.print("[red]⚠️  Using original structure from workspace[/red]")
            self.console.print("[red]Original structures often have missing atoms/residues/hydrogens[/red]")

            if not confirm_with_context(
                self.processor, "Continue with potentially incomplete structure?", default=False,
                module="Metal Site Parameterizer", description="Continue with original structure"
            ):
                self.console.print("[cyan]Recommendation: Run Structure Completeness module first[/cyan]")
                return None
            
            return original_structure
        else:
            self.console.print("[red]❌ No structure found in workspace[/red]")
            return None
        
    # ═══════════════════════════════════════════════════════════════════
    # State persistence (for cross-step data, read back by structure_preprocessor)
    # ═══════════════════════════════════════════════════════════════════

    def _save_workflow_state(self, output_dir: Path):
        """Save step_results to JSON for cross-step coordination."""
        state = {
            "step_results": self._serialize_step_results(),
            "save_time": datetime.now().isoformat(),
        }
        state_file = Path(output_dir) / "workflow_state.json"
        try:
            with open(state_file, 'w') as f:
                json.dump(state, f, indent=2)
        except Exception as e:
            self.console.print(f"[yellow]Warning: Could not save workflow state: {e}[/yellow]")

    def _serialize_step_results(self) -> Dict:
        """Serialize step results for JSON storage."""
        serialized = {}
        for step_key, results in self.step_results.items():
            serialized[step_key] = self._serialize_dict(results)
        return serialized

    def _serialize_dict(self, d: Dict) -> Dict:
        """Recursively serialize dictionary values."""
        serialized = {}
        for key, value in d.items():
            if isinstance(value, Path):
                serialized[key] = str(value)
            elif isinstance(value, dict):
                serialized[key] = self._serialize_dict(value)
            elif isinstance(value, (str, int, float, bool, list)):
                serialized[key] = value
        return serialized


    def _build_small_model_interactive(self, redox_site: RedoxSite, pdb_file: str,
                                       output_dir: Path, interactive: bool) -> Tuple[SmallModelBuilder, List[Tuple[str, int]]]:
        """
        Build small model with interactive residue selection.

        Args:
            redox_site: RedoxSite object
            pdb_file: Path to PDB file
            output_dir: Output directory
            interactive: Enable user prompts

        Returns:
            (SmallModelBuilder, list of selected (chain, resid) tuples)
        """
        small_model = SmallModelBuilder(redox_site, pdb_file, self.console)

        # Get coordinating residues
        coordinating = small_model.get_coordinating_residues()
        all_residues = small_model.get_all_residues_in_site()

        # Show selection interface (Option B format)
        self.console.print("\n[bold cyan]═══ Small Model Residue Selection ═══[/bold cyan]")
        self.console.print("[grey50]Small model: Sidechain-only for QM bonded parameter calculation[/grey50]\n")

        # Get metal center residues
        metal_residues = set()
        for center in redox_site.centers:
            # Find which residue this metal belongs to
            for atom in redox_site.atoms:
                if atom.coords == center.coords:
                    metal_residues.add((atom.chain, atom.resid))
                    break

        # Show coordinating residues
        self.console.print("[yellow]Residues with metal coordination:[/yellow]")
        selected = []
        for chain, resid, resname, coord_atoms in coordinating:
            atoms_str = ", ".join(coord_atoms)
            self.console.print(f"  \\[x] {chain}:{resid} {resname:3s}  ({atoms_str} coordinate to metal)")
            selected.append((chain, resid))

        # Always include metal ions
        for chain, resid in metal_residues:
            if (chain, resid) not in selected:
                # Find residue name
                for atom in redox_site.atoms:
                    if atom.chain == chain and atom.resid == resid:
                        self.console.print(f"  \\[x] {chain}:{resid} {atom.resname:3s}  (metal center)")
                        selected.append((chain, resid))
                        break

        # Show additional residues (excluding metals)
        additional = [(c, r, n) for c, r, n in all_residues
                     if (c, r) not in [(c2, r2) for c2, r2, _, _ in coordinating]
                     and (c, r) not in metal_residues]

        if additional and interactive:
            self.console.print("\n[yellow]Additional residues in RedoxSite:[/yellow]")
            for chain, resid, resname in additional:
                self.console.print(f"  [ ] {chain}:{resid} {resname:3s}  (neighbor)")

            from proprep.utils.prompts import confirm_with_context
            if confirm_with_context(
                self.processor, "\nInclude additional residues in small model?", default=False,
                module="MCPB Step 1", description="Include additional residues in small model"
            ):
                # TODO: Let user select specific additional residues
                # For now, include all
                for chain, resid, resname in additional:
                    selected.append((chain, resid))

        # Bridging-residue choice: the small model must fill 1-residue peptide
        # gaps (else the flanking sidechains' caps collide), but the bridge is
        # only backbone scaffolding. If that intervening residue is charged/bulky
        # (e.g. an ARG between two coordinating CYS) it perturbs the QM and the
        # model's net charge, so offer a neutral GLY bridge — matching the large
        # model's default. Only ask when such a gap actually exists.
        use_gly = True
        if interactive and small_model.has_single_residue_gaps(selected):
            from proprep.utils.prompts import confirm_with_context

            # Name the residue GLY would displace and what that does to the
            # charge. The small model exists for the Seminario force constants,
            # so a bridging residue is there to keep the backbone contiguous --
            # answering no pulls in its full sidechain and its charge.
            gap_note = "Use GLY to bridge small-model 1-residue gaps"
            gap_preview = small_model.preview_gap_residues(selected, max_gap=1)
            if gap_preview:
                self.console.print("\n[bold]Residues in these gaps:[/bold]")
                for chain, resid, resname, charge in gap_preview:
                    suffix = f" (formal charge {charge:+d})" if charge else ""
                    style = "yellow" if charge else "white"
                    self.console.print(
                        f"  {chain}:{resid} [{style}]{resname}[/{style}]{suffix}")

                delta = sum(c for *_rest, c in gap_preview)
                if delta:
                    charged = ", ".join(f"{ch}:{ri} {rn} ({c:+d})"
                                        for ch, ri, rn, c in gap_preview if c)
                    self.console.print(
                        f"  [yellow]⚠ Answering no keeps {charged} in the small "
                        f"model, changing its net charge by {delta:+d} and "
                        f"adding a flexible sidechain to the frequency "
                        f"calculation.[/yellow]")
                    gap_note = (f"Use GLY to bridge small-model 1-residue gaps; "
                                f"keeping {charged} changes the model charge by "
                                f"{delta:+d}")

            use_gly = confirm_with_context(
                self.processor,
                "Bridge 1-residue gaps with neutral GLY (vs the actual PDB residue)?",
                default=True,
                module="MCPB Step 1",
                description=gap_note,
            )

        # Build model
        small_model.build_from_residues(selected, use_gly=use_gly)
        small_model.show_summary()

        return small_model, selected

    def _build_large_model_interactive(self, redox_site: RedoxSite, pdb_file: str,
                                       output_dir: Path, interactive: bool,
                                       small_residues: List[Tuple[str, int]]) -> Tuple[LargeModelBuilder, List[Tuple[str, int]]]:
        """
        Build large model with interactive residue selection and gap filling.

        Args:
            redox_site: RedoxSite object
            pdb_file: Path to PDB file
            output_dir: Output directory
            interactive: Enable user prompts
            small_residues: Residues selected for small model

        Returns:
            (LargeModelBuilder, list of selected (chain, resid) tuples)
        """
        large_model = LargeModelBuilder(redox_site, pdb_file, self.console)

        # Show selection interface
        self.console.print("\n[bold cyan]═══ Large Model Residue Selection ═══[/bold cyan]")
        self.console.print("[grey50]Large model: Full residues for RESP charge fitting[/grey50]\n")

        # Default: same as small model
        selected = list(small_residues)

        # Get coordinating residues
        coordinating = large_model.get_coordinating_residues()

        self.console.print("[yellow]Selected residues (from small model):[/yellow]")
        for chain, resid in selected:
            # Find residue info
            resname = "???"
            for c, r, n, atoms in coordinating:
                if c == chain and r == resid:
                    resname = n
                    atoms_str = ", ".join(atoms)
                    self.console.print(f"  \\[x] {chain}:{resid} {resname:3s}  ({atoms_str})")
                    break

        # Gap filling options
        if interactive:
            from proprep.utils.prompts import confirm_with_context, int_prompt_with_context
            self.console.print("\n[bold]Gap Filling Options:[/bold]")
            if confirm_with_context(
                self.processor, "Fill gaps between selected residues?", default=True,
                module="MCPB Step 1", description="Enable gap filling in large model"
            ):
                max_gap = int_prompt_with_context(
                    self.processor, "Fill gaps up to how many residues?", default=5,
                    module="MCPB Step 1", description="Maximum gap size"
                )

                # Name the residues GLY would displace, and what that does to
                # the model's charge. Substituting a charged residue shifts the
                # total, which then propagates into the QM charge and the RESP
                # constraint; without this the two answers look interchangeable.
                gap_preview = large_model.preview_gap_residues(selected, max_gap)
                gly_note = "Use GLY for gap filling"
                if gap_preview:
                    self.console.print("\n[bold]Residues in these gaps:[/bold]")
                    for chain, resid, resname, charge in gap_preview:
                        if charge:
                            self.console.print(
                                f"  {chain}:{resid} [yellow]{resname}[/yellow] "
                                f"(formal charge {charge:+d})")
                        else:
                            self.console.print(f"  {chain}:{resid} {resname}")

                    delta = sum(c for *_rest, c in gap_preview)
                    if delta:
                        charged = ", ".join(
                            f"{ch}:{ri} {rn} ({c:+d})"
                            for ch, ri, rn, c in gap_preview if c)
                        self.console.print(
                            f"  [yellow]⚠ Answering yes replaces {charged} with "
                            f"neutral GLY, changing the large model's net charge "
                            f"by {-delta:+d}.[/yellow]")
                        gly_note = (f"Use GLY for gap filling; substituting "
                                    f"{charged} changes the model charge by "
                                    f"{-delta:+d}")

                use_gly = confirm_with_context(
                    self.processor, "Use GLY for gap filling (vs actual PDB residues)?", default=True,
                    module="MCPB Step 1", description=gly_note
                )

                # Build with gap filling
                large_model.build_from_residues(selected, max_gap=max_gap, use_gly=use_gly)
            else:
                large_model.build_from_residues(selected, max_gap=0)
        else:
            # Non-interactive: use defaults
            large_model.build_from_residues(selected, max_gap=5, use_gly=True)

        large_model.show_summary()

        return large_model, selected

    def _add_gap_residue_atoms_to_assignments(
        self,
        model: ModelBuilder,
        redox_site: 'RedoxSite',
        preprocessing_types: Dict,
        type_assignments: Dict,
        pdb_file: str
    ):
        """
        Add gap residue atoms to type_assignments from preprocessing data.

        Gap residues (like PHE between two HIE) are included in the large model
        but aren't part of the RedoxSite. Their charges come from preprocessing
        (extracted from prmtop) but were filtered out when converting to
        RedoxSite-only type_assignments.

        Args:
            model: ModelBuilder with model_residues (including gap residues)
            redox_site: Original RedoxSite (to identify which residues are gaps)
            preprocessing_types: Full preprocessing type_assignments (unfiltered)
            type_assignments: Dictionary to update with gap residue atom data
            pdb_file: Path to prepared PDB for atom lookup
        """
        from Bio.PDB import PDBParser
        from .model_builder import CapType

        # Get residue keys in RedoxSite
        site_residue_keys = {(atom.chain, atom.resid) for atom in redox_site.atoms}

        # Load PDB for coordinate lookup
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("prepared", pdb_file)

        # Build residue map from PDB
        pdb_residues = {}
        for model_struct in structure:
            for chain in model_struct:
                for residue in chain:
                    resid = residue.get_id()[1]
                    key = (chain.id, resid)
                    pdb_residues[key] = residue

        added_count = 0
        for model_residue in model.model_residues:
            # Skip capping groups (handled separately), but include gap residues (CapType.FULL)
            if model_residue.cap_type not in (CapType.NONE, CapType.FULL):
                continue

            # Check if this is a gap residue (not in original RedoxSite)
            res_key = (model_residue.chain, model_residue.resid)
            if res_key in site_residue_keys:
                continue  # Already in type_assignments

            # Get atoms from PDB
            pdb_residue = pdb_residues.get(res_key)
            if pdb_residue is None:
                continue

            # Add each atom from the gap residue
            for atom in pdb_residue.get_atoms():
                coords = tuple(round(x, 3) for x in atom.get_coord())

                # Skip if already in type_assignments
                if coords in type_assignments:
                    continue

                # Look up in preprocessing_types
                if coords in preprocessing_types:
                    assignment = preprocessing_types[coords]
                    # Handle both AtomTypeAssignment dataclass and dict
                    if hasattr(assignment, 'charge'):
                        charge = assignment.charge if assignment.charge is not None else 0.0
                        orig_type = assignment.original_type or atom.element
                    else:
                        charge = assignment.get('charge')
                        charge = charge if charge is not None else 0.0
                        orig_type = assignment.get('original_type', atom.element)

                    type_assignments[coords] = {
                        'atom_name': atom.get_id(),
                        'resname': model_residue.resname,
                        'resid': model_residue.resid,
                        'chain': model_residue.chain,
                        'element': atom.element if atom.element else '',
                        'original_type': orig_type,
                        'renamed_type': orig_type,
                        'charge': charge
                    }
                    added_count += 1

        if added_count > 0:
            self.console.print(f"[grey50]  Added {added_count} gap residue atoms to type assignments[/grey50]")

    def _add_capping_atoms_to_assignments(self, model: ModelBuilder, pdb_writer: PDBWriter,
                                          type_assignments: Dict):
        """
        Add capping group atoms to type_assignments for fingerprint generation.

        Capping groups (ACE, NME, GLY) have synthetic coordinates that aren't in
        the original type_assignments. We need to add them with dummy atom types.

        Args:
            model: ModelBuilder with model_residues
            pdb_writer: PDBWriter for getting capping group atoms
            type_assignments: Dictionary to update with capping atom data
        """
        from .model_builder import CapType

        for model_residue in model.model_residues:
            if model_residue.cap_type != CapType.NONE:
                # Get capping group atoms (same as used in XYZ generation)
                cap_atoms_data = pdb_writer._get_capping_group_atoms(model_residue, {})

                # Determine proper residue name for capping group
                if model_residue.cap_type == CapType.ACE:
                    cap_resname = 'ACE'
                elif model_residue.cap_type == CapType.NME:
                    cap_resname = 'NME'
                elif model_residue.cap_type == CapType.GLY:
                    cap_resname = 'GLY'
                else:
                    cap_resname = model_residue.resname

                # Add each atom to type_assignments
                for atom_name, element, x, y, z, charge, serial in cap_atoms_data:
                    coords = (x, y, z)
                    # Create assignment entry with proper atom names and residue names
                    type_assignments[coords] = {
                        'atom_name': atom_name,
                        'resname': cap_resname,
                        'resid': model_residue.resid,
                        'chain': model_residue.chain,
                        'element': element,
                        'original_type': element.upper(),  # Generic: C, H, N, O
                        'renamed_type': element.upper(),   # No renaming for caps
                        'charge': charge
                    }

    def _create_site_from_model(self, model: ModelBuilder, original_site: RedoxSite,
                                pdb_writer: PDBWriter) -> RedoxSite:
        """
        Create a RedoxSite from model that includes capping group atoms.

        This ensures fingerprint files match XYZ files atom-for-atom.

        Args:
            model: ModelBuilder with model_residues (including caps)
            original_site: Original RedoxSite
            pdb_writer: PDBWriter for accessing PDB structure

        Returns:
            New RedoxSite with all atoms from model (including synthetic caps)
        """
        from proprep.structure_prep.comprehensive_redox_detector import RedoxSite, METALS, RedoxSiteAtom
        from .model_builder import CapType

        model_atoms = []

        # Process each residue in model
        for model_residue in model.model_residues:
            if model_residue.cap_type != CapType.NONE:
                # Capping group - need to generate atoms from coordinates
                # Get the capping group atoms using same logic as XYZ writer
                cap_atoms_data = pdb_writer._get_capping_group_atoms(model_residue, {})

                # Determine proper residue name for capping group
                if model_residue.cap_type == CapType.ACE:
                    cap_resname = 'ACE'
                elif model_residue.cap_type == CapType.NME:
                    cap_resname = 'NME'
                elif model_residue.cap_type == CapType.GLY:
                    cap_resname = 'GLY'
                else:
                    cap_resname = model_residue.resname

                # Convert to RedoxSiteAtom objects
                for atom_name, element, x, y, z, charge, serial in cap_atoms_data:
                    # Create atom with proper residue info
                    # Note: ACE is at resid-1, NME is at resid+1 (already in model_residue.resid)
                    atom = RedoxSiteAtom(
                        coords=(x, y, z),
                        element=element,
                        atom_name=atom_name,
                        resname=cap_resname,
                        resid=model_residue.resid,
                        chain=model_residue.chain,
                        properties={'serial_number': serial}
                    )
                    model_atoms.append(atom)
            else:
                # Regular residue from PDB - get atoms from original site
                # ALWAYS update serial numbers from BioPython structure (authoritative source)
                residue = pdb_writer.residue_map.get((model_residue.chain, model_residue.resid))

                for atom in original_site.atoms:
                    if atom.chain == model_residue.chain and atom.resid == model_residue.resid:
                        # Update serial number from BioPython structure (always, to ensure correctness)
                        if residue:
                            # Find matching BioPython atom by name
                            for bio_atom in residue.get_atoms():
                                if bio_atom.get_id() == atom.atom_name:
                                    # Ensure properties dict exists
                                    if not hasattr(atom, 'properties'):
                                        atom.properties = {}
                                    # Always update to correct PDB serial number
                                    atom.properties['serial_number'] = bio_atom.serial_number
                                    break

                        model_atoms.append(atom)

        # Build set of atom coordinates for filtering bonds
        atom_coords = {atom.coords for atom in model_atoms}

        # Filter bonds to only those where both atoms are in the model
        filtered_bonds = []
        for bond in original_site.bonds:
            if bond.atom1_coords in atom_coords and bond.atom2_coords in atom_coords:
                filtered_bonds.append(bond)

        # Filter centers to only those in the model
        filtered_centers = []
        for center in original_site.centers:
            if center.coords in atom_coords:
                filtered_centers.append(center)

        # Create new site with proper constructor
        new_site = RedoxSite(
            site_id=original_site.site_id,
            structure_id=original_site.structure_id
        )

        # Add atoms
        new_site.atoms = model_atoms

        # Add centers
        new_site.centers = filtered_centers

        # Add bonds
        new_site.bonds = filtered_bonds

        return new_site

    def _save_bond_topology(self,
                            bonds: List[Tuple[int, int]],
                            angles: List[Tuple[int, int, int]],
                            serial_to_coords: Dict[int, Tuple[float, float, float]],
                            output_file: Path) -> None:
        """
        Save bond topology to JSON for reuse in Step 2B and Step 3D.

        Args:
            bonds: List of (serial1, serial2) tuples (PDB serial numbers)
            angles: List of (serial1, serial2, serial3) tuples (PDB serial numbers)
            serial_to_coords: Dict mapping PDB serial number → (x, y, z) coordinates
            output_file: Path to output JSON file
        """
        import json

        topology_data = {
            "bonds": bonds,
            "angles": angles,
            "serial_to_coords": {int(k): list(v) for k, v in serial_to_coords.items()}
        }

        with open(output_file, 'w') as f:
            json.dump(topology_data, f, indent=2)

        self.logger.debug(f"Saved bond topology: {len(bonds)} bonds, {len(angles)} angles")


    def _load_bond_topology(self, topology_file: Path) -> Tuple[List, List, Dict]:
        """
        Load bond topology from JSON.

        Args:
            topology_file: Path to bond_topology.json

        Returns:
            Tuple of (bonds, angles, serial_to_coords) where all use PDB serial numbers
        """
        import json

        with open(topology_file) as f:
            data = json.load(f)

        bonds = [tuple(b) for b in data["bonds"]]
        angles = [tuple(a) for a in data["angles"]]
        # Support both old key name (idx_to_coords) and new key name (serial_to_coords)
        coords_key = "serial_to_coords" if "serial_to_coords" in data else "idx_to_coords"
        serial_to_coords = {int(k): tuple(v) for k, v in data[coords_key].items()}

        return bonds, angles, serial_to_coords


