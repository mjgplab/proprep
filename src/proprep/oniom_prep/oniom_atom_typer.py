"""
ONIOM Atom Typer

Atom typing for FULL structures (not just RedoxSite) for ONIOM calculations.
Adapts MCPB 3-tier workflow to handle complete protein structures.
"""

import logging
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass

from Bio.PDB import PDBParser, Structure
from rich.console import Console
from rich.table import Table

from proprep.forcefield_prep.forcefield_data_collector import ForcefieldDataCollector
from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
    float_prompt_with_context,
)
from proprep.utils.file_browser import remap_recorded_index_by_key, annotate_recorded_key
from proprep.forcefield_prep.mcpb.terminal_classifier import TerminalClassifier, TerminalType
from proprep.structure_prep.comprehensive_redox_detector import RedoxSite


@dataclass
class AtomTypeResult:
    """Result of atom typing for a single atom"""
    atom_type: str
    charge: float
    source: str  # "force_field", "metal_manual", "ligand_manual"


class ONIOMAtomTyper:
    """
    3-Tier atom typing for complete protein structures (ONIOM-specific).

    Workflow:
    1. Tier 1: Metal atoms (manual charge/spin input)
    2. Tier 2: Standard organic residues (from force field library)
    3. Tier 3: Non-standard ligands (manual or antechamber)
    """

    def __init__(self,
                 structure: Structure,
                 redox_site: RedoxSite,
                 pdb_file: str,
                 ff_collector: ForcefieldDataCollector,
                 console: Console,
                 logger: Optional[logging.Logger] = None,
                 processor=None):
        """
        Initialize ONIOM atom typer.

        Args:
            structure: BioPython Structure object (full protein)
            redox_site: RedoxSite object (for identifying metals/ligands)
            pdb_file: Path to PDB file (for terminal classification)
            ff_collector: Force field data collector with loaded libraries
            console: Rich console for output
            logger: Optional logger
        """
        self.structure = structure
        self.redox_site = redox_site
        self.pdb_file = pdb_file
        self.ff_collector = ff_collector
        self.console = console
        self.logger = logger or logging.getLogger(__name__)
        self.processor = processor

        # Results storage
        self.atom_types: Dict[Tuple[float, float, float], str] = {}
        self.charges: Dict[Tuple[float, float, float], float] = {}

        # Classification caches
        self.terminal_map: Dict[Tuple[str, int, str], TerminalType] = {}
        self.metal_residues: List[Tuple[str, int, str, str]] = []  # (chain, resid, ins, atom_name)
        self.ligand_atoms: List[Tuple[float, float, float]] = []
        self.standard_residues: set = set()

        # Track custom forcefield files used during atom typing
        self.used_frcmod_files: List[str] = []  # Paths to .frcmod files used

        # Track lib parsers used during atom typing (for connectivity building)
        from proprep.oniom_prep.amber_lib_parser import AmberLibParser
        self.lib_parsers_used: Dict[str, AmberLibParser] = {}  # {resname: AmberLibParser}

        # Track which resnames have been checked for cap atom charge redistribution
        self._cap_checked_resnames: set = set()
        # Store per-atom charge adjustments from cap redistribution {resname: adjustment_per_atom}
        self._cap_charge_adjustments: Dict[str, float] = {}

    def type_all_atoms(self) -> Dict:
        """
        Run complete 3-tier atom typing workflow.

        Returns:
            Dictionary with 'atom_types' and 'charges' mappings
        """
        self.console.print("\n[bold cyan]Assigning atom types to structure...[/bold cyan]\n")

        # Step 1: Classify terminals
        self._classify_terminals()

        # Step 2: Identify metals and ligands from RedoxSite
        self._identify_metals_and_ligands()

        # Step 3: Three-tier atom typing
        self._three_tier_typing()

        self.console.print(f"\n[green]✅ Assigned types for {len(self.atom_types)} atoms[/green]")

        return {
            'atom_types': self.atom_types,
            'charges': self.charges,
            'used_frcmod_files': self.used_frcmod_files
        }

    def _classify_terminals(self):
        """Step 1: Classify terminal residues using MCPB TerminalClassifier"""
        self.console.print("[bold]Classifying terminal residues...[/bold]\n")

        # TerminalClassifier needs a RedoxSite, but we want to classify the whole structure
        # Use the existing RedoxSite - terminals will only be classified for residues in it
        # For residues not in RedoxSite, we'll default to INTERNAL
        terminal_classifier = TerminalClassifier(
            site=self.redox_site,
            pdb_file=self.pdb_file,
            logger=self.logger
        )

        # Run classification with user confirmation
        # (TerminalClassifier handles display via Rich panels/tables)
        self.terminal_map = terminal_classifier.classify_all_residues(require_confirmation=True)

    def _identify_metals_and_ligands(self):
        """Step 2: Identify metal centers from RedoxSite using residue identity"""
        self.metal_residues.clear()
        self.ligand_atoms.clear()

        # Get metal centers from RedoxSite - identify by chain/resid/atom_name
        for center in self.redox_site.centers:
            if center.center_type.value == 'metal_ion':
                # Find the atom in RedoxSite to get its identity
                for atom in self.redox_site.atoms:
                    if atom.coords == center.coords:
                        # Store (chain, resid, insertion_code, atom_name) - normalize strings
                        metal_key = (
                            atom.chain.strip() if atom.chain else '',
                            atom.resid,
                            atom.insertion_code.strip() if atom.insertion_code else '',
                            atom.atom_name
                        )
                        self.metal_residues.append(metal_key)
                        self.logger.debug(f"Found metal: {atom.chain}:{atom.resid} {atom.atom_name} ({atom.resname})")
                        break

        self.logger.debug(f"Found {len(self.metal_residues)} metal atoms")

    def _three_tier_typing(self):
        """
        Step 3: Three-tier atom typing based on RedoxSite membership

        Tier 1: Standard protein (atoms OUTSIDE RedoxSite) - use selected force fields
        Tier 2: RedoxSite atoms (INSIDE RedoxSite) - try FF, then custom params, then manual
        Tier 3: Manual fallback (anything that failed in Tier 1 or Tier 2)
        """
        tier1_atoms = []  # Outside RedoxSite (standard protein)
        tier2_atoms = []  # Inside RedoxSite (special handling)
        tier3_atoms = []  # Manual fallback

        # Build RedoxSite coordinate set for fast lookup
        redox_coords = set()
        for atom in self.redox_site.atoms:
            redox_coords.add(atom.coords)

        self.console.print(f"[grey50]RedoxSite contains {len(redox_coords)} atom coordinates[/grey50]")

        # First pass: classify all atoms into tiers based on RedoxSite membership
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    chain_id = chain.get_id().strip()  # Normalize: strip spaces
                    resid = residue.get_id()[1]
                    insertion = residue.get_id()[2].strip()  # Normalize: strip spaces

                    for atom in residue:
                        coords = tuple(atom.get_coord())
                        atom_name = atom.get_name()

                        # Check if this atom is in RedoxSite (by coordinates)
                        if self._is_in_redoxsite(coords, redox_coords):
                            tier2_atoms.append((atom, residue, chain))
                        else:
                            tier1_atoms.append((atom, residue, chain))

        # Track totals across all tiers
        total_atoms_typed = len(self.atom_types)  # Atoms typed before tier processing

        self.console.print(f"\n[bold]Tier 1:[/bold] Processing {len(tier1_atoms)} atoms outside RedoxSite (standard protein)...")
        tier1_start = len(self.atom_types)
        if tier1_atoms:
            tier1_failed = self._process_tier1_standard_protein(tier1_atoms)
            tier3_atoms.extend(tier1_failed)
        tier1_typed = len(self.atom_types) - tier1_start

        # Print HIS protonation summary right after Tier 1 (where HIS residues live)
        if hasattr(self, "_his_detections") and self._his_detections:
            from collections import Counter
            state_counts = Counter(d[2] for d in self._his_detections)
            self.console.print(f"\n  [bold]HIS protonation auto-detection summary:[/bold]")
            for chain_id, resid, state in self._his_detections:
                self.console.print(f"    HIS {chain_id}:{resid} → {state}", highlight=False)
            summary_parts = [f"{count}×{state}" for state, count in sorted(state_counts.items())]
            self.console.print(f"    Total: {', '.join(summary_parts)}")

        self.console.print(f"\n[bold]Tier 2:[/bold] Processing {len(tier2_atoms)} atoms inside RedoxSite (active region)...")
        tier2_start = len(self.atom_types)
        if tier2_atoms:
            tier2_failed = self._process_tier2_redoxsite(tier2_atoms)
            tier3_atoms.extend(tier2_failed)
        tier2_typed = len(self.atom_types) - tier2_start

        tier3_typed = 0
        if tier3_atoms:
            self.console.print(f"\n[bold]Tier 3:[/bold] Manual entry required for {len(tier3_atoms)} atoms...")
            tier3_start = len(self.atom_types)
            self._process_tier3_manual_fallback(tier3_atoms)
            tier3_typed = len(self.atom_types) - tier3_start

        # Print grand total summary
        total_atoms_typed = len(self.atom_types) - total_atoms_typed
        self.console.print(f"\n[bold cyan]{'=' * 60}[/bold cyan]")
        self.console.print(f"[bold cyan]Atom Typing Summary[/bold cyan]")
        self.console.print(f"[bold cyan]{'=' * 60}[/bold cyan]")
        self.console.print(f"  Tier 1 (Outside RedoxSite): {tier1_typed} atoms typed")
        self.console.print(f"  Tier 2 (Inside RedoxSite):  {tier2_typed} atoms typed")
        if tier3_typed > 0:
            self.console.print(f"  Tier 3 (Manual):            {tier3_typed} atoms typed")
        self.console.print(f"[bold cyan]{'─' * 60}[/bold cyan]")
        self.console.print(f"[bold cyan]  Total:                      {total_atoms_typed} atoms typed[/bold cyan]")
        self.console.print(f"[bold cyan]{'=' * 60}[/bold cyan]\n")

    def _is_in_redoxsite(self, coords: Tuple[float, float, float],
                        redox_coords: set,
                        tolerance: float = 0.001) -> bool:
        """
        Check if coordinates belong to RedoxSite (with floating point tolerance)

        Args:
            coords: Coordinates to check
            redox_coords: Set of RedoxSite atom coordinates
            tolerance: Distance tolerance for coordinate matching

        Returns:
            True if coords match any RedoxSite coordinate within tolerance
        """
        for redox_coord in redox_coords:
            distance = ((coords[0] - redox_coord[0])**2 +
                       (coords[1] - redox_coord[1])**2 +
                       (coords[2] - redox_coord[2])**2)**0.5
            if distance <= tolerance:
                return True
        return False

    def _process_tier1_standard_protein(self, tier1_atoms: List):
        """
        Tier 1: Process standard protein atoms outside RedoxSite

        Fallback chain:
        1. Try user-selected force field
        2. If residue not found, try custom forcefield directory
        3. If still not found, return atoms for Tier 3 manual entry

        This matches Tier 2 logic - custom forcefields work for both inside
        and outside RedoxSite (e.g., transformed terminal residues).
        """
        typed_count = 0
        custom_typed_count = 0
        failed_atoms = []

        # Group atoms by residue for efficient processing
        residue_atoms = {}
        for atom, residue, chain in tier1_atoms:
            res_key = (
                chain.get_id().strip(),
                residue.get_id()[1],
                residue.get_id()[2].strip(),
                residue.get_resname().strip()
            )
            if res_key not in residue_atoms:
                residue_atoms[res_key] = []
            residue_atoms[res_key].append((atom, residue, chain))

        # Process each residue
        for res_key, atoms in residue_atoms.items():
            chain_id, resid, insertion, resname = res_key

            # Try standard force field first
            if self._try_standard_forcefield(atoms, resname):
                typed_count += len(atoms)
            else:
                # Not in standard FF - try custom forcefield params
                custom_result = self._try_custom_forcefield(atoms, resname)

                if custom_result == "typed":
                    custom_typed_count += len(atoms)
                elif custom_result == "failed":
                    # Add to Tier 3 fallback
                    failed_atoms.extend(atoms)

        self.console.print(
            f"  ✓ Typed {typed_count} atoms from standard force field"
        )
        if custom_typed_count > 0:
            self.console.print(
                f"  ✓ Typed {custom_typed_count} atoms from custom forcefield parameters"
            )
        if failed_atoms:
            self.console.print(
                f"  [yellow]⚠ {len(failed_atoms)} atoms need manual entry (Tier 3)[/yellow]"
            )

        return failed_atoms

    def _process_tier2_redoxsite(self, tier2_atoms: List) -> List:
        """
        Tier 2: Process RedoxSite atoms with intelligent fallback chain

        1. Try user-selected force field
        2. If residue not found, check custom forcefield directory (via transformers)
        3. If custom forcefield found, use it with redox/spin state selection
        4. If still not found, return atoms for Tier 3 manual entry

        Returns:
            List of atoms that need manual entry (Tier 3 fallback)
        """
        typed_count = 0
        custom_typed_count = 0
        failed_atoms = []

        # Group atoms by residue for efficient processing
        residue_atoms = {}
        for atom, residue, chain in tier2_atoms:
            res_key = (
                chain.get_id().strip(),
                residue.get_id()[1],
                residue.get_id()[2].strip(),
                residue.get_resname().strip()
            )
            if res_key not in residue_atoms:
                residue_atoms[res_key] = []
            residue_atoms[res_key].append((atom, residue, chain))

        # Process each residue
        for res_key, atoms in residue_atoms.items():
            chain_id, resid, insertion, resname = res_key

            # Try standard force field first
            if self._try_standard_forcefield(atoms, resname):
                typed_count += len(atoms)
            else:
                # Not in standard FF - try custom forcefield params
                custom_result = self._try_custom_forcefield(atoms, resname)

                if custom_result == "typed":
                    custom_typed_count += len(atoms)
                elif custom_result == "failed":
                    # Add to Tier 3 fallback
                    failed_atoms.extend(atoms)

        # Check for cap atoms in non-standard RedoxSite residues
        # (e.g., ORE parameterized with NME cap — cap atoms carry charge
        # that needs redistributing since they're absent from the PDB)
        checked_resnames = set()
        for res_key, atoms in residue_atoms.items():
            resname = res_key[3]
            if resname in checked_resnames:
                continue
            checked_resnames.add(resname)

            if resname not in self._cap_checked_resnames:
                self._cap_checked_resnames.add(resname)
                self._check_for_cap_atoms(resname, atoms)

        # Apply cap charge adjustments to already-assigned charges
        if self._cap_charge_adjustments:
            for res_key, atoms in residue_atoms.items():
                resname = res_key[3]
                charge_adj = self._cap_charge_adjustments.get(resname, 0.0)
                if charge_adj != 0.0:
                    for atom, residue, chain in atoms:
                        coords = tuple(atom.get_coord())
                        if coords in self.charges:
                            self.charges[coords] += charge_adj

        self.console.print(
            f"  ✓ Typed {typed_count} atoms from standard force field"
        )
        if custom_typed_count > 0:
            self.console.print(
                f"  ✓ Typed {custom_typed_count} atoms from custom forcefield parameters"
            )
        if failed_atoms:
            self.console.print(
                f"  [yellow]⚠ {len(failed_atoms)} atoms need manual entry (Tier 3)[/yellow]"
            )

        return failed_atoms

    # Common AMBER residue name aliases: maps non-standard names to valid FF names
    # The user will be prompted to choose for ambiguous cases (e.g. HIS)
    RESIDUE_ALIASES = {
        # Histidine protonation states
        "HIS": ["HIE", "HID", "HIP"],
        # Cysteine variants
        "CYX": ["CYS"],
        # DNA/RNA variants
        "URI": ["U"],
    }

    # N-terminal atom name aliases: PDB names that map to AMBER lib names
    NTERM_ATOM_ALIASES = {
        "H": "H1",  # PDB uses H for first N-H; AMBER uses H1
    }

    def _detect_his_protonation(self, atoms: List) -> str:
        """
        Detect histidine protonation state from PDB hydrogen atoms.

        Args:
            atoms: List of (atom, residue, chain) tuples for one HIS residue

        Returns:
            "HIP" if both HD1 and HE2 present (doubly protonated)
            "HID" if only HD1 present (delta-protonated)
            "HIE" if only HE2 present (epsilon-protonated)
            "" if neither present (ambiguous, needs user input)
        """
        atom_names = {atom.get_name() for atom, _, _ in atoms}
        has_hd1 = "HD1" in atom_names
        has_he2 = "HE2" in atom_names

        if has_hd1 and has_he2:
            return "HIP"
        elif has_hd1:
            return "HID"
        elif has_he2:
            return "HIE"
        return ""

    def _resolve_resname(self, resname: str, atoms: List = None) -> str:
        """
        Resolve ambiguous residue names to valid force field names.

        For HIS residues, detects protonation state from PDB hydrogen atoms
        (HD1/HE2 presence) and resolves per-residue. Only prompts the user
        for genuinely ambiguous cases (no diagnostic hydrogens present).

        Args:
            resname: Residue name from PDB
            atoms: List of (atom, residue, chain) tuples for this residue instance

        Returns the resolved resname, or original if no alias needed.
        """
        if resname not in self.RESIDUE_ALIASES:
            return resname

        aliases = self.RESIDUE_ALIASES[resname]

        # Filter to those actually in the force field
        available = [a for a in aliases if self.ff_collector.has_residue(a)]

        if not available:
            return resname

        if len(available) == 1:
            # Only one option — use it (cache globally)
            cache_key = f"_alias_{resname}"
            if not hasattr(self, cache_key):
                self.console.print(
                    f"[grey50]Mapping {resname} → {available[0]} (only matching FF entry)[/grey50]"
                )
                setattr(self, cache_key, available[0])
            return available[0]

        # For HIS: try to detect protonation from PDB atoms
        if resname == "HIS" and atoms:
            detected = self._detect_his_protonation(atoms)
            if detected and detected in available:
                # Get residue identifier for reporting
                _, residue, chain = atoms[0]
                chain_id = chain.get_id().strip()
                resid = residue.get_id()[1]

                # Track all HIS detections for summary
                if not hasattr(self, "_his_detections"):
                    self._his_detections = []
                self._his_detections.append((chain_id, resid, detected))

                # Report auto-detection (only first time for each unique state)
                report_key = f"_his_reported_{detected}"
                if not hasattr(self, report_key):
                    setattr(self, report_key, True)
                    descriptions = {
                        "HIE": "epsilon-protonated (Nε-H, HE2 present)",
                        "HID": "delta-protonated (Nδ-H, HD1 present)",
                        "HIP": "doubly protonated (HD1 + HE2 present)",
                    }
                    self.console.print(
                        f"  Auto-detected {resname} → [bold]{detected}[/bold]: "
                        f"{descriptions.get(detected, '')}"
                    )
                return detected

        # Fallback: check global cache (from prior user prompt)
        cache_key = f"_alias_{resname}"
        if hasattr(self, cache_key):
            return getattr(self, cache_key)

        # No auto-detection possible — prompt user
        self.console.print(
            f"\n[yellow]Residue '{resname}' is ambiguous in AMBER force fields.[/yellow]"
        )
        self.console.print(
            "[grey50]Could not auto-detect protonation state from PDB hydrogens.[/grey50]"
        )
        self.console.print("Available options:")
        descriptions = {
            "HIE": "epsilon-protonated (Nε-H, most common)",
            "HID": "delta-protonated (Nδ-H)",
            "HIP": "doubly protonated (both N-H, +1 charge)",
        }
        for i, alias in enumerate(available, 1):
            desc = descriptions.get(alias, "")
            desc_str = f" - {desc}" if desc else ""
            self.console.print(f"  {i}. {alias}{desc_str}")

        prot_options_map = {
            str(i): f"{alias} - {descriptions.get(alias, '')}".rstrip(" -")
            for i, alias in enumerate(available, 1)
        }
        choice = prompt_with_context(
            processor=self.processor,
            prompt=f"Select form for remaining ambiguous '{resname}' residues",
            choices=[str(i) for i in range(1, len(available) + 1)],
            default="1",
            module="ONIOM Atom Typer",
            description=f"Select protonation state for {resname}",
            options_map=prot_options_map,
        )
        resolved = available[int(choice) - 1]

        # Cache for remaining ambiguous residues
        setattr(self, cache_key, resolved)
        self.console.print(
            f"[green]✓ Remaining ambiguous '{resname}' residues will be typed as "
            f"'{resolved}'[/green]"
        )
        return resolved

    def _try_standard_forcefield(self, atoms: List, resname: str) -> bool:
        """
        Try to type atoms using standard force field.

        Handles residue aliasing (HIS→HIE/HID/HIP) and N-terminal
        atom name fallbacks (H→H1).

        Returns:
            True if all atoms were typed successfully
        """
        # Resolve ambiguous residue names
        resolved_resname = self._resolve_resname(resname, atoms)

        # Check if residue exists in force field
        if not self.ff_collector.has_residue(resolved_resname):
            return False

        # Track parsers for connectivity building (load standard amino acid .lib files)
        self._ensure_standard_lib_parser_loaded(resolved_resname)

        # Type all atoms in this residue
        success = True
        for atom, residue, chain in atoms:
            coords = tuple(atom.get_coord())
            atom_name = atom.get_name()

            # Get residue key for terminal classification
            res_key = (
                chain.get_id().strip(),
                residue.get_id()[1],
                residue.get_id()[2].strip()
            )
            terminal_type = self.terminal_map.get(res_key, TerminalType.INTERNAL)

            # Build library lookup key
            if terminal_type == TerminalType.NTERM:
                lookup_resname = f"N{resolved_resname}"
            elif terminal_type == TerminalType.CTERM:
                lookup_resname = f"C{resolved_resname}"
            else:
                lookup_resname = resolved_resname

            # Track terminal variant parsers too
            if lookup_resname != resolved_resname:
                self._ensure_standard_lib_parser_loaded(lookup_resname)

            # Try to get from force field library
            atom_def = self.ff_collector.get_atom_type(lookup_resname, atom_name)

            # N-terminal atom name fallback: PDB "H" → AMBER "H1"
            if not atom_def and terminal_type == TerminalType.NTERM:
                alias = self.NTERM_ATOM_ALIASES.get(atom_name)
                if alias:
                    atom_def = self.ff_collector.get_atom_type(lookup_resname, alias)

            if atom_def:
                self.atom_types[coords] = atom_def.atom_type
                self.charges[coords] = atom_def.atom_charge
            else:
                success = False
                self.logger.warning(
                    f"Atom {lookup_resname}-{atom_name} not found in force field"
                )

        return success

    def _ensure_standard_lib_parser_loaded(self, resname: str):
        """
        Load and cache AmberLibParser for standard amino acid .lib files.
        This ensures we have connectivity information for standard residues.

        Args:
            resname: Residue name (e.g., "ALA", "NALA", "CALA")
        """
        from proprep.oniom_prep.amber_lib_parser import AmberLibParser
        import os

        # Skip if already loaded
        if resname in self.lib_parsers_used:
            return

        # Find the standard amino acid .lib file from ff_collector
        # The ff_collector has paths to leaprc files which reference .lib files
        try:
            # Try to find all.lib or amino*.lib in the AMBER dat directory
            amberhome = os.getenv('AMBERHOME', '/opt/homebrew/anaconda3/envs/proprep_env')
            possible_lib_files = [
                f"{amberhome}/dat/leap/lib/amino19.lib",
                f"{amberhome}/dat/leap/lib/amino12.lib",
                f"{amberhome}/dat/leap/lib/all_amino19.lib",
                f"{amberhome}/dat/leap/lib/all_amino12.lib",
                f"{amberhome}/dat/leap/lib/aminoct12.lib",
                f"{amberhome}/dat/leap/lib/aminont12.lib",
            ]

            # Try each possible lib file
            for lib_file in possible_lib_files:
                if os.path.exists(lib_file):
                    try:
                        parser = AmberLibParser(lib_file)
                        parser.parse()

                        # Check if this lib file contains our residue
                        if resname in parser.residues:
                            self.lib_parsers_used[resname] = parser
                            self.logger.debug(f"Loaded standard lib parser for {resname} from {lib_file}")
                            return
                    except Exception as e:
                        self.logger.debug(f"Could not parse {lib_file}: {e}")
                        continue

            self.logger.debug(f"Could not find .lib file containing standard residue '{resname}'")

        except Exception as e:
            self.logger.debug(f"Error loading standard lib parser for {resname}: {e}")

    def _find_lib_files_for_residue(self, resname: str) -> List[Dict]:
        """
        Scan forcefield_params directory for .lib files containing resname

        This works with TRANSFORMED structures (Option 1 approach).
        Scans all .lib files to find which ones contain the target residue name,
        then extracts metadata from directory structure and metadata.json.

        Args:
            resname: Residue name to search for (e.g., "HCR", "DHR", "PHR")

        Returns:
            List of dicts with complete metadata:
            [{
                "lib_path": "/path/to/file.lib",
                "frcmod_path": "/path/to/file.frcmod",
                "cofactor_family": "heme",
                "cofactor_type": "bis_his_c_type",
                "redox_state": "reduced",
                "spin_state": "low_spin",
                "formal_charge": 0,
                "spin_multiplicity": 1,
                "description": "...",
                "set_name": "Guberman_Reduced",
                "is_default": True,
                "version": "1.0",
                "reference": "..."
            }]
        """
        from proprep.forcefield_params import loader as ff_loader
        from proprep.oniom_prep.amber_lib_parser import AmberLibParser
        from pathlib import Path

        results = []
        base_path = ff_loader.get_forcefield_base_path()

        if not base_path.exists():
            self.logger.warning(f"Forcefield params directory not found: {base_path}")
            return results

        # Recursively find all .lib files
        for lib_file in base_path.rglob("*.lib"):
            try:
                # Parse .lib file to check if it contains the residue
                parser = AmberLibParser(str(lib_file))
                parser.parse()

                if resname not in parser.residues:
                    continue  # This .lib doesn't have our residue

                # Extract metadata from directory path
                # Path format: .../cofactor_family/specific_type/redox_state/spin_state/file.lib
                parts = lib_file.relative_to(base_path).parts

                if len(parts) < 5:
                    self.logger.warning(f"Unexpected path structure: {lib_file}")
                    continue

                cofactor_family = parts[0]
                cofactor_type = parts[1]
                redox_state = parts[2]
                spin_state = parts[3]
                lib_filename = parts[4]

                # Find corresponding .frcmod file
                frcmod_file = lib_file.with_suffix('.frcmod')
                frcmod_path = str(frcmod_file) if frcmod_file.exists() else None

                # Load metadata.json for this cofactor type
                cofactor_path = f"{cofactor_family}/{cofactor_type}"
                try:
                    metadata = ff_loader.load_forcefield_metadata(cofactor_path)

                    # Extract state-specific metadata
                    redox_data = metadata.get("redox_states", {}).get(redox_state, {})
                    spin_data = redox_data.get("spin_states", {}).get(spin_state, {})

                    formal_charge = redox_data.get("formal_charge", 0)
                    spin_multiplicity = spin_data.get("spin_multiplicity", 1)
                    description = spin_data.get("description", "")

                    # Check if this lib file is in the forcefield_sets
                    forcefield_sets = spin_data.get("forcefield_sets", {})
                    set_name = None
                    is_default = False
                    version = ""
                    reference = ""

                    # Find which set this lib file belongs to
                    for name, set_info in forcefield_sets.items():
                        if set_info.get("files", {}).get("lib") == lib_filename:
                            set_name = name
                            is_default = set_info.get("is_default", False)
                            version = set_info.get("version", "")
                            reference = set_info.get("reference", "")
                            break

                    if not set_name:
                        # Fallback: use filename as set name
                        set_name = lib_file.stem

                except Exception as e:
                    self.logger.warning(f"Could not load metadata for {cofactor_path}: {e}")
                    # Use defaults
                    formal_charge = 0
                    spin_multiplicity = 1
                    description = ""
                    set_name = lib_file.stem
                    is_default = False
                    version = ""
                    reference = ""

                # Add to results
                results.append({
                    "lib_path": str(lib_file),
                    "frcmod_path": frcmod_path,
                    "cofactor_family": cofactor_family,
                    "cofactor_type": cofactor_type,
                    "redox_state": redox_state,
                    "spin_state": spin_state,
                    "formal_charge": formal_charge,
                    "spin_multiplicity": spin_multiplicity,
                    "description": description,
                    "set_name": set_name,
                    "is_default": is_default,
                    "version": version,
                    "reference": reference
                })

            except Exception as e:
                self.logger.debug(f"Error parsing {lib_file}: {e}")
                continue

        # Sort by: default first, then alphabetically
        results.sort(key=lambda x: (not x["is_default"], x["set_name"]))

        self.logger.debug(f"Found {len(results)} forcefield sets containing residue '{resname}'")
        return results

    def _select_forcefield_set(self, residue_name: str, forcefield_sets: List[Dict]) -> Optional[Dict]:
        """
        Present user with forcefield set choices and let them select

        Args:
            residue_name: Name of residue being typed
            forcefield_sets: List of available forcefield sets

        Returns:
            Selected forcefield set dict, or None if cancelled
        """
        from rich.table import Table
        from rich.prompt import Prompt

        if not forcefield_sets:
            return None

        if len(forcefield_sets) == 1:
            # Auto-select if only one option
            ff_set = forcefield_sets[0]
            self.console.print(
                f"[grey50]Auto-selected forcefield for {residue_name}: "
                f"{ff_set['cofactor_type']} ({ff_set['redox_state']}/{ff_set['spin_state']})[/grey50]"
            )
            return ff_set

        # Display options in Rich table
        self.console.print(f"\n[bold]Found '{residue_name}' in {len(forcefield_sets)} forcefield sets:[/bold]\n")

        table = Table(title=f"Forcefield Options for {residue_name}")
        table.add_column("#", style="cyan", width=3)
        table.add_column("Cofactor Type", style="green")
        table.add_column("Redox State", style="yellow")
        table.add_column("Spin State", style="magenta")
        table.add_column("Charge", style="blue", justify="right")
        table.add_column("Default", style="grey50", width=7)

        for i, ff_set in enumerate(forcefield_sets, 1):
            # Format cofactor type nicely
            cofactor_display = f"{ff_set['cofactor_family']} ({ff_set['cofactor_type']})"

            # Format redox state
            redox_display = f"{ff_set['redox_state']}"

            # Format spin state with multiplicity
            spin_display = f"{ff_set['spin_state']} (2S+1={ff_set['spin_multiplicity']})"

            # Default marker
            default_marker = "✓ Yes" if ff_set['is_default'] else ""

            table.add_row(
                str(i),
                cofactor_display,
                redox_display,
                spin_display,
                str(ff_set['formal_charge']),
                default_marker
            )

        self.console.print(table)

        # Show additional details for each option
        self.console.print("\n[bold]Details:[/bold]")
        for i, ff_set in enumerate(forcefield_sets, 1):
            self.console.print(f"  {i}. {ff_set['set_name']}")
            if ff_set['description']:
                self.console.print(f"     {ff_set['description']}")
            if ff_set['reference']:
                self.console.print(f"     [grey50]Ref: {ff_set['reference']}[/grey50]")

        # Get user choice
        default_idx = next((i for i, fs in enumerate(forcefield_sets, 1) if fs['is_default']), 1)

        # Create options_map for rich context
        options_map = {
            str(i): f"{fs['set_name']} ({fs['redox_state']}/{fs['spin_state']})"
            for i, fs in enumerate(forcefield_sets, 1)
        }

        choice = prompt_with_context(
            processor=self.processor,
            prompt="\nSelect forcefield set",
            choices=[str(i) for i in range(1, len(forcefield_sets) + 1)],
            default=str(default_idx),
            module="ONIOM Atom Typer",
            description="Select forcefield set for metal site",
            options_map=options_map
        )

        # Replay by FF-set identity so a changed/reordered set list re-selects
        # the same set (forcefield_sets are dicts, not file paths).
        _fs_key = lambda fs: f"{fs['set_name']}/{fs['redox_state']}/{fs['spin_state']}"
        choice = remap_recorded_index_by_key(self.processor, forcefield_sets, _fs_key, str(choice))
        selected = forcefield_sets[int(choice) - 1]
        annotate_recorded_key(self.processor, _fs_key(selected))

        self.console.print(
            f"\n[green]✓ Selected: {selected['set_name']} "
            f"({selected['redox_state']}/{selected['spin_state']})[/green]\n"
        )

        return selected

    def _try_custom_forcefield(self, atoms: List, resname: str) -> str:
        """
        Try to type atoms using custom forcefield parameters (Option 1 approach)

        Strategy:
        1. Scan forcefield_params directory for .lib files containing this residue
        2. Present user with available forcefield sets (with metadata)
        3. Parse selected .lib file and extract atom types/charges
        4. Apply to atoms

        This works with TRANSFORMED structures - expects residue names like
        HCR, DHR, PHR (not HEM, HIS).

        Returns:
            "typed" if successful, "failed" if not found or error
        """
        from proprep.oniom_prep.amber_lib_parser import AmberLibParser

        try:
            # Step 1: Find all .lib files containing this residue
            self.console.print(f"\n[grey50]Searching for custom forcefield parameters for {resname}...[/grey50]")
            forcefield_sets = self._find_lib_files_for_residue(resname)

            if not forcefield_sets:
                self.logger.debug(
                    f"No custom forcefield parameters found for residue '{resname}'"
                )
                return "failed"

            # Step 2: Let user select forcefield set
            selected_set = self._select_forcefield_set(resname, forcefield_sets)

            if not selected_set:
                return "failed"

            # Step 3: Parse selected .lib file
            lib_file = selected_set["lib_path"]
            parser = AmberLibParser(lib_file)
            parser.parse()

            # Track lib parser for connectivity building
            if resname not in self.lib_parsers_used:
                self.lib_parsers_used[resname] = parser
                self.logger.debug(f"Tracked lib parser for {resname}: {lib_file}")

            # Track the corresponding .frcmod file for later use
            frcmod_file = selected_set.get("frcmod_path")
            if frcmod_file and frcmod_file not in self.used_frcmod_files:
                self.used_frcmod_files.append(frcmod_file)
                self.logger.debug(f"Tracked frcmod file for ONIOM: {frcmod_file}")

            # Step 4: Apply atom types and charges to atoms
            typed_count = 0
            for atom, residue, chain in atoms:
                coords = tuple(atom.get_coord())
                atom_name = atom.get_name()

                # Get atom type and charge from .lib file
                result = parser.get_atom_type_and_charge(resname, atom_name)

                if result:
                    atom_type, charge = result
                    self.atom_types[coords] = atom_type
                    self.charges[coords] = charge
                    typed_count += 1
                else:
                    self.logger.warning(
                        f"Atom {atom_name} not found in {resname} from {lib_file}"
                    )

            if typed_count == len(atoms):
                # Successfully typed all atoms - suppress individual message,
                # will show in tier summary
                return "typed"
            else:
                self.console.print(
                    f"  [yellow]⚠ Only typed {typed_count}/{len(atoms)} atoms for {resname}[/yellow]"
                )
                return "failed"

        except Exception as e:
            self.logger.error(f"Error in custom forcefield lookup: {e}")
            import traceback
            traceback.print_exc()
            return "failed"

    def _process_tier3_manual_fallback(self, tier3_atoms: List):
        """Tier 3: Manual fallback for atoms that couldn't be typed automatically"""
        if not tier3_atoms:
            return

        # Group atoms by residue
        residue_atoms: Dict[Tuple, List] = {}
        for atom, residue, chain in tier3_atoms:
            res_key = (chain.get_id(), residue.get_id()[1], residue.get_resname())
            if res_key not in residue_atoms:
                residue_atoms[res_key] = []
            residue_atoms[res_key].append((atom, residue, chain))

        # Display summary header
        self.console.print()
        self.console.print("[bold]" + "=" * 60 + "[/bold]")
        self.console.print("[bold]Untyped Residues[/bold]")
        self.console.print("[bold]" + "=" * 60 + "[/bold]\n")

        self.console.print("The following residues could not be typed automatically:\n")

        # Display summary table of untyped residues
        summary_table = Table(show_header=True, header_style="bold")
        summary_table.add_column("Residue", style="cyan")
        summary_table.add_column("Chain", style="green")
        summary_table.add_column("ResID", style="yellow")
        summary_table.add_column("Atom Count", style="magenta")

        unique_resnames = set()
        for res_key, atoms in residue_atoms.items():
            chain_id, resid, resname = res_key
            unique_resnames.add(resname)
            summary_table.add_row(
                resname,
                chain_id,
                str(resid),
                f"{len(atoms)} atoms"
            )

        self.console.print(summary_table)
        self.console.print()

        total_residues = len(residue_atoms)
        total_atoms = len(tier3_atoms)
        self.console.print(f"[bold]Total: {total_residues} residues, {total_atoms} atoms[/bold]\n")

        # Offer to load custom .lib/.frcmod files before manual entry
        self.console.print(
            "You can provide AMBER .lib/.frcmod files to type these residues,\n"
            "or fall back to manual per-atom entry.\n"
        )

        # Try loading user-provided files for each unique untyped residue
        still_untyped = {}
        for res_key, atoms in residue_atoms.items():
            chain_id, resid, resname = res_key
            # Skip if we already typed this resname via a user file
            if resname in self.lib_parsers_used:
                # Apply cached parser
                if self._apply_lib_parser_to_atoms(resname, atoms):
                    continue

            still_untyped[res_key] = atoms

        # For remaining untyped residues, offer to load files
        untyped_resnames = set(rk[2] for rk in still_untyped)
        for resname in sorted(untyped_resnames):
            if self._offer_load_custom_lib(resname):
                # Successfully loaded — apply to all residues with this name
                newly_typed = []
                for res_key, atoms in list(still_untyped.items()):
                    if res_key[2] == resname:
                        if self._apply_lib_parser_to_atoms(resname, atoms):
                            newly_typed.append(res_key)
                for rk in newly_typed:
                    del still_untyped[rk]

        # Manual entry for anything still untyped
        if still_untyped:
            self.console.print(
                f"\n[bold yellow]Manual entry required for "
                f"{sum(len(a) for a in still_untyped.values())} remaining atoms[/bold yellow]\n"
            )

            for res_key, atoms in still_untyped.items():
                chain_id, resid, resname = res_key

                for atom, residue, chain in atoms:
                    coords = tuple(atom.get_coord())
                    element = atom.element.strip() if hasattr(atom, 'element') and atom.element else atom.get_name()[0:2]
                    atom_name = atom.get_name()

                    # Query available atom types for this element
                    self._display_element_atom_types(element)

                    # Prompt for atom type
                    self.console.print(f"\n[yellow]{resname} {atom_name} at {coords}[/yellow]")
                    atom_type = prompt_with_context(
                        processor=self.processor,
                        prompt="  GAFF atom type",
                        default=element.lower(),
                        module="ONIOM Atom Typer",
                        description=f"Enter GAFF atom type for {resname} {atom_name}"
                    )
                    charge_str = prompt_with_context(
                        processor=self.processor,
                        prompt="  Partial charge",
                        default="0.0",
                        module="ONIOM Atom Typer",
                        description=f"Enter partial charge for {resname} {atom_name}"
                    )
                    charge = float(charge_str)

                    self.atom_types[coords] = atom_type
                    self.charges[coords] = charge

                    self.console.print(f"  [green]✓ {resname} {atom_name} parameterized[/green]\n")

    def _offer_load_custom_lib(self, resname: str) -> bool:
        """
        Offer to load a user-provided .lib file (and optional .frcmod) for a residue.

        Also scans the working directory and workspace for matching files.

        Args:
            resname: Residue name to find parameters for

        Returns:
            True if a .lib file was loaded successfully
        """
        import os
        import glob
        from proprep.oniom_prep.amber_lib_parser import AmberLibParser, AmberPrepParser

        self.console.print(f"\n[bold cyan]Residue '{resname}' not found in any force field[/bold cyan]")

        # Auto-scan working directory for .lib and .prep files containing this residue
        cwd = os.getcwd()
        discovered_libs = []

        for pattern in ["*.lib", "*.off", "*.prep"]:
            # Glob CWD itself and subdirectories separately
            # (os.path.join(cwd, "**", pattern) only matches subdirectories)
            for lib_path in glob.glob(os.path.join(cwd, pattern)):
                try:
                    if lib_path.endswith('.prep'):
                        parser = AmberPrepParser(lib_path)
                    else:
                        parser = AmberLibParser(lib_path)
                    parser.parse()
                    if resname in parser.residues:
                        discovered_libs.append(lib_path)
                except Exception:
                    continue
            for lib_path in glob.glob(os.path.join(cwd, "**", pattern), recursive=True):
                if lib_path not in discovered_libs:
                    try:
                        if lib_path.endswith('.prep'):
                            parser = AmberPrepParser(lib_path)
                        else:
                            parser = AmberLibParser(lib_path)
                        parser.parse()
                        if resname in parser.residues:
                            discovered_libs.append(lib_path)
                    except Exception:
                        continue

        if discovered_libs:
            self.console.print(f"\n[green]Found .lib files containing '{resname}' in working directory:[/green]")
            for i, path in enumerate(discovered_libs, 1):
                self.console.print(f"  {i}. {os.path.relpath(path, cwd)}")

            if confirm_with_context(
                self.processor, f"\nLoad discovered file for '{resname}'?", default=True,
                module="ONIOM Atom Typer",
                description=f"Load discovered .lib file for {resname}",
            ):
                if len(discovered_libs) == 1:
                    lib_path = discovered_libs[0]
                else:
                    lib_options_map = {
                        str(i): os.path.relpath(path, cwd)
                        for i, path in enumerate(discovered_libs, 1)
                    }
                    choice = int_prompt_with_context(
                        self.processor, "Select file number", default=1,
                        module="ONIOM Atom Typer",
                        description=f"Select .lib file for {resname}",
                        options_map=lib_options_map,
                    )
                    if 1 <= choice <= len(discovered_libs):
                        lib_path = discovered_libs[choice - 1]
                    else:
                        return False

                return self._load_lib_and_frcmod(resname, lib_path)

        # Manual path entry
        if confirm_with_context(
            self.processor, f"Provide a .lib file path for '{resname}'?", default=False,
            module="ONIOM Atom Typer",
            description=f"Provide .lib file path for {resname}",
        ):
            lib_path = prompt_with_context(
                processor=self.processor,
                prompt=f"Path to .lib file for {resname}",
                module="ONIOM Atom Typer",
                description=f"Path to AMBER .lib file containing {resname}"
            )
            lib_path = os.path.expanduser(lib_path.strip())

            if os.path.exists(lib_path):
                return self._load_lib_and_frcmod(resname, lib_path)
            else:
                self.console.print(f"[red]File not found: {lib_path}[/red]")

        return False

    def _load_lib_and_frcmod(self, resname: str, lib_path: str) -> bool:
        """
        Load a .lib file and optionally its corresponding .frcmod file.

        Args:
            resname: Residue name
            lib_path: Path to .lib file

        Returns:
            True if loaded successfully
        """
        import os
        from proprep.oniom_prep.amber_lib_parser import AmberLibParser, AmberPrepParser

        try:
            if lib_path.endswith('.prep'):
                parser = AmberPrepParser(lib_path)
            else:
                parser = AmberLibParser(lib_path)
            parser.parse()

            if resname not in parser.residues:
                self.console.print(f"[red]Residue '{resname}' not found in {lib_path}[/red]")
                return False

            self.lib_parsers_used[resname] = parser
            file_type = ".prep" if lib_path.endswith('.prep') else ".lib"
            self.console.print(f"[green]✓ Loaded {file_type} file for {resname}: {lib_path}[/green]")

            # Look for corresponding .frcmod file
            # Look for corresponding .frcmod file with multiple naming patterns
            base_dir = os.path.dirname(lib_path)
            base_name = os.path.splitext(os.path.basename(lib_path))[0]
            frcmod_candidates = [
                os.path.splitext(lib_path)[0] + ".frcmod",       # same_name.frcmod
                os.path.join(base_dir, f"frcmod.{base_name}"),   # frcmod.same_name
            ]
            # Also try without the last extension component (e.g., "ore_nme.amber-dyes.gaff2.prep" → "frcmod.ore_nme.amber-dyes.gaff2")
            frcmod_path = None
            for candidate in frcmod_candidates:
                if os.path.exists(candidate):
                    frcmod_path = candidate
                    break

            if not frcmod_path:
                # Broader scan: look for any frcmod file in the same directory
                # that shares a common prefix with the lib/prep file
                for f in os.listdir(base_dir):
                    if f.startswith("frcmod") and base_name.split('.')[0] in f:
                        frcmod_path = os.path.join(base_dir, f)
                        break

            if frcmod_path and os.path.exists(frcmod_path):
                if frcmod_path not in self.used_frcmod_files:
                    self.used_frcmod_files.append(frcmod_path)
                self.console.print(f"[green]✓ Loaded .frcmod file: {frcmod_path}[/green]")
            else:
                # Ask user if they have a .frcmod file
                if confirm_with_context(
                    self.processor, f"Provide a .frcmod file for '{resname}'?", default=False,
                    module="ONIOM Atom Typer",
                    description=f"Provide .frcmod file for {resname}",
                ):
                    frcmod_input = prompt_with_context(
                        processor=self.processor,
                        prompt=f"Path to .frcmod file for {resname}",
                        module="ONIOM Atom Typer",
                        description=f"Path to AMBER .frcmod file for {resname}"
                    )
                    frcmod_input = os.path.expanduser(frcmod_input.strip())
                    if os.path.exists(frcmod_input):
                        if frcmod_input not in self.used_frcmod_files:
                            self.used_frcmod_files.append(frcmod_input)
                        self.console.print(f"[green]✓ Loaded .frcmod file: {frcmod_input}[/green]")
                    else:
                        self.console.print(f"[yellow]Warning: .frcmod file not found: {frcmod_input}[/yellow]")

            return True

        except Exception as e:
            self.console.print(f"[red]Error loading {lib_path}: {e}[/red]")
            self.logger.error(f"Error loading custom .lib file: {e}", exc_info=True)
            return False

    def _apply_lib_parser_to_atoms(self, resname: str, atoms: List) -> bool:
        """
        Apply a cached lib parser to type all atoms in a residue.

        On first call for a given resname, detects cap atoms (atoms in the
        lib/prep file but not in the PDB) and offers to redistribute their
        charge across the matched atoms.

        Args:
            resname: Residue name
            atoms: List of (atom, residue, chain) tuples

        Returns:
            True if all atoms were typed successfully
        """
        parser = self.lib_parsers_used.get(resname)
        if not parser:
            return False

        # On first encounter of this resname, check for cap atoms
        if resname not in self._cap_checked_resnames:
            self._cap_checked_resnames.add(resname)
            self._check_for_cap_atoms(resname, atoms)

        charge_adj = self._cap_charge_adjustments.get(resname, 0.0)

        all_typed = True
        for atom, residue, chain in atoms:
            coords = tuple(atom.get_coord())
            atom_name = atom.get_name()

            result = parser.get_atom_type_and_charge(resname, atom_name)
            if result:
                atom_type, charge = result
                self.atom_types[coords] = atom_type
                self.charges[coords] = charge + charge_adj
            else:
                all_typed = False
                self.logger.warning(f"Atom {atom_name} not found in .lib for {resname}")

        if all_typed:
            self.console.print(
                f"  [green]✓ Typed {len(atoms)} atoms for {resname} from custom .lib file[/green]"
            )

        return all_typed

    def _check_for_cap_atoms(self, resname: str, atoms: List):
        """
        Detect cap atoms in a lib/prep file that are not present in the PDB.

        If cap atoms are found, warns the user and offers to redistribute
        the cap charge across the matched atoms.

        Works with both lib_parsers_used (Tier 3 custom .lib loading) and
        ff_collector (Step 1 force field loading) as data sources.

        Args:
            resname: Residue name
            atoms: List of (atom, residue, chain) tuples for one instance of this residue
        """
        # Try lib parser first, then fall back to ff_collector
        parser = self.lib_parsers_used.get(resname)
        if parser:
            lib_atoms_raw = parser.get_residue_atoms(resname)
            if not lib_atoms_raw:
                return
            # lib parser returns dicts with "name" and "charge" keys
            ff_atom_names = {a["name"] for a in lib_atoms_raw}
            ff_atom_charges = {a["name"]: a["charge"] for a in lib_atoms_raw}
        else:
            # Fall back to ff_collector (e.g., ORE loaded via leaprc in Step 1)
            ff_atom_defs = self.ff_collector.get_residue_atoms(resname)
            if not ff_atom_defs:
                return
            ff_atom_names = {a.atom_name for a in ff_atom_defs}
            ff_atom_charges = {a.atom_name: a.atom_charge for a in ff_atom_defs}

        pdb_atom_names = {atom.get_name() for atom, _, _ in atoms}

        cap_atoms = ff_atom_names - pdb_atom_names
        if not cap_atoms:
            return

        # Compute charges
        cap_charge = sum(
            ff_atom_charges[name] for name in cap_atoms
        )
        matched_charge = sum(
            ff_atom_charges[name] for name in ff_atom_names if name in pdb_atom_names
        )
        total_charge = sum(ff_atom_charges.values())
        nearest_int = round(total_charge)

        self.console.print()
        self.console.print(
            f"[bold yellow]⚠ Cap atom mismatch for residue '{resname}'[/bold yellow]"
        )
        self.console.print(
            f"  The lib/prep file has [cyan]{len(ff_atom_names)}[/cyan] atoms, "
            f"but the PDB has [cyan]{len(pdb_atom_names)}[/cyan]."
        )
        self.console.print(
            f"  [yellow]{len(cap_atoms)}[/yellow] atoms in lib/prep not found in PDB "
            f"(likely parameterization cap): [grey50]{', '.join(sorted(cap_atoms))}[/grey50]"
        )
        self.console.print(
            f"  Charge of cap atoms:     [magenta]{cap_charge:+.4f}[/magenta]"
        )
        self.console.print(
            f"  Charge of matched atoms: [magenta]{matched_charge:+.4f}[/magenta]  "
            f"(expected: [green]{nearest_int:.1f}[/green])"
        )
        self.console.print()
        self.console.print(
            "[grey50]Note: For rigorous calculations, the residue should be re-parameterized "
            "with cap atom charges constrained to zero during charge fitting (e.g., RESP). "
            "Redistributing cap charge is a pragmatic approximation only.[/grey50]"
        )

        if confirm_with_context(
            self.processor,
            f"Redistribute cap charge ({cap_charge:+.4f}) across {len(pdb_atom_names)} matched atoms?",
            default=True,
            module="ONIOM Atom Typer",
            description="Redistribute cap charge across matched atoms",
        ):
            adj = cap_charge / len(pdb_atom_names)
            self._cap_charge_adjustments[resname] = adj
            new_total = matched_charge + cap_charge
            self.console.print(
                f"  [green]✓ Adjusting each atom by {adj:+.6f}  →  "
                f"residue charge = {new_total:.4f}[/green]"
            )
        else:
            self.console.print(
                f"  [yellow]Keeping original charges "
                f"(residue charge = {matched_charge:.4f})[/yellow]"
            )

    def _display_element_atom_types(self, element: str):
        """Display available atom types for an element from force field and GAFF2"""
        # Get atom types for this element
        element_types = self.ff_collector.get_atom_types_for_element(element, include_gaff2=True)

        if not element_types:
            return

        # Display force field table
        if element_types.get('force_field'):
            table = Table(title=f"{element} Atom Types from Selected Force Field", show_header=True)
            table.add_column("Residue", style="cyan")
            table.add_column("Atom Name", style="yellow")
            table.add_column("Atom Type", style="green")
            table.add_column("Charge", style="magenta", justify="right")
            table.add_column("Description", style="grey50")

            for res, name, atype, charge, desc in element_types['force_field'][:10]:  # Limit to 10
                table.add_row(res, name, atype, f"{charge:8.4f}", desc)

            self.console.print(table)
            self.console.print()

        # Display GAFF2 table
        if element_types.get('gaff2'):
            table = Table(title=f"{element} Atom Types from GAFF2 (General Ligands)", show_header=True)
            table.add_column("Residue", style="cyan")
            table.add_column("Atom Name", style="yellow")
            table.add_column("Atom Type", style="green")
            table.add_column("Charge", style="magenta", justify="right")
            table.add_column("Description", style="grey50")

            for res, name, atype, charge, desc in element_types['gaff2']:
                table.add_row(res, name, atype, f"{charge:6.4f}", desc)

            self.console.print(table)
            self.console.print()
