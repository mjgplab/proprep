"""
Structure Completeness Module - Complete Implementation

Detects and repairs missing atoms, missing residues, and alternate locations in PDB structures.
Integrates with ProPrep's modular architecture and workspace system.
"""

import io
import json
import logging
import os
import re
import sys
import urllib.request
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
from Bio.PDB import Atom, Chain, Model, PDBParser, Residue, Structure
from Bio.PDB.PDBIO import PDBIO

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

from proprep.utils.prompts import prompt_with_context, confirm_with_context

# Optional MODELLER import — suppress C-level error messages from
# _modeller.mod_start() when the license key is missing or invalid.
try:
    import ctypes as _ctypes
    _devnull = os.open(os.devnull, os.O_WRONLY)
    _old_stdout = os.dup(1)
    _old_stderr = os.dup(2)
    try:
        os.dup2(_devnull, 1)
        os.dup2(_devnull, 2)
        import modeller
        from modeller.automodel import AutoModel
        HAS_MODELLER = True
    finally:
        _ctypes.CDLL(None).fflush(None)
        os.dup2(_old_stdout, 1)
        os.dup2(_old_stderr, 2)
        os.close(_devnull)
        os.close(_old_stdout)
        os.close(_old_stderr)
except Exception:
    HAS_MODELLER = False
    import sys
    if getattr(sys, 'frozen', False):
        # Running as PyInstaller executable
        print(
            "MODELLER unavailable. Structure repair (missing loops/residues) will be disabled.\n"
            "  MODELLER is free for academic use. Register at:\n"
            "    https://salilab.org/modeller/registration.html\n"
            "  Then save your key (one-time setup):\n"
            "    mkdir -p ~/.proprep && echo 'YOUR_KEY' > ~/.proprep/modeller_key\n"
            "  Or set the environment variable before each session:\n"
            "    export KEY_MODELLER='YOUR_KEY'",
            file=sys.stderr,
        )
    else:
        # Running from source
        print(
            "MODELLER not found. Structure repair (missing loops/residues) will be unavailable.\n"
            "  To install: conda install -c salilab modeller\n"
            "  A free academic license is required: https://salilab.org/modeller/registration.html\n"
            "  After registering, set your license key:\n"
            "    export KEY_MODELLER='your_license_key_here'\n"
            "  Or edit the config.py file shown during conda installation.",
            file=sys.stderr,
        )


# ============================================================================
# DATA STRUCTURES
# ============================================================================

class RepairAction(Enum):
    """Actions for missing residue segments"""
    FILL = "fill"
    PARTIAL_FILL = "partial_fill"
    CAP = "cap"
    SKIP = "skip"
    TER = "ter"          # Leave the gap unfilled but insert a TER record so
                         # tleap treats it as a chain break (no spurious long bond)


@dataclass(frozen=True)
class ResidueIdentity:
    """Immutable residue identifier"""
    chain_id: str
    res_num: int
    res_name: str
    insertion_code: str = " "
    
    def __repr__(self):
        icode = f":{self.insertion_code}" if self.insertion_code.strip() else ""
        return f"{self.res_name}_{self.chain_id}:{self.res_num}{icode}"


@dataclass(frozen=True)
class AtomIdentity:
    """Immutable atom identifier"""
    residue: ResidueIdentity
    atom_name: str
    element: str


@dataclass
class MissingSegment:
    """Represents a contiguous segment of missing residues"""
    chain_id: str
    residues: List[ResidueIdentity]
    start_num: int
    end_num: int
    is_terminal: bool = False
    terminal_type: Optional[str] = None  # 'N' or 'C'
    
    @property
    def length(self) -> int:
        return len(self.residues)
    
    @property
    def is_single(self) -> bool:
        return self.length == 1


@dataclass
class RepairPlan:
    """Plan for repairing a structure"""
    segments_to_fill: List[MissingSegment] = field(default_factory=list)
    segments_to_cap: List[MissingSegment] = field(default_factory=list)
    segments_to_skip: List[MissingSegment] = field(default_factory=list)
    segments_to_ter: List[MissingSegment] = field(default_factory=list)
    mutations: List[Tuple[str, int, str, str]] = field(default_factory=list)

    @property
    def has_fills(self) -> bool:
        return bool(self.segments_to_fill)

    @property
    def has_caps(self) -> bool:
        return bool(self.segments_to_cap)

    @property
    def has_ter(self) -> bool:
        return bool(self.segments_to_ter)

    @property
    def has_mutations(self) -> bool:
        return bool(self.mutations)
    
    @property
    def needs_modeller(self) -> bool:
        return self.has_fills or self.has_mutations


# ============================================================================
# RESIDUE MAPPER - Centralized mapping logic
# ============================================================================

class ResidueMapper:
    """
    Centralized class for all residue number and chain ID mapping.
    
    MULTI-CHAIN SUPPORT:
    - Handles any number of chains (A-Z and beyond)
    - Each chain tracked independently
    - Chain renumbering: MODELLER may rename chains to start from 'A'
    - Residue renumbering: Per-chain, consecutive from 1
    
    Handles the complete mapping chain:
    Original Structure → MODELLER Output → Capped Structure → Final Structure
    """
    
    def __init__(self, console: Console):
        self.console = console
        
        # Core mappings - ALL support multiple chains
        self.chain_mapping: Dict[str, str] = {}  # original_chain → new_chain
        self.residue_mappings: Dict[str, Dict[int, int]] = {}  # new_chain → {original_num → modeller_num}
        self.final_mappings: Dict[str, Dict[int, int]] = {}  # new_chain → {modeller_num → final_num}
    
    def build_modeller_mapping(self, 
                              original_structure: Structure,
                              residues_to_fill: Dict[str, List[ResidueIdentity]],
                              structure_metadata: Dict) -> None:
        """
        Build mapping for MODELLER transformation.
        
        MULTI-CHAIN HANDLING:
        1. Detect all chains in original structure
        2. Predict MODELLER's chain renumbering (A, B, C, ...)
        3. Build per-chain residue mappings independently
        """
        # Detect chain renumbering - HANDLES ALL CHAINS
        original_chains = self._get_chain_order(original_structure)
        self.chain_mapping = self._predict_chain_mapping(original_chains)
        
        self.console.print(f"[cyan]Processing {len(original_chains)} chains: {', '.join(original_chains)}[/cyan]")
        
        # Build residue mappings for EACH chain
        # IMPORTANT: MODELLER numbers residues sequentially across ALL chains, not per-chain
        global_offset = 0
        for original_chain in original_chains:
            new_chain = self.chain_mapping[original_chain]

            # Get all residues that will exist after MODELLER for THIS chain.
            # CRITICAL: preserve PDB file order — MODELLER assigns sequential
            # numbers based on the order residues appear in the file, not by
            # sorted residue number.  Filled residues are standard amino acids
            # that slot into the sequence by position.
            present_ordered = self._get_present_residues_ordered(original_structure, original_chain)
            filling = {r.res_num for r in residues_to_fill.get(original_chain, [])}

            # Insert filled residues at their correct sequence positions.
            # They are always standard residues that precede existing residues
            # with higher numbers, so insert each one just before the first
            # present residue whose number is greater.
            final_residues = list(present_ordered)
            for fill_num in sorted(filling):
                if fill_num not in final_residues:
                    # Find insertion point: before the first element > fill_num
                    insert_idx = 0
                    for i, r in enumerate(final_residues):
                        if r > fill_num:
                            insert_idx = i
                            break
                        insert_idx = i + 1
                    final_residues.insert(insert_idx, fill_num)

            # Determine MODELLER's renumbering behavior
            # MODELLER always renumbers consecutively from 1
            # Key: If final sequence is consecutive 1-N, output matches input (identity mapping)
            #      If final sequence has gaps (e.g., 1-49, 56-100), MODELLER compacts to 1-94 (sequential)

            if final_residues:
                min_res = min(final_residues)
                max_res = max(final_residues)
                num_res = len(final_residues)

                # MODELLER numbers residues globally across all chains:
                # Chain A: 1-N1, Chain B: (N1+1)-(N1+N2), etc.
                if min_res == 1 and max_res == num_res and global_offset == 0:
                    # First chain with consecutive 1-N → Identity mapping
                    mapping = {orig_num: orig_num for orig_num in final_residues}
                    self.console.print(f"[cyan]    Identity mapping: consecutive 1-{num_res}[/cyan]")
                else:
                    # Sequential mapping with global offset
                    # MODELLER renumbers: orig positions → (global_offset+1), (global_offset+2), ...
                    mapping = {orig_num: idx + 1 + global_offset for idx, orig_num in enumerate(final_residues)}
                    start = global_offset + 1
                    end = global_offset + num_res
                    self.console.print(f"[cyan]    Sequential mapping: {min_res}-{max_res} → {start}-{end}[/cyan]")
            else:
                mapping = {}

            self.residue_mappings[new_chain] = mapping

            # Update global offset for next chain
            global_offset += len(final_residues)
        
        self._log_mapping("MODELLER", self.residue_mappings)

    def build_identity_mapping(self, structure: Structure) -> None:
        """
        Build mappings for a structure that will NOT pass through MODELLER
        (e.g. caps-only repairs).

        Chain IDs are preserved exactly (no MODELLER renaming) and the residue
        mapping is the IDENTITY (orig_num → orig_num). It must NOT pre-renumber
        to 1..N here: residue_mappings is the Step-2 ("before capping") map in
        get_final_identity, and without MODELLER nothing renumbers the existing
        residues before _renumber_structure runs. The final 1..N renumbering is
        done authoritatively by _renumber_structure, which records it in
        final_mappings keyed by the capped file's *original* residue numbers
        (Step 3). Applying a 1..N map here too would compose with that Step-3
        map and double-offset every standard residue — e.g. CYS 200 → 118
        (Step 2) → 37 (final_mappings[118]) — landing the redox site on the
        wrong residue's coordinates. (HETATM cofactors were spared only because
        they are excluded from residue_mappings and so skip Step 2.)

        NOTE: build_modeller_mapping keys residue_mappings by the *predicted*
        MODELLER chain name (chains are renamed to A, B, C... by order). When
        MODELLER does not run, the structure keeps its original chain IDs, so
        using that prediction would mis-key the final-number lookup. This
        identity mapping keeps the keys aligned with the real chain IDs.
        """
        self.chain_mapping = {}
        self.residue_mappings = {}

        model = next(iter(structure), None)
        if model is None:
            return

        for chain in model:
            # Identity chain mapping (no renaming without MODELLER).
            self.chain_mapping[chain.id] = chain.id

            # Identity residue mapping: no renumbering happens before
            # _renumber_structure, so Step 2 must be a no-op. _renumber_structure
            # (Step 3, via final_mappings) is the single source of truth for the
            # final 1..N numbering.
            ordered_nums = [r.id[1] for r in chain if r.id[0] == " "]
            self.residue_mappings[chain.id] = {
                orig_num: orig_num for orig_num in ordered_nums
            }

        self._log_mapping("IDENTITY", self.residue_mappings)

    # NOTE: final_mappings (MODELLER-number → final-number, per renamed chain) is
    # now populated authoritatively by CappingHandler._renumber_structure as it
    # renumbers the capped file, rather than predicted up-front. The previous
    # build_capping_mapping() reconstructed it from global cross-chain numbering
    # and mis-keyed MODELLER-renamed chains, which left renamed-chain redox sites
    # unmapped (every atom resolved to a nonexistent global residue number).

    def get_final_identity(self, original_identity: ResidueIdentity) -> ResidueIdentity:
        """
        Map original residue identity to final identity after all transformations.
        
        MULTI-CHAIN: Uses chain_mapping to handle chain renumbering
        """
        # Step 1: Original chain → New chain
        new_chain = self.chain_mapping.get(original_identity.chain_id, original_identity.chain_id)
        
        # Step 2: Original resnum → MODELLER resnum (within new chain)
        modeller_num = self.residue_mappings.get(new_chain, {}).get(
            original_identity.res_num, original_identity.res_num
        )
        
        # Step 3: MODELLER resnum → Final resnum (if capping occurred in this chain)
        final_num = self.final_mappings.get(new_chain, {}).get(modeller_num, modeller_num)
        
        return ResidueIdentity(
            chain_id=new_chain,
            res_num=final_num,
            res_name=original_identity.res_name,
            insertion_code=original_identity.insertion_code
        )
    
    def get_final_atom_identity(self, original_atom: AtomIdentity) -> AtomIdentity:
        """Map original atom identity to final identity"""
        new_residue = self.get_final_identity(original_atom.residue)
        return AtomIdentity(
            residue=new_residue,
            atom_name=original_atom.atom_name,
            element=original_atom.element
        )
    
    def _get_chain_order(self, structure: Structure) -> List[str]:
        """Get ordered list of chain IDs"""
        chains = []
        for model in structure:
            for chain in model:
                if chain.id not in chains:
                    chains.append(chain.id)
        return chains
    
    def _predict_chain_mapping(self, original_chains: List[str]) -> Dict[str, str]:
        """
        Predict how MODELLER will rename chains.
        
        MULTI-CHAIN LOGIC:
        - MODELLER renames chains to start from 'A' consecutively
        - Original B,C,E → MODELLER A,B,C
        - Handles 26+ chains: A-Z, then AA, AB, ...
        """
        if not original_chains:
            return {}
        
        mapping = {}
        for idx, original_chain in enumerate(original_chains):
            # Standard A-Z for first 26 chains
            if idx < 26:
                new_chain = chr(ord('A') + idx)
            else:
                # For 26+, use AA, AB, AC, ...
                first = chr(ord('A') + (idx // 26) - 1)
                second = chr(ord('A') + (idx % 26))
                new_chain = first + second
            
            mapping[original_chain] = new_chain
            
            if original_chain != new_chain:
                self.console.print(
                    f"[yellow]  Chain {original_chain} will be renumbered to {new_chain}[/yellow]"
                )
        
        return mapping
    
    def _get_present_residues(self, structure: Structure, chain_id: str) -> Set[int]:
        """Get set of residue numbers present in structure for a specific chain"""
        present = set()
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                for residue in chain:
                    # Include ALL residues: standard (" ") and HETATMs ("H")
                    # All go through MODELLER and get renumbered
                    present.add(residue.id[1])
        return present

    def _get_present_residues_ordered(self, structure: Structure, chain_id: str) -> list:
        """Get residue numbers in PDB file order (iteration order from BioPython).

        MODELLER assigns sequential numbers based on the order residues appear
        in the PDB file, NOT sorted by residue number.  Using sorted() would
        produce a wrong mapping when HETATM residue IDs are non-contiguous
        (e.g. 801, 804, 802, 803).
        """
        seen = set()
        ordered = []
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                for residue in chain:
                    res_num = residue.id[1]
                    if res_num not in seen:
                        seen.add(res_num)
                        ordered.append(res_num)
        return ordered
    
    def _log_mapping(self, stage: str, mappings: Dict[str, Dict[int, int]]) -> None:
        """Log mapping for debugging"""
        self.console.print(f"[cyan]  {stage} Mapping:[/cyan]")
        for chain_id, mapping in mappings.items():
            if mapping:
                sample = list(mapping.items())[:3]
                more = f"... (+{len(mapping)-3} more)" if len(mapping) > 3 else ""
                self.console.print(f"[cyan]    Chain {chain_id}: {sample}{more}[/cyan]")


# ============================================================================
# STRUCTURE ANALYZER - Detection only
# ============================================================================

class StructureAnalyzer:
    """Pure detection - no repair logic"""
    
    STD_AMINO_ACIDS = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "MSE": "M"
    }
    
    STD_NUCLEOTIDES = {
        "DA": "A", "DT": "T", "DG": "G", "DC": "C",
        "A": "A", "U": "U", "G": "G", "C": "C"
    }
    
    def __init__(self, structure: Structure, metadata: Optional[Any] = None,
                 filter_selections: Optional[Dict] = None, console: Optional[Console] = None):
        self.structure = structure
        self.metadata = metadata
        self.filter_selections = filter_selections
        self.console = console or Console()
        self.ccd_parser = None  # Lazy load
        
        # Results storage
        self.missing_residues_results = {}
        self.missing_atoms_results = {}
        self.alternate_locations_results = {}
    
    def detect_all_issues(self,
                         missing_residue_methods: Optional[List[str]] = None,
                         missing_atom_methods: Optional[List[str]] = None,
                         altloc_methods: Optional[List[str]] = None) -> Dict[str, Any]:
        """Run selected detection methods"""

        # Default to all methods if not specified
        if missing_residue_methods is None:
            missing_residue_methods = ['remark_465', 'seqres_comparison',
                                      'fasta_comparison', 'sequence_gap']
        if missing_atom_methods is None:
            missing_atom_methods = ['remark_470', 'template_comparison']
        if altloc_methods is None:
            altloc_methods = ['altloc_identifier', 'occupancy_value']

        results = {
            'missing_residues': {},
            'missing_atoms': {},
            'alternate_locations': {}
        }

        # Missing residues detection
        for method in missing_residue_methods:
            if method == 'remark_465':
                results['missing_residues']['remark_465'] = self._detect_missing_residues_remark_465()
            elif method == 'seqres_comparison':
                results['missing_residues']['seqres_comparison'] = self._detect_missing_residues_seqres()
            elif method == 'fasta_comparison':
                results['missing_residues']['fasta_comparison'] = self._detect_missing_residues_fasta()
            elif method == 'sequence_gap':
                results['missing_residues']['sequence_gap'] = self._detect_from_gaps()

        # Missing atoms detection
        for method in missing_atom_methods:
            if method == 'remark_470':
                results['missing_atoms']['remark_470'] = self._detect_missing_atoms_remark_470()
            elif method == 'template_comparison':
                results['missing_atoms']['template_comparison'] = self._detect_missing_atoms_templates()

        # Alternate locations detection
        for method in altloc_methods:
            if method == 'altloc_identifier':
                results['alternate_locations']['altloc_identifier'] = self._detect_alternate_locations_altloc()
            elif method == 'occupancy_value':
                results['alternate_locations']['occupancy_value'] = self._detect_alternate_locations_occupancy()

        # Store results
        self.missing_residues_results = results['missing_residues']
        self.missing_atoms_results = results['missing_atoms']
        self.alternate_locations_results = results['alternate_locations']

        return results
    
    # ========================================================================
    # MISSING RESIDUES DETECTION METHODS
    # ========================================================================
    
    def _detect_missing_residues_remark_465(self) -> Dict[str, List[ResidueIdentity]]:
        """Detect missing residues from REMARK 465 records"""
        missing = defaultdict(list)
        
        if not self.metadata or not hasattr(self.metadata, 'missing_res_records'):
            return {}
        
        for chain_id, residues in self.metadata.missing_res_records.items():
            if self._should_analyze_chain(chain_id):
                for res in residues:
                    missing[chain_id].append(
                        ResidueIdentity(chain_id, res['residue_number'], res['residue_name'])
                    )
        
        return dict(missing)
    
    def _detect_missing_residues_seqres(self) -> Dict[str, List[ResidueIdentity]]:
        """Compare atom-derived sequence to SEQRES records"""
        missing = defaultdict(list)

        if not self.metadata or not hasattr(self.metadata, 'seqres_records'):
            return {}

        for model in self.structure:
            for chain in model:
                chain_id = chain.id

                if not self._should_analyze_chain(chain_id):
                    continue

                if chain_id not in self.metadata.seqres_records:
                    continue

                try:
                    # Get SEQRES sequence
                    seqres_residues = self.metadata.seqres_records[chain_id]
                    seqres_seq = self._residue_list_to_sequence(seqres_residues)

                    # Get atom sequence
                    atom_seq = self._get_chain_sequence(chain)

                    # Align and find gaps
                    alignment = self._align_sequences(seqres_seq, atom_seq)
                    if alignment:
                        # Pass chain to get actual PDB residue numbering
                        gaps = self._find_gaps_in_alignment(alignment, seqres_residues, chain)
                        for res_name, res_num in gaps:
                            missing[chain_id].append(
                                ResidueIdentity(chain_id, res_num, res_name)
                            )
                except Exception as e:
                    import traceback
                    self.console.print(f"[red]DEBUG: Error in SEQRES comparison for chain {chain_id}: {e}[/red]")
                    self.console.print(f"[red]DEBUG: Traceback: {traceback.format_exc()}[/red]")
                    raise

        return dict(missing)
    
    def _detect_missing_residues_fasta(self) -> Dict[str, List[ResidueIdentity]]:
        """Compare atom-derived sequence to downloaded FASTA"""
        missing = defaultdict(list)
        
        if not self.metadata or not hasattr(self.metadata, 'header_info'):
            return {}
        
        pdb_id = self.metadata.header_info.get('pdb_id', '').upper()
        if not pdb_id:
            return {}
        
        # Download FASTA
        fasta_sequences = self._download_fasta(pdb_id)
        if not fasta_sequences:
            return {}
        
        # Process each chain
        for model in self.structure:
            for chain in model:
                chain_id = chain.id
                
                if not self._should_analyze_chain(chain_id):
                    continue
                
                if chain_id not in fasta_sequences:
                    continue
                
                # Get sequences
                fasta_seq = fasta_sequences[chain_id]
                atom_seq = self._get_chain_sequence(chain)
                
                # Align and find gaps
                alignment = self._align_sequences(fasta_seq, atom_seq)
                if alignment:
                    # Convert FASTA sequence back to residue names
                    fasta_residues = [self._one_to_three(aa) for aa in fasta_seq]
                    gaps = self._find_gaps_in_alignment(alignment, fasta_residues)
                    
                    for res_name, res_num in gaps:
                        missing[chain_id].append(
                            ResidueIdentity(chain_id, res_num, res_name)
                        )
        
        return dict(missing)
    
    def _detect_from_gaps(self) -> Dict[str, List[ResidueIdentity]]:
        """Detect missing residues from numbering gaps"""
        missing = defaultdict(list)
        
        for model in self.structure:
            for chain in model:
                if not self._should_analyze_chain(chain.id):
                    continue
                
                residues = sorted([r for r in chain if r.id[0] == " "], 
                                key=lambda r: r.id[1])
                
                for i in range(len(residues) - 1):
                    curr_num = residues[i].id[1]
                    next_num = residues[i + 1].id[1]
                    
                    # Skip if insertion codes present
                    if residues[i].id[2] != " " or residues[i + 1].id[2] != " ":
                        continue
                    
                    if next_num - curr_num > 1:
                        # Gap detected
                        for gap_num in range(curr_num + 1, next_num):
                            missing[chain.id].append(
                                ResidueIdentity(chain.id, gap_num, "UNK")
                            )
        
        return dict(missing)
    
    # ========================================================================
    # MISSING ATOMS DETECTION METHODS
    # ========================================================================
    
    def _detect_missing_atoms_remark_470(self) -> Dict[str, List[Tuple[ResidueIdentity, str]]]:
        """Detect missing atoms from REMARK 470 records"""
        missing = defaultdict(list)
        
        if not self.metadata or not hasattr(self.metadata, 'remark_records'):
            return {}
        
        if 470 not in self.metadata.remark_records:
            return {}
        
        current_chain = None
        current_resnum = None
        current_resname = None
        
        for line in self.metadata.remark_records[470]:
            # Look for residue header
            res_match = re.search(r'(\S+)\s+(\S+)\s+(\d+)', line)
            if res_match:
                current_resname = res_match.group(1)
                current_chain = res_match.group(2)
                current_resnum = int(res_match.group(3))
                
                if not self._should_analyze_chain(current_chain):
                    current_chain = None
                    continue
            
            # Look for atom lines
            if current_chain and current_resnum and current_resname:
                atom_match = re.search(r'\s+(\S+)', line)
                if atom_match and not re.search(r'^\s*M RES', line):
                    atom_name = atom_match.group(1)
                    res_identity = ResidueIdentity(current_chain, current_resnum, current_resname)
                    missing[current_chain].append((res_identity, atom_name))
        
        return dict(missing)
    
    def _detect_missing_atoms_templates(self) -> Dict[str, List[Tuple[ResidueIdentity, str]]]:
        """Detect missing atoms by comparing to CCD templates"""
        missing = defaultdict(list)
        
        # Initialize CCD parser
        if self.ccd_parser is None:
            try:
                from proprep.structure_prep.chem_comp_dict_fetcher import CCDParser
                self.ccd_parser = CCDParser(use_cache=True)
            except ImportError:
                self.console.print("[yellow]CCD parser not available, skipping template comparison[/yellow]")
                return {}
        
        # Get unique residue types
        unique_residues = set()
        for model in self.structure:
            for chain in model:
                if not self._should_analyze_chain(chain.id):
                    continue
                for residue in chain:
                    unique_residues.add(residue.resname.strip().upper())
        
        # Pre-fetch templates
        self.console.print(f"[grey50]Fetching templates for {len(unique_residues)} residue types...[/grey50]")
        template_cache = {}
        for res_name in unique_residues:
            template_data = self.ccd_parser.get_residue_data(res_name)
            template_cache[res_name] = template_data
        
        # Check each residue
        for model in self.structure:
            for chain in model:
                if not self._should_analyze_chain(chain.id):
                    continue
                
                for residue in chain:
                    res_name = residue.resname.strip().upper()
                    template_data = template_cache.get(res_name, {})
                    
                    if 'error' in template_data or 'atoms' not in template_data:
                        continue
                    
                    # Get expected atoms
                    template_atoms = {atom['atom_id'] for atom in template_data['atoms']}
                    actual_atoms = {atom.name.strip() for atom in residue}
                    
                    # Find missing (skip hydrogens)
                    missing_atoms = template_atoms - actual_atoms
                    for atom_name in missing_atoms:
                        if not atom_name.startswith('H'):
                            res_identity = ResidueIdentity(
                                chain.id, residue.id[1], res_name, residue.id[2]
                            )
                            missing[chain.id].append((res_identity, atom_name))
        
        return dict(missing)
    
    # ========================================================================
    # ALTERNATE LOCATIONS DETECTION METHODS
    # ========================================================================
    
    def _detect_alternate_locations_altloc(self) -> Dict[str, Dict[str, Set[str]]]:
        """Detect atoms with alternate location identifiers"""
        alt_locations = defaultdict(lambda: defaultdict(set))
        
        for model in self.structure:
            for chain in model:
                if not self._should_analyze_chain(chain.id):
                    continue
                
                for residue in chain:
                    res_key = f"{residue.resname}_{residue.id[1]}"
                    
                    for atom in residue:
                        if hasattr(atom, 'is_disordered') and atom.is_disordered():
                            if hasattr(atom, 'child_dict'):
                                for altloc_id in atom.child_dict.keys():
                                    if altloc_id.strip():
                                        alt_locations[chain.id][res_key].add(altloc_id.strip())
                        else:
                            altloc = atom.altloc.strip() if hasattr(atom, 'altloc') else ""
                            if altloc:
                                alt_locations[chain.id][res_key].add(altloc)
        
        # Filter to only residues with multiple alternates
        filtered = {}
        for chain_id, residues in alt_locations.items():
            filtered_residues = {k: v for k, v in residues.items() if len(v) > 1}
            if filtered_residues:
                filtered[chain_id] = filtered_residues
        
        return filtered
    
    def _detect_alternate_locations_occupancy(self) -> Dict[str, Dict[str, List[Tuple[str, float]]]]:
        """Detect atoms with partial occupancy (indicating alternates)"""
        partial_occ = defaultdict(lambda: defaultdict(list))
        
        for model in self.structure:
            for chain in model:
                if not self._should_analyze_chain(chain.id):
                    continue
                
                for residue in chain:
                    res_key = f"{residue.resname}_{residue.id[1]}"
                    
                    # Group atoms by name
                    atom_groups = defaultdict(list)
                    for atom in residue:
                        atom_groups[atom.name.strip()].append(atom)
                    
                    # Check for atoms with multiple conformations
                    for atom_name, atoms in atom_groups.items():
                        if len(atoms) > 1:  # Multiple conformations
                            for atom in atoms:
                                if hasattr(atom, 'occupancy') and atom.occupancy < 1.0:
                                    altloc = atom.altloc.strip() if hasattr(atom, 'altloc') else 'default'
                                    partial_occ[chain.id][res_key].append((atom_name, atom.occupancy))
        
        return dict(partial_occ)
    
    # ========================================================================
    # HELPER METHODS
    # ========================================================================
    
    def _should_analyze_chain(self, chain_id: str) -> bool:
        """Check if chain should be analyzed based on filter selections"""
        if not self.filter_selections:
            return True
        return str(chain_id) in self.filter_selections
    
    def should_analyze_residue(self, chain_id: str, residue: Residue) -> bool:
        """Check if residue should be analyzed based on filter selections"""
        if not self.filter_selections:
            return True
        
        chain_id = str(chain_id)
        
        if chain_id not in self.filter_selections:
            return False
        
        # Get residue component type
        comp_type = self._get_residue_type(residue)
        
        # Check if this component type is in filter selections
        chain_filters = self.filter_selections.get(chain_id, {})
        if comp_type not in chain_filters:
            return False
        
        # Check if this residue number is in the selected residues
        residue_num = residue.id[1]
        selected_residues = chain_filters.get(comp_type, [])
        
        # Convert selected_residues to list if needed
        if not isinstance(selected_residues, list):
            if isinstance(selected_residues, set):
                selected_residues = list(selected_residues)
            else:
                selected_residues = []
        
        return residue_num in selected_residues
    
    def _get_residue_type(self, residue: Residue) -> str:
        """Get the component type of a residue"""
        res_name = residue.resname.strip().upper()
        
        # Simple classification
        if residue.id[0] == " ":  # Standard residue
            if res_name in self.STD_AMINO_ACIDS:
                return "amino_acid"
            elif res_name in self.STD_NUCLEOTIDES:
                if res_name.startswith("D"):
                    return "dna"
                else:
                    return "rna"
        
        # Other residues
        if res_name in ["HOH", "WAT"]:
            return "water"
        elif res_name in ["NAG", "MAN", "FUC", "BGC"]:
            return "carbohydrate"
        elif res_name in ["FAD", "NAD", "FMN", "HEM", "ATP"]:
            return "cofactor"
        elif res_name in ["CA", "MG", "ZN", "FE", "MN"]:
            return "metal_ion"
        elif residue.id[0] != " ":
            return "ligand"
        
        return "unknown"
    
    def _get_chain_sequence(self, chain: Chain) -> str:
        """Get amino acid sequence from chain"""
        sequence = ""
        for residue in chain:
            if residue.id[0] == " ":
                res_name = residue.resname
                if res_name in self.STD_AMINO_ACIDS:
                    sequence += self.STD_AMINO_ACIDS[res_name]
                elif res_name in self.STD_NUCLEOTIDES:
                    sequence += self.STD_NUCLEOTIDES[res_name]
        return sequence
    
    def _residue_list_to_sequence(self, residue_list: List[str]) -> str:
        """Convert residue name list to single letter sequence"""
        sequence = ""
        for res_name in residue_list:
            if res_name in self.STD_AMINO_ACIDS:
                sequence += self.STD_AMINO_ACIDS[res_name]
            elif res_name in self.STD_NUCLEOTIDES:
                sequence += self.STD_NUCLEOTIDES[res_name]
            else:
                sequence += "X"
        return sequence
    
    def _one_to_three(self, one_letter: str) -> str:
        """Convert one letter code to three letter"""
        for three, one in self.STD_AMINO_ACIDS.items():
            if one == one_letter:
                return three
        return "UNK"
    
    def _download_fasta(self, pdb_id: str) -> Dict[str, str]:
        """Download FASTA sequence from RCSB"""
        pdb_id = pdb_id.lower()
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/download"
        
        try:
            with urllib.request.urlopen(url) as response:
                content = response.read().decode('utf-8')
            
            sequences = {}
            current_chains = []
            current_seq = ""
            
            for line in content.splitlines():
                if line.startswith(">"):
                    # Save previous
                    if current_chains and current_seq:
                        for chain in current_chains:
                            sequences[chain] = current_seq
                    
                    # Parse new header
                    header_match = re.search(r'\|Chains\s+([^|]+)\|', line)
                    if header_match:
                        chains_part = header_match.group(1)
                        chain_matches = re.findall(r'([A-Za-z0-9])(?:\[auth [A-Za-z0-9]\])?', chains_part)
                        current_chains = chain_matches if chain_matches else []
                    else:
                        current_chains = []
                    
                    current_seq = ""
                else:
                    if current_chains:
                        current_seq += line.strip()
            
            # Save last
            if current_chains and current_seq:
                for chain in current_chains:
                    sequences[chain] = current_seq
            
            return sequences
        except Exception as e:
            self.console.print(f"[yellow]Could not download FASTA: {e}[/yellow]")
            return {}
    
    def _align_sequences(self, seq1: str, seq2: str) -> Optional[Tuple[str, str, Dict]]:
        """Align two sequences using Needleman-Wunsch"""
        try:
            from Bio.Align import PairwiseAligner
            
            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5
            
            alignments = aligner.align(seq1, seq2)
            if alignments:
                alignment = alignments[0]

                # Get aligned sequences with gaps using format method
                alignment_str = format(alignment)
                lines = alignment_str.strip().split('\n')

                # The format is:
                # target    0 SEQUENCE... length
                # match       |||||||...
                # query     0 SEQUENCE... length
                # Find the lines with actual sequences
                seq1_aligned = None
                seq2_aligned = None

                for i, line in enumerate(lines):
                    if line.startswith('target'):
                        # Extract sequence from this line (skip 'target', position number, and trailing info)
                        parts = line.split()
                        if len(parts) >= 3:
                            seq1_aligned = parts[2] if seq1_aligned is None else seq1_aligned + parts[2]
                    elif line.startswith('query'):
                        parts = line.split()
                        if len(parts) >= 3:
                            seq2_aligned = parts[2] if seq2_aligned is None else seq2_aligned + parts[2]

                if not seq1_aligned or not seq2_aligned:
                    self.console.print(f"[red]DEBUG: Could not extract aligned sequences from format output[/red]")
                    return None

                # Create mapping
                mapping = {}
                idx1 = idx2 = 0
                for i in range(len(seq1_aligned)):
                    if seq1_aligned[i] != "-":
                        if seq2_aligned[i] != "-":
                            mapping[idx1] = idx2
                            idx2 += 1
                        idx1 += 1
                    elif seq2_aligned[i] != "-":
                        idx2 += 1

                return seq1_aligned, seq2_aligned, mapping
        except ImportError:
            try:
                from Bio import pairwise2
                alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -2, -0.5)
                if alignments:
                    alignment = alignments[0]
                    # Similar processing...
                    return None
            except ImportError:
                pass
        
        return None
    
    def _find_gaps_in_alignment(self, alignment: Tuple[str, str, Dict],
                               reference_residues: List[str],
                               chain=None) -> List[Tuple[str, int]]:
        """Find gaps in alignment and map to residue names/numbers

        Args:
            alignment: Tuple of (seq1_aligned, seq2_aligned, mapping)
            reference_residues: List of residue names from SEQRES
            chain: BioPython Chain object to get actual PDB residue numbering

        Returns:
            List of (residue_name, residue_number) tuples for missing residues
        """
        seq1_aligned, seq2_aligned, mapping = alignment
        gaps = []

        # Build mapping from SEQRES position to actual PDB residue number
        # Strategy: Calculate offset between SEQRES (1-indexed) and PDB numbering
        seqres_to_pdb_num = {}
        pdb_offset = 0

        if chain:
            # Get all present residues from the chain
            present_residues = [r for r in chain if r.id[0] == " "]

            if present_residues:
                # Calculate offset: Find first present residue in SEQRES and compare to PDB number
                first_present_seqres_0idx = None
                temp_idx = 0
                for i in range(len(seq1_aligned)):
                    if seq1_aligned[i] != "-":
                        if seq2_aligned[i] != "-":
                            first_present_seqres_0idx = temp_idx
                            break
                        temp_idx += 1

                if first_present_seqres_0idx is not None:
                    first_present_seqres_1idx = first_present_seqres_0idx + 1
                    first_present_pdb_num = present_residues[0].id[1]
                    pdb_offset = first_present_pdb_num - first_present_seqres_1idx

        # Now collect gaps with their residue numbers using offset
        seq1_idx = 0
        for i in range(len(seq1_aligned)):
            if seq1_aligned[i] != "-" and seq2_aligned[i] == "-":
                # Gap in seq2 (missing residue in structure)
                if seq1_idx < len(reference_residues):
                    res_name = reference_residues[seq1_idx]

                    # Calculate PDB residue number: SEQRES position (1-indexed) + offset
                    seqres_1indexed = seq1_idx + 1
                    res_num = seqres_1indexed + pdb_offset
                    gaps.append((res_name, res_num))

            if seq1_aligned[i] != "-":
                seq1_idx += 1

        return gaps
    
    # ========================================================================
    # VISUALIZATION & REPORTING METHODS
    # ========================================================================
    
    def display_sequence_view(self, chain_id: str, missing_residues: List[ResidueIdentity],
                             missing_atoms_positions: Optional[Set[int]] = None,
                             altloc_positions: Optional[Set[int]] = None,
                             mutation_positions: Optional[Set[int]] = None,
                             nonstandard_mutation_positions: Optional[Set[int]] = None) -> None:
        """Display ASCII sequence view with various issues highlighted"""
        chain = None
        for model in self.structure:
            if chain_id in model:
                chain = model[chain_id]
                break

        if not chain:
            return

        # Get all residue numbers (present + missing)
        present = {r.id[1]: r.resname for r in chain if r.id[0] == " "}
        missing = {r.res_num: r.res_name for r in missing_residues}

        all_nums = sorted(set(present.keys()) | set(missing.keys()))
        if not all_nums:
            return

        min_num, max_num = min(all_nums), max(all_nums)

        # Convert None to empty sets
        missing_atoms_positions = missing_atoms_positions or set()
        altloc_positions = altloc_positions or set()
        mutation_positions = mutation_positions or set()
        nonstandard_mutation_positions = nonstandard_mutation_positions or set()

        self.console.print(f"\n[bold]Chain {chain_id} Sequence[/bold]")

        # Display in blocks of 60 residues
        BLOCK_SIZE = 60
        for start in range(min_num, max_num + 1, BLOCK_SIZE):
            end = min(start + BLOCK_SIZE - 1, max_num)
            block_length = end - start + 1

            # Build sequence line first
            sequence = ""
            for num in range(start, end + 1):
                if num in present:
                    aa = self.STD_AMINO_ACIDS.get(present[num], 'X')
                    sequence += aa
                elif num in missing:
                    aa = self.STD_AMINO_ACIDS.get(missing[num], 'X')
                    sequence += aa  # Plain for ruler calculation
                else:
                    sequence += "-"

            seq_line = f"{start:>4}-{end:<3} {sequence}"

            # Calculate prefix length (everything before the sequence)
            prefix_length = len(seq_line) - len(sequence)

            # Create position ruler with numbers at every 10th position
            ruler = " " * prefix_length
            for i in range(10, block_length + 1, 10):
                number = str(i)
                current_pos = len(ruler) - prefix_length
                spaces_needed = i - current_pos - len(number)
                ruler += " " * spaces_needed + number

            # Pad ruler to match sequence length
            while len(ruler) - prefix_length < len(sequence):
                ruler += " "

            # Print the block
            self.console.print(ruler)

            # Re-build sequence line with proper formatting
            seq_line = f"{start:>4}-{end:<3} "
            for num in range(start, end + 1):
                if num in present:
                    aa = self.STD_AMINO_ACIDS.get(present[num], 'X')
                    # Determine formatting based on issues (priority order)
                    if num in nonstandard_mutation_positions:
                        seq_line += f"[bold bright_red]{aa}[/bold bright_red]"  # Pending non-standard mutation (orange)
                    elif num in mutation_positions:
                        seq_line += f"[bold magenta]{aa}[/bold magenta]"  # Pending standard mutation
                    elif num in altloc_positions:
                        seq_line += f"[bold yellow]{aa}[/bold yellow]"  # Alternate location
                    elif num in missing_atoms_positions:
                        seq_line += f"[bold cyan]{aa}[/bold cyan]"  # Missing atoms
                    else:
                        seq_line += aa  # Normal
                elif num in missing:
                    aa = self.STD_AMINO_ACIDS.get(missing[num], 'X')
                    seq_line += f"[bold green]{aa}[/bold green]"  # Missing residue
                else:
                    seq_line += "-"
            self.console.print(seq_line)
            self.console.print()  # Empty line between blocks

        # Legend
        self.console.print("\n[bold]Legend:[/bold]")
        self.console.print("Normal text: Complete residues")
        self.console.print("[bold green]Green[/bold green]: Missing residues")
        self.console.print("[bold cyan]Cyan[/bold cyan]: Missing atoms")
        self.console.print("[bold yellow]Yellow[/bold yellow]: Alternate locations")
        self.console.print("[bold magenta]Magenta[/bold magenta]: Pending standard mutations")
        self.console.print("[bold bright_red]Orange[/bold bright_red]: Pending non-standard mutations")
        self.console.print("'-': Not in sequence")

    def display_all_chain_sequences(self, missing_by_chain: Dict[str, List[ResidueIdentity]],
                                   pending_mutations: Optional[List[Tuple]] = None,
                                   pending_nonstandard_mutations: Optional[List[Tuple]] = None) -> None:
        """Display sequence views for all chains with missing residues and other issues"""
        for chain_id in sorted(missing_by_chain.keys()):
            # Collect missing atoms positions for this chain
            missing_atoms_positions = set()
            for method_results in self.missing_atoms_results.values():
                if chain_id in method_results:
                    for atom_identity, _ in method_results[chain_id]:
                        missing_atoms_positions.add(atom_identity.res_num)

            # Collect alternate location positions for this chain
            altloc_positions = set()
            for method_results in self.alternate_locations_results.values():
                if chain_id in method_results:
                    for res_key in method_results[chain_id].keys():
                        # Parse res_key which is like "SER_37"
                        parts = res_key.split('_')
                        if len(parts) >= 2:
                            try:
                                pos = int(parts[-1])
                                altloc_positions.add(pos)
                            except ValueError:
                                pass

            # Collect standard mutation positions for this chain
            mutation_positions = set()
            if pending_mutations:
                for mut_chain_id, res_num, _, _ in pending_mutations:
                    if mut_chain_id == chain_id:
                        mutation_positions.add(res_num)

            # Collect non-standard mutation positions for this chain
            nonstandard_mutation_positions = set()
            if pending_nonstandard_mutations:
                for mut_chain_id, res_num, _, _, _ in pending_nonstandard_mutations:
                    if mut_chain_id == chain_id:
                        nonstandard_mutation_positions.add(res_num)

            # Only show if there are any issues in this chain
            if (missing_by_chain[chain_id] or missing_atoms_positions or altloc_positions or
                mutation_positions or nonstandard_mutation_positions):
                self.display_sequence_view(
                    chain_id,
                    missing_by_chain[chain_id],
                    missing_atoms_positions,
                    altloc_positions,
                    mutation_positions,
                    nonstandard_mutation_positions
                )

    def display_missing_residues_report(self) -> None:
        """Display detailed missing residues report"""
        if not self.missing_residues_results:
            self.console.print("[yellow]No missing residues detected[/yellow]")
            return

        for method_name, chain_results in self.missing_residues_results.items():
            if not chain_results:
                continue

            self.console.print(f"\n[bold underline]Missing Residues - {method_name}[/bold underline]")

            for chain_id, residues in chain_results.items():
                if not residues:
                    continue

                sorted_residues = sorted(residues, key=lambda r: r.res_num)

                table = Table(title=f"Chain {chain_id}")
                table.add_column("Position", style="cyan")
                table.add_column("Residue", style="green")
                table.add_column("Type", style="yellow")

                for res in sorted_residues:
                    table.add_row(str(res.res_num), res.res_name, "Missing")

                self.console.print(table)
    
    def display_missing_atoms_report(self) -> None:
        """Display detailed missing atoms report"""
        if not self.missing_atoms_results:
            self.console.print("[yellow]No missing atoms detected[/yellow]")
            return

        for method_name, chain_results in self.missing_atoms_results.items():
            if not chain_results:
                continue

            self.console.print(f"\n[bold underline]Missing Atoms - {method_name}[/bold underline]")

            for chain_id, atoms in chain_results.items():
                if not atoms:
                    continue

                sorted_atoms = sorted(atoms, key=lambda a: (a[0].res_num, a[1]))

                table = Table(title=f"Chain {chain_id}")
                table.add_column("Residue", style="cyan")
                table.add_column("Atom", style="green")

                for res_identity, atom_name in sorted_atoms:
                    res_str = f"{res_identity.res_name} {res_identity.res_num}"
                    table.add_row(res_str, atom_name)

                self.console.print(table)
    
    def display_alternate_locations_report(self) -> None:
        """Display alternate locations report with occupancy information"""
        if not self.alternate_locations_results:
            return

        # Check if there's any actual data to display
        has_data = False
        for method_results in self.alternate_locations_results.values():
            if method_results:
                for chain_residues in method_results.values():
                    if chain_residues:
                        has_data = True
                        break
            if has_data:
                break

        if not has_data:
            return

        self.console.print("\n[bold underline]Alternate Locations[/bold underline]")

        for chain_id, residues in self.alternate_locations_results.get('altloc_identifier', {}).items():
            if not residues:
                continue

            table = Table(title=f"Chain {chain_id}")
            table.add_column("Residue", style="cyan")
            table.add_column("Position", style="yellow")
            table.add_column("Alternates", style="green")
            table.add_column("Occupancies", style="magenta")

            # Get occupancy information from structure
            for res_key, altlocs in residues.items():
                # Parse residue info from key (e.g., "SER_37")
                parts = res_key.split('_')
                if len(parts) == 2:
                    res_name, res_num_str = parts
                    try:
                        res_num = int(res_num_str)
                    except ValueError:
                        continue

                    # Find the residue in the structure to get occupancy
                    occupancies = {}
                    for model in self.structure:
                        if chain_id in model:
                            chain = model[chain_id]
                            for residue in chain:
                                if residue.id[1] == res_num:
                                    for atom in residue:
                                        if hasattr(atom, 'is_disordered') and atom.is_disordered():
                                            if hasattr(atom, 'child_dict'):
                                                for altloc_id, alt_atom in atom.child_dict.items():
                                                    if altloc_id.strip():
                                                        if altloc_id.strip() not in occupancies:
                                                            occupancies[altloc_id.strip()] = []
                                                        if hasattr(alt_atom, 'occupancy'):
                                                            occupancies[altloc_id.strip()].append(alt_atom.occupancy)
                                        else:
                                            altloc = atom.altloc.strip() if hasattr(atom, 'altloc') else ""
                                            if altloc:
                                                if altloc not in occupancies:
                                                    occupancies[altloc] = []
                                                if hasattr(atom, 'occupancy'):
                                                    occupancies[altloc].append(atom.occupancy)
                                    break
                            break

                    # Format occupancy information
                    occ_str = ""
                    for altloc in sorted(altlocs):
                        if altloc in occupancies and occupancies[altloc]:
                            avg_occ = sum(occupancies[altloc]) / len(occupancies[altloc])
                            occ_str += f"{altloc}:{avg_occ:.2f} "
                        else:
                            occ_str += f"{altloc}:? "

                    alt_str = ", ".join(sorted(altlocs))
                    table.add_row(res_name, str(res_num), alt_str, occ_str.strip())

            self.console.print(table)

    def display_mutations_report(self,
                                pending_mutations: Optional[List[Tuple]] = None,
                                pending_nonstandard_mutations: Optional[List[Tuple]] = None) -> None:
        """Display detailed mutations report"""
        has_standard = pending_mutations and len(pending_mutations) > 0
        has_nonstandard = pending_nonstandard_mutations and len(pending_nonstandard_mutations) > 0

        if not has_standard and not has_nonstandard:
            return

        self.console.print("\n[bold underline]Pending Mutations[/bold underline]")

        # Display standard mutations
        if has_standard:
            self.console.print("\n[bold magenta]Standard Mutations[/bold magenta] (Applied by MODELLER)")

            # Group by chain
            mutations_by_chain = {}
            for chain_id, res_num, from_aa, to_aa in pending_mutations:
                if chain_id not in mutations_by_chain:
                    mutations_by_chain[chain_id] = []
                mutations_by_chain[chain_id].append((res_num, from_aa, to_aa))

            for chain_id in sorted(mutations_by_chain.keys()):
                table = Table(title=f"Chain {chain_id}")
                table.add_column("Position", style="cyan")
                table.add_column("From", style="red")
                table.add_column("To", style="green")
                table.add_column("Type", style="yellow")

                for res_num, from_aa, to_aa in sorted(mutations_by_chain[chain_id], key=lambda x: x[0]):
                    # Classify mutation type
                    mutation_type = self._classify_mutation_type(from_aa, to_aa)
                    table.add_row(str(res_num), from_aa, to_aa, mutation_type)

                self.console.print(table)

        # Display non-standard mutations
        if has_nonstandard:
            self.console.print("\n[bold bright_red]Non-Standard Mutations[/bold bright_red] (Atom removal + TLEaP rebuilding)")

            # Group by chain
            ns_mutations_by_chain = {}
            for chain_id, res_num, from_aa, to_aa, atoms_to_keep in pending_nonstandard_mutations:
                if chain_id not in ns_mutations_by_chain:
                    ns_mutations_by_chain[chain_id] = []
                ns_mutations_by_chain[chain_id].append((res_num, from_aa, to_aa, atoms_to_keep))

            for chain_id in sorted(ns_mutations_by_chain.keys()):
                table = Table(title=f"Chain {chain_id}")
                table.add_column("Position", style="cyan")
                table.add_column("From", style="red")
                table.add_column("To", style="green")
                table.add_column("Atoms to Keep", style="yellow")

                for res_num, from_aa, to_aa, atoms_to_keep in sorted(ns_mutations_by_chain[chain_id], key=lambda x: x[0]):
                    atoms_str = ", ".join(atoms_to_keep)
                    table.add_row(str(res_num), from_aa, to_aa, atoms_str)

                self.console.print(table)

    def _classify_mutation_type(self, from_aa: str, to_aa: str) -> str:
        """Classify mutation type for display"""
        # Simple classification
        hydrophobic = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
        polar = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
        charged_pos = {"ARG", "LYS", "HIS"}
        charged_neg = {"ASP", "GLU"}
        small = {"GLY", "ALA"}

        if from_aa == to_aa:
            return "Unchanged"
        elif from_aa in small and to_aa in small:
            return "Small→Small"
        elif from_aa in hydrophobic and to_aa in hydrophobic:
            return "Hydrophobic"
        elif from_aa in polar and to_aa in polar:
            return "Polar"
        elif from_aa in charged_pos and to_aa in charged_pos:
            return "Positive"
        elif from_aa in charged_neg and to_aa in charged_neg:
            return "Negative"
        elif (from_aa in charged_pos or from_aa in charged_neg) and (to_aa in charged_pos or to_aa in charged_neg):
            return "Charge swap"
        elif (from_aa in hydrophobic and to_aa in polar) or (from_aa in polar and to_aa in hydrophobic):
            return "Polar↔Hydrophobic"
        else:
            return "Other"

    def display_summary_report(self,
                               pending_mutations: Optional[List[Tuple]] = None,
                               pending_nonstandard_mutations: Optional[List[Tuple]] = None) -> None:
        """Display comprehensive summary of all issues"""
        table = Table(title="Structure Completeness Summary")
        table.add_column("Issue Type", style="cyan")
        table.add_column("Method", style="yellow")
        table.add_column("Count", style="green")
        table.add_column("Chains", style="magenta")

        # Missing residues - show all methods (even if 0)
        if self.missing_residues_results:
            for method, chain_results in self.missing_residues_results.items():
                total = sum(len(residues) for residues in chain_results.values())
                chains = ", ".join(sorted(chain_results.keys())) if chain_results and total > 0 else "-"
                style = "green" if total == 0 else "red"
                table.add_row("Missing Residues", method, f"[{style}]{total}[/{style}]", chains)

        # Missing atoms - show all methods (even if 0)
        if self.missing_atoms_results:
            for method, chain_results in self.missing_atoms_results.items():
                total = sum(len(atoms) for atoms in chain_results.values())
                chains = ", ".join(sorted(chain_results.keys())) if chain_results and total > 0 else "-"
                style = "green" if total == 0 else "red"
                table.add_row("Missing Atoms", method, f"[{style}]{total}[/{style}]", chains)

        # Alternate locations - show all methods (even if 0)
        if self.alternate_locations_results:
            for method, chain_results in self.alternate_locations_results.items():
                total = sum(len(residues) for residues in chain_results.values())
                chains = ", ".join(sorted(chain_results.keys())) if chain_results and total > 0 else "-"
                style = "green" if total == 0 else "red"
                table.add_row("Alternate Locations", method, f"[{style}]{total}[/{style}]", chains)

        # Pending standard mutations
        if pending_mutations is not None:
            mutation_count = len(pending_mutations)
            if mutation_count > 0:
                # Get unique chains from mutations
                mutation_chains = sorted(set(mut[0] for mut in pending_mutations))
                chains_str = ", ".join(mutation_chains)
                table.add_row("Pending Standard Mutations", "user_specified", f"[yellow]{mutation_count}[/yellow]", chains_str)
            else:
                table.add_row("Pending Standard Mutations", "user_specified", "[green]0[/green]", "-")

        # Pending non-standard mutations
        if pending_nonstandard_mutations is not None:
            ns_mutation_count = len(pending_nonstandard_mutations)
            if ns_mutation_count > 0:
                # Get unique chains from non-standard mutations
                ns_mutation_chains = sorted(set(mut[0] for mut in pending_nonstandard_mutations))
                chains_str = ", ".join(ns_mutation_chains)
                table.add_row("Pending Non-Standard Mutations", "user_specified", f"[yellow]{ns_mutation_count}[/yellow]", chains_str)
            else:
                table.add_row("Pending Non-Standard Mutations", "user_specified", "[green]0[/green]", "-")

        self.console.print(table)

        # Overall status
        has_missing_residues = any(len(residues) > 0 for chain_results in self.missing_residues_results.values() for residues in chain_results.values())
        has_missing_atoms = any(len(atoms) > 0 for chain_results in self.missing_atoms_results.values() for atoms in chain_results.values())
        has_alternate_locations = any(len(residues) > 0 for chain_results in self.alternate_locations_results.values() for residues in chain_results.values())
        has_mutations = pending_mutations and len(pending_mutations) > 0

        has_issues = has_missing_residues or has_missing_atoms or has_alternate_locations or has_mutations

        if not has_issues:
            self.console.print("\n[green]✓ No structural issues found![/green]")
            self.console.print("Structure is complete with:")
            self.console.print("  • No missing residues")
            self.console.print("  • No missing atoms")
            self.console.print("  • No alternate locations")
            self.console.print("  • No pending mutations")


# ============================================================================
# REPAIR ORCHESTRATOR - User interaction and workflow management
# ============================================================================

class RepairOrchestrator:
    """Manages repair planning and execution workflow"""

    def __init__(self, console: Console, processor=None):
        self.console = console
        self.processor = processor
    
    def validate_mutations(self,
                          mutations: List[Tuple[str, int, str, str]],
                          detection_results: Dict[str, Any]) -> Tuple[List, List[str]]:
        """
        Validate mutations against missing residues analysis.
        
        Returns:
            (valid_mutations, warnings)
        """
        valid = []
        warnings = []
        
        # Get missing residues
        missing = self._extract_missing_residues(detection_results)
        missing_positions = defaultdict(set)
        for chain_id, residues in missing.items():
            for res in residues:
                missing_positions[chain_id].add(res.res_num)
        
        # Check each mutation
        for chain_id, res_num, from_aa, to_aa in mutations:
            if res_num in missing_positions.get(chain_id, set()):
                warnings.append(
                    f"Mutation {chain_id}:{res_num}:{from_aa}→{to_aa} is on a missing residue - will be modeled"
                )
            valid.append((chain_id, res_num, from_aa, to_aa))
        
        return valid, warnings
    
    def create_repair_plan(self, 
                        detection_results: Dict[str, Any],
                        structure: Structure,
                        pending_mutations: List[Tuple]) -> RepairPlan:
        """
        Interactive repair planning - creates plan through user interaction.
        
        Returns complete RepairPlan with user decisions.
        """
        plan = RepairPlan()
        
        # Add pending mutations to plan
        plan.mutations = pending_mutations.copy() if pending_mutations else []
        
        # Group missing residues into segments
        missing_residues = self._extract_missing_residues(detection_results)
        segments = self._group_into_segments(missing_residues, structure)
        
        if not segments and not plan.mutations:
            self.console.print("[green]No missing residues to repair or mutations to apply.[/green]")
            self.console.print("[grey50]Note: Missing atoms will be built automatically by TLeaP during topology generation.[/grey50]")
            return plan
        
        # Show sequence visualization FIRST (helps user understand the structure)
        if missing_residues:
            self.console.print("\n[bold cyan]═══ Sequence Visualization ═══[/bold cyan]")
            analyzer = StructureAnalyzer(structure, None, None, self.console)
            analyzer.display_all_chain_sequences(missing_residues, plan.mutations)
        
        # Show overview of segments
        self._display_segments_overview(segments)
        
        # Show mutations table if any
        if plan.mutations:
            self._display_mutations_table(plan.mutations)
        
        # Interactive decision for each segment
        self.console.print("\n[bold]Interactive Repair Planning[/bold]")
        self.console.print("For each missing segment, choose to Fill, Partial fill, Cap, or Skip.")
        self.console.print("Note: Single-residue gaps offer Fill or a TER record (chain break) instead of capping.\n")
        
        for segment in segments:
            action = self._prompt_segment_action(segment)

            # Handle partial fill (returns dict) vs regular actions (returns enum)
            if isinstance(action, dict) and action.get('action') == RepairAction.PARTIAL_FILL:
                # Create a new segment with only the selected residues
                selected_residues = action['selected_residues']
                selected_nums = action['selected_nums']

                partial_segment = MissingSegment(
                    chain_id=segment.chain_id,
                    residues=selected_residues,
                    start_num=min(selected_nums),
                    end_num=max(selected_nums),
                    is_terminal=segment.is_terminal and (min(selected_nums) == segment.start_num or max(selected_nums) == segment.end_num),
                    terminal_type=segment.terminal_type if segment.is_terminal else None
                )
                plan.segments_to_fill.append(partial_segment)

                # Add skipped portions to skip list
                skipped_nums = [r.res_num for r in segment.residues if r.res_num not in selected_nums]
                if skipped_nums:
                    skipped_residues = [r for r in segment.residues if r.res_num in skipped_nums]
                    skipped_segment = MissingSegment(
                        chain_id=segment.chain_id,
                        residues=skipped_residues,
                        start_num=min(skipped_nums),
                        end_num=max(skipped_nums),
                        is_terminal=segment.is_terminal,
                        terminal_type=segment.terminal_type
                    )
                    plan.segments_to_skip.append(skipped_segment)

            elif action == RepairAction.FILL:
                plan.segments_to_fill.append(segment)
            elif action == RepairAction.CAP:
                plan.segments_to_cap.append(segment)
            elif action == RepairAction.TER:
                plan.segments_to_ter.append(segment)
            else:  # SKIP
                plan.segments_to_skip.append(segment)

        # Offer ACE/NME terminal capping BEFORE presenting the summary, so the
        # user's cap choices are part of the plan they review and confirm —
        # rather than being asked *after* approving the plan (the old ordering).
        # Termini already scheduled for capping above are skipped inside, so
        # this never double-offers.
        self._offer_terminal_capping(structure, plan)

        # Show final plan summary
        self._display_plan_summary(plan)

        # Visualize the plan: anchor on residues that exist (flanking
        # for fill/skip, the surviving terminal for cap) so the user
        # can see *where* each action will happen before confirming.
        # Three labels, one per action, so each is independently
        # toggleable in the rep manager.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer

            def _flanking_clauses(segments):
                pairs = set()
                for s in segments:
                    pairs.add((s.chain_id, s.start_num - 1))
                    pairs.add((s.chain_id, s.end_num + 1))
                return [f"(:{c} and {n})" for c, n in sorted(pairs)]

            def _cap_anchor_clauses(segments):
                # Cap segments are terminal — the existing terminal
                # residue is what gets the ACE/NME attached, so anchor
                # there.
                pairs = set()
                for s in segments:
                    if s.terminal_type == 'N':
                        pairs.add((s.chain_id, s.end_num + 1))
                    else:
                        pairs.add((s.chain_id, s.start_num - 1))
                return [f"(:{c} and {n})" for c, n in sorted(pairs)]

            fill_clauses = _flanking_clauses(plan.segments_to_fill)
            cap_clauses = _cap_anchor_clauses(plan.segments_to_cap)
            skip_clauses = _flanking_clauses(plan.segments_to_skip)

            _viewer.unhighlight("fixer_plan_fill")
            _viewer.unhighlight("fixer_plan_cap")
            _viewer.unhighlight("fixer_plan_skip")
            if fill_clauses:
                _viewer.highlight(" or ".join(fill_clauses), style="ball+stick",
                                  color="#33a02c", label="fixer_plan_fill")
            if cap_clauses:
                _viewer.highlight(" or ".join(cap_clauses), style="ball+stick",
                                  color="#ff7f00", label="fixer_plan_cap")
            if skip_clauses:
                _viewer.highlight(" or ".join(skip_clauses), style="ball+stick",
                                  color="#e31a1c", label="fixer_plan_skip")
        except Exception:
            pass

        # Confirm
        if not confirm_with_context(
            processor=self.processor,
            prompt="\nProceed with this repair plan?",
            default=True,
            module="Structure Completeness - Repair",
            description="Confirm repair plan"
        ):
            self.console.print("[yellow]Repair cancelled[/yellow]")

            # Even when the user declines the full plan, an unfilled INTERNAL gap
            # leaves two residues that tLEaP will happily bond across, building a
            # spurious long bond. Offer to drop a TER record at each internal
            # break so the gap reads as a chain end. Terminal gaps need no TER
            # (nothing follows them to bond to).
            internal_segments = [s for s in segments if not s.is_terminal]
            if internal_segments and confirm_with_context(
                processor=self.processor,
                prompt=(
                    f"\nInsert TER records at the {len(internal_segments)} internal "
                    "break(s) so tLEaP won't build long bonds across the gaps?"
                ),
                default=True,
                module="Structure Completeness - Repair",
                description="Insert TER records at declined internal gaps",
            ):
                ter_plan = RepairPlan()
                ter_plan.segments_to_ter = list(internal_segments)
                return ter_plan

            return None

        return plan

    def _enumerate_protein_termini(self, structure: Structure) -> List[Dict[str, Any]]:
        """
        Identify the N- and C-terminus of every protein chain, independent of
        whether any residues are missing. This is what lets capping be offered
        on an otherwise-complete structure.

        Returns one dict per terminus:
            {'chain_id', 'terminal_type' ('N'|'C'), 'res_num', 'res_name',
             'already_capped'}
        """
        protein_resnames = set(StructureAnalyzer.STD_AMINO_ACIDS.keys())
        termini: List[Dict[str, Any]] = []

        # First model only, consistent with the rest of this module.
        model = next(iter(structure), None)
        if model is None:
            return termini

        for chain in model:
            # Standard amino-acid residues only, in sequence order.
            protein_residues = sorted(
                [r for r in chain
                 if r.id[0] == " " and r.resname.strip() in protein_resnames],
                key=lambda r: r.id[1]
            )
            if not protein_residues:
                continue

            # Detect caps already present anywhere in the chain so we never
            # offer to double-cap. ACE caps the N-terminus; NME/NHE the C.
            chain_resnames = {r.resname.strip() for r in chain}
            n_capped = "ACE" in chain_resnames
            c_capped = "NME" in chain_resnames or "NHE" in chain_resnames

            first_res = protein_residues[0]
            last_res = protein_residues[-1]

            termini.append({
                'chain_id': chain.id,
                'terminal_type': 'N',
                'res_num': first_res.id[1],
                'res_name': first_res.resname.strip(),
                'already_capped': n_capped,
            })
            termini.append({
                'chain_id': chain.id,
                'terminal_type': 'C',
                'res_num': last_res.id[1],
                'res_name': last_res.resname.strip(),
                'already_capped': c_capped,
            })

        return termini

    def _offer_terminal_capping(self, structure: Structure, plan: RepairPlan) -> bool:
        """
        Interactively offer to add ACE/NME caps to protein chain termini that
        are NOT already capped and NOT already scheduled for capping by the
        missing-residue plan. Selected caps are appended to
        `plan.segments_to_cap` as synthetic cap-only segments (empty residue
        list; the downstream cap inserter keys off chain_id/terminal_type only).

        Requires no MODELLER — capping is a pure PDB edit.

        Returns True if at least one cap was added to the plan.
        """
        termini = self._enumerate_protein_termini(structure)
        if not termini:
            return False

        # Termini the missing-residue plan already caps — don't double-offer.
        already_planned = {(s.chain_id, s.terminal_type) for s in plan.segments_to_cap}

        # Group by chain so each chain is a single N/C decision.
        by_chain: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(dict)
        for t in termini:
            by_chain[t['chain_id']][t['terminal_type']] = t

        def _available(entry: Optional[Dict[str, Any]], ttype: str, chain_id: str) -> bool:
            if entry is None or entry['already_capped']:
                return False
            return (chain_id, ttype) not in already_planned

        caps_added = 0
        header_shown = False

        for chain_id in sorted(by_chain.keys()):
            ends = by_chain[chain_id]
            n_entry = ends.get('N')
            c_entry = ends.get('C')

            n_avail = _available(n_entry, 'N', chain_id)
            c_avail = _available(c_entry, 'C', chain_id)
            if not n_avail and not c_avail:
                continue

            if not header_shown:
                self.console.print("\n[bold cyan]═══ Terminal Capping (ACE/NME) ═══[/bold cyan]")
                self.console.print(
                    "[grey50]ACE caps the N-terminus, NME caps the C-terminus. Cap atoms get "
                    "placeholder geometry that tLEaP finalizes during topology generation.[/grey50]"
                )
                header_shown = True

            # Show the chain's terminal residues for context.
            def _desc(entry: Optional[Dict[str, Any]]) -> str:
                if entry is None:
                    return "—"
                tag = " [grey50](already capped)[/grey50]" if entry['already_capped'] else ""
                return f"{entry['res_name']}{entry['res_num']}{tag}"

            self.console.print(
                f"\n[bold]Chain {chain_id}[/bold] — N-term: {_desc(n_entry)}, "
                f"C-term: {_desc(c_entry)}"
            )

            # Focus the viewer on the still-cappable terminal residue(s).
            try:
                from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                clauses = []
                if n_avail:
                    clauses.append(f"(:{chain_id} and {n_entry['res_num']})")
                if c_avail:
                    clauses.append(f"(:{chain_id} and {c_entry['res_num']})")
                _viewer.unhighlight("cap_offer_focus")
                if clauses:
                    _viewer.highlight(" or ".join(clauses), style="ball+stick",
                                      color="#ff7f00", label="cap_offer_focus")
            except Exception:
                pass

            # Prompt shape depends on which ends are still available.
            if n_avail and c_avail:
                # options_map suppresses Rich's inline choice list (see
                # prompt_with_context), so print the legend ourselves.
                options_map = {
                    "n": f"N-terminus only — ACE before {n_entry['res_name']}{n_entry['res_num']}",
                    "c": f"C-terminus only — NME after {c_entry['res_name']}{c_entry['res_num']}",
                    "both": "Both termini",
                    "s": "Skip this chain (default)",
                }
                for key in ("n", "c", "both", "s"):
                    self.console.print(f"  [cyan]{key:>4}[/cyan]  {options_map[key]}")
                choice = prompt_with_context(
                    processor=self.processor,
                    prompt=f"Cap which terminus of chain {chain_id}?",
                    choices=["n", "c", "both", "s"],
                    default="s",
                    module="Structure Completeness - Capping",
                    description=f"Terminal capping for chain {chain_id}",
                    options_map=options_map,
                ).lower()
                do_n = choice in ("n", "both")
                do_c = choice in ("c", "both")
            elif n_avail:
                do_n = confirm_with_context(
                    processor=self.processor,
                    prompt=f"Cap the N-terminus of chain {chain_id} "
                           f"(ACE before {n_entry['res_name']}{n_entry['res_num']})?",
                    default=False,
                    module="Structure Completeness - Capping",
                    description=f"N-terminal cap for chain {chain_id}",
                )
                do_c = False
            else:  # c_avail only
                do_c = confirm_with_context(
                    processor=self.processor,
                    prompt=f"Cap the C-terminus of chain {chain_id} "
                           f"(NME after {c_entry['res_name']}{c_entry['res_num']})?",
                    default=False,
                    module="Structure Completeness - Capping",
                    description=f"C-terminal cap for chain {chain_id}",
                )
                do_n = False

            if do_n:
                plan.segments_to_cap.append(MissingSegment(
                    chain_id=chain_id,
                    residues=[],
                    start_num=n_entry['res_num'],
                    end_num=n_entry['res_num'],
                    is_terminal=True,
                    terminal_type='N',
                ))
                caps_added += 1
                self.console.print(
                    f"  [green]✓ Will add ACE cap to chain {chain_id} N-terminus[/green]"
                )
            if do_c:
                plan.segments_to_cap.append(MissingSegment(
                    chain_id=chain_id,
                    residues=[],
                    start_num=c_entry['res_num'],
                    end_num=c_entry['res_num'],
                    is_terminal=True,
                    terminal_type='C',
                ))
                caps_added += 1
                self.console.print(
                    f"  [green]✓ Will add NME cap to chain {chain_id} C-terminus[/green]"
                )

        # Clear the transient viewer focus.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            _viewer.unhighlight("cap_offer_focus")
        except Exception:
            pass

        return caps_added > 0

    def _extract_missing_residues(self, detection_results: Dict) -> Dict[str, List[ResidueIdentity]]:
        """Extract missing residues from detection results"""
        missing = defaultdict(list)
        
        # Use first available method
        for method in ['remark_465', 'seqres_comparison', 'fasta_comparison', 'sequence_gap']:
            if method in detection_results.get('missing_residues', {}):
                for chain_id, residues in detection_results['missing_residues'][method].items():
                    missing[chain_id].extend(residues)
                break
        
        return dict(missing)
    
    def _group_into_segments(self, 
                            missing_residues: Dict[str, List[ResidueIdentity]], 
                            structure: Structure) -> List[MissingSegment]:
        """Group missing residues into contiguous segments"""
        segments = []
        
        for chain_id, residues in missing_residues.items():
            if not residues:
                continue
            
            # Sort by residue number
            sorted_residues = sorted(residues, key=lambda r: r.res_num)
            
            # Get chain bounds
            chain_bounds = self._get_chain_bounds(structure, chain_id)
            
            # Group into contiguous segments
            current_segment = [sorted_residues[0]]
            
            for res in sorted_residues[1:]:
                if res.res_num == current_segment[-1].res_num + 1:
                    current_segment.append(res)
                else:
                    # Finish current segment
                    segment = self._create_segment(current_segment, chain_bounds)
                    segments.append(segment)
                    current_segment = [res]
            
            # Finish last segment
            if current_segment:
                segment = self._create_segment(current_segment, chain_bounds)
                segments.append(segment)
        
        return segments
    
    def _get_chain_bounds(self, structure: Structure, chain_id: str) -> Tuple[int, int]:
        """Get first and last residue numbers in chain"""
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                residues = [r.id[1] for r in chain if r.id[0] == " "]
                if residues:
                    return min(residues), max(residues)
        return 0, 0
    
    def _create_segment(self, residues: List[ResidueIdentity], 
                       chain_bounds: Tuple[int, int]) -> MissingSegment:
        """Create MissingSegment with terminal classification"""
        start_num = residues[0].res_num
        end_num = residues[-1].res_num
        chain_id = residues[0].chain_id
        
        # Determine if terminal
        is_n_terminal = start_num < chain_bounds[0]
        is_c_terminal = end_num > chain_bounds[1]
        
        is_terminal = is_n_terminal or is_c_terminal
        terminal_type = 'N' if is_n_terminal else 'C' if is_c_terminal else None
        
        return MissingSegment(
            chain_id=chain_id,
            residues=residues,
            start_num=start_num,
            end_num=end_num,
            is_terminal=is_terminal,
            terminal_type=terminal_type
        )
    
    def _display_segments_overview(self, segments: List[MissingSegment]) -> None:
        """Display overview table of all segments"""
        if not segments:
            return
        
        self.console.print("\n[bold]Missing Residue Segments[/bold]")
        
        table = Table()
        table.add_column("Chain", style="cyan")
        table.add_column("Segment", style="yellow")
        table.add_column("Length", style="green")
        table.add_column("Type", style="magenta")
        
        for seg in segments:
            seg_str = f"{seg.start_num}-{seg.end_num}" if seg.length > 1 else str(seg.start_num)
            type_str = f"{seg.terminal_type}-terminal" if seg.is_terminal else "Internal"
            table.add_row(seg.chain_id, seg_str, str(seg.length), type_str)
        
        self.console.print(table)
    
    def _display_mutations_table(self, mutations: List[Tuple]) -> None:
        """Display pending mutations table"""
        self.console.print("\n[bold cyan]Pending Mutations[/bold cyan]")
        
        table = Table()
        table.add_column("Chain", style="cyan")
        table.add_column("Position", style="yellow")
        table.add_column("From", style="red")
        table.add_column("To", style="green")
        
        for chain_id, res_num, from_aa, to_aa in mutations:
            table.add_row(chain_id, str(res_num), from_aa, to_aa)
        
        self.console.print(table)
    
    def _prompt_segment_action(self, segment: MissingSegment) -> RepairAction:
        """Prompt user for action on a single segment"""
        # Single-residue gaps: fill it (MODELLER) or just insert a TER record so
        # tLEaP treats the break as a chain end instead of bonding across it.
        # Capping a 1-residue internal gap is not offered: two ACE/NME caps in a
        # single-residue hole make no chemical sense, and TER is the lighter,
        # MODELLER-free way to stop the long-bond artifact.
        if segment.is_single:
            self.console.print(
                f"\n[bold]Single-residue gap: Chain {segment.chain_id}, "
                f"position {segment.start_num} "
                f"({segment.residues[0].res_name})[/bold]"
            )
            if not HAS_MODELLER:
                self.console.print(
                    "  [grey50]MODELLER is unavailable — inserting a TER record "
                    "(chain break) so tLEaP won't build a long bond across the "
                    "gap. Install MODELLER to fill single-residue gaps.[/grey50]"
                )
                return RepairAction.TER

            self.console.print("  \\[f] Fill with MODELLER")
            self.console.print("  \\[t] Insert TER (chain break, no fill)")
            choice = prompt_with_context(
                processor=self.processor,
                prompt="Action",
                choices=["f", "t"],
                default="f",
                module="Structure Completeness - Repair",
                description=f"Repair action for single-residue gap at {segment.chain_id}{segment.start_num}",
                options_map={
                    "f": "Fill with MODELLER",
                    "t": "Insert TER (chain break, no fill)",
                },
            ).lower()
            return RepairAction.FILL if choice == "f" else RepairAction.TER
        
        # Multi-residue segment - ask user
        self.console.print(f"\n[bold]Segment: Chain {segment.chain_id}, residues {segment.start_num}-{segment.end_num}[/bold]")
        self.console.print(f"  Length: {segment.length} residues")
        
        res_names = [r.res_name for r in segment.residues]
        self.console.print(f"  Residues: {', '.join(res_names)}")
        
        type_str = f"{segment.terminal_type}-terminal" if segment.is_terminal else "Internal"
        self.console.print(f"  Type: {type_str}")
        
        # Provide suggestion
        if segment.is_terminal:
            self.console.print("  [grey50]Suggestion: Terminal gaps are often capped with ACE/NME[/grey50]")
        else:
            self.console.print("  [grey50]Suggestion: Internal gaps are usually filled[/grey50]")

        # List the choices explicitly: prompt_with_context suppresses the inline
        # choice list whenever options_map is passed (it assumes the caller
        # printed them), so without this the prompt would show only "Action (f):".
        self.console.print("  \\[f] Fill    \\[p] Partial fill    \\[c] Cap    \\[s] Skip")

        # Focus the viewer on this segment's flanking residues (yellow
        # ball+stick) and halo the chain (blue) so the user sees both
        # the local context for the gap and which chain is being
        # discussed. Single label across iterations — re-firing
        # replaces the prior segment's focus.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            flank = []
            if segment.terminal_type != 'N':
                flank.append(f"(:{segment.chain_id} and {segment.start_num - 1})")
            if segment.terminal_type != 'C':
                flank.append(f"(:{segment.chain_id} and {segment.end_num + 1})")
            _viewer.unhighlight("fixer_segment_focus")
            _viewer.unhighlight("fixer_segment_chain")
            if flank:
                _viewer.highlight(" or ".join(flank), style="ball+stick",
                                  color="#ffff00", label="fixer_segment_focus")
            _viewer.highlight(f":{segment.chain_id}", style="halo",
                              color="#1f78b4", label="fixer_segment_chain")
        except Exception:
            pass

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Action",
            choices=["f", "p", "c", "s"],
            default="f",
            module="Structure Completeness - Repair",
            description=f"Repair action for {segment.chain_id} segment ({segment.length} residues)",
            options_map={"f": "Fill", "p": "Partial fill", "c": "Cap", "s": "Skip"}
        ).lower()

        action_map = {'f': RepairAction.FILL, 'p': RepairAction.PARTIAL_FILL, 'c': RepairAction.CAP, 's': RepairAction.SKIP}
        action = action_map[choice]

        # Handle partial fill
        if action == RepairAction.PARTIAL_FILL:
            return self._prompt_partial_fill(segment)

        # Confirm skip
        if action == RepairAction.SKIP:
            if not confirm_with_context(
                processor=self.processor,
                prompt="Skip this segment? It will remain missing.",
                default=False,
                module="Structure Completeness - Repair",
                description="Confirm skip segment"
            ):
                return self._prompt_segment_action(segment)  # Ask again

        return action

    def _parse_residue_selection(self, selection: str, segment: MissingSegment) -> List[int]:
        """
        Parse residue selection and validate it's a continuous sequence.

        Args:
            selection: User input like "25-27" or "25,26,27"
            segment: The segment being selected from

        Returns:
            List of selected residue numbers

        Raises:
            ValueError: If selection is invalid or not continuous
        """
        # Get available residue numbers from segment
        available_nums = [r.res_num for r in segment.residues]

        selected_nums = []

        # Handle range notation (e.g., "25-27")
        if '-' in selection:
            parts = selection.split('-')
            if len(parts) != 2:
                raise ValueError("Invalid range format. Use: start-end (e.g., 25-27)")

            try:
                start = int(parts[0].strip())
                end = int(parts[1].strip())
            except ValueError:
                raise ValueError("Range values must be integers")

            if start > end:
                raise ValueError(f"Invalid range: {start} > {end}")

            selected_nums = list(range(start, end + 1))

        # Handle comma-separated notation (e.g., "25,26,27")
        elif ',' in selection:
            try:
                selected_nums = [int(x.strip()) for x in selection.split(',')]
            except ValueError:
                raise ValueError("All values must be integers")

        # Single number
        else:
            try:
                selected_nums = [int(selection.strip())]
            except ValueError:
                raise ValueError("Must be a valid integer or range")

        # Validate all selected residues are in the segment
        for num in selected_nums:
            if num not in available_nums:
                raise ValueError(f"Residue {num} is not in this segment ({segment.start_num}-{segment.end_num})")

        # Validate continuous sequence
        if len(selected_nums) > 1:
            sorted_nums = sorted(selected_nums)
            for i in range(len(sorted_nums) - 1):
                if sorted_nums[i + 1] - sorted_nums[i] != 1:
                    raise ValueError(f"Selection must be continuous. Gap found between {sorted_nums[i]} and {sorted_nums[i + 1]}")

        return sorted(selected_nums)

    def _prompt_partial_fill(self, segment: MissingSegment) -> dict:
        """
        Prompt user to select which residues to fill from the segment.

        Returns:
            Dictionary with action type and partial segment details
        """
        self.console.print("\n[bold cyan]Partial Fill Mode[/bold cyan]")
        self.console.print("Select a continuous range of residues to fill from this segment.")
        self.console.print(f"Available residues: {segment.start_num}-{segment.end_num}")
        self.console.print("\n[grey50]Examples:[/grey50]")
        self.console.print(f"  [grey50]{segment.start_num+2}-{segment.end_num}[/grey50] (skip first 2)")
        self.console.print(f"  [grey50]{segment.start_num},{segment.start_num+1}[/grey50] (fill first 2 only)")

        while True:
            selection = prompt_with_context(
                processor=self.processor,
                prompt=f"\n[green]Residues to fill[/green] ({segment.start_num}-{segment.end_num})",
                default=f"{segment.start_num}-{segment.end_num}",
                module="Structure Completeness - Partial Fill",
                description="Select residues for partial fill"
            )

            try:
                selected_nums = self._parse_residue_selection(selection, segment)

                # Show what will be filled
                selected_residues = [r for r in segment.residues if r.res_num in selected_nums]
                selected_names = [r.res_name for r in selected_residues]

                self.console.print(f"\n[green]Selected residues:[/green] {', '.join(selected_names)}")
                self.console.print(f"  Range: {min(selected_nums)}-{max(selected_nums)}")
                self.console.print(f"  Count: {len(selected_nums)} residues")

                # Show what will be skipped
                skipped_nums = [n for n in [r.res_num for r in segment.residues] if n not in selected_nums]
                if skipped_nums:
                    skipped_residues = [r for r in segment.residues if r.res_num in skipped_nums]
                    skipped_names = [r.res_name for r in skipped_residues]
                    self.console.print(f"\n[yellow]Skipped residues:[/yellow] {', '.join(skipped_names)}")
                    if min(skipped_nums) < min(selected_nums):
                        self.console.print(f"  N-terminal skip: {min(skipped_nums)}-{max([n for n in skipped_nums if n < min(selected_nums)])}")
                    if max(skipped_nums) > max(selected_nums):
                        self.console.print(f"  C-terminal skip: {min([n for n in skipped_nums if n > max(selected_nums)])}-{max(skipped_nums)}")

                # Confirm selection
                if confirm_with_context(
                    processor=self.processor,
                    prompt="Confirm partial fill?",
                    default=True,
                    module="Structure Completeness - Partial Fill",
                    description="Confirm partial fill selection"
                ):
                    # Return a dictionary with partial fill details
                    return {
                        'action': RepairAction.PARTIAL_FILL,
                        'selected_residues': selected_residues,
                        'selected_nums': selected_nums
                    }
                else:
                    # Let them try again
                    continue

            except ValueError as e:
                self.console.print(f"[red]Invalid selection:[/red] {e}")
                self.console.print("Please try again.\n")
                continue

    def _display_plan_summary(self, plan: RepairPlan) -> None:
        """Display final repair plan summary"""
        self.console.print("\n[bold]═══ Repair Plan Summary ═══[/bold]")
        
        table = Table(title="Actions to be Performed")
        table.add_column("Action", style="cyan")
        table.add_column("Count", style="yellow")
        table.add_column("Details", style="green")
        
        if plan.segments_to_fill:
            fill_count = sum(s.length for s in plan.segments_to_fill)
            chains = {s.chain_id for s in plan.segments_to_fill}
            table.add_row(
                "Fill (MODELLER)",
                f"{len(plan.segments_to_fill)} segments",
                f"{fill_count} residues, Chains {', '.join(chains)}"
            )
        
        if plan.segments_to_cap:
            chains = {s.chain_id for s in plan.segments_to_cap}
            table.add_row(
                "Cap (ACE/NME)",
                f"{len(plan.segments_to_cap)} segments",
                f"Chains {', '.join(chains)}"
            )
        
        if plan.mutations:
            chains = {m[0] for m in plan.mutations}
            table.add_row(
                "Mutations",
                f"{len(plan.mutations)} mutations",
                f"Chains {', '.join(chains)}"
            )
        
        if plan.segments_to_ter:
            chains = {s.chain_id for s in plan.segments_to_ter}
            table.add_row(
                "TER (chain break)",
                f"{len(plan.segments_to_ter)} segments",
                f"Unfilled; Chains {', '.join(chains)}"
            )

        if plan.segments_to_skip:
            table.add_row(
                "Skip",
                f"{len(plan.segments_to_skip)} segments",
                "Will remain missing"
            )

        self.console.print(table)
        
        # Processing pipeline
        self.console.print("\n[bold]Processing Pipeline:[/bold]")
        if plan.needs_modeller:
            self.console.print("  1. Run MODELLER (fill residues + apply mutations)")
        if plan.has_caps:
            step_num = 2 if plan.needs_modeller else 1
            self.console.print(f"  {step_num}. Add ACE/NME caps")
        if plan.has_ter:
            step_num = 1 + int(plan.needs_modeller) + int(plan.has_caps)
            self.console.print(f"  {step_num}. Insert TER records at unfilled breaks")
        self.console.print("  Final. Synchronize RedoxSites")


# ============================================================================
# MODELLER INTERFACE - MODELLER-specific operations
# ============================================================================

class ModellerInterface:
    """Handles all MODELLER-specific operations"""
    
    def __init__(self, console: Console):
        self.console = console
        self.std_amino_acids = StructureAnalyzer.STD_AMINO_ACIDS
    
    def build_sequences(self,
                       structure: Structure,
                       plan: RepairPlan,
                       structure_metadata: Dict) -> Dict[str, Tuple[str, str]]:
        """
        Build MODELLER template and target sequences.
        
        CRITICAL: Must match original working implementation exactly.
        Returns: {chain_id: (template_seq, target_seq)}
        """
        chain_sequences = {}
        
        # Get present residues from structure
        present_residues = self._get_present_residues(structure)
        
        # Get residues to fill from plan
        residues_to_fill = defaultdict(set)
        for segment in plan.segments_to_fill:
            for res in segment.residues:
                residues_to_fill[res.chain_id].add(res.res_num)
        
        # Get mutations by chain and position
        mutations_by_chain = defaultdict(dict)
        for chain_id, res_num, from_aa, to_aa in plan.mutations:
            mutations_by_chain[chain_id][res_num] = (from_aa, to_aa)
        
        # Build sequences for each chain
        for chain_id in self._get_chain_ids(structure):
            # Get complete residue map from SEQRES with correct PDB numbering
            complete_residues = self._build_residue_map(
                chain_id,
                structure_metadata.get('seqres', {}).get(chain_id, []),
                structure_metadata.get('dbref', {}).get(chain_id, {}),
                structure,
                structure_metadata.get('missing_res_records', {})
            )
            
            # Build template and target sequences
            template_seq, target_seq = self._build_chain_sequences(
                chain_id,
                complete_residues,
                present_residues.get(chain_id, set()),
                residues_to_fill.get(chain_id, set()),
                mutations_by_chain.get(chain_id, {}),
                structure
            )

            # CRITICAL: Add offset gaps for non-standard PDB numbering
            # The offset represents phantom residues that exist in numbering but not in SEQRES
            # E.g., if SEQRES starts at pos 1 but PDB numbering starts at 2, offset=1
            # We've already added gaps for missing SEQRES residues, so only add offset gaps
            chain_obj = None
            for model in structure:
                if chain_id in model:
                    chain_obj = model[chain_id]
                    break

            if chain_obj:
                chain_residues = [r for r in chain_obj if r.id[0] == " "]
                seqres = structure_metadata.get('seqres', {}).get(chain_id, [])

                if chain_residues and seqres:
                    first_pdb_num = chain_residues[0].id[1]
                    # SEQRES represents positions 1-N
                    # If first_pdb_num > len(seqres), something is very wrong
                    # The offset is: first_pdb_num - (number of SEQRES positions before first present)

                    # Count how many SEQRES positions exist before the first present residue
                    # The first present residue is at SEQRES position (first_pdb_num - offset)
                    # Offset = first_pdb_num - first_present_seqres_1indexed
                    # We already calculated this in SEQRES comparison! It's the +1 offset

                    # Simple approach: first SEQRES position should map to PDB position 1
                    # But if PDB starts at position 2, we have 1 phantom position
                    # If PDB starts at 45 and we have 44 SEQRES positions before it, offset is 1

                    # Actually: check if first position in complete_residues matches first_pdb_num
                    if complete_residues:
                        min_complete_pos = min(complete_residues.keys())
                        offset_gaps = min_complete_pos - 1  # Phantom positions before first SEQRES position
                        if offset_gaps > 0:
                            template_seq = "-" * offset_gaps + template_seq
                            target_seq = "-" * offset_gaps + target_seq

            chain_sequences[chain_id] = (template_seq, target_seq)
        
        # CRITICAL: Validate sequences against structure
        self._validate_sequences(chain_sequences, structure)
        
        return chain_sequences
    
    def _validate_sequences(self, chain_sequences: Dict[str, Tuple[str, str]], 
                           structure: Structure) -> None:
        """
        Validate that sequence lengths match actual structure.
        This prevents cryptic MODELLER failures.
        """
        for chain_id, (template_seq, target_seq) in chain_sequences.items():
            # Count what's actually in structure for this chain
            actual_residue_count = 0
            actual_water_count = 0
            actual_other_hetatm_count = 0
            
            for model in structure:
                if chain_id in model:
                    chain = model[chain_id]
                    for residue in chain:
                        if residue.id[0] == " ":  # Standard residue
                            actual_residue_count += 1
                        else:  # HETATM
                            if residue.resname in ['HOH', 'WAT']:
                                actual_water_count += 1
                            else:
                                actual_other_hetatm_count += 1
            
            # Count positions in template
            template_residue_count = sum(1 for c in template_seq if c not in ["-", ".", "w"])
            template_other_hetatm = template_seq.count(".")
            template_water = template_seq.count("w")
            
            # Calculate totals
            actual_total = actual_residue_count + actual_water_count + actual_other_hetatm_count
            template_total = template_residue_count + template_other_hetatm + template_water
            
            if template_total != actual_total:
                self.console.print(f"[red]VALIDATION FAILED for chain {chain_id}:[/red]")
                self.console.print(f"  Structure: {actual_residue_count} residues + "
                                 f"{actual_other_hetatm_count} HETATM + "
                                 f"{actual_water_count} water = {actual_total} total")
                self.console.print(f"  Template: {template_residue_count} residues + "
                                 f"{template_other_hetatm} HETATM + "
                                 f"{template_water} water = {template_total} total")
                raise ValueError(
                    f"Sequence length mismatch for chain {chain_id}: "
                    f"template has {template_total} but structure has {actual_total}"
                )
    
    def create_alignment_file(self, chain_sequences: Dict[str, Tuple[str, str]],
                             output_file: str, structure: Structure) -> None:
        """Create MODELLER alignment file

        Note: We use blank residue ranges (::) to let MODELLER auto-detect from the PDB.
        This is critical because:
        1. Alignment sequences include heteroatoms (.) and waters (W)
        2. PDB numbering may be non-standard (not starting at 1)
        3. MODELLER counts alignment positions differently than PDB residue numbers

        IMPORTANT: Chain order in the alignment MUST match the order chains appear
        in the PDB file (structure order), NOT alphabetical order. MODELLER reads
        chains from the PDB in file order and matches them positionally to alignment
        sequences.
        """
        # Use structure chain order (matching PDB file order) instead of sorted
        chain_ids = [chain.id for model in structure for chain in model
                     if chain.id in chain_sequences]
        first_chain = chain_ids[0] if chain_ids else 'A'
        last_chain = chain_ids[-1] if chain_ids else 'Z'

        with open(output_file, 'w') as f:
            # Template entry - use :: to auto-detect residue range from PDB
            f.write(">P1;template\n")
            f.write(f"structure:input.pdb::{first_chain}::{last_chain}:::: \n")
            
            # Template sequences with chain separators
            for i, chain_id in enumerate(chain_ids):
                template_seq, _ = chain_sequences[chain_id]
                f.write(template_seq)
                if i < len(chain_ids) - 1:
                    f.write("/")
            f.write("*\n\n")
            
            # Target entry
            f.write(">P1;target\n")
            f.write(f"sequence:target:1:{first_chain}:999:{last_chain}:::: \n")
            
            # Target sequences with chain separators
            for i, chain_id in enumerate(chain_ids):
                _, target_seq = chain_sequences[chain_id]
                f.write(target_seq)
                if i < len(chain_ids) - 1:
                    f.write("/")
            f.write("*\n")
    
    def run_modeller(self, work_dir: str, input_pdb: str,
                    alignment_file: str,
                    built_residues: Optional[Set[Tuple[str, int]]] = None
                    ) -> Tuple[bool, str, Optional[Structure]]:
        """Execute MODELLER.

        built_residues: optional set of (chain_id, res_num) in MODELLER
        output numbering identifying residues that MODELLER actually built
        (filled gaps / new residues). When provided, a per-residue DOPE
        assessment of those residues is reported (see _assess_built_region).
        """
        if not HAS_MODELLER:
            return False, "MODELLER not available", None

        import modeller
        from modeller.automodel import AutoModel, assess

        try:
            # Change to work directory so MODELLER outputs there
            original_dir = os.getcwd()
            os.chdir(work_dir)

            # Set up environment
            env = modeller.Environ()
            env.io.atom_files_directory = ['.']
            env.io.hetatm = True
            env.io.water = True

            env.libs.topology.read(file="${LIB}/top_heav.lib")
            env.libs.parameters.read(file="${LIB}/par.lib")

            # Redirect output
            import sys
            from io import StringIO

            old_stdout = sys.stdout
            old_stderr = sys.stderr
            captured_output = StringIO()

            # Freeze all resolved atoms: optimize ONLY the residues MODELLER
            # builds, holding everything else rigid. Without this, AutoModel's
            # default full-model refinement (conjugate gradients + MD/simulated
            # annealing) relaxes EXISTING atoms too, which can drag a metal
            # site's coordinating residue (e.g. a Cys Sγ) out of coordination.
            # A frozen environment lets us fill gaps without perturbing the site
            # — and makes MODELLER safe to run before OR after parameterization.
            _built = built_residues or set()

            class _GapOnlyAutoModel(AutoModel):
                def select_atoms(self):
                    sel = modeller.Selection()
                    added = 0
                    for chain_id, res_num in _built:
                        try:
                            sel.add(self.residues[f'{res_num}:{chain_id}'])
                            added += 1
                        except Exception:
                            # Residue spec didn't resolve (numbering/chain edge
                            # case) — skip it rather than abort the whole repair.
                            pass
                    if added == 0:
                        # Nothing resolved: fall back to default (optimize all)
                        # so MODELLER still has a selection to work on.
                        return modeller.Selection(self)
                    return sel

            model_cls = _GapOnlyAutoModel if _built else AutoModel

            try:
                sys.stdout = captured_output
                sys.stderr = captured_output

                # Create and run model
                mdl = model_cls(env,
                              alnfile=os.path.basename(alignment_file),
                              knowns="template",
                              sequence="target")
                mdl.starting_model = 1
                mdl.ending_model = 1
                # Compute published, calibrated assessment scores directly via
                # the MODELLER API (no stdout scraping): normalized DOPE z-score
                # (Shen & Sali 2006) and GA341 (John & Sali 2003), plus raw DOPE
                # and molpdf for provenance. Results land in mdl.outputs.
                mdl.assess_methods = (assess.normalized_dope, assess.GA341,
                                      assess.DOPE)
                mdl.make()

            finally:
                sys.stdout = old_stdout
                sys.stderr = old_stderr
                # Always change back to original directory
                os.chdir(original_dir)

            modeller_output = captured_output.getvalue()

            # Parse output file (now in work_dir)
            output_pdb = os.path.join(work_dir, "target.B99990001.pdb")
            if not os.path.exists(output_pdb):
                return False, "MODELLER output file not found", None
            
            # Parse structure
            parser = PDBParser(QUIET=True)
            repaired_structure = parser.get_structure("repaired", output_pdb)
            
            # Rigorous, literature-calibrated quality assessment.
            # Pass 1: global scores read straight from the MODELLER API.
            # Pass 2: per-residue DOPE of the residues MODELLER actually built,
            #         relative to the structure's own experimental residues.
            assessment = self._extract_assessment(mdl)
            built_region = None
            if built_residues:
                built_region = self._assess_built_region(
                    env, output_pdb, built_residues
                )
            self._display_quality_assessment(
                assessment, repaired_structure, built_region
            )

            return True, "Success", repaired_structure
            
        except Exception as e:
            return False, f"MODELLER error: {str(e)}", None
    
    def _get_chain_ids(self, structure: Structure) -> List[str]:
        """Get sorted list of chain IDs"""
        chains = []
        for model in structure:
            for chain in model:
                if chain.id not in chains:
                    chains.append(chain.id)
        return sorted(chains)
    
    def _get_present_residues(self, structure: Structure) -> Dict[str, Set[int]]:
        """Get present residue numbers per chain"""
        present = defaultdict(set)
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Include ALL residues: standard and HETATMs
                    present[chain.id].add(residue.id[1])
        return dict(present)
    
    def _build_residue_map(self, chain_id: str, seqres: List[str],
                          dbref: Dict, structure: Structure,
                          missing_res_records: Dict) -> Dict[int, str]:
        """
        Build complete residue map from SEQRES with correct PDB numbering.

        Strategy: SEQRES = missing + present residues (sorted)
        Use REMARK 465 (missing_res_records) + structure to get actual PDB numbering
        """
        residue_map = {}

        if not seqres:
            # No SEQRES data for this chain (common with biological assemblies
            # where SEQRES records only cover asymmetric unit chains).
            # Fall back to building residue map from actual structure residues.
            for model in structure:
                if chain_id in model:
                    chain = model[chain_id]
                    for residue in chain:
                        if residue.id[0] == " ":  # Standard amino acid residue
                            aa_code = self.std_amino_acids.get(residue.resname, 'X')
                            residue_map[residue.id[1]] = aa_code
                    break
            return residue_map

        # Collect all residue numbers that should be in SEQRES
        all_seqres_residue_nums = set()

        # Add missing residues from REMARK 465 (if available)
        if missing_res_records and chain_id in missing_res_records:
            for missing_entry in missing_res_records[chain_id]:
                if isinstance(missing_entry, dict):
                    # Try both 'residue_number' and 'res_num' keys
                    res_num = missing_entry.get('residue_number') or missing_entry.get('res_num')
                    if res_num:
                        all_seqres_residue_nums.add(res_num)
                elif isinstance(missing_entry, int):
                    all_seqres_residue_nums.add(missing_entry)

        # Add present residues from structure (standard residues only, not HETATMs)
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                for residue in chain:
                    if residue.id[0] == " ":  # Standard amino acid residue
                        all_seqres_residue_nums.add(residue.id[1])
                break

        # Sort to get ordered residue numbers
        sorted_residue_nums = sorted(all_seqres_residue_nums)

        # Map SEQRES[i] to PDB residue number
        for i, res_name in enumerate(seqres):
            if i < len(sorted_residue_nums):
                pdb_num = sorted_residue_nums[i]
                aa_code = self.std_amino_acids.get(res_name, 'X')
                residue_map[pdb_num] = aa_code
            else:
                # SEQRES is longer than known residues - shouldn't happen but handle gracefully
                # Fall back to DBREF-based numbering
                pdb_start = dbref.get('pdb_start', dbref.get('seq_begin', 1))
                if isinstance(pdb_start, str):
                    pdb_start = int(pdb_start) if pdb_start.isdigit() else 1
                pdb_num = pdb_start + i
                aa_code = self.std_amino_acids.get(res_name, 'X')
                residue_map[pdb_num] = aa_code

        return residue_map
    
    def _build_chain_sequences(self, chain_id: str, complete_residues: Dict[int, str],
                              present: Set[int], filling: Set[int],
                              mutations: Dict[int, Tuple[str, str]],
                              structure: Structure) -> Tuple[str, str]:
        """Build template and target sequences for a chain"""
        template_seq = ""
        target_seq = ""

        # Process all expected residues
        for pdb_num in sorted(complete_residues.keys()):
            aa = complete_residues[pdb_num]
            
            if pdb_num in present:
                # Present residue
                if pdb_num in mutations:
                    # Mutation: original in template, mutated in target
                    from_aa, to_aa = mutations[pdb_num]
                    from_1 = self.std_amino_acids.get(from_aa, 'X')
                    to_1 = self.std_amino_acids.get(to_aa, 'X')
                    template_seq += from_1
                    target_seq += to_1
                else:
                    # No mutation
                    template_seq += aa
                    target_seq += aa
            
            elif pdb_num in filling:
                # Missing, to fill
                if pdb_num in mutations:
                    # Mutation on missing residue
                    _, to_aa = mutations[pdb_num]
                    to_1 = self.std_amino_acids.get(to_aa, 'X')
                    template_seq += "-"
                    target_seq += to_1
                else:
                    template_seq += "-"
                    target_seq += aa
            
            # else: skipped residue - don't include

        # CRITICAL: Dynamic HETATM counting from actual structure being repaired
        # Count actual HETATM groups in THIS structure (not from metadata)
        water_count = 0
        other_hetatm_count = 0
        
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                for residue in chain:
                    if residue.id[0] != " ":  # HETATM residue
                        if residue.resname in ['HOH', 'WAT']:
                            water_count += 1
                        else:
                            other_hetatm_count += 1
        
        # Add HETATM markers in correct order
        # Non-water HETATM first (ligands, ions, cofactors) using dots
        if other_hetatm_count > 0:
            hetatm_dots = "." * other_hetatm_count
            template_seq += hetatm_dots
            target_seq += hetatm_dots
        
        # Water molecules last using 'w' characters
        if water_count > 0:
            water_chars = "w" * water_count
            template_seq += water_chars
            target_seq += water_chars

        return template_seq, target_seq
    
    # ------------------------------------------------------------------
    # Quality assessment
    #
    # Design note: model quality is reported using MODELLER's own published,
    # calibrated assessment scores rather than an ad-hoc per-residue MOLPDF
    # heuristic. MOLPDF is the optimizer's internal objective function in
    # arbitrary units and is NOT comparable across systems, so it is shown for
    # provenance only. The interpreted metrics are:
    #   • Normalized DOPE z-score — Shen & Sali, Protein Sci. 2006, 15:2507.
    #       z <= -1 indicates a likely-correct (native-like) fold.
    #   • GA341 — John & Sali, Nucleic Acids Res. 2003, 31:3982 (range 0-1;
    #       >= 0.7 indicates a reliable fold).
    # For gap filling, these GLOBAL scores are dominated by the retained
    # experimental coordinates, so the decisive signal is the per-residue
    # assessment of the rebuilt region (_assess_built_region), which compares
    # each built residue to the structure's own experimental residues.
    # ------------------------------------------------------------------

    # Normalized-DOPE z-score above which a built residue is flagged as
    # energetically anomalous relative to its chain's experimental baseline.
    BUILT_RESIDUE_Z_THRESHOLD = 2.0
    # Minimum number of experimental residues needed to form a baseline.
    BUILT_RESIDUE_MIN_BASELINE = 8

    def _extract_assessment(self, mdl) -> Dict:
        """Read calibrated scores from the MODELLER AutoModel outputs.

        Returns a dict with float-or-None values for normalized_dope, ga341,
        dope and molpdf. No stdout parsing: values come straight from
        mdl.outputs (populated by assess_methods), so a MODELLER version
        change cannot silently break scoring.
        """
        result = {
            'normalized_dope': None,
            'ga341': None,
            'dope': None,
            'molpdf': None,
        }

        try:
            output = mdl.outputs[0]
        except (AttributeError, IndexError, TypeError, KeyError):
            return result

        def _as_float(value):
            # GA341 is returned as a list/tuple whose first element is the
            # composite score; everything else is a scalar.
            if isinstance(value, (tuple, list)):
                value = value[0] if value else None
            try:
                return float(value)
            except (TypeError, ValueError):
                return None

        result['normalized_dope'] = _as_float(output.get('Normalized DOPE score'))
        result['ga341'] = _as_float(output.get('GA341 score'))
        result['dope'] = _as_float(output.get('DOPE score'))
        result['molpdf'] = _as_float(output.get('molpdf'))
        return result

    @staticmethod
    def _interpret_normalized_dope(z: float) -> Tuple[str, str]:
        """Return (label, rich_color) for a normalized DOPE z-score."""
        if z <= -1.0:
            return "Native-like fold (z ≤ −1)", "green"
        if z <= 0.0:
            return "Borderline (−1 < z ≤ 0)", "yellow"
        return "Likely incorrect fold (z > 0)", "red"

    @staticmethod
    def _interpret_ga341(score: float) -> Tuple[str, str]:
        """Return (label, rich_color) for a GA341 score."""
        if score >= 0.7:
            return "Reliable fold (≥ 0.7)", "green"
        return "Low reliability (< 0.7)", "yellow"

    def _count_protein(self, structure: Structure) -> Tuple[int, int]:
        """Count standard protein residues and their atoms in a structure."""
        n_res = 0
        n_atoms = 0
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == " ":  # standard (non-HETATM) residue
                        n_res += 1
                        n_atoms += len(list(residue.get_atoms()))
        return n_res, n_atoms

    def _display_quality_assessment(self, assessment: Dict,
                                    repaired_structure: Structure,
                                    built_region: Optional[Dict] = None) -> None:
        """Display calibrated global scores plus the per-residue built-region
        assessment. No invented quality grade or MD-readiness verdict."""
        table = Table(title="MODELLER Model Assessment")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="yellow")
        table.add_column("Interpretation (literature-calibrated)", style="white")

        n_res, n_atoms = self._count_protein(repaired_structure)
        if n_res:
            table.add_row("Protein residues", str(n_res), "")
        if n_atoms:
            table.add_row("Protein atoms", str(n_atoms), "")

        nd = assessment.get('normalized_dope')
        if nd is not None:
            label, color = self._interpret_normalized_dope(nd)
            table.add_row(
                "Normalized DOPE",
                f"{nd:.2f}",
                f"[{color}]{label}[/{color}]  (Shen & Sali 2006)"
            )

        ga = assessment.get('ga341')
        if ga is not None:
            label, color = self._interpret_ga341(ga)
            table.add_row(
                "GA341",
                f"{ga:.3f}",
                f"[{color}]{label}[/{color}]  (John & Sali 2003)"
            )

        dope = assessment.get('dope')
        if dope is not None:
            table.add_row(
                "DOPE (raw)", f"{dope:.1f}",
                "[grey50]non-normalized — provenance only[/grey50]"
            )

        molpdf = assessment.get('molpdf')
        if molpdf is not None:
            table.add_row(
                "molpdf", f"{molpdf:.1f}",
                "[grey50]optimizer-internal — not comparable across systems[/grey50]"
            )

        self.console.print(table)

        if nd is None and ga is None:
            self.console.print(
                "[yellow]Calibrated assessment scores were not computed for "
                "this model.[/yellow]"
            )

        # The global scores above are dominated by retained experimental
        # coordinates when filling gaps; the rebuilt region is assessed
        # separately and is the decisive signal.
        self.console.print(
            "[grey50]Global scores are dominated by retained experimental "
            "coordinates; the rebuilt region is assessed below.[/grey50]"
        )
        self._display_built_region(built_region)

    def _assess_built_region(self, env, model_pdb: str,
                             built_residues: Set[Tuple[str, int]]) -> Dict:
        """Per-residue normalized-DOPE assessment of the residues MODELLER
        built, relative to the experimentally-resolved residues of the SAME
        chain (a within-structure z-score that self-calibrates per system).

        built_residues: (chain_id, res_num) in MODELLER output numbering.
        Returns a result dict; never raises (returns {'error': ...} on failure).
        """
        STANDARD_AA = frozenset({
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
            "TYR", "VAL", "MSE",
        })
        try:
            import io as _io
            import sys as _sys
            import statistics
            from modeller import Selection
            from modeller.scripts import complete_pdb

            # Resolve the model file independently of the current directory.
            model_pdb = os.path.abspath(model_pdb)
            env.io.atom_files_directory = [os.path.dirname(model_pdb) or '.']

            # Suppress MODELLER's banner/log chatter during scoring.
            old_out, old_err = _sys.stdout, _sys.stderr
            _sys.stdout = _sys.stderr = _io.StringIO()
            try:
                mdl = complete_pdb(env, os.path.basename(model_pdb))
                profile = Selection(mdl).get_dope_profile().get_normalized()
            finally:
                _sys.stdout, _sys.stderr = old_out, old_err

            # Per-residue normalized DOPE, aligned position-for-position with
            # mdl.residues (confirmed against the MODELLER EnergyProfile API).
            built = {(c, int(n)) for (c, n) in built_residues}

            # Index the re-parsed model so we can pick a matching strategy.
            model_keys = set()
            for residue in mdl.residues:
                try:
                    model_keys.add((residue.chain.name, int(residue.num)))
                except (ValueError, TypeError):
                    continue
            # Prefer exact (chain, num) matching, but fall back to num-only if
            # the chain keys don't line up — MODELLER chain renaming has been a
            # recurring source of key mismatches in this pipeline, and an empty
            # report would otherwise hide the rebuilt residues entirely.
            if built & model_keys:
                def _is_built(chain, num):
                    return (chain, num) in built
            else:
                built_nums = {n for _, n in built}

                def _is_built(chain, num):
                    return num in built_nums

            per_chain_baseline = defaultdict(list)
            rows = []  # (chain, num, name, energy, is_built)
            for residue, element in zip(mdl.residues, profile):
                if residue.name not in STANDARD_AA:
                    continue
                try:
                    num = int(residue.num)
                except (ValueError, TypeError):
                    continue
                chain = residue.chain.name
                energy = element.energy
                is_built = _is_built(chain, num)
                rows.append((chain, num, residue.name, energy, is_built))
                if not is_built:
                    per_chain_baseline[chain].append(energy)

            results = []
            n_flagged = 0
            for chain, num, name, energy, is_built in rows:
                if not is_built:
                    continue
                baseline = per_chain_baseline.get(chain, [])
                z = None
                note = ""
                if len(baseline) >= self.BUILT_RESIDUE_MIN_BASELINE:
                    mean = statistics.fmean(baseline)
                    sd = statistics.pstdev(baseline)
                    if sd > 0:
                        z = (energy - mean) / sd
                    else:
                        note = "zero-variance baseline"
                else:
                    note = "insufficient experimental baseline"
                flagged = z is not None and z > self.BUILT_RESIDUE_Z_THRESHOLD
                if flagged:
                    n_flagged += 1
                results.append({
                    'chain': chain, 'num': num, 'name': name,
                    'energy': energy, 'z': z, 'flagged': flagged, 'note': note,
                })

            results.sort(key=lambda r: (r['z'] is None, -(r['z'] or 0)))
            return {
                'residues': results,
                'n_built': len(results),
                'n_flagged': n_flagged,
                'threshold': self.BUILT_RESIDUE_Z_THRESHOLD,
            }
        except Exception as e:  # never break a repair over the report
            return {'error': str(e)}

    def _display_built_region(self, built_region: Optional[Dict]) -> None:
        """Render the per-residue built-region assessment."""
        if not built_region:
            return
        if built_region.get('error'):
            self.console.print(
                f"[yellow]Per-residue assessment of the rebuilt region was not "
                f"computed: {built_region['error']}[/yellow]"
            )
            return

        residues = built_region.get('residues', [])
        if not residues:
            return

        threshold = built_region.get('threshold', self.BUILT_RESIDUE_Z_THRESHOLD)
        table = Table(
            title="Rebuilt-Region Assessment (per-residue normalized DOPE)"
        )
        table.add_column("Residue", style="cyan")
        table.add_column("Norm. DOPE", style="yellow", justify="right")
        table.add_column("z vs chain", style="yellow", justify="right")
        table.add_column("", style="white")

        for r in residues:
            if r['z'] is None:
                z_str = "—"
                flag = f"[grey50]{r['note']}[/grey50]" if r['note'] else ""
            elif r['flagged']:
                z_str = f"[red]+{r['z']:.2f}[/red]"
                flag = "[red]⚠ elevated energy[/red]"
            else:
                z_str = f"[green]{r['z']:+.2f}[/green]"
                flag = "[green]✓ within range[/green]"
            table.add_row(
                f"{r['name']} {r['chain']}:{r['num']}",
                f"{r['energy']:.2f}",
                z_str,
                flag,
            )

        self.console.print(table)
        self.console.print(
            f"[grey50]Rule: each rebuilt residue's normalized DOPE energy is "
            f"compared to the same chain's experimentally-resolved residues "
            f"(within-structure z-score); flagged at z > {threshold:.1f}.[/grey50]"
        )
        n_flagged = built_region.get('n_flagged', 0)
        n_built = built_region.get('n_built', 0)
        if n_flagged:
            self.console.print(
                f"[red]{n_flagged}/{n_built} rebuilt residue(s) show elevated "
                f"energy and warrant inspection.[/red]"
            )
        else:
            self.console.print(
                f"[green]No elevated-energy residues among the {n_built} "
                f"rebuilt residue(s).[/green]"
            )

    def visualize_sequence_alignment(self, template_seq: str, target_seq: str, chain_id: str) -> None:
        """
        Visualize the sequence alignment showing matches, insertions, and gaps.
        
        Display format:
        Pos:        10        20        30
                |         |         |
        Template: MAETKVILGSGGSMATYF
                |||||||||  |||||||
        Model:    MAETKVILG--GSMATYF
        
        Legend: | = match, + = insertion, - = deletion
        """
        self.console.print(f"\n[bold]Chain {chain_id} Sequence Alignment:[/bold]")
        
        # Display in blocks of 60 characters
        BLOCK_SIZE = 60
        
        for i in range(0, len(template_seq), BLOCK_SIZE):
            end = min(i + BLOCK_SIZE, len(template_seq))
            
            # Position numbers line
            pos_line = " " * 10  # Initial padding for labels
            for j in range(i, end, 10):
                pos_label = str(j + 10)[:-1]  # Last digit of next position
                spaces = 9 - len(pos_label)
                pos_line += " " * spaces + pos_label + "0"
            self.console.print(pos_line)
            
            # Tick marks line
            tick_line = " " * 10  # Initial padding
            for j in range(i, end):
                tick_line += "|" if (j % 10) == 0 else " "
            self.console.print(tick_line)
            
            # Template sequence
            template_chunk = template_seq[i:end]
            self.console.print(f"Template: {template_chunk}")
            
            # Match line
            match_line = " " * 10  # Consistent padding with sequences
            for t, m in zip(template_chunk, target_seq[i:end]):
                if t == m and t != "-" and m != "-":
                    match_line += "|"  # Matching positions
                elif t == "-" and m != "-":
                    match_line += "+"  # Insertion in model
                elif t != "-" and m == "-":
                    match_line += "-"  # Deletion in model
                else:
                    match_line += " "  # Mismatch or both gaps
            self.console.print(match_line)
            
            # Model/Target sequence
            model_chunk = target_seq[i:end]
            self.console.print(f"Model:    {model_chunk}")
            
            self.console.print()  # Blank line between blocks
        
        # Legend
        self.console.print("[bold]Legend:[/bold]")
        self.console.print("  | = Match")
        self.console.print("  + = Insertion (residue added)")
        self.console.print("  - = Deletion (gap in model)")
        self.console.print()

    def visualize_all_chain_alignments(self, chain_sequences: Dict[str, Tuple[str, str]]) -> None:
        """Visualize alignments for all chains"""
        self.console.print("\n[bold cyan]═══ Sequence Alignments ═══[/bold cyan]")
        
        for chain_id in sorted(chain_sequences.keys()):
            template_seq, target_seq = chain_sequences[chain_id]
            self.visualize_sequence_alignment(template_seq, target_seq, chain_id)

    def display_sequence_summary(self, chain_sequences: Dict[str, Tuple[str, str]]) -> None:
        """Display summary of sequences for all chains"""
        from rich.table import Table
        
        table = Table(title="MODELLER Sequences Summary")
        table.add_column("Chain", style="cyan")
        table.add_column("Length", style="yellow")
        table.add_column("Standard AA", style="green")
        table.add_column("HETATM (.)", style="magenta")
        table.add_column("Water (w)", style="blue")
        table.add_column("Gaps (-)", style="red")
        
        for chain_id in sorted(chain_sequences.keys()):
            template_seq, target_seq = chain_sequences[chain_id]
            
            # Count different types
            aa_count = sum(1 for c in template_seq if c.isalpha() and c not in '.w-/')
            hetatm_count = template_seq.count('.')
            water_count = template_seq.count('w')
            gap_count = template_seq.count('-')
            
            table.add_row(
                chain_id,
                str(len(template_seq)),
                str(aa_count),
                str(hetatm_count),
                str(water_count),
                str(gap_count)
            )
        
        self.console.print()
        self.console.print(table)
        self.console.print()


# ============================================================================
# CAPPING IMPLEMENTATION
# ============================================================================

class CappingHandler:
    """Handles ACE/NME cap insertion"""
    
    def __init__(self, console: Console):
        self.console = console
    
    def add_caps(self, repaired_structure: Structure, plan: RepairPlan,
                mapper: ResidueMapper, output_file: str) -> Structure:
        """Add ACE/NME caps to repaired structure"""
        if not plan.has_caps:
            return repaired_structure
        
        self.console.print("\n[cyan]Adding ACE/NME caps...[/cyan]")
        
        # Save to temp file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            temp_file = tmp.name
        
        try:
            # Save structure
            io = PDBIO()
            io.set_structure(repaired_structure)
            io.save(temp_file)
            
            # Read lines
            with open(temp_file, 'r') as f:
                pdb_lines = f.readlines()
            
            # Build cap insertion plan using mapper
            cap_insertions = self._build_cap_insertion_plan(plan, mapper, repaired_structure)
            
            # Insert caps
            capped_lines = self._insert_caps_with_plan(pdb_lines, cap_insertions, mapper)
            
            # Renumber entire structure
            final_lines = self._renumber_structure(capped_lines, mapper)
            
            # Write final
            with open(output_file, 'w') as f:
                f.writelines(final_lines)
            
            # Parse final structure
            parser = PDBParser(QUIET=True)
            final_structure = parser.get_structure("final", output_file)
            
            self.console.print(f"[green]✓ Caps added successfully[/green]")
            return final_structure
            
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)
    
    def _build_cap_insertion_plan(self, plan: RepairPlan, mapper: ResidueMapper,
                                  repaired_structure: Structure) -> Dict[str, List[Dict]]:
        """
        Build detailed plan for where to insert caps.
        
        Returns: {chain_id: [{'type': 'ACE'|'NME', 'insert_at': modeller_resnum, 
                              'coords_ref': resnum_for_coords}, ...]}
        """
        insertions = defaultdict(list)

        for segment in plan.segments_to_cap:
            # segment.chain_id is the ORIGINAL chain id (the plan is built off the
            # pre-MODELLER structure). MODELLER renames chains (e.g. B,D → A,B),
            # so translate to the renamed id before indexing repaired_structure
            # and keying `insertions`. Without this, a renamed chain's caps are
            # either dropped (id absent) or mis-applied to a coincidentally
            # same-named chain — and `_insert_caps_with_plan` keys by the renamed
            # PDB id, so this is also what keeps the two halves consistent.
            chain_id = mapper.chain_mapping.get(segment.chain_id, segment.chain_id)

            if segment.terminal_type == 'N':
                # N-terminal: ACE before first EXISTING residue (after the gap)
                # repaired_structure already has MODELLER numbering (1-N)
                # Find the first existing residue in the MODELLER-numbered structure
                first_existing_modeller_num = None
                for model in repaired_structure:
                    if chain_id in model:
                        chain = model[chain_id]
                        residues = sorted([r for r in chain if r.id[0] == " "],
                                        key=lambda r: r.id[1])
                        if residues:
                            first_existing_modeller_num = residues[0].id[1]
                            break

                if first_existing_modeller_num is None:
                    self.console.print(
                        f"[yellow]Warning: No existing residues found in chain {chain_id} for N-terminal capping[/yellow]"
                    )
                    continue

                insertions[chain_id].append({
                    'type': 'ACE',
                    'insert_before': first_existing_modeller_num,
                    'ref_residue': first_existing_modeller_num,
                    'segment_start': segment.start_num,
                    'segment_end': segment.end_num
                })
                self.console.print(
                    f"[cyan]  Plan: ACE before {chain_id}:{first_existing_modeller_num} "
                    f"(MODELLER numbering, capping N-terminal gap {segment.start_num}-{segment.end_num})[/cyan]"
                )

            elif segment.terminal_type == 'C':
                # C-terminal: NME after last EXISTING protein residue (before the gap)
                # repaired_structure already has MODELLER numbering (1-N)
                # Find the last existing PROTEIN residue (excluding heteroatoms)
                last_protein_modeller_num = None
                for model in repaired_structure:
                    if chain_id in model:
                        chain = model[chain_id]
                        # Filter for standard amino acids only (hetflag == " ")
                        protein_residues = sorted([r for r in chain if r.id[0] == " "],
                                                key=lambda r: r.id[1])
                        if protein_residues:
                            last_protein_modeller_num = protein_residues[-1].id[1]
                            break

                if last_protein_modeller_num is None:
                    self.console.print(
                        f"[yellow]Warning: No protein residues found in chain {chain_id} for C-terminal capping[/yellow]"
                    )
                    continue

                # Use MODELLER number directly (no mapping needed)
                insertions[chain_id].append({
                    'type': 'NME',
                    'insert_after': last_protein_modeller_num,
                    'ref_residue': last_protein_modeller_num,
                    'segment_start': segment.start_num,
                    'segment_end': segment.end_num
                })
                self.console.print(
                    f"[cyan]  Plan: NME after {chain_id}:{last_protein_modeller_num} "
                    f"(MODELLER numbering, capping C-terminal gap {segment.start_num}-{segment.end_num})[/cyan]"
                )
            
            else:
                # Internal gap the user chose to CAP (not fill): break the chain
                # by placing NME after the residue preceding the gap and ACE
                # before the residue following it. segment.start_num/end_num are
                # ORIGINAL residue numbers, but repaired_structure is in MODELLER
                # numbering, so translate the two gap-bracketing residues through
                # the mapper. (The prior code referenced undefined
                # first_modeller/last_modeller and crashed with a NameError.)
                preceding_res = None
                following_res = None

                orig_to_modeller = mapper.residue_mappings.get(chain_id, {})
                if orig_to_modeller:
                    before = [o for o in orig_to_modeller if o < segment.start_num]
                    after = [o for o in orig_to_modeller if o > segment.end_num]
                    if before:
                        preceding_res = orig_to_modeller[max(before)]
                    if after:
                        following_res = orig_to_modeller[min(after)]
                else:
                    # Caps-only path (identity mapping, MODELLER did not run): the
                    # repaired structure keeps original numbering, so bracket the
                    # gap by scanning it directly.
                    for model in repaired_structure:
                        if chain_id in model:
                            residues = sorted(
                                [r for r in model[chain_id] if r.id[0] == " "],
                                key=lambda r: r.id[1])
                            for res in residues:
                                if res.id[1] < segment.start_num:
                                    preceding_res = res.id[1]
                                elif res.id[1] > segment.end_num and following_res is None:
                                    following_res = res.id[1]
                            break
                
                # NME after preceding residue
                if preceding_res:
                    insertions[chain_id].append({
                        'type': 'NME',
                        'insert_after': preceding_res,
                        'ref_residue': preceding_res,
                        'segment_start': segment.start_num,
                        'segment_end': segment.end_num
                    })
                    self.console.print(
                        f"[cyan]  Plan: NME after {chain_id}:{preceding_res} "
                        f"(internal gap {segment.start_num}-{segment.end_num})[/cyan]"
                    )
                
                # ACE before following residue
                if following_res:
                    insertions[chain_id].append({
                        'type': 'ACE',
                        'insert_before': following_res,
                        'ref_residue': following_res,
                        'segment_start': segment.start_num,
                        'segment_end': segment.end_num
                    })
                    self.console.print(
                        f"[cyan]  Plan: ACE before {chain_id}:{following_res} "
                        f"(internal gap {segment.start_num}-{segment.end_num})[/cyan]"
                    )
        
        return dict(insertions)
    
    def _insert_caps_with_plan(self, pdb_lines: List[str], 
                               cap_insertions: Dict[str, List[Dict]],
                               mapper: ResidueMapper) -> List[str]:
        """Insert caps according to detailed plan"""
        # Parse lines into structure
        lines_by_chain = defaultdict(lambda: defaultdict(list))
        header_lines = []
        footer_lines = []
        
        max_serial = 0
        
        for line in pdb_lines:
            if line.startswith(('ATOM', 'HETATM')):
                chain_id = line[21]
                res_num = int(line[22:26].strip())
                lines_by_chain[chain_id][res_num].append(line)
                
                try:
                    serial = int(line[6:11].strip())
                    max_serial = max(max_serial, serial)
                except:
                    pass
            elif line.startswith(('TER', 'END')):
                footer_lines.append(line)
            else:
                header_lines.append(line)
        
        # Build final structure with caps inserted
        final_lines = header_lines.copy()
        current_serial = max_serial + 1
        
        for chain_id in sorted(lines_by_chain.keys()):
            chain_residues = lines_by_chain[chain_id]
            sorted_resnums = sorted(chain_residues.keys())
            
            # Get cap plan for this chain
            chain_caps = cap_insertions.get(chain_id, [])

            # Residues that get a C-terminal NME bonded to them must lose their
            # OXT: the terminal carboxylate becomes a peptide bond to the cap, so
            # a lingering OXT over-coordinates the carbonyl carbon and breaks
            # tLEaP's mid-chain residue-template match.
            nme_after = {cap.get('insert_after') for cap in chain_caps
                         if cap.get('type') == 'NME' and cap.get('insert_after') is not None}

            # Process each residue, inserting caps at appropriate positions
            for res_num in sorted_resnums:
                # Check if we should insert caps BEFORE this residue
                for cap in chain_caps:
                    if cap.get('insert_before') == res_num:
                        # Insert ACE cap with unique residue number
                        # For N-terminal: use 0
                        # For internal: use res_num - 1 (the gap space before this residue)
                        ace_res_num = 0 if res_num == 1 else res_num - 1
                        ref_lines = chain_residues[res_num]
                        ace_lines = self._create_ace_cap(
                            chain_id, ace_res_num, current_serial, ref_lines[0]
                        )
                        final_lines.extend(ace_lines)
                        current_serial += len(ace_lines)

                        self.console.print(
                            f"[green]  Inserted ACE at {chain_id}:{ace_res_num} (before {res_num})[/green]"
                        )

                # Add the residue itself (drop OXT if an NME caps this residue).
                res_lines = chain_residues[res_num]
                if res_num in nme_after:
                    res_lines = [ln for ln in res_lines
                                 if ln[12:16].strip() != 'OXT']
                final_lines.extend(res_lines)

                # Check if we should insert caps AFTER this residue
                for cap in chain_caps:
                    if cap.get('insert_after') == res_num:
                        # Insert NME cap with unique residue number
                        # Use res_num + 1 (the gap space after this residue, or next sequential for C-terminal)
                        nme_res_num = res_num + 1
                        ref_lines = chain_residues[res_num]
                        nme_lines = self._create_nme_cap(
                            chain_id, nme_res_num, current_serial, ref_lines[-1]
                        )
                        final_lines.extend(nme_lines)
                        current_serial += len(nme_lines)

                        self.console.print(
                            f"[green]  Inserted NME at {chain_id}:{nme_res_num} (after {res_num})[/green]"
                        )
        
        final_lines.extend(footer_lines)
        return final_lines
    
    def _create_ace_cap(self, chain_id: str, res_num: int, serial: int,
                       ref_line: str) -> List[str]:
        """Generate ACE cap HETATM records using reference coordinates"""
        # Extract coordinates from reference line (should be N atom ideally)
        try:
            x = float(ref_line[30:38])
            y = float(ref_line[38:46])
            z = float(ref_line[46:54])
            ref_coords = np.array([x, y, z])
        except:
            # Fallback to default position
            ref_coords = np.array([0.0, 0.0, 0.0])

        # Calculate ACE atom positions relative to reference (the attachment N).
        # Use a trigonal-planar arrangement (~120deg between the carbonyl C's
        # three bonds: to N, O, CH3) rather than a single-axis layout. A
        # collinear layout (CH3, C, N all on one line) makes the C->CH3 and C->N
        # bond vectors antiparallel, so geometry tools that test planarity via a
        # cross product (e.g. PROPKA's SYBYL typing) divide by a zero-length
        # normal and crash. tLEaP rebuilds these coordinates regardless.
        c_coords = ref_coords - np.array([1.32, 0.0, 0.0])   # C->N bond along +x
        o_coords = c_coords + 1.23 * np.array([-0.5, 0.866, 0.0])
        ch3_coords = c_coords + 1.50 * np.array([-0.5, -0.866, 0.0])

        lines = []
        lines.append(
            f"HETATM{serial:5d}  CH3 ACE {chain_id}{res_num:4d}    "
            f"{ch3_coords[0]:8.3f}{ch3_coords[1]:8.3f}{ch3_coords[2]:8.3f}  "
            f"1.00 20.00           C  \n"
        )
        lines.append(
            f"HETATM{serial+1:5d}  C   ACE {chain_id}{res_num:4d}    "
            f"{c_coords[0]:8.3f}{c_coords[1]:8.3f}{c_coords[2]:8.3f}  "
            f"1.00 20.00           C  \n"
        )
        lines.append(
            f"HETATM{serial+2:5d}  O   ACE {chain_id}{res_num:4d}    "
            f"{o_coords[0]:8.3f}{o_coords[1]:8.3f}{o_coords[2]:8.3f}  "
            f"1.00 20.00           O  \n"
        )

        return lines
    
    def _create_nme_cap(self, chain_id: str, res_num: int, serial: int,
                       ref_line: str) -> List[str]:
        """Generate NME cap HETATM records using reference coordinates.

        Only the amide N and H are emitted; the methyl carbon and its three
        hydrogens are intentionally omitted so tLEaP builds them from the
        loaded force field's NME template. Amber force fields disagree on the
        methyl carbon's atom name (some call it CH3, some C), so an explicitly
        placed CH3 would fail to match a template expecting C (and vice versa).
        A missing atom is simply rebuilt under the template's own name, which
        sidesteps the mismatch for whichever protein FF (ff14SB/ff19SB/...) is
        sourced at tLEaP time.
        """
        # Extract coordinates from reference line (should be C atom ideally)
        try:
            x = float(ref_line[30:38])
            y = float(ref_line[38:46])
            z = float(ref_line[46:54])
            ref_coords = np.array([x, y, z])
        except:
            # Fallback to default position
            ref_coords = np.array([10.0, 0.0, 0.0])

        # Place N (bonded to the preceding carbonyl C) and the amide H. Keep H
        # off the N->C attachment axis: a non-collinear second neighbour gives N
        # a well-defined plane, so planarity tests (PROPKA SYBYL typing) don't
        # divide by a zero-length cross product and raise ZeroDivisionError.
        # tLEaP rebuilds these coordinates, plus the omitted methyl group.
        n_coords = ref_coords + np.array([1.32, 0.0, 0.0])   # N->C bond along -x
        h_coords = n_coords + 1.00 * np.array([0.5, 0.866, 0.0])

        lines = []
        lines.append(
            f"HETATM{serial:5d}  N   NME {chain_id}{res_num:4d}    "
            f"{n_coords[0]:8.3f}{n_coords[1]:8.3f}{n_coords[2]:8.3f}  "
            f"1.00 20.00           N  \n"
        )
        lines.append(
            f"HETATM{serial+1:5d}  H   NME {chain_id}{res_num:4d}    "
            f"{h_coords[0]:8.3f}{h_coords[1]:8.3f}{h_coords[2]:8.3f}  "
            f"1.00 20.00           H  \n"
        )

        return lines
    
    def _renumber_structure(self, pdb_lines: List[str], 
                           mapper: ResidueMapper) -> List[str]:
        """Renumber entire structure after cap insertion"""
        renumbered = []
        header = []
        footer = []
        atoms = []
        
        for line in pdb_lines:
            if line.startswith(('ATOM', 'HETATM')):
                atoms.append(line)
            elif line.startswith(('TER', 'END')):
                footer.append(line)
            else:
                header.append(line)
        
        renumbered.extend(header)
        
        # Renumber atoms by chain
        atoms_by_chain = defaultdict(list)
        for line in atoms:
            chain_id = line[21]
            atoms_by_chain[chain_id].append(line)
        
        atom_serial = 1
        for chain_id in sorted(atoms_by_chain.keys()):
            chain_atoms = atoms_by_chain[chain_id]
            
            # Group by residue. Identify a residue by (number, resname,
            # insertion code) — NOT number alone: cap insertion can place an
            # ACE/NME at the same number as an adjacent hetero residue (e.g. an
            # NME and a cofactor both ending up at 622), and grouping by number
            # alone would merge two distinct residues into one slot, leaving
            # duplicate residue numbers in the output. Any resname/icode change
            # therefore starts a new residue.
            residue_groups = []
            current_res = []
            current_key = None

            for line in chain_atoms:
                res_num = int(line[22:26].strip())
                res_key = (res_num, line[17:20], line[26:27])
                if res_key != current_key:
                    if current_res:
                        residue_groups.append((current_key[0], current_res))
                    current_res = [line]
                    current_key = res_key
                else:
                    current_res.append(line)

            if current_res:
                residue_groups.append((current_key[0], current_res))
            
            # Renumber consecutively. Record the authoritative
            # MODELLER-number → final-number map per (renamed) chain as we go.
            # This is the single source of truth for final numbering used by
            # ResidueMapper.get_final_identity during redox-site sync: it is read
            # straight off the file we are writing, so it cannot drift from disk
            # the way a separately-predicted capping map can (the predicted map
            # assumed MODELLER's global cross-chain numbering, but this restarts
            # each chain at 1 — the source of the chain-B mis-mapping).
            new_res_num = 1
            for old_num, res_lines in residue_groups:
                mapper.final_mappings.setdefault(chain_id, {})[old_num] = new_res_num
                for line in res_lines:
                    new_line = f"{line[:6]}{atom_serial:5d}{line[11:22]}{new_res_num:4d}{line[26:]}"
                    renumbered.append(new_line)
                    atom_serial += 1
                new_res_num += 1
        
        renumbered.extend(footer)
        return renumbered


# ============================================================================
# REDOX SITE SYNCHRONIZATION
# ============================================================================

class RedoxSiteSync:
    """Synchronizes RedoxSite objects with repaired structure"""
    
    def __init__(self, console: Console):
        self.console = console
    
    def synchronize_sites(self, 
                         redox_sites: List[Any],
                         final_structure: Structure,
                         mapper: ResidueMapper) -> Dict[str, Any]:
        """
        Update all RedoxSite objects to match final structure.
        
        Returns summary of updates performed.
        """
        if not redox_sites:
            self.console.print("[yellow]No RedoxSite objects to synchronize[/yellow]")
            return {'sites_updated': 0}

        self.console.print("\n[bold cyan]Synchronizing RedoxSite Objects with Repaired Structure[/bold cyan]")
        self.console.print("[grey50]After structure repair, redox site coordinates and residue numbering need to be updated[/grey50]")
        self.console.print("[grey50]to match the repaired structure while preserving site topology and bonds.[/grey50]\n")
        
        summary = {
            'sites_updated': 0,
            'atoms_updated': 0,
            'bonds_updated': 0,
            'centers_updated': 0,
            'details': []
        }
        
        for site in redox_sites:
            site_summary = self._update_site(site, final_structure, mapper)

            if site_summary['atoms_updated'] > 0:
                summary['sites_updated'] += 1
                summary['atoms_updated'] += site_summary['atoms_updated']
                summary['bonds_updated'] += site_summary['bonds_updated']
                summary['centers_updated'] += site_summary['centers_updated']
                summary['details'].append({
                    'site_id': site.site_id,
                    'site_type': getattr(site, 'site_type', 'Unknown'),
                    **site_summary
                })
            else:
                # Debug: Site had no atoms updated
                if site.atoms:
                    first_atom = site.atoms[0]
                    old_identity = AtomIdentity(
                        residue=ResidueIdentity(
                            first_atom.chain, first_atom.resid, first_atom.resname,
                            getattr(first_atom, 'insertion_code', ' ')
                        ),
                        atom_name=first_atom.atom_name,
                        element=first_atom.element
                    )
                    new_identity = mapper.get_final_atom_identity(old_identity)
                    self.console.print(f"[yellow]⚠ Site {site.site_id}: No atoms updated.[/yellow]")
                    self.console.print(f"  Old: {first_atom.chain}:{first_atom.resid}:{first_atom.atom_name}")
                    self.console.print(f"  New: {new_identity.residue.chain_id}:{new_identity.residue.res_num}:{new_identity.atom_name}")
        
        self._display_summary(summary)
        return summary
    
    def _update_site(self, site: Any, final_structure: Structure,
                    mapper: ResidueMapper) -> Dict[str, int]:
        """Update a single RedoxSite object"""
        atoms_updated = 0
        bonds_updated = 0
        centers_updated = 0

        # Display what we're doing
        site_type = getattr(site, 'site_type', 'unknown')
        self.console.print(f"\n[bold cyan]→ Synchronizing {site.site_id}[/bold cyan] [grey50]({site_type})[/grey50]")

        # Show initial composition
        num_atoms = len(site.atoms)
        num_bonds = len(site.bonds) if hasattr(site, 'bonds') else 0
        num_centers = len(site.centers) if hasattr(site, 'centers') else 0

        self.console.print(f"  Site contains: {num_atoms} atoms, {num_bonds} bonds, {num_centers} center(s)")

        # Capture pre-sync state for before→after display.
        # Key by resname only and carry (chain, resid) pairs: MODELLER may
        # rename the chain (e.g. C → B), so keying on (chain, resname) would
        # split a single residue's before/after across two half-empty lines.
        pre_residues = {}
        for atom in site.atoms:
            pre_residues.setdefault(atom.resname, set()).add((atom.chain, atom.resid))
        pre_center_coords = {
            (c.chain, c.resid, c.resname): c.coords
            for c in site.centers
        } if hasattr(site, 'centers') else {}

        # Build old coordinate -> identity mapping before updates
        old_coords_to_identity = {}
        identity_changes = {}

        for atom in site.atoms:
            old_identity = AtomIdentity(
                residue=ResidueIdentity(
                    atom.chain, atom.resid, atom.resname,
                    getattr(atom, 'insertion_code', ' ')
                ),
                atom_name=atom.atom_name,
                element=atom.element
            )
            old_coords_to_identity[atom.coords] = old_identity
        
        # Update atoms
        atoms_not_found = []
        for atom in site.atoms:
            old_identity = AtomIdentity(
                residue=ResidueIdentity(
                    atom.chain, atom.resid, atom.resname,
                    getattr(atom, 'insertion_code', ' ')
                ),
                atom_name=atom.atom_name,
                element=atom.element
            )

            # Get new identity
            new_identity = mapper.get_final_atom_identity(old_identity)
            identity_changes[old_identity] = new_identity

            # Find atom in final structure
            new_coords = self._find_atom_coords(final_structure, new_identity)

            if new_coords:
                old_coords = atom.coords
                atom.coords = new_coords
                atom.chain = new_identity.residue.chain_id
                atom.resid = new_identity.residue.res_num
                atom.resname = new_identity.residue.res_name

                # Update coord_to_pdb mapping
                if hasattr(site, 'coord_to_pdb') and old_coords in site.coord_to_pdb:
                    site.coord_to_pdb[new_coords] = site.coord_to_pdb.pop(old_coords)
                    site.coord_to_pdb[new_coords].update({
                        'chain': new_identity.residue.chain_id,
                        'resid': new_identity.residue.res_num,
                        'resname': new_identity.residue.res_name
                    })

                atoms_updated += 1
            else:
                atoms_not_found.append(f"{old_identity.residue.res_name} {old_identity.residue.chain_id}:{old_identity.residue.res_num} {old_identity.atom_name}")

        # Report atoms not found (indicates potential issues)
        if atoms_not_found:
            self.console.print(f"  [yellow]⚠ {len(atoms_not_found)} atom(s) not found in repaired structure[/yellow]")
            self.console.print(f"    [grey50]This may occur if atoms were removed during repair[/grey50]")
            if len(atoms_not_found) <= 3:
                for atom_desc in atoms_not_found:
                    self.console.print(f"    • {atom_desc}")
            else:
                for atom_desc in atoms_not_found[:3]:
                    self.console.print(f"    • {atom_desc}")
                self.console.print(f"    [grey50]... and {len(atoms_not_found) - 3} more[/grey50]")
        
        # Update bonds
        if hasattr(site, 'bonds'):
            for bond in site.bonds:
                # Capture the pre-repair coordinates as position anchors before
                # we overwrite them: repair may have relabeled a coordinating
                # atom (OE1<->OE2), so we re-resolve each endpoint to the atom
                # physically nearest its old position, not the one that kept its
                # name. coord_to_pdb (updated in the atom loop above) already
                # carries the repaired structure's correct name at each position.
                old1 = bond.atom1_coords
                old2 = bond.atom2_coords
                atom1_identity = old_coords_to_identity.get(old1)
                atom2_identity = old_coords_to_identity.get(old2)

                if atom1_identity and atom2_identity:
                    new_atom1 = identity_changes.get(atom1_identity, atom1_identity)
                    new_atom2 = identity_changes.get(atom2_identity, atom2_identity)

                    new_coords1, _n1 = self._find_bonded_atom_coords(
                        final_structure, new_atom1.residue, old1,
                        new_atom1.element, new_atom1.atom_name)
                    new_coords2, _n2 = self._find_bonded_atom_coords(
                        final_structure, new_atom2.residue, old2,
                        new_atom2.element, new_atom2.atom_name)

                    if new_coords1 and new_coords2:
                        bond.atom1_coords = new_coords1
                        bond.atom2_coords = new_coords2

                        # Recalculate distance
                        bond.distance = float(np.linalg.norm(
                            np.array(new_coords1) - np.array(new_coords2)
                        ))
                        bonds_updated += 1
        
        # Update centers (redox-active atoms or residue centroids)
        if hasattr(site, 'centers'):
            for center in site.centers:
                # Map the residue identity first
                old_res_identity = ResidueIdentity(
                    center.chain, center.resid,
                    getattr(center, 'resname', 'UNK'),
                    getattr(center, 'insertion_code', ' ')
                )
                new_res_identity = mapper.get_final_identity(old_res_identity)

                # If atom_name is None, center is the centroid of the entire residue
                new_coords = None
                if center.atom_name is None:
                    # Find the residue and calculate centroid
                    for model in final_structure:
                        if new_res_identity.chain_id in model:
                            chain = model[new_res_identity.chain_id]
                            for residue in chain:
                                if (residue.id[1] == new_res_identity.res_num and
                                    residue.resname == new_res_identity.res_name and
                                    residue.id[2] == new_res_identity.insertion_code):
                                    # Calculate centroid
                                    coords = [atom.coord for atom in residue]
                                    if coords:
                                        centroid = tuple(round(x, 3) for x in np.mean(coords, axis=0))
                                        new_coords = centroid
                                    break
                            break

                else:
                    # Specific atom - use atom identity mapping
                    old_identity = AtomIdentity(
                        residue=old_res_identity,
                        atom_name=center.atom_name,
                        element=center.element
                    )
                    new_atom_identity = mapper.get_final_atom_identity(old_identity)
                    new_coords = self._find_atom_coords(final_structure, new_atom_identity)
                    new_res_identity = new_atom_identity.residue

                    # Update element if it was None (for specific atom centers only)
                    if new_coords and center.element is None and new_atom_identity.element:
                        center.element = new_atom_identity.element

                if new_coords:
                    old_coords = center.coords
                    center.coords = new_coords
                    center.chain = new_res_identity.chain_id
                    center.resid = new_res_identity.res_num
                    center.resname = new_res_identity.res_name

                    # Update coord_to_pdb
                    if hasattr(site, 'coord_to_pdb') and old_coords in site.coord_to_pdb:
                        site.coord_to_pdb[new_coords] = site.coord_to_pdb.pop(old_coords)
                        site.coord_to_pdb[new_coords].update({
                            'chain': new_res_identity.chain_id,
                            'resid': new_res_identity.res_num,
                            'resname': new_res_identity.res_name
                        })

                    centers_updated += 1

        # Rebuild residue_groups with updated chain/resid info
        if hasattr(site, 'residue_groups'):
            site.residue_groups = {}
            for atom in site.atoms:
                residue_key = (atom.chain, atom.resid, getattr(atom, 'insertion_code', ' '))
                if residue_key not in site.residue_groups:
                    site.residue_groups[residue_key] = []
                site.residue_groups[residue_key].append(atom.coords)

        # Display synchronization results
        if atoms_updated == num_atoms and bonds_updated == num_bonds and centers_updated == num_centers:
            self.console.print(f"  [green]✓ Complete: {atoms_updated} atoms, {bonds_updated} bonds, {centers_updated} centers synchronized[/green]")
        else:
            parts = []
            if atoms_updated > 0:
                parts.append(f"{atoms_updated}/{num_atoms} atoms")
            if bonds_updated > 0:
                parts.append(f"{bonds_updated}/{num_bonds} bonds")
            if centers_updated > 0:
                parts.append(f"{centers_updated}/{num_centers} centers")
            if parts:
                self.console.print(f"  [green]✓ Updated: {', '.join(parts)}[/green]")
            else:
                self.console.print(f"  [yellow]⚠ No updates performed[/yellow]")

        # ── Before→After residue mapping ──
        post_residues = {}
        for atom in site.atoms:
            post_residues.setdefault(atom.resname, set()).add((atom.chain, atom.resid))

        # Show residue ID changes (group by resname; chain may have been renamed)
        for resname in sorted(set(pre_residues) | set(post_residues)):
            old_pairs = pre_residues.get(resname, set())
            new_pairs = post_residues.get(resname, set())
            if old_pairs != new_pairs:
                self.console.print(
                    f"  {resname}: {self._format_chain_residues(old_pairs)} "
                    f"→ {self._format_chain_residues(new_pairs)}"
                )
            else:
                self.console.print(
                    f"  {resname}: {self._format_chain_residues(old_pairs)} (unchanged)"
                )

        # Show bond integrity
        if hasattr(site, 'bonds'):
            coord_keys = set(site.coord_to_pdb.keys()) if hasattr(site, 'coord_to_pdb') else set()
            n_coord = sum(1 for b in site.bonds if b.chemical_type == "coordinate")
            n_cov = len(site.bonds) - n_coord
            broken = []
            for i, bond in enumerate(site.bonds):
                if bond.atom1_coords not in coord_keys or bond.atom2_coords not in coord_keys:
                    a1 = site.coord_to_pdb.get(bond.atom1_coords, {})
                    a2 = site.coord_to_pdb.get(bond.atom2_coords, {})
                    broken.append(f"bond {i} ({bond.chemical_type}): "
                                  f"{a1.get('atom_name','?')}→{a2.get('atom_name','?')}")
            if broken:
                self.console.print(f"  [red]Bonds: {len(broken)} BROKEN coord_to_pdb lookups:[/red]")
                for b in broken:
                    self.console.print(f"    [red]{b}[/red]")
            else:
                self.console.print(f"  Bonds: {n_coord} coordinate + {n_cov} covalent, all OK")

        return {
            'atoms_updated': atoms_updated,
            'bonds_updated': bonds_updated,
            'centers_updated': centers_updated
        }

    @staticmethod
    def _format_chain_residues(pairs: Set[Tuple[str, int]]) -> str:
        """Format a set of (chain, resid) pairs as 'A:1,2,3 B:42,44' for display.

        Grouping per chain keeps the before→after diff readable even when
        MODELLER renames a chain (C → B), so a residue's old and new locations
        appear on the same line rather than as two half-empty entries.
        """
        if not pairs:
            return "(none)"
        by_chain: Dict[str, List[int]] = {}
        for chain, resid in pairs:
            by_chain.setdefault(chain, []).append(resid)
        return " ".join(
            f"{chain}:{','.join(str(r) for r in sorted(ids))}"
            for chain, ids in sorted(by_chain.items())
        )

    def _find_atom_coords(self, structure: Structure,
                         atom_identity: AtomIdentity) -> Optional[Tuple[float, float, float]]:
        """Find atom coordinates in structure"""
        res_id = atom_identity.residue

        for model in structure:
            if res_id.chain_id in model:
                chain = model[res_id.chain_id]

                # Scan ALL residues carrying this number for the atom. Cap
                # insertion can leave two residues sharing a number (e.g. an NME
                # and a following ion both at the same residue number before
                # renumbering merges their group), so bailing on the first match
                # would miss an atom that lives on the second residue.
                for residue in chain:
                    if residue.id[1] == res_id.res_num and atom_identity.atom_name in residue:
                        return tuple(residue[atom_identity.atom_name].get_coord())

        return None

    # Max distance (Å) between a bond's pre-repair coordinate and a candidate
    # atom for the two to be considered the same physical atom. Structure repair
    # nudges a coordinating atom by ~0.1 Å; a symmetry-equivalent partner (the
    # other carboxylate O, the other imidazole N) sits ~2 Å away, so this cleanly
    # separates "same atom, renamed" from "wrong atom".
    _RESNAP_TOL = 1.0

    def _find_bonded_atom_coords(self, structure: Structure,
                                 res_id: "ResidueIdentity",
                                 old_coords: Tuple[float, float, float],
                                 element: Optional[str],
                                 fallback_name: str
                                 ) -> Tuple[Optional[Tuple[float, float, float]], Optional[str]]:
        """Resolve a bonded atom by NEAREST position within its remapped residue.

        Structure repair (MODELLER) can relabel symmetry-equivalent atoms — the
        classic case is a Glu/Asp carboxylate whose OE1/OE2 (OD1/OD2) names swap
        when the sidechain is rebuilt. Carrying the pre-repair atom *name*
        forward would then bind the metal to the wrong (renamed) oxygen and give
        it the generic type. A bond's real intent is the atom at a *position*, so
        resolve to the atom in the mapped residue closest to the pre-repair
        coordinate (restricted to the same element), which tracks the physical
        atom through any renaming. Returns ``(coords, actual_name)``; falls back
        to the name-based lookup if nothing suitable is close enough.
        """
        target_elem = (element or "").strip().upper()
        best_atom = None
        best_d = None
        anchor = np.array(old_coords)
        for model in structure:
            if res_id.chain_id not in model:
                continue
            chain = model[res_id.chain_id]
            for residue in chain:
                if residue.id[1] != res_id.res_num:
                    continue
                for atom in residue:
                    a_elem = (getattr(atom, "element", "") or "").strip().upper()
                    if target_elem and a_elem and a_elem != target_elem:
                        continue
                    d = float(np.linalg.norm(np.array(atom.get_coord()) - anchor))
                    if best_d is None or d < best_d:
                        best_d = d
                        best_atom = atom
        if best_atom is not None and best_d is not None and best_d <= self._RESNAP_TOL:
            return tuple(best_atom.get_coord()), best_atom.name.strip()

        # Nothing close (large rebuild, or element mismatch): keep the old name.
        coords = self._find_atom_coords(
            structure,
            AtomIdentity(residue=res_id, atom_name=fallback_name, element=element),
        )
        return (coords, fallback_name) if coords else (None, None)
    
    def _display_summary(self, summary: Dict) -> None:
        """Display synchronization summary"""
        if summary['sites_updated'] == 0:
            self.console.print("[yellow]No RedoxSite updates needed[/yellow]")
            return
        
        self.console.print(
            f"\n[green]✓ Synchronization Complete[/green]\n"
            f"  Updated {summary['sites_updated']} site(s): "
            f"{summary['atoms_updated']} atoms, {summary['bonds_updated']} bonds, "
            f"{summary['centers_updated']} centers"
        )

        if summary['details']:
            from rich.table import Table
            table = Table(title="Site Synchronization Summary", show_header=True, header_style="bold")
            table.add_column("Site ID", style="cyan", no_wrap=True)
            table.add_column("Type", style="yellow")
            table.add_column("Atoms", style="green", justify="right")
            table.add_column("Bonds", style="magenta", justify="right")
            table.add_column("Centers", style="blue", justify="right")
            
            for detail in summary['details']:
                table.add_row(
                    detail['site_id'],
                    detail['site_type'],
                    str(detail['atoms_updated']),
                    str(detail['bonds_updated']),
                    str(detail['centers_updated'])
                )
            
            self.console.print(table)


# ============================================================================
# MAIN MODULE CLASS
# ============================================================================
# NON-STANDARD MUTATION APPLICATOR
# ============================================================================

class NonStandardMutationApplicator:
    """Applies non-standard residue mutations (atom removal + residue renaming)"""

    def __init__(self, console: Console):
        self.console = console

    def apply_nonstandard_mutations(self,
                                      structure: Structure,
                                      nonstandard_mutations: List[Tuple],
                                      mapper: ResidueMapper,
                                      redox_sites: Optional[List[Any]] = None) -> Tuple[Structure, List[str]]:
        """
        Apply non-standard mutations to structure.

        Args:
            structure: Bio.PDB.Structure object to modify
            nonstandard_mutations: List of (chain_id, res_num, from_aa, to_aa, atoms_to_keep)
            mapper: ResidueMapper for tracking residue renumbering
            redox_sites: Optional list of RedoxSite objects to update

        Returns:
            (modified_structure, warnings_list)
        """
        if not nonstandard_mutations:
            return structure, []

        self.console.print(f"\n[bold cyan]Applying {len(nonstandard_mutations)} non-standard mutation(s)...[/bold cyan]")

        warnings = []
        mutations_applied = []

        for chain_id, res_num, from_aa, to_aa, atoms_to_keep in nonstandard_mutations:
            self.console.print(f"\n[yellow]→ Mutating {chain_id}:{res_num} {from_aa} → {to_aa}[/yellow]")

            # Map the original (chain, residue) through the full transformation
            # chain: MODELLER chain rename (e.g. C → B when chain A precedes it) +
            # residue renumbering + cap insertions.  ResidueMapper keys every map
            # by the *new* chain id, so we must resolve the chain BEFORE indexing
            # the residue maps — looking up by the raw original chain id misses the
            # rename and silently returns the unmapped number against a chain that
            # no longer exists in the repaired structure.
            original_identity = ResidueIdentity(
                chain_id=chain_id, res_num=res_num, res_name=from_aa
            )
            final_identity = mapper.get_final_identity(original_identity)
            mapped_chain = final_identity.chain_id
            mapped_res_num = final_identity.res_num

            if (mapped_chain, mapped_res_num) != (chain_id, res_num):
                self.console.print(
                    f"  [grey50]Mapped residue: {chain_id}:{res_num} → {mapped_chain}:{mapped_res_num}[/grey50]"
                )

            # Find the residue in the repaired structure using the mapped identity
            residue = self._find_residue(structure, mapped_chain, mapped_res_num)

            if not residue:
                warning = (
                    f"Residue {chain_id}:{res_num} "
                    f"(mapped to {mapped_chain}:{mapped_res_num}) not found in structure"
                )
                self.console.print(f"  [red]✗ {warning}[/red]")
                warnings.append(warning)
                continue

            # Verify it matches the expected source residue
            if residue.resname != from_aa:
                warning = f"Residue {chain_id}:{mapped_res_num} is {residue.resname}, expected {from_aa}"
                self.console.print(f"  [red]✗ {warning}[/red]")
                warnings.append(warning)
                continue

            # Apply mutation
            try:
                self._apply_single_mutation(residue, to_aa, atoms_to_keep)
                mutations_applied.append((chain_id, res_num, from_aa, to_aa))
                self.console.print(f"  [green]✓ Mutation applied successfully[/green]")
                self.console.print(f"    Kept {len(atoms_to_keep)} atoms: {', '.join(atoms_to_keep)}")
            except Exception as e:
                warning = f"Failed to apply mutation {chain_id}:{res_num}: {str(e)}"
                self.console.print(f"  [red]✗ {warning}[/red]")
                warnings.append(warning)

        # Update RedoxSite objects if any mutations affected redox sites
        if redox_sites and mutations_applied:
            self._update_redox_sites_for_mutations(redox_sites, mutations_applied, structure, mapper)

        self.console.print(f"\n[green]✓ Applied {len(mutations_applied)}/{len(nonstandard_mutations)} non-standard mutations[/green]")

        return structure, warnings

    def _find_residue(self, structure: Structure, chain_id: str, res_num: int) -> Optional:
        """Find a specific residue in the structure."""
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                for residue in chain:
                    if residue.id[1] == res_num and residue.id[0] == " ":
                        return residue
        return None

    def _apply_single_mutation(self, residue, new_resname: str, atoms_to_keep: List[str]) -> None:
        """
        Apply a single non-standard mutation to a residue.

        Args:
            residue: Bio.PDB.Residue object
            new_resname: New residue name (3-letter code)
            atoms_to_keep: List of atom names to keep
        """
        # Get all current atoms
        all_atoms = list(residue.get_atoms())
        atoms_to_remove = [atom for atom in all_atoms if atom.name not in atoms_to_keep]

        # Remove unwanted atoms
        for atom in atoms_to_remove:
            residue.detach_child(atom.id)

        # Rename the residue
        residue.resname = new_resname

        self.console.print(f"    Removed {len(atoms_to_remove)} atoms")
        self.console.print(f"    Renamed residue to {new_resname}")

    def _update_redox_sites_for_mutations(self,
                                           redox_sites: List[Any],
                                           mutations_applied: List[Tuple],
                                           structure: Structure,
                                           mapper: ResidueMapper) -> None:
        """
        Update RedoxSite objects if mutations affected redox site residues.

        Args:
            redox_sites: List of RedoxSite objects
            mutations_applied: List of (chain_id, res_num, from_aa, to_aa)
            structure: Updated structure
            mapper: ResidueMapper for coordinate tracking
        """
        mutation_positions = {(chain_id, res_num) for chain_id, res_num, _, _ in mutations_applied}

        for site in redox_sites:
            # Check if any site atoms are in mutated residues
            affected_atoms = []
            for atom in site.atoms:
                if (atom.chain, atom.resid) in mutation_positions:
                    affected_atoms.append(atom)

            if affected_atoms:
                self.console.print(f"\n[yellow]⚠ RedoxSite {site.site_id} contains mutated residue(s)[/yellow]")
                self.console.print(f"  Affected atoms: {len(affected_atoms)}")

                # Update atom information (resname changed, some atoms may be missing)
                for atom in site.atoms:
                    if (atom.chain, atom.resid) in mutation_positions:
                        # This runs before synchronize_sites(), so the atom still
                        # carries its ORIGINAL chain/resid — map it through the
                        # MODELLER rename + renumbering to locate it in the
                        # repaired structure (e.g. C:108 → B:<modeller_num>).
                        mapped = mapper.get_final_identity(
                            ResidueIdentity(
                                atom.chain, atom.resid, atom.resname,
                                getattr(atom, 'insertion_code', ' ')
                            )
                        )
                        residue = self._find_residue(
                            structure, mapped.chain_id, mapped.res_num
                        )
                        if residue:
                            # Update resname
                            old_resname = atom.resname
                            atom.resname = residue.resname

                            # Check if atom still exists
                            atom_exists = any(a.name == atom.atom_name for a in residue.get_atoms())
                            if not atom_exists:
                                self.console.print(f"    [red]Warning: Atom {atom.atom_name} no longer exists in {atom.chain}:{atom.resid}[/red]")
                            else:
                                # Update coordinates
                                for res_atom in residue.get_atoms():
                                    if res_atom.name == atom.atom_name:
                                        atom.coords = tuple(res_atom.coord)
                                        break

                self.console.print(f"  [green]✓ RedoxSite {site.site_id} updated[/green]")


# ============================================================================
# MODULE INTERFACE
# ============================================================================

from proprep.utils.module_registry import ProcessingModule, register_module

@register_module
class StructureCompletenessModule(ProcessingModule):
    """
    Complete structure completeness module with detection and repair.
    
    Features:
    - 8 detection methods (4 missing residues, 2 missing atoms, 2 alternate locations)
    - Interactive repair planning with user decisions
    - MODELLER integration for filling residues and applying mutations
    - ACE/NME capping for terminal gaps
    - Complete residue mapping and chain renumbering tracking
    - RedoxSite synchronization
    - Both interactive and batch modes
    """
    
    NAME = "Structure Fixer"
    DESCRIPTION = "Detect and repair missing atoms, residues, and alternate locations"
    VERSION = "2.0.0"
    CATEGORY = "analysis"
    PRIORITY = 3
    
    def initialize(self):
        """Initialize module"""
        self.console = Console()
        self.analyzer = None
        self.results = None
    
    # ========================================================================
    # MODULE INTERFACE
    # ========================================================================
    
    def get_workspace_requirements(self) -> List[str]:
        """Required workspace items - needs experimental structure (not AlphaFold)"""
        return [
            "rcsb_pdb_file | local_pdb_file"
        ]
    
    def get_workspace_outputs(self) -> List[str]:
        """Workspace items this module produces"""
        return [
            "completeness_results",
            "repaired_structure",
            "repaired_pdb_file",
            "mutations_applied",
            "nonstandard_mutations_applied",
        ]
    
    def get_menu_options(self) -> Dict[str, str]:
        """Menu options"""
        modeller_status = "available" if HAS_MODELLER else "unavailable"

        return {
            "analyze": "Analyze structure completeness",
            "process_structure": f"Apply repairs/mutations/caps to structure (MODELLER {modeller_status})",
            "summary": "View structural issues report",
            "export": "Export results to JSON",
        }

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
        has_analysis = self.results is not None
        has_modeller = HAS_MODELLER
        pending_mutations = workspace.get("pending_mutations", [])
        has_mutations = len(pending_mutations) > 0
        has_repaired = workspace.get("repaired_structure") is not None

        # Option 1: Analyze structure - needs an (experimental) structure
        if has_analysis:
            analyze_status, analyze_dep = OptionStatus.COMPLETED, ""
        elif self.can_process(workspace):
            analyze_status, analyze_dep = OptionStatus.AVAILABLE, ""
        else:
            analyze_status = OptionStatus.BLOCKED
            analyze_dep = self.availability_note(workspace) or "Load a structure first"
        options.append(MenuOption(
            key="1",
            description="Analyze structure completeness",
            status=analyze_status,
            dependency_text=analyze_dep,
        ))

        # Option 2: Apply repairs/mutations/caps. Analysis (or pending
        # mutations) is the real prerequisite; MODELLER is only needed for gap
        # filling and mutations — terminal capping works without it, so a
        # missing MODELLER no longer blocks this option.
        if has_repaired:
            # Already repaired - mark as completed
            status = OptionStatus.COMPLETED
            dep_text = ""
        elif not has_analysis and not has_mutations:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to analyze structure first] ○"
        elif not has_modeller:
            status = OptionStatus.READY
            dep_text = "[MODELLER off: capping only] ○"
        else:
            # Ready if we have analysis OR mutations
            status = OptionStatus.READY
            dep_text = ""

        options.append(MenuOption(
            key="2",
            description="Apply repairs/mutations/caps to structure",
            status=status,
            dependency_text=dep_text
        ))

        # Option 3: View report - requires analysis
        if has_analysis:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to analyze structure first] ○"

        options.append(MenuOption(
            key="3",
            description="View structural issues report",
            status=status,
            dependency_text=dep_text
        ))

        # Option 4: Export results - requires analysis
        if has_analysis:
            status = OptionStatus.READY
            dep_text = ""
        else:
            status = OptionStatus.BLOCKED
            dep_text = "[Need to analyze structure first] ○"

        options.append(MenuOption(
            key="4",
            description="Export results to JSON",
            status=status,
            dependency_text=dep_text
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
        has_analysis = self.results is not None
        has_modeller = HAS_MODELLER
        pending_mutations = workspace.get("pending_mutations", [])
        mutations_count = len(pending_mutations)
        has_repaired = workspace.get("repaired_structure") is not None

        # If already repaired, suggest next steps
        if has_repaired:
            return "✓ Structure repaired and saved to workspace. View report (option 3) or press [m] to return to the main menu"

        # Check for pending mutations even without analysis
        if not has_analysis and mutations_count > 0:
            if has_modeller:
                return f"Found {mutations_count} pending mutation{'s' if mutations_count != 1 else ''}. Apply with option 2, or analyze structure first (option 1)"
            else:
                return f"Found {mutations_count} pending mutation{'s' if mutations_count != 1 else ''}, but MODELLER is required to apply them"

        if not has_analysis:
            if not self.can_process(workspace):
                return f"{self.availability_note(workspace) or 'A structure is required'}. Load one via the Structure Loader."
            return "Start by analyzing structure completeness (option 1) to identify missing residues and gaps"
        else:
            # Check if there are issues that need repair
            if self.results:
                # Count issues by category
                missing_residues_count = 0
                missing_atoms_count = 0
                alternate_locations_count = 0

                # Count missing residues
                for method_results in self.results.get('missing_residues', {}).values():
                    missing_residues_count += sum(len(residues) for residues in method_results.values())

                # Count missing atoms
                for method_results in self.results.get('missing_atoms', {}).values():
                    missing_atoms_count += sum(len(atoms) for atoms in method_results.values())

                # Count alternate locations
                for method_results in self.results.get('alternate_locations', {}).values():
                    alternate_locations_count += sum(len(residues) for residues in method_results.values())

                total_issues = missing_residues_count + missing_atoms_count + alternate_locations_count

                if total_issues > 0 or mutations_count > 0:
                    # Build detailed breakdown
                    issue_parts = []
                    if missing_residues_count > 0:
                        issue_parts.append(f"{missing_residues_count} missing residue{'s' if missing_residues_count != 1 else ''}")
                    if missing_atoms_count > 0:
                        issue_parts.append(f"{missing_atoms_count} missing atom{'s' if missing_atoms_count != 1 else ''}")
                    if alternate_locations_count > 0:
                        issue_parts.append(f"{alternate_locations_count} alternate location{'s' if alternate_locations_count != 1 else ''}")
                    if mutations_count > 0:
                        issue_parts.append(f"{mutations_count} pending mutation{'s' if mutations_count != 1 else ''}")

                    issue_summary = ", ".join(issue_parts)

                    if has_modeller:
                        return f"Found {issue_summary}. Apply repairs/mutations with option 2, view details with option 3, or export with option 4"
                    else:
                        return f"Found {issue_summary}. View details with option 3 or export with option 4 (MODELLER needed for repairs)"
                else:
                    return "✓ No structural issues found. Add ACE/NME terminal caps with option 2, view report (option 3), or press [m] to return to the main menu"
            else:
                return "View analysis results (option 3), apply repairs (option 2), or export results (option 4)"
    
    def handle_menu_option(self, option: str) -> bool:
        """Handle menu selection"""
        if option == "analyze":
            self._analyze_structure()
            return True

        elif option == "process_structure":
            # Check prerequisites
            if not self.results:
                self.console.print("[yellow]Analysis required before applying repairs/mutations[/yellow]")
                return False

            # MODELLER is required for gap filling and mutations, but NOT for
            # terminal capping. Let the workflow run either way; it routes to a
            # capping-only path when MODELLER is unavailable.
            if not HAS_MODELLER:
                self.console.print("[yellow]MODELLER not available — gap filling and mutations are disabled.[/yellow]")
                self.console.print("[grey50]Terminal capping (ACE/NME) is still available; install MODELLER for full repairs.[/grey50]")

            self._unified_repair_workflow()
            return True

        elif option == "summary":
            # Enhanced menu should block this, but add fallback check
            if not self.results:
                self.console.print("[yellow]No analysis results available[/yellow]")
                return False

            self._show_summary()
            return True

        elif option == "export":
            # Enhanced menu should block this, but add fallback check
            if not self.results:
                self.console.print("[yellow]No results to export[/yellow]")
                return False

            filename = prompt_with_context(
                processor=self.processor,
                prompt="Output filename",
                default="completeness_results.json",
                module="Structure Completeness",
                description="Enter output filename for completeness results"
            )
            self._export_results(filename)
            return True

        return False
    
    def availability_note(self, workspace):
        """Menu note when unavailable (○). Repair needs an
        experimental structure — AlphaFold models are excluded."""
        return None if self.can_process(workspace) else \
            "Needs an experimental structure (RCSB or local PDB; not AlphaFold)"

    def can_process(self, workspace) -> bool:
        """Check if module can process workspace"""
        # Only process experimental structures (RCSB or local PDB files)
        # AlphaFold structures are complete by design and don't need repair
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)
        status = selector.get_structure_status()
        return status.get("has_experimental", False)
    
    def process(self, workspace) -> Any:
        """Process workspace (batch mode)"""
        # Get structure
        structure, _, _ = self._get_priority_structure(workspace)
        if not structure:
            return workspace
        
        # Get metadata from workspace (rcsb_metadata or local_metadata)
        metadata = workspace.get("rcsb_metadata") or workspace.get("local_metadata")
        filter_selections = workspace.get("filter_selections")
        
        # Run detection with default methods
        self.analyzer = StructureAnalyzer(structure, metadata, filter_selections, self.console)
        self.results = self.analyzer.detect_all_issues()
        
        # Store results
        workspace.set("completeness_results", self.results)
        
        # Check for auto-repair flags
        auto_repair = workspace.get("auto_repair", False)
        
        if auto_repair and HAS_MODELLER:
            # Run automated repair
            self._automated_repair(workspace)
        
        return workspace

    # ========================================================================
    # ANALYSIS METHODS
    # ========================================================================
    
    def _analyze_structure(self):
        """Interactive structure analysis with method selection"""
        workspace = self.processor.workspace
        
        # Get structure
        structure, structure_name, _ = self._get_priority_structure(workspace)
        if not structure:
            self.console.print("[red]No structure available[/red]")
            return
        
        self.console.print(f"[green]Using {structure_name}[/green]")
        
        # Method selection
        self.console.print("\n[bold]Select detection methods:[/bold]")
        
        # Missing residues
        self.console.print("\n[bold]Missing Residues Detection:[/bold]")
        self.console.print("1. REMARK 465 records", highlight=False)
        self.console.print("2. SEQRES comparison", highlight=False)
        self.console.print("3. FASTA comparison", highlight=False)
        self.console.print("4. Sequence gaps", highlight=False)
        self.console.print("5. Skip", highlight=False)

        mr_choice = prompt_with_context(
            processor=self.processor,
            prompt="Choose",
            choices=["1","2","3","4","5"],
            default="1",
            module="Structure Completeness - Analysis",
            description="Select missing residues detection method",
            options_map={
                "1": "REMARK 465 records",
                "2": "SEQRES comparison",
                "3": "FASTA comparison",
                "4": "Sequence gaps",
                "5": "Skip"
            }
        )
        mr_map = {"1": "remark_465", "2": "seqres_comparison", 
                "3": "fasta_comparison", "4": "sequence_gap", "5": None}
        mr_method = [mr_map[mr_choice]] if mr_map[mr_choice] else []
        
        # Missing atoms
        self.console.print("\n[bold]Missing Atoms Detection:[/bold]")
        self.console.print("1. REMARK 470 records", highlight=False)
        self.console.print("2. Template comparison", highlight=False)
        self.console.print("3. Skip", highlight=False)

        ma_choice = prompt_with_context(
            processor=self.processor,
            prompt="Choose",
            choices=["1","2","3"],
            default="1",
            module="Structure Completeness - Analysis",
            description="Select missing atoms detection method",
            options_map={
                "1": "REMARK 470 records",
                "2": "Template comparison",
                "3": "Skip"
            }
        )
        ma_map = {"1": "remark_470", "2": "template_comparison", "3": None}
        ma_method = [ma_map[ma_choice]] if ma_map[ma_choice] else []
        
        # Alternate locations
        self.console.print("\n[bold]Alternate Locations Detection:[/bold]")
        self.console.print("1. AltLoc identifiers", highlight=False)
        self.console.print("2. Occupancy values", highlight=False)
        self.console.print("3. Skip", highlight=False)

        al_choice = prompt_with_context(
            processor=self.processor,
            prompt="Choose",
            choices=["1","2","3"],
            default="1",
            module="Structure Completeness - Analysis",
            description="Select alternate locations detection method",
            options_map={
                "1": "AltLoc identifiers",
                "2": "Occupancy values",
                "3": "Skip"
            }
        )
        al_map = {"1": "altloc_identifier", "2": "occupancy_value", "3": None}
        al_method = [al_map[al_choice]] if al_map[al_choice] else []
        
        # Run detection
        # Get metadata from workspace (rcsb_metadata or local_metadata)
        metadata = workspace.get("rcsb_metadata") or workspace.get("local_metadata")
        filter_selections = workspace.get("filter_selections")

        # Metadata validation (silent - no debug output)
        # Debug statements removed for cleaner user experience

        self.analyzer = StructureAnalyzer(structure, metadata, filter_selections, self.console)
        self.results = self.analyzer.detect_all_issues(
            missing_residue_methods=mr_method,
            missing_atom_methods=ma_method,
            altloc_methods=al_method
        )
        
        # Store results
        workspace.set("completeness_results", self.results)

        # Highlight the residues that *flank* each gap. Missing residues
        # don't exist in the loaded structure, so selecting them would
        # match nothing — instead anchor the visual on the residues that
        # do exist on either side of each contiguous run, so the user
        # can see *where* gaps are even though the gap itself isn't
        # rendered. Missing atoms / altlocs deliberately skipped — they
        # get per-residue treatment in their own picker flow.
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            flanking = set()
            for method_results in self.results.get('missing_residues', {}).values():
                for chain_id, residues in method_results.items():
                    nums = sorted(
                        getattr(r, 'res_num', None) for r in residues
                        if getattr(r, 'res_num', None) is not None
                    )
                    if not nums:
                        continue
                    # Group consecutive runs and anchor on (start-1, end+1)
                    run_start = nums[0]
                    prev = nums[0]
                    for n in nums[1:] + [None]:
                        if n is None or n != prev + 1:
                            flanking.add((chain_id, run_start - 1))
                            flanking.add((chain_id, prev + 1))
                            run_start = n
                        prev = n if n is not None else prev
            if flanking:
                clauses = [f"(:{c} and {n})" for c, n in sorted(flanking)]
                _viewer.highlight(
                    " or ".join(clauses),
                    style="ball+stick",
                    color="#e31a1c",
                    label="fixer_missing_residues",
                )
        except Exception:
            pass

        # Get pending mutations for display
        pending_mutations = workspace.get("pending_mutations", [])
        pending_nonstandard_mutations = workspace.get("pending_nonstandard_mutations", [])

        # Display summary
        self.analyzer.display_summary_report(pending_mutations, pending_nonstandard_mutations)

        # Show sequence visualization automatically
        # Collect all chains that have any issues
        all_chains = set()

        # Get chains with missing residues
        missing_by_chain = {}
        for method, chain_results in self.results['missing_residues'].items():
            for chain_id, residues in chain_results.items():
                if residues:
                    all_chains.add(chain_id)
                    if chain_id not in missing_by_chain:
                        missing_by_chain[chain_id] = []
                    missing_by_chain[chain_id].extend(residues)

        # Get chains with missing atoms
        for method, chain_results in self.results['missing_atoms'].items():
            for chain_id, atoms in chain_results.items():
                if atoms:
                    all_chains.add(chain_id)

        # Get chains with alternate locations
        for method, chain_results in self.results['alternate_locations'].items():
            for chain_id, residues in chain_results.items():
                if residues:
                    all_chains.add(chain_id)

        # Get chains with pending mutations
        if pending_mutations:
            for mut_chain_id, _, _, _ in pending_mutations:
                all_chains.add(mut_chain_id)

        # Get chains with pending non-standard mutations
        if pending_nonstandard_mutations:
            for mut_chain_id, _, _, _, _ in pending_nonstandard_mutations:
                all_chains.add(mut_chain_id)

        # Show sequence for all chains with issues
        if all_chains:
            # Build a dict with all chains (empty list for chains without missing residues)
            all_chains_dict = {}
            for chain_id in sorted(all_chains):
                all_chains_dict[chain_id] = missing_by_chain.get(chain_id, [])
            self.analyzer.display_all_chain_sequences(all_chains_dict, pending_mutations, pending_nonstandard_mutations)

            # Offer detailed reports (optional) - only if there are issues to report
            if confirm_with_context(
                processor=self.processor,
                prompt="\nShow detailed reports?",
                default=False,
                module="Structure Completeness - Analysis",
                description="Show detailed analysis reports"
            ):
                self.analyzer.display_missing_residues_report()
                self.analyzer.display_missing_atoms_report()
                self.analyzer.display_alternate_locations_report()
                self.analyzer.display_mutations_report(pending_mutations, pending_nonstandard_mutations)
                
    def _show_summary(self):
        """Show summary report"""
        if self.analyzer:
            pending_mutations = self.processor.workspace.get("pending_mutations", [])
            pending_nonstandard_mutations = self.processor.workspace.get("pending_nonstandard_mutations", [])
            self.analyzer.display_summary_report(pending_mutations, pending_nonstandard_mutations)
        else:
            self.console.print("[yellow]No analysis results available[/yellow]")

    def _export_results(self, filename: str):
        """Export results to JSON"""
        if not self.results:
            self.console.print("[yellow]No results to export[/yellow]")
            return
        
        try:
            with open(filename, 'w') as f:
                json.dump(self.results, f, indent=2, default=str)
            self.console.print(f"[green]Results exported to {filename}[/green]")
        except Exception as e:
            self.console.print(f"[red]Export failed: {e}[/red]")
    
    # ========================================================================
    # ALTERNATE LOCATION HANDLING
    # ========================================================================

    def _handle_alternate_locations(self, structure: Structure, workspace: Any, save_to_workspace: bool = True) -> Optional[Structure]:
        """
        Handle alternate locations by prompting user to select which to keep.
        Returns the cleaned structure or None if cancelled.

        Args:
            structure: Input structure with alternate locations
            workspace: Workspace object
            save_to_workspace: If True, save cleaned structure to workspace. If False, just return it.
        """
        import copy
        from Bio.PDB import Select, PDBIO

        self.console.print("\n[bold cyan]═══ Alternate Location Selection ═══[/bold cyan]")
        self.console.print("The following residues have alternate locations. Please select which to keep:\n")

        # Best-effort live viewer: launches once, refocuses on the current
        # residue before each prompt. All viewer interactions are wrapped to
        # ensure a viewer failure never blocks the picker flow.
        altloc_viewer = self._setup_altloc_viewer(structure)

        # Collect all alternate locations
        selections = {}  # {(chain_id, res_num): selected_altloc}

        for method, chain_results in self.results.get('alternate_locations', {}).items():
            for chain_id, residues in chain_results.items():
                for res_key, altlocs in residues.items():
                    # Parse residue info
                    parts = res_key.split('_')
                    if len(parts) == 2:
                        res_name, res_num_str = parts
                        try:
                            res_num = int(res_num_str)
                        except ValueError:
                            continue

                        # Get occupancies
                        occupancies = {}
                        for model in structure:
                            if chain_id in model:
                                chain = model[chain_id]
                                for residue in chain:
                                    if residue.id[1] == res_num:
                                        for atom in residue:
                                            if hasattr(atom, 'is_disordered') and atom.is_disordered():
                                                if hasattr(atom, 'child_dict'):
                                                    for altloc_id, alt_atom in atom.child_dict.items():
                                                        if altloc_id.strip():
                                                            if altloc_id.strip() not in occupancies:
                                                                occupancies[altloc_id.strip()] = []
                                                            if hasattr(alt_atom, 'occupancy'):
                                                                occupancies[altloc_id.strip()].append(alt_atom.occupancy)
                                            else:
                                                altloc = atom.altloc.strip() if hasattr(atom, 'altloc') else ""
                                                if altloc:
                                                    if altloc not in occupancies:
                                                        occupancies[altloc] = []
                                                    if hasattr(atom, 'occupancy'):
                                                        occupancies[altloc].append(atom.occupancy)
                                        break
                                break

                        # Display options. With the live viewer open, print each
                        # alternate's viewer colour beside its occupancy so the
                        # red/blue on screen can be read against the numbers.
                        self.console.print(f"[bold]Chain {chain_id}, {res_name} {res_num}:[/bold]")
                        sorted_altlocs = sorted(altlocs)
                        avg_occs: Dict[str, str] = {}
                        for i, altloc in enumerate(sorted_altlocs, 1):
                            avg_occ = "?"
                            if altloc in occupancies and occupancies[altloc]:
                                avg_occ = f"{sum(occupancies[altloc]) / len(occupancies[altloc]):.2f}"
                            avg_occs[altloc] = avg_occ
                            if altloc_viewer is not None:
                                hex_color, color_name = self._altloc_color(altloc, i - 1)
                                self.console.print(
                                    f"  {i}. Alternate {altloc} (occupancy: {avg_occ})  "
                                    f"[{hex_color}]■ {color_name} in viewer[/{hex_color}]"
                                )
                            else:
                                self.console.print(f"  {i}. Alternate {altloc} (occupancy: {avg_occ})")

                        # Refocus the viewer on the current residue before prompting.
                        if altloc_viewer is not None:
                            self._focus_altloc_viewer(
                                altloc_viewer, chain_id, res_name, res_num, sorted_altlocs,
                                occupancies=avg_occs,
                            )

                        # Prompt for selection
                        choices = [str(i) for i in range(1, len(sorted_altlocs) + 1)]
                        choice = prompt_with_context(
                            processor=self.processor,
                            prompt=f"Select alternate to keep",
                            choices=choices,
                            default="1",
                            module="Structure Completeness - Altloc",
                            description=f"Select alternate for {res_name} {chain_id}:{res_num}",
                            options_map={str(i+1): f"Alternate {altloc}" for i, altloc in enumerate(sorted_altlocs)}
                        )
                        selected_altloc = sorted_altlocs[int(choice) - 1]
                        selections[(chain_id, res_num)] = selected_altloc
                        self.console.print(f"[green]✓ Will keep alternate {selected_altloc}[/green]\n")

        if not selections:
            self._teardown_altloc_viewer(altloc_viewer)
            return structure

        # Create a custom selector class
        class AltlocSelector(Select):
            def __init__(self, selections_dict):
                self.selections = selections_dict

            def accept_atom(self, atom):
                residue = atom.get_parent()
                chain = residue.get_parent()
                chain_id = chain.id
                res_num = residue.id[1]

                key = (chain_id, res_num)
                if key in self.selections:
                    selected_altloc = self.selections[key]

                    # Get the altloc for this atom
                    atom_altloc = atom.altloc.strip() if hasattr(atom, 'altloc') else ""

                    # If atom has an altloc, check if it matches the selected one
                    if atom_altloc:
                        return atom_altloc == selected_altloc

                    # If disordered but no altloc string, try child_dict
                    if hasattr(atom, 'is_disordered') and atom.is_disordered():
                        if hasattr(atom, 'child_dict'):
                            return selected_altloc in atom.child_dict

                    # No altloc means this atom is shared across all alternates
                    return True  # Keep atoms with no altloc (shared across alternates)

                return True  # Keep atoms in residues without alternates

        # Save cleaned structure
        self.console.print("[cyan]Removing unselected alternate locations...[/cyan]")

        # Write to temporary file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            tmp_path = tmp.name

        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp_path, AltlocSelector(selections))

        # Re-read the structure to get clean version
        parser = PDBParser(QUIET=True)
        cleaned_structure = parser.get_structure('cleaned', tmp_path)

        # Save to workspace only if requested (i.e., this is the final step)
        if save_to_workspace:
            repaired_pdb_path = "repaired_structure.pdb"
            io.set_structure(cleaned_structure)
            io.save(repaired_pdb_path)

            workspace.set("repaired_structure", cleaned_structure)
            workspace.set("repaired_pdb_file", os.path.abspath(repaired_pdb_path))

            self.console.print(f"[green]✓ Cleaned structure saved to workspace[/green]")
            self.console.print(f"[green]  • repaired_structure (BioPython Structure object)[/green]")
            self.console.print(f"[green]  • repaired_pdb_file: {os.path.abspath(repaired_pdb_path)}[/green]")
        else:
            self.console.print(f"[green]✓ Alternate locations removed[/green]")
            self.console.print(f"[grey50]  (Will be used as input for MODELLER)[/grey50]")

        # Clean up temp file
        os.unlink(tmp_path)

        self._teardown_altloc_viewer(altloc_viewer)
        return cleaned_structure

    # ========================================================================
    # ALT-LOC LIVE VIEWER (best-effort)
    # ========================================================================

    # Per-altloc colors used in the 3D viewer, with the name the prompt prints
    # next to each alternate so the user can tell which occupancy is which.
    _ALTLOC_PALETTE = {
        "A": ("#e74c3c", "red"),
        "B": ("#3498db", "blue"),
        "C": ("#2ecc71", "green"),
        "D": ("#f39c12", "orange"),
        "E": ("#9b59b6", "purple"),
        "F": ("#1abc9c", "teal"),
    }
    _ALTLOC_FALLBACK_PALETTE = [
        ("#e67e22", "dark orange"), ("#34495e", "slate"),
        ("#c0392b", "dark red"), ("#16a085", "sea green"),
    ]

    def _altloc_color(self, alt: str, index: int) -> Tuple[str, str]:
        """(hex, name) for an altloc letter; letters beyond F cycle the fallback palette."""
        return self._ALTLOC_PALETTE.get(
            alt.upper(),
            self._ALTLOC_FALLBACK_PALETTE[index % len(self._ALTLOC_FALLBACK_PALETTE)],
        )
    # Radius (Å) of the environment shell drawn around each alt-loc residue.
    _ALTLOC_ENV_DISTANCE = 5.0

    def _setup_altloc_viewer(self, structure: Structure):
        """Optionally snapshot the all-altloc structure and route the live
        viewer through the coordinator.

        The alt-loc picker auto-runs whenever alternate locations exist; the
        viewer aid that refocuses on each residue is *useful when wanted*
        but should not pop a browser tab in CLI mode unless the user has
        opted in. This method asks before doing anything, so the launch
        decision stays user-initiated rather than being a side-effect of
        the picker firing.

        Returns a small dict ``{pdb_path, prev_labels}`` on success (user
        opted in and the launch succeeded), None when the user declines or
        any step fails. ``prev_labels`` tracks the per-altloc annotation
        labels currently on screen so the next focus call can clear them
        before drawing the new residue's reps. Downstream callers already
        guard ``_focus_altloc_viewer(...)`` with ``if altloc_viewer is not
        None``, so returning None silently disables the per-residue refocus
        without disturbing the picker flow.
        """
        # Ask the user before launching. In CLI mode this prevents an
        # unbidden browser pop; in web-shell mode the iframe is already
        # there so the prompt is mostly a confirmation. Default=False so
        # the silent path is easy for users who don't want the aid.
        if not confirm_with_context(
            processor=self.processor,
            prompt="Launch the 3D viewer to help pick alternate locations?",
            default=False,
            module="Structure Completeness — Alt-Loc Picker",
            description="Optionally launch the structure viewer with per-residue refocus to aid alt-loc selection",
        ):
            return None

        # Show a fixed 5 Å shell of surrounding residues alongside each
        # alt-loc residue. Seeing the immediate neighbours (clashes, H-bond
        # partners, packing) makes it much easier to judge which conformer is
        # the right one; 5 Å is the tight first-contact shell.
        env_distance = self._ALTLOC_ENV_DISTANCE

        try:
            from Bio.PDB import PDBIO
            from Bio.PDB.NeighborSearch import NeighborSearch
            import tempfile
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer

            # Snapshot the input structure (with all altlocs) to a temp PDB
            # the viewer's HTTP server can read. The structure currently
            # in the coordinator may have been written without altlocs,
            # so we need our own copy that preserves the %A/%B records.
            tmp = tempfile.NamedTemporaryFile(
                mode="w", suffix="_altloc_input.pdb", delete=False
            )
            tmp_path = tmp.name
            tmp.close()
            io = PDBIO()
            io.set_structure(structure)
            io.save(tmp_path)

            # Build a NeighborSearch over the first model's atoms once, so each
            # per-residue refocus can resolve its environment shell cheaply.
            # Disordered atoms yield their representative coord, which is fine
            # for proximity detection. Tolerate failure — the picker and the
            # residue-only view must still work without the environment shell.
            neighbor_search = None
            if env_distance > 0:
                try:
                    model = next(iter(structure))
                    atoms = list(model.get_atoms())
                    if atoms:
                        neighbor_search = NeighborSearch(atoms)
                except Exception as e:
                    logger.debug(f"Could not build alt-loc neighbor search: {e}")

            # force=True is correct here because the user just explicitly
            # opted in via the prompt above — this is now a user-initiated
            # view, not an auto-fired workflow waypoint.
            _viewer.show_structure(tmp_path, force=True)
            self.console.print(
                "[grey50]Live 3D viewer is open in your browser; it will refocus "
                "on each residue as you choose.[/grey50]"
            )
            return {
                "pdb_path": tmp_path,
                "prev_labels": [],
                "env_distance": env_distance,
                "neighbor_search": neighbor_search,
            }
        except Exception as e:
            logger.debug(f"Live alt-loc viewer unavailable: {e}")
            return None

    def _focus_altloc_viewer(
        self,
        viewer_state: Dict[str, Any],
        chain_id: str,
        res_name: str,
        res_num: int,
        sorted_altlocs: List[str],
        occupancies: Optional[Dict[str, str]] = None,
    ) -> None:
        """Refocus the coordinator viewer on a single residue's altlocs.

        Clears the previous residue's per-altloc reps, then draws a grey
        licorice scaffold (whole residue) plus one ball+stick rep per
        altloc using the same NGL ``%A`` selector syntax we used in the
        standalone implementation. ``focused=True`` on the scaffold
        triggers the auto-centre — same effect as the old viewer's
        ribbon-hide trick.

        When the user opted into an environment shell (``env_distance``),
        the residues within that radius are drawn as a faint line overlay so
        the surrounding packing/clashes/H-bond partners are visible. The
        shell is a *non-focused* overlay, so the camera still centres on the
        alt-loc residue rather than the whole neighbourhood.
        """
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer

            for stale_label in viewer_state.get("prev_labels", []):
                _viewer.unhighlight(stale_label)

            base = f":{chain_id} and {res_num}"
            new_labels = ["altloc_scaffold"]
            _viewer.highlight(
                base, style="licorice", color="#bdc3c7",
                label="altloc_scaffold", focused=True,
            )
            for i, alt in enumerate(sorted_altlocs):
                color, color_name = self._altloc_color(alt, i)
                label = f"altloc_{alt}"
                new_labels.append(label)
                occ = (occupancies or {}).get(alt)
                display = f"Alt {alt} ({color_name})" + (f", occ {occ}" if occ else "")
                _viewer.highlight(
                    f"{base} and %{alt}", style="ball+stick",
                    color=color, label=label, display_label=display,
                )

            # Draw the surrounding environment shell, if requested. Added after
            # the scaffold/altloc reps so the focused scaffold has already set
            # the camera; this overlay is non-focused and only adds context.
            env_selection = self._altloc_environment_selection(
                viewer_state, chain_id, res_num
            )
            if env_selection:
                _viewer.highlight(
                    env_selection, style="line", color="#7f8c8d",
                    label="altloc_environment", opacity=0.6,
                )
                new_labels.append("altloc_environment")

            viewer_state["prev_labels"] = new_labels
        except Exception as e:
            logger.debug(f"Could not refocus alt-loc viewer: {e}")

    def _altloc_environment_selection(
        self,
        viewer_state: Dict[str, Any],
        chain_id: str,
        res_num: int,
    ) -> Optional[str]:
        """Build an NGL selection for residues near the given alt-loc residue.

        Uses the prebuilt ``NeighborSearch`` (stored in ``viewer_state``) to
        find every residue with an atom within ``env_distance`` Å of any atom
        of the target residue, excludes the target itself, and emits an NGL
        string grouped by chain — e.g. ``(:A and (54 or 55 or 90)) or (:B and (12))``.

        Returns None when no environment was requested, the search is
        unavailable, or nothing falls within the shell (so the caller simply
        skips the overlay).
        """
        env_distance = viewer_state.get("env_distance", 0.0)
        neighbor_search = viewer_state.get("neighbor_search")
        if not env_distance or neighbor_search is None:
            return None

        try:
            # Gather the target residue's atoms across all of its altlocs so
            # the shell is measured from the full residue envelope.
            target_atoms = [
                atom for atom in neighbor_search.atom_list
                if atom.get_parent().id[1] == res_num
                and atom.get_parent().get_parent().id == chain_id
            ]
            if not target_atoms:
                return None

            neighbor_residues = set()
            for atom in target_atoms:
                for residue in neighbor_search.search(
                    atom.coord, env_distance, level="R"
                ):
                    neighbor_residues.add(residue)

            # Group neighbours by chain, dropping the target residue itself.
            by_chain: Dict[str, Set[int]] = defaultdict(set)
            for residue in neighbor_residues:
                r_chain = residue.get_parent().id
                r_num = residue.id[1]
                if r_chain == chain_id and r_num == res_num:
                    continue
                by_chain[r_chain].add(r_num)

            if not by_chain:
                return None

            groups = []
            for r_chain in sorted(by_chain):
                nums = " or ".join(str(n) for n in sorted(by_chain[r_chain]))
                groups.append(f"(:{r_chain} and ({nums}))")
            return " or ".join(groups)
        except Exception as e:
            logger.debug(f"Could not compute alt-loc environment: {e}")
            return None

    def _teardown_altloc_viewer(self, viewer_state: Optional[Dict[str, Any]]) -> None:
        """Clear the per-altloc reps so they don't leak into later hooks.

        Leaves the snapshot structure loaded — Hook 14 (post-MODELLER
        re-detect) will swap to the repaired PDB anyway.
        """
        if not viewer_state:
            return
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer
            for stale_label in viewer_state.get("prev_labels", []):
                _viewer.unhighlight(stale_label)
            viewer_state["prev_labels"] = []
        except Exception as e:
            logger.debug(f"Could not tear down alt-loc viewer: {e}")

    # ========================================================================
    # UNIFIED REPAIR WORKFLOW
    # ========================================================================

    def _unified_repair_workflow(self):
        """Complete repair workflow with all components"""
        workspace = self.processor.workspace

        # Ensure analysis has been run
        if not self.results:
            self.console.print("[yellow]Running structure analysis...[/yellow]")
            self._analyze_structure()
            if not self.results:
                return

        # Get structure
        structure, _, _ = self._get_priority_structure(workspace)
        if not structure:
            self.console.print("[red]No structure available[/red]")
            return

        # Without MODELLER we cannot fill gaps or apply mutations, but adding
        # ACE/NME terminal caps is a pure PDB edit that needs no MODELLER. Offer
        # that capping-only path rather than disabling the whole feature.
        if not HAS_MODELLER:
            self.console.print(
                "[yellow]MODELLER not available — gap filling and mutations are disabled.[/yellow]"
            )
            self.console.print("[grey50]Terminal capping (ACE/NME) is still available.[/grey50]")
            orchestrator = RepairOrchestrator(self.console, processor=self.processor)
            caps_plan = RepairPlan()
            orchestrator._offer_terminal_capping(structure, caps_plan)
            if caps_plan.has_caps:
                self._execute_repair_plan(caps_plan, structure, workspace)
            else:
                self.console.print("[grey50]No terminal caps selected.[/grey50]")
            return

        # Get pending mutations (both standard and non-standard)
        pending_mutations = workspace.get("pending_mutations", [])
        pending_nonstandard_mutations = workspace.get("pending_nonstandard_mutations", [])
        
        # Validate mutations if present
        if pending_mutations:
            orchestrator = RepairOrchestrator(self.console, processor=self.processor)
            valid_mutations, warnings = orchestrator.validate_mutations(
                pending_mutations, self.results
            )
            
            if warnings:
                self.console.print("\n[bold yellow]Mutation Warnings:[/bold yellow]")
                for warning in warnings:
                    self.console.print(f"  [yellow]⚠ {warning}[/yellow]")

                if not confirm_with_context(
                    processor=self.processor,
                    prompt="\nProceed with these mutations?",
                    default=True,
                    module="Structure Completeness - Repair",
                    description="Confirm mutations with warnings"
                ):
                    self.console.print("[yellow]Repair cancelled[/yellow]")
                    return
            
            pending_mutations = valid_mutations
        
        # Check if there's anything to repair with MODELLER
        has_modeller_issues = any([
            any(v for v in self.results.get('missing_residues', {}).values()),
            any(v for v in self.results.get('missing_atoms', {}).values()),
            pending_mutations
        ])

        # Handle alternate locations first (doesn't require MODELLER)
        has_alternate_locations = any(v for v in self.results.get('alternate_locations', {}).values())
        if has_alternate_locations:
            # Only save to workspace if this is the final step (no MODELLER or non-standard mutations needed)
            save_now = not has_modeller_issues and not pending_nonstandard_mutations
            structure = self._handle_alternate_locations(structure, workspace, save_to_workspace=save_now)
            if structure is None:
                return  # User cancelled

        # If no MODELLER issues, check for non-standard mutations
        if not has_modeller_issues:
            if pending_nonstandard_mutations:
                # Apply non-standard mutations directly (no MODELLER needed)
                self._apply_nonstandard_mutations_only(structure, pending_nonstandard_mutations, workspace)
                return

            # No missing residues/atoms/mutations to repair — but the user may
            # still want to add ACE/NME caps to otherwise-complete termini.
            # This turns the old "nothing to do" dead-end into a capping offer.
            orchestrator = RepairOrchestrator(self.console, processor=self.processor)
            caps_plan = RepairPlan()
            orchestrator._offer_terminal_capping(structure, caps_plan)
            if caps_plan.has_caps:
                self._execute_repair_plan(caps_plan, structure, workspace)
            elif has_alternate_locations:
                self.console.print("[green]✓ Structure processing complete (alternate locations resolved)[/green]")
            else:
                self.console.print("[green]No structural issues found[/green]")
            return

        # Create repair plan
        orchestrator = RepairOrchestrator(self.console, processor=self.processor)
        plan = orchestrator.create_repair_plan(
            self.results, structure, pending_mutations
        )

        # create_repair_plan() returns None when the user declines the plan at
        # the "Proceed with this repair plan?" prompt. Without this guard the
        # next line dereferences None and crashes the whole menu/workflow stack
        # ('NoneType' object has no attribute 'needs_modeller') instead of
        # returning cleanly to the menu.
        if plan is None:
            return

        # (Terminal capping is now offered inside create_repair_plan, before the
        # plan summary/confirmation, so caps appear in the plan the user approves.)

        if not plan.needs_modeller and not plan.has_caps and not plan.has_ter:
            if has_alternate_locations:
                # MODELLER turned out to be unnecessary (e.g. the only detected
                # "issue" was missing atoms, which tLEaP builds at topology
                # generation). But we already resolved alternate locations above
                # with save_to_workspace=False, in anticipation of a MODELLER run
                # that isn't happening. Persist the altloc-cleaned structure now
                # so the resolution isn't silently discarded.
                self._save_repaired_structure(structure, workspace)
                self.console.print(
                    "[green]✓ Structure processing complete (alternate locations resolved)[/green]"
                )
            else:
                self.console.print("[yellow]No repairs planned[/yellow]")
            return

        # Execute repair
        self._execute_repair_plan(plan, structure, workspace)

    def _save_repaired_structure(self, structure: Structure, workspace: Any) -> str:
        """Write ``structure`` to repaired_structure.pdb and register it in the
        workspace (``repaired_structure`` object + ``repaired_pdb_file`` path).

        Central helper for the paths that finalize a structure without MODELLER
        (altloc-only resolution, non-standard mutations) so they all persist the
        result identically. Returns the absolute path written.
        """
        output_file = "repaired_structure.pdb"
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file)
        workspace.set("repaired_structure", structure)
        workspace.set("repaired_pdb_file", os.path.abspath(output_file))
        return os.path.abspath(output_file)

    def _apply_nonstandard_mutations_only(self, structure: Structure,
                                           nonstandard_mutations: List[Tuple],
                                           workspace: Any) -> None:
        """
        Apply non-standard mutations without MODELLER (for cases with no structural issues).

        This is a simplified path for when we only have non-standard mutations to apply,
        without any missing residues/atoms that would require MODELLER.
        """
        self.console.print("\n[bold cyan]═══ Non-Standard Mutations ═══[/bold cyan]")

        # Create a minimal mapper (no actual mappings needed since no MODELLER)
        mapper = ResidueMapper(self.console)

        # Apply mutations
        nonstandard_applicator = NonStandardMutationApplicator(self.console)
        redox_sites = workspace.get("detected_redox_sites", [])

        modified_structure, warnings = nonstandard_applicator.apply_nonstandard_mutations(
            structure,
            nonstandard_mutations,
            mapper,
            redox_sites
        )

        if warnings:
            self.console.print(f"\n[yellow]Warnings during non-standard mutation application:[/yellow]")
            for warning in warnings:
                self.console.print(f"  • {warning}")

        # Save modified structure
        output_file = "repaired_structure.pdb"
        io = PDBIO()
        io.set_structure(modified_structure)
        io.save(output_file)

        # Update workspace
        workspace.set("repaired_structure", modified_structure)
        workspace.set("repaired_pdb_file", os.path.abspath(output_file))
        workspace.set("pending_nonstandard_mutations", [])
        workspace.set("nonstandard_mutations_applied", nonstandard_mutations)

        # Update redox sites if they were modified
        if redox_sites:
            workspace.set("detected_redox_sites", redox_sites)

        self.console.print(f"\n[bold green]✓ Non-standard mutations applied![/bold green]")
        self.console.print(f"[cyan]Output: {os.path.abspath(output_file)}[/cyan]")

    def _apply_ter_records(self, pdb_file: str, plan: RepairPlan,
                           mapper: 'ResidueMapper') -> None:
        """Insert a TER record after the residue that precedes each TER-marked
        gap, so tLEaP reads the break as a chain end and does not build a long
        bond across the missing residue(s).

        TER is a file-only concept — Bio.PDB's ``PDBIO`` neither preserves nor
        emits TER — so this post-processes the written PDB text directly, after
        every ``io.save()``. Numbers are translated through ``mapper`` so the
        anchor residue is found under its *final* chain/number (MODELLER output
        numbering when fills ran; identity for the decline / TER-only path).
        """
        if not plan.has_ter:
            return

        # (final_chain, final_resnum) of the residue immediately before each gap.
        anchors = set()
        for seg in plan.segments_to_ter:
            try:
                ident = mapper.get_final_identity(
                    ResidueIdentity(chain_id=seg.chain_id,
                                    res_num=seg.start_num - 1, res_name="")
                )
                anchors.add((ident.chain_id, int(ident.res_num)))
            except Exception:
                anchors.add((seg.chain_id, seg.start_num - 1))

        if not anchors:
            return

        try:
            with open(pdb_file) as fh:
                lines = fh.readlines()
        except OSError:
            return

        out = []
        inserted = 0
        for i, line in enumerate(lines):
            out.append(line)
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21]
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            if (chain, resnum) not in anchors:
                continue
            # Only after the LAST atom of the anchor residue: the next record
            # must not be another atom of the same residue, and must not already
            # be a TER.
            nxt = lines[i + 1] if i + 1 < len(lines) else ""
            if nxt.startswith("TER"):
                continue
            if nxt.startswith(("ATOM", "HETATM")):
                try:
                    if nxt[21] == chain and int(nxt[22:26]) == resnum:
                        continue  # same residue, not the last atom yet
                except ValueError:
                    pass
            out.append("TER\n")
            inserted += 1

        if inserted:
            with open(pdb_file, "w") as fh:
                fh.writelines(out)
            self.console.print(
                f"[grey50]Inserted {inserted} TER record(s) at unfilled chain "
                f"breaks so tLEaP won't bond across the gaps.[/grey50]"
            )

    def _execute_repair_plan(self, plan: RepairPlan, structure: Structure, workspace: Any):
        """Execute the complete repair plan"""
        # Setup
        work_dir = "modeller_work"
        os.makedirs(work_dir, exist_ok=True)
        
        # Initialize components
        mapper = ResidueMapper(self.console)
        modeller_interface = ModellerInterface(self.console)
        capping_handler = CappingHandler(self.console)
        redox_sync = RedoxSiteSync(self.console)
        
        # Parse structure metadata
        structure_metadata = self._parse_structure_metadata(structure, workspace)
        
        # Build MODELLER mapping
        # IMPORTANT: Combine all segments for the same chain (don't overwrite!)
        residues_to_fill = {}
        for seg in plan.segments_to_fill:
            if seg.chain_id not in residues_to_fill:
                residues_to_fill[seg.chain_id] = []
            residues_to_fill[seg.chain_id].extend(seg.residues)

        if plan.needs_modeller:
            mapper.build_modeller_mapping(structure, residues_to_fill, structure_metadata)
        else:
            # Caps-only (or no fills/mutations): MODELLER will not run, so the
            # structure keeps its original chain IDs and numbering. Build an
            # identity mapping keyed by the real chain IDs — build_modeller_mapping
            # would key by predicted (renamed) chains and break the capping lookup.
            mapper.build_identity_mapping(structure)

        # Run MODELLER if needed
        repaired_structure = structure
        if plan.needs_modeller:
            self.console.print("\n[cyan]Running MODELLER...[/cyan]")

            # Save input structure
            input_pdb = os.path.join(work_dir, "input.pdb")
            io = PDBIO()
            io.set_structure(structure)
            io.save(input_pdb)
            
            # Build sequences
            chain_sequences = modeller_interface.build_sequences(
                structure, plan, structure_metadata
            )
            
            self.console.print("\n[bold cyan]═══ Sequence Summary ═══[/bold cyan]")
            modeller_interface.display_sequence_summary(chain_sequences)
            
            # ADD THIS: Display alignments
            modeller_interface.visualize_all_chain_alignments(chain_sequences)
            
            # Ask user to confirm before running MODELLER
            if not confirm_with_context(
                processor=self.processor,
                prompt="\nProceed with MODELLER repair?",
                default=True,
                module="Structure Completeness - Repair",
                description="Confirm MODELLER repair"
            ):
                self.console.print("[yellow]Repair cancelled by user[/yellow]")
                return False
            
            # Create alignment
            aln_file = os.path.join(work_dir, "alignment.ali")
            modeller_interface.create_alignment_file(chain_sequences, aln_file, structure)

            # Resolve the residues MODELLER will build (filled gaps / new
            # residues) from original numbering to MODELLER output numbering,
            # so run_modeller can assess them per-residue. final_mappings is
            # still empty here (caps are added later), so get_final_identity
            # yields pure MODELLER numbering.
            built_residues = set()
            for orig_chain, orig_num in self._find_entirely_new_residues(structure, plan):
                try:
                    ident = mapper.get_final_identity(
                        ResidueIdentity(chain_id=orig_chain, res_num=orig_num,
                                        res_name="")
                    )
                    built_residues.add((ident.chain_id, int(ident.res_num)))
                except Exception:
                    pass

            # Run MODELLER
            if built_residues:
                self.console.print(
                    "[grey50]Optimizing only the rebuilt residues; all resolved "
                    "atoms (including any metal sites) are held fixed.[/grey50]")
            success, message, repaired_structure = modeller_interface.run_modeller(
                work_dir, input_pdb, aln_file, built_residues=built_residues
            )
            
            if not success:
                self.console.print(f"[red]MODELLER failed: {message}[/red]")
                return
        
        # Add caps if needed
        final_output = "repaired_structure.pdb"
        if plan.has_caps:
            # Insert ACE/NME caps and renumber. add_caps() → _renumber_structure()
            # populates mapper.final_mappings authoritatively from the renumbered
            # file itself, so the redox-site synchronizer resolves atoms to their
            # true post-cap positions. (The old separately-predicted capping map
            # assumed MODELLER's global cross-chain numbering and mis-keyed
            # renamed chains, leaving renamed-chain sites unmapped.)
            final_structure = capping_handler.add_caps(
                repaired_structure, plan, mapper, final_output
            )
        else:
            # Save without caps
            io = PDBIO()
            io.set_structure(repaired_structure)
            io.save(final_output)
            final_structure = repaired_structure

        # Apply non-standard mutations (if any)
        nonstandard_mutations = workspace.get("pending_nonstandard_mutations", [])
        if nonstandard_mutations:
            self.console.print("\n[bold cyan]═══ Non-Standard Mutations ═══[/bold cyan]")
            nonstandard_applicator = NonStandardMutationApplicator(self.console)

            # Get redox sites before applying mutations (for update tracking)
            redox_sites = workspace.get("detected_redox_sites", [])

            # Apply mutations
            final_structure, ns_warnings = nonstandard_applicator.apply_nonstandard_mutations(
                final_structure,
                nonstandard_mutations,
                mapper,
                redox_sites
            )

            if ns_warnings:
                self.console.print(f"\n[yellow]Warnings during non-standard mutation application:[/yellow]")
                for warning in ns_warnings:
                    self.console.print(f"  • {warning}")

            # Save updated structure
            io = PDBIO()
            io.set_structure(final_structure)
            io.save(final_output)

            # Clear pending non-standard mutations and mark as applied
            workspace.set("pending_nonstandard_mutations", [])
            # Store applied non-standard mutations separately for tracking
            workspace.set("nonstandard_mutations_applied", nonstandard_mutations)

        # Synchronize RedoxSites
        redox_sites = workspace.get("detected_redox_sites", [])
        if redox_sites:
            sync_summary = redox_sync.synchronize_sites(redox_sites, final_structure, mapper)

            # Save updated redox sites back to workspace
            workspace.set("detected_redox_sites", redox_sites)
            self.console.print(f"[green]✓ {sync_summary['sites_updated']} RedoxSite(s) synchronized with repaired structure[/green]")

            # Check if new residues were added and prompt for re-detection
            entirely_new_residues = self._find_entirely_new_residues(structure, plan)
            if entirely_new_residues:
                choice = self._prompt_redox_redetection(
                    entirely_new_residues, structure,
                    plan=plan, mapper=mapper,
                    final_structure=final_structure,
                    final_pdb_file=os.path.abspath(final_output),
                    redox_sites=redox_sites,
                )
                if choice == "re-detect":
                    self._rerun_redox_detection(final_structure, final_output, workspace)
                    # TER records are a file-level edit; re-detection just
                    # rewrote the file, so stamp them in afterward.
                    self._apply_ter_records(final_output, plan, mapper)
                    # Early return - new detection replaces synchronized sites
                    self.console.print(f"\n[bold green]✓ Structure repair completed![/bold green]")
                    self.console.print(f"[cyan]Output: {final_output}[/cyan]")
                    return

        # Insert any TER records last: PDBIO.save() does not round-trip TER, so
        # this must run after every Bio.PDB write to the output file.
        self._apply_ter_records(final_output, plan, mapper)

        # Update workspace
        workspace.set("repaired_structure", final_structure)
        workspace.set("repaired_pdb_file", os.path.abspath(final_output))
        
        # Clear pending mutations
        if plan.mutations:
            workspace.set("pending_mutations", [])
            workspace.set("mutations_applied", plan.mutations)
        
        self.console.print(f"\n[bold green]✓ Structure repair completed![/bold green]")
        self.console.print(f"[cyan]Output: {final_output}[/cyan]")
    
    def _find_entirely_new_residues(self, original_structure: Structure, plan: RepairPlan) -> List[Tuple[str, int]]:
        """Find residues that were completely missing and added by MODELLER"""
        # Get all residue numbers that existed in original structure
        existing_residues = set()
        for model in original_structure:
            for chain in model:
                for residue in chain:
                    existing_residues.add((chain.id, residue.id[1]))

        # Check which filled residues were entirely new
        entirely_new = []
        for segment in plan.segments_to_fill:
            for residue in segment.residues:
                residue_key = (residue.chain_id, residue.res_num)
                if residue_key not in existing_residues:
                    entirely_new.append(residue_key)

        # Also check mutations that add residues
        for mutation in plan.mutations:
            # Mutations are tuples: (chain_id, res_num, from_aa, to_aa)
            chain_id, res_num, from_aa, to_aa = mutation
            residue_key = (chain_id, res_num)
            if residue_key not in existing_residues:
                entirely_new.append(residue_key)

        return entirely_new

    def _prompt_redox_redetection(self, entirely_new_residues: List[Tuple[str, int]],
                                   original_structure: Structure,
                                   plan: Optional['RepairPlan'] = None,
                                   mapper: Optional['ResidueMapper'] = None,
                                   final_structure: Optional[Structure] = None,
                                   final_pdb_file: Optional[str] = None,
                                   redox_sites: Optional[list] = None) -> str:
        """Prompt user for how to handle RedoxSite updates when new residues were added"""
        from rich.prompt import Prompt

        # Categorize new residues by position
        n_terminal = []
        c_terminal = []
        internal = []

        STANDARD_AMINO_ACIDS = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE"
        }

        # Get protein residue ranges per chain
        chain_protein_ranges = {}
        for model in original_structure:
            for chain in model:
                protein_residues = []
                for residue in chain:
                    if residue.resname in STANDARD_AMINO_ACIDS:
                        protein_residues.append(residue.id[1])

                if protein_residues:
                    chain_protein_ranges[chain.id] = {
                        'min': min(protein_residues),
                        'max': max(protein_residues)
                    }

        # Classify new residues
        for chain_id, res_num in entirely_new_residues:
            if chain_id in chain_protein_ranges:
                protein_range = chain_protein_ranges[chain_id]
                if res_num < protein_range['min']:
                    n_terminal.append(f"{chain_id}:{res_num}")
                elif res_num > protein_range['max']:
                    c_terminal.append(f"{chain_id}:{res_num}")
                else:
                    internal.append(f"{chain_id}:{res_num}")
            else:
                internal.append(f"{chain_id}:{res_num}")

        # Display information
        self.console.print(f"\n[yellow]⚠️  MODELLER added {len(entirely_new_residues)} new residues[/yellow]")
        if n_terminal:
            self.console.print(f"[cyan]N-terminal: {', '.join(n_terminal)}[/cyan]")
        if c_terminal:
            self.console.print(f"[cyan]C-terminal: {', '.join(c_terminal)}[/cyan]")
        if internal:
            self.console.print(f"[cyan]Internal: {', '.join(internal)}[/cyan]")

        self.console.print("These new residues could potentially be part of redox sites.")

        # Switch the viewer to the repaired structure and overlay three
        # separate representations so the user can independently toggle
        # each before deciding whether to re-detect: built residues
        # (MODELLER fills + new mutations), added caps (ACE/NME), and
        # the synchronized redox sites. Built residues need translation
        # via the residue mapper because MODELLER renumbers chains and
        # residues globally; caps are detected by residue name in the
        # final structure (simpler than reconstructing insertion
        # positions); redox sites are already in final numbering after
        # synchronize_sites().
        try:
            from proprep.structure_prep.viewer_coordinator import viewer as _viewer

            if final_pdb_file:
                _viewer.show_structures([final_pdb_file])

            _viewer.unhighlight("fixer_post_modeller_built")
            _viewer.unhighlight("fixer_post_added_caps")
            _viewer.unhighlight("fixer_post_redox_sites")

            if plan and mapper:
                built_pairs = set()
                for segment in plan.segments_to_fill:
                    for r in segment.residues:
                        orig = ResidueIdentity(
                            chain_id=r.chain_id,
                            res_num=r.res_num,
                            res_name=getattr(r, 'res_name', 'XXX'),
                        )
                        final = mapper.get_final_identity(orig)
                        built_pairs.add((final.chain_id, final.res_num))
                for mut in getattr(plan, 'mutations', []) or []:
                    chain_id, res_num, _, to_aa = mut
                    orig = ResidueIdentity(
                        chain_id=chain_id, res_num=res_num, res_name=to_aa,
                    )
                    final = mapper.get_final_identity(orig)
                    built_pairs.add((final.chain_id, final.res_num))
                if built_pairs:
                    clauses = [f"(:{c} and {n})" for c, n in sorted(built_pairs)]
                    _viewer.highlight(
                        " or ".join(clauses), style="ball+stick",
                        color="#1f78b4", label="fixer_post_modeller_built",
                    )

            if final_structure is not None:
                cap_pairs = set()
                for model in final_structure:
                    for chain in model:
                        for residue in chain:
                            if residue.resname in ("ACE", "NME", "NHE"):
                                cap_pairs.add((chain.id, residue.id[1]))
                if cap_pairs:
                    clauses = [f"(:{c} and {n})" for c, n in sorted(cap_pairs)]
                    _viewer.highlight(
                        " or ".join(clauses), style="ball+stick",
                        color="#ff7f00", label="fixer_post_added_caps",
                    )

            if redox_sites:
                redox_pairs = set()
                for site in redox_sites:
                    for c in getattr(site, 'centers', []) or []:
                        chain = getattr(c, 'chain', None)
                        resid = getattr(c, 'resid', None)
                        if chain and resid is not None:
                            redox_pairs.add((chain, resid))
                    for a in getattr(site, 'atoms', []) or []:
                        chain = getattr(a, 'chain', None)
                        resid = getattr(a, 'resid', None)
                        if chain and resid is not None:
                            redox_pairs.add((chain, resid))
                if redox_pairs:
                    clauses = [f"(:{c} and {n})" for c, n in sorted(redox_pairs)]
                    _viewer.highlight(
                        " or ".join(clauses), style="ball+stick",
                        color="#e31a1c", label="fixer_post_redox_sites",
                    )
        except Exception:
            pass

        choice = prompt_with_context(
            processor=self.processor,
            prompt="\nHow would you like to proceed?\n"
            "[bold]re-detect[/bold] - Re-run redox site detection on repaired structure (may include new residues)\n"
            "[bold]proceed[/bold] - Keep updated existing RedoxSite objects only\n"
            "Choose",
            choices=["re-detect", "proceed"],
            default="re-detect",
            module="Structure Completeness - Redox Integration",
            description="Redox site handling after repair",
            options_map={
                "re-detect": "Re-run redox site detection",
                "proceed": "Keep updated existing RedoxSite objects only"
            }
        )

        return choice

    def _rerun_redox_detection(self, final_structure: Structure, final_pdb_file: str, workspace: Any):
        """Re-run redox site detection on the repaired structure"""
        self.console.print("\n[cyan]Re-running redox site detection on repaired structure...[/cyan]")

        # Clear existing RedoxSite objects
        workspace.set("detected_redox_sites", [])

        # Update workspace to use repaired structure
        workspace.set("filtered_structure", final_structure)
        workspace.set("pdb_file", final_pdb_file)

        # Import and run comprehensive redox detector
        try:
            from .comprehensive_redox_detector import ComprehensiveRedoxDetector

            detector = ComprehensiveRedoxDetector()

            # Get processor from workspace if available
            processor = workspace.get("processor")
            if processor:
                detector.set_processor(processor)

            # Run detection
            new_redox_sites = detector.detect_redox_sites(
                structure=final_structure,
                pdb_file=final_pdb_file,
                workspace=workspace
            )

            if new_redox_sites:
                self.console.print(f"[green]✓ Detected {len(new_redox_sites)} redox sites in repaired structure[/green]")
                workspace.set("detected_redox_sites", new_redox_sites)
            else:
                self.console.print("[yellow]No redox sites detected in repaired structure[/yellow]")

        except Exception as e:
            self.console.print(f"[red]Error during re-detection: {e}[/red]")

    # ========================================================================
    # HELPER METHODS
    # ========================================================================

    def _get_priority_structure(self, workspace) -> Tuple[Optional[Structure], str, bool]:
        """Get experimental structure with priority: filtered > interactive selection

        Only works with experimental structures (RCSB or local PDB files).
        AlphaFold structures are complete by design and don't need repair.

        Uses StructureSelector with source_filter="experimental" to exclude
        AlphaFold and other predicted structures.
        """
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console, processor=self.processor)

        # Priority 1: Check for filtered structure first
        filtered = selector.get_structure_by_key(
            "filtered_structure",
            require_exists=False  # BioPython object, not file path
        )
        if filtered:
            return filtered, "filtered structure", False

        # Priority 2: Select from experimental structures only (excludes AlphaFold)
        result = selector.get_structure_object(
            source_filter="experimental",  # Only RCSB and local PDB files
            interactive=True,
            return_key=True,
            silent=False,
        )

        if result is None:
            self.console.print(
                "[yellow]No experimental structures available. "
                "This module only works with RCSB or local PDB files.[/yellow]"
            )
            return None, None, False

        structure, structure_key = result

        # Map structure key to display name
        source_mapping = {
            "rcsb_structure": "RCSB PDB",
            "local_structure": "local"
        }
        source_name = source_mapping.get(structure_key, "selected structure")

        return structure, source_name, True
    
    def _parse_structure_metadata(self, structure: Structure, workspace: Any) -> Dict:
        """Parse structure metadata from various sources"""
        metadata = {}

        # Try to get from workspace (rcsb_metadata or local_metadata)
        ws_metadata = workspace.get("rcsb_metadata") or workspace.get("local_metadata")
        if ws_metadata:
            if hasattr(ws_metadata, 'seqres_records'):
                metadata['seqres'] = ws_metadata.seqres_records
            if hasattr(ws_metadata, 'dbref_records'):
                # Convert dbref_records from list to dict keyed by chain_id
                dbref_list = ws_metadata.dbref_records
                if isinstance(dbref_list, list):
                    # Convert list to dict: {chain_id: dbref_info}
                    dbref_dict = {}
                    for dbref in dbref_list:
                        if isinstance(dbref, dict) and 'chain_id' in dbref:
                            dbref_dict[dbref['chain_id']] = dbref
                    metadata['dbref'] = dbref_dict
                elif isinstance(dbref_list, dict):
                    # Already a dict, use as-is
                    metadata['dbref'] = dbref_list
            if hasattr(ws_metadata, 'missing_res_records'):
                metadata['missing_res_records'] = ws_metadata.missing_res_records

        # Also try to get missing residues from detection results (if we ran detection)
        if hasattr(self, 'results') and self.results and 'missing_residues' in self.results:
            # Use the first available detection method's results
            for method, chain_results in self.results['missing_residues'].items():
                if chain_results:
                    # Convert from ResidueIdentity objects to simple dict format
                    missing_dict = {}
                    for chain_id, residues in chain_results.items():
                        missing_dict[chain_id] = []
                        for res in residues:
                            if hasattr(res, 'res_num'):
                                missing_dict[chain_id].append({
                                    'residue_number': res.res_num,
                                    'residue_name': res.res_name if hasattr(res, 'res_name') else 'UNK'
                                })
                    metadata['missing_res_records'] = missing_dict
                    break  # Use first available method

        # If not available, would need to parse from PDB file
        # For now, return what we have
        return metadata
    
    def _automated_repair(self, workspace: Any):
        """Automated repair for batch mode"""
        # Simple auto-repair strategy:
        # - Fill single gaps
        # - Skip multi-residue gaps
        # - Select highest occupancy for alternates
        
        structure, _, _ = self._get_priority_structure(workspace)
        if not structure:
            return
        
        # Create simplified plan
        plan = RepairPlan()
        
        # Extract single-residue gaps only
        missing = self._extract_missing_residues(self.results)
        for chain_id, residues in missing.items():
            sorted_res = sorted(residues, key=lambda r: r.res_num)
            
            # Group and take single-residue gaps
            current = [sorted_res[0]] if sorted_res else []
            for res in sorted_res[1:]:
                if res.res_num == current[-1].res_num + 1:
                    current.append(res)
                else:
                    if len(current) == 1:
                        seg = MissingSegment(
                            chain_id=chain_id,
                            residues=current,
                            start_num=current[0].res_num,
                            end_num=current[0].res_num
                        )
                        plan.segments_to_fill.append(seg)
                    current = [res]
            
            if current and len(current) == 1:
                seg = MissingSegment(
                    chain_id=chain_id,
                    residues=current,
                    start_num=current[0].res_num,
                    end_num=current[0].res_num
                )
                plan.segments_to_fill.append(seg)
        
        # Execute if there's anything to do
        if plan.segments_to_fill:
            self._execute_repair_plan(plan, structure, workspace)
    
    def _extract_missing_residues(self, results: Dict) -> Dict[str, List[ResidueIdentity]]:
        """Extract missing residues from results"""
        missing = defaultdict(list)
        
        for method in ['remark_465', 'seqres_comparison', 'fasta_comparison', 'sequence_gap']:
            if method in results.get('missing_residues', {}):
                for chain_id, residues in results['missing_residues'][method].items():
                    missing[chain_id].extend(residues)
                break
        
        return dict(missing)
    
    def get_workspace_display(self, workspace_key: str, value: Any, console: Console) -> bool:
        """Custom workspace display for completeness results"""
        if workspace_key == "completeness_results" and isinstance(value, dict):
            # Create summary table
            table = Table(title="Structure Completeness Summary")
            table.add_column("Category", style="cyan")
            table.add_column("Method", style="yellow")
            table.add_column("Issues", style="red")
            table.add_column("Chains", style="magenta")
            
            has_issues = False
            
            # Missing residues
            for method, chain_results in value.get('missing_residues', {}).items():
                total = sum(len(res) for res in chain_results.values())
                if total > 0:
                    has_issues = True
                    chains = ', '.join(chain_results.keys())
                    table.add_row("Missing Residues", method, str(total), chains)
            
            # Missing atoms
            for method, chain_results in value.get('missing_atoms', {}).items():
                total = sum(len(atoms) for atoms in chain_results.values())
                if total > 0:
                    has_issues = True
                    chains = ', '.join(chain_results.keys())
                    table.add_row("Missing Atoms", method, str(total), chains)
            
            # Alternate locations
            for method, chain_results in value.get('alternate_locations', {}).items():
                total = sum(len(res) for res in chain_results.values())
                if total > 0:
                    has_issues = True
                    chains = ', '.join(chain_results.keys())
                    table.add_row("Alternate Locations", method, str(total), chains)
            
            if has_issues:
                console.print(table)
            else:
                console.print("[green]✓ No structural issues detected[/green]")
            
            return True
        
        return False


# ============================================================================
# CLI ENTRY POINT
# ============================================================================

def main():
    """CLI entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze structure completeness")
    parser.add_argument("--pdbid", help="PDB ID to download")
    parser.add_argument("--pdbfile", help="Local PDB file")
    parser.add_argument("--output", help="Output JSON file")
    parser.add_argument("--verbose", "-v", action="store_true")
    
    args = parser.parse_args()
    
    console = Console()
    
    # Get PDB file
    if args.pdbid:
        # Download
        pdb_id = args.pdbid.upper()
        filename = f"{pdb_id}.pdb"
        
        try:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            with urllib.request.urlopen(url) as response:
                with open(filename, 'wb') as f:
                    f.write(response.read())
            console.print(f"[green]Downloaded {filename}[/green]")
        except Exception as e:
            console.print(f"[red]Download failed: {e}[/red]")
            return
    elif args.pdbfile:
        filename = args.pdbfile
    else:
        console.print("[red]Specify --pdbid or --pdbfile[/red]")
        return
    
    # Parse structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", filename)
    
    # Run analysis
    analyzer = StructureAnalyzer(structure, console=console)
    results = analyzer.detect_all_issues()
    
    # Display
    analyzer.display_summary_report()
    
    # Export
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        console.print(f"[green]Results exported to {args.output}[/green]")


if __name__ == "__main__":
    main()