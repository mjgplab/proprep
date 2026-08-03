"""
Terminal Residue Classifier

Multi-method terminal residue detection for accurate AMBER atom typing.
Uses both RedoxSite object and original PDB file for maximum reliability.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set
from enum import Enum
import logging

try:
    from Bio.PDB import PDBParser, Chain, Residue
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

from proprep.structure_prep.comprehensive_redox_detector import (
    RedoxSite, RedoxSiteAtom, SimpleCCDQuerier
)


class TerminalType(Enum):
    """Residue terminal classification"""
    NTERM = "nterm"      # N-terminal (first in chain)
    CTERM = "cterm"      # C-terminal (last in chain)
    INTERNAL = "internal"  # Middle of chain
    ISOLATED = "isolated"  # Single residue (both terminals)


@dataclass
class TerminalDetectionResult:
    """Result of terminal detection analysis"""
    residue_key: Tuple[str, int, str]  # (chain, resid, insertion)
    terminal_type: TerminalType
    confidence: str  # "high", "medium", "low"
    evidence: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


class TerminalClassifier:
    """
    Multi-method terminal detection with user confirmation.

    Uses 6 independent detection methods (priority order):
    1. Protein vs non-protein classification
    2. Backbone atom markers in PDB (OXT, H1/H2/H3)
    3. Peptide bond connectivity with distance tolerance
    4. Sequential atom ID check (C and N are sequential)
    5. Residue numbering in chain (min/max)
    6. Terminal atom markers in RedoxSite (lowest priority)

    All methods vote, results are synthesized with confidence scoring,
    and user confirmation is requested for uncertain cases.
    """

    def __init__(self, site: RedoxSite, pdb_file: str, logger=None, processor=None):
        """
        Initialize terminal classifier.

        Args:
            site: RedoxSite object containing atoms and bonds
            pdb_file: Path to original PDB file
            logger: Optional logger instance
            processor: Optional ProPrep processor for session recording context
        """
        self.site = site
        self.pdb_file = pdb_file
        self.processor = processor
        self.logger = logger or logging.getLogger(__name__)

        # Parse PDB structure
        if BIOPYTHON_AVAILABLE:
            parser = PDBParser(QUIET=True)
            try:
                self.structure = parser.get_structure('structure', pdb_file)
            except Exception as e:
                self.logger.warning(f"Could not parse PDB file: {e}")
                self.structure = None
        else:
            self.logger.warning("BioPython not available - some detection methods disabled")
            self.structure = None

        # Detection results
        self.detection_results: Dict[Tuple, TerminalDetectionResult] = {}

        # CCD querier for residue classification
        self.ccd_querier = SimpleCCDQuerier()

    def classify_all_residues(self,
                             require_confirmation: bool = True) -> Dict[Tuple, TerminalType]:
        """
        Classify all protein residues in the PDB structure.

        For ONIOM atom typing, we need to classify ALL protein residues (not just
        RedoxSite) to assign correct AMBER terminal vs internal atom types.

        Args:
            require_confirmation: If True, prompt user to confirm uncertain cases

        Returns:
            Dictionary mapping residue keys to terminal types
        """
        self.logger.debug("Starting terminal residue classification for full protein...")

        if not self.structure:
            self.logger.warning("No PDB structure available - falling back to RedoxSite residues only")
            # Fallback: use RedoxSite residues
            for res_key in self.site.residue_groups.keys():
                result = self._classify_single_residue(res_key)
                self.detection_results[res_key] = result
        else:
            # Scan entire PDB structure for protein residues
            for model in self.structure:
                for chain in model:
                    for residue in chain.get_residues():
                        # Only process standard protein residues (ATOM records, not HETATM)
                        if residue.get_id()[0] != ' ':
                            continue

                        # Check if it has backbone atoms (protein check)
                        atom_names = {a.get_id() for a in residue.get_atoms()}
                        if not {'N', 'CA', 'C'}.issubset(atom_names):
                            continue

                        # Build residue key
                        chain_id = chain.get_id()
                        resid = residue.get_id()[1]
                        insertion = residue.get_id()[2].strip()
                        res_key = (chain_id, resid, insertion)

                        # Classify this residue
                        result = self._classify_single_residue(res_key)
                        self.detection_results[res_key] = result

        # User confirmation if requested
        if require_confirmation:
            self._confirm_with_user()

        # Return final classifications
        final_classifications = {
            key: result.terminal_type
            for key, result in self.detection_results.items()
        }

        self.logger.debug(f"Classified {len(final_classifications)} residues")
        return final_classifications

    def _classify_single_residue(self,
                                res_key: Tuple[str, int, str]) -> TerminalDetectionResult:
        """
        Apply all detection methods to a single residue.
        Synthesize results into a single classification.

        Args:
            res_key: (chain, resid, insertion_code)

        Returns:
            TerminalDetectionResult with classification and evidence
        """
        chain_id, resid, insertion = res_key

        evidence = []
        warnings = []
        votes = []  # Each method votes for a terminal type

        # Method 1: Check if it's a protein residue
        is_protein, protein_evidence = self._check_is_protein(res_key)
        evidence.extend(protein_evidence)

        if not is_protein:
            # Non-protein (metal, cofactor, etc.) → always INTERNAL
            return TerminalDetectionResult(
                residue_key=res_key,
                terminal_type=TerminalType.INTERNAL,
                confidence="high",
                evidence=evidence + ["Non-protein residue - classified as INTERNAL"],
                warnings=[]
            )

        # Method 2: Check backbone atoms in PDB
        if self.structure:
            backbone_vote, backbone_evidence = self._check_backbone_atoms_pdb(res_key)
            if backbone_vote is not None:
                votes.append(backbone_vote)
            evidence.extend(backbone_evidence)

        # Method 3: Check peptide bond connectivity with distance tolerance
        peptide_vote, peptide_evidence = self._check_peptide_bonds(res_key)
        votes.append(peptide_vote)
        evidence.extend(peptide_evidence)

        # Method 4: Check sequential atom IDs (C and N sequential)
        sequential_vote, sequential_evidence = self._check_sequential_atom_ids(res_key)
        if sequential_vote is not None:
            votes.append(sequential_vote)
        evidence.extend(sequential_evidence)

        # Method 5: Check residue numbering in chain
        if self.structure:
            numbering_vote, numbering_evidence = self._check_residue_numbering(res_key)
            if numbering_vote is not None:
                votes.append(numbering_vote)
            evidence.extend(numbering_evidence)

        # Method 6: Check terminal atom markers in RedoxSite (lowest priority)
        marker_vote, marker_evidence, marker_warnings = self._check_terminal_markers(res_key)
        if marker_vote is not None:
            votes.append(marker_vote)
        evidence.extend(marker_evidence)
        warnings.extend(marker_warnings)

        # Synthesize votes
        final_type, confidence = self._synthesize_votes(votes)

        return TerminalDetectionResult(
            residue_key=res_key,
            terminal_type=final_type,
            confidence=confidence,
            evidence=evidence,
            warnings=warnings
        )

    def _check_is_protein(self, res_key: Tuple) -> Tuple[bool, List[str]]:
        """
        Method 1: Determine if residue is protein.

        Uses PDB structure when available (for full protein classification),
        falls back to RedoxSite if not.

        Returns:
            (is_protein, evidence_list)
        """
        chain_id, resid, insertion = res_key
        evidence = []

        # Try PDB structure first (for full protein classification)
        if self.structure:
            try:
                chain = self.structure[0][chain_id]
                residue_id = (' ', resid, insertion if insertion else ' ')
                residue = chain[residue_id]

                resname = residue.get_resname()
                atom_names = {a.get_id() for a in residue.get_atoms()}

                # Check for backbone atoms
                backbone_atoms = {'N', 'CA', 'C', 'O'}
                has_backbone = backbone_atoms.issubset(atom_names)

                if has_backbone:
                    evidence.append(f"Protein residue with backbone: {resname}")
                    return True, evidence

                # Check against standard amino acids
                is_std_aa = (resname in self.ccd_querier.STANDARD_AMINO_ACIDS or
                            resname in self.ccd_querier.REDOX_ACTIVE_AMINO_ACIDS)

                if is_std_aa:
                    evidence.append(f"{resname} is standard amino acid (missing backbone)")
                    return True, evidence

                # Not protein
                evidence.append(f"{resname} is not a protein residue")
                return False, evidence

            except KeyError:
                # Residue not in PDB, try RedoxSite
                pass
            except Exception as e:
                evidence.append(f"Error reading PDB: {e}")

        # Fallback: Check RedoxSite
        atoms = self.site.get_atoms_by_residue(chain_id, resid, insertion)
        if not atoms:
            evidence.append("Not found in PDB or RedoxSite")
            return False, evidence

        resname = atoms[0].resname
        atom_names = {a.atom_name for a in atoms}

        # Check for backbone atoms
        backbone_atoms = {'N', 'CA', 'C', 'O'}
        has_backbone = backbone_atoms.issubset(atom_names)

        if has_backbone:
            evidence.append(f"Protein residue with backbone: {resname}")
            return True, evidence

        # Check against standard/redox amino acids
        is_std_aa = (resname in self.ccd_querier.STANDARD_AMINO_ACIDS or
                    resname in self.ccd_querier.REDOX_ACTIVE_AMINO_ACIDS)

        if is_std_aa:
            evidence.append(f"{resname} is standard amino acid")
            return True, evidence

        # Not protein
        evidence.append(f"{resname} is not a protein residue")
        return False, evidence

    def _check_backbone_atoms_pdb(self, res_key: Tuple) -> Tuple[Optional[TerminalType], List[str]]:
        """
        Method 2: Check backbone atoms in original PDB.

        Returns:
            (terminal_type_vote, evidence_list)
        """
        chain_id, resid, insertion = res_key
        evidence = []

        if not self.structure:
            return None, ["PDB structure not available"]

        try:
            # Get residue from BioPython structure
            chain = self.structure[0][chain_id]
            residue_id = (' ', resid, insertion if insertion else ' ')
            residue = chain[residue_id]

            atom_names = {atom.get_id() for atom in residue.get_atoms()}

            # Check for terminal markers in PDB
            has_oxt = 'OXT' in atom_names
            n_h_atoms = {a for a in atom_names if a in {'H', 'H1', 'H2', 'H3', 'HN', 'HN1', 'HN2', 'HN3'}}
            has_multiple_n_h = len(n_h_atoms) >= 2

            if has_oxt and has_multiple_n_h:
                evidence.append(f"PDB has both OXT and {len(n_h_atoms)} N-H atoms")
                return TerminalType.ISOLATED, evidence
            elif has_oxt:
                evidence.append("PDB has OXT atom (C-terminal marker)")
                return TerminalType.CTERM, evidence
            elif has_multiple_n_h:
                evidence.append(f"PDB has {len(n_h_atoms)} N-H atoms: {n_h_atoms} (N-terminal marker)")
                return TerminalType.NTERM, evidence
            else:
                evidence.append("PDB has no clear terminal markers in backbone")
                return None, evidence

        except KeyError:
            evidence.append(f"Residue {chain_id}:{resid} not found in PDB structure")
            return None, evidence
        except Exception as e:
            evidence.append(f"Error reading PDB: {e}")
            return None, evidence

    def _check_peptide_bonds(self, res_key: Tuple) -> Tuple[TerminalType, List[str]]:
        """
        Method 3: Check peptide bond connectivity with distance tolerance.

        Uses PDB structure to find adjacent residues (not RedoxSite bonds, since
        RedoxSite may not contain the full protein chain). Checks for peptide bonds
        based on C-N distance (typical ~1.33 Å).

        Returns:
            (terminal_type_vote, evidence_list)
        """
        chain_id, resid, insertion = res_key
        evidence = []

        # Peptide bond distance tolerance (Å)
        PEPTIDE_BOND_DISTANCE_MIN = 1.20  # Minimum C-N distance
        PEPTIDE_BOND_DISTANCE_MAX = 1.50  # Maximum C-N distance

        if not self.structure:
            evidence.append("PDB structure not available for peptide bond check")
            return TerminalType.INTERNAL, evidence

        try:
            # Get this residue from PDB
            chain = self.structure[0][chain_id]
            residue_id = (' ', resid, insertion if insertion else ' ')
            residue = chain[residue_id]

            # Get C and N atoms from this residue
            c_atom = None
            n_atom = None
            for atom in residue.get_atoms():
                if atom.get_id() == 'C':
                    c_atom = atom
                elif atom.get_id() == 'N':
                    n_atom = atom

            if not c_atom or not n_atom:
                evidence.append("Missing C or N backbone atoms in PDB")
                return TerminalType.INTERNAL, evidence

            c_coords = c_atom.get_coord()
            n_coords = n_atom.get_coord()

            # Check for peptide bonds to adjacent residues
            has_prev_bond = False  # N to prev.C
            has_next_bond = False  # C to next.N

            # Search all residues in this chain for potential peptide bonds
            for other_residue in chain.get_residues():
                # Skip this residue
                other_resid = other_residue.get_id()[1]
                if other_resid == resid:
                    continue

                # Get backbone atoms from other residue
                other_c = None
                other_n = None
                for atom in other_residue.get_atoms():
                    if atom.get_id() == 'C':
                        other_c = atom
                    elif atom.get_id() == 'N':
                        other_n = atom

                # Check N-terminal peptide bond (this.N to other.C)
                if other_c and n_atom:
                    other_c_coords = other_c.get_coord()
                    distance = ((n_coords[0] - other_c_coords[0])**2 +
                               (n_coords[1] - other_c_coords[1])**2 +
                               (n_coords[2] - other_c_coords[2])**2)**0.5

                    if PEPTIDE_BOND_DISTANCE_MIN <= distance <= PEPTIDE_BOND_DISTANCE_MAX:
                        has_prev_bond = True
                        evidence.append(f"N→C of res {other_resid} (d={distance:.2f}Å, peptide)")

                # Check C-terminal peptide bond (this.C to other.N)
                if other_n and c_atom:
                    other_n_coords = other_n.get_coord()
                    distance = ((c_coords[0] - other_n_coords[0])**2 +
                               (c_coords[1] - other_n_coords[1])**2 +
                               (c_coords[2] - other_n_coords[2])**2)**0.5

                    if PEPTIDE_BOND_DISTANCE_MIN <= distance <= PEPTIDE_BOND_DISTANCE_MAX:
                        has_next_bond = True
                        evidence.append(f"C→N of res {other_resid} (d={distance:.2f}Å, peptide)")

            # Classify based on connectivity
            if not has_prev_bond and not has_next_bond:
                evidence.append("No peptide bonds found in PDB structure")
                return TerminalType.ISOLATED, evidence
            elif not has_prev_bond:
                evidence.append("No bond to previous residue (N-terminal)")
                return TerminalType.NTERM, evidence
            elif not has_next_bond:
                evidence.append("No bond to next residue (C-terminal)")
                return TerminalType.CTERM, evidence
            else:
                evidence.append("Bonded to both previous and next residues")
                return TerminalType.INTERNAL, evidence

        except KeyError:
            evidence.append(f"Residue {chain_id}:{resid} not found in PDB structure")
            return TerminalType.INTERNAL, evidence
        except Exception as e:
            evidence.append(f"Error checking peptide bonds: {e}")
            return TerminalType.INTERNAL, evidence

    def _check_sequential_atom_ids(self, res_key: Tuple) -> Tuple[Optional[TerminalType], List[str]]:
        """
        Method 4: Check sequential atom IDs/serial numbers.

        In well-formed PDB files, the C atom of one residue and the N atom
        of the next residue should have sequential serial numbers.

        Returns:
            (terminal_type_vote, evidence_list)
        """
        chain_id, resid, insertion = res_key
        evidence = []

        # Get atoms from this residue
        atoms = self.site.get_atoms_by_residue(chain_id, resid, insertion)

        # Find C and N atoms with their serial numbers
        c_atom = None
        n_atom = None
        for atom in atoms:
            if atom.atom_name == 'C':
                c_atom = atom
            elif atom.atom_name == 'N':
                n_atom = atom

        if not c_atom or not n_atom:
            evidence.append("Missing C or N atoms for serial number check")
            return None, evidence

        # Get serial numbers from properties
        c_serial = c_atom.properties.get('serial_number')
        n_serial = n_atom.properties.get('serial_number')

        if c_serial is None or n_serial is None:
            evidence.append("Serial numbers not available in RedoxSite")
            return None, evidence

        # Check for sequential patterns
        has_prev_sequential = False
        has_next_sequential = False

        # Get all atoms from RedoxSite to find adjacent residue backbone atoms
        for other_atom in self.site.atoms:
            other_serial = other_atom.properties.get('serial_number')
            if other_serial is None:
                continue

            # Check if this is a backbone C or N from same chain
            if other_atom.chain != chain_id:
                continue

            # Check for N-terminal: another C has serial = N_serial - 1
            if (other_atom.atom_name == 'C' and
                other_serial == n_serial - 1):
                has_prev_sequential = True
                evidence.append(f"Serial {n_serial}(N) follows {other_serial}(C) of res {other_atom.resid}")

            # Check for C-terminal: another N has serial = C_serial + 1
            if (other_atom.atom_name == 'N' and
                other_serial == c_serial + 1):
                has_next_sequential = True
                evidence.append(f"Serial {other_serial}(N) follows {c_serial}(C) of res {resid}")

        # Only vote if we found evidence
        if not has_prev_sequential and not has_next_sequential:
            return None, evidence
        elif not has_prev_sequential:
            evidence.append("No sequential C before N (likely N-terminal)")
            return TerminalType.NTERM, evidence
        elif not has_next_sequential:
            evidence.append("No sequential N after C (likely C-terminal)")
            return TerminalType.CTERM, evidence
        else:
            evidence.append("Sequential atom IDs suggest internal residue")
            return TerminalType.INTERNAL, evidence

    def _check_residue_numbering(self, res_key: Tuple) -> Tuple[Optional[TerminalType], List[str]]:
        """
        Method 5: Check residue numbering within chain.

        Returns:
            (terminal_type_vote, evidence_list)
        """
        chain_id, resid, insertion = res_key
        evidence = []

        if not self.structure:
            return None, ["PDB structure not available"]

        try:
            chain = self.structure[0][chain_id]

            # Get all protein residues in this chain (ATOM records only)
            protein_resids = []
            for residue in chain.get_residues():
                # Filter for ATOM records (not HETATM)
                if residue.get_id()[0] == ' ':  # Standard residue
                    # Check if it has backbone atoms
                    atom_names = {a.get_id() for a in residue.get_atoms()}
                    if {'N', 'CA', 'C'}.issubset(atom_names):
                        protein_resids.append(residue.get_id()[1])

            if not protein_resids:
                evidence.append("No protein residues found in chain")
                return None, evidence

            protein_resids.sort()
            min_resid = min(protein_resids)
            max_resid = max(protein_resids)

            evidence.append(f"Chain {chain_id} protein residues: {min_resid}-{max_resid}")

            # Check position in chain
            if resid == min_resid and resid == max_resid:
                evidence.append(f"Residue {resid} is the ONLY residue in chain")
                return TerminalType.ISOLATED, evidence
            elif resid == min_resid:
                evidence.append(f"Residue {resid} is FIRST in chain")
                return TerminalType.NTERM, evidence
            elif resid == max_resid:
                evidence.append(f"Residue {resid} is LAST in chain")
                return TerminalType.CTERM, evidence
            else:
                evidence.append(f"Residue {resid} is in middle of chain ({min_resid}-{max_resid})")
                return TerminalType.INTERNAL, evidence

        except Exception as e:
            evidence.append(f"Error checking residue numbering: {e}")
            return None, evidence

    def _check_terminal_markers(self, res_key: Tuple) -> Tuple[Optional[TerminalType], List[str], List[str]]:
        """
        Method 6: Check terminal atom markers in RedoxSite (least reliable).

        Returns:
            (terminal_type_vote, evidence_list, warnings_list)
        """
        chain_id, resid, insertion = res_key
        evidence = []
        warnings = []

        atoms = self.site.get_atoms_by_residue(chain_id, resid, insertion)
        atom_names = {a.atom_name for a in atoms}

        # C-terminal markers
        has_oxt = 'OXT' in atom_names

        # N-terminal markers (various naming conventions)
        n_h_atoms = {a for a in atom_names if a in {'H', 'H1', 'H2', 'H3', 'HN', 'HN1', 'HN2', 'HN3', 'HT1', 'HT2', 'HT3'}}
        has_multiple_n_h = len(n_h_atoms) >= 2

        if not has_oxt and not has_multiple_n_h:
            warnings.append("No terminal markers found (OXT often missing, H atoms often absent)")
            return None, evidence, warnings

        if has_oxt:
            evidence.append("Has OXT atom (C-terminal marker)")

        if has_multiple_n_h:
            evidence.append(f"Has {len(n_h_atoms)} N-H atoms: {n_h_atoms} (N-terminal marker)")
            warnings.append("N-H atoms may be missing if no hydrogens in PDB")

        if has_oxt and has_multiple_n_h:
            return TerminalType.ISOLATED, evidence, warnings
        elif has_oxt:
            return TerminalType.CTERM, evidence, warnings
        elif has_multiple_n_h:
            return TerminalType.NTERM, evidence, warnings
        else:
            return None, evidence, warnings

    def _synthesize_votes(self, votes: List[TerminalType]) -> Tuple[TerminalType, str]:
        """
        Synthesize multiple votes into final classification.

        Priority:
        1. Unanimous votes → high confidence
        2. Majority votes → medium confidence
        3. Tie/conflict → low confidence, use most conservative (INTERNAL)

        Returns:
            (final_terminal_type, confidence_level)
        """
        if not votes:
            return TerminalType.INTERNAL, "low"

        # Count votes
        vote_counts = {}
        for vote in votes:
            vote_counts[vote] = vote_counts.get(vote, 0) + 1

        # Check for unanimous
        if len(vote_counts) == 1:
            return votes[0], "high"

        # Find majority
        max_votes = max(vote_counts.values())
        winners = [vtype for vtype, count in vote_counts.items() if count == max_votes]

        if len(winners) == 1:
            # Clear majority
            confidence = "medium" if max_votes >= len(votes) * 0.6 else "low"
            return winners[0], confidence

        # Tie - use conservative choice
        # Priority: INTERNAL > NTERM/CTERM > ISOLATED
        if TerminalType.INTERNAL in winners:
            return TerminalType.INTERNAL, "low"
        elif TerminalType.NTERM in winners:
            return TerminalType.NTERM, "low"
        elif TerminalType.CTERM in winners:
            return TerminalType.CTERM, "low"
        else:
            return TerminalType.ISOLATED, "low"

    def _confirm_with_user(self):
        """
        Present terminal residue detection results to user.
        Only shows terminals (N-term, C-term, Isolated), not internal residues.
        Uses Rich for consistent ProPrep UX.
        """
        from rich.console import Console
        from rich.panel import Panel
        from rich.table import Table
        from proprep.utils.prompts import confirm_with_context

        console = Console()

        # Educational context
        console.print()
        console.print(Panel(
            "[bold]Why Terminal Detection Matters:[/bold]\n\n"
            "AMBER force fields use different atom types and charges for terminal residues:\n"
            "  • [cyan]N-terminal[/cyan]: NALA, NGLY, etc. (NH3+ group)\n"
            "  • [cyan]C-terminal[/cyan]: CALA, CGLY, etc. (COO- group)\n"
            "  • [cyan]Internal[/cyan]: ALA, GLY, etc. (peptide bond)\n\n"
            "Incorrect classification leads to wrong charges and failed simulations.",
            title="Terminal Residue Classification",
            border_style="cyan"
        ))

        # Filter for only terminal residues (not INTERNAL)
        terminal_results = [
            r for r in self.detection_results.values()
            if r.terminal_type != TerminalType.INTERNAL
        ]

        if not terminal_results:
            console.print("\n[green]✓ No terminal residues detected - all residues are INTERNAL[/green]")
            console.print("[grey50]  (This is expected for structures with capped termini or mid-chain extracts)[/grey50]\n")
            return

        # Organize by confidence
        high_conf = [r for r in terminal_results if r.confidence == "high"]
        med_conf = [r for r in terminal_results if r.confidence == "medium"]
        low_conf = [r for r in terminal_results if r.confidence == "low"]

        # Display high confidence terminals in a table
        if high_conf:
            console.print("\n[bold green]✓ Terminals Detected (High Confidence)[/bold green]\n")

            table = Table(show_header=True, header_style="bold cyan")
            table.add_column("Residue", style="yellow", width=12)
            table.add_column("Type", style="green", width=10)
            table.add_column("Evidence", style="grey50")

            for result in high_conf:
                chain, resid, ins = result.residue_key
                ins_str = ins if ins else ''
                term_type = result.terminal_type.value.upper()

                # Combine top 3 evidence items
                evidence_str = "\n".join(f"• {e}" for e in result.evidence[:3])

                table.add_row(
                    f"{chain}:{resid}{ins_str}",
                    term_type,
                    evidence_str
                )

            console.print(table)

        # Display medium confidence (quick review)
        if med_conf:
            console.print("\n[bold yellow]⚠ Terminals Detected (Medium Confidence - Please Review)[/bold yellow]\n")

            for result in med_conf:
                chain, resid, ins = result.residue_key
                ins_str = ins if ins else ''

                console.print(f"[yellow]{chain}:{resid}{ins_str} → {result.terminal_type.value.upper()}[/yellow]")
                for evidence in result.evidence:
                    console.print(f"  [grey50]• {evidence}[/grey50]")

                # Ask for confirmation
                if not confirm_with_context(
                    self.processor,
                    f"  Accept this classification?",
                    default=True,
                    module="Terminal Classifier",
                    description="Accept automatic terminal classification",
                ):
                    new_type = self._prompt_terminal_type()
                    result.terminal_type = new_type
                    result.evidence.append(f"User override: {new_type.value}")

        # Display low confidence (require review)
        if low_conf:
            console.print("\n[bold red]❌ Uncertain Terminals (Manual Review Required)[/bold red]\n")

            for result in low_conf:
                chain, resid, ins = result.residue_key
                ins_str = ins if ins else ''

                console.print(f"[red]{chain}:{resid}{ins_str} → {result.terminal_type.value.upper()} (uncertain)[/red]")
                console.print("[bold]Evidence:[/bold]")
                for evidence in result.evidence:
                    console.print(f"  [grey50]• {evidence}[/grey50]")

                if result.warnings:
                    console.print("[bold yellow]Warnings:[/bold yellow]")
                    for warning in result.warnings:
                        console.print(f"  [yellow]⚠ {warning}[/yellow]")

                # Require user input
                new_type = self._prompt_terminal_type(default=result.terminal_type)
                if new_type != result.terminal_type:
                    result.terminal_type = new_type
                    result.evidence.append(f"User override: {new_type.value}")
                    result.confidence = "high"
                else:
                    result.evidence.append("User confirmed auto-classification")
                    result.confidence = "high"

        # Summary
        total_residues = len(self.detection_results)
        terminal_count = len(terminal_results)
        internal_count = total_residues - terminal_count

        console.print()
        console.print(Panel(
            f"[bold]Classification Summary[/bold]\n\n"
            f"Total protein residues: {total_residues}\n"
            f"  • {internal_count} INTERNAL residues\n"
            f"  • {terminal_count} TERMINAL residues",
            border_style="green"
        ))
        console.print()

    def _prompt_terminal_type(self, default: TerminalType = None) -> TerminalType:
        """
        Prompt user to select terminal type.

        Args:
            default: Default terminal type if user presses Enter

        Returns:
            Selected TerminalType
        """
        options = {
            '1': TerminalType.NTERM,
            '2': TerminalType.CTERM,
            '3': TerminalType.INTERNAL,
            '4': TerminalType.ISOLATED
        }

        print("  Select terminal type:")
        print("    1. N-terminal")
        print("    2. C-terminal")
        print("    3. Internal (middle of chain)")
        print("    4. Isolated (single residue)")

        if default:
            print(f"    [Press Enter for default: {default.value}]")

        from proprep.utils.prompts import prompt_with_context
        while True:
            choice = prompt_with_context(
                self.processor, "Enter choice (1-4)",
                module="Terminal Classifier",
                description="Select terminal type",
                options_map={
                    "1": "N-terminal", "2": "C-terminal",
                    "3": "Internal (middle of chain)", "4": "Isolated (single residue)",
                },
            ).strip()
            if choice in options:
                return options[choice]
            elif choice == '' and default:
                return default
            else:
                print("  Invalid choice. Please enter 1, 2, 3, or 4.")
