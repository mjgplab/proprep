"""
Terminal Residue Classifier

Multi-method terminal residue detection for accurate AMBER force field parameterization.
Critical for proper atom type assignment (NALA vs ALA vs CALA).
"""

from pathlib import Path
from typing import Dict, Set, Optional, Tuple, List
from enum import Enum
from Bio.PDB import PDBParser, Structure
from rich.console import Console

from proprep.structure_prep.comprehensive_redox_detector import RedoxSite


class TerminalType(Enum):
    """Terminal classification for residues."""
    NTERM = "NTERM"       # N-terminal residue
    CTERM = "CTERM"       # C-terminal residue
    INTERNAL = "INTERNAL" # Internal residue in chain
    ISOLATED = "ISOLATED" # Single residue (not part of chain)


# Standard amino acid residues
STANDARD_AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # Common variants
    'HID', 'HIE', 'HIP', 'CYX', 'CYM'
}


class TerminalClassifier:
    """
    Multi-method terminal detection with 5 independent checks.

    Methods:
    1. Check for complete backbone atoms (N, CA, C, O)
    2. Check terminal atom markers (OXT for C-term, H1/H2/H3 for N-term)
    3. Check peptide bond connectivity (C-N bonds to adjacent residues)
    4. Check residue position in chain (first/last/middle)
    5. Check distance-based peptide bonds (C-N < 1.8 Å)

    All methods vote, results are synthesized with confidence scoring.
    """

    def __init__(self, pdb_file: str, redox_site: RedoxSite, console: Optional[Console] = None):
        """
        Initialize classifier.

        Args:
            pdb_file: Path to PDB file
            redox_site: RedoxSite object with atoms and bonds
            console: Rich console for output
        """
        self.pdb_file = Path(pdb_file)
        self.redox_site = redox_site
        self.console = console or Console()
        self.structure = None
        self.residue_map = {}  # (chain, resid) -> BioPython Residue
        self._parse_structure()

    def _parse_structure(self):
        """Parse PDB structure using BioPython."""
        if not self.pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {self.pdb_file}")

        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('structure', str(self.pdb_file))

        # Build residue map for quick lookup
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    het_flag, resseq, icode = residue.get_id()
                    key = (chain.id, resseq)
                    self.residue_map[key] = residue

    def classify_residue(self, chain_id: str, resid: int, resname: str) -> Tuple[TerminalType, str, List[str]]:
        """
        Classify a single residue using multi-method detection.

        Args:
            chain_id: Chain identifier
            resid: Residue number
            resname: Residue name (3-letter code)

        Returns:
            Tuple of (TerminalType, confidence, evidence_list)
        """
        evidence = []
        votes = []

        # Method 1: Check if protein residue (has backbone atoms)
        is_protein = resname in STANDARD_AMINO_ACIDS
        if not is_protein:
            evidence.append(f"{resname} is not a standard amino acid, no backbone")
            evidence.append("Non-protein residue - classified as INTERNAL")
            return TerminalType.INTERNAL, "high", evidence

        # Get residue from PDB
        residue = self.residue_map.get((chain_id, resid))
        if not residue:
            evidence.append(f"Residue {chain_id}:{resid} not found in PDB")
            return TerminalType.INTERNAL, "low", evidence

        # Method 2: Check complete backbone atoms
        atom_names = {atom.get_id() for atom in residue.get_atoms()}
        backbone_complete = {'N', 'CA', 'C', 'O'}.issubset(atom_names)

        if backbone_complete:
            evidence.append(f"Has complete backbone atoms: {{'N', 'CA', 'C', 'O'}}")
        else:
            evidence.append(f"Missing backbone atoms - has: {atom_names & {'N', 'CA', 'C', 'O'}}")
            votes.append(TerminalType.ISOLATED)

        # Method 3: Check terminal atom markers
        has_oxt = 'OXT' in atom_names
        n_h_atoms = atom_names & {'H', 'H1', 'H2', 'H3', 'HN', 'HN1', 'HN2', 'HN3', 'HT1', 'HT2', 'HT3'}
        has_multiple_n_h = len(n_h_atoms) >= 2

        if has_oxt and has_multiple_n_h:
            evidence.append(f"Has both OXT and {len(n_h_atoms)} N-H atoms")
            votes.append(TerminalType.ISOLATED)
        elif has_oxt:
            evidence.append("Has OXT atom (C-terminal marker)")
            votes.append(TerminalType.CTERM)
        elif has_multiple_n_h:
            evidence.append(f"Has {len(n_h_atoms)} N-H atoms: {n_h_atoms} (N-terminal marker)")
            votes.append(TerminalType.NTERM)

        # Method 4: Check peptide bond connectivity using RedoxSite bonds
        peptide_vote, peptide_evidence = self._check_peptide_bonds_redox_site(chain_id, resid)
        if peptide_vote:
            votes.append(peptide_vote)
            evidence.extend(peptide_evidence)

        # Method 5: Check distance-based peptide bonds in PDB
        distance_vote, distance_evidence = self._check_distance_based_peptide_bonds(chain_id, resid)
        if distance_vote:
            votes.append(distance_vote)
            evidence.extend(distance_evidence)

        # Method 6: Check residue position in chain
        position_vote, position_evidence = self._check_chain_position(chain_id, resid)
        if position_vote:
            votes.append(position_vote)
            evidence.extend(position_evidence)

        # Synthesize votes
        final_type, confidence = self._synthesize_votes(votes)
        return final_type, confidence, evidence

    def _check_peptide_bonds_redox_site(self, chain_id: str, resid: int) -> Tuple[Optional[TerminalType], List[str]]:
        """
        Check peptide bonds using RedoxSite bond information.

        Returns:
            (terminal_type_vote, evidence_list)
        """
        evidence = []

        # Get atoms from RedoxSite
        atoms = [a for a in self.redox_site.atoms
                if a.chain == chain_id and a.resid == resid]

        if not atoms:
            return None, ["Residue not in RedoxSite"]

        # Find C and N coordinates
        c_coords = None
        n_coords = None
        for atom in atoms:
            if atom.atom_name == 'C':
                c_coords = atom.coords
            elif atom.atom_name == 'N':
                n_coords = atom.coords

        if not c_coords or not n_coords:
            return None, ["Missing C or N backbone atoms"]

        # Check bonds
        has_prev_bond = False  # N to prev.C
        has_next_bond = False  # C to next.N

        for bond in self.redox_site.bonds:
            # Skip coordinate bonds
            if bond.chemical_type not in ['covalent', 'unknown']:
                continue

            atom1_info = bond.atom1_residue_info
            atom2_info = bond.atom2_residue_info

            # Check N bonded to prev.C
            if bond.atom1_coords == n_coords or bond.atom2_coords == n_coords:
                other_info = atom2_info if bond.atom1_coords == n_coords else atom1_info
                if (other_info.get('atom_name') == 'C' and
                    other_info.get('chain') == chain_id and
                    other_info.get('resid') == resid - 1):
                    has_prev_bond = True
                    evidence.append(f"N bonded to C of residue {resid-1} (peptide bond)")

            # Check C bonded to next.N
            if bond.atom1_coords == c_coords or bond.atom2_coords == c_coords:
                other_info = atom2_info if bond.atom1_coords == c_coords else atom1_info
                if (other_info.get('atom_name') == 'N' and
                    other_info.get('chain') == chain_id and
                    other_info.get('resid') == resid + 1):
                    has_next_bond = True
                    evidence.append(f"C bonded to N of residue {resid+1} (peptide bond)")

        # Classify
        if not has_prev_bond and not has_next_bond:
            evidence.append("No peptide bonds to adjacent residues")
            return TerminalType.ISOLATED, evidence
        elif not has_prev_bond:
            return TerminalType.NTERM, evidence
        elif not has_next_bond:
            return TerminalType.CTERM, evidence
        else:
            return TerminalType.INTERNAL, evidence

    def _check_distance_based_peptide_bonds(self, chain_id: str, resid: int) -> Tuple[Optional[TerminalType], List[str]]:
        """
        Check peptide bonds using C-N distance in PDB (< 1.8 Å).

        Returns:
            (terminal_type_vote, evidence_list)
        """
        evidence = []

        current_res = self.residue_map.get((chain_id, resid))
        if not current_res:
            return None, []

        # Check previous peptide bond
        has_prev = self._has_distance_bond_to_previous(chain_id, resid)
        if has_prev:
            evidence.append(f"C-N distance to residue {resid-1} < 1.8 Å (peptide bond)")

        # Check next peptide bond
        has_next = self._has_distance_bond_to_next(chain_id, resid)
        if has_next:
            evidence.append(f"C-N distance to residue {resid+1} < 1.8 Å (peptide bond)")

        # Classify
        if not has_prev and not has_next:
            return TerminalType.ISOLATED, evidence
        elif not has_prev:
            return TerminalType.NTERM, evidence
        elif not has_next:
            return TerminalType.CTERM, evidence
        else:
            return TerminalType.INTERNAL, evidence

    def _has_distance_bond_to_previous(self, chain_id: str, resid: int) -> bool:
        """Check if C-N distance to previous residue < 1.8 Å."""
        current_res = self.residue_map.get((chain_id, resid))
        prev_res = self.residue_map.get((chain_id, resid - 1))

        if not current_res or not prev_res:
            return False

        # Check both are standard amino acids
        if (current_res.get_resname() not in STANDARD_AMINO_ACIDS or
            prev_res.get_resname() not in STANDARD_AMINO_ACIDS):
            return False

        # Check atoms exist
        if 'N' not in current_res or 'C' not in prev_res:
            return False

        try:
            n_atom = current_res['N']
            c_atom = prev_res['C']
            distance = n_atom - c_atom  # BioPython distance

            return distance < 1.8
        except:
            return False

    def _has_distance_bond_to_next(self, chain_id: str, resid: int) -> bool:
        """Check if C-N distance to next residue < 1.8 Å."""
        current_res = self.residue_map.get((chain_id, resid))
        next_res = self.residue_map.get((chain_id, resid + 1))

        if not current_res or not next_res:
            return False

        # Check both are standard amino acids
        if (current_res.get_resname() not in STANDARD_AMINO_ACIDS or
            next_res.get_resname() not in STANDARD_AMINO_ACIDS):
            return False

        # Check atoms exist
        if 'C' not in current_res or 'N' not in next_res:
            return False

        try:
            c_atom = current_res['C']
            n_atom = next_res['N']
            distance = c_atom - n_atom  # BioPython distance

            return distance < 1.8
        except:
            return False

    def _check_chain_position(self, chain_id: str, resid: int) -> Tuple[Optional[TerminalType], List[str]]:
        """
        Check residue position in chain (first/last/middle).

        Returns:
            (terminal_type_vote, evidence_list)
        """
        evidence = []

        try:
            # Get chain
            chain = self.structure[0][chain_id]

            # Get all protein residues (standard residues with backbone)
            protein_resids = []
            for residue in chain.get_residues():
                if residue.get_id()[0] == ' ':  # ATOM record
                    atom_names = {a.get_id() for a in residue.get_atoms()}
                    if {'N', 'CA', 'C'}.issubset(atom_names):
                        if residue.get_resname() in STANDARD_AMINO_ACIDS:
                            protein_resids.append(residue.get_id()[1])

            if not protein_resids:
                return None, ["No protein residues in chain"]

            protein_resids.sort()
            min_resid = min(protein_resids)
            max_resid = max(protein_resids)

            # Check position
            if resid == min_resid and resid == max_resid:
                evidence.append(f"Residue {resid} is the ONLY protein residue in chain")
                return TerminalType.ISOLATED, evidence
            elif resid == min_resid:
                evidence.append(f"Residue {resid} is FIRST in chain ({min_resid}-{max_resid})")
                return TerminalType.NTERM, evidence
            elif resid == max_resid:
                evidence.append(f"Residue {resid} is LAST in chain ({min_resid}-{max_resid})")
                return TerminalType.CTERM, evidence
            else:
                evidence.append(f"Residue {resid} is in middle of chain ({min_resid}-{max_resid})")
                return TerminalType.INTERNAL, evidence

        except Exception as e:
            return None, [f"Error checking chain position: {e}"]

    def _synthesize_votes(self, votes: List[TerminalType]) -> Tuple[TerminalType, str]:
        """
        Synthesize votes into final classification.

        Returns:
            (final_type, confidence)
        """
        if not votes:
            return TerminalType.INTERNAL, "low"

        # Count votes
        vote_counts = {}
        for vote in votes:
            vote_counts[vote] = vote_counts.get(vote, 0) + 1

        # Unanimous
        if len(vote_counts) == 1:
            return votes[0], "high"

        # Majority
        max_votes = max(vote_counts.values())
        winners = [vtype for vtype, count in vote_counts.items() if count == max_votes]

        if len(winners) == 1:
            confidence = "high" if max_votes >= len(votes) * 0.7 else "medium"
            return winners[0], confidence

        # Tie - use conservative (INTERNAL)
        if TerminalType.INTERNAL in winners:
            return TerminalType.INTERNAL, "low"
        elif TerminalType.NTERM in winners:
            return TerminalType.NTERM, "low"
        elif TerminalType.CTERM in winners:
            return TerminalType.CTERM, "low"
        else:
            return TerminalType.ISOLATED, "low"

    def get_force_field_residue_name(self, resname: str, terminal_type: TerminalType) -> str:
        """
        Get force field residue name with terminal prefix.

        Examples:
        - ALA + NTERM -> NALA
        - ALA + CTERM -> CALA
        - ALA + INTERNAL -> ALA
        - CU + any -> CU (non-protein)

        Args:
            resname: Base residue name
            terminal_type: Terminal classification

        Returns:
            Force field residue name
        """
        # Non-protein residues don't get terminal prefixes
        if resname not in STANDARD_AMINO_ACIDS:
            return resname

        if terminal_type == TerminalType.NTERM:
            return f"N{resname}"
        elif terminal_type == TerminalType.CTERM:
            return f"C{resname}"
        else:
            # INTERNAL or ISOLATED
            return resname
