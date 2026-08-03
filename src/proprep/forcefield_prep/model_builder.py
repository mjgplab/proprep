"""
Model Builder for QM Calculations

Builds small and large models from RedoxSite for QM calculations.
Based on MCPB.py algorithm from AmberTools.

Small Model: Sidechain-only model for bonded parameter calculation
Large Model: Full residue model for RESP charge fitting

Both models apply MCPB capping rules to prevent dangling bonds.
"""

from pathlib import Path
from typing import List, Set, Tuple, Dict, Optional
from dataclasses import dataclass
from enum import Enum

from Bio.PDB import PDBParser, Structure
from rich.console import Console

from proprep.structure_prep.comprehensive_redox_detector import RedoxSite, RedoxSiteAtom, METALS


class CapType(Enum):
    """Capping group types (MCPB convention)."""
    ACE = "ACE"      # CH3-CO- (N-terminal cap, keeps C=O)
    NME = "NME"      # -NH-CH3 (C-terminal cap, keeps N-H)
    GLY = "GLY"      # Full glycine (backbone only)
    FULL = "FULL"    # Full residue (backbone + sidechain) - for gap residues
    KNH = "KNH"      # Keep N,H backbone atoms
    KCO = "KCO"      # Keep C,O backbone atoms
    ANT = "ANT"      # CH3NH3+ (protonated methylamine)
    ACT = "ACT"      # CH3CO2- (acetate)
    NONE = "NONE"    # Sidechain only (CA→CH3)


@dataclass
class ModelResidue:
    """Residue in a QM model."""
    chain: str
    resid: int
    resname: str
    cap_type: CapType = CapType.NONE
    is_coordinating: bool = False
    coordinating_atoms: List[str] = None  # Atom names that coordinate to metal

    def __post_init__(self):
        if self.coordinating_atoms is None:
            self.coordinating_atoms = []

    def __repr__(self):
        coord_str = f" (coords via {', '.join(self.coordinating_atoms)})" if self.coordinating_atoms else ""
        cap_str = f" [{self.cap_type.value}]" if self.cap_type != CapType.NONE else ""
        return f"{self.chain}:{self.resid} {self.resname}{coord_str}{cap_str}"


class ModelBuilder:
    """
    Base class for building QM models from RedoxSite.

    Implements MCPB capping logic to prevent dangling bonds.
    """

    def __init__(self, redox_site: RedoxSite, pdb_file: str, console: Optional[Console] = None):
        """
        Initialize model builder.

        Args:
            redox_site: RedoxSite object with metal center and coordinating atoms
            pdb_file: Path to PDB file
            console: Rich console for output
        """
        self.redox_site = redox_site
        self.pdb_file = Path(pdb_file)
        self.console = console or Console()

        # Parse PDB structure
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('structure', str(self.pdb_file))

        # Build residue map
        self.residue_map = {}  # (chain, resid) -> BioPython Residue
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    het_flag, resseq, icode = residue.get_id()
                    key = (chain.id, resseq)
                    self.residue_map[key] = residue

        # Model residues
        self.model_residues: List[ModelResidue] = []

    def get_coordinating_residues(self) -> List[Tuple[str, int, str, List[str]]]:
        """
        Get residues that have atoms with coordinate bonds to metal.

        Returns:
            List of (chain, resid, resname, [coordinating_atom_names])
        """
        # Build set of metal center coordinates
        center_coords = {center.coords for center in self.redox_site.centers}

        # Find atoms with coordinate bonds to metals
        coordinating_residues = {}  # (chain, resid) -> (resname, [atom_names])

        for bond in self.redox_site.bonds:
            if bond.chemical_type == 'coordinate':
                # One end is metal, other is ligand
                for coords, info in [(bond.atom1_coords, bond.atom1_residue_info),
                                     (bond.atom2_coords, bond.atom2_residue_info)]:
                    if coords not in center_coords:
                        # This is a ligand atom
                        chain = info.get('chain')
                        resid = info.get('resid')
                        resname = info.get('resname')
                        atom_name = info.get('atom_name')

                        key = (chain, resid)
                        if key not in coordinating_residues:
                            coordinating_residues[key] = (resname, [])
                        coordinating_residues[key][1].append(atom_name)

        # Convert to list
        result = []
        for (chain, resid), (resname, atom_names) in coordinating_residues.items():
            result.append((chain, resid, resname, sorted(set(atom_names))))

        # Sort by chain, then resid
        result.sort(key=lambda x: (x[0], x[1]))

        return result

    def get_all_residues_in_site(self) -> List[Tuple[str, int, str]]:
        """
        Get all unique residues in RedoxSite (coordinating + nearby).

        Returns:
            List of (chain, resid, resname)
        """
        residues = {}  # (chain, resid) -> resname

        for atom in self.redox_site.atoms:
            key = (atom.chain, atom.resid)
            residues[key] = atom.resname

        # Convert to list and sort
        result = [(chain, resid, resname) for (chain, resid), resname in residues.items()]
        result.sort(key=lambda x: (x[0], x[1]))

        return result

    def is_standard_amino_acid(self, resname: str) -> bool:
        """Check if residue is a standard amino acid."""
        STANDARD_AA = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'HID', 'HIE', 'HIP', 'CYX', 'CYM'
        }
        return resname in STANDARD_AA

    def is_metal_ion(self, resname: str) -> bool:
        """Check if residue is a metal ion using comprehensive metal set."""
        return resname.upper() in METALS

    def is_proteogenic(self, chain: str, resid: int) -> bool:
        """
        Check if residue is proteogenic (has backbone atoms N, CA, C, O).

        Only proteogenic residues should get capping groups.
        Organic cofactors/ligands (like HEC) should not be capped.

        Args:
            chain: Chain ID
            resid: Residue number

        Returns:
            True if residue has all backbone atoms (N, CA, C, O)
        """
        residue = self.residue_map.get((chain, resid))
        if residue is None:
            return False

        atom_names = {atom.get_id() for atom in residue.get_atoms()}
        backbone_atoms = {'N', 'CA', 'C', 'O'}

        # Check if all backbone atoms are present
        return backbone_atoms.issubset(atom_names)

    def _get_pdb_resname(self, chain: str, resid: int) -> str:
        """Get residue name from PDB structure."""
        residue = self.residue_map.get((chain, resid))
        if residue is not None:
            return residue.get_resname().strip()
        return "UNK"

    def add_residue(self, chain: str, resid: int, resname: str,
                   cap_type: CapType = CapType.NONE,
                   coordinating_atoms: List[str] = None):
        """
        Add a residue to the model.

        Args:
            chain: Chain ID
            resid: Residue number
            resname: Residue name
            cap_type: Capping group type
            coordinating_atoms: List of atom names that coordinate to metal
        """
        # Guard against adding two residues at the same (chain, resid): a model
        # must have at most one residue per position. Keep the first (the
        # coordinating/gap residue added earlier wins) and skip the duplicate,
        # rather than silently corrupting the model with two overlapping
        # residues at the same number.
        for existing in self.model_residues:
            if existing.chain == chain and existing.resid == resid:
                self.console.print(
                    f"[grey50]  (skipping duplicate at {chain}:{resid}: keeping "
                    f"{existing.resname}, not re-adding {resname})[/grey50]"
                )
                return

        is_coordinating = coordinating_atoms is not None and len(coordinating_atoms) > 0

        residue = ModelResidue(
            chain=chain,
            resid=resid,
            resname=resname,
            cap_type=cap_type,
            is_coordinating=is_coordinating,
            coordinating_atoms=coordinating_atoms or []
        )

        self.model_residues.append(residue)

    def get_residue_atoms_from_pdb(self, chain: str, resid: int) -> List[str]:
        """
        Get atom names for a residue from PDB structure.

        Args:
            chain: Chain ID
            resid: Residue number

        Returns:
            List of atom names
        """
        residue = self.residue_map.get((chain, resid))
        if residue is None:
            return []

        return [atom.get_id() for atom in residue.get_atoms()]

    def show_summary(self):
        """Display summary of model."""
        self.console.print(f"\n[bold]Model Summary:[/bold]")
        self.console.print(f"  Total residues: {len(self.model_residues)}")

        coordinating = [r for r in self.model_residues if r.is_coordinating]
        caps = [r for r in self.model_residues if r.cap_type != CapType.NONE]

        self.console.print(f"  Coordinating residues: {len(coordinating)}")
        self.console.print(f"  Capping groups: {len(caps)}")

        if caps:
            cap_counts = {}
            for r in caps:
                cap_counts[r.cap_type] = cap_counts.get(r.cap_type, 0) + 1

            for cap_type, count in sorted(cap_counts.items(), key=lambda x: x[0].value):
                self.console.print(f"    • {cap_type.value}: {count}")

    def _fill_single_residue_gaps(self, selected_residues: List[Tuple[str, int]],
                                  use_gly: bool = False) -> int:
        """
        Fill 1-residue gaps between coordinating residues.

        When two coordinating amino-acid residues are separated by exactly one
        residue, include that intervening residue to prevent cap collision bugs.
        Only real peptide gaps are filled: a "gap" whose flanks are not both
        amino acids (e.g. between a ligand and a metal) is skipped.

        use_gly mirrors the large model: bridge the gap with a neutral GLY
        (backbone only) instead of the actual PDB residue, so an incidental
        charged/bulky sidechain (e.g. an ARG that merely happens to fall between
        two coordinating cysteines) doesn't perturb the QM or the model's net
        charge. The bridge exists only for backbone continuity, so its sidechain
        is not needed.

        Args:
            selected_residues: List of (chain, resid) tuples for coordinating residues
            use_gly: If True, bridge with GLY; if False, use the actual PDB residue

        Returns:
            Number of residues added
        """
        added_count = 0

        # Group by chain
        by_chain = {}
        for chain, resid in selected_residues:
            if chain not in by_chain:
                by_chain[chain] = []
            by_chain[chain].append(resid)

        # Sort each chain
        for chain in by_chain:
            by_chain[chain].sort()

        # Find and fill 1-residue gaps
        for chain, resids in by_chain.items():
            for i in range(len(resids) - 1):
                resi = resids[i]
                resj = resids[i + 1]
                gap_size = resj - resi - 1

                if gap_size == 1:
                    # Only bridge a real peptide gap: both flanks must sit on a
                    # protein backbone (N, CA, C, O present). Skips ligand→metal
                    # spans like X9E:259→ZN:261 while still handling modified
                    # amino acids that a hardcoded name list would miss.
                    if not (self.is_proteogenic(chain, resi) and self.is_proteogenic(chain, resj)):
                        continue

                    gap_resid = resi + 1
                    residue = self.residue_map.get((chain, gap_resid))

                    if residue:
                        resname = residue.get_resname()
                        if use_gly:
                            self.console.print(
                                f"  • Filling 1-residue gap {chain}:{resi}→{resj}: "
                                f"adding {chain}:{gap_resid} GLY (neutral bridge, was {resname})"
                            )
                            self.add_residue(chain, gap_resid, 'GLY', CapType.GLY, [])
                        else:
                            self.console.print(
                                f"  • Filling 1-residue gap {chain}:{resi}→{resj}: "
                                f"adding {chain}:{gap_resid} {resname}"
                            )
                            self.add_residue(chain, gap_resid, resname, CapType.FULL, [])
                        added_count += 1

        return added_count


    def has_single_residue_gaps(self, selected_residues: List[Tuple[str, int]]) -> bool:
        """True if any two selected residues are separated by exactly one
        proteogenic residue — i.e. _fill_single_residue_gaps would bridge it.

        Used to decide whether asking the user about GLY-vs-actual bridging is
        even meaningful for this site.
        """
        by_chain: Dict[str, List[int]] = {}
        for chain, resid in selected_residues:
            by_chain.setdefault(chain, []).append(resid)
        for chain, resids in by_chain.items():
            resids.sort()
            for a, b in zip(resids, resids[1:]):
                if (b - a - 1) == 1 and self.is_proteogenic(chain, a) and self.is_proteogenic(chain, b):
                    return True
        return False


class SmallModelBuilder(ModelBuilder):
    """
    Builder for small (sidechain) models.

    Small models are used for QM calculation of bonded parameters.
    Converts residues to sidechain-only representation with capping groups.
    """

    def __init__(self, redox_site: RedoxSite, pdb_file: str, console: Optional[Console] = None):
        """Initialize small model builder."""
        super().__init__(redox_site, pdb_file, console)

    def build_from_residues(self, selected_residues: List[Tuple[str, int]],
                            use_gly: bool = False):
        """
        Build small model from selected residues.

        Applies MCPB capping logic based on bonding patterns.

        Args:
            selected_residues: List of (chain, resid) tuples
            use_gly: If True, bridge 1-residue gaps with neutral GLY instead of
                the actual PDB residue (see _fill_single_residue_gaps)
        """
        self.console.print(f"\n[bold cyan]Building small model from {len(selected_residues)} residue(s)...[/bold cyan]")

        # Get coordinating info
        coordinating_info = {}  # (chain, resid) -> [atom_names]
        for chain, resid, resname, atoms in self.get_coordinating_residues():
            coordinating_info[(chain, resid)] = atoms

        # Add selected residues
        self.console.print(f"[grey50]Step 1: Adding residues to model...[/grey50]")
        for chain, resid in selected_residues:
            residue = self.residue_map.get((chain, resid))
            if residue is None:
                continue

            resname = residue.get_resname()
            coord_atoms = coordinating_info.get((chain, resid), [])

            if coord_atoms:
                atoms_str = ", ".join(coord_atoms)
                self.console.print(f"  • {chain}:{resid} {resname:3s} - coordinating atoms: {atoms_str}")
            else:
                self.console.print(f"  • {chain}:{resid} {resname:3s}")

            self.add_residue(chain, resid, resname, CapType.NONE, coord_atoms)

        # Fill 1-residue gaps between coordinating residues
        self.console.print(f"\n[grey50]Step 2: Filling 1-residue gaps between coordinating residues...[/grey50]")
        added_count = self._fill_single_residue_gaps(selected_residues, use_gly=use_gly)
        if added_count == 0:
            self.console.print(f"  • No 1-residue gaps found")

        # Apply MCPB capping logic
        self.console.print(f"\n[grey50]Step 3: Analyzing bonding patterns and applying capping groups...[/grey50]")
        self._apply_small_model_capping(selected_residues)

        # Sort residues by chain and resid
        self.model_residues.sort(key=lambda r: (r.chain, r.resid))

    def _apply_small_model_capping(self, selected_residues: List[Tuple[str, int]]):
        """
        Apply MCPB capping logic for small model.

        Full MCPB implementation based on gene_model_files.py:1700-1823

        Logic:
        - Analyzes which atoms coordinate to metal (N, O, sidechain)
        - Determines if residue is N-term, C-term, or internal
        - Applies appropriate capping based on bonding pattern
        """
        # Build set of all residues currently in model (includes gap-filled residues)
        model_set = {(r.chain, r.resid) for r in self.model_residues}

        # Determine which atoms coordinate to metal for each residue
        bonding_atoms = {}  # (chain, resid) -> list of atom names bonding to metal

        # Build set of metal center coordinates
        center_coords = {center.coords for center in self.redox_site.centers}

        # Find atoms with coordinate bonds to metals
        for bond in self.redox_site.bonds:
            if bond.chemical_type == 'coordinate':
                # One end is metal, other is ligand
                for coords, info in [(bond.atom1_coords, bond.atom1_residue_info),
                                     (bond.atom2_coords, bond.atom2_residue_info)]:
                    if coords not in center_coords:
                        # This is a ligand atom
                        chain = info.get('chain')
                        resid = info.get('resid')
                        atom_name = info.get('atom_name')

                        key = (chain, resid)
                        if key not in bonding_atoms:
                            bonding_atoms[key] = []
                        bonding_atoms[key].append(atom_name)

        # Apply capping rules for each residue
        for residue in list(self.model_residues):
            if residue.cap_type != CapType.NONE:
                continue  # Skip if already capped

            chain = residue.chain
            resid = residue.resid
            resname = residue.resname

            # Get bonding atoms for this residue
            res_bonding_atoms = set(bonding_atoms.get((chain, resid), []))

            # Check if neighbors are in model (includes gap-filled residues)
            prev_in_model = (chain, resid - 1) in model_set
            next_in_model = (chain, resid + 1) in model_set

            # Check if neighbors are caps that require this residue to keep backbone atoms
            prev_is_ace_or_kco = False
            next_is_nme_or_knh = False

            for r in self.model_residues:
                if r.chain == chain and r.resid == resid - 1:
                    if r.cap_type in (CapType.ACE, CapType.KCO):
                        prev_is_ace_or_kco = True
                elif r.chain == chain and r.resid == resid + 1:
                    if r.cap_type in (CapType.NME, CapType.KNH):
                        next_is_nme_or_knh = True

            # Skip metal ions - they don't get capping groups
            if self.is_metal_ion(resname):
                self.console.print(f"  • {chain}:{resid} {resname:3s} - Metal ion, no capping needed")
                continue

            # Skip non-proteogenic residues (organic cofactors/ligands)
            if not self.is_proteogenic(chain, resid):
                self.console.print(f"  • {chain}:{resid} {resname:3s} - Non-proteogenic residue (no backbone), no capping needed")
                continue

            # Determine residue type
            is_nterm = self._is_nterm_residue(chain, resid)
            is_cterm = self._is_cterm_residue(chain, resid)
            is_std = self.is_standard_amino_acid(resname)

            # Apply MCPB capping logic based on bonding pattern and position
            if is_std:
                has_backbone_n = 'N' in res_bonding_atoms or 'N3' in res_bonding_atoms
                has_backbone_o = 'O' in res_bonding_atoms or 'OXT' in res_bonding_atoms
                has_sidechain = bool(res_bonding_atoms - {'N', 'N3', 'O', 'OXT'})

                # Case 1: Both N and O coordinate
                if has_backbone_n and has_backbone_o:
                    # Convert to GLY, cap both neighbors
                    self.console.print(f"  • {chain}:{resid} {resname:3s} - Both N and O coordinate → GLY")
                    residue.cap_type = CapType.GLY
                    if not prev_in_model:
                        self.console.print(f"    ↳ Adding ACE cap at {chain}:{resid-1} (N-terminal)")
                        self.add_residue(chain, resid - 1, 'ACE', CapType.ACE)
                    if not next_in_model:
                        self.console.print(f"    ↳ Adding NME cap at {chain}:{resid+1} (C-terminal)")
                        self.add_residue(chain, resid + 1, 'NME', CapType.NME)

                # Case 2: Only O coordinates
                elif has_backbone_o:
                    if prev_in_model:
                        # The N-terminal neighbour is in the model and its
                        # backbone carbonyl bonds to THIS residue's N. Reducing
                        # to ACE/KCO drops that N, severing the peptide bond and
                        # leaving the neighbour's carbonyl as a dangling acyl
                        # radical (odd-electron fragment) while this residue's CA
                        # gets over-methylated. Keep the backbone instead so the
                        # C(prev)=O…N(this) peptide bond survives.
                        residue.cap_type = CapType.FULL if has_sidechain else CapType.GLY
                        kept = 'FULL' if has_sidechain else 'GLY'
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - O coordinates, bonded to prev residue → keep backbone ({kept})")
                    elif has_sidechain:
                        # Keep C=O, modify rest to sidechain
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - O + sidechain coordinate → KCO (keep C=O)")
                        residue.cap_type = CapType.KCO
                    else:
                        # Only O bonds, convert to ACE (N-terminal side terminates)
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - Only O coordinates → ACE")
                        residue.cap_type = CapType.ACE

                    # Cap next residue
                    if not next_in_model:
                        self.console.print(f"    ↳ Adding NME cap at {chain}:{resid+1} (C-terminal)")
                        self.add_residue(chain, resid + 1, 'NME', CapType.NME)

                # Case 3: Only N coordinates
                elif has_backbone_n:
                    if next_in_model:
                        # Mirror of Case 2: the C-terminal neighbour is in the
                        # model and its backbone N bonds to THIS residue's
                        # carbonyl C. Reducing to NME/KNH drops that C=O,
                        # severing C(this)=O…N(next). Keep the backbone instead.
                        residue.cap_type = CapType.FULL if has_sidechain else CapType.GLY
                        kept = 'FULL' if has_sidechain else 'GLY'
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - N coordinates, bonded to next residue → keep backbone ({kept})")
                    elif has_sidechain:
                        # Keep N-H, modify rest to sidechain
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - N + sidechain coordinate → KNH (keep N-H)")
                        residue.cap_type = CapType.KNH
                    else:
                        # Only N bonds, convert to NME (C-terminal side terminates)
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - Only N coordinates → NME")
                        residue.cap_type = CapType.NME

                    # Cap previous residue
                    if not prev_in_model:
                        self.console.print(f"    ↳ Adding ACE cap at {chain}:{resid-1} (N-terminal)")
                        self.add_residue(chain, resid - 1, 'ACE', CapType.ACE)

                # Case 4: Only sidechain coordinates (or no metal coordination - gap residue)
                else:
                    # Check if we need to keep backbone atoms for neighboring caps
                    if prev_is_ace_or_kco and next_is_nme_or_knh:
                        # Both neighbors are caps - keep full backbone
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - Sidechain only, between caps → GLY")
                        residue.cap_type = CapType.GLY
                    elif prev_is_ace_or_kco:
                        # Previous residue has C=O that bonds to our N
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - Sidechain only, after ACE/KCO → KNH (keep N-H)")
                        residue.cap_type = CapType.KNH
                    elif next_is_nme_or_knh:
                        # Next residue has N-H that bonds to our C=O
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - Sidechain only, before NME/KNH → KCO (keep C=O)")
                        residue.cap_type = CapType.KCO
                    elif prev_in_model and next_in_model:
                        # Middle of connected chain - keep full residue (backbone + sidechain)
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - Sidechain only, middle of chain → FULL (backbone + sidechain)")
                        residue.cap_type = CapType.FULL
                    elif not prev_in_model and next_in_model:
                        # N-terminus of connected chain - add ACE cap, keep full residue for connectivity
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - Sidechain only, N-terminus of chain → ACE + FULL")
                        self.add_residue(chain, resid - 1, 'ACE', CapType.ACE)
                        residue.cap_type = CapType.FULL
                    elif prev_in_model and not next_in_model:
                        # C-terminus of connected chain - keep full residue for connectivity, add NME cap
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - Sidechain only, C-terminus of chain → FULL + NME")
                        residue.cap_type = CapType.FULL
                        self.add_residue(chain, resid + 1, 'NME', CapType.NME)
                    else:
                        # Truly isolated (no neighbors in model)
                        # Backbone is removed (CA→CH3), leaving isolated CH3-sidechain fragment
                        self.console.print(f"  • {chain}:{resid} {resname:3s} - Sidechain only, isolated → no caps (CH3-sidechain)")

            else:
                # Non-standard residue - keep as-is, cap neighbors
                self.console.print(f"  • {chain}:{resid} {resname:3s} - Non-standard residue → keep as-is, cap neighbors")
                if not prev_in_model:
                    self.console.print(f"    ↳ Adding ACE cap at {chain}:{resid-1} (N-terminal)")
                    self.add_residue(chain, resid - 1, 'ACE', CapType.ACE)
                if not next_in_model:
                    self.console.print(f"    ↳ Adding NME cap at {chain}:{resid+1} (C-terminal)")
                    self.add_residue(chain, resid + 1, 'NME', CapType.NME)

    def _is_nterm_residue(self, chain: str, resid: int) -> bool:
        """Check if residue is N-terminal in its chain."""
        # Get all residues in this chain
        chain_residues = [(c, r) for (c, r) in self.residue_map.keys() if c == chain]
        if not chain_residues:
            return False

        min_resid = min(r for c, r in chain_residues)
        return resid == min_resid

    def _is_cterm_residue(self, chain: str, resid: int) -> bool:
        """Check if residue is C-terminal in its chain."""
        # Get all residues in this chain
        chain_residues = [(c, r) for (c, r) in self.residue_map.keys() if c == chain]
        if not chain_residues:
            return False

        max_resid = max(r for c, r in chain_residues)
        return resid == max_resid


class LargeModelBuilder(ModelBuilder):
    """
    Builder for large (full residue) models.

    Large models are used for RESP charge fitting.
    Keeps full residues with optional gap filling.
    """

    def __init__(self, redox_site: RedoxSite, pdb_file: str, console: Optional[Console] = None):
        """Initialize large model builder."""
        super().__init__(redox_site, pdb_file, console)

    def build_from_residues(self, selected_residues: List[Tuple[str, int]],
                           max_gap: int = 5, use_gly: bool = True):
        """
        Build large model from selected residues.

        Args:
            selected_residues: List of (chain, resid) tuples
            max_gap: Maximum gap size to fill (residues)
            use_gly: If True, use GLY for gap filling; if False, use actual PDB residues
        """
        self.console.print(f"\n[bold cyan]Building large model from {len(selected_residues)} residue(s)...[/bold cyan]")

        # Get coordinating info
        coordinating_info = {}  # (chain, resid) -> [atom_names]
        for chain, resid, resname, atoms in self.get_coordinating_residues():
            coordinating_info[(chain, resid)] = atoms

        # Add selected residues
        self.console.print(f"[grey50]Step 1: Adding core residues to model...[/grey50]")
        for chain, resid in selected_residues:
            residue = self.residue_map.get((chain, resid))
            if residue is None:
                continue

            resname = residue.get_resname()
            coord_atoms = coordinating_info.get((chain, resid), [])

            if coord_atoms:
                atoms_str = ", ".join(coord_atoms)
                self.console.print(f"  • {chain}:{resid} {resname:3s} - coordinating atoms: {atoms_str}")
            else:
                self.console.print(f"  • {chain}:{resid} {resname:3s}")

            self.add_residue(chain, resid, resname, CapType.FULL, coord_atoms)

        # Apply gap filling
        if max_gap > 0:
            # _fill_gaps already covers ALL gaps 1..max_gap, including single-
            # residue gaps, so it alone prevents cap collisions. Calling
            # _fill_single_residue_gaps as well would try to re-add the same
            # positions (with the actual residue instead of GLY), producing the
            # "Duplicate residue detected" churn.
            self.console.print(f"\n[grey50]Step 2: Filling gaps ≤{max_gap} residues...[/grey50]")
            gap_count = self._fill_gaps(selected_residues, max_gap, use_gly)
            if gap_count == 0:
                self.console.print(f"  • No gaps found (or all gaps >{max_gap} residues)")
        else:
            # No general gap filling requested, but still bridge single-residue
            # peptide gaps so adjacent caps don't collide.
            gap1_count = self._fill_single_residue_gaps(selected_residues)
            if gap1_count > 0:
                self.console.print(f"  Filled {gap1_count} 1-residue gap(s) between coordinating residues")

        # Apply capping
        self.console.print(f"\n[grey50]Step 3: Applying terminal capping groups...[/grey50]")
        self._apply_large_model_capping(selected_residues)

        # Sort residues by chain and resid
        self.model_residues.sort(key=lambda r: (r.chain, r.resid))

    def _fill_gaps(self, selected_residues: List[Tuple[str, int]],
                   max_gap: int, use_gly: bool) -> int:
        """
        Fill gaps between selected residues.

        Based on MCPB gene_model_files.py:1859-1868

        Returns:
            Number of residues added
        """
        added_count = 0

        # Group by chain
        by_chain = {}
        for chain, resid in selected_residues:
            if chain not in by_chain:
                by_chain[chain] = []
            by_chain[chain].append(resid)

        # Sort each chain
        for chain in by_chain:
            by_chain[chain].sort()

        # Find and fill gaps
        for chain, resids in by_chain.items():
            for i in range(len(resids) - 1):
                resi = resids[i]
                resj = resids[i + 1]
                gap_size = resj - resi - 1

                if 0 < gap_size <= max_gap:
                    # Only fill a real peptide gap: both flanks must sit on a
                    # protein backbone (N, CA, C, O present). Skips ligand→metal/
                    # water spans (e.g. X9E:259→ZN:261) that merely have adjacent
                    # residue numbers, while still handling modified amino acids.
                    if not (self.is_proteogenic(chain, resi) and self.is_proteogenic(chain, resj)):
                        continue

                    # Fill the gap
                    gap_range = list(range(resi + 1, resj))
                    if use_gly:
                        self.console.print(f"  • Filling gap {chain}:{resi}→{resj} ({gap_size} residues) with GLY")
                    else:
                        self.console.print(f"  • Filling gap {chain}:{resi}→{resj} ({gap_size} residues) with PDB residues")

                    for gap_resid in gap_range:
                        if use_gly:
                            # Use GLY
                            self.console.print(f"    ↳ {chain}:{gap_resid} GLY")
                            self.add_residue(chain, gap_resid, 'GLY', CapType.GLY)
                            added_count += 1
                        else:
                            # Use actual PDB residue
                            residue = self.residue_map.get((chain, gap_resid))
                            if residue:
                                resname = residue.get_resname()
                                self.console.print(f"    ↳ {chain}:{gap_resid} {resname}")
                                self.add_residue(chain, gap_resid, resname, CapType.FULL)
                                added_count += 1

        return added_count

    def _apply_large_model_capping(self, selected_residues: List[Tuple[str, int]]):
        """
        Apply MCPB capping logic for large model.

        Based on MCPB gene_model_files.py:1870-1883
        """
        # Get all residues currently in model (including gap-filled)
        model_set = {(r.chain, r.resid) for r in self.model_residues}
        # Track cap positions to detect conflicts
        cap_positions = {}  # (chain, resid) -> cap_type

        # Track caps added
        caps_added = 0

        # Add caps for termini
        for residue in list(self.model_residues):
            chain = residue.chain
            resid = residue.resid
            resname = residue.resname

            # Skip metal ions - they don't get capping groups
            if self.is_metal_ion(resname):
                continue

            # Skip non-proteogenic residues (organic cofactors/ligands)
            if not self.is_proteogenic(chain, resid):
                continue

            # Check if neighbors are in model
            prev_in_model = (chain, resid - 1) in model_set
            next_in_model = (chain, resid + 1) in model_set

            # Add caps for missing neighbors
            if not prev_in_model:
                cap_key = (chain, resid - 1)
                if cap_key in cap_positions and cap_positions[cap_key] != 'ACE':
                    # Conflict: another residue already placed a different cap here
                    # Include the actual residue instead of capping
                    pdb_resname = self._get_pdb_resname(chain, resid - 1)
                    self.console.print(
                        f"  • {chain}:{resid} {resname:3s} - Including {chain}:{resid-1} {pdb_resname} "
                        f"as full residue (bridges to adjacent coordinating residue)")
                    # Replace the existing cap with a full residue
                    self.model_residues = [
                        r for r in self.model_residues
                        if not (r.chain == chain and r.resid == resid - 1)
                    ]
                    self.add_residue(chain, resid - 1, pdb_resname, CapType.FULL)
                    model_set.add(cap_key)
                    cap_positions[cap_key] = 'FULL'
                else:
                    self.console.print(f"  • {chain}:{resid} {resname:3s} - Adding ACE cap at {chain}:{resid-1} (N-terminal)")
                    self.add_residue(chain, resid - 1, 'ACE', CapType.ACE)
                    cap_positions[cap_key] = 'ACE'
                    model_set.add(cap_key)
                    caps_added += 1

            if not next_in_model:
                cap_key = (chain, resid + 1)
                if cap_key in cap_positions and cap_positions[cap_key] != 'NME':
                    # Conflict: another residue already placed a different cap here
                    pdb_resname = self._get_pdb_resname(chain, resid + 1)
                    self.console.print(
                        f"  • {chain}:{resid} {resname:3s} - Including {chain}:{resid+1} {pdb_resname} "
                        f"as full residue (bridges to adjacent coordinating residue)")
                    self.model_residues = [
                        r for r in self.model_residues
                        if not (r.chain == chain and r.resid == resid + 1)
                    ]
                    self.add_residue(chain, resid + 1, pdb_resname, CapType.FULL)
                    model_set.add(cap_key)
                    cap_positions[cap_key] = 'FULL'
                else:
                    self.console.print(f"  • {chain}:{resid} {resname:3s} - Adding NME cap at {chain}:{resid+1} (C-terminal)")
                    self.add_residue(chain, resid + 1, 'NME', CapType.NME)
                    cap_positions[cap_key] = 'NME'
                    model_set.add(cap_key)
                    caps_added += 1

        if caps_added == 0:
            self.console.print(f"  • No terminal caps needed (all residues have neighbors)")
