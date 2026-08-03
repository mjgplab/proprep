"""
ORCA QM/MM Data Structures

Shared data structures for ORCA QM/MM setup.
Maps RedoxSite atoms to prmtop indices by residue ID for correct atom identification.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any

from proprep.structure_prep.comprehensive_redox_detector import RedoxSite


@dataclass
class ORCAQMMMSetup:
    """Complete ORCA QM/MM calculation configuration."""

    # Source data
    redox_site: RedoxSite

    # Atom index sets (0-based, matching prmtop ordering)
    qm_atom_indices: List[int] = field(default_factory=list)
    active_atom_indices: List[int] = field(default_factory=list)

    # File paths
    prmtop_path: str = ""
    rst7_path: str = ""
    pdb_path: str = ""                    # PDB file for ORCA coordinate input
    orca_ff_file: str = ""                # .prms from orca_mm

    # QM settings
    qm_method: str = "B3LYP"
    qm_basis_set: str = "def2-SVP"
    qm_charge: int = 0
    qm_multiplicity: int = 1
    dispersion: str = "D3BJ"

    # QM/MM coupling settings
    embedding: str = "Electrostatic"      # "Electrostatic" or "Mechanical"
    charge_scheme: str = "CS"             # CS, RCD, Z0, Z1, Z2, Z3

    # Job settings
    job_type: str = "Opt"                 # Opt, EnGrad, NumFreq, SP
    n_processors: int = 4
    memory_mb: int = 4000                 # per core
    additional_keywords: str = ""

    # Metadata
    active_region_cutoff: float = 6.0     # Angstroms
    validation_passed: bool = False
    validation_warnings: List[str] = field(default_factory=list)

    # Boundary info (populated during QM region selection)
    boundary_bonds: List[Tuple[int, int]] = field(default_factory=list)  # (qm_idx, mm_idx) pairs

    # Link atom distance overrides (None = use ORCA defaults: C=1.09, N=0.99, O=0.98)
    dist_c_hla: Optional[float] = None
    dist_n_hla: Optional[float] = None
    dist_o_hla: Optional[float] = None

    # Boundary double-counting control
    delete_la_double_counting: bool = True   # ORCA default
    delete_la_bond_double_counting: bool = True  # ORCA default

    # MM water and performance
    rigid_mm_water: bool = False   # ORCA default
    do_nb_fixed_fixed: bool = True   # ORCA default (set False for large systems)

    # Output control
    print_level: int = 1   # ORCA default (0-4)

    # Charge scheme tuning (None = use ORCA defaults)
    scale_cs: Optional[float] = None   # default 0.06
    scale_rcd: Optional[float] = None  # default 0.5

    # Verification mode
    qmmm_setup_only: bool = False   # If True, add QMMMSetup to ! line

    # Use PDB info instead of explicit atom lists
    use_qm_info_from_pdb: bool = False
    use_active_info_from_pdb: bool = False

    # %mm block settings (None = use ORCA defaults)
    coulomb_cutoff: Optional[float] = None    # default 12.0 Angstrom
    lj_cutoff_inner: Optional[float] = None   # default 10.0 Angstrom
    lj_cutoff_outer: Optional[float] = None   # default 12.0 Angstrom
    dielec_const: Optional[float] = None      # default 1.0


def compress_indices(indices: List[int]) -> str:
    """
    Compress a sorted list of 0-based atom indices into ORCA range notation.

    Examples:
        [0, 1, 2, 3, 5, 7, 8, 9] -> "0:3 5 7:9"
        [10] -> "10"
        [0, 2, 4] -> "0 2 4"

    Args:
        indices: Sorted list of 0-based atom indices

    Returns:
        ORCA-formatted range string for use in QMAtoms/ActiveAtoms blocks
    """
    if not indices:
        return ""

    sorted_idx = sorted(set(indices))
    ranges = []
    start = sorted_idx[0]
    end = start

    for i in range(1, len(sorted_idx)):
        if sorted_idx[i] == end + 1:
            end = sorted_idx[i]
        else:
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}:{end}")
            start = sorted_idx[i]
            end = start

    # Final range
    if start == end:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}:{end}")

    return " ".join(ranges)


def _determine_resid_offset(redox_site: RedoxSite, parm) -> int:
    """
    Determine the residue ID offset between RedoxSite and prmtop numbering.

    TLEaP may use 0-based residue numbering while the RedoxSite transformer
    produces 1-based numbers. This checks a sample of atoms to find the
    consistent offset.

    Returns:
        The offset to subtract from RedoxSite resid to get prmtop resid (0 or 1).
    """
    prmtop_residues = set()
    for atom in parm.atoms:
        prmtop_residues.add((atom.residue.name, atom.residue.number))

    hits_at_0 = 0
    hits_at_1 = 0
    for rsa in redox_site.atoms:
        if (rsa.resname, rsa.resid) in prmtop_residues:
            hits_at_0 += 1
        if (rsa.resname, rsa.resid - 1) in prmtop_residues:
            hits_at_1 += 1

    return 1 if hits_at_1 > hits_at_0 else 0


def map_redox_atoms_to_prmtop(redox_site: RedoxSite, parm) -> Tuple[List[int], List[str]]:
    """
    Map RedoxSiteAtom objects to 0-based prmtop atom indices.

    Uses (resname, resid, atom_name) as the matching key. Chain is ignored
    because TLEaP strips chain IDs. The resid offset between RedoxSite and
    prmtop numbering is determined automatically.

    After mapping, completes each matched residue by including ALL atoms
    from the prmtop residue (including hydrogens added by TLEaP that are
    not present in the original PDB/RedoxSite).

    Args:
        redox_site: RedoxSite with atoms to map
        parm: parmed.Structure loaded from prmtop

    Returns:
        Tuple of (sorted 0-based prmtop atom indices, unmapped atom descriptions)
    """
    offset = _determine_resid_offset(redox_site, parm)

    # Build lookup: (resname, resid, atom_name) -> 0-based atom index
    # Chain is omitted because TLEaP strips chain IDs
    prmtop_lookup = {}
    for i, atom in enumerate(parm.atoms):
        key = (atom.residue.name, atom.residue.number, atom.name)
        prmtop_lookup[key] = i

    # Build residue lookup: (resname, resid) -> list of all atom indices
    residue_atoms = {}
    for i, atom in enumerate(parm.atoms):
        res_key = (atom.residue.name, atom.residue.number)
        residue_atoms.setdefault(res_key, []).append(i)

    matched_residues = set()
    unmapped = []

    for rsa in redox_site.atoms:
        prmtop_resid = rsa.resid - offset
        key = (rsa.resname, prmtop_resid, rsa.atom_name)
        idx = prmtop_lookup.get(key)

        if idx is not None:
            matched_residues.add((rsa.resname, prmtop_resid))
        else:
            unmapped.append(f"{rsa.atom_name}:{rsa.resname}{rsa.resid}:{rsa.chain}")

    # Complete residues: include all prmtop atoms (especially TLEaP hydrogens)
    # for every residue that had at least one mapped atom
    qm_indices = set()
    for res_key in matched_residues:
        qm_indices.update(residue_atoms.get(res_key, []))

    return sorted(qm_indices), unmapped
