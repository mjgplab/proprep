"""Cluster cutout for one or two target sites.

A cluster contains:
  1. Bond walk: all residues reachable within `bond_depth` covalent bonds
     of the seed atoms (handles PRN→heme, CYS-CYS disulfides, any bonded
     cofactor without needing residue-name lists).
  2. R-Å radius: any residue with at least one heavy atom within `R` Å of
     a seed heavy atom.

Assumes the structure has been pre-cleaned via `pb_titrate.ions.strip_bulk_ions`,
so any single-atom residue still present is structural. The cluster cutout
does not classify ions itself — it just includes them like any other residue
when they fall within bond walk or R Å.

For pair clusters, takes the union of two single-site cutouts. Since
d_cutoff << 2R is the operating regime (cutoff scales with Debye length;
~10–15 Å vs R ≥ 25 Å), the two single-site spheres overlap and the union
is always contiguous.

# Site abstraction
# -----------------
# `cluster_for_site` and `cluster_for_site_pair` are the canonical entry
# points. They take a `Site` (see `pb_titrate.sites`) which carries an
# *envelope* of structurally coupled residues — necessary for hemes and
# other MCPB groups whose partial charges sum to integer only as a unit.
# Two invariants are enforced:
#   1. the target site's full envelope is unconditionally retained,
#      regardless of cluster_radius;
#   2. for any *other* redox-site envelope (carried on
#      `site.other_redox_envelopes`) that contributes any residue to the
#      cluster, the full envelope is included — partial inclusion would
#      leave non-integer cluster charge.
#
# `cluster_for_target` / `cluster_for_pair` remain as single-residue
# convenience wrappers for callers that don't need envelope handling.
"""
from __future__ import annotations

from typing import Iterable, List, Sequence, Set

import parmed as pmd

from .sites import Site


# ---------------------------------------------------------------------------
# Bond walk (parmed bond graph)
# ---------------------------------------------------------------------------

def bond_walk_from_atoms(seed_atoms: Iterable[pmd.Atom],
                          depth: int = 3) -> Set[int]:
    """BFS over the prmtop bond graph from `seed_atoms`.

    Returns the set of 0-based residue indices reachable within `depth`
    covalent bonds.

    Walks `atom.bond_partners`, which parmed populates from the prmtop's
    BONDS_INC_HYDROGEN and BONDS_WITHOUT_HYDROGEN sections. Coordination
    bonds (e.g., Fe-N for hemes) are included iff they were parameterized
    as bonds (standard for MCPB.py outputs).
    """
    visited_atoms: Set[int] = set()
    visited_residues: Set[int] = set()
    frontier: List[pmd.Atom] = list(seed_atoms)
    for _ in range(depth + 1):
        next_frontier: List[pmd.Atom] = []
        for a in frontier:
            if a.idx in visited_atoms:
                continue
            visited_atoms.add(a.idx)
            visited_residues.add(a.residue.idx)
            for b in a.bond_partners:
                if b.idx not in visited_atoms:
                    next_frontier.append(b)
        frontier = next_frontier
        if not frontier:
            break
    return visited_residues


def bond_walk(start_residue: pmd.Residue, depth: int = 3) -> Set[int]:
    """Single-residue convenience wrapper around `bond_walk_from_atoms`."""
    return bond_walk_from_atoms(start_residue.atoms, depth=depth)


# ---------------------------------------------------------------------------
# R-Å radius
# ---------------------------------------------------------------------------

def _heavy_coords_atoms(atoms: Iterable[pmd.Atom]):
    """Yield (x, y, z) for each heavy atom (atomic_number > 1)."""
    for a in atoms:
        if a.atomic_number > 1:
            yield (a.xx, a.xy, a.xz)


def _heavy_coords(residue: pmd.Residue):
    """Yield (x, y, z) for each heavy atom in a residue."""
    return _heavy_coords_atoms(residue.atoms)


def residues_within_atoms(structure: pmd.Structure,
                            seed_atoms: Iterable[pmd.Atom],
                            R: float,
                            *,
                            exclude_residue_idxs: Set[int] = None) -> Set[int]:
    """Returns 0-based residue indices with any heavy atom within R Å of
    any seed heavy atom. Excludes residues in `exclude_residue_idxs`
    from the search (typically the seed's own residues, since those are
    always retained anyway).

    Uses scipy KDTree if available (much faster on large structures).
    """
    seed_coords = list(_heavy_coords_atoms(seed_atoms))
    if not seed_coords:
        return set()

    excl = exclude_residue_idxs or set()

    try:
        import numpy as np
        from scipy.spatial import cKDTree
        coords = []
        atom_to_res = []
        for r in structure.residues:
            if r.idx in excl:
                continue
            for a in r.atoms:
                if a.atomic_number > 1:
                    coords.append((a.xx, a.xy, a.xz))
                    atom_to_res.append(r.idx)
        if not coords:
            return set()
        tree = cKDTree(np.array(coords))
        out: Set[int] = set()
        for tx, ty, tz in seed_coords:
            for k in tree.query_ball_point([tx, ty, tz], R):
                out.add(atom_to_res[k])
        return out
    except ImportError:
        pass

    # Fallback O(N²) scan
    R2 = R * R
    out: Set[int] = set()
    for r in structure.residues:
        if r.idx in excl:
            continue
        for ax, ay, az in _heavy_coords(r):
            in_range = False
            for tx, ty, tz in seed_coords:
                if (ax-tx)**2 + (ay-ty)**2 + (az-tz)**2 <= R2:
                    in_range = True
                    break
            if in_range:
                out.add(r.idx)
                break
    return out


def residues_within(structure: pmd.Structure,
                     target_residue: pmd.Residue,
                     R: float) -> Set[int]:
    """Single-residue convenience wrapper around `residues_within_atoms`.

    Returns the set including the target residue itself.
    """
    rw = residues_within_atoms(structure, target_residue.atoms, R,
                                exclude_residue_idxs={target_residue.idx})
    rw.add(target_residue.idx)
    return rw


# ---------------------------------------------------------------------------
# Envelope-aware cluster builders (canonical, take a Site)
# ---------------------------------------------------------------------------

def _enforce_other_envelope_completeness(
        keep: Set[int],
        other_envelopes: Sequence[Sequence[int]]) -> Set[int]:
    """All-or-nothing rule for non-target redox-site envelopes.

    For every envelope in `other_envelopes`, if any of its residues is
    already in `keep`, include the entire envelope. This prevents
    fractional cluster charge from partially-included MCPB groups.

    Iterates to a fixed point because pulling in one envelope may put
    its residues inside the radius of yet another envelope (in dense
    multi-cofactor systems).
    """
    keep = set(keep)
    while True:
        added = False
        for env in other_envelopes:
            env_set = set(env)
            if env_set & keep and not env_set.issubset(keep):
                keep |= env_set
                added = True
        if not added:
            break
    return keep


def cluster_for_site(structure: pmd.Structure,
                      site: Site,
                      *,
                      R: float = 25.0,
                      bond_depth: int = 3) -> List[int]:
    """Sorted 0-based residue indices for the cluster around one site.

    Steps:
      1. Bond walk seeded from the *envelope's* atoms (not just the
         titrating residue), depth `bond_depth`.
      2. R-Å heavy-atom search seeded from the envelope's atoms.
      3. Unconditionally include all envelope residues (in case R or
         bond_depth would have clipped any).
      4. All-or-nothing rule for `site.other_redox_envelopes`: any other
         redox-site envelope that contributes at least one residue to
         the cluster gets fully included.

    Re-resolves the envelope by parmed residue index against `structure`,
    so a Site built from one load of a prmtop is valid for any later
    reload of the *same* prmtop file (as commonly happens inside the
    per-state PB loop in `intrinsic.py`).
    """
    envelope_idxs = site.envelope_idxs
    envelope_residues = [structure.residues[i] for i in sorted(envelope_idxs)]
    seed_atoms = [a for r in envelope_residues for a in r.atoms]
    bw = bond_walk_from_atoms(seed_atoms, depth=bond_depth)
    rw = residues_within_atoms(structure, seed_atoms, R,
                                exclude_residue_idxs=envelope_idxs)
    keep = bw | rw | envelope_idxs
    keep = _enforce_other_envelope_completeness(
        keep, site.other_redox_envelopes)
    return sorted(keep)


def cluster_for_site_pair(structure: pmd.Structure,
                           site_i: Site,
                           site_j: Site,
                           *,
                           R: float = 25.0,
                           bond_depth: int = 3) -> List[int]:
    """Union of two single-site clusters, with envelope retention applied
    to both targets.

    `other_redox_envelopes` is taken from both sites and de-duplicated
    by id; the all-or-nothing rule sees the union.
    """
    keep_i = set(cluster_for_site(structure, site_i,
                                    R=R, bond_depth=bond_depth))
    keep_j = set(cluster_for_site(structure, site_j,
                                    R=R, bond_depth=bond_depth))
    keep = keep_i | keep_j

    # Re-apply the all-or-nothing rule against the union of both sites'
    # other-envelopes lists, in case the union pulled new residues in
    # that triggered a partial inclusion.
    seen: List[Sequence[int]] = []
    for env in (*site_i.other_redox_envelopes, *site_j.other_redox_envelopes):
        if env not in seen:
            seen.append(env)
    keep = _enforce_other_envelope_completeness(keep, seen)
    return sorted(keep)


# ---------------------------------------------------------------------------
# Single-residue convenience wrappers (preserve original API)
# ---------------------------------------------------------------------------

def cluster_for_target(structure: pmd.Structure,
                        target_residue: pmd.Residue,
                        *,
                        R: float = 25.0,
                        bond_depth: int = 3) -> List[int]:
    """Cluster around a single residue (no envelope concept).

    Preserved for callers that don't need envelope handling. Equivalent
    to `cluster_for_site` with a Site whose envelope is just the target
    and whose other_redox_envelopes list is empty.
    """
    site = Site(titrating_residue=target_residue,
                envelope_residues=[target_residue])
    return cluster_for_site(structure, site, R=R, bond_depth=bond_depth)


def cluster_for_pair(structure: pmd.Structure,
                      target_i: pmd.Residue,
                      target_j: pmd.Residue,
                      *,
                      R: float = 25.0,
                      bond_depth: int = 3) -> List[int]:
    """Single-residue pair cluster (no envelope concept)."""
    si = Site(titrating_residue=target_i, envelope_residues=[target_i])
    sj = Site(titrating_residue=target_j, envelope_residues=[target_j])
    return cluster_for_site_pair(structure, si, sj,
                                   R=R, bond_depth=bond_depth)


def slice_structure(structure: pmd.Structure,
                     residue_idxs: List[int]) -> pmd.Structure:
    """Reduce ``structure`` to the given 0-based residue indices, in place.

    Strips every other residue and returns the same (now-stripped) object.
    Callers reload the full structure for each PB call and immediately
    reassign, so mutating in place is safe.

    Why ``strip`` rather than parmed mask slicing (``structure[mask]``): the
    mask ``__getitem__`` reduces NTYPES/ATOM_TYPE_INDEX but leaves the
    LENNARD_JONES_ACOEF/BCOEF arrays sized for the ORIGINAL type count, so the
    written prmtop is internally inconsistent (NTYPES no longer matches the LJ
    matrix size). pbsa's ``rdparm2`` then reads past the LJ section and
    segfaults — rarely, depending on byte alignment (observed on AS4-45 in the
    NOS hydride system; other clusters happened to survive the same
    corruption). ``strip`` rebuilds a fully consistent topology and preserves
    residue numbers, which ``apply_state_charges`` relies on to locate the
    titrating residue.
    """
    res_nums = sorted({structure.residues[i].number + 1 for i in residue_idxs})
    if not res_nums:
        raise ValueError("Empty residue selection")
    keep = ",".join(str(n) for n in res_nums)
    structure.strip("!(:%s)" % keep)
    # Collapse the legacy SOLTY/NATYP array. It is unused by modern force
    # fields, but a large distinct atom-type-NAME count (common with custom
    # cofactor atom types) overflows pbsa's static SOLTY array, which aborts
    # rdparm with "a parameter array overflowed". PB (inp=0) never reads it.
    try:
        structure.parm_data["POINTERS"][18] = 1   # NATYP
        structure.parm_data["SOLTY"] = [0.0]
    except (AttributeError, KeyError, TypeError):
        pass  # non-AmberParm or no SOLTY section: nothing to collapse
    return structure


# ---------------------------------------------------------------------------
# Debye length helper (used by coupling.py to set distance cutoff)
# ---------------------------------------------------------------------------

def debye_length(istrng_mM: float, T_K: float = 298.15,
                  epsout: float = 80.0) -> float:
    """Debye screening length κ⁻¹ in Å for `istrng_mM` (1:1 electrolyte).

    Returns ∞ (1e6) if istrng_mM is essentially 0.

        κ⁻¹ = sqrt(εout · ε0 · k_B · T / (2 · N_A · e² · I))

    Numerical: κ⁻¹ [Å] ≈ 3.04 / sqrt(I [mol/L])  at T=298, εout=80
    Generalized for arbitrary T, εout via the formula above.
    """
    if istrng_mM <= 1e-6:
        return 1e6
    I_M = istrng_mM / 1000.0
    base = 3.04 * (epsout / 80.0) ** 0.5 * (T_K / 298.15) ** 0.5
    return base / (I_M ** 0.5)


def recommended_pair_cutoff(istrng_mM: float, *,
                              n_decay: float = 3.0,
                              R_min: float = 6.0) -> float:
    """Return the auto-computed pair-coupling distance cutoff in Å."""
    kappa_inv = debye_length(istrng_mM)
    return max(R_min, n_decay * kappa_inv)
