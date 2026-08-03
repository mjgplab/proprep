"""Charge rebalancing by removal of bulk ions.

After `apply_state.build_production_prmtop` mutates per-residue charges to
match a converged state map, the system net charge typically drifts away
from neutrality by some integer Δq. This module restores the target net
charge by removing single-atom (bulk) ions — never adds.

Why removal-only
----------------
- parmed-only edits: removing a residue is `structure.strip(":N")`. Box,
  coordinates, and everything else stay intact; solvent equilibrates
  around the gap during MD. Adding a new ion would require fabricating
  prmtop atom/charge/mass/type/exclusion entries by hand and finding
  collision-free coordinates — fragile.
- tleap re-solvation is excluded: it reverts custom charges on `loadpdb`,
  destroying the very state map this whole pipeline computed.

Sign convention
---------------
- Δq_protein = +k (state map made the protein more positive)
  → remove k cations  (each cation removal subtracts +1 (or +valence) from
                        the system net charge)
- Δq_protein = −k → remove k anions

Multivalent granularity
-----------------------
Mg²⁺ / bulk Ca²⁺ adjust the net charge in steps of 2. Default ranking
prefers +1 species (Na⁺, K⁺) for fine control; multivalents fall back
only if the +1 pool is exhausted. If the residual cannot be closed
exactly with available species, the function raises rather than leaving a
silent ±1 imbalance.

Distance ranking
----------------
Within a species, the ions farthest from the protein heavy atoms are
removed first — these are the most truly-bulk ions, leaving the diffuse
double layer near the surface mostly intact.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

import parmed as pmd

from .ions import is_single_atom_residue, is_water


# Default species priorities. +1 ions first (fine granularity), then +2.
# Anions are usually only Cl- in practice.
DEFAULT_CATION_PREFERENCE: List[str] = ["Na+", "K+", "Mg2+", "Ca2+"]
DEFAULT_ANION_PREFERENCE:  List[str] = ["Cl-", "Br-", "F-"]


# A residue identifier in the form produced by strip_bulk_ions audit:
# (resnum_1based, resname). Stored as Tuple[int, str] there; we accept
# either order in the keep_ions set.
SiteId = Tuple[str, int]


def _normalize_keep(keep_ions) -> Set[SiteId]:
    """Coerce keep_ions into a normalized {(resname, resnum_1based), ...}.

    Accepts either {(resname, resnum), ...} or {(resnum, resname), ...}
    (the latter is what `strip_bulk_ions` returns in its audit dict).
    """
    if not keep_ions:
        return set()
    out: Set[SiteId] = set()
    for entry in keep_ions:
        a, b = entry
        if isinstance(a, str):
            out.add((a, int(b)))
        else:
            out.add((str(b), int(a)))
    return out


def _ion_charge(residue: pmd.Residue) -> float:
    """Net charge of a single-atom residue (the atom's charge)."""
    return float(residue.atoms[0].charge)


def _protein_heavy_coords(structure: pmd.Structure):
    """Heavy-atom coordinates of every non-ion, non-water residue."""
    coords = []
    for r in structure.residues:
        if is_single_atom_residue(r) or is_water(r):
            continue
        for a in r.atoms:
            if a.atomic_number > 1:
                coords.append((a.xx, a.xy, a.xz))
    return coords


def _rank_ions_by_distance(structure: pmd.Structure,
                            candidate_residues: List[pmd.Residue]
                            ) -> List[Tuple[pmd.Residue, float]]:
    """Return [(residue, dist_to_nearest_protein_atom), ...] sorted
    farthest-first. Uses scipy KDTree if available; else O(N·M)."""
    if not candidate_residues:
        return []
    prot = _protein_heavy_coords(structure)
    if not prot:
        # No protein → just leave order as-is, all distances equal
        return [(r, 0.0) for r in candidate_residues]
    try:
        import numpy as np
        from scipy.spatial import cKDTree
        tree = cKDTree(np.array(prot))
        out = []
        for r in candidate_residues:
            a = r.atoms[0]
            d, _ = tree.query([a.xx, a.xy, a.xz], k=1)
            out.append((r, float(d)))
    except ImportError:
        out = []
        for r in candidate_residues:
            ax, ay, az = r.atoms[0].xx, r.atoms[0].xy, r.atoms[0].xz
            best = float("inf")
            for cx, cy, cz in prot:
                d2 = (ax-cx)**2 + (ay-cy)**2 + (az-cz)**2
                if d2 < best:
                    best = d2
            out.append((r, best ** 0.5))
    out.sort(key=lambda x: -x[1])  # farthest first
    return out


def _bulk_pool(structure: pmd.Structure,
                keep: Set[SiteId],
                sign: int) -> Dict[str, List[pmd.Residue]]:
    """Group bulk single-atom residues by species name.

    `sign` selects: +1 keeps cations (charge > 0.5), −1 keeps anions
    (charge < −0.5). Skips anything in `keep`.
    """
    by_species: Dict[str, List[pmd.Residue]] = {}
    for r in structure.residues:
        if not is_single_atom_residue(r):
            continue
        key = (r.name.strip(), r.number + 1)
        if key in keep:
            continue
        q = _ion_charge(r)
        if (sign > 0 and q > 0.5) or (sign < 0 and q < -0.5):
            by_species.setdefault(r.name.strip(), []).append(r)
    return by_species


def rebalance_ions(structure: pmd.Structure,
                    target_net_q: float = 0.0,
                    *,
                    keep_ions: Optional[Set[SiteId]] = None,
                    cation_preference: Optional[List[str]] = None,
                    anion_preference:  Optional[List[str]] = None,
                    tolerance: float = 0.01,
                    verbose: bool = False,
                    ) -> Dict:
    """Remove bulk ions in-place to bring `structure`'s net charge to
    `target_net_q`. Returns an audit dict.

    Parameters
    ----------
    structure : pmd.Structure
        Mutated in place.
    target_net_q : float
        Desired net system charge after rebalance. Default 0.
    keep_ions : Set[Tuple[str, int]] or Set[Tuple[int, str]], optional
        Single-atom residues to NEVER touch (typically structural ions
        from CryoEM, e.g., the 32 Ca²⁺ at residues 3352–3383 in the bundle).
        Accepts either ordering of the (resname, resnum_1based) tuple to
        play nicely with `strip_bulk_ions`'s audit format.
    cation_preference, anion_preference : List[str], optional
        Species ranking to use when removing. Earlier names tried first.
        Default cation order: ["Na+", "K+", "Mg2+", "Ca2+"] — +1 species
        before +2 for granularity.
    tolerance : float
        Net-charge tolerance in elementary-charge units. After integer
        ion removal, |net_q − target| ≤ tolerance is treated as success.
    verbose : bool
        Print per-species removal counts and before/after net q.

    Returns
    -------
    dict with keys:
        net_q_before, net_q_after, target_net_q, residual,
        removed_per_species (Dict[str, int]),
        removed_residue_ids (List[(resname, resnum_1based)]),
        n_atoms_before, n_atoms_after.

    Raises
    ------
    ValueError if the requested net charge cannot be achieved exactly
    (insufficient ions of needed species, or multivalent granularity
    leaves a non-tolerable residual).
    """
    keep = _normalize_keep(keep_ions)
    cation_pref = list(cation_preference or DEFAULT_CATION_PREFERENCE)
    anion_pref  = list(anion_preference  or DEFAULT_ANION_PREFERENCE)

    n_atoms_before = len(structure.atoms)
    q_before = sum(a.charge for a in structure.atoms)
    delta = target_net_q - q_before

    removed_per_species: Dict[str, int] = {}
    removed_ids: List[SiteId] = []

    # Decide direction: positive delta needs anion removal (each Cl- gone
    # adds +1 to net q); negative delta needs cation removal.
    if abs(delta) <= tolerance:
        info = {
            "net_q_before": q_before, "net_q_after": q_before,
            "target_net_q": target_net_q, "residual": float(delta),
            "removed_per_species": {}, "removed_residue_ids": [],
            "n_atoms_before": n_atoms_before, "n_atoms_after": n_atoms_before,
        }
        if verbose:
            print(f"  rebalance_ions: net q={q_before:+.3f} already within "
                  f"{tolerance} of target {target_net_q:+.3f}; no-op.")
        return info

    if delta > 0:
        # Need to ADD positive charge → REMOVE anions
        pool = _bulk_pool(structure, keep, sign=-1)
        pref = anion_pref
        sign_label = "anion"
    else:
        # Need to ADD negative charge → REMOVE cations
        pool = _bulk_pool(structure, keep, sign=+1)
        pref = cation_pref
        sign_label = "cation"

    # Sort each species's pool by distance (farthest first)
    ranked: Dict[str, List[Tuple[pmd.Residue, float]]] = {}
    for sp in pool:
        ranked[sp] = _rank_ions_by_distance(structure, pool[sp])

    # Greedy removal in preference order. Per-removal Δq = -charge_of_ion.
    # We continue until |residual| ≤ tolerance OR we exhaust the pool.
    remaining = delta
    for species in pref:
        if species not in ranked:
            continue
        for r, _dist in ranked[species]:
            if abs(remaining) <= tolerance:
                break
            ion_q = _ion_charge(r)
            # Removing this ion shifts net q by (-ion_q). Only do so if
            # this brings |remaining| closer to zero (no overshoot beyond
            # tolerance — for multivalents, this prevents flipping sign).
            new_remaining = remaining - (-ion_q)
            if abs(new_remaining) > abs(remaining) + tolerance:
                # Would overshoot — skip this ion (try next species)
                continue
            # Mark for removal
            removed_ids.append((r.name.strip(), r.number + 1))
            removed_per_species[r.name.strip()] = (
                removed_per_species.get(r.name.strip(), 0) + 1)
            remaining = new_remaining
        if abs(remaining) <= tolerance:
            break

    # Apply the removal in one strip call (parmed handles the residue
    # renumbering internally).
    if removed_ids:
        nums = [n for (_n, n) in removed_ids]
        names = sorted({nm for (nm, _) in removed_ids})
        # Build mask by numbers (always exact); name-based mask would
        # over-strip if user keeps some of a species but not others.
        mask = ":" + ",".join(str(n) for n in nums)
        structure.strip(mask)

    q_after = sum(a.charge for a in structure.atoms)
    residual = target_net_q - q_after

    info = {
        "net_q_before":         q_before,
        "net_q_after":          q_after,
        "target_net_q":         target_net_q,
        "residual":             float(residual),
        "removed_per_species":  removed_per_species,
        "removed_residue_ids":  removed_ids,
        "n_atoms_before":       n_atoms_before,
        "n_atoms_after":        len(structure.atoms),
        "sign_label":           sign_label,
    }

    if abs(residual) > tolerance:
        raise ValueError(
            f"rebalance_ions could not reach target net q={target_net_q:+.3f}: "
            f"residual={residual:+.3f} after removing "
            f"{sum(removed_per_species.values())} ions. Pool exhausted or "
            f"multivalent granularity blocks an exact match. "
            f"Removed: {removed_per_species}. "
            f"Consider supplying a +1 ion source or relaxing tolerance.")

    if verbose:
        print(f"  rebalance_ions: net q {q_before:+.3f} → {q_after:+.3f} "
              f"(target {target_net_q:+.3f}; removed "
              f"{sum(removed_per_species.values())} {sign_label}s: "
              f"{removed_per_species})")
    return info
