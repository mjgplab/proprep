"""Generic ion classification for cluster cutout.

Detection is name-agnostic:
  - Single-atom residues are candidate ions (regardless of `r.name`).
  - 3-atom residues with names in WATER_NAMES are candidate waters.

Classification (structural vs bulk) is **user-controlled** via a callable
that takes a single residue and returns True if it should be kept as
structural. Pre-built rule constructors compose into common cases.

Example
-------
    from titrate.ions import (strip_bulk_ions, by_residue_numbers,
                                within_distance_of_protein, combine)

    # Default: any single-atom residue within 4 Å of protein heavy atom kept
    info = strip_bulk_ions(structure)

    # Explicit: keep only specific Ca2+ ions by residue number
    info = strip_bulk_ions(structure,
                             structural_rule=by_residue_numbers([1230, 1245]))

    # Combined: keep all Zn AND any ion within 3 Å of protein
    info = strip_bulk_ions(structure,
                             structural_rule=combine(
                                 by_residue_name({"ZN", "Zn2+"}),
                                 within_distance_of_protein(3.0)))
"""
from __future__ import annotations

from typing import Callable, Iterable, Set

import parmed as pmd


WATER_NAMES = {"WAT", "HOH", "TIP3", "TIP4", "TP3", "TP4", "T3P", "T4P"}


# ---------------------------------------------------------------------------
# Detection (name-agnostic)
# ---------------------------------------------------------------------------

def is_single_atom_residue(r: pmd.Residue) -> bool:
    """True iff the residue has exactly one atom (any element, any name)."""
    return len(r.atoms) == 1


def is_water(r: pmd.Residue) -> bool:
    """True iff the residue is a water (3 atoms with a known water name).

    Detection is by atom count + name list rather than by element so it
    works regardless of model (TIP3P, TIP4P, etc.) and naming convention.
    """
    return r.name.strip() in WATER_NAMES and len(r.atoms) == 3


# ---------------------------------------------------------------------------
# Rule constructors — each returns Callable[[pmd.Residue], bool]
# ---------------------------------------------------------------------------

ResidueRule = Callable[[pmd.Residue], bool]


def within_distance_of_protein(threshold: float = 4.0,
                                 structure: pmd.Structure = None) -> ResidueRule:
    """Keep single-atom residue if within `threshold` Å of any non-ion,
    non-water heavy atom.

    Pass `structure` if you want the rule to precompute once at construction
    time (much faster for repeated application). If `structure` is None,
    the rule scans on each call (slow for large structures).
    """
    if structure is not None:
        # Precompute non-ion/non-water heavy atom coords once
        coords = []
        for r in structure.residues:
            if is_single_atom_residue(r) or is_water(r):
                continue
            for a in r.atoms:
                if a.atomic_number > 1:
                    coords.append((a.xx, a.xy, a.xz))
        try:
            import numpy as np
            from scipy.spatial import cKDTree
            arr = np.array(coords)
            tree = cKDTree(arr)
            t2 = threshold ** 2

            def rule(residue: pmd.Residue) -> bool:
                if not is_single_atom_residue(residue):
                    return False
                a = residue.atoms[0]
                d, _ = tree.query([a.xx, a.xy, a.xz], k=1)
                return d <= threshold
            return rule
        except ImportError:
            pass  # fall through to slow scan

    def rule(residue: pmd.Residue) -> bool:
        if not is_single_atom_residue(residue):
            return False
        if not coords:
            return False
        ax, ay, az = residue.atoms[0].xx, residue.atoms[0].xy, residue.atoms[0].xz
        t2 = threshold ** 2
        for cx, cy, cz in coords:
            if (ax-cx)**2 + (ay-cy)**2 + (az-cz)**2 <= t2:
                return True
        return False
    return rule


def by_residue_numbers(nums: Iterable[int]) -> ResidueRule:
    """Keep single-atom residue if its 1-based residue number is in `nums`."""
    nset: Set[int] = set(int(n) for n in nums)

    def rule(residue: pmd.Residue) -> bool:
        if not is_single_atom_residue(residue):
            return False
        return (residue.number + 1) in nset
    return rule


def by_residue_name(names: Iterable[str]) -> ResidueRule:
    """Keep single-atom residue if its name (stripped) matches any in `names`."""
    nset: Set[str] = set(n.strip() for n in names)

    def rule(residue: pmd.Residue) -> bool:
        if not is_single_atom_residue(residue):
            return False
        return residue.name.strip() in nset
    return rule


def combine(*rules: ResidueRule) -> ResidueRule:
    """Logical OR: keep if any of the rules says keep."""
    def rule(residue: pmd.Residue) -> bool:
        return any(r(residue) for r in rules)
    return rule


def intersect(*rules: ResidueRule) -> ResidueRule:
    """Logical AND: keep only if all of the rules say keep."""
    def rule(residue: pmd.Residue) -> bool:
        return all(r(residue) for r in rules)
    return rule


def negate(inner: ResidueRule) -> ResidueRule:
    """Logical NOT — flip the meaning of `inner`."""
    def rule(residue: pmd.Residue) -> bool:
        return not inner(residue)
    return rule


def keep_none(_residue: pmd.Residue) -> bool:
    """Strip every single-atom residue (no structural ions kept)."""
    return False


def keep_all_ions(_residue: pmd.Residue) -> bool:
    """Keep every single-atom residue (no bulk stripping). Use sparingly."""
    return True


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def strip_bulk_ions(structure: pmd.Structure,
                      *,
                      structural_rule: ResidueRule = None,
                      strip_waters: bool = True,
                      verbose: bool = False) -> dict:
    """Mutate `structure` in place to remove bulk ions and (optionally) waters.

    Returns an audit dict with counts and per-name breakdown.

    Parameters
    ----------
    structure : pmd.Structure
        Loaded prmtop+rst7 (mutated in place).
    structural_rule : Callable[[pmd.Residue], bool], optional
        Returns True for single-atom residues to keep. Default:
        within_distance_of_protein(4.0). Pass `keep_none` to strip all ions
        or `keep_all_ions` to keep all.
    strip_waters : bool
        If True (default), removes all water residues by name+atom-count.
    verbose : bool
        Print a one-line per-name summary of what was stripped vs kept.

    Notes
    -----
    Operates only on monoatomic residues + waters. Anything else
    (proteins, cofactors, lipids, NA bases) is untouched.
    """
    if structural_rule is None:
        structural_rule = within_distance_of_protein(4.0, structure=structure)

    # Audit before
    pre_counts: dict = {}
    for r in structure.residues:
        nm = r.name.strip()
        pre_counts[nm] = pre_counts.get(nm, 0) + 1

    # Decide what to strip. Waters are bulk-stripped by NAME mask
    # (one short string, fast to parse) — enumerating 50k water residue
    # numbers as ":1,2,...,50000" makes parmed's mask parser blow up.
    # Ion residues that fail the structural rule are enumerated by number,
    # since the rule can keep some and strip others within the same name.
    water_names_present: set = set()
    ion_strip_nums: list = []
    kept_ions: list = []
    n_waters_stripped = 0
    for r in structure.residues:
        if is_water(r) and strip_waters:
            water_names_present.add(r.name.strip())
            n_waters_stripped += 1
        elif is_single_atom_residue(r):
            if not structural_rule(r):
                ion_strip_nums.append(r.number + 1)
            else:
                kept_ions.append((r.number + 1, r.name.strip()))

    mask_parts = []
    if water_names_present:
        mask_parts.append(":" + ",".join(sorted(water_names_present)))
    if ion_strip_nums:
        mask_parts.append(":" + ",".join(str(n) for n in ion_strip_nums))
    if mask_parts:
        structure.strip("|".join(mask_parts))
    n_strip_total = n_waters_stripped + len(ion_strip_nums)

    # Audit after
    post_counts: dict = {}
    for r in structure.residues:
        nm = r.name.strip()
        post_counts[nm] = post_counts.get(nm, 0) + 1

    info = {
        "n_atoms_before":     sum(len(r.atoms) for r in structure.residues) + 0,
        "n_atoms_after":      len(structure.atoms),
        "n_residues_stripped": n_strip_total,
        "n_waters_stripped":  n_waters_stripped,
        "n_ions_stripped":    len(ion_strip_nums),
        "kept_structural_ions": kept_ions,
        "pre_counts":  pre_counts,
        "post_counts": post_counts,
    }

    if verbose:
        print(f"  strip_bulk_ions: removed {n_strip_total} residues "
              f"({n_waters_stripped} waters + {len(ion_strip_nums)} ions)")
        print(f"    kept structural ions: {kept_ions or '(none)'}")
        delta = {k: post_counts.get(k, 0) - pre_counts.get(k, 0)
                  for k in pre_counts if pre_counts[k] != post_counts.get(k, 0)}
        if delta:
            print(f"    Δ residue counts: {delta}")
    return info
