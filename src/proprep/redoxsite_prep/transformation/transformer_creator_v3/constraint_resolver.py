"""
Constraint resolver for role assignment.

Takes relationship phrases and resolves each role to a specific residue
in the current RedoxSite. Uses constraint propagation with implicit
elimination.

Algorithm:
  1. For each role, start with candidate set = all residues matching its resname
  2. Apply ALL constraints to narrow candidate sets
  3. Propagate: if any role has exactly 1 candidate, claim it
     → remove that residue from all other roles' candidate sets
  4. Repeat until stable
  5. Report: resolved, ambiguous, or conflicting
"""

import logging
import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set, Tuple

from .data_models import (
    ConnectivityPhrase,
    RelationshipPhrase,
    RoleDefinition,
    SequencePhrase,
    SpatialPhrase,
)

logger = logging.getLogger(__name__)


# Residue identity tuple: (chain, resid)
ResidueKey = Tuple[str, int]


@dataclass
class ResolutionResult:
    """Result of constraint resolution for all roles."""
    resolved: Dict[str, ResidueKey]         # role_label -> (chain, resid)
    ambiguous: Dict[str, List[ResidueKey]]  # role_label -> remaining candidates
    conflicting: List[str]                  # role_labels with 0 candidates

    @property
    def is_complete(self) -> bool:
        return len(self.ambiguous) == 0 and len(self.conflicting) == 0


class ConstraintResolver:
    """
    Resolve roles to specific residues using relationship constraints.

    Args:
        roles: dict of label -> RoleDefinition
        redox_site: the RedoxSite object (for atoms, bonds, coordinates)
    """

    def __init__(self, roles: Dict[str, RoleDefinition], redox_site):
        self.roles = roles
        self.site = redox_site
        # Build lookup structures
        self._atoms_by_residue: Dict[ResidueKey, List[Any]] = {}
        self._residue_resnames: Dict[ResidueKey, str] = {}
        self._build_lookups()

    def _build_lookups(self):
        """Build efficient lookup structures from the RedoxSite."""
        for atom in self.site.atoms:
            key = (atom.chain, atom.resid)
            if key not in self._atoms_by_residue:
                self._atoms_by_residue[key] = []
                self._residue_resnames[key] = atom.resname
            self._atoms_by_residue[key].append(atom)

    def resolve(self) -> ResolutionResult:
        """
        Run constraint resolution.

        Returns ResolutionResult with resolved, ambiguous, and conflicting roles.
        """
        # Step 1: Initialize candidate sets
        candidates: Dict[str, Set[ResidueKey]] = {}
        for label, role in self.roles.items():
            # Candidates are all residues matching the role's resname (or alt_resnames)
            valid_names = {role.resname} | set(role.alt_resnames)
            role_candidates = set()
            for rkey, resname in self._residue_resnames.items():
                if resname in valid_names:
                    role_candidates.add(rkey)
            candidates[label] = role_candidates

        # Auto-resolve centers (should be unique)
        for label, role in self.roles.items():
            if role.is_center:
                center_key = (role.chain, role.resid)
                candidates[label] = {center_key}

        # Step 2: Apply all constraints
        for label, role in self.roles.items():
            for phrase in role.phrases:
                candidates[label] = self._apply_constraint(
                    phrase, label, candidates
                )

        # Step 3: Propagate unique assignments
        candidates = self._propagate(candidates)

        # Step 4: Build result
        resolved = {}
        ambiguous = {}
        conflicting = []

        for label in self.roles:
            cset = candidates[label]
            if len(cset) == 1:
                resolved[label] = next(iter(cset))
            elif len(cset) == 0:
                conflicting.append(label)
            else:
                ambiguous[label] = sorted(cset)

        return ResolutionResult(
            resolved=resolved,
            ambiguous=ambiguous,
            conflicting=conflicting,
        )

    def _propagate(self, candidates: Dict[str, Set[ResidueKey]]) -> Dict[str, Set[ResidueKey]]:
        """
        Propagate unique assignments: if a role has exactly 1 candidate,
        remove that residue from all other roles' candidate sets.
        Repeat until stable.
        """
        changed = True
        while changed:
            changed = False
            for label, cset in candidates.items():
                if len(cset) == 1:
                    claimed = next(iter(cset))
                    for other_label, other_cset in candidates.items():
                        if other_label != label and claimed in other_cset:
                            other_cset.discard(claimed)
                            changed = True
        return candidates

    # ===== CONSTRAINT APPLICATION =====

    def _apply_constraint(
        self,
        phrase: RelationshipPhrase,
        role_label: str,
        candidates: Dict[str, Set[ResidueKey]],
    ) -> Set[ResidueKey]:
        """Apply a single constraint to narrow a role's candidate set."""
        current = candidates[role_label]
        if not current:
            return current

        if isinstance(phrase, ConnectivityPhrase):
            return self._apply_connectivity(phrase, current, candidates)
        elif isinstance(phrase, SpatialPhrase):
            return self._apply_spatial(phrase, current, candidates)
        elif isinstance(phrase, SequencePhrase):
            return self._apply_sequence(phrase, current, candidates)
        return current

    def _apply_connectivity(
        self,
        phrase: ConnectivityPhrase,
        current: Set[ResidueKey],
        candidates: Dict[str, Set[ResidueKey]],
    ) -> Set[ResidueKey]:
        """
        Apply connectivity constraint.
        Keep only candidates that have a matching bond to the target role.
        """
        target_candidates = candidates.get(phrase.target_role, set())
        if not target_candidates:
            return current  # Can't evaluate yet — target not resolved

        result = set()
        for rkey in current:
            if self._has_matching_bond(
                rkey, target_candidates,
                source_atom=phrase.source_atom,
                target_atom=phrase.target_atom,
                bond_type=phrase.bond_type,
            ):
                result.add(rkey)
        return result

    def _has_matching_bond(
        self,
        source_key: ResidueKey,
        target_keys: Set[ResidueKey],
        source_atom: Optional[str] = None,
        target_atom: Optional[str] = None,
        bond_type: Optional[str] = None,
    ) -> bool:
        """Check if source residue has a matching bond to any target residue."""
        for bond in self.site.bonds:
            # Get PDB info for both ends of the bond
            info1 = self.site.coord_to_pdb.get(bond.atom1_coords, {})
            info2 = self.site.coord_to_pdb.get(bond.atom2_coords, {})

            if not info1 or not info2:
                continue

            key1 = (info1.get("chain", ""), info1.get("resid", 0))
            key2 = (info2.get("chain", ""), info2.get("resid", 0))

            # Check both orientations
            for src_info, tgt_info, src_key, tgt_key in [
                (info1, info2, key1, key2),
                (info2, info1, key2, key1),
            ]:
                if src_key != source_key:
                    continue
                if tgt_key not in target_keys:
                    continue

                # Check bond type
                if bond_type and bond.chemical_type != bond_type:
                    continue

                # Check source atom
                if source_atom and src_info.get("atom_name") != source_atom:
                    continue

                # Check target atom
                if target_atom and tgt_info.get("atom_name") != target_atom:
                    continue

                return True

        return False

    def _apply_spatial(
        self,
        phrase: SpatialPhrase,
        current: Set[ResidueKey],
        candidates: Dict[str, Set[ResidueKey]],
    ) -> Set[ResidueKey]:
        """
        Apply spatial constraint.
        For closest/farthest: select the one candidate with min/max distance.
        For within: filter to candidates within the distance threshold.
        """
        target_candidates = candidates.get(phrase.target_role, set())
        if not target_candidates:
            return current

        # Compute distances from each candidate to the target
        distances: Dict[ResidueKey, float] = {}
        for rkey in current:
            dist = self._min_distance(
                rkey, target_candidates,
                source_atom=phrase.source_atom,
                target_atom=phrase.target_atom,
            )
            if dist is not None:
                distances[rkey] = dist

        if not distances:
            return current

        if phrase.comparison == "closest":
            min_dist = min(distances.values())
            return {k for k, d in distances.items() if math.isclose(d, min_dist, abs_tol=0.01)}
        elif phrase.comparison == "farthest":
            max_dist = max(distances.values())
            return {k for k, d in distances.items() if math.isclose(d, max_dist, abs_tol=0.01)}
        elif phrase.comparison == "within":
            return {k for k, d in distances.items() if d <= phrase.distance}

        return current

    def _min_distance(
        self,
        source_key: ResidueKey,
        target_keys: Set[ResidueKey],
        source_atom: Optional[str] = None,
        target_atom: Optional[str] = None,
    ) -> Optional[float]:
        """
        Compute minimum atom-atom distance between source residue and target residues.
        If source_atom/target_atom specified, only use those atoms.
        """
        source_atoms = self._atoms_by_residue.get(source_key, [])
        target_atoms = []
        for tkey in target_keys:
            target_atoms.extend(self._atoms_by_residue.get(tkey, []))

        if not source_atoms or not target_atoms:
            return None

        # Filter to specific atoms if requested
        if source_atom:
            source_atoms = [a for a in source_atoms if a.atom_name == source_atom]
        if target_atom:
            target_atoms = [a for a in target_atoms if a.atom_name == target_atom]

        if not source_atoms or not target_atoms:
            return None

        min_dist = float("inf")
        for sa in source_atoms:
            for ta in target_atoms:
                dx = sa.coords[0] - ta.coords[0]
                dy = sa.coords[1] - ta.coords[1]
                dz = sa.coords[2] - ta.coords[2]
                dist = math.sqrt(dx * dx + dy * dy + dz * dz)
                if dist < min_dist:
                    min_dist = dist

        return min_dist if min_dist < float("inf") else None

    def _apply_sequence(
        self,
        phrase: SequencePhrase,
        current: Set[ResidueKey],
        candidates: Dict[str, Set[ResidueKey]],
    ) -> Set[ResidueKey]:
        """Apply sequence constraint."""
        # Arithmetic: X is +N from Y
        if phrase.offset is not None and phrase.offset_target_role:
            target_cands = candidates.get(phrase.offset_target_role, set())
            result = set()
            for tkey in target_cands:
                target_chain, target_resid = tkey
                expected = (target_chain, target_resid + phrase.offset)
                if expected in current:
                    result.add(expected)
            return result if result else current  # Don't empty if target unresolved

        # Resid comparison: X has a lower/higher resid than Y
        if phrase.resid_comparison and phrase.comparison_target_role:
            target_cands = candidates.get(phrase.comparison_target_role, set())
            if not target_cands:
                return current
            result = set()
            for rkey in current:
                for tkey in target_cands:
                    if phrase.resid_comparison == "lower" and rkey[1] < tkey[1]:
                        result.add(rkey)
                    elif phrase.resid_comparison == "higher" and rkey[1] > tkey[1]:
                        result.add(rkey)
            return result if result else current

        # Chain relation: X is on the same/different chain as Y
        if phrase.chain_relation and phrase.chain_target_role:
            target_cands = candidates.get(phrase.chain_target_role, set())
            target_chains = {tkey[0] for tkey in target_cands}
            if not target_chains:
                return current
            if phrase.chain_relation == "same":
                return {rkey for rkey in current if rkey[0] in target_chains}
            elif phrase.chain_relation == "different":
                return {rkey for rkey in current if rkey[0] not in target_chains}

        return current
