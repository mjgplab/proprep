"""
Layer Classifier - ONIOM Layer Assignment Logic

Assigns atoms from a RedoxSite to ONIOM layers (HIGH/MEDIUM/LOW) based on
user-selected residues with intelligent boundary handling.

Features:
- Residue-based layer selection
- Alpha carbon bridging to prevent peptide bond cutting
- Fragment merging for contiguous QM regions
- Distance-based freezing of distant atoms

© 2024 ProPrep Developer. All rights reserved.
"""

import logging
from typing import Dict, List, Optional, Set, Tuple
from rich.console import Console

from proprep.structure_prep.comprehensive_redox_detector import RedoxSite, RedoxSiteAtom
from .data_structures import ONIOMLayer, FreezeFlag, LayerAssignment
from .structure_utils import calculate_distance, calculate_geometric_center


logger = logging.getLogger(__name__)


class LayerClassifier:
    """
    Assigns ONIOM layer designations to atoms in a RedoxSite.

    Strategy:
    1. User explicitly selects residues for HIGH and/or MEDIUM layers
    2. Add alpha carbons from neighboring residues (X-1 and X+1)
    3. Merge fragments if alpha carbon bridges overlap
    4. Assign freeze flags based on distance or explicit selection
    5. All remaining atoms default to LOW layer
    """

    def __init__(self, redox_site: RedoxSite, logger_instance=None, console=None):
        """
        Initialize layer classifier.

        Args:
            redox_site: RedoxSite containing atoms and bonds
            logger_instance: Optional logger instance
            console: Optional Rich Console for output
        """
        self.redox_site = redox_site
        self.logger = logger_instance or logger
        self.console = console or Console()

        # Storage for results
        self.layer_assignments: Dict[Tuple[float, float, float], LayerAssignment] = {}
        self.residue_to_layer: Dict[Tuple[str, int, str], ONIOMLayer] = {}

    def assign_layers(
        self,
        high_residues: List[Tuple[str, int, str]],
        medium_residues: Optional[List[Tuple[str, int, str]]] = None,
        freeze_distance_cutoff: Optional[float] = None,
        frozen_residues: Optional[List[Tuple[str, int, str]]] = None,
        n_layers: int = 2
    ) -> Dict[Tuple[float, float, float], LayerAssignment]:
        """
        Assign layers to all atoms in RedoxSite.

        Args:
            high_residues: List of (chain, resid, insertion_code) for HIGH layer
            medium_residues: List of residues for MEDIUM layer (3-layer ONIOM only)
            freeze_distance_cutoff: Distance (Å) beyond which atoms are frozen
            frozen_residues: Explicitly frozen residues (overrides distance cutoff)
            n_layers: 2 or 3 (number of ONIOM layers)

        Returns:
            Dict mapping coordinates to LayerAssignment objects

        Raises:
            ValueError: If invalid inputs
        """
        # ========================================
        # STEP 1: Validate inputs
        # ========================================
        if n_layers not in [2, 3]:
            raise ValueError("n_layers must be 2 or 3")

        if n_layers == 2 and medium_residues is not None and len(medium_residues) > 0:
            raise ValueError("Cannot specify medium_residues for 2-layer ONIOM")

        if medium_residues is None:
            medium_residues = []

        if frozen_residues is None:
            frozen_residues = []

        # Validate all specified residues exist in RedoxSite
        all_specified = high_residues + medium_residues + frozen_residues
        for residue_key in all_specified:
            atoms = self.redox_site.get_atoms_by_residue(*residue_key)
            if len(atoms) == 0:
                raise ValueError(f"Residue {residue_key} not found in RedoxSite")

        self.logger.info(f"Starting layer assignment for {len(all_specified)} specified residues")

        # ========================================
        # STEP 2: Build initial layer assignments (user selections)
        # ========================================
        for residue_key in high_residues:
            self._assign_residue_to_layer(residue_key, ONIOMLayer.HIGH, "User selected")

        for residue_key in medium_residues:
            self._assign_residue_to_layer(residue_key, ONIOMLayer.MEDIUM, "User selected")

        self.logger.info(f"Initial assignments: {len(high_residues)} HIGH, {len(medium_residues)} MEDIUM")

        # ========================================
        # STEP 3: Add alpha carbon bridges
        # ========================================
        expanded_high = self._add_alpha_carbon_bridges(high_residues, ONIOMLayer.HIGH)
        expanded_medium = self._add_alpha_carbon_bridges(medium_residues, ONIOMLayer.MEDIUM)

        self.logger.info(f"After CA bridging: HIGH expanded to {len(expanded_high)} residues")
        if expanded_medium:
            self.logger.info(f"After CA bridging: MEDIUM expanded to {len(expanded_medium)} residues")

        # ========================================
        # STEP 4: Merge overlapping fragments
        # ========================================
        merged_high = self._merge_fragments(expanded_high, ONIOMLayer.HIGH)
        merged_medium = self._merge_fragments(expanded_medium, ONIOMLayer.MEDIUM)

        self.logger.info(f"After fragment merging: HIGH has {len(merged_high)} residues")
        if merged_medium:
            self.logger.info(f"After fragment merging: MEDIUM has {len(merged_medium)} residues")

        # ========================================
        # STEP 5: Assign LOW layer to remaining atoms
        # ========================================
        for atom in self.redox_site.atoms:
            if atom.coords not in self.layer_assignments:
                residue_key = (atom.chain, atom.resid, atom.insertion_code)
                assignment = LayerAssignment(
                    coords=atom.coords,
                    layer=ONIOMLayer.LOW,
                    freeze=FreezeFlag.ACTIVE,
                    residue_key=residue_key,
                    atom_name=atom.atom_name,
                    element=atom.element,
                    assignment_reason="Default LOW layer"
                )
                self.layer_assignments[atom.coords] = assignment

        self.logger.info(f"Total atoms assigned: {len(self.layer_assignments)}")

        # ========================================
        # STEP 6: Assign freeze flags
        # ========================================
        if frozen_residues:
            for residue_key in frozen_residues:
                self._freeze_residue(residue_key, "User frozen")

        if freeze_distance_cutoff is not None:
            self._apply_distance_freezing(freeze_distance_cutoff, merged_high)

        self.logger.info("Freeze assignments complete")

        # ========================================
        # STEP 7: Validation and return
        # ========================================
        self._validate_layer_assignments()

        return self.layer_assignments

    def _assign_residue_to_layer(
        self,
        residue_key: Tuple[str, int, str],
        layer: ONIOMLayer,
        reason: str
    ):
        """
        Assign all atoms in a residue to a specific layer.

        Args:
            residue_key: (chain, resid, insertion_code)
            layer: ONIOMLayer to assign
            reason: String explaining assignment reason
        """
        atoms = self.redox_site.get_atoms_by_residue(*residue_key)

        if len(atoms) == 0:
            self.logger.warning(f"Residue {residue_key} has no atoms")
            return

        for atom in atoms:
            assignment = LayerAssignment(
                coords=atom.coords,
                layer=layer,
                freeze=FreezeFlag.ACTIVE,
                residue_key=residue_key,
                atom_name=atom.atom_name,
                element=atom.element,
                assignment_reason=reason
            )
            self.layer_assignments[atom.coords] = assignment

        self.residue_to_layer[residue_key] = layer
        self.logger.debug(f"Assigned {len(atoms)} atoms from {residue_key} to {layer.value}")

    def _add_alpha_carbon_bridges(
        self,
        residue_keys: List[Tuple[str, int, str]],
        layer: ONIOMLayer
    ) -> List[Tuple[str, int, str]]:
        """
        Add alpha carbons from X-1 and X+1 residues to prevent peptide bond cutting.

        Args:
            residue_keys: Initial residue list
            layer: Layer these residues belong to

        Returns:
            Expanded list including bridging residues
        """
        expanded_residues = set(residue_keys)

        for residue_key in residue_keys:
            chain, resid, ins_code = residue_key

            # Try X-1 residue (previous in sequence)
            prev_resid = resid - 1
            prev_key = (chain, prev_resid, ins_code)
            prev_atoms = self.redox_site.get_atoms_by_residue(*prev_key)

            if prev_atoms:
                ca_atom = self._find_ca_atom(prev_atoms)
                if ca_atom and ca_atom.coords not in self.layer_assignments:
                    # Add ONLY the CA atom from X-1
                    assignment = LayerAssignment(
                        coords=ca_atom.coords,
                        layer=layer,
                        freeze=FreezeFlag.ACTIVE,
                        residue_key=prev_key,
                        atom_name="CA",
                        element=ca_atom.element,
                        assignment_reason=f"Alpha carbon bridge for {residue_key}"
                    )
                    self.layer_assignments[ca_atom.coords] = assignment
                    self.logger.debug(f"Added CA bridge from {prev_key}")

            # Try X+1 residue (next in sequence)
            next_resid = resid + 1
            next_key = (chain, next_resid, ins_code)
            next_atoms = self.redox_site.get_atoms_by_residue(*next_key)

            if next_atoms:
                ca_atom = self._find_ca_atom(next_atoms)
                if ca_atom and ca_atom.coords not in self.layer_assignments:
                    assignment = LayerAssignment(
                        coords=ca_atom.coords,
                        layer=layer,
                        freeze=FreezeFlag.ACTIVE,
                        residue_key=next_key,
                        atom_name="CA",
                        element=ca_atom.element,
                        assignment_reason=f"Alpha carbon bridge for {residue_key}"
                    )
                    self.layer_assignments[ca_atom.coords] = assignment
                    self.logger.debug(f"Added CA bridge from {next_key}")

        return list(expanded_residues)

    def _merge_fragments(
        self,
        residue_keys: List[Tuple[str, int, str]],
        layer: ONIOMLayer
    ) -> List[Tuple[str, int, str]]:
        """
        Merge overlapping fragments by including intervening residues.

        Args:
            residue_keys: Residue list (may have gaps)
            layer: Layer these residues belong to

        Returns:
            Merged list with intervening residues included
        """
        if not residue_keys:
            return []

        # Group residues by chain
        chain_groups = {}
        for residue_key in residue_keys:
            chain, resid, ins_code = residue_key
            if chain not in chain_groups:
                chain_groups[chain] = []
            chain_groups[chain].append((resid, ins_code))

        merged_residues = []

        # Process each chain independently
        for chain, resid_list in chain_groups.items():
            # Sort by resid
            sorted_resids = sorted(resid_list, key=lambda x: x[0])

            # Build fragments
            current_fragment = [sorted_resids[0]]

            for i in range(1, len(sorted_resids)):
                prev_resid, prev_ins = current_fragment[-1]
                curr_resid, curr_ins = sorted_resids[i]

                gap = curr_resid - prev_resid

                if gap == 1:
                    # Contiguous
                    current_fragment.append((curr_resid, curr_ins))

                elif gap == 2:
                    # Gap of 1 residue - include middle
                    middle_resid = prev_resid + 1
                    middle_key = (chain, middle_resid, prev_ins)

                    middle_atoms = self.redox_site.get_atoms_by_residue(*middle_key)
                    if middle_atoms:
                        current_fragment.append((middle_resid, prev_ins))
                        self._assign_residue_to_layer(middle_key, layer, "Fragment merge (1-residue gap)")

                    current_fragment.append((curr_resid, curr_ins))

                elif gap == 3:
                    # Gap of 2 residues - include both
                    for offset in [1, 2]:
                        middle_resid = prev_resid + offset
                        middle_key = (chain, middle_resid, prev_ins)

                        middle_atoms = self.redox_site.get_atoms_by_residue(*middle_key)
                        if middle_atoms:
                            current_fragment.append((middle_resid, prev_ins))
                            self._assign_residue_to_layer(middle_key, layer, "Fragment merge (2-residue gap)")

                    current_fragment.append((curr_resid, curr_ins))

                else:
                    # Gap >= 3 - new fragment
                    # Add current fragment to merged list
                    for resid, ins in current_fragment:
                        merged_residues.append((chain, resid, ins))
                    # Start new fragment
                    current_fragment = [(curr_resid, curr_ins)]

            # Add final fragment
            for resid, ins in current_fragment:
                merged_residues.append((chain, resid, ins))

        self.logger.debug(f"Fragment merging: {len(residue_keys)} -> {len(merged_residues)} residues")

        return merged_residues

    def _find_ca_atom(self, atoms: List[RedoxSiteAtom]) -> Optional[RedoxSiteAtom]:
        """
        Find the CA (alpha carbon) atom in a residue.

        Args:
            atoms: List of atoms in the residue

        Returns:
            RedoxSiteAtom for CA, or None if not found
        """
        for atom in atoms:
            if atom.atom_name == "CA":
                return atom
        return None

    def _apply_distance_freezing(
        self,
        cutoff: float,
        active_region_residues: List[Tuple[str, int, str]]
    ):
        """
        Freeze atoms beyond a distance cutoff from the active region.

        Args:
            cutoff: Distance in Angstroms
            active_region_residues: Residues in the active HIGH layer
        """
        # Calculate active region center
        active_coords = []
        for residue_key in active_region_residues:
            atoms = self.redox_site.get_atoms_by_residue(*residue_key)
            for atom in atoms:
                active_coords.append(atom.coords)

        if not active_coords:
            self.logger.warning("No atoms in active region, cannot apply distance freezing")
            return

        center = calculate_geometric_center(active_coords)
        self.logger.info(f"Active region center: {center}")

        # Freeze atoms beyond cutoff
        n_frozen = 0
        for coords, assignment in self.layer_assignments.items():
            # Never freeze HIGH or MEDIUM layer atoms
            if assignment.layer in [ONIOMLayer.HIGH, ONIOMLayer.MEDIUM]:
                continue

            dist = calculate_distance(coords, center)

            if dist > cutoff:
                assignment.freeze = FreezeFlag.FROZEN
                assignment.assignment_reason += f" (Frozen: {dist:.1f}Å from center)"
                n_frozen += 1

        self.logger.info(f"Froze {n_frozen} atoms beyond {cutoff}Å cutoff")

    def _freeze_residue(self, residue_key: Tuple[str, int, str], reason: str):
        """
        Explicitly freeze all atoms in a residue.

        Args:
            residue_key: (chain, resid, insertion_code)
            reason: Explanation for freezing
        """
        atoms = self.redox_site.get_atoms_by_residue(*residue_key)

        for atom in atoms:
            if atom.coords in self.layer_assignments:
                self.layer_assignments[atom.coords].freeze = FreezeFlag.FROZEN
                self.layer_assignments[atom.coords].assignment_reason += f" ({reason})"

        self.logger.debug(f"Froze {len(atoms)} atoms in {residue_key}")

    def _validate_layer_assignments(self):
        """Validate layer assignments for consistency."""
        n_site_atoms = len(self.redox_site.atoms)
        n_assigned = len(self.layer_assignments)

        if n_assigned != n_site_atoms:
            self.logger.warning(
                f"Assignment mismatch: {n_site_atoms} atoms in site, {n_assigned} assigned"
            )

        # Count layers
        n_high = sum(1 for a in self.layer_assignments.values() if a.layer == ONIOMLayer.HIGH)
        n_medium = sum(1 for a in self.layer_assignments.values() if a.layer == ONIOMLayer.MEDIUM)
        n_low = sum(1 for a in self.layer_assignments.values() if a.layer == ONIOMLayer.LOW)

        self.logger.info(f"Layer distribution: HIGH={n_high}, MEDIUM={n_medium}, LOW={n_low}")

        # Warnings
        if n_high < 5:
            self.logger.warning("HIGH layer has very few atoms (<5), check residue selection")

        if n_high > 200:
            self.logger.warning("HIGH layer has many atoms (>200), may be computationally expensive")
