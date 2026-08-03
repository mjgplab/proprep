"""
Connectivity Builder - Build Gaussian connectivity table from parmed topology.

Reads bonds directly from the AMBER prmtop (via parmed), which includes
all covalent bonds including those from RedoxSite bond directives in the
tLEaP input script. Only link atom bonds need to be added separately.
"""

import logging
from typing import Dict, List, Set, Tuple

from .data_structures import ONIOMSetup, ConnectivityEntry


logger = logging.getLogger(__name__)


class ConnectivityBuilder:
    """Build Gaussian ONIOM connectivity table from parmed bond topology."""

    def __init__(self, oniom_setup: ONIOMSetup, parm, logger_instance=None):
        self.oniom_setup = oniom_setup
        self.parm = parm
        self.logger = logger_instance or logger

        self.bonds: Set[Tuple[int, int]] = set()
        self.stats = {
            'prmtop_bonds': 0,
            'link_bonds': 0,
        }

    def build_connectivity(self) -> List[ConnectivityEntry]:
        """
        Build complete connectivity table.

        1. Read all bonds from prmtop (includes RedoxSite coordination bonds)
        2. Add link atom bonds
        3. Convert to Gaussian 1-based ConnectivityEntry format
        """
        # Step 1: All bonds from parmed topology
        # Exclude H-H bonds in TIP3P water (AMBER uses these for SHAKE
        # constraints, but they are not real covalent bonds and cause
        # Gaussian to infer spurious OW-HW-HW angles).
        for bond in self.parm.bonds:
            idx1 = bond.atom1.idx
            idx2 = bond.atom2.idx
            if idx1 in self.oniom_setup.layer_assignments and \
               idx2 in self.oniom_setup.layer_assignments:
                # Skip HW-HW bonds (TIP3P water SHAKE constraint)
                if bond.atom1.type == 'HW' and bond.atom2.type == 'HW':
                    continue
                self._add_bond(idx1, idx2)
                self.stats['prmtop_bonds'] += 1

        # Step 2: Link atom bonds
        self.stats['link_bonds'] = len(self.oniom_setup.link_atoms)

        self.logger.info(
            f"Connectivity: {len(self.bonds)} bonds "
            f"({self.stats['prmtop_bonds']} prmtop, "
            f"{self.stats['link_bonds']} link)"
        )

        # Step 3: Convert to ConnectivityEntry format
        return self._build_connectivity_table()

    def _add_bond(self, idx1: int, idx2: int):
        """Add a bond (canonical order: smaller index first)."""
        self.bonds.add((min(idx1, idx2), max(idx1, idx2)))

    def _build_connectivity_table(self) -> List[ConnectivityEntry]:
        """Convert bonds to Gaussian 1-based ConnectivityEntry list.

        Atoms ordered by their position in layer_assignments iteration
        order (parmed atom order). Link atoms come after real atoms.
        """
        # atom_idx → 1-based Gaussian index
        idx_to_gaussian = {}
        for gaussian_idx, atom_idx in enumerate(
            self.oniom_setup.layer_assignments.keys(), start=1
        ):
            idx_to_gaussian[atom_idx] = gaussian_idx

        # Build adjacency list
        adjacency: Dict[int, List[int]] = {}
        for idx1, idx2 in self.bonds:
            adjacency.setdefault(idx1, []).append(idx2)
            adjacency.setdefault(idx2, []).append(idx1)

        # ConnectivityEntry for each real atom
        entries = []
        for atom_idx in self.oniom_setup.layer_assignments:
            gaussian_idx = idx_to_gaussian[atom_idx]
            neighbors = adjacency.get(atom_idx, [])

            # Upper-triangular: only connections to higher Gaussian indices
            connected = []
            orders = []
            for neighbor_idx in sorted(neighbors):
                g_neighbor = idx_to_gaussian.get(neighbor_idx)
                if g_neighbor and g_neighbor > gaussian_idx:
                    connected.append(g_neighbor)
                    orders.append(1.0)

            atom = self.parm.atoms[atom_idx]
            entries.append(ConnectivityEntry(
                atom_index=gaussian_idx,
                connected_indices=connected,
                bond_orders=orders,
                element=atom.element_name if hasattr(atom, 'element_name') else atom.name[0],
            ))

        # Link atoms are NOT separate entries in Gaussian ONIOM connectivity.
        # They appear on the same line as their QM parent in the atom section.

        return entries
