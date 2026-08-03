"""
Link Atom Placer - QM/MM Boundary Link Atom Generation

Places hydrogen link atoms at bonds crossing the QM/MM boundary in ONIOM calculations.
Uses parmed bond topology to identify boundary bonds.
"""

import logging
from typing import Dict, List, Optional, Tuple
from rich.console import Console

from proprep.structure_prep.comprehensive_redox_detector import RedoxSite
from .data_structures import ONIOMLayer, LinkAtom, LayerAssignment
from .structure_utils import calculate_distance, normalize_vector


logger = logging.getLogger(__name__)

STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN",
}


class LinkAtomPlacer:
    """
    Places hydrogen link atoms at QM/MM boundaries in ONIOM setups.

    Uses parmed bond topology to identify bonds where one atom is in the
    HIGH (or MEDIUM) layer and the other is in the LOW layer.
    """

    def __init__(
        self,
        redox_site: RedoxSite,
        layer_assignments: Dict[int, LayerAssignment],
        boundary_atom_indices: Optional[List[int]] = None,
        parm=None,
        logger_instance=None,
        console=None,
    ):
        self.redox_site = redox_site
        self.layer_assignments = layer_assignments
        self.boundary_atom_indices = boundary_atom_indices or []
        self.parm = parm
        self.logger = logger_instance or logger
        self.console = console or Console()

        self.link_atoms: List[LinkAtom] = []

    def place_link_atoms(
        self,
        custom_scale_factors: Optional[Dict[str, float]] = None,
    ) -> List[LinkAtom]:
        """Identify QM/MM boundary bonds and place link atoms."""
        self.link_atoms.clear()

        scale_factors = self._get_default_scale_factors()
        if custom_scale_factors:
            scale_factors.update(custom_scale_factors)

        boundary_bonds = self._identify_boundary_bonds()

        if not boundary_bonds:
            self.logger.warning("No QM/MM boundary bonds found")
            return []

        self.logger.info(f"Found {len(boundary_bonds)} QM/MM boundary bonds")

        for qm_idx, mm_idx, bond_length in boundary_bonds:
            link_atom = self._place_single_link_atom(
                qm_idx, mm_idx, bond_length, scale_factors)
            if link_atom:
                self.link_atoms.append(link_atom)
                qm_a = self.parm.atoms[qm_idx]
                mm_a = self.parm.atoms[mm_idx]
                self.logger.info(
                    f"Link atom #{len(self.link_atoms)}: "
                    f"{qm_a.residue.name}{qm_a.residue.number}:{qm_a.name} "
                    f"<-> {mm_a.residue.name}{mm_a.residue.number}:{mm_a.name} "
                    f"| {link_atom.boundary_type}"
                )

        self.logger.info(f"Successfully placed {len(self.link_atoms)} link atoms")

        for warning in self._validate_link_atoms():
            self.logger.warning(warning)

        return self.link_atoms

    def _identify_boundary_bonds(self) -> List[Tuple[int, int, float]]:
        """
        Identify bonds crossing QM/MM boundary from parmed topology.

        Only considers bonds involving designated boundary atoms (CA atoms
        identified by _detect_boundary_atoms). Link atoms are placed on
        CA-N or CA-C backbone bonds, not on arbitrary HIGH-LOW crossings
        within cofactors or sidechains.

        Returns list of (qm_atom_idx, mm_atom_idx, bond_eq_length).
        """
        if self.parm is None:
            self.logger.error("No parmed structure available for bond detection")
            return []

        boundary_set = set(self.boundary_atom_indices)
        if not boundary_set:
            self.logger.warning("No boundary atoms designated")
            return []

        boundary = []

        for bond in self.parm.bonds:
            idx1 = bond.atom1.idx
            idx2 = bond.atom2.idx

            # At least one atom must be a designated boundary atom
            if idx1 not in boundary_set and idx2 not in boundary_set:
                continue

            if idx1 not in self.layer_assignments or idx2 not in self.layer_assignments:
                continue

            layer1 = self.layer_assignments[idx1].layer
            layer2 = self.layer_assignments[idx2].layer

            eq_len = bond.type.req if bond.type else calculate_distance(
                (bond.atom1.xx, bond.atom1.xy, bond.atom1.xz),
                (bond.atom2.xx, bond.atom2.xy, bond.atom2.xz),
            )

            if layer1 in (ONIOMLayer.HIGH, ONIOMLayer.MEDIUM) and layer2 == ONIOMLayer.LOW:
                boundary.append((idx1, idx2, eq_len))
            elif layer2 in (ONIOMLayer.HIGH, ONIOMLayer.MEDIUM) and layer1 == ONIOMLayer.LOW:
                boundary.append((idx2, idx1, eq_len))

        return boundary

    def _place_single_link_atom(
        self,
        qm_idx: int,
        mm_idx: int,
        bond_length: float,
        scale_factors: Dict[str, float],
    ) -> Optional[LinkAtom]:
        """Place a single link atom along the QM→MM bond vector."""
        qm_atom = self.parm.atoms[qm_idx]
        mm_atom = self.parm.atoms[mm_idx]

        qm_coords = (qm_atom.xx, qm_atom.xy, qm_atom.xz)
        mm_coords = (mm_atom.xx, mm_atom.xy, mm_atom.xz)

        qm_element = qm_atom.element_name if hasattr(qm_atom, 'element_name') else qm_atom.name[0]
        mm_element = mm_atom.element_name if hasattr(mm_atom, 'element_name') else mm_atom.name[0]

        # Bond vector QM → MM
        bond_vector = (
            mm_coords[0] - qm_coords[0],
            mm_coords[1] - qm_coords[1],
            mm_coords[2] - qm_coords[2],
        )

        try:
            unit_vector = normalize_vector(bond_vector)
        except ValueError as e:
            self.logger.error(f"Failed to normalize bond vector: {e}")
            return None

        # Determine scale factor
        bond_type_key = f"{qm_element}-{mm_element}"
        reverse_key = f"{mm_element}-{qm_element}"
        scale_factor = scale_factors.get(
            bond_type_key,
            scale_factors.get(reverse_key, 0.723)
        )

        # Link atom position
        link_coords = (
            qm_coords[0] + scale_factor * unit_vector[0],
            qm_coords[1] + scale_factor * unit_vector[1],
            qm_coords[2] + scale_factor * unit_vector[2],
        )

        # Classify boundary type
        boundary_type = self._classify_boundary(qm_atom, mm_atom)

        return LinkAtom(
            qm_parent_idx=qm_idx,
            mm_parent_idx=mm_idx,
            coords=link_coords,
            qm_parent_coords=qm_coords,
            mm_parent_coords=mm_coords,
            bond_vector=unit_vector,
            bond_length=bond_length,
            element="H",
            scale_factor=scale_factor,
            boundary_type=boundary_type,
        )

    def _classify_boundary(self, qm_atom, mm_atom) -> str:
        """
        Classify the QM/MM boundary by which atom is on the LOW (mm) side
        and which is on the HIGH (qm) side. Labels match Gaussian's
        recommended fragment caps (Clemente, Apr 2026).
        """
        qm_name, mm_name = qm_atom.name, mm_atom.name

        if qm_name == "CA" and mm_name == "N":
            return "Acetamide cap (N-CA)"
        if qm_name == "CA" and mm_name == "C":
            return "N-methylamide cap (CA-C)"
        if qm_name == "CA" and mm_name == "CB":
            # Second cut of an acetamide / N-methylamide cap (CA in HIGH).
            return "Methyl trim (CA-CB)"
        if qm_name == "CB" and mm_name == "CA":
            # Sidechain-only mode: CA stays LOW.
            return "Sidechain cut (CA-CB)"

        qm_is_protein = qm_atom.residue.name in STANDARD_AMINO_ACIDS
        mm_is_protein = mm_atom.residue.name in STANDARD_AMINO_ACIDS
        if qm_is_protein != mm_is_protein:
            return "Ligand-protein bond"

        return "Other"

    def _validate_link_atoms(self) -> List[str]:
        """Validate quality of link atom placements."""
        warnings = []

        for la in self.link_atoms:
            dist_to_qm = calculate_distance(la.coords, la.qm_parent_coords)
            if dist_to_qm < 0.5:
                warnings.append(
                    f"Link atom too close to QM parent: {dist_to_qm:.2f}A (expected ~0.7-1.0A)")
            elif dist_to_qm > 1.5:
                warnings.append(
                    f"Link atom too far from QM parent: {dist_to_qm:.2f}A (expected ~0.7-1.0A)")

        for i in range(len(self.link_atoms)):
            for j in range(i + 1, len(self.link_atoms)):
                dist = calculate_distance(
                    self.link_atoms[i].coords, self.link_atoms[j].coords)
                if dist < 0.01:
                    warnings.append(f"Duplicate link atoms at {self.link_atoms[i].coords}")

        return warnings

    def _get_default_scale_factors(self) -> Dict[str, float]:
        """Default scale factors for common bond types."""
        return {
            "C-C": 0.723,
            "C-N": 0.70,
            "C-O": 0.68,
            "C-S": 0.73,
            "N-C": 0.70,
            "N-CA": 0.70,
            "CA-CB": 0.723,
            "CB-CG": 0.723,
            "S-C": 0.73,
            "O-C": 0.68,
        }
