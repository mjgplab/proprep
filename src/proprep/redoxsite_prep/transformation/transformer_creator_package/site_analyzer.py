#!/usr/bin/env python3
"""
Site Analyzer

Extracts EXACT transformer requirements from RedoxSite objects in the workspace.
This provides precise, structure-based specifications rather than manual entry.
"""

import logging
from typing import Dict, List, Tuple, Any
from collections import Counter

from .data_models import (
    SiteRequirements,
    CenterRequirements,
    AtomRequirements,
    ResidueRequirement,
    BondRequirements,
    BondGroup,
    CenterTypeEnum
)

logger = logging.getLogger(__name__)


class SiteAnalyzer:
    """Analyze RedoxSite objects to generate exact transformer requirements"""

    def analyze_redox_site(self, redox_site) -> SiteRequirements:
        """
        Extract EXACT requirements from a RedoxSite object.

        No approximations - if the site has 2 HIS, requirement is min=2, max=2.
        Bond patterns are extracted in full detail, not simplified.

        Args:
            redox_site: RedoxSite object from comprehensive_redox_detector

        Returns:
            SiteRequirements with exact specifications
        """
        logger.info(f"Analyzing RedoxSite: {redox_site.site_id}")

        # Extract exact center requirements
        center_reqs = self._analyze_centers(redox_site.centers)

        # Extract exact residue requirements
        atom_reqs = self._analyze_residues(redox_site)

        # Extract complete bond patterns
        bond_reqs = self._analyze_bonds(redox_site.bonds, redox_site)

        return SiteRequirements(
            centers=center_reqs,
            atoms=atom_reqs,
            bonds=bond_reqs
        )

    def _analyze_centers(self, centers: List) -> CenterRequirements:
        """
        Extract EXACT center requirements.

        Count: Exactly how many centers exist
        Types: What CenterTypes are present
        Elements: What specific elements (for metal ions)
        """
        if not centers:
            # Fallback: assume 1 metal ion
            logger.warning("No centers found, defaulting to 1 METAL_ION")
            return CenterRequirements(
                required_count=1,
                center_types=[CenterTypeEnum.METAL_ION],
                elements=None
            )

        # Exact count
        required_count = len(centers)

        # Collect center types (handle both enum and string)
        center_types_set = set()
        for center in centers:
            if hasattr(center, 'center_type'):
                ct = center.center_type
                # Convert to enum if needed
                if hasattr(ct, 'value'):
                    ct_value = ct.value
                else:
                    ct_value = str(ct)

                # Map to our enum
                if ct_value == "metal_ion":
                    center_types_set.add(CenterTypeEnum.METAL_ION)
                elif ct_value == "organometallic_cofactor":
                    center_types_set.add(CenterTypeEnum.ORGANOMETALLIC_COFACTOR)
                elif ct_value == "organic_cofactor":
                    center_types_set.add(CenterTypeEnum.ORGANIC_COFACTOR)
                elif ct_value == "redox_amino_acid":
                    center_types_set.add(CenterTypeEnum.REDOX_AMINO_ACID)

        center_types = list(center_types_set) if center_types_set else [CenterTypeEnum.METAL_ION]

        # Collect elements (for metal ions)
        elements_set = set()
        for center in centers:
            if hasattr(center, 'element') and center.element:
                elements_set.add(center.element.upper())

        elements = list(elements_set) if elements_set else None

        logger.info(f"Centers: count={required_count}, types={[ct.value for ct in center_types]}, elements={elements}")

        return CenterRequirements(
            required_count=required_count,
            center_types=center_types,
            elements=elements
        )

    def _analyze_residues(self, redox_site) -> AtomRequirements:
        """
        Extract EXACT residue requirements.

        For each residue type in the site, set min_count = max_count = actual_count.
        This ensures the transformer matches EXACTLY this site structure.
        """
        # Count residues by name using residue_groups
        residue_counts = Counter()

        for (chain, resid, insertion_code), atom_coords in redox_site.residue_groups.items():
            # Get residue name from first atom in this residue
            if atom_coords:
                first_atom_coord = atom_coords[0]
                atom_info = redox_site.coord_to_pdb.get(first_atom_coord)
                if atom_info and 'resname' in atom_info:
                    resname = atom_info['resname']
                    residue_counts[resname] += 1

        # Build exact requirements (min = max = count)
        required_residues = {}
        for resname, count in residue_counts.items():
            required_residues[resname] = ResidueRequirement(
                min_count=count,
                max_count=count  # Exact match required
            )

        logger.info(f"Residues: {dict(residue_counts)}")

        return AtomRequirements(
            required_residues=required_residues,
            alternative_groups=None  # Could be inferred later if needed
        )

    def _analyze_bonds(self, bonds: List, redox_site) -> BondRequirements:
        """
        Extract COMPLETE bond patterns.

        Do NOT simplify! Extract the full bond structure:
        - All coordinate bonds with specific atom pairs
        - All covalent bonds with specific atom pairs
        - Disulfide bonds
        - Metal-metal bonds

        Organize into bond groups matching the transformer framework format.
        """
        if not bonds:
            logger.warning("No bonds to analyze")
            return BondRequirements(
                required_bond_groups=[],
                require_one_group=False
            )

        # Categorize bonds by chemical type
        coordinate_bonds = []
        covalent_bonds = []
        disulfide_bonds = []
        metal_metal_bonds = []

        for bond in bonds:
            chemical_type = bond.chemical_type

            # Extract atom pair information
            atom1_info = bond.atom1_residue_info
            atom2_info = bond.atom2_residue_info

            # Create bond tuple: ((resname1, atom1), (resname2, atom2))
            bond_tuple = (
                (atom1_info.get('resname', ''), atom1_info.get('atom_name', '')),
                (atom2_info.get('resname', ''), atom2_info.get('atom_name', ''))
            )

            if chemical_type == "coordinate":
                coordinate_bonds.append(bond_tuple)
            elif chemical_type == "covalent":
                covalent_bonds.append(bond_tuple)
            elif chemical_type == "disulfide":
                disulfide_bonds.append(bond_tuple)
            elif chemical_type == "metal-metal":
                metal_metal_bonds.append(bond_tuple)

        # Build bond group
        # All bonds must be present (min_count = total count)
        total_bonds = len(coordinate_bonds) + len(covalent_bonds) + len(disulfide_bonds) + len(metal_metal_bonds)

        bond_types = {}
        if coordinate_bonds:
            bond_types["coordinate"] = coordinate_bonds
        if covalent_bonds:
            bond_types["covalent"] = covalent_bonds
        if disulfide_bonds:
            bond_types["disulfide"] = disulfide_bonds
        if metal_metal_bonds:
            bond_types["metal-metal"] = metal_metal_bonds

        bond_group = BondGroup(
            description=f"Exact bond pattern from {redox_site.site_id}",
            min_count=total_bonds,  # All bonds must match
            bond_types=bond_types
        )

        logger.info(f"Bonds: {total_bonds} total "
                   f"(coordinate={len(coordinate_bonds)}, "
                   f"covalent={len(covalent_bonds)}, "
                   f"disulfide={len(disulfide_bonds)}, "
                   f"metal-metal={len(metal_metal_bonds)})")

        return BondRequirements(
            required_bond_groups=[bond_group],
            require_one_group=True  # This group must match
        )

    def describe_site(self, redox_site) -> str:
        """
        Generate a human-readable description of a RedoxSite.

        Used for displaying site list to user.
        """
        # Count components
        n_centers = len(redox_site.centers)
        n_atoms = len(redox_site.atoms)
        n_bonds = len(redox_site.bonds)

        # Get center info
        center_desc = "No centers"
        if redox_site.centers:
            center = redox_site.centers[0]
            element = getattr(center, 'element', '?')
            resname = getattr(center, 'resname', '?')
            center_desc = f"{element} ({resname})"
            if n_centers > 1:
                center_desc += f" + {n_centers-1} more"

        # Get residue types
        residue_types = set()
        for (chain, resid, insertion), coords in redox_site.residue_groups.items():
            if coords:
                coord = coords[0]
                atom_info = redox_site.coord_to_pdb.get(coord)
                if atom_info:
                    residue_types.add(atom_info.get('resname', '?'))

        residues_desc = ", ".join(sorted(residue_types)[:3])
        if len(residue_types) > 3:
            residues_desc += f" + {len(residue_types)-3} more"

        return f"{center_desc} | Residues: {residues_desc} | {n_atoms} atoms, {n_bonds} bonds"
