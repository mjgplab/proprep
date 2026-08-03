"""
Global Atom Type Registry for Multi-Site Metal Center Parameterization

This module provides coordination of atom type assignment across multiple metal sites
to ensure globally unique M/Y/Z prefixed types. This allows parameters from different
sites to coexist in the same tLEaP session without conflicts.

The Problem:
    When parameterizing Site A and Site B separately, both might generate:
    - Site A: M1, Y1, Y2, Y3...
    - Site B: M1, Y1, Y2, Y3... (CONFLICT!)

The Solution:
    GlobalAtomTypeRegistry analyzes ALL sites first, then assigns:
    - Site A: M1, Y1, Y2, Y3...
    - Site B: M2, Y4, Y5, Y6... (NO CONFLICT!)

Usage:
    # Create registry with all sites to be parameterized
    registry = GlobalAtomTypeRegistry(console)
    registry.register_sites(redox_sites)

    # Assign types across all sites (must be called before get_types_for_site)
    registry.assign_all_types()

    # Get type assignments for a specific site
    site_types = registry.get_types_for_site(site_index)

    # Pass site_types to Step 1 instead of letting it generate its own

Author: ProPrep Development Team
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any, Set
from enum import Enum

from rich.console import Console
from rich.table import Table
from rich.panel import Panel

# Import from structure_prep
from proprep.structure_prep.comprehensive_redox_detector import (
    RedoxSite, RedoxSiteAtom, RedoxCenter, METALS
)


@dataclass
class GlobalTypeAssignment:
    """
    Type assignment for an atom with global context.

    Extends the concept of AtomTypeAssignment to include site-level information.
    """
    # Atom identification (unique across all sites)
    site_index: int
    coords: Tuple[float, float, float]
    chain: str
    resname: str
    resid: int
    atom_name: str
    element: str

    # Type assignment
    original_type: str        # From force field library (to be looked up later)
    global_renamed_type: str  # The globally unique M/Y/Z type

    # Classification
    is_metal: bool
    is_ligand: bool           # Direct metal-coordinating atom

    # For debugging/reporting
    site_id: str = ""


class GlobalAtomTypeRegistry:
    """
    Coordinates atom type assignment across multiple metal sites.

    Ensures globally unique M/Y/Z prefixed atom types so that parameters
    from different sites can coexist in the same tLEaP session.

    Naming scheme (following MCPB conventions):
        M1-M9, MA-MZ: Metal center atoms (35 per prefix)
        Y1-Y9, YA-YZ: Direct metal ligand atoms (35 per prefix)

    For >35 metals, additional prefixes can be used (X, Z, etc.)
    For >35 ligands per metal, the counter continues across sites.
    """

    def __init__(self, console: Console = None, logger=None):
        """
        Initialize the global registry.

        Args:
            console: Rich console for output
            logger: Optional logger instance
        """
        self.console = console or Console()
        self.logger = logger or logging.getLogger(__name__)

        # Registered sites
        self.sites: List[RedoxSite] = []
        self.site_indices: Dict[str, int] = {}  # site_id → index

        # Global counters (continue across sites)
        self.metal_counter = 0      # Next available metal index
        self.ligand_counter = 0     # Next available ligand index

        # Global assignments: (site_index, coords) → GlobalTypeAssignment
        self.global_assignments: Dict[Tuple[int, Tuple[float, float, float]], GlobalTypeAssignment] = {}

        # Per-site tracking for efficient lookup
        self.site_metals: Dict[int, List[Tuple[float, float, float]]] = {}
        self.site_ligands: Dict[int, List[Tuple[float, float, float]]] = {}

        # Assignment complete flag
        self._assigned = False

    def register_sites(self, sites: List[RedoxSite]) -> None:
        """
        Register multiple RedoxSites for global type assignment.

        Args:
            sites: List of RedoxSite objects to parameterize together
        """
        self.sites = sites
        self.site_indices = {site.site_id: idx for idx, site in enumerate(sites)}

        # Initialize per-site tracking
        for idx in range(len(sites)):
            self.site_metals[idx] = []
            self.site_ligands[idx] = []

        self.console.print(f"[cyan]Registered {len(sites)} site(s) for global atom typing[/cyan]")
        for idx, site in enumerate(sites):
            n_centers = len(site.centers)
            n_atoms = len(site.atoms)
            self.console.print(f"  [grey50]Site {idx}: {site.site_id} - {n_centers} center(s), {n_atoms} atoms[/grey50]")

    def register_single_site(self, site: RedoxSite) -> None:
        """
        Convenience method to register a single site.

        Args:
            site: Single RedoxSite object
        """
        self.register_sites([site])

    def assign_all_types(self) -> Dict[Tuple[int, Tuple[float, float, float]], GlobalTypeAssignment]:
        """
        Analyze all registered sites and assign globally unique types.

        This must be called after register_sites() and before get_types_for_site().

        Returns:
            Dictionary of all global assignments
        """
        if not self.sites:
            raise ValueError("No sites registered. Call register_sites() first.")

        self.console.print("\n[bold cyan]Global Atom Type Assignment[/bold cyan]")
        self.console.print(f"[grey50]Analyzing {len(self.sites)} site(s) for global type coordination...[/grey50]")

        # Reset counters and assignments
        self.metal_counter = 0
        self.ligand_counter = 0
        self.global_assignments.clear()

        # Process each site
        for site_idx, site in enumerate(self.sites):
            self._assign_types_for_site(site_idx, site)

        self._assigned = True

        # Display summary
        self._display_assignment_summary()

        return self.global_assignments

    def _assign_types_for_site(self, site_idx: int, site: RedoxSite) -> None:
        """
        Assign types for all atoms in a single site using global counters.

        Args:
            site_idx: Index of the site in self.sites
            site: RedoxSite object
        """
        self.logger.info(f"Assigning types for site {site_idx}: {site.site_id}")

        # Step 1: Identify metals in this site
        # Metals come from two center types:
        # - METAL_ION: the center itself is the metal (e.g., isolated Zn)
        # - ORGANOMETALLIC_COFACTOR: contains embedded metal atoms (e.g., Fe in heme)
        # We check both centers and all site atoms to catch both cases.
        site_metal_coords = set()

        # Check centers for METAL_ION type (center coords ARE the metal)
        for center in site.centers:
            if hasattr(center, 'center_type') and center.center_type.value == 'metal_ion':
                site_metal_coords.add(center.coords)

        # Check all atoms for embedded metals (ORGANOMETALLIC_COFACTOR case)
        # This catches Fe in heme, metals in Fe-S clusters, etc.
        for atom in site.atoms:
            element = atom.element or ""
            if element.upper() in METALS:
                site_metal_coords.add(atom.coords)

        # Step 2: Identify ligands (atoms coordinated to metals)
        site_ligand_coords = set()
        for bond in site.bonds:
            if bond.chemical_type == 'coordinate':
                # One atom is metal, the other is ligand
                if bond.atom1_coords in site_metal_coords:
                    site_ligand_coords.add(bond.atom2_coords)
                elif bond.atom2_coords in site_metal_coords:
                    site_ligand_coords.add(bond.atom1_coords)

        # Store for per-site lookup
        self.site_metals[site_idx] = list(site_metal_coords)
        self.site_ligands[site_idx] = list(site_ligand_coords)

        # Step 3: Assign global types to metals
        for coords in site_metal_coords:
            atom = self._find_atom_by_coords(site, coords)
            if atom:
                global_type = self._generate_systematic_name('M', self.metal_counter)
                self.metal_counter += 1

                self.global_assignments[(site_idx, coords)] = GlobalTypeAssignment(
                    site_index=site_idx,
                    coords=coords,
                    chain=atom.chain,
                    resname=atom.resname,
                    resid=atom.resid,
                    atom_name=atom.atom_name,
                    element=atom.element,
                    original_type="",  # Will be filled by MCPBAtomTyper later
                    global_renamed_type=global_type,
                    is_metal=True,
                    is_ligand=False,
                    site_id=site.site_id
                )

        # Step 4: Assign global types to ligands
        for coords in site_ligand_coords:
            atom = self._find_atom_by_coords(site, coords)
            if atom:
                global_type = self._generate_systematic_name('Y', self.ligand_counter)
                self.ligand_counter += 1

                self.global_assignments[(site_idx, coords)] = GlobalTypeAssignment(
                    site_index=site_idx,
                    coords=coords,
                    chain=atom.chain,
                    resname=atom.resname,
                    resid=atom.resid,
                    atom_name=atom.atom_name,
                    element=atom.element,
                    original_type="",
                    global_renamed_type=global_type,
                    is_metal=False,
                    is_ligand=True,
                    site_id=site.site_id
                )

        # Step 5: Record non-metal/non-ligand atoms (they keep original types)
        for atom in site.atoms:
            if atom.coords not in site_metal_coords and atom.coords not in site_ligand_coords:
                self.global_assignments[(site_idx, atom.coords)] = GlobalTypeAssignment(
                    site_index=site_idx,
                    coords=atom.coords,
                    chain=atom.chain,
                    resname=atom.resname,
                    resid=atom.resid,
                    atom_name=atom.atom_name,
                    element=atom.element,
                    original_type="",
                    global_renamed_type="",  # Empty = use original type from FF
                    is_metal=False,
                    is_ligand=False,
                    site_id=site.site_id
                )

        self.logger.info(f"  Site {site_idx}: {len(site_metal_coords)} metals, {len(site_ligand_coords)} ligands")

    def _find_atom_by_coords(self, site: RedoxSite, coords: Tuple[float, float, float]) -> Optional[RedoxSiteAtom]:
        """Find an atom in the site by its coordinates."""
        for atom in site.atoms:
            if atom.coords == coords:
                return atom
        # Also check centers (metals may only be in centers, not atoms)
        for center in site.centers:
            if center.coords == coords:
                # Create a pseudo-atom from center
                return RedoxSiteAtom(
                    chain=center.chain,
                    resname=center.resname,
                    resid=center.resid,
                    atom_name=center.atom_name or center.element,
                    coords=center.coords,
                    element=center.element or ""
                )
        return None

    def _generate_systematic_name(self, prefix: str, index: int) -> str:
        """
        Generate systematic 2-character atom type name.

        Follows MCPB conventions:
        - Digits first: M1-M9 (9 names)
        - Then letters: MA-MZ (26 names)
        Total: 35 names per prefix

        For additional capacity, switches to secondary prefixes:
        - M: 0-34
        - X: 35-69
        - Additional prefixes as needed

        Args:
            prefix: Primary prefix ('M' for metals, 'Y' for ligands)
            index: 0-based global index

        Returns:
            2-character atom type name
        """
        # Determine which prefix set to use
        names_per_prefix = 35
        prefix_index = index // names_per_prefix
        local_index = index % names_per_prefix

        # Prefix progression for metals: M, X, W, V, U...
        # Prefix progression for ligands: Y, Z, Q, P, O...
        if prefix == 'M':
            prefixes = ['M', 'X', 'W', 'V', 'U', 'T', 'R']
        else:  # 'Y' for ligands
            prefixes = ['Y', 'Z', 'Q', 'P', 'O', 'N', 'K']

        if prefix_index >= len(prefixes):
            raise ValueError(
                f"Exceeded atom type capacity: index {index} requires prefix {prefix_index}, "
                f"but only {len(prefixes)} prefixes available. Consider 3-character types."
            )

        actual_prefix = prefixes[prefix_index]

        # Generate the suffix (1-9, then A-Z)
        if local_index < 9:
            return f"{actual_prefix}{local_index + 1}"
        else:
            letter_index = local_index - 9
            letter = chr(ord('A') + letter_index)
            return f"{actual_prefix}{letter}"

    def get_types_for_site(self, site_index: int) -> Dict[Tuple[float, float, float], GlobalTypeAssignment]:
        """
        Get the global type assignments for a specific site.

        Args:
            site_index: Index of the site (0-based)

        Returns:
            Dictionary mapping coords → GlobalTypeAssignment for that site
        """
        if not self._assigned:
            raise RuntimeError("Types not yet assigned. Call assign_all_types() first.")

        if site_index < 0 or site_index >= len(self.sites):
            raise ValueError(f"Invalid site index: {site_index}. Valid range: 0-{len(self.sites)-1}")

        return {
            coords: assignment
            for (idx, coords), assignment in self.global_assignments.items()
            if idx == site_index
        }

    def get_site_by_id(self, site_id: str) -> Optional[int]:
        """
        Get site index by site_id.

        Args:
            site_id: The site_id string

        Returns:
            Site index or None if not found
        """
        return self.site_indices.get(site_id)

    def get_metal_types_for_site(self, site_index: int) -> List[str]:
        """Get list of metal atom types assigned to a specific site."""
        if not self._assigned:
            raise RuntimeError("Types not yet assigned. Call assign_all_types() first.")

        return [
            assignment.global_renamed_type
            for (idx, coords), assignment in self.global_assignments.items()
            if idx == site_index and assignment.is_metal
        ]

    def get_ligand_types_for_site(self, site_index: int) -> List[str]:
        """Get list of ligand atom types assigned to a specific site."""
        if not self._assigned:
            raise RuntimeError("Types not yet assigned. Call assign_all_types() first.")

        return [
            assignment.global_renamed_type
            for (idx, coords), assignment in self.global_assignments.items()
            if idx == site_index and assignment.is_ligand
        ]

    def _display_assignment_summary(self) -> None:
        """Display a summary of global type assignments."""
        self.console.print("\n[bold green]Global Type Assignment Complete[/bold green]")

        # Create summary table
        table = Table(title="Global Atom Type Summary")
        table.add_column("Site", style="cyan")
        table.add_column("Site ID", style="magenta")
        table.add_column("Metals", style="yellow", justify="right")
        table.add_column("Metal Types", style="yellow")
        table.add_column("Ligands", style="green", justify="right")
        table.add_column("Ligand Types", style="green")

        for site_idx, site in enumerate(self.sites):
            metal_types = self.get_metal_types_for_site(site_idx)
            ligand_types = self.get_ligand_types_for_site(site_idx)

            metal_range = self._format_type_range(metal_types) if metal_types else "-"
            ligand_range = self._format_type_range(ligand_types) if ligand_types else "-"

            table.add_row(
                str(site_idx),
                site.site_id[:20],  # Truncate long IDs
                str(len(metal_types)),
                metal_range,
                str(len(ligand_types)),
                ligand_range
            )

        self.console.print(table)

        # Total summary
        total_metals = self.metal_counter
        total_ligands = self.ligand_counter
        self.console.print(f"\n[grey50]Total unique metal types: {total_metals}[/grey50]")
        self.console.print(f"[grey50]Total unique ligand types: {total_ligands}[/grey50]")
        self.console.print(f"[green]✓ No atom type conflicts between sites[/green]")

    def _format_type_range(self, types: List[str]) -> str:
        """Format a list of types as a range string (e.g., 'M1-M3' or 'M1,M2,M5')."""
        if not types:
            return "-"
        if len(types) == 1:
            return types[0]
        if len(types) <= 3:
            return ",".join(types)
        return f"{types[0]}-{types[-1]}"

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialize registry state to dictionary for saving/resuming.

        Returns:
            Dictionary representation of registry state
        """
        return {
            "metal_counter": self.metal_counter,
            "ligand_counter": self.ligand_counter,
            "site_ids": [site.site_id for site in self.sites],
            "assignments": {
                f"{idx},{coords[0]},{coords[1]},{coords[2]}": {
                    "global_renamed_type": a.global_renamed_type,
                    "is_metal": a.is_metal,
                    "is_ligand": a.is_ligand,
                    "element": a.element,
                    "resname": a.resname,
                    "atom_name": a.atom_name
                }
                for (idx, coords), a in self.global_assignments.items()
            }
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any], sites: List[RedoxSite],
                  console: Console = None) -> 'GlobalAtomTypeRegistry':
        """
        Restore registry from saved dictionary.

        Args:
            data: Dictionary from to_dict()
            sites: List of RedoxSite objects (must match saved site_ids)
            console: Rich console

        Returns:
            Restored GlobalAtomTypeRegistry
        """
        registry = cls(console=console)
        registry.register_sites(sites)

        # Restore counters
        registry.metal_counter = data["metal_counter"]
        registry.ligand_counter = data["ligand_counter"]

        # Restore assignments
        for key_str, assignment_data in data["assignments"].items():
            parts = key_str.split(",")
            idx = int(parts[0])
            coords = (float(parts[1]), float(parts[2]), float(parts[3]))

            # Find the atom in the site
            site = sites[idx]
            atom = registry._find_atom_by_coords(site, coords)

            if atom:
                registry.global_assignments[(idx, coords)] = GlobalTypeAssignment(
                    site_index=idx,
                    coords=coords,
                    chain=atom.chain,
                    resname=assignment_data["resname"],
                    resid=atom.resid,
                    atom_name=assignment_data["atom_name"],
                    element=assignment_data["element"],
                    original_type="",
                    global_renamed_type=assignment_data["global_renamed_type"],
                    is_metal=assignment_data["is_metal"],
                    is_ligand=assignment_data["is_ligand"],
                    site_id=site.site_id
                )

        registry._assigned = True
        return registry
