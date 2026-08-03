"""
Atom Type Counter for MCPB Systematic Naming

This module provides the AtomTypeCounter class which assigns globally unique
M*/Y* atom type names following MCPB conventions.

The key insight: single-site and multi-site parameterization are THE SAME
workflow. The only difference is where the counter starts:

    Single site:  M1, Y1, Y2, Y3, Y4, Y5 (counter starts at 0)
    Multi-site:
        Site 0:   M1, Y1, Y2, Y3, Y4, Y5 (counter starts at 0)
        Site 1:   M2, Y6, Y7, Y8, Y9     (counter continues)
        Site 2:   M3, Y10, Y11, Y12      (counter continues)

This replaces GlobalAtomTypeRegistry which was unnecessarily complex.

Naming Scheme (MCPB conventions):
    Metals:  M1-M9, MA-MZ (35 names), then X1-X9, XA-XZ, W1-W9, WA-WZ...
    Ligands: Y1-Y9, YA-YZ (35 names), then Z1-Z9, ZA-ZZ, Q1-Q9, QA-QZ...

Author: ProPrep Development Team
"""

from typing import Tuple, List, Dict, Optional
from dataclasses import dataclass, field


@dataclass
class AtomTypeCounter:
    """
    Assigns globally unique M*/Y* atom type names.

    Usage for single site:
        counter = AtomTypeCounter()
        metal_type = counter.next_metal_type()  # Returns "M1"
        lig1 = counter.next_ligand_type()       # Returns "Y1"
        lig2 = counter.next_ligand_type()       # Returns "Y2"

    Usage for multi-site:
        counter = AtomTypeCounter()

        # Site 0
        m1 = counter.next_metal_type()          # "M1"
        y1 = counter.next_ligand_type()         # "Y1"
        y2 = counter.next_ligand_type()         # "Y2"

        # Site 1 (counter continues)
        m2 = counter.next_metal_type()          # "M2"
        y3 = counter.next_ligand_type()         # "Y3"

    The counter can also be saved/restored for workflow resume capability.
    """

    # Current counter values
    metal_counter: int = 0
    ligand_counter: int = 0

    # Tracking which types have been assigned
    assigned_metal_types: List[str] = field(default_factory=list)
    assigned_ligand_types: List[str] = field(default_factory=list)

    # Prefix sequences following MCPB conventions
    # Each prefix provides 35 unique names (9 digits + 26 letters)
    METAL_PREFIXES: Tuple[str, ...] = ('M', 'X', 'W', 'V', 'U', 'T', 'R', 'S')
    LIGAND_PREFIXES: Tuple[str, ...] = ('Y', 'Z', 'Q', 'P', 'O', 'N', 'K', 'J')

    NAMES_PER_PREFIX: int = 35  # 9 digits (1-9) + 26 letters (A-Z)

    def next_metal_type(self) -> str:
        """
        Get next available metal type name.

        Returns:
            2-character atom type (e.g., "M1", "M2", ..., "MA", ..., "X1")

        Raises:
            ValueError: If all metal type names exhausted (unlikely: 280 names)
        """
        name = self._generate_name(self.METAL_PREFIXES, self.metal_counter)
        self.metal_counter += 1
        self.assigned_metal_types.append(name)
        return name

    def next_ligand_type(self) -> str:
        """
        Get next available ligand type name.

        Returns:
            2-character atom type (e.g., "Y1", "Y2", ..., "YA", ..., "Z1")

        Raises:
            ValueError: If all ligand type names exhausted (unlikely: 280 names)
        """
        name = self._generate_name(self.LIGAND_PREFIXES, self.ligand_counter)
        self.ligand_counter += 1
        self.assigned_ligand_types.append(name)
        return name

    def _generate_name(self, prefixes: Tuple[str, ...], index: int) -> str:
        """
        Generate a 2-character atom type name from index.

        Naming scheme:
            Index 0-8:   X1, X2, ..., X9   (digits 1-9)
            Index 9-34:  XA, XB, ..., XZ   (letters A-Z)
            Index 35-69: Next prefix with same pattern

        Args:
            prefixes: Tuple of prefix characters to use
            index: 0-based index

        Returns:
            2-character atom type name

        Raises:
            ValueError: If index exceeds capacity
        """
        # Determine which prefix to use
        prefix_index = index // self.NAMES_PER_PREFIX
        local_index = index % self.NAMES_PER_PREFIX

        if prefix_index >= len(prefixes):
            raise ValueError(
                f"Atom type capacity exceeded: index {index} requires "
                f"prefix {prefix_index + 1}, but only {len(prefixes)} prefixes "
                f"available ({prefixes}). Maximum capacity: "
                f"{len(prefixes) * self.NAMES_PER_PREFIX} types."
            )

        prefix = prefixes[prefix_index]

        # Generate suffix
        if local_index < 9:
            # Digits 1-9
            suffix = str(local_index + 1)
        else:
            # Letters A-Z
            letter_index = local_index - 9
            suffix = chr(ord('A') + letter_index)

        return f"{prefix}{suffix}"

    def peek_next_metal_type(self) -> str:
        """Preview next metal type without incrementing counter."""
        return self._generate_name(self.METAL_PREFIXES, self.metal_counter)

    def peek_next_ligand_type(self) -> str:
        """Preview next ligand type without incrementing counter."""
        return self._generate_name(self.LIGAND_PREFIXES, self.ligand_counter)

    def get_statistics(self) -> Dict[str, int]:
        """
        Get current counter statistics.

        Returns:
            Dict with counts and capacity info
        """
        metal_capacity = len(self.METAL_PREFIXES) * self.NAMES_PER_PREFIX
        ligand_capacity = len(self.LIGAND_PREFIXES) * self.NAMES_PER_PREFIX

        return {
            'metals_assigned': self.metal_counter,
            'ligands_assigned': self.ligand_counter,
            'metal_capacity': metal_capacity,
            'ligand_capacity': ligand_capacity,
            'metal_remaining': metal_capacity - self.metal_counter,
            'ligand_remaining': ligand_capacity - self.ligand_counter,
        }

    def reset(self) -> None:
        """Reset counters to initial state."""
        self.metal_counter = 0
        self.ligand_counter = 0
        self.assigned_metal_types.clear()
        self.assigned_ligand_types.clear()

    def to_dict(self) -> Dict:
        """
        Serialize state for saving/resuming workflow.

        Returns:
            Dict suitable for JSON serialization
        """
        return {
            'metal_counter': self.metal_counter,
            'ligand_counter': self.ligand_counter,
            'assigned_metal_types': list(self.assigned_metal_types),
            'assigned_ligand_types': list(self.assigned_ligand_types),
        }

    @classmethod
    def from_dict(cls, data: Dict) -> 'AtomTypeCounter':
        """
        Restore state from saved dictionary.

        Args:
            data: Dict from to_dict()

        Returns:
            AtomTypeCounter with restored state
        """
        counter = cls()
        counter.metal_counter = data.get('metal_counter', 0)
        counter.ligand_counter = data.get('ligand_counter', 0)
        counter.assigned_metal_types = list(data.get('assigned_metal_types', []))
        counter.assigned_ligand_types = list(data.get('assigned_ligand_types', []))
        return counter

    def get_type_range_for_site(
        self,
        n_metals: int,
        n_ligands: int
    ) -> Tuple[List[str], List[str]]:
        """
        Get the type names that WOULD be assigned for a site with
        the given number of metals and ligands.

        This is useful for displaying what types a site will receive
        before actually assigning them.

        Args:
            n_metals: Number of metals in the site
            n_ligands: Number of ligands in the site

        Returns:
            Tuple of (metal_types, ligand_types) lists
        """
        metal_types = []
        ligand_types = []

        # Preview metal types
        for i in range(n_metals):
            metal_types.append(
                self._generate_name(self.METAL_PREFIXES, self.metal_counter + i)
            )

        # Preview ligand types
        for i in range(n_ligands):
            ligand_types.append(
                self._generate_name(self.LIGAND_PREFIXES, self.ligand_counter + i)
            )

        return metal_types, ligand_types


# =============================================================================
# Convenience Functions
# =============================================================================

def is_metal_type(atom_type: str) -> bool:
    """
    Check if an atom type is a metal type (M*, X*, W*, etc.).

    Args:
        atom_type: 2-character atom type

    Returns:
        True if this is a metal type
    """
    if len(atom_type) != 2:
        return False
    return atom_type[0] in AtomTypeCounter.METAL_PREFIXES


def is_ligand_type(atom_type: str) -> bool:
    """
    Check if an atom type is a ligand type (Y*, Z*, Q*, etc.).

    Args:
        atom_type: 2-character atom type

    Returns:
        True if this is a ligand type
    """
    if len(atom_type) != 2:
        return False
    return atom_type[0] in AtomTypeCounter.LIGAND_PREFIXES


def is_mcpb_renamed_type(atom_type: str) -> bool:
    """
    Check if an atom type is an MCPB-renamed type (metal or ligand).

    Args:
        atom_type: Atom type string

    Returns:
        True if this is an M*/Y*/X*/Z*/etc. renamed type
    """
    return is_metal_type(atom_type) or is_ligand_type(atom_type)


# =============================================================================
# Testing
# =============================================================================

def _test_counter():
    """Basic test of AtomTypeCounter functionality."""
    counter = AtomTypeCounter()

    # Test first few metal types
    assert counter.next_metal_type() == "M1"
    assert counter.next_metal_type() == "M2"

    # Test first few ligand types
    assert counter.next_ligand_type() == "Y1"
    assert counter.next_ligand_type() == "Y2"
    assert counter.next_ligand_type() == "Y3"

    # Test statistics
    stats = counter.get_statistics()
    assert stats['metals_assigned'] == 2
    assert stats['ligands_assigned'] == 3

    # Test serialization
    data = counter.to_dict()
    restored = AtomTypeCounter.from_dict(data)
    assert restored.metal_counter == 2
    assert restored.ligand_counter == 3

    # Test continuing after restore
    assert restored.next_metal_type() == "M3"
    assert restored.next_ligand_type() == "Y4"

    # Test digit to letter transition (index 9 = MA)
    counter2 = AtomTypeCounter()
    for _ in range(9):
        counter2.next_metal_type()
    assert counter2.next_metal_type() == "MA"  # 10th metal type

    # Test prefix transition (index 35 = X1)
    counter3 = AtomTypeCounter()
    for _ in range(35):
        counter3.next_metal_type()
    assert counter3.next_metal_type() == "X1"  # 36th metal type

    print("All tests passed!")


if __name__ == "__main__":
    _test_counter()
