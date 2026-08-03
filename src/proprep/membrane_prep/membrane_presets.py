"""
Membrane Composition Presets

Common membrane compositions for quick setup. Each preset defines lipid
composition and ratio, with a description of the biological context.
"""

from dataclasses import dataclass
from typing import List


@dataclass
class MembranePreset:
    """A predefined membrane composition."""
    name: str
    description: str
    lipids: str
    ratio: str
    notes: str = ""


PRESETS: List[MembranePreset] = [
    MembranePreset(
        name="Simple DOPC bilayer",
        description="Symmetric, single-component. Good starting point.",
        lipids="DOPC",
        ratio="1",
    ),
    MembranePreset(
        name="Simple POPC bilayer",
        description="Symmetric, single-component. Widely used model membrane.",
        lipids="POPC",
        ratio="1",
    ),
    MembranePreset(
        name="Mammalian plasma membrane (simple)",
        description="DOPC + cholesterol. Captures basic fluidity modulation.",
        lipids="DOPC:CHL1",
        ratio="2:1",
    ),
    MembranePreset(
        name="Mammalian with sphingomyelin",
        description="DOPC + cholesterol + sphingomyelin. Lipid raft model.",
        lipids="DOPC:CHL1:PSM",
        ratio="2:1:1",
        notes="Requires Lipid21",
    ),
    MembranePreset(
        name="Bacterial inner membrane (E. coli)",
        description="POPE + POPG. Typical Gram-negative inner membrane.",
        lipids="POPE:POPG",
        ratio="3:1",
    ),
    MembranePreset(
        name="Bacterial outer membrane (simple)",
        description="DPPE + DPPG. Simplified outer leaflet composition.",
        lipids="DPPE:DPPG",
        ratio="1:1",
    ),
    MembranePreset(
        name="Mitochondrial inner membrane",
        description="POPC + POPE + cardiolipin. Site of oxidative phosphorylation.",
        lipids="POPC:POPE:TLCL2",
        ratio="4:3:1",
    ),
    MembranePreset(
        name="Thylakoid membrane",
        description="MGDG + DGDG + SQDG + POPG. Photosynthetic membrane.",
        lipids="MGDG:DGDG:SQDG:POPG",
        ratio="5:3:1:1",
    ),
    MembranePreset(
        name="ER membrane (simple)",
        description="DOPC + DOPE + DOPS. Loosely packed, biosynthetic hub.",
        lipids="DOPC:DOPE:DOPS",
        ratio="5:3:1",
    ),
    MembranePreset(
        name="Yeast plasma membrane",
        description="POPC + POPE + ergosterol. Fungal cell membrane.",
        lipids="POPC:POPE:ERG",
        ratio="2:1:1",
    ),
]
