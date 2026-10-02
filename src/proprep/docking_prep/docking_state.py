"""
Every decision of one docking setup, in one object that round-trips through JSON.

It imports nothing from RDKit, so ProPrep starts even where RDKit is missing.

The Molecular Docking menus fill this in step by step; it is saved to the
workspace (key ``docking_state``) after each change so a setup survives leaving
the menu, and it is what a run is built from. Nothing here has a default that
changes a result without being shown: fields start empty or at the value the
menus display as the default, and ``missing`` says what a run still needs.

Residues are Meeko ids ("A:94", "A:105M"); metal keys are "residue|atom"
("A:105|FE") because JSON keys are strings.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Tuple


def metal_key(res_id: str, atom: str) -> str:
    return f"{res_id}|{atom}"


def split_metal_key(key: str) -> Tuple[str, str]:
    res_id, atom = key.split("|", 1)
    return res_id, atom


@dataclass
class ReceptorChoices:
    pdb_path: Optional[str] = None
    chains: List[str] = field(default_factory=list)
    keep_hetero: List[str] = field(default_factory=list)
    hetero_reviewed: bool = False        # menu 3 has been through; until then nothing is kept
    altlocs: Dict[str, str] = field(default_factory=dict)
    altloc_fills: Dict[str, Dict[str, str]] = field(default_factory=dict)   # "A:129" -> {"N": "A"}
    ph: Optional[float] = None
    protonation: Dict[str, str] = field(default_factory=dict)
    protonation_sources: Dict[str, str] = field(default_factory=dict)
    termini: Dict[str, str] = field(default_factory=dict)
    termini_sources: Dict[str, str] = field(default_factory=dict)       # evidence, or "chosen"
    metal_charges: Dict[str, int] = field(default_factory=dict)          # "A:105|FE" -> 2
    metal_types: Dict[str, str] = field(default_factory=dict)            # "A:505|MO" -> "W" or "omit"
    cofactor_edits: Dict[str, List[dict]] = field(default_factory=dict)  # resname -> ChemistryEdit dicts
    flexible: List[str] = field(default_factory=list)

    def build_arguments(self) -> dict:
        """Keyword arguments for receptor_prep.prepare_receptor."""
        from .chemistry_edits import ChemistryEdit
        return dict(
            protonation=dict(self.protonation),
            termini=dict(self.termini),
            metal_charges={split_metal_key(k): v for k, v in self.metal_charges.items()},
            metal_types={split_metal_key(k): v for k, v in self.metal_types.items()},
            altlocs=dict(self.altlocs),
            altloc_fills={k: dict(v) for k, v in self.altloc_fills.items()},
            cofactor_edits={name: [ChemistryEdit.from_dict(e) for e in edits]
                            for name, edits in self.cofactor_edits.items()},
            flexible=list(self.flexible),
        )


@dataclass
class LigandChoices:
    source: Dict[str, object] = field(default_factory=dict)     # LigandBuild.source
    ccd_edits: List[dict] = field(default_factory=list)          # edits to the CCD entry (structure route)
    edits: List[dict] = field(default_factory=list)              # edits to the built ligand
    metal_charges: Dict[str, int] = field(default_factory=dict)  # atom name -> charge
    rigid_bonds: List[List[str]] = field(default_factory=list)
    rotatable_amides: List[List[str]] = field(default_factory=list)
    rigid_macrocycles: Optional[bool] = None

    def edit_objects(self) -> list:
        from .chemistry_edits import ChemistryEdit
        return [ChemistryEdit.from_dict(e) for e in self.edits]


@dataclass
class DockingState:
    receptor: ReceptorChoices = field(default_factory=ReceptorChoices)
    ligand: LigandChoices = field(default_factory=LigandChoices)
    box: Dict[str, object] = field(default_factory=dict)         # SearchBox.to_dict()
    settings: Dict[str, object] = field(default_factory=dict)    # DockingSettings.to_dict()
    output_dir: Optional[str] = None
    last_run: Dict[str, object] = field(default_factory=dict)
    last_campaign: Dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Optional[dict]) -> "DockingState":
        data = dict(data or {})
        return cls(
            receptor=ReceptorChoices(**data.get("receptor", {})),
            ligand=LigandChoices(**data.get("ligand", {})),
            box=dict(data.get("box", {})),
            settings=dict(data.get("settings", {})),
            output_dir=data.get("output_dir"),
            last_run=dict(data.get("last_run", {})),
            last_campaign=dict(data.get("last_campaign", {})),
        )

    def ligand_residue(self) -> Optional[str]:
        """The residue id of a ligand taken from the receptor structure, else None."""
        s = self.ligand.source
        if s.get("kind") != "structure_residue":
            return None
        return f"{s['chain']}:{s['resseq']}{s.get('icode', '')}"

    def unreviewed_hetero(self) -> List[str]:
        """Non-water HETATM residues of the receptor file other than the ligand (menu 3 decides them)."""
        from .receptor_prep import read_atoms
        ligand = self.ligand_residue()
        seen = []
        for atom in read_atoms(self.receptor.pdb_path):
            if (atom.record == "HETATM" and atom.resname not in ("HOH", "WAT", "DOD")
                    and atom.res_id != ligand and atom.res_id not in seen):
                seen.append(atom.res_id)
        return seen

    def undecided_ends(self) -> List[str]:
        """Uncapped chain ends of the receptor without a terminus decision."""
        from .receptor_decisions import chain_ends
        return [e.res_id for e in chain_ends(self.receptor.pdb_path, self.receptor.chains)
                if not e.capped and e.res_id not in self.receptor.termini]

    def missing(self) -> List[str]:
        """What a run still needs, in menu order; empty when it can run."""
        needs = []
        r = self.receptor
        if not r.pdb_path:
            needs.append("a receptor structure (1)")
        elif not r.chains:
            needs.append("receptor chains (1)")
        if not self.ligand.source:
            needs.append("a ligand (2)")
        if r.pdb_path and not r.hetero_reviewed and self.unreviewed_hetero():
            needs.append("which HETATM residues stay (3)")
        if r.pdb_path and r.ph is None:
            needs.append("protonation (4)")
        elif r.pdb_path and r.chains and self.undecided_ends():
            needs.append("chain ends (4)")
        if self.ligand.source and self.ligand.rigid_macrocycles is None:
            needs.append("rotatable bonds reviewed (7)")
        if not self.box:
            needs.append("a search box (8)")
        if not self.settings:
            needs.append("docking settings (9)")
        return needs
