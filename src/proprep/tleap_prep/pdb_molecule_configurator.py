"""
PDB Molecule Configurator

Allows users to configure how PDB chains/segments should be grouped into molecules
for AMBER topology generation. This ensures ATOMS_PER_MOLECULE section is correct
from the start, avoiding the need for ParmEd post-processing.

The key insight: AMBER treats bonded components as single "molecules". For systems
with intermolecular bonds (e.g., protein-heme coordination), each bonded component
must have its atoms grouped contiguously in the PDB file.

Example:
    A protein trimer with coordinated hemes would have 3 molecules:
    - Molecule 1: Protein chains 1-5 + their hemes
    - Molecule 2: Protein chains 6-11 + their hemes
    - Molecule 3: Protein chains 12-16 + their hemes

    The PDB must be ordered: [chains 1-5, hemes for 1-5, chains 6-11, hemes for 6-11, ...]
"""

import json
import logging
import re
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
)

logger = logging.getLogger(__name__)


@dataclass
class ChainSegment:
    """Represents a chain/segment in the PDB file."""
    chain_id: str
    start_residue: int
    end_residue: int
    residue_count: int
    is_protein: bool  # Based on backbone atom presence
    hetero_residues: List[str] = field(default_factory=list)
    protein_residues: List[int] = field(default_factory=list)  # Residue numbers with backbone atoms
    cofactor_residues: Dict[int, str] = field(default_factory=dict)  # res_num -> res_name for non-protein

    def __str__(self):
        seg_type = "Protein" if self.is_protein else "Hetero"
        return f"{seg_type} chain {self.chain_id}: residues {self.start_residue}-{self.end_residue} ({self.residue_count} residues)"

    def get_protein_range_str(self) -> str:
        """Get formatted string of protein residue ranges."""
        if not self.protein_residues:
            return "-"
        return f"{min(self.protein_residues)}-{max(self.protein_residues)}"

    def get_cofactor_str(self) -> str:
        """Get formatted string of cofactor residues."""
        if not self.cofactor_residues:
            return "-"
        # Group by residue name
        by_name = {}
        for res_num, res_name in self.cofactor_residues.items():
            if res_name not in by_name:
                by_name[res_name] = []
            by_name[res_name].append(res_num)

        # Format as "RES(num1,num2,...)"
        parts = []
        for res_name in sorted(by_name.keys()):
            res_nums = sorted(by_name[res_name])
            nums_str = ",".join(str(n) for n in res_nums[:3])
            if len(res_nums) > 3:
                nums_str += f",+{len(res_nums)-3}"
            parts.append(f"{res_name}({nums_str})")

        return ", ".join(parts)


@dataclass
class MoleculeGroup:
    """Represents a molecule composed of multiple chains/segments."""
    name: str
    chain_ids: List[str]
    description: str = ""

    def __str__(self):
        return f"{self.name}: chains {', '.join(self.chain_ids)}"


@dataclass
class MoleculeConfiguration:
    """Complete molecule grouping configuration."""
    pdb_file: str
    molecules: List[MoleculeGroup]
    chain_segments: List[ChainSegment]

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "pdb_file": self.pdb_file,
            "molecules": [asdict(m) for m in self.molecules],
            "chain_segments": [asdict(c) for c in self.chain_segments]
        }

    @classmethod
    def from_dict(cls, data: Dict) -> 'MoleculeConfiguration':
        """Load from dictionary."""
        return cls(
            pdb_file=data["pdb_file"],
            molecules=[MoleculeGroup(**m) for m in data["molecules"]],
            chain_segments=[ChainSegment(**c) for c in data["chain_segments"]]
        )


class PDBReorderMapper:
    """
    Tracks residue number mappings during PDB reordering.

    When reordering PDB for AMBER molecule grouping, residues are renumbered
    sequentially. This mapper tracks the transformations so RedoxSite objects
    can be synchronized.
    """

    def __init__(self, console: Console = None):
        self.console = console or Console()
        # Per-chain mapping: {chain_id: {old_resnum: new_resnum}}
        self.residue_mappings: Dict[str, Dict[int, int]] = {}
        # Chain remapping (if chains are renamed): {old_chain: new_chain}
        self.chain_mapping: Dict[str, str] = {}

    def add_chain_mapping(self, old_chain: str, new_chain: str):
        """Record chain ID mapping if chains are renamed."""
        self.chain_mapping[old_chain] = new_chain

    def add_residue_mapping(self, chain_id: str, old_resnum: int, new_resnum: int):
        """Record residue number mapping for a specific chain."""
        if chain_id not in self.residue_mappings:
            self.residue_mappings[chain_id] = {}
        self.residue_mappings[chain_id][old_resnum] = new_resnum

    def get_new_chain(self, old_chain: str) -> str:
        """Get new chain ID (or return same if no mapping)."""
        return self.chain_mapping.get(old_chain, old_chain)

    def get_new_resnum(self, chain_id: str, old_resnum: int) -> int:
        """Get new residue number (or return same if no mapping)."""
        # Use the chain that exists in mappings (might be mapped chain)
        actual_chain = self.get_new_chain(chain_id)
        return self.residue_mappings.get(actual_chain, {}).get(old_resnum, old_resnum)

    def log_summary(self):
        """Log mapping summary for debugging."""
        if self.chain_mapping:
            self.console.print(f"[cyan]Chain mappings: {self.chain_mapping}[/cyan]")

        for chain_id, mapping in self.residue_mappings.items():
            if mapping:
                old_nums = sorted(mapping.keys())
                new_nums = [mapping[k] for k in old_nums]
                self.console.print(
                    f"[cyan]Chain {chain_id}: {len(mapping)} residues "
                    f"({min(old_nums)}-{max(old_nums)} → {min(new_nums)}-{max(new_nums)})[/cyan]"
                )


class PDBMoleculeConfigurator:
    """Configure molecule grouping for AMBER topology generation."""

    BACKBONE_ATOMS = {"N", "CA", "C", "O"}

    def __init__(self, console: Console = None, processor=None):
        self.console = console or Console()
        self.processor = processor

    def _prompt(self, prompt, default=None, choices=None, description=None, options_map=None):
        """Prompt with optional context recording."""
        return prompt_with_context(
            self.processor,
            f"{prompt} ",
            choices=choices,
            default=default,
            module="tLEaP - Molecule Configurator",
            description=description,
            options_map=options_map,
        )

    def _confirm(self, prompt, default=True, description=None):
        """Confirm with optional context recording."""
        return confirm_with_context(
            self.processor,
            f"{prompt} ",
            default=default,
            module="tLEaP - Molecule Configurator",
            description=description,
        )
        return response

    def synchronize_protonation_results(self, protonation_results: Dict, mapper: PDBReorderMapper) -> Dict[str, int]:
        """
        Synchronize protonation analysis results after PDB reordering and renumbering.

        Updates residue numbers in protonation_results dict to match reordered structure.

        Args:
            protonation_results: Dict mapping residue keys to protonation data
            mapper: PDBReorderMapper with residue number mappings

        Returns:
            Summary dict with count of residues updated
        """
        if not protonation_results:
            return {'residues_updated': 0}

        self.console.print("\n[bold cyan]Synchronizing protonation results with reordered structure...[/bold cyan]")
        self.console.print(f"[grey50]Found {len(protonation_results)} entries in protonation results[/grey50]")

        residues_updated = 0
        updated_results = {}
        sample_updates = []

        for res_key, res_data in protonation_results.items():
            if isinstance(res_data, dict):
                chain = res_data.get('chain', 'A')
                old_resnum = res_data.get('number')

                if old_resnum is not None:
                    # Get new residue number from mapper
                    new_resnum = mapper.get_new_resnum(chain, old_resnum)

                    if new_resnum != old_resnum:
                        # Create updated res_data with new residue number
                        updated_data = res_data.copy()
                        updated_data['number'] = new_resnum

                        # Create new key with updated residue number
                        # Key format is typically "chain:resnum" or similar
                        if ':' in res_key:
                            parts = res_key.split(':')
                            if len(parts) == 2:
                                new_key = f"{chain}:{new_resnum}"
                            else:
                                new_key = res_key  # Keep original key if format is unexpected
                        else:
                            new_key = res_key

                        updated_results[new_key] = updated_data
                        residues_updated += 1

                        # Save first 5 examples for debugging
                        if len(sample_updates) < 5:
                            resname = res_data.get('type', '???')
                            sample_updates.append(f"{chain}:{old_resnum} {resname} → {chain}:{new_resnum}")
                    else:
                        # No change needed
                        updated_results[res_key] = res_data
                else:
                    # No residue number - keep as is
                    updated_results[res_key] = res_data
            else:
                # Not a dict - keep as is
                self.console.print(f"[yellow]Warning: Non-dict entry in protonation results: {res_key}[/yellow]")
                updated_results[res_key] = res_data

        # Update the original dict in-place
        protonation_results.clear()
        protonation_results.update(updated_results)

        self.console.print(f"[green]✓ Updated {residues_updated} residue numbers in protonation results[/green]")

        if sample_updates:
            self.console.print("[grey50]Sample updates:[/grey50]")
            for update in sample_updates:
                self.console.print(f"[grey50]  {update}[/grey50]")

        return {'residues_updated': residues_updated}

    def synchronize_redox_sites(self, redox_sites: List, mapper: PDBReorderMapper) -> Dict[str, int]:
        """
        Synchronize RedoxSite objects after PDB reordering and renumbering.

        Updates ALL RedoxSite attributes to reflect new residue numbers:
        - centers: chain, resid
        - atoms: chain, resid
        - bonds: atom1_residue_info, atom2_residue_info
        - coord_to_pdb: all PDB metadata
        - residue_groups: keys (chain, resid, resname)

        Note: Coordinates are IMMUTABLE and never change.

        Args:
            redox_sites: List of RedoxSite objects to update
            mapper: PDBReorderMapper with residue number mappings

        Returns:
            Summary dict with counts of sites/atoms/bonds updated
        """
        if not redox_sites:
            return {'sites_updated': 0, 'atoms_updated': 0, 'bonds_updated': 0, 'centers_updated': 0}

        self.console.print("\n[bold cyan]Synchronizing RedoxSites with reordered structure...[/bold cyan]")

        sites_updated = 0
        total_atoms_updated = 0
        total_bonds_updated = 0
        total_centers_updated = 0

        for site in redox_sites:
            site_type = getattr(site, 'site_type', 'unknown')
            self.console.print(f"\n[bold cyan]→ {site.site_id}[/bold cyan] [grey50]({site_type})[/grey50]")

            atoms_updated = 0
            bonds_updated = 0
            centers_updated = 0

            # 1. Update CENTERS
            if hasattr(site, 'centers') and site.centers:
                for center in site.centers:
                    old_chain = center.chain
                    old_resid = center.resid

                    new_chain = mapper.get_new_chain(old_chain)
                    new_resid = mapper.get_new_resnum(old_chain, old_resid)

                    if new_chain != old_chain or new_resid != old_resid:
                        center.chain = new_chain
                        center.resid = new_resid
                        centers_updated += 1

            # 2. Update ATOMS
            if hasattr(site, 'atoms') and site.atoms:
                for atom in site.atoms:
                    old_chain = atom.chain
                    old_resid = atom.resid

                    new_chain = mapper.get_new_chain(old_chain)
                    new_resid = mapper.get_new_resnum(old_chain, old_resid)

                    if new_chain != old_chain or new_resid != old_resid:
                        atom.chain = new_chain
                        atom.resid = new_resid
                        atoms_updated += 1

            # 3. Update BONDS (residue_info dicts)
            if hasattr(site, 'bonds') and site.bonds:
                for bond in site.bonds:
                    bond_changed = False

                    # Update atom1_residue_info
                    if hasattr(bond, 'atom1_residue_info') and bond.atom1_residue_info:
                        if 'chain' in bond.atom1_residue_info and 'resid' in bond.atom1_residue_info:
                            old_chain = bond.atom1_residue_info['chain']
                            old_resid = bond.atom1_residue_info['resid']
                            new_chain = mapper.get_new_chain(old_chain)
                            new_resid = mapper.get_new_resnum(old_chain, old_resid)

                            if new_chain != old_chain or new_resid != old_resid:
                                bond.atom1_residue_info['chain'] = new_chain
                                bond.atom1_residue_info['resid'] = new_resid
                                bond_changed = True

                    # Update atom2_residue_info
                    if hasattr(bond, 'atom2_residue_info') and bond.atom2_residue_info:
                        if 'chain' in bond.atom2_residue_info and 'resid' in bond.atom2_residue_info:
                            old_chain = bond.atom2_residue_info['chain']
                            old_resid = bond.atom2_residue_info['resid']
                            new_chain = mapper.get_new_chain(old_chain)
                            new_resid = mapper.get_new_resnum(old_chain, old_resid)

                            if new_chain != old_chain or new_resid != old_resid:
                                bond.atom2_residue_info['chain'] = new_chain
                                bond.atom2_residue_info['resid'] = new_resid
                                bond_changed = True

                    if bond_changed:
                        bonds_updated += 1

            # 4. Rebuild COORD_TO_PDB mapping (values only, coords are keys and never change)
            if hasattr(site, 'coord_to_pdb') and site.coord_to_pdb:
                for coords, pdb_info in site.coord_to_pdb.items():
                    if 'chain' in pdb_info and 'resid' in pdb_info:
                        old_chain = pdb_info['chain']
                        old_resid = pdb_info['resid']
                        new_chain = mapper.get_new_chain(old_chain)
                        new_resid = mapper.get_new_resnum(old_chain, old_resid)

                        if new_chain != old_chain or new_resid != old_resid:
                            pdb_info['chain'] = new_chain
                            pdb_info['resid'] = new_resid

            # 5. Rebuild RESIDUE_GROUPS (keys change from old to new chain/resid)
            if hasattr(site, 'residue_groups') and site.residue_groups:
                new_residue_groups = {}
                for (old_chain, old_resid, resname), atom_coords_list in site.residue_groups.items():
                    new_chain = mapper.get_new_chain(old_chain)
                    new_resid = mapper.get_new_resnum(old_chain, old_resid)
                    new_key = (new_chain, new_resid, resname)
                    new_residue_groups[new_key] = atom_coords_list
                site.residue_groups = new_residue_groups

            # Report updates for this site
            if atoms_updated > 0 or bonds_updated > 0 or centers_updated > 0:
                sites_updated += 1
                total_atoms_updated += atoms_updated
                total_bonds_updated += bonds_updated
                total_centers_updated += centers_updated
                self.console.print(
                    f"  [green]✓[/green] Updated: {centers_updated} centers, "
                    f"{atoms_updated} atoms, {bonds_updated} bonds"
                )
            else:
                self.console.print("  [grey50]No changes needed[/grey50]")

        self.console.print(
            f"\n[green]✓ Synchronized {sites_updated} RedoxSite(s)[/green]"
        )
        self.console.print(
            f"[green]  {total_centers_updated} centers, {total_atoms_updated} atoms, "
            f"{total_bonds_updated} bonds updated[/green]"
        )

        return {
            'sites_updated': sites_updated,
            'atoms_updated': total_atoms_updated,
            'bonds_updated': total_bonds_updated,
            'centers_updated': total_centers_updated
        }

    def analyze_pdb_segments(self, pdb_file: str) -> List[ChainSegment]:
        """
        Analyze PDB file to identify chain segments.

        Uses backbone atom presence (N, CA, C, O) to distinguish protein
        chains from hetero residues, rather than relying on residue names.

        Args:
            pdb_file: Path to PDB file

        Returns:
            List of ChainSegment objects
        """
        segments = []
        current_chain = None
        chain_data = {}  # chain_id -> {residues: set, has_backbone: bool, hetero: set, protein_res: set, cofactor_res: dict}

        with open(pdb_file, 'r') as f:
            for line in f:
                if not line.startswith(('ATOM', 'HETATM')):
                    continue

                chain_id = line[21:22].strip() or ' '
                res_num = int(line[22:26].strip())
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()

                if chain_id not in chain_data:
                    chain_data[chain_id] = {
                        'residues': set(),
                        'has_backbone': False,
                        'hetero': set(),
                        'protein_res': set(),
                        'cofactor_res': {},
                        'residue_backbone': {}  # res_num -> has_backbone
                    }

                chain_data[chain_id]['residues'].add(res_num)

                # Track which residues have backbone atoms
                if res_num not in chain_data[chain_id]['residue_backbone']:
                    chain_data[chain_id]['residue_backbone'][res_num] = False

                # Check for backbone atoms (only from ATOM records, not HETATM)
                # This prevents calcium ions (CA) from being misidentified as backbone
                if line.startswith('ATOM') and atom_name in self.BACKBONE_ATOMS:
                    chain_data[chain_id]['has_backbone'] = True
                    chain_data[chain_id]['residue_backbone'][res_num] = True

                # Track hetero residues
                if line.startswith('HETATM') and res_name not in ('HOH', 'WAT'):
                    chain_data[chain_id]['hetero'].add(res_name)

                # Track residue name for non-protein residues
                if res_num not in chain_data[chain_id]['cofactor_res']:
                    chain_data[chain_id]['cofactor_res'][res_num] = res_name

        # Convert to ChainSegment objects
        for chain_id, data in sorted(chain_data.items()):
            if not data['residues']:
                continue

            # Separate protein and cofactor residues
            protein_res = []
            cofactor_res = {}
            for res_num in data['residues']:
                if data['residue_backbone'].get(res_num, False):
                    protein_res.append(res_num)
                else:
                    # Get residue name from cofactor_res dict
                    res_name = data['cofactor_res'].get(res_num, 'UNK')
                    # Skip waters
                    if res_name not in ('HOH', 'WAT'):
                        cofactor_res[res_num] = res_name

            residues = sorted(data['residues'])
            segments.append(ChainSegment(
                chain_id=chain_id,
                start_residue=residues[0],
                end_residue=residues[-1],
                residue_count=len(residues),
                is_protein=data['has_backbone'],
                hetero_residues=sorted(data['hetero']),
                protein_residues=sorted(protein_res),
                cofactor_residues=cofactor_res
            ))

        return segments

    def configure_molecule_grouping(self, pdb_file: str) -> Optional[MoleculeConfiguration]:
        """
        Interactive configuration of molecule grouping.

        Args:
            pdb_file: Path to PDB file

        Returns:
            MoleculeConfiguration or None if cancelled
        """
        self.console.print(Panel.fit(
            "[bold cyan]PDB Molecule Configuration[/bold cyan]\n\n"
            "Configure how chains should be grouped into molecules for AMBER.\n"
            "This is critical for systems with intermolecular bonds (e.g., metal coordination).",
            border_style="cyan"
        ))

        # Analyze PDB structure
        segments = self.analyze_pdb_segments(pdb_file)

        if not segments:
            self.console.print("[red]No chains found in PDB file[/red]")
            return None

        # Display segments
        self._display_segments(segments)

        # Ask if user wants to configure grouping
        if not self._confirm("Do you want to configure molecule grouping?", default=True, description="Configure molecule grouping"):
            # Default: each chain is its own molecule
            molecules = [
                MoleculeGroup(
                    name=f"Molecule_{i+1}",
                    chain_ids=[seg.chain_id],
                    description=f"Chain {seg.chain_id}"
                )
                for i, seg in enumerate(segments)
            ]
            user_configured = False
        else:
            molecules = self._interactive_molecule_configuration(segments)
            user_configured = True

        if not molecules:
            return None

        # Create configuration
        config = MoleculeConfiguration(
            pdb_file=pdb_file,
            molecules=molecules,
            chain_segments=segments
        )

        # Only display and offer to save if user actively configured grouping
        if user_configured:
            self._display_configuration(config)

            if self._confirm("Save this configuration?", default=True, description="Save configuration"):
                config_file = self._save_configuration(config)
                if config_file:
                    self.console.print(f"[green]Configuration saved to {config_file}[/green]")

        return config

    def _display_segments(self, segments: List[ChainSegment]):
        """Display chain segments in a table."""
        table = Table(title="PDB Chain Segments")
        table.add_column("Chain", style="cyan")
        table.add_column("Type", style="yellow")
        table.add_column("Protein Range", style="green")
        table.add_column("Protein Count", justify="right", style="green")
        table.add_column("Cofactor Residues", style="magenta")

        for seg in segments:
            seg_type = "Protein" if seg.is_protein else "Hetero"
            protein_range = seg.get_protein_range_str()
            protein_count = str(len(seg.protein_residues)) if seg.protein_residues else "-"
            cofactor_str = seg.get_cofactor_str()

            table.add_row(
                seg.chain_id,
                seg_type,
                protein_range,
                protein_count,
                cofactor_str
            )

        self.console.print(table)

    def _interactive_molecule_configuration(self, segments: List[ChainSegment]) -> List[MoleculeGroup]:
        """Interactive molecule grouping configuration."""
        molecules = []
        available_chains = set(seg.chain_id for seg in segments)

        self.console.print("\n[bold]Define molecule groups[/bold]")
        self.console.print("Example: If chains A,B,C bond together with hemes D,E, they form one molecule.\n")

        molecule_num = 1
        while available_chains:
            self.console.print(f"\n[cyan]Available chains: {', '.join(sorted(available_chains))}[/cyan]")

            # Ask for chains in this molecule
            chain_input = self._prompt(
                f"Enter chain IDs for molecule {molecule_num} (comma-separated, or 'done')",
                default="done" if len(molecules) > 0 else None,
                description=f"Chain IDs for molecule {molecule_num}"
            )

            if chain_input.lower() == 'done':
                # Remaining chains each become their own molecule
                for chain_id in sorted(available_chains):
                    molecules.append(MoleculeGroup(
                        name=f"Molecule_{len(molecules) + 1}",
                        chain_ids=[chain_id],
                        description=f"Chain {chain_id}"
                    ))
                break

            # Parse chain IDs
            chain_ids = [c.strip() for c in chain_input.split(',')]
            invalid = [c for c in chain_ids if c not in available_chains]

            if invalid:
                self.console.print(f"[red]Invalid chains: {', '.join(invalid)}[/red]")
                continue

            # Get description
            description = self._prompt(
                "Description for this molecule (optional)",
                default=f"Chains {', '.join(chain_ids)}",
                description="Molecule description"
            )

            molecules.append(MoleculeGroup(
                name=f"Molecule_{molecule_num}",
                chain_ids=chain_ids,
                description=description
            ))

            # Remove used chains
            for chain_id in chain_ids:
                available_chains.remove(chain_id)

            molecule_num += 1

        return molecules

    def _display_configuration(self, config: MoleculeConfiguration):
        """Display final molecule configuration."""
        table = Table(title="Molecule Configuration")
        table.add_column("Molecule", style="cyan")
        table.add_column("Chains", style="yellow")
        table.add_column("Description", style="green")

        for mol in config.molecules:
            table.add_row(
                mol.name,
                ", ".join(mol.chain_ids),
                mol.description
            )

        self.console.print("\n")
        self.console.print(table)

    def _save_configuration(self, config: MoleculeConfiguration) -> Optional[str]:
        """Save configuration to JSON file."""
        pdb_path = Path(config.pdb_file)
        config_file = pdb_path.parent / f"{pdb_path.stem}_molecule_config.json"

        try:
            with open(config_file, 'w') as f:
                json.dump(config.to_dict(), f, indent=2)
            return str(config_file)
        except Exception as e:
            self.console.print(f"[red]Error saving configuration: {e}[/red]")
            return None

    def load_configuration(self, config_file: str) -> Optional[MoleculeConfiguration]:
        """Load configuration from JSON file."""
        try:
            with open(config_file, 'r') as f:
                data = json.load(f)
            return MoleculeConfiguration.from_dict(data)
        except Exception as e:
            self.console.print(f"[red]Error loading configuration: {e}[/red]")
            return None

    def reorder_pdb(self, pdb_file: str, config: MoleculeConfiguration,
                    output_file: str = None) -> Tuple[Optional[str], Optional[PDBReorderMapper]]:
        """
        Reorder PDB file based on molecule configuration with sequential renumbering.

        Renumbers residues sequentially within each chain and atom serials globally.
        Tracks all mappings for RedoxSite synchronization.

        Args:
            pdb_file: Input PDB file path
            config: Molecule configuration
            output_file: Output file path (default: input_reordered.pdb)

        Returns:
            Tuple of (output_file_path, mapper) or (None, None) on error
        """
        if output_file is None:
            pdb_path = Path(pdb_file)
            output_file = str(pdb_path.parent / f"{pdb_path.stem}_reordered.pdb")

        self.console.print(f"[cyan]Reordering and renumbering PDB file...[/cyan]")

        try:
            # Create mapper to track residue number changes
            mapper = PDBReorderMapper(self.console)

            # Read all PDB lines grouped by chain
            chain_atoms = {}  # chain_id -> list of (resnum, line)
            other_lines = []

            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        chain_id = line[21:22].strip() or ' '
                        res_num = int(line[22:26].strip())
                        if chain_id not in chain_atoms:
                            chain_atoms[chain_id] = []
                        chain_atoms[chain_id].append((res_num, line))
                    else:
                        other_lines.append(line)

            # Write reordered and renumbered PDB
            with open(output_file, 'w') as f:
                # Write header lines
                for line in other_lines:
                    if line.startswith(('ATOM', 'HETATM', 'TER', 'END')):
                        continue
                    f.write(line)

                # Global counters - sequential across ALL chains
                atom_serial = 1
                global_resnum = 1  # Sequential across all chains for tLEaP
                total_atoms = 0

                # Process each molecule
                for mol_idx, mol in enumerate(config.molecules):
                    mol_atom_count = 0

                    # Process each chain in this molecule
                    for chain_id in mol.chain_ids:
                        if chain_id not in chain_atoms:
                            continue

                        # Get all atoms for this chain, sorted by residue number
                        atoms_in_chain = chain_atoms[chain_id]

                        # Group by residue
                        residue_groups = {}  # old_resnum -> list of lines
                        for old_resnum, line in atoms_in_chain:
                            if old_resnum not in residue_groups:
                                residue_groups[old_resnum] = []
                            residue_groups[old_resnum].append(line)

                        # Renumber residues sequentially ACROSS ALL CHAINS
                        for old_resnum in sorted(residue_groups.keys()):
                            # Track mapping
                            mapper.add_residue_mapping(chain_id, old_resnum, global_resnum)

                            # Write all atoms in this residue
                            for line in residue_groups[old_resnum]:
                                # Renumber atom serial and residue number
                                new_line = (
                                    f"{line[:6]}{atom_serial:5d}{line[11:22]}"
                                    f"{global_resnum:4d}{line[26:]}"
                                )
                                f.write(new_line)
                                atom_serial += 1
                                mol_atom_count += 1
                                total_atoms += 1

                            global_resnum += 1  # Increment globally across all chains

                    # Add TER record after each molecule
                    if mol_atom_count > 0:
                        f.write("TER\n")

                # Write END record
                f.write("END\n")

            self.console.print(f"[green]✓ Reordered PDB saved to {output_file}[/green]")
            self.console.print(f"[green]  Wrote {total_atoms} atoms in {len(config.molecules)} molecules[/green]")

            # Log mapping summary
            mapper.log_summary()

            return output_file, mapper

        except Exception as e:
            self.console.print(f"[red]✗ Error reordering PDB: {e}[/red]")
            import traceback
            logger.error(f"PDB reordering error: {traceback.format_exc()}")
            return None, None
