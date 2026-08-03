"""
Restraint Mask Generator for ProPrep

This module generates AMBER restraint masks for redox sites based on metal site 
definitions and transformed structures. It analyzes the actual structure to identify
all atoms that are part of redox sites and creates appropriate restraint masks for
molecular dynamics simulations.

Author: ProPrep Development Team
Version: 1.0.0
"""

import logging
import re
from typing import Dict, List, Set, Any, Optional, Tuple
from rich.console import Console

from proprep.utils.prompts import prompt_with_context, confirm_with_context

# Setup logging
logger = logging.getLogger(__name__)


class RestraintMaskGenerator:
    """
    Generator for AMBER restraint masks based on redox site analysis.
    
    This class analyzes metal site definitions and transformed structures to
    identify all atoms that should be restrained during MD minimization.
    It works with transformed structures that have non-standard residue/atom names.
    """
    
    def __init__(self, console: Optional[Console] = None, processor: Optional[Any] = None):
        """
        Initialize the restraint mask generator.

        Args:
            console: Rich console for output (optional)
            processor: Processor object for session recording (optional)
        """
        self.console = console or Console()
        self.processor = processor
        self.metal_sites = []
        self.structure_lines = []
        self.redox_atoms = set()  # Set of (chain, resid, atom_name) tuples
        
    def generate_masks(self, metal_sites: List[Any], pdb_content: str, interactive: bool = True) -> Optional[Dict[str, Any]]:
        """
        Generate restraint masks for redox sites.
        
        Args:
            metal_sites: List of MetalSite objects from workspace
            pdb_content: PDB content string of transformed structure
            interactive: Whether to prompt user for options
            
        Returns:
            Dict containing restraint mask and info, or None if failed
        """
        try:
            from rich.prompt import Prompt, Confirm
            from rich.table import Table
            
            self.metal_sites = metal_sites
            self.structure_lines = pdb_content.strip().split('\n')
            self.redox_atoms = set()
            self.redox_residues = set()  # Set of (chain, resid) tuples
            
            # Step 1: Interactive options
            include_full_residues = False
            restrain_by_residues = False
            
            if interactive:
                self.console.print("\n[bold cyan]===== Restraint Mask Generation Options =====[/bold cyan]")
                
                # First choice: scope
                self.console.print("Choose restraint scope:")
                self.console.print("  1. Metal and coordinating atoms only")
                self.console.print("  2. All atoms in metal and ligand residues")

                scope_choice = prompt_with_context(
                    processor=self.processor,
                    prompt="Select scope",
                    choices=["1", "2"],
                    default="1",
                    module="Restraint Mask Generator",
                    description="Select restraint scope",
                    options_map={
                        "1": "Metal and coordinating atoms only",
                        "2": "All atoms in metal and ligand residues"
                    }
                )
                include_full_residues = (scope_choice == "2")

                # Second choice: granularity
                self.console.print("\nChoose restraint granularity:")
                self.console.print("  1. By residues (restrain entire residues - simpler masks)")
                self.console.print("  2. By individual atoms (select specific atoms - complex masks)")

                granularity_choice = prompt_with_context(
                    processor=self.processor,
                    prompt="Select granularity",
                    choices=["1", "2"],
                    default="1",
                    module="Restraint Mask Generator",
                    description="Select restraint granularity",
                    options_map={
                        "1": "By residues (simpler masks)",
                        "2": "By individual atoms (complex masks)"
                    }
                )
                restrain_by_residues = (granularity_choice == "1")
            
            # Step 2: Analyze metal sites to identify redox atoms
            if include_full_residues:
                self._analyze_metal_sites_full_residues()
            else:
                self._analyze_metal_sites()
            
            # Step 3: Generate restraint mask based on chosen method
            if restrain_by_residues:
                # Residue-based restraints
                restraint_mask = self._create_residue_based_mask(interactive)
                # For residue-based restraints, count atoms in the selected residues
                redox_atom_names = self._count_atoms_in_selected_residues()
            else:
                # Atom-based restraints  
                redox_atom_names = self._parse_structure_for_redox_atoms()
                
                # Interactive editing for atom-based
                if interactive and redox_atom_names:
                    redox_atom_names = self._interactive_edit_atoms(redox_atom_names)
                
                redox_resids = [r['resid'] for r in self._get_residue_info_with_names()]
                restraint_mask = self._create_restraint_mask(redox_atom_names, redox_resids)
            
            if not restraint_mask:
                self.console.print("[yellow]Warning: No redox atoms found for restraint mask[/yellow]")
                # Return backbone-only mask as fallback
                restraint_mask = "@CA,C,O,N&!:WAT"
            
            # Step 4: Compile result info
            method_name = ""
            if restrain_by_residues:
                method_name = "residue_based_" + ("full_residues" if include_full_residues else "coordinating_atoms")
            else:
                method_name = "atom_based_" + ("full_residues" if include_full_residues else "coordinating_atoms")
                
            result = {
                "restraint_mask": restraint_mask,
                "info": {
                    "metal_sites_count": len(self.metal_sites),
                    "included_atoms": list(redox_atom_names) if redox_atom_names else [],
                    "redox_residues": self._get_redox_residues(),
                    "generation_method": method_name,
                    "restraint_by_residues": restrain_by_residues,
                    "user_edited": interactive and (len(redox_atom_names) > 0 or restrain_by_residues)
                }
            }
            
            return result
            
        except Exception as e:
            logger.error(f"Error generating restraint masks: {e}")
            self.console.print(f"[red]Error generating restraint masks: {e}[/red]")
            return None
    
    def _analyze_metal_sites(self):
        """
        Analyze metal sites to identify atoms that should be restrained.
        """
        for site in self.metal_sites:
            try:
                # Add metal atom itself
                if hasattr(site, 'metal_chain') and hasattr(site, 'metal_resid') and hasattr(site, 'metal_name'):
                    if site.metal_chain and site.metal_resid and site.metal_name:
                        self.redox_atoms.add((site.metal_chain, str(site.metal_resid), site.metal_name))
                
                # Add coordinating atoms from ligands
                if hasattr(site, 'ligands') and site.ligands:
                    for ligand in site.ligands:
                        if isinstance(ligand, dict):
                            chain = ligand.get('chain')
                            resid = ligand.get('resid') 
                            atom_name = ligand.get('atom_name')
                            
                            if chain and resid and atom_name:
                                self.redox_atoms.add((chain, str(resid), atom_name))
                
                # For heme sites, add key porphyrin atoms
                if hasattr(site, 'metal_resname') and site.metal_resname:
                    if site.metal_resname.upper() in ['HEM', 'HEME', 'HEC']:
                        heme_chain = site.metal_chain
                        heme_resid = str(site.metal_resid)
                        
                        # Key heme atoms that define the porphyrin ring structure
                        key_heme_atoms = [
                            'NA', 'NB', 'NC', 'ND',  # Nitrogen atoms coordinating to iron
                            'C1A', 'C2A', 'C3A', 'C4A',  # Pyrrole A carbons
                            'C1B', 'C2B', 'C3B', 'C4B',  # Pyrrole B carbons  
                            'C1C', 'C2C', 'C3C', 'C4C',  # Pyrrole C carbons
                            'C1D', 'C2D', 'C3D', 'C4D',  # Pyrrole D carbons
                            'CHA', 'CHB', 'CHC', 'CHD'   # Methine bridge carbons
                        ]
                        
                        for atom in key_heme_atoms:
                            self.redox_atoms.add((heme_chain, heme_resid, atom))
                            
                logger.debug(f"Analyzed metal site {site.site_id}: found {len(self.redox_atoms)} redox atoms")
                
            except Exception as e:
                logger.warning(f"Error analyzing metal site: {e}")
                continue
    
    def _analyze_metal_sites_full_residues(self):
        """
        Analyze metal sites to identify all atoms in metal and ligand residues.
        """
        for site in self.metal_sites:
            try:
                # Add metal residue
                if hasattr(site, 'metal_chain') and hasattr(site, 'metal_resid'):
                    if site.metal_chain and site.metal_resid:
                        self.redox_residues.add((site.metal_chain, str(site.metal_resid)))
                
                # Add coordinating ligand residues
                if hasattr(site, 'ligands') and site.ligands:
                    for ligand in site.ligands:
                        if isinstance(ligand, dict):
                            chain = ligand.get('chain')
                            resid = ligand.get('resid')
                            
                            if chain and resid:
                                self.redox_residues.add((chain, str(resid)))
                                
                logger.debug(f"Analyzed metal site {site.site_id}: found {len(self.redox_residues)} redox residues")
                
            except Exception as e:
                logger.warning(f"Error analyzing metal site: {e}")
                continue
        
        # Now collect all atoms from these residues
        self._collect_atoms_from_residues()
    
    def _collect_atoms_from_residues(self):
        """
        Collect all atoms from the identified redox residues.
        """
        for line in self.structure_lines:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
                
            try:
                # Parse PDB line format
                chain = line[21:22].strip()
                resid = line[22:26].strip()
                atom_name = line[12:16].strip()
                
                # Check if this atom is in our redox residues
                res_key = (chain, resid)
                if res_key in self.redox_residues:
                    self.redox_atoms.add((chain, resid, atom_name))
                    logger.debug(f"Added atom: {atom_name} from residue {chain}:{resid}")
                
            except Exception as e:
                logger.debug(f"Error parsing structure line: {e}")
                continue
    
    def _interactive_edit_atoms(self, atom_names: Set[str]) -> Set[str]:
        """
        Allow user to interactively edit the list of atoms to include.
        
        Args:
            atom_names: Set of atom names found in structure
            
        Returns:
            Modified set of atom names
        """
        from rich.prompt import Prompt, Confirm
        from rich.table import Table
        
        if not atom_names:
            return atom_names
            
        self.console.print(f"\n[bold]Found {len(atom_names)} redox atoms for restraints:[/bold]")
        
        # Group atoms by residue for display
        atoms_by_residue = self._group_atoms_by_residue(atom_names)
        
        # Create a flat list for indexing while preserving visual grouping
        display_atoms = []  # List of (index, atom_name, residue_info)
        index = 1
        
        # Display atoms grouped by residue in a table
        table = Table(title="Redox Atoms for Restraints (Grouped by Residue)")
        table.add_column("#", style="cyan", width=4)
        table.add_column("Residue", style="magenta", width=12)
        table.add_column("Atom Name", style="green")
        table.add_column("Include", style="yellow", width=8)
        
        for (chain, resid), residue_atoms in atoms_by_residue.items():
            residue_label = f"{chain}:{resid}"
            
            # Sort atoms within each residue alphabetically
            sorted_residue_atoms = sorted(residue_atoms)
            
            for i, atom in enumerate(sorted_residue_atoms):
                # Only show residue label on first atom of each residue
                display_residue = residue_label if i == 0 else ""
                table.add_row(str(index), display_residue, atom, "✓")
                display_atoms.append((index, atom, residue_label))
                index += 1
        
        self.console.print(table)

        if not confirm_with_context(
            processor=self.processor,
            prompt="\nWould you like to edit this list?",
            default=False,
            module="Restraint Mask Generator",
            description="Edit restraint atom list"
        ):
            return atom_names
        
        # Allow user to remove atoms
        edited_atoms = set(atom_names)
        
        self.console.print("\n[bold]Edit restraint atoms:[/bold]")
        self.console.print("Enter atom numbers to REMOVE (e.g., '1,3,5' or '1-3,5')")
        self.console.print("Press Enter to keep all atoms")
        
        remove_input = prompt_with_context(
            self.processor, "Atoms to remove",
            module="Restraint Mask Generator",
            description="Atoms to remove from mask",
        ).strip()
        
        if remove_input:
            try:
                # Parse removal indices
                indices_to_remove = self._parse_atom_indices(remove_input, len(display_atoms))
                
                for idx in sorted(indices_to_remove, reverse=True):
                    if 1 <= idx <= len(display_atoms):
                        _, atom_to_remove, residue_info = display_atoms[idx - 1]
                        edited_atoms.discard(atom_to_remove)
                        self.console.print(f"[yellow]Removed: {atom_to_remove} (from {residue_info})[/yellow]")
                
                self.console.print(f"\n[green]Final restraint atoms ({len(edited_atoms)}):[/green] {', '.join(sorted(edited_atoms))}")
                
            except Exception as e:
                self.console.print(f"[red]Error parsing input: {e}[/red]")
                self.console.print("[yellow]Keeping all atoms[/yellow]")
                edited_atoms = atom_names
        
        return edited_atoms
    
    def _group_atoms_by_residue(self, atom_names: Set[str]) -> Dict[Tuple[str, str], List[str]]:
        """
        Group atoms by their residue (chain, resid).
        
        Args:
            atom_names: Set of atom names
            
        Returns:
            Dictionary mapping (chain, resid) to list of atom names
        """
        atoms_by_residue = {}
        
        # Map atom names back to their residues using redox_atoms
        for chain, resid, atom_name in self.redox_atoms:
            if atom_name in atom_names:
                res_key = (chain, resid)
                if res_key not in atoms_by_residue:
                    atoms_by_residue[res_key] = []
                atoms_by_residue[res_key].append(atom_name)
        
        # Sort residues by chain then resid for consistent display
        return dict(sorted(atoms_by_residue.items(), key=lambda x: (x[0][0], int(x[0][1]) if x[0][1].isdigit() else x[0][1])))
    
    
    def _parse_atom_indices(self, input_str: str, max_index: int) -> Set[int]:
        """
        Parse atom indices from user input (supports ranges like '1-3,5,7-9').
        
        Args:
            input_str: User input string
            max_index: Maximum valid index
            
        Returns:
            Set of indices to remove
        """
        indices = set()
        
        for part in input_str.split(','):
            part = part.strip()
            if '-' in part:
                # Range like '1-3'
                try:
                    start, end = map(int, part.split('-', 1))
                    indices.update(range(start, end + 1))
                except ValueError:
                    continue
            else:
                # Single number
                try:
                    indices.add(int(part))
                except ValueError:
                    continue
        
        # Filter valid indices
        return {idx for idx in indices if 1 <= idx <= max_index}
    
    def _create_residue_based_mask(self, interactive: bool = True) -> str:
        """
        Create AMBER restraint mask based on entire residues.
        
        Args:
            interactive: Whether to allow interactive editing
            
        Returns:
            AMBER restraint mask string
        """
        from rich.table import Table
        from rich.prompt import Confirm
        
        # Get residue information with names
        residue_info = self._get_residue_info_with_names()
        
        if interactive and residue_info:
            # Display residues grouped by name and ID
            self.console.print(f"\n[bold]Found {len(residue_info)} redox residues for restraints:[/bold]")
            
            table = Table(title="Redox Residues for Restraints")
            table.add_column("#", style="cyan", width=4)
            table.add_column("Res Name", style="magenta", width=8)
            table.add_column("Chain:ResID", style="green")
            table.add_column("Include", style="yellow", width=8)
            
            for i, res_info in enumerate(residue_info, 1):
                table.add_row(
                    str(i),
                    res_info['resname'],
                    f"{res_info['chain']}:{res_info['resid']}",
                    "✓"
                )
            
            self.console.print(table)

            # Allow removal of entire residues
            if confirm_with_context(
                processor=self.processor,
                prompt="\nWould you like to edit this list?",
                default=False,
                module="Restraint Mask Generator",
                description="Edit restraint residue list"
            ):
                residue_info = self._interactive_edit_residues(residue_info)
        
        # Create mask from selected residues
        if not residue_info:
            return "@CA,C,O,N&!:WAT"
        
        # Build residue-based mask
        residue_specs = []
        for res_info in residue_info:
            residue_specs.append(f":{res_info['resid']}")
        
        # Combine with backbone
        backbone_mask = "@CA,C,O,N&!:WAT"
        residue_mask = "|".join(residue_specs)
        
        combined_mask = f"{backbone_mask}|{residue_mask}"
        return combined_mask
    
    def _get_residue_info_with_names(self) -> List[Dict[str, str]]:
        """
        Get residue information including residue names from PDB structure.
        
        Returns:
            List of dicts with resname, chain, resid info
        """
        residue_info = []
        seen_residues = set()
        
        # Parse structure to get residue names
        for line in self.structure_lines:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
                
            try:
                # Parse PDB line format  
                resname = line[17:20].strip()
                chain = line[21:22].strip()
                resid = line[22:26].strip()
                
                # Check if this is one of our redox residues
                res_key = (chain, resid)
                if res_key in self.redox_residues and res_key not in seen_residues:
                    residue_info.append({
                        'resname': resname,
                        'chain': chain,
                        'resid': resid
                    })
                    seen_residues.add(res_key)
                    
            except Exception as e:
                logger.debug(f"Error parsing structure line for residue names: {e}")
                continue
        
        # Sort by resname then chain then resid
        return sorted(residue_info, key=lambda x: (x['resname'], x['chain'], int(x['resid']) if x['resid'].isdigit() else x['resid']))
    
    def _interactive_edit_residues(self, residue_info: List[Dict[str, str]]) -> List[Dict[str, str]]:
        """
        Allow user to remove entire residues from restraints.
        
        Args:
            residue_info: List of residue information dicts
            
        Returns:
            Modified list of residue information
        """
        self.console.print("\n[bold]Edit residue list:[/bold]")
        self.console.print("Enter residue numbers to REMOVE (e.g., '1,3,5' or '1-3,5')")
        self.console.print("Press Enter to keep all residues")
        
        remove_input = prompt_with_context(
            self.processor, "Residues to remove",
            module="Restraint Mask Generator",
            description="Residues to remove from mask",
        ).strip()
        
        if remove_input:
            try:
                # Parse removal indices
                indices_to_remove = self._parse_atom_indices(remove_input, len(residue_info))
                
                # Remove residues (in reverse order to maintain indices)
                for idx in sorted(indices_to_remove, reverse=True):
                    if 1 <= idx <= len(residue_info):
                        removed = residue_info.pop(idx - 1)
                        self.console.print(f"[yellow]Removed: {removed['resname']} {removed['chain']}:{removed['resid']}[/yellow]")
                
                self.console.print(f"\n[green]Final residues ({len(residue_info)}):[/green]")
                for res in residue_info:
                    self.console.print(f"  {res['resname']} {res['chain']}:{res['resid']}")
                    
            except Exception as e:
                self.console.print(f"[red]Error parsing input: {e}[/red]")
                self.console.print("[yellow]Keeping all residues[/yellow]")
        
        return residue_info
    
    def _parse_structure_for_redox_atoms(self) -> Set[str]:
        """
        Parse the transformed structure to find actual atom names for redox atoms.
        
        Returns:
            Set of atom names found in the structure
        """
        found_atoms = set()
        
        for line in self.structure_lines:
            if not line.startswith(('ATOM  ', 'HETATM')):
                continue
                
            try:
                # Parse PDB line format
                chain = line[21:22].strip()
                resid = line[22:26].strip()
                atom_name = line[12:16].strip()
                
                # Check if this atom is in our redox atoms set
                atom_key = (chain, resid, atom_name)
                if atom_key in self.redox_atoms:
                    found_atoms.add(atom_name)
                    logger.debug(f"Found redox atom: {atom_name} in {chain}:{resid}")
                
            except Exception as e:
                logger.debug(f"Error parsing structure line: {e}")
                continue
        
        return found_atoms
    
    def _create_restraint_mask(self, redox_atom_names: Set[str],
                               redox_resids: Optional[List[str]] = None) -> str:
        """
        Create AMBER restraint mask including backbone and redox atoms.
        
        Args:
            redox_atom_names: Set of redox atom names found in structure
            
        Returns:
            AMBER restraint mask string
        """
        # Start with backbone atoms (excluding water)
        backbone_mask = "@CA,C,O,N&!:WAT"
        
        if not redox_atom_names:
            return backbone_mask
        
        # Sort atom names for consistent output
        sorted_atoms = sorted(redox_atom_names)
        atom_part = "@" + ",".join(sorted_atoms)
        # Scope the redox-atom selection to the redox-site residues so a
        # generic name (CA/CB) doesn't match every same-named atom system-wide.
        if redox_resids:
            res_part = ":" + ",".join(
                str(r) for r in sorted(
                    set(redox_resids),
                    key=lambda x: int(x) if str(x).lstrip("-").isdigit() else 1_000_000,
                )
            )
            redox_mask = f"({res_part})&({atom_part})"
        else:
            redox_mask = atom_part
        
        # Combine backbone and redox atom masks
        combined_mask = f"{backbone_mask}|{redox_mask}"
        
        return combined_mask
    
    def _get_redox_residues(self) -> List[Dict[str, Any]]:
        """
        Get information about redox residues.
        
        Returns:
            List of dictionaries with redox residue info
        """
        redox_residues = []
        
        # Group redox atoms by residue
        residue_atoms = {}
        for chain, resid, atom_name in self.redox_atoms:
            res_key = (chain, resid)
            if res_key not in residue_atoms:
                residue_atoms[res_key] = []
            residue_atoms[res_key].append(atom_name)
        
        # Create residue info
        for (chain, resid), atoms in residue_atoms.items():
            redox_residues.append({
                "chain": chain,
                "resid": resid,
                "atom_count": len(atoms),
                "atoms": sorted(atoms)
            })
        
        return redox_residues

    def _count_atoms_in_selected_residues(self):
        """
        Count atoms in the selected redox residues for residue-based restraints.
        
        Returns:
            set: Set of atom identifiers in the selected residues
        """
        atom_names = set()
        
        # Check if we have a structure to work with
        if not hasattr(self, 'structure') or not self.structure:
            return atom_names
            
        # Get the redox residues that were selected
        redox_residues = self._get_redox_residues()
        
        try:
            # Count atoms in each selected redox residue
            for residue_info in redox_residues:
                chain_id = residue_info["chain"]
                resid = residue_info["resid"]
                
                # Find the residue in the structure
                for model in self.structure:
                    for chain in model:
                        if chain.id == chain_id:
                            for residue in chain:
                                if residue.id[1] == resid:  # residue.id is (' ', resid, icode)
                                    # Add all atoms in this residue
                                    for atom in residue:
                                        atom_names.add(f"{chain_id}:{resid}:{atom.name}")
                                    break
                            break
                    break
                        
        except Exception as e:
            self.console.print(f"[yellow]Warning: Could not count atoms in residues: {e}[/yellow]")
            
        return atom_names