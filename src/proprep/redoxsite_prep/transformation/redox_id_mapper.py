#!/usr/bin/env python3
"""
RedoxSite ID Mapping Module

Handles residue ID space management for RedoxSite transformations.
Ensures transformation products don't conflict with existing residues.

Author: Claude Code Implementation
"""

import logging
from typing import Any, Dict, List, Tuple

from rich.console import Console
from rich.table import Table

logger = logging.getLogger(__name__)

class RedoxSiteIDMapper:
    """Manages residue ID space allocation for RedoxSite transformations"""
    
    def __init__(self, console: Console = None, verbose: bool = True):
        self.console = console or Console()
        self.verbose = verbose
    
    def generate_id_mapping(self, redox_sites: List, space_requirements: Dict[str, Any], 
                          conflicts: List[Dict[str, Any]]) -> Dict[Tuple[str, int], int]:
        """
        Generate ID mapping following the legacy metal site pattern
        
        Args:
            redox_sites: List of RedoxSite objects
            space_requirements: Space analysis results from analyzer
            conflicts: List of detected conflicts
            
        Returns:
            Mapping from (chain, original_id) to new_id
        """
        
        # Always generate mapping to follow legacy pattern, even without conflicts
        # This ensures consistent behavior with the legacy system
        
        # Extract all residues from structure to find safe ranges
        existing_residues = self._extract_existing_residues(redox_sites)
        
        # Process each chain separately and generate mappings
        by_chain = self._group_requirements_by_chain(space_requirements)
        
        id_mapping = {}
        reserved_ranges = {}
        
        for chain, chain_reqs in by_chain.items():
            chain_mapping, chain_reserved = self._generate_chain_mapping(
                chain, chain_reqs, existing_residues.get(chain, set())
            )
            id_mapping.update(chain_mapping)
            reserved_ranges.update(chain_reserved)
        
        # Display mapping table in legacy format
        self._display_legacy_id_mapping_table(id_mapping, space_requirements)
        
        # Display reserved ID ranges
        self._display_reserved_ranges(reserved_ranges, space_requirements)
        
        return id_mapping
    
    def _extract_existing_residues(self, redox_sites: List) -> Dict[str, set]:
        """Extract all existing residue IDs by chain"""
        existing_by_chain = {}
        
        for site in redox_sites:
            for atom in site.atoms:
                chain = atom.chain
                if chain not in existing_by_chain:
                    existing_by_chain[chain] = set()
                existing_by_chain[chain].add(atom.resid)
        
        return existing_by_chain
    
    def _group_requirements_by_chain(self, space_requirements: Dict[str, Any]) -> Dict[str, List]:
        """Group space requirements by chain"""
        by_chain = {}
        
        for site_id, req in space_requirements.items():
            chain = req['anchor_chain']
            if chain not in by_chain:
                by_chain[chain] = []
            by_chain[chain].append((site_id, req))
        
        return by_chain
    
    def _generate_chain_mapping(self, chain: str, chain_reqs: List, 
                              existing_residues: set) -> Tuple[Dict[Tuple[str, int], int], Dict[str, List[int]]]:
        """Generate ID mapping for a single chain following legacy pattern"""
        mapping = {}
        reserved_ranges = {}
        
        # Sort requirements by original anchor ID
        chain_reqs.sort(key=lambda x: x[1]['anchor_id'])
        
        # Find the highest existing residue ID
        max_existing = max(existing_residues) if existing_residues else 0
        
        # Start assigning new IDs after the highest existing ID
        # Leave some buffer space to avoid immediate conflicts
        next_available_id = max_existing + 4
        
        for site_id, req in chain_reqs:
            original_id = req['anchor_id']
            required_count = req['required_count']
            
            # Map original ID to new ID
            new_id = next_available_id
            mapping[(chain, original_id)] = new_id
            
            # Track reserved range for this site
            reserved_range = list(range(new_id, new_id + required_count))
            reserved_ranges[site_id] = reserved_range
            
            # Move to next available block
            next_available_id += required_count
        
        return mapping, reserved_ranges
    
    def _display_legacy_id_mapping_table(self, id_mapping: Dict[Tuple[str, int], int],
                                       space_requirements: Dict[str, Any]):
        """Display ID mapping table in legacy format"""
        if not id_mapping or not self.verbose:
            return

        table = Table(title="Redox Site Residue ID Mapping")
        table.add_column("Redox Site", style="cyan")
        table.add_column("Original ID", style="yellow") 
        table.add_column("New ID", style="green")
        table.add_column("Site Type", style="blue")
        
        # Create site ID to requirement mapping for display
        site_to_req = {}
        for site_id, req in space_requirements.items():
            anchor_key = (req['anchor_chain'], req['anchor_id'])
            site_to_req[anchor_key] = (site_id, req)
        
        # Add rows to table
        for (chain, original_id), new_id in sorted(id_mapping.items()):
            if (chain, original_id) in site_to_req:
                site_id, req = site_to_req[(chain, original_id)]
                # Generate site name in legacy format
                site_name = f"filtered_struct_{site_id.split('_')[-1] if '_' in site_id else site_id}"
                table.add_row(
                    site_name,
                    f"{chain}:{original_id}",
                    str(new_id),
                    req.get('transformer', 'Unknown')
                )
        
        self.console.print(table)
    
    def _display_reserved_ranges(self, reserved_ranges: Dict[str, List[int]],
                               space_requirements: Dict[str, Any]):
        """Display reserved ID ranges in legacy format"""
        if not reserved_ranges or not self.verbose:
            return

        self.console.print("\n[bold]Using new residue numbering:[/bold]")
        
        for site_id, id_range in reserved_ranges.items():
            if site_id in space_requirements:
                req = space_requirements[site_id]
                chain = req['anchor_chain']
                original_id = req['anchor_id']
                
                range_str = f"{id_range[0]}, {id_range[1]}, {id_range[2]}" if len(id_range) >= 3 else ", ".join(map(str, id_range))
                self.console.print(f"residues {range_str} reserved for redox site #{(chain, original_id)}")
    
    def apply_id_mapping(self, pdb_lines: List[str], id_mapping: Dict[Tuple[str, int], int]) -> List[str]:
        """
        Apply ID mapping to PDB structure lines
        
        Args:
            pdb_lines: Original PDB content
            id_mapping: Mapping from (chain, original_id) to new_id
            
        Returns:
            Modified PDB lines with updated residue IDs
        """
        if not id_mapping:
            return pdb_lines
        
        modified_lines = []
        
        for line in pdb_lines:
            if line.startswith(('ATOM', 'HETATM')):
                # Extract current chain and residue ID
                chain = line[21:22]
                try:
                    resid = int(line[22:26].strip())
                except ValueError:
                    # Keep line unchanged if residue ID can't be parsed
                    modified_lines.append(line)
                    continue
                
                # Check if this residue needs mapping
                mapping_key = (chain, resid)
                if mapping_key in id_mapping:
                    new_resid = id_mapping[mapping_key]
                    # Update the residue ID in the line
                    modified_line = line[:22] + f"{new_resid:>4}" + line[26:]
                    modified_lines.append(modified_line)
                else:
                    modified_lines.append(line)
            else:
                modified_lines.append(line)
        
        return modified_lines