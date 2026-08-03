import gzip
import json
import os
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import requests


class CCDParser:
    """Class to retrieve and parse Chemical Component Dictionary data."""
    
    # CCD base URLs
    CCD_BASE_URL = "https://files.wwpdb.org/pub/pdb/data/monomers"
    CCD_MAIN_URL = f"{CCD_BASE_URL}/components.cif.gz"
    CCD_VARIANTS_URL = f"{CCD_BASE_URL}/aa-variants-v1.cif.gz"
    
    # Cache directory
    CACHE_DIR = Path.home() / ".ccd_cache"
    
    def __init__(self, use_cache: bool = True, cache_expiry_days: int = 7):
        """Initialize the CCD parser with cache settings."""
        self.use_cache = use_cache
        self.cache_expiry_days = cache_expiry_days
        
        if use_cache:
            os.makedirs(self.CACHE_DIR, exist_ok=True)
    
    def _download_file(self, url: str, filename: str) -> str:
        """Download a file with caching support."""
        cache_path = self.CACHE_DIR / filename
        
        # Check if we have a cached version that's not expired
        if self.use_cache and cache_path.exists():
            cache_age = time.time() - os.path.getmtime(cache_path)
            if cache_age < self.cache_expiry_days * 24 * 60 * 60:
                # Using cached file silently
                with gzip.open(cache_path, 'rb') as f:
                    return f.read().decode('utf-8')
        
        # Download the file
        print(f"Downloading {url}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        # Decompress the content
        content = gzip.decompress(response.content)
        
        # Cache the content if requested
        if self.use_cache:
            with open(cache_path, 'wb') as f:
                f.write(response.content)
        
        return content.decode('utf-8')
    
    def get_residue_data(self, residue_name: str) -> Dict[str, Any]:
        """Get complete data for a residue from the CCD."""
        residue_name = residue_name.upper()
        
        # Get the main CCD content
        try:
            ccd_content = self._download_file(self.CCD_MAIN_URL, "components.cif.gz")
            residue_block = self._extract_residue_block(ccd_content, residue_name)
            
            if not residue_block:
                return {"error": f"Residue '{residue_name}' not found in the CCD"}
            
            # Parse the residue data
            result = self._parse_residue_data(residue_block, residue_name)
            
            # Check for protonation variants if it's an amino acid
            if result.get('component_type') and ('PEPTIDE' in result.get('component_type', '') or 
                                              'D-PEPTIDE' in result.get('component_type', '')):
                var_content = self._download_file(self.CCD_VARIANTS_URL, "aa-variants-v1.cif.gz")
                variants = self._extract_protonation_variants(var_content, residue_name)
                if variants:
                    result['protonation_variants'] = variants
            
            return result
        
        except requests.exceptions.RequestException as e:
            return {"error": f"Network error: {str(e)}"}
        except Exception as e:
            return {"error": f"Processing error: {str(e)}"}
    
    def _extract_residue_block(self, content: str, residue_name: str) -> Optional[str]:
        """Extract the data block for a specific residue."""
        pattern = rf"data_{residue_name}\s+"
        match = re.search(pattern, content)
        
        if not match:
            return None
        
        # Find the start of the data block
        start_pos = match.start()
        
        # Find the next data block or end of file
        next_match = re.search(r"data_\S+\s+", content[start_pos + len(residue_name) + 6:])
        if next_match:
            end_pos = start_pos + len(residue_name) + 6 + next_match.start()
            return content[start_pos:end_pos]
        else:
            return content[start_pos:]
    
    def _extract_protonation_variants(self, content: str, residue_name: str) -> List[Dict[str, str]]:
        """Extract information about protonation variants for a residue."""
        # Find all variant IDs that include this residue
        pattern = rf"_chem_comp\.id\s+(\S+_{residue_name}_\S+)"
        variant_ids = re.findall(pattern, content)
        
        variants = []
        for variant_id in variant_ids:
            # Extract the data block for this variant
            variant_block = self._extract_residue_block(content, variant_id)
            if variant_block:
                variant_info = self._parse_variant_block(variant_block, variant_id)
                variants.append(variant_info)
        
        return variants
    
    def _parse_variant_block(self, block: str, variant_id: str) -> Dict[str, str]:
        """Parse a protonation variant block."""
        variant = {'id': variant_id}
        
        # Extract basic info
        self._extract_simple_fields(block, variant, '_chem_comp.')
        
        return variant
    
    def _parse_residue_data(self, block: str, residue_name: str) -> Dict[str, Any]:
        """Parse a residue data block into a structured dictionary."""
        result = {'residue_id': residue_name}
        
        # Extract basic information (simple fields)
        self._extract_simple_fields(block, result, '_chem_comp.')
        
        # Extract atom information
        atoms = self._extract_loop_data(block, '_chem_comp_atom.')
        if atoms:
            result['atoms'] = atoms
        
        # Extract bond information
        bonds = self._extract_loop_data(block, '_chem_comp_bond.')
        if bonds:
            result['bonds'] = bonds
        
        # Extract descriptors
        descriptors = self._extract_loop_data(block, '_pdbx_chem_comp_descriptor.')
        if descriptors:
            result['descriptors'] = descriptors
            
            # Extract specific descriptors (SMILES, InChI)
            for desc in descriptors:
                if desc.get('type') == 'SMILES_CANONICAL':
                    result['smiles_canonical'] = desc.get('descriptor')
                elif desc.get('type') == 'SMILES':
                    result['smiles'] = desc.get('descriptor')
                elif desc.get('type') == 'InChI':
                    result['inchi'] = desc.get('descriptor')
                elif desc.get('type') == 'InChIKey':
                    result['inchikey'] = desc.get('descriptor')
        
        # Extract identifiers (systematic names)
        identifiers = self._extract_loop_data(block, '_pdbx_chem_comp_identifier.')
        if identifiers:
            result['identifiers'] = identifiers
            
            # Extract systematic name
            for ident in identifiers:
                if ident.get('type') == 'SYSTEMATIC NAME':
                    result['systematic_name'] = ident.get('identifier')
        
        return result
    
    def _extract_simple_fields(self, block: str, result: Dict[str, Any], prefix: str) -> None:
        """Extract simple key-value fields from a data block."""
        # Find all fields with this prefix
        pattern = rf"{re.escape(prefix)}(\S+)\s+([^\n]+)"
        matches = re.findall(pattern, block)
        
        for field, value in matches:
            # Clean up the value
            value = value.strip()
            
            # Remove quotes if present
            if value.startswith('"') and value.endswith('"'):
                value = value[1:-1]
            
            # Convert "?" to None
            if value == "?":
                value = None
            # Convert numeric values
            elif value.isdigit():
                value = int(value)
            elif re.match(r'^-?\d+\.\d+$', value):
                value = float(value)
            
            # Store in the result dictionary
            result[field] = value
    
    def _extract_loop_data(self, block: str, prefix: str) -> List[Dict[str, Any]]:
        """Extract looped data (like atoms, bonds) from a block."""
        # First, find the loop section that starts with this prefix
        prefix_escaped = re.escape(prefix)
        loop_pattern = rf"loop_\s+{prefix_escaped}\S+"
        loop_match = re.search(loop_pattern, block)
        
        if not loop_match:
            return []
        
        loop_start = loop_match.start()
        
        # Find the next loop or the end of the block
        next_loop = re.search(r"loop_\s+", block[loop_start + 6:])
        if next_loop:
            loop_end = loop_start + 6 + next_loop.start()
            loop_section = block[loop_start:loop_end]
        else:
            loop_section = block[loop_start:]
        
        # Extract the column headers
        header_pattern = rf"{prefix_escaped}(\S+)"
        headers = re.findall(header_pattern, loop_section)
        
        # Find where the data rows start
        header_section = ""
        for header in headers:
            header_line = f"{prefix}{header}"
            header_section += header_line + "\n"
        
        data_start = loop_section.find(header_section) + len(header_section)
        data_section = loop_section[data_start:].strip()
        
        # Parse the data rows
        result = []
        current_row = []
        in_quotes = False
        current_value = ""
        
        # Add a space at the end to ensure the last value is processed
        for char in data_section + " ":
            if char == '"' and not in_quotes:
                in_quotes = True
                current_value += char
            elif char == '"' and in_quotes:
                in_quotes = False
                current_value += char
            elif char.isspace() and not in_quotes:
                if current_value:
                    # Clean the value
                    if current_value.startswith('"') and current_value.endswith('"'):
                        current_value = current_value[1:-1]
                    
                    # Convert "?" to None
                    if current_value == "?":
                        current_value = None
                    
                    current_row.append(current_value)
                    current_value = ""
                
                if char == '\n' and current_row:
                    # End of row, create a dictionary
                    if len(current_row) >= len(headers):
                        row_dict = {}
                        for i, header in enumerate(headers):
                            if i < len(current_row):
                                row_dict[header] = current_row[i]
                        result.append(row_dict)
                    current_row = []
            else:
                current_value += char
        
        # Process any remaining row
        if current_row and len(current_row) >= len(headers):
            row_dict = {}
            for i, header in enumerate(headers):
                if i < len(current_row):
                    row_dict[header] = current_row[i]
            result.append(row_dict)
        
        return result


def get_residue_info(residue_name: str, use_cache: bool = True) -> Dict[str, Any]:
    """Get comprehensive information about a residue from the CCD."""
    parser = CCDParser(use_cache=use_cache)
    return parser.get_residue_data(residue_name)


def display_residue_info(result: Dict[str, Any], verbose: bool = True) -> None:
    """Display residue information in a readable format."""
    if "error" in result:
        print(f"ERROR: {result['error']}")
        return
    
    print("\n" + "="*80)
    print(f"CHEMICAL COMPONENT DICTIONARY: {result.get('residue_id', 'Unknown')}")
    print("="*80)
    
    # Basic information
    print("\n-- BASIC INFORMATION --")
    print(f"Name: {result.get('name', 'N/A')}")
    print(f"Type: {result.get('type', 'N/A')}")
    print(f"PDBx Type: {result.get('pdbx_type', 'N/A')}")
    print(f"Formula: {result.get('formula', 'N/A')}")
    print(f"Formula Weight: {result.get('formula_weight', 'N/A')}")
    print(f"Formal Charge: {result.get('pdbx_formal_charge', 'N/A')}")
    
    if 'one_letter_code' in result:
        print(f"One-letter Code: {result['one_letter_code']}")
    
    if 'three_letter_code' in result:
        print(f"Three-letter Code: {result['three_letter_code']}")
    
    if 'mon_nstd_parent_comp_id' in result and result['mon_nstd_parent_comp_id']:
        print(f"Parent Standard Residue: {result['mon_nstd_parent_comp_id']}")
    
    if 'pdbx_synonyms' in result and result['pdbx_synonyms']:
        print(f"Synonyms: {result['pdbx_synonyms']}")
    
    # Chemical descriptors
    print("\n-- CHEMICAL DESCRIPTORS --")
    print(f"SMILES: {result.get('smiles', 'N/A')}")
    print(f"SMILES Canonical: {result.get('smiles_canonical', 'N/A')}")
    print(f"InChI: {result.get('inchi', 'N/A')}")
    print(f"InChIKey: {result.get('inchikey', 'N/A')}")
    
    if 'systematic_name' in result:
        print(f"Systematic Name: {result['systematic_name']}")
    
    # Protonation variants
    if 'protonation_variants' in result:
        print("\n-- PROTONATION VARIANTS --")
        for i, variant in enumerate(result['protonation_variants']):
            print(f"Variant {i+1}: {variant['id']}")
            print(f"  Name: {variant.get('name', 'N/A')}")
            if 'pdbx_description' in variant:
                print(f"  Description: {variant['pdbx_description']}")
            if 'formula' in variant:
                print(f"  Formula: {variant['formula']}")
            if 'pdbx_formal_charge' in variant:
                print(f"  Formal Charge: {variant['pdbx_formal_charge']}")
            print()
    
    # Atoms
    atoms = result.get('atoms', [])
    atom_count = len(atoms)
    print(f"\n-- ATOM INFORMATION ({atom_count} atoms) --")
    
    if verbose and atoms:
        # Determine which fields to display
        fields = ['atom_id', 'type_symbol', 'charge', 'pdbx_aromatic_flag']
        
        # Add coordinates if they exist
        coord_fields = [field for field in atoms[0].keys() 
                       if 'Cartn' in field and field in atoms[0]]
        
        # Create column list with numbers
        all_fields = fields + coord_fields
        print("Column key:")
        for i, field in enumerate(all_fields, 1):
            print(f"  {i}. {field}")
        
        # Print header with numbers
        header = "  ".join(f"{i:<8}" for i in range(1, len(all_fields) + 1))
        print("\n" + header)
        print("-" * len(header))
        
        # Print atom data
        for atom in atoms:
            line_parts = []
            for field in all_fields:
                value = atom.get(field, 'N/A')
                if value is None:
                    value = 'N/A'
                line_parts.append(f"{value:<8}")
            
            print("  ".join(line_parts))
    
    # Bonds
    bonds = result.get('bonds', [])
    bond_count = len(bonds)
    print(f"\n-- BOND INFORMATION ({bond_count} bonds) --")
    
    if verbose and bonds:
        # Determine which fields to display
        fields = ['atom_id_1', 'atom_id_2', 'value_order', 'pdbx_aromatic_flag']
        
        # Create column list with numbers
        available_fields = [f for f in fields if f in bonds[0]]
        print("Column key:")
        for i, field in enumerate(available_fields, 1):
            print(f"  {i}. {field}")
        
        # Print header with numbers
        header = "  ".join(f"{i:<8}" for i in range(1, len(available_fields) + 1))
        print("\n" + header)
        print("-" * len(header))
        
        # Print bond data
        for bond in bonds:
            line_parts = []
            for field in available_fields:
                value = bond.get(field, 'N/A')
                if value is None:
                    value = 'N/A'
                line_parts.append(f"{value:<8}")
            
            print("  ".join(line_parts))

def to_json(result: Dict[str, Any]) -> str:
    """Convert the result to a JSON string."""
    return json.dumps(result, indent=2)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Retrieve Chemical Component Dictionary information')
    parser.add_argument('residue', help='Residue code (e.g., ALA, ATP)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Show detailed information')
    parser.add_argument('--output', '-o', choices=['display', 'json'], default='display',
                        help='Output format (default: display)')
    parser.add_argument('--no-cache', action='store_true', help='Disable file caching')
    
    args = parser.parse_args()
    
    # Get the residue information
    result = get_residue_info(args.residue, use_cache=not args.no_cache)
    
    # Output the information
    if args.output == 'display':
        display_residue_info(result, verbose=args.verbose)
    elif args.output == 'json':
        print(to_json(result))
