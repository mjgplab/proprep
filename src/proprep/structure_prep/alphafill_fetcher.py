"""
AlphaFill Database Fetcher

Specialist module for retrieving enriched protein structures from AlphaFill Database.
AlphaFill structures are AlphaFold models enriched with transplanted ligands,
cofactors, and metal ions from homologous PDB structures.

Author: ProPrep Development Team
Date: 2025-12-14
"""

import logging
import os
import tempfile
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

import requests
from Bio.PDB import PDBParser, MMCIFParser

logger = logging.getLogger(__name__)


class AlphaFillFetcher:
    """
    Handles retrieval of enriched protein structures from AlphaFill Database.

    AlphaFill enriches AlphaFold predictions with ligands, cofactors, and metal ions
    transplanted from homologous PDB structures. This is particularly valuable for
    redox proteins, metalloproteins, and enzymes requiring cofactors.
    """

    BASE_URL = "https://alphafill.eu/v1/aff"

    def __init__(self, output_dir: Optional[str] = None, timeout: int = 60):
        """
        Initialize the AlphaFill fetcher.

        Args:
            output_dir: Directory to store downloaded structures.
                       If None, uses a temp directory.
            timeout: Request timeout in seconds (default: 60)
        """
        self.output_dir = output_dir or tempfile.mkdtemp(prefix="alphafill_")
        self.timeout = timeout
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        logger.debug(f"AlphaFill fetcher initialized with output_dir: {self.output_dir}")

    def check_availability(self, uniprot_id: str) -> Tuple[bool, Optional[Dict[str, Any]]]:
        """
        Check if AlphaFill structure is available for a UniProt ID.

        Args:
            uniprot_id: UniProt accession ID (e.g., 'P00441')

        Returns:
            Tuple of (available, metadata_dict)
            - available: True if structure exists in AlphaFill DB
            - metadata_dict: Dictionary with basic structure info if available
        """
        logger.info(f"Checking AlphaFill availability for {uniprot_id}")

        try:
            # Check status endpoint
            url = f"{self.BASE_URL}/{uniprot_id}/status"
            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 404:
                logger.debug(f"No AlphaFill structure found for {uniprot_id}")
                return False, None

            response.raise_for_status()
            status_data = response.json()

            # Check if structure is available (not just queued or processing)
            status = status_data.get("status", "unknown")

            if status not in ["finished", "available", "ready"]:
                logger.debug(f"AlphaFill structure for {uniprot_id} has status: {status}")
                return False, None

            # Get basic metadata
            metadata = {
                "uniprot_id": uniprot_id,
                "alphafill_id": status_data.get("id", f"AF-{uniprot_id}-F1"),
                "status": status,
                "has_transplants": status_data.get("has_transplants", False)
            }

            logger.info(f"Found AlphaFill structure for {uniprot_id}: {metadata['alphafill_id']}")
            return True, metadata

        except requests.exceptions.RequestException as e:
            logger.warning(f"Error checking AlphaFill availability for {uniprot_id}: {e}")
            return False, None
        except Exception as e:
            logger.warning(f"Unexpected error checking AlphaFill availability: {e}")
            return False, None

    def retrieve_structure(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Retrieve AlphaFill-enriched structure for a UniProt ID.

        Args:
            uniprot_id: UniProt accession ID (e.g., 'P00441')

        Returns:
            Dictionary containing:
                - pdb_file: Path to downloaded structure file
                - structure: BioPython Structure object
                - uniprot_id: Input UniProt ID
                - alphafill_id: AlphaFill entry ID (e.g., AF-P00441-F1)
                - metadata: Dictionary with structure metadata
                - transplants: List of transplanted ligands/cofactors
            Or None if retrieval failed
        """
        logger.info(f"Retrieving AlphaFill structure for UniProt ID: {uniprot_id}")

        # Check availability first
        available, basic_metadata = self.check_availability(uniprot_id)
        if not available:
            logger.error(f"No AlphaFill structure available for {uniprot_id}")
            return None

        try:
            # Download structure file
            structure_url = f"{self.BASE_URL}/{uniprot_id}"
            logger.debug(f"Downloading structure from {structure_url}")

            response = requests.get(structure_url, timeout=self.timeout)
            response.raise_for_status()

            # Detect file format by inspecting content
            # mmCIF files start with "data_" or contain "loop_"
            # PDB files start with records like "HEADER", "ATOM", etc.
            content_preview = response.content[:500].decode('utf-8', errors='ignore')

            if content_preview.strip().startswith('data_') or 'loop_' in content_preview:
                # mmCIF format
                file_ext = 'cif'
                parser = MMCIFParser(QUIET=True)
                logger.debug("Detected mmCIF format")
            else:
                # PDB format
                file_ext = 'pdb'
                parser = PDBParser(QUIET=True)
                logger.debug("Detected PDB format")

            structure_file = os.path.join(
                self.output_dir,
                f"alphafill_{uniprot_id}.{file_ext}"
            )

            with open(structure_file, 'wb') as f:
                f.write(response.content)

            logger.debug(f"Structure saved to {structure_file}")

            # Parse structure with BioPython
            try:
                structure = parser.get_structure(uniprot_id, structure_file)
            except Exception as parse_error:
                # If parsing fails with one format, try the other
                logger.debug(f"Failed to parse as {file_ext}, trying alternative format: {parse_error}")
                if file_ext == 'pdb':
                    parser = MMCIFParser(QUIET=True)
                    new_ext = 'cif'
                else:
                    parser = PDBParser(QUIET=True)
                    new_ext = 'pdb'

                # Re-save with correct extension
                new_structure_file = os.path.join(
                    self.output_dir,
                    f"alphafill_{uniprot_id}.{new_ext}"
                )
                os.rename(structure_file, new_structure_file)
                structure_file = new_structure_file
                logger.debug(f"Retrying parse with {new_ext} format")
                structure = parser.get_structure(uniprot_id, structure_file)

            # Get transplant metadata
            transplants = self.get_transplant_info(uniprot_id)

            # Build metadata dictionary
            metadata = basic_metadata.copy()
            if transplants:
                metadata['transplant_count'] = len(transplants)
                metadata['transplant_types'] = list(set(t.get('compound_name', 'unknown') for t in transplants))

            result = {
                "pdb_file": structure_file,
                "structure": structure,
                "uniprot_id": uniprot_id,
                "alphafill_id": metadata.get("alphafill_id", f"AF-{uniprot_id}-F1"),
                "metadata": metadata,
                "transplants": transplants or []
            }

            logger.info(f"Successfully retrieved AlphaFill structure for {uniprot_id}")
            if transplants:
                logger.info(f"  Found {len(transplants)} transplanted ligands/cofactors")

            return result

        except requests.exceptions.RequestException as e:
            logger.error(f"Network error retrieving AlphaFill structure for {uniprot_id}: {e}")
            return None
        except Exception as e:
            logger.error(f"Error retrieving AlphaFill structure for {uniprot_id}: {e}")
            logger.exception("Full traceback:")
            return None

    def get_transplant_info(self, uniprot_id: str) -> Optional[List[Dict[str, Any]]]:
        """
        Get detailed information about transplanted ligands and cofactors.

        Args:
            uniprot_id: UniProt accession ID

        Returns:
            List of dictionaries containing transplant information:
                - compound_id: Ligand/cofactor ID (e.g., 'HEM', 'FES', 'ZN')
                - compound_name: Full name of the compound
                - source_pdb: PDB ID from which it was transplanted
                - source_chain: Chain in source PDB
                - confidence: Transplant confidence score (if available)
            Or None if metadata retrieval failed
        """
        logger.debug(f"Fetching transplant metadata for {uniprot_id}")

        try:
            url = f"{self.BASE_URL}/{uniprot_id}/json"
            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 404:
                logger.debug(f"No transplant metadata found for {uniprot_id}")
                return None

            response.raise_for_status()
            metadata = response.json()

            # Parse transplant information
            transplants = []

            # AlphaFill JSON structure may vary; adapt as needed based on actual API response
            # Common fields: transplants, compounds, ligands
            if 'transplants' in metadata:
                for transplant in metadata['transplants']:
                    transplants.append({
                        'compound_id': transplant.get('compound_id', 'unknown'),
                        'compound_name': transplant.get('compound_name', 'Unknown'),
                        'source_pdb': transplant.get('pdb_id', 'unknown'),
                        'source_chain': transplant.get('chain_id', 'unknown'),
                        'confidence': transplant.get('confidence', None),
                        'identity': transplant.get('identity', None),
                        'type': transplant.get('type', 'ligand')
                    })
            elif 'compounds' in metadata:
                # Alternative structure
                for compound in metadata['compounds']:
                    transplants.append({
                        'compound_id': compound.get('id', 'unknown'),
                        'compound_name': compound.get('name', 'Unknown'),
                        'source_pdb': compound.get('source', {}).get('pdb_id', 'unknown'),
                        'source_chain': compound.get('source', {}).get('chain', 'unknown'),
                        'confidence': compound.get('score', None),
                        'type': compound.get('type', 'ligand')
                    })

            logger.debug(f"Found {len(transplants)} transplants for {uniprot_id}")
            return transplants if transplants else None

        except requests.exceptions.RequestException as e:
            logger.warning(f"Error fetching transplant metadata for {uniprot_id}: {e}")
            return None
        except Exception as e:
            logger.warning(f"Unexpected error parsing transplant metadata: {e}")
            return None

    def search_uniprot_by_gene(self, gene_name: str, organism: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Search UniProt for gene name and return potential matches.

        This is a convenience wrapper that uses the same logic as AlphaFoldFetcher
        to maintain consistency in the interface.

        Args:
            gene_name: Gene name (e.g., 'SOD1', 'TP53')
            organism: Optional organism filter (e.g., 'human', 'Homo sapiens')

        Returns:
            List of dictionaries with UniProt entry information
        """
        logger.info(f"Searching UniProt for gene: {gene_name}, organism: {organism}")

        # Build UniProt REST API query
        query_parts = [f"gene:{gene_name}"]
        if organism:
            # Normalize organism name
            if organism.lower() in ['human', 'homo sapiens']:
                query_parts.append("organism_name:Homo sapiens")
            elif organism.lower() in ['mouse', 'mus musculus']:
                query_parts.append("organism_name:Mus musculus")
            elif organism.lower() in ['rat', 'rattus norvegicus']:
                query_parts.append("organism_name:Rattus norvegicus")
            else:
                query_parts.append(f"organism_name:{organism}")

        query = " AND ".join(query_parts)

        # Query UniProt REST API
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            "query": query,
            "format": "json",
            "size": 10
        }

        try:
            response = requests.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()

            results = []
            for entry in data.get("results", []):
                uniprot_id = entry.get("primaryAccession")
                protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown")
                organism_name = entry.get("organism", {}).get("scientificName", "Unknown")
                gene_names = [g.get("geneName", {}).get("value") for g in entry.get("genes", []) if "geneName" in g]
                reviewed = entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"
                length = entry.get("sequence", {}).get("length", 0)

                results.append({
                    "uniprot_id": uniprot_id,
                    "protein_name": protein_name,
                    "organism": organism_name,
                    "gene_names": gene_names,
                    "reviewed": reviewed,
                    "length": length
                })

            logger.info(f"Found {len(results)} UniProt entries for gene {gene_name}")
            return results

        except Exception as e:
            logger.error(f"Error searching UniProt: {e}")
            return []
