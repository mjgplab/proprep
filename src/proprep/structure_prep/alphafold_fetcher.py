"""
AlphaFold Database Fetcher

Specialist module for retrieving protein structures from the AlphaFold Database.
Uses BioPython's alphafold_db module to query and download structures.

Author: ProPrep Development Team
Date: 2025-11-08
"""

import logging
import os
import tempfile
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.alphafold_db import get_predictions, get_structural_models_for

logger = logging.getLogger(__name__)


class AlphaFoldFetcher:
    """
    Handles retrieval of protein structures from AlphaFold Database.

    Uses BioPython's alphafold_db module to access pre-computed AlphaFold
    predictions via the EBI AlphaFold Database API.
    """

    def __init__(self, output_dir: Optional[str] = None):
        """
        Initialize the AlphaFold fetcher.

        Args:
            output_dir: Directory to store downloaded structures.
                       If None, uses a temp directory.
        """
        self.output_dir = output_dir or tempfile.mkdtemp(prefix="alphafold_")
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        logger.debug(f"AlphaFold fetcher initialized with output_dir: {self.output_dir}")

    def check_availability(self, uniprot_id: str) -> Tuple[bool, Optional[Dict[str, Any]]]:
        """
        Check if AlphaFold structure is available for a UniProt ID.

        Args:
            uniprot_id: UniProt accession ID (e.g., 'P00520')

        Returns:
            Tuple of (available, metadata_dict)
            - available: True if structure exists in AlphaFold DB
            - metadata_dict: Dictionary with structure info if available
        """
        try:
            predictions = list(get_predictions(uniprot_id))

            if not predictions:
                logger.debug(f"No AlphaFold predictions found for {uniprot_id}")
                return False, None

            # Get first prediction (usually there's only one)
            prediction = predictions[0]

            # Extract metadata
            metadata = {
                "uniprot_id": uniprot_id,
                "alphafold_id": prediction.get("entryId", f"AF-{uniprot_id}-F1"),
                "organism": prediction.get("organismScientificName", "Unknown"),
                "protein_name": prediction.get("uniprotDescription", "Unknown"),
                "length": prediction.get("sequenceLength", 0),
                "version": prediction.get("latestVersion", "unknown"),
                "created_date": prediction.get("created", "unknown")
            }

            logger.info(f"Found AlphaFold structure for {uniprot_id}: {metadata['alphafold_id']}")
            return True, metadata

        except Exception as e:
            logger.warning(f"Error checking AlphaFold availability for {uniprot_id}: {e}")
            return False, None

    def retrieve_structure(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Retrieve AlphaFold structure for a UniProt ID.

        Args:
            uniprot_id: UniProt accession ID (e.g., 'P00520')

        Returns:
            Dictionary containing:
                - pdb_file: Path to downloaded PDB file
                - structure: BioPython Structure object
                - uniprot_id: Input UniProt ID
                - alphafold_id: AlphaFold entry ID
                - metadata: Dictionary with structure metadata
                - confidence: Confidence metrics (if available)
            Or None if retrieval failed
        """
        logger.info(f"Retrieving AlphaFold structure for UniProt ID: {uniprot_id}")

        # Check availability first
        available, metadata = self.check_availability(uniprot_id)
        if not available:
            logger.error(f"No AlphaFold structure available for {uniprot_id}")
            return None

        try:
            # Download structure using BioPython
            mmcif_parser = MMCIFParser(QUIET=True)

            structures = list(get_structural_models_for(
                uniprot_id,
                mmcif_parser=mmcif_parser,
                directory=self.output_dir
            ))

            if not structures:
                logger.error(f"Failed to download structure for {uniprot_id}")
                return None

            # Get first structure
            structure = structures[0]

            # Find the downloaded CIF file
            alphafold_id = metadata.get("alphafold_id", f"AF-{uniprot_id}-F1")
            cif_file = None
            for file in Path(self.output_dir).glob(f"*{uniprot_id}*.cif"):
                cif_file = str(file)
                break

            if not cif_file:
                logger.warning(f"CIF file not found, searching for PDB file")
                for file in Path(self.output_dir).glob(f"*{uniprot_id}*.pdb"):
                    cif_file = str(file)
                    break

            # Extract confidence metrics from structure (pLDDT scores in B-factor column)
            confidence = self._extract_confidence_metrics(structure)

            result = {
                "pdb_file": cif_file,
                "structure": structure,
                "uniprot_id": uniprot_id,
                "alphafold_id": alphafold_id,
                "metadata": metadata,
                "confidence": confidence
            }

            logger.info(f"Successfully retrieved {alphafold_id}")
            logger.debug(f"Structure file: {cif_file}")
            logger.debug(f"Mean pLDDT: {confidence.get('mean_plddt', 'N/A')}")

            return result

        except Exception as e:
            logger.error(f"Error retrieving AlphaFold structure for {uniprot_id}: {e}")
            return None

    def _extract_confidence_metrics(self, structure) -> Dict[str, Any]:
        """
        Extract confidence metrics from AlphaFold structure.

        AlphaFold stores pLDDT scores in the B-factor column.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary with confidence metrics
        """
        plddt_scores = []

        try:
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            # pLDDT is stored in B-factor
                            plddt_scores.append(atom.get_bfactor())

            if plddt_scores:
                mean_plddt = sum(plddt_scores) / len(plddt_scores)

                # Categorize confidence regions
                very_high = sum(1 for s in plddt_scores if s > 90)
                confident = sum(1 for s in plddt_scores if 70 <= s <= 90)
                low = sum(1 for s in plddt_scores if 50 <= s < 70)
                very_low = sum(1 for s in plddt_scores if s < 50)

                return {
                    "mean_plddt": round(mean_plddt, 2),
                    "residue_count": len(plddt_scores),
                    "very_high_confidence": very_high,
                    "confident": confident,
                    "low_confidence": low,
                    "very_low_confidence": very_low,
                    "quality_category": self._categorize_quality(mean_plddt)
                }
        except Exception as e:
            logger.warning(f"Could not extract confidence metrics: {e}")

        return {"mean_plddt": None}

    def _categorize_quality(self, mean_plddt: float) -> str:
        """Categorize overall model quality based on mean pLDDT."""
        if mean_plddt > 90:
            return "Very High"
        elif mean_plddt > 70:
            return "Confident"
        elif mean_plddt > 50:
            return "Low"
        else:
            return "Very Low"

    def search_uniprot_by_gene(self, gene_name: str, organism: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Search UniProt for gene name and return potential matches.

        Args:
            gene_name: Gene name (e.g., 'SRC', 'TP53')
            organism: Optional organism filter (e.g., 'human', 'Homo sapiens')

        Returns:
            List of dictionaries with UniProt entry information
        """
        import requests

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
            "size": 10  # Limit to top 10 results
        }

        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            results = []
            for entry in data.get("results", []):
                # Extract relevant information
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

    def retrieve_from_blast_results(self, blast_results: List[Dict[str, Any]], max_structures: int = 10) -> Dict[str, Dict[str, Any]]:
        """
        Retrieve AlphaFold structures for BLAST hits.

        Args:
            blast_results: List of BLAST hit dictionaries (must contain 'accession' or 'uniprot_id')
            max_structures: Maximum number of structures to retrieve

        Returns:
            Dictionary mapping UniProt IDs to structure data
        """
        logger.info(f"Retrieving AlphaFold structures for {len(blast_results)} BLAST hits")

        structures = {}
        retrieved_count = 0

        for idx, hit in enumerate(blast_results):
            if retrieved_count >= max_structures:
                break

            # Extract UniProt ID from BLAST result
            uniprot_id = hit.get("uniprot_id") or hit.get("accession")

            if not uniprot_id:
                logger.debug(f"BLAST hit {idx+1} has no UniProt ID, skipping")
                continue

            # Check availability
            available, metadata = self.check_availability(uniprot_id)

            if available:
                # Store availability info for now, retrieve later if user selects
                structures[uniprot_id] = {
                    "available": True,
                    "metadata": metadata,
                    "blast_rank": idx + 1,
                    "blast_identity": hit.get("identity", 0),
                    "blast_evalue": hit.get("evalue", 0)
                }
                retrieved_count += 1
            else:
                structures[uniprot_id] = {
                    "available": False,
                    "blast_rank": idx + 1
                }

        logger.info(f"Found {retrieved_count} AlphaFold structures available out of {len(blast_results)} BLAST hits")
        return structures
