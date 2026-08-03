"""
UniProt Search Module - UniProt Database API Integration

This module provides search and data retrieval functionality for the UniProt
protein database, allowing users to find proteins by gene name, protein name,
or UniProt ID, and retrieve associated PDB structures along with rich metadata.

Author: ProPrep Development Team
"""

import logging
import requests
from typing import List, Dict, Optional, Any
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class UniProtEntry:
    """Represents a UniProt protein entry with metadata."""
    accession: str  # Primary UniProt accession (e.g., P00441)
    entry_name: str  # Entry name (e.g., SODC_HUMAN)
    protein_name: str  # Recommended protein name
    gene_name: Optional[str]  # Primary gene name
    organism_scientific: str  # Scientific organism name
    organism_common: Optional[str]  # Common organism name
    sequence: str  # Amino acid sequence
    sequence_length: int  # Number of residues
    molecular_weight: float  # Molecular weight in Da

    # Functional information
    function: Optional[str]  # Protein function description
    catalytic_activity: List[Dict[str, str]]  # EC numbers and reactions
    cofactors: List[str]  # Required cofactors (metals, etc.)

    # Structural features
    metal_binding_sites: List[Dict[str, Any]]  # Metal binding sites
    disulfide_bonds: List[Dict[str, Any]]  # Disulfide bonds
    active_sites: List[Dict[str, Any]]  # Active sites
    binding_sites: List[Dict[str, Any]]  # Other binding sites
    modifications: List[Dict[str, Any]]  # Post-translational modifications

    # Clinical/biological context
    subcellular_location: List[str]  # Cellular localization
    disease_associations: List[Dict[str, str]]  # Associated diseases

    # Cross-references
    pdb_ids: List[Dict[str, str]]  # Associated PDB structures with metadata
    cross_references: Dict[str, int]  # Counts of other database links


class UniProtSearcher:
    """
    Search the UniProt protein database using the REST API.

    Provides methods to search by gene name, protein name, or UniProt ID,
    and retrieve comprehensive protein metadata including associated PDB structures.
    """

    SEARCH_API_URL = "https://rest.uniprot.org/uniprotkb/search"
    ENTRY_API_URL = "https://rest.uniprot.org/uniprotkb"

    # Common organism mappings
    ORGANISM_MAP = {
        "human": "9606",
        "mouse": "10090",
        "rat": "10116",
        "yeast": "559292",  # S. cerevisiae
        "e.coli": "83333",
        "ecoli": "83333",
        "fly": "7227",  # D. melanogaster
        "worm": "6239",  # C. elegans
        "zebrafish": "7955",
        "arabidopsis": "3702",
    }

    def __init__(self, timeout: int = 30):
        """
        Initialize UniProt searcher.

        Args:
            timeout: HTTP request timeout in seconds
        """
        self.timeout = timeout

    def normalize_organism(self, organism: str) -> str:
        """
        Normalize organism name to taxonomy ID.

        Args:
            organism: Common or scientific name, or taxonomy ID

        Returns:
            Taxonomy ID as string
        """
        if not organism:
            return organism

        # Check if already a taxonomy ID
        if organism.isdigit():
            return organism

        organism_lower = organism.lower().strip()
        return self.ORGANISM_MAP.get(organism_lower, organism)

    def search(
        self,
        query: str,
        organism: Optional[str] = None,
        max_results: int = 10
    ) -> List[str]:
        """
        Search UniProt by gene name, protein name, or UniProt ID.

        Args:
            query: Search query (gene name, protein name, or UniProt ID)
            organism: Optional organism filter (common name, scientific name, or taxonomy ID)
            max_results: Maximum number of results to return

        Returns:
            List of UniProt accession IDs matching the search
        """
        logger.debug(f"Searching UniProt for '{query}'" +
                    (f" in {organism}" if organism else ""))

        # Build search query
        # Try multiple search strategies for better results
        search_queries = []

        # Check if it looks like a UniProt accession
        if len(query) in [6, 10] and query.replace('_', '').replace('-', '').isalnum():
            search_queries.append(f"accession:{query}")

        # Gene name search
        search_queries.append(f"gene:{query}")

        # Protein name search
        search_queries.append(f"protein_name:{query}")

        # Combine with OR
        combined_query = " OR ".join([f"({q})" for q in search_queries])

        # Prefer reviewed entries (Swiss-Prot) over unreviewed (TrEMBL)
        combined_query = f"({combined_query}) AND reviewed:true"

        # Add organism filter if provided
        if organism:
            org_id = self.normalize_organism(organism)
            if org_id.isdigit():
                combined_query = f"({combined_query}) AND organism_id:{org_id}"
            else:
                combined_query = f"({combined_query}) AND organism_name:{org_id}"

        params = {
            "query": combined_query,
            "format": "json",
            "size": max_results
        }

        try:
            response = requests.get(
                self.SEARCH_API_URL,
                params=params,
                timeout=self.timeout
            )
            response.raise_for_status()

            data = response.json()
            results = data.get("results", [])

            accessions = [entry["primaryAccession"] for entry in results]
            logger.debug(f"Found {len(accessions)} UniProt entries")

            return accessions

        except requests.exceptions.RequestException as e:
            logger.error(f"UniProt search failed: {e}")
            return []

    def get_entry(self, accession: str) -> Optional[UniProtEntry]:
        """
        Retrieve detailed information for a UniProt entry.

        Args:
            accession: UniProt accession ID

        Returns:
            UniProtEntry object with full metadata, or None if not found
        """
        logger.debug(f"Fetching UniProt entry {accession}")

        try:
            url = f"{self.ENTRY_API_URL}/{accession}.json"
            response = requests.get(url, timeout=self.timeout)
            response.raise_for_status()

            data = response.json()

            # Extract all metadata
            entry = self._parse_entry(data)
            return entry

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to fetch UniProt entry {accession}: {e}")
            return None

    def _parse_entry(self, data: Dict[str, Any]) -> UniProtEntry:
        """Parse UniProt JSON data into UniProtEntry object."""

        # Basic info
        accession = data.get("primaryAccession", "")
        entry_name = data.get("uniProtkbId", "")

        # Protein names
        protein_name = "Unknown"
        if "proteinDescription" in data:
            desc = data["proteinDescription"]
            if "recommendedName" in desc:
                protein_name = desc["recommendedName"]["fullName"]["value"]

        # Gene names
        gene_name = None
        if "genes" in data and data["genes"]:
            gene_data = data["genes"][0]
            if "geneName" in gene_data:
                gene_name = gene_data["geneName"]["value"]

        # Organism
        organism_scientific = ""
        organism_common = None
        if "organism" in data:
            org = data["organism"]
            organism_scientific = org.get("scientificName", "")
            organism_common = org.get("commonName")

        # Sequence
        sequence = ""
        sequence_length = 0
        molecular_weight = 0.0
        if "sequence" in data:
            seq = data["sequence"]
            sequence = seq.get("value", "")
            sequence_length = seq.get("length", 0)
            molecular_weight = seq.get("molWeight", 0.0)

        # Function
        function = None
        catalytic_activity = []
        cofactors = []
        subcellular_location = []
        disease_associations = []

        if "comments" in data:
            for comment in data["comments"]:
                comment_type = comment.get("commentType")

                if comment_type == "FUNCTION":
                    texts = comment.get("texts", [])
                    if texts:
                        function = texts[0]["value"]

                elif comment_type == "CATALYTIC_ACTIVITY":
                    reaction = comment.get("reaction", {})
                    if reaction:
                        catalytic_activity.append({
                            "name": reaction.get("name", ""),
                            "ec": reaction.get("ecNumber", "")
                        })

                elif comment_type == "COFACTOR":
                    if "cofactors" in comment:
                        for cofactor in comment["cofactors"]:
                            cofactors.append(cofactor.get("name", ""))

                elif comment_type == "SUBCELLULAR_LOCATION":
                    if "subcellularLocations" in comment:
                        for loc in comment["subcellularLocations"]:
                            if "location" in loc:
                                subcellular_location.append(loc["location"]["value"])

                elif comment_type == "DISEASE":
                    disease = comment.get("disease", {})
                    if disease:
                        disease_associations.append({
                            "id": disease.get("diseaseId", ""),
                            "name": disease.get("diseaseAccession", ""),
                            "description": disease.get("description", "")
                        })

        # Features
        metal_binding_sites = []
        disulfide_bonds = []
        active_sites = []
        binding_sites = []
        modifications = []

        if "features" in data:
            for feature in data["features"]:
                ftype = feature.get("type")
                location = feature.get("location", {})
                start = location.get("start", {}).get("value")
                end = location.get("end", {}).get("value", start)
                description = feature.get("description", "")

                feature_data = {
                    "start": start,
                    "end": end,
                    "description": description
                }

                if ftype == "Metal binding":
                    ligand = feature.get("ligand", {})
                    feature_data["ligand"] = ligand.get("name", "")
                    metal_binding_sites.append(feature_data)

                elif ftype == "Disulfide bond":
                    disulfide_bonds.append(feature_data)

                elif ftype == "Active site":
                    active_sites.append(feature_data)

                elif ftype == "Binding site":
                    binding_sites.append(feature_data)

                elif ftype == "Modified residue":
                    modifications.append(feature_data)

        # PDB cross-references
        pdb_ids = []
        cross_references = {}

        if "uniProtKBCrossReferences" in data:
            for ref in data["uniProtKBCrossReferences"]:
                db = ref["database"]
                cross_references[db] = cross_references.get(db, 0) + 1

                if db == "PDB":
                    pdb_data = {
                        "id": ref["id"],
                        "method": "",
                        "resolution": "",
                        "chains": ""
                    }

                    # Extract properties
                    if "properties" in ref:
                        for prop in ref["properties"]:
                            key = prop["key"]
                            value = prop["value"]
                            if key == "Method":
                                pdb_data["method"] = value
                            elif key == "Resolution":
                                pdb_data["resolution"] = value
                            elif key == "Chains":
                                pdb_data["chains"] = value

                    pdb_ids.append(pdb_data)

        return UniProtEntry(
            accession=accession,
            entry_name=entry_name,
            protein_name=protein_name,
            gene_name=gene_name,
            organism_scientific=organism_scientific,
            organism_common=organism_common,
            sequence=sequence,
            sequence_length=sequence_length,
            molecular_weight=molecular_weight,
            function=function,
            catalytic_activity=catalytic_activity,
            cofactors=cofactors,
            metal_binding_sites=metal_binding_sites,
            disulfide_bonds=disulfide_bonds,
            active_sites=active_sites,
            binding_sites=binding_sites,
            modifications=modifications,
            subcellular_location=subcellular_location,
            disease_associations=disease_associations,
            pdb_ids=pdb_ids,
            cross_references=cross_references
        )

    def search_and_get(
        self,
        query: str,
        organism: Optional[str] = None
    ) -> Optional[UniProtEntry]:
        """
        One-step search and retrieve operation.

        Searches for a protein and returns the first matching entry with full metadata.

        Args:
            query: Search query (gene name, protein name, or UniProt ID)
            organism: Optional organism filter

        Returns:
            UniProtEntry object for the first match, or None if no results
        """
        # Search for matches
        accessions = self.search(query, organism, max_results=1)

        if not accessions:
            return None

        # Get detailed entry for first match
        return self.get_entry(accessions[0])
