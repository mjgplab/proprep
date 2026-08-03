"""
PDB Search Module - RCSB PDB Search API Integration

This module provides search functionality for the RCSB Protein Data Bank,
allowing users to find structures by protein name, gene name, organism, and
other criteria without needing to know the exact PDB ID.

Author: ProPrep Development Team
"""

import logging
import requests
from typing import List, Dict, Optional, Any
from dataclasses import dataclass
from datetime import datetime

logger = logging.getLogger(__name__)


@dataclass
class PDBSearchResult:
    """Represents a single PDB search result."""
    pdb_id: str
    title: str
    experimental_method: str
    resolution: Optional[float]
    release_date: str
    organism: Optional[str]
    chains: List[str]
    has_ligands: bool
    ligand_info: List[Dict[str, str]]  # List of ligands with id, name, formula

    def __str__(self):
        res_str = f"{self.resolution:.2f}" if self.resolution else "N/A"
        return f"{self.pdb_id}: {self.title} ({self.experimental_method}, {res_str} Å)"


class PDBSearcher:
    """
    Search the RCSB Protein Data Bank using the REST API.

    Provides methods to search by protein/gene name, filter results by
    experimental method, resolution, organism, and presence of specific
    ligands or metals.
    """

    SEARCH_API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
    DATA_API_URL = "https://data.rcsb.org/rest/v1/core/entry"

    # Common organism name mappings
    ORGANISM_MAP = {
        "human": "Homo sapiens",
        "mouse": "Mus musculus",
        "rat": "Rattus norvegicus",
        "yeast": "Saccharomyces cerevisiae",
        "e.coli": "Escherichia coli",
        "ecoli": "Escherichia coli",
        "fly": "Drosophila melanogaster",
        "worm": "Caenorhabditis elegans",
        "zebrafish": "Danio rerio",
        "arabidopsis": "Arabidopsis thaliana",
    }

    # Common metal ions relevant for redox sites
    METAL_IONS = {
        "zinc": "ZN",
        "zn": "ZN",
        "iron": "FE",
        "fe": "FE",
        "copper": "CU",
        "cu": "CU",
        "manganese": "MN",
        "mn": "MN",
        "nickel": "NI",
        "ni": "NI",
        "cobalt": "CO",
        "co": "CO",
        "magnesium": "MG",
        "mg": "MG",
        "calcium": "CA",
        "ca": "CA",
    }

    def __init__(self, timeout: int = 30):
        """
        Initialize PDB searcher.

        Args:
            timeout: HTTP request timeout in seconds
        """
        self.timeout = timeout

    def normalize_organism(self, organism: str) -> str:
        """
        Normalize organism name to scientific name.

        Args:
            organism: Common or scientific name

        Returns:
            Scientific name if mapping exists, otherwise original name
        """
        if not organism:
            return organism

        organism_lower = organism.lower().strip()
        return self.ORGANISM_MAP.get(organism_lower, organism)

    def search_by_name(
        self,
        protein_name: str,
        organism: Optional[str] = None,
        max_results: int = 100
    ) -> tuple[List[str], int]:
        """
        Search PDB by protein/gene name.

        Args:
            protein_name: Protein or gene name to search for
            organism: Optional organism filter (common or scientific name)
            max_results: Maximum number of results to return

        Returns:
            Tuple of (list of PDB IDs, total count of all matching structures)
        """
        logger.debug(f"Searching PDB for '{protein_name}'" +
                    (f" in {organism}" if organism else ""))

        # Normalize organism name
        if organism:
            organism = self.normalize_organism(organism)

        # Build query - search in structure title
        # Note: PDB structures may use either gene names (e.g., "SOD1") or
        # full protein names (e.g., "superoxide dismutase") in titles.
        # Users may need to try both types of search terms.
        query_node = {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "struct.title",
                "operator": "contains_words",
                "value": protein_name
            }
        }

        # Add organism filter if provided
        if organism:
            organism_node = {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.scientific_name",
                    "operator": "exact_match",
                    "value": organism
                }
            }
            query_node = {
                "type": "group",
                "logical_operator": "and",
                "nodes": [query_node, organism_node]
            }

        # Construct full query
        query = {
            "query": query_node,
            "return_type": "entry",
            "request_options": {
                "results_content_type": ["experimental"],
                "return_all_hits": False,
                "scoring_strategy": "combined"
                # Sort by relevance score (most relevant first)
                # No explicit sort means results ordered by relevance
            }
        }

        if max_results:
            query["request_options"]["paginate"] = {
                "start": 0,
                "rows": max_results
            }

        try:
            response = requests.post(
                self.SEARCH_API_URL,
                json=query,
                timeout=self.timeout,
                headers={"Content-Type": "application/json"}
            )
            response.raise_for_status()

            results = response.json()
            logger.debug(f"API response type: {type(results)}, keys: {list(results.keys()) if isinstance(results, dict) else 'N/A'}")

            # Check if results is a dict with result_set
            if isinstance(results, dict) and "result_set" in results:
                result_set = results["result_set"]
                total_count = results.get("total_count", len(result_set))
                logger.debug(f"result_set type: {type(result_set)}, length: {len(result_set) if isinstance(result_set, list) else 'N/A'}")
                if result_set:
                    logger.debug(f"First element type: {type(result_set[0])}, value: {result_set[0]}")
                pdb_ids = [hit["identifier"] for hit in result_set]
            else:
                logger.error(f"Unexpected API response format: {type(results)}, content: {str(results)[:200]}")
                return [], 0

            logger.debug(f"Found {len(pdb_ids)} of {total_count} total structures")
            return pdb_ids, total_count

        except requests.exceptions.RequestException as e:
            logger.error(f"PDB search failed: {e}")
            return [], 0
        except (KeyError, ValueError, TypeError) as e:
            import traceback
            logger.error(f"Error parsing PDB search results: {e}, response type: {type(results) if 'results' in locals() else 'undefined'}")
            logger.debug(f"Response content: {str(results)[:500] if 'results' in locals() else 'N/A'}")
            logger.debug(f"Full traceback: {traceback.format_exc()}")
            return [], 0

    def apply_filters(
        self,
        pdb_ids: List[str],
        max_resolution: Optional[float] = None,
        methods: Optional[List[str]] = None,
        has_metal: Optional[str] = None,
        min_date: Optional[str] = None
    ) -> List[str]:
        """
        Filter PDB IDs by additional criteria.

        Args:
            pdb_ids: List of PDB IDs to filter
            max_resolution: Maximum resolution in Angstroms (X-ray/Cryo-EM only)
            methods: List of experimental methods (e.g., ["X-RAY DIFFRACTION", "SOLUTION NMR"])
            has_metal: Metal ion that must be present (e.g., "ZN", "FE")
            min_date: Minimum release date in YYYY-MM-DD format

        Returns:
            Filtered list of PDB IDs
        """
        if not pdb_ids:
            return []

        # Build filter query nodes
        filter_nodes = []

        # Resolution filter
        if max_resolution is not None:
            filter_nodes.append({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal",
                    "value": max_resolution
                }
            })

        # Experimental method filter
        if methods:
            method_nodes = []
            for method in methods:
                method_nodes.append({
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": method
                    }
                })

            if len(method_nodes) == 1:
                filter_nodes.append(method_nodes[0])
            else:
                filter_nodes.append({
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": method_nodes
                })

        # Metal ligand filter
        if has_metal:
            metal_code = self.METAL_IONS.get(has_metal.lower(), has_metal.upper())
            filter_nodes.append({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_nonpolymer_entity_container_identifiers.chem_comp_monomers",
                    "operator": "exact_match",
                    "value": metal_code
                }
            })

        # Date filter
        if min_date:
            filter_nodes.append({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_accession_info.initial_release_date",
                    "operator": "greater_or_equal",
                    "value": min_date
                }
            })

        # PDB ID filter (must be in original list)
        id_nodes = [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_container_identifiers.entry_id",
                    "operator": "exact_match",
                    "value": pdb_id
                }
            }
            for pdb_id in pdb_ids
        ]

        id_filter = {
            "type": "group",
            "logical_operator": "or",
            "nodes": id_nodes
        }

        # Combine all filters
        if filter_nodes:
            filter_nodes.append(id_filter)
            query_node = {
                "type": "group",
                "logical_operator": "and",
                "nodes": filter_nodes
            }
        else:
            query_node = id_filter

        query = {
            "query": query_node,
            "return_type": "entry",
            "request_options": {
                "return_all_hits": True
            }
        }

        try:
            response = requests.post(
                self.SEARCH_API_URL,
                json=query,
                timeout=self.timeout,
                headers={"Content-Type": "application/json"}
            )
            response.raise_for_status()

            results = response.json()
            filtered_ids = [hit["identifier"] for hit in results.get("result_set", [])]

            logger.debug(f"Filtered from {len(pdb_ids)} to {len(filtered_ids)} structures")
            return filtered_ids

        except requests.exceptions.RequestException as e:
            logger.error(f"PDB filter failed: {e}")
            return pdb_ids  # Return original list if filtering fails

    def get_structure_details(self, pdb_ids: List[str]) -> List[PDBSearchResult]:
        """
        Fetch detailed information for a list of PDB IDs.

        Args:
            pdb_ids: List of PDB IDs

        Returns:
            List of PDBSearchResult objects with detailed information
        """
        results = []

        for pdb_id in pdb_ids:
            try:
                # Fetch entry data
                url = f"{self.DATA_API_URL}/{pdb_id}"
                response = requests.get(url, timeout=self.timeout)
                response.raise_for_status()
                data = response.json()

                # Extract basic info
                struct_info = data.get("struct", {})
                if isinstance(struct_info, dict):
                    title = struct_info.get("title", "No title")
                else:
                    title = "No title"
                    logger.warning(f"Unexpected struct format for {pdb_id}: {type(struct_info)}")

                # Experimental method
                exptl_list = data.get("exptl", [])
                if isinstance(exptl_list, list) and len(exptl_list) > 0:
                    exptl = exptl_list[0]
                    if isinstance(exptl, dict):
                        method = exptl.get("method", "Unknown")
                    else:
                        method = "Unknown"
                else:
                    method = "Unknown"

                # Resolution (may not exist for NMR)
                resolution = None
                entry_info = data.get("rcsb_entry_info", {})
                if isinstance(entry_info, dict):
                    res_value = entry_info.get("resolution_combined")
                    # resolution_combined is a list with one element
                    if isinstance(res_value, list) and len(res_value) > 0:
                        resolution = res_value[0]
                    elif isinstance(res_value, (int, float)):
                        resolution = float(res_value)

                # Release date
                release_date = "Unknown"
                accession_info = data.get("rcsb_accession_info", {})
                if isinstance(accession_info, dict):
                    release_date = accession_info.get("initial_release_date", "Unknown")

                # Organism (from first entity)
                organism = None
                organisms_list = data.get("rcsb_entity_source_organism", [])
                if isinstance(organisms_list, list) and len(organisms_list) > 0:
                    org_data = organisms_list[0]
                    if isinstance(org_data, dict):
                        organism = org_data.get("scientific_name")

                # Chain IDs
                chains = []
                struct_asym = data.get("struct_asym", [])
                if isinstance(struct_asym, list):
                    chains = [asym.get("id", "") for asym in struct_asym if isinstance(asym, dict)]

                # Ligand information
                ligand_info = []
                has_ligands = False
                nonpoly_entities = data.get("nonpolymer_entities", [])
                if isinstance(nonpoly_entities, list):
                    for entity in nonpoly_entities:
                        if not isinstance(entity, dict):
                            continue
                        chem_comp = entity.get("rcsb_nonpolymer_entity_container_identifiers", {})
                        if isinstance(chem_comp, dict):
                            comp_id = chem_comp.get("comp_id")
                            if comp_id:
                                has_ligands = True
                                nonpoly_entity = entity.get("rcsb_nonpolymer_entity", {})
                                if isinstance(nonpoly_entity, dict):
                                    ligand_info.append({
                                        "id": comp_id,
                                        "name": nonpoly_entity.get("pdbx_description", comp_id),
                                        "formula": str(nonpoly_entity.get("formula_weight", ""))
                                    })

                result = PDBSearchResult(
                    pdb_id=pdb_id.upper(),
                    title=title,
                    experimental_method=method,
                    resolution=resolution,
                    release_date=release_date,
                    organism=organism,
                    chains=chains,
                    has_ligands=has_ligands,
                    ligand_info=ligand_info
                )

                results.append(result)

            except Exception as e:
                logger.warning(f"Failed to fetch details for {pdb_id}: {e}")
                # Add minimal result
                results.append(PDBSearchResult(
                    pdb_id=pdb_id.upper(),
                    title="[Details unavailable]",
                    experimental_method="Unknown",
                    resolution=None,
                    release_date="Unknown",
                    organism=None,
                    chains=[],
                    has_ligands=False,
                    ligand_info=[]
                ))

        return results

    def search_and_filter(
        self,
        protein_name: str,
        organism: Optional[str] = None,
        max_resolution: Optional[float] = None,
        methods: Optional[List[str]] = None,
        has_metal: Optional[str] = None,
        max_results: int = 100
    ) -> tuple[List[PDBSearchResult], int]:
        """
        One-step search and filter operation.

        Args:
            protein_name: Protein or gene name to search for
            organism: Optional organism filter
            max_resolution: Maximum resolution in Angstroms
            methods: List of experimental methods to include
            has_metal: Metal ion that must be present
            max_results: Maximum number of results

        Returns:
            Tuple of (list of PDBSearchResult objects sorted by resolution, total count)
        """
        # Initial search
        pdb_ids, total_count = self.search_by_name(protein_name, organism, max_results=max_results * 2)

        if not pdb_ids:
            return [], 0

        # Apply filters if any
        if max_resolution or methods or has_metal:
            pdb_ids = self.apply_filters(
                pdb_ids,
                max_resolution=max_resolution,
                methods=methods,
                has_metal=has_metal
            )

        # Limit to max_results
        pdb_ids = pdb_ids[:max_results]

        # Fetch detailed information
        results = self.get_structure_details(pdb_ids)

        # Sort by resolution (best first), putting NMR/no-resolution at end
        results.sort(key=lambda r: (r.resolution is None, r.resolution or 999.0))

        return results, total_count
