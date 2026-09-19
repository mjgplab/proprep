"""
PDB Search Module - RCSB PDB Search API Integration

This module provides search functionality for the RCSB Protein Data Bank,
allowing users to find structures by protein name, gene name, organism, the
ligand they contain, and other criteria without needing to know the exact PDB ID.

Author: ProPrep Development Team
"""

import logging
import re
import requests
from typing import List, Dict, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime

logger = logging.getLogger(__name__)


# The elements that are NOT metals: hydrogen and its isotopes as the PDB writes them, the
# nonmetals, the metalloids, the halogens and the noble gases, plus X, the PDB's symbol for an
# atom of unknown element. Every other element is a metal. The set is written this way round
# because it is closed: no new nonmetal will be found, whereas lists of metals go stale (the
# PDB holds americium, curium and californium ligands) and lists of metal-containing ligand
# codes are out of date every week. Metalloids count as not metal: arsenic is in the
# cacodylate buffer of thousands of entries, boron in hundreds of inhibitors.
NONMETAL_ELEMENTS = frozenset({
    "H", "D", "T", "He",
    "B", "C", "N", "O", "F", "Ne",
    "Si", "P", "S", "Cl", "Ar",
    "Ge", "As", "Se", "Br", "Kr",
    "Sb", "Te", "I", "Xe",
    "At", "Rn",
    "X",
})


def formula_elements(formula: str) -> List[str]:
    """Element symbols in a PDB formula such as "C34 H32 Fe N4 O4", as written, in order."""
    symbols = []
    for token in (formula or "").split():
        match = re.fullmatch(r"([A-Za-z]{1,2})\d*", token)
        if match:
            symbols.append(match.group(1))
    return symbols


def metals_in_formula(formula: str) -> List[str]:
    """The metal elements in a PDB formula, as title-case symbols."""
    return [symbol.title() for symbol in formula_elements(formula)
            if symbol.title() not in NONMETAL_ELEMENTS]


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
    matched_ligands: List[str] = field(default_factory=list)  # searched-for ligand codes present
    ligand_rscc: Optional[float] = None  # best real-space correlation among matched ligand copies

    def metal_components(self, element: Optional[str] = None) -> List[str]:
        """The non-polymer components that contain a metal, read from their formulas.

        Each is its PDB code, which for an ion carries the deposited oxidation state
        (FE2, FE), followed by the metal in brackets when the component is more than
        the ion itself: "HEM (Fe)", "SF4 (Fe)", "ZN". These are the metal-containing
        components, not the groups that coordinate the metal.

        Args:
            element: a metal's symbol, to keep only the components containing that metal
        """
        wanted = element.title() if element else None
        components = []
        for ligand in self.ligand_info:
            formula = ligand.get("formula", "")
            metals = metals_in_formula(formula)
            if not ((wanted in metals) if wanted else metals):
                continue
            code = ligand.get("id", "")
            bare_ion = len(formula_elements(formula)) == 1
            components.append(code if bare_ion else f"{code} ({', '.join(metals)})")
        return sorted(components)

    def metals(self) -> List[str]:
        """The metal elements present in this entry's ligands."""
        found: List[str] = []
        for ligand in self.ligand_info:
            for metal in metals_in_formula(ligand.get("formula", "")):
                if metal not in found:
                    found.append(metal)
        return found

    def __str__(self):
        res_str = f"{self.resolution:.2f}" if self.resolution else "N/A"
        return f"{self.pdb_id}: {self.title} ({self.experimental_method}, {res_str} Å)"


@dataclass
class LigandDefinition:
    """A ligand (chemical component) defined in the PDB, identified by its code."""
    comp_id: str
    name: str
    formula: str
    smiles: Optional[str] = None
    relation: str = ""  # how the search arrived at this ligand
    similarity: Optional[float] = None  # RCSB fingerprint similarity to the ligand searched from


class PDBSearchError(Exception):
    """The RCSB service could not be reached, or it rejected the request."""


class PDBSearcher:
    """
    Search the RCSB Protein Data Bank using the REST API.

    Provides methods to search by protein/gene name, filter results by
    experimental method, resolution and organism, and to find the structures
    containing given ligands.
    """

    SEARCH_API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
    GRAPHQL_API_URL = "https://data.rcsb.org/graphql"

    # Entries per GraphQL details request
    DETAILS_BATCH_SIZE = 100

    # The most hits the Search API returns for one request
    MAX_ROWS = 10000

    # Search attribute holding the code of each ligand in an entry
    LIGAND_ATTRIBUTE = "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id"

    # Ligands listed for a name lookup; the number found is always reported
    NAME_MATCH_ROWS = 50

    # A related-ligand search compares the formula of every ligand containing
    # the same non-hydrogen atoms. Past this many (single atoms, small
    # fragments) the comparison is skipped, and the caller is told.
    RELATED_FORMULA_MAX = 500

    _LIGAND_QUERY = """
    query($ids: [String!]!) {
      chem_comps(comp_ids: $ids) {
        chem_comp { id name formula }
        rcsb_chem_comp_descriptor { SMILES_stereo SMILES }
      }
    }
    """

    _ENTRY_DETAILS_QUERY = """
    query($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        struct { title }
        exptl { method }
        rcsb_entry_info { resolution_combined }
        rcsb_accession_info { initial_release_date }
        polymer_entities {
          rcsb_entity_source_organism { scientific_name }
          rcsb_polymer_entity_container_identifiers { auth_asym_ids }
        }
        nonpolymer_entities {
          nonpolymer_comp { chem_comp { id name formula } }
          nonpolymer_entity_instances {
            rcsb_nonpolymer_instance_validation_score { RSCC }
          }
        }
      }
    }
    """

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

    def count_by_name(self, protein_name: str, organism: Optional[str] = None) -> int:
        """
        Count the structures a title search would find, without fetching them.

        Raises:
            PDBSearchError: the service could not be reached
        """
        nodes = [{
            "type": "terminal",
            "service": "text",
            "parameters": {"attribute": "struct.title", "operator": "contains_words", "value": protein_name}
        }]
        if organism:
            nodes.append({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.scientific_name",
                    "operator": "exact_match",
                    "value": self.normalize_organism(organism)
                }
            })
        return self._post_search({
            "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
            "return_type": "entry",
            "request_options": {"results_content_type": ["experimental"], "return_counts": True}
        }).get("total_count", 0)

    def apply_filters(
        self,
        pdb_ids: List[str],
        max_resolution: Optional[float] = None,
        methods: Optional[List[str]] = None,
        min_date: Optional[str] = None
    ) -> List[str]:
        """
        Filter PDB IDs by additional criteria.

        Args:
            pdb_ids: List of PDB IDs to filter
            max_resolution: Maximum resolution in Angstroms (X-ray/Cryo-EM only)
            methods: List of experimental methods (e.g., ["X-RAY DIFFRACTION", "SOLUTION NMR"])
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

    def get_structure_details(
        self,
        pdb_ids: List[str],
        ligand_ids: Optional[List[str]] = None
    ) -> List[PDBSearchResult]:
        """
        Fetch detailed information for a list of PDB IDs.

        The entry record of the REST Data API carries no organism, chain or
        ligand information (those live on the polymer and nonpolymer entity
        records), so the details come from the GraphQL Data API, which returns
        an entry together with its entities, many entries per request.

        Args:
            pdb_ids: List of PDB IDs
            ligand_ids: Ligand codes of interest. When given, each result lists
                which of them the entry contains (``matched_ligands``) and the
                best real-space correlation among their copies (``ligand_rscc``).

        Returns:
            List of PDBSearchResult objects, in the order of ``pdb_ids``
        """
        wanted = {code.upper() for code in (ligand_ids or [])}
        found: Dict[str, PDBSearchResult] = {}

        for start in range(0, len(pdb_ids), self.DETAILS_BATCH_SIZE):
            batch = pdb_ids[start:start + self.DETAILS_BATCH_SIZE]
            try:
                response = requests.post(
                    self.GRAPHQL_API_URL,
                    json={"query": self._ENTRY_DETAILS_QUERY, "variables": {"ids": batch}},
                    timeout=self.timeout
                )
                response.raise_for_status()
                entries = (response.json().get("data") or {}).get("entries") or []
            except (requests.exceptions.RequestException, ValueError) as e:
                logger.warning(f"Failed to fetch details for {len(batch)} entries: {e}")
                continue

            for entry in entries:
                if isinstance(entry, dict) and entry.get("rcsb_id"):
                    found[entry["rcsb_id"].upper()] = self._result_from_entry(entry, wanted)

        results = []
        for pdb_id in pdb_ids:
            result = found.get(pdb_id.upper())
            if result is None:
                logger.warning(f"No details returned for {pdb_id}")
                result = PDBSearchResult(
                    pdb_id=pdb_id.upper(),
                    title="[Details unavailable]",
                    experimental_method="Unknown",
                    resolution=None,
                    release_date="Unknown",
                    organism=None,
                    chains=[],
                    has_ligands=False,
                    ligand_info=[]
                )
            results.append(result)

        return results

    @staticmethod
    def _result_from_entry(entry: Dict[str, Any], wanted: set) -> PDBSearchResult:
        """Build a PDBSearchResult from one GraphQL ``entries`` element."""
        title = (entry.get("struct") or {}).get("title") or "No title"

        exptl = entry.get("exptl") or []
        method = (exptl[0] or {}).get("method") or "Unknown" if exptl else "Unknown"

        # resolution_combined is a list (absent for NMR)
        res_list = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
        resolution = res_list[0] if res_list else None

        release_date = (entry.get("rcsb_accession_info") or {}).get("initial_release_date") or "Unknown"

        organism = None
        chains: List[str] = []
        for entity in entry.get("polymer_entities") or []:
            for source in entity.get("rcsb_entity_source_organism") or []:
                if organism is None and source.get("scientific_name"):
                    organism = source["scientific_name"]
            ids = entity.get("rcsb_polymer_entity_container_identifiers") or {}
            for chain in ids.get("auth_asym_ids") or []:
                if chain not in chains:
                    chains.append(chain)

        ligand_info = []
        matched: List[str] = []
        rscc_values: List[float] = []
        for entity in entry.get("nonpolymer_entities") or []:
            comp = ((entity.get("nonpolymer_comp") or {}).get("chem_comp")) or {}
            comp_id = comp.get("id")
            if not comp_id:
                continue
            ligand_info.append({
                "id": comp_id,
                "name": comp.get("name") or comp_id,
                "formula": comp.get("formula") or ""
            })
            if comp_id.upper() in wanted:
                if comp_id not in matched:
                    matched.append(comp_id)
                for instance in entity.get("nonpolymer_entity_instances") or []:
                    for score in instance.get("rcsb_nonpolymer_instance_validation_score") or []:
                        if score.get("RSCC") is not None:
                            rscc_values.append(score["RSCC"])

        return PDBSearchResult(
            pdb_id=entry["rcsb_id"].upper(),
            title=title,
            experimental_method=method,
            resolution=resolution,
            release_date=release_date,
            organism=organism,
            chains=chains,
            has_ligands=bool(ligand_info),
            ligand_info=ligand_info,
            matched_ligands=matched,
            ligand_rscc=max(rscc_values) if rscc_values else None
        )

    def search_and_filter(
        self,
        protein_name: str,
        organism: Optional[str] = None,
        max_resolution: Optional[float] = None,
        methods: Optional[List[str]] = None,
        max_results: int = 100
    ) -> tuple[List[PDBSearchResult], int]:
        """
        One-step search and filter operation.

        Args:
            protein_name: Protein or gene name to search for
            organism: Optional organism filter
            max_resolution: Maximum resolution in Angstroms
            methods: List of experimental methods to include
            max_results: Maximum number of results

        Returns:
            Tuple of (list of PDBSearchResult objects sorted by resolution, total count)
        """
        # Initial search. Filters remove hits, so fetch more only when there are filters
        filtering = bool(max_resolution or methods)
        fetch = min(max_results * 2, self.MAX_ROWS) if filtering else max_results
        pdb_ids, total_count = self.search_by_name(protein_name, organism, max_results=fetch)

        if not pdb_ids:
            return [], 0

        # Apply filters if any
        if max_resolution or methods:
            pdb_ids = self.apply_filters(
                pdb_ids,
                max_resolution=max_resolution,
                methods=methods
            )

        # Limit to max_results
        pdb_ids = pdb_ids[:max_results]

        # Fetch detailed information
        results = self.get_structure_details(pdb_ids)

        # Sort by resolution (best first), putting NMR/no-resolution at end
        results.sort(key=lambda r: (r.resolution is None, r.resolution or 999.0))

        return results, total_count

    # ------------------------------------------------------------------
    # Search by ligand
    # ------------------------------------------------------------------

    def _post_search(self, query: Dict[str, Any]) -> Dict[str, Any]:
        """POST a Search API query; no hits (HTTP 204) gives an empty dict.

        Raises:
            PDBSearchError: the service could not be reached or rejected the query
        """
        try:
            response = requests.post(self.SEARCH_API_URL, json=query, timeout=self.timeout)
            if response.status_code == 204:
                return {}
            response.raise_for_status()
            return response.json()
        except (requests.exceptions.RequestException, ValueError) as e:
            raise PDBSearchError(str(e)) from e

    def _ligand_node(self, comp_ids: List[str]) -> Dict[str, Any]:
        """Query node matching entries that contain any of the ligand codes."""
        return {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": self.LIGAND_ATTRIBUTE,
                "operator": "in",
                "value": [code.upper() for code in comp_ids]
            }
        }

    def get_ligand_definitions(self, comp_ids: List[str]) -> List[LigandDefinition]:
        """
        Look ligand codes up in the PDB chemical component dictionary.

        Returns:
            Definitions in the order asked for; codes the PDB does not define are absent.

        Raises:
            PDBSearchError: the service could not be reached
        """
        codes = list(dict.fromkeys(code.upper() for code in comp_ids if code))
        found: Dict[str, LigandDefinition] = {}
        for start in range(0, len(codes), self.DETAILS_BATCH_SIZE):
            batch = codes[start:start + self.DETAILS_BATCH_SIZE]
            try:
                response = requests.post(
                    self.GRAPHQL_API_URL,
                    json={"query": self._LIGAND_QUERY, "variables": {"ids": batch}},
                    timeout=self.timeout
                )
                response.raise_for_status()
                comps = (response.json().get("data") or {}).get("chem_comps") or []
            except (requests.exceptions.RequestException, ValueError) as e:
                raise PDBSearchError(str(e)) from e
            for comp in comps:
                chem = (comp or {}).get("chem_comp") or {}
                if not chem.get("id"):
                    continue
                descriptor = comp.get("rcsb_chem_comp_descriptor") or {}
                found[chem["id"].upper()] = LigandDefinition(
                    comp_id=chem["id"].upper(),
                    name=chem.get("name") or "",
                    formula=chem.get("formula") or "",
                    smiles=descriptor.get("SMILES_stereo") or descriptor.get("SMILES")
                )
        return [found[code] for code in codes if code in found]

    def lookup_ligands(self, text: str) -> tuple[List[LigandDefinition], int]:
        """
        Find ligands from what the user typed: a ligand code, or words of a name.

        The text is tried both ways. Names and synonyms match on whole words
        ("biopterin" does not match "tetrahydrobiopterin").

        Returns:
            Tuple of (ligands, number of name matches in the PDB). At most
            NAME_MATCH_ROWS name matches are listed.

        Raises:
            PDBSearchError: the service could not be reached
        """
        text = text.strip()
        by_code: List[LigandDefinition] = []
        if re.fullmatch(r"[A-Za-z0-9]{1,5}", text):
            by_code = self.get_ligand_definitions([text])
            for ligand in by_code:
                ligand.relation = "code"

        name_nodes = [
            {
                "type": "terminal",
                "service": "text_chem",
                "parameters": {"attribute": attribute, "operator": "contains_phrase", "value": text}
            }
            for attribute in ("chem_comp.name", "rcsb_chem_comp_synonyms.name")
        ]
        hits = self._post_search({
            "query": {"type": "group", "logical_operator": "or", "nodes": name_nodes},
            "return_type": "mol_definition",
            "request_options": {"paginate": {"start": 0, "rows": self.NAME_MATCH_ROWS}}
        })
        name_total = hits.get("total_count", 0)
        # mol_definition also returns BIRD molecules (PRD_...), which are not ligand codes
        name_ids = [hit["identifier"] for hit in hits.get("result_set", [])
                    if not hit["identifier"].startswith("PRD_")]
        known = {ligand.comp_id for ligand in by_code}
        by_name = self.get_ligand_definitions([i for i in name_ids if i.upper() not in known])
        for ligand in by_name:
            ligand.relation = "name"

        return by_code + by_name, name_total

    @staticmethod
    def split_formula(formula: str) -> tuple[Dict[str, int], int]:
        """Split a formula such as "C9 H15 N5 O3" into ({non-hydrogen atoms}, hydrogen count)."""
        heavy: Dict[str, int] = {}
        hydrogens = 0
        for token in (formula or "").split():
            match = re.fullmatch(r"([A-Za-z]{1,2})(\d*)", token)
            if not match:
                continue
            element, count = match.group(1), int(match.group(2) or 1)
            if element.upper() in ("H", "D"):
                hydrogens += count
            else:
                heavy[element] = heavy.get(element, 0) + count
        return heavy, hydrogens

    def find_related_ligands(self, seed: LigandDefinition) -> tuple[List[LigandDefinition], List[str]]:
        """
        Find ligands related to one ligand, in two ways.

        1. The same non-hydrogen atoms with any number of hydrogens: other
           redox or protonation states, tautomers and stereoisomers.
        2. RCSB fingerprint similarity to the ligand's SMILES: analogs. The
           score, and which ligands are similar enough to be returned, are RCSB's.

        Nothing is excluded here; a ligand that only resembles the seed is
        labelled as such and the choice is left to the user.

        Returns:
            Tuple of (related ligands without the seed, notes on anything skipped)

        Raises:
            PDBSearchError: the service could not be reached
        """
        notes: List[str] = []
        seed_heavy, seed_hydrogens = self.split_formula(seed.formula)

        same_atoms_ids: List[str] = []
        if seed_heavy:
            formula_node = {
                "type": "terminal",
                "service": "chemical",
                "parameters": {
                    "type": "formula",
                    # Element symbols exactly as the PDB writes them: the match is case-sensitive
                    "value": " ".join(f"{el}{n}" for el, n in seed_heavy.items()),
                    "match_subset": True
                }
            }
            count = self._post_search({
                "query": formula_node,
                "return_type": "mol_definition",
                "request_options": {"return_counts": True}
            }).get("total_count", 0)
            if count > self.RELATED_FORMULA_MAX:
                notes.append(
                    f"{count} ligands contain the atoms of {seed.comp_id}; the search for the same "
                    f"atoms with a different number of hydrogens is skipped above {self.RELATED_FORMULA_MAX}."
                )
            elif count:
                hits = self._post_search({
                    "query": formula_node,
                    "return_type": "mol_definition",
                    "request_options": {"return_all_hits": True}
                })
                same_atoms_ids = [hit["identifier"] for hit in hits.get("result_set", [])]

        scores: Dict[str, float] = {}
        if seed.smiles:
            hits = self._post_search({
                "query": {
                    "type": "terminal",
                    "service": "chemical",
                    "parameters": {
                        "type": "descriptor",
                        "descriptor_type": "SMILES",
                        "value": seed.smiles,
                        "match_type": "fingerprint-similarity"
                    }
                },
                "return_type": "mol_definition",
                "request_options": {"return_all_hits": True}
            })
            scores = {hit["identifier"].upper(): hit.get("score") for hit in hits.get("result_set", [])}
        else:
            notes.append(f"{seed.comp_id} has no SMILES in the PDB, so no similarity search was made.")

        candidates = [i for i in dict.fromkeys(same_atoms_ids + list(scores))
                      if i.upper() != seed.comp_id and not i.startswith("PRD_")]
        related = self.get_ligand_definitions(candidates)

        for ligand in related:
            heavy, hydrogens = self.split_formula(ligand.formula)
            ligand.similarity = scores.get(ligand.comp_id)
            if heavy == seed_heavy:
                difference = hydrogens - seed_hydrogens
                ligand.relation = f"same atoms, {difference:+d} H" if difference else "same formula"
            elif ligand.comp_id in scores:
                ligand.relation = "similar structure"
        # The formula search returns supersets (extra elements); only ligands that are
        # the same atoms or that RCSB scored as similar are related
        related = [ligand for ligand in related if ligand.relation]
        related.sort(key=lambda l: (l.relation == "similar structure", -(l.similarity or 0.0), l.comp_id))
        return related, notes

    def count_entries_by_ligand(self, comp_ids: List[str]) -> tuple[int, Dict[str, int]]:
        """
        Count the experimental entries containing each ligand code.

        Returns:
            Tuple of (entries containing any of the codes, {code: entries containing it})

        Raises:
            PDBSearchError: the service could not be reached
        """
        codes = [code.upper() for code in comp_ids]
        if not codes:
            return 0, {}
        hits = self._post_search({
            "query": self._ligand_node(codes),
            "return_type": "entry",
            "request_options": {
                "results_content_type": ["experimental"],
                "paginate": {"start": 0, "rows": 1},
                "facets": [{
                    "name": "by_ligand",
                    "aggregation_type": "terms",
                    "attribute": self.LIGAND_ATTRIBUTE,
                    "min_interval_population": 1,
                    "max_num_intervals": 65000
                }]
            }
        })
        per_code = {code: 0 for code in codes}
        for facet in hits.get("facets", []):
            for bucket in facet.get("buckets", []):
                if bucket.get("label") in per_code:
                    per_code[bucket["label"]] = bucket.get("population", 0)
        return hits.get("total_count", 0), per_code

    def search_by_ligands(
        self,
        comp_ids: List[str],
        max_results: int
    ) -> tuple[List[PDBSearchResult], int]:
        """
        Find the experimental entries containing any of the ligand codes.

        Entries come best resolution first (entries with no resolution last),
        so ``max_results`` keeps the best-resolved ones.

        Returns:
            Tuple of (results with ``matched_ligands`` and ``ligand_rscc`` filled, total count)

        Raises:
            PDBSearchError: the service could not be reached
        """
        hits = self._post_search({
            "query": self._ligand_node(comp_ids),
            "return_type": "entry",
            "request_options": {
                "results_content_type": ["experimental"],
                "paginate": {"start": 0, "rows": max_results},
                "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}]
            }
        })
        pdb_ids = [hit["identifier"] for hit in hits.get("result_set", [])]
        return self.get_structure_details(pdb_ids, ligand_ids=comp_ids), hits.get("total_count", 0)
