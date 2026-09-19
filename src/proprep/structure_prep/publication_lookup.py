"""
Publication Lookup - Europe PMC

Looks up the article a structure was published in: its abstract, and whether a
free full text exists and where. The citation itself (authors, title, journal,
PubMed ID, DOI) is in the PDB file's JRNL records and needs no network.

ProPrep does not download articles. Most are behind a subscription, and the
PDF services of PubMed Central and of the publishers refuse scripted downloads;
the links given here open in a browser, where a reader's own access applies.
"""

import html
import logging
import re
from dataclasses import dataclass, field
from typing import List, Optional

import requests

logger = logging.getLogger(__name__)

EUROPE_PMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


class PublicationLookupError(Exception):
    """Europe PMC could not be reached, or it rejected the request."""


@dataclass
class Publication:
    """What Europe PMC holds on one article."""
    title: str
    abstract: str
    pmid: Optional[str]
    pmcid: Optional[str]
    doi: Optional[str]
    open_access: bool  # licensed for reuse, not merely free to read
    free_full_text_urls: List[str] = field(default_factory=list)


def _plain_text(marked_up: str) -> str:
    """Europe PMC text without its markup (<h4>, <i>, <sub>, entities)."""
    text = re.sub(r"</h4>", ": ", marked_up or "")
    text = re.sub(r"<[^>]+>", "", text)
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def pubmed_url(pmid: str) -> str:
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"


def doi_url(doi: str) -> str:
    return f"https://doi.org/{doi}"


def lookup_publication(
    pmid: Optional[str] = None,
    doi: Optional[str] = None,
    timeout: int = 15
) -> Optional[Publication]:
    """
    Find an article in Europe PMC by PubMed ID, or by DOI when there is no PubMed ID.

    Returns:
        The publication, or None if Europe PMC has no such article

    Raises:
        PublicationLookupError: Europe PMC could not be reached
    """
    if pmid:
        query = f"EXT_ID:{pmid} AND SRC:MED"
    elif doi:
        query = f'DOI:"{doi}"'
    else:
        return None

    try:
        response = requests.get(
            EUROPE_PMC_SEARCH_URL,
            params={"query": query, "format": "json", "resultType": "core", "pageSize": 1},
            timeout=timeout
        )
        response.raise_for_status()
        hits = (response.json().get("resultList") or {}).get("result") or []
    except (requests.exceptions.RequestException, ValueError) as e:
        raise PublicationLookupError(str(e)) from e

    if not hits:
        return None
    hit = hits[0]

    free_urls = []
    for link in (hit.get("fullTextUrlList") or {}).get("fullTextUrl") or []:
        if link.get("availability") in ("Free", "Open access") and link.get("url") not in free_urls:
            free_urls.append(link.get("url"))

    return Publication(
        title=_plain_text(hit.get("title")),
        abstract=_plain_text(hit.get("abstractText")),
        pmid=hit.get("pmid"),
        pmcid=hit.get("pmcid"),
        doi=hit.get("doi"),
        open_access=hit.get("isOpenAccess") == "Y",
        free_full_text_urls=free_urls
    )
