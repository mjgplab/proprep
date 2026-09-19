"""The article behind a structure: citation from the file's JRNL records on load; abstract and
free-full-text availability from Europe PMC in the metadata viewer. ProPrep downloads no article."""

import io
import os
import types

import pytest
import requests
from rich.console import Console

import proprep.structure_prep.pdb_loader as pdb_loader
import proprep.structure_prep.publication_lookup as pl
from proprep.structure_prep.pdb_loader import PDBMetadataExtractor
from proprep.structure_prep.publication_lookup import Publication, PublicationLookupError, lookup_publication
from proprep.structure_prep.structure_loader import StructureLoaderModule

PUBLISHED = """\
HEADER    OXIDOREDUCTASE                          02-DEC-04   2BF7
TITLE     LEISHMANIA MAJOR PTERIDINE REDUCTASE 1 IN COMPLEX WITH NADP AND BIOPTERIN
JRNL        AUTH   A.W.SCHUETTELKOPF,L.W.HARDY,S.M.BEVERLEY,W.N.HUNTER
JRNL        TITL   STRUCTURES OF LEISHMANIA MAJOR PTERIDINE REDUCTASE
JRNL        TITL 2 COMPLEXES [REVEAL] THE ACTIVE SITE FEATURES
JRNL        REF    J.MOL.BIOL.                   V. 352   105 2005
JRNL        REFN                   ISSN 0022-2836
JRNL        PMID   16055151
JRNL        DOI    10.1016/J.JMB.2005.06.076
REMARK   2 RESOLUTION.    2.40 ANGSTROMS.
"""
UNPUBLISHED = """\
HEADER    TRANSFERASE                             05-MAR-18   5ZFK
JRNL        AUTH   K.H.KILLIVALAVAN
JRNL        TITL   UDP GLUCOSE ALPHA TETRAHYDROBIOPTERIN GLYCOSYLTRANSFERASE
JRNL        REF    TO BE PUBLISHED
JRNL        REFN
"""


def _pdb(tmp_path, text):
    path = tmp_path / "x.pdb"
    path.write_text(text)
    return str(path)


def _shown_on_load(tmp_path, text):
    console = Console(width=200, record=True, file=io.StringIO(), force_terminal=False)
    m = StructureLoaderModule()
    m.processor = types.SimpleNamespace(console=console)
    m.initialize()
    m.console = console
    m._show_primary_publication(PDBMetadataExtractor(_pdb(tmp_path, text)))
    return console.export_text()


def test_load_shows_the_publication_from_the_file_without_any_network(tmp_path, monkeypatch):
    def no_network(*a, **k):
        raise AssertionError("loading a structure must not look anything up")
    monkeypatch.setattr(pl.requests, "get", no_network)

    out = _shown_on_load(tmp_path, PUBLISHED)

    assert "STRUCTURES OF LEISHMANIA MAJOR PTERIDINE REDUCTASE COMPLEXES [REVEAL] THE ACTIVE SITE FEATURES" in out
    assert "A.W.SCHUETTELKOPF" in out and "J.MOL.BIOL. V. 352 105 2005" in out
    assert "https://pubmed.ncbi.nlm.nih.gov/16055151/" in out
    assert "https://doi.org/10.1016/J.JMB.2005.06.076" in out


def test_load_says_when_the_structure_was_never_published(tmp_path):
    out = _shown_on_load(tmp_path, UNPUBLISHED)
    assert "Not published" in out and "pubmed" not in out


def test_load_shows_nothing_for_a_file_without_a_citation(tmp_path):
    assert _shown_on_load(tmp_path, "HEADER    X\nATOM      1  N   ALA A   1       0.0 0.0 0.0\n") == ""


class _Response:
    def __init__(self, payload):
        self.payload = payload

    def raise_for_status(self):
        pass

    def json(self):
        return self.payload


def test_lookup_reads_abstract_and_free_links_and_strips_markup(monkeypatch):
    asked = {}

    def get(url, params=None, timeout=None):
        asked.update(params)
        return _Response({"resultList": {"result": [{
            "title": "The <i>folE</i> riboswitch.", "pmid": "38769061", "pmcid": "PMC11317127", "doi": "10.1093/nar/gkae377",
            "abstractText": "<h4>Background</h4>Fe<sup>2+</sup> &amp; pterin.", "isOpenAccess": "Y",
            "fullTextUrlList": {"fullTextUrl": [
                {"availability": "Subscription required", "url": "https://doi.org/x"},
                {"availability": "Open access", "url": "https://europepmc.org/articles/PMC11317127"},
                {"availability": "Free", "url": "https://europepmc.org/articles/PMC11317127"}]}}]}})
    monkeypatch.setattr(pl.requests, "get", get)

    publication = lookup_publication(pmid="38769061", doi="ignored when there is a PubMed ID")

    assert asked["query"] == "EXT_ID:38769061 AND SRC:MED"
    assert publication.title == "The folE riboswitch." and publication.abstract == "Background: Fe2+ & pterin."
    assert publication.open_access and publication.free_full_text_urls == ["https://europepmc.org/articles/PMC11317127"]


def test_lookup_distinguishes_unknown_article_from_unreachable_service(monkeypatch):
    monkeypatch.setattr(pl.requests, "get", lambda *a, **k: _Response({"resultList": {"result": []}}))
    assert lookup_publication(doi="10.1/none") is None
    assert lookup_publication() is None

    def down(*a, **k):
        raise requests.exceptions.ConnectionError("no route to host")
    monkeypatch.setattr(pl.requests, "get", down)
    with pytest.raises(PublicationLookupError):
        lookup_publication(pmid="1")


def test_metadata_viewer_shows_abstract_and_availability(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(pdb_loader, "lookup_publication", lambda pmid=None, doi=None: Publication(
        title="t", abstract="Ternary complexes with cofactor and biopterin.", pmid=pmid, pmcid=None, doi=doi, open_access=False))
    PDBMetadataExtractor(_pdb(tmp_path, PUBLISHED)).display_citations()
    out = capsys.readouterr().out
    assert "Ternary complexes with cofactor and biopterin." in out
    assert "no free version known to Europe PMC" in out and "https://doi.org/10.1016/J.JMB.2005.06.076" in out


def test_metadata_viewer_offline_still_shows_the_citation(tmp_path, monkeypatch, capsys):
    def down(pmid=None, doi=None):
        raise PublicationLookupError("no route to host")
    monkeypatch.setattr(pdb_loader, "lookup_publication", down)
    PDBMetadataExtractor(_pdb(tmp_path, PUBLISHED)).display_citations()
    out = capsys.readouterr().out
    assert "PubMed ID: 16055151" in out and "Europe PMC could not be reached (no route to host)" in out


def test_metadata_viewer_makes_no_lookup_for_an_unpublished_structure(tmp_path, monkeypatch, capsys):
    def never(pmid=None, doi=None):
        raise AssertionError("nothing to look up")
    monkeypatch.setattr(pdb_loader, "lookup_publication", never)
    PDBMetadataExtractor(_pdb(tmp_path, UNPUBLISHED)).display_citations()
    assert "TO BE PUBLISHED" in capsys.readouterr().out


@pytest.mark.skipif(not os.environ.get("PROPREP_NETWORK_TESTS"), reason="set PROPREP_NETWORK_TESTS=1 to query Europe PMC")
def test_live_europe_pmc():
    paywalled = lookup_publication(pmid="16055151")             # 2BF7, J. Mol. Biol. 2005
    assert "biopterin" in paywalled.abstract and not paywalled.open_access and paywalled.free_full_text_urls == []
    open_access = lookup_publication(pmid="38769061")           # 8XZN, Nucleic Acids Res. 2024
    assert open_access.open_access and open_access.pmcid == "PMC11317127" and open_access.free_full_text_urls
    assert lookup_publication(doi="10.1093/nar/gkae377").pmid == "38769061"
