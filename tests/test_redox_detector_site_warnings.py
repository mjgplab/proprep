#!/usr/bin/env python3
"""
Redox Site Detector: the template results table must expose two facts the
count-based search used to discard.

1. How far each added residue was from the search boundary. A count search
   always returns the requested number of residues, so a mono-His heme at the
   end of a truncated multiheme chain silently picks up the next-nearest His
   (9YUQ: His 143 at 5.1 Å for heme 302, versus 2.0 Å for its real ligand).
2. Whether a residue belongs to more than one site. One side chain cannot
   coordinate two centres; this is a definition error.

Run with: pytest tests/test_redox_detector_site_warnings.py
"""

import io
import sys
from pathlib import Path

import pytest
from rich.console import Console

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep.structure_prep.comprehensive_redox_detector import (  # noqa: E402
    ComprehensiveRedoxDetector,
    RedoxSite,
    RedoxSiteAtom,
    SearchResult,
    _farthest_search_residue,
    _record_residue_search_distances,
    _shared_residue_notes,
)


def _atom(chain, resname, resid, name, x, element="N", icode=""):
    return RedoxSiteAtom(chain, resname, resid, name, (x, 0.0, 0.0), element, insertion_code=icode)


def _heme_site(number, heme_resid, his_resids, his_icode=""):
    site = RedoxSite(f"site_{number}", "test")
    site.site_type = "heme_bis_his_c_type"
    site.add_atom(_atom("N", "HEC", heme_resid, "FE", float(heme_resid), "FE"))
    for r in his_resids:
        site.add_atom(_atom("N", "HIS", r, "NE2", float(r), icode=his_icode))
    return site


def _search_result(residue_distances):
    atoms = []
    for (chain, resid, icode), dists in residue_distances.items():
        for i, d in enumerate(dists):
            atoms.append({"chain": chain, "resname": "HIS", "resid": resid, "atom_name": f"N{i}",
                          "element": "N", "coords": (0.0, 0.0, float(i)), "distance": d,
                          "insertion_code": icode})
    return SearchResult(detected_atoms=atoms, detected_residues=[], search_parameters=None,
                        boundary_coords=[], total_atoms_found=len(atoms))


def test_record_keeps_min_distance_per_selected_residue_only():
    site = _heme_site(2, 302, [])
    result = _search_result({("N", 190, " "): [4.16, 2.05], ("N", 143, " "): [7.75, 5.08],
                             ("N", 999, " "): [1.0]})
    _record_residue_search_distances(site, result, [("N", 190, " "), ("N", 143, "")])
    assert site.residue_search_distances == {("N", 190, " "): 2.05, ("N", 143, " "): 5.08}


def test_farthest_residue_is_labelled_with_its_resname():
    site = _heme_site(2, 302, [190, 143])
    site.residue_search_distances = {("N", 190, " "): 2.05, ("N", 143, " "): 5.08}
    assert _farthest_search_residue(site) == ("HIS N:143", 5.08)
    assert _farthest_search_residue(_heme_site(3, 303, [])) is None


def test_shared_residue_is_reported_on_both_sites_regardless_of_icode_spelling():
    site1 = _heme_site(1, 301, [143, 216], his_icode="")
    site2 = _heme_site(2, 302, [143, 190], his_icode=" ")
    site3 = _heme_site(3, 303, [27, 103])
    notes = _shared_residue_notes([site1, site2, site3])
    assert notes == {"site_1": ["HIS N:143 also in site 2"],
                     "site_2": ["HIS N:143 also in site 1"]}


def test_results_table_shows_distance_and_shared_residue():
    buf = io.StringIO()
    detector = ComprehensiveRedoxDetector.__new__(ComprehensiveRedoxDetector)
    detector.console = Console(file=buf, width=160, force_terminal=False)
    detector.processor = None
    site1 = _heme_site(1, 301, [143, 216])
    site1.residue_search_distances = {("N", 143, " "): 2.03, ("N", 216, " "): 1.99}
    detector.final_sites = [site1]
    site2 = _heme_site(2, 302, [143, 190])
    site2.residue_search_distances = {("N", 190, " "): 2.05, ("N", 143, " "): 5.08}
    site3 = _heme_site(3, 303, [27, 103])
    site3.residue_search_distances = {("N", 27, " "): 2.0, ("N", 103, " "): 2.1}
    results = [{"site_idx": 2, "site": site2, "status": "success", "error": None},
               {"site_idx": 3, "site": site3, "status": "success", "error": None}]

    detector._review_template_results("heme_bis_his_c_type", results, structure=None)

    out = buf.getvalue()
    assert "HIS N:143 5.1 Å" in out
    assert "HIS N:143 also in site 2" in out and "HIS N:143 also in site 1" in out
    assert "2 site(s) share a residue" in out
    assert "site 1, site 2" in out
    assert detector.final_sites == [site1, site2, site3]
