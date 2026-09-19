"""Structure Loader: search the PDB by the ligand a structure contains.

The searcher is tested against canned RCSB responses; the loader flow against a
stand-in searcher. `PROPREP_NETWORK_TESTS=1` adds a run against the live RCSB
services (the flow was developed against them: tetrahydrobiopterin, H4B, and
its other redox states BHS / HBI / H2B / BIO)."""

import io
import json
import os
import sys
import types

import pytest
import requests
from rich.console import Console

import proprep.structure_prep.pdb_searcher as ps
import proprep.structure_prep.structure_loader as sl
import proprep.utils.prompts as prompts
from proprep.structure_prep.pdb_searcher import (
    LigandDefinition, PDBSearcher, PDBSearchError, PDBSearchResult, metals_in_formula,
)
from proprep.structure_prep.structure_loader import StructureLoaderModule

MODULE = "Structure Loader - Ligand Search"


# ---------------------------------------------------------------------------
# Searcher, against canned RCSB responses
# ---------------------------------------------------------------------------

class _Response:
    def __init__(self, payload=None, status=200):
        self.payload, self.status_code = payload, status

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.exceptions.HTTPError(f"{self.status_code}")

    def json(self):
        return self.payload


FORMULAS = {
    "H4B": "C9 H15 N5 O3", "BHS": "C9 H15 N5 O3", "HBI": "C9 H13 N5 O3", "BIO": "C9 H11 N5 O3",
    "HM9": "C9 H12 N5 O3 P", "4AB": "C9 H16 N6 O2",
}


def _fake_rcsb(monkeypatch, formula_hits, similar):
    """Route requests.post by URL and query shape to canned RCSB answers."""
    def post(url, json=None, timeout=None, **_):
        if url == PDBSearcher.GRAPHQL_API_URL:
            ids = json["variables"]["ids"]
            return _Response({"data": {"chem_comps": [
                {"chem_comp": {"id": i, "name": f"name of {i}", "formula": FORMULAS[i]},
                 "rcsb_chem_comp_descriptor": {"SMILES_stereo": "C"}}
                for i in ids if i in FORMULAS]}})
        params = json["query"]["parameters"]
        if params.get("type") == "formula":
            assert params["value"] == "C9 N5 O3" and params["match_subset"] is True
            if json["request_options"].get("return_counts"):
                return _Response({"total_count": len(formula_hits)})
            return _Response({"result_set": [{"identifier": i} for i in formula_hits]})
        if params.get("type") == "descriptor":
            return _Response({"result_set": [{"identifier": i, "score": s} for i, s in similar.items()]})
        raise AssertionError(f"unexpected query {json}")
    monkeypatch.setattr(ps.requests, "post", post)


def test_split_formula_keeps_element_case_and_counts_hydrogens():
    assert PDBSearcher.split_formula("C34 H32 Fe N4 O4") == ({"C": 34, "Fe": 1, "N": 4, "O": 4}, 32)
    assert PDBSearcher.split_formula("Zn") == ({"Zn": 1}, 0)
    assert PDBSearcher.split_formula("H4 N") == ({"N": 1}, 4)
    assert PDBSearcher.split_formula("") == ({}, 0)


def test_related_ligands_are_labelled_and_supersets_dropped(monkeypatch):
    # HM9 has the atoms of H4B plus phosphorus: the subset search returns it, it is not related
    _fake_rcsb(monkeypatch, formula_hits=["H4B", "BHS", "HBI", "BIO", "HM9"],
               similar={"H4B": 1.0, "BHS": 1.0, "HBI": 1.0, "4AB": 0.79})
    seed = LigandDefinition("H4B", "5,6,7,8-TETRAHYDROBIOPTERIN", FORMULAS["H4B"], smiles="C")

    related, notes = PDBSearcher().find_related_ligands(seed)

    assert notes == []
    assert {l.comp_id: l.relation for l in related} == {
        "BHS": "same formula", "HBI": "same atoms, -2 H", "BIO": "same atoms, -4 H",
        "4AB": "similar structure"}
    assert [l.comp_id for l in related][-1] == "4AB"          # lookalikes after the same-atom ligands
    assert {l.comp_id: l.similarity for l in related}["BIO"] is None


def test_related_ligand_formula_search_is_skipped_and_reported_when_too_broad(monkeypatch):
    _fake_rcsb(monkeypatch, formula_hits=["X"] * (PDBSearcher.RELATED_FORMULA_MAX + 1), similar={})
    seed = LigandDefinition("H4B", "", FORMULAS["H4B"], smiles="C")
    related, notes = PDBSearcher().find_related_ligands(seed)
    assert related == []
    assert len(notes) == 1 and str(PDBSearcher.RELATED_FORMULA_MAX) in notes[0]


def test_no_hits_is_empty_but_an_unreachable_service_raises(monkeypatch):
    monkeypatch.setattr(ps.requests, "post", lambda *a, **k: _Response(status=204))
    assert PDBSearcher()._post_search({}) == {}

    def down(*a, **k):
        raise requests.exceptions.ConnectionError("no route to host")
    monkeypatch.setattr(ps.requests, "post", down)
    with pytest.raises(PDBSearchError):
        PDBSearcher()._post_search({})
    with pytest.raises(PDBSearchError):
        PDBSearcher().get_ligand_definitions(["H4B"])


def test_details_come_back_in_the_requested_order_with_matched_ligands(monkeypatch):
    def entry(pdb_id, comps):
        return {"rcsb_id": pdb_id, "struct": {"title": f"title {pdb_id}"}, "exptl": [{"method": "X-RAY DIFFRACTION"}],
                "rcsb_entry_info": {"resolution_combined": [1.5]},
                "rcsb_accession_info": {"initial_release_date": "2002-05-22T00:00:00Z"},
                "polymer_entities": [{"rcsb_entity_source_organism": [{"scientific_name": "Homo sapiens"}],
                                      "rcsb_polymer_entity_container_identifiers": {"auth_asym_ids": ["A", "B"]}}],
                "nonpolymer_entities": [
                    {"nonpolymer_comp": {"chem_comp": {"id": c, "name": c, "formula": ""}},
                     "nonpolymer_entity_instances": [
                         {"rcsb_nonpolymer_instance_validation_score": [{"RSCC": r}]} for r in rscc]}
                    for c, rscc in comps.items()]}
    payload = {"data": {"entries": [entry("2G6H", {"H4B": [0.91, 0.96], "HEM": [0.99]}), entry("1J8U", {"FE2": []})]}}
    monkeypatch.setattr(ps.requests, "post", lambda *a, **k: _Response(payload))

    first, second, missing = PDBSearcher().get_structure_details(["1j8u", "2G6H", "9ZZZ"], ligand_ids=["h4b"])

    assert (first.pdb_id, first.matched_ligands, first.ligand_rscc) == ("1J8U", [], None)
    assert (second.pdb_id, second.matched_ligands, second.ligand_rscc) == ("2G6H", ["H4B"], 0.96)
    assert second.organism == "Homo sapiens" and second.chains == ["A", "B"]
    assert [l["id"] for l in second.ligand_info] == ["H4B", "HEM"] and second.has_ligands
    assert missing.title == "[Details unavailable]"


# ---------------------------------------------------------------------------
# Loader flow, against a stand-in searcher
# ---------------------------------------------------------------------------

def _result(pdb_id, resolution, ligand):
    return PDBSearchResult(pdb_id=pdb_id, title=f"title of {pdb_id}", experimental_method="X-RAY DIFFRACTION",
                           resolution=resolution, release_date="2002", organism=None, chains=["A"],
                           has_ligands=True, ligand_info=[{"id": ligand}], matched_ligands=[ligand], ligand_rscc=0.9)


class _Searcher:
    """What RCSB knows, in miniature."""
    NAME_MATCH_ROWS = 50
    MAX_ROWS = 10000
    known = {"H4B": "5,6,7,8-TETRAHYDROBIOPTERIN", "HBI": "7,8-DIHYDROBIOPTERIN", "HEM": "HEME"}
    entries = [_result("1LTZ", 1.4, "HBI"), _result("1J8U", 1.5, "H4B"), _result("2G6H", 2.0, "H4B"),
               _result("1A6M", 1.0, "HEM")]
    calls = []

    def get_ligand_definitions(self, codes):
        return [LigandDefinition(c, self.known[c], "C9 H15 N5 O3") for c in codes if c in self.known]

    def lookup_ligands(self, text):
        found = self.get_ligand_definitions([text.upper()])
        for ligand in found:
            ligand.relation = "code"
        return found, 0

    def find_related_ligands(self, seed):
        return [LigandDefinition("HBI", self.known["HBI"], "C9 H13 N5 O3", relation="same atoms, -2 H", similarity=1.0)], []

    def count_entries_by_ligand(self, codes):
        per_code = {c: sum(1 for e in self.entries if c in e.matched_ligands) for c in codes}
        return sum(per_code.values()), per_code

    def count_by_name(self, protein_name, organism=None):
        return len(self.entries)

    def search_and_filter(self, protein_name, organism=None, max_results=None):
        type(self).calls.append((protein_name, organism, max_results))
        return self.entries[:max_results], len(self.entries)

    def search_by_ligands(self, codes, max_results):
        type(self).calls.append((list(codes), max_results))
        hits = [e for e in self.entries if set(e.matched_ligands) & set(codes)]
        return hits[:max_results], len(hits)


def _loader(console=None, session_manager=None):
    m = StructureLoaderModule()
    m.processor = types.SimpleNamespace(
        console=console or Console(file=open(os.devnull, "w")), session_manager=session_manager)
    if session_manager is None:
        del m.processor.session_manager
    m.initialize()
    if console is not None:
        m.console = console          # initialize() makes its own
    loaded = []
    m._load_selected_structures = lambda results, workspace: loaded.extend(r.pdb_id for r in results) or workspace
    return m, loaded


def _scripted_prompts(monkeypatch, answers):
    """Answer prompts in order; keep every call for the recording-contract checks."""
    answers, seen = iter(answers), []

    def fake_prompt(processor, prompt, choices=None, default=None, module=None,
                    description=None, options_map=None, show_choices=None):
        seen.append(dict(prompt=prompt, choices=choices, default=default, module=module,
                         description=description, options_map=options_map))
        return next(answers)
    monkeypatch.setattr(sl, "prompt_with_context", fake_prompt)
    monkeypatch.setattr(prompts, "prompt_with_context", fake_prompt)   # the retry helpers call it from there
    return seen


def test_flow_from_code_through_related_ligands_to_a_loaded_structure(monkeypatch):
    monkeypatch.setattr(sl, "PDBSearcher", _Searcher)
    _Searcher.calls.clear()
    seen = _scripted_prompts(monkeypatch, ["h4b", "h4b", "H4B, hbi", "2", "2"])
    m, loaded = _loader()

    m._search_pdb_by_ligand_interactive({})

    assert _Searcher.calls == [(["H4B", "HBI"], 2)]      # the codes typed, and the limit typed
    assert loaded == ["1J8U"]                            # row 2 of the two best-resolved
    assert seen[2]["default"] == "H4B"                   # the code the user entered is the offered default
    for call in seen:                                    # what the session editor shows
        assert call["module"] in (MODULE, "Structure Loader - PDB Search")
        assert call["description"]


def test_prompt_text_does_not_depend_on_the_ligand(monkeypatch):
    """Replay matches prompt text exactly, so no prompt may carry search state."""
    monkeypatch.setattr(sl, "PDBSearcher", _Searcher)
    texts = []
    for ligand in ("H4B", "HEM"):
        seen = _scripted_prompts(monkeypatch, [ligand, "", ligand, "", "1"])
        _loader()[0]._search_pdb_by_ligand_interactive({})
        texts.append([call["prompt"] for call in seen])
    assert texts[0] == texts[1]


def test_an_unknown_code_is_reported_and_asked_again_never_dropped(monkeypatch):
    monkeypatch.setattr(sl, "PDBSearcher", _Searcher)
    _Searcher.calls.clear()
    console = Console(width=200, record=True, file=io.StringIO(), force_terminal=False)
    _scripted_prompts(monkeypatch, ["H4B", "", "H4B,QQQ", "", "H4B", "", "1"])
    m, loaded = _loader(console)

    m._search_pdb_by_ligand_interactive({})

    out = console.export_text()
    assert "Not a PDB ligand code: QQQ" in out
    assert "Enter at least one ligand code" in out
    assert _Searcher.calls == [(["H4B"], 100)]           # searched only once every code was valid


def test_back_at_the_codes_prompt_searches_nothing(monkeypatch):
    monkeypatch.setattr(sl, "PDBSearcher", _Searcher)
    _Searcher.calls.clear()
    _scripted_prompts(monkeypatch, ["H4B", "", "back"])
    m, loaded = _loader()
    workspace = {"kept": True}
    assert m._search_pdb_by_ligand_interactive(workspace) is workspace
    assert _Searcher.calls == [] and loaded == []


def test_an_unreachable_service_is_said_so_not_reported_as_no_results(monkeypatch):
    class Down(_Searcher):
        def lookup_ligands(self, text):
            raise PDBSearchError("no route to host")
    monkeypatch.setattr(sl, "PDBSearcher", Down)
    console = Console(width=200, record=True, file=io.StringIO(), force_terminal=False)
    _scripted_prompts(monkeypatch, ["H4B"])
    m, _ = _loader(console)
    workspace = {"kept": True}

    assert m._search_pdb_by_ligand_interactive(workspace) is workspace
    out = console.export_text()
    assert "could not be reached" in out and "no route to host" in out
    assert "No PDB ligand" not in out


def test_rcsb_menu_offers_the_search_and_keeps_the_old_labels(monkeypatch):
    seen = _scripted_prompts(monkeypatch, ["4"])
    m, _ = _loader()
    m._search_pdb_by_ligand_interactive = lambda workspace: "searched"

    assert m._load_from_pdb_database({}) == "searched"
    menu = seen[0]
    assert list(menu["options_map"]) == menu["choices"]
    # Recorded sessions find "Back" (4 before this option existed) by its label
    assert menu["options_map"]["5"] == "Back to source selection"
    assert [menu["options_map"][k] for k in "123"] == [
        "Enter PDB ID directly", "Search by entry title keyword", "Search by gene/protein (UniProt)"]


def test_recorded_selection_replays_by_pdb_id_when_the_results_have_moved(monkeypatch, tmp_path):
    """The PDB grows: on replay a better-resolved entry has pushed the chosen one down a row."""
    from proprep.utils.session_recorder import SessionManager
    monkeypatch.setattr(sl, "PDBSearcher", _Searcher)
    typed = ["H4B", "", "H4B", "", "1"]                  # row 1 = 1J8U when recorded
    recorded = ["H4B", "", "H4B", "100", "1"]            # Enter at the limit prompt records its default

    def run(stdin_text, start):
        monkeypatch.setattr(sys, "stdin", io.StringIO(stdin_text))
        manager = SessionManager()
        m, loaded = _loader(Console(width=200, file=io.StringIO(), force_terminal=False), manager)
        start(manager)
        try:
            m._search_pdb_by_ligand_interactive({})
        finally:
            consumed = manager.replayer.interaction_index if manager.replayer else None
            manager.stop()
        return consumed, loaded

    log = str(tmp_path / "session.json")
    _, loaded = run("\n".join(typed) + "\n", lambda sm: sm.start_recording(log, {}))
    assert loaded == ["1J8U"]
    interactions = json.load(open(log))["interactions"]
    assert [i["response"] for i in interactions] == recorded
    assert interactions[-1]["context"]["selected_key"] == "1J8U"

    monkeypatch.setattr(_Searcher, "entries", [_result("9NEW", 1.1, "H4B")] + _Searcher.entries)
    consumed, loaded = run("", lambda sm: sm.start_replay(log))
    assert consumed == len(interactions)                 # no live prompt
    assert loaded == ["1J8U"]                            # not 9NEW, which now holds row 1


# ---------------------------------------------------------------------------
# Title search: how many structures are listed is the user's choice
# ---------------------------------------------------------------------------

def test_title_search_reports_the_count_and_asks_how_many_to_list(monkeypatch):
    monkeypatch.setattr(sl, "PDBSearcher", _Searcher)
    _Searcher.calls.clear()
    console = Console(width=200, record=True, file=io.StringIO(), force_terminal=False)
    seen = _scripted_prompts(monkeypatch, ["hydroxylase", "", "3", "back"])
    m, _ = _loader(console)

    m._search_pdb_interactive({})

    assert _Searcher.calls == [("hydroxylase", None, 3)]          # the typed limit reaches the search
    assert "4 structures match" in console.export_text()          # said before the limit is asked
    limit = seen[2]
    assert limit["prompt"] == ("Maximum number of structures to list, most relevant first "
                               "(range: 1-10000, default: 100)")
    assert limit["default"] == "100" and limit["module"] == "Structure Loader - PDB Search" and limit["description"]


def test_title_search_offline_is_said_so(monkeypatch):
    class Down(_Searcher):
        def count_by_name(self, protein_name, organism=None):
            raise PDBSearchError("no route to host")
    monkeypatch.setattr(sl, "PDBSearcher", Down)
    console = Console(width=200, record=True, file=io.StringIO(), force_terminal=False)
    _scripted_prompts(monkeypatch, ["hydroxylase", ""])
    m, _ = _loader(console)

    m._search_pdb_interactive({})

    out = console.export_text()
    assert "could not be reached" in out and "No structures found" not in out


def test_search_and_filter_fetches_what_was_asked_unless_filtering(monkeypatch):
    fetched = []
    searcher = PDBSearcher()
    monkeypatch.setattr(searcher, "search_by_name", lambda n, o, max_results: fetched.append(max_results) or ([], 0))
    searcher.search_and_filter("kinase", max_results=7000)
    searcher.search_and_filter("kinase", max_results=7000, max_resolution=2.0)
    assert fetched == [7000, PDBSearcher.MAX_ROWS]                # never more than RCSB returns


# ---------------------------------------------------------------------------
# Metals: read from the ligands' formulas, never from a list of ligand codes
# ---------------------------------------------------------------------------

def _with_ligands(pdb_id, ligands):
    result = _result(pdb_id, 2.0, "H4B")
    result.ligand_info = [{"id": code, "formula": formula} for code, formula in ligands.items()]
    return result


def test_a_metal_is_any_element_that_is_not_a_nonmetal():
    assert metals_in_formula("C34 H32 Fe N4 O4") == ["Fe"]          # heme
    assert metals_in_formula("Fe4 S4") == ["Fe"]                    # iron-sulfur cluster
    assert metals_in_formula("C62 H89 Co N13 O14 P") == ["Co"]      # cobalamin
    assert metals_in_formula("Cm") == ["Cm"]                        # curium: on no list of metals ProPrep has
    assert metals_in_formula("C2 H6 As O2") == []                   # cacodylate: a metalloid is not a metal
    assert metals_in_formula("C5 H11 N O2 Se") == []
    assert metals_in_formula("D2 O") == [] and metals_in_formula("X") == [] and metals_in_formula("") == []
    assert metals_in_formula("ZN") == ["Zn"]                        # whatever the case


def test_metal_components_of_a_result():
    nos = _with_ligands("2G6H", {"ARG": "C6 H15 N4 O2", "H4B": "C9 H15 N5 O3", "HEM": "C34 H32 Fe N4 O4",
                                 "ZN": "Zn", "CAC": "C2 H6 As O2"})
    # the code, which for an ion carries the deposited oxidation state; the metal named when it is inside a cofactor
    assert nos.metal_components() == ["HEM (Fe)", "ZN"] and nos.metals() == ["Fe", "Zn"]
    assert nos.metal_components("fe") == ["HEM (Fe)"] and nos.metal_components("Cu") == []
    assert _with_ligands("1J8U", {"FE2": "Fe"}).metal_components() == ["FE2"]
    assert _with_ligands("X", {"SF4": "Fe4 S4", "FCO": "C3 Fe N2 O", "CUZ": "Cu4 S"}).metal_components() == [
        "CUZ (Cu)", "FCO (Fe)", "SF4 (Fe)"]


def _filter_metal(monkeypatch, results, typed):
    console = Console(width=250, record=True, file=io.StringIO(), force_terminal=False)
    seen = _scripted_prompts(monkeypatch, ["3", typed])
    m, _ = _loader(console)
    return [r.pdb_id for r in m._apply_search_filters(results)], console.export_text(), seen


def test_metal_filter_finds_iron_whatever_ligand_holds_it(monkeypatch):
    results = [_with_ligands("1LTZ", {"FE": "Fe"}), _with_ligands("1J8U", {"FE2": "Fe"}),
               _with_ligands("2G6H", {"HEM": "C34 H32 Fe N4 O4", "ZN": "Zn"}),
               _with_ligands("2BFP", {"NAP": "C21 H28 N7 O17 P3"}), _with_ligands("1B66", {"ZN": "Zn"})]

    kept, _, seen = _filter_metal(monkeypatch, results, "FE")       # as an old recording typed it
    assert kept == ["1LTZ", "1J8U", "2G6H"]
    assert seen[1]["prompt"] == "Metal type (e.g., Zn, Fe, Cu) or leave blank for any metal"   # unchanged, for replay

    kept, out, _ = _filter_metal(monkeypatch, results, "")
    assert kept == ["1LTZ", "1J8U", "2G6H", "1B66"] and "(Fe, Zn)" in out


def test_metal_filter_reports_an_absent_metal_and_keeps_the_results(monkeypatch):
    results = [_with_ligands("2G6H", {"HEM": "C34 H32 Fe N4 O4"}), _with_ligands("2BFP", {"NAP": "C21 H28 N7 O17 P3"})]
    kept, out, _ = _filter_metal(monkeypatch, results, "hem")       # a ligand code, not an element
    assert kept == ["2G6H", "2BFP"]
    assert "contains 'Hem'" in out and "Metals present: Fe" in out and "unchanged" in out

    kept, out, _ = _filter_metal(monkeypatch, [results[1]], "")
    assert kept == ["2BFP"] and "unchanged" in out


def test_results_table_has_a_metals_column_and_no_marker():
    console = Console(width=200, record=True, file=io.StringIO(), force_terminal=False)
    m, _ = _loader(console)
    m._display_results_table([_with_ligands("2G6H", {"HEM": "C34 H32 Fe N4 O4", "ZN": "Zn", "ACT": "C2 H3 O2"})], 0, 1, 1, 1, 1)
    out = console.export_text()
    assert "HEM (Fe)  ZN" in out and "[M" not in out
    assert "Metal ligands" not in out          # to a coordination chemist that means the groups binding the metal


# ---------------------------------------------------------------------------
# Live RCSB services
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not os.environ.get("PROPREP_NETWORK_TESTS"), reason="set PROPREP_NETWORK_TESTS=1 to query RCSB")
def test_live_biopterin_redox_states():
    searcher = PDBSearcher()

    by_name, _ = searcher.lookup_ligands("tetrahydrobiopterin")
    assert {"H4B", "BHS"} <= {l.comp_id for l in by_name}
    assert searcher.lookup_ligands("h4b")[0][0].relation == "code"

    seed = searcher.get_ligand_definitions(["H4B"])[0]
    related, notes = searcher.find_related_ligands(seed)
    relation = {l.comp_id: l.relation for l in related}
    assert relation["BHS"] == "same formula"
    assert relation["HBI"] == relation["H2B"] == "same atoms, -2 H"
    assert relation["BIO"] == "same atoms, -4 H"
    assert relation["4AB"] == "similar structure"        # 4-amino analog: a lookalike, labelled as one

    codes = ["H4B", "BHS", "HBI", "H2B", "BIO"]
    total, per_code = searcher.count_entries_by_ligand(codes)
    assert per_code["H4B"] > 700 and per_code["HBI"] >= 28 and total <= sum(per_code.values())

    assert searcher.count_by_name("phenylalanine hydroxylase", "human") == \
        searcher.search_by_name("phenylalanine hydroxylase", "human", max_results=1)[1] > 0

    results, total = searcher.search_by_ligands(["HBI", "BIO"], 10)
    assert len(results) == 10 and total >= 34
    resolutions = [r.resolution for r in results]
    assert resolutions == sorted(resolutions)
    assert all(set(r.matched_ligands) & {"HBI", "BIO"} for r in results)
    assert results[0].pdb_id == "1LTZ" and results[0].ligand_rscc is not None
