"""
Molecular Docking menu 3: chain ends and the protonation review.

A chain end is a real terminus (charged NH3+ or COO-) or a chain break (left
neutral, as if the chain went on). The evidence is read from the file:
REMARK 465 (residues not modelled), else SEQRES (every residue modelled
means both ends are the chain's own), else an OXT on the last residue. The
code used to skip every proposal whenever REMARK 465 was absent, which is
exactly the case of a complete chain (1SDU): Meeko then built Pro 1 as a
neutral end and Phe 99, which carries an OXT, as a charged one. An end with
no evidence is now asked about, with only the two choices that fit a first
or a last residue; a run needs every uncapped end decided.

The review said "pdb2pqr assigned 32 titratable residues" over a table of 3:
it now says how many are in their usual charged form and not listed, that
any residue can be changed, and 'a' lists all of them.

Run with: pytest tests/test_docking_termini_and_protonation_review.py
"""

import os
import shutil

import pytest

from proprep.docking_prep import receptor_decisions as rd
from tests.test_docking_module import _Processor, _scripted

DATA = os.path.join(os.path.dirname(__file__), "data", "docking")


def _copy(tmp_path, name, keep=lambda line: True):
    path = tmp_path / name.replace(".txt", "")
    with open(os.path.join(DATA, name)) as handle:
        path.write_text("".join(line for line in handle if keep(line)))
    return str(path)


# ------------------------------------------------------------- evidence ---

def test_a_complete_chain_has_real_termini_from_seqres(tmp_path):
    proposals = rd.propose_termini(_copy(tmp_path, "2FDN.pdb.txt"), ["A"])
    assert proposals["A:1"].value == "N-terminus" and proposals["A:55"].value == "C-terminus"
    assert proposals["A:1"].source == "SEQRES: all 55 residues modelled, none missing"


def test_without_seqres_an_oxt_marks_the_c_terminus_and_the_n_end_is_open(tmp_path):
    path = _copy(tmp_path, "2FDN.pdb.txt", keep=lambda line: not line.startswith("SEQRES"))
    proposals = rd.propose_termini(path, ["A"])
    assert proposals["A:55"].value == "C-terminus" and proposals["A:55"].source == "OXT modelled on this residue"
    assert "A:1" not in proposals


def test_remark_465_still_wins_and_a_capped_end_needs_nothing():
    assert rd.propose_termini(os.path.join(DATA, "2CBA.pdb.txt"), ["A"])["A:3"].value == "chain break"
    ends = {e.res_id: e for e in rd.chain_ends(os.path.join(DATA, "1HRC.pdb.txt"), ["A"])}
    assert ends["A:1"].capped and not ends["A:104"].capped


# ------------------------------------------------------------------ menu ---

@pytest.fixture
def module(tmp_path):
    from proprep.docking_prep.docking_module import MolecularDockingModule
    m = MolecularDockingModule()
    m.processor = _Processor(tmp_path)
    m.initialize()
    return m


def test_an_end_without_evidence_is_asked_with_only_the_choices_that_fit(module, tmp_path, monkeypatch):
    from proprep.docking_prep.docking_menus_receptor import _decide_open_ends
    r = module.state.receptor
    r.pdb_path, r.chains, r.ph = _copy(tmp_path, "2D5M_zinc_excerpt.pdb.txt"), ["A"], 7.0
    assert "chain ends (4)" in module.state.missing()
    asked = []
    import proprep.utils.prompts as prompts
    def answer(processor, prompt, *a, **k):
        asked.append((prompt, k))
        return {"A:39": "b", "A:150": "c"}[prompt.split()[0]]
    monkeypatch.setattr(prompts, "prompt_with_context", answer)
    _decide_open_ends(module)
    assert r.termini == {"A:39": "chain break", "A:150": "C-terminus"}
    assert r.termini_sources == {"A:39": "chosen", "A:150": "chosen"}
    first, last = asked
    assert first[1]["choices"] == ["n", "b"] and last[1]["choices"] == ["c", "b"]
    assert "default" not in first[1] or first[1]["default"] is None          # no answer is assumed
    assert "the first residue of chain A" in first[0]
    assert "chain ends (4)" not in module.state.missing()


def test_the_review_counts_what_it_hides_and_any_residue_can_be_changed(module, tmp_path, monkeypatch):
    pytest.importorskip("meeko", reason="Meeko is required for the template families")
    from proprep.docking_prep.docking_menus_receptor import _residues, _review_protonation
    r = module.state.receptor
    r.pdb_path, r.chains, r.ph = _copy(tmp_path, "2FDN.pdb.txt"), ["A"], 7.0
    r.protonation = {"A:27": "ASP", "A:28": "ASP", "A:6": "GLU"}
    r.protonation_sources = {k: "pdb2pqr/PROPKA at pH 7.0" for k in r.protonation}
    r.termini = {"A:1": "N-terminus", "A:55": "C-terminus"}
    r.termini_sources = {"A:1": "SEQRES: all 55 residues modelled, none missing",
                         "A:55": "SEQRES: all 55 residues modelled, none missing"}
    asked = []
    _scripted(monkeypatch, ["a", "A:27", "ASH", ""], asked)
    _review_protonation(module, _residues(r.pdb_path))
    shown = module.processor.out.getvalue()
    assert "Of the 3 titratable residues, 3 are in their usual charged form at pH 7.0 (Asp-, Glu-, Lys+, ...) " \
           "and 0 are not, or are His (which has no single usual form). Only those 0 are listed ('a' lists all)." \
        in shown.replace("\n", " ").replace("  ", " ")
    assert "Type any residue to change it, listed or not." in shown
    assert "Protonation: every titratable residue" in shown                  # 'a' lists them all
    assert r.protonation["A:27"] == "ASH" and r.protonation_sources["A:27"] == "chosen"
    assert "N-terminus: charged NH3+ (the chain starts here)" in shown
    assert "SEQRES: all 55 residues modelled, none missing" in shown


def test_the_review_is_drawn_in_the_viewer_and_a_residue_can_be_clicked(module, tmp_path, monkeypatch):
    """Each listed residue is drawn with its template as a label (side chain), each chain end with
    what it is built as (its terminal atoms); 'p' picks the residue to change by clicking it."""
    pytest.importorskip("meeko", reason="Meeko is required for the template families")
    from unittest.mock import MagicMock
    import proprep.structure_prep.viewer_coordinator as vc
    from proprep.docking_prep import docking_ui
    from proprep.docking_prep.docking_menus_receptor import _residues, _review_protonation
    fake = MagicMock()
    fake.current_structures.return_value = []
    monkeypatch.setattr(vc, "viewer", fake)
    r = module.state.receptor
    r.pdb_path, r.chains, r.ph = _copy(tmp_path, "2FDN.pdb.txt"), ["A"], 7.0
    r.protonation = {"A:27": "ASH", "A:6": "GLU"}
    r.protonation_sources = {"A:27": "pdb2pqr/PROPKA at pH 7.0", "A:6": "pdb2pqr/PROPKA at pH 7.0"}
    r.termini = {"A:1": "N-terminus", "A:55": "C-terminus"}
    monkeypatch.setattr(docking_ui, "pick", lambda m, kind, prompt, index: {
        "atom": {"index": 5, "name": "CD", "resname": "GLU", "resno": 6, "chain": "A", "inscode": ""}})
    _scripted(monkeypatch, ["p", "GLH", ""], [])
    _review_protonation(module, _residues(r.pdb_path))
    assert r.protonation["A:6"] == "GLH" and r.protonation_sources["A:6"] == "chosen"
    first = fake.replace_annotations.call_args_list[0].args[1]
    rows = {e["label"]: e for e in first}
    assert rows["dock_prot_A:27"]["selection"] == "(:A and 27) and sidechainAttached"
    assert rows["dock_prot_A:27"]["text"] == "ASH"
    assert rows["dock_prot_end_A:1"]["selection"] == "(:A and 1) and .N"
    assert rows["dock_prot_end_A:1"]["text"] == "ALA1 NH3+"
    assert rows["dock_prot_end_A:55"]["text"].endswith("COO-")
    assert "dock_prot_A:6" not in rows                                     # GLU is its usual form: not listed
    second = {e["label"]: e for e in fake.replace_annotations.call_args_list[1].args[1]}
    assert second["dock_prot_A:6"]["text"] == "GLH"                         # changed, so now listed
    assert fake.replace_annotations.call_args_list[-1].args == ("dock_prot_", [])   # cleared on leaving
