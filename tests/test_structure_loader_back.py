"""Structure Loader pickers say 'back', not 'cancel' (which read as
cancelling the structure just loaded). The old words still work so
recorded sessions replay."""

import os

from rich.console import Console

import proprep.structure_prep.structure_loader as sl
from proprep.structure_prep.structure_loader import StructureLoaderModule


class _Workspace(dict):
    def has(self, k):
        return k in self

    def set(self, k, v):
        self[k] = v


class _Processor:
    def __init__(self):
        self.console = Console(file=open(os.devnull, "w"))
        self.workspace = _Workspace()

    def _get_workspace(self):
        return self.workspace


def _module():
    m = StructureLoaderModule()
    m.processor = _Processor()
    m.initialize()
    return m


def _run_picker(monkeypatch, method, answer):
    """Drive the AlphaFold/AlphaFill entry picker: gene-name prompt, a
    stubbed one-hit search, then ``answer`` at the entry prompt."""
    m = _module()
    seen = {}
    answers = iter(["abl1", "", answer])  # gene, organism, entry

    def fake_prompt(processor, prompt, **kw):
        if "Select entry" in prompt:
            seen["prompt"] = prompt
            seen["choices"] = kw.get("choices")
            seen["show_choices"] = kw.get("show_choices")
        return next(answers)

    class _Fetcher:
        def search_uniprot_by_gene(self, *a, **k):
            return [{"uniprot_id": "P00520", "protein_name": "x", "organism": "y",
                     "length": 1, "gene_names": "ABL1", "reviewed": True}]

        def retrieve_structure(self, *a, **k):
            raise AssertionError("should not download after 'back'")

    monkeypatch.setattr(sl, "prompt_with_context", fake_prompt)
    monkeypatch.setattr(m, "_get_alphafold_fetcher", lambda: _Fetcher())
    ws = m.processor.workspace
    getattr(m, method)(ws)
    return seen, ws


def test_alphafold_picker_offers_back_and_accepts_legacy_c(monkeypatch):
    for answer in ("b", "c"):
        seen, ws = _run_picker(monkeypatch, "_retrieve_by_gene_name", answer)
        assert "'b' to go back" in seen["prompt"]
        assert "cancel" not in seen["prompt"]
        assert seen["choices"][-2:] == ["b", "c"]
        assert seen["show_choices"] is False
        assert "alphafold_pdb_file" not in ws


def test_no_cancel_wording_in_loader_prompts():
    src = open(sl.__file__).read()
    assert "'c' to cancel" not in src
    assert "Type 'cancel'" not in src
    # Numbered menu items and their echoes: "6. Cancel", '"6": "Cancel"',
    # "[grey50]Cancelled[/grey50]".
    assert '. Cancel"' not in src
    assert '"Cancel"' not in src
    assert "Cancelled[/" not in src
