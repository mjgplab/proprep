"""MODELLER holds the resolved atoms fixed, or says that it cannot and asks.

ProPrep has MODELLER refine only the residues it rebuilds. When none of those
residues could be found in MODELLER's model, the code fell back to refining
EVERY atom, without a word, under a console line saying "all resolved atoms
(including any metal sites) are held fixed". With no rebuilt residue identified
at all, a plain AutoModel did the same. On a 6R2Q fragment the unfrozen run
left 0 of 487 resolved heavy atoms in place; the frozen one left all 487.

The run is now refused unless the user allows it, and after every run the
atoms still in place are counted from the two files, by position: MODELLER
renames every chain, renumbers every residue and exchanges the names of
equivalent atoms (OD1/OD2, NH1/NH2, ...) without moving them, so names cannot
be used.
"""

import io
import shutil
import types
from pathlib import Path

import pytest
from rich.console import Console

from proprep.structure_prep import structure_completeness as sc

FRAGMENT = Path(__file__).parent.parent / "data" / "6r2q_chainB_400_420.txt"
ONE_LETTER = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I",
                  LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
_NAMES = {int(l[22:26]): l[17:20] for l in FRAGMENT.read_text().splitlines() if l.startswith("ATOM")}
SEQUENCE = "".join(ONE_LETTER[_NAMES[n]] for n in sorted(_NAMES))      # chain B 400-420, read from the fixture
GAP = range(409, 412)                            # cut out, for MODELLER to rebuild
BUILT = {("A", n - 399) for n in GAP}            # MODELLER's numbering: one chain A, from 1

needs_modeller = pytest.mark.skipif(not sc.HAS_MODELLER, reason="MODELLER not installed or not licensed")


def _interface():
    out = io.StringIO()
    interface = sc.ModellerInterface(Console(file=out, width=300))
    for name in ("_extract_assessment", "_assess_built_region", "_display_quality_assessment"):
        setattr(interface, name, lambda *a, **k: None)
    return interface, out


def _work_dir(tmp_path):
    lines = [l for l in FRAGMENT.read_text().splitlines(keepends=True) if l.startswith("ATOM")]
    (tmp_path / "template.pdb").write_text("".join(l for l in lines if int(l[22:26]) not in GAP) + "END\n")
    gapped = "".join("-" if 400 + i in GAP else aa for i, aa in enumerate(SEQUENCE))
    (tmp_path / "aln.pir").write_text(
        f">P1;template\nstructureX:template:400:B:420:B::::\n{gapped}*\n"
        f">P1;target\nsequence:target:::::::0.00: 0.00\n{SEQUENCE}*\n")
    return str(tmp_path)


# ── With real MODELLER ──────────────────────────────────────────────────

@needs_modeller
def test_frozen_run_leaves_every_resolved_atom_where_it_was(tmp_path):
    interface, out = _interface()
    ok, message, structure = interface.run_modeller(_work_dir(tmp_path), "template.pdb", "aln.pir",
                                                    built_residues=BUILT)
    report = interface.freeze_report
    assert ok and structure is not None and report["frozen"] is True
    assert report["resolved"] == len(BUILT) and report["unresolved"] == []
    assert report["atoms_in_place"] == report["atoms_resolved"] > 100
    assert "exactly where they were" in out.getvalue()
    # MODELLER did rename and renumber: by name nothing would have been found.
    assert {l[21] for l in (tmp_path / "target.B99990001.pdb").read_text().splitlines() if l.startswith("ATOM")} == {"A"}


@needs_modeller
def test_rebuilt_residues_that_are_not_in_the_model_stop_the_run(tmp_path):
    """Original numbering passed where MODELLER's was expected: nothing resolves."""
    interface, _ = _interface()
    wrong = {("B", n) for n in GAP}
    ok, message, structure = interface.run_modeller(_work_dir(tmp_path), "template.pdb", "aln.pir",
                                                    built_residues=wrong)
    assert (ok, message, structure) == (False, sc.FREEZE_NOT_APPLICABLE, None)
    assert "none of the 3 rebuilt residues was found" in interface.freeze_report["not_applicable"]
    assert "B:409" in interface.freeze_report["not_applicable"]
    assert not (tmp_path / "target.B99990001.pdb").exists()          # stopped before refining anything


@needs_modeller
def test_allowed_unfrozen_run_says_how_many_atoms_moved(tmp_path):
    interface, out = _interface()
    ok, _, _ = interface.run_modeller(_work_dir(tmp_path), "template.pdb", "aln.pir",
                                      built_residues={("B", n) for n in GAP}, allow_unfrozen=True)
    report = interface.freeze_report
    assert ok and report["frozen"] is False
    assert report["atoms_in_place"] < report["atoms_resolved"]
    assert "Refined without holding the resolved atoms" in out.getvalue()


@needs_modeller
def test_a_rebuilt_residue_that_is_not_found_is_named(tmp_path):
    interface, out = _interface()
    ok, _, _ = interface.run_modeller(_work_dir(tmp_path), "template.pdb", "aln.pir",
                                      built_residues=(BUILT - {("A", 12)}) | {("A", 999)})
    assert ok and interface.freeze_report["unresolved"] == ["A:999"]
    assert "A:999" in out.getvalue() and "not refined" in out.getvalue()
    assert interface.freeze_report["atoms_in_place"] == interface.freeze_report["atoms_resolved"]


# ── Without MODELLER ────────────────────────────────────────────────────

def test_no_rebuilt_residue_identified_is_refused_before_modeller_starts(monkeypatch, tmp_path):
    monkeypatch.setattr(sc, "HAS_MODELLER", True)
    interface, _ = _interface()
    assert interface.run_modeller(str(tmp_path), "template.pdb", "aln.pir", built_residues=set()) == \
        (False, sc.FREEZE_NOT_APPLICABLE, None)
    assert interface.freeze_report["not_applicable"] == "no rebuilt residue was identified"


def _module(monkeypatch, answer):
    module = sc.StructureCompletenessModule.__new__(sc.StructureCompletenessModule)
    out = io.StringIO()
    module.processor = types.SimpleNamespace(console=Console(file=out, width=300))
    if not isinstance(getattr(type(module), "console", None), property):
        module.console = module.processor.console
    asked = []

    def confirm(processor, prompt, **kwargs):
        asked.append((prompt, kwargs.get("default")))
        return answer

    monkeypatch.setattr(sc, "confirm_with_context", confirm)
    calls = []

    def run_modeller(work_dir, input_pdb, aln_file, built_residues=None, allow_unfrozen=False):
        calls.append(allow_unfrozen)
        return (True, "Success", "structure") if allow_unfrozen else (False, sc.FREEZE_NOT_APPLICABLE, None)

    interface = types.SimpleNamespace(run_modeller=run_modeller,
                                      freeze_report={"not_applicable": "none of the 4 rebuilt residues was found"})
    return module, interface, asked, calls, out


def test_the_user_is_told_and_asked_and_the_default_is_to_stop(monkeypatch):
    module, interface, asked, calls, out = _module(monkeypatch, answer=False)
    result = module._run_modeller_holding_resolved_atoms(interface, ".", "input.pdb", "aln.ali", {("A", 1)})
    assert result == (False, sc.REPAIR_DECLINED, None)
    assert asked == [("Run MODELLER anyway, letting every atom move?", False)]
    assert calls == [False]                                   # never run unfrozen
    said = out.getvalue()
    assert "cannot be held fixed" in said and "none of the 4 rebuilt residues" in said
    assert "metal sites" in said and "The structure is unchanged" in said


def test_unfrozen_run_happens_only_when_chosen(monkeypatch):
    module, interface, _, calls, _ = _module(monkeypatch, answer=True)
    result = module._run_modeller_holding_resolved_atoms(interface, ".", "input.pdb", "aln.ali", {("A", 1)})
    assert result == (True, "Success", "structure") and calls == [False, True]


def test_an_ordinary_modeller_failure_is_passed_on_without_asking(monkeypatch):
    module, _, asked, _, _ = _module(monkeypatch, answer=True)
    failing = types.SimpleNamespace(run_modeller=lambda *a, **k: (False, "MODELLER error: x", None), freeze_report={})
    assert module._run_modeller_holding_resolved_atoms(failing, ".", "i", "a", {("A", 1)})[1] == "MODELLER error: x"
    assert asked == []
