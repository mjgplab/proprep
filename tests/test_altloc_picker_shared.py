"""
One alternate-location picker for the Structure Fixer and Molecular Docking.

Docking used to ask "Alternate location to use for all of them", then
"Residue to choose individually", which read as one question asked twice and
applied one letter to every residue whatever its occupancies. It now runs the
Structure Fixer's picker residue by residue: occupancies, atoms covered, a
partial alternate flagged and completed from another, and the same viewer
colours. These tests pin that both callers ask the same question, that a
partial alternate keeps its residue whole in both, and that docking's later
steps (pdb2pqr's input, the ligand built from a residue) honour the fills.

Run with: pytest tests/test_altloc_picker_shared.py
"""

import io
import os
import shutil

import pytest
from rich.console import Console

from proprep.structure_prep import altloc_picker

DATA = os.path.join(os.path.dirname(__file__), "data", "docking")

# LEU 129 as PDB 8ZST models it: A and B complete, C the side chain only.
LEU_129 = """\
ATOM      1  N  AARG A 128      10.000  10.000  10.000  0.60 10.00           N
ATOM      2  CA AARG A 128      11.000  10.000  10.000  0.60 10.00           C
ATOM      3  C   ARG A 128      12.000  10.000  10.000  1.00 10.00           C
ATOM      4  O   ARG A 128      12.500  10.000  10.000  1.00 10.00           O
ATOM      5  N  ALEU A 129      13.000  10.000  10.000  0.60 10.00           N
ATOM      6  CA ALEU A 129      14.000  10.000  10.000  0.60 10.00           C
ATOM      7  C  ALEU A 129      15.000  10.000  10.000  0.60 10.00           C
ATOM      8  O  ALEU A 129      15.500  10.000  10.000  0.60 10.00           O
ATOM      9  CB ALEU A 129      14.000  11.000  10.000  0.60 10.00           C
ATOM     10  N  BLEU A 129      13.100  10.100  10.000  0.30 10.00           N
ATOM     11  CA BLEU A 129      14.100  10.100  10.000  0.30 10.00           C
ATOM     12  C  BLEU A 129      15.100  10.100  10.000  0.30 10.00           C
ATOM     13  O  BLEU A 129      15.600  10.100  10.000  0.30 10.00           O
ATOM     14  CB BLEU A 129      14.100  11.100  10.000  0.30 10.00           C
ATOM     15  CB CLEU A 129      14.200  11.200  10.000  0.10 10.00           C
END
"""


class _Workspace(dict):
    def set(self, key, value):
        self[key] = value


class _Processor:
    def __init__(self, folder):
        self.out = io.StringIO()
        self.console = Console(file=self.out, width=160)
        self.workspace = _Workspace(working_directory=str(folder))

    def _get_workspace(self):
        return self.workspace


def _scripted(monkeypatch, answers):
    """Answer prompts in order; returns the list of (prompt, kwargs) asked."""
    import proprep.utils.prompts as prompts
    asked = []

    def answer(processor=None, prompt=None, *args, **kwargs):
        asked.append((prompt, kwargs))
        assert answers, f"no scripted answer left for: {prompt}"
        value = answers.pop(0)
        return str(kwargs["default"]) if value == "" and kwargs.get("default") is not None else value

    monkeypatch.setattr(prompts, "prompt_with_context", answer)
    monkeypatch.setattr(prompts, "confirm_with_context",
                        lambda processor=None, prompt=None, **k: answer(processor, prompt, **{**k, "default": None}).startswith("y"))
    return asked


# ----------------------------------------------------------- shared model ---

def test_a_partial_alternate_is_flagged_and_completed():
    residues = altloc_picker.residues_from_records(
        (l[21], int(l[22:26]), l[26].strip(), l[17:20], l[12:16].strip(), l[16].strip(), float(l[54:60]))
        for l in LEU_129.splitlines() if l.startswith("ATOM"))
    leu = residues[("A", 129, "")]
    assert leu.letters == ["A", "B", "C"]
    assert leu.missing("C") == ["C", "CA", "N", "O"] and leu.missing("A") == []
    assert leu.fills("C") == {"N": "A", "CA": "A", "C": "A", "O": "A"}
    assert leu.average_occupancy("B") == "0.30"
    assert residues[("A", 128, "")].letters == ["A"]          # a lone label: nothing to choose


def test_keeps_atom():
    assert altloc_picker.keeps_atom("CB", "", "C", {})
    assert altloc_picker.keeps_atom("CB", "C", "C", {})
    assert altloc_picker.keeps_atom("N", "A", "C", {"N": "A"})
    assert not altloc_picker.keeps_atom("N", "B", "C", {"N": "A"})
    assert not altloc_picker.keeps_atom("N", "A", "B", None)


# ------------------------------------------------------- Structure Fixer ---

def test_structure_fixer_picks_per_residue_through_the_shared_picker(tmp_path, monkeypatch):
    from Bio.PDB import PDBParser
    from proprep.structure_prep.structure_completeness import StructureCompletenessModule
    pdb = tmp_path / "in.pdb"
    pdb.write_text(LEU_129)
    structure = PDBParser(QUIET=True).get_structure("t", str(pdb))
    fixer = StructureCompletenessModule.__new__(StructureCompletenessModule)
    fixer.processor = _Processor(tmp_path)
    fixer.console = fixer.processor.console
    fixer.results = {"alternate_locations": {"altloc_identifier": {"A": {"LEU_129": {"A", "B", "C"}}}}}
    asked = _scripted(monkeypatch, ["n", "3"])

    cleaned = fixer._handle_alternate_locations(structure, fixer.processor.workspace, save_to_workspace=False)

    assert [p for p, _ in asked] == [altloc_picker.VIEWER_PROMPT, "Select alternate to keep"]
    assert asked[1][1]["module"] == "Structure Completeness - Altloc"            # replays keep matching
    assert asked[1][1]["description"] == "Select alternate for LEU A:129"
    leu = next(r for r in cleaned.get_residues() if r.id[1] == 129)
    assert {a.get_name() for a in leu} == {"N", "CA", "C", "O", "CB"}
    assert leu["CB"].coord[0] == pytest.approx(14.2) and leu["N"].coord[0] == pytest.approx(13.0)
    shown = fixer.processor.out.getvalue()
    assert "3. Alternate C (occupancy: 0.10, 1 of 5 atoms)" in shown
    assert "keeping C from A, CA from A, N from A, O from A" in shown


# --------------------------------------------------------------- docking ---

class _Docking:
    """The parts of MolecularDockingModule the receptor menus use."""

    def __init__(self, folder, pdb, chains):
        from proprep.docking_prep.docking_state import DockingState
        self.processor = _Processor(folder)
        self.console = self.processor.console
        self.state = DockingState()
        self.state.receptor.pdb_path, self.state.receptor.chains = str(pdb), chains


def test_docking_asks_each_residue_and_offers_the_earlier_choice_again(tmp_path, monkeypatch):
    from proprep.docking_prep.docking_menus_receptor import choose_altlocs
    pdb = tmp_path / "2CBA.pdb"
    shutil.copy(os.path.join(DATA, "2CBA.pdb.txt"), pdb)
    m = _Docking(tmp_path, pdb, ["A"])
    asked = _scripted(monkeypatch, ["n", "2", "", "", ""])
    choose_altlocs(m, revisit=True)
    assert [p for p, _ in asked] == [altloc_picker.VIEWER_PROMPT] + ["Select alternate to keep"] * 4
    assert [k["description"] for _, k in asked[1:]] == [
        "Select alternate for HIS A:4", "Select alternate for GLU A:14",
        "Select alternate for HIS A:64", "Select alternate for GLN A:136"]
    assert all(k["module"] == "Molecular Docking" for _, k in asked)
    assert m.state.receptor.altlocs == {"A:4": "B", "A:14": "A", "A:64": "A", "A:136": "A"}
    shown = m.processor.out.getvalue()
    assert "Alternate location to use for all of them" not in shown
    assert "4 residue(s) have alternate locations; choose one for each." in shown

    # Menu 2 does not ask again about residues already decided ...
    asked.clear()
    choose_altlocs(m)
    assert asked == []
    # ... menu 1 does, with the earlier choice as the default, and a change clears protonation.
    m.state.receptor.ph, m.state.receptor.protonation = 7.0, {"A:4": "HID"}
    asked = _scripted(monkeypatch, ["n", "", "2", "", ""])
    choose_altlocs(m, revisit=True)
    assert [k["default"] for _, k in asked[1:]] == ["2", "1", "1", "1"]
    assert m.state.receptor.altlocs["A:14"] == "B"
    assert m.state.receptor.ph is None and m.state.receptor.protonation == {}
    assert "Changed alternates (A:14) move atoms" in m.processor.out.getvalue()


def test_docking_keeps_a_partial_alternate_residue_whole(tmp_path, monkeypatch):
    from proprep.docking_prep.docking_menus_receptor import choose_altlocs
    from proprep.docking_prep.receptor_decisions import _protein_only_copy
    from proprep.docking_prep.receptor_prep import read_atoms
    pdb = tmp_path / "leu.pdb"
    pdb.write_text(LEU_129)
    m = _Docking(tmp_path, pdb, ["A"])
    asked = _scripted(monkeypatch, ["n", "3"])
    choose_altlocs(m)
    r = m.state.receptor
    assert len(asked) == 2                                       # ARG 128 has one letter: kept, not asked
    assert r.altlocs == {"A:128": "A", "A:129": "C"}
    assert r.altloc_fills == {"A:129": {"N": "A", "CA": "A", "C": "A", "O": "A"}}
    assert "Only one alternate location is present, so it is kept: A:128 (A)." in m.processor.out.getvalue()

    out = tmp_path / "protein.pdb"
    _protein_only_copy(str(pdb), ["A"], r.altlocs, str(out), r.altloc_fills)
    leu = {a.name: a.xyz for a in read_atoms(str(out)) if a.resseq == 129}
    assert set(leu) == {"N", "CA", "C", "O", "CB"}
    assert leu["CB"][0] == pytest.approx(14.2) and leu["N"][0] == pytest.approx(13.0)
    assert all(a.altloc == "" for a in read_atoms(str(out)))
    # the saved state round-trips the fills, and the receptor build receives them
    from proprep.docking_prep.docking_state import DockingState
    again = DockingState.from_dict(m.state.to_dict())
    assert again.receptor.build_arguments()["altloc_fills"] == r.altloc_fills


def test_a_ligand_alternate_that_lacks_atoms_takes_them_from_another(tmp_path):
    Chem = pytest.importorskip("rdkit.Chem", reason="RDKit is required for docking ligand intake")
    pytest.importorskip("gemmi", reason="gemmi is required to read CCD entries")
    from proprep.docking_prep import ligand_sources as ls
    from proprep.docking_prep.ccd_chemistry import component_from_block
    with open(os.path.join(DATA, "ccd_BTN.cif.txt")) as handle:
        btn = component_from_block(handle.read(), "BTN")
    with open(os.path.join(DATA, "1STP_BTN.pdb.txt")) as handle:
        lines = [l for l in handle.read().splitlines() if l.startswith("HETATM")]
    tail = {"C11", "O11", "O12", "C10"}                          # alternate B models only the ring end
    records = []
    for line in lines:
        records.append(line[:16] + "A" + line[17:])
        if line[12:16].strip() not in tail:
            records.append(line[:16] + "B" + line[17:30] + f"{float(line[30:38]) + 1.0:8.3f}" + line[38:])
    path = tmp_path / "altloc.pdb"
    path.write_text("\n".join(records) + "\nEND\n")

    residue = altloc_picker.residues_from_records(
        (l[21], int(l[22:26]), l[26].strip(), l[17:20], l[12:16].strip(), l[16].strip(), float(l[54:60]))
        for l in records)[("A", 300, "")]
    plan = residue.fills("B")
    assert set(plan) == tail and set(plan.values()) == {"A"}
    with pytest.raises(ValueError, match="missing from the structure"):
        ls.from_structure_residue(str(path), "A", 300, "", "BTN", btn, altloc="B")
    build = ls.from_structure_residue(str(path), "A", 300, "", "BTN", btn, altloc="B", altloc_fills=plan)
    assert build.source["altloc_fills"] == plan
    assert "which B does not model" in build.notes[0]
