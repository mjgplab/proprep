"""
Regression test: an unfilled INTERNAL gap can be closed with a TER record
instead of a MODELLER fill or ACE/NME caps.

Declining the repair plan (or picking "TER" on a single-residue gap, or having
no MODELLER) used to leave the two gap-bracketing residues bonded across the
break by tLEaP, producing a spurious long bond through a metal site. The repair
now offers to drop a TER record after the residue preceding each internal gap.

TER is a file-only edit (Bio.PDB's PDBIO does not round-trip it), so
_apply_ter_records post-processes the written PDB text. The anchor residue is
found under its *final* numbering via the ResidueMapper.
"""

import proprep.structure_prep.structure_completeness as sc
from proprep.structure_prep.structure_completeness import (
    ResidueMapper, RepairPlan, MissingSegment,
)
from rich.console import Console


def _module():
    for obj in vars(sc).values():
        if isinstance(obj, type) and "_apply_ter_records" in obj.__dict__:
            inst = obj.__new__(obj)
            inst.console = Console()
            return inst
    raise AssertionError("could not find the class owning _apply_ter_records")


def _atom(serial, name, resname, chain, resnum):
    return (
        f"ATOM  {serial:>5} {name:<4} {resname:<3} {chain}{resnum:>4}"
        f"      0.000   0.000   0.000  1.00  0.00           C\n"
    )


def _write_pdb(path, rows):
    # rows: list of (serial, name, resname, chain, resnum)
    with open(path, "w") as fh:
        for i, r in enumerate(rows, 1):
            fh.write(_atom(i, *r))
        fh.write("END\n")


def _identity_mapper():
    m = ResidueMapper(Console())
    m.chain_mapping = {}
    m.residue_mappings = {}
    m.final_mappings = {}
    return m


def test_ter_inserted_after_residue_before_internal_gap(tmp_path):
    pdb = tmp_path / "s.pdb"
    # Chain A: residues 5, 6 present, 7 missing (the gap), 8 present.
    _write_pdb(pdb, [
        ("N", "ALA", "A", 5), ("CA", "ALA", "A", 5),
        ("N", "GLY", "A", 6), ("CA", "GLY", "A", 6),
        ("N", "SER", "A", 8), ("CA", "SER", "A", 8),
    ])
    seg = MissingSegment(chain_id="A", residues=[], start_num=7, end_num=7,
                         is_terminal=False, terminal_type=None)
    plan = RepairPlan(segments_to_ter=[seg])

    _module()._apply_ter_records(str(pdb), plan, _identity_mapper())

    lines = pdb.read_text().splitlines()
    ter_idx = [i for i, l in enumerate(lines) if l.startswith("TER")]
    assert len(ter_idx) == 1
    # TER must sit right after residue 6's LAST atom and before residue 8.
    assert lines[ter_idx[0] - 1][22:26].strip() == "6"
    assert lines[ter_idx[0] + 1][22:26].strip() == "8"


def test_no_ter_when_plan_has_none(tmp_path):
    pdb = tmp_path / "s.pdb"
    _write_pdb(pdb, [("N", "ALA", "A", 5), ("N", "GLY", "A", 6)])
    before = pdb.read_text()

    _module()._apply_ter_records(str(pdb), RepairPlan(), _identity_mapper())

    assert pdb.read_text() == before  # untouched


def test_ter_anchor_translated_through_modeller_numbering(tmp_path):
    pdb = tmp_path / "s.pdb"
    # MODELLER compacted the gap out: original 166,167 stay; 199 -> 168.
    # The file is in MODELLER (final) numbering.
    _write_pdb(pdb, [
        ("N", "ALA", "A", 166), ("N", "GLY", "A", 167),
        ("N", "SER", "A", 168),
    ])
    # Original gap started at 199 (so anchor = original 198). But 198 was part
    # of the gap; the residue actually before it in the file is 167. Use a gap
    # whose preceding residue (start_num-1) maps to a real final number.
    seg = MissingSegment(chain_id="A", residues=[], start_num=199, end_num=210,
                         is_terminal=False, terminal_type=None)
    plan = RepairPlan(segments_to_ter=[seg])

    mapper = ResidueMapper(Console())
    mapper.chain_mapping = {"A": "A"}
    # original -> modeller: 198 -> 167 (residue before the gap)
    mapper.residue_mappings = {"A": {166: 166, 167: 167, 198: 167, 199: 168}}
    mapper.final_mappings = {}

    _module()._apply_ter_records(str(pdb), plan, mapper)

    lines = pdb.read_text().splitlines()
    ter_idx = [i for i, l in enumerate(lines) if l.startswith("TER")]
    assert len(ter_idx) == 1
    assert lines[ter_idx[0] - 1][22:26].strip() == "167"


def test_ter_not_duplicated_if_already_present(tmp_path):
    pdb = tmp_path / "s.pdb"
    with open(pdb, "w") as fh:
        fh.write(_atom(1, "N", "ALA", "A", 5))
        fh.write(_atom(2, "N", "GLY", "A", 6))
        fh.write("TER\n")
        fh.write(_atom(3, "N", "SER", "A", 8))
    seg = MissingSegment(chain_id="A", residues=[], start_num=7, end_num=7,
                         is_terminal=False, terminal_type=None)
    plan = RepairPlan(segments_to_ter=[seg])

    _module()._apply_ter_records(str(pdb), plan, _identity_mapper())

    assert pdb.read_text().count("TER") == 1
