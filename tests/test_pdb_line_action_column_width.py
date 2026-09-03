#!/usr/bin/env python3
"""
The transformation executor edits PDB lines in place by column. Every field it
touches (chain, residue id, insertion code, atom name) has a fixed width, so an
action must never change the length of the line.

Regression: ``change_insertion_code: ""`` (used by the c-type heme transformers
when they move Cys/His side chains into the heme residue) was spliced in as an
empty string, deleting column 27 and shifting x/y/z left by one. A His claimed
by two sites went through the splice twice and its y field became unparseable.

Run with: pytest tests/test_pdb_line_action_column_width.py
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (  # noqa: E402
    TransformationExecutor,
)

LINE = "ATOM    849  CB  HIS N 143      13.576   0.094 -91.896  1.00 40.32           C  "
FAR_LINE = "ATOM   1201  SG  CYS N 215A    -74.177 -56.938-106.695  1.00 40.32           S  "

MIGRATE = {
    "change_residue_name": "HEC",
    "change_residue_id": 301,
    "change_chain_id": "N",
    "change_insertion_code": "",
    "rename_atoms": {"CB": "CB2"},
    "convert_to_hetatm": True,
}


@pytest.fixture
def executor():
    return TransformationExecutor(verbose=False)


def _fields(line):
    return {
        "name": line[12:16], "resname": line[17:20], "chain": line[21],
        "resid": line[22:26], "icode": line[26], "x": line[30:38],
        "y": line[38:46], "z": line[46:54], "element": line[76:78],
    }


def test_empty_insertion_code_keeps_line_width(executor):
    out = executor._apply_action_to_pdb_line(LINE, MIGRATE)
    assert len(out) == len(LINE)
    f = _fields(out)
    assert out.startswith("HETATM")
    assert f["name"] == " CB2" and f["resname"] == "HEC"
    assert f["chain"] == "N" and f["resid"] == " 301" and f["icode"] == " "
    assert (f["x"], f["y"], f["z"]) == ("  13.576", "   0.094", " -91.896")
    assert f["element"] == " C"


def test_second_pass_over_same_atom_still_parses(executor):
    once = executor._apply_action_to_pdb_line(LINE, MIGRATE)
    twice = executor._apply_action_to_pdb_line(once, MIGRATE)
    assert twice == once
    assert float(twice[38:46]) == pytest.approx(0.094)


def test_clearing_a_real_insertion_code_does_not_drop_a_sign(executor):
    out = executor._apply_action_to_pdb_line(FAR_LINE, MIGRATE)
    assert len(out) == len(FAR_LINE)
    assert out[26] == " "
    assert float(out[46:54]) == pytest.approx(-106.695)


def test_explicit_insertion_code_and_empty_chain_are_one_column(executor):
    out = executor._apply_action_to_pdb_line(
        LINE, {"change_insertion_code": "B", "change_chain_id": ""})
    assert len(out) == len(LINE)
    assert out[26] == "B" and out[21] == " "
    assert float(out[30:38]) == pytest.approx(13.576)
