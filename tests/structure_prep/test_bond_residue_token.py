"""The bond-definition prompt accepts residue IDs (bare resid, resname+resid,
chain:resid) so the user can read IDs off the viewer instead of translating to
table rows, while still accepting the old 1-based row numbers for replaying
older session logs.
"""

import pytest

from proprep.structure_prep.comprehensive_redox_detector import (
    resolve_bond_residue_token as resolve,
)


# (chain, resname, resid, insertion_code) — mirrors the bond-definition table.
SITE = [
    ("A", "E4Z", 201, ""), ("A", "MN", 202, ""), ("A", "MN", 203, ""),
    ("A", "HIS", 41, ""), ("A", "ASP", 108, ""), ("A", "ILE", 120, ""),
    ("A", "GLU", 119, ""), ("A", "GLU", 80, ""), ("A", "HOH", 301, ""),
    ("A", "HOH", 302, ""),
]


@pytest.mark.parametrize("token,expected", [
    ("202", 1),          # bare resid
    ("80", 7),
    ("MN202", 1),        # resname+resid
    ("GLU80", 7),
    ("E4Z201", 0),       # resname containing a digit
    ("A:203", 2),        # chain:resid
    ("a:203", 2),        # case-insensitive chain
    ("mn202", 1),        # case-insensitive resname
])
def test_residue_id_forms(token, expected):
    assert resolve(token, SITE) == expected


@pytest.mark.parametrize("token,expected", [
    ("1", 0),  # row 1 (no residue has resid 1)
    ("2", 1),  # row 2
    ("10", 9),
])
def test_row_number_fallback(token, expected):
    # Backward compatibility for old session logs that recorded table rows.
    assert resolve(token, SITE) == expected


def test_residue_id_takes_precedence_over_row():
    # A site where a residue's resid collides with a valid row number: the
    # residue ID wins (new behaviour). Resid 2 exists at index 2.
    site = [("A", "GLY", 1, ""), ("A", "CYS", 2, ""), ("A", "ZN", 3, "")]
    assert resolve("2", site) == 1  # resid 2 -> index 1, not row 2 -> index 1
    assert resolve("3", site) == 2  # resid 3 -> index 2, not row 3 (out of range anyway)


def test_no_match_and_bad_input_report_error():
    errs = []
    assert resolve("MN80", SITE, on_error=errs.append) is None   # resname/resid mismatch
    assert resolve("999", SITE, on_error=errs.append) is None    # no resid, not a row
    assert resolve("abc", SITE, on_error=errs.append) is None    # unparseable
    assert resolve("", SITE, on_error=errs.append) is None       # empty
    assert len(errs) == 3  # empty token returns silently


def test_ambiguous_resid_across_chains_needs_chain():
    site = [("A", "HIS", 50, ""), ("B", "HIS", 50, "")]
    errs = []
    assert resolve("50", site, on_error=errs.append) is None
    assert "ambiguous" in errs[0].lower()
    # Qualifying with chain resolves it.
    assert resolve("A:50", site) == 0
    assert resolve("B:50", site) == 1
