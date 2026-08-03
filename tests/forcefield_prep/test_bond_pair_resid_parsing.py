"""
The bond-definition prompt accepts RESIDUE IDs (the numbers shown in the table
and 3D viewer), not just table row indices — so the user can type what they read
off the viewer without converting to a row number.

_parse_bond_residue_pairs interprets each token as a residue ID first, falling
back to a 1-based table index only when the token is not a residue ID in the
site. Post-tLEaP residues are globally renumbered and unique, so residue IDs
resolve unambiguously.
"""

from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor


class _Console:
    def print(self, *a, **k):
        pass


def _p():
    obj = object.__new__(StructurePreprocessor)
    obj.console = _Console()
    return obj


# 7K0W-like site: (chain, resname, resid, icode). Blank chains, unique resids.
SITE = [
    ("", "HID", 45, ""),   # index 0
    ("", "GLU", 65, ""),   # index 1
    ("", "ASP", 93, ""),   # index 2
    ("", "E4Z", 182, ""),  # index 3
    ("", "MN", 185, ""),   # index 4
    ("", "MN", 186, ""),   # index 5
]


def test_residue_ids_map_to_indices():
    # "185-45 185-182" == bond MN185 to HID45 and to E4Z182
    result = _p()._parse_bond_residue_pairs("185-45 185-182", SITE)
    # source index 4 (resid 185) -> targets index 0 (45), index 3 (182)
    assert result == {4: [0, 3]}


def test_grouping_by_source_residue_id():
    result = _p()._parse_bond_residue_pairs("185-45 185-65 186-93", SITE)
    assert result == {4: [0, 1], 5: [2]}


def test_comma_and_space_separators():
    assert _p()._parse_bond_residue_pairs("185-45, 186-93", SITE) == {4: [0], 5: [2]}


def test_unknown_residue_id_is_rejected():
    # 999 is neither a residue ID nor a valid index -> pair dropped.
    result = _p()._parse_bond_residue_pairs("999-45", SITE)
    assert result == {}


def test_low_number_prefers_residue_id_over_index():
    # A site whose residue IDs collide with the 1..N index range: resid wins.
    site = [
        ("", "MN", 5, ""),   # index 0, resid 5
        ("", "HID", 2, ""),  # index 1, resid 2
        ("", "GLU", 9, ""),  # index 2, resid 9
    ]
    # "5-2" should mean resid 5 (index 0) -> resid 2 (index 1), NOT index 5/2.
    result = _p()._parse_bond_residue_pairs("5-2", site)
    assert result == {0: [1]}


def test_index_fallback_when_token_is_not_a_resid():
    # Tokens that are not residue IDs but are valid table indices still work
    # (legacy '1-2' behavior). Here resids are 45.. so 1 and 2 are indices.
    result = _p()._parse_bond_residue_pairs("1-2", SITE)
    assert result == {0: [1]}
