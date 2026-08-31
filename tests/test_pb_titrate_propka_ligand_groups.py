"""parse_propka_groups must recover ligand groups and must not confuse a
C-terminal residue's terminus group with its side-chain group.

ProPKA's .pka file identifies ligand groups by residue name plus ATOM name and
never writes their residue number:

    PRN  CG A     4.54       4.50                OCO

so parse_propka_summary drops every heme propionate. On a multi-heme cytochrome
that is 48 real PB sites. The Python API keeps the owning atom, so residue
identity is recoverable there.

The collision this guards against: a C-terminal Lys carries BOTH a "LYS"
side-chain group and a "COO" C-terminus group on the same residue. Keying on
group type alone would overwrite the side-chain pKa (~11.5) with the terminus
pKa (~2). Observed on all six C-terminal Lys of the 9YUQ filament.
"""
import pytest

from proprep.pb_titrate.propka_compare import (
    _GROUP_TYPE_RESNAMES, parse_propka_groups)


class _Atom:
    def __init__(self, res_name, res_num, name):
        self.res_name, self.res_num, self.name = res_name, res_num, name


class _Group:
    def __init__(self, gtype, res_name, res_num, atom_name, pka, used=True):
        self.type, self.pka_value, self._used = gtype, pka, used
        self.atom = _Atom(res_name, res_num, atom_name)

    def use_in_calculations(self):
        return self._used


class _Conf:
    def __init__(self, groups):
        self.groups = groups


class _Mol:
    def __init__(self, groups):
        self.conformations = {"1A": _Conf(groups), "AVR": _Conf(groups)}


GROUPS = [
    _Group("COO", "ASP", 122, "CG", 8.47),
    _Group("COO", "GLU", 419, "CD", 7.53),
    _Group("LYS", "LYS", 1150, "NZ", 11.83),
    _Group("TYR", "TYR", 55, "OH", 16.97),
    _Group("HIS", "HIS", 32, "ND1", 6.41),
    _Group("OCO", "PRN", 1201, "CG", 7.26),      # heme propionate
    _Group("OCO", "PRN", 1213, "CG", 3.21),
    # C-terminal Lys: side chain AND terminus, same residue
    _Group("LYS", "LYS", 198, "NZ", 10.53),
    _Group("COO", "LYS", 198, "C", 2.10),
    # not PB sites
    _Group("ARG", "ARG", 45, "CZ", 12.5),
    _Group("N+", "VAL", 1, "N", 7.87),
    _Group("NAR", "HCO", 1196, "NA", -5.18),     # Fe-coordinated pyrrole N
    # excluded by ProPKA itself
    _Group("COO", "ASP", 999, "CG", 3.0, used=False),
]


@pytest.fixture
def parsed(monkeypatch):
    mod = type("m", (), {"single": staticmethod(lambda *a, **k: _Mol(GROUPS))})
    monkeypatch.setitem(__import__("sys").modules, "propka.run", mod)
    return parse_propka_groups("ignored.pdb")


def test_ligand_carboxylates_recovered(parsed):
    """The bug: heme propionates were absent entirely."""
    assert parsed[("PRN", 1201)] == 7.26
    assert parsed[("PRN", 1213)] == 3.21


def test_standard_residues_kept(parsed):
    assert parsed[("ASP", 122)] == 8.47
    assert parsed[("GLU", 419)] == 7.53
    assert parsed[("LYS", 1150)] == 11.83
    assert parsed[("TYR", 55)] == 16.97
    assert parsed[("HIS", 32)] == 6.41


def test_c_terminus_does_not_overwrite_lysine_side_chain(parsed):
    """COO on a LYS residue is the terminus, not a side chain."""
    assert parsed[("LYS", 198)] == 10.53


def test_non_pb_group_types_skipped(parsed):
    for key in (("ARG", 45), ("VAL", 1), ("HCO", 1196)):
        assert key not in parsed
    for gtype in ("ARG", "N+", "NAR"):
        assert gtype not in _GROUP_TYPE_RESNAMES


def test_groups_propka_excluded_are_skipped(parsed):
    assert ("ASP", 999) not in parsed


def test_missing_propka_returns_empty(monkeypatch):
    """No propka installed must degrade to {}, not raise."""
    import builtins
    real = builtins.__import__

    def fake(name, *a, **k):
        if name == "propka.run":
            raise ImportError("no propka")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake)
    assert parse_propka_groups("ignored.pdb") == {}
