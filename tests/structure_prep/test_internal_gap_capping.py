"""
Regression test: capping (not filling) an INTERNAL gap must not crash.

The internal-gap branch of CappingHandler._build_cap_insertion_plan referenced
undefined names (first_modeller / last_modeller) and an uninitialized
following_res, so choosing 'c' (cap) on an internal gap raised NameError. The
fix translates the gap-bracketing residues through the ResidueMapper (MODELLER
path) or scans the structure directly (caps-only/identity path).

An internal gap chosen for capping splits the chain: NME after the residue
preceding the gap, ACE before the residue following it.
"""

from rich.console import Console

import proprep.structure_prep.structure_completeness as sc
from proprep.structure_prep.structure_completeness import (
    ResidueMapper, RepairPlan, MissingSegment,
)


def _capping_handler():
    for obj in vars(sc).values():
        if isinstance(obj, type) and "_build_cap_insertion_plan" in obj.__dict__:
            inst = obj.__new__(obj)
            inst.console = Console()
            return inst
    raise AssertionError("could not find the class owning _build_cap_insertion_plan")


class _Res:
    def __init__(self, num):
        self.id = (" ", num, " ")


class _Chain:
    def __init__(self, nums):
        self._r = [_Res(n) for n in nums]

    def __iter__(self):
        return iter(self._r)


class _Model:
    def __init__(self, chains):
        self._c = chains

    def __contains__(self, k):
        return k in self._c

    def __getitem__(self, k):
        return self._c[k]


def _internal_gap_plan():
    seg = MissingSegment(chain_id="A", residues=[], start_num=168, end_num=198,
                         is_terminal=False, terminal_type=None)
    return RepairPlan(segments_to_cap=[seg])


def test_internal_gap_cap_modeller_path():
    """MODELLER ran: bracket residues come from the mapper (compacted numbers)."""
    inst = _capping_handler()
    mapper = ResidueMapper(Console())
    mapper.chain_mapping = {"A": "A"}
    # original -> modeller; gap 168-198 absent, so 199 compacts to 168
    mapper.residue_mappings = {"A": {166: 166, 167: 167, 199: 168, 200: 169}}

    out = inst._build_cap_insertion_plan(_internal_gap_plan(), mapper,
                                         repaired_structure=[])
    ins = out["A"]
    nme = next(d for d in ins if d["type"] == "NME")
    ace = next(d for d in ins if d["type"] == "ACE")
    assert nme["insert_after"] == 167     # residue before the gap
    assert ace["insert_before"] == 168    # residue after the gap (modeller-numbered)


def test_internal_gap_cap_caps_only_path():
    """No MODELLER (identity map): scan the structure in original numbering."""
    inst = _capping_handler()
    mapper = ResidueMapper(Console())
    mapper.chain_mapping = {}
    mapper.residue_mappings = {}
    repaired = [_Model({"A": _Chain([166, 167, 199, 200])})]

    out = inst._build_cap_insertion_plan(_internal_gap_plan(), mapper,
                                         repaired_structure=repaired)
    ins = out["A"]
    nme = next(d for d in ins if d["type"] == "NME")
    ace = next(d for d in ins if d["type"] == "ACE")
    assert nme["insert_after"] == 167
    assert ace["insert_before"] == 199
