#!/usr/bin/env python3
"""
ONIOM Preparer: the suggested model-system charge must be a sum of per-fragment
formal charges, not one rounding of pooled partial charges.

Whole AMBER residues sum to exact integers. A side chain cut at CA-CB does not
(ff19SB amino19.lib: ASP -0.858, GLU -0.882, HIP +0.943, LYS +1.026), so four
trimmed carboxylates pool to -3.43, which a single round() turns into -3.

Run with: pytest tests/test_oniom_model_charge.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep.oniom_prep.oniom_qmmm_preparator import suggested_model_charge  # noqa: E402

# Side-chain partial-charge sums from amino19.lib, spread over a few atoms.
SIDECHAIN = {"ASP": -0.858, "GLU": -0.882, "HIP": 0.943, "LYS": 1.026, "SER": 0.055}


def _fragments(resnames):
    """Build a charges dict and one atom-index group per residue."""
    charges, groups, idx = {}, [], 0
    for name in resnames:
        total = SIDECHAIN[name]
        group = []
        for k in range(3):
            charges[idx] = total / 3
            group.append(idx)
            idx += 1
        groups.append(group)
    return charges, groups


def test_four_trimmed_carboxylates_are_minus_four():
    charges, groups = _fragments(["ASP", "ASP", "ASP", "ASP"])
    assert round(sum(charges.values())) == -3          # what pooling gives
    assert suggested_model_charge(charges, [], groups) == -4


def test_mixed_trimmed_side_chains():
    charges, groups = _fragments(["ASP", "GLU", "HIP", "HIP", "SER"])
    assert suggested_model_charge(charges, [], groups) == 0
    charges, groups = _fragments(["ASP", "ASP", "GLU", "LYS"])
    assert suggested_model_charge(charges, [], groups) == -2


def test_whole_residues_and_trimmed_fragments_add():
    # A whole residue with an exact -1 (e.g. a heme propionate) plus two ASP side chains.
    charges = {0: -0.4, 1: -0.6}
    trimmed_charges, groups = _fragments(["ASP", "ASP"])
    charges.update({k + 10: v for k, v in trimmed_charges.items()})
    groups = [[i + 10 for i in g] for g in groups]
    assert suggested_model_charge(charges, [[0, 1]], groups) == -3


def test_missing_atoms_count_as_zero():
    assert suggested_model_charge({}, [[1, 2]], [[3]]) == 0
