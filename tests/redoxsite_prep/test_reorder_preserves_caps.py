"""
Internal caps must stay in the chain when the PDB is reordered for tLEaP.

_reorder_for_tleap binned every ACE into one list and every NME into another,
then wrote `all ACE / all protein / all NME` per chain. That assumes a cap can
only ever be terminal.

An internal cap breaks it. Capping an unfilled gap puts an NME after the last
residue before the gap and an ACE before the first residue after it, mid-chain.
Binning moved that pair to opposite ends: on a 4UHX run the ACE landed before
residue 1 and tLEaP bonded it to the N-terminus 70 A away, while the gap it had
been guarding was left open for a 20.8 A peptide bond. Both were fully
parameterized -- a C-N bond has parameters whatever the distance -- so the
build "succeeded".

    ACE1  C -> MET2  N :  70.13 A
    THR168 C -> PRO169 N : 20.79 A

Caps are the only residues that are part of the peptide chain while recorded as
HETATM, which is why a reorder keyed on record type mishandles them.
"""

import math

import pytest

from proprep.redoxsite_prep.transformation.redox_transformation_manager import (
    RedoxTransformationManager,
)


def _atom(record, serial, name, resname, chain, resid, xyz):
    x, y, z = xyz
    return (f"{record:<6s}{serial:>5d} {name:^4s} {resname:>3s} {chain}{resid:>4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n")


def _capped_gap_structure():
    """A chain broken by a capped gap: ...166 167 [NME 168][ACE 169] 170 171...

    The two fragments are far apart, as an unfilled gap's flanks are.
    """
    lines = []
    serial = 1
    for resid, x in ((166, 0.0), (167, 3.8)):
        for name, dx in (("N", 0.0), ("CA", 1.4), ("C", 2.4), ("O", 2.9)):
            lines.append(_atom("ATOM", serial, name, "ALA", "A", resid, (x + dx, 0.0, 0.0)))
            serial += 1
    lines.append(_atom("HETATM", serial, "N", "NME", "A", 168, (7.5, 0.0, 0.0))); serial += 1
    lines.append(_atom("HETATM", serial, "CH3", "NME", "A", 168, (8.5, 0.0, 0.0))); serial += 1
    # The far fragment, 30 A away.
    lines.append(_atom("HETATM", serial, "CH3", "ACE", "A", 169, (30.0, 0.0, 0.0))); serial += 1
    lines.append(_atom("HETATM", serial, "C", "ACE", "A", 169, (31.0, 0.0, 0.0))); serial += 1
    lines.append(_atom("HETATM", serial, "O", "ACE", "A", 169, (31.4, 1.0, 0.0))); serial += 1
    for resid, x in ((170, 32.3), (171, 36.1)):
        for name, dx in (("N", 0.0), ("CA", 1.4), ("C", 2.4), ("O", 2.9)):
            lines.append(_atom("ATOM", serial, name, "ALA", "A", resid, (x + dx, 0.0, 0.0)))
            serial += 1
    lines.append(_atom("HETATM", serial, "FE1", "FES", "A", 300, (50.0, 0.0, 0.0)))
    lines.append("END\n")
    return lines


def _reorder(lines):
    manager = RedoxTransformationManager.__new__(RedoxTransformationManager)
    return manager._reorder_for_tleap(lines)


def _residue_order(lines):
    order, seen = [], set()
    for line in lines:
        if line[:6] in ("ATOM  ", "HETATM"):
            key = (int(line[22:26]), line[17:20].strip())
            if key not in seen:
                seen.add(key)
                order.append(key)
    return order


# --------------------------------------------------------------------------- #
# the caps stay put
# --------------------------------------------------------------------------- #

def test_an_internal_cap_pair_keeps_its_position():
    order = _residue_order(_reorder(_capped_gap_structure()))

    assert order[:6] == [(166, "ALA"), (167, "ALA"), (168, "NME"),
                         (169, "ACE"), (170, "ALA"), (171, "ALA")]


def test_the_ace_does_not_migrate_to_the_chain_start():
    """The reported failure: ACE at residue 1, bonded to the N-terminus."""
    order = _residue_order(_reorder(_capped_gap_structure()))

    assert order[0] != (169, "ACE")


def test_the_nme_does_not_migrate_to_the_chain_end():
    order = _residue_order(_reorder(_capped_gap_structure()))
    protein = [r for r in order if r[1] != "FES"]

    assert protein[-1] != (168, "NME")


def test_a_real_hetatm_still_moves_after_the_protein():
    order = _residue_order(_reorder(_capped_gap_structure()))

    assert order[-1] == (300, "FES")


# --------------------------------------------------------------------------- #
# the break is closed off
# --------------------------------------------------------------------------- #

def test_a_ter_follows_the_nme():
    """Without it tLEaP bonds straight across the gap the caps bracket."""
    out = _reorder(_capped_gap_structure())

    nme_last = max(i for i, l in enumerate(out)
                   if l[:6] == "HETATM" and l[17:20].strip() == "NME")

    assert out[nme_last + 1].startswith("TER")


def test_no_bondable_gap_survives():
    """
    The property that matters: no C->N pair left adjacent across a break.
    tLEaP only refuses when it cannot find PARAMETERS, and a C-N bond has them
    at any distance -- so a 70 A peptide bond builds silently.
    """
    out = _reorder(_capped_gap_structure())

    coords, prev, spans = {}, None, []
    for line in out:
        if line.startswith("TER"):
            prev = None
            continue
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        key = (int(line[22:26]), line[17:20].strip())
        coords.setdefault(key, {})[line[12:16].strip()] = (
            float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if prev and prev != key:
            c, n = coords.get(prev, {}).get("C"), coords.get(key, {}).get("N")
            if c and n:
                spans.append((prev, key, math.dist(c, n)))
        prev = key

    assert not [s for s in spans if s[2] > 2.0], f"bondable gap: {spans}"


def test_a_terminal_cap_is_unaffected():
    """The case the old binning handled; it must keep working."""
    lines = [
        _atom("HETATM", 1, "CH3", "ACE", "A", 0, (0.0, 0.0, 0.0)),
        _atom("HETATM", 2, "C", "ACE", "A", 0, (1.0, 0.0, 0.0)),
        _atom("ATOM", 3, "N", "ALA", "A", 1, (2.4, 0.0, 0.0)),
        _atom("ATOM", 4, "C", "ALA", "A", 1, (3.4, 0.0, 0.0)),
        _atom("HETATM", 5, "N", "NME", "A", 2, (4.8, 0.0, 0.0)),
        "END\n",
    ]

    assert _residue_order(_reorder(lines)) == [(0, "ACE"), (1, "ALA"), (2, "NME")]


def test_existing_ter_lines_are_rebuilt_not_duplicated():
    lines = _capped_gap_structure()
    lines.insert(8, "TER\n")

    out = _reorder(lines)

    # One after the NME, one closing the chain.
    assert sum(1 for l in out if l.startswith("TER")) == 2


def test_chains_are_kept_separate():
    lines = [
        _atom("ATOM", 1, "N", "ALA", "A", 1, (0.0, 0.0, 0.0)),
        _atom("ATOM", 2, "N", "ALA", "B", 1, (20.0, 0.0, 0.0)),
        "END\n",
    ]

    out = _reorder(lines)

    assert sum(1 for l in out if l.startswith("TER")) == 2
