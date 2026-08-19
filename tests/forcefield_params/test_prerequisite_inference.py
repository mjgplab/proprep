"""
What a deposited parameter set requires, read from the types it uses.

A set is only loadable if the force fields defining its atom types are sourced.
metadata declares that as prerequisites.leaprc_groups and the Topology
Generator enforces what it finds -- but almost nothing wrote it. Of nineteen
library entries only four declared anything: the two built-in flavins, zinc/cys4
and two hemes. Every MCPB-generated metal site declared nothing, and so did an
imported set, so the requirement existed only in the depositor's head.

The types are in the files at deposit time, so it can be read rather than
guessed: a GAFF set uses c3/os/p5, an MCPB set uses N/CT/XC alongside its
M*/Y*. This is also what replaces the old prerequisites panel, which asserted
that every cofactor needs a protein force field and pattern-matched residue
NAMES to decide GAFF2 was needed -- while flavin/fad declared exactly that in
metadata, unread.
"""

import pytest

from proprep.forcefield_params.prerequisites import (
    GAFF2_LEAPRC, PROTEIN_LEAPRCS, atom_types_used, infer_leaprc_groups,
    type_sources,
)


available = pytest.mark.skipif(
    not any(type_sources()), reason="Amber parameter files not locatable")


def _lib(tmp_path, name, atom_types):
    path = tmp_path / f"{name}.lib"
    rows = "\n".join(f' "A{i}" "{t}" 0 1 131072 {i} 6 0.0'
                     for i, t in enumerate(atom_types, 1))
    path.write_text(
        f"!!index array str\n"
        f"!entry.{name}.unit.atoms table  str name  str type\n{rows}\n"
        f"!entry.{name}.unit.atomspertinfo table\n")
    return path


def _mol2(tmp_path, name, atom_types):
    path = tmp_path / f"{name}.mol2"
    rows = "\n".join(f"  {i} A{i} 0.0 0.0 0.0 {t} 1 LIG 0.0"
                     for i, t in enumerate(atom_types, 1))
    path.write_text(f"@<TRIPOS>MOLECULE\n{name}\n@<TRIPOS>ATOM\n{rows}\n@<TRIPOS>BOND\n")
    return path


# --------------------------------------------------------------------------- #
# reading the types a set uses
# --------------------------------------------------------------------------- #

def test_types_are_read_from_a_library(tmp_path):
    lib = _lib(tmp_path, "FAD", ["c3", "os", "p5"])

    assert atom_types_used([lib]) == {"c3", "os", "p5"}


def test_types_are_read_from_a_mol2(tmp_path):
    mol2 = _mol2(tmp_path, "LIG", ["ca", "ha"])

    assert atom_types_used([mol2]) == {"ca", "ha"}


def test_a_missing_file_contributes_nothing(tmp_path):
    assert atom_types_used([tmp_path / "nope.lib"]) == set()


# --------------------------------------------------------------------------- #
# what those types imply
# --------------------------------------------------------------------------- #

@available
def test_a_gaff_set_requires_gaff2(tmp_path):
    """The reported FAD import: all GAFF2 types, empty MASS section."""
    lib = _lib(tmp_path, "FAD", ["c", "c3", "ca", "os", "p5", "ho"])

    assert infer_leaprc_groups([lib]) == [[GAFF2_LEAPRC]]


@available
def test_an_mcpb_set_requires_a_protein_forcefield(tmp_path):
    """Its coordinating residues keep N/CT/XC alongside the new M*/Y*."""
    lib = _lib(tmp_path, "CM1", ["N", "H", "XC", "CT", "Y1", "C", "O"])

    assert infer_leaprc_groups([lib]) == [list(PROTEIN_LEAPRCS)]


@available
def test_a_set_using_both_requires_both(tmp_path):
    """
    4hux_moco: HO on the hydroxo plus the pterin's GAFF types. Its tleap input
    loads both, and nothing declared either.
    """
    lib = _lib(tmp_path, "MS1", ["HO", "M3", "YA", "c6", "cd"])

    groups = infer_leaprc_groups([lib])

    assert [GAFF2_LEAPRC] in groups
    assert list(PROTEIN_LEAPRCS) in groups
    assert len(groups) == 2


@available
def test_custom_types_alone_imply_nothing(tmp_path):
    """M*/Y* are defined by the set's own addAtomTypes, not by a leaprc."""
    lib = _lib(tmp_path, "FS1", ["M1", "M2", "Y5", "Y6"])

    assert infer_leaprc_groups([lib]) == []


@available
def test_a_type_defined_by_both_does_not_vote(tmp_path):
    """
    Types shared between GAFF and the protein sets discriminate nothing;
    counting them would make every set claim to need everything.
    """
    gaff_only, protein_only = type_sources()

    assert not (gaff_only & protein_only)


def test_no_files_means_no_declaration(tmp_path):
    """Empty is the schema's "nothing declared", which consumers handle."""
    assert infer_leaprc_groups([]) == []


# --------------------------------------------------------------------------- #
# it agrees with the entries that declared prerequisites by hand
# --------------------------------------------------------------------------- #

@available
def test_it_reproduces_the_hand_written_flavin_declaration(tmp_path):
    """flavin/fad declares [['leaprc.gaff2']]; inference must not contradict it."""
    lib = _lib(tmp_path, "FAD", ["c3", "os", "oh", "ca", "nd"])

    assert infer_leaprc_groups([lib]) == [[GAFF2_LEAPRC]]


@available
def test_the_protein_group_matches_the_hand_written_one():
    """
    zinc/cys4 declares ff14SB / ff19SB / ff14SBonlysc / constph / conste as
    alternatives. Inference offers the same set, so a picker satisfying one
    satisfies the other.
    """
    from proprep.forcefield_params.loader import get_prerequisite_leaprc_groups

    try:
        declared = get_prerequisite_leaprc_groups("zinc/cys4")
    except Exception:
        pytest.skip("zinc/cys4 not available")

    assert declared, "zinc/cys4 should declare a protein-FF group"
    assert set(declared[0]) == set(PROTEIN_LEAPRCS)
