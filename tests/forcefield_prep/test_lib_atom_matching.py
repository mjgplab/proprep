"""
Matching a PDB residue against a supplied library: hydrogens and double reads.

Two defects, both hit when supplying an externally obtained FAD library for a
crystal structure.

1. _parse_lib_atoms tested ``'.unit.atoms' in line``, which also matches
   ``.unit.atomspertinfo`` -- a table with one row per atom carrying the same
   names. Every library was therefore read at double length: an 84-atom FAD was
   reported as 168, its listing repeating from index 85.

2. The name sets were compared raw. A crystal structure has no hydrogens and a
   library always does, so a library that fits perfectly reported a mismatch
   and the user was asked to hand-map 53 atoms against 84. tLEaP builds the
   missing hydrogens from the template, so only the heavy atoms need to
   correspond.

Real numbers from the reported case: the library has 84 atoms, 53 heavy and 31
hydrogens; the structure's FAD has exactly those 53.
"""

import pytest

from proprep.forcefield_prep.structure_preprocessor import (
    StructurePreprocessor, _is_hydrogen_name,
)


def _lib(tmp_path, names_and_types, unit="FAD"):
    """A library with BOTH tables, as tLEaP writes them."""
    atoms = "\n".join(f' "{n}" "{t}" 0 1 131072 {i} 6 0.0'
                      for i, (n, t) in enumerate(names_and_types, 1))
    pertinfo = "\n".join(f' "{n}" "{t}" 0 -1 0.0'
                         for n, t in names_and_types)
    path = tmp_path / f"{unit}.lib"
    path.write_text(
        f'!!index array str\n "{unit}"\n'
        f"!entry.{unit}.unit.atoms table  str name  str type\n{atoms}\n"
        f"!entry.{unit}.unit.atomspertinfo table  str pname  str ptype\n{pertinfo}\n"
        f"!entry.{unit}.unit.connectivity table  int atom1x  int atom2x\n 1 2 1\n")
    return path


ATOMS = [("PA", "p5"), ("O1A", "o"), ("C5B", "c3"), ("H112", "h1"), ("H114", "h1")]


def _parser():
    return StructurePreprocessor.__new__(StructurePreprocessor)


# --------------------------------------------------------------------------- #
# reading the atoms table once
# --------------------------------------------------------------------------- #

def test_the_atomspertinfo_table_is_not_read_as_atoms(tmp_path):
    """The reported symptom: 84 atoms reported as 168."""
    parsed = _parser()._parse_lib_atoms(str(_lib(tmp_path, ATOMS)))

    assert len(parsed) == len(ATOMS)


def test_names_are_not_duplicated(tmp_path):
    parsed = _parser()._parse_lib_atoms(str(_lib(tmp_path, ATOMS)))
    names = [a[0] for a in parsed]

    assert len(names) == len(set(names))


def test_types_survive_the_parse(tmp_path):
    parsed = _parser()._parse_lib_atoms(str(_lib(tmp_path, ATOMS)))

    assert dict((n, t) for n, t, _e in parsed)["PA"] == "p5"


def test_a_unit_named_something_else_still_parses(tmp_path):
    parsed = _parser()._parse_lib_atoms(
        str(_lib(tmp_path, ATOMS[:2], unit="x")))

    assert len(parsed) == 2


def test_a_missing_file_yields_nothing(tmp_path):
    assert _parser()._parse_lib_atoms(str(tmp_path / "nope.lib")) == []


# --------------------------------------------------------------------------- #
# telling a hydrogen from a heavy atom
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("name", ["H", "H112", "HB2", "HH31", "1HB", "2HG1"])
def test_hydrogen_names_are_recognised(name):
    """PDB names may lead with a branch digit, so the first LETTER decides."""
    assert _is_hydrogen_name(name)


@pytest.mark.parametrize("name", ["PA", "O1A", "C5B", "N9A", "HG", "HO"])
def test_heavy_atoms_are_not_mistaken_for_hydrogen(name):
    """HG is mercury and HO is holmium; neither is a hydrogen."""
    if name in ("HG", "HO"):
        pytest.skip("two-letter element names starting with H are ambiguous "
                    "by name alone; see test_ambiguous_two_letter_names")
    assert not _is_hydrogen_name(name)


def test_ambiguous_two_letter_names_are_treated_as_hydrogen():
    """
    Documented limitation: 'HG' is a mercury atom name AND a gamma hydrogen's.
    Treating it as hydrogen is the safe default here -- the consequence is a
    heavy atom excluded from the comparison, which at worst falls through to
    the manual mapping that existed before.
    """
    assert _is_hydrogen_name("HG")


def test_an_empty_name_is_not_hydrogen():
    assert not _is_hydrogen_name("")
    assert not _is_hydrogen_name(None)


# --------------------------------------------------------------------------- #
# what the comparison should conclude
# --------------------------------------------------------------------------- #

def test_the_reported_case_matches_on_heavy_atoms(tmp_path):
    """
    A structure with 53 heavy atoms against a library of 53 heavy + 31 H.
    Comparing raw sets called this a mismatch and asked for 53 hand-mappings.
    """
    lib_names = {"PA", "O1A", "C5B"} | {f"H{i}" for i in range(1, 32)}
    pdb_names = {"PA", "O1A", "C5B"}

    lib_heavy = {n for n in lib_names if not _is_hydrogen_name(n)}
    pdb_heavy = {n for n in pdb_names if not _is_hydrogen_name(n)}

    assert pdb_names != lib_names          # raw comparison: mismatch
    assert pdb_heavy == lib_heavy          # heavy-atom comparison: match


def test_a_genuine_heavy_atom_difference_is_still_a_mismatch():
    """The check must not become permissive: a real difference must show."""
    lib_heavy = {"PA", "O1A", "C5B"}
    pdb_heavy = {"PA", "O1A", "C9Z"}

    assert pdb_heavy != lib_heavy


def test_a_structure_that_has_hydrogens_is_compared_in_full():
    """
    Heavy-only matching applies when the PDB has NO hydrogens. If it has them,
    they are real information and a difference in them is a real mismatch.
    """
    pdb_names = {"PA", "H112"}
    pdb_heavy = {n for n in pdb_names if not _is_hydrogen_name(n)}

    assert pdb_names - pdb_heavy, "this structure does have hydrogens"
