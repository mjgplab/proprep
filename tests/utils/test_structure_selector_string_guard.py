"""
Regression test for StructureSelector.get_structure_object string-guard.

A *_structure workspace key is supposed to hold a parsed BioPython Structure
object. The redox transformation manager, however, stores a *path string* under
transformed_structure (the modAA parameterizer reads it as a path). Before the
guard, get_structure_object mapped transformed_pdb_file -> transformed_structure,
found that string, and returned it unvalidated. Callers then iterated the string
character-by-character and crashed on `chain.id` with:

    'str' object has no attribute 'id'

Reproduced live on a structure with FES/MTE/MOS/FAD redox sites: the FES metal
units classified fine (metal test never touches the structure), then MTE — a
non-metal cofactor — fell through to the peptide-backbone test, which loads the
structure and crashed.

These tests pin the boundary: a string-valued *_structure key is skipped, and a
genuine Structure object is still returned.
"""

from Bio.PDB.StructureBuilder import StructureBuilder

from proprep.utils.structure_selector import StructureSelector
from proprep.forcefield_prep.forcefield_parameterizer import ForcefieldParameterizer


class _FakeWorkspace:
    """Minimal workspace: get_structure_object only ever calls .get()."""

    def __init__(self, data):
        self._data = data

    def get(self, key, default=None):
        return self._data.get(key, default)


def _real_structure():
    sb = StructureBuilder()
    sb.init_structure("t")
    sb.init_model(0)
    sb.init_chain("A")
    sb.init_seg(" ")
    sb.init_residue("ALA", " ", 1, " ")
    sb.init_atom("CA", (0.0, 0.0, 0.0), 0.0, 1.0, " ", "CA", 1, "C")
    return sb.get_structure()


def _selector(data):
    return StructureSelector(_FakeWorkspace(data), processor=None)


def test_string_valued_structure_key_is_skipped():
    # transformed_pdb_file is a valid path candidate; its transformed_structure
    # sibling holds a PATH STRING (the real-world bug). The selector must not
    # hand that string back.
    selector = _selector({
        "transformed_pdb_file": "/nonexistent/transformed.pdb",
        "transformed_structure": "/nonexistent/transformed.pdb",
    })

    result = selector.get_structure_object(
        priority_override=ForcefieldParameterizer._PRIORITY_KEYS,
        silent=True,
    )

    # No genuine Structure object anywhere -> None, and never a str.
    assert result is None
    assert not isinstance(result, str)


def test_real_structure_object_is_still_returned():
    # Guard must not over-reject: a genuine Structure object under a *_structure
    # key is returned unchanged.
    structure = _real_structure()
    selector = _selector({
        "repaired_pdb_file": "/nonexistent/repaired.pdb",
        "repaired_structure": structure,
    })

    result = selector.get_structure_object(
        priority_override=ForcefieldParameterizer._PRIORITY_KEYS,
        silent=True,
    )

    assert result is structure


def test_string_key_skipped_but_lower_priority_object_wins():
    # A string transformed_structure at top priority must not shadow a real
    # object available at lower priority: the selector skips the string and
    # falls through to the genuine Structure.
    structure = _real_structure()
    selector = _selector({
        "transformed_pdb_file": "/nonexistent/transformed.pdb",
        "transformed_structure": "/nonexistent/transformed.pdb",  # string (bug)
        "filtered_pdb_file": "/nonexistent/filtered.pdb",
        "filtered_structure": structure,                          # real object
    })

    result = selector.get_structure_object(
        priority_override=ForcefieldParameterizer._PRIORITY_KEYS,
        silent=True,
    )

    assert result is structure
