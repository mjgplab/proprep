"""
The exported redox-site PDB must put atom names and elements in the right columns.

_export_to_pdb_single / _export_to_pdb_separate built their ATOM records with
two off-by-one column errors:

    f"ATOM  {n:5d}  {name:<4s}..."          two spaces -> name began in col 14
    f"...  1.00 99.00           {el:>2s}"   11 spaces  -> element in cols 78-79

The name one is not cosmetic. A two-letter element must START in column 13; a
one-letter element is right-justified into column 14 so the rest of the field
can hold a remoteness indicator. Writing every name from column 14 emitted
`` FE1``, which resolves to element F, and the element field that would have
corrected the guess sat one column past where readers look. BioPython read the
iron in these files as fluorine and the molybdenum as X.
"""

import warnings

import pytest

warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser  # noqa: E402

from proprep.structure_prep.comprehensive_redox_detector import _pdb_atom_line  # noqa: E402


class _Atom:
    def __init__(self, atom_name, element, resname="LIG", chain="A", resid=1):
        self.atom_name = atom_name
        self.element = element
        self.resname = resname
        self.chain = chain
        self.resid = resid
        self.coords = (-46.078, -17.593, -47.380)


@pytest.mark.parametrize("name,element,expected_name_field", [
    ("FE1", "FE", "FE1 "),    # two-letter element starts in column 13
    ("FE2", "Fe", "FE2 "),
    ("MO",  "MO", "MO  "),
    ("ZN",  "ZN", "ZN  "),
    ("S1",  "S",  " S1 "),    # one-letter element indented to column 14
    ("PA",  "P",  " PA "),    # FAD phosphate, NOT protactinium
    ("CA",  "C",  " CA "),    # alpha carbon
    ("CA",  "CA", "CA  "),    # calcium ion — same name, different column
    ("HH31", "H", "HH31"),    # a 4-character name has nowhere to indent
])
def test_atom_name_column_justification(name, element, expected_name_field):
    line = _pdb_atom_line(1, _Atom(name, element))

    assert line[12:16] == expected_name_field, repr(line)


def test_fixed_columns():
    line = _pdb_atom_line(7, _Atom("FE1", "FE", resname="FES", chain="A", resid=3001))

    assert line[0:6] == "ATOM  "
    assert line[6:11] == "    7"          # serial, 7-11
    assert line[16] == " "                # altLoc, 17
    assert line[17:20] == "FES"           # resName, 18-20
    assert line[21] == "A"                # chainID, 22
    assert line[22:26] == "3001"          # resSeq, 23-26
    assert line[30:38].strip() == "-46.078"
    assert line[54:60] == "  1.00"        # occupancy, 55-60
    assert line[60:66] == " 99.00"        # tempFactor, 61-66
    assert line[76:78] == "Fe"            # element, 77-78
    assert len(line.rstrip("\n")) == 78


def test_biopython_reads_back_every_element(tmp_path):
    """The symptom that started this: iron parsed as fluorine."""
    atoms = [
        _Atom("FE1", "FE", "FES", resid=1),
        _Atom("MO", "MO", "MOS", resid=2),
        _Atom("PA", "P", "FAD", resid=3),
        _Atom("CA", "C", "ALA", resid=4),
        _Atom("CA", "CA", "CA", resid=5),
        _Atom("HH31", "H", "ACE", resid=6),
    ]
    pdb = tmp_path / "sites.pdb"
    pdb.write_text("".join(_pdb_atom_line(i, a) for i, a in enumerate(atoms, 1)) + "END\n")

    parsed = list(PDBParser(QUIET=True).get_structure("x", str(pdb)).get_atoms())

    assert len(parsed) == len(atoms)
    for written, got in zip(atoms, parsed):
        assert got.element.upper() == written.element.upper(), (
            f"{written.atom_name} written as {written.element}, read as {got.element}")


def test_long_fields_cannot_shift_later_columns():
    """A 4-char name and a 3-char resname must not push the element out."""
    line = _pdb_atom_line(99999, _Atom("HH31", "H", resname="NME", resid=9999))

    assert line[17:20] == "NME"
    assert line[76:78] == " H"
    assert len(line.rstrip("\n")) == 78


def test_missing_element_leaves_the_field_blank_not_shifted():
    line = _pdb_atom_line(1, _Atom("C1B", ""))

    assert line[12:16] == " C1B"
    assert line[76:78] == "  "
    assert len(line.rstrip("\n")) == 78
