"""
Metal detection must not read a PDB atom NAME as an element symbol.

Triage classified FAD as an organometallic cofactor. FAD contains no metal: its
two phosphate atoms are named PA and PB, and `_has_embedded_metal` consulted the
atom name even when the element field was populated, so "PA".title() == "Pa"
matched protactinium and "PB".title() == "Pb" matched lead.

The name fallback exists for PDBs with a blank element column, so it is kept —
but gated on the field actually being blank, and taught the PDB column
convention: a two-letter element starts in column 13 ("FE1 " is iron), while a
one-letter element is right-justified into column 14 (" PA " is phosphorus,
" CA " an alpha carbon, "CA  " a calcium ion).
"""

import warnings

import pytest

warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser  # noqa: E402

from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor  # noqa: E402


def _residue(tmp_path, lines, name="t.pdb"):
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\nEND\n")
    s = PDBParser(QUIET=True).get_structure("x", str(p))
    return next(s.get_residues())


def _sp():
    return StructurePreprocessor.__new__(StructurePreprocessor)


# Columns:      1-6    7-11 13-16 18-20 22 23-26                            77-78
FAD_WITH_ELEMENT = [
    "HETATM10001  PA  FAD A3006     -64.995 -38.702 -39.942  1.00 99.00           P",
    "HETATM10002  PB  FAD A3006     -63.495 -38.633 -39.779  1.00 99.00           P",
    "HETATM10003  O1A FAD A3006     -65.603 -38.813 -41.320  1.00 99.00           O",
    "HETATM10004  C5B FAD A3006     -65.553 -39.902 -39.021  1.00 99.00           C",
]

# Same residue with the element column stripped, as older PDBs have it.
FAD_NO_ELEMENT = [line[:76] for line in FAD_WITH_ELEMENT]

HEME_NO_ELEMENT = [
    "HETATM 5001 FE   HEM A 400      10.000  10.000  10.000  1.00 20.00",
    "HETATM 5002  NA  HEM A 400      11.000  10.000  10.000  1.00 20.00",
    "HETATM 5003  C1A HEM A 400      12.000  10.000  10.000  1.00 20.00",
]


def test_fad_is_not_organometallic(tmp_path):
    """The reported bug: PA/PB read as protactinium/lead."""
    res = _residue(tmp_path, FAD_WITH_ELEMENT)

    assert _sp()._has_embedded_metal(res) is False


def test_fad_is_not_organometallic_without_an_element_column(tmp_path):
    """Even on the fallback path, the column convention says PA is phosphorus."""
    res = _residue(tmp_path, FAD_NO_ELEMENT)

    assert _sp()._has_embedded_metal(res) is False


def test_heme_iron_is_still_found_without_an_element_column(tmp_path):
    """The fallback exists for this case and has to keep working.

    Note NA here is a porphyrin nitrogen, not sodium — the same convention
    that saves FAD saves this too.
    """
    res = _residue(tmp_path, HEME_NO_ELEMENT)

    assert _sp()._has_embedded_metal(res) is True


def test_fe2s2_still_reads_as_a_metal_cluster(tmp_path):
    res = _residue(tmp_path, [
        "HETATM10001 FE1  FES A3001     -46.078 -17.593 -47.380  1.00 52.71          FE",
        "HETATM10002 FE2  FES A3001     -43.100 -18.216 -47.517  1.00 49.26          FE",
        "HETATM10003  S1  FES A3001     -44.425 -16.836 -48.664  1.00 49.85           S",
        "HETATM10004  S2  FES A3001     -44.700 -19.000 -46.400  1.00 49.85           S",
    ])
    sp = _sp()

    assert sp._has_embedded_metal(res) is True
    assert sp._is_pure_metal_cluster(res) is True


def test_element_field_wins_over_the_name(tmp_path):
    """A stated non-metal element is not second-guessed by the name."""
    res = _residue(tmp_path, [
        "HETATM 7001  PA  XYZ A 500      10.000  10.000  10.000  1.00 20.00           P",
    ])

    # Single atom: goes through _is_isolated_metal.
    assert _sp()._is_isolated_metal(res) is False


def test_calcium_ion_is_still_a_metal(tmp_path):
    """CA in column 13 is calcium; the convention has to cut both ways."""
    res = _residue(tmp_path, [
        "HETATM 8001 CA    CA A 600      10.000  10.000  10.000  1.00 20.00",
    ])

    assert _sp()._is_isolated_metal(res) is True


@pytest.mark.parametrize("fullname,expected", [
    (" PA ", "P"),    # FAD phosphate — phosphorus, not protactinium
    (" PB ", "P"),    # FAD phosphate — phosphorus, not lead
    (" CA ", "C"),    # alpha carbon, not calcium
    ("CA  ", "Ca"),   # calcium ion
    ("FE1 ", "Fe"),   # cluster iron
    ("MO  ", "Mo"),   # molybdenum
    (" N  ", "N"),
    (" C1A", "C"),
])
def test_column_convention(tmp_path, fullname, expected):
    class _Atom:
        def __init__(self, full):
            self.fullname = full
            self.name = full.strip()

    assert StructurePreprocessor._element_from_atom_name(_Atom(fullname)) == expected
