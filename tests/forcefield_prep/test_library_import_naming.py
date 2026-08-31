"""Import-wizard naming and the atom-type check.

Two silent wrongs: the residue name was taken from the .lib FILENAME rather
than the unit inside it (a GDP.lib holding unit "gdp" recorded "GDP", which no
structure carries), and when Amber's parameter files could not be found every
declared atom type was listed as new with nothing saying why.
"""

from unittest.mock import patch

from rich.console import Console

from proprep.forcefield_prep import library_promotion as lp
from proprep.forcefield_prep.library_promotion import (
    _prompt_imported_atom_types, library_unit_names,
)

LIB = '''!!index array str
 "gdp"
!entry.gdp.unit.atoms table  str name  str type  int typex
 "PB" "P" 0
'''

TWO_UNIT_LIB = '''!!index array str
 "AAA"
 "BBB"
!entry.AAA.unit.atoms table  str name  str type  int typex
'''

FRCMOD = """polyphosphate mods
MASS
CT 12.01
O3 16.00

BOND
"""


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


# --------------------------------------------------------------------------
# the unit name, not the filename
# --------------------------------------------------------------------------

def test_the_unit_name_is_read_from_the_file_not_the_filename(tmp_path):
    assert library_unit_names(_write(tmp_path, "GDP.lib", LIB)) == ["gdp"]


def test_every_unit_is_returned_in_file_order(tmp_path):
    assert library_unit_names(_write(tmp_path, "two.lib", TWO_UNIT_LIB)) == ["AAA", "BBB"]


def test_an_unreadable_library_yields_nothing(tmp_path):
    assert library_unit_names(tmp_path / "does_not_exist.lib") == []


def test_parsing_stops_at_the_first_entry_block(tmp_path):
    """Atom names must not be mistaken for unit names."""
    assert "PB" not in library_unit_names(_write(tmp_path, "GDP.lib", LIB))


# --------------------------------------------------------------------------
# saying so when the types cannot be checked
# --------------------------------------------------------------------------

def _prompt_types(tmp_path, known):
    console = Console(record=True, width=100, force_terminal=False)
    frcmod = _write(tmp_path, "phos.frcmod", FRCMOD)
    with patch.object(lp, "_known_atom_types", return_value=known), \
         patch.object(lp, "prompt_with_context", return_value=""):
        _prompt_imported_atom_types(console, None, frcmod)
    return console.export_text()


def test_it_says_when_it_could_not_check_the_types(tmp_path):
    out = _prompt_types(tmp_path, set())
    assert "unverified" in out.lower() or "could not check" in out.lower()
    assert "AMBERHOME" in out


def test_the_warning_admits_standard_types_are_being_listed(tmp_path):
    out = _prompt_types(tmp_path, set())
    assert "CT" in out


def test_no_warning_when_the_parameter_files_were_found(tmp_path):
    out = _prompt_types(tmp_path, {"CT", "OS", "P", "H1", "O2"})
    assert "unverified" not in out.lower()
    assert "O3" in out          # the one genuinely new type is still asked about
