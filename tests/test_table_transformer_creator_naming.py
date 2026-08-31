"""Addressing and naming in the interactive transformer editor.

Two defects made the editor unusable for the case it exists to serve --
rewriting a PDB residue so it matches an imported force field:

  * every residue was unaddressable when the site carried BioPython's space
    insertion code, because the lookup compared it against ``''``;
  * new names were force-uppercased, so a lib whose tLEaP unit is lowercase
    (``gdp``) could never be matched by a rename.
"""

import pytest

from proprep.structure_prep.comprehensive_redox_detector import (
    RedoxSite, RedoxSiteAtom,
)
from proprep.redoxsite_prep.transformation.table_transformer_creator import (
    RecipeBuilder, RecipeError, _Structure, apply_command,
)

# The prep-derived GDP library uses the old ff94 asterisk names; the PDB uses
# primes. These are the atoms that have to be renamed for the two to meet.
PRIMED = ["O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]


def _site(insertion_code):
    """A one-residue GDP site whose atoms carry the given insertion code."""
    site = RedoxSite(site_id="site_1", structure_id="6OIM")
    site.atoms = [
        RedoxSiteAtom(
            chain="A", resname="GDP", resid=302, atom_name=name,
            coords=(float(i), 0.0, 0.0), element=name[0],
            insertion_code=insertion_code,
        )
        for i, name in enumerate(["PB"] + PRIMED)
    ]
    site.bonds = []
    return site


def _builder(insertion_code):
    site = _site(insertion_code)
    return RecipeBuilder(_Structure.from_redox_site(site), site)


# --------------------------------------------------------------------------
# insertion-code normalization
# --------------------------------------------------------------------------

@pytest.mark.parametrize("insertion_code", [" ", "", None])
def test_residue_is_addressable_whatever_the_absent_icode_looks_like(insertion_code):
    """BioPython's ' ', an empty string and None all mean "no insertion code"."""
    builder = _builder(insertion_code)
    assert builder.structure.resname_of("A", 302) == "GDP"
    msg = apply_command(builder, ["rename_atom", "A", "302", "O5'", "O5*"])
    assert "O5*" in msg


def test_space_icode_does_not_leak_into_the_structure():
    assert all(a.icode == "" for a in _builder(" ").structure.atoms)


def test_a_real_insertion_code_is_preserved():
    builder = _builder("B")
    assert builder.structure.atoms[0].icode == "B"
    # ...and is not silently addressable as the no-icode residue.
    assert builder.structure.resname_of("A", 302) is None


# --------------------------------------------------------------------------
# case handling
# --------------------------------------------------------------------------

def test_residue_rename_preserves_lowercase():
    """tLEaP unit names are case-sensitive: 'gdp' must not become 'GDP'."""
    builder = _builder(" ")
    apply_command(builder, ["rename_res", "A", "302", "gdp"])
    assert builder.structure.resname_of("A", 302) == "gdp"
    op = builder.operations[-1]
    assert op["action"]["change_residue_name"] == "gdp"


def test_atom_rename_preserves_case_of_the_new_name():
    builder = _builder(" ")
    apply_command(builder, ["rename_atom", "A", "302", "PB", "Pb1"])
    assert {a.name for a in builder.structure.find("A", 302)} >= {"Pb1"}


def test_old_atom_name_may_be_typed_in_any_case():
    builder = _builder(" ")
    apply_command(builder, ["rename_atom", "A", "302", "o5'", "O5*"])
    # The replay applier matches the PDB spelling exactly, so the recorded key
    # must be the structure's name, not what the user typed.
    assert builder.operations[-1]["action"]["rename_atoms"] == {"O5'": "O5*"}


def test_unknown_atom_still_raises():
    builder = _builder(" ")
    with pytest.raises(RecipeError, match="not found"):
        apply_command(builder, ["rename_atom", "A", "302", "O9'", "O9*"])


# --------------------------------------------------------------------------
# the whole job the editor exists for
# --------------------------------------------------------------------------

def test_prime_to_asterisk_sweep_matches_the_ff94_library():
    builder = _builder(" ")
    for name in PRIMED:
        apply_command(builder, ["rename_atom", "A", "302", name, name[:-1] + "*"])
    apply_command(builder, ["rename_res", "A", "302", "gdp"])

    names = {a.name for a in builder.structure.find("A", 302)}
    assert names == {"PB"} | {n[:-1] + "*" for n in PRIMED}
    assert builder.structure.resname_of("A", 302) == "gdp"

    mappings = {}
    for op in builder.operations:
        mappings.update(op["action"].get("rename_atoms", {}))
    assert mappings == {n: n[:-1] + "*" for n in PRIMED}
