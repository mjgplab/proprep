"""Three names are asked for across the import and transformer flows --
library entry, parameter set, template -- and the prompts alone do not say
which level each one sits at. Each flow explains its own before asking.
"""

from rich.console import Console

from proprep.forcefield_prep.library_promotion import _print_naming_help


def _render(family="small_molecules", is_metal=False):
    console = Console(record=True, width=100, force_terminal=False)
    _print_naming_help(console, family, is_metal)
    return console.export_text()


def test_both_levels_are_named_and_distinguished():
    out = _render()
    assert "Entry name" in out and "Set name" in out
    assert "which molecule" in out
    assert "whose numbers" in out


def test_the_entry_path_shows_where_it_lands():
    assert "specialized_residues/small_molecules/" in _render()


def test_it_says_several_sets_can_share_one_entry():
    out = _render()
    assert "several" in out
    assert "default" in out


def test_metal_sites_are_described_as_sites_not_molecules():
    out = _render(family="metal_sites", is_metal=True)
    assert "which site" in out
    assert "molecule" not in out
