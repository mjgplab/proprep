"""
GLY gap filling can change the large model's net charge without saying so.

Two runs on the same 4UHX Fe2S2 site produced different models:

    Aug 15  gap filling = actual PDB residues  GLY115 PHE116 ARG150  104 atoms  charge -1
    Aug 16  gap filling = GLY                  GLY115 GLY116 GLY150   89 atoms  charge -2

The Fe2S2 motif is C114-x-x-C117 ... C149-x-C151, so the gap at A:150 sits
between two coordinating cysteines -- and in this protein it is an arginine.
Answering "use GLY" replaced a +1 residue with a neutral one, and the model
charge moved by exactly that unit. Nothing in the prompt said so; the two
answers read as interchangeable formatting choices.

``preview_gap_residues`` reports what the gaps actually contain, so the charge
consequence is visible before the answer is given.
"""

from types import SimpleNamespace

import pytest

from proprep.forcefield_prep.model_builder import LargeModelBuilder


class _Residue:
    def __init__(self, resname):
        self._resname = resname

    def get_resname(self):
        return self._resname


def _builder(residues, proteogenic=None):
    """A LargeModelBuilder with only the state preview_gap_residues touches."""
    builder = LargeModelBuilder.__new__(LargeModelBuilder)
    builder.residue_map = {k: _Residue(v) for k, v in residues.items()}
    allowed = proteogenic if proteogenic is not None else set(residues)
    builder.is_proteogenic = lambda chain, resid: (chain, resid) in allowed
    return builder


# The reported site: two Cys pairs, one bridged by GLY/PHE, one by ARG.
SITE = {
    ('A', 114): 'CYS', ('A', 115): 'GLY', ('A', 116): 'PHE', ('A', 117): 'CYS',
    ('A', 149): 'CYS', ('A', 150): 'ARG', ('A', 151): 'CYS',
}
SELECTED = [('A', 114), ('A', 117), ('A', 149), ('A', 151)]


def test_the_reported_gaps_are_reported():
    preview = _builder(SITE).preview_gap_residues(SELECTED, max_gap=5)

    assert [(c, r, n) for c, r, n, _ in preview] == [
        ('A', 115, 'GLY'), ('A', 116, 'PHE'), ('A', 150, 'ARG'),
    ]


def test_the_arginine_is_flagged_as_charged():
    preview = _builder(SITE).preview_gap_residues(SELECTED, max_gap=5)

    charges = {resname: charge for _c, _r, resname, charge in preview}
    assert charges['ARG'] == +1
    assert charges['PHE'] == 0
    assert charges['GLY'] == 0


def test_the_charge_delta_matches_the_reported_run():
    """-1 with the real residues, -2 with GLY: one unit, from the ARG."""
    preview = _builder(SITE).preview_gap_residues(SELECTED, max_gap=5)

    assert sum(charge for *_r, charge in preview) == +1


@pytest.mark.parametrize("resname,charge", [
    ('ARG', +1), ('LYS', +1), ('HIP', +1),
    ('ASP', -1), ('GLU', -1),
    ('PHE', 0), ('GLY', 0), ('HIE', 0), ('HID', 0), ('ASH', 0), ('GLH', 0),
])
def test_formal_charges(resname, charge):
    residues = {('A', 1): 'CYS', ('A', 2): resname, ('A', 3): 'CYS'}
    preview = _builder(residues).preview_gap_residues(
        [('A', 1), ('A', 3)], max_gap=5)

    assert preview[0][3] == charge


def test_neutral_histidine_is_not_flagged():
    """HIE/HID carry no charge; only HIP does."""
    residues = {('A', 1): 'CYS', ('A', 2): 'HIE', ('A', 3): 'CYS'}

    preview = _builder(residues).preview_gap_residues([('A', 1), ('A', 3)], max_gap=5)

    assert preview[0][3] == 0


# --------------------------------------------------------------------------- #
# it must agree with what _fill_gaps would actually do
# --------------------------------------------------------------------------- #

def test_a_gap_wider_than_max_gap_is_not_previewed():
    preview = _builder(SITE).preview_gap_residues(SELECTED, max_gap=1)

    # A:114->117 is a 2-residue gap; A:149->151 is a 1-residue gap.
    assert [(c, r) for c, r, *_ in preview] == [('A', 150)]


def test_a_nonproteogenic_flank_is_not_bridged():
    """The ligand->metal span guard in _fill_gaps."""
    residues = {('A', 259): 'X9E', ('A', 260): 'HOH', ('A', 261): 'ZN'}
    builder = _builder(residues, proteogenic=set())

    assert builder.preview_gap_residues([('A', 259), ('A', 261)], max_gap=5) == []


def test_adjacent_residues_have_no_gap():
    assert _builder(SITE).preview_gap_residues(
        [('A', 114), ('A', 115)], max_gap=5) == []


def test_gaps_are_found_per_chain():
    residues = {('A', 1): 'CYS', ('A', 2): 'ARG', ('A', 3): 'CYS',
                ('B', 1): 'CYS', ('B', 2): 'GLU', ('B', 3): 'CYS'}

    preview = _builder(residues).preview_gap_residues(
        [('A', 1), ('A', 3), ('B', 1), ('B', 3)], max_gap=5)

    assert sorted((c, r, n) for c, r, n, _ in preview) == [
        ('A', 2, 'ARG'), ('B', 2, 'GLU')]
    assert sum(ch for *_r, ch in preview) == 0     # +1 and -1 cancel


def test_unsorted_selection_still_finds_the_gap():
    """Selections arrive in coordination order, not residue order."""
    preview = _builder(SITE).preview_gap_residues(
        [('A', 151), ('A', 114), ('A', 149), ('A', 117)], max_gap=5)

    assert ('A', 150, 'ARG', 1) in preview


def test_a_residue_missing_from_the_structure_is_skipped():
    residues = {('A', 1): 'CYS', ('A', 3): 'CYS'}   # A:2 absent

    assert _builder(residues).preview_gap_residues(
        [('A', 1), ('A', 3)], max_gap=5) == []


# --------------------------------------------------------------------------- #
# both model builders need it
# --------------------------------------------------------------------------- #

def test_the_small_model_builder_has_it_too():
    """
    The prompt that surprised the user was the SMALL model's. The preview
    lives on the shared base so both builders offer it; the small model fills
    only 1-residue gaps, so it asks with max_gap=1.
    """
    from proprep.forcefield_prep.model_builder import ModelBuilder, SmallModelBuilder

    assert hasattr(SmallModelBuilder, "preview_gap_residues")
    assert SmallModelBuilder.preview_gap_residues is ModelBuilder.preview_gap_residues
    assert LargeModelBuilder.preview_gap_residues is ModelBuilder.preview_gap_residues


def test_the_small_model_sees_only_single_residue_gaps():
    """max_gap=1 matches _fill_single_residue_gaps, which is all it bridges."""
    from proprep.forcefield_prep.model_builder import SmallModelBuilder

    builder = SmallModelBuilder.__new__(SmallModelBuilder)
    builder.residue_map = {k: _Residue(v) for k, v in SITE.items()}
    builder.is_proteogenic = lambda chain, resid: (chain, resid) in SITE

    preview = builder.preview_gap_residues(SELECTED, max_gap=1)

    # A:114->117 is a 2-residue gap and is not bridged; A:149->151 is.
    assert [(c, r, n) for c, r, n, _ in preview] == [('A', 150, 'ARG')]
