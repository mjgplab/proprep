"""
Regression test for the small-model capping bug that broke a peptide backbone.

7K0W's small model had GLU104 (sidechain-coordinating, kept FULL) adjacent to
ILE105 (backbone-O-coordinating). ILE105 was reduced to an ACE cap, which DROPS
its backbone N — severing the GLU104(C=O)…ILE105(N) peptide bond. That left
GLU104's carbonyl as a dangling acyl radical and over-methylated ILE105's CA,
producing an ODD-electron small model (587 e-) that Gaussian rejected against
multiplicity 11 ("combination of multiplicity 11 and 587 electrons is impossible").

Fix: when a backbone-O-coordinating residue's N-terminal neighbour is in the
model, keep its backbone (GLY/FULL) instead of ACE so the peptide bond survives.
Case 3 (backbone-N coordinating with C-terminal neighbour in model) is the mirror.
"""

from proprep.forcefield_prep.model_builder import (
    SmallModelBuilder, ModelResidue, CapType,
)


class _Console:
    def print(self, *a, **k):
        pass


class _Atom:
    def __init__(self, name):
        self._n = name
    def get_id(self):
        return self._n


class _Res:
    def __init__(self, names):
        self._atoms = [_Atom(n) for n in names]
    def get_atoms(self):
        return self._atoms


class _Bond:
    """A coordinate bond: atom1 = metal center, atom2 = ligand atom."""
    def __init__(self, metal_coords, lig_coords, chain, resid, atom_name):
        self.chemical_type = 'coordinate'
        self.atom1_coords = metal_coords
        self.atom1_residue_info = {'chain': 'A', 'resid': 185, 'atom_name': 'MN'}
        self.atom2_coords = lig_coords
        self.atom2_residue_info = {'chain': chain, 'resid': resid, 'atom_name': atom_name}


class _Center:
    def __init__(self, coords):
        self.coords = coords


class _Site:
    def __init__(self, centers, bonds):
        self.centers = centers
        self.bonds = bonds


def _builder(model_residues, bonds, residue_map):
    b = object.__new__(SmallModelBuilder)
    b.console = _Console()
    b.model_residues = model_residues
    b.residue_map = residue_map
    b.redox_site = _Site([_Center((0.0, 0.0, 0.0))], bonds)
    # Terminal-position helpers aren't relevant to these cases.
    b._is_nterm_residue = lambda chain, resid: False
    b._is_cterm_residue = lambda chain, resid: False
    return b


AA_BACKBONE = ['N', 'CA', 'C', 'O']


def test_backbone_O_residue_keeps_N_when_prev_in_model():
    # GLU104 (sidechain OE2 coordinates) — FULL, and ILE105 (backbone O
    # coordinates) adjacent. ILE105 must NOT become ACE (which drops N).
    glu = ModelResidue('A', 104, 'GLU', CapType.NONE)
    ile = ModelResidue('A', 105, 'ILE', CapType.NONE)
    bonds = [
        _Bond((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), 'A', 104, 'OE2'),
        _Bond((0.0, 0.0, 0.0), (2.0, 2.0, 2.0), 'A', 105, 'O'),
    ]
    residue_map = {
        ('A', 104): _Res(AA_BACKBONE + ['CB', 'CG', 'CD', 'OE1', 'OE2']),
        ('A', 105): _Res(AA_BACKBONE + ['CB']),
    }
    b = _builder([glu, ile], bonds, residue_map)
    b._apply_small_model_capping([('A', 104), ('A', 105)])

    caps = {(r.chain, r.resid): r.cap_type for r in b.model_residues}
    # The peptide bond is preserved: ILE105 keeps its backbone (GLY), not ACE.
    assert caps[('A', 105)] == CapType.GLY, caps
    assert caps[('A', 105)] != CapType.ACE
    # GLU104 stays FULL and its N-terminal side is capped by an ACE at 103.
    assert caps[('A', 104)] == CapType.FULL
    assert caps.get(('A', 103)) == CapType.ACE
    # C-terminal side of ILE105 (residue 106 absent) gets an NME cap.
    assert caps.get(('A', 106)) == CapType.NME


def test_isolated_backbone_O_residue_still_becomes_ACE():
    # No N-terminal neighbour in the model -> legacy ACE behavior preserved.
    ile = ModelResidue('A', 105, 'ILE', CapType.NONE)
    bonds = [_Bond((0.0, 0.0, 0.0), (2.0, 2.0, 2.0), 'A', 105, 'O')]
    residue_map = {('A', 105): _Res(AA_BACKBONE + ['CB'])}
    b = _builder([ile], bonds, residue_map)
    b._apply_small_model_capping([('A', 105)])

    caps = {(r.chain, r.resid): r.cap_type for r in b.model_residues}
    assert caps[('A', 105)] == CapType.ACE, caps
    assert caps.get(('A', 106)) == CapType.NME


def test_backbone_N_residue_keeps_CO_when_next_in_model():
    # Mirror case: residue coordinates via backbone N, C-terminal neighbour in
    # model -> must keep C=O (GLY), not reduce to NME.
    r1 = ModelResidue('A', 50, 'ALA', CapType.NONE)   # backbone-N coordinating
    r2 = ModelResidue('A', 51, 'ALA', CapType.NONE)   # in-model neighbour
    bonds = [
        _Bond((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), 'A', 50, 'N'),
        _Bond((0.0, 0.0, 0.0), (3.0, 3.0, 3.0), 'A', 51, 'N'),
    ]
    residue_map = {
        ('A', 50): _Res(AA_BACKBONE + ['CB']),
        ('A', 51): _Res(AA_BACKBONE + ['CB']),
    }
    b = _builder([r1, r2], bonds, residue_map)
    b._apply_small_model_capping([('A', 50), ('A', 51)])

    caps = {(r.chain, r.resid): r.cap_type for r in b.model_residues}
    assert caps[('A', 50)] == CapType.GLY, caps
    assert caps[('A', 50)] != CapType.NME
