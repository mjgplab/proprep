"""
Scaffolding residues must be constrained to their own charge, not to zero.

The large model's non-metal-site residues -- ACE/NME caps and any residue added
to bridge a gap -- get a RESP group constraint and no mol2 file. That target was
hardcoded to 0.00000, which is right for a cap and for a GLY bridge and wrong
for a charged residue.

4UHX site 1 keeps the real ARG150 as the bridge between the coordinating CYS149
and CYS151. The ESP was computed with its guanidinium present, but pinning the
ARG at zero left its +1 nowhere to go except the free atoms, and RESP put it on
the nearest ones:

    CYM114  -0.4770   remote from ARG
    CYM117  -0.4899   remote from ARG
    CYM149  -0.2748   adjacent to ARG
    CYM151  +0.2344   adjacent to ARG      <- a cysteine cannot be +0.23
    FES1311 +0.0074
    Total   -1.0000   <- the whole model charge, scaffolding forced to 0

With the ARG at its real +1 the site sums to -2, which is [Fe2S2(SCys)4]2-.

The target is derived from the force field rather than a table. Only whole
residues are ever added as bridges and a whole residue carries an integer
charge, so a non-integral sum did not come from one -- it is a cap, whose
synthetic atom names only partially resolve -- and 0 is correct for those.
"""

import pytest

from proprep.forcefield_prep.mcpb.resp_input_generator import RESPInputGenerator


# ff14SB ARG, which sums to exactly +1.
ARG = {
    'N': -0.3479, 'H': 0.2747, 'CA': -0.2637, 'HA': 0.1560,
    'CB': -0.0007, 'HB2': 0.0327, 'HB3': 0.0327,
    'CG': 0.0390, 'HG2': 0.0285, 'HG3': 0.0285,
    'CD': 0.0486, 'HD2': 0.0687, 'HD3': 0.0687,
    'NE': -0.5295, 'HE': 0.3456, 'CZ': 0.8076,
    'NH1': -0.8627, 'HH11': 0.4478, 'HH12': 0.4478,
    'NH2': -0.8627, 'HH21': 0.4478, 'HH22': 0.4478,
    'C': 0.7341, 'O': -0.5894,
}
# ff14SB GLY, which sums to 0.
GLY = {'N': -0.4157, 'H': 0.2719, 'CA': -0.0252, 'HA2': 0.0698,
       'HA3': 0.0698, 'C': 0.5973, 'O': -0.5679}
# ff14SB GLU, which sums to -1.
GLU = {'N': -0.5163, 'H': 0.2936, 'CA': 0.0397, 'HA': 0.1105,
       'CB': 0.0560, 'HB2': -0.0173, 'HB3': -0.0173,
       'CG': 0.0136, 'HG2': -0.0425, 'HG3': -0.0425,
       'CD': 0.8054, 'OE1': -0.8188, 'OE2': -0.8188,
       'C': 0.5366, 'O': -0.5819}


def _atoms(resname, charges, chain='A', resid=150):
    return [{'chain': chain, 'resid': resid, 'resname': resname,
             'atom_name': name, 'element': 'C'} for name in charges]


def _lookup(resname, charges):
    return {(resname, name): q for name, q in charges.items()}


def _target(resname, charges, lookup=None):
    gen = RESPInputGenerator()
    atoms = _atoms(resname, charges)
    return gen._scaffolding_group_charge(
        list(range(1, len(atoms) + 1)), atoms, None,
        lookup if lookup is not None else _lookup(resname, charges))


# --------------------------------------------------------------------------- #
# the reported case
# --------------------------------------------------------------------------- #

def test_an_arginine_bridge_is_constrained_to_plus_one():
    assert _target('ARG', ARG) == pytest.approx(1.0)


def test_the_arg_sum_really_is_integral():
    """Guards the premise: ff14SB ARG sums to +1, so the integer test passes."""
    assert sum(ARG.values()) == pytest.approx(1.0, abs=1e-3)


def test_a_glycine_bridge_is_still_zero():
    """The default answer to the gap prompt must be unaffected."""
    assert _target('GLY', GLY) == pytest.approx(0.0)


def test_a_glutamate_bridge_is_constrained_to_minus_one():
    assert _target('GLU', GLU) == pytest.approx(-1.0)


# --------------------------------------------------------------------------- #
# caps must stay at zero
# --------------------------------------------------------------------------- #

def test_an_unresolvable_residue_is_zero():
    """A cap's synthetic atom names have no library entry."""
    assert _target('ACE', {'CH3': 0.0, 'HH31': 0.0}, lookup={}) == 0.0


def test_a_partially_resolvable_residue_is_zero():
    """
    The dangerous case: some atoms resolve and some do not, so a naive sum
    would produce an arbitrary fractional target.
    """
    partial = {k: v for k, v in list(ARG.items())[:6]}
    assert _target('ARG', ARG, lookup=_lookup('ARG', partial)) == 0.0


def test_a_nonintegral_sum_is_rejected():
    """Only whole residues are added as bridges, so this is not one."""
    assert _target('XXX', {'A': 0.3, 'B': 0.2}) == 0.0


def test_a_sum_just_off_an_integer_is_accepted():
    """Rounding in library charges must not disqualify a real residue."""
    nudged = dict(ARG)
    nudged['O'] = ARG['O'] + 0.004

    assert _target('ARG', nudged) == pytest.approx(1.0)


# --------------------------------------------------------------------------- #
# the constraint reaches the file
# --------------------------------------------------------------------------- #

def _written_groups(tmp_path, groups):
    gen = RESPInputGenerator()
    atoms = [{'chain': 'A', 'resid': 1, 'resname': 'X',
              'atom_name': f'A{i}', 'element': 'C'} for i in range(1, 5)]
    out = tmp_path / "resp1.in"
    gen._write_respin1(out, atoms, -1, {}, groups)
    return out.read_text()


def test_a_nonzero_group_target_is_written(tmp_path):
    text = _written_groups(tmp_path, [([1, 2], 1.0)])

    assert "1.00000" in text
    assert "0.00000" not in text


def test_a_zero_group_target_is_still_written(tmp_path):
    text = _written_groups(tmp_path, [([1, 2], 0.0)])

    assert "0.00000" in text


def test_a_negative_group_target_is_written(tmp_path):
    assert "-1.00000" in _written_groups(tmp_path, [([1, 2], -1.0)])


def test_stage_two_writes_the_same_targets(tmp_path):
    """Both stages carry the constraint; only stage 1 was checked above."""
    gen = RESPInputGenerator()
    atoms = [{'chain': 'A', 'resid': 1, 'resname': 'X',
              'atom_name': f'A{i}', 'element': 'C'} for i in range(1, 5)]
    out = tmp_path / "resp2.in"
    gen._write_respin2(out, atoms, -1, {}, [], [([1, 2], 1.0)], set())

    assert "1.00000" in out.read_text()


# --------------------------------------------------------------------------- #
# the charge lookup the target is derived from
# --------------------------------------------------------------------------- #

def test_the_charge_lookup_skips_atoms_with_no_charge():
    """A withheld cluster atom has charge None; it must not enter the lookup."""
    lookup = RESPInputGenerator._build_charge_lookup({
        (1.0, 2.0, 3.0): {'resname': 'FES', 'atom_name': 'S1', 'charge': None},
        (4.0, 5.0, 6.0): {'resname': 'CYM', 'atom_name': 'SG', 'charge': -0.88},
    })

    assert ('FES', 'S1') not in lookup
    assert lookup[('CYM', 'SG')] == pytest.approx(-0.88)


def test_the_charge_lookup_tolerates_an_empty_input():
    assert RESPInputGenerator._build_charge_lookup(None) == {}
    assert RESPInputGenerator._build_charge_lookup({}) == {}
