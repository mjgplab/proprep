"""
Mol2 generation: the right site's fingerprint, and metal bonds that exist.

Two defects behind ~160 tleap errors on a two-site 4UHX build.

1. Every atom in site 1's deposited library carried the type ``XX``:

       FS1.lib:  "FE1" "XX" ...      should be M1, M2, Y5, Y6
       MS1.lib:  "MO"  "M3" ...      site 2 was correct

   Step 3D locates standard.fingerprint through step_results["step_1"], which
   is restored from mcpb_step_results -- one workspace key shared by every
   site. With no per-site workflow_state.json, whichever site wrote step_1
   last supplied the fingerprint for both. The two sites' PDB serial ranges do
   not overlap (1764-2323/20396-20399 vs 20362-20404), so every lookup missed
   and fell to the 'XX' default. XX is not an Amber type, so tleap reported a
   missing parameter for every bond and angle in the residue.

2. _COVALENT_RADII had no metals, so Mo fell to the 0.77 default. Against the
   reported coordinates:

       Mo-O2  1.72 A  <  1.88 A cutoff   bonded
       Mo-O1  1.96 A  >  1.88 A cutoff   MISSED
       Mo-S   2.33 A  >  2.27 A cutoff   MISSED

   The MOS library held a 5-atom residue with 2 of its 4 bonds. Fe-S at
   ~2.2-2.3 A was passing only because it sat just under that accidental
   cutoff.

Adding metal radii then required excluding metal-metal pairs: Fe...Fe in
Fe2S2 is ~2.7 A and clears any radius cutoff, but MCPB does not bond them.
"""

import pytest

from proprep.forcefield_prep.mcpb.mol2_writer import Mol2Writer


# MOS as reported: Mo with a sulfido, an oxo and a hydroxo.
MOS = [
    {'serial': 20400, 'atom_name': 'MO', 'element': 'MO',
     'x': -29.569, 'y': -17.017, 'z': -41.812},
    {'serial': 20401, 'atom_name': 'S', 'element': 'S',
     'x': -29.451, 'y': -17.649, 'z': -39.574},
    {'serial': 20402, 'atom_name': 'O1', 'element': 'O',
     'x': -27.666, 'y': -16.557, 'z': -41.673},
    {'serial': 20403, 'atom_name': 'O2', 'element': 'O',
     'x': -29.607, 'y': -18.617, 'z': -42.453},
    {'serial': 20404, 'atom_name': 'H1', 'element': 'H',
     'x': -27.532, 'y': -15.655, 'z': -41.974},
]

# FES as reported: two irons bridged by two sulfides.
FES = [
    {'serial': 20396, 'atom_name': 'FE1', 'element': 'FE',
     'x': -46.078, 'y': -17.593, 'z': -47.380},
    {'serial': 20397, 'atom_name': 'FE2', 'element': 'FE',
     'x': -43.100, 'y': -18.216, 'z': -47.517},
    {'serial': 20398, 'atom_name': 'S1', 'element': 'S',
     'x': -44.425, 'y': -16.836, 'z': -48.664},
    {'serial': 20399, 'atom_name': 'S2', 'element': 'S',
     'x': -44.818, 'y': -19.213, 'z': -46.534},
]


def _bonds(atoms):
    writer = Mol2Writer()
    ordered = sorted(atoms, key=lambda a: a['serial'])
    ids = {a['serial']: i + 1 for i, a in enumerate(ordered)}
    names = {ids[a['serial']]: a['atom_name'] for a in ordered}
    return {tuple(sorted((names[i], names[j])))
            for i, j in writer._derive_intraresidue_bonds(atoms, ids)}


# --------------------------------------------------------------------------- #
# metal bond perception
# --------------------------------------------------------------------------- #

def test_the_molybdenum_cofactor_gets_all_four_bonds():
    assert _bonds(MOS) == {
        ('MO', 'S'), ('MO', 'O1'), ('MO', 'O2'), ('H1', 'O1'),
    }


@pytest.mark.parametrize("pair", [('MO', 'O1'), ('MO', 'S')])
def test_the_bonds_that_were_missed_are_now_found(pair):
    assert tuple(sorted(pair)) in _bonds(MOS)


def test_the_iron_sulfur_cluster_keeps_its_four_fe_s_bonds():
    assert _bonds(FES) == {
        ('FE1', 'S1'), ('FE1', 'S2'), ('FE2', 'S1'), ('FE2', 'S2'),
    }


def test_no_metal_metal_bond_is_invented():
    """
    Fe...Fe is ~2.7 A and clears any radius-based cutoff, but MCPB bridges the
    irons through the sulfides. Bonding them would need an M1-M2 force
    constant that Seminario never derives, and would fabricate a 3-membered
    ring for angle and dihedral generation to walk.
    """
    assert ('FE1', 'FE2') not in _bonds(FES)


def test_molybdenum_is_in_the_radius_table():
    """The default is carbon-like and silently under-bonds every metal."""
    assert Mol2Writer._COVALENT_RADII['MO'] == pytest.approx(1.54)
    assert Mol2Writer._COVALENT_RADII['FE'] == pytest.approx(1.32)


def test_geminal_hydrogens_are_still_not_bonded():
    """The property the radius approach was chosen for; metals must not break it."""
    hh = [
        {'serial': 1, 'atom_name': 'HB2', 'element': 'H', 'x': 0.0, 'y': 0.0, 'z': 0.0},
        {'serial': 2, 'atom_name': 'HB3', 'element': 'H', 'x': 1.78, 'y': 0.0, 'z': 0.0},
    ]

    assert _bonds(hh) == set()


# --------------------------------------------------------------------------- #
# the fingerprint must describe this model
# --------------------------------------------------------------------------- #

FINGERPRINT = """114-CYM-SG 1771 SH -> Y1
1311-FES-FE1 20396 FE -> M1
LINK 20396-FE1 1771-SG
"""


def _pdb_line(record, serial, name, resname, resid, x, y, z, element):
    """A column-correct PDB line; hand-typed ones drift and misparse."""
    return (f"{record:<6s}{serial:>5d} {name:^4s} {resname:>3s} A{resid:>4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n")


PDB = (
    _pdb_line("ATOM", 1771, "SG", "CYM", 114, -42.093, -17.443, -53.093, "S")
    + _pdb_line("HETATM", 20396, "FE1", "FES", 1311, -46.078, -17.593, -47.380, "FE")
)


def _write(tmp_path, fingerprint=FINGERPRINT, pdb=PDB):
    fp = tmp_path / "standard.fingerprint"
    fp.write_text(fingerprint)
    pdb_file = tmp_path / "large.pdb"
    pdb_file.write_text(pdb)
    topo = tmp_path / "bond_topology.json"
    topo.write_text('{"bonds": []}')
    return fp, pdb_file, topo


def test_a_fingerprint_from_another_site_is_refused(tmp_path):
    """
    The reported failure. Serials from site 2 against site 1's model: the old
    code wrote XX for every atom and the error surfaced only at the tleap
    build, after the libraries had been deposited.
    """
    other = "1310-MTE-N1 20362 nc -> nc\n1312-MOS-MO 20400 MO -> M3\n"
    fp, pdb_file, topo = _write(tmp_path, fingerprint=other)

    with pytest.raises(ValueError, match="different model"):
        Mol2Writer().write_mol2_files(
            output_dir=str(tmp_path / "out"), pdb_file=str(pdb_file),
            fingerprint_file=str(fp), resp_charges=[0.0, 0.0],
            bond_topology_file=str(topo))


def test_the_message_names_the_placeholder_and_the_files(tmp_path):
    other = "1312-MOS-MO 20400 MO -> M3\n"
    fp, pdb_file, topo = _write(tmp_path, fingerprint=other)

    with pytest.raises(ValueError) as exc:
        Mol2Writer().write_mol2_files(
            output_dir=str(tmp_path / "out"), pdb_file=str(pdb_file),
            fingerprint_file=str(fp), resp_charges=[0.0, 0.0],
            bond_topology_file=str(topo))

    message = str(exc.value)
    assert "XX" in message
    assert "standard.fingerprint" in message and "large.pdb" in message


def test_the_matching_fingerprint_is_accepted(tmp_path):
    fp, pdb_file, topo = _write(tmp_path)

    files = Mol2Writer().write_mol2_files(
        output_dir=str(tmp_path / "out"), pdb_file=str(pdb_file),
        fingerprint_file=str(fp), resp_charges=[-0.88, 0.83],
        bond_topology_file=str(topo))

    assert files
    written = (tmp_path / "out" / "FES1311.mol2").read_text()
    assert " M1 " in written and "XX" not in written


def test_a_partial_overlap_is_allowed(tmp_path):
    """
    A fingerprint covers only the site's atoms while large.pdb also holds caps,
    so overlap is normally partial. Only a total mismatch is an error.
    """
    padded = PDB + _pdb_line(
        "ATOM", 2261, "CH3", "ACE", 148, -40.0, -10.0, -40.0, "C")
    fp, pdb_file, topo = _write(tmp_path, pdb=padded)

    files = Mol2Writer().write_mol2_files(
        output_dir=str(tmp_path / "out"), pdb_file=str(pdb_file),
        fingerprint_file=str(fp), resp_charges=[-0.88, 0.83, 0.0],
        bond_topology_file=str(topo))

    assert files


def test_an_empty_fingerprint_is_not_treated_as_a_mismatch(tmp_path):
    """Nothing to compare is a different problem, reported elsewhere."""
    fp, pdb_file, topo = _write(tmp_path, fingerprint="")

    Mol2Writer().write_mol2_files(
        output_dir=str(tmp_path / "out"), pdb_file=str(pdb_file),
        fingerprint_file=str(fp), resp_charges=[0.0, 0.0],
        bond_topology_file=str(topo))
