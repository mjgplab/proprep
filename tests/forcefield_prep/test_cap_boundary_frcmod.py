"""
Regression test for the ACE/NME cap-boundary bonded-term bug.

When a metal site coordinates through a residue's backbone carbonyl O, that O is
retyped (O -> Y6). In the QM small model the residue's C-terminal neighbour is
truncated to an NME cap, so the peptide bond ILE:C - NME:N crosses a cap boundary.

The bonded topology used to be borrowed only from the metal-free preprocessing
prmtop, filtered to the *typed* fingerprint atoms. Cap atoms are never typed, so
the C-N bond was dropped and the boundary angle (Y6-C-N) and dihedral
(N-C-Y6-M1) never enumerated -> tleap aborted with "No angle/torsion parameters"
on the assembled topology (which HAS those bonds). The earlier distance-union
"fix" reused the same fingerprint filter and was inert for the cap boundary.

The real fix: perceive the cap-boundary bond by distance over the small-model
coords and type the cap atom with its standard ff14SB type, so the boundary
terms enumerate and resolve to ordinary force-field values.
"""

from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager


def _bare():
    return MetalSiteWorkflowManager.__new__(MetalSiteWorkflowManager)


def _write_small_pdb(tmp_path):
    # A backbone fragment: CA-C(=O) where O ligates the metal (fingerprint atoms),
    # plus an NME cap whose N sits at a real peptide-bond distance (~1.33 A) from C,
    # and an ACE cap whose C sits at ~1.33 A from a downstream N (fingerprint).
    # Coordinates are chosen so only the true C-N / C-N boundary pairs bond.
    lines = [
        # serial  name res  resid      x       y       z            elem
        ("ATOM", 10, "CA", "ILE", 105, 0.000, 0.000, 0.000, "C"),
        ("ATOM", 11, "C",  "ILE", 105, 1.520, 0.000, 0.000, "C"),
        ("ATOM", 12, "O",  "ILE", 105, 2.100, 1.080, 0.000, "O"),   # -> Y6 (ligates)
        ("ATOM", 20, "N",  "NME", 106, 2.100, -1.150, 0.000, "N"),  # cap, ~1.33 A from C
        ("ATOM", 21, "H",  "NME", 106, 1.600, -2.000, 0.000, "H"),  # cap H, far from any fp atom
        ("ATOM", 30, "N",  "GLU", 104, -1.330, 0.000, 0.000, "N"),  # fingerprint N
        ("ATOM", 40, "C",  "ACE", 103, -2.150, 1.080, 0.000, "C"),  # cap, ~1.33 A from GLU104:N
    ]
    text = ""
    for rec, serial, name, res, resid, x, y, z, elem in lines:
        text += (
            f"{rec:<6}{serial:>5} {name:<4} {res:>3} A{resid:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2}\n"
        )
    p = tmp_path / "small.pdb"
    p.write_text(text)
    return p


def test_perceives_and_types_cap_boundary_bonds(tmp_path):
    pdb = _write_small_pdb(tmp_path)
    # Fingerprint = the typed atoms only (CA, C, O of ILE105 and N of GLU104).
    fp_serials = {10, 11, 12, 30}

    inst = _bare()
    bonds, caps = inst._perceive_cap_boundary_bonds(str(pdb), fp_serials)

    bondset = {frozenset(b) for b in bonds}
    # NME:N bonds to the carbonyl C, ACE:C bonds to the downstream N.
    assert frozenset((20, 11)) in bondset, bonds
    assert frozenset((40, 30)) in bondset, bonds
    # The cap H (21) is not near any fingerprint atom -> no spurious boundary bond.
    assert not any(21 in b for b in bonds), bonds

    # Cap atoms are typed with their standard ff14SB types, renamed=False path.
    assert caps[20]["renamed_type"] == "N"
    assert caps[40]["renamed_type"] == "C"


def test_no_caps_when_none_present(tmp_path):
    # A small model with no ACE/NME residues yields no boundary bonds.
    text = (
        "ATOM     10 CA   ILE A 105       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM     11 C    ILE A 105       1.520   0.000   0.000  1.00  0.00           C\n"
    )
    p = tmp_path / "small.pdb"
    p.write_text(text)
    inst = _bare()
    bonds, caps = inst._perceive_cap_boundary_bonds(str(p), {10, 11})
    assert bonds == []
    assert caps == {}
