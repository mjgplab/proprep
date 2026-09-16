"""Would a term written under shared atom types reach other residues?

The penalty table cannot answer that: a covalent-linkage torsion such as
CT-S-CT-CT carries no penalty (parmchk2 matches it through X-CT-S-X) yet no
standard residue contains it under ff14SB/ff19SB (MET and CYS use 2C), so a
refit is local. The ARG guanidinium H-N2-CA-N2, also unpenalized, is defined
explicitly by the force field, so a refit would rewrite every arginine.
"""

import os
import shutil

import pytest

from proprep.forcefield_prep import shared_type_terms as st

PARM = """title
C  12.01
CT 12.01

C   H   HO

C -CT  317.0  1.522

C -CT-N  63.0  110.1

X -CT-S -X    3    1.00          0.0             3.
CT-S -2C-CT   1    0.50          0.0             3.
H -N2-CA-N2   1    2.40        180.0             2.

CA-CA-CA-CT   1.1  180.0  2.0
"""

LIB = """!!index array str
 "AAA"
 "BBB"
!entry.AAA.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg
 "N" "N" 0 1 131072 1 7 0.0
 "CA" "XC" 0 1 131072 2 6 0.0
 "CB" "2C" 0 1 131072 3 6 0.0
 "SG" "S" 0 1 131072 4 16 0.0
 "C" "C" 0 1 131072 5 6 0.0
!entry.AAA.unit.connect array int
 1
 5
!entry.AAA.unit.connectivity table  int atom1x  int atom2x  int flags
 1 2 1
 2 3 1
 3 4 1
 2 5 1
!entry.BBB.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg
 "N" "N" 0 1 131072 1 7 0.0
 "CA" "XC" 0 1 131072 2 6 0.0
 "C" "C" 0 1 131072 3 6 0.0
!entry.BBB.unit.connect array int
 1
 3
!entry.BBB.unit.connectivity table  int atom1x  int atom2x  int flags
 1 2 1
 2 3 1
"""


def test_specific_dihedrals_ignore_wildcards_bonds_angles_and_impropers(tmp_path):
    p = tmp_path / "parm.dat"; p.write_text(PARM)
    got = st.parm_specific_dihedrals([str(p)])
    assert st.canonical(("CT", "S", "2C", "CT")) in got
    assert st.canonical(("H", "N2", "CA", "N2")) in got
    assert all("X" not in q for q in got)
    assert len(got) == 2                       # the improper (no IDIVF) is not a dihedral


def test_standard_residue_quads_include_intra_and_cross_residue_paths(tmp_path):
    p = tmp_path / "a.lib"; p.write_text(LIB)
    got = st.standard_residue_type_quads([str(p)])
    assert st.canonical(("N", "XC", "2C", "S")) in got          # within AAA
    assert st.canonical(("S", "2C", "XC", "C")) in got
    assert st.canonical(("2C", "XC", "C", "N")) in got          # AAA tail -> BBB head
    assert st.canonical(("XC", "C", "N", "XC")) in got          # spans the peptide bond
    assert st.canonical(("C", "N", "XC", "C")) in got           # into the next residue
    assert st.canonical(("CT", "S", "CT", "CT")) not in got


def test_conflict_reasons(tmp_path, monkeypatch):
    monkeypatch.setattr(st, "_protein_ff_facts",
                        lambda: ({st.canonical(("H", "N2", "CA", "N2"))},
                                 {st.canonical(("HC", "CT", "CT", "N"))}))
    assert "defines" in st.shared_type_dihedral_conflict(("N2", "CA", "N2", "H"))
    assert "standard residue contains" in st.shared_type_dihedral_conflict(("N", "CT", "CT", "HC"))
    assert st.shared_type_dihedral_conflict(("CT", "S", "CT", "CT")) is None


_AH = os.environ.get("AMBERHOME") or os.path.dirname(os.path.dirname(shutil.which("tleap") or "/x/y"))


@pytest.mark.skipif(not os.path.isfile(f"{_AH}/dat/leap/lib/amino19.lib"), reason="Amber data not available")
def test_real_amber_data_localizes_the_linkage_and_flags_the_corrections(monkeypatch):
    monkeypatch.setenv("AMBERHOME", _AH)
    st.reset_cache()
    try:
        assert st.amber_data_available()
        assert st.shared_type_dihedral_conflict(("CT", "S", "CT", "CT")) is None       # linkage: local
        assert st.shared_type_dihedral_conflict(("2C", "2C", "S", "CT"))               # methionine
        assert st.shared_type_dihedral_conflict(("H", "N2", "CA", "N2"))               # arginine
        assert st.shared_type_dihedral_conflict(("CA", "C", "OH", "HO"))               # tyrosine
        assert st.shared_type_dihedral_conflict(("2C", "XC", "C", "N"))                # backbone
    finally:
        st.reset_cache()
