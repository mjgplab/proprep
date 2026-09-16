"""The transformer editor shows what the linked library expects.

After an import, the site-match line called a residue-NAME match "100% of
its atoms", and the editor never showed the library's atom names, so an
O5'/O5* difference was invisible until tLEaP failed. Now the match line
reports names and atoms separately, and the editor prints a three-way
comparison at startup and on 'lib <chain> <resid>'.
"""

from types import SimpleNamespace

from rich.console import Console

from proprep.forcefield_prep import library_promotion as lp
from proprep.forcefield_prep.forcefield_parameterizer import _residue_atom_overlap
from proprep.redoxsite_prep.transformation import table_transformer_creator as tc
from proprep.structure_prep.comprehensive_redox_detector import RedoxSite, RedoxSiteAtom

PRIMED = ["PB", "O1B", "O5'", "C5'", "C4'", "O4'", "C1'", "N9"]
STARRED = ["PB", "O1B", "O5*", "C5*", "C4*", "O4*", "C1*", "N9", "H8", "HO2'"]

LIB = ('!!index array str\n "GDP"\n'
       '!entry.GDP.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n'
       + "".join(f' "{n}" "X" 0 1 131072 {i} 6 0.0\n' for i, n in enumerate(STARRED, 1))
       + '!entry.GDP.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg\n'
       ' "PB" "P" 0 -1 0.0\n')


def _site():
    site = RedoxSite(site_id="site_2", structure_id="6OIM")
    site.atoms = [
        RedoxSiteAtom(chain="A", resname="GDP", resid=302, atom_name=n,
                      coords=(float(i), 0.0, 0.0), element=n[0], insertion_code="")
        for i, n in enumerate(PRIMED)
    ] + [RedoxSiteAtom(chain="A", resname="MG", resid=301, atom_name="MG",
                        coords=(9.0, 0.0, 0.0), element="MG", insertion_code="")]
    site.bonds = []
    return site


# --------------------------------------------------------------------------
# library atom table
# --------------------------------------------------------------------------

def test_library_unit_atoms_keeps_order_spelling_and_hydrogens(tmp_path):
    lib = tmp_path / "GDP.lib"
    lib.write_text(LIB)
    assert lp.library_unit_atoms(lib) == {"GDP": STARRED}


# --------------------------------------------------------------------------
# the site-match line after an import
# --------------------------------------------------------------------------

def test_overlap_counts_heavy_atoms_and_names_the_misses():
    site = SimpleNamespace(atoms=[SimpleNamespace(resname="GDP", atom_name=n) for n in PRIMED + ["H8"]])
    lib_names = {n.upper() for n in STARRED if not n.startswith("H")}
    hit, total, unmatched = _residue_atom_overlap(site, "GDP", lib_names)
    assert (hit, total) == (3, 8)                       # PB O1B N9 match; H8 not counted
    assert unmatched == ["O5'", "C5'", "C4'", "O4'", "C1'"]


def test_match_message_separates_name_match_from_atom_match(tmp_path, monkeypatch):
    from proprep.forcefield_prep import forcefield_parameterizer as fp
    out = Console(record=True, width=200, quiet=False, file=open(tmp_path / "o.txt", "w"))
    module = fp.ForcefieldParameterizer.__new__(fp.ForcefieldParameterizer)
    module.console = out
    ws = SimpleNamespace(get=lambda k, d=None: [_site()] if k == "detected_redox_sites" else d)
    module.processor = SimpleNamespace(_get_workspace=lambda: ws)
    monkeypatch.setattr(fp, "_imported_library_atom_names",
                        lambda result: {n.upper() for n in STARRED if not n.startswith("H")})
    monkeypatch.setattr(fp, "confirm_with_context", lambda *a, **k: False)
    module._offer_transformer_for_import({"residue_name": "GDP",
                                          "library_path": "/x/specialized_residues/small_molecules/GDP"})
    text = out.export_text()
    assert "100%" not in text
    assert "same name as the library" in text
    assert "3 of 8 heavy atoms have a name in the library" in text
    assert "O5'" in text and "C1'" in text


# --------------------------------------------------------------------------
# the editor
# --------------------------------------------------------------------------

class _Proc:
    def __init__(self):
        self.console = Console(record=True, width=120, file=open("/dev/null", "w"))


def _creator(units, answers, monkeypatch):
    proc = _Proc()
    creator = tc.TableTransformerCreator(
        proc, _site(), forcefield_default={"path": "small_molecules/GDP",
                                           "redox_state": "single_state", "spin_state": "default"})
    creator._library_units = units
    it = iter(answers)
    monkeypatch.setattr(tc, "prompt_with_context", lambda *a, **k: next(it))
    monkeypatch.setattr(tc, "confirm_with_context", lambda *a, **k: True)
    return creator, proc


def test_startup_comparison_names_the_primes_and_the_stars(monkeypatch):
    creator, proc = _creator({"GDP": STARRED}, ["quit"], monkeypatch)
    assert creator.create() is None
    text = proc.console.export_text()
    assert "Linked library: small_molecules/GDP [single_state/default]" in text
    assert "A 302 GDP (8 heavy atoms) vs library unit GDP (8 heavy atoms)" in text
    assert "in both (3)" in text
    assert "in the structure only (5)" in text and "O5' C5' C4' O4' C1'" in text
    assert "in the library only (5)" in text and "O5* C5* C4* O4* C1*" in text
    assert "(2 library hydrogens not compared)" in text


def test_lib_command_lists_units_and_compares_a_residue(monkeypatch):
    creator, proc = _creator({"GDP": STARRED}, ["lib", "lib A 301", "lib B 9", "quit"], monkeypatch)
    creator.create()
    text = proc.console.export_text()
    assert "unit GDP: 8 heavy atoms, 2 hydrogens" in text
    assert "hydrogens: H8 HO2'" in text
    assert "A 301 MG (1 heavy atoms) vs library unit GDP (closest unit by atom names)" in text
    assert "No residue B 9 in this site" in text


def test_lib_command_explains_when_no_library_is_linked(monkeypatch):
    creator, proc = _creator({}, ["lib", "quit"], monkeypatch)
    creator.create()
    assert "No library is linked yet" in proc.console.export_text()


def test_comparison_reports_case_only_differences():
    structure = tc._Structure.from_redox_site(_site())
    units = {"GDP": ["pb", "O1B", "O5*"]}
    text = tc.render_library_comparison(structure, ("A", 302, ""), units)
    assert "case differs (1) -- tLEaP is case-sensitive:" in text
    assert "PB -> pb" in text
    assert "in the library only (1)" in text and "O5*" in text
