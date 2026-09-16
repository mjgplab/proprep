"""Importing parameters shipped as prep + frcmod (the Bryce database shape).

The library, the loader and the transformers all want an OFF/lib; a prep on
its own had no consumer, so the import wizard refused it. It now converts the
prep with tLEaP at import time and deposits the resulting lib.
"""

import os
import shutil
from pathlib import Path

import pytest
from rich.console import Console

from proprep.forcefield_prep import library_promotion as lp

QUIET = Console(quiet=True)

SINGLE_PREP = """    0    0    2

heme from somewhere

HEM   INT  0
CORRECT  OMIT DU  BEG
  0.0000
   1  DUMM  DU    M    0  -1  -2     0.000     0.000     0.000   0.00000

DONE
STOP
"""

MULTI_PREP = """    1    1    2
db94.dat
ALANINE

 ALA  INT     1
 CORR OMIT DU   BEG
   0.00000
   1  DUMM  DU    M    0  -1  -2     0.000     0.000     0.000  0.00000

DONE
GLYCINE

 GLY  XYZ     1
 CORR OMIT DU   BEG
   0.00000

DONE
STOP
"""


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


# --------------------------------------------------------------------------
# the residue name comes from inside the prep
# --------------------------------------------------------------------------

def test_residue_names_are_read_from_the_int_xyz_line(tmp_path):
    assert lp.prep_residue_names(_write(tmp_path, "bryce_heme.prep", SINGLE_PREP)) == ["HEM"]


def test_every_residue_of_a_multi_residue_prep_is_returned_in_order(tmp_path):
    assert lp.prep_residue_names(_write(tmp_path, "amino.in", MULTI_PREP)) == ["ALA", "GLY"]


def test_a_missing_or_empty_prep_yields_nothing(tmp_path):
    assert lp.prep_residue_names(tmp_path / "nope.prep") == []
    assert lp.prep_residue_names(_write(tmp_path, "junk.prep", "no residue here\n")) == []


# --------------------------------------------------------------------------
# refusals that must not reach tleap
# --------------------------------------------------------------------------

def test_conversion_refuses_a_digit_leading_residue_name(tmp_path):
    prep = _write(tmp_path, "x.prep", SINGLE_PREP.replace("HEM   INT", "9E2   INT"))
    lib, why = lp.prep_to_lib(prep, tmp_path)
    assert lib is None and "9E2" in why and "number" in why


def test_conversion_refuses_a_prep_without_a_residue(tmp_path):
    lib, why = lp.prep_to_lib(_write(tmp_path, "junk.prep", "nothing\n"), tmp_path)
    assert lib is None and "no residue definition" in why


def test_conversion_says_when_tleap_is_missing(tmp_path, monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: None)
    lib, why = lp.prep_to_lib(_write(tmp_path, "h.prep", SINGLE_PREP), tmp_path)
    assert lib is None and "tleap" in why


# --------------------------------------------------------------------------
# the wizard accepts a prep where it used to insist on a lib
# --------------------------------------------------------------------------

def test_wizard_converts_a_prep_and_deposits_the_lib(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    frcmod = _write(tmp_path, "heme.frcmod", "remark\nMASS\n\nBOND\n\n")
    prep = _write(tmp_path, "heme.prep", SINGLE_PREP)
    answers = iter([frcmod, prep])
    monkeypatch.setattr(lp, "_prompt_existing_file",
                        lambda console, processor, label, **kw: next(answers, None))
    monkeypatch.setattr(lp, "_prompt_category", lambda console, processor: "1")
    monkeypatch.setattr(lp, "_prompt_imported_atom_types", lambda console, processor, f: [])
    converted = {}

    def fake_convert(prep_path, out_dir, frcmod_path=None, lib_name=None):
        converted.update(prep=prep_path, frcmod=frcmod_path, out_dir=out_dir)
        lib = Path(out_dir) / "HEM.lib"
        lib.write_text('!!index array str\n "HEM"\n!entry.HEM.unit.atoms table  str name  str type  int typex\n "FE" "FE" 0\n')
        return str(lib), ""

    monkeypatch.setattr(lp, "prep_to_lib", fake_convert)
    captured = {}

    def fake_build(console, processor, **kw):
        captured.update(kw)
        return None            # stop before writing to the real library

    monkeypatch.setattr(lp, "_build_request_interactively", fake_build)

    assert lp.run_import_wizard(QUIET, None) is None
    assert converted["prep"] == prep and converted["frcmod"] == frcmod
    assert captured["lib_file"].endswith("HEM.lib")
    assert captured["prep_file"] == prep                  # the prep rides along
    assert captured["residue_name"] == "HEM"              # unit name, not filename


def test_wizard_stops_cleanly_when_the_prep_cannot_be_converted(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    frcmod = _write(tmp_path, "x.frcmod", "remark\nMASS\n\n")
    prep = _write(tmp_path, "x.prep", SINGLE_PREP.replace("HEM   INT", "9E2   INT"))
    answers = iter([frcmod, prep])
    monkeypatch.setattr(lp, "_prompt_existing_file",
                        lambda console, processor, label, **kw: next(answers, None))
    built = []
    monkeypatch.setattr(lp, "_build_request_interactively",
                        lambda *a, **k: built.append(1))
    assert lp.run_import_wizard(QUIET, None) is None
    assert not built


def test_wizard_still_takes_a_lib_and_then_an_optional_prep(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    frcmod = _write(tmp_path, "x.frcmod", "remark\nMASS\n\n")
    lib = _write(tmp_path, "x.lib", '!!index array str\n "XYZ"\n!entry.XYZ.unit.atoms table  str name  str type  int typex\n')
    answers = iter([frcmod, lib, None])
    monkeypatch.setattr(lp, "_prompt_existing_file",
                        lambda console, processor, label, **kw: next(answers, None))
    monkeypatch.setattr(lp, "_prompt_category", lambda console, processor: "1")
    monkeypatch.setattr(lp, "_prompt_imported_atom_types", lambda console, processor, f: [])
    monkeypatch.setattr(lp, "prep_to_lib", lambda *a, **k: pytest.fail("no conversion for a lib"))
    captured = {}
    monkeypatch.setattr(lp, "_build_request_interactively",
                        lambda console, processor, **kw: captured.update(kw))
    lp.run_import_wizard(QUIET, None)
    assert captured["lib_file"] == lib and captured["prep_file"] is None


# --------------------------------------------------------------------------
# real tleap, on Amber's own prep files (only where AmberTools is on PATH)
# --------------------------------------------------------------------------

def _amber_dat():
    for root in (os.environ.get("AMBERHOME"), os.path.dirname(os.path.dirname(shutil.which("tleap") or "/x/y"))):
        if root and Path(root, "dat", "leap", "prep", "chcl3.in").is_file():
            return Path(root, "dat", "leap")
    return None


@pytest.mark.skipif(not (shutil.which("tleap") and _amber_dat()), reason="tleap / Amber data not available")
def test_real_tleap_converts_single_and_multi_residue_preps(tmp_path):
    dat = _amber_dat()
    lib, why = lp.prep_to_lib(dat / "prep" / "chcl3.in", tmp_path, frcmod_path=dat / "parm" / "frcmod.chcl3")
    assert lib, why
    assert lp.library_unit_names(lib) == ["CL3"]
    assert "CL3" in lp.library_atom_names(lib) or lp.library_atom_names(lib)   # atoms table present
    lib2, why2 = lp.prep_to_lib(dat / "prep" / "all_amino03.in", tmp_path)
    assert lib2, why2
    units = lp.library_unit_names(lib2)
    # saveOff writes the index alphabetically, so compare as sets
    assert len(units) == 33
    assert set(units) == set(lp.prep_residue_names(dat / "prep" / "all_amino03.in"))
