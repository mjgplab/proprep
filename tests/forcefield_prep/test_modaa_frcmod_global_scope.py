"""A modified residue's frcmod must not change the standard residues.

The modified-AA route runs parmchk2 with -a Y, so its frcmod carries a copy
of every bonded term of the residue keyed by the standard amber types. tLEaP
loads that file after the protein leaprc and gives it priority for every
residue in the system. That is harmless exactly when the copies equal what
the leaprc loads, and only penalized / ATTN lines are ever modified.

Two ways that broke:
  * the ff19SB path handed parmchk2 parm19.dat without frcmod.ff19SB, so the
    copies carried the older parm19 values and replaced ff19SB's multi-term
    ARG guanidinium and TYR hydroxyl torsions in every ARG and TYR;
  * Seminario's "all" scope refined the directly matched (standard) terms
    from the residue's Hessian and wrote them under the shared types.
"""

import os
import re
import shutil
import subprocess

import pytest

from proprep.forcefield_prep import modified_amino_acid_parameterizer as M

AH = os.environ.get("AMBERHOME") or os.path.dirname(os.path.dirname(shutil.which("tleap") or "/x/y"))
_HAVE = shutil.which("tleap") and shutil.which("parmchk2") and os.path.isfile(f"{AH}/dat/leap/prep/all_amino03.in")


# --------------------------------------------------------------------------
# parmchk2 is given the ff19SB corrections
# --------------------------------------------------------------------------

def test_run_parmchk2_passes_afrc(tmp_path, monkeypatch):
    seen = {}

    def fake_run(cmd, **kw):
        seen["cmd"] = cmd
        open(cmd[cmd.index("-o") + 1], "w").write("x")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(subprocess, "run", fake_run)
    prep = tmp_path / "x.prep"; prep.write_text("p")
    M.run_parmchk2(str(prep), str(tmp_path / "x.frcmod"), parm_dat_file="/p/parm19.dat",
                   afrc_file="/p/frcmod.ff19SB")
    cmd = seen["cmd"]
    assert cmd[cmd.index("-p") + 1] == "/p/parm19.dat"
    assert cmd[cmd.index("-afrc") + 1] == "/p/frcmod.ff19SB"
    assert "-a" in cmd and cmd[cmd.index("-a") + 1] == "Y"


def test_ff19sb_choice_supplies_the_correction_file(tmp_path, monkeypatch):
    """The FF chooser inside generate_bonded_parameters passes frcmod.ff19SB along."""
    amber = tmp_path / "amber"
    (amber / "dat" / "leap" / "parm").mkdir(parents=True)
    (amber / "dat" / "leap" / "parm" / "parm19.dat").write_text("parm")
    (amber / "dat" / "leap" / "parm" / "frcmod.ff19SB").write_text("frc")
    monkeypatch.setenv("AMBERHOME", str(amber))
    monkeypatch.chdir(tmp_path)
    (tmp_path / "cs1.prep").write_text("prep")
    seen = {}
    monkeypatch.setattr(M, "run_parmchk2",
                        lambda *a, **k: seen.update(k) or {"success": False, "error": "stop here"})
    answers = iter(["1", "CS1", "1"])   # prep #1, residue symbol, ff19SB
    monkeypatch.setattr(M, "prompt_with_context", lambda *a, **k: next(answers))
    M.generate_bonded_parameters("CS1", None)
    assert seen["parm_dat_file"].endswith("parm19.dat")
    assert seen["afrc_file"].endswith("frcmod.ff19SB")


# --------------------------------------------------------------------------
# Seminario's 'all' scope is refused with an explanation, and re-prompted
# --------------------------------------------------------------------------

def test_seminario_all_scope_is_explained_and_reprompted(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    m = object.__new__(M.ModifiedAAWorkflowManager)
    m.console = M._console
    m.processor = None
    m.amino_acid = "CS1"
    m.conformers = ["xtal"]
    m._from_structure_constraint_mode = "restrain"
    (tmp_path / "cs1_xtal_opt.fchk").write_text("fchk")
    frcmod = tmp_path / "cs1.frcmod"
    frcmod.write_text("remark\nMASS\n\nBOND\nCT-N   337.00   1.449\n"
                      "CT-S1  230.00   1.810       same as CT-S , penalty score=  2.0\n\n"
                      "ANGLE\nC -CT-N    63.000     110.100\n\nDIHE\n\nIMPROPER\n\nNONBON\n\nEND\n")
    answers = iter(["a", "a", "b"])
    monkeypatch.setattr(M, "prompt_with_context", lambda *a, **k: next(answers))
    monkeypatch.setattr(m, "_ac_to_mol2", lambda rs: "cs1.mol2")
    got = {}
    import proprep.forcefield_prep.seminario_refinement as sr
    monkeypatch.setattr(sr, "run_seminario_refinement_workflow",
                        lambda **kw: got.update(kw) or {"refinement_success": False, "message": "stub"})
    m._maybe_seminario("CS1", str(frcmod), "cs1.prep")
    # after two refused 'a's the third answer 'b' ran, with ONLY the by-analogy term
    assert [p[0] for p in got["selected_params"]] == ["CT-S1"]
    assert next(answers, "consumed") == "consumed"


# --------------------------------------------------------------------------
# regression: a ProPrep-style frcmod loaded into an ff19SB peptide changes nothing
# --------------------------------------------------------------------------

def _split_prep(dest):
    lines = open(f"{AH}/dat/leap/prep/all_amino03.in").read().splitlines()
    head, blocks, cur = lines[:2], [], []
    for ln in lines[2:]:
        if ln.strip() == "STOP":
            break
        cur.append(ln)
        if ln.strip() == "DONE":
            blocks.append(cur); cur = []
    out = {}
    for b in blocks:
        m = next((re.match(r"\s*(\S+)\s+INT", l) for l in b if re.match(r"\s*(\S+)\s+INT", l)), None)
        if m and m.group(1) in ("ARG", "TYR", "LYS"):
            p = dest / f"prep_{m.group(1)}.in"
            p.write_text("\n".join(head + b + ["STOP"]) + "\n")
            out[m.group(1)] = str(p)
    return out


def _terms(parm7):
    import parmed
    p = parmed.load_file(parm7)
    out = []
    for t in p.bonds:
        out.append(("b", t.atom1.idx, t.atom2.idx, round(t.type.k, 4), round(t.type.req, 4)))
    for t in p.angles:
        out.append(("a", t.atom1.idx, t.atom2.idx, t.atom3.idx, round(t.type.k, 4), round(t.type.theteq, 4)))
    for t in p.dihedrals:
        v = (round(t.type.phi_k, 4), t.type.per, round(t.type.phase, 2))
        if t.improper:
            # tLEaP may list an improper's two leading atoms in either order
            # once a frcmod is loaded; the term is the same (same central
            # atom, same outer atoms, same parameters), so canonicalize.
            out.append(("i", t.atom3.idx, tuple(sorted((t.atom1.idx, t.atom2.idx, t.atom4.idx))), v))
        else:
            out.append(("d", (t.atom1.idx, t.atom2.idx, t.atom3.idx, t.atom4.idx), v))
    return sorted(out, key=str)


@pytest.mark.skipif(not _HAVE, reason="tleap/parmchk2/Amber data not available")
def test_ff19sb_frcmod_does_not_alter_standard_residues(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("AMBERHOME", AH)
    preps = _split_prep(tmp_path)
    frcmods = []
    for name, prep in preps.items():
        out = f"{name}.frcmod"
        # exactly what generate_bonded_parameters does for the ff19SB choice
        r = M.run_parmchk2(prep, out, parm_dat_file=f"{AH}/dat/leap/parm/parm19.dat",
                           afrc_file=f"{AH}/dat/leap/parm/frcmod.ff19SB")
        assert r["success"], r
        frcmods.append(out)
    script = ["source leaprc.protein.ff19SB", "m = sequence { ACE ARG TYR LYS NME }",
              "saveamberparm m base.parm7 base.rst7"]
    script += [f'loadamberparams "{f}"' for f in frcmods]
    script += ["saveamberparm m over.parm7 over.rst7", "quit"]
    (tmp_path / "t.in").write_text("\n".join(script) + "\n")
    subprocess.run(["tleap", "-f", "t.in"], capture_output=True, text=True, check=True)
    assert _terms("base.parm7") == _terms("over.parm7")

    # and the control: without the correction file the ARG/TYR torsions DO change
    bad = []
    for name, prep in preps.items():
        out = f"{name}_bad.frcmod"
        assert M.run_parmchk2(prep, out, parm_dat_file=f"{AH}/dat/leap/parm/parm19.dat")["success"]
        bad.append(out)
    script = ["source leaprc.protein.ff19SB", "m = sequence { ACE ARG TYR LYS NME }"]
    script += [f'loadamberparams "{f}"' for f in bad] + ["saveamberparm m bad.parm7 bad.rst7", "quit"]
    (tmp_path / "t2.in").write_text("\n".join(script) + "\n")
    subprocess.run(["tleap", "-f", "t2.in"], capture_output=True, text=True, check=True)
    assert _terms("base.parm7") != _terms("bad.parm7")
