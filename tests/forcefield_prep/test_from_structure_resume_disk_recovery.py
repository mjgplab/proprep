"""
Regression tests: the Route B (from-structure) steps 8->9->10 recover their
inputs from disk when resumed.

The workflow's ``step_results`` is per-process and in-memory only — it is not
persisted to workflow_state.json. A step completed in a PRIOR session (step 6
in particular pauses for external QM) therefore leaves no in-memory result. When
the user resumes and runs a later step, that step used to fail with "Missing ESP
/ AC / prep / frcmod file" even though the artifact was on disk the whole time.

Each step now falls back to the artifact's canonical on-disk name (deterministic
for the from-structure route), while still preferring an in-memory result when
present. These tests pin both behaviors on a synthetic run directory.
"""

import os

import pytest

from proprep.forcefield_prep import modified_amino_acid_parameterizer as M


def _mgr(tmp_path, step_results=None):
    m = object.__new__(M.ModifiedAAWorkflowManager)
    m.console = M._console
    m.processor = None
    m.amino_acid = "CS1"
    m.conformer_mode = "from_structure"
    m.conformers = ["xtal"]
    m.starting_pdb = str(tmp_path / "MOV_CYS_capped_H.pdb")
    m.source_residues = [{"name": "CYS"}, {"name": "MOV"}]
    m.step_results = step_results if step_results is not None else {}
    return m


def _seed_run_dir(tmp_path, charge=0):
    """Write the canonical from-structure artifacts a resumed run must recover."""
    (tmp_path / "cs1_combined.esp").write_text("esp\n")
    (tmp_path / "CS1.ac").write_text("ac\n")
    (tmp_path / "cs1.prep").write_text("prep\n")
    (tmp_path / "cs1.frcmod").write_text("frcmod\n")
    # intermediates that must NOT be mistaken for the final frcmod
    (tmp_path / "cs1_temp.frcmod").write_text("temp\n")
    (tmp_path / "cs1_gaff.frcmod").write_text("gaff\n")
    (tmp_path / "MOV_CYS_capped_H.pdb").write_text("REMARK\n")
    (tmp_path / "cs1_xtal_opt.gjf").write_text(
        f"%chk=cs1_xtal_opt.chk\n%mem=8GB\n#p B3LYP/6-31+G(d) Opt Freq\n\n"
        f"title\n\n{charge} 1\nC 0.0 0.0 0.0\n")


@pytest.fixture
def run_dir(tmp_path, monkeypatch):
    _seed_run_dir(tmp_path)
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_step8_recovers_esp_ac_and_charge_on_resume(run_dir, monkeypatch):
    cap = {}
    monkeypatch.setattr(
        M, "generate_and_run_residuegen",
        lambda a, ac, esp, nc, processor=None: cap.update(ac=ac, esp=esp, nc=nc)
        or {"success": True, "prep_file": "cs1.prep", "residue_symbol": "CS1"})

    r = _mgr(run_dir)._run_step_8()

    assert r["success"] is True
    assert cap["esp"] == "cs1_combined.esp"
    assert cap["ac"] == "CS1.ac"
    assert cap["nc"] == 0


def test_step8_recovers_nonzero_charge_from_opt_gjf(tmp_path, monkeypatch):
    _seed_run_dir(tmp_path, charge=-2)
    monkeypatch.chdir(tmp_path)
    cap = {}
    monkeypatch.setattr(
        M, "generate_and_run_residuegen",
        lambda a, ac, esp, nc, processor=None: cap.update(nc=nc)
        or {"success": True, "prep_file": "cs1.prep", "residue_symbol": "CS1"})

    _mgr(tmp_path)._run_step_8()
    assert cap["nc"] == -2


def test_step9_recovers_residue_symbol_on_resume(run_dir, monkeypatch):
    cap = {}
    monkeypatch.setattr(
        M, "generate_bonded_parameters",
        lambda rs, proc, standalone_use=True: cap.update(symbol=rs)
        or {"success": True, "final_frcmod": "cs1.frcmod", "prep_file": "cs1.prep"})
    m = _mgr(run_dir)
    monkeypatch.setattr(m, "_maybe_seminario", lambda *a, **k: None)
    monkeypatch.setattr(m, "_maybe_torsion_refit", lambda *a, **k: None)

    r = m._run_step_9_from_structure(interactive=False)
    assert r["success"] is True
    assert cap["symbol"] == "CS1"


def test_step10_recovers_prep_and_final_frcmod_on_resume(run_dir, monkeypatch):
    cap = {}
    monkeypatch.setattr(
        M, "integrate_modaa_from_structure",
        lambda **kw: cap.update(prep=kw["prep_file"], frcmod=kw["frcmod_file"])
        or {"success": True})

    r = _mgr(run_dir)._run_step_10_from_structure()

    assert r["success"] is True
    assert cap["prep"] == "cs1.prep"
    # the FINAL combined frcmod, never an intermediate
    assert cap["frcmod"] == "cs1.frcmod"


def test_step10_prefers_in_memory_results_over_disk(run_dir, monkeypatch):
    cap = {}
    monkeypatch.setattr(
        M, "integrate_modaa_from_structure",
        lambda **kw: cap.update(frcmod=kw["frcmod_file"]) or {"success": True})

    # A forward run has step 9's result in memory; it must win over the disk name.
    m = _mgr(run_dir, step_results={
        "step_8": {"prep_file": "cs1.prep"},
        "step_9": {"prep_file": "cs1.prep", "frcmod_file": "cs1_temp.frcmod"}})
    m._run_step_10_from_structure()
    assert cap["frcmod"] == "cs1_temp.frcmod"


def test_final_frcmod_recovery_ignores_intermediates(run_dir):
    # _temp / _gaff frcmods exist alongside; the resolver must pick only the final.
    got = _mgr(run_dir)._from_structure_frcmod_on_disk()
    assert got == "cs1.frcmod"


def test_ac_on_disk_helper_finds_and_misses(run_dir):
    m = _mgr(run_dir)
    assert m._from_structure_ac_on_disk() == "CS1.ac"
    os.remove(run_dir / "CS1.ac")
    assert m._from_structure_ac_on_disk() is None


def test_ac_to_mol2_recovers_ac_from_disk_on_resume(run_dir, monkeypatch):
    # Seminario runs inside step 9 and reaches back to step 7's AC. On resume
    # step_results is empty, so the connectivity MOL2 must be built from the
    # on-disk CS1.ac rather than aborting with "Could not build a MOL2".
    seen = {}

    def fake_run(cmd, **kw):
        seen["cmd"] = cmd
        # antechamber would write the MOL2; simulate that so the path resolves.
        out = cmd[cmd.index("-o") + 1]
        (run_dir / out).write_text("@<TRIPOS>MOLECULE\n")
        class _R:  # stand-in for CompletedProcess
            returncode = 0
        return _R()

    monkeypatch.setattr(M.subprocess, "run", fake_run)

    m = _mgr(run_dir)  # empty step_results = resumed run
    mol2 = m._ac_to_mol2("CS1")

    assert mol2 == "cs1.mol2"
    # It fed the recovered on-disk AC to antechamber, with amber atom types.
    assert seen["cmd"][seen["cmd"].index("-i") + 1] == "CS1.ac"
    assert seen["cmd"][seen["cmd"].index("-at") + 1] == "amber"


def test_ac_to_mol2_returns_none_when_no_ac_anywhere(run_dir, monkeypatch):
    os.remove(run_dir / "CS1.ac")
    called = {"run": False}
    monkeypatch.setattr(M.subprocess, "run",
                        lambda *a, **k: called.update(run=True))
    m = _mgr(run_dir)
    assert m._ac_to_mol2("CS1") is None
    assert called["run"] is False  # bailed before invoking antechamber


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-v"]))
