"""The shared torsion-refinement engine, and the two routes that call it.

Three defects motivated the engine and are pinned here:

  1. paramfit's WRITE_FRCMOD emits only the fitted terms, so the CREST route
     handed tLEaP an incomplete frcmod. The merge must produce a complete file
     and replace multi-term dihedrals whole, in either orientation (the old
     merge reversed the CHARACTER string, so ``c3-c3-os-c`` never matched
     its reverse).
  2. Each dihedral was fitted alone against the sm-6 prmtop. The engine pools
     every scan into one fit against a topology built from the current frcmod.
  3. Waiting for Gaussian ticked sm-7 green. Both routes now pause.
"""

import os
import shutil
import subprocess

import pytest
from rich.console import Console

from proprep.forcefield_prep import torsion_refinement as tr

QUIET = Console(quiet=True)


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------

def _gjf(path, *, flags=False, route="B3LYP/6-31G(d) Opt=ModRedundant Freq IOp(7/33=1) SCRF=(PCM,Solvent=Water)",
         tail=("D 1 2 3 4 F", "D 2 3 4 5 S 12 30.0"), chk=True):
    """A 5-atom Gaussian input in the style of either route."""
    with open(path, "w") as f:
        if chk:
            f.write("%chk=x.chk\n")
        f.write("%mem=8GB\n%nprocshared=4\n")
        f.write(f"#p {route}\n\n")
        f.write("title line\n\n0 1\n")
        for k, el in enumerate(["C", "C", "O", "C", "H"]):
            if flags:
                f.write(f"{el:<2s} {0:2d}  {k * 1.0:12.6f} {0.0:12.6f} {0.0:12.6f}\n")
            else:
                f.write(f"{el:<2s}    {k * 1.0:12.6f} {0.0:12.6f} {0.0:12.6f}\n")
        f.write("\n")
        for ln in tail:
            f.write(ln + "\n")
        f.write("\n")
    return path


def _mol2(path):
    """C1-C2-O3-C4 chain plus H5 on C4, typed c3 c3 os c3 hc."""
    atoms = [("C1", "c3"), ("C2", "c3"), ("O3", "os"), ("C4", "c3"), ("H5", "hc")]
    with open(path, "w") as f:
        f.write("@<TRIPOS>MOLECULE\nmol\n 5 4 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
        for i, (name, t) in enumerate(atoms, 1):
            f.write(f"{i:7d} {name:<4s} {i*1.0:10.4f} {0.0:10.4f} {0.0:10.4f} {t:<5s} 1 MOL {0.0:8.4f}\n")
        f.write("@<TRIPOS>BOND\n")
        for i, (a, b) in enumerate([(1, 2), (2, 3), (3, 4), (4, 5)], 1):
            f.write(f"{i:6d} {a:5d} {b:5d} 1\n")
    return path


def _frcmod(path, dihe_lines):
    with open(path, "w") as f:
        f.write("remark\nMASS\nc3 12.01\nos 16.00\n\nBOND\nc3-os  301.5  1.439\n\nANGLE\n"
                "c3-c3-os  67.7  108.4\n\nDIHE\n")
        for ln in dihe_lines:
            f.write(ln + "\n")
        f.write("\nIMPROPER\nc -c3-n -o   10.5  180.0  2.0\n\nNONBON\n\nEND\n")
    return path


def _scan(label, quad, energies, n_atoms=5, route="b3lyp/6-31g(d) opt=modredundant nosymm"):
    s = tr.TorsionScan(label=label, gjf=f"{label}.gjf", log=f"{label}.log",
                       atom_indices=(1, 2, 3, 4), quad=quad, route=route, n_atoms=n_atoms)
    s.geometries = [[(6, float(k), 0.0, 0.0)] * n_atoms for k in range(len(energies))]
    s.energies = list(energies)
    return s


# --------------------------------------------------------------------------
# quads and frcmod lines
# --------------------------------------------------------------------------

def test_quad_key_is_orientation_independent():
    assert tr.quad_matches(("c3", "c3", "os", "c"), ("c", "os", "c3", "c3"))
    assert not tr.quad_matches(("c3", "c3", "os", "c"), ("c3", "os", "c3", "c"))


def test_frcmod_dihe_types_reads_two_character_types_in_fixed_columns():
    line = "c3-c3-os-c    1    1.150       180.000           2.000    penalty score=  8.0"
    assert tr.frcmod_dihe_types(line) == ("c3", "c3", "os", "c")
    # the reversed orientation matches as TYPES, which the old character-reversal missed
    assert tr.quad_matches(tr.frcmod_dihe_types("c -os-c3-c3   1  0.0  0.0  1.0"),
                           ("c3", "c3", "os", "c"))


def test_quad_from_name_rejects_non_dihedrals():
    assert tr.quad_from_name("c3-os") is None
    assert tr.quad_from_name("c3-c3-os-c ") == ("c3", "c3", "os", "c")


# --------------------------------------------------------------------------
# deriving a scan input from an existing input
# --------------------------------------------------------------------------

def test_derive_scan_input_keeps_model_and_replaces_scan_line(tmp_path):
    template = _gjf(tmp_path / "opt.gjf", flags=True)
    out = tr.derive_scan_input(str(template), str(tmp_path / "scan.gjf"),
                               (1, 2, 3, 4), 24, 15.0,
                               coords=[(0, 0, k) for k in range(5)])
    text = open(out).read()
    route = [ln for ln in text.splitlines() if ln.startswith("#p")][0].lower()
    assert "freq" not in route and "iop(" not in route
    assert "opt=modredundant" in route and "nosymm" in route
    assert "scrf=(pcm,solvent=water)" in route          # solvent carried forward
    assert "%chk" not in text                            # no Hessian needed for a scan
    assert "%mem=8GB" in text
    assert "D 1 2 3 4 F" not in text                     # cannot freeze the scanned torsion
    assert "D 2 3 4 5 S 12 30.0" not in text             # old scan line dropped
    assert "D 1 2 3 4 S 24 15.0" in text
    atoms = [ln for ln in text.splitlines() if ln.startswith(("C ", "O ", "H "))]
    assert len(atoms) == 5
    assert atoms[4].split()[1] == "0"                    # freeze flag preserved
    assert float(atoms[4].split()[4]) == 4.0             # coords replaced
    assert tr.read_scan_line(out) == ((1, 2, 3, 4), 24, 15.0)


def test_derive_scan_input_adds_opt_when_template_is_single_point(tmp_path):
    template = _gjf(tmp_path / "sp.gjf", route="HF/6-31G* SP pop=mk iop(6/33=2)", tail=(), chk=False)
    out = tr.derive_scan_input(str(template), str(tmp_path / "scan.gjf"), (2, 3, 4, 5), 6, 60.0)
    route = [ln for ln in open(out).read().splitlines() if ln.startswith("#p")][0]
    assert route.startswith("#p Opt=ModRedundant HF/6-31G* SP")
    assert "pop=" not in route.lower()


def test_derive_scan_input_rejects_wrong_coordinate_count(tmp_path):
    template = _gjf(tmp_path / "opt.gjf")
    with pytest.raises(ValueError):
        tr.derive_scan_input(str(template), str(tmp_path / "s.gjf"), (1, 2, 3, 4), 24, 15.0,
                             coords=[(0, 0, 0)] * 4)


def test_scan_from_gjf_resolves_types_through_the_mol2(tmp_path):
    gjf = _gjf(tmp_path / "s.gjf")
    mol2 = _mol2(tmp_path / "m.mol2")
    scan = tr.scan_from_gjf(str(gjf), str(tmp_path / "s.log"), "s", mol2_file=str(mol2))
    assert scan.atom_indices == (2, 3, 4, 5) and scan.n_steps == 12 and scan.step_size == 30.0
    assert scan.quad == ("c3", "os", "c3", "hc")
    assert not scan.log_exists


def test_atoms_for_quad_matches_either_orientation(tmp_path):
    mol2 = _mol2(tmp_path / "m.mol2")
    assert tr.atoms_for_quad(("c3", "c3", "os", "c3"), str(mol2)) == (1, 2, 3, 4)
    assert tr.atoms_for_quad(("c3", "os", "c3", "c3"), str(mol2)) == (4, 3, 2, 1)
    assert tr.atoms_for_quad(("n", "c3", "os", "c3"), str(mol2)) is None


# --------------------------------------------------------------------------
# merging paramfit's fitted-terms-only output into a complete frcmod
# --------------------------------------------------------------------------

def test_merge_replaces_multi_term_dihedral_whole_and_keeps_everything_else(tmp_path):
    original = _frcmod(tmp_path / "orig.frcmod", [
        "c3-c3-os-c3   1    0.383       0.000          -3.000    penalty score= 12.0",
        "c3-c3-os-c3   1    0.100       180.000         2.000",
        "hc-c3-c3-os   1    0.156       0.000           3.000",
    ])
    fitted = tmp_path / "fitted.frcmod"
    fitted.write_text("Generated frcmod with paramfit\nDIHE\n"
                      "c3-os-c3-c3   1    0.911       0.000          -3.000\n"
                      "c3-os-c3-c3   1    0.222       180.000         2.000\n\n")
    out = tr.merge_fitted_dihedrals(str(original), str(fitted), [("c3", "c3", "os", "c3")],
                                    str(tmp_path / "merged.frcmod"), QUIET)
    text = open(out).read()
    assert "MASS\nc3 12.01" in text and "c3-os  301.5" in text    # complete file
    assert "0.383" not in text and "0.100" not in text            # both old terms gone
    assert "0.911" in text and "0.222" in text                    # both fitted terms in
    assert "hc-c3-c3-os   1    0.156" in text                     # unrelated term kept
    assert "IMPROPER\nc -c3-n -o" in text
    dihe = text.split("DIHE\n")[1].split("\nIMPROPER")[0]
    assert dihe.count("c3-os-c3-c3") == 2


def test_merge_leaves_frcmod_unchanged_when_nothing_was_fitted(tmp_path):
    original = _frcmod(tmp_path / "orig.frcmod", ["hc-c3-c3-os   1    0.156       0.000           3.000"])
    (tmp_path / "empty.frcmod").write_text("DIHE\n\n")
    out = tr.merge_fitted_dihedrals(str(original), str(tmp_path / "empty.frcmod"),
                                    [("a", "b", "c", "d")], str(tmp_path / "m.frcmod"), QUIET)
    assert open(out).read() == open(original).read()


# --------------------------------------------------------------------------
# pooling
# --------------------------------------------------------------------------

def test_check_compatible_refuses_mixed_levels_of_theory():
    a = _scan("a", ("c3", "c3", "os", "c3"), [-1.0, -1.1, -1.0, -0.9, -1.0])
    b = _scan("b", ("hc", "c3", "c3", "os"), [-1.0, -1.1, -1.0, -0.9, -1.0],
              route="wb97xd/def2tzvp opt=modredundant nosymm")
    ok, why = tr.check_compatible([a, b])
    assert not ok and "levels of theory" in why
    b.route = "b3lyp/6-31g(d) opt=modredundant nosymm geom=printinputorient"
    assert tr.check_compatible([a, b])[0]   # job-control tokens do not count


def test_check_compatible_refuses_different_atom_counts():
    a = _scan("a", ("c3", "c3", "os", "c3"), [-1.0] * 5)
    b = _scan("b", ("c3", "c3", "os", "c3"), [-1.0] * 5, n_atoms=6)
    assert not tr.check_compatible([a, b])[0]


def test_write_joint_inputs_pools_frames_in_order(tmp_path):
    a = _scan("a", ("c3", "c3", "os", "c3"), [-1.0, -1.1, -1.0, -0.9, -1.0])
    b = _scan("b", ("hc", "c3", "c3", "os"), [-2.0, -2.1, -2.0])
    mdcrd, energies, n, bounds = tr.write_joint_inputs([a, b], "mol", str(tmp_path), QUIET)
    assert n == 8 and bounds == [("a", 0, 5), ("b", 5, 8)]
    assert [float(x) for x in open(energies).read().split()] == a.energies + b.energies
    lines = open(mdcrd).read().splitlines()
    assert len(lines) == 1 + 8 * 2          # title + 15 coords per frame = 2 lines each


def test_per_scan_residuals_slices_by_boundary(tmp_path):
    f = tmp_path / "e.dat"
    f.write_text("# idx amber+K qm\n" + "".join(f"{i} {10.0 + r} 10.0\n"
                                                for i, r in enumerate([1, -1, 1, 3, 0, 0])))
    rows = tr.per_scan_residuals(str(f), [("a", 0, 3), ("b", 3, 6)])
    assert rows[0][0] == "a" and abs(rows[0][1] - 1.0) < 1e-9 and rows[0][2] == 1.0
    assert rows[1][0] == "b" and rows[1][2] == 3.0


def test_parse_selection():
    assert tr.parse_selection("none", [1, 2, 3]) == []
    assert tr.parse_selection("all", [2, 5]) == [2, 5]
    assert tr.parse_selection("1, 3-5, 9", [1, 2, 3, 4, 5]) == [1, 3, 4, 5]
    assert tr.parse_selection("x", [1]) is None


# --------------------------------------------------------------------------
# the joint fit: one K, one selection with every quad, one fit, one merge
# --------------------------------------------------------------------------

def test_run_joint_torsion_fit_fits_all_quads_at_once(tmp_path, monkeypatch):
    from proprep.forcefield_prep import paramfit_refinement as pf
    calls = {}
    monkeypatch.setattr(pf, "run_paramfit_k_fitting",
                        lambda prmtop, mdcrd, en, n, console: calls.update(k=(prmtop, n)) or 3.5)

    def fake_set(prmtop, params_file, selected, console, force_constants_only=True):
        calls["selected"] = [s[0] for s in selected]
        calls["fc_only"] = force_constants_only
        open(params_file, "w").write("x")
        return True

    def fake_fit(prmtop, mdcrd, en, params, n, k, out, console):
        calls["fit"] = (prmtop, n, k)
        with open(out, "w") as f:
            f.write("DIHE\nc3-c3-os-c3   1    0.5   0.0   3.0\nhc-c3-c3-os   1    0.2   0.0   3.0\n")
        with open(out.replace(".frcmod", "_energies.dat"), "w") as f:
            f.write("".join(f"{i} 1.0 1.0\n" for i in range(n)))
        return True

    monkeypatch.setattr(pf, "run_paramfit_set_params_automated", fake_set)
    monkeypatch.setattr(pf, "run_paramfit_parameter_fitting", fake_fit)
    monkeypatch.setattr(pf, "check_fitted_parameters", lambda f, c: True)

    frcmod = _frcmod(tmp_path / "cur.frcmod", [
        "c3-c3-os-c3   1    0.383       0.000           3.000    penalty score= 12.0",
        "hc-c3-c3-os   1    0.156       0.000           3.000    penalty score=  4.0",
        "c3-c3-c3-c3   1    0.180       0.000           3.000",
    ])
    a = _scan("a", ("c3", "c3", "os", "c3"), [-1.0, -1.1, -1.0, -0.9, -1.0])
    a2 = _scan("a2", ("c3", "os", "c3", "c3"), [-1.0, -1.1, -1.0, -0.9, -1.0])   # same quad, reversed
    b = _scan("b", ("hc", "c3", "c3", "os"), [-1.0, -1.1, -1.0, -0.9, -1.0])
    res = tr.run_joint_torsion_fit([a, a2, b], mol_name="mol", prmtop=str(tmp_path / "fit.parm7"),
                                   frcmod=str(frcmod), work_dir=str(tmp_path / "w"), console=QUIET)
    assert res["refinement_success"]
    assert res["quads"] == ["c3-c3-os-c3", "hc-c3-c3-os"]           # de-duplicated
    assert calls["k"][1] == 15 and calls["fit"][1] == 15               # all 3 scans pooled
    assert calls["fit"][2] == 3.5                                     # K passed through
    assert "c3-c3-os-c3" in calls["selected"] and "c3-os-c3-c3" in calls["selected"]
    assert "hc-c3-c3-os" in calls["selected"] and "so-3c-3c-ch" not in calls["selected"]
    assert calls["fc_only"] is False                                   # barrier + phase + n
    merged = open(res["merged_frcmod"]).read()
    assert "0.383" not in merged and "0.156" not in merged
    assert "0.5   0.0   3.0" in merged and "c3-c3-c3-c3   1    0.180" in merged
    assert [r[0] for r in res["residuals"]] == ["a", "a2", "b"]


# --------------------------------------------------------------------------
# tleap-backed topology (only where AmberTools is on PATH)
# --------------------------------------------------------------------------

_HAVE_AMBER = all(shutil.which(x) for x in ("tleap", "antechamber", "parmchk2"))


@pytest.mark.skipif(not _HAVE_AMBER, reason="AmberTools binaries not on PATH")
def test_build_prmtop_keeps_mol2_atom_order(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    pdb = tmp_path / "meoh.pdb"
    pdb.write_text(
        "HETATM    1  C1  MOH     1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  MOH     1       1.420   0.000   0.000  1.00  0.00           O\n"
        "HETATM    3  H1  MOH     1      -0.390   1.030   0.000  1.00  0.00           H\n"
        "HETATM    4  H2  MOH     1      -0.390  -0.515   0.892  1.00  0.00           H\n"
        "HETATM    5  H3  MOH     1      -0.390  -0.515  -0.892  1.00  0.00           H\n"
        "HETATM    6  H4  MOH     1       1.740   0.900   0.000  1.00  0.00           H\n"
        "END\n")
    subprocess.run(["antechamber", "-i", str(pdb), "-fi", "pdb", "-o", "meoh.mol2", "-fo", "mol2",
                    "-at", "gaff2", "-pf", "y"], check=True, capture_output=True)
    subprocess.run(["parmchk2", "-i", "meoh.mol2", "-f", "mol2", "-o", "meoh.frcmod", "-s", "2"],
                   check=True, capture_output=True)
    prmtop = tr.build_gaff2_prmtop("meoh.mol2", "meoh.frcmod", str(tmp_path / "fit"), QUIET)
    assert prmtop and os.path.exists(prmtop)
    import parmed
    names = [a.name for a in parmed.load_file(prmtop).atoms]
    assert names == ["C1", "O1", "H1", "H2", "H3", "H4"]
    # a Gaussian input listing the same atoms passes the order guard
    gjf = tmp_path / "meoh.gjf"
    gjf.write_text("#p HF/6-31G* opt\n\nt\n\n0 1\n" + "".join(
        f"{el}  0.0 0.0 {k}.0\n" for k, el in enumerate("COHHHH")) + "\n")
    assert tr.check_atom_order("meoh.mol2", str(gjf)) == (True, "")
    bad = tmp_path / "bad.gjf"
    bad.write_text("#p HF/6-31G* opt\n\nt\n\n0 1\n" + "".join(
        f"{el}  0.0 0.0 {k}.0\n" for k, el in enumerate("OCHHHH")) + "\n")
    assert not tr.check_atom_order("meoh.mol2", str(bad))[0]


# --------------------------------------------------------------------------
# small-molecule step sm-7
# --------------------------------------------------------------------------

def _runner(tmp_path):
    from proprep.forcefield_prep.small_molecule_parameterizer import SmallMolWorkflowRunner
    r = SmallMolWorkflowRunner.__new__(SmallMolWorkflowRunner)
    r.console = QUIET
    r.processor = None
    r.interactive = False
    r.mol_name = "LIG"
    r.mol_name_lower = "lig"
    r.net_charge, r.multiplicity = 0, 1
    r.results = {"parameter_files": {}}
    r.penalties = []
    r.refinement_config = {}
    r.gaussian_settings = None
    r.gaussian_file = r.gaussian_log_file = None
    r.mol2_file = r.frcmod_file = r.chk_file = r.prmtop_file = r.current_frcmod = None
    r.mol_pdb_file = None
    r.success_files = []
    return r


def test_refinement_selection_survives_a_fresh_session(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    r = _runner(tmp_path)
    r.frcmod_file = _frcmod(tmp_path / "lig.frcmod", [
        "c3-c3-os-c3   1    0.383       0.000           3.000    penalty score= 12.0",
        "hc-c3-c3-os   1    0.156       0.000           3.000    penalty score=  4.0",
    ])
    from proprep.forcefield_prep.small_molecule_parameterizer import analyze_frcmod_penalties
    r.penalties = analyze_frcmod_penalties(r.frcmod_file, QUIET)
    dihe = [i for i, p in enumerate(r.penalties) if p[3] == "DIHE"]
    r.refinement_config = {"bonds_angles": [], "dihedrals": dihe, "dihedral_method": "pes"}
    r._save_refinement_selection()

    fresh = _runner(tmp_path)
    fresh.frcmod_file = "lig.frcmod"
    fresh._load_refinement_selection()
    assert fresh.refinement_config["dihedral_method"] == "pes"
    assert sorted(fresh.penalties[i][0] for i in fresh.refinement_config["dihedrals"]) == \
        sorted(r.penalties[i][0] for i in dihe)


def test_sm7_writes_scans_and_reports_pending_when_no_logs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    r = _runner(tmp_path)
    r.mol2_file = str(_mol2(tmp_path / "lig.mol2"))
    r.frcmod_file = r.current_frcmod = str(_frcmod(tmp_path / "lig.frcmod", [
        "c3-c3-os-c3   1    0.383       0.000           3.000    penalty score= 12.0"]))
    r.gaussian_file = str(_gjf(tmp_path / "lig.gjf", tail=(), chk=True))
    outcome = r._refine_dihedrals_pes([("c3-c3-os-c3", 12.0, "penalty", "DIHE")])
    assert outcome == "pending"
    gjf = tmp_path / "lig_pes_scans" / "lig_c3_c3_os_c3_scan.gjf"
    assert gjf.exists()
    assert tr.read_scan_line(str(gjf)) == ((1, 2, 3, 4), 24, 15.0)


def test_sm7_returns_checkpoint_while_scans_are_outstanding(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    r = _runner(tmp_path)
    r.frcmod_file = str(_frcmod(tmp_path / "lig.frcmod", [
        "c3-c3-os-c3   1    0.383       0.000           3.000    penalty score= 12.0"]))
    from proprep.forcefield_prep.small_molecule_parameterizer import analyze_frcmod_penalties
    r.penalties = analyze_frcmod_penalties(r.frcmod_file, QUIET)
    r.refinement_config = {"bonds_angles": [], "dihedrals": [0], "dihedral_method": "pes"}
    monkeypatch.setattr(r, "_refine_dihedrals_pes", lambda params: "pending")
    assert r._step_refinement() == {"checkpoint": True}


def test_sm7_step_is_declared_a_checkpoint():
    from proprep.forcefield_prep.small_molecule_parameterizer import SMALL_MOL_WORKFLOW_STEPS
    step = next(s for s in SMALL_MOL_WORKFLOW_STEPS if s.id == "sm-7")
    assert step.checkpoint and step.optional


def test_sm7_crest_result_is_merged_not_used_raw(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    r = _runner(tmp_path)
    r.mol2_file = str(_mol2(tmp_path / "lig.mol2"))
    r.frcmod_file = r.current_frcmod = str(_frcmod(tmp_path / "lig.frcmod", [
        "c3-c3-os-c3   1    0.383       0.000           3.000    penalty score= 12.0",
        "c3-c3-c3-c3   1    0.180       0.000           3.000"]))
    monkeypatch.setattr(r, "_fit_prmtop", lambda d: str(tmp_path / "fit.parm7"))
    fitted = tmp_path / "lig_paramfit" / "lig_fitted.frcmod"
    fitted.parent.mkdir()
    fitted.write_text("DIHE\nc3-c3-os-c3   1    0.9   0.0   3.0\n\n")
    import proprep.forcefield_prep.paramfit_refinement as pf
    monkeypatch.setattr(pf, "run_paramfit_refinement_workflow",
                        lambda **kw: {"refinement_success": True, "fitted_frcmod": str(fitted), "message": ""})
    assert r._refine_dihedrals_crest([("c3-c3-os-c3", 12.0, "penalty", "DIHE")]) == "done"
    merged = open(r.current_frcmod).read()
    assert "MASS" in merged and "c3-c3-c3-c3   1    0.180" in merged and "0.9   0.0   3.0" in merged
    assert "0.383" not in merged


# --------------------------------------------------------------------------
# modified amino acid step 9 (both routes)
# --------------------------------------------------------------------------

def _manager(tmp_path, mode="from_structure"):
    from proprep.forcefield_prep.modified_amino_acid_parameterizer import ModifiedAAWorkflowManager
    m = ModifiedAAWorkflowManager.__new__(ModifiedAAWorkflowManager)
    m.console = QUIET
    m.processor = None
    m.amino_acid = "CS1"
    m.conformer_mode = mode
    m.conformers = ["xtal"] if mode == "from_structure" else ["ahelix", "bsheet"]
    m.step_results = {}
    m._conformer_pdb_map = {}
    return m


def test_modaa_refinement_is_skipped_in_batch_mode(tmp_path):
    assert _manager(tmp_path)._maybe_dihedral_refinement("CS1", "x.frcmod", "x.prep",
                                                         interactive=False) == (None, False)


def test_modaa_pools_existing_scan_and_added_scan_then_pauses(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    m = _manager(tmp_path)
    mol2 = _mol2(tmp_path / "cs1.mol2")
    monkeypatch.setattr(m, "_ac_to_mol2", lambda rs: str(mol2))
    # the step-3 linkage scan (log present) drives c3-c3-os-c3 = atoms 1-2-3-4
    _gjf(tmp_path / "cs1_xtal_scan.gjf", flags=True, tail=("D 1 2 3 4 S 24 15.0",))
    (tmp_path / "cs1_xtal_scan.log").write_text("")
    frcmod = _frcmod(tmp_path / "cs1.frcmod", [
        "c3-c3-os-c3   1    0.383       0.000           3.000    penalty score= 12.0",
        "c3-os-c3-hc   1    0.156       0.000           3.000    penalty score=  4.0",
    ])
    import proprep.forcefield_prep.modified_amino_acid_parameterizer as M
    answers = iter(["all"])
    monkeypatch.setattr(M, "confirm_with_context", lambda *a, **k: True)
    monkeypatch.setattr(M, "prompt_with_context", lambda *a, **k: next(answers, "15.0"))
    monkeypatch.setattr(M, "int_prompt_with_context", lambda *a, **k: 12)
    collected = []
    monkeypatch.setattr(tr, "collect_scan", lambda sc, console, min_points=5: collected.append(sc.label) or False)

    refit, pending = m._maybe_dihedral_refinement("CS1", str(frcmod), "cs1.prep", interactive=True)
    assert refit is None and pending is True
    added = tmp_path / "torsion_refit" / "cs1_xtal_scan_c3_os_c3_hc.gjf"
    assert added.exists()
    assert tr.read_scan_line(str(added)) == ((2, 3, 4, 5), 12, 15.0)
    # the added scan was derived from the step-3 scan input: same restraint style, new S line
    assert "D 1 2 3 4 S 24 15.0" not in added.read_text()
    assert "step-3 scan, xtal" in collected     # the existing scan was pooled, not just the new one


def test_modaa_step9_pauses_both_routes_when_scans_are_pending(tmp_path, monkeypatch):
    import proprep.forcefield_prep.modified_amino_acid_parameterizer as M
    monkeypatch.setattr(M, "generate_bonded_parameters",
                        lambda rs, proc, standalone_use=True: {"success": True, "final_frcmod": "f", "prep_file": "p"})
    for mode in ("from_structure", "de_novo"):
        m = _manager(tmp_path, mode)
        m.step_results["step_8"] = {"residue_symbol": "CS1"}
        monkeypatch.setattr(m, "_maybe_seminario", lambda *a, **k: None, raising=False)
        monkeypatch.setattr(m, "_maybe_dihedral_refinement", lambda *a, **k: (None, True))
        r = m._run_step_9(interactive=True)
        assert r["status"] == M.ModifiedAAWorkflowManager.PAUSE_STATUS, mode
        assert m._checklist_aa_9_parmchk2() == {"checkpoint": True}, mode
        assert "step_9" not in m.step_results


def test_modaa_step9_uses_the_refit_frcmod_on_both_routes(tmp_path, monkeypatch):
    import proprep.forcefield_prep.modified_amino_acid_parameterizer as M
    monkeypatch.setattr(M, "generate_bonded_parameters",
                        lambda rs, proc, standalone_use=True: {"success": True, "final_frcmod": "f", "prep_file": "p"})
    for mode in ("from_structure", "de_novo"):
        m = _manager(tmp_path, mode)
        m.step_results["step_8"] = {"residue_symbol": "CS1"}
        monkeypatch.setattr(m, "_maybe_seminario", lambda *a, **k: None, raising=False)
        monkeypatch.setattr(m, "_maybe_dihedral_refinement", lambda *a, **k: ("refined.frcmod", False))
        assert m._run_step_9(interactive=True)["frcmod_file"] == "refined.frcmod", mode


def test_aa9_is_declared_a_checkpoint_in_both_step_lists():
    from proprep.forcefield_prep.modified_amino_acid_parameterizer import (
        MODIFIED_AA_STEPS, MODIFIED_AA_FROM_STRUCTURE_STEPS)
    for steps in (MODIFIED_AA_STEPS, MODIFIED_AA_FROM_STRUCTURE_STEPS):
        assert next(s for s in steps if s.id == "aa-9").checkpoint


# --------------------------------------------------------------------------
# real paramfit, end to end (only where AmberTools is on PATH; ~1 s)
# --------------------------------------------------------------------------

_HAVE_PARAMFIT = _HAVE_AMBER and shutil.which("paramfit")


@pytest.mark.skipif(not _HAVE_PARAMFIT, reason="paramfit/AmberTools not on PATH")
def test_real_paramfit_joint_fit_keeps_periodicity_and_builds(tmp_path, monkeypatch):
    """K fit → pty-driven SET_PARAMS → genetic fit → merge, on a synthetic
    methanol H-C-O-H scan. Pins two things a mocked fit cannot: paramfit's
    prompts are answered so that periodicity stays an integer, and the
    merged frcmod is complete enough for tleap to build from."""
    import math
    import numpy as np
    monkeypatch.chdir(tmp_path)
    base = np.array([[0.0, 0.0, 0.0], [1.42, 0.0, 0.0], [-0.39, 1.03, 0.0],
                     [-0.39, -0.515, 0.892], [-0.39, -0.515, -0.892], [1.74, 0.90, 0.0]])
    els = ["C", "O", "H", "H", "H", "H"]
    with open("meoh.pdb", "w") as f:
        for i, (e, (x, y, z)) in enumerate(zip(els, base), 1):
            f.write(f"HETATM{i:5d} {e}{i:<3d} MOH     1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {e:>2s}\n")
        f.write("END\n")
    subprocess.run(["antechamber", "-i", "meoh.pdb", "-fi", "pdb", "-o", "meoh.mol2", "-fo", "mol2",
                    "-at", "gaff2", "-c", "gas", "-pf", "y"], check=True, capture_output=True)
    subprocess.run(["parmchk2", "-i", "meoh.mol2", "-f", "mol2", "-o", "meoh.frcmod", "-s", "2", "-a", "Y"],
                   check=True, capture_output=True)

    def rot(p, o, d, ang):
        d = d / np.linalg.norm(d)
        v = p - o
        return o + v * math.cos(ang) + np.cross(d, v) * math.sin(ang) + d * np.dot(d, v) * (1 - math.cos(ang))

    frames, energies = [], []
    for k in range(12):
        ang = math.radians(30 * k)
        g = base.copy()
        g[5] = rot(base[5], base[1], base[1] - base[0], ang)
        frames.append([({"C": 6, "O": 8, "H": 1}[e], *map(float, xyz)) for e, xyz in zip(els, g)])
        energies.append(-115.0 + 0.6 * (1 + math.cos(3 * ang)) / 627.5095)
    with open("scan.gjf", "w") as f:
        f.write("%mem=2GB\n#p HF/6-31G* Opt=ModRedundant NoSymm\n\nscan\n\n0 1\n")
        for e, xyz in zip(els, base):
            f.write(f"{e}  {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n")
        f.write("\nD 3 1 2 6 S 12 30.0\n\n")
    scan = tr.scan_from_gjf("scan.gjf", "scan.log", "h-c-o-h", mol2_file="meoh.mol2")
    scan.geometries, scan.energies, scan.n_atoms = frames, energies, 6
    scan.route = "hf/6-31g* opt=modredundant nosymm"
    assert scan.quad == ("h1", "c3", "oh", "ho")

    prmtop = tr.build_gaff2_prmtop("meoh.mol2", "meoh.frcmod", str(tmp_path / "fit"), QUIET)
    res = tr.run_joint_torsion_fit([scan], mol_name="meoh", prmtop=prmtop, frcmod="meoh.frcmod",
                                   work_dir=str(tmp_path), console=QUIET, fit_equilibrium=True)
    assert res["refinement_success"], res["message"]
    merged = open(res["merged_frcmod"]).read()
    line = next(ln for ln in merged.splitlines() if tr.frcmod_dihe_types(ln) == ("h1", "c3", "oh", "ho"))
    idivf, barrier, phase, period = line[11:].split()[:4]
    assert float(period) == 3.0                      # periodicity held, not fitted
    assert float(barrier) > 0
    assert "MASS" in merged and "NONBON" in merged   # complete file, not paramfit's fragment
    assert tr.build_gaff2_prmtop("meoh.mol2", res["merged_frcmod"], str(tmp_path / "fit2"), QUIET)
    assert res["residuals"] and res["residuals"][0][1] < 1.0


def _refit_run(tmp_path, monkeypatch, conflict):
    """Run the modAA refinement with one existing (unpenalized) scan and a stubbed collision check."""
    monkeypatch.chdir(tmp_path)
    m = _manager(tmp_path)
    mol2 = _mol2(tmp_path / "cs1.mol2")
    monkeypatch.setattr(m, "_ac_to_mol2", lambda rs: str(mol2))
    _gjf(tmp_path / "cs1_xtal_scan.gjf", flags=True, tail=("D 1 2 3 4 S 24 15.0",))
    (tmp_path / "cs1_xtal_scan.log").write_text("")
    frcmod = _frcmod(tmp_path / "cs1.frcmod", [
        "c3-c3-os-c3   1    0.383       0.000           3.000",          # no penalty: wildcard-covered
        "c3-os-c3-hc   1    0.156       0.000           3.000    penalty score=  4.0",
    ])
    import proprep.forcefield_prep.modified_amino_acid_parameterizer as M
    from proprep.forcefield_prep import shared_type_terms as st
    monkeypatch.setattr(M, "confirm_with_context", lambda *a, **k: True)
    monkeypatch.setattr(M, "prompt_with_context", lambda *a, **k: "none")
    monkeypatch.setattr(st, "shared_type_dihedral_conflict", lambda q: conflict)
    monkeypatch.setattr(st, "amber_data_available", lambda: True)
    pooled = []
    monkeypatch.setattr(tr, "collect_scan", lambda sc, console, min_points=5: pooled.append(sc.label) or False)
    return m._maybe_dihedral_refinement("CS1", str(frcmod), "cs1.prep", interactive=True), pooled


def test_modaa_does_not_refit_a_quad_shared_with_the_protein_ff(tmp_path, monkeypatch):
    """A quad the force field defines, or a standard residue contains, is shown but not fitted:
    the refit would be written under shared types and reach every residue using it."""
    (refit, pending), pooled = _refit_run(tmp_path, monkeypatch,
                                          conflict="the protein force field defines it explicitly")
    assert pooled == []
    assert (refit, pending) == (None, False)


def test_modaa_refits_an_unpenalized_quad_that_is_local(tmp_path, monkeypatch):
    """No penalty does NOT mean shared: a wildcard-covered quad absent from every
    standard residue (the covalent-linkage case) is local and may be refit."""
    (refit, pending), pooled = _refit_run(tmp_path, monkeypatch, conflict=None)
    assert pooled == ["step-3 scan, xtal"]      # pooled; then pending because collect_scan stub says no data
    assert pending is True
