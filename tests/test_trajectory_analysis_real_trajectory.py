"""Every MD Manager trajectory analysis runs on a real solvated trajectory.

Seen in 1.20.0: of the 22 analyses, 5 ran. The rest called pytraj functions
that do not exist (pt.cluster.hierarchical, pt.radial, pt.calc_chi), passed
keywords pytraj does not take, read tuples as arrays, or called pytraj on bare
Frames; nine error handlers then died on an undefined `logger`. Others ran but
were wrong: RMSD/PCA/RMSF rotated the shared trajectory in place (changing
every later imaged analysis), "average structure" was the last frame, DSSP on
the default C-alpha selection saw no H-bonds, PCA's variance always summed to
100%, and phi of one residue was paired with psi of its neighbour.

None of that can be caught with a mocked analyzer, so these tests use the
solvated trpzip2 trajectory that ships with pytraj, and check results against
an independent calculation rather than only checking that nothing raises.
"""

from __future__ import annotations

import ast
import io
import types
from pathlib import Path

import numpy as np
import pytest
from rich.console import Console
from scipy.spatial.distance import cdist

pt = pytest.importorskip("pytraj")

from proprep.md_prep import molecular_dynamics_manager as mdm_module  # noqa: E402
from proprep.utils import prompts as prompts_module  # noqa: E402
from proprep.md_prep.molecular_dynamics_manager import MolecularDynamicsManager as MDM  # noqa: E402
from proprep.md_prep.trajectory_analyzer import TrajectoryAnalyzer  # noqa: E402

TOP, NC = pt.datafiles.tz2_ortho_parm7, pt.datafiles.tz2_ortho_nc
PROTEIN = "1-12"  # SWTWENGKWTWK; residue 13 is the NHE cap (no CA), the rest is water


@pytest.fixture()
def analyzer():
    return TrajectoryAnalyzer(TOP, NC)


@pytest.fixture(scope="module")
def reference():
    return pt.iterload(NC, top=TOP)  # immutable, independent of the analyzer


# ---------------------------------------------------------------- analyzer

def test_protein_residues_come_from_the_topology(analyzer):
    assert analyzer.get_protein_residue_range() == PROTEIN
    assert TrajectoryAnalyzer._residue_ids_to_range([1, 2, 3, 7, 9, 10]) == "1-3,7,9-10"


def test_non_string_mask_raises_instead_of_segfaulting(analyzer):
    with pytest.raises(TypeError):
        analyzer.calculate_contact_frequency(mask={"mask": "@CA", "description": "C-alpha atoms"})


def test_no_analysis_moves_the_shared_trajectory(analyzer):
    before = analyzer.traj.xyz.copy()
    analyzer.calculate_rmsd("@CA", 0)
    analyzer.calculate_rmsd("@CA", -1)
    analyzer.calculate_rmsf("@CA")
    analyzer.calculate_bfactors("@CA")
    analyzer.calculate_pca("@CA", 3)
    analyzer.calculate_clustering("@CA", 3, "kmeans")
    analyzer.calculate_pairwise_rmsd("@CA")
    assert np.array_equal(before, analyzer.traj.xyz)


def test_rmsd_matches_pytraj_and_average_reference_is_not_the_last_frame(analyzer, reference):
    to_first = analyzer.calculate_rmsd("@CA", 0)
    assert np.allclose(to_first, pt.rmsd(reference, "@CA", ref=0))

    to_average = analyzer.calculate_rmsd("@CA", -1)
    to_last = pt.rmsd(reference, "@CA", ref=reference.n_frames - 1)
    assert to_average.min() > 0.01          # a frame is never identical to the average
    assert not np.allclose(to_average, to_last)


def test_bfactors_follow_from_fitted_rmsf(analyzer):
    result = analyzer.calculate_bfactors("@CA")
    assert result["residue_ids"] == list(range(1, 13))
    assert np.allclose(result["bfactors"], 8 * np.pi ** 2 / 3 * result["rmsf"] ** 2)
    assert np.allclose(np.sort(result["rmsf"]), np.sort(analyzer.calculate_rmsf("@CA")), rtol=0.35)


def test_contacts_native_count_and_q(analyzer, reference):
    result = analyzer.calculate_contacts("@CA", "@CA", 6.0, 0)
    ca = reference.top.select("@CA")
    first = cdist(reference[0].xyz[ca], reference[0].xyz[ca])
    assert result["n_native_contacts"] == int(((first <= 6.0) & ~np.eye(len(ca), dtype=bool)).sum())
    assert result["q_values"][0] == 1.0
    assert all(0.0 <= q <= 1.0 for q in result["q_values"])


def test_contact_frequency_agrees_with_contact_map(analyzer):
    contacts = analyzer.calculate_contacts("@CA", "@CA", 6.0, 0)
    frequency = analyzer.calculate_contact_frequency("@CA", 6.0)
    assert np.isclose(sum(frequency["frequencies"]), contacts["contact_frequency"].sum())


def test_salt_bridges_find_the_glu_lys_pairs(analyzer, reference):
    result = analyzer.calculate_salt_bridges(distance_cutoff=30.0)  # wide: every GLU-LYS pair forms
    assert {b["pair"] for b in result["salt_bridges"]} == {"GLU5-LYS8", "GLU5-LYS12"}
    glu5_lys8 = next(b for b in result["salt_bridges"] if b["pair"] == "GLU5-LYS8")
    direct = np.min([pt.distance(reference, f":5@{a} :8@NZ") for a in ("OE1", "OE2")], axis=0)
    assert np.allclose(glu5_lys8["distances"], direct, atol=1e-3)


def test_pca_variance_is_relative_to_all_motion(analyzer, reference):
    result = analyzer.calculate_pca("@CA", 3)
    _, (all_eigenvalues, _) = pt.pca(reference["@CA"], mask="*", n_vecs=36)
    assert np.allclose(result["variance_explained"], all_eigenvalues[:3] / all_eigenvalues.sum() * 100, rtol=1e-3)
    assert result["cumulative_variance"][-1] < 100.0
    assert len(result["projections"]["PC1"]) == reference.n_frames


@pytest.mark.parametrize("algorithm", ["kmeans", "hierarchical"])
def test_cluster_representative_is_the_medoid_of_its_own_cluster(analyzer, reference, algorithm):
    result = analyzer.calculate_clustering("@CA", 3, algorithm)
    assignments, representatives = result["assignments"], result["representatives"]
    assert sum(result["populations"].values()) == reference.n_frames

    pairwise = pt.pairwise_rmsd(reference, mask="@CA")
    for cluster_id in range(result["n_clusters"]):
        members = np.where(assignments == cluster_id)[0]
        medoid = members[np.argmin(pairwise[np.ix_(members, members)].sum(axis=1))]
        assert assignments[representatives[cluster_id]] == cluster_id
        assert representatives[cluster_id] == medoid


def test_pairwise_rmsd_matrix(analyzer, reference):
    result = analyzer.calculate_pairwise_rmsd("@CA")
    matrix = result["matrix"]
    assert matrix.shape == (reference.n_frames, reference.n_frames)
    assert np.allclose(matrix, matrix.T) and np.allclose(np.diag(matrix), 0.0)
    assert np.allclose(matrix[0], pt.rmsd(reference, "@CA", ref=0), atol=1e-3)
    assert np.allclose(analyzer.calculate_pairwise_rmsd("@CA", subsample=2)["matrix"], matrix[::2, ::2], atol=1e-3)


def test_sasa(analyzer, reference):
    sasa = analyzer.calculate_sasa(f":{PROTEIN}", 1.4)
    assert np.allclose(sasa, pt.molsurf(reference, f":{PROTEIN}", probe=1.4))


def test_dssp_widens_an_atom_selection_to_whole_residues(analyzer):
    from_ca = analyzer.calculate_dssp("@CA")
    assert from_ca["mask"] == f":{PROTEIN}"
    assert np.array_equal(from_ca["per_residue"], analyzer.calculate_dssp("*")["per_residue"])
    assert np.mean(from_ca["sheet_pct"]) > 0            # trpzip2 is a beta hairpin
    total = np.array(from_ca["helix_pct"]) + from_ca["sheet_pct"] + from_ca["turn_pct"] + from_ca["coil_pct"]
    assert np.allclose(total, 100.0)


def test_dssp_comes_from_the_cpptraj_program_never_from_pytraj(analyzer, reference, monkeypatch):
    """Found by the 1.21.0 installer test on Rocky: pt.dssp corrupts the heap on the
    linux-64 pytraj, and a LATER pytraj call aborts the whole process (DSSP then
    other analyses: 11 of 20 sessions died mid-run; without DSSP, 0 of 20). The
    cpptraj executable runs the same algorithm: 40 of 40 clean, identical codes.
    cpptraj sees a stripped copy numbered from 1, so labels must come from the
    real topology, and a partial selection is assigned with its partners present."""
    real_dssp = pt.dssp
    monkeypatch.setattr(pt, "dssp", lambda *a, **k: pytest.fail("pt.dssp must not run in the ProPrep process"))

    part = analyzer.calculate_dssp(":3-8", label="part")
    assert part["residue_labels"] == ["THR:3", "TRP:4", "GLU:5", "ASN:6", "GLY:7", "LYS:8"] and part["mask"] == ":3-8"

    whole = analyzer.calculate_dssp("*", label="whole")
    assert whole["per_residue"].shape == (reference.n_frames, 12)
    assert np.array_equal(part["per_residue"], whole["per_residue"][:, 2:8])       # same residues, same answer
    # the same codes pytraj gives (it is safe to ask it here: macOS never showed the fault, and this is a test process)
    assert np.array_equal(whole["per_residue"], real_dssp(reference, mask=f":{PROTEIN}")[1])


def test_dssp_reports_a_cpptraj_failure(analyzer, monkeypatch):
    import subprocess
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 1, stdout="", stderr="Error: boom"))
    with pytest.raises(RuntimeError, match="cpptraj secstruct failed"):
        analyzer.calculate_dssp("@CA")


def test_ramachandran_pairs_phi_and_psi_of_the_same_residue(analyzer, reference):
    result = analyzer.calculate_ramachandran()
    assert result["residues"] == list(range(2, 13))     # residue 1 has no phi
    row = result["residues"].index(6)
    assert np.allclose(result["phi"][row], pt.calc_phi(reference, resrange="6")[0].values)
    assert np.allclose(result["psi"][row], pt.calc_psi(reference, resrange="6")[0].values)


def test_chi1_and_omega(analyzer):
    chi1 = analyzer.calculate_ramachandran(None, "chi1")
    assert chi1["chi"].shape == (len(chi1["residues"]), 10) and 7 not in chi1["residues"]  # GLY has no chi1
    omega = analyzer.calculate_ramachandran(None, "omega")
    assert abs(np.abs(omega["omega"]).mean() - 180.0) < 20.0


def test_rdf(analyzer, reference):
    distances, g = analyzer.calculate_water_radial_distribution(f":{PROTEIN}&!@H=", ":WAT@O", 8.0, 0.5)
    expected = pt.rdf(reference, solvent_mask=":WAT@O", solute_mask=f":{PROTEIN}&!@H=", maximum=8.0, bin_spacing=0.5)
    assert np.allclose(distances, expected[0]) and np.allclose(g, expected[1])
    with pytest.raises(ValueError):
        analyzer.calculate_water_radial_distribution(":NOPE", ":WAT@O")


@pytest.mark.parametrize("shell_width, max_distance, n_shells", [(2.0, 8.0, 4), (3.0, 9.0, 3), (5.0, 5.0, 1)])
def test_water_shells_are_imaged_cumulative_counts(analyzer, reference, shell_width, max_distance, n_shells):
    result = analyzer.calculate_water_shells(None, shell_width, max_distance)
    assert result["solute_mask"] == f":{PROTEIN}&!@H=" and result["n_shells"] == n_shells
    populations = result["populations"]
    assert (populations >= 0).all()
    within = pt.watershell(reference, solute_mask=f":{PROTEIN}&!@H=", solvent_mask=":WAT@O",
                           lower=max_distance / 2, upper=max_distance, dtype="ndarray")[1]
    assert np.allclose(populations.sum(axis=1), within)


def test_density_map_accounts_for_every_position(analyzer):
    result = analyzer.calculate_density_map(":WAT@O", 2.0)
    assert np.isclose(result["density"].sum(), 1.0)
    assert result["density_slice"].shape == result["grid_dims"][:2]


def test_density_grid_limit_is_honoured_and_reported(analyzer):
    result = analyzer.calculate_density_map(":WAT@O", 1.0, max_grid_points=10, n_peaks=3)
    assert max(result["grid_dims"]) <= 10
    assert result["requested_spacing"] == 1.0 and result["grid_spacing"] > 1.0
    assert len(result["peaks"]) == 3 and np.isclose(result["density"].sum(), 1.0)


def test_density_solute_frame_differs_from_lab_frame_and_works_on_a_copy(analyzer):
    before = analyzer.traj.xyz.copy()
    lab = analyzer.calculate_density_map(":WAT@O", 1.0, label="lab")
    solute = analyzer.calculate_density_map(":WAT@O", 1.0, fit_mask="@CA,C,N", label="solute")
    assert solute["fit_mask"] == "@CA,C,N" and lab["fit_mask"] is None
    assert np.isclose(solute["density"].sum(), 1.0)
    assert not np.allclose(lab["grid_min"], solute["grid_min"])     # the fit really moved the waters it grids
    assert np.array_equal(before, analyzer.traj.xyz)                # ...but on a copy, not the shared trajectory


def test_exposed_method_choices_change_the_result(analyzer):
    """Each new prompt has to reach pytraj, not just be recorded."""
    assert not np.allclose(analyzer.calculate_pairwise_rmsd("@CA", metric="rms")["matrix"],
                           analyzer.calculate_pairwise_rmsd("@CA", metric="dme")["matrix"])
    fitted, unfitted = analyzer.calculate_pca("@CA", 2, fit=True), analyzer.calculate_pca("@CA", 2, fit=False)
    assert not np.allclose(fitted["eigenvalues"], unfitted["eigenvalues"])
    assert unfitted["cumulative_variance"][-1] < 100.0
    on_backbone = analyzer.calculate_rmsf("@CA", alignment_mask="@CA,C,N", reference=0)
    on_one_residue = analyzer.calculate_rmsf("@CA", alignment_mask=":6", reference=9, label="other")
    assert not np.allclose(on_backbone, on_one_residue)
    assert len(analyzer.calculate_contacts("@CA", "@CA", 6.0, 0, persistence_threshold=0.0)["persistent_contacts"]) \
        >= len(analyzer.calculate_contacts("@CA", "@CA", 6.0, 0, persistence_threshold=0.9)["persistent_contacts"])
    with pytest.raises(ValueError):
        analyzer.calculate_rmsf("@CA", alignment_mask=":NOPE")
    with pytest.raises(ValueError):
        analyzer.calculate_clustering("@CA", 3, "hierarchical", linkage="ward")


def test_vector_length_is_the_distance(analyzer, reference):
    result = analyzer.calculate_vector_analysis(":1@CA", ":12@CA")
    assert np.allclose(result["magnitudes"], pt.distance(reference, ":1@CA :12@CA"), atol=1e-3)
    assert ((0 <= result["polar_angles"]) & (result["polar_angles"] <= 180)).all()
    with pytest.raises(ValueError):
        analyzer.calculate_vector_analysis(":1@CA", ":400@CA")


def test_imaged_analysis_does_not_depend_on_what_ran_before(analyzer):
    fresh = analyzer.calculate_water_shells(None, 2.0, 6.0)["populations"]
    analyzer.calculate_rmsd("@CA", 0)
    analyzer.calculate_pca("@CA", 2)
    analyzer.calculate_rmsf("@CA")
    assert np.array_equal(fresh, analyzer.calculate_water_shells(None, 2.0, 6.0, label="again")["populations"])


# ---------------------------------------------------------------- MD Manager menu

MASKS = [("first atom", ":1@CA"), ("second atom", ":5@CA"), ("third atom", ":6@CA"), ("fourth atom", ":12@CA"),
         ("solute mask", f":{PROTEIN}&!@H="), ("vector start", ":1@CA"), ("vector end", ":12@CA"),
         ("additional", "n"), ("another", "n")]

MENU = [
    ("_analyze_rmsd", []), ("_analyze_rmsd", [("select reference frame", "2")]),
    ("_analyze_rmsd", [("select region", "4")]), ("_analyze_rmsf", []), ("_analyze_contacts", []),
    ("_analyze_salt_bridges", []), ("_analyze_sasa", [("select region", "4")]), ("_analyze_dssp", []),
    ("_analyze_ramachandran", []), ("_analyze_ramachandran", [("dihedral type", "2")]),
    ("_analyze_ramachandran", [("dihedral type", "3")]), ("_analyze_bfactors", []), ("_analyze_distance", []),
    ("_analyze_angle", []), ("_analyze_dihedral", []), ("_analyze_vector", []), ("_analyze_pca", []),
    ("_analyze_clustering", [("select algorithm", "1")]), ("_analyze_clustering", [("select algorithm", "2")]),
    ("_analyze_autocorrelation", []), ("_analyze_pairwise_rmsd", []), ("_analyze_hbonds", []),
    ("_analyze_water_rdf", []), ("_analyze_water_shells", []),
    ("_analyze_water_shells", [("select solute", "2"), ("solute mask", f":{PROTEIN}")]),
    ("_analyze_density_maps", []), ("_analyze_radius_of_gyration", []),
    ("_analyze_contact_frequency_per_residue", []),
    ("_analyze_image_distance", []), ("_analyze_image_distance", [("select the solute", "2")]),
]


def _drive(monkeypatch, analyzer, method, answers, asked=None):
    """Run one menu analysis with scripted answers; returns what reached the console."""
    def prompt(processor, prompt, *args, **kwargs):
        if asked is not None:
            asked.append((processor, prompt, kwargs))
        for fragment, answer in list(answers) + MASKS:
            if fragment in str(prompt).lower():
                return answer
        assert kwargs.get("default") is not None, f"unscripted prompt without a default: {prompt!r}"
        return kwargs["default"]

    # the retry helpers reach prompt_with_context through their own module
    monkeypatch.setattr(mdm_module, "prompt_with_context", prompt)
    monkeypatch.setattr(prompts_module, "prompt_with_context", prompt)
    console = Console(width=110, record=True, file=io.StringIO(), force_terminal=False)
    m = object.__new__(MDM); m.processor = types.SimpleNamespace(console=console)
    getattr(m, method)(analyzer)
    return m, console.export_text()


MENU_IDS = [f"{m}-{'-'.join(a for _, a in ans) or 'defaults'}" for m, ans in MENU]


@pytest.mark.parametrize("method, answers", MENU, ids=MENU_IDS)
def test_menu_analysis_runs_clean_end_to_end(monkeypatch, analyzer, method, answers):
    """Analysis AND its display code; these swallow exceptions, so read the console."""
    before = analyzer.traj.xyz.copy()

    _, out = _drive(monkeypatch, analyzer, method, answers)

    failures = [line.strip() for line in out.splitlines() if "Error" in line or "failed" in line.lower()]
    assert not failures, "\n".join(failures)
    assert np.array_equal(before, analyzer.traj.xyz), "analysis moved the shared trajectory"


# ---------------------------------------------------------------- session recording

@pytest.mark.parametrize("method, answers", MENU, ids=MENU_IDS)
def test_every_prompt_carries_the_session_recording_contract(monkeypatch, analyzer, method, answers):
    """processor + module + description, and an options_map naming every choice."""
    asked = []
    m, _ = _drive(monkeypatch, analyzer, method, answers, asked)
    assert asked
    for processor, prompt, kwargs in asked:
        assert processor is m.processor, prompt
        assert kwargs.get("module", "").startswith("MD Manager - "), prompt
        assert kwargs.get("description"), prompt
        if kwargs.get("choices"):
            assert set(kwargs["choices"]) == set(kwargs.get("options_map") or {}), prompt


def _is_constant(node) -> bool:
    if isinstance(node, ast.Constant):
        return True
    if isinstance(node, (ast.List, ast.Tuple)):
        return all(_is_constant(e) for e in node.elts)
    if isinstance(node, ast.Dict):
        return all(_is_constant(k) for k in node.keys) and all(_is_constant(v) for v in node.values)
    return False


def test_analysis_prompt_text_choices_and_labels_never_depend_on_the_trajectory():
    """Replay uses a recorded answer only if the prompt TEXT matches exactly, and
    resolves a recorded choice by its options_map LABEL. Anything computed from the
    trajectory (a frame count, a residue range) in either would break replay on
    another system, so all three must be literals in the source."""
    source = Path(mdm_module.__file__).read_text()
    tree = ast.parse(source)
    manager = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "MolecularDynamicsManager")
    first = next(f.lineno for f in manager.body if getattr(f, "name", "") == "_prompt_frame_interval")
    prompting = {"prompt_with_context", "prompt_float_with_retry", "prompt_int_with_retry"}
    forwarding = {"_prompt_frame_index": 1, "_prompt_angle_box": 0}    # self.<name>(...): index of the text

    checked, offenders = 0, []
    for function in (f for f in manager.body if isinstance(f, ast.FunctionDef) and f.lineno >= first):
        # a name passed as choices/options_map is judged by what the function assigned to it
        assigned = {t.id: a.value for a in ast.walk(function) if isinstance(a, ast.Assign)
                    for t in a.targets if isinstance(t, ast.Name)}
        for call in (c for c in ast.walk(function) if isinstance(c, ast.Call)):
            name = getattr(call.func, "id", None) or getattr(call.func, "attr", None)
            if name in prompting:
                text = call.args[1] if len(call.args) > 1 else next((k.value for k in call.keywords if k.arg == "prompt"), None)
            elif name in forwarding and len(call.args) > forwarding[name]:
                text = call.args[forwarding[name]]
            else:
                continue
            if isinstance(text, ast.Name) and text.id in [a.arg for a in function.args.args]:
                continue                    # a helper passing its own `prompt` parameter through
            checked += 1
            parts = [("text", text)] + [(k.arg, k.value) for k in call.keywords if k.arg in ("choices", "options_map")]
            for what, node in parts:
                if isinstance(node, ast.Name):
                    node = assigned.get(node.id, node)
                if not _is_constant(node):
                    offenders.append(f"line {call.lineno}: {what} is computed")

    assert checked > 60
    assert not offenders, "\n".join(offenders)


def test_frame_interval_is_asked_when_the_file_stores_no_times(monkeypatch, analyzer):
    assert not analyzer.has_frame_times and analyzer.time_label == "Frame"
    assert analyzer.frame_times == [float(i) for i in range(10)]        # frame numbers, not an assumed time step

    answers = iter(["soon", "-1", "2.5"])                               # two rejected, then accepted
    monkeypatch.setattr(mdm_module, "prompt_with_context", lambda *a, **k: next(answers))
    console = Console(width=110, record=True, file=io.StringIO(), force_terminal=False)
    m = object.__new__(MDM); m.processor = types.SimpleNamespace(console=console)
    m._prompt_frame_interval(analyzer)

    out = console.export_text()
    assert "'soon' is not a positive number" in out and "'-1' is not a positive number" in out
    assert analyzer.has_frame_times and analyzer.time_label == "Time (ps)"
    assert analyzer.frame_times[:3] == [0.0, 2.5, 5.0]


def test_a_rejected_number_is_reported_not_silently_replaced(monkeypatch, analyzer):
    """Seen in 1.20.0: '4,5' typed as a cutoff ran the analysis at 4.5 without a word."""
    typed = iter(["4,5", "6.0"])

    def prompt(processor, prompt, *args, **kwargs):
        if "distance cutoff" in str(prompt).lower():
            return next(typed)
        return kwargs["default"]

    monkeypatch.setattr(mdm_module, "prompt_with_context", prompt)
    monkeypatch.setattr(prompts_module, "prompt_with_context", prompt)
    console = Console(width=110, record=True, file=io.StringIO(), force_terminal=False)
    m = object.__new__(MDM); m.processor = types.SimpleNamespace(console=console)
    m._analyze_salt_bridges(analyzer)

    out = console.export_text()
    assert "Please enter a valid number" in out
    assert "Distance cutoff: 6.0 Å" in out


def test_region_menu_offers_the_real_protein_range(monkeypatch, analyzer):
    monkeypatch.setattr(mdm_module, "prompt_with_context", lambda *args, **kwargs: "4")
    console = Console(width=110, record=True, file=io.StringIO(), force_terminal=False)
    m = object.__new__(MDM); m.processor = types.SimpleNamespace(console=console)
    region = m._get_analysis_region_selection(analyzer, "RMSD")
    assert region["mask"] == f":{PROTEIN}&!@H="       # not "80% of all residues", which is mostly water
    assert f"Protein residues {PROTEIN}" in console.export_text()


# ---------------------------------------------------------------- the real recorder

RECORDED_RUNS = [
    # region, algorithm=hierarchical, clusters, metric=distance RMSD, linkage=complete
    ("_analyze_clustering", ["2", "2", "3", "3", "3"],
     {"Select region": "C-alpha atoms (@CA)", "Select algorithm": "Hierarchical",
      "Select distance metric": "Distance RMSD", "Select linkage": "Complete linkage"}),
    # selection, spacing, grid limit, frame of reference=solute, fit mask, peaks
    ("_analyze_density_maps", [":WAT@O", "2.0", "40", "1", "@CA", "3"],
     {"Select frame of reference": "Solute frame (autoimage, then fit on a mask)"}),
    # region, reference frame, max lag (Enter = half the trajectory)
    ("_analyze_autocorrelation", ["1", "2", ""], {"Select region": "Backbone atoms (@C,CA,N,O&!:WAT)"}),
]


@pytest.mark.parametrize("method, typed, labels", RECORDED_RUNS, ids=[r[0] for r in RECORDED_RUNS])
def test_recorded_session_replays_without_a_single_live_prompt(monkeypatch, tmp_path, method, typed, labels):
    """Through ProPrep's own SessionManager: record with typed answers, replay with none."""
    import json
    import sys
    from proprep.utils.session_recorder import SessionManager

    def run(stdin_text, start):
        monkeypatch.setattr(sys, "stdin", io.StringIO(stdin_text))
        console = Console(width=110, record=True, file=io.StringIO(), force_terminal=False)
        m = object.__new__(MDM)
        m.processor = types.SimpleNamespace(console=console, session_manager=SessionManager())
        start(m.processor.session_manager)
        analyzer = TrajectoryAnalyzer(TOP, NC)
        try:
            getattr(m, method)(analyzer)       # an unanswerable prompt reads the empty stdin and raises EOFError
        finally:
            manager = m.processor.session_manager
            consumed = manager.replayer.interaction_index if manager.replayer else None   # stop() drops the replayer
            manager.stop()
        out = console.export_text()
        assert "Error" not in out, out
        return consumed, analyzer, out

    log = str(tmp_path / "session.json")
    _, recorded_analyzer, recorded_out = run("\n".join(typed) + "\n", lambda sm: sm.start_recording(log, {}))

    interactions = json.load(open(log))["interactions"]
    assert [i["response"] for i in interactions] == typed
    for interaction in interactions:                          # what the session editor shows
        assert interaction["context"]["module"].startswith("MD Manager - ")
        assert interaction["context"]["description"]
    recorded_labels = {i["prompt"]: i["context"].get("option_label") for i in interactions}
    for prompt, label in labels.items():
        assert recorded_labels[prompt] == label

    consumed, replayed_analyzer, replayed_out = run("", lambda sm: sm.start_replay(log))

    assert consumed == len(interactions)      # every answer came from the log
    assert replayed_out == recorded_out
    for key, value in recorded_analyzer.data.items():
        assert replayed_analyzer.data[key] == value, key
