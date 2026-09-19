"""Structure viewer: trajectory entry point, superpose fix, depth cue, movies.

Seen in 1.20.0:
- The only way to a playing trajectory was MD Manager > Analyze > pick a
  simulation > "Select analysis type" > 4, a prompt that defaults to 3. mgp,
  who asked for the feature, could not find it. It belongs in the Structure
  Viewer: pick a topology and a trajectory from disk.
- The Trajectory panel's "superpose" box scrambled the molecule. NGL's
  Superposition.transform() copies a Float32Array (what Trajectory hands it for
  every frame) into a 4-wide matrix using ONE index for the 3-wide source and
  the 4-wide destination. With NGL 2.5.0, fitting a rigid tumbled 8-atom chain
  onto itself gave RMSD 15.6 A and a 25.6 A "bond"; patched, RMSD 0.000 A.
- Save Screenshot wrote the default Dark background (30,30,30) as (3,3,3):
  NGL's makeImage emits an opaque background in linear light.

Ported from MDViz (Rocky, /opt/lab/viz_server): movie export and the depth-cue
toggle. Checked in headless Chrome against a 300-frame chignolin trajectory:
bonds stay 1.59 A under superpose, fog goes 50/100 <-> 100/100 and is saved
with the scene, and a 23-frame H.264 MP4 (ffprobe: High, 30 fps) downloads with
a (30,30,30) background. That run also found that NGL 2.5.0's
loadFrame(array, cb) never calls cb, which hung every interpolated movie frame.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from rich.console import Console

from proprep.md_prep import trajectory_view
from proprep.structure_prep.interactive_structure_viewer import InteractiveStructureViewer

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = (ROOT / "src/proprep/structure_prep/templates/ngl_viewer.html").read_text()


class _WS(dict):
    def set(self, k, v): self[k] = v


class _Proc:
    def __init__(self, ws):
        self._ws = ws
        self.console = Console(width=200, record=True)
        self.workspace = ws

    def _get_workspace(self): return self._ws


def _viewer(tmp_path, structures=()):
    v = InteractiveStructureViewer.__new__(InteractiveStructureViewer)
    v.selected_structures = list(structures)
    v.available_annotations = {}
    v.annotation_config = {}
    v.viewer_config = {}
    v.shape_config = {}
    v.scene_override = None
    v._last_saved_scene = None
    v.processor = _Proc(_WS(output_dir=str(tmp_path)))
    return v


# ---------------------------------------------------------------- entry point

def _no_structures(v, monkeypatch):
    monkeypatch.setattr(v, "_has_structures", lambda workspace: False)


def test_trajectory_option_is_gated_on_what_the_structure_loader_loaded(tmp_path, monkeypatch):
    from proprep.utils.enhanced_menu import OptionStatus

    v = _viewer(tmp_path); _no_structures(v, monkeypatch)
    ws = v.processor.workspace
    options = v.get_menu_options()
    assert list(options)[:2] == ["launch", "trajectory"]

    def status(option):
        enhanced = v.get_enhanced_menu_options(ws)
        assert [o.key for o in enhanced] == [str(i) for i in range(1, len(options) + 1)]
        return dict(zip(options, enhanced))[option]

    # nothing loaded: blocked, and the note says where to load it
    assert status("trajectory").status == OptionStatus.BLOCKED
    assert "Structure Loader" in status("trajectory").dependency_text
    assert not v.can_process(ws) and "topology + trajectory" in v.availability_note(ws)

    ws.set("parm7_file", "sys.prmtop")                       # a topology alone (the rst7 route) is not a trajectory
    assert status("trajectory").status == OptionStatus.BLOCKED

    ws.set("trajectory_files", ["prod1.nc", "prod2.nc"])     # what Structure Loader > Load AMBER files stores
    assert status("trajectory").status == OptionStatus.AVAILABLE
    assert status("launch").status == OptionStatus.BLOCKED   # still no PDB structure to launch
    assert status("load_scene").status == OptionStatus.AVAILABLE
    assert v.can_process(ws) and v.availability_note(ws) is None


def test_an_option_added_by_other_work_gets_the_structure_gate_not_a_keyerror(tmp_path, monkeypatch):
    """feat/viewer-electron-density adds a "density" option to this same menu; merged
    with a per-option table that did not name it, the menu died with KeyError."""
    from proprep.utils.enhanced_menu import OptionStatus

    v = _viewer(tmp_path); _no_structures(v, monkeypatch)
    options = dict(v.get_menu_options(), density="Launch viewer with electron density")
    monkeypatch.setattr(v, "get_menu_options", lambda: options)

    enhanced = dict(zip(options, v.get_enhanced_menu_options(v.processor.workspace)))
    assert enhanced["density"].status == OptionStatus.BLOCKED and enhanced["density"].key == str(len(options))
    assert "Load a structure first" in enhanced["density"].dependency_text


def test_trajectory_workflow_plays_the_loaded_files_without_asking_for_any(tmp_path, monkeypatch):
    prmtop = tmp_path / "run" / "sys.prmtop"; prmtop.parent.mkdir(); prmtop.write_text("x")
    segments = [tmp_path / "run" / "prod1.nc", tmp_path / "run" / "prod2.dcd"]
    for f in segments: f.write_text("x")
    shown = []
    import proprep.utils.file_browser as fb
    monkeypatch.setattr(fb, "file_browser", lambda *a, **k: pytest.fail("the Structure Loader already chose the files"))
    monkeypatch.setattr(trajectory_view, "prepare_and_show",
                        lambda processor, console, top, traj, out_dir, module: shown.append((top, traj, out_dir, module)) or True)

    v = _viewer(tmp_path)
    v.processor.workspace.set("parm7_file", str(prmtop))
    v.processor.workspace.set("trajectory_files", [str(f) for f in segments])
    assert v._view_trajectory_workflow() is True

    # every segment, in order, written to the project directory
    assert shown == [(str(prmtop), [str(f) for f in segments], str(tmp_path / "viewer"), "Structure Viewer - Trajectory")]
    assert "Segment 2:" in v.console.export_text()


def test_trajectory_workflow_says_where_to_load_when_nothing_is_loaded(tmp_path, monkeypatch):
    monkeypatch.setattr(trajectory_view, "prepare_and_show", lambda *a, **k: pytest.fail("nothing is loaded"))
    v = _viewer(tmp_path)
    assert v._view_trajectory_workflow() is False
    assert "Structure Loader > Load AMBER topology & coordinate files" in v.console.export_text()

    v.processor.workspace.set("parm7_file", str(tmp_path / "gone.prmtop"))
    v.processor.workspace.set("trajectory_files", [str(tmp_path / "gone.nc")])
    assert v._view_trajectory_workflow() is False
    assert "no longer exist" in v.console.export_text()


def test_cpptraj_reads_every_loaded_segment_in_order(tmp_path):
    script = trajectory_view.cpptraj_input("p.prmtop", ["a.nc", "b.dcd", "c.xtc"], "v.pdb", "v.nc", strip_solvent=False)
    assert script.splitlines()[:5] == ["parm p.prmtop", "trajin a.nc", "trajin b.dcd", "trajin c.xtc", "autoimage"]
    assert trajectory_view.cpptraj_input("p", "one.nc", "a", "b", strip_solvent=False).splitlines()[1] == "trajin one.nc"


def test_prepare_and_show_asks_about_solvent_then_opens_the_viewer(tmp_path, monkeypatch):
    import proprep.utils.prompts as prompts
    import proprep.structure_prep.viewer_coordinator as coordinator
    asked, opened = [], []
    monkeypatch.setattr(prompts, "confirm_with_context",
                        lambda processor, prompt, **kw: asked.append((prompt, kw["module"], kw["description"])) or False)
    monkeypatch.setattr(trajectory_view, "write_view_files",
                        lambda top, traj, out, strip_solvent: (Path(out) / "v.pdb", Path(out) / "v.nc"))
    monkeypatch.setattr(trajectory_view, "frame_count", lambda nc: 300)
    monkeypatch.setattr(coordinator.viewer, "show_trajectory",
                        lambda pdb, nc, show_waters, force: opened.append((pdb, nc, show_waters, force)))
    console = Console(width=120, record=True)

    assert trajectory_view.prepare_and_show(object(), console, "s.prmtop", "p.nc", str(tmp_path), module="X - Trajectory")

    assert asked == [("Strip water and ions from the viewed trajectory? (no = keep them visible)",
                      "X - Trajectory", "Strip solvent for viewing")]
    assert opened == [(str(tmp_path / "v.pdb"), str(tmp_path / "v.nc"), True, True)]   # kept solvent -> waters shown
    assert "300 frames" in console.export_text()


def test_md_manager_topology_prompt_text_does_not_carry_the_file_count():
    """Replay matches prompt text exactly, so a count in it breaks replay on another directory."""
    source = (ROOT / "src/proprep/md_prep/molecular_dynamics_manager.py").read_text()
    assert re.findall(r'f"Select topology file[^"]*\{[^"]*"', source) == []
    for text in ('"Select topology file",', '"Select topology file, or \'n\' to browse elsewhere"',
                 '"Select topology file, \'browse\' to pick manually, or \'cancel\'"', '"Select topology file, or \'cancel\'"'):
        assert text in source, text


# ---------------------------------------------------------------- superpose fix

def _script_block(marker: str) -> str:
    return next(b for b in re.findall(r"<script>(.*?)</script>", TEMPLATE, flags=re.S) if marker in b)


def test_superposition_patch_is_installed_before_the_viewer_code():
    assert TEMPLATE.index("patchSuperposition") < TEMPLATE.index("new NGL.Stage(")
    assert "typeof NGL === 'undefined'" in _script_block("patchSuperposition")     # CDN down: no exception


@pytest.mark.skipif(shutil.which("node") is None, reason="needs node to run the template's JavaScript")
def test_superposition_patch_applies_the_fit_matrix_atom_by_atom(tmp_path):
    """The template's own patch block, run against a stand-in NGL whose stock
    transform has NGL's bug, on a rotation + translation with a known answer."""
    rng = np.random.default_rng(7)
    xyz = rng.normal(scale=6.0, size=(9, 3))
    a, b = np.radians(70), np.radians(25)
    rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    rx = np.array([[1, 0, 0], [0, np.cos(b), -np.sin(b)], [0, np.sin(b), np.cos(b)]])
    rotation, shift = rx @ rz, np.array([12.0, -7.0, 3.0])
    matrix = np.eye(4); matrix[:3, :3] = rotation; matrix[:3, 3] = shift
    expected = xyz @ rotation.T + shift

    script = tmp_path / "patch.js"
    script.write_text(
        "const calls = [];\n"
        "globalThis.NGL = { Superposition: function () {} };\n"
        "NGL.Superposition.prototype.transform = function (c) { calls.push('stock'); };\n"
        + _script_block("patchSuperposition") +
        f"\nconst sp = new NGL.Superposition(); sp.transformationMatrix = {{ elements: {json.dumps(matrix.T.ravel().tolist())} }};\n"
        f"const coords = new Float32Array({json.dumps(xyz.ravel().tolist())});\n"
        "sp.transform(coords); sp.transform({ atomCount: 3 });\n"
        "console.log(JSON.stringify({ coords: Array.from(coords), calls }));\n")
    out = json.loads(subprocess.run(["node", str(script)], capture_output=True, text=True, check=True).stdout)

    assert np.allclose(np.array(out["coords"]).reshape(-1, 3), expected, atol=1e-4)
    assert out["calls"] == ["stock"]            # a Structure argument still takes NGL's own path


# ---------------------------------------------------------------- depth cue

def test_depth_cue_toggles_ngl_fog_and_has_a_button():
    assert 'id="fog-btn"' in TEMPLATE and 'onclick="toggleDepthCue()"' in TEMPLATE
    assert "on: { fogNear: 50, fogFar: 100 }, off: { fogNear: 100, fogFar: 100 }" in TEMPLATE
    assert "depth_cue: depthCue," in TEMPLATE                                       # saved with a scene
    assert "applyDepthCue(cfg.depth_cue)" in TEMPLATE                               # and restored from one


def test_depth_cue_round_trips_through_a_scene_file(tmp_path, monkeypatch):
    pdb = tmp_path / "prot.pdb"; pdb.write_text("END\n")
    payload = {"name": "flat", "representations": {"0": []}, "camera": None,
               "camera_type": "orthographic", "background": "#ffffff", "depth_cue": False}
    path = _viewer(tmp_path, [str(pdb)])._save_scene_payload(payload)["path"]
    assert json.loads(Path(path).read_text())["depth_cue"] is False

    fresh = _viewer(tmp_path)
    monkeypatch.setattr(fresh, "_launch_viewer", lambda open_browser=True: True)
    assert fresh.load_scene(path) is True
    assert fresh._build_viewer_config()["depth_cue"] is False


# ---------------------------------------------------------------- movie + screenshot

def test_movie_panel_appears_with_the_trajectory_and_records_with_m():
    for element in ("movie-section", "mv-start", "mv-end", "mv-stride", "mv-interp", "mv-fps", "mv-scale",
                    "mv-quality", "mv-aa", "mv-record", "movie-progress", "movie-status", "rec-shield", "rec-badge"):
        assert f'id="{element}"' in TEMPLATE, element
    assert "document.getElementById('movie-section').classList.remove('hidden');" in TEMPLATE
    assert "if (e.key === 'm') { e.preventDefault(); toggleMovie(); }" in TEMPLATE
    # while recording every key but Esc is inert
    keys = TEMPLATE[TEMPLATE.index("document.addEventListener('keydown'"):]
    assert keys.index("if (movie) {") < keys.index("if (isTypingInField() || e.metaKey")


MP4_MUXER = ROOT / "src/proprep/structure_prep/templates/mp4-muxer-5.2.2.min.js"


def test_mp4_writer_ships_with_proprep_and_is_loaded_only_when_recording():
    """Movie export works offline: the writer is vendored, not fetched from a CDN."""
    import hashlib
    assert hashlib.sha256(MP4_MUXER.read_bytes()).hexdigest() == \
        "23696c3869c78ae0743a9231e6940b357aaddf563a6bc7cc5b3c58303c0b5ef5"     # npm mp4-muxer@5.2.2, as MDViz runs it
    licence = MP4_MUXER.with_name("mp4-muxer-5.2.2.LICENSE.txt").read_text()
    assert licence.startswith("MIT License") and "Vanilagy" in licence            # MIT: the notice travels with the file

    assert "const MP4_MUXER_URL = '/vendor/mp4-muxer.min.js';" in TEMPLATE
    assert "mp4-muxer" not in "".join(re.findall(r"<script src=[^>]+>", TEMPLATE))     # the viewer loads without it
    assert "cdn.jsdelivr.net/npm/mp4-muxer" not in TEMPLATE
    record = TEMPLATE[TEMPLATE.index("async function recordMovie()"):]
    assert "await loadMp4Muxer();" in record[:record.index("function toggleMovie()")]


def test_viewer_server_serves_the_vendored_mp4_writer(tmp_path):
    import urllib.error
    import urllib.request

    import http.server
    import threading
    from proprep.structure_prep import viewer_server as vs
    handler = vs.ViewerHTTPRequestHandler
    handler.template_path = str(MP4_MUXER.with_name("ngl_viewer.html"))
    httpd = http.server.HTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    try:
        base = f"http://127.0.0.1:{httpd.server_address[1]}"
        with urllib.request.urlopen(base + "/vendor/mp4-muxer.min.js?v=1") as r:
            assert r.headers["Content-Type"] == "application/javascript"
            assert r.read() == MP4_MUXER.read_bytes()
        with pytest.raises(urllib.error.HTTPError):                     # one fixed file, not a directory
            urllib.request.urlopen(base + "/vendor/../ngl_viewer.html")
    finally:
        httpd.shutdown(); httpd.server_close()


def test_interpolated_movie_frames_wait_on_the_cache_not_on_ngls_callback():
    """NGL 2.5.0: loadFrame(array, cb) never calls cb, and setFrameInterpolated
    passes an array for any frame that is not cached yet."""
    record = TEMPLATE[TEMPLATE.index("async function recordMovie()"):]
    ensure = record.index("await ensureFrameCached(traj, from); await ensureFrameCached(traj, to);")
    assert ensure < record.index("traj.setFrameInterpolated(")
    assert "if (traj.frameCache[i]) return Promise.resolve();" in TEMPLATE


def test_images_are_rendered_transparent_and_painted_on_the_real_background():
    assert "transparent: false" not in TEMPLATE
    helper = TEMPLATE[TEMPLATE.index("async function renderOnBackground("):TEMPLATE.index("function captureScreenshot()")]
    assert "transparent: true" in helper and "BACKGROUNDS[backgroundIndex].color" in helper
    screenshot = TEMPLATE[TEMPLATE.index("function captureScreenshot()"):]
    assert "renderOnBackground(" in screenshot[:400]
    capture = TEMPLATE[TEMPLATE.index("const capture = async () => {"):]
    assert "transparent: true" in capture[:500] and "ctx.fillStyle = BACKGROUNDS[backgroundIndex].color;" in capture[:1400]


@pytest.mark.skipif(shutil.which("node") is None, reason="needs node to syntax-check the template's JavaScript")
def test_template_javascript_parses(tmp_path):
    for i, block in enumerate(re.findall(r"<script>(.*?)</script>", TEMPLATE, flags=re.S)):
        f = tmp_path / f"block{i}.js"; f.write_text(block)
        result = subprocess.run(["node", "--check", str(f)], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr


def test_h_hides_the_control_panel_even_without_a_trajectory():
    """As in MDViz. The key handler used to return early when no trajectory was
    loaded, which would have made h work only in trajectory views."""
    keys = TEMPLATE[TEMPLATE.index("document.addEventListener('keydown'"):]
    keys = keys[:keys.index("});")]
    h = keys.index("if (e.key === 'h' || e.key === 'H') { e.preventDefault(); toggleSidebar(); return; }")
    assert keys.index("if (movie) {") < h                                     # inert while a movie records
    assert keys.index("if (isTypingInField() || e.metaKey || e.ctrlKey || e.altKey) return;") < h   # an h typed in a box is a letter
    assert h < keys.index("if (!trajState) return;")                          # but it needs no trajectory
    assert "fills the window (h)" in TEMPLATE                                 # the button says so
