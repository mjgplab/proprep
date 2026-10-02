"""Structure viewer: the coordinate axes toggle.

mgp oriented a structure onto its principal axes (Structure Orientation
module) and saw nothing in the viewer to show it: NGL frames every loaded
structure the same way whatever its coordinates, so the rotated structure
looked like the original. The module did drop axis markers, but as four
unlabelled spheres 15 A out along each axis, which sit inside any protein of
ordinary size. The viewer now draws the frame itself, from an "Axes" button
under View Controls: labelled arrows from the origin, each reaching just past
the far edge of the structure along its axis, in colours chosen by computed
contrast for each background. The orientation module asks for them through
the viewer config, and a scene keeps the setting.

Looked at in headless Chrome (NGL 2.5.0) against crambin, raw and after
centring + principal-axis rotation: arrows and letters drawn on the dark,
white and black backgrounds, the component gone after toggling off, the
config request applied once at load and not re-applied by a later refresh,
no page errors.
"""

from __future__ import annotations

import json
import re
from pathlib import Path


from proprep.structure_prep import viewer_coordinator as vc_mod
from proprep.structure_prep.interactive_structure_viewer import InteractiveStructureViewer
from proprep.structure_prep.structure_orientation import StructureOrientationModule

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = (ROOT / "src/proprep/structure_prep/templates/ngl_viewer.html").read_text()


def _luminance(hex_colour):
    rgb = [int(hex_colour.lstrip("#")[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    lin = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def _contrast(a, b):
    hi, lo = sorted((_luminance(a), _luminance(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


def _viewer(tmp_path, structures=()):
    v = InteractiveStructureViewer.__new__(InteractiveStructureViewer)
    v.selected_structures = list(structures)
    v.available_annotations = {}
    v.annotation_config = {}
    v.viewer_config = {}
    v.shape_config = {}
    v.trajectory_files = {}
    v.density_by_file = {}
    v.scene_override = None
    v.processor = None
    v._scene_dir = lambda: str(tmp_path)
    return v


# ---------------------------------------------------------------- the page

def test_axes_button_toggles_a_labelled_arrow_component():
    assert 'id="axes-btn"' in TEMPLATE and 'onclick="toggleAxes()"' in TEMPLATE
    assert ">Axes: Off</button>" in TEMPLATE                      # off until asked for
    assert "shape.addArrow([0, 0, 0], tip, rgb, shaft, name + ' axis')" in TEMPLATE
    assert "new NGL.TextBuffer({" in TEMPLATE and "showBorder: true" in TEMPLATE
    assert "stage.removeComponent(axesComponent)" in TEMPLATE      # off takes it away
    assert "Math.abs(box.min[k]), Math.abs(box.max[k])" in TEMPLATE  # reach past the structure


def test_axes_colours_contrast_with_every_background_and_with_each_other():
    """The user is visually impaired: hue alone is not enough, luminance must
    separate the axes too, and each set must clear its backgrounds."""
    light = re.search(r"light: \{ X: '(#[0-9a-f]{6})', Y: '(#[0-9a-f]{6})', Z: '(#[0-9a-f]{6})'", TEMPLATE).groups()
    dark = re.search(r"dark:\s+\{ X: '(#[0-9a-f]{6})', Y: '(#[0-9a-f]{6})', Z: '(#[0-9a-f]{6})'", TEMPLATE).groups()
    for colours, backgrounds in ((light, ("#1e1e1e", "#000000")), (dark, ("#ffffff",))):
        for c in colours:
            for bg in backgrounds:
                assert _contrast(c, bg) >= 5.9, (c, bg, _contrast(c, bg))
        lums = sorted(_luminance(c) for c in colours)
        for lo, hi in zip(lums, lums[1:]):
            assert (hi + 0.05) / (lo + 0.05) >= 1.3, colours
    assert "densityLuminance(parseInt(BACKGROUNDS[backgroundIndex].color.slice(1), 16)) < 0.5" in TEMPLATE
    assert TEMPLATE.count("rebuildAxes();") >= 2                   # recoloured when the background changes


def test_axes_setting_is_saved_with_a_scene_and_asked_for_through_the_config():
    assert "axes: axesOn," in TEMPLATE                              # saved
    assert "applyAxes(cfg.axes)" in TEMPLATE                        # restored
    assert TEMPLATE.count("applyAxesFromConfig(config);") == 2      # first load + live refresh
    assert "if (wanted === axesFromConfig)" in TEMPLATE             # a repeat does not undo the button


def test_axes_setting_round_trips_through_a_scene_file(tmp_path, monkeypatch):
    pdb = tmp_path / "prot.pdb"; pdb.write_text("END\n")
    payload = {"name": "framed", "representations": {"0": []}, "camera": None,
               "camera_type": "orthographic", "background": "#ffffff", "depth_cue": True, "axes": True}
    path = _viewer(tmp_path, [str(pdb)])._save_scene_payload(payload)["path"]
    assert json.loads(Path(path).read_text())["axes"] is True

    fresh = _viewer(tmp_path)
    monkeypatch.setattr(fresh, "_launch_viewer", lambda open_browser=True: True)
    assert fresh.load_scene(path) is True
    assert fresh._build_viewer_config()["axes"] is True


def test_viewer_config_request_reaches_the_page_and_a_scene_wins(tmp_path):
    pdb = tmp_path / "prot.pdb"; pdb.write_text("END\n")
    v = _viewer(tmp_path, [str(pdb)])
    assert "axes" not in v._build_viewer_config()                   # nothing asked: nothing sent
    v.viewer_config = {"axes": True}
    assert v._build_viewer_config()["axes"] is True
    v.scene_override = {"_for": [str(pdb)], "representations": {}, "axes": False, "scene_id": "s@1"}
    assert v._build_viewer_config()["axes"] is False


# ---------------------------------------------------------------- coordinator + orientation module

class _Viewer:
    def __init__(self, structures):
        self.selected_structures = structures
        self.viewer_config = {}
        self.pushed = 0
        self.launched = 0

    def update_annotations(self, *a, **k):
        self.pushed += 1
        return True

    def _launch_viewer(self, open_browser=True):
        self.launched += 1


def _coordinator(monkeypatch, structures, running):
    c = vc_mod.ViewerCoordinator()
    fake = _Viewer(structures)
    monkeypatch.setattr(c, "_ensure_viewer", lambda: fake)
    monkeypatch.setattr(c, "is_running", lambda: running)
    monkeypatch.setattr(c, "_owns_live_server", lambda v: running)   # the running server is this viewer's
    monkeypatch.setattr(vc_mod, "_is_web_shell_mode", lambda: False)
    return c, fake


def test_show_axes_sets_the_request_and_pushes_it_to_a_running_viewer(monkeypatch):
    c, fake = _coordinator(monkeypatch, ["/x/prot.pdb"], running=True)
    c.show_axes(True)
    assert fake.viewer_config == {"axes": True} and fake.pushed == 1 and fake.launched == 0
    c.show_axes(False)
    assert fake.viewer_config == {"axes": False} and fake.pushed == 2


def test_show_axes_is_silent_in_cli_mode_without_a_viewer_unless_forced(monkeypatch):
    c, fake = _coordinator(monkeypatch, ["/x/prot.pdb"], running=False)
    c.show_axes(True)
    assert fake.viewer_config == {"axes": True} and fake.launched == 0
    c.show_axes(True, force=True)
    assert fake.launched == 1


def test_show_axes_needs_a_structure(monkeypatch):
    c, fake = _coordinator(monkeypatch, [], running=True)
    c.show_axes(True)
    assert fake.viewer_config == {} and fake.pushed == 0


def test_orientation_module_asks_for_the_axes_instead_of_spheres(monkeypatch):
    calls = []

    class _Coord:
        def show_structure(self, path, *, force=False):
            calls.append(("structure", path, force))

        def show_axes(self, on=True, *, force=False):
            calls.append(("axes", on, force))

        def show_sphere(self, *a, **k):
            raise AssertionError("the sphere markers are gone")

    monkeypatch.setattr(vc_mod, "viewer", _Coord())
    mod = StructureOrientationModule.__new__(StructureOrientationModule)
    mod._show_orientation_view("/x/prot.pdb")
    assert calls == [("structure", "/x/prot.pdb", False), ("axes", True, False)]
    assert not hasattr(StructureOrientationModule, "_AXIS_MARKERS")
