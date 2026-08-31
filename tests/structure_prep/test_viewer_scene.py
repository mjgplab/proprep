"""Scene save/load for the Structure Viewer: the page's representation and camera
state round-trips through a JSON file in the project directory and back into
the served config."""
import io
import json
import os

import pytest
from rich.console import Console

from proprep.structure_prep.interactive_structure_viewer import InteractiveStructureViewer
from proprep.structure_prep.viewer_server import ViewerHTTPRequestHandler, ViewerServer


class _WS(dict):
    def set(self, k, v): self[k] = v


class _Proc:
    def __init__(self, ws):
        self.console = Console(file=io.StringIO(), force_terminal=False)
        self._ws = ws
    def _get_workspace(self): return self._ws


def _viewer(tmp_path, structures):
    v = InteractiveStructureViewer.__new__(InteractiveStructureViewer)
    v.selected_structures = structures
    v.available_annotations = {}
    v.annotation_config = {}
    v.viewer_config = {}
    v.shape_config = {}
    v.scene_override = None
    v._last_saved_scene = None
    v.processor = _Proc(_WS(output_dir=str(tmp_path)))
    return v


PAYLOAD = {
    "name": "my figure",
    "representations": {"0": [
        {"id": "default_protein", "label": "Protein", "selection": "protein", "style": "cartoon",
         "color": "#0078d4", "visible": True, "opacity": 0.4},
        {"id": "user_1", "label": "Heme", "selection": "[HEC]", "style": "licorice", "color": "element", "visible": True},
    ]},
    "camera": {"orientation": [float(i) for i in range(16)]},
    "camera_type": "perspective",
    "background": "#ffffff",
}


def test_save_writes_relative_paths_and_records_workspace(tmp_path):
    pdb = tmp_path / "sub" / "prot.pdb"; pdb.parent.mkdir(); pdb.write_text("END\n")
    v = _viewer(tmp_path, [str(pdb)])
    result = v._save_scene_payload(PAYLOAD)
    assert result["ok"]
    path = result["path"]
    assert path == str(tmp_path / "my_figure.scene.json")        # sanitized name, project dir
    scene = json.loads(open(path).read())
    assert scene["format"] == "proprep-scene"
    assert scene["structures"] == [{"index": 0, "name": "prot", "path": os.path.join("sub", "prot.pdb")}]
    assert scene["representations"]["0"][1]["id"] == "user_1"     # hand-added rep persisted
    assert scene["camera"]["orientation"][5] == 5.0 and scene["background"] == "#ffffff"
    assert v.processor._ws["viewer_scenes"] == [path]


def test_load_resolves_paths_and_config_carries_scene(tmp_path, monkeypatch):
    pdb = tmp_path / "prot.pdb"; pdb.write_text("END\n")
    v = _viewer(tmp_path, [str(pdb)])
    path = v._save_scene_payload(PAYLOAD)["path"]

    fresh = _viewer(tmp_path, [])
    launched = []
    monkeypatch.setattr(fresh, "_launch_viewer", lambda open_browser=True: launched.append(open_browser) or True)
    assert fresh.load_scene(path) is True
    assert fresh.selected_structures == [str(pdb)]
    assert launched == [True]

    cfg = fresh._build_viewer_config()
    reps = cfg["structures"][0]["representations"]
    assert [r["id"] for r in reps] == ["default_protein", "user_1"]   # scene reps, not defaults
    assert reps[0]["opacity"] == 0.4
    assert cfg["camera"]["orientation"] == [float(i) for i in range(16)]
    assert cfg["camera_type"] == "perspective" and cfg["background"] == "#ffffff"
    assert cfg["scene_id"].startswith("my_figure.scene.json@")


def test_load_refuses_missing_structures(tmp_path):
    pdb = tmp_path / "prot.pdb"; pdb.write_text("END\n")
    v = _viewer(tmp_path, [str(pdb)])
    path = v._save_scene_payload(PAYLOAD)["path"]
    pdb.unlink()
    fresh = _viewer(tmp_path, [])
    assert fresh.load_scene(path) is False
    assert fresh.scene_override is None


def test_scene_override_is_dropped_when_structures_change(tmp_path):
    pdb = tmp_path / "prot.pdb"; pdb.write_text("END\n")
    v = _viewer(tmp_path, [str(pdb)])
    v.scene_override = {"_for": [str(pdb)], "representations": {0: [{"id": "x", "label": "x", "selection": "all",
                        "style": "line", "color": "red", "visible": True}]}, "camera": None}
    assert v._build_viewer_config()["structures"][0]["representations"][0]["id"] == "x"
    v.selected_structures = [str(tmp_path / "other.pdb")]
    ids = [r["id"] for r in v._build_viewer_config()["structures"][0]["representations"]]
    assert "x" not in ids and "default_protein" in ids


def test_server_scene_request_rides_on_version_and_post_hits_sink(tmp_path):
    seen = []
    server = ViewerServer(config={"structures": []}, structure_files=[], port=8799,
                          scene_sink=lambda payload: seen.append(payload) or {"ok": True, "path": "x"})
    token = server.request_scene("fig1")
    assert ViewerHTTPRequestHandler.scene_request == {"token": token, "name": "fig1"}
    server.clear_scene_request()
    assert ViewerHTTPRequestHandler.scene_request is None
    # drive do_POST without a socket
    handler = ViewerHTTPRequestHandler.__new__(ViewerHTTPRequestHandler)
    body = json.dumps({"name": "fig1", "representations": {}}).encode()
    handler.path = "/scene"; handler.headers = {"Content-Length": str(len(body))}
    handler.rfile = io.BytesIO(body); handler.wfile = io.BytesIO()
    sent = {}
    handler.send_response = lambda code: sent.setdefault("code", code)
    handler.send_header = lambda k, v: None
    handler.end_headers = lambda: None
    handler.do_POST()
    assert sent["code"] == 200 and seen[0]["name"] == "fig1"
    assert json.loads(handler.wfile.getvalue()) == {"ok": True, "path": "x"}
