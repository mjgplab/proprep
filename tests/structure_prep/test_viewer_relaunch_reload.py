"""A viewer tab left open in the browser follows ProPrep's relaunches.

Showing another structure list (a docked pose with the receptor) or re-reading
edited files restarts the viewer server. In the terminal the browser shows the
plain viewer page, and it never changed: the new server restarted its config
version at 1, which an untouched page already held, and the page rebuilt only
representations, never structures. Each server start now carries an instance
id (in /config and /version) and the page reloads when it changes, keeping the
camera when the files are the same (checked in headless Chrome with NGL). A
relaunch also reclaims the previous server's port: asking for 8765 again,
with 8765 held elsewhere, moved the server and stranded the tab.
"""

import json
import socket
import urllib.request
from unittest.mock import MagicMock

from proprep.structure_prep import interactive_structure_viewer as isv
from proprep.structure_prep.viewer_server import ViewerServer


def _free_port():
    with socket.socket() as s:
        s.bind(("localhost", 0))
        return s.getsockname()[1]


def _get(port, path):
    with urllib.request.urlopen(f"http://localhost:{port}{path}", timeout=5) as r:
        return json.loads(r.read())


def test_every_server_start_has_its_own_instance_id(tmp_path):
    ids = []
    for _ in range(2):
        server = ViewerServer(config={"structures": []}, structure_files=[], port=_free_port())
        assert server.start(open_browser=False)
        try:
            version, config = _get(server.port, "/version"), _get(server.port, "/config")
            assert version["instance"] and version["instance"] == config["_instance_id"]
            ids.append(version["instance"])
        finally:
            server.stop()
    assert ids[0] != ids[1]


def test_the_page_reloads_on_a_new_instance_and_keeps_the_camera():
    with open(isv.__file__.replace("interactive_structure_viewer.py", "templates/ngl_viewer.html")) as handle:
        page = handle.read()
    assert "vinfo.instance !== loadedInstance" in page and "location.reload()" in page
    assert "saveOrientationForReload()" in page and "restoreOrientationAfterReload();" in page


def test_a_relaunch_reclaims_the_previous_port(monkeypatch):
    previous = MagicMock(port=8799)
    previous.is_running.return_value = True
    monkeypatch.setattr(isv, "_active_viewer_server", previous)
    ports = []

    class FakeServer:
        def __init__(self, **kwargs):
            ports.append(kwargs["port"])
            self.port = kwargs["port"]
        def start(self, open_browser=True):
            return False                     # stop before printing or binding anything

    monkeypatch.setattr("proprep.structure_prep.viewer_server.ViewerServer", FakeServer)
    v = isv.InteractiveStructureViewer.__new__(isv.InteractiveStructureViewer)
    v.selected_structures, v.available_annotations, v.annotation_config = ["a.pdb"], {}, {}
    v.viewer_config, v.shape_config, v.trajectory_files, v.density_by_file = {}, {}, {}, {}
    v.scene_override, v.processor = None, None
    v._scene_dir = lambda: "."
    monkeypatch.delenv("PROPREP_BATCH", raising=False)
    v._launch_viewer(open_browser=False)
    assert ports == [8799]
    previous.stop.assert_called_once()
