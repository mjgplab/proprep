"""Viewer picks for the command line: click an atom or a bond, and ProPrep gets it.

Picks never left the page: the measurement and density tools use clicks
inside the browser only. A command-line tool (Molecular Docking: fix a bond
order, set a charge, choose a rotatable bond) can now ask the open viewer
for one click. The request rides on /version like a scene request, the page
shows it in a high-contrast banner with Cancel (and Escape), and the click
comes back as POST /pick. Picking is never the only way: every caller also
takes typed input, as the user needs (no pointer-precision-only controls).
"""

import io
import json
import threading
from pathlib import Path

import pytest

from proprep.structure_prep import viewer_coordinator as vc_mod
from proprep.structure_prep.viewer_server import ViewerHTTPRequestHandler, ViewerServer

ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = (ROOT / "src/proprep/structure_prep/templates/ngl_viewer.html").read_text()


def _post(path, payload):
    handler = ViewerHTTPRequestHandler.__new__(ViewerHTTPRequestHandler)
    body = json.dumps(payload).encode()
    handler.path = path
    handler.headers = {"Content-Length": str(len(body))}
    handler.rfile = io.BytesIO(body)
    handler.wfile = io.BytesIO()
    sent = {}
    handler.send_response = lambda code: sent.setdefault("code", code)
    handler.send_header = lambda k, v: None
    handler.end_headers = lambda: None
    handler.do_POST()
    return sent["code"], json.loads(handler.wfile.getvalue())


def _version():
    handler = ViewerHTTPRequestHandler.__new__(ViewerHTTPRequestHandler)
    handler.wfile = io.BytesIO()
    handler.send_response = lambda code: None
    handler.send_header = lambda k, v: None
    handler.end_headers = lambda: None
    handler.serve_version()
    return json.loads(handler.wfile.getvalue())


@pytest.fixture
def server():
    s = ViewerServer(config={"structures": []}, structure_files=[], port=8798)
    yield s
    s.clear_pick_request()


def test_pick_request_rides_on_version_and_the_click_comes_back(server):
    token = server.request_pick("bond", "Click the bond to make double", structure=1)
    assert _version()["pick_request"] == {"token": token, "kind": "bond",
                                          "prompt": "Click the bond to make double", "structure": 1}
    click = {"token": token, "kind": "bond", "structure": 1,
             "atom1": {"index": 3, "name": "C3"}, "atom2": {"index": 4, "name": "C4"}}
    waiting = {}
    thread = threading.Thread(target=lambda: waiting.setdefault("r", server.wait_for_pick(token, 5)))
    thread.start()
    assert _post("/pick", click) == (200, {"ok": True})
    thread.join()
    assert waiting["r"] == click
    assert "pick_request" not in _version()                  # answered requests are cleared


def test_a_click_for_another_request_is_refused(server):
    token = server.request_pick("atom", "Click an atom")
    code, body = _post("/pick", {"token": token + 1, "atom": {"index": 0}})
    assert code == 409 and body["ok"] is False
    assert server.wait_for_pick(token, 0.3) is None           # nothing arrived; request cleared


def test_cancel_comes_back_and_bad_kinds_are_refused(server):
    token = server.request_pick("atom", "Click an atom")
    _post("/pick", {"token": token, "cancelled": True})
    assert server.wait_for_pick(token, 1) == {"token": token, "cancelled": True}
    with pytest.raises(ValueError):
        server.request_pick("residue", "Click a residue")


def test_coordinator_pick_returns_none_without_a_viewer_and_on_cancel(monkeypatch, server):
    from proprep.structure_prep import interactive_structure_viewer as isv
    c = vc_mod.ViewerCoordinator()
    monkeypatch.setattr(isv, "_active_viewer_server", None, raising=False)
    assert c.pick("atom", "Click an atom") is None
    monkeypatch.setattr(isv, "_active_viewer_server", server, raising=False)
    monkeypatch.setattr(c, "is_running", lambda: True)
    monkeypatch.setattr(server, "wait_for_pick", lambda token, timeout: {"token": token, "cancelled": True})
    assert c.pick("atom", "Click an atom") is None
    monkeypatch.setattr(server, "wait_for_pick", lambda token, timeout: {"token": token, "atom": {"index": 2}})
    assert c.pick("atom", "Click an atom")["atom"] == {"index": 2}

    def interrupted(token, timeout):
        raise KeyboardInterrupt
    monkeypatch.setattr(server, "wait_for_pick", interrupted)
    assert c.pick("atom", "Click an atom") is None            # Ctrl-C in the terminal stops waiting


def test_page_banner_is_high_contrast_and_cancellable():
    assert 'id="pick-banner" role="status" aria-live="assertive"' in TEMPLATE
    assert "background:#000000; color:#ffffff;" in TEMPLATE          # 21:1
    assert "if (pickRequest && e.key === 'Escape') { e.preventDefault(); cancelPick(); return; }" in TEMPLATE
    assert "enterPickMode(vinfo.pick_request);" in TEMPLATE
    assert "exitPickMode();          // the command line stopped waiting" in TEMPLATE


def test_the_answering_click_does_not_also_measure_or_pick_density():
    assert TEMPLATE.count("if (pickRequest || clickTakenByPick)") == 2
    assert "stage.signals.clicked.add(onPickRequestClick);" in TEMPLATE
    assert "if (!pickingProxy.bond) {" in TEMPLATE                   # a bond pick needs a bond


def test_web_shell_forwards_the_pick():
    from proprep.web import server as web
    assert "/pick" in web._VIEWER_PROXY_BASE_PATHS
    assert any(getattr(route, "path", None) == "/pick" and "POST" in getattr(route, "methods", set())
               for route in web.app.routes)
