"""A user-requested view does not open a second browser tab when one is open.

An explicit view (force=True; the redox detector's view commands, a docking
pose asked for) used to relaunch the viewer with a new tab every time,
because an open tab could not follow a relaunch. Since an open tab reloads
itself on a new server instance, a new tab opened beside it was a duplicate
(seen by mgp with the redox detector). A tab counts as open when it polled
the server within ViewerServer.PAGE_OPEN_WINDOW seconds.
"""

import time
from unittest.mock import MagicMock

import pytest

from proprep.structure_prep import viewer_coordinator as vc_mod
from proprep.structure_prep.viewer_server import ViewerHTTPRequestHandler, ViewerServer


@pytest.fixture
def coordinator(monkeypatch):
    v = MagicMock()
    v.selected_structures = ["a.pdb"]
    v.annotation_config, v.viewer_config, v.shape_config, v.trajectory_files = {}, {}, {}, {}
    c = vc_mod.ViewerCoordinator()
    c._ensure_viewer = lambda: v
    c.is_running = lambda: True
    c._owns_live_server = lambda v: True          # the coordinator's own server
    monkeypatch.setattr(vc_mod, "_is_web_shell_mode", lambda: False)
    return c, v


def _page(monkeypatch, seconds_ago):
    monkeypatch.setattr(ViewerHTTPRequestHandler, "browser_opened_at", 0.0)
    monkeypatch.setattr(ViewerHTTPRequestHandler, "last_page_poll",
                        0.0 if seconds_ago is None else time.monotonic() - seconds_ago)


def test_page_open_follows_the_last_poll(monkeypatch):
    _page(monkeypatch, None)
    assert not ViewerServer.page_open()
    _page(monkeypatch, 1.0)
    assert ViewerServer.page_open()
    _page(monkeypatch, ViewerServer.PAGE_OPEN_WINDOW + 1)
    assert not ViewerServer.page_open()


@pytest.mark.parametrize("call", [
    lambda c: c.show_structure("b.pdb", force=True),
    lambda c: c.show_structures(["b.pdb", "pose.sdf"], force=True),
    lambda c: c.show_trajectory("b.pdb", "b.nc", force=True),
])
def test_a_forced_view_with_a_tab_open_opens_no_new_tab(coordinator, monkeypatch, call):
    c, v = coordinator
    _page(monkeypatch, 0.5)
    call(c)
    assert v._launch_viewer.call_args.kwargs["open_browser"] is False      # relaunched; the tab reloads itself


@pytest.mark.parametrize("call", [
    lambda c: c.show_structure("b.pdb", force=True),
    lambda c: c.show_structures(["b.pdb", "pose.sdf"], force=True),
])
def test_a_forced_view_with_no_tab_open_opens_one(coordinator, monkeypatch, call):
    c, v = coordinator
    _page(monkeypatch, None)
    call(c)
    assert v._launch_viewer.call_args.kwargs["open_browser"] is True


def test_the_same_structure_asked_for_again_with_a_tab_open_does_nothing(coordinator, monkeypatch):
    c, v = coordinator
    _page(monkeypatch, 0.5)
    c.show_structure("a.pdb", force=True)
    assert not v._launch_viewer.called
    _page(monkeypatch, None)                       # the user closed the tab
    c.show_structure("a.pdb", force=True)
    assert v._launch_viewer.call_args.kwargs["open_browser"] is True


def test_the_version_poll_records_the_page(monkeypatch):
    _page(monkeypatch, None)
    handler = ViewerHTTPRequestHandler.__new__(ViewerHTTPRequestHandler)
    handler.wfile = MagicMock()
    handler.send_response = handler.send_header = handler.end_headers = lambda *a, **k: None
    handler.serve_version()
    assert ViewerServer.page_open()


def test_a_server_started_is_stopped_by_the_next_launch_even_if_reporting_failed(monkeypatch, tmp_path):
    """The active server was recorded only after the success message; anything raising in
    between left it running unrecorded, and the next launch moved to another port."""
    import socket
    import webbrowser
    from proprep.structure_prep import interactive_structure_viewer as isv
    monkeypatch.setattr(webbrowser, "open", lambda *a, **k: True)
    monkeypatch.delenv("PROPREP_BATCH", raising=False)
    monkeypatch.setattr(isv, "_active_viewer_server", None)
    with socket.socket() as s:
        s.bind(("localhost", 0))
        free = s.getsockname()[1]

    v = isv.InteractiveStructureViewer.__new__(isv.InteractiveStructureViewer)
    v.selected_structures, v.available_annotations, v.annotation_config = [], {}, {}
    v.viewer_config, v.shape_config, v.trajectory_files, v.density_by_file = {}, {}, {}, {}
    v.scene_override, v.processor = None, None
    v._scene_dir = lambda: str(tmp_path)

    import proprep.structure_prep.viewer_server as vs
    real_init = vs.ViewerServer.__init__
    monkeypatch.setattr(vs.ViewerServer, "__init__", lambda self, *a, **k: real_init(self, *a, **{**k, "port": free}))
    monkeypatch.setattr(v, "_report_headless_access", lambda *a: (_ for _ in ()).throw(RuntimeError("boom")))
    assert v._launch_viewer(open_browser=False) is False          # the report raised after the start
    first = isv._active_viewer_server
    assert first is not None and first.is_running()                # but the server is recorded
    monkeypatch.setattr(v, "_report_headless_access", lambda *a: None)
    try:
        assert v._launch_viewer(open_browser=False) is True
        assert not first.is_running()                              # the next launch stopped it
        assert isv._active_viewer_server.port == free              # and reclaimed its port
    finally:
        isv._stop_active_viewer_server()


def test_a_tab_just_opened_counts_as_open_until_it_polls(monkeypatch):
    """A replayed session asks for views back to back: the tab opened for the first is
    still loading (no poll yet) when the second comes, and must not get a sibling."""
    now = time.monotonic()
    monkeypatch.setattr(ViewerHTTPRequestHandler, "last_page_poll", 0.0)
    monkeypatch.setattr(ViewerHTTPRequestHandler, "browser_opened_at", now - 2)
    assert ViewerServer.page_open()                                    # loading
    monkeypatch.setattr(ViewerHTTPRequestHandler, "browser_opened_at", now - ViewerServer.PAGE_LOAD_GRACE - 1)
    assert not ViewerServer.page_open()                                # never polled: gone
    monkeypatch.setattr(ViewerHTTPRequestHandler, "browser_opened_at", now - 20)
    monkeypatch.setattr(ViewerHTTPRequestHandler, "last_page_poll", now - 15)
    assert not ViewerServer.page_open()                                # it polled, then stopped: closed


def test_the_coordinator_takes_over_a_server_another_viewer_started(monkeypatch):
    """The main menu's Structure Viewer is its own viewer object. After it launched, the
    coordinator saw a server running and its structure unchanged and did nothing, and its
    highlights went to a server it never started: no module's views appeared (mgp, with the
    Redox Site Detector). It now relaunches its own viewer on that server's port, keeping
    its annotations, and the open tab reloads with them."""
    from proprep.structure_prep import interactive_structure_viewer as isv
    other_server = MagicMock()
    monkeypatch.setattr(isv, "_active_viewer_server", other_server)    # started by the module's viewer
    v = MagicMock()
    v.selected_structures = ["a.pdb"]
    v.annotation_config, v.viewer_config, v.shape_config, v.trajectory_files = {"site": {}}, {}, {}, {}
    v.server = None                                                    # the coordinator's viewer started nothing
    c = vc_mod.ViewerCoordinator()
    c._ensure_viewer = lambda: v
    c.is_running = lambda: True
    monkeypatch.setattr(vc_mod, "_is_web_shell_mode", lambda: False)
    _page(monkeypatch, 0.5)

    assert not c._owns_live_server(v)
    c.show_structure("a.pdb")                                          # same structure, not its server
    assert v._launch_viewer.call_args.kwargs == {"open_browser": False}
    assert v.annotation_config == {"site": {}}                         # kept
    v.reset_mock()
    c.highlight(":A and 5", label="x")
    assert v._launch_viewer.called and not v.update_annotations.called  # took over rather than update nothing

    v.server = other_server                                            # once it owns the live server
    v.reset_mock()
    c.highlight(":A and 6", label="y")
    assert v.update_annotations.called and not v._launch_viewer.called
