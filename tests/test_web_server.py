"""HTTP/websocket behaviour of proprep.web.server in local and hosted mode.

Uses Starlette's TestClient (real ASGI app, in-process) with a fake PTY
session so nothing forks ProPrep.
"""

from __future__ import annotations

import asyncio
import json
import os
import zipfile
from io import BytesIO

import pytest
from starlette.testclient import TestClient

from proprep.web import server
from proprep.web.pty_session import PtyExit
from proprep.web.seats import SeatInfo, load_or_create_seats


class FakeSession:
    instances: list["FakeSession"] = []

    def __init__(self, *, cwd, shell_url=None, transcript=True, cols=120, rows=32):
        self.cwd, self.shell_url = cwd, shell_url
        self.transcript_path = None
        self.inputs: list[bytes] = []
        self.closed = False
        self.queue: asyncio.Queue = asyncio.Queue()
        FakeSession.instances.append(self)

    def start(self):
        pass

    async def read_output(self):
        return await self.queue.get()

    def write_input(self, data):
        self.inputs.append(data)

    def resize(self, cols, rows):
        self.size = (cols, rows)

    async def close(self, *, grace: float = 5.0):
        self.closed = True

    def emit(self, data: bytes):
        self.queue.put_nowait(data)

    def exit(self, code=0):
        self.queue.put_nowait(PtyExit(code))


# The announce endpoints only accept loopback callers; TestClient's default
# peer address is "testclient", so fixtures present as 127.0.0.1.
LOOPBACK = ("127.0.0.1", 50000)


@pytest.fixture(autouse=True)
def _clean(monkeypatch):
    FakeSession.instances.clear()
    monkeypatch.delenv("PROPREP_WEB_SEATS_FILE", raising=False)
    monkeypatch.delenv("PROPREP_WEB_BIND", raising=False)
    monkeypatch.setenv("PROPREP_WEB_NO_TRANSCRIPT", "1")
    yield
    FakeSession.instances.clear()


@pytest.fixture
def local_client(tmp_path, monkeypatch):
    monkeypatch.setenv("PROPREP_WEB_BIND", "http://127.0.0.1:8123")
    server.configure_registry(
        seats=[SeatInfo("local", None, str(tmp_path / "proj"))], local=True, session_factory=FakeSession,
    )
    with TestClient(server.app, client=LOOPBACK) as client:
        yield client


@pytest.fixture
def hosted(tmp_path, monkeypatch):
    seats = load_or_create_seats(tmp_path / "seats", 2)
    monkeypatch.setenv("PROPREP_WEB_SEATS_FILE", str(tmp_path / "seats" / "seats.json"))
    monkeypatch.setenv("PROPREP_WEB_BIND", "http://127.0.0.1:8123")
    server.configure_registry(seats=seats, local=False, session_factory=FakeSession)
    with TestClient(server.app, client=LOOPBACK) as client:
        yield client, seats


def recv_until_bytes(ws, limit=10):
    """Return (control_frames, first_binary_frame)."""
    controls = []
    for _ in range(limit):
        msg = ws.receive()
        if "bytes" in msg and msg["bytes"] is not None:
            return controls, msg["bytes"]
        if "text" in msg and msg["text"] is not None:
            controls.append(json.loads(msg["text"]))
    return controls, None


# --- local mode ----------------------------------------------------------------

def test_local_index_and_seat_info_need_no_token(local_client):
    assert local_client.get("/").status_code == 200
    info = local_client.get("/seat-info").json()
    assert info == {"name": "local", "hosted": False, "alive": False, "cwd": info["cwd"]}
    assert local_client.get("/healthz").text == "ok"


def test_local_terminal_attach_input_and_replay_on_reconnect(local_client):
    with local_client.websocket_connect("/ws/term") as ws:
        first = json.loads(ws.receive_text())
        assert first["type"] == "attached" and first["replay"] is False and first["alive"] is True
        sess = FakeSession.instances[-1]
        # child announces to the loopback bind address, never the Host header
        assert sess.shell_url == "http://127.0.0.1:8123"
        sess.emit(b"prompt> ")
        assert ws.receive_bytes() == b"prompt> "
        ws.send_bytes(b"1\n")
        ws.send_text(json.dumps({"type": "resize", "cols": 90, "rows": 30}))
    # socket closed; the child must still be alive
    assert not sess.closed
    assert sess.inputs == [b"1\n"] and sess.size == (90, 30)

    with local_client.websocket_connect("/ws/term") as ws:
        controls, replay = recv_until_bytes(ws)
        assert controls[0]["type"] == "attached" and controls[0]["replay"] is True
        assert replay == b"prompt> "
        assert FakeSession.instances[-1] is sess, "reconnect re-attaches, does not restart"


def test_local_exit_then_reconnect_starts_fresh_session(local_client):
    with local_client.websocket_connect("/ws/term") as ws:
        ws.receive_text()
        sess = FakeSession.instances[-1]
        sess.emit(b"bye\n")
        assert ws.receive_bytes() == b"bye\n"
        sess.exit(0)
        assert json.loads(ws.receive_text()) == {"type": "exit", "returncode": 0}
    with local_client.websocket_connect("/ws/term") as ws:
        first = json.loads(ws.receive_text())
        assert first["replay"] is False and first["alive"] is True
        assert FakeSession.instances[-1] is not sess


def test_local_arms_auto_shutdown_only_after_last_socket(local_client, monkeypatch):
    kills = []
    monkeypatch.setattr(server.os, "kill", lambda pid, sig: kills.append((pid, sig)))
    monkeypatch.setattr(server, "_SHUTDOWN_GRACE_SECONDS", 0.05)
    with local_client.websocket_connect("/ws/term") as ws:
        ws.receive_text()
    assert server._shutdown_task is not None
    # give the grace timer a chance to fire inside the app's loop
    import time
    for _ in range(50):
        if kills:
            break
        time.sleep(0.02)
        local_client.get("/healthz")
    assert kills and kills[0][1] == server.signal.SIGINT


def test_local_announce_and_proxy_503_without_viewer(local_client):
    assert local_client.get("/viewer").status_code == 503
    r = local_client.post("/internal/viewer-announce", json={"port": 9999})
    assert r.status_code == 204
    assert server.registry().default.viewer_port == 9999
    assert local_client.post("/internal/viewer-announce", json={"port": 0}).status_code == 400


def test_announce_rejects_non_loopback_callers(tmp_path):
    server.configure_registry(
        seats=[SeatInfo("local", None, str(tmp_path))], local=True, session_factory=FakeSession,
    )
    with TestClient(server.app, client=("203.0.113.9", 4444)) as remote:
        assert remote.post("/internal/viewer-announce", json={"port": 9999}).status_code == 403
    assert server.registry().default.viewer_port is None


# --- hosted mode ---------------------------------------------------------------

def test_hosted_requires_token_and_sets_cookie(hosted):
    client, seats = hosted
    assert client.get("/").status_code == 403
    assert client.get("/?seat=bogus").status_code == 403
    assert client.get("/seat-info").status_code == 403
    r = client.get(f"/?seat={seats[1].token}", follow_redirects=False)
    assert r.status_code == 303 and r.headers["location"] == "/"
    assert client.cookies.get(server.SEAT_COOKIE) == seats[1].token
    # cookie now carries the seat for everything else
    assert client.get("/").status_code == 200
    info = client.get("/seat-info").json()
    assert info["name"] == "02" and info["hosted"] is True


def test_hosted_websocket_rejects_missing_token(hosted):
    client, seats = hosted
    from starlette.websockets import WebSocketDisconnect
    with pytest.raises(WebSocketDisconnect) as exc:
        with client.websocket_connect("/ws/term") as ws:
            ws.receive_text()
    assert exc.value.code == 4003


def test_hosted_seats_are_isolated_and_announce_per_seat(hosted):
    client, seats = hosted
    a = TestClient(server.app, client=LOOPBACK)
    b = TestClient(server.app, client=LOOPBACK)
    a.get(f"/?seat={seats[0].token}")
    b.get(f"/?seat={seats[1].token}")
    with a.websocket_connect("/ws/term") as wa, b.websocket_connect("/ws/term") as wb:
        wa.receive_text(); wb.receive_text()
        sa, sb = FakeSession.instances[-2], FakeSession.instances[-1]
        assert sa.cwd.endswith("seat_01") and sb.cwd.endswith("seat_02")
        assert sa.shell_url == f"http://127.0.0.1:8123/seat/{seats[0].token}"
        sa.emit(b"only A")
        assert wa.receive_bytes() == b"only A"
        # seat-scoped announce routes the viewer proxy per seat
        r = a.post(f"/seat/{seats[0].token}/internal/viewer-announce", json={"port": 7001})
        assert r.status_code == 204
        assert server.registry().lookup(seats[0].token).viewer_port == 7001
        assert server.registry().lookup(seats[1].token).viewer_port is None
        assert b.get("/viewer").status_code == 503      # B has no viewer yet
        assert a.post("/seat/nope/internal/viewer-announce", json={"port": 1}).status_code == 404
    assert not sa.closed and not sb.closed


def test_hosted_never_arms_auto_shutdown(hosted, monkeypatch):
    client, seats = hosted
    client.get(f"/?seat={seats[0].token}")
    with client.websocket_connect("/ws/term") as ws:
        ws.receive_text()
    assert server._shutdown_task is None


def test_hosted_second_tab_displaces_first(hosted):
    client, seats = hosted
    client.get(f"/?seat={seats[0].token}")
    with client.websocket_connect("/ws/term") as w1:
        w1.receive_text()
        with client.websocket_connect("/ws/term") as w2:
            w2.receive_text()
            assert json.loads(w1.receive_text()) == {"type": "displaced"}


def test_download_zips_seat_directory(hosted):
    client, seats = hosted
    client.get(f"/?seat={seats[0].token}")
    root = seats[0].cwd
    os.makedirs(os.path.join(root, "sub"), exist_ok=True)
    with open(os.path.join(root, "a.txt"), "w") as fh:
        fh.write("alpha")
    with open(os.path.join(root, "sub", "b.pdb"), "w") as fh:
        fh.write("ATOM")
    r = client.get("/download")
    assert r.status_code == 200
    assert r.headers["content-type"] == "application/zip"
    assert "proprep_01_" in r.headers["content-disposition"]
    with zipfile.ZipFile(BytesIO(r.content)) as zf:
        assert sorted(zf.namelist()) == ["a.txt", "sub/b.pdb"]
        assert zf.read("a.txt") == b"alpha"


def test_configure_registry_from_seats_file_env(tmp_path, monkeypatch):
    seats = load_or_create_seats(tmp_path / "seats", 1)
    monkeypatch.setenv("PROPREP_WEB_SEATS_FILE", str(tmp_path / "seats" / "seats.json"))
    reg = server.configure_registry(session_factory=FakeSession)
    assert reg.local is False and reg.lookup(seats[0].token).name == "01"
    monkeypatch.delenv("PROPREP_WEB_SEATS_FILE")
    reg = server.configure_registry(session_factory=FakeSession)
    assert reg.local is True and reg.default.name == "local"
