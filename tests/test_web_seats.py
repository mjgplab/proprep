"""Seat/session lifecycle (proprep.web.seats) with a fake PTY session."""

from __future__ import annotations

import asyncio
import json

import pytest

from proprep.web.pty_session import PtyExit
from proprep.web.seats import (
    Scrollback,
    Seat,
    SeatInfo,
    SeatRegistry,
    load_or_create_seats,
)


class FakeSession:
    """Stands in for PtySession: output is whatever the test feeds it."""

    instances: list["FakeSession"] = []

    def __init__(self, *, cwd, shell_url=None, transcript=True, cols=120, rows=32):
        self.cwd = cwd
        self.shell_url = shell_url
        self.transcript = transcript
        self.cols, self.rows = cols, rows
        self.transcript_path = None
        self.inputs: list[bytes] = []
        self.closed = False
        self.queue: asyncio.Queue = asyncio.Queue()
        FakeSession.instances.append(self)

    def start(self):
        self.started = True

    async def read_output(self):
        return await self.queue.get()

    def write_input(self, data: bytes):
        self.inputs.append(data)

    def resize(self, cols, rows):
        self.cols, self.rows = cols, rows

    async def close(self, *, grace: float = 5.0):
        self.closed = True

    # test helpers
    def emit(self, data: bytes):
        self.queue.put_nowait(data)

    def exit(self, code: int = 0):
        self.queue.put_nowait(PtyExit(code))


class FakeClient:
    """Records what a seat sends to an attached websocket."""

    def __init__(self):
        self.frames: list = []
        self.closed_with = None

    async def send_bytes(self, data: bytes):
        self.frames.append(data)

    async def send_text(self, data: str):
        self.frames.append(json.loads(data))

    async def close(self, code: int = 1000, reason=None):
        self.closed_with = (code, reason)


async def settle():
    for _ in range(5):
        await asyncio.sleep(0)


@pytest.fixture(autouse=True)
def _reset_instances():
    FakeSession.instances.clear()
    yield
    FakeSession.instances.clear()


# --- Scrollback ------------------------------------------------------------

def test_scrollback_bounded_and_cuts_at_newline():
    sb = Scrollback(limit=20)
    sb.append(b"first line\n")
    sb.append(b"second line\n")
    sb.append(b"third\n")
    snap = sb.snapshot()
    assert len(snap) <= 20
    assert not snap.startswith(b"ine")          # never starts mid-line
    assert snap.endswith(b"third\n")


def test_scrollback_clear():
    sb = Scrollback()
    sb.append(b"abc")
    sb.clear()
    assert len(sb) == 0 and sb.snapshot() == b""


# --- Seat --------------------------------------------------------------------

def test_attach_replays_then_streams_live(tmp_path):
    async def run():
        seat = Seat(SeatInfo("01", "tok", str(tmp_path / "s1")), session_factory=FakeSession)
        assert seat.ensure_session("http://127.0.0.1:1/seat/tok") is True
        sess = FakeSession.instances[-1]
        assert sess.shell_url == "http://127.0.0.1:1/seat/tok"
        sess.emit(b"hello ")
        sess.emit(b"world\n")
        await settle()

        c1 = FakeClient()
        await seat.attach(c1)
        assert c1.frames[0]["type"] == "attached"
        assert c1.frames[0]["replay"] is True and c1.frames[0]["alive"] is True
        assert c1.frames[1] == b"hello world\n"

        sess.emit(b"more\n")
        await settle()
        assert c1.frames[-1] == b"more\n"

        seat.write_input(b"q\n")
        assert sess.inputs == [b"q\n"]
        seat.resize(100, 40)
        assert (sess.cols, sess.rows) == (100, 40)
        await seat.close()
        assert sess.closed
    asyncio.run(run())


def test_detach_keeps_session_and_reattach_replays(tmp_path):
    async def run():
        seat = Seat(SeatInfo("01", "tok", str(tmp_path / "s1")), session_factory=FakeSession)
        seat.ensure_session(None)
        sess = FakeSession.instances[-1]
        c1 = FakeClient()
        await seat.attach(c1)
        sess.emit(b"before drop\n")
        await settle()
        seat.detach(c1)                       # network blip
        assert seat.alive and not sess.closed
        sess.emit(b"while away\n")            # output while nobody is attached
        await settle()
        assert seat.ensure_session(None) is False   # still the same child
        c2 = FakeClient()
        await seat.attach(c2)
        assert c2.frames[0]["replay"] is True
        assert c2.frames[1] == b"before drop\nwhile away\n"
        await seat.close()
    asyncio.run(run())


def test_second_client_displaces_first(tmp_path):
    async def run():
        seat = Seat(SeatInfo("01", "tok", str(tmp_path / "s1")), session_factory=FakeSession)
        seat.ensure_session(None)
        sess = FakeSession.instances[-1]
        c1, c2 = FakeClient(), FakeClient()
        await seat.attach(c1)
        await seat.attach(c2)
        assert c1.frames[-1] == {"type": "displaced"}
        assert c1.closed_with[0] == 4001
        sess.emit(b"x")
        await settle()
        assert c2.frames[-1] == b"x" and b"x" not in c1.frames
        await seat.close()
    asyncio.run(run())


def test_exit_notifies_client_and_next_ensure_starts_fresh(tmp_path):
    async def run():
        seat = Seat(SeatInfo("01", "tok", str(tmp_path / "s1")), session_factory=FakeSession)
        seat.ensure_session(None)
        first = FakeSession.instances[-1]
        c1 = FakeClient()
        await seat.attach(c1)
        first.emit(b"bye\n")
        first.exit(3)
        await settle()
        assert seat.alive is False and seat.exit_code == 3
        assert c1.frames[-1] == {"type": "exit", "returncode": 3}
        # A late attach reports the dead session honestly...
        c2 = FakeClient()
        await seat.attach(c2)
        assert c2.frames[0]["alive"] is False and c2.frames[0]["returncode"] == 3
        # ...and a new ensure starts a fresh child with an empty scrollback.
        assert seat.ensure_session(None) is True
        assert FakeSession.instances[-1] is not first
        assert seat.scrollback.snapshot() == b"" and seat.alive
        await seat.close()
    asyncio.run(run())


def test_output_to_dead_client_detaches_it(tmp_path):
    class BrokenClient(FakeClient):
        async def send_bytes(self, data):
            raise RuntimeError("socket gone")

    async def run():
        seat = Seat(SeatInfo("01", "tok", str(tmp_path / "s1")), session_factory=FakeSession)
        seat.ensure_session(None)
        sess = FakeSession.instances[-1]
        bad = BrokenClient()
        await seat.attach(bad)
        sess.emit(b"x")
        await settle()
        assert seat.attached is None
        assert seat.scrollback.snapshot() == b"x"   # still buffered for the next client
        await seat.close()
    asyncio.run(run())


def test_ensure_session_creates_cwd(tmp_path):
    async def run():
        cwd = tmp_path / "deep" / "seat"
        seat = Seat(SeatInfo("01", "tok", str(cwd)), session_factory=FakeSession)
        seat.ensure_session(None)
        assert cwd.is_dir()
        await seat.close()
    asyncio.run(run())


# --- Registry + persistence --------------------------------------------------

def test_load_or_create_seats_is_stable_across_restarts(tmp_path):
    first = load_or_create_seats(tmp_path / "seats", 3)
    assert [s.name for s in first] == ["01", "02", "03"]
    assert len({s.token for s in first}) == 3
    assert all((tmp_path / "seats" / f"seat_{s.name}").is_dir() for s in first)
    again = load_or_create_seats(tmp_path / "seats", 3)
    assert [(s.name, s.token) for s in again] == [(s.name, s.token) for s in first]
    grown = load_or_create_seats(tmp_path / "seats", 5)
    assert [s.token for s in grown[:3]] == [s.token for s in first]
    assert len(grown) == 5
    shrunk = load_or_create_seats(tmp_path / "seats", 2)
    assert len(shrunk) == 5, "handed-out links must survive a smaller restart"
    data = json.loads((tmp_path / "seats" / "seats.json").read_text())
    assert len(data["seats"]) == 5


def test_registry_lookup_hosted_and_local(tmp_path):
    seats = load_or_create_seats(tmp_path / "seats", 2)
    reg = SeatRegistry(seats, local=False, session_factory=FakeSession)
    assert reg.lookup(seats[1].token).name == "02"
    assert reg.lookup("nope") is None
    assert reg.lookup(None) is None
    assert reg.lookup("") is None

    local = SeatRegistry([SeatInfo("local", None, str(tmp_path))], local=True, session_factory=FakeSession)
    assert local.lookup(None) is local.default
    assert local.lookup("anything") is local.default
    with pytest.raises(ValueError):
        SeatRegistry(seats, local=True, session_factory=FakeSession)
