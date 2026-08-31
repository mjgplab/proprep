"""Seats: PTY sessions that outlive the websockets attached to them.

Phase A tied a :class:`~proprep.web.pty_session.PtySession` to one terminal
websocket: the child was forked when the socket opened and killed when it
closed. That is right for a laptop (``127.0.0.1``, close the tab and the
shell goes away) and wrong for a hosted workshop, where a Wi-Fi blip or a
page reload must not destroy an attendee's session.

Phase B introduces the *seat*: a long-lived holder of one PTY session, its
recent output (scrollback for replay on reconnect), the viewer port its
child announced, and whichever websocket is currently attached. Websockets
come and go; the seat and its child stay until the child exits or the seat
is closed.

Two deployment shapes share this code:

* **Local mode** (default ``proprep-web``): one implicit seat with no token,
  working directory = the launch directory. The server still auto-shuts
  down shortly after the last browser leaves, so behaviour on a laptop is
  unchanged apart from reloads now *reattaching* instead of restarting.
* **Hosted mode** (``proprep-web --seats N``): N seats, each with its own
  working directory and a random URL token persisted in ``seats.json`` so
  links survive restarts. The token is the only access control, which is
  what a one-off workshop wants: no accounts, one link per attendee.
"""

from __future__ import annotations

import asyncio
import hmac
import json
import logging
import secrets
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Awaitable, Callable, Optional, Protocol

from proprep.web.pty_session import PtyExit, PtySession

logger = logging.getLogger(__name__)

SEATS_FILE_NAME = "seats.json"
DEFAULT_SCROLLBACK_BYTES = 256 * 1024


# ---------------------------------------------------------------------------
# Session protocol (so tests can substitute a fake for the real PTY)
# ---------------------------------------------------------------------------

class SessionLike(Protocol):
    transcript_path: Optional[Path]

    def start(self) -> None: ...
    async def read_output(self) -> bytes | PtyExit: ...
    def write_input(self, data: bytes) -> None: ...
    def resize(self, cols: int, rows: int) -> None: ...
    async def close(self, *, grace: float = 5.0) -> None: ...


SessionFactory = Callable[..., SessionLike]


# ---------------------------------------------------------------------------
# Scrollback
# ---------------------------------------------------------------------------

class Scrollback:
    """Bounded byte buffer of recent PTY output, replayed on (re)attach.

    Raw bytes (ANSI included) so xterm.js re-renders colours and layout.
    When the head is trimmed we cut forward to the next newline so a replay
    does not start inside an escape sequence.
    """

    def __init__(self, limit: int = DEFAULT_SCROLLBACK_BYTES) -> None:
        self._buf = bytearray()
        self._limit = limit

    def append(self, data: bytes) -> None:
        self._buf += data
        excess = len(self._buf) - self._limit
        if excess > 0:
            nl = self._buf.find(b"\n", excess)
            cut = nl + 1 if nl != -1 else excess
            del self._buf[:cut]

    def snapshot(self) -> bytes:
        return bytes(self._buf)

    def clear(self) -> None:
        self._buf.clear()

    def __len__(self) -> int:
        return len(self._buf)


# ---------------------------------------------------------------------------
# Attached-client protocol: the seat only needs to send to it
# ---------------------------------------------------------------------------

class ClientLike(Protocol):
    async def send_bytes(self, data: bytes) -> None: ...
    async def send_text(self, data: str) -> None: ...
    async def close(self, code: int = 1000, reason: Optional[str] = None) -> None: ...


# ---------------------------------------------------------------------------
# Seat
# ---------------------------------------------------------------------------

@dataclass
class SeatInfo:
    name: str
    token: Optional[str]
    cwd: str

    def to_json(self) -> dict[str, Any]:
        return {"name": self.name, "token": self.token, "cwd": self.cwd}


class Seat:
    """One long-lived PTY session plus the state needed to re-attach to it."""

    def __init__(
        self,
        info: SeatInfo,
        *,
        session_factory: SessionFactory = PtySession,
        transcript: bool = True,
        scrollback_bytes: int = DEFAULT_SCROLLBACK_BYTES,
    ) -> None:
        self.info = info
        self._session_factory = session_factory
        self._transcript = transcript
        self.session: Optional[SessionLike] = None
        self.scrollback = Scrollback(scrollback_bytes)
        self.attached: Optional[ClientLike] = None
        self.control_clients: set[Any] = set()
        self.viewer_port: Optional[int] = None
        self.exit_code: Optional[int] = None
        self.last_detach: Optional[float] = None
        self.cols = 120
        self.rows = 32
        self._pump_task: Optional[asyncio.Task] = None
        # Serialises "append to scrollback + forward to client" against
        # "snapshot scrollback + send replay + attach", so a chunk can never
        # slip between the snapshot and the replay.
        self._out_lock = asyncio.Lock()

    # ----- identity ---------------------------------------------------------

    @property
    def name(self) -> str:
        return self.info.name

    @property
    def token(self) -> Optional[str]:
        return self.info.token

    @property
    def cwd(self) -> str:
        return self.info.cwd

    @property
    def alive(self) -> bool:
        return self.session is not None and self.exit_code is None

    # ----- lifecycle --------------------------------------------------------

    def ensure_session(self, shell_url: Optional[str]) -> bool:
        """Start the child if there is none (or the previous one exited).

        Returns True if a new child was started. A fresh child starts with an
        empty scrollback: the previous session's output belongs to the
        transcript, not to the new terminal.
        """
        if self.alive:
            return False
        if self._pump_task is not None and not self._pump_task.done():
            # The old pump is still winding down after an exit; let it finish
            # before we replace the session it is reading from.
            self._pump_task.cancel()
        self.scrollback.clear()
        self.exit_code = None
        self.viewer_port = None
        Path(self.cwd).mkdir(parents=True, exist_ok=True)
        session = self._session_factory(
            cwd=self.cwd, shell_url=shell_url, transcript=self._transcript,
            cols=self.cols, rows=self.rows,
        )
        session.start()
        self.session = session
        self._pump_task = asyncio.create_task(self._pump(session), name=f"seat-pump-{self.name}")
        logger.info("Seat %s: started session (cwd=%s)", self.name, self.cwd)
        return True

    async def _pump(self, session: SessionLike) -> None:
        """Drain the PTY forever, buffering and forwarding to the attached client."""
        while True:
            chunk = await session.read_output()
            if isinstance(chunk, PtyExit):
                self.exit_code = chunk.returncode
                logger.info("Seat %s: child exited with %d", self.name, chunk.returncode)
                client = self.attached
                if client is not None:
                    try:
                        await client.send_text(json.dumps({"type": "exit", "returncode": chunk.returncode}))
                    except Exception:
                        pass
                return
            async with self._out_lock:
                self.scrollback.append(chunk)
                client = self.attached
                if client is not None:
                    try:
                        await client.send_bytes(chunk)
                    except Exception:
                        # The socket is gone; the handler's finally will
                        # detach, but don't keep hammering it meanwhile.
                        if self.attached is client:
                            self.attached = None

    async def attach(self, client: ClientLike) -> None:
        """Make ``client`` the live terminal for this seat.

        Any previously attached client is closed (one connection per seat,
        as in Phase A). The new client receives an ``attached`` control
        frame, then the scrollback replay, then live output.
        """
        old = self.attached
        if old is not None and old is not client:
            self.attached = None
            try:
                await old.send_text(json.dumps({"type": "displaced"}))
            except Exception:
                pass
            try:
                await old.close(code=4001, reason="displaced by a newer connection")
            except Exception:
                pass
        async with self._out_lock:
            replay = self.scrollback.snapshot()
            await client.send_text(json.dumps({
                "type": "attached",
                "seat": self.name,
                "replay": len(replay) > 0,
                "alive": self.alive,
                "returncode": self.exit_code,
            }))
            if replay:
                await client.send_bytes(replay)
            self.attached = client
        self.last_detach = None

    def detach(self, client: ClientLike) -> None:
        if self.attached is client:
            self.attached = None
            self.last_detach = time.monotonic()

    # ----- I/O passthrough --------------------------------------------------

    def write_input(self, data: bytes) -> None:
        if self.session is not None and self.exit_code is None:
            self.session.write_input(data)

    def resize(self, cols: int, rows: int) -> None:
        self.cols, self.rows = max(1, int(cols)), max(1, int(rows))
        if self.session is not None and self.exit_code is None:
            self.session.resize(self.cols, self.rows)

    async def close(self) -> None:
        """Kill the child (if any) and drop the attached client."""
        if self._pump_task is not None:
            self._pump_task.cancel()
            try:
                await self._pump_task
            except (asyncio.CancelledError, Exception):
                pass
            self._pump_task = None
        if self.session is not None:
            try:
                await self.session.close()
            except Exception:
                pass
            self.session = None
        self.attached = None


# ---------------------------------------------------------------------------
# Registry + persistence
# ---------------------------------------------------------------------------

def new_token() -> str:
    # 12 bytes -> 16 URL-safe chars: unguessable, still typable from a slide.
    return secrets.token_urlsafe(12)


def load_or_create_seats(seats_dir: Path, count: int) -> list[SeatInfo]:
    """Return ``count`` seats under ``seats_dir``, reusing ``seats.json``.

    Existing seats keep their tokens (links stay valid across restarts);
    missing ones are appended. ``count`` smaller than the file's length
    keeps the extra seats too: a link that was handed out must not die
    because the instructor restarted with a smaller number.
    """
    seats_dir = Path(seats_dir).expanduser().resolve()
    seats_dir.mkdir(parents=True, exist_ok=True)
    path = seats_dir / SEATS_FILE_NAME
    seats: list[SeatInfo] = []
    if path.is_file():
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            for entry in data.get("seats", []):
                seats.append(SeatInfo(str(entry["name"]), str(entry["token"]), str(entry["cwd"])))
        except (OSError, ValueError, KeyError, TypeError) as exc:
            logger.warning("Ignoring unreadable %s: %s", path, exc)
            seats = []
    width = max(2, len(str(count)))
    existing = {s.name for s in seats}
    n = 1
    while len(seats) < count:
        name = f"{n:0{width}d}"
        n += 1
        if name in existing:
            continue
        seats.append(SeatInfo(name, new_token(), str(seats_dir / f"seat_{name}")))
    for s in seats:
        Path(s.cwd).mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"seats": [s.to_json() for s in seats]}, indent=2) + "\n", encoding="utf-8")
    try:
        path.chmod(0o600)  # tokens are credentials
    except OSError:
        pass
    return seats


class SeatRegistry:
    """All seats of one server process.

    ``local`` (no tokens) resolves every request to the single default seat;
    hosted mode resolves by token with a constant-time compare and returns
    None for anything else.
    """

    def __init__(
        self,
        seats: list[SeatInfo],
        *,
        local: bool,
        session_factory: SessionFactory = PtySession,
        transcript: bool = True,
    ) -> None:
        if local and len(seats) != 1:
            raise ValueError("local mode has exactly one seat")
        self.local = local
        self._seats = [
            Seat(info, session_factory=session_factory, transcript=transcript) for info in seats
        ]
        self._by_token = {s.token: s for s in self._seats if s.token}

    @property
    def seats(self) -> list[Seat]:
        return list(self._seats)

    @property
    def default(self) -> Seat:
        return self._seats[0]

    def lookup(self, token: Optional[str]) -> Optional[Seat]:
        if self.local:
            return self.default
        if not token:
            return None
        for candidate, seat in self._by_token.items():
            if hmac.compare_digest(candidate, token):
                return seat
        return None

    async def close_all(self) -> None:
        for seat in self._seats:
            await seat.close()
