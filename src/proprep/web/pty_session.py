"""PTY-backed ProPrep subprocess.

Each browser tab gets one PtySession. The child execs ``python -m proprep``
inside a PTY so Rich's ANSI rendering survives end-to-end into xterm.js.
The subprocess sees ``PROPREP_WEB_SHELL=1`` so later phases can detect the
shell and route the structure viewer's HTTP server through the same web
process instead of opening a new browser tab.
"""

from __future__ import annotations

import asyncio
import datetime
import errno
import fcntl
import logging
import os
import pty
import re
import signal
import struct
import sys
import termios
from dataclasses import dataclass
from pathlib import Path
from typing import IO

logger = logging.getLogger(__name__)

# Bytes-per-read from the master fd. xterm.js handles short writes fine, so
# small chunks reduce input-to-render latency. 4 KiB is a reasonable balance.
READ_CHUNK = 4096

# Strip ANSI control sequences so the on-disk transcript is plain, greppable
# text rather than raw Rich output. Covers CSI (colors, cursor moves), OSC
# (window-title) and lone two-byte escapes — everything ProPrep's menus and
# result blocks emit. Bare carriage returns from progress redraws are left
# intact; they are rare in ProPrep's line-oriented output.
_ANSI_RE = re.compile(
    rb"\x1b\[[0-?]*[ -/]*[@-~]"            # CSI ... final byte
    rb"|\x1b\][^\x07\x1b]*(?:\x07|\x1b\\)"  # OSC ... BEL or ST terminated
    rb"|\x1b[@-OQ-WYZ\\]"                   # two-byte Fe escape, EXCLUDING the
                                            # string introducers P/X/]/^/_ (those
                                            # start OSC/DCS/SOS/PM/APC strings and
                                            # must not be eaten as a 2-char escape,
                                            # esp. mid-stream before the body/ST)
)

def _strip_ansi_stream(carry: bytes, data: bytes) -> tuple[bytes, bytes]:
    """Strip ANSI escapes from one stream chunk, boundary-safe.

    ``.sub()`` per chunk cannot see an escape that straddles a read boundary, so
    it would emit both halves verbatim and they would reassemble in the file.
    This prepends the previous chunk's held-back tail, strips complete escapes,
    then holds back any still-incomplete trailing escape.

    After stripping, ``_ANSI_RE`` has removed every *complete* escape, so any ESC
    still present must belong to a single trailing escape that the read boundary
    cut short. Everything from the first remaining ESC to the end is that partial
    (an unfinished CSI, or an OSC whose ST terminator — itself ESC + ``\\`` — was
    split), so we carry it whole into the next chunk. Any leftover at EOF is an
    incomplete escape and is correctly dropped.

    Returns ``(bytes_to_write, carry_for_next_chunk)``.
    """
    clean = _ANSI_RE.sub(b"", carry + data)
    idx = clean.find(0x1b)
    if idx == -1:
        return clean, b""
    return clean[:idx], clean[idx:]


def _default_transcript_path() -> Path:
    """Return a fresh, unique transcript path under ~/.proprep/web_sessions/.

    Package dirs are read-only for us (see project conventions), so session
    logs live under the user's ~/.proprep. Microsecond precision keeps two
    tabs opened in the same second from colliding.
    """
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    log_dir = Path.home() / ".proprep" / "web_sessions"
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir / f"session_{stamp}.log"


@dataclass
class PtyExit:
    """Sentinel sent on the output queue when the child terminates."""
    returncode: int


class PtySession:
    """A single PTY-backed ProPrep child.

    Producer side (PTY -> websocket): a thread runs a blocking ``os.read``
    loop on the master fd and pushes bytes onto an asyncio queue via
    ``loop.call_soon_threadsafe``. (Asyncio's ``add_reader`` would also
    work, but on macOS with kqueue the readiness signaling for PTY masters
    has historically been quirky; a dedicated thread is more portable.)

    Consumer side (websocket -> PTY): the websocket task calls
    ``write_input`` / ``resize`` directly; ``os.write`` is non-blocking
    enough on a master fd that it does not need its own thread.
    """

    def __init__(
        self,
        *,
        cols: int = 120,
        rows: int = 32,
        cwd: str | None = None,
        shell_url: str | None = None,
        transcript: bool = True,
    ) -> None:
        self._cols = cols
        self._rows = rows
        self._cwd = cwd or os.getcwd()
        self._shell_url = shell_url
        self._transcript_enabled = transcript
        self._transcript_file: IO[bytes] | None = None
        # Populated by start() when a transcript is opened; read by the server
        # so it can tell the user where the session log landed.
        self.transcript_path: Path | None = None
        self._pid: int | None = None
        self._master_fd: int | None = None
        self._loop: asyncio.AbstractEventLoop | None = None
        self._out_queue: asyncio.Queue[bytes | PtyExit] | None = None
        self._reader_thread: "asyncio.Future | None" = None
        self._closed = False

    # ----- lifecycle --------------------------------------------------------

    def start(self) -> None:
        """Fork the PTY child. Must be called from within an asyncio loop."""
        if self._pid is not None:
            raise RuntimeError("PtySession already started")
        self._loop = asyncio.get_running_loop()
        self._out_queue = asyncio.Queue(maxsize=256)

        pid, master_fd = pty.fork()
        if pid == 0:
            # Child: replace this process with `python -m proprep`.
            try:
                os.chdir(self._cwd)
            except OSError:
                pass
            env = os.environ.copy()
            env["PROPREP_WEB_SHELL"] = "1"
            if self._shell_url:
                # ViewerServer.start() POSTs {port:N} to
                # ``$PROPREP_WEB_SHELL_URL/internal/viewer-announce`` once
                # it's listening, so the shell can iframe the viewer.
                env["PROPREP_WEB_SHELL_URL"] = self._shell_url
            env["TERM"] = env.get("TERM", "xterm-256color")
            env["COLORTERM"] = env.get("COLORTERM", "truecolor")
            # Force unbuffered stdio so prompts/output reach the master fd
            # immediately.
            env["PYTHONUNBUFFERED"] = "1"
            try:
                # ``python -m proprep.main`` matches the existing if-__main__
                # guard in proprep/main.py. (The package has no __main__.py,
                # so ``python -m proprep`` would silently no-op.)
                os.execvpe(sys.executable, [sys.executable, "-m", "proprep.main"], env)
            except Exception as exc:  # pragma: no cover - exec failure path
                # Write to fd 2 (the slave's stderr) so the parent sees it.
                os.write(2, f"proprep.web: failed to exec child: {exc}\n".encode())
                os._exit(127)

        # Parent
        self._pid = pid
        self._master_fd = master_fd
        # Apply initial window size before any output is produced.
        self._apply_winsize(self._cols, self._rows)

        # Open the plain-text session transcript before any output flows, so
        # the reader loop can tee into it from the very first byte. A failure
        # here is non-fatal: the shell still works, just without a log.
        if self._transcript_enabled:
            try:
                path = _default_transcript_path()
                self._transcript_file = open(path, "wb")
                self.transcript_path = path
                header = (
                    f"# ProPrep web-shell session transcript\n"
                    f"# started: {datetime.datetime.now().isoformat(timespec='seconds')}\n"
                    f"# cwd:     {self._cwd}\n\n"
                ).encode()
                self._transcript_file.write(header)
                self._transcript_file.flush()
            except OSError as exc:
                logger.warning("Could not open session transcript: %s", exc)
                self._transcript_file = None
                self.transcript_path = None

        # Spawn the reader in a worker thread.
        self._reader_thread = self._loop.run_in_executor(None, self._reader_loop)
        logger.info("PtySession started: pid=%d cols=%d rows=%d", pid, self._cols, self._rows)

    async def close(self, *, grace: float = 5.0) -> None:
        """Terminate the child and clean up. Idempotent."""
        if self._closed:
            return
        self._closed = True

        pid = self._pid
        if pid is not None:
            for sig, label in ((signal.SIGTERM, "SIGTERM"), (signal.SIGKILL, "SIGKILL")):
                try:
                    os.kill(pid, sig)
                except ProcessLookupError:
                    break
                # Wait up to `grace` seconds (only on the SIGTERM pass).
                if sig is signal.SIGTERM:
                    deadline = asyncio.get_event_loop().time() + grace
                    while asyncio.get_event_loop().time() < deadline:
                        try:
                            wpid, _ = os.waitpid(pid, os.WNOHANG)
                        except ChildProcessError:
                            break
                        if wpid == pid:
                            break
                        await asyncio.sleep(0.05)
                    else:
                        continue  # SIGTERM didn't take, fall through to KILL
                    break

        if self._master_fd is not None:
            try:
                os.close(self._master_fd)
            except OSError:
                pass
            self._master_fd = None

        if self._reader_thread is not None:
            try:
                await asyncio.wait_for(self._reader_thread, timeout=1.0)
            except (asyncio.TimeoutError, Exception):
                pass

        # Best-effort reap (in case SIGKILL happened above).
        if pid is not None:
            try:
                os.waitpid(pid, os.WNOHANG)
            except ChildProcessError:
                pass

        # Flush and close the transcript last, so any final bytes the reader
        # loop wrote are on disk before we let go of the handle.
        if self._transcript_file is not None:
            try:
                self._transcript_file.flush()
                self._transcript_file.close()
            except OSError:
                pass
            self._transcript_file = None

    # ----- I/O --------------------------------------------------------------

    async def read_output(self) -> bytes | PtyExit:
        """Yield the next chunk of bytes from the child, or a PtyExit sentinel."""
        assert self._out_queue is not None, "PtySession not started"
        return await self._out_queue.get()

    def write_input(self, data: bytes) -> None:
        if self._master_fd is None:
            return
        try:
            os.write(self._master_fd, data)
        except OSError as exc:
            logger.debug("write_input failed: %s", exc)

    def resize(self, cols: int, rows: int) -> None:
        cols = max(1, int(cols))
        rows = max(1, int(rows))
        self._cols, self._rows = cols, rows
        self._apply_winsize(cols, rows)

    # ----- internal ---------------------------------------------------------

    def _apply_winsize(self, cols: int, rows: int) -> None:
        if self._master_fd is None:
            return
        try:
            # struct: rows, cols, xpixel, ypixel
            packed = struct.pack("HHHH", rows, cols, 0, 0)
            fcntl.ioctl(self._master_fd, termios.TIOCSWINSZ, packed)
        except OSError as exc:
            logger.debug("TIOCSWINSZ failed: %s", exc)

    def _reader_loop(self) -> None:
        fd = self._master_fd
        loop = self._loop
        queue = self._out_queue
        tfile = self._transcript_file
        ansi_carry = b""  # trailing partial escape held across read boundaries
        assert fd is not None and loop is not None and queue is not None

        returncode = 0
        try:
            while True:
                try:
                    data = os.read(fd, READ_CHUNK)
                except OSError as exc:
                    if exc.errno in (errno.EIO,):
                        # Child closed the slave end -> EOF on Linux.
                        break
                    if exc.errno == errno.EBADF:
                        break
                    raise
                if not data:
                    break
                # Tee an ANSI-stripped copy to the transcript. This is the one
                # chokepoint every output byte passes through, so nothing the
                # browser sees can escape the log. Flush per chunk so the file
                # survives a crash mid-session. Never let a log error kill the
                # shell: the transcript is a convenience, the terminal is not.
                if tfile is not None:
                    try:
                        chunk, ansi_carry = _strip_ansi_stream(ansi_carry, data)
                        tfile.write(chunk)
                        tfile.flush()
                    except OSError as exc:
                        logger.debug("transcript write failed: %s", exc)
                        tfile = None
                # Push bytes to the asyncio queue. put_nowait so a stalled
                # consumer can't block the reader thread; if the queue is
                # full we drop oldest. (xterm.js keeps up easily; the bound
                # exists only as a safety valve.)
                fut = asyncio.run_coroutine_threadsafe(queue.put(data), loop)
                try:
                    fut.result(timeout=2.0)
                except Exception:
                    logger.warning("PtySession output queue stalled; dropping chunk")
        finally:
            # Try to capture the child's exit code without blocking.
            if self._pid is not None:
                try:
                    _wpid, status = os.waitpid(self._pid, os.WNOHANG)
                    if _wpid == self._pid and os.WIFEXITED(status):
                        returncode = os.WEXITSTATUS(status)
                except ChildProcessError:
                    pass
            try:
                asyncio.run_coroutine_threadsafe(queue.put(PtyExit(returncode)), loop).result(timeout=1.0)
            except Exception:
                pass
