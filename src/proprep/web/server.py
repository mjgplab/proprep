"""FastAPI app for the ProPrep web shell.

Routes:
    GET  /                 -> static index.html (xterm.js shell); ``?seat=<token>``
                              selects a seat in hosted mode and sets the seat cookie
    GET  /healthz          -> liveness probe
    GET  /seat-info        -> {name, hosted, alive} for the topbar
    GET  /shell-theme      -> launch-time UI theme
    GET  /download         -> zip of the seat's working directory
    GET  /static/*         -> static assets
    WS   /ws/term          -> bidirectional terminal channel (attaches to the seat)
    WS   /ws/control       -> page-side event channel (viewer announcements)
    POST /internal/viewer-announce
    POST /seat/{token}/internal/viewer-announce
                           -> the PTY child's viewer server reports its port
    GET  /viewer, /config, /version, /structure/{path}
                           -> reverse proxy to the seat's announced viewer port

The terminal websocket protocol is intentionally tiny:

    Client -> Server:
        - text frame "{json}" interpreted as a control message
            {"type": "resize", "cols": int, "rows": int}
        - binary frame: raw stdin bytes for the PTY child

    Server -> Client:
        - text frame {"type": "attached", "seat": str, "replay": bool,
                      "alive": bool, "returncode": int|null}
          sent first, followed by the scrollback replay as one binary frame
        - binary frame: raw stdout/stderr bytes from the PTY child
        - text frame {"type": "exit", "returncode": int}
        - text frame {"type": "displaced"}  (a newer connection took the seat)

Sessions belong to seats (see ``seats.py``), not to websockets: a dropped
connection leaves the child running and the next connection re-attaches
with a replay. One terminal websocket per seat at a time; a newer one
displaces the older.

Seat resolution: hosted mode requires a token, taken from ``?seat=`` on the
page load (which sets the ``proprep_seat`` cookie) or from that cookie on
every later request, so the viewer iframe's own fetches and the websocket
upgrades are seat-scoped without the page knowing about tokens. Local mode
has one seat and ignores tokens.

Configuration arrives through the environment, populated by ``__main__``:
    PROPREP_WEB_SEATS_FILE   path to seats.json -> hosted mode
    PROPREP_WEB_BIND         http://127.0.0.1:<port>, for the child's announce URL
    PROPREP_WEB_NO_TRANSCRIPT, PROPREP_WEB_THEME, ... (see __main__)
"""

from __future__ import annotations

import asyncio
import datetime
import json
import logging
import os
import signal
import tempfile
import zipfile
from pathlib import Path
from typing import Any, Optional

import httpx
from fastapi import FastAPI, Request, WebSocket, WebSocketDisconnect
from fastapi.responses import FileResponse, PlainTextResponse, RedirectResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from starlette.background import BackgroundTask

from proprep.web.seats import Seat, SeatInfo, SeatRegistry, SessionFactory, load_or_create_seats
from proprep.web.pty_session import PtySession

logger = logging.getLogger(__name__)

_STATIC_DIR = Path(__file__).resolve().parent / "static"
SEAT_COOKIE = "proprep_seat"

# Local mode only: auto-shutdown when the browser is gone. The grace window
# is long enough that a page reload (which closes and reopens the websocket)
# doesn't trip it, but short enough that closing the tab feels instant.
_SHUTDOWN_GRACE_SECONDS = 2.0

_registry: Optional[SeatRegistry] = None
_active_terms = 0
_shutdown_task: Optional[asyncio.Task] = None
_first_connect_seen = False
_proxy_client: Optional[httpx.AsyncClient] = None


def _seats_file() -> Optional[Path]:
    raw = os.environ.get("PROPREP_WEB_SEATS_FILE", "").strip()
    return Path(raw) if raw else None


def _hosted() -> bool:
    return _seats_file() is not None


def configure_registry(
    *,
    seats: Optional[list[SeatInfo]] = None,
    local: Optional[bool] = None,
    session_factory: SessionFactory = PtySession,
) -> SeatRegistry:
    """Build (or rebuild) the process-wide registry.

    Called at app startup from the environment; tests call it directly with
    a fake session factory.
    """
    global _registry, _active_terms, _shutdown_task, _first_connect_seen
    transcript = os.environ.get("PROPREP_WEB_NO_TRANSCRIPT", "0") != "1"
    if seats is None:
        seats_file = _seats_file()
        if seats_file is not None:
            data = json.loads(seats_file.read_text(encoding="utf-8"))
            seats = [SeatInfo(str(e["name"]), str(e["token"]), str(e["cwd"])) for e in data["seats"]]
            local = False
        else:
            seats = [SeatInfo("local", None, os.getcwd())]
            local = True
    if local is None:
        local = all(s.token is None for s in seats)
    _registry = SeatRegistry(seats, local=local, session_factory=session_factory, transcript=transcript)
    _active_terms = 0
    _shutdown_task = None
    _first_connect_seen = False
    return _registry


def registry() -> SeatRegistry:
    if _registry is None:
        return configure_registry()
    return _registry


app = FastAPI(title="ProPrep Web Shell", docs_url=None, redoc_url=None)
app.mount("/static", StaticFiles(directory=str(_STATIC_DIR)), name="static")


@app.on_event("startup")
async def _startup() -> None:
    global _proxy_client
    # Single shared client; HTTP/1.1 keep-alive against localhost is fine.
    _proxy_client = httpx.AsyncClient(timeout=httpx.Timeout(10.0, connect=2.0))
    registry()


@app.on_event("shutdown")
async def _shutdown() -> None:
    global _proxy_client
    if _proxy_client is not None:
        await _proxy_client.aclose()
        _proxy_client = None
    if _registry is not None:
        await _registry.close_all()


# ---------------------------------------------------------------------------
# Seat resolution
# ---------------------------------------------------------------------------

def _token_from(scope_obj: Any) -> Optional[str]:
    """``?seat=`` wins over the cookie so a new link can replace a stale one."""
    token = scope_obj.query_params.get("seat")
    if token:
        return token
    return scope_obj.cookies.get(SEAT_COOKIE)


def _resolve_seat(scope_obj: Any) -> Optional[Seat]:
    return registry().lookup(_token_from(scope_obj))


def _forbidden() -> Response:
    return Response(
        "This ProPrep link is not valid. Ask the workshop host for your seat link.",
        status_code=403,
        media_type="text/plain",
    )


def _bind_base() -> Optional[str]:
    raw = os.environ.get("PROPREP_WEB_BIND", "").strip()
    return raw.rstrip("/") or None


def _shell_url_for(seat: Seat, ws: WebSocket) -> str:
    """Absolute URL the PTY child POSTs viewer announcements to.

    Preferred: the loopback bind address the launcher recorded, with a
    per-seat path so the announcement identifies its seat. This never
    leaves the machine, so proxies, TLS and auth in front of the shell are
    irrelevant to it. Fallback (someone ran uvicorn by hand): the Host
    header, as in Phase A.
    """
    base = _bind_base()
    if base is None:
        host_header = ws.headers.get("host")
        if host_header:
            base = f"http://{host_header}"
        else:
            base = f"http://{ws.url.hostname or '127.0.0.1'}:{ws.url.port or 8000}"
    if seat.token:
        return f"{base}/seat/{seat.token}"
    return base


# ---------------------------------------------------------------------------
# Plain routes
# ---------------------------------------------------------------------------

@app.get("/healthz", response_class=PlainTextResponse)
async def healthz() -> str:
    return "ok"


@app.get("/")
async def index(request: Request) -> Response:
    seat = _resolve_seat(request)
    if seat is None:
        return _forbidden()
    token = request.query_params.get("seat")
    if token and not registry().local:
        # Strip the token from the address bar and remember it in a cookie so
        # the iframe's fetches and the websockets are seat-scoped.
        resp: Response = RedirectResponse("/", status_code=303)
        resp.set_cookie(SEAT_COOKIE, token, httponly=True, samesite="lax", path="/")
        return resp
    return FileResponse(_STATIC_DIR / "index.html", media_type="text/html")


@app.get("/seat-info")
async def seat_info(request: Request) -> Response:
    seat = _resolve_seat(request)
    if seat is None:
        return _forbidden()
    return Response(
        json.dumps({
            "name": seat.name,
            "hosted": not registry().local,
            "alive": seat.alive,
            "cwd": seat.cwd,
        }),
        media_type="application/json",
    )


@app.get("/shell-theme")
async def shell_theme() -> dict:
    """Launch-time UI theme for the terminal pane (set by proprep-web args).

    Distinct from the proxied ``/config`` route (which belongs to the viewer).
    ``fontSize`` is None unless ``--font-size`` was passed; the browser then
    keeps its own remembered size instead of being overridden on every load.
    """
    try:
        font_size = int(os.environ.get("PROPREP_WEB_FONT_SIZE", "") or 0) or None
    except ValueError:
        font_size = None
    return {
        "theme": os.environ.get("PROPREP_WEB_THEME", "dark"),
        "highContrast": os.environ.get("PROPREP_WEB_HIGH_CONTRAST", "0") == "1",
        "fontSize": font_size,
    }


def _zip_directory(root: Path, dest: Path) -> int:
    """Zip ``root`` into ``dest``; returns the number of files written."""
    count = 0
    with zipfile.ZipFile(dest, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(root.rglob("*")):
            if path.is_symlink() or not path.is_file():
                continue
            zf.write(path, path.relative_to(root).as_posix())
            count += 1
    return count


@app.get("/download")
async def download(request: Request) -> Response:
    """The seat's working directory as a zip, for attendees to take home."""
    seat = _resolve_seat(request)
    if seat is None:
        return _forbidden()
    root = Path(seat.cwd)
    if not root.is_dir():
        return Response("project directory not found", status_code=404, media_type="text/plain")
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"proprep_{seat.name}_{stamp}.zip"
    tmp = tempfile.NamedTemporaryFile(prefix="proprep_dl_", suffix=".zip", delete=False)
    tmp.close()
    dest = Path(tmp.name)
    try:
        await asyncio.get_running_loop().run_in_executor(None, _zip_directory, root, dest)
    except Exception as exc:
        dest.unlink(missing_ok=True)
        logger.exception("download: zipping %s failed", root)
        return Response(f"could not package project: {exc}", status_code=500, media_type="text/plain")
    return FileResponse(
        dest,
        media_type="application/zip",
        filename=filename,
        background=BackgroundTask(lambda: dest.unlink(missing_ok=True)),
    )


# ---------------------------------------------------------------------------
# Terminal websocket
# ---------------------------------------------------------------------------

async def _shutdown_after_grace() -> None:
    """Local mode: sleep the grace period; if still no terminals, exit."""
    try:
        await asyncio.sleep(_SHUTDOWN_GRACE_SECONDS)
    except asyncio.CancelledError:
        return
    if _active_terms == 0:
        logger.info("No browser reconnected within %.1fs; shutting down.", _SHUTDOWN_GRACE_SECONDS)
        # SIGINT is what uvicorn already handles for clean shutdown
        # (sets server.should_exit, drains, then exits).
        os.kill(os.getpid(), signal.SIGINT)


@app.websocket("/ws/term")
async def ws_term(ws: WebSocket) -> None:
    global _active_terms, _shutdown_task, _first_connect_seen
    seat = _resolve_seat(ws)
    if seat is None:
        await ws.close(code=4003, reason="invalid seat")
        return
    if _shutdown_task is not None and not _shutdown_task.done():
        _shutdown_task.cancel()
        _shutdown_task = None
    _active_terms += 1
    _first_connect_seen = True

    await ws.accept()
    try:
        started = seat.ensure_session(_shell_url_for(seat, ws))
    except Exception as exc:
        logger.exception("Seat %s: failed to start PTY session: %s", seat.name, exc)
        await ws.send_text(json.dumps({"type": "error", "message": f"failed to start: {exc}"}))
        await ws.close()
        _active_terms -= 1
        return
    if started and seat.session is not None and getattr(seat.session, "transcript_path", None):
        logger.info("Seat %s: session transcript: %s", seat.name, seat.session.transcript_path)
    logger.info("Seat %s: terminal websocket %s", seat.name, "started" if started else "re-attached")

    try:
        await seat.attach(ws)
        while True:
            try:
                msg = await ws.receive()
            except WebSocketDisconnect:
                break
            if msg.get("type") == "websocket.disconnect":
                break
            data = msg.get("bytes")
            if data:
                seat.write_input(data)
                continue
            text = msg.get("text")
            if not text:
                continue
            try:
                control = json.loads(text)
            except json.JSONDecodeError:
                logger.debug("Ignoring non-JSON text frame: %r", text[:64])
                continue
            kind = control.get("type")
            if kind == "resize":
                cols = int(control.get("cols") or 0) or 80
                rows = int(control.get("rows") or 0) or 24
                seat.resize(cols, rows)
            else:
                logger.debug("Unknown control type: %s", kind)
    finally:
        seat.detach(ws)
        try:
            await ws.close()
        except Exception:
            pass
        _active_terms -= 1
        if registry().local and _active_terms == 0 and _first_connect_seen:
            # Arm the grace timer. A fresh websocket within the window will
            # cancel it. Hosted seats never shut the server down.
            _shutdown_task = asyncio.create_task(_shutdown_after_grace())


# ---------------------------------------------------------------------------
# Viewer integration: announce + reverse proxy + control channel
# ---------------------------------------------------------------------------

# Path prefixes the NGL viewer page expects at its origin. Because the
# iframe shares the parent's origin (we proxy here), the page's relative
# fetch('/config') etc. land on these routes and we forward them to the
# seat's announced viewer port.
_VIEWER_PROXY_BASE_PATHS = ("/viewer", "/config", "/version", "/structure", "/scene")

_LOOPBACK = {"127.0.0.1", "::1", "localhost"}


async def _broadcast_control(seat: Seat, message: dict[str, Any]) -> None:
    """Send a JSON message to every control client of ``seat``.

    Failures eject the client from the set so a stuck socket can't wedge
    future broadcasts.
    """
    payload = json.dumps(message)
    dead: list[WebSocket] = []
    for client in list(seat.control_clients):
        try:
            await client.send_text(payload)
        except Exception:
            dead.append(client)
    for client in dead:
        seat.control_clients.discard(client)


async def _handle_announce(seat: Seat, req: Request) -> Response:
    """The PTY child's ViewerServer calls this when it starts, so we can
    proxy to it and tell the page to source the iframe.

    Body: ``{"port": <int>}``. Loopback-only: the child runs on this host
    and announces to the bind address, never through a proxy.
    """
    client_host = req.client.host if req.client else None
    if client_host not in _LOOPBACK:
        return Response("forbidden", status_code=403)
    try:
        body = await req.json()
        port = int(body["port"])
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return Response("bad request", status_code=400)
    if not (1 <= port <= 65535):
        return Response("bad port", status_code=400)
    seat.viewer_port = port
    logger.info("Seat %s: viewer server announced on port %d", seat.name, port)
    await _broadcast_control(seat, {"type": "viewer_announce", "port": port})
    return Response(status_code=204)


@app.post("/internal/viewer-announce")
async def viewer_announce(req: Request) -> Response:
    """Local mode (and Phase A children without a seat path)."""
    reg = registry()
    seat = reg.default if reg.local else _resolve_seat(req)
    if seat is None:
        return Response("unknown seat", status_code=404)
    return await _handle_announce(seat, req)


@app.post("/seat/{token}/internal/viewer-announce")
async def viewer_announce_seat(token: str, req: Request) -> Response:
    seat = registry().lookup(token)
    if seat is None:
        return Response("unknown seat", status_code=404)
    return await _handle_announce(seat, req)


@app.websocket("/ws/control")
async def ws_control(ws: WebSocket) -> None:
    """Page-side control channel: pushes ``viewer_announce`` to the page."""
    seat = _resolve_seat(ws)
    if seat is None:
        await ws.close(code=4003, reason="invalid seat")
        return
    await ws.accept()
    seat.control_clients.add(ws)
    # Replay any cached state so a freshly-loaded page catches up.
    if seat.viewer_port is not None:
        try:
            await ws.send_text(json.dumps({"type": "viewer_announce", "port": seat.viewer_port}))
        except Exception:
            pass
    try:
        while True:
            try:
                msg = await ws.receive()
            except WebSocketDisconnect:
                break
            if msg.get("type") == "websocket.disconnect":
                break
            # Page->server messages aren't used yet; ignore.
    finally:
        seat.control_clients.discard(ws)
        try:
            await ws.close()
        except Exception:
            pass


async def _proxy_to_viewer(
    path: str,
    request: Request,
    *,
    forward_query: bool = True,
) -> Response:
    """Proxy a single GET to the seat's announced viewer port.

    Streams the upstream body to keep large PDB structures memory-cheap.
    Headers we don't proxy: hop-by-hop (Connection, Transfer-Encoding,
    Keep-Alive), and Content-Length (httpx may decompress and a stale
    length would mismatch).

    ``forward_query=False`` is used for the ``/viewer`` HTML page because
    the upstream HTTP handler does an exact-string path match; a cache-
    buster query like ``?_v=8765&t=...`` would make it 404.
    """
    seat = _resolve_seat(request)
    if seat is None:
        return _forbidden()
    if seat.viewer_port is None or _proxy_client is None:
        logger.info("Seat %s: proxy %s -> 503 (no viewer announced yet)", seat.name, path)
        return Response(
            "viewer not running — launch the structure viewer from the CLI",
            status_code=503,
            media_type="text/plain",
        )
    url = f"http://127.0.0.1:{seat.viewer_port}{path}"
    if forward_query and request.url.query:
        url = f"{url}?{request.url.query}"
    try:
        upstream = await _proxy_client.send(
            _proxy_client.build_request("GET", url),
            stream=True,
        )
    except httpx.RequestError as exc:
        logger.warning("Proxy upstream error for %s: %s", url, exc)
        return Response(
            f"viewer unreachable at {url}: {exc}",
            status_code=502,
            media_type="text/plain",
        )

    skip = {"connection", "keep-alive", "transfer-encoding", "content-length", "content-encoding"}
    headers = {k: v for k, v in upstream.headers.items() if k.lower() not in skip}

    async def body_iter():
        try:
            async for chunk in upstream.aiter_raw():
                yield chunk
        finally:
            await upstream.aclose()

    return StreamingResponse(
        body_iter(),
        status_code=upstream.status_code,
        headers=headers,
        media_type=upstream.headers.get("content-type"),
    )


@app.get("/viewer")
async def proxy_viewer_root(request: Request) -> Response:
    # Drop query — see _proxy_to_viewer's docstring.
    return await _proxy_to_viewer("/viewer", request, forward_query=False)


@app.get("/config")
async def proxy_viewer_config(request: Request) -> Response:
    return await _proxy_to_viewer("/config", request)


@app.get("/version")
async def proxy_viewer_version(request: Request) -> Response:
    return await _proxy_to_viewer("/version", request)


@app.get("/structure/{rest:path}")
async def proxy_viewer_structure(rest: str, request: Request) -> Response:
    return await _proxy_to_viewer(f"/structure/{rest}", request)


@app.post("/scene")
async def proxy_viewer_scene(request: Request) -> Response:
    """Forward the page's scene POST to the seat's viewer (the GET proxy
    above streams; this is a small JSON body, sent whole)."""
    seat = _resolve_seat(request)
    if seat is None:
        return _forbidden()
    if seat.viewer_port is None or _proxy_client is None:
        return Response("viewer not running", status_code=503, media_type="text/plain")
    url = f"http://127.0.0.1:{seat.viewer_port}/scene"
    body = await request.body()
    try:
        upstream = await _proxy_client.post(
            url, content=body, headers={"content-type": "application/json"})
    except httpx.RequestError as exc:
        return Response(f"viewer unreachable at {url}: {exc}", status_code=502, media_type="text/plain")
    return Response(upstream.content, status_code=upstream.status_code,
                    media_type=upstream.headers.get("content-type", "application/json"))
