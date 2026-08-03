"""FastAPI app for the ProPrep web shell.

Routes:
    GET  /            -> static index.html (xterm.js shell)
    GET  /healthz     -> liveness probe
    GET  /static/*    -> static assets
    WS   /ws/term     -> bidirectional terminal channel

The terminal websocket protocol is intentionally tiny:

    Client -> Server:
        - text frame "{json}" interpreted as a control message
            {"type": "resize", "cols": int, "rows": int}
        - binary frame: raw stdin bytes for the PTY child

    Server -> Client:
        - binary frame: raw stdout/stderr bytes from the PTY child
        - text frame "{json}" for control messages
            {"type": "exit", "returncode": int}

Phase A is single-user, single-connection. We accept one terminal websocket
at a time; if a second one connects we close the first.
"""

from __future__ import annotations

import asyncio
import json
import logging
import os
import signal
from pathlib import Path
from typing import Any

import httpx
from fastapi import FastAPI, Request, WebSocket, WebSocketDisconnect
from fastapi.responses import FileResponse, PlainTextResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles

from proprep.web.pty_session import PtyExit, PtySession

logger = logging.getLogger(__name__)

_STATIC_DIR = Path(__file__).resolve().parent / "static"

# Auto-shutdown when the browser is gone. The grace window is long enough
# that a page reload (which closes and reopens the websocket) doesn't trip
# it, but short enough that closing the tab feels instant.
_SHUTDOWN_GRACE_SECONDS = 2.0

# Module-globals are sufficient here — proprep.web is single-user, single-
# process. ``_active_terms`` counts open terminal websockets;
# ``_shutdown_task`` is the pending grace-period timer (cancelled when a
# new websocket connects in time).
_active_terms = 0
_shutdown_task: asyncio.Task | None = None
_first_connect_seen = False

# Reverse-proxy state. The PTY child's ViewerServer POSTs its port to
# ``/internal/viewer-announce`` once it's up; we stash it here and the
# /viewer, /config, /version, /structure/{path} routes proxy to that port.
# ``_control_clients`` is the broadcast set for the page-side control
# websocket (used to push viewer-up announcements to the iframe).
_viewer_port: int | None = None
_control_clients: set[WebSocket] = set()
_proxy_client: httpx.AsyncClient | None = None


app = FastAPI(title="ProPrep Web Shell", docs_url=None, redoc_url=None)
app.mount("/static", StaticFiles(directory=str(_STATIC_DIR)), name="static")


@app.on_event("startup")
async def _startup() -> None:
    global _proxy_client
    # Single shared client; HTTP/1.1 keep-alive against localhost is fine.
    _proxy_client = httpx.AsyncClient(timeout=httpx.Timeout(10.0, connect=2.0))


@app.on_event("shutdown")
async def _shutdown() -> None:
    global _proxy_client
    if _proxy_client is not None:
        await _proxy_client.aclose()
        _proxy_client = None


@app.get("/healthz", response_class=PlainTextResponse)
async def healthz() -> str:
    return "ok"


@app.get("/")
async def index() -> FileResponse:
    return FileResponse(_STATIC_DIR / "index.html", media_type="text/html")


@app.get("/shell-theme")
async def shell_theme() -> dict:
    """Launch-time UI theme for the terminal pane (set by proprep-web args).

    Distinct from the proxied ``/config`` route (which belongs to the viewer).
    The browser (web/static/app.js) fetches this once and applies the xterm
    theme + contrast. Read live from the environment the launcher populated.
    """
    return {
        "theme": os.environ.get("PROPREP_WEB_THEME", "dark"),
        "highContrast": os.environ.get("PROPREP_WEB_HIGH_CONTRAST", "0") == "1",
    }


async def _shutdown_after_grace() -> None:
    """Sleep ``_SHUTDOWN_GRACE_SECONDS``; if still no terminals, exit."""
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
    # Cancel any pending auto-shutdown — a (possibly reloaded) browser is back.
    if _shutdown_task is not None and not _shutdown_task.done():
        _shutdown_task.cancel()
        _shutdown_task = None
    _active_terms += 1
    _first_connect_seen = True

    await ws.accept()
    # The PTY child needs an absolute URL to POST viewer announcements to.
    # The Host header is the most reliable source — it's literally what the
    # browser thinks our address is. Falls back to the parsed ws.url if a
    # client somehow doesn't send Host (shouldn't happen in practice).
    host_header = ws.headers.get("host")
    if host_header:
        shell_origin = f"http://{host_header}"
    else:
        shell_origin = f"http://{ws.url.hostname or '127.0.0.1'}:{ws.url.port or 8000}"
    logger.info("Terminal websocket connected; shell_origin=%s", shell_origin)
    # Session transcript is on by default; PROPREP_WEB_NO_TRANSCRIPT=1 (set by
    # proprep-web --no-transcript) disables it.
    want_transcript = os.environ.get("PROPREP_WEB_NO_TRANSCRIPT", "0") != "1"
    session = PtySession(shell_url=shell_origin, transcript=want_transcript)
    try:
        session.start()
    except Exception as exc:
        logger.exception("Failed to start PTY session: %s", exc)
        await ws.send_text(json.dumps({"type": "error", "message": f"failed to start: {exc}"}))
        await ws.close()
        _active_terms -= 1
        return
    if session.transcript_path is not None:
        logger.info("Session transcript: %s", session.transcript_path)

    async def pump_pty_to_ws() -> None:
        while True:
            chunk = await session.read_output()
            if isinstance(chunk, PtyExit):
                try:
                    await ws.send_text(json.dumps({"type": "exit", "returncode": chunk.returncode}))
                except Exception:
                    pass
                return
            try:
                await ws.send_bytes(chunk)
            except Exception:
                return

    async def pump_ws_to_pty() -> None:
        while True:
            try:
                msg = await ws.receive()
            except WebSocketDisconnect:
                return
            if msg.get("type") == "websocket.disconnect":
                return
            # Binary stdin
            data = msg.get("bytes")
            if data:
                session.write_input(data)
                continue
            # Text control frame
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
                session.resize(cols, rows)
            else:
                logger.debug("Unknown control type: %s", kind)

    pump_out = asyncio.create_task(pump_pty_to_ws())
    pump_in = asyncio.create_task(pump_ws_to_pty())
    try:
        done, pending = await asyncio.wait({pump_out, pump_in}, return_when=asyncio.FIRST_COMPLETED)
        for task in pending:
            task.cancel()
        for task in pending:
            try:
                await task
            except (asyncio.CancelledError, Exception):
                pass
    finally:
        await session.close()
        try:
            await ws.close()
        except Exception:
            pass
        _active_terms -= 1
        if _active_terms == 0 and _first_connect_seen:
            # Arm the grace timer. A fresh websocket within the window will
            # cancel it.
            _shutdown_task = asyncio.create_task(_shutdown_after_grace())


# ---------------------------------------------------------------------------
# Viewer integration: announce + reverse proxy + control channel
# ---------------------------------------------------------------------------

# Path prefixes the NGL viewer page expects at its origin. Because the
# iframe shares the parent's origin (we proxy here), the page's relative
# fetch('/config') etc. land on these routes and we forward them to the
# announced viewer port. ``_VIEWER_PROXY_BASE_PATHS`` is documentation; the
# routes below define the actual handlers.
_VIEWER_PROXY_BASE_PATHS = ("/viewer", "/config", "/version", "/structure")


async def _broadcast_control(message: dict[str, Any]) -> None:
    """Send a JSON message to every connected control client.

    Failures eject the client from the set so a stuck socket can't wedge
    future broadcasts.
    """
    payload = json.dumps(message)
    dead: list[WebSocket] = []
    for ws in list(_control_clients):
        try:
            await ws.send_text(payload)
        except Exception:
            dead.append(ws)
    for ws in dead:
        _control_clients.discard(ws)


@app.post("/internal/viewer-announce")
async def viewer_announce(req: Request) -> Response:
    """Internal endpoint: the PTY child's ViewerServer calls this when it
    starts, so we can proxy to it and tell the page to source the iframe.

    Body: ``{"port": <int>}``. Loopback-only by design (we bind 127.0.0.1).
    """
    global _viewer_port
    try:
        body = await req.json()
        port = int(body["port"])
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return Response("bad request", status_code=400)
    if not (1 <= port <= 65535):
        return Response("bad port", status_code=400)
    _viewer_port = port
    logger.info("Viewer server announced on port %d", port)
    await _broadcast_control({"type": "viewer_announce", "port": port})
    return Response(status_code=204)


@app.websocket("/ws/control")
async def ws_control(ws: WebSocket) -> None:
    """Page-side control channel.

    Used to push events to the parent page (currently just
    ``viewer_announce`` so the page knows when to source the iframe).
    Phase B will reuse this channel for structured-prompt events.
    """
    await ws.accept()
    _control_clients.add(ws)
    # Replay any cached state so a freshly-loaded page catches up.
    if _viewer_port is not None:
        try:
            await ws.send_text(json.dumps({"type": "viewer_announce", "port": _viewer_port}))
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
        _control_clients.discard(ws)
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
    """Proxy a single GET to the announced viewer port.

    Streams the upstream body to keep large PDB structures memory-cheap.
    Headers we don't proxy: hop-by-hop (Connection, Transfer-Encoding,
    Keep-Alive), and Content-Length (httpx may decompress and a stale
    length would mismatch).

    ``forward_query=False`` is used for the ``/viewer`` HTML page because
    the upstream HTTP handler does an exact-string path match — a cache-
    buster query like ``?_v=8765&t=...`` (used in the iframe src so a
    relaunch forces the browser to refetch) would make ``self.path`` no
    longer equal ``"/viewer"`` and the handler would 404.
    """
    if _viewer_port is None or _proxy_client is None:
        logger.info("Proxy %s -> 503 (no viewer announced yet)", path)
        return Response(
            "viewer not running — launch the structure viewer from the CLI",
            status_code=503,
            media_type="text/plain",
        )
    url = f"http://127.0.0.1:{_viewer_port}{path}"
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
