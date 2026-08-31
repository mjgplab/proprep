"""
Interactive Structure Viewer - HTTP Server

Simple HTTP server to serve the NGL viewer and structure files.
Uses Python's built-in http.server for minimal dependencies.

Author: ProPrep Development Team
Date: 2025-11-14
"""

import os
import getpass
import json
import logging
import socket
import sys
import threading
import webbrowser
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


# Text-mode browsers that Python's ``webbrowser`` registers as fallbacks
# whenever TERM is set, regardless of whether a display exists. Handing the
# NGL viewer page to one of these is never useful (the page is a WebGL app)
# and is actively harmful: ``GenericBrowser.open()`` blocks on the child, so
# a text browser either takes over the CLI's terminal or, as seen with
# lynx 2.8.9, segfaults on the page.
_CONSOLE_BROWSERS = {"www-browser", "links", "links2", "elinks", "lynx", "w3m"}


def gui_browser_unavailable_reason() -> Optional[str]:
    """Return None if a real browser can be opened here, else why it cannot.

    An explicitly set ``BROWSER`` environment variable always wins — that is
    the user telling us what to launch, including deliberately neutralising
    the launch with something like ``BROWSER=/bin/true``.
    """
    if os.environ.get("BROWSER"):
        return None

    if sys.platform.startswith("win") or sys.platform == "darwin":
        return None

    if not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")):
        return "no DISPLAY or WAYLAND_DISPLAY"

    try:
        browser = webbrowser.get()
    except webbrowser.Error:
        return "no browser registered"

    name = os.path.basename(getattr(browser, "name", "") or "")
    if name in _CONSOLE_BROWSERS:
        return f"only a text-mode browser ({name}) is available"

    return None


def ssh_forward_hint(port: int) -> Optional[str]:
    """An ``ssh -L`` command that tunnels ``port`` to the caller's machine.

    Returns None when this process is not on the far end of an SSH
    connection, in which case we have no idea how the user reaches this host
    and should not guess. SSH_CONNECTION is
    "<client ip> <client port> <server ip> <server port>", so the port the
    user actually dialled is field 4 (sshd may not be on 22).
    """
    conn = os.environ.get("SSH_CONNECTION", "").split()
    if len(conn) != 4:
        return None

    server_ip, server_port = conn[2], conn[3]
    try:
        user = getpass.getuser()
    except Exception:
        user = os.environ.get("USER", "")

    port_opt = "" if server_port == "22" else f"-p {server_port} "
    account = f"{user}@{server_ip}" if user else server_ip
    return f"ssh {port_opt}-L {port}:localhost:{port} {account}"


class ViewerHTTPRequestHandler(SimpleHTTPRequestHandler):
    """Custom HTTP request handler for the structure viewer."""

    # Class variables set by ViewerServer.
    # config_version monotonically increments each time the server's config
    # is replaced, so the browser can poll /version cheaply and only re-fetch
    # /config when something actually changed.
    config = {}
    config_version = 0
    template_path = ""
    structure_files = []
    # Scene save/load. ``scene_sink`` is a callable(payload) -> dict set by the
    # viewer object; the page POSTs its state to /scene and the sink writes the
    # file. ``scene_request`` rides on /version so the CLI can ask the page to
    # save without a push channel: the page sees a new token and POSTs.
    scene_sink = None
    scene_request = None

    def do_GET(self):
        """Handle GET requests."""
        if self.path == "/" or self.path == "/viewer":
            # Serve the HTML viewer template
            self.serve_html()
        elif self.path == "/config":
            # Serve the viewer configuration JSON
            self.serve_config()
        elif self.path == "/version":
            # Cheap poll endpoint: returns current config version
            self.serve_version()
        elif self.path.startswith("/structure/"):
            # Serve PDB structure files
            self.serve_structure()
        elif self.path == "/favicon.ico":
            # Silently ignore favicon requests (browsers always request this)
            self.send_response(204)  # No Content
            self.end_headers()
        else:
            # Default handler for other files
            super().do_GET()

    # Errors that mean "the client hung up", not "the server failed".
    _CLIENT_GONE = (BrokenPipeError, ConnectionResetError, ConnectionAbortedError)

    def _handle_serve_error(self, what: str, exc: Exception) -> None:
        """Report a failure while serving, distinguishing a vanished client.

        The browser polls /version every 1.5 s and never aborts, so any stall
        long enough to outlast its patience (a compiled extension such as
        MODELLER holding the GIL for minutes, with this single-threaded server
        unable to drain its backlog) leaves a queue of requests whose sockets
        are closed by the time they are answered. That is ordinary, and logging
        it at ERROR made a normal disconnect look like a fault.

        It also must not be answered: send_error would write a 500 to the same
        closed socket and raise again inside the handler.
        """
        if isinstance(exc, self._CLIENT_GONE):
            logger.debug("Client disconnected while serving %s: %s", what, exc)
            return
        logger.error(f"Error serving {what}: {exc}")
        try:
            self.send_error(500, f"Internal server error: {exc}")
        except self._CLIENT_GONE:
            logger.debug("Client also gone before the 500 for %s could be sent", what)

    def serve_html(self):
        """Serve the NGL viewer HTML template."""
        try:
            with open(self.template_path, 'r') as f:
                content = f.read()

            self.send_response(200)
            self.send_header('Content-type', 'text/html; charset=utf-8')
            self.send_header('Content-Length', len(content.encode('utf-8')))
            self.end_headers()
            self.wfile.write(content.encode('utf-8'))

        except FileNotFoundError:
            self.send_error(404, f"Template not found: {self.template_path}")
        except Exception as e:
            self._handle_serve_error("HTML", e)

    def serve_config(self):
        """Serve the viewer configuration as JSON."""
        try:
            payload = dict(self.config)
            payload["_config_version"] = self.config_version
            config_json = json.dumps(payload, indent=2)

            self.send_response(200)
            self.send_header('Content-type', 'application/json; charset=utf-8')
            self.send_header('Content-Length', len(config_json.encode('utf-8')))
            self.send_header('Cache-Control', 'no-store')
            self.end_headers()
            self.wfile.write(config_json.encode('utf-8'))

        except Exception as e:
            self._handle_serve_error("config", e)

    def do_POST(self):
        """POST /scene: the page hands over its representation and camera state."""
        if self.path != "/scene":
            self.send_error(404, "Unknown endpoint")
            return
        try:
            length = int(self.headers.get('Content-Length') or 0)
            raw = self.rfile.read(length) if length else b""
            payload = json.loads(raw.decode("utf-8") or "{}")
            sink = type(self).scene_sink
            if sink is None:
                result = {"ok": False, "error": "no scene handler registered"}
                status = 503
            else:
                result = sink(payload) or {"ok": True}
                status = 200 if result.get("ok", True) else 400
            body = json.dumps(result).encode("utf-8")
            self.send_response(status)
            self.send_header('Content-type', 'application/json; charset=utf-8')
            self.send_header('Content-Length', len(body))
            self.end_headers()
            self.wfile.write(body)
        except Exception as e:
            self._handle_serve_error("scene", e)

    def serve_version(self):
        """Tiny endpoint for the browser's poll loop."""
        try:
            info = {"version": self.config_version}
            if self.scene_request:
                info["scene_request"] = self.scene_request
            body = json.dumps(info).encode("utf-8")
            self.send_response(200)
            self.send_header('Content-type', 'application/json; charset=utf-8')
            self.send_header('Content-Length', len(body))
            self.send_header('Cache-Control', 'no-store')
            self.end_headers()
            self.wfile.write(body)
        except Exception as e:
            self._handle_serve_error("version", e)

    def serve_structure(self):
        """Serve a PDB structure file."""
        try:
            # Extract structure index from path: /structure/0 -> index 0
            path_parts = self.path.split('/')
            index = int(path_parts[2])

            if index < 0 or index >= len(self.structure_files):
                self.send_error(404, f"Structure index {index} not found")
                return

            structure_file = self.structure_files[index]

            if not os.path.exists(structure_file):
                self.send_error(404, f"Structure file not found: {structure_file}")
                return

            with open(structure_file, 'r') as f:
                content = f.read()

            self.send_response(200)
            self.send_header('Content-type', 'chemical/x-pdb; charset=utf-8')
            self.send_header('Content-Disposition', f'inline; filename="{os.path.basename(structure_file)}"')
            self.send_header('Content-Length', len(content.encode('utf-8')))
            # The URL /structure/<index> is stable, but the file behind an index
            # is rewritten in place (e.g. orientation overwrites *_oriented.pdb)
            # or rebound to a new path on a server relaunch. Without no-store the
            # browser serves the previously-cached body, so NGL re-renders the
            # OLD coordinates while no-store overlays (axis markers from /config)
            # update — the structure looks unchanged. Match /config and /version.
            self.send_header('Cache-Control', 'no-store')
            self.send_header('Pragma', 'no-cache')
            self.end_headers()
            self.wfile.write(content.encode('utf-8'))

        except ValueError:
            self.send_error(400, "Invalid structure index")
        except Exception as e:
            self._handle_serve_error("structure", e)

    def log_message(self, format, *args):
        """Override to suppress routine HTTP request logging."""
        # Suppress routine GET request logs to avoid cluttering console
        # Errors are still logged via send_error() which uses logger.error()
        pass

    def handle_one_request(self):
        """Override to catch BrokenPipeError from closed browser connections."""
        try:
            super().handle_one_request()
        except BrokenPipeError:
            # Browser closed connection before we could respond (e.g., favicon request)
            # This is normal behavior - just ignore it
            pass
        except ConnectionResetError:
            # Connection reset by browser - also normal
            pass


class ViewerServer:
    """
    HTTP server for the Interactive Structure Viewer.

    Serves the NGL viewer HTML template and provides endpoints for:
    - /viewer - Main viewer page
    - /config - Configuration JSON
    - /structure/{index} - PDB structure files
    """

    def __init__(self, config: Dict, structure_files: List[str], port: int = 8765,
                 scene_sink=None):
        """
        Initialize the viewer server.

        Args:
            config: Viewer configuration dict
            structure_files: List of PDB file paths
            port: Port to run server on (default: 8765)
            scene_sink: callable(payload) -> dict that persists a scene the
                page POSTs to /scene (see InteractiveStructureViewer)
        """
        self.config = config
        self.structure_files = structure_files
        self._scene_token = 0
        self.port = port
        self.server = None
        self.thread = None
        self.template_path = self._get_template_path()
        self._config_lock = threading.Lock()

        # Set by start(): whether a browser tab was actually opened, and if
        # not, the reason — callers print a port-forwarding hint instead.
        self.browser_opened = False
        self.headless_reason: Optional[str] = None

        # Set class variables for the request handler
        ViewerHTTPRequestHandler.config = self.config
        ViewerHTTPRequestHandler.config_version = 1
        ViewerHTTPRequestHandler.structure_files = self.structure_files
        ViewerHTTPRequestHandler.template_path = self.template_path
        ViewerHTTPRequestHandler.scene_sink = scene_sink
        ViewerHTTPRequestHandler.scene_request = None

    def request_scene(self, name: str) -> int:
        """Ask the open page to POST its current scene under ``name``.

        The page polls /version; a new token there makes it collect its state
        and POST /scene, which lands in ``scene_sink``. Returns the token so the
        caller can wait for that particular save."""
        self._scene_token += 1
        ViewerHTTPRequestHandler.scene_request = {"token": self._scene_token, "name": name}
        return self._scene_token

    def clear_scene_request(self) -> None:
        ViewerHTTPRequestHandler.scene_request = None

    def update_config(self, new_config: Dict) -> int:
        """Replace the served config and bump the version counter.

        The browser polls ``/version`` and re-fetches ``/config`` whenever the
        version changes, so this is how a long-running viewer is told to
        re-render with new annotations / focus / representations.

        Returns the new version number.
        """
        with self._config_lock:
            self.config = new_config
            ViewerHTTPRequestHandler.config = new_config
            ViewerHTTPRequestHandler.config_version += 1
            return ViewerHTTPRequestHandler.config_version

    def _get_template_path(self) -> str:
        """Get the path to the ngl_viewer.html template."""
        from proprep.utils.paths import get_package_dir
        template_path = get_package_dir() / "structure_prep" / "templates" / "ngl_viewer.html"

        if not template_path.exists():
            raise FileNotFoundError(f"Viewer template not found: {template_path}")

        return str(template_path)

    def find_available_port(self, start_port: int = 8765, max_attempts: int = 10) -> int:
        """
        Find an available port starting from start_port.

        Args:
            start_port: Port to start searching from
            max_attempts: Maximum number of ports to try

        Returns:
            Available port number

        Raises:
            RuntimeError: If no available port found
        """
        for port in range(start_port, start_port + max_attempts):
            try:
                sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                sock.bind(('localhost', port))
                sock.close()
                return port
            except OSError:
                continue

        raise RuntimeError(f"Could not find available port in range {start_port}-{start_port + max_attempts - 1}")

    def start(self, open_browser: bool = True) -> bool:
        """
        Start the HTTP server in a background thread.

        Args:
            open_browser: Whether to automatically open the browser. The
                ``PROPREP_WEB_SHELL`` env var overrides this to False — in
                that mode the parent web shell hosts the viewer in an
                iframe and we must not pop a separate tab.

        Returns:
            True if server started successfully
        """
        try:
            # Find available port if the default is in use
            try:
                self.server = HTTPServer(('localhost', self.port), ViewerHTTPRequestHandler)
            except OSError:
                logger.debug(f"Port {self.port} in use, finding alternative...")
                self.port = self.find_available_port(self.port)
                self.server = HTTPServer(('localhost', self.port), ViewerHTTPRequestHandler)

            # Start server in background thread
            self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
            self.thread.start()

            logger.debug(f"Viewer server started on http://localhost:{self.port}")

            in_web_shell = bool(os.environ.get("PROPREP_WEB_SHELL"))
            if open_browser and not in_web_shell:
                self.headless_reason = gui_browser_unavailable_reason()
                if self.headless_reason:
                    logger.debug(
                        "Not opening a browser (%s); viewer is served on port %d",
                        self.headless_reason,
                        self.port,
                    )
                else:
                    url = f"http://localhost:{self.port}/viewer"
                    webbrowser.open(url)
                    self.browser_opened = True
                    logger.debug(f"Opened browser at {url}")

            if in_web_shell:
                # Tell the parent uvicorn process which port we're on so it
                # can reverse-proxy /viewer/* to us and source the iframe.
                self._announce_to_web_shell()

            return True

        except Exception as e:
            logger.error(f"Failed to start viewer server: {e}")
            return False

    def _announce_to_web_shell(self) -> None:
        """Best-effort POST {port:N} to the parent web shell.

        Quietly tolerated on failure — the viewer still works via the
        direct URL even if the parent never hears about it. The shell does
        get an audible info log either way so a missing announcement is
        diagnosable.
        """
        shell_url = os.environ.get("PROPREP_WEB_SHELL_URL")
        if not shell_url:
            logger.info(
                "PROPREP_WEB_SHELL is set but PROPREP_WEB_SHELL_URL is empty; "
                "viewer running on http://localhost:%d/viewer (no parent to notify)",
                self.port,
            )
            return

        # Wait briefly for our own listen socket to accept before we tell
        # anyone about it. ``HTTPServer.__init__`` binds the socket but the
        # ``serve_forever`` thread may not have entered ``accept()`` yet —
        # without this guard the parent can hit a tiny window of 502s on
        # the first /version poll right after a relaunch.
        if not self._wait_until_accepting(timeout=1.0):
            logger.warning("Viewer port %d not accepting after 1.0s; announcing anyway", self.port)

        endpoint = shell_url.rstrip("/") + "/internal/viewer-announce"
        body = json.dumps({"port": self.port}).encode("utf-8")
        try:
            import urllib.request
            req = urllib.request.Request(
                endpoint,
                data=body,
                headers={"Content-Type": "application/json"},
                method="POST",
            )
            urllib.request.urlopen(req, timeout=2.0).read()
            logger.debug("Announced viewer (port %d) to web shell at %s", self.port, shell_url)
        except Exception as exc:
            logger.warning(
                "Could not announce viewer port %d to web shell at %s: %s",
                self.port, shell_url, exc,
            )

    def _wait_until_accepting(self, timeout: float = 1.0) -> bool:
        """Probe localhost:self.port until a TCP connect succeeds, or timeout."""
        import time
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            try:
                with socket.create_connection(("127.0.0.1", self.port), timeout=0.2):
                    return True
            except OSError:
                time.sleep(0.02)
        return False

    def stop(self):
        """Stop the HTTP server."""
        if self.server:
            logger.debug("Stopping viewer server")
            self.server.shutdown()
            self.server.server_close()
            self.server = None
            logger.debug("Viewer server stopped")

    def is_running(self) -> bool:
        """Check if server is running."""
        return self.server is not None and self.thread is not None and self.thread.is_alive()

    def get_url(self) -> str:
        """Get the viewer URL."""
        return f"http://localhost:{self.port}/viewer"
