"""Entry point: ``python -m proprep.web``.

Boots the FastAPI app on localhost (single-user, local launch) and opens
the default browser. POSIX-only (uses pty.fork in pty_session).

Multiple concurrent instances are supported: by default the launcher
scans upward from the requested port for the first free one, so the
second instance lands on 8001, the third on 8002, and so on. Pass
``--strict-port`` to disable the scan and fail if the requested port
is already taken.
"""

import argparse
import logging
import os
import socket
import sys
import threading
import time
import webbrowser
from typing import Optional


PORT_SCAN_LIMIT = 20  # how many ports to try above the requested start


def _find_free_port(host: str, start_port: int,
                     max_attempts: int = PORT_SCAN_LIMIT) -> int:
    """Return the first free port in [start_port, start_port + max_attempts).

    Raises OSError if none are free. There is a tiny TOCTOU window between
    the probe close() and uvicorn's bind() — acceptable for a local-dev
    launcher (and the same race exists with port=0).
    """
    last_err: Optional[OSError] = None
    for offset in range(max_attempts):
        port = start_port + offset
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind((host, port))
                return port
            except OSError as exc:
                last_err = exc
                continue
    raise OSError(
        f"No free port in {start_port}..{start_port + max_attempts - 1} "
        f"on {host}: last error {last_err}")


def main() -> int:
    parser = argparse.ArgumentParser(prog="proprep.web", description="ProPrep Web Shell")
    parser.add_argument("--host", default="127.0.0.1", help="Bind host (default: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=8000,
                          help="Starting port for the auto-scan (default: 8000). "
                               "If busy, the launcher tries 8001, 8002, ...")
    parser.add_argument("--strict-port", action="store_true",
                          help="Fail instead of scanning if the requested port is busy")
    parser.add_argument("--no-browser", action="store_true", help="Do not auto-open the browser")
    _THEME_CHOICES = ["dark", "light", "black", "white"]
    parser.add_argument("theme_pos", nargs="?", choices=_THEME_CHOICES, default=None,
                          metavar="THEME",
                          help="Terminal pane color scheme as a positional word, "
                               "e.g. 'proprep-web white'. Same choices as --theme.")
    parser.add_argument("--theme", choices=_THEME_CHOICES, default=None,
                          help="Terminal pane color scheme (default: dark). "
                               "'black'/'white' are aliases for dark/light. "
                               "The structure viewer pane always stays dark.")
    parser.add_argument("--font-size", type=int, default=None, metavar="PX",
                          help="Terminal pane font size in pixels (default: 13, "
                               "range 8-32). Useful for screenshots and figures. "
                               "Also adjustable live with Ctrl/Cmd +/- in the page.")
    parser.add_argument("--high-contrast", action="store_true",
                          help="Boost terminal text contrast for accessibility "
                               "(xterm minimumContrastRatio=7, WCAG AAA).")
    parser.add_argument("--no-transcript", action="store_true",
                          help="Do not write a plain-text session log. By default "
                               "each session is teed to "
                               "~/.proprep/web_sessions/session_<timestamp>.log.")
    parser.add_argument("--log-level", default="info", help="uvicorn log level")
    args = parser.parse_args()

    # Accept both 'proprep-web white' (positional) and '--theme white'; the
    # explicit flag wins if both are given, else the positional, else dark.
    # 'black'/'white' are user-friendly aliases for the canonical dark/light.
    chosen = args.theme or args.theme_pos or "dark"
    theme = {"black": "dark", "white": "light"}.get(chosen, chosen)
    # Hand the launch-time UI choices to the server process (same process, so
    # env is the simplest channel) — it serves them at GET /shell-theme and the
    # browser applies them. See web/static/app.js.
    os.environ["PROPREP_WEB_THEME"] = theme
    os.environ["PROPREP_WEB_HIGH_CONTRAST"] = "1" if args.high_contrast else "0"
    # Empty means "not specified": the browser then keeps whatever size the
    # user last set in-page (localStorage), falling back to the 13px default.
    if args.font_size is not None:
        if not 8 <= args.font_size <= 32:
            print(f"proprep-web: --font-size must be 8-32 (got {args.font_size})",
                  file=sys.stderr)
            return 2
        os.environ["PROPREP_WEB_FONT_SIZE"] = str(args.font_size)
    else:
        os.environ["PROPREP_WEB_FONT_SIZE"] = ""
    os.environ["PROPREP_WEB_NO_TRANSCRIPT"] = "1" if args.no_transcript else "0"

    if sys.platform.startswith("win"):
        print("proprep.web requires POSIX (pty module). Windows is not supported.", file=sys.stderr)
        return 2

    try:
        import uvicorn  # noqa: F401
    except ImportError:
        print(
            "proprep.web requires the 'web' extras. Install with:\n"
            "    pip install 'proprep[web]'",
            file=sys.stderr,
        )
        return 2

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    requested_port = args.port
    if args.strict_port:
        port = requested_port
    else:
        try:
            port = _find_free_port(args.host, requested_port)
        except OSError as exc:
            print(f"proprep-web: {exc}", file=sys.stderr)
            return 2
        if port != requested_port:
            print(f"  (port {requested_port} was busy; using {port})")

    url = f"http://{args.host}:{port}/"
    print(f"ProPrep Web Shell listening on {url}")
    print("  - Terminal pane runs the unmodified ProPrep CLI")
    print("  - Structure Viewer will dock on the right when launched from the CLI")
    if args.no_transcript:
        print("  - Session transcript: disabled (--no-transcript)")
    else:
        print("  - Session transcript: ~/.proprep/web_sessions/session_<timestamp>.log")
    print()

    if not args.no_browser:
        threading.Thread(target=lambda: (time.sleep(0.6), webbrowser.open(url)), daemon=True).start()

    import uvicorn
    uvicorn.run(
        "proprep.web.server:app",
        host=args.host,
        port=port,
        log_level=args.log_level,
        reload=False,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
