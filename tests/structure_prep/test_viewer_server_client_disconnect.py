"""
A browser that hung up is not a server error.

Running MODELLER logged, mid-run:

    ERROR - viewer_server.py:105 - Error serving version: [Errno 32] Broken pipe

The viewer page polls /version every 1.5 s and never aborts, and the server is a
single-threaded HTTPServer. A compiled extension that holds the GIL for minutes
(MODELLER ran ~2.5 min here) starves the server thread, so polls pile up on the
listen backlog and their sockets are closed by the time each one is answered.
Ordinary, and already handled: handle_one_request has caught BrokenPipeError
since it was written, with the comment "This is normal behavior - just ignore
it".

What defeated that net was the per-handler `except Exception`, which caught the
BrokenPipeError first — logged it at ERROR and then called send_error, writing a
500 to the same closed socket.
"""

import logging

import pytest

from proprep.structure_prep.viewer_server import ViewerHTTPRequestHandler


class _Handler(ViewerHTTPRequestHandler):
    """Bare handler: no socket, no BaseHTTPRequestHandler __init__."""

    def __init__(self):  # noqa: D107 - deliberately skips the base __init__
        self.sent_errors = []
        self.raise_on_send_error = None

    def send_error(self, code, message=None, explain=None):
        self.sent_errors.append((code, message))
        if self.raise_on_send_error:
            raise self.raise_on_send_error


@pytest.mark.parametrize("exc", [
    BrokenPipeError(32, "Broken pipe"),
    ConnectionResetError(54, "Connection reset by peer"),
    ConnectionAbortedError(53, "Software caused connection abort"),
])
def test_client_disconnect_is_not_an_error(exc, caplog):
    h = _Handler()
    with caplog.at_level(logging.DEBUG, logger="proprep.structure_prep.viewer_server"):
        h._handle_serve_error("version", exc)

    assert h.sent_errors == [], "must not write a 500 to a closed socket"
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING], \
        "a normal disconnect was logged as a fault"
    assert any("disconnected" in r.getMessage() for r in caplog.records)


def test_a_real_failure_is_still_reported(caplog):
    h = _Handler()
    with caplog.at_level(logging.DEBUG, logger="proprep.structure_prep.viewer_server"):
        h._handle_serve_error("config", ValueError("bad json"))

    assert h.sent_errors == [(500, "Internal server error: bad json")]
    assert any(r.levelno == logging.ERROR for r in caplog.records)


def test_a_client_that_vanishes_mid_500_does_not_raise(caplog):
    """The 500 for a real failure can itself hit a closed socket."""
    h = _Handler()
    h.raise_on_send_error = BrokenPipeError(32, "Broken pipe")

    with caplog.at_level(logging.DEBUG, logger="proprep.structure_prep.viewer_server"):
        h._handle_serve_error("structure", ValueError("boom"))  # must not raise

    assert h.sent_errors == [(500, "Internal server error: boom")]


def test_serve_version_survives_a_dead_socket(caplog):
    """The reported path end to end."""
    class _DeadSocketHandler(_Handler):
        def send_response(self, *a, **k):
            pass

        def send_header(self, *a, **k):
            pass

        def end_headers(self):
            pass

        @property
        def wfile(self):
            class _W:
                def write(self, _data):
                    raise BrokenPipeError(32, "Broken pipe")
            return _W()

    h = _DeadSocketHandler()
    with caplog.at_level(logging.DEBUG, logger="proprep.structure_prep.viewer_server"):
        h.serve_version()   # must not raise

    assert h.sent_errors == []
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


def test_outer_net_still_swallows_a_disconnect():
    """handle_one_request's existing guard must remain in place."""
    class _Boom(ViewerHTTPRequestHandler):
        def __init__(self):
            pass

        def _super_call(self):
            raise BrokenPipeError(32, "Broken pipe")

    import http.server

    calls = []

    def fake_super(self):
        calls.append(1)
        raise BrokenPipeError(32, "Broken pipe")

    original = http.server.BaseHTTPRequestHandler.handle_one_request
    http.server.BaseHTTPRequestHandler.handle_one_request = fake_super
    try:
        _Boom().handle_one_request()   # must not raise
    finally:
        http.server.BaseHTTPRequestHandler.handle_one_request = original

    assert calls == [1]
