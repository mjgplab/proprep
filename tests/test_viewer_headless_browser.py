"""Never hand the NGL viewer page to a text-mode browser.

Python's ``webbrowser`` registers lynx/w3m/links as fallbacks whenever TERM
is set, even with no display at all. On a headless workstation that meant
``webbrowser.open()`` ran lynx against the WebGL viewer page, on the CLI's
own terminal, blocking until it crashed. These tests pin the guard that
decides when a browser may be opened, and the ssh -L hint printed instead.
"""

import os
import sys
from unittest import mock

import pytest

from proprep.structure_prep.viewer_server import (
    ViewerServer,
    gui_browser_unavailable_reason,
    ssh_forward_hint,
)


def _reason(env, platform="linux", browser_name=None):
    ctx = [
        mock.patch.dict(os.environ, env, clear=True),
        mock.patch.object(sys, "platform", platform),
    ]
    if browser_name is not None:
        fake = mock.Mock()
        fake.name = browser_name
        ctx.append(mock.patch("webbrowser.get", return_value=fake))
    with ctx[0], ctx[1]:
        if browser_name is None:
            return gui_browser_unavailable_reason()
        with ctx[2]:
            return gui_browser_unavailable_reason()


class TestGuiBrowserAvailability:
    def test_headless_ssh_session_is_refused(self):
        # TERM set, no DISPLAY: exactly the case that found lynx.
        assert _reason({"TERM": "xterm-256color"}) == "no DISPLAY or WAYLAND_DISPLAY"

    def test_text_browser_refused_even_with_display(self):
        reason = _reason({"DISPLAY": ":0"}, browser_name="/usr/bin/lynx")
        assert reason is not None and "lynx" in reason

    @pytest.mark.parametrize("env", [{"DISPLAY": ":0"}, {"WAYLAND_DISPLAY": "wayland-0"}])
    def test_gui_browser_allowed(self, env):
        assert _reason(env, browser_name="firefox") is None

    def test_explicit_browser_env_var_wins(self):
        # The documented escape hatch: BROWSER=/bin/true, or a wrapper that
        # forwards the URL to the user's own machine.
        assert _reason({"BROWSER": "/bin/true", "TERM": "xterm"}) is None

    @pytest.mark.parametrize("platform", ["darwin", "win32"])
    def test_gui_platforms_never_blocked(self, platform):
        assert _reason({}, platform=platform) is None


class TestSshForwardHint:
    def test_uses_the_port_sshd_was_reached_on(self):
        env = {"SSH_CONNECTION": "192.168.1.5 51234 24.47.30.228 2222", "USER": "mjgp"}
        with mock.patch.dict(os.environ, env, clear=True):
            assert ssh_forward_hint(8765) == (
                "ssh -p 2222 -L 8765:localhost:8765 mjgp@24.47.30.228"
            )

    def test_default_port_omits_the_flag(self):
        env = {"SSH_CONNECTION": "10.0.0.2 4444 10.0.0.9 22", "USER": "mjgp"}
        with mock.patch.dict(os.environ, env, clear=True):
            assert ssh_forward_hint(8765) == "ssh -L 8765:localhost:8765 mjgp@10.0.0.9"

    def test_no_guess_when_not_over_ssh(self):
        with mock.patch.dict(os.environ, {"USER": "mjgp"}, clear=True):
            assert ssh_forward_hint(8765) is None


class TestServerStart:
    def _start(self, env, browser_name=None, port=8899):
        opened = []
        fake = mock.Mock()
        fake.name = browser_name or "firefox"
        with mock.patch.dict(os.environ, env, clear=True), \
                mock.patch.object(sys, "platform", "linux"), \
                mock.patch("webbrowser.get", return_value=fake), \
                mock.patch("webbrowser.open", side_effect=opened.append):
            server = ViewerServer(config={"structures": []}, structure_files=[], port=port)
            try:
                assert server.start(open_browser=True) is True
                return server, opened
            finally:
                server.stop()

    def test_headless_start_serves_but_opens_nothing(self):
        server, opened = self._start({"TERM": "xterm"})
        assert opened == []
        assert server.browser_opened is False
        assert server.headless_reason == "no DISPLAY or WAYLAND_DISPLAY"

    def test_desktop_start_still_opens_a_tab(self):
        server, opened = self._start({"DISPLAY": ":0"}, port=8898)
        assert opened == ["http://localhost:8898/viewer"]
        assert server.browser_opened is True
        assert server.headless_reason is None
