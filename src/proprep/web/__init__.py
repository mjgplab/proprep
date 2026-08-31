"""ProPrep Web Shell.

A browser front-end that hosts the unmodified ProPrep CLI in an xterm.js
terminal pane and docks the live NGL structure viewer beside it.

Phase A: a PTY-based terminal multiplexer with no module changes.
Phase B: sessions live in *seats* that outlive their websockets (reconnect
with replay), a hosted mode serving N token-addressed seats from one
process for workshops, loopback-only viewer announcements, and a project
download. See ``seats.py`` and ``server.py``.

Launch with: ``proprep-web`` (local) or ``proprep-web --seats N`` (hosted).
POSIX-only.
"""
