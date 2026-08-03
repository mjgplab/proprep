"""ProPrep Web Shell.

A browser front-end that hosts the unmodified ProPrep CLI in an xterm.js
terminal pane and (in later phases) docks the live NGL structure viewer
beside it. Phase A is a PTY-based terminal multiplexer with no module
changes; later phases promote individual prompts to web widgets.

Launch with: ``python -m proprep.web`` (POSIX-only).
"""
