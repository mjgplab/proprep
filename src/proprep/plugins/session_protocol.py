"""
SessionRecorder Protocol — the cross-package session-recording contract.

Both ProPrep's and etanalyze's existing ``SessionRecorder`` classes
satisfy this Protocol structurally; nothing has to change in their
recorder implementations. Defined here (rather than in
``proprep/utils/session_recorder.py``) so plugins can ``from
proprep.plugins import SessionRecorder`` without dragging in ProPrep's
concrete recorder + interceptor code.

When ProPrep launches a plugin it can pass its recorder into the
plugin via ``HostContext.session_recorder``; the plugin then routes
its own interactions through the host recorder so a single session
file captures the whole run. Plugins that don't want to participate
ignore the field; standalone plugin runs (no host) leave it ``None``.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Protocol, runtime_checkable


@runtime_checkable
class SessionRecorder(Protocol):
    """Structural type for a session-recording sink.

    The three methods plus the three attributes below are all the
    plugin and the host need to know about each other's recorder.
    Anything richer (replay, hybrid, file-format details) stays
    internal to whichever package owns the concrete implementation.
    """

    recording: bool
    """Whether the recorder is currently writing events."""

    record_file: Optional[str]
    """Path to the on-disk session file, when set."""

    session_data: Dict[str, Any]
    """The in-memory session dict (mirrors the on-disk JSON)."""

    def record_interaction(
        self,
        interaction_type: str,
        prompt: str,
        response: str,
        choices: Optional[List[str]] = None,
        context: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Append one interaction to the session.

        ``context`` is an extensible dict; the host adds ``module`` /
        ``description`` / ``options_map`` / ``option_label`` post-hoc
        via the prompt-utility wrappers. Plugins that route writes
        through the host recorder can also stamp a ``plugin_scope``
        key into ``context`` so a future replayer can distinguish
        host-origin from plugin-origin events.
        """

    def start_recording(
        self, metadata: Optional[Dict[str, Any]] = None
    ) -> None:
        """Begin a session. Idempotent for already-active recorders."""

    def stop_recording(self) -> None:
        """Finalise the session file. Idempotent."""
