"""
HostContext — what ProPrep hands to a plugin's ``launch`` method.

Every field is optional except ``working_dir``. Plugins must tolerate
missing keys in ``seed_state`` (the host workspace may not yet hold
the value the plugin would like).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from logging import Logger
from pathlib import Path
from typing import Any, Dict, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    # Imported under TYPE_CHECKING so HostContext doesn't drag the
    # rich.console module into every importer of this file.
    from rich.console import Console
    from .session_protocol import SessionRecorder


@dataclass
class HostContext:
    """Bundle ProPrep hands to a plugin at launch time."""

    working_dir: Path
    """Process CWD when the plugin was invoked. Plugins should treat
    this as an immutable hint; if they ``cd`` they take responsibility
    for restoring it."""

    seed_state: Dict[str, Any] = field(default_factory=dict)
    """Shallow copy of host workspace values for the keys the plugin
    declared via ``PluginMetadata.consumes_workspace_keys``. Any key
    the host doesn't have is silently omitted — never present as
    ``None`` — so plugins disambiguate "absent" from "explicitly
    null"."""

    session_recorder: Optional["SessionRecorder"] = None
    """Host's active session recorder, when one is running. The
    plugin should write its interactions through this recorder so a
    single session file captures the run; plugins that wrap it in a
    scope-annotating proxy take responsibility for the annotation
    convention."""

    console: Optional["Console"] = None
    """Shared Rich console. Falls back to ``None`` only when the host
    is in a non-TTY context (rare); plugins should create their own
    if absent."""

    logger: Optional[Logger] = None
    """Per-plugin child logger named ``proprep.plugin.<name>``. Users
    can filter all plugin output via the parent ``proprep.plugin``
    logger."""

    api_version: int = 1
    """Mirror of ``PLUGIN_API_VERSION`` at the time the host built the
    context. Plugins can read this to behave conservatively when the
    host is older than they expect."""
