"""
ProPrep plugin infrastructure.

External packages register here via the ``proprep.plugins`` entry-point
group; ProPrep discovers them at startup, validates the contract, and
splices their declared stage/tool into the menu. Plugins that don't
match the current ``PLUGIN_API_VERSION`` are dropped with a warning.

Public surface:

* :class:`ProPrepPlugin` — the structural contract every plugin
  satisfies (``get_metadata``, ``is_available``, ``launch``).
* :class:`PluginMetadata` — the dataclass plugins return from
  ``get_metadata()`` describing where they live in the menu and which
  workspace keys they want pulled from the host.
* :class:`HostContext` — the package ProPrep hands to ``launch``: a
  working directory, a seed-state dict (declarative pull of host
  workspace keys), an optional shared session recorder, console, and
  child logger.
* :class:`SessionRecorder` — runtime-checkable Protocol both ProPrep's
  and any plugin's recorder satisfy structurally; lets the host pass
  its recorder into the plugin so a single session file captures
  everything that happened.
* :func:`discover_plugins` — entry-point scan, returning the validated
  ``{name: plugin}`` mapping.

Design rationale lives in
``etanalyze/docs/PROPREP_PLUGIN_INTEGRATION.md``.
"""

from .discovery import discover_plugins
from .host import HostContext
from .metadata import PLUGIN_API_VERSION, PluginMetadata
from .protocol import ProPrepPlugin
from .session_protocol import SessionRecorder

__all__ = [
    "PLUGIN_API_VERSION",
    "HostContext",
    "PluginMetadata",
    "ProPrepPlugin",
    "SessionRecorder",
    "discover_plugins",
]
