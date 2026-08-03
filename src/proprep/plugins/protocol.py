"""
ProPrepPlugin — the structural contract every plugin satisfies.

Three methods, no inheritance required: any class with the right
signatures is a plugin. The Protocol is ``runtime_checkable`` so
``isinstance(obj, ProPrepPlugin)`` works in tests, but the production
discovery path doesn't rely on it (plugins that *fail* the contract
silently surface as exceptions during ``get_metadata`` / ``launch``).
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from .host import HostContext
from .metadata import PluginMetadata


@runtime_checkable
class ProPrepPlugin(Protocol):
    """The plugin contract.

    Plugins typically expose a module-level ``PLUGIN`` instance:

        # mypkg/proprep_plugin.py
        class MyPlugin:
            def get_metadata(self): ...
            def is_available(self): ...
            def launch(self, host): ...
        PLUGIN = MyPlugin()

    and register it in ``pyproject.toml``:

        [project.entry-points."proprep.plugins"]
        mypkg = "mypkg.proprep_plugin:PLUGIN"
    """

    def get_metadata(self) -> PluginMetadata:
        """Return the plugin's static metadata. Called once per
        ProPrep startup. Cheap; pure; no side effects."""

    def is_available(self) -> bool:
        """Self-check: do the plugin's runtime dependencies resolve
        in this env? Called once per ProPrep startup before menu
        splicing. Cheap (import probes, file-existence checks); a
        return value of ``False`` quietly drops the plugin from the
        menu so users don't see a tool that would crash on launch."""

    def launch(self, host: HostContext) -> None:
        """Run the plugin in-process. Returns when the user exits
        the plugin's UI; exceptions propagate but ProPrep's dispatch
        wrapper catches them and prints a panel rather than killing
        the host."""
