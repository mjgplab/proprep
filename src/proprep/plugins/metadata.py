"""
PluginMetadata — what every plugin tells ProPrep about itself.

Returned from ``ProPrepPlugin.get_metadata()`` at discovery time.
ProPrep uses it to (a) splice the plugin into the workflow menu under
the declared stage, (b) validate the API version against the host's
``PLUGIN_API_VERSION``, and (c) build the seed_state dict by pulling
the declared workspace keys from its own workspace.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional


# Bumped only on intentional breaking changes to the plugin API. When
# this number changes, every plugin author has to rev their adapter;
# discover_plugins drops mismatched plugins with a warning naming the
# expected vs. observed versions.
PLUGIN_API_VERSION = 1


@dataclass
class PluginMetadata:
    """Where the plugin lives in the menu + what it wants from the host."""

    name: str
    """Stable identifier; appears in the menu sentinel key
    (``__plugin::<name>``) and as the logger child name
    (``proprep.plugin.<name>``)."""

    version: str
    """Plugin's own version string. Surfaced in the menu / errors."""

    description: str
    """One-liner shown next to the menu entry."""

    stage: str
    """Workflow stage the plugin's tool lives under. If the stage
    doesn't exist yet in ``WorkflowStateManager.STAGES``, ProPrep
    extends the list."""

    stage_name: str
    """Human-readable stage label (e.g. ``"Analyze"``)."""

    stage_description: str = ""
    """Optional longer description used in workflow guides."""

    stage_order: int = 999
    """Display order of the stage; lower comes first. Native ProPrep
    stages occupy 10, 20, 30, 40, 50; plugin defaults at 999 so they
    sort after natives unless they choose otherwise."""

    tool_order: int = 0
    """Order among siblings under the same stage. Ties broken
    alphabetically by plugin name."""

    api_version: int = PLUGIN_API_VERSION
    """Plugin's declared API version. ``discover_plugins`` drops
    plugins where this doesn't match the host's ``PLUGIN_API_VERSION``."""

    consumes_workspace_keys: List[str] = field(default_factory=list)
    """Declarative pull list. ProPrep's ``build_plugin_seed`` reads
    each key from its workspace and copies the value verbatim into
    ``HostContext.seed_state`` if present. Adding a new shared concept
    is a one-line edit on the plugin side; ProPrep's seed builder
    stays generic."""

    display_name: Optional[str] = None
    """Human-friendly label rendered in the menu list. When unset,
    falls back to ``name`` — but ``name`` doubles as the sentinel
    identifier (lowercase, no spaces) so plugins with mixed-case
    branding (``"ETAnalyze"``, ``"VMD-render"``) should set this
    field explicitly."""

    stage_shortcut: Optional[str] = None
    """Single-character shortcut for jumping to the plugin's stage
    in the workflow menu's nav footer (alongside ProPrep's native
    [l]oad / [d]etect / etc.). Lowercase letter; collisions with
    existing shortcuts (or with another plugin's) are silently
    dropped at registration time and logged as a warning, so a
    plugin's stage remains reachable via [n]ext / [b]ack even when
    its shortcut is unavailable."""
