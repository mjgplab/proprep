"""
Plugin discovery via the ``proprep.plugins`` entry-point group.

Called once at ProPrep startup. Plugins that:

* fail to import,
* don't expose the right shape,
* declare an ``api_version`` mismatching the host's ``PLUGIN_API_VERSION``,
* return ``False`` from ``is_available()``,

are silently dropped (with a warning logged via the
``proprep.plugins`` logger). ProPrep continues either way — a
broken plugin never breaks the host.
"""

from __future__ import annotations

import logging
from importlib.metadata import entry_points
from typing import Dict

from .metadata import PLUGIN_API_VERSION
from .protocol import ProPrepPlugin

logger = logging.getLogger(__name__)


PLUGIN_ENTRY_POINT_GROUP = "proprep.plugins"


def discover_plugins() -> Dict[str, ProPrepPlugin]:
    """Scan installed packages for ProPrep plugins.

    Returns ``{plugin_name: plugin_instance}``. Empty dict when no
    entry-points are registered (the common case for a fresh ProPrep
    install with no plugins). Never raises.
    """
    found: Dict[str, ProPrepPlugin] = {}

    try:
        eps = entry_points(group=PLUGIN_ENTRY_POINT_GROUP)
    except Exception:  # noqa: BLE001 — defensive: importlib edge cases
        # entry_points() can raise on certain mis-installed packages
        # (PEP 517 builds with broken metadata). Treat as "no plugins".
        logger.warning(
            "Could not enumerate %r entry-points; "
            "no plugins will be loaded.",
            PLUGIN_ENTRY_POINT_GROUP,
            exc_info=True,
        )
        return found

    for ep in eps:
        try:
            plugin = ep.load()
        except Exception:  # noqa: BLE001 — plugin import failures are isolated
            logger.warning(
                "Plugin %r failed to import; skipping.",
                ep.name, exc_info=True,
            )
            continue

        try:
            meta = plugin.get_metadata()
        except Exception:  # noqa: BLE001
            logger.warning(
                "Plugin %r raised during get_metadata(); skipping.",
                ep.name, exc_info=True,
            )
            continue

        if meta.api_version != PLUGIN_API_VERSION:
            logger.warning(
                "Plugin %r declares api_version=%d but ProPrep expects "
                "%d; skipping. Upgrade %s%s.",
                meta.name, meta.api_version, PLUGIN_API_VERSION,
                "the plugin" if meta.api_version < PLUGIN_API_VERSION else "ProPrep",
                "",
            )
            continue

        try:
            available = plugin.is_available()
        except Exception:  # noqa: BLE001
            logger.warning(
                "Plugin %r raised during is_available(); skipping.",
                meta.name, exc_info=True,
            )
            continue

        if not available:
            logger.info(
                "Plugin %r reported is_available()=False; not surfacing.",
                meta.name,
            )
            continue

        if meta.name in found:
            logger.warning(
                "Plugin name collision on %r; later registration wins.",
                meta.name,
            )
        found[meta.name] = plugin

    return found
