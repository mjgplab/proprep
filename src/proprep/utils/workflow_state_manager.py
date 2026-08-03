"""
Workflow State Manager for ProPrep

Manages the current stage in the guided workflow mode and handles
stage transitions.
"""

from typing import Literal, Optional, List
from pathlib import Path
import json
from datetime import datetime


# Define workflow stages
# WorkflowStage = Literal["load", "detect", "compare", "fix", "simulate", "analyze"]  # DISABLED for pre-analysis release
WorkflowStage = Literal["load", "detect", "compare", "fix", "simulate"]  # See docs/developer/ANALYSIS_MODULES_REENABLE.md


class WorkflowStateManager:
    """Manages workflow progression state."""

    # Define the ordered workflow stages
    # STAGES = ["load", "detect", "compare", "fix", "simulate", "analyze"]  # DISABLED for pre-analysis release
    STAGES = ["load", "detect", "compare", "fix", "simulate"]  # See docs/developer/ANALYSIS_MODULES_REENABLE.md

    # Map stages to display names
    # STAGE_NAMES = {  # DISABLED for pre-analysis release
    STAGE_NAMES = {
        "load": "Load",
        "detect": "Detect",
        "compare": "Compare",
        "fix": "Fix",
        "simulate": "Simulate",
        # "analyze": "Analyze"  # See docs/developer/ANALYSIS_MODULES_REENABLE.md
    }

    # Display ordering for native stages — keeps a stable position
    # when plugin-contributed stages get spliced in via
    # ``register_plugin_stage``. Plugins choose their own order
    # value (``PluginMetadata.stage_order``); native stages occupy
    # 10/20/30/40/50 so a plugin defaulting to 999 lands after them.
    STAGE_ORDER = {
        "load": 10,
        "detect": 20,
        "compare": 30,
        "fix": 40,
        "simulate": 50,
    }

    @classmethod
    def register_plugin_stage(cls, meta) -> None:
        """Splice a plugin-declared stage into ``STAGES`` /
        ``STAGE_NAMES`` / ``STAGE_ORDER`` if it isn't already present.

        Idempotent: re-registering the same stage from another plugin
        (two plugins targeting the same new stage) leaves the existing
        entry alone — the stage's ``stage_name`` and ``stage_order``
        are taken from whichever plugin registered first. The second
        plugin's tool still gets added under the same stage by the
        menu's ``register_plugin_tool``.

        Class-level state mutation is intentional — there's only one
        ProPrep run per process, and the dict's lifetime matches the
        process. ``meta`` is duck-typed as ``PluginMetadata`` to keep
        ``utils/`` importable without dragging in the plugin package.
        """
        stage = meta.stage
        if stage in cls.STAGES:
            return  # already registered — first plugin wins
        # Insert in stage_order position so menus render in the
        # intended sequence even when plugins land mid-list.
        cls.STAGE_ORDER[stage] = meta.stage_order
        cls.STAGE_NAMES[stage] = meta.stage_name
        # Re-sort STAGES by STAGE_ORDER so the listing reflects the
        # new entry's position. Stable sort preserves prior ordering
        # for equal keys.
        cls.STAGES = sorted(
            list(cls.STAGES) + [stage],
            key=lambda s: cls.STAGE_ORDER.get(s, 999),
        )

    def __init__(self, workspace=None):
        """
        Initialize workflow state manager.

        Args:
            workspace: Optional workspace object to sync state with
        """
        self.workspace = workspace
        self._current_stage = "load"  # Always start at load
        self._completed_stages = set()
        self._skipped_stages = set()

        # Load state from workspace if available
        if workspace:
            self._load_from_workspace()

    def _load_from_workspace(self):
        """Load workflow state from workspace."""
        if self.workspace:
            self._current_stage = self.workspace.get("workflow_current_stage", "load")
            self._completed_stages = set(self.workspace.get("workflow_completed_stages", []))
            self._skipped_stages = set(self.workspace.get("workflow_skipped_stages", []))

    def _save_to_workspace(self):
        """Save workflow state to workspace."""
        if self.workspace:
            self.workspace.set("workflow_current_stage", self._current_stage)
            self.workspace.set("workflow_completed_stages", list(self._completed_stages))
            self.workspace.set("workflow_skipped_stages", list(self._skipped_stages))

    def get_current_stage(self) -> WorkflowStage:
        """Get the current workflow stage."""
        return self._current_stage

    def get_stage_index(self, stage: WorkflowStage) -> int:
        """Get the index of a stage in the workflow."""
        return self.STAGES.index(stage)

    def get_current_stage_index(self) -> int:
        """Get the index of the current stage."""
        return self.get_stage_index(self._current_stage)

    def is_stage_completed(self, stage: WorkflowStage) -> bool:
        """Check if a stage has been completed."""
        return stage in self._completed_stages

    def is_stage_skipped(self, stage: WorkflowStage) -> bool:
        """Check if a stage has been skipped."""
        return stage in self._skipped_stages

    def is_stage_current(self, stage: WorkflowStage) -> bool:
        """Check if a stage is the current one."""
        return stage == self._current_stage

    def is_stage_future(self, stage: WorkflowStage) -> bool:
        """Check if a stage is in the future (not yet reached)."""
        current_idx = self.get_current_stage_index()
        stage_idx = self.get_stage_index(stage)
        return stage_idx > current_idx

    def mark_completed(self, stage: Optional[WorkflowStage] = None):
        """
        Mark a stage as completed.

        Args:
            stage: Stage to mark as completed. If None, marks current stage.
        """
        if stage is None:
            stage = self._current_stage

        self._completed_stages.add(stage)
        # Remove from skipped if it was there
        self._skipped_stages.discard(stage)
        self._save_to_workspace()

    def mark_skipped(self, stage: Optional[WorkflowStage] = None):
        """
        Mark a stage as skipped.

        Args:
            stage: Stage to mark as skipped. If None, marks current stage.
        """
        if stage is None:
            stage = self._current_stage

        self._skipped_stages.add(stage)
        # Don't mark as completed if skipped
        self._completed_stages.discard(stage)
        self._save_to_workspace()

    def advance_stage(self, auto_advance: bool = False):
        """
        Move to the next workflow stage.

        Args:
            auto_advance: If True, automatically mark current as completed.
                         Only happens for 'load' stage.
        """
        current_idx = self.get_current_stage_index()

        # Auto-advance only works for load stage
        if auto_advance and self._current_stage == "load":
            self.mark_completed(self._current_stage)

        # Move to next stage if not at the end
        if current_idx < len(self.STAGES) - 1:
            self._current_stage = self.STAGES[current_idx + 1]
            self._save_to_workspace()

    def go_back(self):
        """
        Go back to the previous workflow stage.
        Marks the current stage as neither completed nor skipped.
        """
        current_idx = self.get_current_stage_index()

        # Revert current stage (remove from completed/skipped)
        self._completed_stages.discard(self._current_stage)
        self._skipped_stages.discard(self._current_stage)

        # Move to previous stage if not at the beginning
        if current_idx > 0:
            self._current_stage = self.STAGES[current_idx - 1]
            self._save_to_workspace()

    def can_go_back(self) -> bool:
        """Check if we can go back to a previous stage."""
        return self.get_current_stage_index() > 0

    def can_advance(self) -> bool:
        """Check if we can advance to the next stage."""
        return self.get_current_stage_index() < len(self.STAGES) - 1

    def get_stage_status(self, stage: WorkflowStage) -> Literal["completed", "current", "skipped", "future"]:
        """
        Get the status of a stage.

        Returns:
            "completed", "current", "skipped", or "future"
        """
        if self.is_stage_completed(stage):
            return "completed"
        elif self.is_stage_current(stage):
            return "current"
        elif self.is_stage_skipped(stage):
            return "skipped"
        else:
            return "future"

    def get_stage_symbol(self, stage: WorkflowStage) -> str:
        """
        Get the display symbol for a stage based on its status.

        Returns:
            Symbol: ✓ (completed), ▶ (current), ⊘ (skipped), ○ (future)
        """
        status = self.get_stage_status(stage)
        symbols = {
            "completed": "✓",
            "current": "▶",
            "skipped": "⊘",
            "future": "○"
        }
        return symbols[status]

    def get_stage_color(self, stage: WorkflowStage) -> str:
        """
        Get the display color for a stage based on its status.

        Returns:
            Rich color name: green, yellow, or grey50

        ``grey50`` (a concrete mid-grey) rather than the ``dim``
        attribute: ``dim`` lowers contrast *relative to the
        background*, so it fades to near-invisible on a light/white
        terminal, whereas ``grey50`` reads against both light and dark
        backgrounds.
        """
        status = self.get_stage_status(stage)
        # ``current`` uses ``dark_violet`` (#8700d7), not ``yellow`` or an
        # orange: yellow is near-invisible on white (~1.x:1 contrast), and
        # any orange collides semantically with ProPrep's unmet-prereq ⚠
        # warnings (also orange). Violet is a distinct hue from green
        # (completed), grey50 (future/skipped), AND the warning orange, so
        # the "you are here" stage reads as active rather than as an alert,
        # while clearing the bold-text contrast floor on white and dark.
        colors = {
            "completed": "green",
            "current": "dark_violet",
            "skipped": "grey50",
            "future": "grey50"
        }
        return colors[status]

    def get_all_stages(self) -> List[WorkflowStage]:
        """Get list of all workflow stages in order."""
        return self.STAGES.copy()

    def jump_to_stage(self, target: WorkflowStage):
        """
        Jump directly to a target stage, auto-skipping intermediate stages.

        Stages between current and target (when jumping forward) are marked
        as skipped. When jumping backward, intermediate stages are cleared.
        """
        target_idx = self.get_stage_index(target)
        current_idx = self.get_current_stage_index()

        if target_idx > current_idx:
            # Jumping forward — mark intermediate stages as skipped
            for i in range(current_idx, target_idx):
                stage = self.STAGES[i]
                if stage not in self._completed_stages:
                    self._skipped_stages.add(stage)
        elif target_idx < current_idx:
            # Jumping backward — clear intermediate stages
            for i in range(target_idx, current_idx + 1):
                stage = self.STAGES[i]
                self._completed_stages.discard(stage)
                self._skipped_stages.discard(stage)

        self._current_stage = target
        self._save_to_workspace()

    def reset(self):
        """Reset workflow state to beginning."""
        self._current_stage = "load"
        self._completed_stages.clear()
        self._skipped_stages.clear()
        self._save_to_workspace()
