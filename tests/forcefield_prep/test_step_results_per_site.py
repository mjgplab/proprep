"""
MCPB step_results belong to one site, not to the workspace.

They used to be restored in __init__ from `mcpb_step_results`, a single
workspace key every site shared, so step_1 belonged to whichever site wrote it
last. Three separate failures came from that:

  - site 1 was typed from site 2's standard.fingerprint. The PDB serial ranges
    do not overlap, so every atom became the 'XX' placeholder and tleap
    reported ~150 missing parameters.
  - site 1's RESP was fitted against site 2's charge constraint; a -1 model
    against -3 pushed a metal to +5.6.
  - "Step-1 records charge -3 for this site but its ESP was computed at -1".

Each was patched at its consumer. This partitions the storage instead, so a
consumer cannot read another site's step_1 in the first place.

The site is not known at construction -- callers assign provided_redox_site
immediately afterwards -- so that setter is where the restore happens.
"""

import pytest

from proprep.forcefield_prep.metal_site_parameterizer import (
    DEFAULT_SITE_KEY, MetalSiteWorkflowManager,
)


class _Workspace:
    def __init__(self, data=None):
        self.data = data or {}

    def get(self, key, default=None):
        return self.data.get(key, default)

    def set(self, key, value):
        self.data[key] = value


class _Processor:
    def __init__(self, workspace):
        self._workspace = workspace

    def _get_workspace(self):
        return self._workspace


class _Site:
    def __init__(self, site_id):
        self.site_id = site_id


def _manager(processor, site_id=None):
    manager = MetalSiteWorkflowManager(processor=processor)
    if site_id is not None:
        manager.provided_redox_site = _Site(site_id)
    return manager


def _run(processor, site_id, **step_1):
    manager = _manager(processor, site_id)
    manager.step_results["step_1"] = step_1
    manager._save_step_results()
    return manager


# --------------------------------------------------------------------------- #
# the reported failure
# --------------------------------------------------------------------------- #

def test_a_site_does_not_inherit_the_site_that_ran_last():
    workspace = _Workspace()
    processor = _Processor(workspace)

    _run(processor, "site_1", charge=-1, fingerprint="/run/site_1/standard.fingerprint")
    _run(processor, "site_2", charge=-3, fingerprint="/run/site_2/standard.fingerprint")

    revisited = _manager(processor, "site_1")

    assert revisited.step_results["step_1"]["charge"] == -1
    assert "site_1" in revisited.step_results["step_1"]["fingerprint"]


def test_saving_one_site_does_not_erase_another():
    workspace = _Workspace()
    processor = _Processor(workspace)

    _run(processor, "site_1", charge=-1)
    _run(processor, "site_2", charge=-3)

    stored = workspace.data["mcpb_step_results"]
    assert sorted(stored) == ["site_1", "site_2"]
    assert stored["site_1"]["step_1"]["charge"] == -1
    assert stored["site_2"]["step_1"]["charge"] == -3


def test_results_are_stored_under_the_site_id():
    workspace = _Workspace()
    _run(_Processor(workspace), "site_7", charge=0)

    assert "site_7" in workspace.data["mcpb_step_results"]


def test_a_site_with_no_saved_results_starts_empty():
    workspace = _Workspace()
    processor = _Processor(workspace)
    _run(processor, "site_1", charge=-1)

    fresh = _manager(processor, "site_2")

    assert fresh.step_results == {}


def test_reassigning_the_site_reloads_its_results():
    """One manager reused across sites must not carry results over."""
    workspace = _Workspace()
    processor = _Processor(workspace)
    _run(processor, "site_1", charge=-1)
    _run(processor, "site_2", charge=-3)

    manager = _manager(processor, "site_1")
    assert manager.step_results["step_1"]["charge"] == -1

    manager.provided_redox_site = _Site("site_2")
    assert manager.step_results["step_1"]["charge"] == -3


# --------------------------------------------------------------------------- #
# the standalone workflow and older workspaces
# --------------------------------------------------------------------------- #

def test_the_standalone_workflow_still_round_trips():
    """No site assigned: results go to a default slot and come back."""
    workspace = _Workspace()
    processor = _Processor(workspace)

    manager = _manager(processor)
    manager.step_results["step_1"] = {"charge": -2}
    manager._save_step_results()

    assert DEFAULT_SITE_KEY in workspace.data["mcpb_step_results"]
    assert _manager(processor).step_results["step_1"]["charge"] == -2


def test_a_workspace_written_before_partitioning_is_readable():
    """The legacy flat shape has step keys at the top level."""
    workspace = _Workspace({"mcpb_step_results": {"step_1": {"charge": -2}}})

    manager = _manager(_Processor(workspace))

    assert manager.step_results["step_1"]["charge"] == -2


def test_a_labelled_site_does_not_adopt_legacy_unlabelled_results():
    """
    Legacy results belong to SOME site but there is no telling which. Handing
    them to a named site would reintroduce exactly the bug being fixed.
    """
    workspace = _Workspace({"mcpb_step_results": {"step_1": {"charge": -3}}})

    manager = _manager(_Processor(workspace), "site_1")

    assert manager.step_results == {}


def test_saving_replaces_the_legacy_shape_rather_than_nesting_inside_it():
    workspace = _Workspace({"mcpb_step_results": {"step_1": {"charge": -3}}})
    processor = _Processor(workspace)

    _run(processor, "site_1", charge=-1)

    stored = workspace.data["mcpb_step_results"]
    assert "step_1" not in stored
    assert stored["site_1"]["step_1"]["charge"] == -1


# --------------------------------------------------------------------------- #
# the accessor other code uses
# --------------------------------------------------------------------------- #

def test_provided_redox_site_reads_back():
    """Call sites use getattr(self, 'provided_redox_site', None)."""
    manager = _manager(_Processor(_Workspace()), "site_1")

    assert getattr(manager, "provided_redox_site", None).site_id == "site_1"


def test_it_is_none_before_assignment():
    assert getattr(_manager(_Processor(_Workspace())),
                   "provided_redox_site", None) is None


def test_a_dict_site_is_accepted():
    """Resumed sessions hand dict-form sites around."""
    workspace = _Workspace()
    manager = MetalSiteWorkflowManager(processor=_Processor(workspace))
    manager.provided_redox_site = {"site_id": "site_4"}
    manager.step_results["step_1"] = {"charge": 0}
    manager._save_step_results()

    assert "site_4" in workspace.data["mcpb_step_results"]


def test_no_processor_is_not_an_error():
    manager = MetalSiteWorkflowManager()
    manager.provided_redox_site = _Site("site_1")

    assert manager.step_results == {}
