"""Regression test for PDB Filter structure-object resolution.

Bug: after generating a biological assembly (e.g. 1RAC -> 12 chains), the
Biological Assembly Generator stores only the assembly FILE PATH in the
workspace (``biological_assembly_pdb_file``), never a parsed
``biological_assembly_structure`` object. When PDB Filter selected the
assembly, ``_get_structure_from_workspace`` matched the path but found no
cached structure object, then SILENTLY fell back to a priority lookup that
returned the original 4-chain ``rcsb_structure``. The worker prefers a
passed structure over re-parsing the file, so the interface/BSA and
composition analyses ran on 4 chains (the asymmetric unit) instead of 12.

Fix: when the selected file matches a registered key but has no cached
structure object, the resolver returns ``None`` so the worker parses THAT
file, rather than substituting a structure for a different file.
"""
import types

from rich.console import Console

from proprep.structure_prep.pdb_filter import PDBFilterModule


def _module():
    """A PDBFilterModule instance with just enough wiring to exercise
    ``_get_structure_from_workspace`` (``console`` is a property backed by
    ``processor.console``); __init__ is bypassed deliberately."""
    m = object.__new__(PDBFilterModule)
    m.processor = types.SimpleNamespace(console=Console())
    return m


def test_selected_file_without_cached_structure_returns_none():
    """Selecting the biological assembly (path present, no structure object,
    a different structure cached) must NOT return the wrong cached object."""
    m = _module()
    rcsb_struct = object()  # the 4-chain asymmetric-unit object
    workspace = {
        "biological_assembly_pdb_file": "/tmp/assembly_12.pdb",
        "rcsb_pdb_file": "/tmp/orig_4.pdb",
        "rcsb_structure": rcsb_struct,
        # deliberately NO "biological_assembly_structure"
    }

    result = m._get_structure_from_workspace(
        workspace, pdb_file="/tmp/assembly_12.pdb", silent=True
    )

    assert result is None, (
        "resolver must return None so the worker parses the selected "
        f"assembly file; got {result!r} (would analyze the wrong structure)"
    )
    assert result is not rcsb_struct


def test_selected_file_with_cached_structure_is_returned():
    """Control: when the selected file DOES have a cached structure object,
    the resolver returns it unchanged."""
    m = _module()
    asm_struct = object()
    workspace = {
        "biological_assembly_pdb_file": "/tmp/assembly_12.pdb",
        "biological_assembly_structure": asm_struct,
        "rcsb_pdb_file": "/tmp/orig_4.pdb",
        "rcsb_structure": object(),
    }

    result = m._get_structure_from_workspace(
        workspace, pdb_file="/tmp/assembly_12.pdb", silent=True
    )

    assert result is asm_struct
