"""Depositing a second set must not destroy the first one's files.

Sets are metadata KEYS, not directories, so every set in one (redox, spin)
state shares a single folder. The copy was a bare shutil.copy2, so two sets
whose files happened to share a basename silently overwrote each other while
both metadata leaves still pointed at that name.
"""

import pytest

from proprep.forcefield_params.user_library import (
    PromotionRequest, _copy_into, _reserved_basenames, _rollback,
)


def _request(set_name="set_b"):
    return PromotionRequest(
        family="small_molecules", type_name="gdp", set_name=set_name,
        frcmod_src="x.frcmod", lib_srcs=["x.lib"],
        redox_state="single_state", spin_state="default",
    )


METADATA = {
    "redox_states": {
        "single_state": {
            "spin_states": {
                "default": {
                    "forcefield_sets": {
                        "set_a": {"files": {"frcmod": "GDP.frcmod", "lib": "GDP.lib"}},
                        "set_b": {"files": {"frcmod": "own.frcmod", "lib": "own.lib"}},
                    }
                }
            }
        }
    }
}


# --------------------------------------------------------------------------
# which names are spoken for
# --------------------------------------------------------------------------

def test_other_sets_files_are_reserved():
    assert _reserved_basenames(METADATA, _request(), "set_b") == {"GDP.frcmod", "GDP.lib"}


def test_a_sets_own_files_are_not_reserved_against_itself():
    """Re-depositing the same set must still overwrite its own files."""
    reserved = _reserved_basenames(METADATA, _request(), "set_a")
    assert "GDP.lib" not in reserved
    assert "own.lib" in reserved


def test_list_valued_file_entries_are_covered():
    metadata = {"redox_states": {"single_state": {"spin_states": {"default": {
        "forcefield_sets": {"set_a": {"files": {"lib": ["a.lib", "b.lib"]}}}}}}}}
    assert _reserved_basenames(metadata, _request(), "set_b") == {"a.lib", "b.lib"}


# --------------------------------------------------------------------------
# the copy itself
# --------------------------------------------------------------------------

def test_a_reserved_name_with_different_content_is_written_alongside(tmp_path):
    dest = tmp_path / "state"
    dest.mkdir()
    (dest / "GDP.lib").write_text("set A's parameters")
    src = tmp_path / "GDP.lib"
    src.write_text("set B's parameters")

    created, replaced = [], []
    name = _copy_into(str(src), dest, created, replaced, reserved={"GDP.lib"})

    assert name == "GDP_2.lib"
    assert (dest / "GDP.lib").read_text() == "set A's parameters"
    assert (dest / "GDP_2.lib").read_text() == "set B's parameters"


def test_an_identical_file_is_shared_rather_than_duplicated(tmp_path):
    dest = tmp_path / "state"
    dest.mkdir()
    (dest / "GDP.lib").write_text("same bytes")
    src = tmp_path / "GDP.lib"
    src.write_text("same bytes")

    created, replaced = [], []
    name = _copy_into(str(src), dest, created, replaced, reserved={"GDP.lib"})

    assert name == "GDP.lib"
    assert not (dest / "GDP_2.lib").exists()


def test_an_unreserved_name_is_written_as_is(tmp_path):
    dest = tmp_path / "state"
    dest.mkdir()
    src = tmp_path / "new.lib"
    src.write_text("fresh")

    created, replaced = [], []
    assert _copy_into(str(src), dest, created, replaced) == "new.lib"
    assert created == [dest / "new.lib"]


# --------------------------------------------------------------------------
# rollback
# --------------------------------------------------------------------------

def test_rollback_restores_a_file_the_deposit_overwrote(tmp_path):
    dest = tmp_path / "state"
    dest.mkdir()
    victim = dest / "GDP.lib"
    victim.write_text("the parameters that were already there")
    src = tmp_path / "GDP.lib"
    src.write_text("the replacement")

    created, replaced = [], []
    _copy_into(str(src), dest, created, replaced)          # deliberate overwrite
    assert victim.read_text() == "the replacement"
    assert created == []                                    # nothing NEW was made

    _rollback(tmp_path / "metadata.json", None, created, [], replaced)
    assert victim.read_text() == "the parameters that were already there"


def test_rollback_still_removes_files_the_deposit_created(tmp_path):
    dest = tmp_path / "state"
    dest.mkdir()
    src = tmp_path / "new.lib"
    src.write_text("fresh")

    created, replaced = [], []
    _copy_into(str(src), dest, created, replaced)
    _rollback(tmp_path / "metadata.json", None, created, [], replaced)
    assert not (dest / "new.lib").exists()
