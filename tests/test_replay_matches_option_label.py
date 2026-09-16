"""A recorded menu answer replays by its option label, not its bare number.

Menus get reordered between releases (1.19.1 moved the MD engine menu to
sander, sander.MPI, pmemd, pmemd.MPI, pmemd.cuda; pmemd went from 2 to 3).
The recording stores both the typed key and the label it meant, so replay
re-resolves the key against the current options_map.
"""

from __future__ import annotations

import types

from proprep.utils.prompts import resolve_replayed_key

OLD_MAP = {"1": "sander (single CPU)", "2": "pmemd (single CPU, optimized)",
           "3": "pmemd.MPI (multi-CPU)", "4": "pmemd.cuda (GPU acceleration)"}
NEW_MAP = {"1": "sander (single CPU)", "2": "sander.MPI (multi-CPU, AmberTools)",
           "3": "pmemd (single CPU, optimized)", "4": "pmemd.MPI (multi-CPU)",
           "5": "pmemd.cuda (GPU acceleration)"}


def _processor(replaying, label):
    replayer = types.SimpleNamespace(
        replaying=replaying,
        last_returned_interaction={"context": {"option_label": label}} if label else None,
    )
    return types.SimpleNamespace(session_manager=types.SimpleNamespace(replayer=replayer, recorder=None))


def test_recorded_pmemd_replays_as_pmemd_after_reorder(capsys):
    proc = _processor(True, OLD_MAP["2"])
    assert resolve_replayed_key(proc, "2", NEW_MAP) == "3"
    assert "now means a different option" in capsys.readouterr().out


def test_unchanged_key_passes_through(capsys):
    proc = _processor(True, OLD_MAP["1"])
    assert resolve_replayed_key(proc, "1", NEW_MAP) == "1"
    assert capsys.readouterr().out == ""


def test_live_input_is_untouched():
    proc = _processor(False, OLD_MAP["2"])
    assert resolve_replayed_key(proc, "2", NEW_MAP) == "2"


def test_no_recorded_label_keeps_the_key():
    proc = _processor(True, None)
    assert resolve_replayed_key(proc, "2", NEW_MAP) == "2"


def test_ambiguous_label_keeps_the_key():
    proc = _processor(True, "same")
    assert resolve_replayed_key(proc, "1", {"1": "other", "2": "same", "3": "same"}) == "1"


def test_no_options_map_or_processor():
    assert resolve_replayed_key(None, "2", NEW_MAP) == "2"
    assert resolve_replayed_key(_processor(True, "x"), "2", None) == "2"
