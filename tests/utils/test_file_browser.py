#!/usr/bin/env python3
"""Tests for the shared interactive file browser (proprep.utils.file_browser).

Covers:
- _parse_indices: bare-N, comma lists, inclusive ranges, dedup, range
  normalization, single-mode rejection of multi input, out-of-range.
- file_browser interactive flow (prompt_with_context monkeypatched to feed a
  scripted sequence of commands): bare-N file select, directory navigation,
  multi-select + ranges, q cancel, skip, invalid re-prompt.
- Filename-based session replay: the recorder remembers basenames, so when
  files are added/removed between record and replay the right file is still
  selected, and a vanished file falls through to live input.

Run with: pytest tests/utils/test_file_browser.py
"""

import sys
from pathlib import Path

import pytest
from rich.console import Console

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from proprep.utils import file_browser as fb
from proprep.utils.file_browser import file_browser, SKIP, _parse_indices


# --------------------------------------------------------------------------- #
# _parse_indices (pure)
# --------------------------------------------------------------------------- #

class TestParseIndices:
    def test_single(self):
        assert _parse_indices("3", count=5, multi=False) == [3]

    def test_single_rejects_multi(self):
        with pytest.raises(ValueError):
            _parse_indices("1,2", count=5, multi=False)

    def test_multi_commas(self):
        assert _parse_indices("1,3,5", count=5, multi=True) == [1, 3, 5]

    def test_multi_whitespace(self):
        assert _parse_indices("1 3 5", count=5, multi=True) == [1, 3, 5]

    def test_range(self):
        assert _parse_indices("1-5", count=5, multi=True) == [1, 2, 3, 4, 5]

    def test_mixed_range_and_list(self):
        assert _parse_indices("1-3,5,7", count=9, multi=True) == [1, 2, 3, 5, 7]

    def test_reversed_range_normalizes(self):
        assert _parse_indices("5-3", count=9, multi=True) == [3, 4, 5]

    def test_dedup_preserves_order(self):
        assert _parse_indices("3,1,3,2", count=9, multi=True) == [3, 1, 2]

    def test_overlapping_ranges_dedup(self):
        assert _parse_indices("1-3,2-4", count=9, multi=True) == [1, 2, 3, 4]

    def test_out_of_range(self):
        with pytest.raises(ValueError):
            _parse_indices("7", count=5, multi=True)

    def test_non_numeric(self):
        with pytest.raises(ValueError):
            _parse_indices("x", count=5, multi=True)

    def test_empty(self):
        with pytest.raises(ValueError):
            _parse_indices("   ", count=5, multi=True)


# --------------------------------------------------------------------------- #
# Helpers for interactive tests
# --------------------------------------------------------------------------- #

class _ScriptedPrompt:
    """Stand-in for prompt_with_context that returns scripted responses.

    Records the (module, description) of each call for assertions.
    """

    def __init__(self, responses):
        self._responses = list(responses)
        self.calls = []

    def __call__(self, processor, prompt, **kwargs):
        self.calls.append(kwargs)
        if not self._responses:
            raise AssertionError("scripted prompt ran out of responses")
        return self._responses.pop(0)


@pytest.fixture
def tree(tmp_path):
    """A small directory tree: 3 pdb files + a subdir with one more."""
    (tmp_path / "alpha.pdb").write_text("A")
    (tmp_path / "beta.pdb").write_text("B")
    (tmp_path / "gamma.pdb").write_text("C")
    (tmp_path / "notes.txt").write_text("ignore me")
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "delta.pdb").write_text("D")
    return tmp_path


def _run(monkeypatch, responses, tree, **kwargs):
    scripted = _ScriptedPrompt(responses)
    monkeypatch.setattr(fb, "prompt_with_context", scripted)
    result = file_browser(
        directory=str(tree),
        extensions=[".pdb"],
        console=Console(file=open("/dev/null", "w")),
        **kwargs,
    )
    return result, scripted


# --------------------------------------------------------------------------- #
# Interactive flow
# --------------------------------------------------------------------------- #

class TestInteractive:
    def test_bare_number_selects_file(self, monkeypatch, tree):
        # Items at root: [1] .. (parent), [2] sub/, [3] alpha, [4] beta, [5] gamma
        result, _ = _run(monkeypatch, ["3"], tree)
        assert Path(result).name == "alpha.pdb"

    def test_q_cancels(self, monkeypatch, tree):
        result, _ = _run(monkeypatch, ["q"], tree)
        assert result is None

    def test_extension_filter_hides_txt(self, monkeypatch, tree):
        # Only .pdb listed → indices never reach notes.txt. Selecting the last
        # file index lands on gamma.pdb, proving notes.txt isn't in the list.
        result, _ = _run(monkeypatch, ["5"], tree)
        assert Path(result).name == "gamma.pdb"

    def test_navigate_into_subdir_then_select(self, monkeypatch, tree):
        # [2] sub/ -> navigate; in sub: [1] .. , [2] delta.pdb
        result, _ = _run(monkeypatch, ["2", "2"], tree)
        assert Path(result).name == "delta.pdb"

    def test_parent_shortcut(self, monkeypatch, tree):
        # Go into sub, '..' back up, then pick alpha (index 3 at root).
        result, _ = _run(monkeypatch, ["2", "..", "3"], tree)
        assert Path(result).name == "alpha.pdb"

    def test_invalid_then_valid(self, monkeypatch, tree):
        result, _ = _run(monkeypatch, ["99", "3"], tree)
        assert Path(result).name == "alpha.pdb"

    def test_multi_select_list(self, monkeypatch, tree):
        result, _ = _run(monkeypatch, ["3,4,5"], tree, multi=True)
        assert [Path(p).name for p in result] == ["alpha.pdb", "beta.pdb", "gamma.pdb"]

    def test_multi_select_range(self, monkeypatch, tree):
        result, _ = _run(monkeypatch, ["3-5"], tree, multi=True)
        assert [Path(p).name for p in result] == ["alpha.pdb", "beta.pdb", "gamma.pdb"]

    def test_skip_returns_sentinel(self, monkeypatch, tree):
        result, _ = _run(monkeypatch, ["skip"], tree, optional=True)
        assert result is SKIP

    def test_selecting_directory_index_navigates_not_selects(self, monkeypatch, tree):
        # Picking [2] (sub/) must navigate, not return it as a file; then pick
        # delta inside. If it had wrongly returned the dir, the second response
        # would be unused and the result would be the directory path.
        result, scripted = _run(monkeypatch, ["2", "2"], tree)
        assert Path(result).name == "delta.pdb"
        assert scripted._responses == []  # both responses consumed


# --------------------------------------------------------------------------- #
# Filename-based session replay
# --------------------------------------------------------------------------- #

class _FakeSessionManager:
    def __init__(self, replay_ctx=None):
        self._replay_ctx = replay_ctx
        self.annotations = []

    def get_last_replayed_interaction(self):
        if self._replay_ctx is None:
            return None
        return {"prompt": "\nSelection", "context": self._replay_ctx}

    def is_recording(self):
        return self._replay_ctx is None  # recording when not replaying

    def annotate_last_recorded(self, extra):
        self.annotations.append(extra)


class _FakeProcessor:
    def __init__(self, replay_ctx=None):
        self.session_manager = _FakeSessionManager(replay_ctx)


def _browse(monkeypatch, tree, processor, live_response, **kwargs):
    scripted = _ScriptedPrompt([live_response])
    monkeypatch.setattr(fb, "prompt_with_context", scripted)
    return file_browser(
        directory=str(tree),
        extensions=[".pdb"],
        console=Console(file=open("/dev/null", "w")),
        processor=processor,
        **kwargs,
    )


class TestReplay:
    def test_record_annotates_basename(self, monkeypatch, tree):
        proc = _FakeProcessor(replay_ctx=None)  # recording
        result = _browse(monkeypatch, tree, proc, "3")  # alpha.pdb
        assert Path(result).name == "alpha.pdb"
        assert proc.session_manager.annotations[-1] == {"selected_files": ["alpha.pdb"]}

    def test_replay_resolves_by_name_after_index_shift(self, monkeypatch, tree):
        # Recorded "beta.pdb". Now insert a file that sorts before it so beta's
        # index shifts. The recorded numeric response is deliberately WRONG (1).
        (tree / "aardvark.pdb").write_text("Z")
        proc = _FakeProcessor(replay_ctx={"selected_files": ["beta.pdb"]})
        result = _browse(monkeypatch, tree, proc, "1")  # stale index — ignored
        assert Path(result).name == "beta.pdb"

    def test_replay_backcompat_single_key(self, monkeypatch, tree):
        # Old recordings stored the singular 'selected_file' key.
        proc = _FakeProcessor(replay_ctx={"selected_file": "gamma.pdb"})
        result = _browse(monkeypatch, tree, proc, "1")
        assert Path(result).name == "gamma.pdb"

    def test_replay_missing_file_falls_through_to_live(self, monkeypatch, tree):
        # Recorded file no longer exists -> resolver falls through and the live
        # numeric response (3 = alpha) is used instead of crashing/mis-picking.
        proc = _FakeProcessor(replay_ctx={"selected_files": ["deleted.pdb"]})
        result = _browse(monkeypatch, tree, proc, "3")
        assert Path(result).name == "alpha.pdb"

    def test_replay_multi_by_name(self, monkeypatch, tree):
        proc = _FakeProcessor(replay_ctx={"selected_files": ["gamma.pdb", "alpha.pdb"]})
        result = _browse(monkeypatch, tree, proc, "1", multi=True)
        assert [Path(p).name for p in result] == ["gamma.pdb", "alpha.pdb"]


class TestExtensions:
    def test_path_factory_single(self, monkeypatch, tree):
        result, _ = _run(monkeypatch, ["3"], tree, path_factory=Path)
        assert isinstance(result, Path) and result.name == "alpha.pdb"

    def test_path_factory_multi(self, monkeypatch, tree):
        result, _ = _run(monkeypatch, ["3,4"], tree, multi=True, path_factory=Path)
        assert all(isinstance(p, Path) for p in result)
        assert [p.name for p in result] == ["alpha.pdb", "beta.pdb"]

    def test_extra_command_returns_value(self, monkeypatch, tree):
        sentinel = ("FOUND", "deep/x.pdb")
        calls = []

        def handler(current_dir):
            calls.append(current_dir)
            return sentinel

        result, _ = _run(
            monkeypatch, ["find"], tree,
            extra_commands={"find": ("Search recursively", handler)},
        )
        assert result is sentinel          # returned verbatim, not wrapped
        assert len(calls) == 1             # handler got the current dir

    def test_extra_command_none_keeps_browsing(self, monkeypatch, tree):
        # find returns None (nothing found) -> keep browsing -> then pick alpha.
        result, _ = _run(
            monkeypatch, ["find", "3"], tree,
            extra_commands={"find": ("Search recursively", lambda d: None)},
        )
        assert Path(result).name == "alpha.pdb"

    def test_extra_command_case_insensitive(self, monkeypatch, tree):
        result, _ = _run(
            monkeypatch, ["FIND"], tree,
            extra_commands={"find": ("Search", lambda d: "hit")},
        )
        assert result == "hit"


class TestRemapRecordedIndex:
    """The precomputed-list remap primitive used by Tier-1 index prompts."""

    def _paths(self, *names):
        return [Path("/work") / n for n in names]

    def test_not_replaying_passes_through(self):
        proc = _FakeProcessor(replay_ctx=None)  # recording
        assert fb.remap_recorded_index(proc, self._paths("a.prmtop", "b.prmtop"), "2") == "2"

    def test_replay_remaps_to_new_index(self):
        # Recorded b.prmtop at old index 1; list now has it at index 2.
        proc = _FakeProcessor(replay_ctx={"selected_file": "b.prmtop"})
        paths = self._paths("a.prmtop", "b.prmtop")
        assert fb.remap_recorded_index(proc, paths, "1") == "2"

    def test_replay_missing_file_falls_back_to_literal(self):
        proc = _FakeProcessor(replay_ctx={"selected_file": "gone.prmtop"})
        paths = self._paths("a.prmtop", "b.prmtop")
        assert fb.remap_recorded_index(proc, paths, "1") == "1"

    def test_verbatim_token_passes_through(self):
        # 'cancel'/'browse' carry no basename → returned unchanged on replay.
        proc = _FakeProcessor(replay_ctx={"option_label": "whatever"})
        assert fb.remap_recorded_index(proc, self._paths("a.prmtop"), "cancel") == "cancel"

    def test_annotate_then_remap_roundtrip(self):
        rec = _FakeProcessor(replay_ctx=None)
        paths = self._paths("x.prmtop", "y.prmtop")
        fb.annotate_selected_path(rec, paths[1])
        assert rec.session_manager.annotations[-1] == {"selected_file": "y.prmtop"}
        # Replay with that annotation, list reordered → resolves y at new index 1.
        rep = _FakeProcessor(replay_ctx={"selected_file": "y.prmtop"})
        assert fb.remap_recorded_index(rep, self._paths("y.prmtop", "x.prmtop"), "2") == "1"


class TestRemapRecordedIndexByKey:
    """The name/ID-keyed remap primitive for non-file option lists (Tier 2B)."""

    def _sets(self, *names):
        return [{"set_name": n, "extra": 1} for n in names]

    KEY = staticmethod(lambda fs: fs["set_name"])

    def test_not_replaying_passes_through(self):
        proc = _FakeProcessor(replay_ctx=None)
        items = self._sets("ff14SB", "ff19SB")
        assert fb.remap_recorded_index_by_key(proc, items, self.KEY, "2") == "2"

    def test_replay_remaps_by_key_after_reorder(self):
        # Recorded ff19SB at old index 2; now first → resolves to index 1.
        proc = _FakeProcessor(replay_ctx={"selected_key": "ff19SB"})
        items = self._sets("ff19SB", "ff14SB")
        assert fb.remap_recorded_index_by_key(proc, items, self.KEY, "2") == "1"

    def test_replay_missing_key_falls_back(self):
        proc = _FakeProcessor(replay_ctx={"selected_key": "gone"})
        items = self._sets("ff14SB", "ff19SB")
        assert fb.remap_recorded_index_by_key(proc, items, self.KEY, "1") == "1"

    def test_key_fn_raising_is_skipped(self):
        # A malformed item whose key_fn raises must not crash the resolve.
        proc = _FakeProcessor(replay_ctx={"selected_key": "ff19SB"})
        items = [None, {"set_name": "ff19SB"}]  # key_fn(None) raises
        assert fb.remap_recorded_index_by_key(proc, items, self.KEY, "9") == "2"

    def test_annotate_key_then_remap_roundtrip(self):
        rec = _FakeProcessor(replay_ctx=None)
        fb.annotate_recorded_key(rec, "ff19SB")
        assert rec.session_manager.annotations[-1] == {"selected_key": "ff19SB"}
        rep = _FakeProcessor(replay_ctx={"selected_key": "ff19SB"})
        items = self._sets("ff14SB", "ff19SB")
        assert fb.remap_recorded_index_by_key(rep, items, self.KEY, "1") == "2"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
