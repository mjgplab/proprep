"""Shared interactive file/directory browser.

Historically ProPrep grew ~11 near-duplicate file pickers across modules
(structure loader, MD manager, force-field prep, ORCA/QM-MM prep, frame
extractor, ...). They drifted in selection syntax (``select N`` vs bare ``N``),
cancel keyword (``exit`` vs ``q``), entry coloring, help text, multi-select
grammar, and — most importantly — whether they were captured by the session
recorder at all (a few used raw ``input()``).

This module provides a single ``file_browser()`` that all of them can call so
the UX is uniform:

* **bare ``N`` selection** — a number on a *file* selects it, a number on a
  *directory* navigates into it (no ``select`` prefix);
* **``q`` cancels** (returns ``None``);
* **multi-select** with comma lists and inclusive ranges (``1-5,7,9``) when
  ``multi=True``;
* one entry format (``[NN] 📁 dir`` / ``[NN] 📄 file (detail)``) and one help
  block;
* **filename-based session record/replay** — the recorder remembers the
  *basename(s)* chosen, not the index, so files being added/removed between
  runs can't make replay pick the wrong file.

Returns ``str`` (single mode), ``list[str]`` (``multi=True``), the ``SKIP``
sentinel (``optional=True`` and the user skipped), or ``None`` (cancelled).
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Callable, List, Optional, Sequence, Union

from rich.console import Console

from proprep.utils.prompts import prompt_with_context


# Sentinel returned when ``optional=True`` and the user chooses to skip. It is
# distinct from ``None`` (cancel) so callers can tell "proceed without a file"
# apart from "abort the whole operation".
class _Skip:
    __slots__ = ()

    def __repr__(self) -> str:  # pragma: no cover - debug aid only
        return "<file_browser.SKIP>"


SKIP = _Skip()


def _parse_indices(raw: str, count: int, multi: bool) -> List[int]:
    """Parse a selection string into a list of 1-based indices.

    In ``multi`` mode tokens may be bare numbers or inclusive ranges (``2-4``),
    separated by commas and/or whitespace, e.g. ``"1-5, 7 9"``. In single mode
    only one bare integer is accepted. Returns the (de-duplicated, order-
    preserving) indices. Raises ``ValueError`` on any malformed or out-of-range
    token so the caller can re-prompt without a partial selection.
    """
    tokens = [t for t in re.split(r"[,\s]+", raw.strip()) if t]
    if not tokens:
        raise ValueError("empty selection")
    if not multi and len(tokens) != 1:
        raise ValueError("single selection only")

    nums: List[int] = []
    for tok in tokens:
        m = re.fullmatch(r"(\d+)\s*-\s*(\d+)", tok) if multi else None
        if m:
            lo, hi = int(m.group(1)), int(m.group(2))
            if lo > hi:
                lo, hi = hi, lo
            span = range(lo, hi + 1)
        else:
            span = (int(tok),)  # raises ValueError on non-numeric token
        for n in span:
            if not (1 <= n <= count):
                raise ValueError(f"index {n} out of range 1-{count}")
            if n not in nums:
                nums.append(n)
    return nums


def _build_items(current_dir: str, extensions: Optional[Sequence[str]],
                 entry_detail: Optional[Callable[[str], str]]):
    """List the parent link, subdirectories and matching files in a directory.

    Returns a list of ``(display, abspath, kind)`` tuples where ``kind`` is
    ``"dir"`` or ``"file"``. Files are filtered by ``extensions`` (case-
    insensitive; ``None`` means all files). ``entry_detail`` optionally renders
    a trailing annotation (size, date, ...) for each file.
    """
    items: List[tuple] = []

    parent = os.path.dirname(current_dir)
    if parent and parent != current_dir:
        items.append((".. (parent)", parent, "dir"))

    try:
        names = sorted(os.listdir(current_dir))
    except (PermissionError, FileNotFoundError):
        return items

    exts = None
    if extensions is not None:
        exts = {e.lower() for e in extensions}

    subdirs = [n for n in names
               if not n.startswith(".") and os.path.isdir(os.path.join(current_dir, n))]
    for name in subdirs:
        items.append((name + "/", os.path.join(current_dir, name), "dir"))

    files = []
    for name in names:
        path = os.path.join(current_dir, name)
        if not os.path.isfile(path):
            continue
        if exts is not None and os.path.splitext(name)[1].lower() not in exts:
            continue
        files.append((name, path))

    for name, path in files:
        display = name
        if entry_detail is not None:
            try:
                extra = entry_detail(path)
            except Exception:
                extra = ""
            if extra:
                display = f"{name}  ({extra})"
        items.append((display, path, "file"))

    return items


def default_size_detail(path: str) -> str:
    """``entry_detail`` helper: human-readable file size."""
    try:
        size = os.path.getsize(path)
    except OSError:
        return ""
    if size < 1024:
        return f"{size} B"
    if size < 1024 * 1024:
        return f"{size / 1024:.1f} KB"
    return f"{size / (1024 * 1024):.1f} MB"


def _recorded_context(processor):
    """Return the just-replayed interaction's context dict, or None.

    ``None`` means the previous prompt was answered by live input (record path)
    rather than replayed, so there is nothing to re-resolve.
    """
    if not (processor and hasattr(processor, "session_manager")):
        return None
    replayed = processor.session_manager.get_last_replayed_interaction()
    if replayed is None:
        return None
    return replayed.get("context") or {}


def _annotate(processor, extra) -> None:
    """Attach resolver hints (chosen basenames / nav target) to the recording."""
    if not (processor and hasattr(processor, "session_manager")):
        return
    sm = processor.session_manager
    if sm.is_recording():
        sm.annotate_last_recorded(extra)


def annotate_selected_path(processor, path) -> None:
    """Record a file choice by basename on the just-recorded interaction.

    Pair with :func:`remap_recorded_index` to make a "pick one of N precomputed
    files by 1-based index" prompt replay-robust: the recorded basename lets
    replay re-resolve the file's position even if the candidate list changed.
    """
    _annotate(processor, {"selected_file": os.path.basename(str(path))})


def remap_recorded_index(processor, paths, live_response: str) -> str:
    """Replay-resolve an index prompt over a precomputed file list, by name.

    For prompts where the caller prints a numbered list of `paths` and reads a
    1-based index via ``prompt_with_context`` (then calls
    :func:`annotate_selected_path` on the pick). On a *record* run, or when the
    prompt was answered live, returns ``live_response`` unchanged. On *replay*,
    if the recorded interaction tagged a basename, returns the index string that
    points at that same file in the CURRENT `paths` (so an added/removed/reordered
    candidate can't mis-pick). If the recorded file is gone, falls back to the
    recorded literal response. Non-file verbatim answers (e.g. 'browse'/'cancel')
    carry no basename, so they pass through unchanged.
    """
    ctx = _recorded_context(processor)
    if not ctx:
        return live_response
    name = ctx.get("selected_file")
    if not name and ctx.get("selected_files"):
        name = ctx["selected_files"][0]
    if name:
        for i, p in enumerate(paths, 1):
            if os.path.basename(str(p)) == name:
                return str(i)
    return live_response


def annotate_recorded_key(processor, key) -> None:
    """Record a non-file choice by a stable semantic key.

    The keyed sibling of :func:`annotate_selected_path`, for "pick one of N
    named items by index" prompts where the items are NOT file paths (force-field
    sets, saved cluster profiles/run plans, search hits, ...). `key` is whatever
    stable identifier the caller's ``key_fn`` returns for the chosen item
    (a set name, profile name, accession, ...).
    """
    _annotate(processor, {"selected_key": str(key)})


def remap_recorded_index_by_key(processor, items, key_fn, live_response: str) -> str:
    """Replay-resolve an index prompt over a list of non-file items, by key.

    Like :func:`remap_recorded_index` but for lists whose elements are not file
    paths — `key_fn(item)` extracts a stable identifier (e.g. ``lambda fs:
    fs['set_name']``). On replay, if the recorded interaction tagged a
    ``selected_key``, returns the index string of the item whose key matches in
    the CURRENT `items` list (so the named item is re-found even if the list
    changed). Passes ``live_response`` through on record runs, when no key was
    recorded, or when the recorded item is gone.
    """
    ctx = _recorded_context(processor)
    if not ctx:
        return live_response
    key = ctx.get("selected_key")
    if key is None:
        return live_response
    for i, item in enumerate(items, 1):
        try:
            if str(key_fn(item)) == key:
                return str(i)
        except Exception:
            continue
    return live_response


def file_browser(
    directory: str = ".",
    extensions: Optional[Sequence[str]] = None,
    *,
    console: Optional[Console] = None,
    processor=None,
    multi: bool = False,
    label: str = "file",
    entry_detail: Optional[Callable[[str], str]] = None,
    allow_path_jump: bool = False,
    optional: bool = False,
    path_factory: Optional[Callable[[str], object]] = None,
    extra_commands: Optional[dict] = None,
    module: str = "File Browser",
):
    """Interactively browse for one or more files.

    Args:
        directory: Starting directory.
        extensions: Allowed file extensions (e.g. ``[".pdb"]``); ``None`` = all.
        console: Rich console (created if omitted).
        processor: ProPrep processor, for session recording + prompt context.
        multi: Allow selecting several files at once (``N,M`` and ``N-M``).
            Changes the return type to ``list``.
        label: Noun used in prompts/messages (e.g. ``"PDB file"``).
        entry_detail: Optional ``callable(path) -> str`` annotation per file
            (size, date, ...). See :func:`default_size_detail`.
        allow_path_jump: Accept a typed path — a file selects it, a directory
            navigates to it.
        optional: Offer a ``skip`` command that returns :data:`SKIP`.
        path_factory: Optional ``callable(str) -> object`` applied to each
            selected path before returning (e.g. ``pathlib.Path``). Defaults to
            returning the path string unchanged. Not applied to values returned
            by ``extra_commands`` handlers.
        extra_commands: Optional ``{keyword: (help_text, handler)}`` of custom
            commands. ``keyword`` is matched case-insensitively against the
            whole typed command. ``handler(current_dir_str)`` returns a value to
            return verbatim from the browser (short-circuiting selection), or
            ``None`` to keep browsing. Used for e.g. recursive ``find``.
        module: Module label for session-recording context.

    Returns:
        - ``multi=False``: the selected path (``str`` or ``path_factory`` type),
          or ``None`` if cancelled.
        - ``multi=True``: ``list`` of selected paths, or ``None``.
        - :data:`SKIP` if ``optional`` and the user skipped.
        - whatever an ``extra_commands`` handler returns.
    """
    wrap = path_factory or (lambda p: p)
    extra_commands = extra_commands or {}
    if console is None:
        console = Console()

    current_dir = os.path.abspath(directory)

    while True:
        items = _build_items(current_dir, extensions, entry_detail)
        file_idx = {  # basename -> 1-based position, for filename-based replay
            os.path.basename(p): i
            for i, (_d, p, kind) in enumerate(items, 1) if kind == "file"
        }

        console.print(f"\n[bold]Current directory:[/bold] [cyan]{current_dir}[/cyan]")
        if items:
            for i, (display, _path, kind) in enumerate(items, 1):
                icon = "📁" if kind == "dir" else "📄"
                console.print(f"  [{i:2}] {icon} {display}", highlight=False)
        else:
            console.print(f"  [grey50]No directories or {label}s here[/grey50]")

        # One unified, context-sensitive help block.
        console.print("\n[bold]Commands:[/bold]")
        if multi:
            console.print("  [cyan]N[/cyan]          - open dir N, or select file N")
            console.print("  [cyan]N,M  N-M[/cyan]   - select multiple files (e.g. 1-3,5)")
        else:
            console.print("  [cyan]N[/cyan]          - open dir N, or select file N")
        console.print("  [cyan]..[/cyan]         - go up one level")
        if allow_path_jump:
            console.print("  [cyan]<path>[/cyan]     - jump to a file or directory path")
        for keyword, (help_text, _handler) in extra_commands.items():
            console.print(f"  [cyan]{keyword}[/cyan]".ljust(24) + f" - {help_text}")
        if optional:
            console.print(f"  [cyan]skip[/cyan]       - continue without a {label}")
        console.print("  [cyan]q[/cyan]          - cancel")

        options_map = {str(i): d for i, (d, _p, _k) in enumerate(items, 1)}
        options_map["q"] = "Cancel"
        for keyword, (help_text, _handler) in extra_commands.items():
            options_map[keyword] = help_text
        if optional:
            options_map["skip"] = f"Skip ({label})"

        raw = (prompt_with_context(
            processor,
            "\nSelection",
            default="q",
            module=module,
            description=f"Browse for {label}",
            options_map=options_map,
        ) or "").strip()

        # --- Session replay: re-resolve by filename, ignoring stale indices ---
        ctx = _recorded_context(processor)
        if ctx:
            remap = _resolve_from_recording(ctx, items, file_idx, multi, console, label)
            if remap is _FALL_THROUGH:
                pass  # recorded target gone; fall through to live handling below
            elif remap is _NAVIGATE:
                current_dir = ctx.get("_nav_dir_resolved")  # set by resolver
                continue
            else:
                # resolved file path(s) — apply path_factory like a live pick
                return [wrap(p) for p in remap] if isinstance(remap, list) else wrap(remap)

        cmd = raw.lower()

        if cmd == "q":
            return None
        if optional and cmd == "skip":
            return SKIP
        if cmd in extra_commands:
            result = extra_commands[cmd][1](current_dir)
            if result is not None:
                return result  # handler's value is returned verbatim (no wrap)
            continue
        if raw == "..":
            parent = os.path.dirname(current_dir)
            if parent and parent != current_dir:
                _annotate(processor, {"nav_target": ".."})
                current_dir = parent
            else:
                console.print("[yellow]Already at the filesystem root.[/yellow]")
            continue

        # Typed path jump (opt-in).
        if allow_path_jump and (os.sep in raw or raw.startswith("~")):
            target = os.path.abspath(os.path.expanduser(raw))
            if os.path.isfile(target):
                if extensions is not None and \
                        os.path.splitext(target)[1].lower() not in {e.lower() for e in extensions}:
                    console.print(f"[yellow]{Path(target).name} is not a {label}.[/yellow]")
                    continue
                _annotate(processor, {"selected_files": [os.path.basename(target)]})
                return [wrap(target)] if multi else wrap(target)
            if os.path.isdir(target):
                _annotate(processor, {"nav_target": os.path.basename(target.rstrip(os.sep)) or target})
                current_dir = target
                continue
            console.print(f"[yellow]Path not found: {raw}[/yellow]")
            continue

        # Numeric selection / navigation.
        try:
            indices = _parse_indices(raw, len(items), multi)
        except ValueError:
            usage = ("[red]Enter a number, N,M / N-M for multiple, '..', or 'q'.[/red]"
                     if multi else
                     "[red]Enter a single number, '..', or 'q'.[/red]")
            console.print(usage)
            continue

        # A lone index on a directory means "navigate into it".
        if len(indices) == 1:
            _display, path, kind = items[indices[0] - 1]
            if kind == "dir":
                _annotate(processor, {"nav_target": os.path.basename(path.rstrip(os.sep)) or ".."})
                current_dir = path
                continue

        # Otherwise every index must point at a file.
        chosen: List[str] = []
        bad = False
        for n in indices:
            _display, path, kind = items[n - 1]
            if kind != "file":
                console.print(f"[red]Item {n} is a directory, not a {label}.[/red]")
                bad = True
                break
            if path not in chosen:
                chosen.append(path)
        if bad or not chosen:
            continue

        for p in chosen:
            console.print(f"[green]Selected: {Path(p).name}[/green]")
        _annotate(processor, {"selected_files": [os.path.basename(p) for p in chosen]})
        return [wrap(p) for p in chosen] if multi else wrap(chosen[0])


# Resolver outcomes for the replay path.
_FALL_THROUGH = object()
_NAVIGATE = object()


def _resolve_from_recording(ctx, items, file_idx, multi, console, label):
    """Re-resolve a replayed selection against the *current* directory listing.

    Returns the resolved file path(s) for a file selection, ``_NAVIGATE`` (with
    ``ctx['_nav_dir_resolved']`` set) for a directory move, or ``_FALL_THROUGH``
    if the recorded target no longer exists (so the live numeric handling runs).
    """
    # File selection — match by basename, never by the recorded index.
    recorded_files = ctx.get("selected_files")
    if recorded_files is None and ctx.get("selected_file"):
        recorded_files = [ctx["selected_file"]]  # back-compat: old single-file key
    if recorded_files:
        resolved: List[str] = []
        for name in recorded_files:
            pos = file_idx.get(name)
            if pos is None:
                console.print(
                    f"[red]Recorded {label} '{name}' is not in this directory; "
                    f"replay can't resolve it — falling through to live input.[/red]"
                )
                return _FALL_THROUGH
            resolved.append(items[pos - 1][1])
        return resolved if multi else resolved[0]

    # Directory navigation — match the recorded target name among the dirs.
    nav = ctx.get("nav_target")
    if nav:
        if nav == "..":
            ctx["_nav_dir_resolved"] = next(
                (p for d, p, k in items if k == "dir" and d.startswith("..")),
                None,
            )
            if ctx["_nav_dir_resolved"]:
                return _NAVIGATE
            return _FALL_THROUGH
        for display, path, kind in items:
            if kind == "dir" and (os.path.basename(path.rstrip(os.sep)) == nav
                                  or display.rstrip("/") == nav):
                ctx["_nav_dir_resolved"] = path
                return _NAVIGATE
        console.print(
            f"[red]Recorded directory '{nav}' is gone; falling through to live input.[/red]"
        )
        return _FALL_THROUGH

    return _FALL_THROUGH
