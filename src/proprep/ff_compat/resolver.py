"""
Tool 2 — collision avoidance: rewrite one or more selected FF sets with
renamed atom-type codes so a set of conflicting selections can coexist in
one tleap session.

Operates entirely on filesystem inputs and outputs; the bundled FF files
are never mutated. Rewrites land under
``<workspace>/.ff_resolved/<set_id_sanitized>/``.

Standalone CLI: ``proprep.ff_compat.resolve``. Programmatic entry point:
``resolve()`` in this module.

Scope: this resolves *atom-type* collisions (same 2-letter code declared
with different parameters). Residue-template collisions (same residue
name in two libs with different atom charges) are out of scope here; the
existing ``generate_unique_residue_names`` pipeline handles the
user-MCPB-derived path, and a follow-up will extend the resolver for the
bundled-author-variants case (e.g., Guberman vs Henriques c-heme).
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import re
import string
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

from proprep.ff_compat.matrix import discover_all_sets
from proprep.ff_compat.parser import FFSignature, parse_set


# Type-spec column widths per section. Each type is padded to 2 chars; hyphens
# separate types. So a 4-type DIHE spec is 11 chars (2*4 + 3*1).
_SECTION_TYPE_COUNT = {"MASS": 1, "NONB": 1, "BOND": 2, "ANGL": 3, "DIHE": 4, "IMPR": 4}
_SECTION_SPEC_WIDTH = {sec: n * 2 + (n - 1) for sec, n in _SECTION_TYPE_COUNT.items()}

_FRCMOD_SECTION_ALIASES = {
    "MASS": "MASS",
    "BOND": "BOND",
    "ANGL": "ANGL", "ANGLE": "ANGL",
    "DIHE": "DIHE",
    "IMPR": "IMPR",
    "NONB": "NONB", "NONBON": "NONB",
}

# Reserved literals that must never appear as a renamed atom-type.
_RESERVED_TOKENS = {"X"}


# ---------------------------------------------------------------------------
# Free-code pool
# ---------------------------------------------------------------------------

def _generate_all_codes() -> List[str]:
    """Enumerate 2-character alphanumeric codes, deterministic order.

    Pool size: 62 * 62 = 3844 codes. Tleap-compatible per the standard
    AMBER convention (alphanumeric type names of width <= 4, here width 2).
    Ordering prefers letter-letter combinations first (most natural-looking
    AMBER-style codes), then letter-digit, then digit-letter, then digit-digit.
    """
    letters = string.ascii_uppercase + string.ascii_lowercase
    digits = string.digits
    pool: List[str] = []
    for a in letters:
        for b in letters:
            pool.append(a + b)
    for a in letters:
        for b in digits:
            pool.append(a + b)
    for a in digits:
        for b in letters:
            pool.append(a + b)
    for a in digits:
        for b in digits:
            pool.append(a + b)
    return pool


def allocate_free_codes(
    types_to_rename: List[str],
    forbidden: Set[str],
    hint_keep_first_char: bool = True,
) -> Dict[str, str]:
    """Assign each type in ``types_to_rename`` a fresh 2-char code not in ``forbidden``.

    If ``hint_keep_first_char`` is True (default), the allocator prefers
    new codes whose first character matches the old type's first character —
    keeps debug output more readable (e.g., M7 -> M9, not M7 -> Aa).
    Falls back to the unrestricted pool if a same-first-char code is unavailable.

    Args:
        types_to_rename: old type names to assign new codes to.
        forbidden: codes that must not be used (already in use by any set
            in the resolution scope, plus the running renamed-codes list).
        hint_keep_first_char: prefer new codes sharing first char with old.

    Returns:
        Mapping ``{old_type: new_type}``.
    """
    pool = _generate_all_codes()
    forbidden = set(forbidden) | _RESERVED_TOKENS
    mapping: Dict[str, str] = {}

    for old in sorted(types_to_rename):
        new_type: Optional[str] = None
        if hint_keep_first_char and old:
            first = old[0]
            for candidate in pool:
                if candidate[0] == first and candidate not in forbidden and candidate != old:
                    new_type = candidate
                    break
        if new_type is None:
            for candidate in pool:
                if candidate not in forbidden and candidate != old:
                    new_type = candidate
                    break
        if new_type is None:
            raise RuntimeError(
                f"Free 2-char-code pool exhausted while renaming {old!r}; "
                f"forbidden set size = {len(forbidden)}"
            )
        mapping[old] = new_type
        forbidden.add(new_type)

    return mapping


# ---------------------------------------------------------------------------
# Frcmod rewriting
# ---------------------------------------------------------------------------

def _detect_section(line: str) -> Optional[str]:
    """Return canonical section name if ``line`` is a section header, else None."""
    stripped = line.strip()
    if not stripped:
        return None
    head = stripped.split(maxsplit=1)[0].upper()
    return _FRCMOD_SECTION_ALIASES.get(head)


def _is_data_line(line: str, section: Optional[str]) -> bool:
    """Heuristic: data lines have content after stripping and aren't headers."""
    if section is None:
        return False
    stripped = line.strip()
    if not stripped:
        return False
    if stripped.startswith(("#", "REMARK")):
        return False
    # Treat any line whose first non-space character is alphanumeric as data.
    return stripped[0].isalnum()


def _rewrite_type_spec_slice(slice_str: str, section: str,
                             rename_map: Dict[str, str]) -> str:
    """Rewrite the type-spec column substring, preserving width.

    The slice has the form e.g. ``"C -CX"`` (BOND, 5 chars). We split by '-',
    strip each token, rename, repad to 2 chars, rejoin. Result length matches
    input slice length because:
      - each type is rendered as 2 chars (padded with space on the right),
      - hyphen positions are fixed.
    """
    n = _SECTION_TYPE_COUNT[section]
    tokens = slice_str.split("-")
    if len(tokens) != n:
        # Malformed; bail out leaving the slice unchanged.
        return slice_str
    new_tokens = []
    for tok in tokens:
        old = tok.strip()
        new = rename_map.get(old, old)
        new_tokens.append(new.ljust(2))
    return "-".join(new_tokens)


def _rewrite_frcmod_data_line(line: str, section: str,
                              rename_map: Dict[str, str]) -> str:
    """Rewrite a single data line in a section, replacing renamed type tokens."""
    if section not in _SECTION_SPEC_WIDTH:
        return line

    spec_width = _SECTION_SPEC_WIDTH[section]
    # Locate the start of the type spec: skip leading whitespace.
    leading_ws_len = len(line) - len(line.lstrip(" "))
    type_spec_start = leading_ws_len
    type_spec_end = type_spec_start + spec_width

    if type_spec_end > len(line):
        # Line shorter than expected; safest to leave alone.
        return line

    old_spec = line[type_spec_start:type_spec_end]
    new_spec = _rewrite_type_spec_slice(old_spec, section, rename_map)
    if new_spec == old_spec:
        return line

    # Splice the rewritten spec back in, preserving the trailing content
    # (numeric values + comments) byte-for-byte.
    return line[:type_spec_start] + new_spec + line[type_spec_end:]


def rewrite_frcmod(src_path: Path, dst_path: Path,
                   rename_map: Dict[str, str],
                   header_comments: Optional[List[str]] = None) -> None:
    """Read frcmod from ``src_path``, apply rename_map, write to ``dst_path``.

    ``header_comments`` are emitted as ``# ...`` lines on the source title
    line (joined with ``;``) so tleap still treats line 1 as the file title
    and doesn't try to parse our metadata as parameters.
    """
    text = Path(src_path).read_text()
    lines = text.splitlines()
    out_lines: List[str] = []
    section: Optional[str] = None
    title_emitted = False

    for line in lines:
        new_section = _detect_section(line)
        if not title_emitted and new_section is None:
            # First non-section line is the frcmod title. Append a "rewrite tag"
            # so a human reading the file knows it came from the resolver.
            if header_comments:
                tag = "  [rewritten by proprep.ff_compat.resolve; see manifest]"
                out_lines.append(line + tag)
            else:
                out_lines.append(line)
            title_emitted = True
            continue
        if new_section is not None:
            section = new_section
            out_lines.append(line)
            continue
        if _is_data_line(line, section):
            out_lines.append(_rewrite_frcmod_data_line(line, section, rename_map))
        else:
            out_lines.append(line)

    dst_path.parent.mkdir(parents=True, exist_ok=True)
    dst_path.write_text("\n".join(out_lines) + ("\n" if text.endswith("\n") else ""))


# ---------------------------------------------------------------------------
# OFF library rewriting
# ---------------------------------------------------------------------------

_LIB_ATOMS_HEADER_RE = re.compile(r"^!entry\.(\S+)\.unit\.atoms\b")
_LIB_ATOMSPERTINFO_HEADER_RE = re.compile(r"^!entry\.(\S+)\.unit\.atomspertinfo\b")


def rewrite_lib(src_path: Path, dst_path: Path,
                rename_map: Dict[str, str]) -> None:
    """Rewrite the atom-type column of every atoms / atomspertinfo table.

    OFF format: in the ``!entry.<unit>.unit.atoms`` and
    ``!entry.<unit>.unit.atomspertinfo`` tables, column 2 is the atom type
    (quoted). All other columns are passed through unchanged.
    """
    text = Path(src_path).read_text()
    out_lines: List[str] = []
    in_atoms_table = False

    for line in text.splitlines():
        if _LIB_ATOMS_HEADER_RE.match(line.strip()) or _LIB_ATOMSPERTINFO_HEADER_RE.match(line.strip()):
            in_atoms_table = True
            out_lines.append(line)
            continue
        if line.strip().startswith("!"):
            in_atoms_table = False
            out_lines.append(line)
            continue
        if not in_atoms_table or not line.strip():
            out_lines.append(line)
            continue

        # Atom row: ` "atom_name" "atom_type" 0 1 131072 1 6 -0.035100`
        # Pertinfo row: ` "atom_name" "atom_type" 0 -1 0.0`
        # Rewrite by preserving the quoted-token structure exactly.
        match = re.match(r'^(\s*)"([^"]*)"\s+"([^"]*)"(\s+.*)?$', line)
        if not match:
            out_lines.append(line)
            continue
        leading, atom_name, atom_type, rest = match.groups()
        new_type = rename_map.get(atom_type, atom_type)
        if new_type == atom_type:
            out_lines.append(line)
        else:
            out_lines.append(f'{leading}"{atom_name}" "{new_type}"{rest or ""}')

    dst_path.parent.mkdir(parents=True, exist_ok=True)
    dst_path.write_text("\n".join(out_lines) + ("\n" if text.endswith("\n") else ""))


# ---------------------------------------------------------------------------
# addAtomTypes emission for renamed types
# ---------------------------------------------------------------------------

def _format_add_atom_types_entry(type_code: str, element: str, hybridization: str) -> str:
    """One tleap addAtomTypes entry line, matching the metadata.json bundled style."""
    return f'{{ "{type_code}" "{element}" "{hybridization}" }}'


def _build_add_atom_types_entries(
    rename_map: Dict[str, str],
    signature: "FFSignature",
) -> Tuple[List[str], List[str]]:
    """Produce addAtomTypes lines for every renamed code in a set.

    Returns ``(entries, warnings)`` where:
      - ``entries`` is the list of formatted lines for tleap; one per new code
        whose source type had a known element. Sorted by the new code.
      - ``warnings`` is the list of human-readable strings naming the renames
        we couldn't synthesize an entry for (no lib usage to infer element).
        These are propagated to the caller so the resolver CLI can surface
        them, and so the manifest preserves them for the picker integration.
    """
    entries: List[Tuple[str, str]] = []
    warnings: List[str] = []
    for old_type, new_type in sorted(rename_map.items(), key=lambda kv: kv[1]):
        decl = signature.infer_type_declaration(old_type)
        if decl is None:
            warnings.append(
                f"could not infer element+hybridization for {old_type!r} -> {new_type!r} "
                f"(type has no lib atom in this set)"
            )
            continue
        element, hyb = decl
        entries.append((new_type, _format_add_atom_types_entry(new_type, element, hyb)))
    entries.sort(key=lambda kv: kv[0])
    return [line for _, line in entries], warnings


# ---------------------------------------------------------------------------
# Resolution orchestration
# ---------------------------------------------------------------------------

@dataclass
class SetOutput:
    set_id: str
    was_renamed: bool
    frcmod_path: Path
    lib_paths: List[Path]
    rename_map: Dict[str, str] = field(default_factory=dict)
    # tleap addAtomTypes entry lines declaring the renamed types' element
    # + hybridization. Empty list when ``was_renamed`` is False.
    add_atom_types_entries: List[str] = field(default_factory=list)
    # Filesystem path to the per-set addatomtypes.txt sidecar (the generator
    # ingests this alongside the rewritten frcmod/lib).
    add_atom_types_path: Optional[Path] = None


@dataclass
class ResolveResult:
    workspace_dir: Path
    manifest_path: Path
    set_outputs: Dict[str, SetOutput] = field(default_factory=dict)
    # Aggregated addAtomTypes lines across every renamed set, ready to drop
    # into a tleap ``addAtomTypes { ... }`` block.
    add_atom_types_block: str = ""


def _sanitize_set_id(set_id: str) -> str:
    """Turn a set_id into a filesystem-safe directory name."""
    return set_id.replace("/", "_").replace(":", "__")


def _file_hash(path: Path) -> str:
    h = hashlib.sha256()
    h.update(Path(path).read_bytes())
    return h.hexdigest()[:12]


def _selection_hash(set_ids: List[str], descriptors: Dict[str, Dict],
                    rename_choices: Dict[str, str]) -> str:
    """Stable hash of the resolution inputs for idempotence."""
    parts: List[str] = []
    for sid in sorted(set_ids):
        d = descriptors[sid]
        files = [d["frcmod"], *d["libs"]]
        for f in sorted(files):
            parts.append(f"{sid}|{f}|{_file_hash(f)}")
    for pair, target in sorted(rename_choices.items()):
        parts.append(f"rename|{pair}|{target}")
    h = hashlib.sha256("\n".join(parts).encode()).hexdigest()
    return h[:16]


def _index_descriptors(set_ids: List[str]) -> Dict[str, Dict]:
    """Look up filesystem details for each selected set via the matrix discovery API."""
    all_descriptors = {d["set_id"]: d for d in discover_all_sets(include_user=True)}
    missing = [sid for sid in set_ids if sid not in all_descriptors]
    if missing:
        raise ValueError(f"Set IDs not found in bundled or user FF discovery: {missing}")
    return {sid: all_descriptors[sid] for sid in set_ids}


def resolve(
    set_ids: List[str],
    rename_choices: Dict[str, str],
    workspace: Path,
    prompt_for_pair: Optional[Callable[[str, str], str]] = None,
    hint_keep_first_char: bool = True,
) -> ResolveResult:
    """Resolve atom-type conflicts among ``set_ids`` by renaming one side per pair.

    Args:
        set_ids: every FF set the user is planning to load into one tleap session.
        rename_choices: maps pair_key ``"A|B"`` -> set_id to rename. When a
            conflicting pair has no entry in this dict, ``prompt_for_pair`` is
            invoked. If neither is provided for a conflicting pair, raises.
        workspace: directory under which ``.ff_resolved/`` is created.
        prompt_for_pair: optional callable ``(a_id, b_id) -> set_id_to_rename``.
            Used by the interactive CLI; in batch/test mode pass None and use
            ``rename_choices`` exclusively.
        hint_keep_first_char: see ``allocate_free_codes``.

    Returns:
        ResolveResult with paths to the resolved frcmod + lib(s) per set and
        the manifest.
    """
    if len(set_ids) < 2:
        raise ValueError("resolve() needs at least 2 set_ids to compare")

    descriptors = _index_descriptors(set_ids)

    # Parse every selected set into a signature.
    signatures: Dict[str, FFSignature] = {}
    for sid in set_ids:
        d = descriptors[sid]
        signatures[sid] = parse_set(sid, d["frcmod"], d["libs"])

    # Pairwise compare to find conflicts; resolve which side to rename for each.
    from proprep.ff_compat.matrix import compare_signatures
    sorted_ids = sorted(set_ids)
    rename_targets: Dict[str, Set[str]] = {sid: set() for sid in sorted_ids}

    for i, a_id in enumerate(sorted_ids):
        for b_id in sorted_ids[i + 1:]:
            report = compare_signatures(signatures[a_id], signatures[b_id])
            if not report:
                continue
            pair = f"{a_id}|{b_id}"
            target = rename_choices.get(pair)
            if target is None and prompt_for_pair is not None:
                target = prompt_for_pair(a_id, b_id)
            if target is None:
                raise ValueError(
                    f"Conflict between {a_id} and {b_id} has no rename target; "
                    f"provide one in rename_choices['{pair}'] or pass prompt_for_pair."
                )
            if target not in (a_id, b_id):
                raise ValueError(
                    f"rename target {target!r} for pair {pair!r} must be one of the pair members"
                )
            # Collect the conflicting types in the target side.
            for entry in report["entries"]:
                if entry["section"] == "LIB":
                    continue  # library-side conflicts handled separately
                key = entry["key"]
                if isinstance(key, str):
                    rename_targets[target].add(key)
                else:
                    for t in key:
                        if t != "X":
                            rename_targets[target].add(t)

    # Compute global forbidden set: every type already declared by any set in
    # the resolution scope, plus reserved tokens.
    all_used_types: Set[str] = set()
    for sig in signatures.values():
        all_used_types |= sig.all_types()

    # Allocate codes for each set's rename pool. We compose per-set allocations
    # sequentially so a later set's renames don't accidentally re-pick a
    # code that a prior set just claimed.
    forbidden = set(all_used_types)
    per_set_renames: Dict[str, Dict[str, str]] = {}
    for sid in sorted_ids:
        renames = rename_targets.get(sid) or set()
        if not renames:
            per_set_renames[sid] = {}
            continue
        mapping = allocate_free_codes(
            types_to_rename=list(renames),
            forbidden=forbidden,
            hint_keep_first_char=hint_keep_first_char,
        )
        forbidden |= set(mapping.values())
        per_set_renames[sid] = mapping

    # Materialize the workspace.
    selection_hash = _selection_hash(sorted_ids, descriptors, rename_choices)
    workspace = Path(workspace).expanduser().resolve()
    out_root = workspace / ".ff_resolved" / selection_hash
    out_root.mkdir(parents=True, exist_ok=True)

    set_outputs: Dict[str, SetOutput] = {}
    inference_warnings: List[str] = []
    aggregated_entries: List[str] = []
    for sid in sorted_ids:
        d = descriptors[sid]
        renames = per_set_renames[sid]
        if not renames:
            set_outputs[sid] = SetOutput(
                set_id=sid,
                was_renamed=False,
                frcmod_path=Path(d["frcmod"]),
                lib_paths=[Path(p) for p in d["libs"]],
                rename_map={},
            )
            continue

        set_dir = out_root / _sanitize_set_id(sid)
        set_dir.mkdir(parents=True, exist_ok=True)
        new_frcmod = set_dir / Path(d["frcmod"]).name
        rename_summary = ", ".join(f"{k}->{v}" for k, v in sorted(renames.items()))
        rewrite_frcmod(
            Path(d["frcmod"]), new_frcmod, renames,
            header_comments=[
                f"FF compat rewrite: set_id = {sid}",
                f"Source frcmod: {d['frcmod']}",
                f"Renamed types: {rename_summary}",
                f"Manifest: {out_root / 'manifest.json'}",
            ],
        )

        new_libs: List[Path] = []
        for lib in d["libs"]:
            new_lib = set_dir / Path(lib).name
            rewrite_lib(Path(lib), new_lib, renames)
            new_libs.append(new_lib)

        # Emit addAtomTypes lines for each renamed code.
        entries, warns = _build_add_atom_types_entries(renames, signatures[sid])
        add_path: Optional[Path] = None
        if entries:
            add_path = set_dir / "addatomtypes.txt"
            add_path.write_text(
                "# FF compat addAtomTypes for renamed codes\n"
                f"# set_id = {sid}\n"
                + "\n".join(entries) + "\n"
            )
            aggregated_entries.extend(entries)
        for w in warns:
            inference_warnings.append(f"[{sid}] {w}")

        set_outputs[sid] = SetOutput(
            set_id=sid,
            was_renamed=True,
            frcmod_path=new_frcmod,
            lib_paths=new_libs,
            rename_map=renames,
            add_atom_types_entries=entries,
            add_atom_types_path=add_path,
        )

    manifest = {
        "schema_version": 1,
        "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
        "selection_hash": selection_hash,
        "selected_set_ids": sorted_ids,
        "rename_choices": rename_choices,
        "outputs": {
            sid: {
                "was_renamed": so.was_renamed,
                "frcmod": str(so.frcmod_path),
                "libs": [str(p) for p in so.lib_paths],
                "rename_map": so.rename_map,
                "add_atom_types_entries": so.add_atom_types_entries,
                "add_atom_types_path": str(so.add_atom_types_path) if so.add_atom_types_path else None,
            }
            for sid, so in set_outputs.items()
        },
        "inference_warnings": inference_warnings,
    }
    manifest_path = out_root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))

    block = ("addAtomTypes {\n"
             + "\n".join(f"    {line}" for line in aggregated_entries)
             + "\n}") if aggregated_entries else ""

    return ResolveResult(
        workspace_dir=out_root,
        manifest_path=manifest_path,
        set_outputs=set_outputs,
        add_atom_types_block=block,
    )
