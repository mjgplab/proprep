"""
CLI: resolve atom-type collisions among a selection of FF sets.

    python -m proprep.ff_compat.resolve \\
        --select <set_id> --select <set_id> [--select ...] \\
        --workspace <path> \\
        [--rename <set_id_to_rename>] [--rename ...] \\
        [--no-prompt]

For every pair of conflicting selections, the user must choose which side
to rename. Two ways to specify:

    --rename <set_id>     # this set will be renamed wherever it appears in a conflicting pair
    (interactive)         # if no --rename covers a conflicting pair, prompt the user

With ``--no-prompt`` the tool errors out instead of prompting, suitable
for batch use / CI.

Outputs land at ``<workspace>/.ff_resolved/<selection_hash>/``; the
manifest path is printed on stdout.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from proprep.ff_compat.resolver import resolve


def _interactive_prompt(a_id: str, b_id: str) -> str:
    print(f"\nConflict between:\n  [a] {a_id}\n  [b] {b_id}")
    while True:
        choice = input("Which set's atom types should be renamed? [a/b]: ").strip().lower()
        if choice in ("a", "1"):
            return a_id
        if choice in ("b", "2"):
            return b_id
        print("Please type 'a' or 'b'.")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Resolve FF atom-type collisions.")
    parser.add_argument("--select", action="append", required=True, metavar="SET_ID",
                        help="A set_id (cofactor_path:redox_state:spin_state:set_name) "
                             "to include in the selection. Specify at least twice.")
    parser.add_argument("--rename", action="append", default=[], metavar="SET_ID",
                        help="A set_id that should be renamed wherever it appears in a "
                             "conflicting pair. Repeat for multiple sets.")
    parser.add_argument("--workspace", type=Path, required=True,
                        help="Directory under which .ff_resolved/ is created.")
    parser.add_argument("--no-prompt", action="store_true",
                        help="Error out when a conflicting pair has no rename target "
                             "specified via --rename. Default: prompt interactively.")
    args = parser.parse_args(argv)

    if len(args.select) < 2:
        parser.error("--select must be given at least twice (two sets to compare)")

    # Pre-fill rename_choices from --rename flags. For each pair (sorted), if
    # either side is in --rename, that side wins. If both are listed, the
    # earlier one (alphabetical) wins. We populate the entire choice map first.
    from proprep.ff_compat.matrix import compare_signatures
    from proprep.ff_compat.parser import parse_set
    from proprep.ff_compat.resolver import _index_descriptors

    sorted_ids = sorted(set(args.select))
    rename_set = set(args.rename)

    # Pre-resolve the conflict map so the user sees what's about to happen.
    descriptors = _index_descriptors(sorted_ids)
    sigs = {sid: parse_set(sid, descriptors[sid]["frcmod"], descriptors[sid]["libs"])
            for sid in sorted_ids}
    rename_choices = {}
    for i, a in enumerate(sorted_ids):
        for b in sorted_ids[i + 1:]:
            report = compare_signatures(sigs[a], sigs[b])
            if not report:
                continue
            pair_key = f"{a}|{b}"
            chosen = None
            if a in rename_set:
                chosen = a
            elif b in rename_set:
                chosen = b
            if chosen is not None:
                rename_choices[pair_key] = chosen

    prompt = None if args.no_prompt else _interactive_prompt
    result = resolve(
        set_ids=sorted_ids,
        rename_choices=rename_choices,
        workspace=args.workspace,
        prompt_for_pair=prompt,
    )

    print(f"\nResolved selection ({len(sorted_ids)} sets):")
    for sid, out in result.set_outputs.items():
        if out.was_renamed:
            renames = ", ".join(f"{k}->{v}" for k, v in sorted(out.rename_map.items()))
            print(f"  [renamed]  {sid}")
            print(f"             {renames}")
            if out.add_atom_types_path:
                print(f"             addAtomTypes: {out.add_atom_types_path}")
        else:
            print(f"  [original] {sid}")
    if result.add_atom_types_block:
        print("\nAggregated addAtomTypes block (drop into tleap script):")
        print(result.add_atom_types_block)
    print(f"\nManifest:  {result.manifest_path}")
    print(f"Workspace: {result.workspace_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
