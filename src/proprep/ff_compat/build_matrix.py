"""
CLI: build the FF compatibility matrix.

    python -m proprep.ff_compat.build_matrix [--bundled-only] [--out PATH] [-v]

Default output: forcefield_params/specialized_residues/compatibility.json (bundled
matrix, committed to the repo).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from proprep.forcefield_params.loader import get_forcefield_base_path
from proprep.ff_compat.matrix import build_matrix, summarize_matrix, write_matrix


def main(argv: list | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build FF compatibility matrix.")
    parser.add_argument("--bundled-only", action="store_true",
                        help="Scan only bundled FF sets (skip ~/.proprep/forcefield_params/).")
    parser.add_argument("--out", type=Path, default=None,
                        help="Output path. Default: <bundled specialized_residues>/compatibility.json")
    parser.add_argument("--detailed", action="store_true",
                        help="Include the full conflict entry list per pair. Default "
                             "produces a summary (counts + shared types only). The full "
                             "form is large (multi-MB) and is intended for the resolver "
                             "tool; the summary form is what's committed to the repo.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print discovered sets and per-set type counts.")
    args = parser.parse_args(argv)

    matrix = build_matrix(include_user=not args.bundled_only,
                          verbose=args.verbose,
                          detailed=args.detailed)

    if args.out is None:
        out = get_forcefield_base_path() / "compatibility.json"
    else:
        out = args.out

    write_matrix(matrix, out)
    print(summarize_matrix(matrix))
    print(f"\nWrote: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
