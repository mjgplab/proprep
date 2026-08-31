#!/bin/bash
# ============================================================
# Export the "clean" public snapshot of a ProPrep release
# ============================================================
# Produces the tree that is published to github.com/mjgplab/proprep and
# uploaded to Zenodo (see docs/RELEASE_PROCEDURE.md, steps 5 to 7):
#
#   the files tracked at tag vX.Y.Z
#   minus  docs/  examples/  prototypes/  tools/
#   with   README.public.md renamed to README.md
#
# These rules were derived file-for-file from the 1.16.0 mirror and this
# script is checked against that listing in tests/test_release_scripts.py.
#
#   bash scripts/export_public_snapshot.sh vX.Y.Z [OUTDIR]
#
# Writes OUTDIR/proprep-vX.Y.Z-clean/ and OUTDIR/proprep-vX.Y.Z-clean.zip
# (OUTDIR defaults to releases/public/, which is gitignored).
# ============================================================
set -euo pipefail

TAG="${1:?usage: export_public_snapshot.sh vX.Y.Z [outdir]}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${2:-$ROOT/releases/public}"
VER="${TAG#v}"
NAME="proprep-${TAG}-clean"
SNAP="$OUT/$NAME"

case "$VER" in
    [0-9]*.[0-9]*.[0-9]*) ;;
    *) echo "error: '$TAG' does not look like vX.Y.Z" >&2; exit 2 ;;
esac
if ! git -C "$ROOT" rev-parse -q --verify "refs/tags/$TAG" >/dev/null; then
    echo "error: tag $TAG does not exist in $ROOT (tag the release first)" >&2
    exit 2
fi

rm -rf "$SNAP"
mkdir -p "$SNAP"
# git archive exports tracked files only (never gitignored ones such as
# CLAUDE.md, build/, dist/, leap.log).
git -C "$ROOT" archive --format=tar "$TAG" | tar -x -C "$SNAP"

# Exclusions: developer-only material that is not part of the public release.
rm -rf "$SNAP/docs" "$SNAP/examples" "$SNAP/prototypes" "$SNAP/tools"
# The public README is the one users see on GitHub.
mv "$SNAP/README.public.md" "$SNAP/README.md"

# Sanity: the exported tree must pin the version its tag claims.
if ! grep -q "^version = \"$VER\"" "$SNAP/pyproject.toml"; then
    echo "error: pyproject.toml at $TAG does not say version $VER" >&2
    exit 1
fi

rm -f "$OUT/$NAME.zip"
( cd "$OUT" && zip -qr "$NAME.zip" "$NAME" )

echo "tree: $SNAP  ($(find "$SNAP" -type f | wc -l | tr -d ' ') files)"
echo "zip:  $OUT/$NAME.zip"
