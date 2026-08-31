#!/bin/bash
# ============================================================
# ProPrep release version-lockstep check
# ============================================================
# The ProPrep version is pinned in SEVEN independent files that must
# always agree (see docs/RELEASE_PROCEDURE.md, "The version is pinned in
# SEVEN files"). The release bump is done by hand (in practice by a
# Claude session), and history shows which pins get forgotten:
#   - install_proprep.sh (missed in 1.12.0)
#   - update_proprep_in_ambertools.sh (sat at 1.14.0 through 1.16.0)
#   - constructor/construct.yaml (added 2026-08-28, two lines)
#
# Run this BEFORE `conda build` and again BEFORE the public snapshot on
# every release. It exits non-zero if any pin disagrees or is missing.
#
#   bash scripts/check_version_lockstep.sh          # check they agree
#   bash scripts/check_version_lockstep.sh 1.17.0   # also assert == 1.17.0
#   make check-version VERSION=1.17.0               # same, via Makefile
# ============================================================
set -u

# Resolve repo root from this script's location so it runs from anywhere.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 2

SEMVER='[0-9]+\.[0-9]+\.[0-9]+'

# label|file|sed program extracting the pin (targets the pin line, not any
# semver-looking string in the file).
PINS=(
  "pyproject.toml|pyproject.toml|s/^version = \"($SEMVER)\".*/\1/p"
  "setup.py|setup.py|s/.*version=\"($SEMVER)\".*/\1/p"
  "recipe/meta.yaml|recipe/meta.yaml|s/.*set version = \"($SEMVER)\".*/\1/p"
  "install_proprep.sh|install_proprep.sh|s/^PROPREP_VERSION=\"($SEMVER)\".*/\1/p"
  "update_proprep_in_ambertools.sh|update_proprep_in_ambertools.sh|s/^PROPREP_VERSION=\"($SEMVER)\".*/\1/p"
  "CITATION.cff|CITATION.cff|s/^version: ($SEMVER).*/\1/p"
  "constructor/construct.yaml (version:)|constructor/construct.yaml|s/^version: ($SEMVER).*/\1/p"
  "constructor/construct.yaml (proprep=)|constructor/construct.yaml|s/^ *- proprep=($SEMVER).*/\1/p"
)

fail=0
first=""
for entry in "${PINS[@]}"; do
    label="${entry%%|*}"
    rest="${entry#*|}"
    file="${rest%%|*}"
    prog="${rest#*|}"
    if [ ! -f "$file" ]; then
        value="<file missing>"
        fail=1
    else
        value="$(sed -nE "$prog" "$file" | head -1)"
        if [ -z "$value" ]; then
            value="<not found>"
            fail=1
        fi
    fi
    printf '%-42s %s\n' "$label" "$value"
    if [ -z "$first" ] && [ "$value" != "<not found>" ] && [ "$value" != "<file missing>" ]; then
        first="$value"
    elif [ -n "$first" ] && [ "$value" != "$first" ]; then
        fail=1
    fi
done

if [ "$fail" -ne 0 ]; then
    echo
    echo "ERROR: version pins are NOT in lockstep (or a pin was not found)."
    echo "       Bump every file listed above to the same X.Y.Z before"
    echo "       building, publishing, or exporting the public snapshot."
    echo "       See docs/RELEASE_PROCEDURE.md."
    exit 1
fi

# Optional: assert the agreed version equals an expected value.
if [ "$#" -ge 1 ]; then
    if [ "$first" != "$1" ]; then
        echo
        echo "ERROR: pins agree on $first but the release is supposed to be $1."
        exit 1
    fi
fi

echo
echo "OK: all version pins agree on $first"
