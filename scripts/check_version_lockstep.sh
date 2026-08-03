#!/bin/bash
# ============================================================
# ProPrep release version-lockstep check
# ============================================================
# The ProPrep version is pinned in FOUR independent files that must
# always agree. The release-bump step is done by hand (in practice by
# a Claude session) and install_proprep.sh — which lives at the repo
# root, apart from the packaging files, and must also be copied to the
# separate mjgplab/proprep dist-repo — is the one most often forgotten.
#
# Run this BEFORE `conda build` on every release. It exits non-zero if
# the four pins disagree, so a missed file fails loudly instead of
# silently shipping an installer that fetches the wrong version.
#
#   bash scripts/check_version_lockstep.sh          # check they agree
#   bash scripts/check_version_lockstep.sh 1.12.0   # also assert == 1.12.0
# ============================================================
set -u

# Resolve repo root from this script's location so it runs from anywhere.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 2

# Each pin: "label|file|sed-extract-expression". The sed programs target
# the specific pin line per file, not any semver-looking string.
pyproject=$(sed -nE 's/^version = "([0-9]+\.[0-9]+\.[0-9]+)".*/\1/p'         pyproject.toml       | head -1)
setup=$(sed -nE 's/.*version="([0-9]+\.[0-9]+\.[0-9]+)".*/\1/p'             setup.py             | head -1)
recipe=$(sed -nE 's/.*set version = "([0-9]+\.[0-9]+\.[0-9]+)".*/\1/p'      recipe/meta.yaml     | head -1)
installer=$(sed -nE 's/^PROPREP_VERSION="([0-9]+\.[0-9]+\.[0-9]+)".*/\1/p'  install_proprep.sh   | head -1)

printf '%-28s %s\n' "pyproject.toml"       "${pyproject:-<not found>}"
printf '%-28s %s\n' "setup.py"             "${setup:-<not found>}"
printf '%-28s %s\n' "recipe/meta.yaml"     "${recipe:-<not found>}"
printf '%-28s %s\n' "install_proprep.sh"   "${installer:-<not found>}"

fail=0
for v in "$pyproject" "$setup" "$recipe" "$installer"; do
    [ -z "$v" ] && fail=1
done

# All four must be identical.
if [ "$pyproject" != "$setup" ] || [ "$pyproject" != "$recipe" ] || [ "$pyproject" != "$installer" ]; then
    fail=1
fi

if [ "$fail" -ne 0 ]; then
    echo
    echo "ERROR: version pins are NOT in lockstep (or a pin was not found)."
    echo "       Bump all four to the same X.Y.Z before building/releasing."
    echo "       Reminder: install_proprep.sh must ALSO be copied to the"
    echo "       mjgplab/proprep dist-repo and pushed after publishing."
    exit 1
fi

# Optional: assert the agreed version equals an expected value.
if [ "$#" -ge 1 ]; then
    if [ "$pyproject" != "$1" ]; then
        echo
        echo "ERROR: pins agree at $pyproject but you expected $1."
        exit 1
    fi
    echo
    echo "OK: all four pins agree at $pyproject (== expected $1)."
else
    echo
    echo "OK: all four version pins agree at $pyproject."
fi
