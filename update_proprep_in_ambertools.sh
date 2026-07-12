#!/bin/bash
# ============================================================
# Update ProPrep inside an existing (conda) AmberTools env
# ============================================================
# For users who installed AmberTools via conda (dacase::ambertools-dac)
# and want the current ProPrep in that SAME environment, rather than a
# separate ProPrep env.
#
# What it does:
#   1. Installs the standalone ProPrep conda package into your AmberTools env.
#   2. Force-reinstalls it so its files win over the older copy AmberTools
#      bundles (the two packages ship the same paths).
#   3. Deletes the stale bundled metadata that would otherwise make
#      `proprep --version` misreport the old version.
#   4. Verifies the result.
#
# Usage:
#   bash update_proprep_in_ambertools.sh [ENV_NAME]
#
# If ENV_NAME is omitted, it uses the active env, or asks.
# ============================================================

set -euo pipefail

PROPREP_VERSION="1.14.0"
CHANNELS="-c mjgplab -c dacase -c salilab -c bioconda -c conda-forge"

# --- Pick the target environment ---------------------------------
ENV_NAME="${1:-${CONDA_DEFAULT_ENV:-}}"
if [ -z "$ENV_NAME" ] || [ "$ENV_NAME" = "base" ]; then
    echo "Which conda environment contains your AmberTools install?"
    read -r -p "Environment name: " ENV_NAME </dev/tty
fi

# --- macOS: feed conda the true OS version (harmless elsewhere) --
if [ "$(uname)" = "Darwin" ] && [ -z "${CONDA_OVERRIDE_OSX:-}" ]; then
    export CONDA_OVERRIDE_OSX="$(sw_vers -productVersion)"
fi

echo "============================================================"
echo "  Updating ProPrep to ${PROPREP_VERSION} in env: ${ENV_NAME}"
echo "============================================================"

# --- Guard: the env must exist -----------------------------------
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "ERROR: conda env '${ENV_NAME}' not found. Run 'conda env list' to see your envs." >&2
    exit 1
fi
ENV_PREFIX=$(conda env list | awk -v e="$ENV_NAME" '$1==e {print $NF}')

# --- Guard: this must be a CONDA AmberTools env (Case A) ----------
# If AmberTools was built from source or installed another way there is no
# 'ambertools-dac' conda package, and installing ProPrep here would pull in a
# second ~600 MB copy of AmberTools. Refuse rather than do that.
if ! conda list -n "$ENV_NAME" 2>/dev/null | grep -qi '^ambertools-dac'; then
    echo "" >&2
    echo "This env has no conda 'ambertools-dac' package, so AmberTools was probably" >&2
    echo "built from source or installed another way. This script is only for conda" >&2
    echo "AmberTools (dacase::ambertools-dac). Running it here would install a second," >&2
    echo "parallel copy of AmberTools. Stopping. (Ask for the source-build recipe.)"  >&2
    exit 1
fi

# --- 1. Install ProPrep (+ its extra deps) into the env ----------
echo ""
echo "[1/3] Installing ProPrep ${PROPREP_VERSION} ..."
conda install -n "$ENV_NAME" $CHANNELS "proprep=${PROPREP_VERSION}" -y

# --- 2. Force ProPrep's files to win over the bundled copy --------
echo ""
echo "[2/3] Making ProPrep authoritative over the AmberTools-bundled copy ..."
conda install -n "$ENV_NAME" $CHANNELS "proprep=${PROPREP_VERSION}" --force-reinstall -y

# --- 3. Remove the stale bundled *.egg-info (version shadow) ------
echo ""
echo "[3/3] Removing stale bundled metadata ..."
find "$ENV_PREFIX"/lib/python*/site-packages -maxdepth 1 \
     -name 'proprep-[0-9]*.egg-info' -exec rm -rf {} + 2>/dev/null || true

# --- Verify ------------------------------------------------------
echo ""
echo "Verifying ..."
REPORTED=$(conda run -n "$ENV_NAME" proprep --version 2>&1 | head -1 || true)
if echo "$REPORTED" | grep -qF "$PROPREP_VERSION"; then
    echo ""
    echo "  Success: ${REPORTED}"
    echo "  ProPrep ${PROPREP_VERSION} is now active in the '${ENV_NAME}' env."
    echo ""
    echo "  Use it with:"
    echo "    conda activate ${ENV_NAME}"
    echo "    proprep"
else
    echo ""
    echo "  WARNING: expected ${PROPREP_VERSION} but got:" >&2
    echo "    ${REPORTED:-<no version reported>}" >&2
    echo "  Check: conda list -n ${ENV_NAME} proprep   (should read 1.14.0 / mjgplab)" >&2
    exit 1
fi
