#!/bin/bash
# ============================================================
# ProPrep Installer
# ============================================================
# This script installs ProPrep and all its dependencies into
# a self-contained conda environment.
#
# Prerequisites:
#   - Miniforge or Anaconda/Miniconda installed
#     (https://github.com/conda-forge/miniforge#download)
#   - A MODELLER license key (free for academics)
#     Register at: https://salilab.org/modeller/registration.html
#
# Usage:
#   bash install_proprep.sh
#
# After installation:
#   conda activate ProPrep
#   source $CONDA_PREFIX/amber.sh
#   proprep
# ============================================================

# Version of ProPrep to install/verify. Bump this ONE line each release;
# every install command and the post-install check below read from it.
PROPREP_VERSION="1.14.0"

show_help() {
cat << HELPTEXT
ProPrep Installation Guide
===========================

Quick Install (Recommended)
----------------------------
Run this script:

    bash install_proprep.sh

The script handles everything: creates a conda environment, installs ProPrep
with all dependencies (including AmberTools and MODELLER), and verifies the
installation.

Manual Install
--------------

1. Install Conda

   If you don't have conda, install Miniforge:
   https://github.com/conda-forge/miniforge#download

   Follow the defaults during installation. Restart your terminal when done.

2. Get a MODELLER License Key

   MODELLER is used for homology modeling and requires a free academic license.
   Register at: https://salilab.org/modeller/registration.html

   You'll receive a license key by email. Set it as an environment variable
   before installing:

       export KEY_MODELLER="your_key_here"

   To make this permanent, add the line above to your ~/.bashrc or ~/.zshrc.

3. Create Environment and Install

       conda create --name ProPrep python=3.12 -y
       conda install -n ProPrep -c mjgplab -c dacase -c salilab -c bioconda -c conda-forge proprep=${PROPREP_VERSION} -y
       conda run -n ProPrep pip install tmtools

Usage
-----
Each time you want to use ProPrep:

    conda activate ProPrep
    source \$CONDA_PREFIX/amber.sh
    proprep

For the browser-based UI (web shell), use:

    proprep-web

Troubleshooting
---------------

"conda: command not found"
    Restart your terminal after installing Miniforge, or run:
    source ~/.bashrc (Linux) or source ~/.zshrc (macOS)

MODELLER license warning
    If you skipped the MODELLER key during install, homology modeling won't
    work. Set the key and reinstall MODELLER:
        export KEY_MODELLER="your_key_here"
        conda install -n ProPrep -c salilab modeller --force-reinstall -y

FreeSASA not available
    Run: conda install -n ProPrep -c conda-forge freesasa

Permission errors on macOS
    If macOS blocks the install, open System Settings > Privacy & Security
    and allow the blocked items.
HELPTEXT
}

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    show_help
    exit 0
fi

set -e

# Old conda on macOS misreads the OS as 10.16 (the Big Sur version-compat
# shim), so any dependency requiring __osx >=11.0 looks unsatisfiable and the
# solve fails. Feed conda the TRUE product version so detection is correct.
# Safe by construction: on a genuinely pre-11 macOS this passes the real low
# value and the solve still (correctly) fails rather than forcing a bad env.
if [[ "$(uname)" == "Darwin" && -z "$CONDA_OVERRIDE_OSX" ]]; then
    export CONDA_OVERRIDE_OSX="$(sw_vers -productVersion)"
fi

ENV_NAME="ProPrep"
PYTHON_VERSION="3.12"

echo "============================================================"
echo "  ProPrep Installer"
echo "============================================================"
echo ""

# Check for conda
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda is not installed or not in your PATH."
    echo ""
    echo "Please install Miniforge first:"
    echo "  https://github.com/conda-forge/miniforge#download"
    echo ""
    echo "After installing, restart your terminal and run this script again."
    exit 1
fi

# Check for MODELLER license key
echo "MODELLER requires a license key (free for academic use)."
echo "Register at: https://salilab.org/modeller/registration.html"
echo ""
if [ -z "$KEY_MODELLER" ]; then
    printf "Enter your MODELLER license key (or press Enter to skip): "
    modeller_key=""
    while IFS= read -r -s -n 1 char; do
        # Empty char means Enter/newline was pressed -> done
        [ -z "$char" ] && break
        if [ "$char" = $'\177' ] || [ "$char" = $'\b' ]; then
            # Backspace/Delete: drop last char and erase one star on screen
            if [ -n "$modeller_key" ]; then
                modeller_key="${modeller_key%?}"
                printf '\b \b'
            fi
        else
            modeller_key+="$char"
            printf '*'
        fi
    done
    echo ""
    if [ -n "$modeller_key" ]; then
        export KEY_MODELLER="$modeller_key"
    else
        echo ""
        echo "NOTE: Skipping MODELLER key. Homology modeling will not be"
        echo "available until you configure the key manually. See the"
        echo "install guide for instructions."
        echo ""
    fi
else
    echo "Found MODELLER key in environment variable."
fi

# Check if environment already exists
echo ""
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Environment '$ENV_NAME' already exists."
    echo ""
    echo "  1. Update ProPrep (keep existing environment)"
    echo "  2. Full reinstall (remove and recreate environment)"
    echo "  3. Cancel"
    echo ""
    # Read the choice from the controlling terminal, not stdin. Under the
    # documented `curl -fsSL ... | bash` one-liner, bash's stdin IS the piped
    # script, so a plain `read` hits EOF and silently defaults to option 1
    # (update-in-place) -- which skips the orphan-purge in option 2 and leaves
    # a stale version shadow in place. Reading from /dev/tty lets the user
    # actually choose. A truly headless run (no tty) still falls back to 1.
    if [ -r /dev/tty ]; then
        read -p "Choose [1/2/3] (1): " install_choice </dev/tty
    else
        read -p "Choose [1/2/3] (1): " install_choice
    fi
    install_choice=${install_choice:-1}

    if [[ "$install_choice" == "1" ]]; then
        echo ""
        echo "Updating ProPrep..."
        conda install -n "$ENV_NAME" -c mjgplab -c dacase -c salilab -c bioconda -c conda-forge "proprep=${PROPREP_VERSION}" -y
        conda run -n "$ENV_NAME" pip install --upgrade tmtools

        # Purge any stale pip-installed proprep metadata that would shadow the
        # conda version. A leftover proprep-*.egg-info (from an old
        # `pip install`) is returned first by importlib.metadata, so
        # `proprep --version` reports the ghost version and the assertion below
        # fails even though conda installed the right package. `conda install`
        # does not remove these; option 2 nukes the whole prefix, but the
        # update path has to clear them surgically. Targets .egg-info only so
        # conda's own proprep-<ver>.dist-info is never touched.
        ENV_PREFIX=$(conda env list | awk -v e="$ENV_NAME" '$1==e {print $NF}')
        if [ -n "$ENV_PREFIX" ] && [ -d "$ENV_PREFIX" ]; then
            rm -rf "$ENV_PREFIX"/lib/python*/site-packages/proprep-*.egg-info
        fi

        echo ""
        echo "Verifying installation..."
        VERSION=$(conda run -n "$ENV_NAME" proprep --version 2>&1 || true)
        if echo "$VERSION" | grep -qF "$PROPREP_VERSION"; then
            echo "  ProPrep updated successfully: $VERSION"
        else
            echo ""
            echo "  ERROR: expected ProPrep $PROPREP_VERSION but conda resolved:"
            echo "    ${VERSION:-<no version reported>}"
            echo "  An unmet dependency (or incorrect OS detection) likely forced"
            echo "  conda to a different version. See the conflict messages above."
            exit 1
        fi

        echo ""
        echo "============================================================"
        echo "  Update complete!"
        echo "============================================================"
        echo ""
        echo "  To use ProPrep:"
        echo ""
        echo "    conda activate $ENV_NAME"
        echo "    source \$CONDA_PREFIX/amber.sh"
        echo "    proprep              # interactive CLI"
        echo "    proprep-web          # browser-based UI"
        echo ""
        echo "============================================================"
        exit 0

    elif [[ "$install_choice" == "2" ]]; then
        # Capture the env path BEFORE removal so we can purge whatever conda
        # leaves behind. `conda env remove` only deletes conda-tracked packages;
        # pip-installed artifacts (e.g. a stale proprep .egg-info from an old
        # `pip install`) survive in site-packages and then shadow the real
        # version metadata on the next "fresh" install. Deleting the whole env
        # directory clears them so the reinstall is genuinely clean.
        ENV_PREFIX=$(conda env list | awk -v e="$ENV_NAME" '$1==e {print $NF}')
        conda env remove -n "$ENV_NAME" -y
        if [ -n "$ENV_PREFIX" ] && [ -d "$ENV_PREFIX" ]; then
            echo "Purging leftover files in $ENV_PREFIX ..."
            rm -rf "$ENV_PREFIX"
        fi
    else
        echo "Cancelled."
        exit 0
    fi
fi

# Create conda environment
echo ""
echo "Creating conda environment '$ENV_NAME' with Python $PYTHON_VERSION..."
conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y

# Install ProPrep and dependencies
echo ""
echo "Installing ProPrep and dependencies..."
echo "(This may take several minutes)"
echo ""
conda install -n "$ENV_NAME" -c mjgplab -c dacase -c salilab -c bioconda -c conda-forge "proprep=${PROPREP_VERSION}" -y

# Install PyPI-only dependencies
echo ""
echo "Installing additional Python packages..."
conda run -n "$ENV_NAME" pip install tmtools

# Verify installation
echo ""
echo "Verifying installation..."
VERSION=$(conda run -n "$ENV_NAME" proprep --version 2>&1 || true)
if echo "$VERSION" | grep -qF "$PROPREP_VERSION"; then
    echo "  ProPrep installed successfully: $VERSION"
else
    echo ""
    echo "  ERROR: expected ProPrep $PROPREP_VERSION but conda resolved:"
    echo "    ${VERSION:-<no version reported>}"
    echo "  An unmet dependency (or incorrect OS detection) likely forced"
    echo "  conda to a different version. See the conflict messages above."
    echo ""
    echo "  On macOS, if conda misreports your OS version (e.g. 10.16), try:"
    echo "    SYSTEM_VERSION_COMPAT=0 bash install_proprep.sh"
    exit 1
fi

echo ""
echo "============================================================"
echo "  Installation complete!"
echo "============================================================"
echo ""
echo "  To use ProPrep:"
echo ""
echo "    conda activate $ENV_NAME"
echo "    source \$CONDA_PREFIX/amber.sh"
echo "    proprep              # interactive CLI"
echo "    proprep-web          # browser-based UI"
echo ""
if [ -z "$KEY_MODELLER" ]; then
    echo "  To configure MODELLER later:"
    echo "    1. Get a key at https://salilab.org/modeller/registration.html"
    echo "    2. Add this to your ~/.bashrc or ~/.zshrc:"
    echo "       export KEY_MODELLER=\"your_key_here\""
    echo ""
fi
echo "============================================================"
