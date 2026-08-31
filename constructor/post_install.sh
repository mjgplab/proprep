#!/bin/bash
# Runs on the attendee's machine after the environment is laid down.
# constructor exports $PREFIX (install location).

# ambertools-dac vendors an OLD proprep (1.0.0) into the same
# site-packages/proprep path as the standalone package, plus a stale
# proprep-1.0.0-py3.12.egg-info that sorts ahead of our dist-info and makes
# `proprep --version` lie. Package install order is topological (proprep
# depends on ambertools-dac, so our files land last), but the vendored
# metadata must still go. Mirrors pin_proprep() in install_proprep.sh.
rm -rf "$PREFIX"/lib/python*/site-packages/proprep-[0-9]*.egg-info

# Verify the standalone package is the one on disk: proprep.web only exists
# in the standalone releases, never in the vendored 1.0.0 copy.
if ! "$PREFIX/bin/python" -c "import proprep.web" >/dev/null 2>&1; then
    echo "WARNING: the AmberTools-vendored ProPrep appears to have overwritten the standalone package." >&2
    echo "         Re-run: $PREFIX/bin/conda install -c mjgplab proprep --force-reinstall" >&2
    exit 1
fi

"$PREFIX/bin/proprep" --version || true
