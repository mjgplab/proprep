"""PyInstaller runtime hook: configure MODELLER inside the bundle.

This hook runs before the main script. It pre-loads a patched version of
modeller.config into sys.modules so that when modeller.__init__ does
``from modeller import config``, it gets our version instead of the
conda-installed config.py that hard-codes 'XXXX'.

License key resolution order:
  1. KEY_MODELLER environment variable
  2. ~/.proprep/modeller_key file (first non-blank line)
  3. Falls back to 'XXXX' (import will fail gracefully; C-level error
     output is suppressed by the importing modules)

The install directory points to the MODELLER data files bundled inside
the executable (residue topologies, parameter libraries, etc.).
"""

import glob
import os
import sys
import types


def _find_bundled_modeller_dir():
    """Find the MODELLER data directory bundled inside the executable."""
    base = getattr(sys, '_MEIPASS', '')
    if base:
        pattern = os.path.join(base, 'modeller-*')
        matches = glob.glob(pattern)
        if matches:
            return matches[0]
    return ''


def _read_license_key():
    """Read the MODELLER license key from env var or ~/.proprep/modeller_key."""
    # 1. Environment variable takes precedence
    key = os.environ.get('KEY_MODELLER', '').strip()
    if key:
        return key

    # 2. Check persistent key file
    key_file = os.path.join(os.path.expanduser('~'), '.proprep', 'modeller_key')
    try:
        with open(key_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    return line
    except (OSError, IOError):
        pass

    # 3. Fall back to invalid placeholder
    return 'XXXX'


config = types.ModuleType('modeller.config')
config.license = _read_license_key()
config.install_dir = _find_bundled_modeller_dir()
config.__file__ = __file__

sys.modules['modeller.config'] = config
