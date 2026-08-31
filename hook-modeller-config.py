"""PyInstaller runtime hook: configure MODELLER inside the bundle.

Runs before the entry script. Delegates to ``proprep.utils.modeller_license``,
which pre-loads a stand-in ``modeller.config`` into sys.modules so that when
``modeller/__init__`` does ``from modeller import config`` it gets the runtime
key (``KEY_MODELLER`` or ``~/.proprep/modeller_key``) and the install_dir of
the MODELLER data bundled under ``sys._MEIPASS/modeller-<ver>/``.

The real ``modeller.config`` is excluded from the bundle by proprep.spec
because it carries the build machine's license key; the stand-in is therefore
always installed in frozen mode (with 'XXXX' when no runtime key exists, in
which case the import fails gracefully and structure repair is disabled).

The inline fallback below is only reached if the proprep package itself cannot
be imported from the bundle, and mirrors the original hook behaviour.
"""

try:
    from proprep.utils.modeller_license import configure_modeller_license
    configure_modeller_license()
except Exception:  # pragma: no cover - defensive fallback for a broken bundle
    import glob as _glob
    import os as _os
    import sys as _sys
    import types as _types

    def _fallback_key():
        key = _os.environ.get('KEY_MODELLER', '').strip()
        if key:
            return key
        try:
            with open(_os.path.join(_os.path.expanduser('~'), '.proprep', 'modeller_key')) as fh:
                for line in fh:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        return line
        except OSError:
            pass
        return 'XXXX'

    _matches = _glob.glob(_os.path.join(getattr(_sys, '_MEIPASS', ''), 'modeller-*'))
    _config = _types.ModuleType('modeller.config')
    _config.license = _fallback_key()
    _config.install_dir = _matches[0] if _matches else ''
    _config.__file__ = __file__
    _sys.modules['modeller.config'] = _config
