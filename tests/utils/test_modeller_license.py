"""Runtime MODELLER key lookup (proprep.utils.modeller_license)."""

import sys
import types

import pytest

from proprep.utils import modeller_license as ml


@pytest.fixture
def clean_modules(monkeypatch):
    """Snapshot/restore sys.modules and purge any modeller entries."""
    saved = dict(sys.modules)
    for name in list(sys.modules):
        if name == "modeller" or name.startswith("modeller."):
            del sys.modules[name]
    yield
    sys.modules.clear()
    sys.modules.update(saved)


@pytest.fixture
def no_key(monkeypatch, tmp_path):
    monkeypatch.delenv(ml.KEY_ENV_VAR, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))  # Path.home() -> empty dir
    monkeypatch.delattr(sys, "frozen", raising=False)
    monkeypatch.delattr(sys, "_MEIPASS", raising=False)
    return tmp_path


@pytest.fixture
def fake_modeller(monkeypatch, tmp_path):
    """A stand-in 'modeller' package with a baked config.py, ahead on sys.path."""
    pkg = tmp_path / "site" / "modeller"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text(
        "from modeller import config\n"
        "LICENSE_SEEN = config.license\n"
        "INSTALL_DIR_SEEN = config.install_dir\n"
    )
    (pkg / "config.py").write_text("install_dir = r'/opt/fake/modeller-10.8'\nlicense = 'BAKED'\n")
    monkeypatch.syspath_prepend(str(tmp_path / "site"))
    return pkg


# --- key sources ---------------------------------------------------------

def test_env_var_wins_over_file(no_key, monkeypatch):
    keyfile = ml.key_file_path()
    keyfile.parent.mkdir(parents=True)
    keyfile.write_text("FILEKEY\n")
    monkeypatch.setenv(ml.KEY_ENV_VAR, "ENVKEY")
    assert ml.read_runtime_key() == "ENVKEY"
    assert ml.runtime_key_source() == "$KEY_MODELLER"


def test_key_file_skips_comments_and_blanks(no_key):
    keyfile = ml.key_file_path()
    keyfile.parent.mkdir(parents=True)
    keyfile.write_text("# my modeller key\n\n   FILEKEY  \nIGNORED\n")
    assert ml.read_runtime_key() == "FILEKEY"
    assert ml.runtime_key_source() == str(keyfile)


def test_no_sources_gives_none(no_key):
    assert ml.read_runtime_key() is None
    assert ml.runtime_key_source() is None


# --- configure_modeller_license ------------------------------------------

def test_noop_without_runtime_key_keeps_baked_config(no_key, fake_modeller, clean_modules):
    """install_proprep.sh users (key baked at install) must see no change."""
    assert ml.configure_modeller_license() == ml.STATUS_NO_RUNTIME_KEY
    assert "modeller.config" not in sys.modules
    import modeller  # noqa: F401  (the fake)
    assert modeller.LICENSE_SEEN == "BAKED"


def test_runtime_key_overrides_baked_config(no_key, fake_modeller, clean_modules, monkeypatch):
    monkeypatch.setenv(ml.KEY_ENV_VAR, "RUNTIME")
    assert ml.configure_modeller_license() == ml.STATUS_CONFIGURED
    stand_in = sys.modules["modeller.config"]
    assert stand_in.license == "RUNTIME"
    # install_dir is taken from the installed config, not invented
    assert stand_in.install_dir == "/opt/fake/modeller-10.8"
    import modeller
    assert modeller.LICENSE_SEEN == "RUNTIME"
    assert modeller.INSTALL_DIR_SEEN == "/opt/fake/modeller-10.8"


def test_idempotent(no_key, fake_modeller, clean_modules, monkeypatch):
    monkeypatch.setenv(ml.KEY_ENV_VAR, "RUNTIME")
    assert ml.configure_modeller_license() == ml.STATUS_CONFIGURED
    assert ml.configure_modeller_license() == ml.STATUS_CONFIGURED


def test_too_late_if_modeller_already_imported(no_key, fake_modeller, clean_modules, monkeypatch):
    import modeller  # noqa: F401  baked import happens first
    monkeypatch.setenv(ml.KEY_ENV_VAR, "RUNTIME")
    assert ml.configure_modeller_license() == ml.STATUS_ALREADY_IMPORTED
    assert not getattr(sys.modules.get("modeller.config"), "__proprep_runtime_key__", False)


def test_no_modeller_installed(no_key, clean_modules, monkeypatch, tmp_path):
    monkeypatch.setenv(ml.KEY_ENV_VAR, "RUNTIME")
    monkeypatch.setattr(ml.importlib.util, "find_spec", lambda name: None)
    assert ml.configure_modeller_license() == ml.STATUS_NO_MODELLER
    assert "modeller.config" not in sys.modules


def test_frozen_always_seeds_from_meipass(no_key, clean_modules, monkeypatch, tmp_path):
    """In a PyInstaller bundle the real config is excluded, so we always seed."""
    bundle = tmp_path / "bundle"
    (bundle / "modeller-10.8").mkdir(parents=True)
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setattr(sys, "_MEIPASS", str(bundle), raising=False)
    assert ml.configure_modeller_license() == ml.STATUS_CONFIGURED
    stand_in = sys.modules["modeller.config"]
    assert stand_in.license == ml.PLACEHOLDER_KEY
    assert stand_in.install_dir == str(bundle / "modeller-10.8")
    monkeypatch.setenv(ml.KEY_ENV_VAR, "RUNTIME")
    assert ml.configure_modeller_license() == ml.STATUS_CONFIGURED
    assert sys.modules["modeller.config"].license == "RUNTIME"
