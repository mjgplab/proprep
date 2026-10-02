#!/usr/bin/env python3
"""
Startup notice of a newer ProPrep release (proprep.utils.update_check).

Only copies built by the mjgplab conda recipe may check GitHub; the source
tree, and so the AmberTools build, must keep proprep._distribution.CHANNEL =
None. The recipe's own sed lines are run here against a copy of the file, so
a reworded line fails this test rather than a release build.

Run with: pytest tests/test_update_check.py
"""

import json
import re
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from proprep import _distribution  # noqa: E402
from proprep.utils import update_check  # noqa: E402
from proprep.utils.update_check import UpdateCheck, notice_lines, parse_version  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RELEASE = {
    "tag": "v1.23.0",
    "url": "https://github.com/mjgplab/proprep/releases/tag/v1.23.0",
    "assets": ["ProPrep-1.23.0-Linux-x86_64.sh", "ProPrep-1.23.0-MacOSX-arm64.sh"],
}
NOW = datetime(2026, 9, 27, 12, 0, tzinfo=timezone.utc)


# ---- the AmberTools switch ----

def test_source_tree_is_not_our_build():
    assert _distribution.CHANNEL is None
    assert update_check.is_our_build() is False
    assert update_check.start("1.0.0") is None


def test_recipe_marks_only_its_own_build(tmp_path):
    recipe = (ROOT / "recipe" / "meta.yaml").read_text()
    commands = [c for c in re.findall(r"^\s*- (.+)$", recipe, re.M)
                if "_distribution.py" in c]
    assert len(commands) == 3, commands
    (tmp_path / "src" / "proprep").mkdir(parents=True)
    target = tmp_path / "src" / "proprep" / "_distribution.py"
    target.write_text((ROOT / "src" / "proprep" / "_distribution.py").read_text())
    for command in commands:
        subprocess.run(command, shell=True, cwd=tmp_path, check=True)
    namespace = {}
    exec(target.read_text(), namespace)
    assert namespace["CHANNEL"] == "mjgplab"
    assert sorted(p.name for p in target.parent.iterdir()) == ["_distribution.py"]


def test_start_respects_the_preference(monkeypatch):
    monkeypatch.setattr(_distribution, "CHANNEL", "mjgplab")

    class Off:
        def get_update_check_enabled(self):
            return False

    assert update_check.is_our_build() is True
    assert update_check.start("1.0.0", settings=Off()) is None


# ---- what the notice says ----

def test_parse_version():
    assert parse_version("v1.22.0") == (1, 22, 0)
    assert parse_version("1.9.2") == (1, 9, 2)
    assert parse_version("unknown") is None
    assert parse_version("1.23.0rc1") is None
    assert parse_version(None) is None


def test_versions_compare_numerically_not_as_text():
    assert notice_lines({"tag": "v1.10.0", "assets": []}, "1.9.2", prefix="/x/ProPrep", by_installer=False)


@pytest.mark.parametrize("current", ["1.23.0", "1.24.0", "unknown"])
def test_nothing_to_say_when_current_or_newer(current):
    assert notice_lines(RELEASE, current, prefix="/x/ProPrep", by_installer=False) == []


def test_installer_copy_gets_the_installer_command():
    lines = notice_lines(RELEASE, "1.22.0", prefix="/Users/a/ProPrep", by_installer=True,
                         plat="MacOSX-arm64")
    assert lines[0] == "ProPrep 1.23.0 is available (you have 1.22.0)."
    assert "  bash ProPrep-1.23.0-MacOSX-arm64.sh -b -u -p /Users/a/ProPrep" in lines
    assert lines[-1] == "Release page: " + RELEASE["url"]


def test_installer_copy_waits_for_its_installer():
    # The GitHub release goes up before the installers are built and attached.
    assert notice_lines(RELEASE, "1.22.0", prefix="/p", by_installer=True, plat="MacOSX-x86_64") == []


def test_conda_env_from_the_install_script():
    lines = notice_lines(RELEASE, "1.22.0", prefix="/opt/conda/envs/ProPrep", by_installer=False)
    assert "  " + update_check.INSTALL_SCRIPT_CMD in lines


def test_conda_env_inside_ambertools():
    lines = notice_lines(RELEASE, "1.22.0", prefix="/opt/conda/envs/amber26", by_installer=False)
    assert lines[2].endswith("update_proprep_in_ambertools.sh | bash -s -- amber26")


def test_installer_platform_names_match_the_release_assets():
    assert update_check.installer_platform("Darwin", "arm64") == "MacOSX-arm64"
    assert update_check.installer_platform("Darwin", "x86_64") == "MacOSX-x86_64"
    assert update_check.installer_platform("Linux", "x86_64") == "Linux-x86_64"
    assert update_check.installer_platform("Linux", "aarch64") is None


def test_installer_marker(tmp_path):
    assert update_check.installed_by_installer(str(tmp_path)) is False
    (tmp_path / update_check.INSTALLER_MARKER).write_text("x")
    assert update_check.installed_by_installer(str(tmp_path)) is True
    post_install = (ROOT / "constructor" / "post_install.sh").read_text()
    assert f'"$PREFIX/{update_check.INSTALLER_MARKER}"' in post_install


# ---- asking GitHub at most once a day ----

def _check(tmp_path, fetch, now=NOW, current="1.22.0"):
    return UpdateCheck(current, state_dir=tmp_path, fetch=fetch, now=now,
                       prefix="/x/ProPrep", by_installer=False).begin()


def _state(tmp_path):
    return json.loads((tmp_path / update_check.STATE_FILE).read_text())


def test_first_check_fetches_saves_and_explains_once(tmp_path):
    calls = []
    lines = _check(tmp_path, lambda: calls.append(1) or RELEASE).lines()
    assert calls == [1]
    assert lines[0].startswith("ProPrep checks GitHub once a day")
    assert lines[1] == "ProPrep 1.23.0 is available (you have 1.22.0)."
    state = _state(tmp_path)
    assert state["release"] == RELEASE and state["explained"] is True

    # Next launch the same day: no request, no explanation, same notice.
    lines = _check(tmp_path, lambda: pytest.fail("fetched twice in a day"),
                   now=NOW + timedelta(hours=3)).lines()
    assert lines[0] == "ProPrep 1.23.0 is available (you have 1.22.0)."


def test_asks_again_after_a_day(tmp_path):
    _check(tmp_path, lambda: RELEASE).lines()
    calls = []
    _check(tmp_path, lambda: calls.append(1) or RELEASE, now=NOW + timedelta(hours=25)).lines()
    assert calls == [1]


def test_offline_says_nothing_and_asks_again_next_launch(tmp_path):
    def offline():
        raise OSError("no network")

    lines = _check(tmp_path, offline).lines()
    assert lines == [lines[0]] and lines[0].startswith("ProPrep checks GitHub")
    assert "checked_at" not in _state(tmp_path)


def test_rate_limited_waits_a_day_and_keeps_the_last_answer(tmp_path):
    _check(tmp_path, lambda: RELEASE).lines()
    later = NOW + timedelta(hours=25)
    lines = _check(tmp_path, lambda: {"status": 403}, now=later).lines()
    assert lines[0] == "ProPrep 1.23.0 is available (you have 1.22.0)."
    state = _state(tmp_path)
    assert state["checked_at"] == later.isoformat() and state["release"] == RELEASE


def test_slow_github_does_not_hold_startup(tmp_path):
    import threading
    release_gate = threading.Event()

    def slow():
        release_gate.wait(5)
        return RELEASE

    check = _check(tmp_path, slow)
    lines = check.lines(wait_s=0.05)
    assert lines == [lines[0]]  # only the one-time explanation
    release_gate.set()
    check._thread.join(5)
    # The late answer is saved and shown at the next launch.
    assert _state(tmp_path)["release"] == RELEASE


def test_corrupt_state_file_is_ignored(tmp_path):
    (tmp_path / update_check.STATE_FILE).write_text("{not json")
    lines = _check(tmp_path, lambda: RELEASE).lines()
    assert "ProPrep 1.23.0 is available (you have 1.22.0)." in lines


# ---- Preferences ----

def test_setting_defaults_on_and_persists(tmp_path):
    from proprep.utils.settings_manager import SettingsManager

    assert SettingsManager(tmp_path).get_update_check_enabled() is True
    SettingsManager(tmp_path).set_update_check_enabled(False)
    assert SettingsManager(tmp_path).get_update_check_enabled() is False


def test_preferences_option_is_constant_and_recorded(monkeypatch, tmp_path):
    """The label and choices must not depend on state (session replay)."""
    from rich.console import Console
    from proprep.application import menu_commands
    from proprep.utils import settings_manager

    monkeypatch.setattr(settings_manager.Path, "home", lambda: tmp_path)
    monkeypatch.setattr(_distribution, "CHANNEL", "mjgplab")

    class Proc:
        console = Console(file=open("/dev/null", "w"))

    cmd = menu_commands.PreferencesMenuCommand.__new__(menu_commands.PreferencesMenuCommand)
    cmd.processor = Proc()
    labels = []
    for enabled in (True, False):
        settings_manager.SettingsManager().set_update_check_enabled(enabled)
        cmd._show_menu()
        labels.append(cmd.options["4"])
    assert labels[0] == labels[1] == ("Check for new ProPrep releases at startup", "toggle_update_check")

    calls = []

    def fake_prompt(processor, prompt, **kw):
        calls.append((prompt, kw))
        return "2"

    monkeypatch.setattr(menu_commands, "prompt_with_context", fake_prompt)
    settings_manager.SettingsManager().set_update_check_enabled(True)
    cmd._toggle_update_check()
    assert settings_manager.SettingsManager().get_update_check_enabled() is False
    prompt, kw = calls[0]
    assert prompt == "\nCheck for new releases at startup"
    assert kw["module"] == "Preferences" and kw["description"]
    assert list(kw["options_map"]) == kw["choices"] == ["1", "2"]

    # A copy that is not ours explains and does not prompt.
    monkeypatch.setattr(_distribution, "CHANNEL", None)
    calls.clear()
    cmd._toggle_update_check()
    assert calls == []
