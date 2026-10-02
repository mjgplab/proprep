"""Startup notice when a newer ProPrep release is on GitHub.

At startup ProPrep asks GitHub for the latest release of mjgplab/proprep and,
when it is newer than the running version, prints a few lines under the
banner: the new version, the one update command that fits how this copy was
installed, and a link to the release page. Nothing is downloaded or changed.

Rules, in the order they are applied:

- Only copies built by the mjgplab conda recipe check at all
  (proprep._distribution.CHANNEL == "mjgplab"). The AmberTools build, source
  checkouts and installs from the public source never contact GitHub.
- The user can turn the check off in Preferences; the first time it runs it
  says so, once.
- GitHub is asked at most once per CHECK_INTERVAL; the answer is kept in
  ~/.proprep/update_check.json. At a workshop many machines share one
  address, and GitHub answers only about 60 unauthenticated requests per
  hour per address.
- The request runs in a background thread. Startup waits at most
  STARTUP_WAIT_S for it after the banner; an answer that arrives later is
  still saved and shown at the next launch. Offline, nothing is printed.
- The request carries no information about the user: it is a plain GET of
  the public releases endpoint.
- Installer users are told only once the release has their platform's
  installer attached (the installers are built after the GitHub release is
  published, docs/RELEASE_PROCEDURE.md steps 6 and 8); conda users are told
  as soon as the release is published, since the conda package goes up
  first (step 4).
"""

import json
import platform
import re
import sys
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Callable, List, Optional

REPO = "mjgplab/proprep"
LATEST_RELEASE_URL = f"https://api.github.com/repos/{REPO}/releases/latest"
CHECK_INTERVAL = timedelta(hours=24)
REQUEST_TIMEOUT_S = 5
STARTUP_WAIT_S = 1.5

# Written into the install prefix by constructor/post_install.sh, so a copy
# installed by a ProPrep-<version>-<platform>.sh installer can say so.
INSTALLER_MARKER = ".proprep_installer"
STATE_FILE = "update_check.json"
# install_proprep.sh's environment name (ENV_NAME="ProPrep").
INSTALL_SCRIPT_ENV = "ProPrep"
INSTALL_SCRIPT_CMD = "curl -fsSL https://raw.githubusercontent.com/mjgplab/proprep/main/install_proprep.sh | bash"
AMBERTOOLS_ENV_CMD = ("curl -fsSL https://raw.githubusercontent.com/mjgplab/proprep/main/"
                      "update_proprep_in_ambertools.sh | bash -s -- {env}")
PREFERENCES_HINT = "Preferences (p) > Check for new ProPrep releases at startup"


def is_our_build() -> bool:
    """True only for a copy built by the mjgplab conda recipe."""
    try:
        from proprep._distribution import CHANNEL
    except Exception:
        return False
    return CHANNEL == "mjgplab"


def parse_version(text) -> Optional[tuple]:
    """(major, minor, patch) from '1.22.0' or 'v1.22.0'; None for anything else."""
    match = re.fullmatch(r"v?(\d+)\.(\d+)\.(\d+)", str(text or "").strip())
    return tuple(int(p) for p in match.groups()) if match else None


def installer_platform(system: Optional[str] = None, machine: Optional[str] = None) -> Optional[str]:
    """The platform part of the installer file names, or None where none is built."""
    system = system or platform.system()
    machine = machine or platform.machine()
    if system == "Darwin":
        return "MacOSX-arm64" if machine == "arm64" else "MacOSX-x86_64"
    if system == "Linux" and machine in ("x86_64", "AMD64"):
        return "Linux-x86_64"
    return None


def installed_by_installer(prefix: Optional[str] = None) -> bool:
    return (Path(prefix or sys.prefix) / INSTALLER_MARKER).is_file()


def _fetch_latest_release() -> dict:
    """GitHub's latest (non-draft, non-prerelease) release as a small dict.

    Raises on a network failure. An HTTP error status returns {"status": code}
    so the caller can wait out a rate limit instead of retrying at every launch.
    """
    import requests

    response = requests.get(
        LATEST_RELEASE_URL,
        headers={"Accept": "application/vnd.github+json"},
        timeout=REQUEST_TIMEOUT_S,
    )
    if response.status_code != 200:
        return {"status": response.status_code}
    data = response.json()
    return {
        "tag": data.get("tag_name", ""),
        "url": data.get("html_url", ""),
        "assets": [a.get("name", "") for a in data.get("assets", [])],
    }


def notice_lines(release: Optional[dict], current_version: str,
                 prefix: Optional[str] = None, by_installer: Optional[bool] = None,
                 plat: Optional[str] = None) -> List[str]:
    """The lines to print for this release, or [] when there is nothing to say."""
    if not release or "tag" not in release:
        return []
    latest, current = parse_version(release["tag"]), parse_version(current_version)
    if latest is None or current is None or latest <= current:
        return []
    new = ".".join(map(str, latest))
    prefix = str(prefix or sys.prefix)
    by_installer = installed_by_installer(prefix) if by_installer is None else by_installer
    url = release.get("url") or f"https://github.com/{REPO}/releases"

    lines = [f"ProPrep {new} is available (you have {current_version})."]
    if by_installer:
        plat = plat or installer_platform()
        asset = f"ProPrep-{new}-{plat}.sh" if plat else None
        if asset is None or asset not in release.get("assets", []):
            return []  # the installer for this platform is not attached yet
        lines += [f"To update, download {asset} from the release page and run:",
                  f"  bash {asset} -b -u -p {prefix}"]
    else:
        env = Path(prefix).name
        command = INSTALL_SCRIPT_CMD if env == INSTALL_SCRIPT_ENV else AMBERTOOLS_ENV_CMD.format(env=env)
        lines += ["To update, run:", f"  {command}"]
    lines.append(f"Release page: {url}")
    return lines


class UpdateCheck:
    """One startup check. Build it with start(); show its result with show()."""

    def __init__(self, current_version: str, state_dir: Optional[Path] = None,
                 fetch: Callable[[], dict] = _fetch_latest_release,
                 now: Optional[datetime] = None, prefix: Optional[str] = None,
                 by_installer: Optional[bool] = None, plat: Optional[str] = None):
        self.current_version = current_version
        self.state_path = Path(state_dir or Path.home() / ".proprep") / STATE_FILE
        self._fetch = fetch
        self._now = now or datetime.now(timezone.utc)
        self._prefix, self._by_installer, self._plat = prefix, by_installer, plat
        self._release: Optional[dict] = None
        self._state = self._load_state()
        self._thread: Optional[threading.Thread] = None

    def _load_state(self) -> dict:
        try:
            state = json.loads(self.state_path.read_text())
            return state if isinstance(state, dict) else {}
        except Exception:
            return {}

    def _save_state(self) -> None:
        try:
            self.state_path.parent.mkdir(parents=True, exist_ok=True)
            self.state_path.write_text(json.dumps(self._state, indent=2))
        except OSError:
            pass

    def _checked_recently(self) -> bool:
        try:
            checked = datetime.fromisoformat(self._state["checked_at"])
        except (KeyError, TypeError, ValueError):
            return False
        return self._now - checked < CHECK_INTERVAL

    def _run(self) -> None:
        try:
            answer = self._fetch()
        except Exception:
            return  # offline or unreachable: say nothing, ask again next launch
        self._state["checked_at"] = self._now.isoformat()
        if "tag" in answer:
            self._state["release"] = answer
            self._release = answer
        self._save_state()

    def begin(self) -> "UpdateCheck":
        # The saved answer is shown when it is recent, and also when GitHub is
        # too slow to answer before the banner is done.
        self._release = self._state.get("release")
        if not self._checked_recently():
            self._thread = threading.Thread(target=self._run, name="proprep-update-check", daemon=True)
            self._thread.start()
        return self

    def lines(self, wait_s: float = STARTUP_WAIT_S) -> List[str]:
        if self._thread is not None:
            self._thread.join(timeout=wait_s)
        out = []
        if not self._state.get("explained"):
            out.append("ProPrep checks GitHub once a day for a newer release. "
                       f"To turn this off: {PREFERENCES_HINT}.")
            self._state["explained"] = True
            self._save_state()
        return out + notice_lines(self._release, self.current_version, self._prefix,
                                  self._by_installer, self._plat)

    def show(self, console, wait_s: float = STARTUP_WAIT_S) -> None:
        try:
            lines = self.lines(wait_s)
        except Exception:
            return
        for line in lines:
            bold = line.startswith("ProPrep ") and " is available " in line
            # soft_wrap: Rich would otherwise break a long command at the
            # terminal width, and the copied command would fail.
            console.print(line, style="bold" if bold else None, highlight=False,
                          markup=False, soft_wrap=True)


def start(current_version: str, settings=None) -> Optional[UpdateCheck]:
    """Begin the startup check, or return None when this copy must not check.

    Never raises: a failure here must not stop ProPrep from starting.
    """
    try:
        if not is_our_build():
            return None
        if settings is None:
            from proprep.utils.settings_manager import SettingsManager
            settings = SettingsManager()
        if not settings.get_update_check_enabled():
            return None
        return UpdateCheck(current_version).begin()
    except Exception:
        return None
