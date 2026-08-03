"""The restraint manager is the single authority over restraint keywords.

An imported .mdin used to be copied byte-for-byte into the run directory, which
silently discarded every restraint configured in the MD Manager's Step 3. The
simulation ran unrestrained and said nothing. Builtin templates were unaffected,
because their `template_path` is a logical string that does not resolve on disk,
so the bug only bit users who imported their own input -- the case least likely
to be noticed and most likely to matter.

These tests pin the two halves of the fix: restraints reach a custom mdin, and
overriding a hand-written value is announced rather than done quietly.
"""

import tempfile
import types
from pathlib import Path

import pytest
from rich.console import Console

from proprep.md_prep.molecular_dynamics_manager import MolecularDynamicsManager


class _Mgr(MolecularDynamicsManager):
    """Bare instance: `console` is a read-only property on the real class."""

    def __init__(self):
        self._console = Console(width=100, quiet=True)

    @property
    def console(self):
        return self._console


@pytest.fixture
def mgr():
    return _Mgr()


IMPORTED_MDIN = """Imported heating step
 &cntrl
  imin=0, ntx=1, irest=0,
  nstlim=5000, dt=0.002,
  ntr=1, restraintmask=':1-58', restraint_wt=5.0,
  ntpr=100,
 /
"""

NO_RESTRAINT_MDIN = """Imported production step
 &cntrl
  imin=0, nstlim=500000, dt=0.002, ntpr=1000,
 /
"""


def _cfg(restraints):
    return types.SimpleNamespace(restraints=restraints)


def _apply(mgr, text, cfg):
    with tempfile.TemporaryDirectory() as d:
        return mgr._apply_configured_restraints(text, cfg, Path(d))


# ── Restraints reach a custom mdin ──────────────────────────────────────

def test_restraints_are_applied_to_imported_mdin(mgr):
    out = _apply(mgr, NO_RESTRAINT_MDIN,
                 _cfg({'restraintmask': {'mask': '@CA,C,N', 'weight': 10.0}}))
    assert "restraintmask='@CA,C,N'" in out
    assert "restraint_wt=10.0" in out
    assert "ntr=1" in out


def test_manager_values_replace_the_files_own(mgr):
    """The file's mask and weight must not survive alongside the new ones."""
    out = _apply(mgr, IMPORTED_MDIN,
                 _cfg({'restraintmask': {'mask': '@CA', 'weight': 10.0}}))
    assert "@CA" in out
    assert ":1-58" not in out, "stale restraintmask survived"
    assert "5.0" not in out, "stale restraint_wt survived"


def test_no_restraints_configured_leaves_content_untouched(mgr):
    """Absent restraints, a user's input must pass through byte-identical."""
    assert _apply(mgr, IMPORTED_MDIN, _cfg(None)) == IMPORTED_MDIN
    assert _apply(mgr, IMPORTED_MDIN, _cfg({})) == IMPORTED_MDIN


# ── Overriding is announced, not silent ─────────────────────────────────

def test_override_of_hand_written_values_is_reported(mgr):
    hits = mgr._warn_if_restraints_overwrite_input(
        IMPORTED_MDIN, _cfg({'restraintmask': {'mask': '@CA', 'weight': 10.0}}))
    assert set(hits) == {"ntr", "restraintmask", "restraint_wt"}


def test_nothing_reported_when_the_input_sets_no_restraints(mgr):
    hits = mgr._warn_if_restraints_overwrite_input(
        NO_RESTRAINT_MDIN, _cfg({'restraintmask': {'mask': '@CA', 'weight': 10.0}}))
    assert hits == []


def test_nothing_reported_when_no_restraints_are_configured(mgr):
    assert mgr._warn_if_restraints_overwrite_input(IMPORTED_MDIN, _cfg(None)) == []


def test_disang_alone_does_not_claim_an_override(mgr):
    """DISANG layers on top via its own file; it does not replace ntr/restraintmask."""
    hits = mgr._warn_if_restraints_overwrite_input(
        IMPORTED_MDIN, _cfg({'disang': {'file': 'rst.disang'}}))
    assert hits == []


def test_commented_out_keywords_are_not_counted(mgr):
    text = " &cntrl\n  ntpr=100,   ! ntr=1, restraintmask=':1-58'\n /\n"
    hits = mgr._warn_if_restraints_overwrite_input(
        text, _cfg({'restraintmask': {'mask': '@CA', 'weight': 1.0}}))
    assert hits == []


def test_substring_keywords_do_not_false_match(mgr):
    """'my_ntr=' and 'nntr=' must not register as ntr."""
    text = " &cntrl\n  my_ntr=1, xrestraint_wt=3.0,\n /\n"
    hits = mgr._warn_if_restraints_overwrite_input(
        text, _cfg({'restraintmask': {'mask': '@CA', 'weight': 1.0}}))
    assert hits == []


# ── The seam itself ─────────────────────────────────────────────────────

def test_execution_path_no_longer_copies_mdin_verbatim():
    """The copy2 shortcut bypassed the restraint pass entirely.

    Guarding on source text rather than behaviour because the surrounding
    function needs a fully wired manager to call.
    """
    src = (Path(__file__).parent.parent
           / "src/proprep/md_prep/molecular_dynamics_manager.py").read_text()
    assert "shutil.copy2(sim_config.mdin_path, mdin_file)" in src, (
        "expected the guarded fallback copy to still exist"
    )
    # It must now be reachable only from the read-failure fallback, which says so.
    idx = src.index("shutil.copy2(sim_config.mdin_path, mdin_file)")
    preceding = src[max(0, idx - 600):idx]
    assert "could not read" in preceding, (
        "copy2 of an imported mdin must only happen when reading it failed"
    )
