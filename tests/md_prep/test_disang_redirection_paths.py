"""Regression: nmropt reads each restraint-file redirection line into a FIXED
80-char buffer, so 'KEYWORD=' + path must fit in 80 chars or sander silently
truncates the filename. A deep sim_dir used to clip both DUMPAVE and LISTOUT to
the same stem (they then collide). The fix emits all three redirections RELATIVE
to sim_dir, since sander is always launched with cwd=sim_dir.
"""

from pathlib import Path

from proprep.md_prep.molecular_dynamics_manager import MolecularDynamicsManager


TEMPLATE = """&cntrl
  imin=1, maxcyc=1000,
/
"""

DISANG_CONFIG = {
    "file": "/Users/someone/project/restraints.disang",
    "dumpave_file": "restraints_dump.txt",
    "listout_file": "restraints_violations.txt",
    "dump_freq": 500,
}


def _apply(sim_dir):
    mgr = MolecularDynamicsManager.__new__(MolecularDynamicsManager)
    return mgr._apply_disang_to_template(TEMPLATE, DISANG_CONFIG, Path(sim_dir))


def _redir_lines(content):
    return {
        ln.split("=", 1)[0]: ln.split("=", 1)[1]
        for ln in content.splitlines()
        if ln.startswith(("DISANG=", "DUMPAVE=", "LISTOUT="))
    }


def test_redirections_are_relative_to_sim_dir():
    deep = "/Users/someone/project/simulations/batch_20260702_162630/step1"
    redir = _redir_lines(_apply(deep))
    # DUMPAVE/LISTOUT live IN sim_dir -> bare basenames.
    assert redir["DUMPAVE"] == "restraints_dump.txt"
    assert redir["LISTOUT"] == "restraints_violations.txt"
    # DISANG is above sim_dir -> a short '../…' hop, not an absolute path.
    assert redir["DISANG"] == "../../../restraints.disang"
    assert not redir["DISANG"].startswith("/")


def test_every_redirection_line_fits_the_80_char_nmropt_buffer():
    # Pathological depth that would overflow if we emitted absolute paths.
    deep = "/Users/someone/really/deeply/nested/project/simulations/" \
           "batch_20260702_162630/step1"
    content = _apply(deep)
    for ln in content.splitlines():
        if ln.startswith(("DISANG=", "DUMPAVE=", "LISTOUT=")):
            assert len(ln) <= 80, f"redirection line exceeds nmropt buffer: {ln!r} ({len(ln)})"


def test_dumpave_and_listout_do_not_collide():
    # The original bug truncated both to the same stem.
    redir = _redir_lines(_apply("/a/b/c/step1"))
    assert redir["DUMPAVE"] != redir["LISTOUT"]


def test_redirection_block_has_no_comment_line():
    # nmropt's redirection reader has no comment support (separate prior fix).
    content = _apply("/a/b/c/step1")
    block = content.split("&wt type='END'")[-1]
    for ln in block.splitlines():
        assert not ln.strip().startswith("!"), f"comment leaked into nmropt block: {ln!r}"
