"""Membrane Builder: a failed tLEaP hydrogen pass is explained, and never papered over.

Seen on 6R2Q. tLEaP stopped on "Could not find bond parameter for atom types:
SH - S". The builder printed "tLEaP did not produce topology files" with no
reason, because it looked for one in stderr and tLEaP writes its errors to
stdout. It then dropped ``--notprotonate --nottrim`` without asking and let
packmol-memgen protonate the protein: the bilayer was built around a protein
with all 20 hemes removed (0 of 20 Fe atoms left) and with one more ASH and one
more GLH than the Protonation State Analyzer had assigned.
"""

import io
import types

from rich.console import Console

from proprep.membrane_prep import membrane_builder as mb
from proprep.membrane_prep import packmol_runner

# Verbatim from the 6R2Q run (paths shortened): two occurrences of one error,
# then another, as tLEaP prints them on stdout.
TLEAP_STDOUT = """\
Checking Unit.

/opt/env/bin/teLeap: Error!
Could not find bond parameter for atom types: SH - S
        for atom SG at position 70.176000, 10.988000, 5.768000
        and atom SG at position 70.511000, 11.709000, 3.888000.

/opt/env/bin/teLeap: Error!
Could not find bond parameter for atom types: SH - S
        for atom SG at position 7.333000, -2.671000, 25.631000
        and atom SG at position 6.338000, -3.396000, 27.211000.

/opt/env/bin/teLeap: Error!
Could not find angle parameter for atom types: HS - SH - S
        for atom HG at position 69.934429, 12.111198, 6.438051,
/opt/env/bin/teLeap: Warning!
Parameter file was not saved.
"""


def test_tleap_errors_are_read_from_stdout_once_each():
    assert mb.MembraneBuilderModule._tleap_error_messages(TLEAP_STDOUT) == [
        "Could not find bond parameter for atom types: SH - S",
        "Could not find angle parameter for atom types: HS - SH - S",
    ]
    assert mb.MembraneBuilderModule._tleap_error_messages("") == []
    assert mb.MembraneBuilderModule._tleap_error_messages(None) == []


def _builder_with_failed_hydrogen_pass(monkeypatch, tmp_path, answers):
    """A builder whose tLEaP pass fails; returns it, its output and the packmol calls."""
    pdb = tmp_path / "protein.pdb"
    pdb.write_text("END\n")
    module = mb.MembraneBuilderModule()
    out = io.StringIO()
    module.processor = types.SimpleNamespace(console=Console(file=out, width=300))
    module.config.protein_pdb = str(pdb)
    module.config.skip_protonation = True

    asked, packmol_calls = [], []

    def confirm(processor, prompt, **kwargs):
        asked.append((prompt.strip(), kwargs.get("default")))
        return answers.pop(0)

    def run_packmol(args, work_dir, console, **kwargs):
        packmol_calls.append(list(args))
        return packmol_runner.PackmolResult(success=False, error_message="stop here")

    monkeypatch.setattr(mb, "confirm_with_context", confirm)
    monkeypatch.setattr(module, "_show_review", lambda workspace: None)
    monkeypatch.setattr(module, "_run_pre_tleap_hydrogen_pass", lambda workspace, work_dir: None)
    monkeypatch.setattr(packmol_runner, "find_packmol_memgen", lambda: "/x/packmol-memgen")
    monkeypatch.setattr(packmol_runner, "run_packmol_memgen", run_packmol)
    return module, out, asked, packmol_calls


def test_failed_hydrogen_pass_stops_the_build_by_default(monkeypatch, tmp_path):
    module, out, asked, packmol_calls = _builder_with_failed_hydrogen_pass(
        monkeypatch, tmp_path, answers=[True, False])

    assert module._run_build(workspace=None) is False

    # Asked, with stopping as the default, and told what the alternative costs.
    assert asked[1] == ("Build anyway with packmol-memgen's own protonation?", False)
    said = out.getvalue()
    assert "removes residues it does not" in said and "protonation" in said
    assert packmol_calls == []
    assert module.config.skip_protonation is True


def test_fallback_runs_only_when_chosen(monkeypatch, tmp_path):
    module, _, _, packmol_calls = _builder_with_failed_hydrogen_pass(
        monkeypatch, tmp_path, answers=[True, True])

    module._run_build(workspace=None)

    assert len(packmol_calls) == 1
    assert "--notprotonate" not in packmol_calls[0] and "--nottrim" not in packmol_calls[0]
