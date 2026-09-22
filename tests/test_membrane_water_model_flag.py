"""Membrane Builder tells packmol-memgen the water model it shows the user.

The review panel said "Water model: TIP3P" and packmol-memgen then printed
"Water model was not set. Using tip3p for ff14SB": ``--ffwat`` was never
passed, so packmol-memgen guessed from ``--ffprot``. packmol-memgen uses the
water model only in its own ``--parametrize`` step, which ProPrep does not
run, so the packed system was the same; the message was simply untrue to what
the user had chosen. Checked against the real tool when written: no message
with ``--ffwat``, and ``--ffwat fb4`` stops it with "invalid choice".
"""

import re
import subprocess

import pytest

from proprep.membrane_prep.membrane_config import PACKMOL_MEMGEN_WATER_MODELS, MembraneConfig
from proprep.membrane_prep.packmol_runner import find_packmol_memgen


def _ffwat(config):
    args = config.to_cli_args()
    return args[args.index("--ffwat") + 1] if "--ffwat" in args else None


def test_auto_water_model_is_passed_as_the_model_it_resolves_to():
    assert _ffwat(MembraneConfig(ffprot="ff14SB")) == "tip3p"
    assert _ffwat(MembraneConfig(ffprot="ff19SB")) == "opc"


def test_chosen_water_model_is_passed():
    assert _ffwat(MembraneConfig(ffprot="ff14SB", water_model="opc")) == "opc"


@pytest.mark.parametrize("model", ["fb4", "opc3pol", "tip4pd-a99SBdisp"])
def test_a_model_packmol_memgen_does_not_know_is_left_out(model):
    """It comes from the force-field selection; passing it would stop packmol-memgen."""
    config = MembraneConfig(water_model=model)
    assert config.packmol_accepts_water_model is False
    assert _ffwat(config) is None
    assert config.effective_water_model == model     # ProPrep still builds with it


def test_the_list_matches_the_installed_packmol_memgen():
    """argparse names the accepted values when given one it does not know."""
    exe = find_packmol_memgen()
    if exe is None:
        pytest.skip("packmol-memgen not installed")
    proc = subprocess.run([exe, "--ffwat", "not-a-water-model"], capture_output=True, text=True, timeout=120)
    match = re.search(r"--ffwat: invalid choice: 'not-a-water-model' \(choose from ([^)]*)\)", proc.stderr)
    assert match, proc.stderr[-500:]
    accepted = {c.strip().strip("'") for c in match.group(1).split(",")}
    assert accepted == set(PACKMOL_MEMGEN_WATER_MODELS)
