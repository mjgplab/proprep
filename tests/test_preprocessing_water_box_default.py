"""Regression: a water FF selected during MCPB preprocessing records box='none'
(the metal-free preprocessing build adds no box). When the production topology
later emits an explicit solvation box, that placeholder must resolve to the
model's real AMBER solvent unit, otherwise tleap gets `solvateBox mol none ...`
and fails. Covers `_default_water_box` and the box passthrough in
`_build_standard_forcefield_section`.
"""

from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator


def test_default_water_box_maps_models():
    f = TLeapInputGenerator._default_water_box
    assert f({"name": "tip3p", "leaprc": "leaprc.water.tip3p", "box": "none"}) == "TIP3PBOX"
    assert f({"name": "opc"}) == "OPCBOX"
    assert f({"name": "opc3"}) == "OPC3BOX"        # opc3 before opc
    assert f({"name": "tip4pew"}) == "TIP4PEWBOX"  # tip4pew before tip4p
    assert f({"name": "tip4p"}) == "TIP4PBOX"
    assert f({"name": "spce"}) == "SPCBOX"


def test_default_water_box_falls_back_to_tip3p():
    assert TLeapInputGenerator._default_water_box({"name": "weird"}) == "TIP3PBOX"
    assert TLeapInputGenerator._default_water_box({}) == "TIP3PBOX"


def test_build_ff_section_resolves_none_box():
    # box='none' from preprocessing must resolve to a real box; an explicit box
    # must pass through unchanged.
    gen = TLeapInputGenerator.__new__(TLeapInputGenerator)

    def fake_ws(key, default=None):
        return {
            "protein": {"name": "ff14SB", "leaprc": "leaprc.protein.ff14SB"},
            "water": {"name": "tip3p", "leaprc": "leaprc.water.tip3p", "box": "none"},
        }
    gen.get_from_workspace = fake_ws
    _, water_box = gen._build_standard_forcefield_section()
    assert water_box == "TIP3PBOX"

    def fake_ws2(key, default=None):
        return {"water": {"name": "opc", "leaprc": "leaprc.water.opc", "box": "OPCBOX"}}
    gen.get_from_workspace = fake_ws2
    _, water_box = gen._build_standard_forcefield_section()
    assert water_box == "OPCBOX"
