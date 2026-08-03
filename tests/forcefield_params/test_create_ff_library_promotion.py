"""Regression tests for create_ff_library after it was refactored to deposit
through the merge-safe promote_state writer.

The original implementation json.dump'd a fresh metadata.json per call, so a
second redox/spin state for the same site silently erased the first. These
tests pin the fix (states coexist), the consistent metal_sites/<type> path,
and overwrite-this-leaf-only semantics on a re-run.

The fingerprint/mol2/tleap helpers are monkeypatched so the test is hermetic
(no real MCPB file formats or tleap needed); we exercise the new staging +
delegation + merge logic, not the parsers (those are covered elsewhere).
"""

import shutil
from pathlib import Path

import pytest

from proprep.forcefield_params import loader
from proprep.forcefield_prep.mcpb import integration_utils


@pytest.fixture
def user_lib(tmp_path, monkeypatch):
    base = tmp_path / "userlib"
    monkeypatch.setattr(loader, "get_user_forcefield_base_path", lambda: base)
    return base


@pytest.fixture
def fake_mcpb(monkeypatch, tmp_path):
    """Stub the fingerprint parser, mol2 renamer, and mol2→lib conversion."""
    monkeypatch.setattr(integration_utils, "parse_fingerprint",
                        lambda fp, asn: {
                            "atom_type_entries": ['{ "M1" "Fe" "sp3" }'],
                            "residues": {(90, "HID"): []},
                        })

    def fake_rename(src, new_name, out_path):
        shutil.copy2(src, out_path)
        return out_path
    monkeypatch.setattr(integration_utils, "rename_mol2_residue", fake_rename)

    # Force the mol2 fallback (no tleap): conversion returns no lib files.
    monkeypatch.setattr(integration_utils, "_convert_mol2_to_lib",
                        lambda mol2s, frcmods, types, outdir: [])

    src = tmp_path / "src"; src.mkdir()
    mol2 = src / "HID90.mol2"; mol2.write_text("@<TRIPOS>MOLECULE\n")
    frcmod = src / "site_bonded.frcmod"; frcmod.write_text("MASS\n")
    return {"mol2": str(mol2), "frcmod": str(frcmod)}


def _call(site_type, redox, spin, fake_mcpb):
    return integration_utils.create_ff_library(
        site_type=site_type,
        description=f"{site_type} {redox}",
        mol2_files=[fake_mcpb["mol2"]],
        frcmod_files=[fake_mcpb["frcmod"]],
        fingerprint_path="ignored",
        assignments_path=None,
        residue_name_map={(90, "HID"): "HD1"},
        redox_state=redox,
        spin_state=spin,
    )


def test_deposits_under_metal_sites(user_lib, fake_mcpb):
    result = _call("zinc_his3", "oxidized", "low_spin", fake_mcpb)
    assert "/metal_sites/zinc_his3" in result["library_path"]
    assert "/mcpb/" not in result["library_path"]
    sets = loader.discover_forcefield_files("metal_sites/zinc_his3",
                                            "oxidized", "low_spin")
    assert [s["name"] for s in sets] == ["mcpb_generated"]
    # atom types are recorded for metal sites
    meta = loader.load_forcefield_metadata("metal_sites/zinc_his3")
    spin = meta["redox_states"]["oxidized"]["spin_states"]["low_spin"]
    assert spin["atom_types"] == ['{ "M1" "Fe" "sp3" }']
    assert spin["residue_name"] == {"HID90": "HD1"}


def test_second_state_does_not_clobber_first(user_lib, fake_mcpb):
    """The original bug: a 2nd state erased the 1st. Now both must survive."""
    _call("zinc_his3", "oxidized", "low_spin", fake_mcpb)
    _call("zinc_his3", "reduced", "low_spin", fake_mcpb)

    meta = loader.load_forcefield_metadata("metal_sites/zinc_his3")
    assert set(meta["redox_states"]) == {"oxidized", "reduced"}
    assert loader.discover_forcefield_files("metal_sites/zinc_his3",
                                            "oxidized", "low_spin")
    assert loader.discover_forcefield_files("metal_sites/zinc_his3",
                                            "reduced", "low_spin")


def test_rerun_same_state_overwrites_only_that_leaf(user_lib, fake_mcpb):
    _call("zinc_his3", "oxidized", "low_spin", fake_mcpb)
    _call("zinc_his3", "reduced", "low_spin", fake_mcpb)
    # Re-run oxidized with a new description → replaces just that leaf.
    integration_utils.create_ff_library(
        site_type="zinc_his3", description="oxidized REDONE",
        mol2_files=[fake_mcpb["mol2"]], frcmod_files=[fake_mcpb["frcmod"]],
        fingerprint_path="ignored", assignments_path=None,
        residue_name_map={(90, "HID"): "HD1"},
        redox_state="oxidized", spin_state="low_spin",
    )
    meta = loader.load_forcefield_metadata("metal_sites/zinc_his3")
    # reduced still there
    assert "reduced" in meta["redox_states"]
    ox_set = (meta["redox_states"]["oxidized"]["spin_states"]["low_spin"]
              ["forcefield_sets"]["mcpb_generated"])
    assert ox_set["description"] == "oxidized REDONE"
