"""Tests for the user-library writer (forcefield_params/user_library.py).

These exercise the upsert/merge semantics, the collision policies, and the
round-trip-validate-or-rollback guarantee. They are hermetic: the loader and
the writer are both redirected to a tmp dir by monkeypatching
``get_user_forcefield_base_path``, and no external tools (tleap, antechamber)
are invoked — the writer copies pre-made artifact files and the loader reads
JSON + checks file existence only.
"""

import json
from pathlib import Path

import pytest

from proprep.forcefield_params import loader, user_library
from proprep.forcefield_params.user_library import (
    PromotionRequest,
    promote_state,
    LibraryCollisionError,
    UserLibraryError,
)


@pytest.fixture
def user_lib(tmp_path, monkeypatch):
    """Redirect the user library base to a tmp dir for writer AND loader."""
    base = tmp_path / "userlib"
    monkeypatch.setattr(loader, "get_user_forcefield_base_path", lambda: base)
    return base


@pytest.fixture
def artifacts(tmp_path):
    """A pair of plausible source .frcmod / .lib files in a workspace dir."""
    ws = tmp_path / "workspace"
    ws.mkdir()
    frcmod = ws / "lig.frcmod"
    frcmod.write_text("remark goes here\nMASS\n\nBOND\n")
    lib = ws / "lig.lib"
    lib.write_text("!!index array str\n \"LIG\"\n")
    return frcmod, lib


def _small_molecule_request(frcmod, lib, **overrides):
    """A minimal valid single-state (small-molecule) promotion request."""
    kwargs = dict(
        family="small_molecules",
        type_name="LIG",
        set_name="Guberman_LIG_RESP",
        frcmod_src=str(frcmod),
        lib_srcs=[str(lib)],
        residue_meta={"description": "Test ligand", "references": ["ref A"]},
        spin_meta={"residue_name": "LIG", "atom_types": []},
        set_meta={"description": "GAFF + AM1-BCC", "version": "1.0"},
    )
    kwargs.update(overrides)
    return PromotionRequest(**kwargs)


# --------------------------------------------------------------------------- #
# create
# --------------------------------------------------------------------------- #

def test_create_new_residue_is_discoverable(user_lib, artifacts):
    frcmod, lib = artifacts
    result = promote_state(_small_molecule_request(frcmod, lib))

    # Files landed at <family>/<type>/single_state/default/
    state_dir = Path(result["state_dir"])
    assert state_dir == user_lib / "small_molecules" / "LIG" / "single_state" / "default"
    assert (state_dir / "lig.frcmod").is_file()
    assert (state_dir / "lig.lib").is_file()

    # metadata.json is valid per the loader and the set is discoverable.
    meta = loader.load_forcefield_metadata("small_molecules/LIG")
    assert meta["cofactor_type"] == "LIG"
    assert meta["source"] == "user"
    sets = loader.discover_forcefield_files("small_molecules/LIG",
                                            "single_state", "default")
    assert [s["name"] for s in sets] == ["Guberman_LIG_RESP"]
    assert sets[0]["frcmod"].endswith("lig.frcmod")

    # And it shows up in the unified cross-library discovery.
    types = loader.get_available_cofactor_types()
    assert "small_molecules/LIG" in types


def test_single_lib_is_stored_as_string_not_list(user_lib, artifacts):
    frcmod, lib = artifacts
    result = promote_state(_small_molecule_request(frcmod, lib))
    meta = json.loads(Path(result["metadata_path"]).read_text())
    files = (meta["redox_states"]["single_state"]["spin_states"]["default"]
             ["forcefield_sets"]["Guberman_LIG_RESP"]["files"])
    assert files["lib"] == "lig.lib"  # bare string, not ["lig.lib"]


# --------------------------------------------------------------------------- #
# merge / upsert
# --------------------------------------------------------------------------- #

def test_add_second_state_preserves_first(user_lib, tmp_path):
    """One-state-at-a-time: a second redox state must not clobber the first."""
    ws = tmp_path / "ws"
    ws.mkdir()
    fr_ox = ws / "ox.frcmod"; fr_ox.write_text("MASS\n")
    lb_ox = ws / "ox.lib"; lb_ox.write_text("\"HEME_OX\"\n")
    fr_red = ws / "red.frcmod"; fr_red.write_text("MASS\n")
    lb_red = ws / "red.lib"; lb_red.write_text("\"HEME_RED\"\n")

    promote_state(PromotionRequest(
        family="heme", type_name="my_heme", set_name="set_ox",
        frcmod_src=str(fr_ox), lib_srcs=[str(lb_ox)],
        redox_state="oxidized", spin_state="low_spin",
        spin_meta={"residue_name": "HOX", "atom_types": []},
    ))
    promote_state(PromotionRequest(
        family="heme", type_name="my_heme", set_name="set_red",
        frcmod_src=str(fr_red), lib_srcs=[str(lb_red)],
        redox_state="reduced", spin_state="low_spin",
        spin_meta={"residue_name": "HRD", "atom_types": []},
    ))

    meta = loader.load_forcefield_metadata("heme/my_heme")
    assert set(meta["redox_states"]) == {"oxidized", "reduced"}
    assert loader.discover_forcefield_files("heme/my_heme", "oxidized", "low_spin")
    assert loader.discover_forcefield_files("heme/my_heme", "reduced", "low_spin")


def test_add_second_set_same_state(user_lib, artifacts):
    frcmod, lib = artifacts
    promote_state(_small_molecule_request(frcmod, lib))
    promote_state(_small_molecule_request(frcmod, lib, set_name="Alt_RESP"))

    sets = loader.discover_forcefield_files("small_molecules/LIG",
                                            "single_state", "default")
    names = {s["name"] for s in sets}
    assert names == {"Guberman_LIG_RESP", "Alt_RESP"}
    # Exactly one default after the second (default) set demotes the first.
    assert sum(1 for s in sets if s["is_default"]) == 1


# --------------------------------------------------------------------------- #
# collision policies
# --------------------------------------------------------------------------- #

def test_collision_error_is_default_and_leaves_state_intact(user_lib, artifacts):
    frcmod, lib = artifacts
    promote_state(_small_molecule_request(frcmod, lib))
    before = (user_lib / "small_molecules" / "LIG" / "metadata.json").read_text()

    with pytest.raises(LibraryCollisionError):
        promote_state(_small_molecule_request(frcmod, lib))  # same set_name

    after = (user_lib / "small_molecules" / "LIG" / "metadata.json").read_text()
    assert before == after  # nothing changed


def test_collision_version_bump(user_lib, artifacts):
    frcmod, lib = artifacts
    promote_state(_small_molecule_request(frcmod, lib))
    result = promote_state(_small_molecule_request(
        frcmod, lib, on_collision=user_library.ON_COLLISION_VERSION))
    assert result["set_name"] == "Guberman_LIG_RESP_v2"

    sets = loader.discover_forcefield_files("small_molecules/LIG",
                                            "single_state", "default")
    assert {s["name"] for s in sets} == {"Guberman_LIG_RESP",
                                         "Guberman_LIG_RESP_v2"}


def test_collision_overwrite(user_lib, artifacts):
    frcmod, lib = artifacts
    promote_state(_small_molecule_request(frcmod, lib))
    result = promote_state(_small_molecule_request(
        frcmod, lib, on_collision=user_library.ON_COLLISION_OVERWRITE,
        set_meta={"description": "v2 charges", "version": "2.0"}))
    assert result["set_name"] == "Guberman_LIG_RESP"

    meta = json.loads((user_lib / "small_molecules" / "LIG"
                       / "metadata.json").read_text())
    block = (meta["redox_states"]["single_state"]["spin_states"]["default"]
             ["forcefield_sets"]["Guberman_LIG_RESP"])
    assert block["version"] == "2.0"


# --------------------------------------------------------------------------- #
# failure handling
# --------------------------------------------------------------------------- #

def test_missing_source_raises(user_lib, artifacts):
    frcmod, lib = artifacts
    req = _small_molecule_request(frcmod, lib, frcmod_src="/nope/missing.frcmod")
    with pytest.raises(FileNotFoundError):
        promote_state(req)
    assert not (user_lib / "small_molecules").exists()


def test_roundtrip_failure_rolls_back_new_residue(user_lib, artifacts):
    """A leaf that won't validate (missing required residue_name) must leave
    NO trace — no metadata.json, no copied files, no empty dirs."""
    frcmod, lib = artifacts
    bad = _small_molecule_request(frcmod, lib, spin_meta={"atom_types": []})
    with pytest.raises(UserLibraryError):
        promote_state(bad)
    assert not (user_lib / "small_molecules").exists()


def test_roundtrip_failure_restores_existing_metadata(user_lib, artifacts):
    """A failed SECOND promotion must restore the pre-existing metadata."""
    frcmod, lib = artifacts
    promote_state(_small_molecule_request(frcmod, lib))
    good = (user_lib / "small_molecules" / "LIG" / "metadata.json").read_text()

    bad = _small_molecule_request(
        frcmod, lib, set_name="Broken",
        redox_state="single_state", spin_state="other",
        spin_meta={"atom_types": []},  # no residue_name → invalid
    )
    with pytest.raises(UserLibraryError):
        promote_state(bad)

    assert (user_lib / "small_molecules" / "LIG"
            / "metadata.json").read_text() == good
    # The half-written 'other' spin dir must be gone.
    assert not (user_lib / "small_molecules" / "LIG" / "single_state"
                / "other").exists()


# --------------------------------------------------------------------------- #
# metal-site shape (multi-lib + atom types)
# --------------------------------------------------------------------------- #

def test_metal_site_multi_lib_and_atom_types(user_lib, tmp_path):
    ws = tmp_path / "ws"; ws.mkdir()
    fr = ws / "site.frcmod"; fr.write_text("MASS\n")
    l1 = ws / "ZN1.lib"; l1.write_text("\"ZN1\"\n")
    l2 = ws / "CY1.lib"; l2.write_text("\"CY1\"\n")

    result = promote_state(PromotionRequest(
        family="zinc", type_name="my_zn_his3",
        set_name="mcpb_generated",
        frcmod_src=str(fr), lib_srcs=[str(l1), str(l2)],
        spin_meta={
            "residue_name": {"ZN90": "ZN1", "CYS5": "CY1"},
            "atom_types": ['{ "ZM" "Zn" "sp3" }', '{ "SZ" "S" "sp3" }'],
        },
        state_meta={"formal_charge": -1},
    ))

    meta = json.loads(Path(result["metadata_path"]).read_text())
    spin = meta["redox_states"]["single_state"]["spin_states"]["default"]
    assert spin["atom_types"] == ['{ "ZM" "Zn" "sp3" }', '{ "SZ" "S" "sp3" }']
    assert spin["residue_name"] == {"ZN90": "ZN1", "CYS5": "CY1"}
    # Multiple libs stored as a list.
    files = spin["forcefield_sets"]["mcpb_generated"]["files"]
    assert files["lib"] == ["ZN1.lib", "CY1.lib"]
    # Loader returns the lib list form.
    sets = loader.discover_forcefield_files("zinc/my_zn_his3",
                                            "single_state", "default")
    assert isinstance(sets[0]["lib"], list) and len(sets[0]["lib"]) == 2


def test_single_frcmod_is_stored_as_string_not_list(user_lib, artifacts):
    """No extra frcmods → files.frcmod stays a bare string (loader too)."""
    frcmod, lib = artifacts
    result = promote_state(_small_molecule_request(frcmod, lib))
    files = (json.loads(Path(result["metadata_path"]).read_text())
             ["redox_states"]["single_state"]["spin_states"]["default"]
             ["forcefield_sets"]["Guberman_LIG_RESP"]["files"])
    assert files["frcmod"] == "lig.frcmod"  # bare string, not ["lig.frcmod"]
    sets = loader.discover_forcefield_files("small_molecules/LIG",
                                            "single_state", "default")
    assert isinstance(sets[0]["frcmod"], str)


def test_extra_frcmods_recorded_as_list_and_loaded(user_lib, tmp_path):
    """A metal site with a bonded frcmod plus a ligand GAFF frcmod records BOTH
    in files.frcmod (a list) and the loader returns absolute paths to each, so
    the tleap generator emits a loadamberparams for every one."""
    ws = tmp_path / "ws"; ws.mkdir()
    bonded = ws / "site_1_bonded.frcmod"; bonded.write_text("MASS\nM1 54.9\n")
    ligand = ws / "e4z.frcmod"; ligand.write_text("MASS\nc 12.0\n")
    l1 = ws / "MN1.lib"; l1.write_text("\"MN1\"\n")
    l2 = ws / "EZ1.lib"; l2.write_text("\"EZ1\"\n")

    result = promote_state(PromotionRequest(
        family="metal_sites", type_name="my_site",
        set_name="mcpb_generated",
        frcmod_src=str(bonded), lib_srcs=[str(l1), str(l2)],
        extra_frcmod_srcs=[str(ligand)],
        spin_meta={"residue_name": {"MN90": "MN1", "LIG5": "EZ1"},
                   "atom_types": ['{ "M1" "Mn" "sp3" }']},
    ))

    # Both frcmods copied into the state dir.
    state_dir = Path(result["state_dir"])
    assert (state_dir / "site_1_bonded.frcmod").is_file()
    assert (state_dir / "e4z.frcmod").is_file()

    # Recorded as a list, primary first.
    files = (json.loads(Path(result["metadata_path"]).read_text())
             ["redox_states"]["single_state"]["spin_states"]["default"]
             ["forcefield_sets"]["mcpb_generated"]["files"])
    assert files["frcmod"] == ["site_1_bonded.frcmod", "e4z.frcmod"]

    # Loader hands back both absolute paths for the generator to load.
    sets = loader.discover_forcefield_files("metal_sites/my_site",
                                            "single_state", "default")
    frcmods = sets[0]["frcmod"]
    assert isinstance(frcmods, list) and len(frcmods) == 2
    assert [Path(p).name for p in frcmods] == ["site_1_bonded.frcmod", "e4z.frcmod"]
    assert all(Path(p).is_absolute() and Path(p).is_file() for p in frcmods)
