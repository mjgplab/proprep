"""
Per-site file resolution, and a library frcmod named for its entry.

``step_results`` is restored from ``mcpb_step_results``, a single workspace key
every site shares, so ``step_1`` belongs to whichever site wrote it last. The
checklist sets ``step_3a`` per site, so that is what identifies the site being
processed.

Reading step_1 unconditionally cross-wired site 1 to site 2's
standard.fingerprint on a 4UHX build. Their PDB serial ranges do not overlap,
so every atom-type lookup missed and the deposited libraries were written with
the ``XX`` placeholder -- ~150 tleap "could not find parameter" errors, none of
them visible until the build.

Separately: the bonded frcmod was named for the working directory's ``site_N``,
which carries no meaning once deposited and means something different in every
structure. It now takes the library entry's name.
"""

import pytest

from proprep.forcefield_prep.metal_site_parameterizer import MetalSiteWorkflowManager
from proprep.forcefield_prep.mcpb.integration_utils import _slug


# --------------------------------------------------------------------------- #
# which site's files
# --------------------------------------------------------------------------- #

def _manager(step_results):
    manager = MetalSiteWorkflowManager.__new__(MetalSiteWorkflowManager)
    manager.step_results = step_results
    return manager


def test_step_3a_identifies_the_site():
    """step_3a is set per site by the checklist; step_1 is shared."""
    manager = _manager({
        "step_1": {"output_files": {
            "standard_fingerprint": "/run/site_2/models/standard.fingerprint"}},
        "step_3a": {"output_dir": "/run/site_1/models"},
    })

    assert str(manager._site_models_dir()) == "/run/site_1/models"


def test_the_shared_step_1_does_not_win():
    """The reported failure: site 1 processed with site 2's fingerprint."""
    manager = _manager({
        "step_1": {"output_files": {
            "standard_fingerprint": "/run/site_2/models/standard.fingerprint"}},
        "step_3a": {"output_dir": "/run/site_1/models"},
    })

    assert "site_2" not in str(manager._site_models_dir())


def test_step_1_is_used_when_there_is_no_step_3a():
    """The standalone workflow, where step_1 is this site's."""
    manager = _manager({"step_1": {"output_files": {
        "standard_fingerprint": "/run/site_1/models/standard.fingerprint"}}})

    assert str(manager._site_models_dir()) == "/run/site_1/models"


def test_large_pdb_is_an_accepted_fallback():
    manager = _manager({"step_1": {"output_files": {
        "large_pdb": "/run/site_1/models/large.pdb"}}})

    assert str(manager._site_models_dir()) == "/run/site_1/models"


def test_an_empty_step_3a_falls_through():
    manager = _manager({
        "step_1": {"output_files": {"large_pdb": "/run/site_1/models/large.pdb"}},
        "step_3a": {},
    })

    assert str(manager._site_models_dir()) == "/run/site_1/models"


def test_no_usable_path_raises_rather_than_guessing():
    with pytest.raises(ValueError, match="models directory"):
        _manager({"step_1": {"output_files": {}}})._site_models_dir()


# --------------------------------------------------------------------------- #
# the deposited frcmod's name
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("entry,expected", [
    ("4hux_fe2s2", "4hux_fe2s2"),
    ("4hux_moco", "4hux_moco"),
    ("my site 1", "my_site_1"),
    ("Fe4S4/oxidized", "Fe4S4_oxidized"),
    ("__odd__", "odd"),
])
def test_entry_names_become_usable_filename_stems(entry, expected):
    assert _slug(entry) == expected


def test_an_empty_entry_name_yields_nothing_to_use():
    """The caller falls back to the source stem rather than writing '_bonded'."""
    assert _slug("") == ""
    assert _slug(None) == ""


def _deposit(tmp_path, monkeypatch, site_type, sources, extras=()):
    """Run create_ff_library's staging and report the staged frcmod names."""
    from proprep.forcefield_prep.mcpb import integration_utils as iu

    staged_names = []

    def fake_convert(renamed_mol2, staged_frcmod, atom_type_entries, staging):
        staged_names.extend(sorted(p.name for p in map(__import__('pathlib').Path,
                                                       staged_frcmod)))
        raise RuntimeError("stop after staging")

    monkeypatch.setattr(iu, "_convert_mol2_to_lib", fake_convert)

    mol2 = tmp_path / "FES1311.mol2"
    mol2.write_text("@<TRIPOS>MOLECULE\nFES\n")
    fingerprint = tmp_path / "standard.fingerprint"
    fingerprint.write_text("1311-FES-FE1 20396 FE -> M1\n")
    frcmods = []
    for name in sources:
        p = tmp_path / name
        p.write_text("MASS\n")
        frcmods.append(str(p))
    extra_paths = []
    for name in extras:
        p = tmp_path / name
        p.write_text("MASS\n")
        extra_paths.append(str(p))

    with pytest.raises(RuntimeError, match="stop after staging"):
        iu.create_ff_library(
            site_type=site_type,
            description="test entry",
            mol2_files=[str(mol2)],
            frcmod_files=frcmods,
            fingerprint_path=str(fingerprint),
            assignments_path=None,
            residue_name_map={(1311, "FES"): "FS1"},
            redox_state="oxidized",
            spin_state="high_spin",
            extra_frcmod_files=extra_paths or None,
            base_dir=tmp_path / "lib",
        )
    return staged_names


def test_the_bonded_frcmod_is_named_for_the_entry(tmp_path, monkeypatch):
    names = _deposit(tmp_path, monkeypatch, "4hux_fe2s2", ["site_1_bonded.frcmod"])

    assert "4hux_fe2s2_bonded.frcmod" in names
    assert "site_1_bonded.frcmod" not in names


def test_a_ligand_frcmod_keeps_its_own_name(tmp_path, monkeypatch):
    """Those identify the ligand and are shared between entries."""
    names = _deposit(tmp_path, monkeypatch, "4hux_moco",
                     ["site_2_bonded.frcmod"], extras=["mte.frcmod"])

    assert "4hux_moco_bonded.frcmod" in names
    assert "mte.frcmod" in names


def test_two_bonded_frcmods_do_not_collide(tmp_path, monkeypatch):
    names = _deposit(tmp_path, monkeypatch, "entry",
                     ["site_1_bonded.frcmod", "site_2_bonded.frcmod"])

    assert sorted(names) == ["entry_2_bonded.frcmod", "entry_bonded.frcmod"]
