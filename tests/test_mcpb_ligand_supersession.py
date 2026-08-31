"""Regression: a metal site's organic ligand must not reach tLEaP twice.

On 7K0W the baloxavir ligand E4Z is parameterized standalone by the
small-molecule parameterizer (``small_molecule_params_E4Z/e4z.{lib,frcmod}``)
and then again as part of the di-Mn site, which renames the residue E4Z -> EZ1
and deposits its own copy of the GAFF frcmod into the site library. Before this
fix both registrations survived into the template, so tLEaP saw:

    loadoff       .../small_molecule_params_E4Z/e4z.lib      <- unit E4Z
    loadoff       .../metal_sites/DiMn_E4Z/.../EZ1.lib       <- unit EZ1
    loadamberparams .../small_molecule_params_E4Z/e4z.frcmod
    loadamberparams .../metal_sites/DiMn_E4Z/.../e4z.frcmod

The standalone lib names a unit that no longer occurs in the renamed PDB, and
the frcmod is loaded from two paths. mcpb-4 must drop the superseded halves.
"""

import os

from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor


SMALL = "/demo/7K0W_baloxavir_SI/small_molecule_params_E4Z"
SITE = "/home/u/.proprep/forcefield_params/specialized_residues/metal_sites/DiMn_E4Z/default/high_spin"


class _WS:
    def __init__(self, data):
        self._d = dict(data)

    def get(self, key, default=None):
        return self._d.get(key, default)


def _prep(organic=None, orgmet=None):
    p = StructurePreprocessor.__new__(StructurePreprocessor)
    p.workspace = _WS({
        "preprocessing_organic_ff": organic or {},
        "preprocessing_organometallic_ff": orgmet or {},
    })
    return p


def test_standalone_ligand_files_are_superseded():
    p = _prep(organic={
        "E4Z": {"lib_file": f"{SMALL}/e4z.lib",
                "mol2_file": f"{SMALL}/e4z.mol2",
                "frcmod_file": f"{SMALL}/e4z.frcmod"},
        # A ligand NOT in the site must be left completely alone.
        "LIG": {"lib_file": "/demo/other/lig.lib",
                "frcmod_file": "/demo/other/lig.frcmod"},
    })

    superseded = p._superseded_component_files(
        {"E4Z", "GLU", "HIS"}, {"E4Z": f"{SMALL}/e4z.frcmod"})

    assert os.path.realpath(f"{SMALL}/e4z.lib") in superseded
    assert os.path.realpath(f"{SMALL}/e4z.mol2") in superseded
    assert os.path.realpath(f"{SMALL}/e4z.frcmod") in superseded
    assert os.path.realpath("/demo/other/lig.lib") not in superseded
    assert os.path.realpath("/demo/other/lig.frcmod") not in superseded


def test_ondisk_fallback_frcmod_is_superseded_without_registry_entry():
    # When the organic step cached, the registry may carry no entry at all and
    # _collect_ligand_frcmods resolves the frcmod by globbing. That path still
    # gets copied into the site library, so it must still be superseded.
    p = _prep(organic={})
    superseded = p._superseded_component_files(
        {"E4Z"}, {"E4Z": f"{SMALL}/e4z.frcmod"})
    assert os.path.realpath(f"{SMALL}/e4z.frcmod") in superseded


def test_merge_drops_superseded_and_appends_site_library():
    existing_libs = [f"{SMALL}/e4z.lib", "/demo/other/lig.lib"]
    site_libs = [f"{SITE}/EZ1.lib", f"{SITE}/MN1.lib", f"{SITE}/HD1.lib"]
    superseded = {os.path.realpath(f"{SMALL}/e4z.lib")}

    merged = StructurePreprocessor._merge_ff_file_list(
        existing_libs, site_libs, superseded)

    assert f"{SMALL}/e4z.lib" not in merged      # superseded by EZ1.lib
    assert "/demo/other/lig.lib" in merged        # unrelated ligand survives
    assert merged[-3:] == site_libs


def test_merge_drops_superseded_frcmod_keeping_library_copy():
    existing = [f"{SMALL}/e4z.frcmod"]
    additions = [f"{SITE}/site_1_bonded.frcmod", f"{SITE}/e4z.frcmod"]
    superseded = {os.path.realpath(f"{SMALL}/e4z.frcmod")}

    merged = StructurePreprocessor._merge_ff_file_list(
        existing, additions, superseded)

    # Exactly one e4z.frcmod, and it is the site-library copy.
    assert merged == [f"{SITE}/site_1_bonded.frcmod", f"{SITE}/e4z.frcmod"]


def test_merge_dedupes_on_repeated_mcpb4_runs():
    additions = [f"{SITE}/EZ1.lib", f"{SITE}/MN1.lib"]
    first = StructurePreprocessor._merge_ff_file_list([], additions, set())
    second = StructurePreprocessor._merge_ff_file_list(first, additions, set())
    assert second == first


def test_merge_dedupes_relative_and_absolute_spellings(tmp_path):
    f = tmp_path / "e4z.lib"
    f.write_text("!!index\n")
    rel = os.path.relpath(str(f), os.getcwd())

    merged = StructurePreprocessor._merge_ff_file_list([str(f)], [rel], set())
    assert len(merged) == 1


def test_merge_tolerates_empty_and_none_entries():
    merged = StructurePreprocessor._merge_ff_file_list(
        None, [f"{SITE}/EZ1.lib", None, ""], set())
    assert merged == [f"{SITE}/EZ1.lib"]
