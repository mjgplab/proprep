"""
mcpb-4 (checklist Step 15) deposits ONE library entry per metal site.

Step 15 used to flatten every ``site_*`` directory into a single
``create_ff_library`` call under one prompted site_type/redox/spin. A protein
with two metal sites therefore got one merged library entry, which could not be
reused on a structure carrying only one of the sites, and one merged reuse
transformer that matched neither site: ``evaluate_redox_site`` counts the rename
table's resnames against a SINGLE redox site and requires ``met == total``, so a
table spanning two sites fails on both.

Residue NAMING still runs over the union of all sites — that is what keeps
site 2's Cys from colliding with site 1's — but the deposit is partitioned.

These tests drive the real handler with stubbed prompts and a stubbed
``create_ff_library``, so the partitioning is exercised rather than replicated.
"""

import io

import pytest
from rich.console import Console

from proprep.forcefield_prep import structure_preprocessor as sp_mod
from proprep.forcefield_prep.mcpb import integration_utils
from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor


class _WS:
    """Minimal workspace: the handler only get()s and set()s."""

    def __init__(self, d=None):
        self.d = dict(d or {})

    def get(self, k, default=None):
        return self.d.get(k, default)

    def set(self, k, v):
        self.d[k] = v


def _make_site(root, index, fingerprint_lines, mol2_stems):
    """Lay out one ``site_N`` directory the way mcpb-1..3 leave it."""
    site = root / f"site_{index}"
    models = site / "models"
    models.mkdir(parents=True)
    (models / "standard.fingerprint").write_text("\n".join(fingerprint_lines) + "\n")
    for stem in mol2_stems:
        (models / f"{stem}.mol2").write_text("@<TRIPOS>MOLECULE\n")
    bonded = site / "bonded_params"
    bonded.mkdir()
    (bonded / f"site{index}_bonded.frcmod").write_text("MASS\n")
    return site


def _preprocessor(tmp_path, workspace=None):
    obj = StructurePreprocessor.__new__(StructurePreprocessor)
    obj._output_dir = tmp_path
    obj._final_pdb = None            # section F (PDB rename) is skipped
    obj.redox_sites = []
    obj.workspace = workspace if workspace is not None else _WS()
    obj.console = Console(file=io.StringIO(), width=200)
    obj.processor = None
    obj._collect_ligand_frcmods = lambda resnames: {}
    obj._collect_ligand_atom_maps = lambda resnames: {}
    return obj


@pytest.fixture
def two_sites(tmp_path):
    # Site 1: a Mn center with an Asp ligand.
    _make_site(tmp_path, 1, [
        "202-MN-MN     1  MN -> M1",
        "108-ASP-OD1   2  o  -> Y1",
        "108-ASP-CG    3  c  -> c",
    ], ["MN202", "ASP108"])
    # Site 2: an independent Mn center with a Glu ligand.
    _make_site(tmp_path, 2, [
        "203-MN-MN     1  MN -> M2",
        "245-GLU-OE1   2  o  -> Y2",
        "245-GLU-CD    3  c  -> c",
    ], ["MN203", "GLU245"])
    return tmp_path


def _install_stubs(monkeypatch, site_types, redox="oxidized", spin="high_spin"):
    """Stub the prompts and the library writer; return the recorded calls."""
    calls = []
    type_queue = list(site_types)

    def fake_prompt(processor, prompt, **kwargs):
        if "site type" in prompt:
            return type_queue.pop(0)
        if "Redox state" in prompt:
            return redox
        if "Spin state" in prompt:
            return spin
        raise AssertionError(f"unexpected prompt: {prompt}")

    def fake_confirm(processor, prompt, **kwargs):
        return True

    def fake_create(**kwargs):
        calls.append(kwargs)
        n = len(calls)
        return {
            "library_path": f"/lib/specialized_residues/metal_sites/{kwargs['site_type']}",
            "metadata_path": f"/lib/metadata_{n}.json",
            "renamed_mol2_files": [f"/lib/{v}.lib" for v in kwargs["residue_name_map"].values()],
            "frcmod_files": list(kwargs["frcmod_files"]),
            "atom_type_entries": list(kwargs["atom_type_entries"] or []),
        }

    monkeypatch.setattr(sp_mod, "prompt_with_context", fake_prompt)
    monkeypatch.setattr(sp_mod, "confirm_with_context", fake_confirm)
    monkeypatch.setattr(integration_utils, "create_ff_library", fake_create)
    return calls


def test_each_site_gets_its_own_library_entry(two_sites, monkeypatch):
    calls = _install_stubs(monkeypatch, ["mn_asp", "mn_glu"])
    pre = _preprocessor(two_sites)

    result = pre._checklist_mcpb_4_integration()

    # One deposit per site, not one merged deposit.
    assert len(calls) == 2, f"expected 2 library entries, got {len(calls)}"
    assert [c["site_type"] for c in calls] == ["mn_asp", "mn_glu"]

    first, second = calls

    # Each entry carries only its OWN residues.
    assert set(first["residue_name_map"]) == {(202, "MN"), (108, "ASP")}
    assert set(second["residue_name_map"]) == {(203, "MN"), (245, "GLU")}

    # ...its own mol2 files...
    assert all("site_1" in m for m in first["mol2_files"])
    assert all("site_2" in m for m in second["mol2_files"])

    # ...its own bonded frcmod...
    assert len(first["frcmod_files"]) == 1 and "site1_bonded" in first["frcmod_files"][0]
    assert len(second["frcmod_files"]) == 1 and "site2_bonded" in second["frcmod_files"][0]

    # ...and only its own M*/Y* atom types.
    joined_first = " ".join(first["atom_type_entries"])
    joined_second = " ".join(second["atom_type_entries"])
    assert '"M1"' in joined_first and '"Y1"' in joined_first
    assert "M2" not in joined_first and "Y2" not in joined_first
    assert '"M2"' in joined_second and '"Y2"' in joined_second
    assert "M1" not in joined_second and "Y1" not in joined_second

    assert "2 library entries" in result["summary"]


def test_residue_names_stay_unique_across_sites(two_sites, monkeypatch):
    """Naming still runs over the union — partitioning must not reintroduce
    the collision that a per-site naming pass would cause."""
    calls = _install_stubs(monkeypatch, ["mn_asp", "mn_glu"])
    pre = _preprocessor(two_sites)

    pre._checklist_mcpb_4_integration()

    all_names = [n for c in calls for n in c["residue_name_map"].values()]
    assert len(all_names) == len(set(all_names)), f"collision across sites: {all_names}"


def test_duplicate_identity_is_rejected_and_reprompted(two_sites, monkeypatch):
    """Two sites cannot share one library key: their residue names differ, and
    promote_state overwrites a repeated key, so the first site would be lost."""
    # Site 2 is first offered the same identity as site 1, then a distinct one.
    calls = _install_stubs(monkeypatch, ["mn_site", "mn_site", "mn_site_b"])
    pre = _preprocessor(two_sites)

    pre._checklist_mcpb_4_integration()

    assert [c["site_type"] for c in calls] == ["mn_site", "mn_site_b"]


def test_workspace_registers_every_site_for_tleap(two_sites, monkeypatch):
    """tLEaP runs once over the whole structure, so every site's files and
    types are still registered together — deduped."""
    _install_stubs(monkeypatch, ["mn_asp", "mn_glu"])
    ws = _WS()
    pre = _preprocessor(two_sites, workspace=ws)

    pre._checklist_mcpb_4_integration()

    libs = ws.get("preprocessing_lib_files")
    frcmods = ws.get("preprocessing_frcmod_files")
    types = " ".join(ws.get("preprocessing_atom_types"))

    assert len(libs) == 4                      # 2 residues per site
    assert len(libs) == len(set(libs))
    assert len(frcmods) == 2
    for t in ("M1", "Y1", "M2", "Y2"):
        assert f'"{t}"' in types, f"{t} missing from tLEaP atom types: {types}"


def test_single_site_still_deposits_once(tmp_path, monkeypatch):
    """The single-site path is unchanged: one prompt set, one entry."""
    _make_site(tmp_path, 1, [
        "202-ZN-ZN     1  ZN -> M1",
        "108-CYS-SG    2  s  -> Y1",
    ], ["ZN202", "CYS108"])
    calls = _install_stubs(monkeypatch, ["zinc_cys"])
    pre = _preprocessor(tmp_path)

    result = pre._checklist_mcpb_4_integration()

    assert len(calls) == 1
    assert calls[0]["site_type"] == "zinc_cys"
    assert set(calls[0]["residue_name_map"]) == {(202, "ZN"), (108, "CYS")}
    assert "1 library entry" in result["summary"]


def test_site_dirs_ordered_numerically(tmp_path, monkeypatch):
    """site_10 must sort after site_2 — the numeric suffix indexes the
    redox-site list, so a lexical sort would mislabel the sites."""
    for idx, (resid, resname) in enumerate(
            [(1, "MN"), (2, "MN"), (10, "MN")], start=1):
        n = {1: 1, 2: 2, 3: 10}[idx]
        _make_site(tmp_path, n, [
            f"{resid}-{resname}-MN   1  MN -> M{n}",
        ], [f"{resname}{resid}"])

    calls = _install_stubs(monkeypatch, ["a", "b", "c"])
    pre = _preprocessor(tmp_path)

    pre._checklist_mcpb_4_integration()

    ordered = [c["mol2_files"][0] for c in calls]
    assert "site_1/" in ordered[0]
    assert "site_2/" in ordered[1]
    assert "site_10/" in ordered[2]


_PDB_LINES = [
    "ATOM      1  N   ASP A 108      11.104  13.207  10.000  1.00  0.00           N",
    "ATOM      2  OD1 ASP A 108      12.104  13.207  10.000  1.00  0.00           O",
    "ATOM      3  N   GLU A 245      21.104  13.207  10.000  1.00  0.00           N",
    "ATOM      4  OE1 GLU A 245      22.104  13.207  10.000  1.00  0.00           O",
    "HETATM    5 MN    MN A 202      12.500  13.500  10.000  1.00  0.00          MN",
    "HETATM    6 MN    MN A 203      22.500  13.500  10.000  1.00  0.00          MN",
    "END",
]


def test_one_reuse_transformer_per_site(two_sites, monkeypatch):
    """The merged rename table could match neither site.

    AutoRenameTransformerBase.evaluate_redox_site counts the table's resnames
    against ONE redox site and requires met == total, so a table holding both
    sites' residues fails on both. Each site must get its own transformer.
    """
    _install_stubs(monkeypatch, ["mn_asp", "mn_glu"])

    emitted = []

    def fake_emit(rename_table, **kwargs):
        emitted.append((rename_table, kwargs))
        return f"/transformers/{kwargs['name']}.py"

    monkeypatch.setattr(
        "proprep.redoxsite_prep.transformation.auto_rename.emit_rename_transformer",
        fake_emit,
    )

    pdb = two_sites / "prepared.pdb"
    pdb.write_text("\n".join(_PDB_LINES) + "\n")
    pre = _preprocessor(two_sites, workspace=_WS({"prepared_pdb": str(pdb)}))

    pre._checklist_mcpb_4_integration()

    assert len(emitted) == 2, f"expected one transformer per site, got {len(emitted)}"

    first_resnames = {e["resname"] for e in emitted[0][0]}
    second_resnames = {e["resname"] for e in emitted[1][0]}

    # Each table holds only its own site's residues — never the union.
    assert first_resnames == {"MN", "ASP"}
    assert second_resnames == {"MN", "GLU"}

    # Each transformer points at its own site's deposited library.
    assert emitted[0][1]["site_types"] == ["mn_asp"]
    assert emitted[1][1]["site_types"] == ["mn_glu"]


def test_one_failed_deposit_does_not_lose_the_others(two_sites, monkeypatch):
    """A per-site deposit can now fail on its own; the others must survive."""
    calls = []

    def fake_prompt(processor, prompt, **kwargs):
        if "site type" in prompt:
            return f"site_{len(calls)}"
        return "default"

    def fake_create(**kwargs):
        calls.append(kwargs)
        if len(calls) == 1:
            raise ValueError("create_ff_library: no frcmod files to deposit")
        return {
            "library_path": "/lib/specialized_residues/metal_sites/ok",
            "metadata_path": "/lib/metadata.json",
            "renamed_mol2_files": ["/lib/A.lib"],
            "frcmod_files": ["/lib/a.frcmod"],
            "atom_type_entries": [],
        }

    monkeypatch.setattr(sp_mod, "prompt_with_context", fake_prompt)
    monkeypatch.setattr(sp_mod, "confirm_with_context",
                        lambda processor, prompt, **kw: True)
    monkeypatch.setattr(integration_utils, "create_ff_library", fake_create)

    pre = _preprocessor(two_sites)
    result = pre._checklist_mcpb_4_integration()

    assert len(calls) == 2, "second site was not attempted after the first failed"
    assert "1 library entry" in result["summary"]
    assert "not deposited" in result["summary"]
