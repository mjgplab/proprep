"""Redox sites imported from JSON still match their transformers.

Seen on 6R2Q: sites detected earlier were imported from
``6R2Q_noH_redox_sites.json`` and, in the Redox Site Preparer, both disulfide
sites showed "Compatible Transformers: None" and "disulfide ✗ Not compatible",
2/3 requirements met: "Both centers must have is_disulfide_bonded property
(expected 2, got 0)". The exporter wrote seven fields per center and no
``properties``; the importer built centers with none. Everything detection had
learned about a center was lost in the round trip. Treating such a site as
no_transformation leaves two CYS that tLEaP cannot bond.

The file below is that site as 1.21.0 exported it (C:112-C:118 of 6R2Q).
"""

import io
import json

import numpy as np
from rich.console import Console

from proprep.redoxsite_prep.transformation.transformers.disulfide import DisulfideTransformer
from proprep.structure_prep import comprehensive_redox_detector as crd

SG = {112: [70.176, 10.988, 5.768], 118: [70.511, 11.709, 3.888]}


def _center(resid, properties=None):
    data = {"chain": "C", "resname": "CYS", "resid": resid, "atom_name": "SG", "element": "S",
            "center_type": "redox_amino_acid", "coordinates": SG[resid]}
    if properties is not None:
        data["properties"] = properties
    return data


def _site_file(tmp_path, centers, treatment="bonded", chemical_type="disulfide"):
    site = {
        "site_id": "site_1", "structure_id": "temp", "site_type": "disulfide", "centers": centers,
        "atoms": [dict(_center(r), **{"center_type": None}) for r in SG],
        "bonds": [{"atom1": {"chain": "C", "resname": "CYS", "resid": 112, "atom_name": "SG"},
                   "atom2": {"chain": "C", "resname": "CYS", "resid": 118, "atom_name": "SG"},
                   "bond_type": "interresidue", "chemical_type": chemical_type, "distance": 2.0412,
                   "atom1_element": "S", "atom2_element": "S",
                   "atom1_coordinates": SG[112], "atom2_coordinates": SG[118], "treatment": treatment}],
    }
    for atom in site["atoms"]:
        atom.pop("center_type")
    path = tmp_path / "sites.json"
    path.write_text(json.dumps({"source_pdb": "x.pdb", "transformer_mappings": {"disulfide": "disulfide"},
                                "sites": [site]}))
    return str(path)


def _import(path):
    sites, _ = crd._import_from_json(path)
    return sites[0]


def test_a_file_written_without_properties_still_matches_the_disulfide_transformer(tmp_path):
    site = _import(_site_file(tmp_path, [_center(112), _center(118)]))
    evaluation = DisulfideTransformer.evaluate_redox_site(site)
    assert evaluation.is_valid and evaluation.requirements_met == 3
    first, second = site.centers
    assert first.properties == {"is_disulfide_bonded": True, "disulfide_partner_chain": "C",
                                "disulfide_partner_res": 118, "disulfide_bond_distance": 2.04}
    assert second.properties["disulfide_partner_res"] == 112


def test_only_a_bonded_disulfide_bond_says_so(tmp_path):
    """A pair the user declined to bond, or a bond of another kind, is not made a disulfide."""
    for kwargs in ({"treatment": "nonbonded"}, {"chemical_type": "covalent"}):
        site = _import(_site_file(tmp_path, [_center(112), _center(118)], **kwargs))
        assert all("is_disulfide_bonded" not in c.properties for c in site.centers)
        assert not DisulfideTransformer.evaluate_redox_site(site).is_valid


def test_what_the_file_says_about_a_center_is_not_overridden(tmp_path):
    rejected = {"is_disulfide_bonded": False, "disulfide_rejected": True}
    site = _import(_site_file(tmp_path, [_center(112, rejected), _center(118, dict(rejected))]))
    assert [c.properties["is_disulfide_bonded"] for c in site.centers] == [False, False]


def test_properties_survive_export_and_import(tmp_path, monkeypatch):
    site = _import(_site_file(tmp_path, [_center(112), _center(118)]))
    live = site.centers[0].properties
    live.update({"from_ssbond_record": np.bool_(True), "disulfide_bond_distance": np.float64(2.04),
                 "neighbours": (1, 2), "not_for_json": object()})
    site.atoms[0].properties["serial_number"] = np.int64(8059)

    monkeypatch.chdir(tmp_path)
    crd._export_to_json([site], str(tmp_path / "roundtrip.pdb"), console=Console(file=io.StringIO()),
                        transformer_mappings={"disulfide": "disulfide"})
    written = next(tmp_path.glob("roundtrip*redox_sites*.json"))
    raw = json.loads(written.read_text())["sites"][0]
    assert raw["centers"][0]["properties"]["from_ssbond_record"] is True
    assert "not_for_json" not in raw["centers"][0]["properties"]
    assert raw["atoms"][0]["properties"] == {"serial_number": 8059}
    assert "properties" not in raw["atoms"][1]                       # nothing to say, nothing written

    again = _import(str(written))
    assert again.centers[0].properties == {
        "is_disulfide_bonded": True, "disulfide_partner_chain": "C", "disulfide_partner_res": 118,
        "disulfide_bond_distance": 2.04, "from_ssbond_record": True, "neighbours": [1, 2]}
    assert again.atoms[0].properties == {"serial_number": 8059}
    assert DisulfideTransformer.evaluate_redox_site(again).is_valid
