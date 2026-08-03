#!/usr/bin/env python3
"""Envelope-retention unit test for cluster cutout.

No PBSA, no force-field math — just cluster construction. Validates the
two invariants required for multi-residue MCPB groups:

  1. The target site's envelope is unconditionally retained, even when
     `cluster_radius` is so small it would otherwise clip envelope members.
  2. Any other detected_redox_sites entry that contributes residues to the
     cluster has its full envelope included (all-or-nothing rule); partial
     inclusion is forbidden because it leaves fractional charge.

Uses BPTI as a stand-in for a multi-residue MCPB system: we construct
synthetic "redox-site envelopes" by grouping a handful of BPTI residues
into a fake envelope and check that cluster_for_site honors them.
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import List, Tuple

import parmed as pmd

from proprep.pb_titrate.cluster import (
    cluster_for_site, cluster_for_site_pair, cluster_for_target,
)
from proprep.pb_titrate.sites import Site, discover_sites


HERE = Path(__file__).resolve().parent
PRMTOP = HERE.parent / "bpti_asp3" / "bpti.prmtop"
RST7   = HERE.parent / "bpti_asp3" / "bpti.rst7"


def _residue_centroid(s: pmd.Structure, idx: int) -> Tuple[float, float, float]:
    r = s.residues[idx]
    xs = [a.xx for a in r.atoms]
    ys = [a.xy for a in r.atoms]
    zs = [a.xz for a in r.atoms]
    return sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)


def _distance(a, b) -> float:
    return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2) ** 0.5


def _make_fake_redox_site(s: pmd.Structure, idxs: List[int],
                           site_id: str) -> SimpleNamespace:
    """Build an object that quacks like comprehensive_redox_detector.RedoxSite
    enough for `discover_sites` and `resolve_envelope` to consume.

    Only `residue_groups` and `site_id` are required.
    """
    rg = {}
    for i in idxs:
        r = s.residues[i]
        key = (str(getattr(r, "chain", "") or ""),
               int(r.number) + 1,
               str(getattr(r, "insertion_code", "") or ""))
        rg[key] = []  # value unused by pb_titrate.sites
    return SimpleNamespace(site_id=site_id, residue_groups=rg)


def test_single_residue_path_unchanged():
    """`cluster_for_target` (the legacy single-residue API) must give the
    same set of residues with or without the new Site abstraction."""
    s = pmd.load_file(str(PRMTOP), str(RST7))
    target = next(r for r in s.residues
                  if r.name == "AS4" and r.number+1 == 3)

    legacy = set(cluster_for_target(s, target, R=15.0, bond_depth=3))

    # Same answer via the new Site path with a single-residue envelope
    site = Site(titrating_residue=target, envelope_residues=[target])
    new = set(cluster_for_site(s, site, R=15.0, bond_depth=3))

    assert legacy == new, (
        f"Single-residue path drift: legacy={len(legacy)} new={len(new)} "
        f"diff={legacy ^ new}")
    print(f"  ✓ single-residue path: legacy and new agree on {len(legacy)} residues")


def test_envelope_retained_under_tight_radius():
    """An envelope containing residues far from the titrating one must
    still be retained when cluster_radius is too small to reach them."""
    s = pmd.load_file(str(PRMTOP), str(RST7))

    # AS4-3 is at idx 2. Pick two residues that are spatially distant from
    # AS4-3 to form a synthetic envelope. Verify that even at R=2 Å (which
    # would normally find no other residue), all envelope members are kept.
    as4_3 = s.residues[2]
    far_indices = []
    as4_centroid = _residue_centroid(s, 2)
    distances = sorted(
        ((i, _distance(as4_centroid, _residue_centroid(s, i)))
         for i in range(len(s.residues)) if i != 2),
        key=lambda x: -x[1]
    )
    far_indices = [distances[0][0], distances[1][0]]  # two farthest residues

    envelope_residues = [as4_3] + [s.residues[i] for i in far_indices]
    site = Site(titrating_residue=as4_3, envelope_residues=envelope_residues)

    keep = set(cluster_for_site(s, site, R=2.0, bond_depth=0))
    expected = {2, *far_indices}
    assert expected.issubset(keep), (
        f"Envelope dropped at tight R: expected {expected} ⊆ keep but "
        f"got keep={keep} (missing {expected - keep})")
    print(f"  ✓ envelope retained at R=2 Å, bond_depth=0: "
          f"{len(envelope_residues)} envelope residues all in cluster of {len(keep)}")


def test_partial_overlap_triggers_full_inclusion():
    """If another redox site has any residue inside the cluster, its full
    envelope must be included (all-or-nothing).

    Geometry in BPTI: PRO-2 (idx 1) sits ~5 Å from AS4-3; LYS-15 (idx 14)
    sits ~29 Å away. With R=8 Å around AS4-3, PRO-2 is in range and LYS-15
    is not — i.e. partial overlap of the envelope {1, 14}. The rule must
    pull idx 14 into the cluster.
    """
    s = pmd.load_file(str(PRMTOP), str(RST7))
    as4_3 = s.residues[2]

    # Synthetic "other redox site" envelope: one near residue, one far.
    other_envelope = [1, 14]
    R_test = 8.0

    # Sanity: confirm that without the rule, exactly one of {1, 14} ends up
    # inside the radius.
    site_no_rule = Site(
        titrating_residue=as4_3,
        envelope_residues=[as4_3],
        other_redox_envelopes=[],
    )
    keep_no_rule = set(cluster_for_site(s, site_no_rule, R=R_test, bond_depth=0))
    in_range = {x for x in (1, 14) if x in keep_no_rule}
    out_of_range = {x for x in (1, 14) if x not in keep_no_rule}
    assert len(in_range) == 1 and len(out_of_range) == 1, (
        f"Geometry assumption failed: at R={R_test} Å neither one nor both "
        f"of {{1,14}} are in range — got in={in_range}, out={out_of_range}. "
        f"Cluster size {len(keep_no_rule)}. Re-pick test geometry.")

    # With the all-or-nothing rule, both must end up in the cluster.
    site_with_rule = Site(
        titrating_residue=as4_3,
        envelope_residues=[as4_3],
        other_redox_envelopes=[other_envelope],
    )
    keep_with_rule = set(cluster_for_site(s, site_with_rule, R=R_test, bond_depth=0))
    assert {1, 14}.issubset(keep_with_rule), (
        f"Partial overlap not fixed: keep_with_rule had {keep_with_rule & {1,14}}, "
        f"missing {set([1,14]) - keep_with_rule}")
    print(f"  ✓ partial overlap detected (idx {in_range} inside R={R_test}, "
          f"{out_of_range} outside); rule pulled in the full envelope "
          f"{{1, 14}}; cluster grew {len(keep_no_rule)} → {len(keep_with_rule)}")


def test_discover_sites_finds_titratables():
    """`discover_sites(structure)` with no redox sites must return one Site
    per AS4/GL4/HIP/LYS/CYS/TYR/PRN residue, each with a single-residue
    envelope."""
    s = pmd.load_file(str(PRMTOP), str(RST7))
    sites = discover_sites(s)
    names = [(site.resname, site.resnum) for site in sites]
    expected = {("AS4", 3), ("TYR", 10), ("LYS", 15), ("TYR", 21),
                ("TYR", 23), ("LYS", 26), ("TYR", 35), ("LYS", 41),
                ("LYS", 46)}
    assert set(names) == expected, f"discovery mismatch: got {names}"
    for site in sites:
        assert len(site.envelope_residues) == 1
        assert site.redox_site_id is None
        assert site.other_redox_envelopes == []
    print(f"  ✓ discover_sites found {len(sites)} single-residue sites in BPTI")


def test_discover_sites_uses_redox_envelope():
    """`discover_sites(structure, redox_sites)` must wrap a titratable
    residue inside a synthetic redox-site envelope into a multi-residue
    Site, with the other redox sites surfaced on `other_redox_envelopes`."""
    s = pmd.load_file(str(PRMTOP), str(RST7))

    # Fake redox-site #1: bundles AS4-3 (idx 2) + neighbors at idx 1, 3, 4.
    rs1 = _make_fake_redox_site(s, [1, 2, 3, 4], site_id="fake_heme_1")
    # Fake redox-site #2: bundles LYS-15 (idx 14) + neighbors idx 13, 16.
    rs2 = _make_fake_redox_site(s, [13, 14, 16], site_id="fake_heme_2")

    sites = discover_sites(s, redox_sites=[rs1, rs2])

    # Find the AS4-3 site
    as4_site = next(site for site in sites if site.resname == "AS4")
    assert as4_site.redox_site_id == "fake_heme_1"
    assert as4_site.envelope_idxs == {1, 2, 3, 4}
    # The OTHER redox site (rs2) must appear in other_redox_envelopes
    assert [13, 14, 16] in as4_site.other_redox_envelopes, (
        f"rs2 envelope not in as4_site.other_redox_envelopes: "
        f"{as4_site.other_redox_envelopes}")
    # rs1 (the owner) must NOT appear in other_redox_envelopes
    rs1_idxs = [1, 2, 3, 4]
    assert rs1_idxs not in as4_site.other_redox_envelopes

    # Find the LYS-15 site (idx 14)
    lys15_site = next(s_ for s_ in sites if s_.resname == "LYS" and s_.resnum == 15)
    assert lys15_site.redox_site_id == "fake_heme_2"
    assert lys15_site.envelope_idxs == {13, 14, 16}
    assert [1, 2, 3, 4] in lys15_site.other_redox_envelopes

    # Other titratable residues should still be single-residue, but their
    # other_redox_envelopes should list BOTH fake redox sites.
    tyr_10 = next(s_ for s_ in sites if s_.resname == "TYR" and s_.resnum == 10)
    assert tyr_10.redox_site_id is None
    assert tyr_10.envelope_idxs == {9}
    assert sorted(tyr_10.other_redox_envelopes) == sorted([[1,2,3,4], [13,14,16]])

    print(f"  ✓ discover_sites built multi-residue envelopes for AS4-3 "
          f"(envelope={sorted(as4_site.envelope_idxs)}) and LYS-15 "
          f"(envelope={sorted(lys15_site.envelope_idxs)})")


def main():
    print("=== Envelope retention unit test ===\n")
    print("[1] Single-residue path (legacy API parity)")
    test_single_residue_path_unchanged()
    print("\n[2] Envelope retention under tight radius")
    test_envelope_retained_under_tight_radius()
    print("\n[3] Partial-overlap → all-or-nothing inclusion")
    test_partial_overlap_triggers_full_inclusion()
    print("\n[4] discover_sites: BPTI without redox sites")
    test_discover_sites_finds_titratables()
    print("\n[5] discover_sites: BPTI with synthetic redox envelopes")
    test_discover_sites_uses_redox_envelope()
    print("\n=== ALL TESTS PASSED ===")


if __name__ == "__main__":
    main()
