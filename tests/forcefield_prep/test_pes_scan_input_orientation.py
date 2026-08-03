"""
Regression tests for parse_pes_scan_log on a NoSymm relaxed dihedral scan.

The modified-AA from-structure route writes its scan with NoSymm (so the atom
order stays stable for the per-point ESP job), which makes Gaussian print
"Input orientation:" and never "Standard orientation:". Two bugs followed:

  1. Geometry extraction keyed only on "Standard orientation:", so a NoSymm scan
     yielded ZERO geometries — success=True but no points — and step 4 then
     failed with "Could not parse the relaxed scan log".

  2. Energies were taken from the "Summary of Optimized Potential Surface Scan"
     section, whose column-batched layout (an index header row "1 2 3 4 5" above
     each block of five) parsed as spurious (angle=1, energy=5) pairs — one
     garbage point per block, and an absurd energy span.

Both are exercised here on a synthetic 3-point, 3-atom NoSymm scan log.
"""

from rich.console import Console

from proprep.forcefield_prep.pes_scan_refinement import parse_pes_scan_log

HARTREE_TO_KCAL = 627.509


def _input_orientation_block(coords):
    """A Gaussian 'Input orientation:' block for the given (Z, x, y, z) rows."""
    head = (
        "                          Input orientation:                          \n"
        " ---------------------------------------------------------------------\n"
        " Center     Atomic      Atomic             Coordinates (Angstroms)\n"
        " Number     Number       Type             X           Y           Z\n"
        " ---------------------------------------------------------------------\n"
    )
    body = "".join(
        f" {i:6d} {z:10d}           0     {x:11.6f} {y:11.6f} {zc:11.6f}\n"
        for i, (z, x, y, zc) in enumerate(coords, 1)
    )
    tail = " ---------------------------------------------------------------------\n"
    return head + body + tail


def _scan_point(coords, energy):
    """One relaxed-scan point: an SCF energy, its geometry, its stationary point."""
    return (
        f" SCF Done:  E(RB3LYP) = {energy:.9f}     A.U. after   9 cycles\n"
        + _input_orientation_block(coords)
        + "    -- Stationary point found.\n"
    )


def _synthetic_nosymm_scan_log(energies):
    """A minimal NoSymm relaxed-scan log: no 'Standard orientation:' anywhere,
    plus a column-batched summary section whose index header would mis-parse."""
    natoms_line = " NAtoms=      3 NQM=        3 NQMF=       0\n"
    # Three atoms (C, O, H); a distinct x per point so points are not identical.
    points = "".join(
        _scan_point([(6, 0.0 + k, 0.0, 0.0), (8, 1.2, 0.0, 0.0), (1, 2.0, 0.0, 0.0)], e)
        for k, e in enumerate(energies)
    )
    # The trap: an index header "1 2 3" over an Eigenvalues row, then geometry
    # (R/D) rows. The old parser skipped the Eigenvalues row and read the header
    # as (angle=1, energy=3).
    summary = (
        " Summary of Optimized Potential Surface Scan (add -0.0 to energies):\n"
        "                           1         2         3\n"
        "     Eigenvalues --    "
        + "  ".join(f"{e:.5f}" for e in energies) + "\n"
        "           R1           1.20000   1.20100   1.20200\n"
        "           D1           0.00000  15.00000  30.00000\n"
    )
    return (natoms_line + points + summary
            + " Normal termination of Gaussian 16.\n")


def test_nosymm_scan_geometries_are_parsed(tmp_path):
    energies = [-382.100000, -382.090000, -382.095000]
    log = tmp_path / "cs1_xtal_scan.log"
    log.write_text(_synthetic_nosymm_scan_log(energies))

    r = parse_pes_scan_log(str(log), Console())

    # Bug 1: geometries came back empty despite a completed scan.
    assert r["success"] is True
    assert len(r["geometries"]) == 3
    assert len(r["energies"]) == 3
    assert r["n_atoms"] == 3
    assert r["elements"] == ["C", "O", "H"]


def test_nosymm_scan_energies_are_physical_not_summary_indices(tmp_path):
    energies = [-382.100000, -382.090000, -382.095000]
    log = tmp_path / "cs1_xtal_scan.log"
    log.write_text(_synthetic_nosymm_scan_log(energies))

    r = parse_pes_scan_log(str(log), Console())

    # Bug 2: the real (paired) energies, not the summary's index-header garbage.
    assert r["energies"] == energies
    span_kcal = (max(r["energies"]) - min(r["energies"])) * HARTREE_TO_KCAL
    # ~6 kcal/mol here; the index-header bug produced an absurd span (>1000).
    assert span_kcal < 50.0
    # The lowest-energy point is point 0, not a summary artifact.
    assert min(range(len(r["energies"])), key=lambda i: r["energies"][i]) == 0


def test_standard_orientation_still_parsed(tmp_path):
    """A conventional (reoriented) log with Standard orientation must still work —
    adding Input orientation must not regress the Standard-orientation path."""
    energies = [-100.0, -100.2]
    block = (
        "                         Standard orientation:                        \n"
        " ---------------------------------------------------------------------\n"
        " Center     Atomic      Atomic             Coordinates (Angstroms)\n"
        " Number     Number       Type             X           Y           Z\n"
        " ---------------------------------------------------------------------\n"
        "      1          6           0     0.000000    0.000000    0.000000\n"
        "      2          8           0     1.200000    0.000000    0.000000\n"
        " ---------------------------------------------------------------------\n"
    )
    log = tmp_path / "route_a_scan.log"
    log.write_text(
        " NAtoms=      2 NQM=        2\n"
        + "".join(
            f" SCF Done:  E(RB3LYP) = {e:.9f}     A.U. after   9 cycles\n"
            + block + "    -- Stationary point found.\n"
            for e in energies)
        + " Normal termination\n")

    r = parse_pes_scan_log(str(log), Console())
    assert r["success"] is True
    assert len(r["geometries"]) == 2
    assert r["energies"] == energies


if __name__ == "__main__":  # pragma: no cover
    import pytest
    raise SystemExit(pytest.main([__file__, "-v"]))
