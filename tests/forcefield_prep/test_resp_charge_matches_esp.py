"""
RESP must be constrained to the charge its ESP was computed with.

A two-site run fitted site 1 against the wrong total charge:

    site 1  Gaussian: Charge = -1  Multiplicity = 11
            resp1.in:  -3  104
            result:    max charge +5.638, mean |q| 1.13

MetalSiteWorkflowManager.__init__ restores step_results from ONE workspace key,
``mcpb_step_results``, that every site shares. Whichever site ran last leaves its
step_1 there, so the mcpb-3 handler read site 2's charge (-3, multiplicity 1)
while processing site 1. Site 2 looked fine only because -3 was its own charge.

RESP does not fail on a mismatched constraint — it fits point charges to the
supplied potential subject to whatever total it is given, and distributes the
difference over the atoms. So the charge is now read from the site's OWN
Gaussian artifacts, which cannot be another site's.
"""

from pathlib import Path

import pytest

from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor

_esp = StructurePreprocessor._esp_charge_multiplicity


GJF = """%Chk=large_resp.chk
%Mem=48GB
%NProcShared=24
# B3LYP/def2SVP Opt Pop(MK,ReadRadii) IOp(6/33=2)

MCPB Large Model - Opt + ESP for RESP charges

-1 11
Fe   -1   0.000   0.000   0.000
S    -1   1.900   0.000   0.000

Fe 1.383

Fe 1.383

"""

LOG = """ Entering Gaussian System
 Charge = -1 Multiplicity =11
 Standard orientation:
 Normal termination of Gaussian 16
"""


def test_reads_the_charge_gaussian_actually_used(tmp_path):
    (tmp_path / "large_resp.log").write_text(LOG)

    assert _esp(tmp_path) == (-1, 11)


def test_falls_back_to_the_input_when_there_is_no_log(tmp_path):
    (tmp_path / "large_resp.gjf").write_text(GJF)

    assert _esp(tmp_path) == (-1, 11)


def test_the_log_wins_over_the_input(tmp_path):
    """The .gjf is what was asked for; the log is what ran."""
    (tmp_path / "large_resp.gjf").write_text(GJF.replace("-1 11", "-3 1"))
    (tmp_path / "large_resp.log").write_text(LOG)

    assert _esp(tmp_path) == (-1, 11)


def test_nothing_to_read_reports_unknown(tmp_path):
    """None must be distinguishable from a real charge of 0."""
    assert _esp(tmp_path) == (None, 1)


def test_a_malformed_charge_line_is_not_guessed(tmp_path):
    (tmp_path / "large_resp.gjf").write_text(GJF.replace("-1 11", "not a charge"))

    charge, _mult = _esp(tmp_path)
    assert charge is None


@pytest.mark.parametrize("charge,mult", [(-3, 1), (0, 1), (2, 6), (-1, 11)])
def test_round_trips_the_pairs_a_metal_site_uses(tmp_path, charge, mult):
    (tmp_path / "large_resp.log").write_text(
        f" Charge = {charge} Multiplicity ={mult}\n")

    assert _esp(tmp_path) == (charge, mult)


def test_two_sites_do_not_share_a_charge(tmp_path):
    """The reported failure: one dict, two sites, one charge."""
    site1 = tmp_path / "site_1" / "models"
    site2 = tmp_path / "site_2" / "models"
    site1.mkdir(parents=True)
    site2.mkdir(parents=True)
    (site1 / "large_resp.log").write_text(" Charge = -1 Multiplicity =11\n")
    (site2 / "large_resp.log").write_text(" Charge = -3 Multiplicity = 1\n")

    assert _esp(site1) == (-1, 11)
    assert _esp(site2) == (-3, 1)


def test_reading_is_per_directory_not_per_process(tmp_path):
    """Repeated calls must not cache the first site's answer."""
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    (a / "large_resp.log").write_text(" Charge = -1 Multiplicity =11\n")
    (b / "large_resp.log").write_text(" Charge = -3 Multiplicity = 1\n")

    assert [_esp(a), _esp(b), _esp(a)] == [(-1, 11), (-3, 1), (-1, 11)]
