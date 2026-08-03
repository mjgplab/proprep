"""
The large-model RESP route must not emit the contradictory 'Opt SP'.

Gaussian's Opt and SP are mutually exclusive job types; emitting both is
ambiguous. The large model is a geometry optimization whose Merz-Kollman ESP
(Pop=MK) is evaluated at the optimized geometry, so the route carries 'Opt'
(not 'Opt SP', and not 'SP').
"""

import re

from proprep.forcefield_prep.mcpb.qm_interface import (
    QMInterface, QMSoftware, QMCalculationMode,
)


def _route_line(gjf_text):
    for line in gjf_text.splitlines():
        if line.startswith("#"):
            return line
    return ""


def test_large_resp_route_has_opt_not_sp(tmp_path):
    pdb = tmp_path / "large.pdb"
    pdb.write_text(
        "ATOM      1  MN  MN  A 185      0.000   0.000   0.000  1.00  0.00          MN\n"
        "ATOM      2  O   WAT A 183      2.100   0.000   0.000  1.00  0.00           O\n"
        "END\n"
    )
    qm = QMInterface(QMSoftware.GAUSSIAN, QMCalculationMode.GUIDED)
    esp_keywords = "Pop(MK,ReadRadii) IOp(6/33=2,6/41=10,6/42=17)"
    qm.generate_input_files(
        pdb_file=pdb,
        output_dir=tmp_path,
        charge=0,
        multiplicity=11,
        functional="HF",
        basis_set="6-31G*",
        job_type="Opt",
        additional_keywords=esp_keywords,
        title_card="MCPB Large Model - Opt + ESP for RESP charges",
        output_name="large_resp",
        include_seminario_iop=False,
    )
    route = _route_line((tmp_path / "large_resp.gjf").read_text())

    # 'Opt' is present as a standalone keyword...
    assert re.search(r"(?<!\w)Opt(?!\w)", route), route
    # ...and there is no standalone 'SP' job keyword.
    assert not re.search(r"(?<!\w)SP(?!\w)", route), f"unexpected SP in route: {route}"
    # ESP fitting keywords survive.
    assert "Pop(MK,ReadRadii)" in route
    assert "IOp(6/33=2,6/41=10,6/42=17)" in route
