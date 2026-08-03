#!/usr/bin/env python3
"""Apply mean-field state map to BPTI prmtop and verify net charge."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))

from proprep.pb_titrate.solvers.mean_field import make_site, solve
from proprep.pb_titrate.apply_state import build_production_prmtop, state_map_from_sites


BPTI_PRMTOP = HERE.parent / "bpti_asp3" / "bpti.prmtop"
BPTI_RST7   = HERE.parent / "bpti_asp3" / "bpti.rst7"


def main():
    sites = [
        make_site("AS4",  3, pH=7.0),
        make_site("TYR", 10, pH=7.0),
        make_site("LYS", 15, pH=7.0),
    ]
    result = solve(
        sites=sites,
        cluster_prmtop=BPTI_PRMTOP, cluster_rst7=BPTI_RST7,
        pH=7.0, max_iter=3, workers=3,
        work_root=HERE / "mf_work",
        pbsa_params={"epsin": 4.0, "space": 0.5, "nfocus": 2,
                     "bcopt": 5, "istrng": 150},
    )

    state_map = state_map_from_sites(result["sites"])
    info = build_production_prmtop(
        BPTI_PRMTOP, BPTI_RST7,
        state_map=state_map,
        output_prefix=HERE / "bpti_production",
    )

    print(f"\n=== Production prmtop built ===")
    print(f"  Input  net q: {info['net_q_in']:+.4f} ({info['n_atoms_in']} atoms)")
    print(f"  Output net q: {info['net_q_out']:+.4f} ({info['n_atoms_out']} atoms)")
    print(f"  Δq:           {info['delta_q']:+.4f}")
    for (rn, ri), chk in info["site_check"].items():
        print(f"  {rn}-{ri:<3}  actual={chk['actual']:+.4f}  "
              f"expected={chk['expected']:+.4f}  Δ={chk['delta']:+.4f}")
    print(f"\n  {info['output_prmtop']}")
    print(f"  {info['output_rst7']}")


if __name__ == "__main__":
    main()
