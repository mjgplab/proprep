#!/usr/bin/env python3
"""Coupling-precompute + MC test on BPTI 3-site.

Compares MC result to the mean-field result. For weakly-coupled sites
(BPTI is weakly coupled — surface-exposed residues, no strong pairs)
both should agree closely.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))

from proprep.pb_titrate.coupling import (
    make_coupling_site, compute_self_energies, compute_pair_couplings,
    coupling_summary,
)
from proprep.pb_titrate.solvers.monte_carlo import (
    run_mcmc, assign_states_from_mc, chemistry_label,
)
from proprep.pb_titrate.pb_backend import RT_LN10

PRM = HERE.parent / "bpti_asp3" / "bpti.prmtop"
RST = HERE.parent / "bpti_asp3" / "bpti.rst7"


def main():
    sites = [
        make_coupling_site("AS4",  3),
        make_coupling_site("TYR", 10),
        make_coupling_site("LYS", 15),
    ]
    print(f"=== BPTI 3-site coupling + MC test ===")
    for s in sites:
        ref_label = chemistry_label(s, s.ref_chem)
        chems = [chemistry_label(s, c) for c in s.chemistries]
        print(f"  {s.resname}-{s.resnum}: chemistries={chems}, ref={ref_label}")

    pbsa_params = {"epsin": 4.0, "space": 0.5, "nfocus": 2,
                    "bcopt": 5, "istrng": 150}

    print(f"\n--- Step 1: self-energies (one PB per site per chem) ---")
    se = compute_self_energies(
        sites, PRM, RST,
        workers=3, work_root=HERE / "coupling_work",
        pbsa_params=pbsa_params,
    )
    for site in sites:
        for chem in site.chemistries:
            ddg = se.intrinsic_ddG[(site.key, chem)]
            label = chemistry_label(site, chem)
            print(f"  {site.resname}-{site.resnum:<3}  {label}: ΔΔG_PB={ddg:+9.3f} kcal/mol")

    print(f"\n--- Step 2: pair couplings ---")
    W = compute_pair_couplings(
        sites, PRM, RST,
        workers=3, work_root=HERE / "coupling_work",
        pbsa_params=pbsa_params,
    )
    print(coupling_summary(W))

    print(f"\n--- Step 3: Metropolis MC ---")
    result = run_mcmc(se, W, pH=7.0,
                       n_steps=200_000, n_equil=20_000, seed=42)
    print(f"  Acceptance: {result.acceptance_ratio*100:.1f}%")
    print(f"  Marginals (post-equilibration):")
    for site in sites:
        marg = result.marginals[site.key]
        labeled = ", ".join(f"{chemistry_label(site, c)}={p:.4f}"
                             for c, p in sorted(marg.items()))
        print(f"    {site.resname}-{site.resnum:<3}  {labeled}")

    assign = assign_states_from_mc(result)
    print(f"\n  Dominant chemistry per site:")
    for site in sites:
        chem = assign[site.key]
        print(f"    {site.resname}-{site.resnum:<3}  {chemistry_label(site, chem)}")

    print(f"\n  Mean-field reference (from earlier run):")
    print(f"    AS4-3   DEPROT")
    print(f"    TYR-10  PROT")
    print(f"    LYS-15  PROT")


if __name__ == "__main__":
    main()
