"""WorkflowChecklist for the PB Titrate refinement step.

Invoked from the Topology Generator sub-menu (option `pb_titrate`). Eight
steps mirror the original titrate run.py orchestrator but bind to ProPrep's
WorkflowChecklist, prompts, and session recorder.

  pbt-1  Configure PB parameters
  pbt-2  Discover & confirm titratable sites
  pbt-3  Single-site intrinsic pKa
  pbt-4  Coupling precompute (optional)
  pbt-5  Solve (mean-field / Monte Carlo / exact enumeration)
  pbt-6  Review & override per-site recommendations
  pbt-7  Apply to prmtop (charges + ion rebalance)
  pbt-8  Persist recommendations + CSV report

Workspace contract
==================

Reads:
  prmtop, rst7                — from Topology Generator
  output_dir                  — base directory for pb_titrate_work/
  detected_redox_sites        — for envelope-aware multi-residue sites
  propka_pka_values           — optional, used as informational priors
  simulation_pH               — optional, default for pH prompt

Writes:
  pb_titrate_params           — configured PB parameters dict
  pb_titrate_sites            — serialized site list (resname/resnum/envelope)
  pb_titrate_self_energies    — per-site state energies pickled at output_dir
  pb_titrate_coupling         — coupling matrix pickled at output_dir
  pb_titrate_solver_result    — solver output (marginals, dominant states)
  pb_titrate_state_map        — final per-site assigned state (after review)
  prmtop_titrated, rst7_titrated  — output topology with charges baked in
  titrate_recommendations     — per-site recommendation dict consumed by
                                 the cpin generation step's "Use titrate
                                 recommendations" option
  titrate_report_csv          — human-readable per-site CSV report
"""
from __future__ import annotations

import csv
import glob
import os
import pickle
import traceback
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from proprep.utils.workflow_checklist import WorkflowChecklist, WorkflowStep
from proprep.utils.prompts import (
    prompt_with_context, confirm_with_context,
    int_prompt_with_context, float_prompt_with_context,
)


MODULE_NAME = "PB Titrate"


# ============================================================================
# Educational content (shown via [i] keybind in the checklist)
# ============================================================================

PB_TITRATE_OVERVIEW = """\
[bold]What this module does[/bold]
PB Titrate predicts the dominant protonation state of every titratable site
in your prepared topology at a chosen pH. It returns per-site recommendations
that can be (a) baked directly into a fixed-charge prmtop for production MD,
or (b) handed to constant-pH MD as initial states (consumed by the cpin
generation step in the Topology Generator).

The pipeline is the [italic]titrate[/italic] tool, originally written as a standalone
ProPrep-adjacent project and then vendored into ProPrep as the [italic]pb_titrate[/italic]
package. Multi-residue redox cofactors (hemes, FeS clusters) are handled
through a [italic]Site[/italic] abstraction so their envelopes stay integer-charge during
PB calculations.

[bold]Theory: Bashford thermodynamic cycle[/bold]
For each titratable site we compare two protonation free energies:

  ΔG[italic]_protein[/italic]  = E_DEPROT − E_PROT, with the site embedded in the protein cluster
  ΔG[italic]_model[/italic]    = the same difference for an ACE–X–NME tripeptide in solvent

  ΔpKa = (ΔG[italic]_protein[/italic] − ΔG[italic]_model[/italic]) / (RT ln10)
  pKa[italic]_int[/italic] = pKa[italic]_model[/italic] + ΔpKa

The "model" leg has the same Born and grid artefacts as the "protein" leg, so
they cancel in the difference. This is the four-call cycle implemented in
[italic]intrinsic.py[/italic]: cluster_PROT, cluster_DEPROT, model_PROT, model_DEPROT.

When other titratable sites are also non-reference, their charges contribute a
pair interaction W[italic]_ij[/italic] to the site's effective pKa. Step 4 precomputes the
W[italic]_ij[/italic] matrix; step 5's solver (mean-field / MC / enumeration) folds it in.

PB calculations use AmberTools' [italic]pbsa[/italic] backend with the linearized PB equation
(npbopt=0); ionic strength enters through Debye-Hückel screening of the
boundary potential.

[bold]Ionizable residues handled[/bold]
  AS4 (Asp), GL4 (Glu), HIP (His; HID/HIE/HIP triad), LYS, CYS, TYR
  PRN (heme propionate); other redox-cofactor envelopes via the Site abstraction

[bold]Implemented solvers (step 5)[/bold]
  • [bold]Mean-field[/bold] iteration — each site sees the [italic]average[/italic] field of every other
    site. Iterates until no site changes its assigned state, or max_iter
    is reached. Fast and deterministic; cheapest path when coupling is weak.
    Live-PB variant runs without a precomputed W matrix.
  • [bold]Monte Carlo (Metropolis)[/bold] — single-site flips against the precomputed
    energy function (self-energies + W). Defaults: 100 000 steps, 10 000
    equilibration. Correct in the sampling limit; handles correlated flips.
  • [bold]Exact enumeration[/bold] — sums all chemistry-state combinations exactly.
    Cost is the product of per-site chemistry counts, so it's only practical
    when that product is small.

  Auto-recommendation ([italic]select_solver.recommend()[/italic]) uses three quantities:
  the largest |W[italic]_ij[/italic]| in the matrix, the count of "near-pKa" sites (any two
  chemistries within 1 kT at the target pH), and the total state-combination
  count. Rules: ≤16 total combinations → enumerate; max |W| > 2 kT with at
  least one near-pKa site, or ≥5 near-pKa sites → monte_carlo; otherwise
  → mean_field. You can override the recommendation.

[bold]Multi-residue envelopes[/bold]
For heme propionates and redox cofactors, the titrating residue sits inside an
envelope of co-bonded residues that must move together to keep the cluster's
net charge integer (a hard PB requirement). The Site object carries the
envelope, and per-site PB clusters are built around the union of the active
site's envelope plus the [italic]other-redox[/italic] envelopes held at their reference state.

[bold]The 8-step pipeline[/bold]
  1. Configure PB parameters (εin, ionic strength, grid, target pH, cluster radius)
  2. Discover sites (titratable residues + envelopes from detected_redox_sites)
  3. Single-site pKa (4-state Bashford cycle per site — the heaviest step)
  4. Coupling precompute (optional; needed for MC / enumerate)
  5. Solve (mean-field / MC / enumeration — auto-pick or override)
  6. Review & override (per-site table; flag near-50/50 sites)
  7. Apply to prmtop (bake state charges, optionally rebalance ions)
  8. Persist recommendations (workspace dict + CSV; cpin step consumes these)

[bold]Defaults[/bold]
  εin = 4, ionic strength = 150 mM, grid spacing = 0.5 Å, focusing levels = 2,
  parallel workers = 4, target pH = 7.0, cluster radius = 25 Å. Higher εin
  damps both self-energy and pair coupling; finer grid tightens accuracy at
  cubic cost; cluster radius caps each per-site PB job size.

[bold]Constant-pH MD handoff: a note on net charge[/bold]
If you forward these recommendations to Amber's discrete-state CpHMD via the
cpin step, two non-neutralities arise:
  • [italic]Static t=0 offset[/italic]: tleap+leaprc.constph neutralized the system for
    cpinutil's default starting state (state 0). Setting initial states to PB
    Titrate's recommendations changes the net charge at t=0 by Δq = (sum of
    recommended charges) − (sum of state-0 charges).
  • [italic]Dynamic during MD[/italic]: every accepted Metropolis flip changes net charge
    by ±1 (or ±2 for HIP↔HID).
Amber's discrete-state CpHMD (Swails et al. 2014) handles this not with a
co-titrating ion mechanism but by a [italic]hybrid scheme[/italic]: regular MD runs in
explicit solvent with PME, but the MC accept/reject decision is evaluated
in GB implicit solvent (where neutrality is irrelevant). After an accepted
move, the solute is held fixed for [italic]ntrelax[/italic] solvent-only steps so waters
reorganize around the new charge. Net charge during explicit propagation
is corrected by PME's uniform background charge — small offsets (±a few
electrons) are tolerated; this is an accepted limitation of the method.
There is no co-titrating-ion or [italic]cpein[/italic] file in the discrete-state pipeline.
Strict unit-cell neutrality with PME requires [italic]continuous[/italic] CpHMD
(λ-dynamics, with titratable water/ions) — a different code path that does
not consume cpinutil's cpin file.

For implicit-solvent CpHMD ([italic]igb[/italic]), there's no Ewald constraint at all and
the offset is a non-issue.

[bold]Background reading on PB-pKa methodology[/bold]
The pb_titrate code itself is not a citation of any of these papers, but they
are the canonical references for the underlying methods:
  Bashford & Karplus, Biochemistry 29, 10219 (1990) — thermodynamic cycle
  Yang, Gunner, Sampogna, Sharp & Honig, Proteins 15, 252 (1993) — FDPB pKas
  Antosiewicz, McCammon & Gilson, J. Mol. Biol. 238, 415 (1994) — predictions
  Ullmann & Knapp, Eur. Biophys. J. 28, 533 (1999) — review of solvers
  Swails, York & Roitberg, JCTC 10, 1341 (2014) — Amber discrete-state explicit-solvent CpHMD
"""


# ============================================================================
# Step definitions
# ============================================================================

PB_TITRATE_STEPS: List[WorkflowStep] = [
    WorkflowStep(
        id="pbt-1", name="Configure PB Parameters",
        description="εin, ionic strength, grid spacing, workers, target pH",
        handler="_checklist_pbt_1_configure",
        section="PB Titrate",
        info=("[bold]Purpose[/bold]\n"
              "Pin the physics knobs that drive every downstream PB call.\n\n"
              "[bold]Parameters[/bold]\n"
              "  • [italic]εin[/italic] — interior protein dielectric. The default is 4, "
              "matching AMBER's constph convention. Higher εin damps both "
              "self-energy and pair coupling.\n"
              "  • [italic]Ionic strength[/italic] (mM) — sets Debye screening. The "
              "Debye length at 25 °C, εout=80 is κ⁻¹ ≈ 3.04 / √I (Å, with I in "
              "mol/L); equivalently κ ≈ 0.328 √I Å⁻¹. 150 mM is physiological. "
              "PB enters this through the linearized Poisson-Boltzmann equation "
              "(npbopt=0).\n"
              "  • [italic]Grid spacing[/italic] (Å) — finer = more accurate but cubic "
              "in cost. 0.5 Å is the production default.\n"
              "  • [italic]Focusing levels[/italic] (nfocus) — number of nested grids "
              "from coarse-box → site-centred fine grid. 2 is standard.\n"
              "  • [italic]bcopt[/italic] — outer boundary condition. Per the AmberTools "
              "manual, valid values are 1 (zero-potential / conductor), 5 "
              "(Debye-Hückel from grid charges; default), 6 (DH with charge-"
              "singularity removal), and 10 (periodic). pb_titrate uses 5.\n"
              "  • [italic]Workers[/italic] — parallel pbsa processes. ≤8 is generally "
              "safe; the binding constraint is RAM, not CPU count.\n"
              "  • [italic]Target pH[/italic] — the pH at which states will be "
              "recommended. Comes from simulation_pH if previously set.\n"
              "  • [italic]Cluster radius[/italic] (Å) — per-site PB cluster cutoff "
              "around the titrating residue. Empty = full structure (slow). "
              "25 Å is the default."),
    ),
    WorkflowStep(
        id="pbt-m", name="Minimize Structure",
        description="Relieve bad contacts before PB (sander/pmemd); optional but recommended",
        handler="_checklist_pbt_m_minimize",
        section="PB Titrate",
        optional=True,
        checkpoint=True,
        checkpoint_message=("Runs an Amber energy minimization locally and "
                            "blocks until it finishes. Resumable."),
        info=("[bold]Purpose[/bold]\n"
              "Run a short energy minimization with an Amber MD engine to "
              "relieve steric clashes [italic]before[/italic] any PB calculation. PB pKa "
              "shifts are very sensitive to bad contacts: an overlapping atom "
              "pair produces a spurious reaction-field energy that contaminates "
              "every leg of the Bashford cycle (pbt-3/pbt-4). A structure fresh "
              "from tLEaP — especially one that was never equilibrated — should "
              "almost always be minimized here first.\n\n"
              "[bold]What it does[/bold]\n"
              "Minimizes the full [italic]solvated[/italic] system (waters/ions included, "
              "constant-volume PBC when a box is present), then hands the "
              "relaxed coordinates to pbt-2, which strips solvent for the PB "
              "calls and keeps the solvated copy for state apply (pbt-7). The "
              "topology (prmtop) is unchanged; only the coordinates (rst7) "
              "move.\n\n"
              "[bold]Defaults[/bold]\n"
              "Mirror the bundled aqueous-metalloprotein minimization template "
              "(maxcyc=5000, ncyc=1000, cut=10 Å, backbone+metal restraints at "
              "10 kcal/mol·Å² on @CA,C,O,N). The prompts let you override max "
              "cycles, the SD→CG switch, the cutoff, and the restraint mask / "
              "weight, and pick the engine (sander / pmemd / pmemd.MPI).\n\n"
              "[bold]Restraints[/bold]\n"
              "Restraining heavy backbone atoms lets clashing side chains and "
              "hydrogens relax without dragging the structure away from the "
              "built/experimental conformation you intend to titrate. Set "
              "ntr=no for a fully free minimization if you want maximal strain "
              "relief.\n\n"
              "[bold]Skipping[/bold]\n"
              "Optional: skip it if your input was already equilibrated. The "
              "step runs locally and blocks until the engine exits; on a "
              "cluster, minimize beforehand and skip this step."),
    ),
    WorkflowStep(
        id="pbt-2", name="Discover Sites",
        description="Find titratable residues; envelope multi-residue sites from detected_redox_sites",
        handler="_checklist_pbt_2_discover",
        section="PB Titrate",
        dependencies=["pbt-1", "pbt-m"],
        info=("[bold]Purpose[/bold]\n"
              "Walk the prmtop, identify every titratable residue, and group "
              "co-bonded residues into [italic]envelopes[/italic] when they belong to "
              "a multi-residue redox cofactor.\n\n"
              "[bold]Recognised residues[/bold]\n"
              "  AS4, GL4 (constph aspartate / glutamate)\n"
              "  HIP (histidine; HID ⇌ HIE ⇌ HIP triad)\n"
              "  LYS, CYS, TYR\n"
              "  PRN (heme propionate)\n\n"
              "[bold]Why envelopes matter[/bold]\n"
              "PB requires every cluster to carry an integer net charge. A "
              "lone PRN protonation flips charge by ±1, but a heme also has "
              "Fe redox states that change the [italic]cofactor[/italic] charge. The Site "
              "abstraction binds the titrating residue to its envelope of "
              "co-bonded residues (Cys-S–Fe links, propionate–porphyrin "
              "covalency, etc.) so they move as one chemical entity.\n\n"
              "[bold]Inputs read[/bold]\n"
              "  detected_redox_sites — supplies envelope topology for redox "
              "cofactors. If empty, only standard amino acids titrate.\n\n"
              "[bold]Exclusion[/bold]\n"
              "Use the prompt to drop sites you know aren't physically "
              "relevant (e.g., a buried residue you've already locked from "
              "external pKa data). Format: AS4:34,LYS:120."),
    ),
    WorkflowStep(
        id="pbt-3", name="Single-Site pKa",
        description="4-state Bashford cycle per site (PB calls, parallel-friendly)",
        handler="_checklist_pbt_3_single_site",
        section="PB Titrate",
        dependencies=["pbt-2"],
        checkpoint=True,
        checkpoint_message=("Long-running PB calculations — one set of "
                              "calls per site. Resumable."),
        info=("[bold]Purpose[/bold]\n"
              "Compute the [italic]intrinsic[/italic] pKa of every site — the pKa it would "
              "have in the protein if all other titratable sites were held at "
              "their model-compound reference state. This is the first piece "
              "of the system energy model that step 5 will eventually solve "
              "(see step 4's info for the full equation).\n\n"
              "[bold]Bashford 4-state cycle (per site)[/bold]\n"
              "Four PB calls give the energies needed to close the thermodynamic "
              "box:\n"
              "  (a) protonated,   in protein cluster        → E_PROT,protein\n"
              "  (b) deprotonated, in protein cluster        → E_DEPROT,protein\n"
              "  (c) protonated,   in solvent (model compound) → E_PROT,model\n"
              "  (d) deprotonated, in solvent (model compound) → E_DEPROT,model\n"
              "  ΔΔG = (b − a) − (d − c)\n"
              "  ΔpKa = ΔΔG / (RT ln10)\n"
              "  pKa  = pKa[italic]_model[/italic] + ΔpKa\n"
              "Differencing protein and model legs cancels Born and grid "
              "discretisation artefacts that are common to both.\n\n"
              "[bold]What this produces for the energy model[/bold]\n"
              "ΔE[italic]_self,i[/italic](s[italic]_i[/italic]) = [(E[italic]_protein[/italic] − E[italic]_model[/italic])][italic]_state s_i[/italic] — "
              "the protein-environment shift of site [italic]i[/italic] in chemistry s[italic]_i[/italic], "
              "with every other titratable site at its reference (closest-to-"
              "neutral) chemistry. Step 5 uses these as the diagonal "
              "(per-site) part of the system free energy.\n\n"
              "[bold]Multi-state residues[/bold]\n"
              "HIS uses three chemistries (HID, HIE, HIP); the routine "
              "computes pairwise transfer energies and stores all of them. "
              "PRN/heme cofactors emit a multi-state energy table the solver "
              "consumes directly.\n\n"
              "[bold]Cost[/bold]\n"
              "Roughly 4 × N[italic]_sites[/italic] PB calls (more for multi-state residues), "
              "parallelised across [italic]workers[/italic] processes.\n\n"
              "[bold]Resumable[/bold]\n"
              "Per-site work directories under pb_titrate_work/single_site/; "
              "results pickled to single_site_pka.pkl. Failed sites are "
              "logged but don't abort the run."),
    ),
    WorkflowStep(
        id="pbt-4", name="Coupling Precompute",
        description="Pair coupling matrix W_ij (optional; for coupled solvers)",
        handler="_checklist_pbt_4_coupling",
        section="PB Titrate",
        dependencies=["pbt-3"],
        optional=True,
        checkpoint=True,
        checkpoint_message=("Pair PB calls within per-pair clusters; "
                              "Debye-derived distance cutoff bounds work."),
        info=("[bold]Purpose[/bold]\n"
              "Build the pair-interaction matrix W[italic]_ij[/italic] giving the free-energy "
              "cost of flipping site [italic]i[/italic]'s state when site [italic]j[/italic] is also non-"
              "reference. Together with the self-energies from step 3, this "
              "completes the analytical system free-energy model that step 5 "
              "evaluates.\n\n"
              "[bold]Full system free-energy model[/bold]\n"
              "Given a configuration s = (s[italic]_1[/italic], …, s[italic]_N[/italic]) of chemistries:\n"
              "  [bold]G(s) = Σ[/bold][italic]_i[/italic] ΔE[italic]_self,i[/italic](s[italic]_i[/italic])\n"
              "       + [bold]Σ[/bold][italic]_{i<j}[/italic] W[italic]_ij[/italic](s[italic]_i[/italic], s[italic]_j[/italic])\n"
              "       + RT ln10 · [bold]Σ[/bold][italic]_i[/italic] (n[italic]_i[/italic](s[italic]_i[/italic]) · pH − pKa[italic]_corr[/italic](s[italic]_i[/italic]))\n"
              "ΔE[italic]_self,i[/italic] comes from step 3 (Bashford 4-state cycle: protein-"
              "vs-model environmental shift). W[italic]_ij[/italic] comes from this step. "
              "The pH-dependent term is the Henderson-Hasselbalch contribution: "
              "n[italic]_i[/italic] is the number of titratable protons in state s[italic]_i[/italic] and "
              "pKa[italic]_corr[/italic] is the calibration constant from cpinutil's "
              "residue templates.\n\n"
              "[bold]Why couplings? Why not just self-energies?[/bold]\n"
              "ΔE[italic]_self,i[/italic] is computed with [italic]every other site held at reference[/italic]. "
              "It captures site [italic]i[/italic]'s individual preference but not how site [italic]i[/italic]'s "
              "state changes site [italic]j[/italic]'s energy. Two charged residues 8 Å "
              "apart electrostatically interact; flipping one shifts the "
              "other's preference. W[italic]_ij[/italic] is the cross-term that captures this.\n\n"
              "[bold]Second-difference formula for W_ij[/bold]\n"
              "For each chemistry pair (s[italic]_i[/italic], s[italic]_j[/italic]), the coupling is the "
              "joint-perturbation energy minus the two single-site "
              "perturbations plus the reference (all-at-ref) energy:\n"
              "  W[italic]_ij[/italic](s, t) = E(i=s, j=t)\n"
              "             − E(i=s, j=ref)\n"
              "             − E(i=ref, j=t)\n"
              "             + E(i=ref, j=ref)\n"
              "Each of those four E(...) is one PB call inside the same "
              "cluster. The second difference cancels Born terms, grid "
              "artefacts, and any cluster-specific contributions that are "
              "shared by all four energies, leaving only the joint i↔j "
              "interaction. In linearized PB this has a clean physical "
              "meaning — the work to charge site [italic]j[/italic] against the field site "
              "[italic]i[/italic] produces — so like-charged neighbors give positive W "
              "(penalty for both being charged) and opposite-charged "
              "neighbors give negative W (favorable to flip both).\n\n"
              "[bold]Why 'single' PB calls appear in the progress log[/bold]\n"
              "The four-energy formula needs the two singles E(i=s, j=ref) "
              "and E(i=ref, j=t) in addition to the pair E(i=s, j=t). With a "
              "cluster cutoff, the singles must be re-computed [italic]inside each "
              "pair's cluster[/italic] (different cluster = different PB problem; the "
              "formula's cancellation only works when all four energies share "
              "the same cluster boundary). Without a cutoff (full structure), "
              "the singles are the same number step 3 already computed and "
              "could be reused — there's a known efficiency gap on that path "
              "(~10% of Phase 2's PB calls are redundant in the no-cutoff "
              "case).\n\n"
              "[bold]Distance cutoff[/bold]\n"
              "Pairs with site centroids beyond a Debye-derived cutoff "
              "contribute negligibly and are skipped: "
              "[italic]recommended_pair_cutoff(I) = max(6 Å, 3·κ⁻¹)[/italic]. At I=150 mM "
              "this is ≈ 24 Å; at I=650 mM it tightens to ≈ 11 Å. Skipped "
              "pairs get W_ij = 0.\n\n"
              "[bold]When to skip[/bold]\n"
              "Skip if you're going to run [italic]mean-field[/italic] anyway. The live-PB "
              "mean-field path doesn't need a precomputed W matrix.\n\n"
              "[bold]When to run[/bold]\n"
              "Multi-heme proteins, salt bridges, or any system where you "
              "expect cooperative protonation (e.g., carboxyl–carboxyl pairs "
              "sharing a proton). Required for Monte Carlo and exact "
              "enumeration solvers.\n\n"
              "[bold]Cost[/bold]\n"
              "With a cluster cutoff: 4 PB calls per kept pair within that "
              "pair's own cluster, parallelised across [italic]workers[/italic] processes. "
              "Without a cutoff: 1 reference + (Σ_i n[italic]_nonref(i)[/italic]) singles + "
              "(Σ_{i<j} n[italic]_nonref(i)[/italic]·n[italic]_nonref(j)[/italic]) pairs, all sharing one "
              "cluster."),
    ),
    WorkflowStep(
        id="pbt-5", name="Solve",
        description="Auto-pick mean-field / Monte Carlo / exact enumeration",
        handler="_checklist_pbt_5_solve",
        section="PB Titrate",
        dependencies=["pbt-3"],
        info=("[bold]Purpose[/bold]\n"
              "Evaluate the system free-energy model built in steps 3+4 to "
              "decide each site's assigned state at the target pH. All four "
              "solvers operate on the same model; they differ only in how "
              "they handle correlations between sites.\n\n"
              "[bold]The energy model[/bold]\n"
              "For any configuration s = (s[italic]_1[/italic], …, s[italic]_N[/italic]):\n"
              "  G(s) = Σ[italic]_i[/italic] ΔE[italic]_self,i[/italic](s[italic]_i[/italic])         [from step 3]\n"
              "       + Σ[italic]_{i<j}[/italic] W[italic]_ij[/italic](s[italic]_i[/italic], s[italic]_j[/italic])    [from step 4]\n"
              "       + RT ln10 · Σ[italic]_i[/italic] (n[italic]_i[/italic] · pH − pKa[italic]_corr,i[/italic])   [Henderson-Hasselbalch]\n"
              "Boltzmann population:\n"
              "  P(s) = exp(−G(s)/kT) / Z,  with Z = Σ[italic]_s[/italic] exp(−G(s)/kT)\n"
              "Per-site marginal in chemistry c:\n"
              "  P[italic]_i[/italic](c) = Σ[italic]_{s : s_i = c}[/italic] P(s)\n"
              "The solver's job is to estimate the marginals (or pick the "
              "dominant state per site).\n\n"
              "[bold]Mean-field iteration[/bold]\n"
              "Each site sees the [italic]average[/italic] field of every other site, then "
              "picks the lowest-free-energy chemistry given that field. The "
              "iteration repeats until [italic]no site changes its assigned state[/italic], "
              "or max_iter is reached (default 5). Equivalent to factoring "
              "P(s) ≈ Π[italic]_i[/italic] P[italic]_i[/italic](s[italic]_i[/italic]) — sites assumed independent.\n"
              "If you ran step 4 (precomputed W matrix), step 5 currently uses "
              "a one-shot adapter that picks the dominant chemistry per site "
              "from the precomputed self-energies and pH term — it does not "
              "iterate against the static W matrix. The full iterative "
              "mean-field path runs without step 4 by issuing live PB calls "
              "every iteration ([italic]solvers/mean_field.py[/italic]).\n"
              "  + Cheap, deterministic\n"
              "  − Misses correlations a static W matrix could capture\n\n"
              "[bold]Cluster mean-field[/bold]\n"
              "Threshold the coupling graph at |W| > [italic]threshold_kT[/italic] · kT "
              "(default 1 kT) and take connected components as 'clusters'. "
              "Within each cluster, [italic]enumerate[/italic] every joint chemistry-state "
              "combination exactly. Between clusters, treat coupling at MF "
              "level. Iterate per-cluster marginals to self-consistency. "
              "Frustrated pairs (where coupling overrides individual "
              "preference) all end up inside a single cluster, so the "
              "anti-correlated joint distribution is captured exactly.\n"
              "  + Exact within strongly-coupled subsystems\n"
              "  + Deterministic; converges in a few iterations\n"
              "  + Faster than MC on most protein systems\n"
              "  − Requires the coupling matrix from step 4\n"
              "  − Aggressive thresholds can produce a single huge cluster "
              "(safety bound: max_cluster_states = 65 536 = 2¹⁶)\n\n"
              "[bold]Monte Carlo (Metropolis)[/bold]\n"
              "Single-site-flip MC against the precomputed G(s). At each "
              "step, propose a random chemistry for a random site; accept "
              "with probability min(1, exp(−ΔG/kT)). Defaults: 100 000 steps, "
              "10 000 equilibration. Reports acceptance ratio + per-site "
              "lockup. Requires the coupling matrix from step 4.\n"
              "  + Exact in the sampling limit; handles correlated flips\n"
              "  − Stochastic; tightly coupled clusters need long runs\n\n"
              "[bold]Exact enumeration[/bold]\n"
              "Sums G(s) over every configuration s exactly. Cost is "
              "n[italic]_total[/italic] = Π[italic]_i[/italic] n[italic]_chemistries(i)[/italic]. For N binary sites that is "
              "2^N; multi-state residues like HIS multiply by their own "
              "chemistry count. Feasible up to n[italic]_total[/italic] ≈ 100,000 on "
              "commodity hardware. Requires step 4.\n"
              "  + Gold-standard answer when feasible — no sampling noise, "
              "no convergence question, no acceptance ratio\n"
              "  − Exponential cost in n_total\n\n"
              "[bold]Auto-recommend[/bold] ([italic]select_solver.recommend()[/italic])\n"
              "  • [italic]n_total[/italic] ≤ ENUMERATE_MAX_STATES (10⁵) → [bold]enumerate[/bold] "
              "(preferred whenever feasible — strictly better than MC)\n"
              "  • frustrated pair exists → [bold]cluster_mean_field[/bold]\n"
              "  • max |W[italic]_ij[/italic]| > 2 kT and ≥1 near-pKa site → [bold]monte_carlo[/bold]\n"
              "  • ≥5 near-pKa sites → [bold]monte_carlo[/bold]\n"
              "  • otherwise → [bold]mean_field[/bold]\n"
              "A site is \"near-pKa\" when two of its chemistries are within "
              "1 kT at the target pH. A pair is \"frustrated\" when |W| "
              "exceeds the smaller site's individual preference margin "
              "(so coupling can flip a site against its own intrinsic "
              "preference). You can override the recommendation."),
    ),
    WorkflowStep(
        id="pbt-6", name="Review & Override",
        description="Per-site table; near-50/50 sites flagged; allow user overrides",
        handler="_checklist_pbt_6_review",
        section="PB Titrate",
        dependencies=["pbt-5"],
        info=("[bold]Purpose[/bold]\n"
              "Human-in-the-loop quality control before you commit a state "
              "assignment to a prmtop or cpin.\n\n"
              "[bold]What to look for[/bold]\n"
              "  • [italic]Near-degenerate sites[/italic] — fractions in (0.4, 0.6). The "
              "calculation can't decide; experimental pKa data or chemical "
              "intuition should override.\n"
              "  • [italic]Sites near catalytic residues[/italic] — even small numerical "
              "errors can flip an active-site residue. Cross-check against "
              "PROPKA, H++, or experimental pKa values.\n"
              "  • [italic]Buried CYS / TYR[/italic] — often anomalously low pKa if "
              "stabilised by a positive partner; if you don't expect "
              "deprotonation in the experimental conditions, override.\n"
              "  • [italic]Heme propionates (PRN)[/italic] — protonation depends on "
              "neighboring oxidation states; revisit after confirming Fe "
              "redox assignment.\n\n"
              "[bold]Override format[/bold]\n"
              "  RESNAME:RESNUM=STATE_NAME\n"
              "State names come from [italic]cpinutil.py --describe[/italic] and are literal "
              "strings like \"STATE 0\", \"STATE 1\" — not residue codes. AS4 "
              "has 5 states, HIP has 3. Run [italic]cpinutil.py --describe AS4[/italic] (or GL4, "
              "HIP, etc.) to see the per-state Prot Cnt and pKa Corr columns "
              "and pick the right one. Example:\n"
              "  AS4:34=STATE 0\n"
              "  HIP:120=STATE 2\n"
              "Empty input ends the override loop."),
    ),
    WorkflowStep(
        id="pbt-7", name="Apply to Prmtop",
        description="Bake state charges into prmtop; rebalance ions",
        handler="_checklist_pbt_7_apply",
        section="PB Titrate",
        dependencies=["pbt-6"],
        info=("[bold]Purpose[/bold]\n"
              "Translate per-site state assignments into a [italic]fixed-charge[/italic] "
              "production prmtop suitable for vanilla MD. (For constant-pH "
              "MD, you keep the constph residues and use the recommendations "
              "as initial states via the cpin step instead.)\n\n"
              "[bold]Mechanics[/bold]\n"
              "  • Look up canonical partial charges for each (resname, "
              "state) pair from the residues registry\n"
              "  • Overwrite atom CHARGE entries in the prmtop\n"
              "  • Optionally [italic]rebalance bulk ions[/italic] to keep net charge integer\n\n"
              "[bold]Ion rebalance[/bold]\n"
              "If the per-residue state changes shift the system's net charge "
              "by Δq, the rebalancer [italic]removes[/italic] Δq bulk counter-ions to restore "
              "neutrality — it never adds new ions (avoiding fabricated "
              "prmtop entries and collision-free placement). Within a species, "
              "ions farthest from the protein heavy atoms are removed first. "
              "Default ranking prefers monovalent species (Na⁺/K⁺) for fine "
              "control; multivalents fall back only if the +1 pool is "
              "exhausted. If the residual cannot be closed exactly with the "
              "available species, it raises rather than silently leaving an "
              "imbalance. Without rebalance, PME runs on a non-neutral cell "
              "and you'll see drift.\n\n"
              "[bold]Outputs[/bold]\n"
              "  prmtop_titrated, rst7_titrated — written next to "
              "pb_titrate_work/, with `_titrated` suffix on the original "
              "prmtop stem. Original files are not overwritten."),
    ),
    WorkflowStep(
        id="pbt-pdb", name="Write modern-FF PDB",
        description="Bake recommended states into a PDB (standard names) for a modern-FF rebuild",
        handler="_checklist_pbt_pdb_write",
        section="PB Titrate",
        optional=True,
        dependencies=["pbt-5"],
        info=("[bold]Purpose[/bold]\n"
              "For [italic]conventional[/italic] fixed-charge MD (not constant-pH), bake the "
              "recommended protonation states into a PDB using [italic]standard[/italic] Amber "
              "residue names, then rebuild the topology with a modern protein "
              "force field via the Topology Generator.\n\n"
              "[bold]Why not just use the titrated prmtop (step 7)?[/bold]\n"
              "PB-titrate works on a constant-pH topology sourced from "
              "oldff/leaprc.ff10, which lacks the ff14SB/ff19SB backbone "
              "dihedral and CMAP corrections. Step 7 bakes charges into that "
              "ff10 prmtop — correct for [italic]constant-pH[/italic] MD (with a cpin). For "
              "[italic]fixed-protonation[/italic] MD you want a modern FF, which means "
              "rebuilding from a PDB. This step produces that PDB.\n\n"
              "[bold]What it does[/bold]\n"
              "Renames every titratable residue in the minimized, "
              "solvent-stripped structure (Ca²⁺ + cofactors retained) to the "
              "standard residue encoding its state:\n"
              "  AS4→ASP/ASH, GL4→GLU/GLH, HIP→HID/HIE/HIP, LYS→LYS/LYN, "
              "CYS→CYS/CYM.\n"
              "Metal-coordinating residues excluded in step 2 are written at "
              "their fixed metal-bound state (deprotonated carboxylate → ASP/"
              "GLU). Writes pb_titrate_work/<stem>_renamed.pdb and sets the "
              "workspace key [italic]pb_rename_pdb_file[/italic], which the Topology "
              "Generator auto-selects with highest priority.\n\n"
              "[bold]Then[/bold]\n"
              "Re-enter the Topology Generator, choose a modern protein FF "
              "(ff19SB → OPC water, or ff14SB → TIP3P) and explicit solvation; "
              "tleap re-solvates fresh with your minimized coordinates and the "
              "baked-in protonation states.\n\n"
              "[bold]Tyrosinate caveat[/bold]\n"
              "Neither ff14SB nor ff19SB has a deprotonated-tyrosine residue. "
              "A Tyr predicted deprotonated is written as [italic]neutral[/italic] TYR with a "
              "prominent warning; supply custom TYM parameters if you need a "
              "true tyrosinate."),
    ),
    WorkflowStep(
        id="pbt-8", name="Persist Recommendations",
        description="Workspace dict + CSV report (consumed by cpin step)",
        handler="_checklist_pbt_8_persist",
        section="PB Titrate",
        # Persisting consumes the per-site state map produced by review
        # (pbt-6) — it does NOT need the titrated prmtop (pbt-7), which is an
        # independent branch (fixed-charge MD vs. cpin handoff). Depending on
        # pbt-7 wrongly warned when users went review → persist for the cpin
        # path without applying to prmtop.
        dependencies=["pbt-6"],
        info=("[bold]Purpose[/bold]\n"
              "Stash a structured record of every recommendation so it can "
              "(a) be consumed automatically by the cpin generation step, "
              "and (b) be inspected/audited later by humans.\n\n"
              "[bold]Workspace key[/bold]\n"
              "  titrate_recommendations: dict keyed by (resname, resnum) "
              "with state_id, state_name, prot_count, pka_corr, net_charge.\n\n"
              "[bold]How cpin uses it[/bold]\n"
              "Topology Generator option 6 (cpin generation) offers \"Use "
              "PB Titrate recommendations\" which translates this dict into "
              "the [italic]-states[/italic] argument list for cpinutil, so constant-pH MD "
              "starts from sensible per-site initial states rather than the "
              "default (all-protonated or all-deprotonated).\n\n"
              "[bold]CSV report[/bold]\n"
              "  pb_titrate_work/titrate_recommendations.csv — one row per "
              "site, suitable for spreadsheet review and inclusion in "
              "supplementary information.\n\n"
              "[bold]CpHMD handoff caveat[/bold]\n"
              "If you push these states into a cpin and run discrete-state "
              "CpHMD, the system will be non-neutral at t=0 by Δq = (sum of "
              "recommended charges) − (sum of cpinutil-default state-0 "
              "charges), and will fluctuate further as MC moves are accepted. "
              "Amber's discrete-state CpHMD (Swails 2014) does the MC step "
              "in GB and tolerates the residual non-neutrality during "
              "explicit-solvent propagation via PME's background-charge "
              "correction — there is no [italic]cpein[/italic]/co-titrating-ion mechanism. "
              "Strict neutrality with PME requires continuous (λ-dynamics) "
              "CpHMD, which does not read this cpin. See the overview "
              "([italic]i[/italic] from the menu) for details."),
    ),
    WorkflowStep(
        id="pbt-cmp", name="Compare pKa with ProPKA",
        description="Re-run ProPKA on the minimized structure; PB-vs-ProPKA table + CSV",
        handler="_checklist_pbt_cmp_propka",
        section="PB Titrate",
        optional=True,
        dependencies=["pbt-3"],
        info=("[bold]Purpose[/bold]\n"
              "Cross-check the PB single-site pKas against ProPKA as an "
              "independent empirical predictor.\n\n"
              "[bold]Why it re-runs ProPKA[/bold]\n"
              "PB-titrate computes pKas on the [italic]minimized[/italic] structure. The "
              "Protonation State Analyzer's earlier ProPKA run (if any) was "
              "typically on the [italic]unminimized[/italic] input, so comparing to it "
              "would mix the method difference (PB vs ProPKA) with a coordinate "
              "difference. To keep the comparison honest, this step runs ProPKA "
              "on the exact minimized, solvent-stripped structure PB used.\n\n"
              "[bold]What it does[/bold]\n"
              "Writes a standard-named, hydrogen-free PDB of the minimized "
              "structure (ProPKA reconstructs protons from heavy atoms), runs "
              "propka3 on it, and lines the two predictions up per residue:\n"
              "  pb_titrate_work/pka_vs_propka.csv — resname, resnum, PB pKa, "
              "ProPKA pKa, Δ(PB−ProPKA).\n\n"
              "[bold]Notes[/bold]\n"
              "  • Matching is by residue type + number; chains are ignored.\n"
              "  • ProPKA omits non-titratable groups (e.g. disulfide-bonded "
              "CYS), which then show no ProPKA value.\n"
              "  • Metal-coordinating residues excluded in step 2 won't appear "
              "(they aren't PB sites); their large method gap is the metal "
              "effect PB captures and ProPKA doesn't.\n"
              "  • Requires propka3 on PATH (ships with the Protonation State "
              "Analyzer's dependencies); if absent, the step reports the PB "
              "pKas alone."),
    ),
]


# ============================================================================
# Workflow executor
# ============================================================================

class PBTitrateWorkflow:
    """Handler-bearing executor for the WorkflowChecklist.

    Each step is a `_checklist_pbt_N_xxx` method returning a result dict
    with at least `{"summary": str}`. Raising any exception marks the
    step as failed in the checklist state.
    """

    def __init__(self, processor):
        self.processor = processor
        self.console = processor.console

    # -- workspace helpers --------------------------------------------------

    def _ws(self):
        return self.processor.workspace

    def _params(self) -> Dict[str, Any]:
        return self._ws().get("pb_titrate_params") or {}

    def _output_dir(self) -> Path:
        out = Path(self._ws().get("output_dir") or ".") / "pb_titrate_work"
        out.mkdir(parents=True, exist_ok=True)
        return out

    # Water residue names — always stripped before PB. PB treats solvent
    # implicitly via the εout/Debye-Hückel model; explicit waters in a cluster
    # are stuck at εin=4 with their tip3p charges and produce wildly wrong
    # reaction-field energies. Water names are unambiguous, so a name match is
    # safe here (unlike ions — see _classify_ions_for_strip).
    WATER_RESNAMES = frozenset({
        "WAT", "HOH", "T3P", "TIP3", "T4P", "TIP4", "T5P", "TIP5",
        "OPC", "OPC3", "SOL", "SPC", "SPCE",
    })

    def _redox_site_reskeys(self) -> set:
        """(chain, resid_1based, insertion_code) keys for every residue in any
        detected RedoxSite. These are user-declared structural components — a
        structural Ca²⁺/Zn²⁺ site, a heme envelope, etc. — and must never be
        stripped as 'bulk', even though by element/name they look like a
        free ion."""
        from .sites import _residue_keys_in_redox_site
        keys = set()
        for site in (self._ws().get("detected_redox_sites", []) or []):
            try:
                keys |= _residue_keys_in_redox_site(site)
            except Exception:
                continue
        return keys

    def _classify_ions_for_strip(self, structure, redox_keys, cutoff):
        """Classify each monatomic ion residue as 'keep' (structural) or
        'strip' (bulk).

        Bulk neutralizing ions and structural ions can share an element AND a
        residue name (Amber writes calcium as ``CA`` for both), so a name mask
        cannot separate them. Instead each single-atom ion is kept iff it is
        either (a) part of a detected RedoxSite — the user's explicit
        structural declaration — or (b) coordinated: within ``cutoff`` Å of a
        non-solvent, non-ion heavy atom (the backstop for structural ions that
        weren't declared, e.g. EF-hand Ca²⁺). Everything else is bulk.

        Returns a list of dicts {resname, resnum, chain, atom_idx, decision,
        reason}; empty if the structure has no monatomic ions.
        """
        import numpy as np
        from .sites import _parmed_reskey

        ion_residues = [r for r in structure.residues
                        if len(r.atoms) == 1
                        and (r.name or "").strip() not in self.WATER_RESNAMES]
        if not ion_residues:
            return []

        coords = np.asarray(structure.coordinates, dtype=float)
        ion_atom_idxs = {r.atoms[0].idx for r in ion_residues}
        coord_idxs = [a.idx for a in structure.atoms
                      if a.idx not in ion_atom_idxs
                      and (a.residue.name or "").strip() not in self.WATER_RESNAMES
                      and a.atomic_number > 1]
        coord_xyz = coords[coord_idxs] if coord_idxs else None
        coord_atoms = [structure.atoms[i] for i in coord_idxs]

        out = []
        for r in ion_residues:
            a = r.atoms[0]
            entry = {
                "resname": (r.name or "").strip(),
                "resnum": r.number + 1,
                "chain": str(getattr(r, "chain", "") or ""),
                "atom_idx": a.idx,
            }
            if _parmed_reskey(r) in redox_keys:
                entry["decision"], entry["reason"] = "keep", "RedoxSite"
            elif coord_xyz is not None and len(coord_xyz):
                d = np.linalg.norm(coord_xyz - coords[a.idx], axis=1)
                j = int(np.argmin(d))
                if d[j] <= cutoff:
                    ca = coord_atoms[j]
                    entry["decision"] = "keep"
                    entry["reason"] = (
                        f"coordinated: {(ca.residue.name or '').strip()}"
                        f"{ca.residue.number + 1}-{ca.name} @ {d[j]:.2f}Å")
                else:
                    entry["decision"], entry["reason"] = "strip", "bulk (not coordinated)"
            else:
                entry["decision"], entry["reason"] = "strip", "bulk (not coordinated)"
            out.append(entry)
        return out

    def _confirm_ion_classification(self, classified) -> None:
        """Show the keep/strip split and let the user toggle any ion by residue
        number before stripping. Mutates ``classified`` entries' 'decision' in
        place."""
        console = self.console
        keep = [e for e in classified if e["decision"] == "keep"]
        bulk = [e for e in classified if e["decision"] == "strip"]

        def _label(e):
            tag = f":{e['chain']}" if e["chain"] else ""
            return f"{e['resname']} {e['resnum']}{tag}  [grey50]{e['reason']}[/grey50]"

        console.print("\n[bold]Ion classification for PB strip[/bold]")
        console.print("  [green]KEEP (structural):[/green]"
                      + ("" if keep else " [grey50]none[/grey50]"))
        for e in keep:
            console.print(f"    {_label(e)}")
        console.print("  [yellow]STRIP (bulk):[/yellow]"
                      + ("" if bulk else " [grey50]none[/grey50]"))
        for e in bulk:
            console.print(f"    {_label(e)}")

        resp = prompt_with_context(
            processor=self.processor,
            prompt=("Residue numbers to reclassify (toggle KEEP↔STRIP), "
                    "space/comma-separated; Enter to accept"),
            default="",
            module="PB Titrate",
            description="Confirm/override ion strip classification",
        )
        toks = (resp or "").replace(",", " ").split()
        if not toks:
            return
        by_num = {}
        for e in classified:
            by_num.setdefault(e["resnum"], []).append(e)
        for t in toks:
            try:
                num = int(t)
            except ValueError:
                console.print(f"  [yellow]ignored '{t}' (not a residue number)[/yellow]")
                continue
            matches = by_num.get(num)
            if not matches:
                console.print(f"  [yellow]residue {num} is not a monatomic ion — ignored[/yellow]")
                continue
            for e in matches:
                e["decision"] = "strip" if e["decision"] == "keep" else "keep"
                tag = f":{e['chain']}" if e["chain"] else ""
                console.print(f"  [cyan]toggled {e['resname']} {num}{tag} "
                              f"→ {e['decision'].upper()}[/cyan]")

    def _make_dry_prmtop(self, prmtop: Path, rst7: Path
                           ) -> Tuple[Path, Path, int, int]:
        """Strip waters + bulk ions, write a 'dry' prmtop to output_dir.

        Waters are stripped by name. Monatomic ions are classified structural
        (kept) vs bulk (stripped) by RedoxSite membership + a coordination
        backstop — a residue-name mask can't separate a structural Ca²⁺ from a
        bulk Ca²⁺ counter-ion, so kept structural ions stay in the PB calc and
        in the protonation-renamed PDB while bulk ions are removed. The original
        solvated prmtop is preserved for step 7 (state application + ion
        rebalance), which still needs the bulk ions.

        Returns (dry_prmtop, dry_rst7, n_residues_removed, n_atoms_removed). If
        nothing needed stripping, returns the original paths and zero counts.
        """
        import parmed as pmd
        from .metal_coordination import DEFAULT_COORDINATION_CUTOFF
        s = pmd.load_file(str(prmtop), str(rst7))
        n_atoms_before = len(s.atoms)
        n_res_before = len(s.residues)

        waters = sorted({(r.name or "").strip() for r in s.residues
                         if (r.name or "").strip() in self.WATER_RESNAMES})
        classified = self._classify_ions_for_strip(
            s, self._redox_site_reskeys(), DEFAULT_COORDINATION_CUTOFF)

        if not waters and not classified:
            return prmtop, rst7, 0, 0

        if classified:
            self._confirm_ion_classification(classified)
        bulk = [e for e in classified if e["decision"] == "strip"]

        if not waters and not bulk:
            return prmtop, rst7, 0, 0

        # Strip bulk ions first by ORIGINAL atom index (a 1-based @-mask),
        # before the water strip shifts indices; then waters by name.
        if bulk:
            s.strip("@" + ",".join(str(e["atom_idx"] + 1) for e in bulk))
        if waters:
            s.strip(":" + ",".join(waters))

        n_removed = n_atoms_before - len(s.atoms)
        n_res_removed = n_res_before - len(s.residues)
        out_dir = self._output_dir() / "dry"
        out_dir.mkdir(parents=True, exist_ok=True)
        dry_prm = out_dir / (prmtop.stem + "_dry.prmtop")
        dry_rst = out_dir / (prmtop.stem + "_dry.rst7")
        s.save(str(dry_prm), overwrite=True)
        s.save(str(dry_rst), overwrite=True)
        return dry_prm, dry_rst, n_res_removed, n_removed

    def _resolve_prmtop_rst7(self) -> Tuple[Path, Path]:
        """Find a prmtop+rst7 pair in the workspace's output_dir.

        For single-state runs, expects exactly one *.prmtop in output_dir
        (or a sibling subdirectory). For batch-microstate runs, prompts the
        user to pick one — the workflow runs per-microstate.
        """
        ws = self._ws()
        # Canonical workspace keys: the Topology Generator writes
        # parm7_file/rst7_file and PB Titrate's own apply step (pbt-7)
        # promotes to the same pair; the older prmtop/prmtop_file aliases are
        # kept for back-compat. Reading parm7_file here is what stops this
        # from silently falling through to the glob below when PB Titrate is
        # launched straight after topology generation.
        prmtop = (ws.get("prmtop") or ws.get("parm7_file")
                  or ws.get("prmtop_file"))
        rst7   = ws.get("rst7") or ws.get("rst7_file")
        if prmtop and rst7 and Path(prmtop).exists() and Path(rst7).exists():
            return Path(prmtop), Path(rst7)

        # Glob fallback
        out = ws.get("output_dir", ".")
        prmtops = sorted(glob.glob(os.path.join(out, "**/*.prmtop"),
                                    recursive=True))
        if not prmtops:
            raise FileNotFoundError(
                f"No *.prmtop found under {out}. Run Topology Generator first.")

        if len(prmtops) == 1:
            prmtop = Path(prmtops[0])
        else:
            self.console.print(f"\n  Found {len(prmtops)} prmtop files:")
            options_map = {}
            for i, p in enumerate(prmtops):
                self.console.print(f"    {i}: {p}")
                options_map[str(i)] = p
            choice = int_prompt_with_context(
                self.processor,
                "Which prmtop to titrate?",
                default=0,
                module=MODULE_NAME,
                description="Pick prmtop for PB titration",
                options_map=options_map,
            )
            prmtop = Path(prmtops[choice])

        rst7 = prmtop.with_suffix(".rst7")
        if not rst7.exists():
            raise FileNotFoundError(
                f"Found {prmtop} but no matching {rst7}")
        return prmtop, rst7

    def _minimized_rst7_for(self, prmtop: Path) -> Optional[Path]:
        """Return the minimized rst7 for ``prmtop`` if pbt-m produced one.

        ``pb_titrate_minimized`` maps an absolute prmtop path to its minimized
        rst7. Keying by prmtop (rather than a single scalar) keeps this
        multi-state-ready: a future per-microstate pbt-m just adds one entry
        per microstate topology.
        """
        mapping = self._ws().get("pb_titrate_minimized") or {}
        hit = mapping.get(str(Path(prmtop).resolve()))
        if hit and Path(hit).exists():
            return Path(hit)
        return None

    def _make_pbsa_kwargs(self) -> Dict[str, Any]:
        p = self._params()
        return {
            "epsin":  p["epsin"],
            "istrng": p["istrng"],
            "space":  p["space"],
            "nfocus": p["nfocus"],
            "bcopt":  p["bcopt"],
        }

    # -- step 1: configure --------------------------------------------------

    def _checklist_pbt_1_configure(self) -> Dict[str, Any]:
        ws = self._ws()
        cur = self._params()

        self.console.print("\n[bold cyan]Configure PB calculation parameters[/bold cyan]")
        self.console.print(
            "[grey50]Each prompt has a short explanation of what it controls and a "
            "suggested default in brackets. Press Enter to accept the default.[/grey50]\n")

        # --- Interior dielectric ---
        self.console.print(
            "[bold]Interior dielectric (εin)[/bold] — models polarizability of the\n"
            "protein interior. Lower values treat the protein as rigid; higher\n"
            "values approximate dynamic averaging of conformations. The right\n"
            "value depends on your system (rigidity of the relevant region,\n"
            "burial depth of the titratable groups). The current default is\n"
            "shown in brackets.")
        epsin  = float_prompt_with_context(
            self.processor, "  εin",
            default=cur.get("epsin", 4.0),
            module=MODULE_NAME, description="Interior dielectric")

        # --- Ionic strength ---
        self.console.print(
            "\n[bold]Ionic strength (mM)[/bold] — monovalent salt concentration;\n"
            "controls Debye screening of long-range electrostatics. Set this to\n"
            "match the salt condition of the simulation you intend to run.\n"
            "0 = vacuum / no-salt limit.")
        istrng_default = cur.get("istrng", ws.get("ionic_strength_mM", 150.0))
        istrng = float_prompt_with_context(
            self.processor, "  ionic strength (mM)",
            default=istrng_default,
            module=MODULE_NAME, description="Ionic strength in mM")

        # --- Grid spacing ---
        self.console.print(
            "\n[bold]PB grid spacing (Å)[/bold] — finite-difference resolution.\n"
            "Smaller spacing = more accurate but slower (PB cost grows as\n"
            "~1/spacing^3). Convergence depends on the system; do a sweep\n"
            "of a few values on one site if accuracy matters for your case.")
        space  = float_prompt_with_context(
            self.processor, "  grid spacing (Å)",
            default=cur.get("space", 0.5),
            module=MODULE_NAME, description="PB grid spacing")

        # --- Focusing levels ---
        self.console.print(
            "\n[bold]PB focusing levels[/bold] — number of nested grids; each\n"
            "step focuses tighter around the target. More levels improves the\n"
            "potential near the target at marginal extra cost; the right number\n"
            "depends on system size and how deeply buried the target sits.")
        nfocus = int_prompt_with_context(
            self.processor, "  focusing levels",
            default=cur.get("nfocus", 2),
            module=MODULE_NAME, description="PB focusing levels")

        # --- Boundary condition ---
        self.console.print(
            "\n[bold]Boundary condition (bcopt)[/bold] — how pbsa sets the\n"
            "potential at the grid edge:\n"
            "  5 = Debye-Hückel from grid charges  (physical, stable)\n"
            "  6 = DH with charge-singularity removal  (more accurate, slower)\n"
            "  1 = Zero-potential / conductor  (cheap; testing only)")
        bcopt = int(prompt_with_context(
            self.processor, "  bcopt",
            choices=["5", "6", "1"],
            default=str(cur.get("bcopt", 5)),
            module=MODULE_NAME, description="PB boundary condition (pbsa bcopt)",
            options_map={
                "5": "Debye-Hückel from grid charges",
                "6": "DH with charge-singularity removal",
                "1": "Zero-potential (conductor)",
            }))

        # --- Workers ---
        self.console.print(
            "\n[bold]Parallel pbsa workers[/bold] — concurrent pbsa processes.\n"
            "Each worker handles one PB call (one site × one state). The right\n"
            "number depends on your hardware: enough to keep cores busy without\n"
            "thermal throttling on your specific machine.")
        workers = int_prompt_with_context(
            self.processor, "  workers",
            default=cur.get("workers", 4),
            module=MODULE_NAME, description="Parallel pbsa workers")

        # --- Target pH ---
        self.console.print(
            "\n[bold]Target simulation pH[/bold] — pH at which to predict\n"
            "protonation states. Match this to your intended simulation pH.")
        pH = float_prompt_with_context(
            self.processor, "  pH",
            default=cur.get("pH", ws.get("simulation_pH", 7.0)),
            module=MODULE_NAME, description="Target pH for state recommendations")

        # --- Cluster cutoff (Y/N then optional radius) ---
        self.console.print(
            "\n[bold]Cluster cutoff[/bold] — restrict each PB call to a sub-cluster\n"
            "around the target site instead of using the full structure. Trades\n"
            "accuracy for speed; for very large systems the full-structure path\n"
            "may be impractical. Off = use full structure.\n"
            "\n"
            "How the cutoff is applied (it is NOT a single sphere):\n"
            "  1. Take every heavy atom in the site's envelope (the residue\n"
            "     itself, or every residue of the detected redox site it\n"
            "     belongs to).\n"
            "  2. Form the union of R-Å balls around each of those seed\n"
            "     heavy atoms.\n"
            "  3. Any residue with at least one heavy atom inside that union\n"
            "     is kept in full (all of its atoms, including those farther\n"
            "     than R).\n"
            "  4. Plus a covalent bond walk seeded from the same atoms.\n"
            "  5. Plus the all-or-nothing rule: any other detected redox site\n"
            "     that contributes any residue to the cluster gets fully\n"
            "     included (no partial inclusion).\n"
            "Effective cluster extent ≈ R + the envelope's own diameter,\n"
            "often a bit more after residue promotion and bond walk.")
        prev_cluster_radius = cur.get("cluster_radius")
        use_cutoff_default = (prev_cluster_radius is not None
                              if "cluster_radius" in cur else True)
        use_cutoff = confirm_with_context(
            self.processor,
            "  Use a cluster cutoff (instead of the full structure)?",
            default=use_cutoff_default,
            module=MODULE_NAME,
            description="Use cluster cutoff sphere instead of full structure")
        if use_cutoff:
            radius_default = float(prev_cluster_radius) if prev_cluster_radius else 25.0
            cluster_radius = float_prompt_with_context(
                self.processor, "  cutoff radius (Å)",
                default=radius_default,
                module=MODULE_NAME, description="Cluster cutoff radius (Å)")
        else:
            cluster_radius = None

        # --- Targeted coupling (subset restriction) ---
        self.console.print(
            "\n[bold]Targeted coupling[/bold] — when mean-field (step 5) does not\n"
            "converge, a small frustrated subset of titratable sites keeps\n"
            "flipping. Computing the full pairwise coupling matrix (step 4) over\n"
            "every site to resolve them is wasteful: the cost grows as the square\n"
            "of the site count. Targeted mode restricts the coupling + Monte-Carlo\n"
            "work to the unconverged sites plus every titratable site whose nearest\n"
            "heavy atom lies within a chosen radius of them. All other sites are\n"
            "held fixed at their converged mean-field state. Leave off for the\n"
            "normal full run.")
        coupling_targeted = confirm_with_context(
            self.processor,
            "  Restrict coupling to unconverged sites + near neighbors?",
            default=bool(cur.get("coupling_targeted", False)),
            module=MODULE_NAME,
            description="Restrict coupling to unconverged sites + near neighbors")
        if coupling_targeted:
            self.console.print(
                "\n  [bold]Subgraph radius (Å)[/bold] — a titratable site joins the\n"
                "  coupling subgraph if its nearest heavy atom comes within this\n"
                "  distance of the nearest heavy atom of any unconverged site.\n"
                "  This same value also becomes the pair-distance cutoff inside the\n"
                "  subgraph, replacing the ionic-strength auto value (which can be\n"
                "  tens of Å at low salt).")
            coupling_subgraph_radius = float_prompt_with_context(
                self.processor, "  subgraph radius (Å)",
                default=float(cur.get("coupling_subgraph_radius", 12.0)),
                module=MODULE_NAME,
                description="Targeted-coupling subgraph radius (Å)")
            coupling_distance_cutoff = coupling_subgraph_radius
        else:
            coupling_subgraph_radius = cur.get("coupling_subgraph_radius", 12.0)
            coupling_distance_cutoff = cur.get("coupling_distance_cutoff")

        params = {
            "epsin": epsin, "istrng": istrng, "space": space,
            "nfocus": nfocus, "bcopt": bcopt,
            "workers": workers, "pH": pH,
            "cluster_radius": cluster_radius,
            "coupling_targeted": coupling_targeted,
            "coupling_subgraph_radius": coupling_subgraph_radius,
            "coupling_distance_cutoff": coupling_distance_cutoff,
        }
        ws.set("pb_titrate_params", params)

        self.console.print(f"\n[green]✓ PB parameters saved[/green]")
        radius_summary = f"R={cluster_radius:.1f} Å" if cluster_radius else "R=full structure"
        targeted_summary = (f", targeted coupling @ {coupling_subgraph_radius:.0f} Å"
                            if coupling_targeted else "")
        return {"summary": (f"εin={epsin}, I={istrng:.0f} mM, pH={pH}, "
                              f"{radius_summary}, {workers} workers"
                              f"{targeted_summary}")}

    # -- step m: minimize structure -----------------------------------------

    def _checklist_pbt_m_minimize(self) -> Dict[str, Any]:
        from .minimize import (
            MinimizationParams, SUPPORTED_ENGINES, run_minimization,
            DEFAULT_MAXCYC, DEFAULT_NCYC, DEFAULT_CUT, DEFAULT_RESTRAINT_WT,
            DEFAULT_RESTRAINTMASK, DEFAULT_ENGINE, DEFAULT_MPI_TASKS,
        )

        ws = self._ws()
        # Minimize the base (un-minimized) structure — never a prior pbt-m
        # output. _resolve_prmtop_rst7 returns the base pair; the minimized
        # substitution lives in pbt-2, not the resolver.
        prmtop, rst7 = self._resolve_prmtop_rst7()

        self.console.print("\n[bold cyan]Minimize structure before PB[/bold cyan]")
        self.console.print(f"  prmtop: {prmtop}")
        self.console.print(f"  rst7:   {rst7}")

        # Reuse a prior minimization if one already exists for this prmtop
        # (step run earlier, or being re-run): a minimization can be long, and
        # the result is deterministic for fixed inputs. Look in the workspace
        # map first, then fall back to the on-disk min.rst7 (covers a fresh
        # session where the map wasn't persisted).
        prmtop_key = str(Path(prmtop).resolve())
        mapping = dict(ws.get("pb_titrate_minimized") or {})
        work_dir = self._output_dir() / "minimization" / prmtop.stem
        prior_rst7 = None
        mapped = mapping.get(prmtop_key)
        if mapped and Path(mapped).exists():
            prior_rst7 = Path(mapped)
        elif (work_dir / "min.rst7").exists():
            prior_rst7 = work_dir / "min.rst7"
        if prior_rst7 is not None:
            self.console.print(
                f"\n  [green]Found an existing minimized structure:[/green] "
                f"{prior_rst7}")
            mdout = prior_rst7.parent / "min.out"
            if mdout.exists():
                from .minimize import final_energy
                e = final_energy(mdout)
                if e is not None:
                    self.console.print(
                        f"  [grey50]Final energy: {e:.2f} kcal/mol[/grey50]")
            if confirm_with_context(
                    self.processor,
                    "  Reuse it and skip re-running minimization?",
                    default=True, module=MODULE_NAME,
                    description="Reuse existing minimized structure"):
                mapping[prmtop_key] = str(Path(prior_rst7).resolve())
                ws.set("pb_titrate_minimized", mapping)
                self.console.print(
                    "  [green]Using existing minimized coordinates; "
                    "minimization not re-run.[/green]")
                return {"summary": ("reused existing minimized structure "
                                    f"({prior_rst7.name})")}
            self.console.print("  [grey50]Re-running minimization...[/grey50]")

        self.console.print(
            "\n[grey50]PB pKa shifts are sensitive to steric clashes; a short "
            "minimization relieves bad contacts before any PB call. Defaults "
            "mirror the bundled aqueous-metalloprotein minimization template; "
            "press Enter to accept each.[/grey50]\n")

        engine_options = {
            "sander": "Serial CPU (always available; slowest)",
            "pmemd": "Optimized CPU (recommended)",
            "pmemd.MPI": "Multi-core CPU via mpirun (parallel)",
        }
        self.console.print("  [bold]Engine options:[/bold]")
        for name in SUPPORTED_ENGINES:
            marker = " [grey50](default)[/grey50]" if name == DEFAULT_ENGINE else ""
            self.console.print(
                f"    [cyan]{name}[/cyan] — {engine_options[name]}{marker}")
        engine = prompt_with_context(
            self.processor, "  Engine",
            choices=list(SUPPORTED_ENGINES),
            default=DEFAULT_ENGINE,
            module=MODULE_NAME, description="MD engine for minimization",
            options_map=engine_options)
        mpi_tasks = DEFAULT_MPI_TASKS
        if engine == "pmemd.MPI":
            mpi_tasks = int_prompt_with_context(
                self.processor, "  MPI tasks (-np)",
                default=DEFAULT_MPI_TASKS,
                module=MODULE_NAME, description="MPI task count for pmemd.MPI")

        maxcyc = int_prompt_with_context(
            self.processor, "  Max minimization cycles (maxcyc)",
            default=DEFAULT_MAXCYC,
            module=MODULE_NAME, description="Maximum minimization cycles")
        ncyc = int_prompt_with_context(
            self.processor, "  Steepest-descent → CG switch (ncyc)",
            default=DEFAULT_NCYC,
            module=MODULE_NAME, description="SD-to-CG switch cycle")
        cut = float_prompt_with_context(
            self.processor, "  Nonbonded cutoff Å (cut)",
            default=DEFAULT_CUT,
            module=MODULE_NAME, description="Nonbonded cutoff in Angstroms")

        ntr = confirm_with_context(
            self.processor,
            "  Restrain atoms during minimization? (relieves clashes "
            "without large conformational drift)",
            default=True,
            module=MODULE_NAME, description="Enable positional restraints")
        restraintmask = DEFAULT_RESTRAINTMASK
        restraint_wt = DEFAULT_RESTRAINT_WT
        if ntr:
            restraintmask = prompt_with_context(
                self.processor, "  Restraint mask (AMBER syntax)",
                default=DEFAULT_RESTRAINTMASK,
                module=MODULE_NAME, description="Atoms to restrain")
            restraint_wt = float_prompt_with_context(
                self.processor, "  Restraint force constant (kcal/mol·Å²)",
                default=DEFAULT_RESTRAINT_WT,
                module=MODULE_NAME, description="Positional restraint weight")

        params = MinimizationParams(
            maxcyc=maxcyc, ncyc=ncyc, cut=cut, ntr=ntr,
            restraint_wt=restraint_wt, restraintmask=restraintmask,
            engine=engine, mpi_tasks=mpi_tasks)

        self.console.print()
        minimized_rst7 = run_minimization(
            prmtop, rst7, work_dir, params, console=self.console)

        # Record in the prmtop-keyed map so pbt-2 picks up the minimized
        # coordinates (multi-state-ready: one entry per topology).
        mapping = dict(ws.get("pb_titrate_minimized") or {})
        mapping[str(Path(prmtop).resolve())] = str(Path(minimized_rst7).resolve())
        ws.set("pb_titrate_minimized", mapping)

        restr = (f"restrained ({restraintmask} @ {restraint_wt})"
                 if ntr else "unrestrained")
        return {"summary": (f"{engine}, {maxcyc} cyc, {restr} → "
                            f"{minimized_rst7.name}")}

    # -- step 2: discover sites ---------------------------------------------

    def _checklist_pbt_2_discover(self) -> Dict[str, Any]:
        import parmed as pmd
        from .sites import discover_sites

        ws = self._ws()
        original_prmtop, original_rst7 = self._resolve_prmtop_rst7()

        # If pbt-m minimized this structure, titrate the minimized
        # coordinates. The prmtop is unchanged by minimization; only the rst7
        # moves. This flows through to step 7 too, because the solvated rst7
        # cached below is what apply_state reads.
        minimized = self._minimized_rst7_for(original_prmtop)
        if minimized is not None:
            self.console.print(
                "  [green]Using minimized coordinates from the Minimize "
                "Structure step (pbt-m).[/green]")
            original_rst7 = minimized

        self.console.print(f"\n[bold cyan]Discovering titratable sites[/bold cyan]")
        self.console.print(f"  prmtop: {original_prmtop}")
        self.console.print(f"  rst7:   {original_rst7}")

        # If the prmtop has explicit waters or bulk ions, build a stripped
        # 'dry' version for PB calculations. PB treats solvent and ions
        # implicitly through εout and Debye screening — explicit waters at
        # εin=4 in a PB cluster compete with that description and produce
        # wildly wrong reaction-field energies (saw a -9 pKa where the
        # correct value was ~3 because TIP3P waters were inside the
        # cluster). The original prmtop is preserved and used in step 7
        # for state application + ion rebalance.
        dry_prmtop, dry_rst7, n_res_stripped, n_atoms_stripped = (
            self._make_dry_prmtop(original_prmtop, original_rst7))
        if n_res_stripped > 0:
            self.console.print(
                f"  [cyan]Stripped {n_res_stripped} solvent/ion residue(s) "
                f"({n_atoms_stripped} atoms) for PB calcs.[/cyan]")
            self.console.print(f"  dry prmtop: {dry_prmtop}")
            self.console.print(
                f"  [grey50](original solvated prmtop kept for step 7's state "
                f"application + ion rebalance.)[/grey50]")

        structure = pmd.load_file(str(dry_prmtop), str(dry_rst7))
        redox_sites = ws.get("detected_redox_sites") or []
        sites, terminal_excluded = discover_sites(structure, redox_sites)

        # Surface excluded terminal residues clearly — Amber's constph has
        # no multi-state libraries for terminal AAs, so titrating them
        # would silently corrupt charges (see CLYS/CASP gap in
        # aminoct10.lib / aminont10.lib).
        if terminal_excluded:
            self.console.print(
                f"\n  [yellow]Excluded {len(terminal_excluded)} terminal "
                f"residue(s) — constph has no multi-state libraries for "
                f"chain termini, so these are held at their prmtop "
                f"protonation state:[/yellow]")
            for site, kind in terminal_excluded:
                marker_atoms = sorted(
                    {a.name for a in site.titrating_residue.atoms}
                    & {"H1", "H2", "H3", "OXT"})
                self.console.print(
                    f"    [yellow]{site.resname}-{site.resnum}[/yellow] "
                    f"({kind}; detected via "
                    f"{', '.join(marker_atoms) or '?'})")

        if not sites:
            raise RuntimeError(
                "No titratable internal residues "
                "(AS4/GL4/HIP/LYS/CYS/TYR/PRN) found in prmtop. Was the "
                "topology built with leaprc.constph? "
                f"({len(terminal_excluded)} terminal residue(s) were "
                f"excluded — constph cannot titrate those.)")

        counts = Counter(s.resname for s in sites)
        self.console.print(f"\n  Found {len(sites)} titratable sites:")
        for name, n in sorted(counts.items()):
            self.console.print(f"    {name}: {n}")

        multi = [s for s in sites if s.is_multi_residue]
        if multi:
            self.console.print(f"\n  Multi-residue sites (envelope > 1):")
            for s in multi[:10]:
                self.console.print(
                    f"    {s.resname}-{s.resnum}: envelope = "
                    f"{len(s.envelope_residues)} residues "
                    f"(redox site: {s.redox_site_id})")
            if len(multi) > 10:
                self.console.print(f"    ... and {len(multi)-10} more")

        # Optional exclusion
        exclude_str = prompt_with_context(
            self.processor,
            "Exclude any sites? (RESNAME:RESNUM,RESNAME:RESNUM,... or empty)",
            default="",
            module=MODULE_NAME, description="Sites to exclude from titration")
        excluded_count = 0
        if exclude_str.strip():
            excludes = set()
            for tok in exclude_str.split(","):
                tok = tok.strip()
                if not tok: continue
                if ":" not in tok:
                    self.console.print(f"  [yellow]skipping malformed token "
                                          f"{tok!r}[/yellow]")
                    continue
                rn, rn_num = tok.split(":", 1)
                excludes.add((rn.strip(), int(rn_num.strip())))
            kept_sites = [s for s in sites if (s.resname, s.resnum) not in excludes]
            excluded_count = len(sites) - len(kept_sites)
            sites = kept_sites
            self.console.print(f"  After exclusion: {len(sites)} sites "
                                  f"({excluded_count} excluded)")

        # Metal-coordination handling. A titrating side chain in inner-sphere
        # contact with a metal cation is metal-locked (not pH-titratable in the
        # holo state) and breaks fixed-background PB. Detect geometrically,
        # drop them from the titratable set, and record the charge state to
        # hold them at in the background (consumed by pbt-3). Restoring their
        # real charge also de-biases the environment every other site sees.
        from .metal_coordination import (
            find_metal_coordinated_sites, classify_fix,
            DEFAULT_COORDINATION_CUTOFF)
        metal_fixed: List[Dict[str, Any]] = []
        coordinated = find_metal_coordinated_sites(
            structure, sites, cutoff=DEFAULT_COORDINATION_CUTOFF)
        if coordinated:
            self.console.print(
                f"\n  [bold yellow]{len(coordinated)} titratable site(s) "
                f"directly coordinate a metal ion[/bold yellow] "
                f"(within {DEFAULT_COORDINATION_CUTOFF} Å):")
            for rec in coordinated:
                m = rec["metals"][0]
                extra = (f" (+{rec['n_metals']-1} more metal)"
                         if rec["n_metals"] > 1 else "")
                self.console.print(
                    f"    {rec['resname']}-{rec['resnum']}  "
                    f"{rec['coord_atom']}→{m['element']}"
                    f"{m['resname']}{m['resnum']} {m['distance']} Å{extra}")
            self.console.print(
                "  [grey50]These aren't pH-titratable while metal-bound and "
                "corrupt PB pKas. Recommended: exclude them and fix them at "
                "their metal-bound charge.[/grey50]")
            if confirm_with_context(
                    self.processor,
                    "  Exclude these metal-coordinating sites and fix their "
                    "background charge?",
                    default=True, module=MODULE_NAME,
                    description="Exclude metal-coordinating sites"):
                drop = set()
                for rec in coordinated:
                    key = (rec["resname"], rec["resnum"])
                    cls = classify_fix(rec["resname"], rec["n_metals"])
                    fix_mode = cls["fix_mode"]
                    if cls["needs_prompt"]:
                        # Let the user confirm the chemistry for the ambiguous
                        # cases (His protonation / bridging histidinate / rare
                        # coordinators). Options depend on what is representable.
                        opts = []
                        if cls["fix_mode"] == "neutral":
                            opts = [("neutral", "neutral imidazole (HID/HIE)"),
                                    ("deprotonated", "deprotonated (-1) if a state exists"),
                                    ("keep", "keep prmtop charges (exclude only)")]
                            default_choice = "neutral"
                        else:
                            opts = [("neutral", "neutral form"),
                                    ("keep", "keep prmtop charges (exclude only)")]
                            default_choice = "keep"
                        self.console.print(
                            f"\n    {rec['resname']}-{rec['resnum']}: "
                            f"{cls['note']}")
                        choice = prompt_with_context(
                            self.processor,
                            f"      Fix {rec['resname']}-{rec['resnum']} as",
                            choices=[o[0] for o in opts],
                            default=default_choice,
                            module=MODULE_NAME,
                            description="Metal coordinator fix state",
                            options_map={o[0]: o[1] for o in opts})
                        fix_mode = None if choice == "keep" else choice
                    drop.add(key)
                    if fix_mode is not None:
                        metal_fixed.append({
                            "resname": rec["resname"], "resnum": rec["resnum"],
                            "fix_mode": fix_mode,
                            "metal": (f"{rec['metals'][0]['element']}"
                                      f"{rec['metals'][0]['resnum']}"),
                            "distance": rec["min_distance"],
                        })
                before = len(sites)
                sites = [s for s in sites if (s.resname, s.resnum) not in drop]
                self.console.print(
                    f"  [green]Excluded {before - len(sites)} metal-coordinating "
                    f"site(s); {len(metal_fixed)} fixed at metal-bound charge "
                    f"in the background.[/green]")
        ws.set("pb_titrate_metal_fixed", metal_fixed)

        # Persist (serializable form). Indices stable across same-prmtop reloads.
        serialized = [
            {"resname": s.resname, "resnum": s.resnum,
             "envelope_idxs": sorted(s.envelope_idxs),
             "redox_site_id": s.redox_site_id,
             "other_redox_envelopes": [sorted(env)
                                          for env in s.other_redox_envelopes]}
            for s in sites
        ]
        ws.set("pb_titrate_sites", serialized)
        # `pb_titrate_prmtop` / `pb_titrate_rst7` are the PB-ready (dry)
        # files. `pb_titrate_prmtop_solvated` / `pb_titrate_rst7_solvated`
        # are the originals — used by step 7 for state apply + rebalance.
        # Stored as Path (not str) so the MD manager's str-only auto-
        # discovery scanner doesn't offer these intermediates as candidate
        # production prmtops.
        ws.set("pb_titrate_prmtop", Path(dry_prmtop).resolve())
        ws.set("pb_titrate_rst7", Path(dry_rst7).resolve())
        ws.set("pb_titrate_prmtop_solvated", Path(original_prmtop).resolve())
        ws.set("pb_titrate_rst7_solvated", Path(original_rst7).resolve())
        ws.set("pb_titrate_terminal_excluded",
                [{"resname": s.resname, "resnum": s.resnum, "kind": kind}
                  for s, kind in terminal_excluded])

        term_str = (f", {len(terminal_excluded)} terminal-excluded"
                     if terminal_excluded else "")
        return {"summary": (f"{len(sites)} sites ({len(multi)} multi-residue, "
                              f"{excluded_count} user-excluded{term_str})")}

    # -- step 3: single-site pKa --------------------------------------------

    def _site_objects(self) -> List["Site"]:
        """Reconstitute Site objects from the serialized workspace form,
        bound to a freshly loaded structure."""
        import parmed as pmd
        from .sites import Site

        ws = self._ws()
        s = pmd.load_file(str(ws.get("pb_titrate_prmtop")),
                            str(ws.get("pb_titrate_rst7")))
        out = []
        for entry in ws.get("pb_titrate_sites") or []:
            envelope = [s.residues[i] for i in entry["envelope_idxs"]]
            target = next(r for r in envelope
                          if r.name == entry["resname"]
                          and r.number + 1 == entry["resnum"])
            out.append(Site(
                titrating_residue=target,
                envelope_residues=envelope,
                redox_site_id=entry["redox_site_id"],
                other_redox_envelopes=entry["other_redox_envelopes"],
            ))
        return out

    def _checklist_pbt_3_single_site(self) -> Dict[str, Any]:
        from .intrinsic import compute_pka_multistate
        from .residues import get_residue

        ws = self._ws()
        params = self._params()
        sites = self._site_objects()
        prmtop = Path(ws.get("pb_titrate_prmtop"))
        rst7   = Path(ws.get("pb_titrate_rst7"))

        # Build the all-neutral background (MEAD/Karlsberg convention):
        # every non-target titratable is held in its charge-closest-to-zero
        # form so the cluster's reference charge distribution is physical.
        # Without this, the prmtop-default STATE 0 charges leave AS4/GL4/PRN
        # at -1 and LYS/HIP at +1 simultaneously, which gives huge spurious
        # pKa shifts (especially for buried propionates near other
        # carboxylates).
        background = {
            (s.resname, s.resnum): get_residue(s.resname).neutral_state
            for s in sites
        }
        # Metal-coordinating residues (excluded from titration in pbt-2) are
        # held in the background at their metal-bound charge — not neutral —
        # so the cluster's reference charge distribution is physical and the
        # +metal field every other site feels is correct.
        metal_fixed = ws.get("pb_titrate_metal_fixed") or []
        n_metal_fixed = 0
        for entry in metal_fixed:
            res = get_residue(entry["resname"])
            mode = entry["fix_mode"]
            state = (res.deprot_state if mode == "deprotonated"
                     else res.neutral_state)
            background[(entry["resname"], entry["resnum"])] = state
            n_metal_fixed += 1
        if n_metal_fixed:
            self.console.print(
                f"  [grey50]{n_metal_fixed} metal-coordinating residue(s) held "
                f"at their metal-bound charge in the background.[/grey50]")

        self.console.print(f"\n[bold cyan]Single-site intrinsic pKa per site[/bold cyan]")
        self.console.print(f"  {len(sites)} sites at pH {params['pH']}; "
                              f"cluster_radius={params.get('cluster_radius')}")
        self.console.print(f"  background: all-neutral "
                              f"(non-target titratables held at their "
                              f"closest-to-neutral state)")
        if not confirm_with_context(
                self.processor,
                f"Run {len(sites)} single-site PB calculations now?",
                default=True,
                module=MODULE_NAME,
                description="Confirm single-site pKa run"):
            raise RuntimeError("User cancelled single-site PB run")

        pbsa_kw = self._make_pbsa_kwargs()
        wd_root = self._output_dir() / "single_site"
        # Model compounds are per-resname and shared across all sites. Cache
        # them under the run's work dir (writable; never the installed package,
        # which may be read-only). build_model prefers any prebuilt compound
        # shipped with the package and only builds the rest here.
        model_cache_dir = self._output_dir() / "model_compounds"
        results: Dict[Tuple[str, int], Dict[str, Any]] = {}
        failed: List[Tuple[str, int, str]] = []

        # Resume support. Each completed site is saved to the aggregate pickle
        # immediately (not just at the end), so an interrupted run can pick up
        # where it left off. Sites already present are skipped; previously
        # failed sites are retried (e.g. after a pbsa fix). At the PB-leg level
        # run_pbsa also reuses any parseable {label}.out on disk, so even a site
        # that wasn't recorded but whose outputs exist re-reads instead of
        # recomputing.
        pkl = self._output_dir() / "single_site_pka.pkl"
        recompute_all = False
        if pkl.exists():
            try:
                with pkl.open("rb") as f:
                    prior = pickle.load(f)
                results = dict(prior.get("results", {}))
            except Exception:
                results = {}
            if results:
                self.console.print(
                    f"\n  [green]Found {len(results)} previously-computed "
                    f"site(s) in {pkl.name}.[/green]")
                recompute_all = confirm_with_context(
                    self.processor,
                    "  Recompute ALL sites from scratch? "
                    "(no = reuse completed, run only missing/failed)",
                    default=False, module=MODULE_NAME,
                    description="Recompute all single-site pKas vs resume")
                if recompute_all:
                    results = {}

        def _save() -> None:
            with pkl.open("wb") as f:
                pickle.dump({"results": results, "failed": failed}, f)

        for i, site in enumerate(sites, 1):
            key = (site.resname, site.resnum)
            wd = wd_root / f"{site.resname}{site.resnum:04d}"
            if key in results and not recompute_all:
                pka = results[key].get("pKa_eff")
                self.console.print(
                    f"  [{i}/{len(sites)}] {site.resname}-{site.resnum}... "
                    f"[grey50]cached[/grey50] pKa="
                    f"{pka if pka is None else f'{pka:.2f}'}")
                continue
            self.console.print(f"  [{i}/{len(sites)}] {site.resname}-{site.resnum}...",
                                  end="")
            try:
                r = compute_pka_multistate(
                    cluster_prmtop=prmtop, cluster_rst7=rst7,
                    target_resname=site.resname,
                    target_resnum_1based=site.resnum,
                    pH=params["pH"], work_dir=wd,
                    model_cache_dir=model_cache_dir,
                    cluster_radius=params.get("cluster_radius"),
                    background=background,
                    site=site,
                    **pbsa_kw)
                results[key] = r
                # drop any stale failure record for this site
                failed = [t for t in failed if (t[0], t[1]) != key]
                pka = r.get("pKa_eff")
                self.console.print(f"  pKa={pka if pka is None else f'{pka:.2f}'}")
            except Exception as exc:
                self.console.print(f" [red]FAILED: {exc}[/red]")
                failed = [t for t in failed if (t[0], t[1]) != key]
                failed.append((site.resname, site.resnum, str(exc)))
            _save()  # persist progress after every site (resume-safe)

        ws.set("pb_titrate_self_energies", str(pkl))

        # Compact per-site pKa map for downstream consumers (cpin generation,
        # session-replay UI). Keyed by (resname, resnum) → pKa_eff (float) or
        # None if compute_pka_multistate couldn't extract one.
        pka_map: Dict[Tuple[str, int], Optional[float]] = {}
        for key, r in results.items():
            pka_map[key] = r.get("pKa_eff")
        ws.set("pb_titrate_pka", pka_map)

        # Human-readable CSV of the per-site pKas (the pickle is the authoritative
        # record but isn't directly viewable). Sorted by residue number.
        csv_path = self._output_dir() / "single_site_pka.csv"
        with csv_path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["resname", "resnum", "pKa_eff", "dominant_state_idx"])
            for (rn, num) in sorted(results, key=lambda k: (k[1], k[0])):
                r = results[(rn, num)]
                pka = r.get("pKa_eff")
                w.writerow([rn, num,
                            "" if pka is None else f"{pka:.4f}",
                            r.get("dominant_state_idx")])
        ws.set("pb_titrate_pka_csv", str(csv_path))
        self.console.print(f"  [green]pKa table written:[/green] {csv_path}")

        n_ok = len(results)
        if n_ok == 0:
            raise RuntimeError(
                f"All {len(sites)} single-site runs failed; aborting workflow.")
        return {"summary": f"{n_ok}/{len(sites)} sites; {len(failed)} failed"}

    # -- step 4: coupling precompute (optional) -----------------------------

    def _load_unconverged_sites(self) -> List[Tuple[str, int]]:
        """Seed set for targeted coupling: the titratable sites mean-field
        (step 5) could not settle, persisted as ``pb_titrate_unconverged_sites``
        (a list of ``[resname, resnum]``). Empty when the last solve converged
        or never ran."""
        raw = self._ws().get("pb_titrate_unconverged_sites") or []
        return [(r[0], int(r[1])) for r in raw]

    def _restrict_to_subgraph(self, coupling_sites, seed_keys, radius,
                              prmtop, rst7):
        """Keep only the seed sites plus every coupling site whose nearest
        heavy atom lies within ``radius`` Å of the nearest heavy atom of any
        seed site (closest-approach, not centroid).

        Returns the filtered CouplingSite list. If ``seed_keys`` is empty the
        full list is returned unchanged (nothing to target)."""
        import parmed as pmd

        seed = {(rn, int(num)) for rn, num in seed_keys}
        if not seed:
            return coupling_sites

        struct = pmd.load_file(str(prmtop), str(rst7))

        # Heavy-atom coordinates per (resname, resnum). resnum is 1-based to
        # match CouplingSite.resnum / the rest of pb_titrate (r.number + 1).
        heavy: Dict[Tuple[str, int], List[Tuple[float, float, float]]] = {}
        for r in struct.residues:
            key = (r.name, r.number + 1)
            coords = [(a.xx, a.xy, a.xz) for a in r.atoms if a.atomic_number > 1]
            if coords:
                heavy[key] = coords

        seed_coords: List[Tuple[float, float, float]] = []
        for k in seed:
            seed_coords.extend(heavy.get(k, []))
        if not seed_coords:
            self.console.print(
                "  [yellow]Targeted coupling: none of the unconverged sites were "
                "found in the structure; running the full site set.[/yellow]")
            return coupling_sites

        r2 = float(radius) * float(radius)

        def within(key) -> bool:
            for (ax, ay, az) in heavy.get(key, []):
                for (bx, by, bz) in seed_coords:
                    if (ax-bx)**2 + (ay-by)**2 + (az-bz)**2 <= r2:
                        return True
            return False

        kept = [s for s in coupling_sites
                if (s.resname, s.resnum) in seed or within((s.resname, s.resnum))]
        return kept

    def _checklist_pbt_4_coupling(self) -> Dict[str, Any]:
        from .coupling import (
            make_coupling_site, compute_self_energies,
            compute_pair_couplings, coupling_summary,
            estimate_self_energy_pb_calls, estimate_pair_pb_calls,
        )

        ws = self._ws()
        params = self._params()
        prmtop = Path(ws.get("pb_titrate_prmtop"))
        rst7   = Path(ws.get("pb_titrate_rst7"))

        self.console.print(f"\n[bold cyan]Coupling matrix precompute[/bold cyan]")
        if not confirm_with_context(
                self.processor,
                "Run pairwise coupling? (Optional — only needed for "
                "Monte Carlo / exact enumeration. Skip for mean-field.)",
                default=True,
                module=MODULE_NAME,
                description="Run coupling precompute"):
            return {"summary": "skipped (mean-field-only)"}

        sites_serialized = ws.get("pb_titrate_sites") or []
        coupling_sites = [make_coupling_site(s["resname"], s["resnum"])
                            for s in sites_serialized]

        # Targeted coupling: restrict to the unconverged (flipper) sites plus
        # their near neighbors when mean-field (step 5) did not settle. This
        # turns an all-sites N² coupling run into a small subgraph. The sites
        # dropped from the subgraph are frozen at their converged mean-field
        # state in the PB background (see the extra_background block below).
        targeted = bool(params.get("coupling_targeted", False))
        frozen_outside: List[Tuple[str, int]] = []
        if targeted:
            seed = self._load_unconverged_sites()
            if not seed:
                self.console.print(
                    "  [yellow]Targeted coupling is on but no unconverged sites "
                    "are recorded (run mean-field at step 5 first, or it "
                    "converged). Running the full site set.[/yellow]")
            else:
                radius = float(params.get("coupling_subgraph_radius", 12.0))
                keep_keys = {(s.resname, s.resnum) for s in
                             self._restrict_to_subgraph(
                                 coupling_sites, seed, radius, prmtop, rst7)}
                n_before = len(coupling_sites)
                frozen_outside = [(s.resname, s.resnum) for s in coupling_sites
                                  if (s.resname, s.resnum) not in keep_keys]
                coupling_sites = [s for s in coupling_sites
                                  if (s.resname, s.resnum) in keep_keys]
                self.console.print(
                    f"  [cyan]Targeted coupling: {len(coupling_sites)}/{n_before} "
                    f"sites (seed {len(seed)} unconverged + neighbors within "
                    f"{radius:.0f} Å); {len(frozen_outside)} other titratable "
                    f"site(s) frozen at their converged state.[/cyan]")

        site_envelope_objs = self._site_objects()
        site_envelopes = {(s.resname, s.resnum): s for s in site_envelope_objs}

        wd = self._output_dir() / "coupling"
        pbsa_kw = self._make_pbsa_kwargs()

        # ---- Phase 1: self-energies -------------------------------------
        # Reuse step 3's per-site state_energies when available. Step 3 ran
        # the same PB calculations (same all-neutral background, same
        # cluster radius, same site envelopes, same pbsa params), so its
        # cached state_energies dicts substitute directly here and the
        # worker pool is spared from redoing them.
        precomputed: Dict[Tuple[str, int], Dict[int, Dict[str, Any]]] = {}
        step3_pkl = ws.get("pb_titrate_self_energies")
        if step3_pkl and Path(step3_pkl).exists():
            try:
                with open(step3_pkl, "rb") as f:
                    step3_data = pickle.load(f)
                for key, result in (step3_data.get("results") or {}).items():
                    se_dict = result.get("state_energies") if isinstance(result, dict) else None
                    if se_dict:
                        precomputed[key] = se_dict
            except Exception as exc:
                self.console.print(
                    f"  [yellow]Could not load step 3 cache "
                    f"({Path(step3_pkl).name}): {exc}; recomputing all "
                    f"self-energies.[/yellow]")
                precomputed = {}

        precomputed_keys = {s.key for s in coupling_sites if s.key in precomputed}
        n_cached = len(precomputed_keys)
        n_self_pb = estimate_self_energy_pb_calls(
            coupling_sites, precomputed_keys=precomputed_keys)
        if n_cached > 0:
            self.console.print(
                f"\n[bold]Phase 1/2 — self-energies[/bold]: "
                f"{len(coupling_sites)} sites "
                f"([green]{n_cached} reused from step 3[/green], "
                f"{len(coupling_sites) - n_cached} fresh PB), "
                f"~{n_self_pb} PB calls ({params['workers']} workers)")
        else:
            self.console.print(
                f"\n[bold]Phase 1/2 — self-energies[/bold]: "
                f"{len(coupling_sites)} sites, ~{n_self_pb} PB calls "
                f"({params['workers']} workers)")

        def _self_progress(done: int, total: int, label: str):
            self.console.print(f"  [{done}/{total}] {label}")

        # Metal-coordinating residues are excluded from titration but must be
        # held at their fixed metal-bound charge in every self-energy PB call
        # (mirrors the pbt-3 background). They are not CouplingSites, so they
        # only enter as extra background. (Pair couplings W_ij are a second
        # difference and cancel any fixed background, so they need no analog.)
        from .residues import get_residue
        extra_background = {}
        for entry in (ws.get("pb_titrate_metal_fixed") or []):
            res = get_residue(entry["resname"])
            state = (res.deprot_state if entry["fix_mode"] == "deprotonated"
                     else res.neutral_state)
            extra_background[(entry["resname"], entry["resnum"])] = state
        if extra_background:
            self.console.print(
                f"  [grey50]{len(extra_background)} metal-coordinating residue(s) "
                f"held at metal-bound charge in the coupling background.[/grey50]")

        # Targeted mode: freeze every titratable site outside the subgraph at
        # its converged mean-field state, so the subgraph's self-energies see
        # the correct fixed environment (not the default reference state).
        if frozen_outside:
            state_map = ws.get("pb_titrate_state_map") or {}
            n_frozen = 0
            n_missing = 0
            for key in frozen_outside:
                name = state_map.get(key)
                if name is None:
                    name = state_map.get(f"{key[0]}-{key[1]}")  # defensive
                if name is None:
                    n_missing += 1
                    continue
                res = get_residue(key[0])
                state = next((st for st in res.states if st.name == name), None)
                if state is None:
                    n_missing += 1
                    continue
                extra_background[key] = state
                n_frozen += 1
            self.console.print(
                f"  [grey50]{n_frozen} out-of-subgraph titratable site(s) frozen at "
                f"their converged state in the background"
                + (f"; {n_missing} had no recorded state (left at reference)."
                   if n_missing else ".") + "[/grey50]")

        se = compute_self_energies(
            coupling_sites, prmtop, rst7,
            workers=params["workers"], work_root=wd,
            pbsa_params=pbsa_kw,
            cluster_radius=params.get("cluster_radius"),
            site_envelopes=site_envelopes,
            model_cache_dir=self._output_dir() / "model_compounds",
            extra_background=extra_background or None,
            precomputed_state_energies=precomputed if precomputed else None,
            on_progress=_self_progress)
        self.console.print(f"  ✓ self-energies computed")

        # ---- Phase 2: pair couplings ------------------------------------
        # In targeted mode the configured subgraph radius doubles as the
        # pair-distance cutoff, overriding the ionic-strength auto value
        # (which is 3·κ⁻¹ and can be tens of Å at low salt). None = auto.
        pair_cutoff = (params.get("coupling_distance_cutoff")
                       if targeted else None)
        pair_plan = estimate_pair_pb_calls(
            coupling_sites, prmtop, rst7,
            cluster_radius=params.get("cluster_radius"),
            distance_cutoff=pair_cutoff,
            pbsa_params=pbsa_kw)
        cutoff_str = (f"{pair_plan['distance_cutoff_A']:.1f} Å"
                       if pair_plan['distance_cutoff_A'] != float("inf")
                       else "∞ (full structure)")
        self.console.print(
            f"\n[bold]Phase 2/2 — pair couplings[/bold]: "
            f"{pair_plan['n_pairs_kept']}/{pair_plan['n_pairs_total']} "
            f"pairs within {cutoff_str} "
            f"({pair_plan['n_pairs_skipped']} skipped beyond cutoff), "
            f"~{pair_plan['n_pb_calls']} PB calls")

        # Live PB-call counter; print every Nth call to avoid log spam.
        stride = max(1, pair_plan['n_pb_calls'] // 50)
        def _pair_progress(done: int, total: int, label: str):
            if done == total or done % stride == 0:
                self.console.print(f"  [{done}/{total}] {label}")

        W = compute_pair_couplings(
            coupling_sites, prmtop, rst7,
            workers=params["workers"], work_root=wd,
            pbsa_params=pbsa_kw,
            cluster_radius=params.get("cluster_radius"),
            distance_cutoff=pair_cutoff,
            site_envelopes=site_envelopes,
            on_progress=_pair_progress)

        self.console.print(f"\n{coupling_summary(W)}")

        pkl = self._output_dir() / "coupling.pkl"
        with pkl.open("wb") as f:
            pickle.dump({"sites": coupling_sites, "self_energies": se,
                         "coupling": W}, f)
        ws.set("pb_titrate_coupling", str(pkl))

        return {"summary": f"coupling matrix for {len(coupling_sites)} sites"}

    # -- step 5: solve ------------------------------------------------------

    def _checklist_pbt_5_solve(self) -> Dict[str, Any]:
        ws = self._ws()
        params = self._params()

        self.console.print(f"\n[bold cyan]Solve for state populations[/bold cyan]")

        coupling_pkl = ws.get("pb_titrate_coupling")
        if not coupling_pkl or not Path(coupling_pkl).exists():
            # Mean-field via live PB iteration is the path when coupling
            # wasn't precomputed.
            return self._solve_mean_field_live()

        # Coupling precomputed — pick a solver against the static matrix
        with open(coupling_pkl, "rb") as f:
            data = pickle.load(f)
        from .select_solver import (recommend, total_states,
                                     ENUMERATE_MAX_STATES, CMF_MAX_CLUSTER_STATES,
                                     cluster_feasibility_table)
        rec = recommend(data["self_energies"], data["coupling"], pH=params["pH"])
        n_total = total_states(data["sites"])
        self.console.print(f"  recommended: [bold]{rec.solver}[/bold]")
        self.console.print(f"  reason: {rec.reason}")

        # cluster_mean_field feasibility table — show the threshold↔cluster-size
        # tradeoff up front so the user isn't surprised by a crash mid-run (the
        # failure mode: strongly-coupled sites percolate into one huge cluster).
        cmf_feas = cluster_feasibility_table(
            data["sites"], data["coupling"])
        cmf_ever_feasible = any(r["feasible"] for r in cmf_feas)
        cmf_feasible_at_1kT = next(
            (r["feasible"] for r in cmf_feas if r["threshold_kT"] == 1.0), False)

        # List the available solvers explicitly. The prompt itself only
        # echoes the default; without this block the user sees just
        # "Solver to run (auto): " and has no idea what the alternatives are.
        enum_feasible = n_total <= ENUMERATE_MAX_STATES
        enum_label = (f"Exact enumeration ({n_total:,} total states)"
                      if enum_feasible
                      else f"Exact enumeration (UNFEASIBLE: {n_total:,} > "
                           f"{ENUMERATE_MAX_STATES:,} state limit)")
        if cmf_feasible_at_1kT:
            cmf_label = "exact within tightly coupled clusters, MF between"
        elif cmf_ever_feasible:
            thr = next(r["threshold_kT"] for r in cmf_feas if r["feasible"])
            cmf_label = (f"[yellow]needs threshold ≥ {thr:.0f} kT to be feasible "
                         f"(demotes coupling — MC is better)[/yellow]")
        else:
            cmf_label = ("[red]INFEASIBLE: clusters too large to enumerate at "
                         "any threshold — use Monte Carlo[/red]")
        self.console.print("\n  Available solvers:")
        self.console.print(f"    [cyan]auto[/cyan]                "
                              f"→ run the recommendation above ({rec.solver})")
        self.console.print(f"    [cyan]mean_field[/cyan]          "
                              f"→ self-consistent iteration; per-site marginals only "
                              f"(no joint correlations)")
        self.console.print(f"    [cyan]cluster_mean_field[/cyan]  "
                              f"→ {cmf_label}")
        self.console.print(f"    [cyan]monte_carlo[/cyan]         "
                              f"→ Metropolis sampling; gives joint distribution; "
                              f"built-in convergence check available")
        self.console.print(f"    [cyan]enumerate[/cyan]           "
                              f"→ {enum_label}\n")

        choice = prompt_with_context(
            self.processor,
            "Solver to run",
            choices=["auto", "mean_field", "cluster_mean_field",
                     "monte_carlo", "enumerate"],
            default="auto",
            module=MODULE_NAME,
            description="Solver selection",
            options_map={"auto": rec.solver,
                          "mean_field": "Mean-field iteration",
                          "cluster_mean_field":
                              "Cluster mean-field (exact within tightly "
                              "coupled clusters, MF between)",
                          "monte_carlo": "Monte Carlo (Metropolis)",
                          "enumerate":
                              f"Exact enumeration (n_total={n_total:,}, "
                              f"limit {ENUMERATE_MAX_STATES:,})"})
        solver = rec.solver if choice == "auto" else choice

        if solver == "monte_carlo":
            state_map = self._run_mc_solver(data, params, n_total)
        elif solver == "enumerate":
            from .solvers.enumerate import enumerate_all
            result = enumerate_all(data["self_energies"], data["coupling"],
                                    pH=params["pH"])
            state_map = self._state_map_from_marginals(
                data["sites"], result.marginals)
            self.console.print(f"  ✓ enumerated {result.n_states_enumerated} states")
        elif solver == "cluster_mean_field":
            state_map = self._run_cmf_solver(data, params, n_total, cmf_feas)
            if state_map is None:
                # CMF was infeasible and the user opted to fall back to MC.
                self.console.print(
                    "\n  [cyan]Falling back to Monte Carlo.[/cyan]")
                state_map = self._run_mc_solver(data, params, n_total)
                solver = "monte_carlo"
        else:
            # mean_field on the precomputed matrix
            state_map = self._solve_mean_field_static(data, params["pH"])

        # What did coupling change? Compare the coupled solve to the mean-field
        # baseline (the per-site states already in the workspace from pbt-5's
        # mean-field run, if any) over the solved sites.
        if solver != "mean_field":
            self._report_coupling_diff(ws, state_map)

        # The coupled solvers run only over the sites in coupling.pkl. In
        # targeted mode that is the flipper subgraph (a subset), so MERGE the
        # solved assignments over the existing full map rather than replacing
        # it — otherwise the out-of-subgraph sites (held at their converged
        # mean-field state) would be silently dropped. In a full (non-targeted)
        # run coupling.pkl covers every site, so the merge is a no-op overwrite.
        existing = ws.get("pb_titrate_state_map") or {}
        merged = dict(existing)
        merged.update(state_map)
        ws.set("pb_titrate_state_map", merged)
        ws.set("pb_titrate_solver_result", solver)
        n_solved = len(state_map)
        n_total_map = len(merged)
        if n_total_map > n_solved:
            self.console.print(
                f"  [grey50]Merged {n_solved} coupled-solver assignment(s) into the "
                f"full {n_total_map}-site map ({n_total_map - n_solved} sites kept "
                f"at their mean-field state).[/grey50]")
        return {"summary": (f"solver={solver}, {n_solved} solved, "
                            f"{n_total_map} total states assigned")}

    def _solve_mean_field_live(self) -> Dict[str, Any]:
        """Mean-field iteration via live PB calls (no precomputed matrix)."""
        from .solvers.mean_field import make_site, solve

        ws = self._ws()
        params = self._params()
        prmtop = Path(ws.get("pb_titrate_prmtop"))
        rst7   = Path(ws.get("pb_titrate_rst7"))

        sites_serialized = ws.get("pb_titrate_sites") or []
        sites = [make_site(s["resname"], s["resnum"], pH=params["pH"])
                  for s in sites_serialized]
        max_iter = int_prompt_with_context(
            self.processor, "Max mean-field iterations",
            default=5, module=MODULE_NAME,
            description="Max mean-field iterations")

        # Cluster cutout + fixed-charge background, same as the single-site
        # (pbt-3) and coupling (pbt-4) paths. Without these the worker ran PB on
        # the FULL structure (NATYP overflow -> "rdparm: a parameter array
        # overflowed"); these keep each PB call to a single-site cutout.
        from .residues import get_residue
        site_envelopes = {(s.resname, s.resnum): s
                          for s in self._site_objects()}
        extra_background = {}
        for entry in (ws.get("pb_titrate_metal_fixed") or []):
            res = get_residue(entry["resname"])
            state = (res.deprot_state if entry["fix_mode"] == "deprotonated"
                     else res.neutral_state)
            extra_background[(entry["resname"], entry["resnum"])] = state

        def _run(n_iter):
            # solve() mutates `sites` in place (each carries its .state), so a
            # follow-up call with the same list continues from the current
            # assignment rather than restarting.
            return solve(
                sites=sites,
                cluster_prmtop=prmtop,
                cluster_rst7=rst7,
                pH=params["pH"], max_iter=n_iter,
                workers=params["workers"],
                work_root=self._output_dir() / "mean_field",
                pbsa_params=self._make_pbsa_kwargs(),
                cluster_radius=params.get("cluster_radius"),
                site_envelopes=site_envelopes,
                extra_background=extra_background or None,
                model_cache_dir=self._output_dir() / "model_compounds",
            )

        result = _run(max_iter)
        total_iters = result.get("iterations", max_iter)

        # Accumulate the flipper set across every _run() call. solve() returns
        # the union of changed sites over its OWN last few sweeps; a period-2
        # oscillator shows only ~half its members in any single sweep, and each
        # _run() resets that window, so we union across calls too. Over-
        # inclusion is harmless (a few extra coupled sites in the subgraph);
        # under-inclusion would leave a flipper uncoupled. Empty when converged.
        flipper_union: set = set()
        if not result.get("converged", False):
            flipper_union.update(
                tuple(k) for k in result.get("flipper_sites",
                                             result.get("last_changed_sites", [])))

        # If mean-field hit max_iter without converging, offer to continue (it
        # resumes from the current state) or accept the partial result. A
        # genuinely frustrated/strongly-coupled subset oscillates forever no
        # matter how many iterations are added — the right tool for those is the
        # coupling matrix (step 4) + a coupled solver (Monte Carlo / cluster
        # mean-field) at step 5, which capture correlated flips mean-field can't.
        while not result.get("converged", False):
            n_left = result.get("last_n_changes", 0)
            flippers = result.get("flipper_sites",
                                  result.get("last_changed_sites", []))
            self.console.print(
                f"\n  [yellow]Mean-field did NOT converge after "
                f"{total_iters} iteration(s): {n_left} site(s) still "
                f"changing.[/yellow]")
            if flippers:
                shown = ", ".join(f"{rn}-{num}" for rn, num in flippers[:20])
                more = f" (+{len(flippers)-20} more)" if len(flippers) > 20 else ""
                self.console.print(f"  [grey50]Still flipping: {shown}{more}[/grey50]")
            self.console.print(
                "  [grey50]Typically a strongly-coupled subset. Adding iterations "
                "helps only if the count is still decaying; a frustrated cluster "
                "oscillates forever. For those, accept here and run Coupling "
                "Precompute (step 4) + a coupled solver (Monte Carlo / cluster "
                "mean-field) at step 5.[/grey50]")
            more_iter = int_prompt_with_context(
                self.processor,
                "  Run how many more mean-field iterations? (0 = accept current)",
                default=0, module=MODULE_NAME,
                description="Additional mean-field iterations")
            if more_iter <= 0:
                break
            result = _run(more_iter)
            total_iters += result.get("iterations", more_iter)
            if result.get("converged", False):
                flipper_union.clear()  # settled — no flippers to record
            else:
                flipper_union.update(
                    tuple(k) for k in result.get("flipper_sites",
                                                 result.get("last_changed_sites", [])))

        # solve() returns an updated `sites` list with .state set on each
        state_map = {(s.resname, s.resnum): s.state.name
                      for s in result["sites"]}
        ws.set("pb_titrate_state_map", state_map)
        ws.set("pb_titrate_solver_result", "mean_field_live")

        # Persist the still-flipping (unconverged) sites so a later coupling run
        # (step 4) can restrict itself to them in targeted mode. This is the
        # UNION over every sweep window across all _run() calls, so period-2
        # oscillators are fully captured. Cleared to [] on convergence. Stored
        # as [[resname, resnum], ...], ordered by residue number.
        if result.get("converged", False) or not flipper_union:
            ws.set("pb_titrate_unconverged_sites", [])
        else:
            unconverged = [[rn, int(num)] for rn, num in
                           sorted(flipper_union, key=lambda k: (k[1], k[0]))]
            ws.set("pb_titrate_unconverged_sites", unconverged)
            self.console.print(
                f"  [grey50]Recorded {len(unconverged)} unconverged site(s) "
                f"(union across all sweeps) for targeted coupling — enable it "
                f"in step 1, then run step 4.[/grey50]")

        conv = "converged" if result.get("converged") else "NOT converged"
        return {"summary": (f"mean_field (live), {len(state_map)} states "
                              f"assigned, {conv} after {total_iters} iter")}

    def _solve_mean_field_static(self, data, pH: float,
                                    max_iter: int = 50) -> Dict[Tuple[str, int], str]:
        """Mean-field iteration on a precomputed self-energy + coupling matrix.

        For each site i, given the currently assigned chemistry of every
        other site, compute G_eff(i, chem) = ddG(i, chem) + ph_term(i, chem)
        + Σ_{j≠i} W((i, chem), (j, chem_j)) and pick the lowest-G_eff
        chemistry. Repeat until no site changes (or max_iter is hit).

        Initial state for every site is its ref_chem (= the residue's
        neutral form), keeping this consistent with single-site, coupling,
        and live mean-field.
        """
        from .pb_backend import RT_LN10
        sites = data["sites"]
        self_e = data["self_energies"]
        W = data["coupling"]

        # Current chemistry per site; start at the neutral reference.
        current: Dict[Tuple[str, int], Tuple[int, float]] = {
            site.key: site.ref_chem for site in sites
        }

        for it in range(max_iter):
            changed = 0
            for site in sites:
                best_chem = None
                best_g = float("inf")
                for chem in site.chemistries:
                    ddg = self_e.intrinsic_ddG.get((site.key, chem), 0.0)
                    ph_term = RT_LN10 * (chem[0] * pH - chem[1])
                    pair_term = 0.0
                    for other in sites:
                        if other.key == site.key:
                            continue
                        chem_j = current[other.key]
                        pair_term += W.W.get((site.key, chem), {}).get(
                            (other.key, chem_j), 0.0)
                    g_eff = ddg + ph_term + pair_term
                    if g_eff < best_g:
                        best_g = g_eff
                        best_chem = chem
                if best_chem != current[site.key]:
                    current[site.key] = best_chem
                    changed += 1
            if changed == 0:
                break

        out: Dict[Tuple[str, int], str] = {}
        for site in sites:
            state = site.chem_to_state[current[site.key]]
            out[(site.resname, site.resnum)] = state.name
        return out

    # -- pbt-5 solver helpers (richer UX) -----------------------------------

    AMBIGUOUS_MARGINAL = 0.60   # dominant-state probability below this = unsure

    def _run_mc_solver(self, data, params, n_total: int
                       ) -> Dict[Tuple[str, int], str]:
        """Run Monte Carlo, print a per-site dominant-state + confidence table
        (flagging ambiguous sites), then optionally run a multi-seed
        convergence check and adopt its consensus. Returns the state map."""
        from .solvers.monte_carlo import run_mcmc, check_convergence

        # Default to a generous chain — MC is pure arithmetic on the
        # precomputed matrix (no PB calls), so 1e6 steps run in seconds.
        n_steps = int_prompt_with_context(
            self.processor, "MC steps",
            default=1_000_000, module=MODULE_NAME,
            description="Monte Carlo steps")
        n_equil = int_prompt_with_context(
            self.processor, "Equilibration steps",
            default=n_steps // 10, module=MODULE_NAME,
            description="MC equilibration")
        result = run_mcmc(data["self_energies"], data["coupling"],
                            pH=params["pH"], n_steps=n_steps, n_equil=n_equil)

        # Acceptance ratio — reframed. It reflects landscape ruggedness (how
        # many proposed flips are uphill), NOT result quality. A low ratio on a
        # strongly-coupled system is expected and does not mean the marginals
        # are wrong; convergence is what matters, checked separately below.
        acc = result.acceptance_ratio * 100.0
        self.console.print(
            f"\n  MC: {n_steps:,} steps, acceptance [bold]{acc:.1f}%[/bold] "
            f"[grey50](landscape ruggedness — low is normal for strong coupling; "
            f"NOT a quality score)[/grey50]")

        self._print_marginal_table(result.sites, result.marginals,
                                    title="MC dominant states")

        # Built-in convergence check — the automated "run again, compare"
        # test the user previously had to do by hand.
        if confirm_with_context(
                self.processor,
                "Run a convergence check (independent replicate chains, "
                "compare per-site states)?",
                default=True, module=MODULE_NAME,
                description="Run MC convergence check"):
            n_rep = int_prompt_with_context(
                self.processor, "  Number of replicate chains",
                default=4, module=MODULE_NAME,
                description="MC convergence replicate count")

            def _cprog(done, total, seed, ratio):
                self.console.print(
                    f"    replicate {done}/{total} (seed {seed}): "
                    f"acceptance {ratio*100:.1f}%")

            conv = check_convergence(
                data["self_energies"], data["coupling"],
                pH=params["pH"], n_steps=n_steps, n_equil=n_equil,
                n_replicates=n_rep, on_progress=_cprog)

            if conv.all_agree:
                self.console.print(
                    f"  [green]✓ Converged: all {n_rep} replicate chains agree "
                    f"on the dominant state for every site.[/green]")
            else:
                self.console.print(
                    f"  [yellow]⚠ {len(conv.unstable_sites)} site(s) disagree "
                    f"across replicate chains — undersampled or near-"
                    f"degenerate:[/yellow]")
                for key in sorted(conv.unstable_sites, key=lambda k: k[1]):
                    site = next(s for s in result.sites if s.key == key)
                    votes = conv.per_site_dominant[key]
                    labels = [site.residue.state_label(site.chem_to_state[c])
                              for c in votes]
                    self.console.print(
                        f"      {key[0]}-{key[1]}: votes {labels} "
                        f"(agreement {conv.per_site_agreement[key]*100:.0f}%, "
                        f"mean p={conv.mean_marginal[key]:.2f})")
                self.console.print(
                    "  [grey50]Consensus (modal) state adopted for these; review "
                    "and set them by hand in the next step if it matters.[/grey50]")

            # Adopt the consensus map (identical to the single-run map when all
            # chains agree; resolves ties by majority vote otherwise).
            return {(s.resname, s.resnum):
                    s.chem_to_state[conv.consensus[s.key]].name
                    for s in result.sites}

        return self._state_map_from_marginals(result.sites, result.marginals)

    def _run_cmf_solver(self, data, params, n_total: int, cmf_feas
                        ) -> Optional[Dict[Tuple[str, int], str]]:
        """Run cluster mean-field. On the 'cluster too large to enumerate'
        failure, show the threshold↔feasibility table and let the user raise
        the threshold or fall back to Monte Carlo (returns None to signal
        fallback)."""
        from .solvers.cluster_mean_field import solve as cmf_solve

        threshold = float_prompt_with_context(
            self.processor,
            "Cluster threshold (kT) — pairs with |W| above this merge into "
            "the same cluster (higher = smaller clusters)",
            default=1.0, module=MODULE_NAME,
            description="Cluster MF coupling-graph edge threshold (kT)")
        while True:
            try:
                result = cmf_solve(
                    data["self_energies"], data["coupling"],
                    pH=params["pH"], threshold_kT=threshold)
            except ValueError as exc:
                self.console.print(f"\n  [red]Cluster mean-field cannot run: "
                                   f"{exc}[/red]")
                self.console.print(
                    "\n  [bold]Threshold → largest-cluster feasibility:[/bold]")
                for r in cmf_feas:
                    mark = ("[green]feasible[/green]" if r["feasible"]
                            else "[red]too big[/red]")
                    self.console.print(
                        f"    {r['threshold_kT']:>4.0f} kT  →  "
                        f"{r['largest_cluster_sites']:>3d} sites, "
                        f"{r['largest_joint_states']:.1e} states   {mark}")
                self.console.print(
                    "  [grey50]Raising the threshold breaks big clusters apart but "
                    "demotes those couplings to mean-field level — if you have "
                    "to raise it far, Monte Carlo is the better tool.[/grey50]")
                feasible = [r for r in cmf_feas if r["feasible"]]
                if feasible and confirm_with_context(
                        self.processor,
                        f"Retry at a higher threshold "
                        f"(≥ {feasible[0]['threshold_kT']:.0f} kT needed)?",
                        default=False, module=MODULE_NAME,
                        description="Retry cluster MF at higher threshold"):
                    threshold = float_prompt_with_context(
                        self.processor, "  new threshold (kT)",
                        default=float(feasible[0]['threshold_kT']),
                        module=MODULE_NAME,
                        description="Retry cluster MF threshold (kT)")
                    continue
                return None   # caller falls back to Monte Carlo
            break

        n_singletons = sum(1 for c in result.cluster_sizes if c == 1)
        n_multi = len(result.cluster_sizes) - n_singletons
        big_clusters = sorted(result.cluster_sizes, reverse=True)[:5]
        self.console.print(
            f"  ✓ cluster MF: {n_multi} multi-site cluster(s) "
            f"(sizes {big_clusters}), {n_singletons} singletons; "
            f"converged in {result.n_iterations} iterations "
            f"({'yes' if result.converged else 'NO — hit max_iter'})")
        self._print_marginal_table(result.sites, result.marginals,
                                    title="Cluster mean-field dominant states")
        return self._state_map_from_marginals(result.sites, result.marginals)

    def _print_marginal_table(self, sites, marginals, *, title: str) -> None:
        """Per multi-state site: dominant state label + marginal probability,
        flagging ambiguous (near-50/50) sites. Single-chemistry sites are
        omitted (nothing to decide)."""
        rows = []
        for site in sites:
            if site.n_chemistries() <= 1:
                continue
            marg = marginals[site.key]
            best = max(marg, key=marg.get)
            p = marg[best]
            label = site.residue.state_label(site.chem_to_state[best])
            rows.append((site.resname, site.resnum, label, p,
                         p < self.AMBIGUOUS_MARGINAL))
        if not rows:
            return
        n_amb = sum(1 for r in rows if r[4])
        self.console.print(f"\n  [bold]{title}[/bold] "
                           f"({len(rows)} multi-state sites, "
                           f"{n_amb} ambiguous):")
        for rn, num, label, p, amb in sorted(rows, key=lambda r: r[1]):
            flag = " [yellow]← ambiguous[/yellow]" if amb else ""
            color = "yellow" if amb else "white"
            self.console.print(
                f"      {rn}-{num}: [{color}]{label}[/{color}] "
                f"p={p:.2f}{flag}")
        if n_amb:
            self.console.print(
                f"  [grey50]{n_amb} site(s) below p={self.AMBIGUOUS_MARGINAL:.2f} "
                f"are genuinely split — consider setting them by hand in the "
                f"review step.[/grey50]")

    def _report_coupling_diff(self, ws, state_map) -> None:
        """Show how the coupled solve differs from the mean-field baseline
        already in the workspace, over the sites the coupled solver touched.
        Makes the value of running coupling visible."""
        baseline = ws.get("pb_titrate_state_map") or {}
        if not baseline:
            return
        changed = []
        for key, new_name in state_map.items():
            old = baseline.get(key)
            if old is not None and old != new_name:
                changed.append((key, old, new_name))
        n = len(state_map)
        if not changed:
            self.console.print(
                f"  [grey50]Coupling vs mean-field: identical over all {n} "
                f"solved site(s) — coupling confirmed the mean-field states."
                f"[/grey50]")
            return
        self.console.print(
            f"  [cyan]Coupling changed {len(changed)}/{n} solved site(s) vs "
            f"the mean-field baseline:[/cyan]")
        from .residues import get_residue

        def _labeled(resname, state_name):
            """Append the human label (HIP/HID/PROT/DEPROT…) to a STATE k name."""
            try:
                res = get_residue(resname)
                lbl = res.state_label(res.state_by_name(state_name))
                return f"{state_name} ({lbl})"
            except Exception:
                return state_name

        for (rn, num), old, new in sorted(changed, key=lambda r: r[0][1]):
            self.console.print(
                f"      {rn}-{num}: {_labeled(rn, old)} → {_labeled(rn, new)}")

    def _state_map_from_marginals(self, sites, marginals
                                    ) -> Dict[Tuple[str, int], str]:
        """Pick the dominant state per site from MC/enumerate marginals."""
        out: Dict[Tuple[str, int], str] = {}
        for site in sites:
            marg = marginals[site.key]
            best_chem = max(marg, key=marg.get)
            state = site.chem_to_state[best_chem]
            out[(site.resname, site.resnum)] = state.name
        return out

    def _report_mc_result(self, result, n_total: int) -> None:
        """Print acceptance ratio with interpretation, list locked sites,
        and suggest exact enumeration when feasible.

        Empirical bands for single-flip Metropolis MC on protonation states
        (no theoretical optimum analogous to the 0.234 result for continuous
        random-walk MH; depends on coupling strength and how many sites
        sit far from their pKa midpoint):

           > 35%   excellent — chain explores freely
          15–35%   healthy — typical for well-conditioned systems
           5–15%   borderline — run longer or check marginal convergence
           < 5%    low — strong coupling and/or sites far from midpoint;
                   marginals on dominant states still reliable, but rare-
                   state populations are noisy
        """
        from .select_solver import ENUMERATE_MAX_STATES

        ratio = result.acceptance_ratio
        ratio_pct = ratio * 100.0
        if ratio >= 0.35:
            band, color = "excellent", "green"
            interpretation = "chain exploring freely."
        elif ratio >= 0.15:
            band, color = "healthy", "green"
            interpretation = "typical for well-conditioned systems."
        elif ratio >= 0.05:
            band, color = "borderline", "yellow"
            interpretation = ("run longer or check marginal convergence "
                              "with block averages.")
        else:
            band, color = "low", "red"
            interpretation = (
                "strong coupling and/or many sites far from their titration "
                "midpoint at this pH. Marginals on dominant states are "
                "reliable; rare-state populations may be noisy.")

        self.console.print(
            f"  ✓ MC: acceptance [{color}]{ratio_pct:.1f}% ({band})[/{color}] "
            f"— {interpretation}")

        # Surface "locked" sites that never visited a second chemistry
        # during post-equilibration sampling. These are usually fine
        # (their state is so heavily preferred there's nothing to sample),
        # but seeing them helps the user calibrate trust in the result.
        locked = []
        for site in result.sites:
            marg = result.marginals[site.key]
            n_visited = sum(1 for v in marg.values() if v > 0.0)
            if n_visited <= 1 and len(site.chemistries) > 1:
                locked.append(f"{site.resname}-{site.resnum}")
        if locked:
            n_with_options = sum(1 for s in result.sites
                                 if len(s.chemistries) > 1)
            self.console.print(
                f"    [grey50]({len(locked)}/{n_with_options} multi-state sites "
                f"never flipped during sampling — fine if they're far from "
                f"pKa, watch out if you expected ambiguity.)[/grey50]")
            if len(locked) <= 8:
                self.console.print(f"    [grey50]locked: {', '.join(locked)}[/grey50]")

        # Suggest exact enumeration when it would be feasible — strictly
        # better than MC (no acceptance ratio, no convergence question, no
        # rare-state noise) when the joint state space is tractable.
        if n_total <= ENUMERATE_MAX_STATES:
            self.console.print(
                f"  [grey50]Tip: with {n_total:,} total joint states, "
                f"[cyan]enumerate[/cyan] would give the exact Boltzmann "
                f"distribution in seconds (no sampling noise). Re-run "
                f"step 5 and pick [cyan]enumerate[/cyan] if you want a "
                f"noise-free reference.[/grey50]")

    # -- step 6: review & override ------------------------------------------

    def _checklist_pbt_6_review(self) -> Dict[str, Any]:
        from .residues import get_residue

        ws = self._ws()
        state_map = ws.get("pb_titrate_state_map") or {}

        self.console.print(f"\n[bold cyan]Review recommended states[/bold cyan]")
        self.console.print(f"  {len(state_map)} sites; review and optionally override.")
        self.console.print(f"  [grey50]The cpin step consumes the literal "
                              f"\"STATE k\" name; the label in parentheses is "
                              f"a human-readable hint.[/grey50]\n")

        # Print table with both the cpin state name and a chemistry label.
        for (rn, num), state_name in sorted(state_map.items()):
            res = get_residue(rn)
            try:
                state = res.state_by_name(state_name)
                label = res.state_label(state)
                q = state.net_charge
                q_str = f"q={q:+.0f}" if abs(q - round(q)) < 0.01 else f"q={q:+.2f}"
            except (KeyError, ValueError):
                label = "?"
                q_str = ""
            self.console.print(
                f"    {rn}-{num:<5d}  →  {state_name:<8} "
                f"([cyan]{label}[/cyan]"
                + (f", {q_str}" if q_str else "")
                + ")")

        # Allow override
        if not confirm_with_context(
                self.processor,
                "Override any recommendations?",
                default=False,
                module=MODULE_NAME,
                description="Override recommendations"):
            return {"summary": f"{len(state_map)} states confirmed"}

        n_overrides = 0
        while True:
            line = prompt_with_context(
                self.processor,
                "Override (RESNAME:RESNUM=STATE_NAME), or empty to finish",
                default="",
                module=MODULE_NAME,
                description="Per-site state override")
            if not line.strip():
                break
            try:
                lhs, state_name = line.split("=", 1)
                rn, num = lhs.split(":", 1)
                key = (rn.strip(), int(num.strip()))
                state_name = state_name.strip()
                if key not in state_map:
                    self.console.print(f"  [yellow]not in state map: {key}[/yellow]")
                    continue
                # Validate the state name exists for this residue.
                try:
                    res = get_residue(key[0])
                    new_state = res.state_by_name(state_name)
                    label = res.state_label(new_state)
                except (KeyError, ValueError) as exc:
                    self.console.print(f"  [red]invalid state {state_name!r} "
                                          f"for {key[0]}: {exc}[/red]")
                    continue
                state_map[key] = state_name
                n_overrides += 1
                self.console.print(f"  ✓ {rn}-{num} → {state_name} "
                                      f"([cyan]{label}[/cyan])")
            except Exception as exc:
                self.console.print(f"  [red]parse error: {exc}[/red]")

        ws.set("pb_titrate_state_map", state_map)
        return {"summary": f"{len(state_map)} states ({n_overrides} overridden)"}

    # -- step 7: apply to prmtop --------------------------------------------

    def _checklist_pbt_7_apply(self) -> Dict[str, Any]:
        import parmed as pmd
        from .apply_state import build_production_prmtop, apply_state_map
        from .residues import get_residue

        ws = self._ws()
        state_map_named = ws.get("pb_titrate_state_map") or {}
        # Use the solvated prmtop for state apply if step 2 saved it —
        # waters/ions are needed here for ion rebalancing. Fall back to
        # the dry version (= original) if the input wasn't solvated.
        prmtop_str = (ws.get("pb_titrate_prmtop_solvated")
                       or ws.get("pb_titrate_prmtop"))
        rst7_str   = (ws.get("pb_titrate_rst7_solvated")
                       or ws.get("pb_titrate_rst7"))
        prmtop = Path(prmtop_str)
        rst7   = Path(rst7_str)

        self.console.print(f"\n[bold cyan]Applying state map to prmtop[/bold cyan]")

        # Translate (resname, resnum) → state_name into
        # (resname, resnum) → ProtonationState
        state_map_obj = {}
        for (rn, num), name in state_map_named.items():
            res = get_residue(rn)
            state = next((s for s in res.states if s.name == name), None)
            if state is None:
                raise RuntimeError(
                    f"Unknown state {name!r} for residue {rn}; "
                    f"valid: {[s.name for s in res.states]}")
            state_map_obj[(rn, num)] = state

        # Dry-run: load + apply locally to report the net charges *before*
        # asking the user to commit to rebalance or no-rebalance.
        s_preview = pmd.load_file(str(prmtop), str(rst7))
        q_in = sum(a.charge for a in s_preview.atoms)
        # Expected delta from the state map: sum of (target.net_charge -
        # STATE-0 net charge) for each titratable that changes.
        expected_delta = 0.0
        for (rn, num), st in state_map_obj.items():
            res = get_residue(rn)
            expected_delta += st.net_charge - res.states[0].net_charge
        apply_state_map(s_preview, state_map_obj)
        q_after_apply = sum(a.charge for a in s_preview.atoms)

        # Probe the bulk-ion pool before prompting.
        from .rebalance import _bulk_pool, _normalize_keep
        bulk_cations = _bulk_pool(s_preview, _normalize_keep(None), sign=+1)
        bulk_anions  = _bulk_pool(s_preview, _normalize_keep(None), sign=-1)
        n_cations = sum(len(v) for v in bulk_cations.values())
        n_anions  = sum(len(v) for v in bulk_anions.values())

        self.console.print(
            f"\n  Net charge of input prmtop:               {q_in:+.4f}")
        self.console.print(
            f"  Expected Δq from state map:               {expected_delta:+.4f}")
        self.console.print(
            f"  Net charge after applying state map:      {q_after_apply:+.4f}")
        residual_to_zero = -q_after_apply
        self.console.print(
            f"  Residual charge to neutralise:            {residual_to_zero:+.4f}")
        self.console.print(
            f"  Bulk ions in system: {n_cations} cation(s), {n_anions} anion(s)")
        if abs(q_after_apply - round(q_after_apply)) > 0.01:
            self.console.print(
                f"  [yellow]⚠ Net charge is not integer — state delta "
                f"semantics broken or input prmtop has fractional charges.["
                f"/yellow]")

        # Decide whether to even offer rebalance: only if the pool can plausibly
        # cover the residual.
        can_rebalance = True
        if abs(residual_to_zero) <= 0.01:
            can_rebalance = False
            self.console.print(
                f"  [green]✓ Already neutral; no rebalance needed.[/green]")
        elif residual_to_zero > 0 and n_anions == 0:
            can_rebalance = False
            self.console.print(
                f"  [yellow]No bulk anions in the system — cannot rebalance "
                f"a +{residual_to_zero:.2f} residual by ion removal.[/yellow]")
        elif residual_to_zero < 0 and n_cations == 0:
            can_rebalance = False
            self.console.print(
                f"  [yellow]No bulk cations in the system — cannot rebalance "
                f"a {residual_to_zero:.2f} residual by ion removal.[/yellow]")

        if can_rebalance:
            do_rebalance = confirm_with_context(
                self.processor,
                "Rebalance bulk ions to keep net charge integer?",
                default=True,
                module=MODULE_NAME,
                description="Ion rebalance after state apply")
        else:
            self.console.print(
                f"  [grey50]The unrebalanced prmtop will run under PME with a "
                f"non-zero monopole, which neutralises against a uniform "
                f"plasma. This is fine for many uses but introduces a small "
                f"finite-size correction (Hub & van der Vegt 2014). To rebalance "
                f"properly, solvate first or supply ions of the appropriate "
                f"sign.[/grey50]")
            do_rebalance = False
            if abs(residual_to_zero) > 0.01:
                if not confirm_with_context(
                        self.processor,
                        f"Write the prmtop with non-zero net charge "
                        f"({q_after_apply:+.2f})?",
                        default=True,
                        module=MODULE_NAME,
                        description="Confirm writing unrebalanced prmtop"):
                    raise RuntimeError("User cancelled — solvate / add ions, "
                                          "then re-run step 7.")

        out_prefix = Path(self._output_dir().parent) / (prmtop.stem + "_titrated")
        try:
            info = build_production_prmtop(
                input_prmtop=prmtop, input_rst7=rst7,
                state_map=state_map_obj,
                output_prefix=out_prefix,
                rebalance=do_rebalance,
            )
        except ValueError as exc:
            # Rebalance couldn't reach target — fall back to writing the
            # unbalanced prmtop with a clear message rather than crashing.
            self.console.print(f"  [yellow]rebalance failed: {exc}[/yellow]")
            self.console.print(f"  [yellow]Falling back to writing the "
                                  f"unrebalanced prmtop.[/yellow]")
            info = build_production_prmtop(
                input_prmtop=prmtop, input_rst7=rst7,
                state_map=state_map_obj,
                output_prefix=out_prefix,
                rebalance=False,
            )

        self.console.print(f"\n  net q in:                {info['net_q_in']:+.4f}")
        self.console.print(f"  net q after state apply: {info['net_q_after_apply']:+.4f}")
        self.console.print(f"  net q out:               {info['net_q_out']:+.4f}")
        self.console.print(f"  Δ (out − in):            {info['delta_q']:+.4f}")
        if info.get("rebalance"):
            r = info["rebalance"]
            removed = r.get("removed_per_species", {})
            if removed:
                rem_str = ", ".join(f"{n}× {sp}" for sp, n in sorted(removed.items()))
                self.console.print(f"  ions removed:            {rem_str}")
            else:
                self.console.print(f"  ions removed:            none")
        self.console.print(f"\n  prmtop: {info['output_prmtop']}")
        self.console.print(f"  rst7:   {info['output_rst7']}")

        # Descriptive aliases (kept for backward compatibility with the
        # docstring contract).
        out_prmtop_abs = str(Path(info["output_prmtop"]).resolve())
        out_rst7_abs   = str(Path(info["output_rst7"]).resolve())
        ws.set("prmtop_titrated", out_prmtop_abs)
        ws.set("rst7_titrated",   out_rst7_abs)

        # Promote to the canonical workspace keys that the rest of ProPrep
        # consumes (parm7_file / rst7_file are read by the MD manager,
        # structure_preprocessor, metal_site_parameterizer, etc.). After
        # state apply, the titrated prmtop *is* the production topology.
        ws.set("parm7_file", out_prmtop_abs)
        ws.set("rst7_file",  out_rst7_abs)

        # Append to md_structure_pairs so the MD manager's pair picker
        # offers the titrated pair under a clear name.
        existing_pairs = ws.get("md_structure_pairs") or []
        new_entry = {
            "name": Path(out_prmtop_abs).stem,
            "prmtop": out_prmtop_abs,
            "rst7":   out_rst7_abs,
        }
        if not any(p.get("prmtop") == out_prmtop_abs
                    and p.get("rst7") == out_rst7_abs
                    for p in existing_pairs):
            existing_pairs.append(new_entry)
            ws.set("md_structure_pairs", existing_pairs)

        self.console.print()
        self.console.print(
            f"  [bold cyan]Workspace updated[/bold cyan] — the active "
            f"prmtop/rst7 (workspace keys [italic]parm7_file[/italic] / "
            f"[italic]rst7_file[/italic]) now point to the titrated outputs. "
            f"Downstream modules (MD manager, etc.) will pick up the "
            f"titrated topology automatically.")
        if existing_pairs[-1] is new_entry:
            self.console.print(
                f"  [grey50]Added to md_structure_pairs as "
                f"\"{new_entry['name']}\".[/grey50]")

        return {"summary": (f"prmtop written, Δq={info['delta_q']:+.4f}, "
                              f"rebalance={do_rebalance}")}

    # -- step 8: persist recommendations ------------------------------------

    def _checklist_pbt_8_persist(self) -> Dict[str, Any]:
        from .residues import get_residue

        ws = self._ws()
        state_map_named = ws.get("pb_titrate_state_map") or {}

        self.console.print(f"\n[bold cyan]Persisting recommendations[/bold cyan]")

        # Build the recommendations dict: (resname, resnum) → {state_id, state_name, ...}
        # The cpin generation step (Topology Generator option 6) consumes
        # this to construct cpinutil's `-states` argument list.
        recommendations: Dict[Tuple[str, int], Dict[str, Any]] = {}
        for (rn, num), state_name in state_map_named.items():
            res = get_residue(rn)
            state_idx = next(
                (i for i, s in enumerate(res.states) if s.name == state_name),
                None)
            if state_idx is None:
                continue
            state = res.states[state_idx]
            recommendations[(rn, num)] = {
                "state_id":   state_idx,
                "state_name": state_name,
                "prot_count": state.prot_count,
                "pka_corr":   state.pka_corr,
                "net_charge": state.net_charge,
            }

        ws.set("titrate_recommendations", recommendations)

        # Human-readable CSV
        csv_path = Path(self._output_dir()) / "titrate_recommendations.csv"
        with csv_path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["resname", "resnum", "state_id", "state_name",
                        "prot_count", "pka_corr", "net_charge"])
            for (rn, num), rec in sorted(recommendations.items()):
                w.writerow([rn, num, rec["state_id"], rec["state_name"],
                            rec["prot_count"], rec["pka_corr"],
                            rec["net_charge"]])
        ws.set("titrate_report_csv", str(csv_path))

        self.console.print(f"  ✓ {len(recommendations)} site recommendations stored")
        self.console.print(f"  ✓ CSV: {csv_path}")
        self.console.print(
            f"\n[grey50]These will be used as initial states by the cpin "
            f"generation step (Topology Generator option 6) when the "
            f"\"Use titrate recommendations\" option is selected.[/grey50]")
        return {"summary": f"{len(recommendations)} recommendations persisted"}

    # -- step pdb: write modern-FF renamed PDB (optional) -------------------

    def _checklist_pbt_pdb_write(self) -> Dict[str, Any]:
        from .pdb_rename import build_rename_map, write_renamed_pdb
        import parmed as pmd

        ws = self._ws()
        state_map_named = ws.get("pb_titrate_state_map") or {}
        if not state_map_named:
            raise RuntimeError(
                "No recommended states found — run Solve (step 5) "
                "(and optionally Review, step 6) first.")
        metal_fixed = ws.get("pb_titrate_metal_fixed") or []

        # Source = the dry (minimized, solvent-stripped) structure: Ca²⁺ and
        # cofactors retained, coordinates = pre-PB minimization, no waters so
        # tleap re-solvates fresh with the modern water model.
        dry_prmtop = Path(ws.get("pb_titrate_prmtop"))
        dry_rst7 = Path(ws.get("pb_titrate_rst7"))

        self.console.print(
            "\n[bold cyan]Write modern-FF PDB (renamed protonation states)"
            "[/bold cyan]")
        self.console.print(f"  source: {dry_prmtop.name} (minimized, dry)")

        structure = pmd.load_file(str(dry_prmtop), str(dry_rst7))
        rename_map, report = build_rename_map(
            structure, state_map_named, metal_fixed)

        # Summarize the mapping by source and surface every warning.
        from collections import Counter
        by_source = Counter(r["source"] for r in report)
        self.console.print(
            f"  renaming {len(report)} titratable residue(s): "
            f"{by_source.get('recommended',0)} from recommendations, "
            f"{by_source.get('metal-fixed',0)} metal-fixed, "
            f"{by_source.get('fallback',0)} fallback")

        warnings_list = [r for r in report if r["warning"]]
        for r in warnings_list:
            self.console.print(
                f"  [yellow]⚠ {r['warning']}[/yellow]")

        # Build notes (e.g. heme propionate PRP/PRD → FixedpH set). Not a
        # problem — show each distinct note once so it isn't buried per-residue.
        notes_seen = set()
        for r in report:
            note = r.get("note")
            if note and note not in notes_seen:
                notes_seen.add(note)
                heme = [x for x in report if x.get("note") == note]
                resns = ", ".join(f"{x['resname']}-{x['resnum']}→{x['new_resname']}"
                                  for x in heme)
                self.console.print(f"  [cyan]ⓘ {note}[/cyan]")
                self.console.print(f"    [grey50]({resns})[/grey50]")

        fallbacks = [r for r in report if r["source"] == "fallback"]
        if fallbacks:
            self.console.print(
                f"  [yellow]{len(fallbacks)} residue(s) had no PB "
                f"recommendation and used a standard-pH default — review "
                f"above.[/yellow]")

        out_pdb = self._output_dir() / (dry_prmtop.stem.replace("_dry", "")
                                        + "_renamed.pdb")
        n = write_renamed_pdb(dry_prmtop, dry_rst7, out_pdb, rename_map)
        ws.set("pb_rename_pdb_file", str(out_pdb.resolve()))

        # Also drop a CSV of the full rename map for the record.
        import csv as _csv
        map_csv = self._output_dir() / "pb_rename_map.csv"
        with map_csv.open("w", newline="") as f:
            w = _csv.writer(f)
            w.writerow(["resname", "resnum", "new_resname", "source",
                        "warning", "note"])
            for r in sorted(report, key=lambda x: x["resnum"]):
                w.writerow([r["resname"], r["resnum"], r["new_resname"],
                            r["source"], r["warning"] or "", r.get("note") or ""])

        self.console.print(f"\n  [green]Renamed {n} residue(s).[/green]")
        self.console.print(f"  [green]PDB:[/green] {out_pdb}")
        self.console.print(f"  [green]map:[/green] {map_csv}")
        self.console.print(
            "  [grey50]Hydrogens stripped — the PB topology's H come from a "
            "different FF; tleap re-adds them per the FF you select.[/grey50]")
        self.console.print(
            "\n  [grey50]Re-enter the Topology Generator: it will auto-select "
            "this PDB (highest priority). Choose a modern protein FF "
            "(ff19SB→OPC, or ff14SB→TIP3P) and explicit solvation; tleap "
            "re-solvates with your minimized coordinates and these "
            "protonation states.[/grey50]")
        n_warn = len(warnings_list)
        return {"summary": (f"{n} residues renamed → {out_pdb.name}"
                            + (f"; {n_warn} warning(s)" if n_warn else ""))}

    # -- step cmp: compare with ProPKA (optional, on-demand) ----------------

    def _compute_effective_pka(self, ws, pb_pka):
        """Effective (coupling-aware) pKa per site for the ProPKA comparison.

        Titrates the coupled system over a pH grid using the precomputed
        coupling matrix (coupling.pkl; NO PB calls) and reads each coupled
        site's effective pKa as its protonation-fraction 0.5 crossing. Sites
        absent from the coupling matrix (uncoupled, or frozen out of the
        subgraph in targeted mode) keep their single-site PB pKa — that is
        their honest effective value. Returns {(resname,resnum): float|str|None};
        returns the single-site map unchanged if coupling.pkl is absent.
        """
        import pickle as _pickle
        # Default: everyone at their single-site value.
        effective = dict(pb_pka)

        coupling_pkl = ws.get("pb_titrate_coupling")
        if not coupling_pkl or not Path(coupling_pkl).exists():
            self.console.print(
                "  [grey50]No coupling matrix on disk — effective pKa = single-site "
                "(run Coupling Precompute + Solve to get coupled values).[/grey50]")
            return effective

        try:
            with open(coupling_pkl, "rb") as f:
                data = _pickle.load(f)
            self_e = data["self_energies"]
            W = data["coupling"]
        except Exception as exc:
            self.console.print(
                f"  [yellow]Could not load coupling matrix ({exc}); effective "
                f"pKa = single-site.[/yellow]")
            return effective

        from .effective_pka import coupled_pka_titration
        n_sites = len(self_e.sites)
        self.console.print(
            f"  Computing effective (coupled) pKa: titrating {n_sites} coupled "
            f"site(s) over pH 2–12 (Monte Carlo, no PB calls)…")

        def _prog(done, total, pH):
            if done == total or done % 5 == 0:
                self.console.print(f"    [{done}/{total}] pH {pH:.1f}")

        pka, _curves = coupled_pka_titration(self_e, W, on_progress=_prog)
        # Overlay coupled values onto the single-site defaults.
        for key, val in pka.items():
            effective[key] = val
        return effective

    def _checklist_pbt_cmp_propka(self) -> Dict[str, Any]:
        from . import propka_compare as pc

        ws = self._ws()
        pb_pka = dict(ws.get("pb_titrate_pka") or {})
        if not pb_pka:
            pkl = self._output_dir() / "single_site_pka.pkl"
            if pkl.exists():
                with pkl.open("rb") as f:
                    res = pickle.load(f).get("results", {})
                pb_pka = {k: v.get("pKa_eff") for k, v in res.items()}
        if not pb_pka:
            raise RuntimeError(
                "No PB pKas found — run Single-Site pKa (step 4) first.")

        prmtop = Path(ws.get("pb_titrate_prmtop"))
        rst7 = Path(ws.get("pb_titrate_rst7"))

        self.console.print("\n[bold cyan]Compare PB pKas with ProPKA[/bold cyan]")
        self.console.print(
            f"  Re-running ProPKA on the same minimized structure PB used: "
            f"{prmtop.name}")

        exe = pc.find_propka()
        if exe is None:
            self.console.print(
                "  [yellow]propka3 not found on PATH — showing PB pKas only.\n"
                "  Install ProPKA to enable the side-by-side comparison.[/yellow]")
            self._print_pb_only(pb_pka)
            return {"summary": f"{len(pb_pka)} PB pKas (ProPKA unavailable)"}

        cmp_dir = self._output_dir() / "propka_compare"
        pdb_path = cmp_dir / "min_for_propka.pdb"
        try:
            pc.write_propka_pdb(prmtop, rst7, pdb_path)
            pka_path = pc.run_propka(pdb_path, cmp_dir)
        except (FileNotFoundError, RuntimeError) as exc:
            self.console.print(f"  [yellow]ProPKA step failed: {exc}[/yellow]")
            self.console.print("  [yellow]Showing PB pKas only.[/yellow]")
            self._print_pb_only(pb_pka)
            return {"summary": f"{len(pb_pka)} PB pKas (ProPKA failed)"}

        propka = pc.parse_propka_summary(pka_path)
        # The .pka file identifies LIGAND groups by residue name plus ATOM name
        # and never writes their residue number ("PRN  CG A  4.54  4.50  OCO"),
        # so the file parse silently drops every heme propionate even though
        # ProPKA does assess them. Recover those from ProPKA's Python API,
        # which keeps the owning atom. Additive only: keys already parsed from
        # the file win, so the well-tested standard-residue path is unchanged.
        try:
            api_groups = pc.parse_propka_groups(pdb_path)
        except Exception as exc:                     # noqa: BLE001
            api_groups = {}
            self.console.print(
                f"  [yellow]ProPKA API pass skipped ({type(exc).__name__}: {exc}); "
                f"ligand groups such as heme propionates will be missing.[/yellow]")
        recovered = {k: v for k, v in api_groups.items() if k not in propka}
        if recovered:
            propka.update(recovered)
            by_name = sorted({rn for rn, _ in recovered})
            self.console.print(
                f"  [grey50]{len(recovered)} ligand group(s) recovered via the "
                f"ProPKA API ({', '.join(by_name)}) - the .pka file records no "
                f"residue number for these.[/grey50]")

        # Effective (coupling-aware) pKa: titrate the coupled system over a pH
        # grid using the precomputed coupling matrix (no PB calls). Sites in the
        # coupling matrix get a coupled pKa; sites NOT in it (uncoupled / frozen
        # in targeted mode) keep their single-site value (that IS their honest
        # effective pKa). Skipped gracefully if coupling.pkl is absent.
        effective_pka = self._compute_effective_pka(ws, pb_pka)

        rows = pc.build_comparison(pb_pka, propka, effective_pka=effective_pka)

        csv_path = self._output_dir() / "pka_vs_propka.csv"
        with csv_path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["resname", "resnum", "pb_pka_single_site",
                        "pb_pka_effective", "propka_pka",
                        "delta_single_minus_propka",
                        "delta_effective_minus_propka"])
            for r in rows:
                eff = r.get("effective_pka")
                eff_str = (f"{eff:.4f}" if isinstance(eff, (int, float))
                           else ("" if eff is None else str(eff)))
                w.writerow([
                    r["resname"], r["resnum"],
                    "" if r["pb_pka"] is None else f"{r['pb_pka']:.4f}",
                    eff_str,
                    "" if r["propka_pka"] is None else f"{r['propka_pka']:.4f}",
                    "" if r["delta"] is None else f"{r['delta']:.4f}",
                    "" if r.get("delta_eff") is None else f"{r['delta_eff']:.4f}",
                ])
        ws.set("pb_titrate_pka_vs_propka_csv", str(csv_path))

        self._print_comparison_table(rows)
        n_matched = sum(1 for r in rows if r["propka_pka"] is not None)

        # Aggregate agreement statistics over residues with both predictions.
        # Compare ProPKA against the EFFECTIVE (coupling-aware) PB pKa when
        # available; summary_statistics falls back to single-site if no
        # effective value exists, and reports which basis it used.
        pH = self._params().get("pH", 7.0)
        stats = pc.summary_statistics(rows, pH, use="effective")
        if stats:
            r = stats["pearson_r"]
            r_str = f"{r:.3f}" if r is not None else "n/a"
            basis_lbl = ("EFFECTIVE / coupled PB pKa"
                         if stats.get("basis") == "effective"
                         else "single-site PB pKa")
            self.console.print(
                f"\n  [bold]Summary[/bold] (ProPKA vs {basis_lbl}; "
                f"n={stats['n']} with both values):")
            self.console.print(
                f"    mean Δ (PB−ProPKA)        {stats['mean_delta']:+.2f}"
                f"   [grey50](PB = {stats.get('basis','?')})[/grey50]")
            self.console.print(
                f"    median Δ                  {stats['median_delta']:+.2f}")
            self.console.print(
                f"    mean |Δ| (MAE)            {stats['mae']:.2f}")
            self.console.print(
                f"    RMSD                      {stats['rmsd']:.2f}")
            self.console.print(
                f"    Pearson r                 {r_str}")
            self.console.print(
                f"    within 2 pH units         "
                f"{stats['frac_within_2']*100:.0f}%")
            self.console.print(
                f"    protonation call @ pH {stats['pH']:.1f}    "
                f"{stats['binary_agreement']*100:.0f}% "
                f"({stats['binary_n_agree']}/{stats['n']})")
            # Explain the denominator: sites whose effective pKa is censored
            # (titration midpoint outside the pH window) have no numeric pKa and
            # are excluded from these statistics.
            n_locked = sum(1 for r in rows
                           if isinstance(r.get("effective_pka"), str))
            n_no_propka = sum(1 for r in rows if r["propka_pka"] is None)
            excl_bits = []
            if n_locked:
                excl_bits.append(f"{n_locked} with a censored (locked) "
                                 f"effective pKa")
            if n_no_propka:
                excl_bits.append(f"{n_no_propka} with no ProPKA value")
            if excl_bits:
                self.console.print(
                    f"    [grey50]n={stats['n']} of {len(rows)} sites; excluded "
                    f"{' and '.join(excl_bits)}.[/grey50]")

        self.console.print(f"\n  [green]CSV:[/green] {csv_path}")
        self.console.print(
            f"  [grey50]{n_matched}/{len(rows)} sites matched a ProPKA value "
            f"(unmatched = ProPKA omits non-titratable groups, e.g. "
            f"disulfide CYS).[/grey50]")
        if stats:
            return {"summary": (
                f"PB vs ProPKA: {n_matched}/{len(rows)} matched → "
                f"{csv_path.name}; RMSD={stats['rmsd']:.2f}, "
                f"call agree {stats['binary_agreement']*100:.0f}%")}
        return {"summary": (f"PB vs ProPKA: {n_matched}/{len(rows)} matched "
                            f"→ {csv_path.name}")}

    def _print_pb_only(self, pb_pka: Dict[Tuple[str, int], Optional[float]]) -> None:
        from rich.table import Table
        t = Table(title="PB single-site pKas")
        t.add_column("res")
        t.add_column("num", justify="right")
        t.add_column("PB pKa", justify="right")
        for (rn, num) in sorted(pb_pka, key=lambda k: (k[1], k[0])):
            v = pb_pka[(rn, num)]
            t.add_row(rn, str(num), "n/a" if v is None else f"{v:.2f}")
        self.console.print(t)

    def _print_comparison_table(self, rows: List[Dict[str, Any]]) -> None:
        from rich.table import Table
        pH = self._params().get("pH", 7.0)
        # Show the effective column only if any row carries an effective value.
        has_eff = any(r.get("effective_pka") is not None for r in rows)
        t = Table(title="PB vs ProPKA  (both on the minimized structure)")
        t.add_column("res")
        t.add_column("num", justify="right")
        t.add_column("PB pKa\n(single)", justify="right")
        if has_eff:
            t.add_column("PB pKa\n(coupled)", justify="right")
        t.add_column("ProPKA", justify="right")
        t.add_column("Δ (PB−ProPKA)", justify="right")

        def _protonated(pka, pH):
            """Protonation call: protonated when pKa > pH. Returns True/False for
            a numeric pKa; for a sentinel ('< x' below pH → deprotonated,
            '> x' above pH → protonated) resolve when the bound settles it,
            else None (undetermined)."""
            if isinstance(pka, (int, float)):
                return pka > pH
            if isinstance(pka, str):
                s = pka.strip()
                if s.startswith("<"):
                    try:
                        return False if float(s[1:]) <= pH else None
                    except ValueError:
                        return None
                if s.startswith(">"):
                    try:
                        return True if float(s[1:]) >= pH else None
                    except ValueError:
                        return None
            return None

        for r in rows:
            pb = "n/a" if r["pb_pka"] is None else f"{r['pb_pka']:.2f}"
            pp = "—" if r["propka_pka"] is None else f"{r['propka_pka']:.2f}"
            eff = r.get("effective_pka")
            eff_str = (f"{eff:.2f}" if isinstance(eff, (int, float))
                       else ("—" if eff is None else str(eff)))  # sentinel as-is

            # The PB pKa that drives the comparison: coupled when available
            # (float or sentinel), else single-site.
            pb_for_cmp = eff if (has_eff and eff is not None) else r["pb_pka"]

            # Δ column: only defined when BOTH the chosen PB pKa and ProPKA are
            # numeric. For a censored (sentinel) effective pKa, Δ is undefined —
            # show "—" rather than silently falling back to the single-site Δ
            # (which would mislead: a >12 site can AGREE with ProPKA while its
            # single-site Δ looks large).
            if isinstance(pb_for_cmp, (int, float)) and r["propka_pka"] is not None:
                d_val = pb_for_cmp - r["propka_pka"]
                d = f"{d_val:+.2f}"
            else:
                d_val = None
                d = "—"

            # Protonation-call disagreement at the system pH (the scientifically
            # important flag) → magenta. Large numeric Δ (≥3) but same call →
            # yellow (as before). Disagreement takes precedence.
            call_pb = _protonated(pb_for_cmp, pH)
            call_pp = _protonated(r["propka_pka"], pH)
            if (call_pb is not None and call_pp is not None
                    and call_pb != call_pp):
                style = "magenta"
            elif d_val is not None and abs(d_val) >= 3.0:
                style = "yellow"
            else:
                style = None

            if has_eff:
                t.add_row(r["resname"], str(r["resnum"]), pb, eff_str, pp, d,
                          style=style)
            else:
                t.add_row(r["resname"], str(r["resnum"]), pb, pp, d, style=style)
        self.console.print(t)
        if has_eff:
            self.console.print(
                "  [grey50]Δ uses the coupled PB pKa; '—' Δ = effective pKa is "
                "censored (titration midpoint outside the pH window). "
                "[magenta]magenta[/magenta] = PB and ProPKA disagree on "
                "protonation at the target pH; [yellow]yellow[/yellow] = same "
                "call but |Δ| ≥ 3.[/grey50]")


# ============================================================================
# Entry point
# ============================================================================

def run_pb_titrate_workflow(processor) -> bool:
    """Entry point invoked from `TLeapInputGenerator.run_pb_titrate`."""
    import warnings as _warnings
    from parmed.exceptions import AmberWarning as _AmberWarning

    workflow = PBTitrateWorkflow(processor)
    output_dir = workflow._output_dir()
    checklist = WorkflowChecklist(
        steps=PB_TITRATE_STEPS,
        executor=workflow,
        processor=processor,
        workflow_name="PB Titrate",
        console=processor.console,
        state_dir=output_dir,
        info_text=PB_TITRATE_OVERVIEW,
    )
    # Suppress parmed's "Molecule atoms are not contiguous! Attempting to
    # reorder to fix." warning for the duration of the workflow. parmed
    # emits this on every load and save when MOLECULE_PTR doesn't match
    # the bond-derived molecule layout (which happens whenever we slice a
    # cluster cutout that splits a heme's residues across non-contiguous
    # index ranges, and on every load of a prmtop that already had a
    # mismatch). The reorder is correct and there's nothing for the user
    # to act on — it's pure noise during step 3/4 (hundreds of save/load
    # cycles per run). Other AmberWarnings are still emitted normally.
    with _warnings.catch_warnings():
        _warnings.filterwarnings(
            "ignore",
            message=r".*Molecule atoms are not contiguous.*",
            category=_AmberWarning)
        return checklist.run()
