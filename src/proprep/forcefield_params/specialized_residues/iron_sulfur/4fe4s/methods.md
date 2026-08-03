# Fe4S4 cluster force-field methods

Authoritative build-method reference for the `iron_sulfur/4fe4s` parameter sets.
Companion to `README.md` (naming scheme) and the fuller narrative in
`proprep/docs/fe4s4_methodology.md` / `fe4s4_set_comparison.md`.

## 1. QM

Model built from PDB 1HIP. Two nested models per broken-symmetry (BS) spin-layer
assignment:

- **Small model** (SCH3-capped cluster) → geometry + Hessian for the Seminario
  bond/angle force constants.
- **Large model** (ACE-CYS-NME cap on each of the 4 ligating Cys, 98 atoms) →
  Merz-Kollman (MK) ESP for RESP charges.

Both at **UB3LYP/6-31G(d,p) + GD3 + C-PCM(diethyl ether, eps=4)**, Gaussian 16,
`int=ultrafine`. Large-model geometry optimizations freeze the protein backbone
(`opt=modredundant`) to hold the crystallographic scaffold; consequently the
large-model frequency step shows ~12-14 small imaginary modes (60-90 cm^-1) for
**every** assignment — an expected artifact of the constrained optimization, not
a defect.

Redox states / BS guesses (see README for the residue-code table):
- oxidized [Fe4S4]^3+ : 6 guesses (guess1a..guess3b) → FO1..FO6
- reduced HiPIP [Fe4S4]^2+ : 3 guesses → FR1..FR3
- reduced Fd [Fe4S4]^1+ : 6 guesses → FD1..FD6

## 2. Per-guess parameter sets (FD1..FD6, FO1..FO6, FR1..FR3)

Built with **MCPB.py** (Li & Merz):
- **Bonds/angles**: Seminario method from the small-model Hessian (`-s 2`).
- **Charges**: 2-stage RESP from the large-model MK ESP with ff14SB CYM backbone
  charges frozen (`-s 3c`). The 4 Fe are kept **distinct** (atom types X1-X4 /
  layer-specific), reflecting each Fe's ferric / mixed-valence role in that BS
  solution (from Mulliken spin populations).
- **Fe Lennard-Jones**: formal-state-matched TIP3P IOD (Li-Merz): Fe3+ Rmin/2
  1.386 / eps 0.01357; Fe2+ 1.409 / 0.01721; Fe2.5+ (mixed-valence) = the
  arithmetic mean, 1.3975 / 0.01539.

The 2-stage RESP for one guess is reproducible from its deposited inputs:
```
espf = <mklog basename>.esp                 # get_esp_from_gau(mklog)
resp -O -i resp1.in -o resp1.out -p resp1.pch -t resp1.chg -e $espf -s resp1_calc.esp
resp -O -i resp2.in -o resp2.out -p resp2.pch -q resp1.chg -t resp2.chg -e $espf -s resp2_calc.esp
```
(`resp2_calc.esp` is the `-s` OUTPUT, not an input.) Verified to reproduce the
deposited per-guess `resp2.chg` to max|diff| = 0.000000.

## 3. Averaged deployment (FD0 / FO0 / FR0)

A single transferable set per redox state with equivalenced cluster atom types
(FS/SB/SC), for studies driven by the integer redox-state change rather than a
specific Fe labeling. Two distinct averaging operations:

- **Bond/angle constants + composition-weighted Fe LJ**: arithmetic mean across
  the BS guesses (`average_parameters.py`). Fe LJ is composition-weighted between
  Fe2+/Fe3+ TIP3P IOD by the cluster's formal Fe3+:Fe2+ ratio; for FD (0.25 Fe3+
  / 0.75 Fe2+): Rmin/2 = 1.4033, eps = 0.01630. These are small-model / IOD
  quantities and are **not** affected by any large-model charge issue.

- **Charges**: a **joint multi-conformer group-equivalence RESP** over the
  concatenated MK ESP of all N guesses, with the 4 Fe, 4 bridging S, and 4 Cys SG
  each equivalenced to one charge across all conformers (so all 4*N instances of each
  share a single value), ff14SB backbone charge-constrained and ACE/NME caps neutral
  per conformer. This is a joint least-squares fit, NOT an arithmetic mean of the
  per-guess charges (which would differ by ~0.03-0.05 e on FS/SB — for FD, mean gives
  FS +0.737 vs the joint fit's +0.698). The fit uses the **standard 2-stage RESP
  protocol — identical to the per-guess FD1-6/FO1-6/FR1-3 sets** (stage 1 qwt=0.0005
  all free; stage 2 qwt=0.001, `iqopt=2`, heavy atoms frozen, methylene HB2=HB3
  equivalenced) — applied jointly with the cross-conformer core equivalencing above.
  See 3.1 for why 2-stage (not the single-stage originally deposited).

### 3.1 Reconstruction note (2026-07-06)

The original `multiconformer_resp.py`, `average_parameters.py`, and the per-guess
lib builder `regenerate_v3.py` are **not present on disk**; only their outputs and
the method description survived. The joint RESP driver was reconstructed
(`HiPIP/Parameterization/multiconformer_resp_recon.py`) from this method statement
and the exact `resp` multi-molecule input format (AmberTools `resp.F` subroutine
`mult_mol`): N molecule decks (`nmol=N`), each ending in a blank-line terminator,
with Fe2/3/4->Fe1 / S->S1 / SG->SG1 set by IVARY in the atom block and the
per-molecule Fe1/S1/SG1 linked by inter-molecular equivalencing cards. Per-molecule
charge-constraint atom references use molecule index 1 (local); `lagrange` supplies
the molecule offset via `ifirst`.

**What the original deposit actually did, and why it was replaced.** Reproducing the
deposited FD0 exactly showed it came from a **single-stage** group-equivalence joint
RESP at qwt ~= 0.001 with the CM hydrogens fit **freely** (no methylene equivalence).
That single stage reproduces every deposited FD0 atom to <=1.1e-4 e (FS +0.697527 vs
+0.697535, ... HB2 +0.061805, HB3 +0.068189), which confirmed the mechanism — but it
also revealed a **defect**: the deposited averaged sets carry **inequivalent HB2 !=
HB3** methylene charges (FD0 0.0619/0.0683; FO0 -0.013/-0.026; FR0 0.033/0.026). Two
protons on one carbon have no basis for different charges; the split is an artifact of
skipping the standard RESP stage 2, whose sole job is methylene equivalencing. It also
made the averaged sets the only members of their families *not* fit by the 2-stage
protocol used for every per-guess set.

**Corrected method (2026-07-06): standard 2-stage joint fit.** All three averaged sets
(FD0, FO0, FR0) were rebuilt with the per-guess 2-stage RESP protocol run jointly
(stage 1 qwt=0.0005 all-free with core equivalenced across conformers; stage 2
qwt=0.001 `iqopt=2` freezing heavy atoms and equivalencing methylene HB2=HB3). The
machinery is validated: with cross-conformer core equivalencing switched **off**, the
joint fit reproduces all six deposited per-guess `resp2.chg` to max|diff| = 1e-6, so
switching it on yields a faithful properly-averaged set. Effect vs the deposit: HB2 now
= HB3 exactly; the core is essentially unchanged by the protocol (the collapse fix, not
the protocol, moves FD0's Fe); side chains shift only slightly (FD0 CB +0.003 ->
-0.010). System charge is exact per state. The single-stage reconstruction driver
(`multiconformer_resp_recon.py`) is retained for provenance; the deployed fit uses the
2-stage joint generator.

## 4. guess3a (FD5) collapse fix + averaged-set 2-stage rebuild (2026-07-06)

**Problem.** The reduced-Fd `guess3a` **large-model** BS optimization had collapsed
to a non-BS delocalized state: SCF E = -10389.30 Ha (~+80 kcal/mol above its BS
siblings at -10389.434..-10389.440), <S**2> = 7.38 vs ~8.5, Fe spins
-3.90/+1.76/+1.94/+2.05. Its MK ESP and therefore its RESP charges were invalid
(mean Fe +0.648, outlier Fe2 +0.478). This poisoned **FD5** (= guess3a directly)
and **FD0** (guess3a is 1 of the 6 conformers in the joint fit). The small-model
Seminario bonded parameters were never affected. All other guesses (FO*, FR*,
FD1-4, FD6) are clean in both models.

**Re-convergence.** The collapse is at the SCF-convergence level, not the geometry.
guess3a was re-run by seeding the SCF from the guess3b converged orbitals with the
alpha/beta MO blocks swapped (a mirror-image BS state, `swap_ab.py`), then a
`CalcFC` constrained optimization. Result: the correct Fe2,Fe3-alpha / Fe1,Fe4-beta
layer state, Fe spins -3.75/+3.86/+3.86/-3.74, <S**2> = 8.54, E = -10389.4341
(degenerate with guess3b), stable wavefunction, freq in family (14 small imaginary
modes). (An independent level-shift re-run reached the same basin but a lopsided,
higher-energy solution; the spin-swap result was used.)

**FD5 re-fit.** MCPB `-s 3c` 2-stage RESP re-run on the new MK ESP. Corrected Fe
0.8743/0.6033/0.6820/0.7541 (mean +0.728); bridging S more negative; changes
localized to the core (RMS change over all 98 atoms 0.026; backbone untouched).

**FD0 re-derivation.** FD0 was rebuilt with the corrected 2-stage joint protocol (3.1)
and guess3a replaced by the corrected conformer (guess1a, 1b, 2a, 2b,
**guess3a_refit_B**, 3b). So FD0 carries **both** fixes at once: the guess3a collapse
correction and the single-stage -> 2-stage method switch.

| type | deposited (single-stage, collapsed g3a) | corrected (2-stage, fixed g3a) |
|------|-----------------------------------------|--------------------------------|
| FS (Fe)      | +0.697535 | **+0.723004** |
| SB (bridge S)| -0.796580 | **-0.817572** |
| SC (Cys SG)  | -0.685209 | **-0.686771** |
| CB           | +0.002784 | **-0.010345** |
| HB2          | +0.061872 | **+0.070192** |
| HB3          | +0.068296 | **+0.070192** (= HB2) |

The Fe lift (+0.698 -> +0.723) and more-negative bridging S (-0.797 -> -0.818) are the
proper ionic [Fe4S4] character the collapse washed out (the 2-stage protocol alone
moves the FD core <0.005 e; the collapse fix supplies the rest). The side chain barely
moves apart from HB2/HB3 merging to their mean. System charge (cluster + 4 Cys) =
-3.000000 exact.

**FO0 / FR0 re-derivation.** Neither oxidized nor reduced-HiPIP had a BS collapse, so
these were rebuilt with the 2-stage protocol switch **only** (original guesses). Cores
move <0.007 e (FO0 Fe +0.7078 -> +0.7146; FR0 Fe +0.6825 -> +0.6862); side chains are
within ~0.005 e of the deposit; HB2 now = HB3. System charges = -1.000004 (FO0, matching
the original deposit's rounding) and -2.000000 (FR0) exact.

**Validation.** All three corrected libs build under `tleap` (ff14SB + the unchanged
frcmod) with Errors = 0 (the 2 per-unit warnings are the expected non-integral cluster
/ CM fragment charges — they sum to the state charge only in combination). Backups:
`Guberman_FD{0,5}_RESP.lib.bak_collapsed_guess3a` (true pre-collapse-fix original) and
`Guberman_F{D,O,R}0_RESP.lib.bak_pre_2stage` (state immediately before the method
switch).

## 5. Provenance of inputs

Per-guess MCPB dirs and ESP grids: `HiPIP/Parameterization/reduced_fd_mcpb/guess*`
(+ `guess3a_refit_B` for the corrected guess3a). QM logs:
`HiPIP/Parameterization/reduced_fd_large_cpcm/` (originals) and `.../fd_g3a_B/`
(the corrected guess3a chain). RESP toolchain: AmberTools 25 `resp`.
