# Membrane systems and molecular docking (scope s8)

Everything below was driven through `proprep --full-menu` (ProPrep 1.23.0) on
2026-11-21 unless a scenario is marked SOURCE-ONLY.

Sessions used (transcripts are the primary evidence):

| session | cwd / project | transcript | ProPrep session log |
|---|---|---|---|
| `s8-mem1` | `runs/s8/mem1/mem2rh1` | `<repo>/manual/runs/_sessions/s8-mem1/transcript.txt` | `<repo>/manual/runs/s8/mem1/mem2rh1/proprep_session_20261121_123055.json` |
| `s8-mem2` | `runs/s8/mem2/kcsa` | `<repo>/manual/runs/_sessions/s8-mem2/transcript.txt` | `<repo>/manual/runs/s8/mem2/kcsa/proprep_session_20261121_124621.json` |
| `s8-mem3` | `runs/s8/mem3/solvonly` | `<repo>/manual/runs/_sessions/s8-mem3/transcript.txt` | `<repo>/manual/runs/s8/mem3/solvonly/proprep_session_20261121_133500.json` |
| `s8-dock1` | `runs/s8/dock1/dock1hsg` | `<repo>/manual/runs/_sessions/s8-dock1/transcript.txt` | `<repo>/manual/runs/s8/dock1/dock1hsg/proprep_session_20261121_124112.json` |

Environment facts observed in these runs (useful for the manual's "what ProPrep
calls" box): the external programs are taken from `$AMBERHOME/bin` first, and on
this machine `AMBERHOME=/opt/homebrew/Caskroom/miniforge/base/envs/etanalyze_env`,
so the screens show `packmol-memgen` and `teLeap` running out of **etanalyze_env**,
not the ProPrep env. `vina`, `vina_split` and `autogrid4` were found and used from
the ProPrep env, so all three scoring functions (vina / vinardo / ad4) work here.
MEMEMBED and PPM3 both ran (they ship with packmol-memgen).

---

## Feature: 4g. Membrane Builder

**What it is for (plain English):** it puts your protein into a lipid bilayer —
or builds a bare patch of membrane, or just a water box — by calling
`packmol-memgen`. It tilts/centres the protein so its membrane-spanning part sits
in the bilayer plane, packs lipids, water and ions around it, and writes one
packed PDB (`bilayer.pdb`) that the Topology Generator (5a) can parametrise.
It does **not** minimise or run MD; that is 5b.

**Needs before you start:**
- `packmol-memgen` on the path or in `$AMBERHOME/bin` (installed here).
- For a protein system: a structure loaded with 1a (or a local PDB chosen inside
  the tool). The protein should contain **only residues tLEaP can type** — see
  Scenario 4: ions/waters left in the file made the hydrogen-adding step fail.
- Nothing else; an empty bilayer or water-only box needs no structure at all.
- Expect minutes of packing: ~2 min for a 40x40 Å empty patch, ~8.5 min for a
  4-chain channel in a 70x70 Å POPC patch (this laptop).

**Where to find it:** Main menu `4g` (full menu) → `1. Build membrane-protein
system`. It then shows a single dashboard screen with nine numbered settings
groups plus `r` (review), `g` (generate & run), `b` (back). In workflow mode it
appears in "STRUCTURE PREPARATION" next to 4f Structure Orientator. No CLI flag;
it is menu-driven only.

The dashboard (exact layout):

```
  SYSTEM
    1   Protein Selection        none
    2   Lipid Composition        Required
    3   Solvent & Ions           TIP3P, 0.15 M K+/Cl-
  GEOMETRY
    4   Box & Membrane Dims      15.0 Å pad, 17.5 Å water
    5   Protein Orientation      Auto (memembed)
    6   Specialized Geometry     Curvature, CompEL...
  PACKING
    7   PACKMOL Settings         tol=2.0, nloop=20/100
  PARAMETRIZATION
    8   Force Fields             ff14SB / lipid21
    9   Custom Params & Solutes  Solutes: none
    r   Review Full Configuration
    g   Generate & Run
```

### Scenario 1: Empty (protein-free) bilayer, POPC:cholesterol 3:1 — RUN

- Starting point: fresh session `s8-mem1`, project name `mem2rh1`, nothing loaded.
- Steps:
  1. Prompt `Project name (.):` → typed `mem2rh1` (creates
     `runs/s8/mem1/mem2rh1/`, where every output lands).
  2. Prompt `Enter your choice (1a):` → typed `4g`.
  3. Prompt `Enter your choice (1):` → typed `1` ("Build membrane-protein system").
  4. Prompt `Enter choice (2):` → typed `1` (Protein Selection).
  5. Prompt `Select (3):` → typed `3` ("No protein (empty bilayer)").
     Screen: `No protein — building empty bilayer.` /
     `Empty bilayers need an explicit x/y patch size (packmol-memgen --distxy_fix).`
  6. Prompt `Membrane patch size (Å, x and y) (70.0):` → typed `abc` first, which
     answered `Please enter a valid number`, then `40` (small = fast).
  7. `2` Lipid Composition → `1` Build custom composition → `1` Symmetric →
     lipid `POPC`, ratio `3` → lipid `chl1` (lower case accepted, upper-cased
     internally), ratio `1` → `d` (done). Each entry is echoed with its full
     chemical name and formal charge, and a running table shows lipid/ratio/charge.
  8. `3` Solvent & Ions → `3` salt concentration → `-1` rejected
     (`Value must be at least 0.0`), then `0.1`; `4` ion types → cation `Na+`,
     anion Enter (keeps `Cl-`); `6` manual charge delta → `0`.
  9. `4` Box & Membrane Dims → `2` water layer → `12`.
  10. `7` PACKMOL Settings → `9` time limit → `0.25` (h) as an unattended-run guard.
  11. `8` Force Fields → `1` guided walkthrough: protein `1` (None — no protein),
      modified AA Enter (none), DNA/RNA/carbohydrate `1` (None), lipid `2`
      (lipid21), small molecule `1`, water `6` (TIP3P), divalent ions `1`.
  12. `r` to review, then `g` → `Proceed with build? [y/n] (y):` → `y`.
- What the user sees: the review panel ends with the literal command, which is
  the single most useful line for a reader who knows packmol-memgen:
  `packmol-memgen -p None -l POPC:CHL1 -r 3:1 --dist_wat 12.0 --distxy_fix 40.0 --salt --saltcon 0.1 --salt_override --salt_c Na+ --ffprot ff14SB --fflip lipid21 --ffwat tip3p -o bilayer.pdb --overwrite`.
  During packing ProPrep rewrites packmol's progress in plain English, e.g.
  `loop 92 of 100 | objective 3.54e+03 (lowest so far 895) | ... | 1.6 s per loop | at most 13 s more if every loop is needed`,
  with a one-off explanation that the objective jumping up is normal and that
  "lowest so far" is the number to watch. Warnings seen (from packmol-memgen):
  `Experimental value for POPC volume not available in parm file. Using estimated 1227 A^3 instead...` and
  `Sterols and/or PI are used, and the concentration is too high or the system too small to add proper side padding...`
  At the end: `packmol used all its loops and stopped short of a perfect packing... the close contacts left in it are for the minimization to relax`
  then
  ```
  Build successful!
    Output PDB: .../mem2rh1/bilayer.pdb
    Box: 40.0 × 40.0 × 70.0 Å
    Water molecules: 1,280
    Ions: 2 Na+, 2 Cl-
  Workspace updated. Next steps:
    1. Run Topology Generator to parametrize the system
    2. Run Molecular Dynamics Manager for minimization/equilibration
  ```
- Result / files produced (all in `<repo>/manual/runs/s8/mem1/mem2rh1/`):
  - `bilayer.pdb` — the packed system, 10,776 atoms (3,840 WAT atoms, lipids split
    into `PC`/`PA`/`OL` residues by `charmmlipid2amber.py`, 1,036 `CHL` atoms, 2 Na+, 2 Cl-).
    **It has no CRYST1 line** — the box is set later by tLEaP (`setBox mol vdw`).
  - `bilayer.pdb_FORCED` — packmol's last (imperfect) structure, kept by packmol-memgen.
  - `packmol.inp`, `packmol.log`, `packmol-memgen.log` — inputs/logs; the memgen log
    repeats the exact execution line and the box bounds (x -20..20, y -20..20, z -35..35).
- Time: packing 12:39:03 → 12:41 ≈ 2 min (20 loops per component, then 100 all-together).
- Problems: none fatal. The sterol warning appears because 25% cholesterol in a
  40 Å patch is dense; the manual should say a bigger patch avoids it.

### Scenario 2: Empty bilayer → 5a Topology Generator (how the files feed on) — RUN

- Starting point: straight after Scenario 1, same session (workspace now holds
  "membrane system" + the packed PDB + the force fields chosen in 4g).
- Steps: `m` → `5a` → `2` (Generate tLEaP input for single state) → template
  shown → `Edit the tLEaP template? [y/n] (y)` → Enter (opens **nano**; exited
  with Ctrl-X then `n` for "don't save") → `Enter output filename for tLEaP input
  (bilayer_tleap.in):` Enter → `View the tLEaP input file? [y/n] (y)` → `n` →
  `4` (Generate prmtop/rst7) → `Run bilayer_tleap.in? [y/n] (y)` → Enter.
- What the user sees: 5a does **not** re-ask anything about force fields or
  solvation for a membrane system —
  `Membrane system — using forcefields selected in membrane builder: lipids: lipid21 / water: TIP3P`,
  `Solvation skipped — system already solvated by membrane builder`,
  `Membrane system: vdW box + neutralization pass`, and for a protein-free system
  `Empty bilayer — no bond definitions needed.` The generated template is:
  ```
  source leaprc.lipid21
  source leaprc.water.tip3p
  mol = loadpdb .../bilayer.pdb
  setBox mol vdw
  addionsrand mol Na+ 0
  addionsrand mol Cl- 0
  check mol
  saveamberparm mol bilayer.prmtop bilayer.rst7
  ```
  tLEaP finished with `Errors = 0; Warnings = 132; Notes = 1` — 89 "Name change in
  pdb file residue N" warnings (lipid atom/residue renaming), one
  `addIonsRand: 1st Ion & target unit have charges of the same sign`, and ~40
  `Close contact of 0.8–1.5 angstroms between nonbonded atoms ...` (exactly the
  overlaps the message during packing said minimisation would fix). ProPrep then
  runs a ParmEd pass with an explanatory panel (validity check → rediscover
  molecules → optional 12-6-4) and reports `✓ Molecule definitions are correct`.
- Result / files produced: `bilayer_tleap.in`, `bilayer.prmtop` (3.5 MB),
  `bilayer.rst7` (box line `46.608 48.076 73.712 90 90 90` — note tLEaP's vdW box
  is bigger than packmol's 40x40x70 request), `bilayer_parmed.pdb`,
  `bilayer_parmed_replay.py` (re-runnable ParmEd script), `leap.log`.
- Time: ~40 s for tLEaP + ParmEd on 10.8k atoms.
- Problems: the 132 warnings scroll for several screens and look alarming to a
  beginner; they are normal for a packed bilayer.

### Scenario 3: `undo` inside the force-field walkthrough — RUN

- Starting point: mid-way through the Force Fields walkthrough (Scenario 1 step 11).
- Steps: at `Select rna forcefield (2):` typed `undo` → a table of the last
  interactions appeared (numbered 49-60 with module/question/answer) →
  `Rewind to which answer? (60)` → typed `60` →
  `[r]e-answer it live, or [c]hange it and replay what followed? ([x] cancel) (r):`
  (typing `y` here gives `Received 'y', which is not one of the available options`)
  → typed `r`.
- What the user sees: `[Rewinding: question 60 will be asked again; 1 answer from
  there on are set aside in a backup of the session log]`, then ProPrep **replays
  the whole session** from the start with `[REPLAY: ...]` markers against each
  previously-answered prompt, prints `[Replay complete - now recording new
  interactions]`, and stops at the DNA force-field question again.
- Result: a backup session log
  `proprep_session_20261121_123055.json.bak-20261121_123729`; the replay is fast
  (seconds) but the screen scroll is long.
- Problem worth flagging: the rewind silently discarded a setting made *before*
  the rewind point in a different menu — the PACKMOL time limit (0.25 h) was back
  to `none` in a later review.

### Scenario 4: Membrane protein 1BL8 (KCSA) as deposited — FAILED (and how to get past it)

- Starting point: session `s8-mem2`, project `kcsa`, 1BL8 downloaded with
  `1a → 1 → 1 → 1 → 1BL8 → 1` (asymmetric unit). On load ProPrep warns
  `Non-standard residue numbering detected ... Chain A: First residue is ALA #23`
  for all four chains.
- Steps: `4g → 1`; `1` Protein Selection shows `Current: filtered.../1BL8.pdb`
  already picked from the workspace (option 1 = "Use workspace structure");
  `5` Orientation → `2` (PPM3) with `6 Compare with the OPM database: Yes`;
  `2` Lipid Composition → `3` raw string → typed `POPC:FOO` →
  `Unknown lipid: 'FOO'` (rejected, nothing set) → `2 → 1` custom, `2` Asymmetric
  → upper `POPC` ratio 1 → `d`, lower `POPE` ratio 1 → `d` (shows `POPC//POPE`);
  `4` box → padding `10`, water `12`; `8` Force Fields walkthrough with protein
  `3` (ff14SB), lipid lipid21, water TIP3P; `g` → `y`.
- What the user sees, in order:
  ```
  tLEaP hydrogen addition
  Adding hydrogens via tLEaP using AMBER library templates
  tLEaP script: .../pretleap_hydrogen.in
  tLEaP did not produce topology files
  Full log: .../leap.log
  tLEaP could not add the hydrogens. packmol-memgen can add them itself
  (pdb2pqr/reduce), but it then removes residues it does not recognise
  (cofactors such as hemes, metal sites) and assigns protonation states of its
  own, replacing the ones in this structure.
  Build anyway with packmol-memgen's own protonation? [y/n] (n):
  ```
  Answering `y` continued: PPM3 ran
  (`PPM3: ppm3in.pdb emin= -79.9 thickn= 33.6+- 1.5 tilt= 0.+- 1.` /
  `Placed from 2824 atoms PPM3 kept (fit RMSD 0.0005 Å); all 2824 atoms moved as one body`),
  then the OPM cross-check, which is a genuinely nice screen:
  ```
  Asking the OPM database for 1BL8...
  OPM has 1BL8. Against its placement (2824 atoms paired with OPM's file):
    same way up as in OPM
    membrane normal: 0.1° away from OPM's
    protein centre along the normal: z = +5.2 Å here, +4.9 Å in OPM (difference +0.3 Å)
    OPM's hydrophobic thickness: 32.4 Å
  ```
  Packing then ran ~5 min and finished, but ProPrep **rejected the result**:
  ```
  Removed 4 spurious TER record(s) from output PDB
  The protein in bilayer.pdb is not the one packmol-memgen was given: atom 6 is
  ALA HA in residue 1 of the file, where the protein given has LEU N in residue 2.
  The Topology Generator's bond directives would not apply to it.
  The build is not recorded.
  ```
- Root cause (from `leap.log`): the pre-tLEaP hydrogen script loads the raw PDB,
  which still has 3 K+ ions and 1 crystal water in the selectivity filter; tLEaP
  logs `Bond: maximum coordination exceeded on .R<K 411>.A<K 1>` and
  `ATOMS NOT BONDED: .R<K 411>.A<K 1> .R<WAT 414>.A<O 1>` and then stops without
  writing a prmtop. packmol-memgen's own protonation (`pdb2pqr`) then trimmed and
  renumbered residues, so ProPrep's safety check (packed protein must match the
  protein it handed over) failed and the build was discarded.
- Files left behind: `bilayer.pdb` + `bilayer.pdb_FORCED` exist on disk but are
  **not** registered in the workspace, `protein_oriented.pdb`,
  `_ppm3_orientation/` (ppm3in.pdb, ppm3inout.pdb, ppm.inp, ppm3.log, res.lib,
  datapar*, fort.*), `opm_1bl8.pdb` (downloaded OPM reference), `leap.log`.
- Time: ~6.5 min wasted (orientation ~1 min + packing ~5 min) before the rejection.
- How to get past it: strip the non-protein residues first (Scenario 5).

### Scenario 5: Same protein after 4a PDB Filter — RUN (full success, end to end)

- Starting point: same session, 1BL8 loaded, 4g has just refused to record a build.
- Steps:
  1. `m → 4a → 1` (Filter PDB structure). It first prints a buried-surface-area
     matrix and a "Chain Topology" panel (4 chains, cyclic), asks
     `Detect redox-active sites before filtering? [y/n] (y)` → `n`,
     `Would you like to view the structure in the interactive 3D viewer? [y/n] (n)` → Enter,
     a residue-composition table, `Do you want to proceed with filtering? [y/n] (y)` → Enter.
  2. One prompt per **component**, in file order:
     `Chain A - Protein` → `r`, `Chain A - POTASSIUM ION` → `d`,
     `Chain B - Protein` → `r`, `Chain B - Water` → `d`,
     `Chain C - Protein` → `r`, `Chain D - Protein` → `r`.
     (Each offers `[r] Retain entire component / [s] Select specific residues /
     [d] Discard entire component`.)
  3. Summary panel lists only `Amino_acid: 97 residues` per chain →
     `Do you want to apply these filters? [y/n] (y)` → `y` →
     `Filtered structure automatically saved to: filtered_structure.pdb`.
  4. `m → 4g → 1`. **The membrane configuration had been reset** by the new
     structure: `Protein: filtered_structure.pdb`, `Lipids: not set`,
     `Orientation: Auto (memembed)`, `Box padding: 15.0 Å`.
  5. `2 → 3` raw lipid string → `POPC` → `Ratio string (or Enter to skip) ():` → `4`
     (so POPC only; the ratio for a single lipid is cosmetic).
  6. `4` box → padding `8`, water layer `10`.
  7. `5` Orientation → `1` (MEMEMBED) — `Set: Automatic orientation via MEMEMBED`.
  8. `g` → review shows
     `packmol-memgen -p filtered_structure.pdb --notprotonate --nottrim -l POPC -r 4 --dist 8.0 --dist_wat 10.0 --salt --salt_override --ffprot ff14SB --fflip lipid21 --ffwat tip3p -o bilayer.pdb --overwrite`
     → `y`.
- What the user sees: this time `Hydrogens added via tLEaP → protein_with_h.pdb`,
  then `Orienting protein (MEMEMBED...)`, then
  `Preprocessing protein... Protein charge: 4 / Protein mass: 39710 Da /
  Estimated protein volume: 34206 A^3`, nine packing components, and
  ```
  Build successful!
    Output PDB: bilayer.pdb
    Box: 69.9 × 69.9 × 78.5 Å
    Water molecules: 5,142
    Ions: 9 K+, 13 Cl-
  ```
  Note the `--notprotonate --nottrim` flags: because ProPrep added the hydrogens
  itself, packmol-memgen is told to leave the protein alone — which is why the
  verification passes this time.
- Then 5a: `m → 5a → 2`: template now sources `leaprc.protein.ff14SB` +
  `leaprc.lipid21` + `leaprc.water.tip3p`, loads `bilayer.pdb`, `setBox mol vdw`,
  `addionsrand mol K+ 0`, `addionsrand mol Cl- 0`. Because no redox sites were
  detected it asks `No bond definitions found. Run comprehensive RedoxSite
  detection first.` → `Do you want to continue without bond definitions? [y/n] (n)`
  → `y` → `Continuing with empty bond definitions.` → filename Enter → `n` to
  viewing → `4` → `1. Generate missing only`.
  tLEaP: `Errors = 0; Warnings = 820; Notes = 1` (581 close-contact lines).
  ParmEd flagged one real chemistry question:
  `Detected two cysteine residues whose sulfur atoms are within 3 Angstroms. Rename CYS to CYX ... and use the 'bond' command in tleap`.
- Result / files produced (`<repo>/manual/runs/s8/mem2/kcsa/`):
  `filtered_structure.pdb`, `pretleap_hydrogen.in/.prmtop/.rst7`,
  `protein_with_h.pdb`, `protein_with_hin_memembed.log` (MEMEMBED's own log),
  `bilayer.pdb` (37,282 atoms: 15,426 WAT, POPC as PC/PA/OL, K+/Cl-),
  `bilayer.pdb_FORCED`, `bilayer_tleap.in`, `bilayer.prmtop` (11.6 MB),
  `bilayer.rst7`, `bilayer_parmed.pdb` (`CRYST1 78.980 78.390 83.691`),
  `bilayer_parmed_replay.py`, `packmol*.log`, `leap.log`.
- Time: hydrogens + MEMEMBED ~1 min, packing 13:12:45 → 13:21:16 (8.5 min, 4.1-5.5 s
  per all-together loop), tLEaP + ParmEd ~1.5 min. ~11 min total.
- Problems: (a) PDB Filter's residue-selection prompt crashed on a range —
  typing `1-97` at `Select residues (comma-separated indices, 'all', or 'none')`
  produced `Error executing PDB Filter.filter_pdb_structure: invalid literal for
  int() with base 10: '1-97'` and dumped the user back to the filter menu (losing
  the selections made so far). Use comma-separated numbers or `all`.
  (b) After filtering, 4g silently forgot the lipid/box/orientation choices.

### Scenario 6: Water-only box (no lipids) around a soluble protein — RUN

- Starting point: session `s8-mem3`, project `solvonly`, `1UBQ` loaded via 1a.
- Steps: `4g → 1` → `2` Lipid Composition → `4` (No lipids, solvate only). Panel
  explains it builds "the same conventions as a membrane build ... Useful for
  preparing a matched non-membrane reference of the same system", then asks
  `Build a solvent-only system (no lipids)? [y/n] (y)` → `y` and
  `Use a cubic box? [y/n] (n)` → `n`. Then `4 → 1` padding `10`, `g`, `y`.
- What the user sees: review shows `MEMBRANE / Solvate only: yes (no lipids, water
  + ions only)` and the command
  `packmol-memgen -p protein_with_h.pdb --notprotonate --nottrim --dist 10.0 --solvate --salt --salt_override ...`.
  Build output: `Box: 59.5 × 59.5 × 66.0 Å`, `Water molecules: 7,884`,
  `Ions: 21 K+, 21 Cl-`.
- Result / files: `<repo>/manual/runs/s8/mem3/solvonly/bilayer.pdb`
  (+ `protein_with_h.pdb`, `pretleap_hydrogen.*`, `packmol*`). The output is still
  called `bilayer.pdb` even with no lipids — confusing, worth a note.
- Time: ~1 min.
- Problems: the review still prints `FORCE FIELDS ... Lipid: lipid21` and
  `leaprc needed: leaprc.lipid21` for a system with no lipids.

### Scenario 7: Every remaining settings screen exercised (values + resulting flags) — RUN

Driven in `s8-mem1`/`s8-mem3` purely to read the screens and confirm which
`packmol-memgen` flag each choice becomes (the `r` review prints the command, so
this is the cheapest way to document them):

- Specialized Geometry (`6`): curvature `0.0`; Gaussian deformation C `10` D `20`
  H `0`; channel plug radius `8`; head plane `18` / tail plane `4`; APL offset
  `1.1` / lipid segment multiplier `1.0`; self-assembly toggled on
  (`Self-assembly mode: enabled`); double bilayer `y` → `Charge imbalance
  (electron charge units) (0):` `20` → `Imbalance ion: [1] Cation [2] Anion` `1`;
  PBC toggled on (`PBC: enabled`). Resulting command line (reviewed, not built):
  `packmol-memgen -p None --apl_offset 1.1 --solvate --cubic --pbc --xygauss 10.0 20.0 0.0 --channel_plug 8.0 --headplane 18.0 --tailplane 4.0 --self_assembly --double --charge_imbalance 20 ...`
- PACKMOL Settings (`7`): `7` Troubleshooting asks three yes/no questions —
  `Randomize badly-placed molecules?` (`--movebadrandom`),
  `Add short-range overlap penalty?` (`--short_penalty`),
  `Skip individual packing steps (packall)?`; `8` Output options asks
  `Save intermediate PDB snapshots?` (`--traj`),
  `Create optimization function plot?` (`--plot`), `Use random seed?`;
  `9` time limit prints "the build shows packmol's progress loop by loop and
  Ctrl-C stops it cleanly, so a limit is only a safeguard for an unattended run".
- Orientation (`5`): `4` MEMEMBED options sub-menu (algorithm GA / Grid / Direct /
  GA×5, beta-barrel mode, keep ligands after alignment); `5` N-terminus
  orientation → typed `in,out` → `--n_ter in --n_ter out`; `3` Pre-oriented →
  `Set: Pre-oriented (alignment skipped)` → `--preoriented`;
  `6` toggles the OPM cross-check.
- Solvent & Ions (`1` water model): list is `Auto (currently TIP3P from ff14SB)`,
  TIP3P, OPC, OPC3, SPCE, SPCEB, TIP4PEW, TIP4PD, FB3; choosing OPC showed
  `Water model: OPC` with the "(auto: ...)" suffix gone.
- Lipid library browser (`2 → 2`): `packmol-memgen includes 4055 lipid types.`
  Browse by category lists 10 categories with charges and counts (PC 73, PE 73,
  PS 50 charge -1, PG 64, PA 50, PI 1323 charges -7..-1, SM 7, Sterols 5,
  Cardiolipins 2401, SIRAH coarse-grain 9); each category prints a table of
  name / full name / APL (Å²) / charge. Search for `cholesterol` → `CHL1, APL 19`;
  search for `zzzq` → `No lipids matching 'zzzq'`.
- Custom Params & Solutes (`9`): `1 Add solute` → `Solute PDB file path:` →
  typed `nofile.pdb` → `File not found: nofile.pdb` (returns to the sub-menu).
  Remaining prompts (SOURCE-ONLY, read from
  `<repo>/src/proprep/membrane_prep/membrane_builder.py`
  `_add_solute`): formal charge, `Count (integer) or concentration (e.g. '0.05M')`,
  `Place in membrane (vs. in solution)?`, optional minimum cylindrical distance
  from the protein, then `1 Provide frcmod + lib files` / `2 Standard residue`
  with the warning that tLEaP cannot build a topology for a non-standard solute
  without parameters.
- Protein Selection (`1 → 2` "Load a different PDB file"): opens a file browser
  (`N` to pick, `..` up, a typed path, `q` cancel); picking
  `protein_with_h.pdb` set `Protein: protein_with_h.pdb`.

### Options reference — 4g Membrane Builder

| Where | Choice | What it does |
|---|---|---|
| Main dashboard | `1` Protein Selection | pick the solute |
| | `2` Lipid Composition | what the bilayer is made of (required unless solvate-only) |
| | `3` Solvent & Ions | water model, salt, counterions, charge fixes |
| | `4` Box & Membrane Dims | padding, water layer, leaflet width, fixed sizes, lipid counting |
| | `5` Protein Orientation | MEMEMBED / PPM3 / pre-oriented, OPM cross-check |
| | `6` Specialized Geometry | curvature, plug, CompEL, self-assembly, PBC |
| | `7` PACKMOL Settings | tolerance, loops, troubleshooting, outputs, time limit |
| | `8` Force Fields | one guided walkthrough that also fixes 5a's choices |
| | `9` Custom Params & Solutes | extra molecules in the box + their frcmod/lib |
| | `r` / `g` / `b` | review (prints the equivalent command) / build / back |
| 1 Protein Selection | `1` Use workspace structure | the structure loaded in 1a |
| | `2` Load a different PDB file | file browser |
| | `3` No protein (empty bilayer) | then asks the x/y patch size → `--distxy_fix` |
| 2 Lipid Composition | `1` Build custom composition | symmetric or asymmetric; lipid+ratio per leaflet, validated against the library |
| | `2` Search / browse lipid library | search by name/keyword or browse 10 categories |
| | `3` Enter raw lipid string | `LIPID1:LIPID2`, `UPPER//LOWER`, `UPPER///LOWER` (double bilayer); unknown names rejected |
| | `4` No lipids (solvate only) | water+ions only → `--solvate` (+ optional `--cubic`) |
| 3 Solvent & Ions | `1` water model | auto / tip3p / opc / opc3 / spce / spceb / tip4pew / tip4pd / fb3 → `--ffwat` |
| | `2` toggle salt | `--salt` on/off |
| | `3` salt concentration | 0.0-5.0 M → `--saltcon` |
| | `4` ion types | cation / anion strings → `--salt_c` / `--salt_a` |
| | `5` toggle counterions | neutralisation on/off |
| | `6` manual charge adjustment | integer added on top of the auto-detected redox correction |
| 4 Box & Dims | `1` box padding | XY padding → `--dist` |
| | `2` water layer | Z water thickness → `--dist_wat` |
| | `3` leaflet width | per-leaflet thickness (default 23 Å) |
| | `4` fix XY | `--distxy_fix` (required for an empty bilayer) |
| | `5` fix X Y Z | fully fixed box |
| | `6` lipid counting | APL-based or volume-based |
| | `7` clear fixed dims | back to automatic |
| 5 Orientation | `1` MEMEMBED | default automatic orientation |
| | `2` PPM3 | the OPM method; ProPrep moves the whole structure into the PPM3 frame (`protein_oriented.pdb`) |
| | `3` Pre-oriented | `--preoriented`, use for an OPM download or a 4f-oriented file |
| | `4` MEMEMBED options | algorithm (GA/Grid/Direct/GA×5), beta-barrel mode, keep ligands |
| | `5` N-terminus orientation | per-chain `in`/`out` → `--n_ter` |
| | `6` OPM cross-check | compares the placement with OPM's (internet; reports only) |
| 6 Specialized | `1`-`8` | curvature, `--xygauss`, `--channel_plug`, `--headplane`/`--tailplane`, `--apl_offset`/lipid offset, `--self_assembly`, `--double` + `--charge_imbalance` (CompEL, Kutzner 2011), `--pbc` |
| 7 PACKMOL | `1`-`9` | clash tolerance (2.0 Å), protein atom radius (1.5), nloop (20), nloop_all (100), GENCAN iterations (20), move fraction (0.05), troubleshooting trio, output trio, time limit in hours |
| 8 Force Fields | `1` guided walkthrough | protein (10 options incl. constant pH/redox ff10), modified AA, DNA, RNA, carbohydrate, lipid (lipid21/lipid17), small molecule (GAFF2/GAFF), water (11), divalent ion sets (13). Choices are reused by 5a without re-asking |
| 9 Solutes | `1` add / `2` remove | PDB path, formal charge, count or molar concentration, in-membrane vs in-solution, minimum distance from protein, frcmod+lib or "standard residue" |

---

## Feature: 5d. Molecular Docking

**What it is for (plain English):** it predicts how a small molecule sits in your
protein's pocket, and how strongly, using AutoDock Vina. ProPrep prepares both
sides for you (Meeko/pdb2pqr — no force-field parameters needed), lets you define
the search box from the crystal ligand, from residues, or by typing numbers, runs
the docking, and reports the poses with their scores and (when the ligand came
from the structure) the RMSD back to the crystal pose. It can also run one
receptor against a whole ligand library ("campaign").

**Needs before you start:** a structure loaded with 1a (the menu is greyed out
with `⚠ Needs a loaded structure` otherwise — selecting it just re-prints the
warning); `vina` on the path; `autogrid4` as well if you want AD4 scoring (both
present here). For a library campaign, load the SDF/mol2/SMILES file first with
`1a → 4`. No AMBER parameters, no tLEaP, no MD.

**Where to find it:** Main menu `5d` → `1. Dock a ligand into the structure`,
then a single dashboard of nine numbered setup steps plus `r` run, `c` campaign,
`v` view results, `b` back. The dashboard shows what is still missing on the `r`
line, e.g. `r Run docking (needs a receptor structure (1), a ligand (2), a search
box (8), docking settings (9))`.

### Scenario 8: 1HSG redock of MK1 — receptor preparation path — RUN

- Starting point: session `s8-dock1`, project `dock1hsg`; `5d` first tried with
  nothing loaded → `⚠ Needs a loaded structure` and the menu re-printed.
  Then 1HSG downloaded (`1a → 1 → 1 → 1 → 1HSG → 1`; the metadata screen shows
  title, method, 2.0 Å resolution and the JBC reference).
- Steps:
  1. `5d → 1`, then `1` Receptor. Screen: `Current structure: .../1HSG.pdb` →
     `Use this structure as the receptor? [y/n] (y)` → `y` →
     `Chains: A (99 residues), B (99 residues)` →
     `Chains to include (e.g. A or A,B; Enter for all) (all):` → typed `C` →
     `No chain C.` (re-asked) → Enter (= all).
  2. `4` Protonation → `pH for protonation states (range: 0.0-14.0, default: 7.0)`
     → Enter → `Running pdb2pqr with PROPKA at pH 7.0 ...` →
     `pdb2pqr assigned 34 titratable residues.` Then the clever summary:
     "Of the 34 titratable residues, 31 are in their usual charged form at pH 7.0
     (Asp-, Glu-, Lys+, ...) and 3 are not, or are His ... Only those 3 are
     listed ('a' lists all)". Table: `A:69 HIS HID`, `B:25 ASP ASH`,
     `B:69 HIS HID`, each with `from = pdb2pqr/PROPKA at pH 7.0`. A second table
     lists the chain ends with their justification
     (`A:1 PRO first N-terminus: charged NH3+`, `from SEQRES: all 99 residues
     modelled, none missing`).
     At `Residue or chain end to change ...` I typed `A:69` → `Template for A:69
     HIS [HID/HIE/HIP]` → `HIE` (table now shows `from = chosen`); `a` → every
     titratable residue listed; `Z:1` → `No residue Z:1.`; `A:1` →
     `A:1 PRO, the first residue of chain A: [n] N-terminus/[b] chain break (n):`;
     Enter to finish.
  3. `5` Metals/cofactors → `No metals are kept in the receptor; no charges to set.`
     `No cofactors are kept (only standard residues and ions); no chemistry to review.`
  4. `6` Flexible → `Residues whose side chains move during docking (e.g.
     A:45,A:88; 'none' for a rigid receptor; 'keep' to leave as is; 'p' to pick in
     the viewer) (keep):` → typed `A:50` (dashboard then showed `6 Flexible: A:50`),
     later set back with `none`.
- Result: no files yet; the setup is saved in the session state (and later copied
  into every run's `docking.json`).
- Problems: nothing blocking here, but see Scenario 10 — the chain-end decision
  made in step 2 is what breaks or fixes the run.

### Scenario 9: Choosing the ligand from the structure, and the rotatable-bond review — RUN

- Steps:
  1. `2` Ligand → `Ligand source: [s] residue of the structure, [m] SMILES,
     [f] file (s):` → `s` → table `1 | B:902 | MK1 | 45 atoms` →
     `Ligand number (range: 1-1, default: 1)` → `1`.
  2. Screen: the CCD name of MK1 (indinavir), then
     `Ligand: 45 heavy atoms, 47 hydrogens, net charge +0; charged atoms: none`,
     `SMILES: CC(C)(C)NC(=O)C1CN(Cc2cccnc2)CCN1CC(O)CC(Cc1ccccc1)C(=O)NC1c2ccccc2CC1O`,
     and provenance lines (`45 heavy atoms from MK1 B902 in 1HSG.pdb`,
     `Bond orders, formal charges and hydrogens from CCD MK1 (net charge +0)`).
     The edit prompt (`remove NAME | charge NAME VALUE | bond NAME1 NAME2 ORDER |
     add-h NAME; 'p' to pick in the viewer; 'drop' ...`) was left empty.
     Picking the ligand this way **takes it out of the receptor** automatically
     (`3 Kept HETATM: no other HETATM residues in the file; MK1 B:902 is the ligand (2)`).
  3. `7` Rotatable → table of 16 bonds with `kind` (single / amide / conjugated
     single) and `now` (rotates / rigid) and `14 rotatable bonds.` Typing `3`
     flipped the amide `C3-N2` to rotating (15), typing `3` again flipped it back;
     `c` ("hold every conjugated single bond rigid") changed nothing for MK1.
     Enter → `7 Rotatable: reviewed`.
- Result / files: `docking/ligand_preview.sdf` and `docking/ccd_preview_MK1.sdf`
  in the project folder.

### Scenario 10: First real run — FAILED on the C-terminal residue, then fixed — RUN (both halves)

- Starting point: receptor + ligand + box (`8 → l` around the ligand, padding 8 →
  `Box centre (11.63, 22.47, 5.86), edges 30.9 x 26.2 x 24.6 A (19925 A^3)`) +
  settings (`9` → scoring `ad4`, everything else default: seed 42,
  exhaustiveness 8, 9 poses, energy range 3.0, min RMSD 1.0, spacing 0.375, CPU 0).
- First attempt: `r` → `Building the ligand and receptor and docking; this takes
  seconds to minutes ...` then a Meeko error:
  ```
  No template matched for residue_key='A:99' ... CPHE heavy_miss=1 ...
  Run docking failed:
  Error: Creation of data structure for receptor failed.
  Details:
  - Template matching failed for: ['A:99', 'B:99']
  These residues can be ignored with option --delete_bad_res ...
  ```
  (a Python traceback is printed above the message). Cause: 1HSG contains **no
  OXT atoms**, so the charged C-terminus template `CPHE` is missing one heavy
  atom — but ProPrep had proposed "C-terminus: charged COO-" because SEQRES says
  all 99 residues are modelled.
- Fix that worked, entirely inside the app: `4` Protonation → Enter (pH 7.0) →
  typed `A:99` → `A:99 PHE, the last residue of chain A: [c] C-terminus/[b] chain
  break (c):` → `b`; same for `B:99`; Enter.
- Second attempt: `r` succeeded.
  ```
  Done in 59.0 s. Output: .../docking/docking_20261121_130216
  - AutoGrid maps: 82x68x64 points at 0.375 A, types C NA A OA N HD; distance-dependent dielectric.
  - Vina: WARNING: At low exhaustiveness, it may be impossible to utilize all CPUs.
  - 8 poses within 3.0 kcal/mol of the best, ad4 scoring, seed 42, exhaustiveness 8.
  ```
  Pose table (total / inter / intra / torsions / RMSD to crystal): best
  `-12.88 … 4.05 Å`, and pose 5 `-11.24 … 1.78 Å` — i.e. **AD4's top-ranked pose
  was not the crystal pose** for this ligand/box.
- Result / files (`runs/s8/dock1/dock1hsg/docking/docking_20261121_130216/`):
  `docking.json` (the complete setup, including every protonation assignment),
  `ligand.pdbqt`, `receptor_rigid.pdbqt`, `poses.sdf` + `pose_1..8.sdf`,
  `vina_output.txt`, and the AutoGrid products `receptor.{A,C,HD,N,NA,OA,d,e}.map`,
  `receptor.gpf`, `receptor.glg`, `receptor.maps.fld`, `receptor.maps.xyz`.
- Problem to flag for the manual: "C-terminus vs chain break" is a question a
  biologist cannot answer without knowing about OXT atoms, and the error message
  (Meeko's) does not point at the fix. ProPrep's own default is the one that fails
  on this very common kind of file.

### Scenario 11: Vinardo and Vina redocks (small, fast settings) — RUN

- Steps: `9` → scoring `vinardo`, seed `7`, exhaustiveness `4`, poses `5`,
  then Enter for energy range / min RMSD / spacing / CPU → `r`. Repeated with
  scoring `vina`.
- Results:
  - vinardo: `Done in 29.1 s`, best pose `-9.07 kcal/mol`, **RMSD to crystal 0.41 Å**
    (folder `docking_20261121_131537`).
  - vina: `Done in 17.4 s`, best pose `-11.56 kcal/mol`, **RMSD 0.46 Å**
    (folder `docking_20261121_131645`).
  - For reference, `vina_output.txt` holds Vina's own table
    (`mode | affinity | dist from best mode`).
- Reading: both empirical functions put the top pose essentially on the crystal
  pose; AD4 (which needs maps and took 3.5x longer) did not. Scores are **not**
  comparable between functions — the menu says so before asking.

### Scenario 12: Flexible side chains — RUN

- Steps: ligand back to `MK1 B:902`, `6` Flexible → typed `A:25,A:50`
  (catalytic Asp25 + Ile50 of the flap), box `l` padding 8, `7` Enter, `r`.
- Result: `Done in 26.9 s`, best pose `-11.43`, RMSD 0.50 Å
  (`docking_20261121_133128`). Extra files compared with a rigid run:
  `receptor_flex.pdbqt` and `flexible_residues_pose_1..5.pdb` (the side chains as
  docked, one file per pose).

### Scenario 13: Ligand from SMILES, and the box trap — RUN

- Steps: `2` Ligand → `m` → `SMILES:` → `CC(=O)Nc1ccc(O)cc1` (paracetamol) →
  `Seed for the 3D build (RDKit ETKDG) (42):` Enter →
  `Minimise the 3D build with MMFF94?` Enter.
- What the user sees: `Ligand: 11 heavy atoms, 9 hydrogens, net charge +0`,
  `3D coordinates from RDKit ETKDGv3, seed 42.`,
  `Atom names generated for all 20 atoms (element + number).`,
  `Geometry minimised with MMFF94 (converged).`
  Chemistry edits tried at the edit prompt: `remove ZZ9` →
  `'remove ZZ9: no atom named ZZ9.'`; `remove H9` →
  `1 chemistry edit(s) applied: remove H9`; `drop` → `Dropped: remove H9`.
- Switching the ligand **cleared** `3 Kept HETATM` (MK1 is a HETATM again) and
  `8 Search box` (it had been defined from the old ligand).
- `3` Kept HETATM screen: table of HETATM residues with `kind`/`receptor` columns
  and a water line: `Waters: 127 in the file, 0 kept. ... 'w' keeps or leaves out
  all 127 at once`. Typing `1` toggled MK1 keep/leave-out; `w` kept all 127 waters
  (listing the first eight ids), `w` again dropped them; `B:312` kept that single
  water, typing it again dropped it. I left MK1 **kept** for this run, which makes
  it part of the rigid receptor — and `5` then offered its chemistry for review:
  `MK1 ...: 92 atoms, net formal charge +0 ... The entry is consistent. Press
  Enter to use it as written`.
- `7` Rotatable for paracetamol showed only 3 bonds (`C2-N1 amide rigid`,
  `N1-C3 conjugated single rotates`, `C6-O2 conjugated single rotates`); `c` made
  both conjugated bonds rigid (`0 rotatable bonds.`), `2,3` restored them.
- `8 → l` printed the warning `A ligand built from SMILES sits wherever RDKit put
  it, not in the receptor; a box around it is rarely what you want.` and produced
  `Box centre (-0.23, 0.2, -0.17)` — in empty space. Running it anyway:
  `Done in 0.7 s`, five poses around `-4.15 kcal/mol`, **no RMSD column** and no
  warning that the box does not overlap the protein. Re-doing the box as
  `8 → c` with centre `11.6 22.5 5.9` and edges `22 22 22` gave a plausible
  `-4.74` for a fragment squeezed beside the retained MK1.
- Problem to flag: a box that contains no receptor atoms still "succeeds". The
  only hint is the sub-second runtime and the weak score.

### Scenario 14: Search-box definition, all three modes — RUN

- `l` around the ligand, padding 8 → `edges 30.9 x 26.2 x 24.6 A (19925 A^3);
  the ligand as placed, padding 8 A`.
- `r` around residues → `Residues to enclose (e.g. A:43,A:88,A:128)` → typed
  `A:25,B:25,A:50` → padding Enter (8) →
  `Box centre (13.57, 23.63, 2.72), edges 26.9 x 30.1 x 24.3 A; residues A:25, A:50, B:25, padding 8 A`.
  A bad id (`A:999`) → `Not in the structure: A:999.` and the box is left unset.
- `c` typed → six prompts (`Box centre x/y/z (A)` default 0.0, `Box edge x/y/z (A)`
  default 20.0, minimum 1.0) → `22.0 x 22.0 x 22.0 A, typed`.

### Scenario 15: Library + docking campaign — RUN

- Library prepared as a plain SMILES file (three molecules, names in column 2):
  `<repo>/manual/runs/s8/ligands/mini_library.smi`.
- Loading it: `1a → 4` ("Load a ligand file or library (SDF/mol2/SMILES) for
  docking") → file browser → pasted the absolute path →
  `mini_library.smi: 3 molecule(s) / 0: paracetamol / 1: benzoic_acid /
  2: not_a_real_drug` → `Stored for Molecular Docking: ...`.
- Single-ligand "from file" path (`2 → f`): `Molecule number in the file (0 for
  the first) (0):` → `1` → **`The ligand could not be built: Unsupported ligand
  file type '.smi'; use .sdf, .mol or .mol2.`** — but the dashboard still showed
  `2 Ligand: mini_library.smi record 1`. So .smi works for a campaign and for the
  loader, but not for a single ligand.
- Campaign (`c`) prompts, in order: library summary (`3 molecules.` +
  "Protonation and formal charges are taken as written in the library: ProPrep
  does not change a ligand's protonation for pH"), then
  `Ligands whose stereocentres are not specified: [s] skip, [d] dock as built (s):`
  → `d`; `Records with 2D coordinates: [b] build 3D with RDKit, [s] skip (b)`;
  `Seed for 3D builds (RDKit ETKDG) (42)`; `Minimise 3D builds with MMFF94?
  [mmff94/none]`; `Hold every conjugated single bond rigid (ester C(=O)-O,
  polyene single bonds)? [y/n] (n)`; `Let Meeko open macrocycles so the ring
  flexes? [y/n] (y)`; `Building the receptor ...`; `3 ligand(s) to dock.`;
  `Start docking? (Ctrl-C stops; finished ligands are kept) [y/n] (y)` → `y`.
- Output:
  ```
  3 molecules in mini_library.smi; 0 already done.
  Prepared 3; skipped 0 (reasons in results.csv).
  [1/3] paracetamol: -5.34 kcal/mol (1 s)
  [2/3] benzoic_acid: -5.16 kcal/mol (1 s)
  [3/3] not_a_real_drug: -4.70 kcal/mol (1 s)
  Campaign: 3 docked, 0 skipped, 0 failed.
  ```
  plus a "Top 3 by score" table with `score (kcal/mol)`, `per heavy atom`,
  `heavy atoms`, `net charge` (ranking by score-per-heavy-atom is offered as a
  column, a nice touch against size bias).
- Result / files (`docking/campaign_mini_library_20261121_132752/`):
  `results.csv` (index, name, status, best_score, best_inter,
  score_per_heavy_atom, poses, heavy_atoms, rotatable_bonds, net_charge, seconds,
  notes), `top_poses.sdf`, `manifest.json` (library + receptor SHA-256, so an
  interrupted campaign can be resumed only if nothing changed),
  `receptor_rigid.pdbqt`, and `ligands/00001_paracetamol/…` per ligand
  (`ligand.pdbqt`, `poses.sdf`, `receptor_rigid.pdbqt`, `vina_output.txt`).
- Resume path (SOURCE-ONLY, `docking_prep/docking_menus_batch.py`): if a
  campaign with exactly the same manifest is found, it prints
  `An unfinished campaign made with exactly these choices is in … (N of M done).`
  and asks `Resume it?`.

### Scenario 16: Viewing and exporting results — RUN

- `v` with both a run and a campaign present first asks which:
  `The run: docking_20261121_132205` / `The campaign: campaign_mini_library_... (newer)`.
- Run view: `Run folder: …`, a `pose | total | inter | RMSD to crystal` table,
  then `Pose to show in the viewer (0 for none)` (0 = just the table; a number
  opens the browser viewer with that pose on the receptor and the box drawn).
  When the viewer is closed, the run itself prints
  `The viewer is not open; 'v' shows any pose in it.`
- Campaign view: the same summary panel + `Top 20 by score`, then
  `Rank of the ligand to show in the viewer (0 for none)`; choosing a rank writes
  `best_pose.sdf` next to that ligand's `poses.sdf`.
- Exporting is simply the files: `poses.sdf` (every pose, best first, each record
  carrying `rank`, `energy … (kcal/mol)` and `RMSD to reference (A)` as SD tags),
  `pose_N.sdf` per pose, `top_poses.sdf` for a campaign, `results.csv` for the
  table, and the `.pdbqt` files for anyone who wants to re-run Vina by hand.

### Options reference — 5d Molecular Docking

| Where | Choice | What it does |
|---|---|---|
| Dashboard | `1` Receptor | confirm the loaded structure, pick chains, choose alternate locations; warns about residues missing heavy atoms ("Rebuild them with the Structure Fixer first") |
| | `2` Ligand | `s` a HETATM residue of the structure (bond orders from the CCD), `m` SMILES (ETKDG seed + optional MMFF94), `f` a record of a loaded .sdf/.mol/.mol2 file; then free-form chemistry edits (`remove`, `charge`, `bond`, `add-h`, `p` pick, `drop`) |
| | `3` Kept HETATM | per-residue keep/leave-out, `w` for all waters at once, `p` to click one |
| | `4` Protonation | pH → pdb2pqr/PROPKA; review table of unusual residues + all His; per-residue template (HID/HIE/HIP, ASH, GLH …); chain ends as real terminus or chain break |
| | `5` Metals & cofactors | formal charges for kept metals, CCD chemistry review per kept cofactor |
| | `6` Flexible | residue ids whose side chains move (`none` / `keep` / `p`) |
| | `7` Rotatable | toggle individual bonds, `c` to freeze all conjugated single bonds, macrocycle question when a big ring is present |
| | `8` Search box | `l` around the ligand + padding, `r` around residues + padding, `c` typed centre and edges |
| | `9` Settings | scoring `vina` / `vinardo` / `ad4` (with a printed explanation of all three), seed (not 0), exhaustiveness (default 8), number of poses (9), energy range (3.0 kcal/mol), minimum RMSD between poses (1.0 Å), grid spacing (0.375 Å), CPU cores (0 = all) |
| | `r` Run docking | prepares both sides, runs Vina (+ AutoGrid for ad4), prints the pose table |
| | `c` Campaign | one receptor/box/settings against the loaded library, with its own stereo/2D/seed/minimisation/conjugation/macrocycle policies, resume support and `results.csv` |
| | `v` View results | pose tables for the last run and/or campaign, optional viewer |
| Scoring guide (printed) | vina | empirical, no charges/electrostatics, fast |
| | vinardo | re-parameterised Vina terms, same speed |
| | ad4 | AutoDock4 force field with Gasteiger charges and desolvation; needs AutoGrid maps per receptor+box, so slower; the screen cites ProPrep's own carbonic-anhydrase validation |

---

## Notes for the manual writer

**Jargon that needs a plain-English gloss:** bilayer / leaflet, APL (area per
lipid), `distxy_fix` vs padding, "pre-oriented", MEMEMBED, PPM3, OPM, CompEL /
charge imbalance, self-assembly packing, PACKMOL "loops" and "objective",
`prmtop`/`rst7`, PDBQT, Meeko, exhaustiveness, scoring function, RMSD to crystal
pose, HID/HIE/HIP and ASH/GLH, OXT.

**Recommended beginner path (membrane):** load the structure (1a) → **remove
ions, waters and anything exotic with 4a** → 4g: option 1 (use workspace
structure), option 2 → raw string with a single lipid such as `POPC`, option 4 →
padding 10-15 Å and water 15-17.5 Å, option 5 → MEMEMBED (or PPM3 if you want the
OPM comparison), option 8 → the guided walkthrough (ff14SB/ff19SB + lipid21 +
matching water), then `r` to read the review and `g` to build. Afterwards go
straight to 5a option 2 then option 4; 5a will not re-ask the force-field or
solvation questions.

**Recommended beginner path (docking):** 1a to load the complex → 5d option 1
(Enter for all chains) → option 2 → `s` and pick the crystal ligand → option 4
(Enter for pH 7, then **check the chain-end rows**) → option 8 → `l` with the
default 8 Å padding → option 9 → `vina`, exhaustiveness 8, 9 poses → `r`. If the
crystal ligand is present, the RMSD column tells you immediately whether the
settings are good enough (<2 Å = the method reproduces reality here).

**Common mistakes seen:**
1. Building a membrane around a file that still has metal ions or crystal waters —
   the tLEaP hydrogen step dies and the fallback produces a system ProPrep refuses
   to record (Scenario 4). Filter first.
2. Accepting "C-terminus: charged COO-" for a PDB with no OXT atoms — docking
   fails with a Meeko template error (Scenario 10). Answer `b` (chain break) for
   those ends, or add OXT with the Structure Fixer.
3. Using box mode `l` with a SMILES ligand — the box lands at the origin and the
   docking "succeeds" meaninglessly (Scenario 13).
4. Leaving the crystal ligand kept in the receptor (option 3) while docking
   something else into the same pocket — the site is blocked.
5. Typing a residue *range* (`1-97`) in 4a's residue selector — it crashes the
   tool (see below).
6. Choosing a big exhaustiveness/pose count or `ad4` for a first try; vina at
   exhaustiveness 4-8 answered in 17-30 s here, `ad4` took 59 s plus map building.

**Things that surprised me / look like bugs (described, not fixed):**
- *Config reset:* selecting a different protein in 4g (or filtering the structure
  and coming back) silently clears the lipid composition, box and orientation
  choices. Nothing says so.
- *`undo` side effect:* rewinding inside the force-field walkthrough also reverted
  a PACKMOL time limit set earlier in a different sub-menu.
- *Solvate-only wording:* with "No lipids" the output file is still
  `bilayer.pdb`, and the review still lists `Lipid: lipid21` /
  `leaprc needed: leaprc.lipid21` for a lipid-free system.
- *Topology Generator status line:* after the membrane build, 5a's suggestion read
  `✓ Topology files generated (1 prmtop files)` **before** anything was generated —
  it was counting the internal `pretleap_hydrogen.prmtop` written by the
  hydrogen-adding step, which also appears as a second row ("already have topology
  files") in the tLEaP file table at option 4. Users may think they are done.
- *4a PDB Filter crash:* `Select residues (comma-separated indices, 'all', or
  'none')` → `1-97` raises `invalid literal for int() with base 10: '1-97'`,
  logs a traceback and aborts the whole filter run.
- *Docking state after a failed ligand build:* the `.smi` rejection still left
  `2 Ligand: mini_library.smi record 1` on the dashboard.
- *No box/receptor overlap check* in docking (see mistake 3).
- *Program provenance:* screens show `packmol-memgen` and `teLeap` coming from
  `$AMBERHOME` (a *different* conda env than ProPrep's). If the manual tells
  people to check their install, it should mention `AMBERHOME` explicitly.
- *AD4 ranking:* in the one case tried, AD4's best-scoring pose was 4.05 Å from
  the crystal pose while its 5th pose was 1.78 Å; vina/vinardo got it right at
  rank 1. Worth a sentence that the top score is not automatically the answer.

**Nice touches worth keeping in the manual:** the review screen's "Equivalent
command" line (self-documenting, and copy-pasteable for cluster runs); the
plain-English packmol progress with a time estimate and the note that leftover
clashes are the minimiser's job; the OPM cross-check report; the protonation
summary that hides the 31 boring residues and explains why; the RMSD-to-crystal
column as a built-in sanity check; `results.csv`/`manifest.json` for campaigns.
