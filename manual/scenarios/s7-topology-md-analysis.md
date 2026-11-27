# Topology generation, molecular dynamics, and analysis (scope s7)

ProPrep 1.23.0, always started with `proprep --full-menu` (plus one `proprep --analysis` run).
Test proteins: **1UBQ** (ubiquitin; raw PDB with crystal waters), **1CRN** (crambin; 3 disulfides),
**1HRC** (cytochrome c; heme c), and a CYX-renamed crambin PDB loaded from a local file.

Sessions (transcripts are in `manual/runs/_sessions/<name>/transcript.txt`):

| Session | Working folder | What it covered | Session JSON |
|---|---|---|---|
| `s7-main` | `runs/s7/main/ubq` | 1UBQ: 5a bond editor, tLEaP in/out, tLEaP run, implicit-solvent failure, cpin, PB Titrate menu, the whole 5b setup, run, monitor, analysis, templates library, SLURM | `runs/s7/main/ubq/proprep_session_20261121_122842.json` |
| `s7-crn` | `runs/s7/crn/crn` | 1CRN: 4e first, rectangular per-axis box, two salts, disulfide failures, 2a+4d batch microstates, 5a option 3 | `runs/s7/crn/crn/proprep_session_20261121_132821.json` |
| `s7-crn2` | `runs/s7/crn2` | CYX crambin loaded as a local PDB, Constant pH (ff10) force field, disulfide bonds, cpin with pKa filter | `runs/s7/crn2/proprep_session_20261121_140524.json` |
| `s7-pair` | `runs/s7/pair/pair` | Loading an existing prmtop/rst7 pair, pmemd failure, the "step from scratch" wizard, a sander minimization | `runs/s7/pair/pair/proprep_session_20261121_134811.json` |
| `s7-hrc` | `runs/s7/hrc/hrc` | 1HRC heme straight into 5a (fails) | `runs/s7/hrc/hrc/proprep_session_20261121_135939.json` |
| `s7-ana` | `runs/s7/main/ubq` | `proprep --analysis` | `runs/s7/main/ubq/proprep_session_20261121_135814.json` |

**Environment note found while testing:** this machine already has `AMBERHOME` set to a different
conda env (`etanalyze_env`). ProPrep keeps an AMBERHOME that is already set (see
`utils/amber_env.py`), so every tleap, sander, cpptraj, ParmEd and pytraj call in these tests used
**etanalyze_env's AmberTools**, not the ProPrep env's. Every tLEaP message shows
`.../envs/etanalyze_env/bin/teLeap`. The manual should say: "activate the ProPrep environment (or
unset AMBERHOME) before starting ProPrep if you have another Amber installed."

---

## Feature: 5a. Topology Generator

**What it is for (plain English):** it turns your prepared protein (PDB) into the two files an Amber
simulation needs. The *topology* (`.prmtop`) describes every atom, bond and force-field parameter.
The *coordinates* (`.rst7`) give every atom's position. Along the way it adds a water box and salt
ions. It does this by writing an input script for Amber's **tLEaP** program and running it.

**Needs before you start:** a loaded structure (1a), a packed membrane (4g), or a prmtop/rst7 pair.
In practice you also need:
- crystal waters removed if you want implicit solvent;
- disulfide cysteines named **CYX** *and* SG–SG bond definitions;
- every non-standard residue (heme, metals, ligands) prepared through 2a → 4d (and 2c if needed).

Protonation (4e) is optional; the generator never asked for it. If 4e was run, its
`structure_with_md_names.pdb` is used automatically. AmberTools `tleap` and ParmEd must be
installed (both are bundled).

**Where to find it:** full menu key `5a` (section "5. SIMULATION SETUP & EXECUTION"). In Workflow
Mode it is in the **Simulate** stage. There is no command-line flag.

Menu (as shown with a structure loaded):
```
1. Edit redox site bond definitions
2. Generate tLEaP input for single state
3. Generate tLEaP inputs for all redox microstates   (needs batch microstate metadata from 4d)
4. Generate prmtop/rst7 files from tLEaP input files (needs a tLEaP input)
5. Refine protonation states via PB (PBSA)            (needs a topology)
6. Generate cpin/cein file for constant pH or redox MD (needs a topology)
[b]ack · [m]main · [x]exit
```
Disabled items have a ○ and a ⚠ reason. Choosing a disabled item just shows the menu again, with
no message. ● marks a step that is done.

### Scenario 1: Opening 5a with nothing loaded — RUN
- Starting point: new project `ubq`, nothing loaded.
- Steps: `5a` at "Enter your choice (1a):"
- What the user sees: the header reads "⚠ Needs a structure, a packed membrane, or a parm7/rst7 pair".
  Only option 1 is enabled; option 2 shows "⚠ Load a structure first". Typing `2` just shows the menu again.
- Session log: `runs/_sessions/s7-main/transcript.txt`

### Scenario 2: Bond editor (option 1) on a protein with no detected redox sites — RUN
- Starting point: 1UBQ loaded from RCSB (1a → 1 → 1 → 1 → `1UBQ` → 1 asymmetric unit → 6 Done).
- Steps:
  1. `1` → "No bond definitions found. Gathering from modules... No bond definitions found. Run comprehensive RedoxSite detection first."
  2. "Do you want to continue without bond definitions? [y/n] (n):" → Enter (n) → "Please run RedoxSite detection and try again." The **tLEaP Bond Editor Dashboard** still opens.
  3. `h` → help screen (filters d/m/p/o/all; a add; `#` edit; `x#` delete; s save; q quit). Bond format: `bond mol.ResidueID.AtomName mol.ResidueID.AtomName`, e.g. `bond mol.123.SG mol.456.SG`.
  4. `a` → "Enter tLEaP bond command:" → Enter on an empty line gives "Cancelled".
  5. `a` → `bond mol.1.N mol.76.C` → "Select category: 1 Covalent, 2 Disulfide (SG-SG), 3 Coordination, 4 Metal-Metal, 5 Peptide Backbone, 6 Other" (default 6) → `5` → "Added bond to Peptide Backbone".
  6. `1` (edit bond 1) → Enter keeps it → "No changes made".
  7. `x1` → "Confirm deletion? (n)" → `y` → "Bond deleted".
  8. `q` → "✓ Saved 0 bond definitions to workspace". **Quitting saves; it does not discard.**
- Problems: the editor accepts any text as a "bond". In session s7-crn a stray `n` became bond #2,
  literally `n`, with no validation (see Scenario 10). Watch out for the "Press Enter to continue"
  pause after each add: anything typed there is swallowed.
- Session log: `runs/_sessions/s7-main/transcript.txt`

### Scenario 3: Single-state tLEaP input, first pass (ff19SB / TIP3P, octahedron, 150 mM NaCl, edit in nano) — RUN, ends in a trap
- Starting point: 1UBQ raw RCSB file (it still contains 58 crystal waters), no protonation step.
- Steps (prompt → answer):
  1. `2`. An info panel explains tLEaP and the 5 stages.
  2. "Select protein forcefield (2):" → Enter (2 = **ff19SB**).
  3. "Select modified amino acids ():" → Enter (none).
  4. "Select dna forcefield (2):" → Enter. **Default is OL24, not None**, even for a protein-only system. The same goes for RNA (OL3), carbohydrates (GLYCAM_06j) and small molecules (GAFF2). All of them get loaded unless you type `1`.
  5. RNA → Enter; carbohydrate → Enter; lipid → `1` (None); small molecule → Enter (GAFF2).
  6. "Select water model (2):" → `6` (**TIP3P**). No warning appears, even though the list itself says ff19SB should be used with OPC.
  7. "Select divalent+ ion parameters (6):" → Enter. The recommended default follows the chosen water model (here "12-6-4 TIP3P").
  8. A summary of the selected force fields appears, then "No transformer info … This is normal if no transformations were performed."
  9. "Select solvent model [1/2] (1):" → Enter (explicit).
  10. "Select box shape [1/2] (1):" → Enter (truncated octahedron, "~29% fewer waters").
  11. "Buffer distance (Angstroms) (10.0):" → see Scenario 6 for invalid entries, then `10`.
  12. Salt #1: "Select cation (1)" Enter (Na+); "Select anion (1)" Enter (Cl−); "NaCl concentration" `5` (Custom) → "Enter concentration in mM (150)" → `150`.
  13. "Add another salt? (n)" → Enter.
  14. ProPrep runs a **tLEaP info pass** ("Net charge: +0.0, Water count: 4792"). It then shows a large panel, "Ion Calculation — 150 mM NaCl, SPLIT Method", with the formula and a caveat that NPT later raises the effective concentration to about 158–165 mM. Result: 13 Na+ + 13 Cl−.
  15. "tLEaP template saved". The template is printed with line numbers.
  16. "Edit the tLEaP template? [y/n] (y):" → Enter. The template opens in the terminal editor; on this Mac `nano` is really **UW PICO 5.09**. I typed a comment line, then Ctrl-X, `y`, Enter. ProPrep said "Template updated successfully!"
  17. "Using user-loaded RCSB PDB PDB file: …/1UBQ.pdb". STEP 2: "Molecule contiguity … repaired by ParmEd after the build". STEP 3: "No bond definitions found…"
  18. "Do you want to continue without bond definitions? [y/n] (n):" → Enter (n) → "Please run RedoxSite detection and try again." **You are dropped back at the 5a menu. No input file is written.**
- Problems: **the default answer `n` throws away the whole run.** For an ordinary protein with no
  special bonds the user must type `y`. The force-field and solvation choices are remembered (see
  Scenario 4). The template edit is kept too.
- Session log: `runs/_sessions/s7-main/transcript.txt` (from "Generate tLEaP Input for Single State")

### Scenario 4: Single-state input, second pass (re-use previous settings) — RUN
- Steps:
  1. `2`. It prints "Previous settings: protein: {'name': 'ff19SB', 'leaprc': … }". This is a **raw Python dictionary** dumped on screen, which is confusing to read. It ends with "Solvation: truncated octahedron, 10.0 A buffer, 150 mM NaCl".
  2. "Reconfigure all settings (forcefields, solvation, template)? [y/n] (n):" → Enter (keep).
  3. The template is shown, with my earlier comment line still at the top.
  4. "Edit the tLEaP template? (y)" → `n`.
  5. "Do you want to continue without bond definitions? (n)" → `y` → "Continuing with empty bond definitions."
  6. "Enter output filename for tLEaP input (1UBQ_tleap.in):" → Enter → "Successfully wrote tLEaP input file to: 1UBQ_tleap.in".
  7. "View the tLEaP input file? (y)" → Enter shows a boxed listing of the final file.
- Result: `runs/s7/main/ubq/1UBQ_tleap.in`. Its contents: `source leaprc.protein.ff19SB / DNA.OL24 / RNA.OL3 / GLYCAM_06j-1 / gaff2 / water.tip3p`, `loadamberparams frcmod.ions1lm_1264_tip3p` + `frcmod.ions234lm_1264_tip3p`, `mol = loadpdb <abs path>/1UBQ.pdb`, `solvateoct mol TIP3PBOX 10.00`, `addionsrand mol Na+ 13 Cl- 13`, `check mol`, `saveamberparm mol 1UBQ.prmtop 1UBQ.rst7`, `quit`. The placeholder sections are filled with "# No custom atom types needed" and similar comments.
- Note: output is named **`.prmtop`**, not `.parm7`.

### Scenario 5: Run tLEaP (option 4) on 1UBQ — RUN (success)
- Steps: `4` → a status table (tLEaP file | PDB files "All exist" | expected outputs | "prmtop ✗ rst7 ✗") → "Run 1UBQ_tleap.in? (y)" → Enter.
- What the user sees:
  - A list of 57 lines "Missing TER record between HETATM HOH 77 and HOH 78…", then "✓ Inserted 57 TER record(s) in 1UBQ.pdb". **ProPrep edits your PDB in place** and keeps `1UBQ.pdb.ter_backup`.
  - "Running: tleap -s -f 1UBQ_tleap.in ✓ tLEaP completed successfully". Summary: "Errors = 0; Warnings = 13; Notes = 1". All 13 warnings are "Close contact of 1.1–1.5 Å between nonbonded atoms H and H1", which is normal for freshly added hydrogens.
  - An info panel titled "About to run: ParmEd validation" (validity check, rediscover molecules, apply the 12-6-4 ion correction). Then: "✓ Initial validity check passed", "✓ Molecule definitions are correct", "Applying 12-6-4 LJ potential (water model: TIP3P)… Detected free ions: Cl-, Na+", "✓ Added LENNARD_JONES_CCOEF", "✓ Wrote 1UBQ_parmed.pdb".
  - A validation table showing prmtop ✓ and rst7 ✓, followed by "All expected topology files were generated successfully!"
- Files produced (`runs/s7/main/ubq/`):
  - `1UBQ.prmtop` (2.9 MB): the topology, overwritten by ParmEd with the 12-6-4 terms.
  - `1UBQ.rst7` (574 kB): solvated coordinates. 15,729 atoms: 76 residues, 4,824 waters, 26 ions.
  - `1UBQ_parmed.pdb`: a PDB in the same atom order as the topology.
  - `1UBQ_parmed_replay.py`: a script to redo the ParmEd pass by hand.
  - `leap.log`: the full tLEaP log.
  - `1UBQ.pdb.ter_backup`: the PDB as it was before TER records were inserted.
- Time: a few seconds.

### Scenario 6: Invalid answers at solvation prompts — RUN
- Buffer distance `abc` → "Please enter a number". `-5` → "Value should be between 5.0 and 25.0 Angstroms".
- Custom salt `5000` mM → "Concentration should be between 0 and 1000 mM".
- In every case the prompt simply repeats.

### Scenario 7: Implicit solvent (GB) with ff14SB on raw 1UBQ — FAILED
- Steps: `2` → "Reconfigure? (n)" → `y` → protein `3` (ff14SB) → `1` (None) for modified AA, DNA, RNA, carbohydrates, lipids and small molecules → water `1` (None = implicit). The ion menu default becomes "13 Default only (Recommended — no explicit water model)" → Enter. It then reports "No water forcefield selected → using implicit solvent (GB/PB continuum model)" and "Implicit solvent selected - no solvation section in template". No GB radii (`set default PBradii …`) are added. Then: edit `n`, continue without bonds `y`, filename `1UBQ_gb_tleap.in`, view `n`. Next, `4` → "All topology files already exist. 1 Regenerate all / 2 Select specific / 3 Skip (3)" → `2` → pick `1` (1UBQ_gb_tleap.in).
- What the user sees: "✗ tLEaP failed (exit code: 31)", "Errors = 1; Warnings = 121", warnings "Unknown residue: HOH … -no luck", then "Results: 0/1 files processed successfully".
- The real cause is only in `leap.log`: "FATAL: Atom .R<HOH 120>.A<O 1> does not have a type … Failed to generate parameters". **The FATAL line is not shown in ProPrep's summary.**
- **Data loss:** both input files write to the same `1UBQ.prmtop`. The failed run left `1UBQ.prmtop` at **0 bytes**, destroying the good explicit-solvent topology from Scenario 5. ProPrep never warns that two inputs share an output name.
- How to get past it: strip crystal waters first (4a PDB Filter), and give implicit-solvent runs their own output names by editing the template.

### Scenario 8: Option 4 when outputs already exist — RUN
- With several inputs: "Run options: All topology files already exist. 1. Regenerate all (N files) 2. Select specific files 3. Skip (default 3)". Option 2 prints a numbered list (✓ = topology exists, ○ = needs generation) and accepts "comma-separated numbers, 'all', or 'skip'".
- With one input: "Topology already exists. Regenerate X? [y/n] (n)".
- A 0-byte prmtop left by a failed run still counts as "✓ exists": the table shows "prmtop ✓ rst7 ✗", yet the summary says "already have topology files".

### Scenario 9: 1CRN after protonation (4e), ff19SB/OPC, rectangular per-axis box, KCl + MgCl2 — template RUN, tLEaP FAILED
- Starting point: project `crn`, 1CRN from RCSB.
  - 4e → 1 (PROPKA; probability threshold 0.5; net charge −0.14).
  - 4e → 3 "Set MD residue names": constant pH? `n`; set all standard? `n`; accept recommended per type; output `structure_with_md_names.pdb`.
  - **All 6 cysteines stay `CYS`**; none becomes CYX even though 1CRN has three SSBOND records.
- 5a → 1 (bond editor; "continue without bonds?" `y`; empty; `q`) → 2:
  - protein Enter (ff19SB); modified AA Enter; DNA, RNA, carbohydrates, lipids, small molecules all `1`; water Enter (**OPC**, recommended for ff19SB); ions Enter ("12-6-4 OPC (Recommended — matches OPC)" → loads `frcmod.ionslm_1264_opc`).
  - Explicit; box `2` Rectangular → "Buffer mode: 1 Uniform (Recommended) / 2 Per-axis" → `2` → X `10`, Y `10`, Z `15`.
  - Salt #1: K+ (`2`), Cl−, 150 mM. "Add another salt?" `y` → Salt #2: Mg2+ (`6`), Cl−, Custom `5` → `10` mM.
  - "Charge Neutralization — Choose which salt's ions are used to balance the protein charge via addions: 1 KCl / 2 MgCl2 (1)" → Enter.
  - Info pass: Q = +0, 4,164 waters. Panel "Ion Calculation — Multi-Salt": 11 K+ + 11 Cl−; 1 Mg2+ + 2 Cl−; total 1 Mg2+, 11 K+, 13 Cl−; "Charge balance verified".
  - Template: `solvateBox mol OPCBOX { 10 10 15 }`, `addionsrand mol K+ 11 Cl- 11`, `addionsrand mol MG 1 Cl- 2`.
  - Edit `n`; continue without bonds `y`; filename default `structure_with_md_names_tleap.in`; view `n`.
- 5a → 4 → Enter. "✓ TER records valid" → "✓ tLEaP completed successfully" (sic) → "Errors = 42": "Could not find bond parameter for atom types: SH - SH", "No torsion terms for 2C-SH-SH-HS", and so on. tLEaP read the PDB's CONECT records and bonded the SG atoms of CYS residues that still carry their HG hydrogens. Then: "There are missing parameters. Parameter file was not saved." → "✗ ParmEd validation failed: Could not identify file format" (a Python traceback is printed).
- **Misleading messages:** "Results: 1/1 files processed successfully", an empty validation table, "Some expected files are missing", and then the menu suggestion "✓ Topology files generated (1 prmtop files)". The prmtop is 0 bytes.

### Scenario 10: 1CRN — add the three disulfides by hand in the bond editor (still CYS names) — FAILED
- 5a → 1 → for each pair: `a` → `bond mol.3.SG mol.40.SG` → category `2` (Disulfide) → Enter. The same for 4–32 and 16–26. Then `d` (filter) shows the 3 bonds; `s` → "✓ Saved 3 bond definitions … Disulfide: 3 bonds".
- 5a → 2 → keep settings → filename `crn_ss_tleap.in`. The file now contains "# Disulfide bonds (CYS/CYX SG-SG)" followed by the three `bond` lines.
- 5a → 4 → "Select specific" → 1 → the same SH–SH parameter errors (60 errors). "Results: 0/1".
- Lesson for the manual: bond lines alone are not enough. The residues must also be named **CYX**.
  4e did not rename them; 4d with the `disulfide` transformer does (Scenario 11).
- Bond-editor bug seen here: typing `n` at "New bond command" saved a bond literally called `n`.

### Scenario 11: Multi-microstate batch (option 3) after 2a + 4d batch on 1CRN — FAILED (silently produces nothing)
- 2a → 2 (detect). Defaults were accepted; "Found 3 disulfide bond(s)". Per-site refinement: `finish` at "Select residues to add"; for sites 2 and 3 the SG–SG bonds were defined (`1-2` → atom `6` SG → target `1`). Then "5. Continue".
- 4d → 3 (batch): auto-assign failed ("No mapping found"), so the suggested `1-3:"disulfide"` was accepted (Enter, then accept `y`). "Enter selection (required)" → `single:2`, giving 4 microstates: OLOLOL, RLOLOL, OLRLOL, OLOLRL. Generate all.
  - Output: `transformed_microstate_00{1-4}.pdb`, `*_redox_sites_updated.json`, `microstate_metadata.json`.
  - **Possible bug (4d, outside this scope):** microstate 002 (site 1 "reduced") still names Cys3 and Cys40 as CYX.
- 5a now shows "2. ○ … ⚠ Batch microstates already generated" (so single-state is blocked) and "3. ●".
- 5a → 3 → template shown with `# PDB_FILE_SECTION` / `# SOLVATION_SECTION` placeholders → "Edit? (n)" → STEP 2/3/4 → "Found 4 microstates to process" → "Forcefield for UNK (disulfide, oxidized…): No forcefield files needed for disulfide" (twice) → **back at the menu. No `*_tleap.in` files written, no error message.**
- Cause (source, `tleap_input_generator.py` ≈ line 2814): `if not selected_forcefields: return False`. A batch where every site is disulfide-only has no force-field files to pick, so it aborts silently.

### Scenario 12: CYX crambin as a local PDB with Constant pH (ff10) — RUN (first silently wrong, then correct)
- Starting point: new session in `runs/s7/crn2` (project name Enter = current folder). 1a → 1 → 4 (local) → `crn_cyx.pdb` (a copy of microstate 001, all CYX).
- 5a → 2:
  - protein `8` "Constant pH (ff10)". The template sources `leaprc.constph`.
  - Water list now reads "6 TIP3P (Recommended — matches Constant pH (ff10)) … 3-point, traditional (not recommended)". Two contradictory labels on one line.
  - Defaults: TIP3P, 12-6-4 TIP3P, octahedron 10 Å, 150 mM NaCl (8 + 8).
  - **First attempt, no bond definitions:** tLEaP "Errors = 0; Warnings = 6", all "Illegal CONECT record in pdb file" (the CONECT serial numbers no longer match). ParmEd: "⚠ Detected CYX residue with a Sulfur atom not bonded to another CYX Sulfur! Did you forget to make the disulfide bond?". It still ends "Results: 1/1 files processed successfully". **The topology is built without the disulfides. Only the yellow ParmEd line hints at it.**
  - **Second attempt:** 5a → 1 → add the 3 `bond mol.X.SG mol.Y.SG` lines (category 2) → `q` → 2 (keep settings, edit n) → 4 → "Regenerate? (n)" → `y`. Result: "Errors = 0", "✓ Initial validity check passed", no CYX warning.
- Files: `runs/s7/crn2/crn_cyx.prmtop`, `crn_cyx.rst7`, `crn_cyx_tleap.in`, `crn_cyx_parmed.pdb`, `crn_cyx_parmed_replay.py`, `leap.log`.

### Scenario 13: Heme protein 1HRC straight into 5a — FAILED (expected)
- Project `hrc`, 1HRC from RCSB, no other preparation. 5a → 2 → ff19SB / OPC / 12-6-4 OPC, octahedron 10 Å, salt `4` "Neutralize only (0 mM)".
- The info pass said "Net charge: +8.0, Water count: 5703", giving `addions mol Cl- 0` ("places 8 Cl- near positive charges"). It did **not** warn that the heme (HEC) is unknown.
- Continue without bonds `y`; 4 → Enter → "✓ Inserted 124 TER record(s)" → "✗ tLEaP failed (exit code: 31)", "Errors = 159; Warnings = 23".
  - Warning 1: "Unknown residue: HEC number: 105 type: Terminal/last".
  - leap.log: "FATAL: Atom .R<HEC 105>.A<FE 1> does not have a type" for every heme atom.
  - Errors: "Could not find bond parameter for atom types: SH - " etc. These come from the heme-thioether links to Cys14/Cys17, made through CONECT records.
- The menu then still says "✓ Topology files generated (1 prmtop files)". The prmtop is 0 bytes.
- The intended route is SOURCE-ONLY here (and belongs to the redox scopes): 2a detect → 4d with a heme transformer. Bundled parameter sets are in `forcefield_params/specialized_residues/heme/`: `bis_his_b_type`, `bis_his_c_type`, `cys_axial_b_type`, `his_met_axial_c_type`; the last one fits cytochrome c. Back in 5a, the "Select Forcefield Parameters" step then fills `# ATOM_TYPES_SECTION` / `# FORCEFIELD_PARAMETERS_SECTION`. I started 2a on 1HRC, but its site-refinement dialogue (an "ACE A:0" site and boundary-atom questions) is long and is another scope's job, so I stopped.

### Scenario 14: Starting from an existing prmtop/rst7 pair — RUN
- New project `pair`. 1a → `2` "Load AMBER topology & coordinate files" → a file browser (N, `..`, `<path>`, q).
  - I typed the folder path `…/batch_20261121_125742/shared` → `2` (1UBQ.prmtop).
  - "Coordinate source: [r]st7/inpcrd or [t]rajectory (.nc/.mdcrd/…) (r)" → `r`.
  - The browser **goes back to the project folder** instead of the topology's folder, so I typed the full rst7 path.
  - "✓ Loaded topology + rst7 pair into workspace. Use the Molecular Dynamics Manager or QM/MM Preparator to proceed."
- Main menu: 5a and 5b now show ✓. Inside 5a, though, options 2–6 are all still disabled ("Load a structure first", "Need to generate topology files first"). **PB refinement and cpin generation cannot use a loaded pair.** The pair is usable only from 5b (Scenario 5b-21).

### Scenario 15: cpin for a normal ff19SB topology (option 6, explicit) — RUN (allowed, but chemically incomplete)
- 1UBQ: 5a → 6.
  - "Step 1: Simulation Type: 1 Implicit (GB) / 2 Explicit (requires custom radii in prmtop)" → `2`.
  - "✓ Found 8 titratable residues… 7 LYS, 1 TYR. pKa values: all from cpinutil --describe (residue-type defaults; run PB Titrate to get per-site values)".
  - Selection menu: 1 all / 2 by type / 3 by number / 4 exclude / 5 pKa range → Enter (all) → "Confirm selection? (y)".
  - Initial states: 1 cpinutil defaults / 2 by target pH → Enter → filename `1UBQ.cpin`.
- Runs `cpinutil.py -p 1UBQ.prmtop -op 1UBQ_cpin.prmtop -o ./1UBQ.cpin -system system -resnums 6 11 27 29 33 48 59 63`. ParmWarning: "change RADII: :AS4,GL4,PRN… matches no atoms".
- Files: `1UBQ.cpin`, `1UBQ_cpin.prmtop`, `1UBQ_cpin_generate.sh`.
- **Side effect:** "Updated md_structure_pairs: 1UBQ.prmtop → 1UBQ_cpin.prmtop". From then on, 5b offers "Use 1UBQ (cpin-ready)" by default and asks "Enable Constant pH MD for production steps? (y)".
- Problem: the topology was not built with the Constant-pH force field, so ASP/GLU/HIS cannot titrate. ProPrep did not warn about this.

### Scenario 16: cpin for a Constant-pH topology with pKa filter and target-pH states — RUN
- crn2 (Scenario 12 topology): 5a → 6 → `2` explicit → "✓ Found 2 titratable residues… 2 TYR". ASP43 and GLU23 are missing because they were never renamed AS4/GL4; 4e → 3 → "constant pH? y" is presumably needed.
- `5` (pKa range): "Minimum pKa (blank = none)" Enter, "Maximum…" Enter → "Filtered selection: 2 residues". Confirm.
- Initial states `2` → "Target pH (7.0)" `7.0` → "TYR (pKa 9.6): pH 7.0 < pKa → State 0 = OH (protonated, neutral)".
- Filename `crn_cyx.cpin`. ProPrep says "ProPrep applies these to the production steps automatically: -cpin -cpout -cprestrt; icnstph=2, solvph, ntcnstph".
- Files: `runs/s7/crn2/crn_cyx.cpin`, `crn_cyx_cpin.prmtop`, `crn_cyx_cpin_generate.sh`.

### Scenario 17: PB Titrate (option 5) — menu RUN, calculations SOURCE-ONLY
- 5a → 5 opens a checklist:
  1. Configure PB Parameters
  2. Minimize Structure (optional)
  3. Discover Sites
  4. Single-Site pKa
  5. Coupling Precompute (optional)
  6. Solve (mean-field / Monte Carlo / exact)
  7. Review & Override
  8. Apply to Prmtop
  9. Write modern-FF PDB (optional)
  10. Persist Recommendations
  11. Compare pKa with ProPKA (optional)
- Keys: `n` next, `v` view structure, `s` save & quit, `q` quit, `#` run step, `#s` skip, `i`/`#i` info.
- `i` shows the defaults: εin = 4, ionic strength 150 mM, grid 0.5 Å, 2 focusing levels, 4 workers, pH 7.0, cluster radius 25 Å. It also explains the net-charge issues in constant-pH MD and lists literature.
- `q` → "Quit without saving? (Progress will be lost) (n)". Enter means *no* and returns to the checklist; `q` again then `y` leaves.
- Opening it creates `pb_titrate_work/` in the project. I did not run the PB steps: they are long-running and big enough for their own scope.

### Options reference (5a)

| Prompt | Choices (default in **bold**) | Notes |
|---|---|---|
| Protein FF | 1 None, **2 ff19SB** (use OPC), 3 ff14SB, 4 ff15ipq (needs SPC/Eb), 5 fb15 (TIP3P-FB/TIP4P-FB), 6 ff03.r1, 7 ff99SBildn, 8 Constant pH (ff10), 9 Constant Redox (ff10), 10 Constant pH + Redox (ff10) | 8–10 load `leaprc.constph` etc. |
| Modified AA (multi) | **none**, 2 Modified AA ff19SB, 3 Modified ff14SB, 4 Phospho ff19SB, 5 Phospho ff14SB, 6 Phospho ff99SB, 7 Mimetic (ff15ipq), 8 Fluorinated (ff15ipq) | comma-separated |
| DNA FF | 1 None, **2 OL24**, 3 OL21, 4 OL15, 5 bsc1, 6 tumuc1 | note: only OL24/OL21 are compatible with ff19SB |
| RNA FF | 1 None, **2 OL3**, 3 LJbb, 4 ROC, 5 YIL, 6 Shaw, 7 modrna08 | |
| Carbohydrate FF | 1 None, **2 GLYCAM_06j**, 3 GLYCAM_06EPb | |
| Lipid FF | 1 None, **2 lipid21**, 3 lipid17 | |
| Small molecule FF | 1 None, **2 GAFF2**, 3 GAFF | |
| Water | 1 None (implicit), 2 OPC, 3 OPC3, 4 OPC3-pol, 5 TIP4P-Ew, 6 TIP3P, 7 SPC/E, 8 SPC/Eb, 9 TIP5P, 10 TIP3P-FB, 11 TIP4P-FB | default = recommended for the protein FF (OPC for ff19SB, TIP3P for ff10 constant pH; after a previous run it is the last choice) |
| Divalent+ ions | 1 None, 2 12-6-4 OPC, 3 12-6 OPC, 4 12-6 OPC3, 5 IOD OPC, 6 12-6-4 TIP3P, 7 12-6-4 SPC/E, 8 12-6-4 TIP4P-Ew, 9 12-6 TIP3P, 10 IOD TIP3P, 11 12-6 SPC/E, 12 Joung-Cheatham TIP3P, 13 Default only | default follows the water model; 12-6-4 sets trigger a ParmEd add12_6_4 pass |
| Solvent model | **1 Explicit**, 2 Implicit (GB/PB) | choosing water "None" forces implicit |
| Box shape | **1 Truncated octahedron**, 2 Rectangular | octahedron: one uniform buffer only |
| Buffer mode (rectangular only) | **1 Uniform**, 2 Per-axis (X, Y, Z) | range 5–25 Å, default 10 |
| Cation | **1 Na+**, 2 K+, 3 Li+, 4 Cs+, 5 Rb+, 6 Mg2+, 7 Ca2+, 8 Zn2+ | divalent ions give a 1:2 salt |
| Anion | **1 Cl−**, 2 Br−, 3 F−, 4 I− | |
| Concentration | **1 150 mM**, 2 100 mM, 3 50 mM, 4 Neutralize only (0 mM), 5 Custom (0–1000 mM) | SPLIT method (Machado & Pantano 2020) |
| Add another salt? | **n** / y | with ≥2 salts: "Use which salt for neutralization" |
| Reconfigure all settings? (re-run) | **n** keep / y | |
| Edit the tLEaP template? | **y** (opens $EDITOR/nano) / n | edits persist for later runs |
| Continue without bond definitions? | **n** (aborts!) / y | type y for normal proteins |
| Output filename | **`<pdbstem>_tleap.in`** | output prmtop/rst7 names come from the template (`<stem>.prmtop`) |
| View the tLEaP input file? | **y** / n | |
| Option 4 run | Run X? **y**; or Regenerate all / Select specific / **Skip** | |
| Option 6 type | **1 Implicit (GB)**, 2 Explicit (custom radii) | |
| Option 6 residues | **1 all**, 2 by type, 3 by number, 4 exclude, 5 pKa range | |
| Option 6 initial states | **1 cpinutil defaults**, 2 by target pH | |
| Bond editor | d / m / p / o / all filters; a add; `#` edit; x# delete; s save; q quit (saves); h help | categories: 1 Covalent, 2 Disulfide, 3 Coordination, 4 Metal-Metal, 5 Peptide, **6 Other** |

---

## Feature: 5b. Molecular Dynamics Manager

**What it is for (plain English):** after you have a topology, this tool sets up and runs the actual
simulation. Usually that means four phases: relax clashes (minimization), warm up to body or room
temperature (heating), let the water box settle (equilibration), then collect data (production). It
writes the Amber input files, runs them on your computer or writes cluster job scripts, shows live
progress, and analyses the results (RMSD, RMSF, H-bonds and more).

**Needs before you start:** a prmtop + rst7 pair, either made in 5a or loaded in 1a → 2. Also an MD
engine: here only **sander / sander.MPI** (pmemd and pmemd.cuda are not installed). Analysis needs
pytraj and cpptraj (bundled).

**Where to find it:** full-menu key `5b`; Workflow Mode → Simulate stage. `proprep --analysis`
jumps straight to the analysis browser.

Menu:
```
1. Setup and configure simulations
2. Execute simulations          (needs setup)
3. Monitor simulation progress  (needs a running simulation)
4. Analyze completed simulations
5. Manage Templates & Protocols
6. Manage cluster profiles (HPC configurations)
7. Manage run plans (per-protocol resource assignments)
```
Setup is a 6-step wizard with a breadcrumb bar:
`[Files] → Protocol → Restraints → Queue → Hardware → Review`.

### Scenario 1: Setup step 1, structure files — RUN
- 5b → 1. "Structure Pairs: 1 selected — 1. 1UBQ_cpin.prmtop + 1UBQ.rst7" (picked automatically after Scenario 5a-15). Options: 1 Use 1UBQ (cpin-ready), 2 Create custom pair, 3 Browse, 4 Remove pair; n Next; h Help; x Exit.
- `h` explains topology vs coordinates and the match types ("exact match", "modified" for constant pH).
- `2` → "Select Topology File: 1. 1UBQ.prmtop (workspace) 2. 1UBQ.prmtop (workspace) 3. 1UBQ_cpin.prmtop (directory)". **The same file is listed twice.** `1` → coordinates (also listed twice) `1` → "✓ Added custom pair: 1UBQ_1UBQ".
- `4` → choose `1` (1UBQ_cpin) → removed. Then `n`.
- Typing Enter at "Enter choice:" (it has no default) gives "Received '', which is not one of the available options".

### Scenario 2: Predefined protocol "Workshop (CPU, under an hour)", edited down to a few minutes — RUN
- "Protocol Source: 1 predefined / 2 load saved / 3 build custom from catalog / 4 create step from scratch (wizard); h overview; b back; x exit" → `1`.
- The protocol list is in the options reference. `4` → description: "measured 41 min with sander.MPI on 4 cores … ~75 min serial". Steps: min 1000 cycles; heat 20 ps; NPT 100 ps; NVT production 300 ps; total 420 ps.
- Commands: `p <n>` preview, `u` use as-is, `s` select specific steps, `e` edit, `b` back.
- `p 1` → the full annotated mdin (imin=1, ntmin=1, ncyc=500, maxcyc=1000, ntr=1, restraint_wt=10, restraintmask='@CA,C,O,N&!:WAT').
- `e` → "Protocol name (Protocol_1UBQ_1UBQ_20261121_1243)" Enter → the workflow editor. Its commands: `edit <n>`, `name <n>`, `files <n>`, `add before/after <n>`, `remove <n>`, `move <a> to <b>`, `preview`, `info`, `save`, `cancel`, `help`.
  - `edit 1` → shows name, type, description, Input `1UBQ.rst7` / Output `001_minimize.rst`, and the parameter list. Sub-menu: `params` (opens the text editor), `name`, `desc`, `files`, `back`. In `params` I set `ncyc=100, maxcyc=300`. (My first cursor edit deleted imin/ntmin by mistake. **ProPrep saved an mdin without imin and gave no warning**; I fixed it.)
  - `edit 2` → replaced the &cntrl/&wt block: nstlim=1000 (2 ps), ntwx=50, ntpr=50, ramp 0→300 K over 750 steps then hold.
  - `remove 3` → "Remove 'NPT Equilibration'? (n)" `y`.
  - `files 3` → "Input coordinate file (step_2.rst)", "Output (step_3.rst)". Editing showed a confusing naming mismatch: the edit screen shows `001_minimize.rst`/`002_heating.rst`, but `preview` shows the chain as `step_1.rst → step_2.rst`. When I typed `002_heating.rst` as step 3's input, the preview chain broke (step 3 read `002_heating.rst`, which no step writes). Setting it back to `step_2.rst`/`step_3.rst` fixed the preview. (Removing a step renumbers the chain automatically; `files` overrides it.) The run itself later used `step1/step1.rst7` and similar, so these names are cosmetic, but they mislead.
  - `edit 3` params → nstlim=500, ntwr=250, ntwx=25, ntpr=25.
  - `info` → asks for workflow name and description. **It edits them; it does not "show info".**
  - `preview` → "File Flow" table with `pmemd -i step_N.in …` commands (it says pmemd even though sander was chosen later), then "Preview template content for a step? (n)". That y/n prompt rejects other commands ("Received 'files 3'; please enter y or n").
  - `save`.
- "Step 2 Complete — 1 protocol(s) configured" (table: Structure | Protocol | Source = predefined). Then "Constant pH MD Available — CPIN file 1UBQ.cpin … Enable Constant pH MD for production steps? (y)" → `n`. Then "✓ Added 3 simulations to queue".
- The saved protocol goes to **`~/.proprep/md_workflows/user/custom/`**, which is shared by all projects. I deleted it afterwards (see Scenario 18).

### Scenario 3: Restraints step — positional mask generator and per-step weights — RUN
- Step 3 table: # | Step | Mask | Weight | DISANG. Options: (e) edit positional, (g) generate mask, (d) DISANG, (r) GROUP, (s) continue; (b) back, (x) exit.
- `g`:
  - "Extracting PDB from topology using cpptraj… ✓ 1UBQ_extracted.pdb". Numbering panel: "Chain A: residues 1-102 (+4824 water)". The 26 ions are counted as residues.
  - "No redox sites detected — (a) General structural restraints" → `a` → atom types: (b) backbone @CA,C,O,N, (h) heavy !@H=, (c) CA, (u) custom → `c`.
  - "Apply to which chains? (all)" → `A:1-70` → "Exclude water? (y)" `y`.
  - Preview "FINAL COMBINED MASK: @CA&:1-70&!:WAT". (a) accept / (e) start over / (c) cancel → `a`. Saved as `redox_restraint_mask`.
- Then "Edit Positional Restraints — Select steps (e.g. 2-4, all) (all)" → `2`. Presets are shown (`none` clears) → Enter keeps the generated mask → "Step 2 — System Heating [1.0]" → `2.0` → (s) accept.

### Scenario 4: DISANG restraint, interactive entry — FAILED
- `d` → "MD RESTRAINT CONFIGURATION" menu:
  - 1 Import DISANG file; 2 Configure interactively; 3 Display; 4 Export; 5 Mask; 6 GROUP; 7 Display GROUP; 8 Export GROUP; `d` done; `x` exit.
- `2` → `a` (add) → selection mode (s sequence / n non-standard / r redox sites / **m manual**) → `m` → `A:1 A:76` → "Residue A:1 not found". `1 76` → "Invalid format '1'. Use 'CHAIN:RESID'". `:1 :76` → "Residue :1 not found".
- Cause: the PDB extracted from the topology has **no chain IDs** ("Structure chains: 1 chains ( )"), so no manual entry can match. Sequence view (`s`) lists all 4,926 residues, mostly water.

### Scenario 5: DISANG restraint, import from file, export, apply to a step — RUN (with an error message)
- Wrote `runs/s7/main/ubq/s7_test.disang` (`&rst iat=5,1226, r1=30, r2=35, r3=40, r4=45, rk2=5, rk3=5, /`; atoms 5 and 1226 are the CA of Met1 and Gly76). The importer only lists `*.disang` files in the project folder, so a file named `.RST` was not found.
- `1` → "Suggested: s7_test.disang" Enter → **"✗ Error importing restraints: 'MDRestraintModule' object has no attribute 'get_from_workspace'"**. Yet `3` (Display) then shows "1 Distance … Flat: [35.0-40.0], k=5.0/5.0 Active" with an empty Atoms column.
- `4` → "DISANG filename (restraints.disang)" Enter → "✓ … 1 active restraints written". The file is valid.
- `d` → "Select steps for DISANG integration … (all)" → `3` → "Restraint dump frequency (500)" → `25` → "✓ DISANG applied to 1 step(s)". The table shows `restraints.disang` on step 3.
- In the run: step 3's mdin got `nmropt=1`, `&wt type='DUMPFREQ', istep1=25`, `DISANG=../../../restraints.disang`, `DUMPAVE=restraints_dump.txt`, `LISTOUT=restraints_violations.txt`. sander read it ("Restraints will be read from file"). The dump file shows a CA–CA distance of 37.3–37.6 Å. It works despite the error message.

### Scenario 6: GROUP restraints — RUN (configured, not applied)
- `r` → info panel ("cannot be combined with restraintmask in the same step") → "Continue? (y)":
  - "How many restraint groups? (1)" (1–10) → title "Group 1" → force constant (10.0) → `5.0`.
  - "How many FIND criteria? (1)" → `CA * * *` (format: atom_name atom_type tree_name residue_name; `*` = wildcard).
  - "Residue ranges" (≤ 7) → `1-76`.
- Preview of the GROUP block (`Group 1 / 5.0 / FIND / CA * * * / SEARCH / RES 1 76 / END / END`). "Save? (y)" → `y`. A developer log line leaks onto the screen: "WARNING - This method is deprecated. Use set() instead."
- "Apply GROUP to which steps? (all)" → `none`.

### Scenario 7: Queue management — RUN
- Step 4 table (# | Structure | Simulation | Step | Status ACTIVE). Options: k keep, r reorder, t toggle; b, x.
- `t` → commands `hold n`, `active n`, `hold 1,2`, `active 1-4`, `hold all`, `active all`, `b`. `hold 2` → HOLD; `active 2` → ACTIVE.
- `r` → "Enter new order (e.g. '1,3,2,4' or '6-19,1-5')" → `1,2`.
- With only one simulation: "✓ Only one simulation — proceeding to hardware configuration", and the prompt accepts only `n`, `b`, `x`.

### Scenario 8: Hardware configuration — RUN
- Panel explaining CPU for minimization and NPT, GPU for NVT/production. "Available: 11 CPUs, 0 GPUs". Options: 1 recommended, 2 per-simulation, 3 all the same, p load cluster profile, l load run plan, s save run plan; n next, h why CPU vs GPU, b, x.
- `1` (recommended) → minimization `sander.MPI`, heating and production **`pmemd.cuda`**, even with 0 GPUs and pmemd not installed. "Number of MPI processes (1-11) (5)" → `4`. "No GPUs detected on this machine. Defaulting to GPU 0. This is expected if setting up on a login node…"
- `3` (all the same) → list "1 sander, 2 sander.MPI, 3 pmemd, 4 pmemd.MPI, 5 pmemd.cuda". **This list does not mark the missing engines** → `1` → "3x sander".
- `2` (individual) → **here missing engines *are* marked "(not installed)"**. sim 1 → `2` sander.MPI, 4 processes; sim 2 → sander.MPI, 4; sim 3 → Enter (keep sander).
- `p` → "1. TAMU ACES (tamu-aces, bundled)" (`v1` shows the JSON, `c` cancels) → `1` → profile card: SLURM; modules GCC/13.2.0, OpenMPI/4.1.8, Amber/24; resource classes cpu-mpi (pmemd.MPI, 96 tasks), gpu-h100 (pmemd.cuda.DPFP), gpu-a30, gpu-debug. "Validation warnings: Missing required user fields: account". "Use 'tamu-aces' as the active cluster profile? (y)" → `n` → `c`.
- `h` → full explanation of the GPU fixed-precision (SPFP) force clipping and the fixed PME grid.

### Scenario 9: Review & Save, "Write MDIN only" — RUN
- Review header: Protocol, Cluster profile "(none — manual mode)", Run plan. Summary table: # | System | Simulation | Ensemble | Sim Time | Engine.
  - **Sim Time showed the original 1000 cycles / 20 ps / 300 ps, not my edited values.** The written files had the edited values.
- Options: 1 Save setup to workspace (execute later), 2 Write MDIN only, 3 Write files and run locally now, 4 Write MDIN + SLURM scripts (manual wizard), 5 Write and submit via sbatch (needs a profile); b, x.
- `2` → "Created run directory …/simulations/batch_20261121_125651", "✓ Copied structure files", "✓ Created stepN/simulation.mdin" → "Add extended production run cycles? (n)" → `n` → "✓ Created master run script (run_workflow.sh) ✓ README".
- Files: `batch_20261121_125651/{1UBQ.prmtop, 1UBQ.rst7, README.txt, run_workflow.sh, step1..3/simulation.mdin}`. `run_workflow.sh` chains `mpirun -np 4 sander.MPI -O -i simulation.mdin -p ../1UBQ.prmtop -c ../1UBQ.rst7 -o step1.mdout -r step1.rst7 -x step1.nc -ref ../1UBQ.rst7 -inf step1.mdinfo`, then step2 from `../step1/step1.rst7`, and so on. It skips a step whose rst7 already exists.
- Note: positional-restraint rewriting is text-based. On my one-line &cntrl it left the fragment `ntr=1,  !:WAT',` (harmless, because `!` starts a comment, but ugly).

### Scenario 10: Execute the queue (minimization → heating → short production) with sander — RUN (success)
- 5b → 2 → "Execute Options: 1 Execute queue (3), 2 Manage queue; 3 Monitor; 4 Analyze; b". `1` → queue list → "Proceed with execution of 3 active simulations? (y)" → `y`.
- **Execution makes a fresh run folder** (`batch_20261121_125742`); it does not reuse the "MDIN only" one. Each step gets its own folder (`01_1UBQ_1UBQ_Energy_Minimization/`, `02_…_System_Heating/`, `03_…_NVT_Production_MD/`, plus `shared/`). Each folder holds `simulation.mdin`, `run_simulation.sh`, `config.json`, `simulation.stdout/stderr`, `mdinfo`, `*.mdout`, `*.rst7`, `*.nc`.
- Each step's input coordinates are copied in as `1UBQ.rst7`, even though for steps 2 and 3 it is really the previous step's restart (checked by md5). The naming is confusing but the chaining is correct.
- Screen: "Starting simulation in background: mpirun -np 4 sander.MPI … ✓ started (PID 4080)". Steps 2 and 3 show "⏸️ Queued (depends on …) Will auto-start when dependency completes", then "Started background protocol manager". The Execute menu then shows "▶ Running: Energy Minimization ⏸ Next: System Heating (+1 more pending)" and later "✓ Finished" for all 3.
- Timings for 15,729 atoms on this Mac:

  | Step | Engine | Work | Wall time | Rate |
  |---|---|---|---|---|
  | Minimization | sander.MPI ×4 | 300 cycles | 6.0 s | |
  | Heating | sander.MPI ×4 | 1,000 steps (2 ps) | 22.2 s | 22 ms/step, ≈7.7 ns/day |
  | Production | sander serial | 500 steps (1 ps) | 33.2 s | |

- Results: final heating temperature 291.7 K; production average 294.5 K. Heating trajectory 20 frames, production 20 frames.

### Scenario 11: Monitoring a running simulation — RUN
- Execute menu → `3` (or 5b → 3) → "Currently tracked simulations: 1UBQ_1UBQ_System Heating (running, runtime 0:00:10)" plus "Other recent simulations".
- Monitor menu:
  - Data: 1 status, 2 energy summary, 3 temperature summary, 4 pressure summary, 5 timing.
  - ASCII plots: 6 total energy, 7 temperature, 8 pressure, 9 kinetic energy, 10 potential energy, 11 RMS gradient.
  - 12 periodic-image distance.
  - r refresh, b back.
- `1` → table: Step 1,000, Time 2.000 ps, Temperature 291.68 K, Total Energy −40,661.52.
- `7` → ASCII temperature ramp from 0 to 291.68 K (21 points).
- `5` → "No timing information available" (seen just after the step finished).
- `12`:
  - Explanation, then solute choice: 1 protein residues 1–76, 2 everything not solvent or single-atom ion, 3 specific residues, 4 custom mask → `1`.
  - "Analyse every how many frames (1)" → Enter. The cutoff is read from the mdout (10 Å). "Nonbonded cutoff to compare with (Enter = shown)" → Enter.
  - Result: truncated octahedron, edges 61.99 Å; closest approach 29.28 Å (frame 17, ARG 74 HH22 ↔ image of LYS 63 HZ2); mean 30.01 ± 0.35 Å; 0 of 20 frames below the cutoff; "never comes within the 10 Å cutoff…". Plus an ASCII plot.

### Scenario 12: Energetics analysis (mdout) — RUN
- 5b → 4 → "Available simulations: 1. 03_…Production (completed 12:58) 2. 02_…Heating 3. 01_…Minimization" → "Select simulation (1-3) or 'browse'".
- Production → "Select Analysis Type: 1 Energetics (mdout), 2 Trajectory (structural), 3 Combined (recommended), 4 View trajectory in the structure viewer" → `1`.
  - "Simulation Type Detection … Detected: NVT Equilibration/Production (confidence: medium) Evidence: ntp=0, ntb=1". Choose 1 Minimization, 2 Heating NVT, 3 Heating NPT, 4 NVT Eq/Prod, 5 NPT Eq/Prod → Enter.
  - **The "NVT Equilibration/Production Metrics" table was empty.** Plots: temperature vs time, total energy vs time, temperature histogram (mean 294.50 K).
- Minimization (selecting `3`, no trajectory, so it goes straight to energetics): detected with high confidence. Table: Total Steps 300; Initial −26,515 → Final −60,537 kcal/mol; Final RMS gradient 1.1249; GMAX 52.15; "Monotonic Decrease: Yes". Plots of energy and RMS gradient.
- Heating with `3` Combined: **misdetected as "NVT Equilibration/Production"**; I corrected it to `2` Heating (NVT). Table: total 2.00 ps, initial 0.0 K, final 291.7 K, "Avg Final Temp (last 20%) 282.7 ± 7.0 K (target 300 K)", "Energy Drift +23,836 kcal/mol" (an energy drift during heating is expected). Then it continues into trajectory analysis.

### Scenario 13: Trajectory analysis (pytraj) on the 1 ps production — RUN
- Analysis type `2` → file browser (N select, `N,M`/`N-M` multi-select, `..`, `find`, `q`). `find` searched the step folder and found 1 file → `1`.
  - "Trajectory Information: 20 frames, 15,729 atoms, 4,926 residues, 2.0–3.0 ps, periodic".
- The menu has 23 analyses (see options reference). Answer with a comma-separated list. Region menus are shared: 1 backbone `@C,CA,N,O&!:WAT`, **2 C-alpha**, 3 all atoms, 4 protein residues 1–76 no-H, 5 specific residues, 6 custom mask.
- What I ran (prompts → my answers → key result):
  - **1 RMSD**: region `1` backbone; reference 1 frame 0 / 2 average / 3 specific → frame 0. Final 0.65 Å. "Analyze additional regions? (n)".
  - Export: "1 CSV, 2 JSON, 3 both, 4 don't (4)" → `3` → `runs/s7/main/ubq/03_1UBQ_1UBQ_NVT_Production_MD_trajectory_analysis.csv` and `.json`. Later exports **overwrite the same file names**; the CSV holds every analysis run so far in the session.
  - **2 RMSF**: CA; alignment mask (@CA,C,N); reference frame 0. Mean 0.28 Å; top 10%: residues 10, 15, 46, 61, 65, 69, 75, 76 (0.72 Å); per-residue plot.
  - **9 Distance**: `:1@CA` to `:76@CA`, description "N-C terminus". Mean 38.44 ± 0.70 Å, histogram.
  - **21 Radius of gyration**: region 4 → 11.90 ± 0.03 Å.
  - **6 DSSP**: region 4 → α-helix 19.5%, β-sheet 33.9%, turn 19.0%, coil 27.6%; time plots.
  - **3 Contacts**: contact type (1 = CA), cutoff 4.5 Å, persistence 0.5, native reference 1 first / 2 last / 3 specific → 162 native contacts, 158 persistent, Q-value table.
  - **4 Salt bridges**: cutoff 4.0 Å (relaxed 5.0, strict 3.2), persistence 0.5 → 11 acidic / 11 basic residues, 8 salt bridges, all persistent; distance trace.
  - **5 SASA**: CA region → mean 2,940 Å² (SASA on CA atoms only is not meaningful; the default region is poor for SASA).
  - **7 Ramachandran**: type 1 φ/ψ, 2 χ1, 3 ω; residue range (all); α and β region boxes are editable → 17 residues (23%) α, 23 (31%) β.
  - **8 B-factors**: B = 8π²/3 · RMSF², mean 2.22 Å², interpretation ranges shown.
  - **13 PCA**: CA, 3 components, fit (1 to average / 2 none) → PC1 54.7%, first 3 = 83.2%; scree, projections, PC1–PC2 plot.
  - **14 Clustering**: CA; 1 K-means / 2 hierarchical; 5 clusters; metric 1 best-fit RMSD / 2 no-fit / 3 distance RMSD; seed 1 → 5 clusters (25%–15%).
  - **16 Pairwise RMSD**: subsampling 1 none … 5 custom; metric → 190 pairs, mean 0.40 Å.
  - **17 H-bonds**: donor/acceptor masks (blank = all), 3.0 Å, 135° → 77 unique pairs, 36 ± 4 per frame; top list, e.g. ASP52 OD2 – LYS27 NZ 100%.
  - **15 Autocorrelation** (of an RMSD series): region, reference frame, max lag (≤ half the frames).
  - **18 Water RDF**: "Enter solute mask (e.g. :42-45)" — I left it blank and it was **accepted silently**, giving a "Solute:" (empty) RDF to `:WAT@O`, first peak 2.75 Å. Type a mask such as `:1-76`.
  - **19 Water shells**: solute 1 protein / 2 custom; shell width 2.0 Å; max 10 Å → 5 shells (cpptraj watershell, 602 solute heavy atoms, 4,824 waters).
  - **22 Contact frequency**: region and cutoff prompts.
  - **11 Dihedral**: `:10@N :10@CA :10@CB :10@CG` on Gly10 (no CB) → **pytraj "Failed to setup action" plus a full Python traceback, and you are thrown out to the MD Manager menu** (the loaded trajectory is lost). With `:4@N :4@CA :4@CB :4@CG` (Phe4 χ1) it worked: mean −51.7°.
  - **12 Vector**: start `:23-34@CA`, end `:23@CA` → magnitude 8.43 Å, polar angle 134.9°.
- `0` → "Returning to MD Manager menu…".

### Scenario 14: Analyses 10 (Angle) and 20 (Density maps) cannot be chosen — FAILED (bug)
- Typing `15,18,19,20,22` immediately printed "Returning to MD Manager menu…".
- Source (`molecular_dynamics_manager.py`, `_get_analysis_selection_menu`): `if '0' in selection: return None`. This is a substring test, so **any answer containing a 0 (10, 20) exits the menu**. Angle (10) and Density maps (20) can never run, alone or in a list.

### Scenario 15: Trajectory file browser limits — RUN
- From a step folder, `..` then `find` lists the `.nc` files in all step folders. Its prompt "Select file (1-2) or 'cancel'" accepts **one** file only; `2,1` was rejected and dropped back to the browser.
- Multi-select (`N,M`) only works inside one folder. To join heating and production trajectories, use `proprep --analysis` (Scenario 22), whose browser can stage files across folders.

### Scenario 16: Custom protocol from the catalog, constant-pH production, extended cycles, SLURM scripts (manual wizard) — RUN (with a SLURM bug)
- Setup → protocol source `3` → "Template Catalog": all steps of all 5 built-in protocols, numbered 1–35. Commands: `<n>` / `<n>,<n>` add, `p <n>` preview, `c <n> <n>` compare, `r` remove last, `d` done, `b`.
  - **User-imported templates are not listed here.**
- `18,21` (Workshop minimization and production). `c 18 22` → side-by-side parameter table (maxcyc 1000 vs 500, restraint_wt 10 vs 50, cut 10 vs 12…). `d` → name `s7_custom_test`.
- Constant pH `y`: "Target pH (7.0)", "ntcnstph (100)", "ntrelax (200)", a note on saltcon, "saltcon (0.1)" → "✓ Constant pH MD enabled icnstph=2, solvph=7.0, ntcnstph=100, ntrelax=200, saltcon=0.1 … applied to production steps only".
- Restraints `s`; queue toggles (Scenario 7); hardware `3` → sander; review `4`:
  - "✓ Staged CPIN into step2/ (1UBQ.cpin) ✓ Staged modified topology into step2/".
  - "Add extended production run cycles? (n)" → `y` → "Number of additional production cycles (10)" → `2`.
  - SLURM wizard, all defaults: job name `amber_simulation`, partition `gpu`, nodes 1, tasks/node 1, CPUs/task 4, "Use GPU? (y)", GPUs 1, GPU type blank, memory 16G, time 02:00:00, email n, "Load modules? (y)", modules `amber/22`, conda n, advanced n → summary → "Use this configuration? (y)".
- Files:
  - `simulations/batch_20261121_132159/` (local layout; `run_workflow.sh` has an extended-production loop for 2 cycles chaining `.rst7` and `.cprestrt`).
  - `simulations/slurm_scripts/{step1.sh, step2.sh, submit_workflow.sh, slurm_config.json, step1/, step2/1UBQ.cpin, step2/1UBQ_cpin.prmtop}`.
  - `submit_workflow.sh` uses `sbatch --parsable --dependency=afterok:$STEP1_JOBID`.
- **Bugs and oddities in the SLURM output:**
  - `slurm_scripts/step2.sh` runs `sander … -p ../1UBQ.prmtop … -cpin 1UBQ.cpin`, i.e. the **non-cpin topology** with a cpin file. The local `run_workflow.sh` correctly uses `1UBQ_cpin.prmtop`.
  - It requests `--gres=gpu:1` on partition `gpu` for the CPU-only engine `sander`.
  - `step2.sh` also contains a literal `#SBATCH --dependency=afterok:$STEP1_JOBID`, which is invalid if step2.sh is submitted by itself.
- Option 5 (sbatch now) — SOURCE-ONLY: needs a loaded cluster profile and a SLURM system.

### Scenario 17: Templates & Protocols library — RUN
- 5b → 5. Import: 1 .mdin template, 2 .json protocol. Templates: 3 statistics, 4 preview, 5 delete user templates. Protocols: 6 statistics, 7 delete custom. b back.
- `3` → 35 built-in templates (Minimization 8, Heating 5, Equilibration 10, Production 4, Other 8). Library locations: base = the package `md_templates`; User = `~/.proprep/md_templates/user`.
- `4` → a numbered list of all 35 templates with descriptions; `0` cancels.
- `1` → "No .mdin files found in current directory. Enter path to .mdin file:" → `import_test/s7_quick_min.mdin` → name (s7_quick_min) → description → "Detected simulation type: minimization ✓ … Saved to: **<project>/md_templates/user/custom/s7_quick_min.mdin**". That is a different place from the "User" location the statistics screen reports.
- `6` → "Custom Protocols: Protocol_1UBQ_1UBQ_20261121_1243 (**? steps**)".
- `7` → pick `1,2` → "Associated templates to delete (orphaned): 3" → "Delete 2 protocols and 3 orphaned templates? (n)" `y` → "✓ Successfully deleted 2 protocols". Statistics then show 0 custom, but **the JSON files were still in `~/.proprep/md_workflows/user/custom/`.** I removed my own leftover files by hand so other users of `~/.proprep` are not affected.

### Scenario 18: "Load a saved protocol" — RUN (always empty)
- Protocol source `2` → "No saved protocols found. 1. Import protocol from file (JSON)". This happened even while the library (Scenario 17) listed a custom protocol. Saved protocols do not show up here.

### Scenario 19: Cluster profiles (6) and run plans (7) menus — RUN (menus only)
- 6: lists "TAMU ACES tamu-aces (bundled)". Actions: 1 show details, 2 fill required fields (quick), 3 edit full JSON in $EDITOR, 4 create from template, 5 export, 6 import, 7 delete user profile.
  - `4` → "Name for new profile (e.g. 'my-cluster')" → empty → "Invalid name (cannot be empty or start with '_')".
  - The template (`md_prep/cluster_profiles/_template.json`) documents: scheduler (slurm only), account, qos, email, modules, conda_env, pre/post commands, extra #SBATCH directives, resource_classes (mode cpu/gpu, partition, gres, memory, binary, mpi_launcher — use `srun`, not `mpirun`), and defaults.
- 7: "Available run plans: (none — save one from Step 4)". In fact run plans are saved in **Step 5 (Hardware)**. Actions: 1 show, 2 export, 3 import, 4 delete, 5 remap to another cluster.

### Scenario 20: Engine not installed (pmemd) — FAILED (expected)
- Session `s7-pair`, loaded pair (Scenario 5a-14). Setup → predefined 4 → `s` (select steps) `1` → name `s7_minonly` → restraints s → queue n → hardware: "Current engine: **pmemd**" (the default!) → n → review `3` "Write files and run locally now" → "✓ Configuration saved, switching to execution…" → 1 → y.
- Screen: "Error starting simulation: [Errno 2] No such file or directory: 'pmemd'", "✗ Failed to start", "Simulations started: 0".
- The execute menu's "Manage queue" cannot change the engine; you must go back through Setup → Hardware.
- **Beginners will hit this: the default engine is pmemd, which is not part of the free install.**

### Scenario 21: "Create a step from scratch" (AMBER configuration wizard), then run it with sander — RUN
- Protocol source `4` → overview of 13 sections → "Educational Content Level: 1 Terse / **2 Compact** / 3 Verbose" → `1` → "View parameter reference guide? (n)".
- Wizard prompts seen (all defaults except the type):
  - Step 1: type 1 Minimization (imin=1), 2 MD, 3 Trajectory analysis (imin=5), 4 MD analysis (imin=6), 5 Socket server (imin=7), 6 Custom → `1`.
  - 2 NRESPA? (n); 3 boundary conditions; 4 input/restart.
  - 5 Minimization: maxcyc (10000), method (1), ncyc (1000), drms (0.0001), step size.
  - 9 **Pressure control** (asked even for a minimization): pres0 1.0, barostat 2 (Monte Carlo), mcbarint 100, taup 2.0, compressibility 44.6.
  - 10: cutoff 10 Å, force switching n, dielectric n, nsnb 0, PME y, manual PME n, external field n.
  - 11: ntpr 100, ntwr 1000, wrap n, file formats.
  - 12: varying conditions, PMEMD-specific n, energy decomposition n, belly n, reference guide n.
  - Step name `Custom_minimization_1UBQ`, description, protocol name, "Add more steps? (y)".
- The generated parameters include **ntb=2, ntp=1, barostat=2 (constant pressure) for an energy minimization**. sander accepted it, warning "mcbarint is greater than the number of steps. This is effectively constant volume."
- I edited maxcyc to 200 (`edit 1` → params), saved, and set hardware `2` → `1` sander. Run locally → finished in 11.5 s, final energy −5.6566E+04, RMS 0.98.
- Files: `runs/s7/pair/pair/simulations/batch_20261121_135635/01_1UBQ_Custom_minimization_1UBQ/`.

### Scenario 22: Not run — SOURCE-ONLY
- **Analysis option 4, "View trajectory in the structure viewer"**: not run, because it opens a web page in the browser (utility `ua`, another scope).
- **Plugins**: ProPrep discovers external packages through the `proprep.plugins` entry-point group and adds them to the menu (`plugins/`, `menu_commands._dispatch_plugin_in_menu`). None are installed in the ProPrep env (`entry_points(group='proprep.plugins')` is empty), so no plugin items appear. The analysis menu has no plugin hook; its 23 analyses are built in.
- **sbatch submission (review option 5)** and **"Save setup to workspace (execute later)" (option 1)**: not exercised.

### Options reference (5b)

| Where | Choices |
|---|---|
| Built-in protocols | **Membrane Protein Minimization** (2: restrained min, unrestrained min); **Basic Equilibration** (6: min 10 kcal backbone+metal; heat 0→300 K 2 ns; NPT density 1 ns CPU; NPT relax 3 ns GPU; NVT eq 4 ns; NVT production 100 ns); **Staged Equilibration 310 K** (9: SD+CG min 50 kcal; heat 0→310 K; NPT density CPU; NPT relax GPU; release 10/5/2 kcal; unrestrained settle 1 ns; production 50 ns NPT); **Workshop CPU** (4: min 1000; heat 20 ps; NPT 100 ps; NVT 300 ps); **York Lab tutorial** (14: solvent min, 3× NPT, solvent min, heat, NPT solvent eq 5 ns, solute min, solute heat 25 kcal, eq 25/10/5/2 kcal, production 100 ns) |
| Protocol commands | `p <n>`, `u` use as-is, `s` select steps, `e` edit, `b` |
| Workflow editor | `edit/name/files <n>`, `add before/after <n>`, `remove <n>`, `move a to b`, `preview`, `info` (edits name/description), `save`, `cancel`, `help`; per step: `params` (editor), `name`, `desc`, `files`, `back` |
| Catalog builder | `<n>` / `<n>,<n>` add, `p <n>`, `c <a> <b>` compare, `r` remove last, `d` done, `b` |
| Constant pH (if a cpin exists) | enable y/n; Target pH 7.0; ntcnstph 100; ntrelax 200; saltcon 0.1 M |
| Restraints | (e) edit masks/weights per step (`none` clears); (g) generator: (b) backbone, (h) heavy atoms, (c) CA, (u) custom × chains/ranges (`A:10-50`) × exclude water; (d) DISANG: import / interactive (s, n, r, m modes; 2 residues = distance, 3 = angle, 4 = torsion) / display / export, dump frequency; (r) GROUP: 1–10 groups, k, 1–5 FIND criteria, ≤ 7 residue ranges |
| Queue | k keep, r reorder, t toggle (`hold`/`active` n, lists, ranges, `all`) |
| Hardware | 1 recommended, 2 per simulation (engine + MPI count / GPU id), 3 bulk engine, p cluster profile, l load run plan, s save run plan, h help. Engines: sander, sander.MPI, pmemd, pmemd.MPI, pmemd.cuda |
| Review | 1 save to workspace, 2 MDIN only, 3 write + run locally, 4 MDIN + SLURM (manual wizard), 5 write + sbatch (needs profile); then extended production cycles y/n and count |
| SLURM wizard | job name, partition, nodes, tasks/node, CPUs/task, GPU y/n, count, type, memory, time, email, modules, conda env, advanced |
| Execute menu | 1 execute queue, 2 manage queue (1 protocols, 2 steps, 3 active/hold, 4 clear queue), 3 monitor, 4 analyze |
| Monitor | 1 status, 2 energy, 3 temperature, 4 pressure, 5 timing, 6–11 ASCII plots (Etot, T, P, KE, PE, RMS gradient), 12 periodic image distance, r refresh |
| Analysis type | 1 energetics, 2 trajectory, 3 combined, 4 structure viewer; type override 1 Min, 2 Heat NVT, 3 Heat NPT, 4 NVT eq/prod, 5 NPT eq/prod |
| Trajectory analyses | 1 RMSD, 2 RMSF, 3 Contacts, 4 Salt bridges, 5 SASA, 6 DSSP, 7 Ramachandran, 8 B-factors, 9 Distance, 10 Angle (unreachable), 11 Dihedral, 12 Vector, 13 PCA, 14 Clustering, 15 Autocorrelation, 16 Pairwise RMSD, 17 H-bonds, 18 Water RDF, 19 Water shells, 20 Density maps (unreachable), 21 Rg, 22 Contact frequency, 23 Periodic image distance, 0 exit |
| Export | 1 CSV, 2 JSON, 3 both, 4 none (default) → `<project>/<simname>_trajectory_analysis.csv/json` |

---

## Feature: `proprep --analysis` (command-line shortcut)

**What it is for:** to open the MD analysis browser directly, without going through the menus, for
example to look at a simulation you ran earlier.

**Needs before you start:** run it inside a folder that contains a `simulations/` tree (or browse to
the files). It cannot be combined with `--pdbid`/`--pdbfile`.

### Scenario 1: Analyse two joined trajectories from the command line — RUN
- `cd runs/s7/main/ubq; proprep --analysis`.
  - Asked "Select menu mode [1/2]", then **"Save this as your default preference? (n)"**. Answer `n`, or add `--full-menu` to skip both.
  - "Using current directory …", "Jumping to Simulation Analysis…", then the list of 3 completed simulations.
- `browse` → a richer browser: `select N`, `select N,M`, `add N` (stage files and keep browsing other folders), `done`, `clear`, `cd N`/`cd path`/`cd ..`, `find`, `q`.
- `find` → 16 files (mdouts, `.nc`, many prmtops, including the 0-byte `1UBQ.prmtop` from Scenario 5a-7) → "Select file(s) (comma-separated)" → `5,4,14` (heating .nc, production .nc, `shared/1UBQ.prmtop`).
  - "Loading 2 Trajectory Segments… Concatenating… ✓ Loaded 40 total frames" (0.1–3.0 ps).
- `21` Rg (region 4) → 11.84 ± 0.06 Å → export 4 → `0`. This lands in the **ProPrep main menu**, not back in the analysis browser.
- Session log: `runs/_sessions/s7-ana/transcript.txt`

---

## Notes for the manual writer

**Jargon to explain:** topology (prmtop/parm7) vs coordinates (rst7/inpcrd); force field; water model
(3-point vs 4-point); truncated octahedron; buffer; neutralizing counter-ions vs bulk salt; 12-6-4
ion model; SPLIT ion counting; tLEaP; ParmEd; GB/implicit solvent; disulfide/CYX; CONECT/TER records;
minimization / heating / NVT / NPT / production; restraint mask (Amber mask syntax `:res@atom`);
DISANG / NMR restraints; GROUP restraints; sander vs pmemd vs pmemd.cuda; MPI; SLURM; constant-pH MD
(cpin, AS4/GL4/HIP); RMSD/RMSF/Rg/DSSP/PCA.

**Recommended beginner path that worked end to end** (1UBQ, about 2 minutes of compute):
1. 1a load.
2. (Optional: remove waters with 4a; run 4e.)
3. 5a → 2. Accept defaults (ff19SB + OPC is the recommended pair; I used TIP3P), but type `1` for unneeded DNA/RNA/sugar/lipid force fields. Edit template `n`. **Answer `y` to "continue without bond definitions".**
4. 5a → 4.
5. 5b → 1 → n → 1 → Workshop protocol → `e` to shorten it.
6. **Hardware: set every step to `sander` (or `sander.MPI`).**
7. Review 3 (run locally) → 5b → 3 to monitor → 5b → 4 to analyse.

**Things that surprised me or look like bugs (described, not fixed):**
1. The 5a default answer `n` to "continue without bond definitions?" aborts the whole input generation.
2. Failed tLEaP runs truncate the target prmtop to 0 bytes and still print "✓ Topology files generated". Two input files can share one output name with no warning. The real tLEaP FATAL reason is not shown on screen.
3. Implicit solvent is offered without removing crystal waters (fatal) and without setting GB radii.
4. Disulfide proteins fail unless the cysteines are CYX **and** bond lines are present. 4e does not rename them. CONECT records are silently "Illegal" after transformation, so a CYX build can "succeed" without disulfides; only a yellow ParmEd line warns.
5. 5a option 3 (microstate batch) silently writes nothing when no site needs force-field files (disulfide-only).
6. Analyses 10 (Angle) and 20 (Density maps) can never be selected: any answer containing `0` exits.
7. DISANG manual entry cannot find residues (no chain IDs in the topology-extracted PDB). DISANG import prints an AttributeError but works.
8. The default MD engine is pmemd (not installed). "Recommended" picks pmemd.cuda with 0 GPUs. Only the per-simulation list marks engines "(not installed)".
9. The SLURM step for constant-pH production uses the non-cpin prmtop and asks for a GPU for sander.
10. The review table shows the template's original step lengths, not edited ones. Preview shows `pmemd` commands and `step_N.rst` names that differ from the real run.
11. "Load a saved protocol" never lists saved protocols. Deleting custom protocols leaves their JSON files on disk. Imported templates are saved in the project but the library reports `~/.proprep`. "? steps" is shown for every custom protocol.
12. The wizard's default minimization contains constant-pressure settings (ntb=2, ntp=1).
13. Heating is auto-detected as "NVT production". The production metrics table is empty.
14. An invalid atom mask in the dihedral analysis throws a traceback and exits the analysis (loaded trajectory lost). An empty RDF solute mask is accepted.
15. cpin generation on a non-constant-pH topology is allowed with no warning, and silently switches 5b's default structure pair to the `_cpin.prmtop`.
16. Raw Python dictionaries are shown as "Previous settings". INFO/WARNING log lines leak into the screen (trajectory_analyzer, mpsa_workspace).
17. The TER-record fixer edits the user's PDB in place (a backup is kept).
18. With another `AMBERHOME` set, ProPrep silently uses that Amber installation.

**Driver note (for anyone re-running):** several prompts (`Enter …:` with no default) are only
printed after the next key arrives, so a screen can look stuck after an answer. Pressing Enter or
looking at the transcript reveals the waiting prompt.
