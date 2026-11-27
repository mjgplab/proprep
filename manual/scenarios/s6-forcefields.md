# Force fields and parameterization (scope s6)

Recorded against ProPrep 1.23.0, `--full-menu`, 2026-11-21. All sessions ran with an
isolated home (`HOME=manual/runs/s6/home`), so the user library written here is
`manual/runs/s6/home/.proprep/...`, not the shared `~/.proprep`.

**Important environment note (affects every scenario below).** The shell that launched
the first sessions had `AMBERHOME=/opt/homebrew/Caskroom/miniforge/base/envs/etanalyze_env`
(a different AmberTools install). ProPrep uses whatever `AMBERHOME` is set to: the
Force Field Explorer read *that* install's force-field files, and the metal workflow
ran `.../etanalyze_env/bin/teLeap`. From session `s6-metal-3` onwards I started ProPrep with
`AMBERHOME=<ProPrep env>` and `PATH=<ProPrep env>/bin:$PATH`. That is what `conda activate ProPrep`
normally sets up. The force-field lists differ slightly between the two installs (see 2b,
Scenario 1). **Manual writer: tell users to always `conda activate ProPrep` and not to
set `AMBERHOME` by hand to some other Amber installation.**

---

## Feature: 2b. Force Field Explorer

**What it is for (plain English):** A read-only browser for the AMBER force-field
files on your computer. It shows the atom types, residue templates (atoms, types and
charges), and bond/angle/dihedral/van-der-Waals numbers in a force field. It also
shows ProPrep's built-in cofactor parameter sets and your own deposited parameters. It
can also check every atom in your loaded structure against the chosen force field and
report which atoms the force field does not cover. Use it to answer "does AMBER already
know this residue?" before you start parameterizing.

**Needs before you start:** Nothing (AmberTools must be installed, `AMBERHOME` set).
Option 7, "Type Structure Atoms", also uses a loaded structure. Without one it opens a file browser.

**Where to find it:** Main menu `2b` (section "2. REDOX SITE DETECTION"). In Workflow mode it is
one of the "REDOX SITE DETECTION" tools ("Detect redox-active sites, explore the parameter
space of AMBER FFs ..."). No command-line flag.

### Scenario 1: Browse every force-field category   — RUN
- Starting point: fresh project `ffe` (runs/s6/explorer), nothing loaded.
- Steps:
  1. Prompt "Enter your choice (1a):" -> typed `2b`.
  2. Screen "Step 1: Select Force Field Categories" (table Index / Category / Available Files):
     1 Protein (9), 2 Modified Amino Acids (8), 3 DNA (5), 4 RNA (6), 5 Water (12),
     6 GAFF (General) (2), 7 Lipid (1), 8 Carbohydrate (2), 9 Dye (1), 10 Constant pH/Redox (2),
     11 Other (4; 6 with the ProPrep AmberTools), 12 Redox Cofactors (built-in) (69), 13 Custom/User-Generated ("User files").
     Hints printed: "Protein simulations: 1,5", "Protein with ligands: 1,5,6", "DNA simulations: 3,5".
  3. Prompt "Select categories (comma-separated) (1,5):" -> typed `1,2,3,4,5,6,7,8,9,10,11`.
  4. Screen "Step 2: Select Force Fields" lists every leaprc with a friendly name,
     description, ★ = recommended, ⚠ = add-on that needs a base protein FF. Full list seen:
     - Protein: fb15, ff03.r1, protein.ff03ua, ff14SB, protein.ff14SBonlysc, ff15ipq
       ("REQUIRES SPC/Eb water"), protein.ff15ipq-vac, ★ff19SB ("use with OPC water"), protein.ff19ipq.
       (With the ProPrep env AmberTools the list has protein.ff12pol and no protein.ff19ipq.)
     - Modified Amino Acids (⚠ add-ons): Fluorinated AA (ff15ipq), Mimetic residues (ff15ipq),
       Phosphorylated AA (ff99SB / ff14SB / ff19SB), phosfb18, Modified AA (ff14SB),
       Modified AA (ff19SB) "Selenomethionine, MTSL spin labels, etc."
     - DNA: OL15, OL21, ★OL24, bsc1, tumuc1. RNA: LJbb, ★OL3, ROC, Shaw, YIL, modrna08.
     - Water: fb3, fb3mod, fb4, ★OPC, OPC3, OPC3-pol, SPC/E, SPC/Eb, TIP3P ("traditional (not
       recommended)"), tip4pd, tip4pd-a99SBdisp, TIP4P-Ew.
     - GAFF: GAFF ("less accurate"), ★GAFF2. Lipid: ★lipid21. Carbohydrate: GLYCAM_06EPb, ★GLYCAM_06j.
       Dye: amberdyes. Constant pH/Redox: Constant Redox (ff10), Constant pH (ff10). Other: ffAM1, ffPM3, music, xFPchromophores.
  5. Prompt "Select force field(s) (comma-separated) (4,37,42):" -> typed `8,32,42` (ff19SB + OPC + GAFF2).
- What the user sees: "Loading: leaprc.protein.ff19SB, leaprc.water.opc, leaprc.gaff2" …
  "Loaded: 1536 atoms, 182 residues, 1524 bonds, 10272 angles, 1748 dihedrals, 1358 VDW params", then the Explorer menu
  (1 Overview, 2 Browse Atom Types, 3 Browse Residues, 4 Search Parameters, 5 Reverse Lookup,
  6 Load Different Force Field, 7 Type Structure Atoms, b Back).
- Problems: the suggested default `(4,37,42)` = ff14SB + **TIP3P** + GAFF2. That is not the ★-recommended
  ff19SB + OPC, and TIP3P is labelled "not recommended" on the same screen. Typing a number out of
  range or text is silently ignored (only valid numbers are kept).
- Session log: runs/_sessions/s6-ffe-1/transcript.txt; session file runs/s6/explorer/ffe/proprep_session_20261121_122935.json

### Scenario 2: Overview of what was loaded   — RUN
- Steps: "Enter your choice (1):" -> `1`.
- What the user sees: table "Force Field Statistics" (Atom definitions 1536, Residues 182, Residue aliases 5,
  Mass 227, Bond 1524, Angle 10272, Dihedral 1748, VDW 1358, Improper 3). Then "Loaded Files": Leaprc
  (3), Lib/Off (amino19.lib, aminoct12.lib, aminont12.lib, atomic_ions.lib, solvents.lib), Parm.dat
  (parm19.dat, gaff2.dat), Frcmod (frcmod.ff19SB, frcmod.opc, frcmod.ionslm_126_opc).
- Problems: "Improper parameters 3" for ff19SB+GAFF2 looks too low (possible parsing gap; describe, not verified).

### Scenario 3: Browse atom types (filter and full list)   — RUN (bug found)
- Steps: `2` -> prompt "Enter atom type filter (e.g., 'CT'), Enter for all, 'b' to go back ():" -> `CX`.
  Result: 2 types, CX "protein C-alpha (new to ff10)" and lowercase GAFF `cx` "Sp3 carbons in triangle systems".
  The filter is case-insensitive and matches any part of the name.
- Then Enter (all): "Atom Types (1362 total) — Page 1/69", navigation "[n]ext [q]uit".
- **Bug:** the first pages are filled with junk "types" such as `-0.00006` with descriptions like
  "5.13593 6.22093 6.57327 …". These are numbers from the ff19SB CMAP grid, read in as atom
  types. It happens with both AmberTools installs (1270 "types" with ff19SB+phosaa19SB+OPC). Real
  types are only reachable by filtering. Advise users to type a filter.

### Scenario 4: Browse residues and view one template   — RUN
- Steps: `3` -> "Enter residue filter (e.g., 'ALA'), Enter for all, 'b' to go back ():" -> `HI`.
  Table of 12 matches (CHID/CHIE/CHIP from aminoct12.lib, HID/HIE/HIP from amino19.lib, NHID…; HIS shown as
  "alias:HIE"). "[#] select [q]uit" -> `5` shows "Residue: HID (17 atoms)" with Atom Name / Type / Charge / Source
  (e.g. CA XC 0.0188, ND1 NA -0.3811, HD1 H 0.3649).
- Gotcha observed: after a residue detail is shown you are back at the *filter* prompt. Typing `q` there
  does not leave. It searches for residues containing "q" (got QSPCFWBOX). Only `b` leaves.
  I accidentally typed `4` and `a` as filters too.

### Scenario 5: Search parameters (bond/angle/dihedral/VDW/mass, plus bad input)   — RUN
- Steps: `4` -> menu a Bond (2 types), b Angle (3), c Dihedral (4), d VDW (1), e Mass (1), q Back.
  - `a` then "Enter 2 atom types (e.g., 'CT N' or 'CT-N')" -> `XC N`: "Bond: XC-N  Force Const. 337.0  Eq Length 1.4490  parm19.dat".
  - `b` -> `N XC C`: angle C-XC-N, 63.0, 110.10°, parm19.dat (order of the outer types does not matter).
  - `c` -> `X XC C X`: matched X-C-XC-X (wildcard, barrier 0.000). `c` -> `C N XC C`: 4 Fourier terms, all
    barrier 0.000 from frcmod.ff19SB (this is expected for ff19SB, where CMAP carries the backbone torsions).
  - `d` -> `ca`: GAFF2 VDW radius 1.8606, depth 0.0988.
  - `e` -> `ZZ`: "No mass parameter found for 'ZZ'".
  - `a` -> `foo`: "Expected 2 atom types, got 1".
- Case matters: `CA` (protein) and `ca` (GAFF aromatic) are different types.

### Scenario 6: Reverse lookup of an atom type   — RUN
- Steps: `5` -> "Enter atom type to search for (e.g., 'CT')" -> `XC`.
  Shows Mass 12.010 and VDW, then sub-menu "1. Bonds (9) 2. Angles (66) 3. Dihedrals (92) b. Back" -> `1` lists
  all 9 bonds (C-XC, C8-XC, XC-2C, … with source file). `b` returns.

### Scenario 7: Built-in redox cofactor parameter sets   — RUN
- Steps: `6` (Load Different Force Field) -> category `12`.
- What the user sees: panel "Built-in Redox Cofactor Parameters … Bundled .frcmod / .lib sets for hemes,
  iron-sulfur clusters, and other redox-active cofactors". Then a 69-row table: Cofactor / Redox / Spin-Variant / Set / Notes.
  The families are:
  - flavin/fad and flavin/fmn (5 redox states each: oxidized, semiquinone anionic/neutral, hydroquinone anionic/neutral),
  - heme/bis_his_b_type (ox/red, low spin, CM5 or RESP), heme/bis_his_c_type (Henriques sets, fixed pH / constant pH / constant E),
    heme/cys_axial_b_type (high spin), heme/his_met_axial_c_type (CM5, CM5_FeLocal, RESPmc, RESPsc; fixed or constant pH),
  - iron_sulfur/4fe4s (oxidized FO0–FO6, reduced_fd FD0–FD6, reduced_hipip FR0–FR3),
  - nicotinamide/nad and nadp (ox/red), pterin/biopterin (6 states), zinc/cys4 (Guberman_ZNC_RESP).
  The Notes column is very long and wraps over 15+ lines per row. The table is about 750 lines of terminal output.
- Prompt "Select set(s) (comma/range or 'none') (none):" -> `12,69` (bis-His b-heme oxidized RESP + Zn(Cys)4).
  "✓ Loaded 89 atom definitions from Guberman_HBO_RESP.lib … ✓ Loaded 11 atom definitions from Guberman_ZNC_RESP.lib".
- Then Browse Residues shows CYZ, FHO, HBO, RHO, ZNC. Residue CYZ shows SG type SZ charge -0.8185, ZNC = ZN type ZM charge 0.9617.
  Atom types: M1 (Fe), SZ ("S(γ) of Cys-thiolate bound to Zn — forked from ff14SB SH"), Y1–Y6, ZM ("Zn metal (Li et al. JCTC 2013 …)").
  Search bond `ZM SZ` -> 32.3 kcal/mol/Å², 2.4310 Å. Angle `SZ ZM SZ` -> 26.9, 109.33°.
- Problems: the header line "Loaded:" is empty when only cofactor sets are loaded (cosmetic).
  `SZ ZN` / `ZN` find nothing because the metal *type* is ZM, while ZN is the atom *name*. Users must search by type.

### Scenario 8: Custom/User-Generated files (load your own .lib/.frcmod)   — RUN
- Steps: `6` -> `13`. Panel "Custom Force Field File Loader … Supported formats: .lib, .off, .frcmod, .mol2, .prep".
  A file browser opens in the project folder. Commands: N, N,M / N-M (multi-select), `..`, `q`. This browser has no "<path>" command.
  First time: `q` (cancel) -> back to main menu.
  Second time (session s6-maa-1): `..`, `..`, `2` (bcc/), `2-3` selected mk1_bcc.frcmod + mk1_bcc.lib ->
  "✓ Loaded parameters from mk1_bcc.frcmod ✓ Loaded 92 atom definitions from mk1_bcc.lib" -> "Load another file? (y)" -> `n`.
- Result: residue MK1 browsable (C1 type c6 charge 0.1648, …).

### Scenario 9: Type Structure Atoms (atom-typing report)   — RUN
- 9a, no structure loaded (s6-ffe-1): `7` opens a file browser ("Current directory: … [1] .. (parent)"). Cancelled with `q`.
- 9b, structure 3AJ4 loaded, FF = ff19SB + Phosphorylated AA (ff19SB) + OPC (s6-maa-1):
  `7` -> "Using workspace structure: … Use this structure? [y/n] (y)" -> `y` (it used the H-stripped 3AJ4_noH.pdb).
  Report: "Total atoms 2155, Typed 2142 (99.4%), Untyped 13 (0.6%), Unique atom types 27, Total charge (typed) -234.2566".
  "Untyped residues: A:EDO112 (4 atoms), A:SEP111 (1 atoms), B:EDO111, B:EDO112". There is a note explaining that per-residue charge
  sums are non-integer for heavy-atom-only structures. Paginated "Per-Residue Charge Summary" -> `q`.
  Then "View [a]ll atoms, [t]yped only, [u]ntyped only, [e]xport report, [d]one" -> `u` lists the 13 untyped atoms with coordinates
  (SEP OXT, all EDO atoms). `e` -> "Export path (3AJ4_noH_atom_typing.txt)" -> Enter -> "Report exported to: 3AJ4_noH_atom_typing.txt". `d` done.
- Result: runs/s6/modaa/maa/3AJ4_noH_atom_typing.txt (text report).
- Useful finding for the manual: phosphoserine (SEP) **is already typed** by the "Phosphorylated AA (ff19SB)"
  add-on. Only its terminal OXT was untyped, because here SEP is a free amino acid. So the user did not
  need to parameterize SEP. The Explorer is the place to check this first.
- 9c (s6-sm-1, 1HSG with ff19SB+GAFF2+MK1 library set): the session's control process died while at
  "Use this structure?". Not repeated with 1HSG (9b covers the feature).

### Scenario 10: Your deposited parameters appear in the Explorer   — RUN
- After the library import (see "Local parameter library" below), category 12's count rose from 69 to 71.
  Two new rows were listed: "69 small_molecules/MK1 single_state default MK1_AM1BCC Indinavir (MK1) from 1HSG, GAFF2 + AM1-BCC charges"
  and "70 … MK1_AM1BCC_v2 default — …". zinc/cys4 moved to row 71.
- Selecting `69`: "✓ Loaded parameters from mk1_bcc.frcmod ✓ Loaded 92 atom definitions from mk1_bcc.lib Loaded 2 bundled file(s)".
  Browse Residues -> MK1 (92 atoms, mk1_bcc.lib).
- Confusing: user entries are listed under a category called "Redox Cofactors (built-in)" and the message says "bundled",
  although they are neither redox cofactors nor built in. Row numbers of the real built-in sets also shift.

### Options reference (2b)
| Prompt | Options | Meaning |
|---|---|---|
| Select categories | numbers, comma-separated (default `1,5`) | 1–11 = AmberTools leaprc families; 12 = ProPrep bundled cofactor sets **plus your user library**; 13 = browse for your own .lib/.off/.frcmod/.mol2/.prep |
| Select force field(s) | numbers (default ff14SB+TIP3P+GAFF2 when those categories are chosen) | ★ recommended, ⚠ add-on that must be combined with a protein FF |
| Select set(s) (cat. 12) | numbers, ranges `2-4`, `none` | loads matching .frcmod + .lib |
| Explorer menu | 1 Overview, 2 Atom types, 3 Residues, 4 Search, 5 Reverse lookup, 6 Load different FF, 7 Type structure atoms, b Back | |
| Atom/residue filter | text (substring, case-insensitive), Enter = all, `b` = back | `q` is treated as a filter, not "quit" |
| Paginated tables | `n` next, `q` quit, `#` select (residues) | |
| Search type | a bond, b angle, c dihedral, d VDW, e mass, q back | types separated by space or `-`; `X` = wildcard in dihedrals |
| Typing report view | a all, t typed, u untyped, e export, d done | export writes `<pdb>_atom_typing.txt` |

---

## Feature: 2c. Force Field Parameterizer

**What it is for (plain English):** Makes AMBER parameters for anything in your structure that the standard protein
force field does not know: a drug or ligand ("small molecule"), a chemically modified amino acid, or a metal ion
together with the residues holding it ("metal site"). It finds these residues, sorts them into the three kinds, and
runs a step-by-step checklist for each. The checklist writes the quantum-chemistry (Gaussian) input files and
then turns Gaussian's results into charges and force constants.

**Needs before you start:** A loaded structure (menu 1a). **Gaussian 16 (`g16`, `formchk`) for every route that
derives charges.** ProPrep only writes the Gaussian inputs. The user must run them, then come back and resume.
For metal sites, run the Redox Site Detector (2a) first. ProPrep asks for this itself in step 11 if you have not.
AmberTools (antechamber, parmchk2, tleap, cpptraj, reduce, propka) is used and is installed.
Optional refinements need CREST/xtb (not installed) or Gaussian PES scans.

**Where to find it:** Main menu `2c` (shows "⚠ Needs a loaded structure" until one is loaded). Workflow mode:
"REDOX SITE DETECTION" stage. Sub-menu:
1 Analyze and configure residue classifications · 2 Import existing parameters into your library ·
3 Parameterize residues (new or resume) · 4 View parameterization status · 5 Display detailed help information.

### Scenario 1: Help screen   — RUN
- `5` shows "Forcefield Parameterization Help". It describes the 8 small-molecule steps (two-step Gaussian:
  B3LYP/6-31+G(d) opt+freq with IOp(7/33=1), then HF/6-31G(d) Merz-Kollman ESP via --link1--). It also covers the refinement options
  (Seminario; PES scan of 24 points at 15°; CREST), the modified-AA 10 steps (ACE-XXX-NME, α-helix φ/ψ=-60/-45 and β-sheet -135/135, HF/6-31G*,
  residuegen multi-conformer RESP) and output files.
- Problem: the Metal Site section is **out of date**. It says "Step 2: Metal center parameter building 🚧 Coming soon"
  and so on, but the real metal checklist has 15 steps through Force Field Integration (see below).

### Scenario 2: Analyze non-standard residues (1HSG)   — RUN
- Starting point: 1HSG loaded from RCSB (asymmetric unit). `2c` -> `1` -> sub-menu (1 Analyze, 2 Change classification,
  3 Configure classification settings, 4 Update modified amino acid mapping table, 5 Return) -> `1`.
- What the user sees: panel "Analysis Methodology" with a 6-rule order: user classification > metal content → metal_site (MCPB) >
  peptide backbone → modified_amino_acid > nucleic acid/carbohydrate → unsupported > small organic (GAFF2) > unknown.
  "Step 1: Extracting from RedoxSites – No RedoxSites detected", "Step 2: Scanning Structure ✓ Found 1 … MK1 (B:902) - 45 atoms",
  "Downloading https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz..." (first time only), "✓ MK1 (B:902): small molecule".
  Table: "MK1 (B:902) | Small Molecule | Organic small molecule | Small Molecule Parameterizer available".
- Other structures analyzed: 1ZNI gave 3 Zn metal sites + 3 **CL as "Unknown — Auto-classification could not decide"**.
  3AJ4 gave SEP "modified amino acid — (peptide backbone)", "Modified AA Parameterizer (de-novo)", and 3× EDO small molecule.

### Scenario 3: Classification settings and reclassification   — RUN
- `3` (settings): table "Small Molecule Atom Range 2-200 atoms" (the only configurable item), Metal detection, CCD, Backbone detection always on.
  "Modify small molecule atom count range? [y/n] (n)" -> `n`.
- Inconsistency: this screen lists "Priority 2: RedoxSite Classification", but the analysis screen says units are classified
  "by STRUCTURAL CONTENT … not by any label carried over from redox detection".
- `2` (change classification): table of residues -> "Enter residue number to reclassify (or 'q' to cancel)" -> `1` -> choices
  1 Modified amino acid, 2 Small molecule (← Current), 3 Metal site, 4 Unknown -> `2`.
  "Successfully reclassified MK1 from 'Small molecule' to 'Small molecule' … ✓ Classification updated and saved to workspace".

### Scenario 4: Modified amino acid mapping table   — RUN
- `4` -> panel "Manage Classification Mappings" -> actions `view_aa / edit_aa / view_user / edit_user / clear_user / back`.
- `view_aa`: 148 built-in mappings (e.g. MSE→MET, SEP→SER, TPO→THR, PTR→TYR, CSO→CYS, AIB→ALA).
  Oddity: CSD (3-sulfinoalanine, a cysteine derivative) maps to ALA.
- `edit_aa` -> "Modified amino acid code (3-letter):" `XYZ` -> numbered list of 20 parents -> "Choose parent amino acid [1/…/20]".
  Typing `ALA` is rejected ("not one of the available options"); `1` accepted.
  Saved to `~/.proprep/modified_aa_map.json` (here runs/s6/home/.proprep/modified_aa_map.json).

### Scenario 5: Small molecule MK1 (indinavir, 1HSG) — steps 1–3   — RUN
- Steps:
  1. `2c` -> `3` "Select a residue to parameterize" (table) -> "Enter site number(s) to parameterize (or 'q' to cancel)" -> `1`.
     Note: "Select the sites … one (e.g. '3') or several (e.g. '1,2,3'). Each goes to its own parameterizer."
  2. Panel "Residue Naming" -> "Output residue name (MK1):" -> Enter (rules: 1–3 characters, starts with a letter, alphanumeric).
  3. Checklist "Small Molecule Parameterization": 1 Extract Coordinates, 2 Hydrogen Addition, 3 Gaussian Input, 4 RESP Charges,
     5 Bonded Parameters (parmchk2), 6 AMBER Topology Files (tLEaP), 7 Parameter Refinement (optional), 8 Force Field Integration.
     Controls: "[n]ext | [v]iew structure | [s]ave & quit | [q]uit | [#] to run/re-run step | [#s] skip optional step".
  4. `n` -> Step 1: "✓ Extracted coordinates to mk1.pdb" (45 atoms).
  5. `n` -> Step 2: "Hydrogen atoms: 0 (0.0%) Low hydrogen content detected … typically ~40-60%" -> "Would you like to add hydrogens using reduce? (y)" -> `y`.
     Panel explains NOFLIP vs FLIP, -OH, -NUCLEAR, -Keep. Prompts: "Mode [NOFLIP/FLIP] (NOFLIP)" Enter; "Add hydrogens to OH and SH groups? (y)" Enter;
     "Use nuclear X-H distances? (n)" Enter; "Keep original bond lengths …? (n)" Enter.
     "Running reduce: reduce -NOFLIP -OH mk1.pdb > mk1_H.pdb ✓ Added 47 hydrogen atoms" (correct for neutral indinavir C36H47N5O4).
     "Launch interactive structure viewer? (y)" -> `n`; "Edit H atoms? [add/remove/done] (done)" -> `done`.
  6. `n` -> Step 3: panel "Gaussian Calculation Overview" (keywords and the anion/SCRF warning). Prompts:
     "Net charge for MK1 (0)" -> `0`; "Spin multiplicity (1)" -> (I typed `n` by mistake: "Invalid multiplicity, using 1");
     "Memory allocation (10GB)" Enter; "Number of processors (4)" Enter; "Customize calculation keywords? (n)" -> `n`.
     "✓ Generated Gaussian input: mk1.gjf … Charge: 0, Multiplicity: 1, Memory: 10GB, Processors: 4, Atoms: 92".
- Result / files (runs/s6/ligand/hsg/small_molecule_params_MK1/):
  - mk1.pdb – ligand heavy atoms extracted from 1HSG; mk1_H.pdb – with reduce hydrogens;
  - mk1.gjf – two-job Gaussian input: job 1 `# opt freq b3lyp/6-31+G(d) nosym int=ultrafine IOp(7/33=1)`, charge/mult `0 1`, 92 atoms;
    `--link1--` job 2 `#P HF/6-31G(d) Pop=mk IOp(6/33=2,6/41=10,6/42=10) nosym int=ultrafine guess=read geom=allcheck` (%oldchk=mk1.chk, %chk=mk1_esp.chk);
  - workflow_state.json – resumable checklist state.
- What the user must do with Gaussian: run `g16 mk1.gjf` in that folder so that `mk1.log` appears (plus `mk1.chk`;
  `formchk mk1.chk mk1.fchk` is needed later for Seminario refinement). Then resume (Scenario 7).
- If the charge is negative, ProPrep shows an "Anion Solvation Warning" and asks "Add implicit solvation (SCRF) to Gaussian keywords?
  (y)" and "Solvent name (Water)" (source; seen live in the metal route).
- Session: runs/_sessions/s6-sm-1/transcript.txt; runs/s6/ligand/hsg/proprep_session_20261121_123644.json

### Scenario 6: Small molecule step 4 without Gaussian, and forcing later steps   — FAILED (expected) / bug
- `n` -> Step 4: "✗ Step 4 failed: Gaussian output not found: mk1.log  Run Gaussian first: g16 mk1.gjf".
  The user also sees a raw Python **traceback** ("RuntimeError: Gaussian output not found …"). It looks like a crash but is not one.
- `5` (jump to step 5) -> "Step sm-5 requires these steps first: sm-4  Run anyway? (may fail) [y/n] (n)".
  Typing `7s` or `s` here gives "please enter y or n". With `y`: "Running parmchk2... ✗ Step 5 failed: sequence item 2: expected str instance, NoneType found"
  plus a TypeError traceback. **Bug-ish:** "Run anyway" without step 4 gives a meaningless error instead of "no mol2 yet".
- `7s` -> "Skipped Parameter Refinement." (step 7 marked ○).
- There is **no AM1-BCC (or any non-Gaussian) charge option** in this route. Step 4 only runs `antechamber -c resp` on the Gaussian log.

### Scenario 7: Save & quit, check status, resume, bad name, quit without saving   — RUN
- `s` -> "State saved. You can resume later." back at the 2c menu.
- `4` (status) -> "Parameterization Status Summary / No parameterizations found in workspace". This is misleading, because a
  half-finished MK1 job exists on disk.
- `3` -> `1` -> naming prompt -> typed `9E2` -> "'9E2' is not a valid residue name for tLEaP. Use 1-3 characters … (e.g. X9E)" -> `MK1`
  -> panel "Resume Available … Started / Last updated / Working directory … Progress: 3/8 steps completed" ->
  "Resume from saved state? (y)" -> `y` -> checklist reappears at "[n]ext (6. AMBER Topology Files)". Answering `n` deletes the saved state and starts fresh (source).
- `q` -> "Quit without saving? (Progress will be lost) [y/n] (n)" -> `y` -> back to 2c menu.

### Scenario 8: A simpler ligand – ethylene glycol (EDO, 3AJ4), and several copies of one ligand   — RUN
- 3AJ4 has three EDO molecules. Each is listed as a separate unit (Instances = 1).
- `3` -> `2,3` (EDO A:112 and EDO B:111). EDO A:112: name Enter; step 1 "Extracted 4 atoms"; step 2 reduce "Added 6 hydrogens (10 atoms total)";
  step 3 charge `0`, multiplicity `1`, defaults -> "Generated edo.gjf (charge=0, mult=1)". `s` save & quit.
- ProPrep then automatically launched the second unit, EDO B:111, **in the same folder** small_molecule_params_EDO. After the name prompt it offered
  "Resume Available … Progress: 3/8 steps completed", i.e. the *first* copy's state.
  **Hazard:** answering `n` there would delete the first copy's work. For identical ligands, parameterize one copy only.
  (The 2c "Every metal site needs parameters…" message says the same for metals: reuse the transformer emitted at the end instead of re-deriving.)
- Files: runs/s6/modaa/maa/small_molecule_params_EDO/{edo.pdb, edo_H.pdb, edo.gjf, workflow_state.json}.

### Scenario 9: Small molecule steps 5–8 (after Gaussian)   — SOURCE-ONLY (Gaussian not installed)
- Step 4 runs `antechamber -i mk1.log -fi gout -o mk1.mol2 -fo mol2 -c resp -nc <charge> …`, checks the frequencies, writes
  mk1_atom_name_mapping.json.
- Step 5 runs `parmchk2 -i mk1.mol2 -f mol2 -o mk1.frcmod -s 2 -a Y`, shows a penalty table, and asks "Parameters to refine" (numbers, `all`,
  Enter/`none` to skip). Bonds/angles are refined by Seminario (needs mk1.fchk: "formchk mk1.chk mk1.fchk").
  Dihedrals use "Dihedral method": 1 PES scan, 2 CREST, 3 skip (default 1 if ≤2 dihedrals, else 2). CREST/xtb are not installed.
- Step 6 asks "tLEaP mode": `gas` / `aqueous` / `both` (default both). Aqueous asks "Select water model" and "Box size (Å)".
  Writes mk1_gas.parm7, mk1.lib, etc.
- Step 7 (optional): Seminario writes mk1_seminario.frcmod. PES asks "Number of scan steps (24 = 15° increments for 360°)" and "Have you completed the PES scan calculations?".
- Step 8 "Force Field Integration" registers lib + frcmod for the Topology Generator, then offers the library deposit
  ("Reuse these parameters? … Add to your library (n)"), see the library section below.

### Scenario 10: AM1-BCC charges for MK1 (done outside ProPrep, as a user would)   — RUN (outside the app)
Because the app has no AM1-BCC option, I made AM1-BCC parameters with the installed AmberTools from ProPrep's own
hydrogenated ligand. Then I deposited them with 2c option 2 (next section). Commands (runs/s6/bcc/):
```
antechamber -i mk1_H.pdb -fi pdb -o mk1_bcc.mol2 -fo mol2 -c bcc -nc 0 -at gaff2 -rn MK1 -s 0 -pf y   # ~4 min (sqm AM1, 92 atoms)
parmchk2 -i mk1_bcc.mol2 -f mol2 -o mk1_bcc.frcmod -s 2
tleap: source leaprc.gaff2 / MK1 = loadmol2 mk1_bcc.mol2 / loadamberparams mk1_bcc.frcmod / saveoff MK1 mk1_bcc.lib
```
Result: mk1_bcc.mol2 (GAFF2 types + AM1-BCC charges), mk1_bcc.frcmod (only DIHE entries, all penalty 0), mk1_bcc.lib (unit MK1).
tleap "Errors = 0; Warnings = 0".

### Scenario 11: Modified amino acid – phosphoserine SEP (3AJ4), de-novo route, steps 1–4   — RUN
- How the structure was found: 1a -> RCSB -> "4. Search by ligand (code or name)" -> `SEP` -> SEP "PHOSPHOSERINE C3 H8 N O6 P, 38 structures" ->
  Enter, Enter, Enter -> list sorted by resolution -> `1` (3AJ4, 1.00 Å) -> "Download 3AJ4? (y)" -> "Remove hydrogen atoms? (y)" `y` (1804 H removed).
- Steps: `2c` -> `1` -> `1` (analysis) -> `5` -> `3` -> `1` (SEP).
  "Parameterizing Modified Amino Acid: SEP  Full Name: PHOSPHOSERINE  ⚠ Note: The CCD may have an entry for this residue code that differs from
  your intended compound (e.g. 'TYO' …). Parent Standard Residue: SER" -> "Use SER as parent residue? (y)" -> `y`.
  "IMPORTANT WORKFLOW INFORMATION … will use the parent amino acid (SER) structure as a starting point. During the workflow, **you will need to
  modify this structure to create the modified amino acid (SEP)**." Output dir modified_aa_params_SEP.
  Checklist (10 steps): 1 Initial Structure (tleap), 2 Backbone Angles (cpptraj), 3 Gaussian Input, 4 PES Scan Results, 5 Extract Structures,
  6 ESP Generation, 7 AC File, 8 RESP Charges (residuegen), 9 Bonded Parameters (parmchk2 + GAFF2 patching; optional Seminario/torsion scan),
  10 Force Field Integration.
  - `n` step 1: "✓ Generated ace-ser-nme.parm7 / .rst7 / .pdb".
  - `n` step 2: "✓ ser.ahelix.mol2, ser.bsheet.mol2".
  - `n` step 3: "Enter memory for Gaussian job (8GB)", "number of processors (4)", "Gaussian keywords (HF/6-31G* opt=(modredundant) freq iop(6/33=2,6/42=6) pop=mk)",
    "title (Dihedral scan for amino acid)", "Enter charge for the molecule (0)" -> `-2`, "Enter multiplicity (1)" -> `1`,
    "Auto-generate a relaxed sidechain torsion scan …? (n)" -> `n`. ProPrep finds the φ atoms [5,7,9,16] and ψ atoms [7,9,16,18], then "✓ Generated ser_ahelix.gjf, ser_bsheet.gjf".
  - `n` step 4: panel "Expected files: ser_ahelix_pes.log, ser_bsheet_pes.log … Simple: save output as ser_ahelix.log … PES scan: …_pes.log …
    auto-detects which you chose based on filenames" -> "Proceed (y), Pause to run Gaussian (p), or Quit (q)? (p)" -> `p` ->
    "⏸ Step 4 paused: Waiting on an external calculation". `s` save & quit.
- Result files (runs/s6/modaa/maa/modified_aa_params_SEP/): tleap.in, leap.log, ace-ser-nme.*, cpptraj.in, ser.ahelix.mol2, ser.bsheet.mol2,
  ser_ahelix.gjf, ser_bsheet.gjf (each 25 lines of atoms = **ACE-SER-NME, 22 atoms, charge line `-2 1`, with φ/ψ frozen: `D 5 7 9 16 F`, `D 7 9 16 18 F`**).
- **Important for the manual:** the generated Gaussian input is still plain *serine* (no phosphate) but carries the charge you typed (−2).
  The user must open each .gjf in GaussView (or similar), build the phosphate onto OG, and keep the atom order for the frozen dihedrals.
  Only then run Gaussian. Running it unchanged would compute a serine dianion. Also: SEP is already in AMBER's phosaa19SB library (2b Scenario 9b),
  so this whole route was unnecessary for SEP. Users should check the Explorer first.
- Steps 5–10 SOURCE-ONLY: extract optimized conformers, HF ESP, AC file from lowest-energy structure, residuegen multi-conformer RESP,
  parmchk2, write sep.prep / sep.frcmod / sep.lib, register for Topology Generator, offer library deposit.
  There is also a "from-structure (conjugate)" route when a RedoxSite bundles an amino acid with a covalent partner, and a protonation-change
  procedure (scan neutral, re-optimize charged with dihedral/atom freezing) described in Help. Neither was run.
- Session: runs/_sessions/s6-maa-1/transcript.txt; runs/s6/modaa/maa/proprep_session_20261121_143651.json

### Scenario 12: Metal site – 1ZNI Zn B:31, without Redox Site Detector, chloride ions kept   — FAILED at step 10
- Steps: 1ZNI loaded; 2c analyze -> 3 Zn metal sites + 3 CL "Unknown". `3` -> `1` (ZN B:31) ->
  "Metal Site Parameterization: 1 site selected … will run: ZN (B:31) not run: ZN (B:32), ZN (D:31) … Every metal site needs parameters before the
  topology build will succeed. For a site equivalent to one being run, apply the reuse transformer emitted at the end … No redox sites detected yet --
  will be detected during preprocessing (Step 0f)" -> "Proceed …? (y)" `y`.
- 15-step checklist "Metal Site Parameterization":
  STRUCTURE PREPARATION 1 Structure Filtering (optional), 2 Structure Completeness (optional), 3 Structure Triage;
  COMPONENT PARAMETERIZATION 4 Protein: FF + Protonation, 5 Water: Model Selection, 6 Organic Small Molecules, 7 Organometallic Small Molecules,
  8 Metal Clusters, 9 Isolated Metal Ions; STRUCTURE ASSEMBLY 10 Structure Recombination, 11 Redox Site Sync;
  METAL SITE PARAMETERIZATION 12 MCPB Atom Typing, 13 Bonded Parameters, 14 RESP Charges, 15 Force Field Integration.
- Step 1 = the PDB Filter: chain interface (buried surface area) map, "Select chains … (all)", then "[r] Retain / [s] Select / [d] Discard" for every
  chain component (11 prompts); kept everything.
- Step 2 = Structure Fixer: chose `5` skip missing-residue detection, `1` REMARK 470, `1` AltLoc identifiers. 4 alternate locations found
  (chains B, C, D). "Repair …? (y)" -> "MODELLER not available — gap filling and mutations are disabled. Terminal capping (ACE/NME) is still available".
  For each chain "Cap which terminus of chain X? (s)" -> `s`. "No repaired structure produced, using current". The altlocs were not resolved.
- Step 3 Triage: "102 protein, 3 organic, 3 isolated metal, 103 water". **The 3 chloride ions are filed under "B Organic small molecules"**.
- Step 4: protein FF menu (ff19SB ★, ff14SB, ff15ipq, fb15, ff03.r1, ff99SBildn, Constant pH/Redox ff10 variants) `1`;
  modified-AA add-ons `1` none. Then the protonation analysis: pH 7.0, threshold 0.5, PROPKA `y`, pKa table, "Is this for constant pH simulation? (n)",
  "Set all residues to standard …? (n)", per-type "Use recommended names for all CYS/GLU/HIS/LYS/TYR residues? (y)".
  All CYS were recommended as CYS (pKa 99.99), although insulin's cysteines are disulfide-bonded.
  His B10, which binds Zn through NE2, was recommended as HIE. A metal-bound NE2 should carry no H, so HID would be right.
  structure_with_md_names.pdb was written.
- Step 5 water: OPC (Recommended) … TIP3P-FB/TIP4P-FB; chose `1` OPC.
- Step 6 (chlorides as "organic residues"): "[1] I have parameters (lib/frcmod files) [2] I need to generate parameters".
  - CL B:33: `1`, "FF residue name in parameter files (CL)" -> `Cl-`. The lib file browser here has **no path entry**, so I copied AMBER's atomic_ions.lib into the project
    and picked it. Frcmod: `skip`. **Bug:** "Atom name mismatch … heavy atoms only in the lib: AG". It reads the *first* unit of the multi-ion library (Ag+),
    not the `Cl-` unit named. Mapping prompt -> `all`.
  - CL B:34: `2` launched a nested Small Molecule Parameterizer for a single chlorine atom (Gaussian needed); `s` -> "Parameterization saved (incomplete) for CL".
  - CL B:35: same as B:33.
- Steps 7–9: "No organometallic residues", "No pure metal clusters", "3 metals will be removed before tLEaP, reinserted after".
- Step 10: "Residue CL → AG … Configured tLEaP with 1 lib files". Shows the tLEaP template and asks "Edit the tLEaP template? (y)" -> `n`.
  "No bond definitions found. Run comprehensive RedoxSite detection first. Do you want to continue without bond definitions? (n)" -> `y`.
  Then: output name, "View the tLEaP input file? (y)", "Run …? (y)". Result: "FATAL: Atom .R<AG 82>.A<CL 2> does not have a type", "✗ Step 10 failed: tLEaP failed".
- Workaround that worked (Scenario 13): **discard the chloride ions in step 1**. They are crystallization ions; tLEaP can add counter-ions later.
- Session: runs/_sessions/s6-metal-2/transcript.txt; runs/s6/metal/zni/

### Scenario 13: Metal site – 1ZNI, chlorides discarded, detection done inside step 11   — RUN to step 13 (QM inputs empty: bug)
- Same as 12, but in step 1 "Chain B - CHLORIDE ION" -> `d`. His B10 set to HID with "Use recommended names for all HIS? (y)" -> `n` ->
  "Customize HIS names ()" syntax "Single name … Overrides: 94,119:HIE … Chain-specific: A:94:HIE" -> `10:HID`.
- Step 10: tLEaP ran with the ProPrep teLeap: "✓ tLEaP completed successfully … Errors = 0; Warnings = 25". ParmEd validation panel explains
  validity check / rediscover molecules / optional 12-6-4. Warning: "two cysteine residues whose sulfur atoms are within 3 Angstroms. Rename CYS to CYX …"
  (disulfides **not** made). "Inserted isolated metal: ZN as residue 206/207/208 … ✓ Created prepared_structure.pdb ✓ Collected 1988 atoms (types + charges)".
- Step 11 "Redox Site Sync": "No existing redox sites found … Tip: Run Redox Site Detector before Force Field Parameterizer. Run interactive redox site detection now? (y)" -> `y`.
  The full Redox Site Detector runs. It only offers the *filtered* or *RCSB* structure, not prepared_structure.pdb. Findings:
  6 disulfides; ZN B:31 "No nearby atoms found" (it sits on the crystal 3-fold axis, so its His ligands are in symmetry copies);
  ZN B:32 – HIS B:10 NE2 1.9 Å; ZN D:31 – HIS D:10 NE2 2.1 Å + HOH D:32 O 2.3 Å. I defined SG–SG bonds for each disulfide and Zn–N/O bonds
  (about 6 typed answers per bond) -> "Created 9 fully refined redox sites".
- Step 12 "MCPB Atom Typing" loops over all metal sites (7, 8, 9), not only the one chosen:
  - Skips disulfide sites with the confusing message "organic cofactor (no metal) — parameterized by the Small Molecule Parameterizer, not MCPB".
  - Per metal: "Enter formal charge for ZN … (e.g., +2, +3, -1)" -> `+2` (Enter alone = "Invalid charge format", repeated); "Enter unpaired electrons … (e.g., 0, 1, 5)" -> `0`.
  - Sites 8 and 9: "Cluster Core Atom Charges — N atom(s) belong to a withheld cluster residue, so they carry no prmtop charge …" followed by one
    "Formal charge for HIS D:10 atom N (N) …" prompt **per ligand atom** (10–11 prompts). These are ordinary His/water ligands, not a cluster. Answered 0.
  - Site 8 warnings: "Atom CG in HIS:10 … not found in PDB file" ×5 (alternate conformations).
  - Small/large model questions: "Fill gaps between selected residues? (y)", "Fill gaps up to how many residues? (5)", "Use GLY for gap filling (vs actual PDB residues)? (y)".
  - **Bug:** every model reports "Total residues: 0 … ✓ Wrote PDB file: small.pdb (0 atoms, charge=0.00)" and the same for large.pdb.
    "Suggested charge from PDB: 0" is offered even after +2 was entered for Zn.
  - Gaussian prompts, small model: "Memory allocation (GB) (4)", "Number of processors (4)", "Total charge of small model (0)" (typed 2 for site 9),
    "Spin multiplicity (1)", "Use recommended method and route (B3LYP/6-31G*)? (y)", freezing options 1 Full / 2 heavy atoms / 3 backbone / 4 custom.
    "✅ Generated Gaussian input: small_freq.gjf", then Next-steps panel "cd …/site_N/models; g16 small_freq.gjf; formchk small_freq.chk small_freq.fchk".
    Large model: "Use same resources …? (y)", "Total charge of large model (0)", "Spin multiplicity", "Add implicit solvation (SCRF) to the large model? (n/y)",
    "Use the same method as the small model? (y)", freezing, "Use standard ESP route keywords? (y)" (panel explains Pop(MK,ReadRadii), IOp(6/33=2,6/41=10,6/42=17))
    -> large_resp.gjf. "✓ Step 12 completed".
- Step 13: "Site 7/8/9: Gaussian output not found Expected: …/small_freq.fchk and small_freq.log … Checkpoint: Gaussian geometry optimization required.
  ⏸ Step 13 paused". Jumping to `14`: "Step mcpb-3 requires these steps first: mcpb-2  Run anyway? (n)".
- Files (runs/s6/metal2/zni/metal_site_params_ZN_B_31/): metal_free.pdb/.prmtop/.rst7/_tleap.in, metal_free_parmed_replay.py, tleap_output.pdb,
  prepared_structure.pdb, site_7|8|9/models/{small.pdb, large.pdb, small_freq.gjf, large_resp.gjf, standard.fingerprint, large.fingerprint, atom_type_assignments.json}.
  **The .gjf files contain only the route, title and `charge mult` line, with no atoms.** For example site_9 small_freq.gjf has
  `# B3LYP/6-31G* Opt Freq Geom=PrintInputOrient Integral=(Grid=UltraFine) IOp(7/33=1)` / `2 1` / (nothing).
  large_resp.gjf ends with the metal radius "Zn 1.373" twice.
- Likely cause (observation, not fixed): prepared_structure.pdb has no chain IDs and renumbered residues (e.g. His B10 → residue 31, Zn → 206+),
  while the sites use the original chain:number labels, so no residues are matched.
- Session: runs/_sessions/s6-metal-3/transcript.txt; runs/s6/metal2/zni/proprep_session_20261121_130836.json

### Scenario 14: Metal site – 1ZNI with Redox Site Detector (2a) run first   — FAILED at step 10 (bug)
- Followed the in-app tip: 2a first (Zn only + 6 disulfides with bonds, site ZN D:31 with His D:10 and HOH D:32), then 2c -> parameterize ZN D:31.
- New prompt in step 4 when redox sites exist: a table of redox-site residues and "Enter indices to exclude (e.g., 1,3,5-7 or 'all') ()". Its note:
  "Residues coordinated to metal centers typically exhibit unusually shifted pKa values". In step 1 filtering, chains holding redox sites offer
  extra choices "[c] Keep redox site residues [f] Keep redox site + and - N flanking residues [o] Override redox protection and discard (⚠️ breaks redox sites)".
- Step 10 writes bond commands with bare PDB residue numbers: "disulfide: 3 bonds" (6 were defined), and `bond mol.7.SG mol.7.SG` (A:7–B:7 became a self-bond).
  Also `bond mol.31.ZN mol.10.NE2` to a zinc that was removed before tLEaP. tLEaP: "Fatal Error! bond: Argument #2 is of type String must be of type: [atom]".
- Session: runs/_sessions/s6-metal-6/transcript.txt (resumed with `proprep --resume-session …/proprep_session_20261121_133619.json`, see notes).

### Scenario 15: Metal site – rubredoxin Fe(III)(Cys)4 (1IRO), 2a first   — FAILED at step 10 (same bug)
- 1IRO: one chain, Fe(III) bound by Cys 6, 9, 39, 42, no disulfides. 2a: "Fe … CYS A:42 SG 2.2 Å, A:9 2.2, A:39 2.3, A:6 2.3". A pairwise SG–SG table
  is shown when none are within 2.5 Å ("Specify disulfide bonds as index pairs … or 'none'") -> `none`. Four Fe–SG bonds defined.
- In 2c, steps 1 and 2 were skipped with `1s`, `2s` ("Skipped Structure Filtering.").
- Step 4 CYS table has a "Redox site" column and the message "a metal-bound cysteine is normally CYM". PROPKA recommended CYM only for Cys6.
  I typed `n` then `CYM` -> "CYM: 4 residues".
- Step 10 tLEaP input contains `bond mol.55.FE mol.42.SG` … while FE has been removed from metal_free.pdb -> "Fatal Error! bond: Argument #1 …".
- Workaround tried: re-ran 2a with "2. Skip bond definition" for the site, resumed, re-ran step 10 ("Reconfigure all settings …? (n)").
  The old Fe–SG bond lines were still written and it failed again.
- Session: runs/_sessions/s6-iro-2/transcript.txt; runs/s6/metal4/iro/

### Scenario 16: Metal site – rubredoxin Fe(III), detection only at step 11   — RUN to step 13 (QM inputs empty: bug)
- No 2a beforehand. Steps 1–2 skipped, step 4 ff19SB + CYM ×4, OPC, steps 6–9 automatic. Step 10: "Do you want to continue without bond definitions?" `y`
  -> tLEaP OK, "✓ Reinserted 1 metal(s) ✓ Collected 1221 atoms".
- Step 11: detection run on the H-stripped structure, Fe site with 4 Fe–SG bonds -> "site_1: 25 atoms, 4 bonds".
- Step 12: "Enter formal charge for FE A:55 …" -> `+3`; "Enter unpaired electrons …" -> `5` (high-spin d5).
  Then 24 "Formal charge for CYS A:6 atom N …" prompts (answered 0, and -1 for each SG).
  Renaming: FE → M1, SG → Y1–Y4 (standard.fingerprint lines like "55-FE-FE 1222 FE -> M1", "6-CYS-SG 106 S -> Y1", "LINK 1222-FE 616-SG").
  Small model "Suggested charge from PDB: 0" (wrong; correct is Fe³⁺ + 4 thiolates = −1) -> typed `-1`; "Suggested multiplicity (2S+1): 6" (correct) -> Enter.
  "Anion Solvation Warning (charge = -1) … Add implicit solvation (SCRF)? (y)" -> `n`.
  "Use recommended method and route (B3LYP/6-31G*)? (y)" -> `n` -> "DFT functional (B3LYP)", "Basis set (6-31G*)",
  "Full route line (# B3LYP/6-31G* Opt Freq Geom=PrintInputOrient Integral=(Grid=UltraFine) IOp(7/33=1))".
  The source also has GenECP/ECP prompts for heavy metals (not shown for Fe with 6-31G*).
  Large model: charge `-1`, multiplicity 6, SCRF default **y** here ("Solvent name (Gaussian keyword) (Water)").
- Again "Total residues: 0" and both .gjf files have only `-1 6` and no atoms (large_resp.gjf: `… Pop(MK,ReadRadii) … SCRF=(Solvent=Water)`, "Fe 1.400" twice).
  prepared_structure.pdb again has no chain letter and Fe renumbered to residue 164.
- Step 13 paused waiting for site_1/models/small_freq.fchk + small_freq.log.
- Files: runs/s6/metal5/iro/metal_site_params_FE_A_55/site_1/models/*. Session: runs/_sessions/s6-iro-3/transcript.txt; runs/s6/metal5/iro/proprep_session_20261121_142210.json

### Scenario 17: Metal steps 13–15 and the bonded-vs-nonbonded choices   — SOURCE-ONLY
- What the user must do with Gaussian (from the on-screen Next Steps): in each `site_N/models/` run `g16 small_freq.gjf` then
  `formchk small_freq.chk small_freq.fchk`, and `g16 large_resp.gjf`. Then resume the checklist (2c -> 3 -> same site -> "Resume from saved state? y").
- Step 13 (MCPB "step 2"): can pick existing outputs or "Generate new Gaussian input (if Step 1 input was modified)". It asks for the log/fchk paths
  ("Enter path to existing fchk file or type 'convert' to generate from .chk"), "Seminario frequency scaling factor", and writes a pre-frcmod and then
  a Seminario frcmod (bonds/angles to the metal).
- Step 14 (MCPB "step 3a–3d"): ESP extraction from large_resp.log, RESP input generation and fit ("charge_fit" folder; "Select backbone restraint level"),
  then mol2 files.
- Step 15: combines everything into tLEaP-ready lib/frcmod, emits a reuse transformer for equivalent sites, and offers the library deposit. For metals it also asks
  "Redox state label (as_parameterized)", "Spin state label (default)" and "Formal charge of the site".
- Metal models available: ProPrep implements the **bonded (MCPB.py-style) model** for the site it parameterizes. The other metals in the structure are treated
  as **isolated ions**: removed before tLEaP and reinserted (step 9). For a metal–water contact, the preprocessor's bond editor can instead keep the
  water **restrained/nonbonded** ("Metal-water contact: model the water as RESTRAINED … held near the metal by an MD distance restraint instead of a bonded MCPB ligand? (n)").
  This prompt did not appear in the Redox Site Detector bond editor I used. The nonbonded **12-6-4** Li/Merz ion option belongs to the topology/ion choices
  (ParmEd panel: "12-6-4 LJ correction — not requested for this system"), not to this parameterizer.

### Scenario 18: NMR zinc finger 1ZNF (aborted)   — RUN (stopped early)
- 1ZNF loads 37 NMR models (74 HETATM ZN lines in the file) and asks "This structure contains 7844 hydrogen atoms … Remove hydrogen atoms? (y)".
  Abandoned in favour of crystal structures, because every model's Zn would appear. Session: runs/_sessions/s6-metal-1/transcript.txt.

### Options reference (2c)
| Where | Prompt / choice | Options and meaning |
|---|---|---|
| 2c menu | 1–5 | 1 analyze/configure, 2 import into library, 3 parameterize (new or resume), 4 status, 5 help; `b` back, `m` main, `x` exit |
| Analyze sub-menu | 1–5 | 1 analyze, 2 change classification (1 modified AA, 2 small molecule, 3 metal site, 4 unknown), 3 settings (only small-molecule atom range, default 2–200), 4 mapping table (`view_aa/edit_aa/view_user/edit_user/clear_user/back`), 5 return |
| Parameterize | site numbers | one (`3`) or several (`1,2,3`); each runs its own checklist one after another |
| Any checklist | action | `n` next, `#` run/re-run step, `#s` skip optional step, `v` view structure, `s` save & quit (resumable), `q` quit (asks "Quit without saving?") |
| Small mol. naming | Output residue name | 1–3 alphanumeric, starts with a letter (`9E2` rejected) |
| Small mol. H addition | add H? / mode / -OH / nuclear / keep | default y / NOFLIP (FLIP flips Asn/Gln/His, protein-oriented) / y / n (1.0 Å X-ray vs 1.08 Å nuclear) / n |
| Small mol. Gaussian | net charge, multiplicity (2S+1), SCRF (anions only, default y), solvent (Water), memory 10GB, processors 4, customize keywords (n) | opt: `opt freq b3lyp/6-31+G(d) nosym int=ultrafine IOp(7/33=1)`; ESP: `HF/6-31G(d) Pop=mk IOp(6/33=2,6/41=10,6/42=10)` |
| Small mol. charges | (none) | RESP from Gaussian only; **no AM1-BCC option** |
| Small mol. refinement | parameters to refine; dihedral method | numbers/`all`/`none`; 1 PES scan, 2 CREST (needs crest/xtb), 3 skip |
| Small mol. tLEaP | tLEaP mode | gas / aqueous / both (default) |
| Modified AA | parent residue; memory 8GB; processors 4; keywords `HF/6-31G* opt=(modredundant) freq iop(6/33=2,6/42=6) pop=mk`; charge; multiplicity; auto sidechain scan (n); step 4 y/p/q | simple (2 conformers: α-helix, β-sheet) vs PES (`*_pes.log`) workflow chosen by log file names |
| Metal step 4 | protein FF; modified-AA add-on; exclude redox residues; pH; threshold; PROPKA; constant pH; per-type names | metal-bound Cys → CYM, metal-bound His → HID/HIE depending on which N binds |
| Metal step 5 | water model | OPC ★, OPC3, OPC3-pol, TIP4P-Ew, TIP3P, SPC/E, SPC/Eb, TIP5P, TIP3P-FB, TIP4P-FB |
| Metal step 6/7 | organic: 1 have lib/frcmod, 2 generate; organometallic: 1 params incl. metal, 2 organic part only, 3 generate organic part | |
| Metal step 12 | metal formal charge; unpaired electrons; ligand-atom formal charges; gap filling (y, 5, GLY); small/large charge & multiplicity; SCRF; method/route (B3LYP/6-31G*, or custom functional/basis/route, GenECP); freezing 1–4; standard ESP route | multiplicity suggested from unpaired electrons (5 → 6); charge suggestion is unreliable (always 0 here) |

---

## Feature: Local parameter library (deposit, find, reuse)

**What it is for (plain English):** A personal shelf of finished parameters (frcmod + lib) kept in your home folder.
Once something is deposited, every future ProPrep project can find it in the Force Field Explorer and the Topology Generator.
You do not have to parameterize the same ligand or metal site twice.

**Needs before you start:** A finished .frcmod and its .lib/.off (optional .prep, which is converted to a lib).
They can come from ProPrep or from anywhere else.

**Where to find it:** 2c option 2 "Import existing parameters into your library" (always available, no structure needed).
Also offered automatically at the final "Force Field Integration" step of each parameterizer ("Reuse these parameters? … Add to your library").
Stored in `~/.proprep/forcefield_params/specialized_residues/<small_molecules|modified_aa|metal_sites>/<ENTRY>/metadata.json` + files.

### Scenario 1: Import AM1-BCC MK1 parameters   — RUN
- Steps (s6-sm-1, at 2c menu):
  1. `2` -> panel "Import parameters … Point ProPrep at a finished .frcmod and its .lib/.off (and an optional .prep)". File browser (this one *does* accept `<path>`).
  2. `2` (bcc_params/) -> `2` mk1_bcc.frcmod -> `2` mk1_bcc.lib -> prep: `skip`.
  3. "What kind of parameters are these? 1 Small molecule / ligand … 2 Modified amino acid … 3 Metal site …" "Category (1)" -> `1`.
  4. "The library defines unit MK1; the filename says mk1_bcc. tLEaP matches the unit name and the match is case-sensitive, so MK1 is recorded as the residue name."
     Panel "Two names, two levels": Entry name = which molecule (one folder each); Set name = whose numbers (several sets per entry, exactly one default).
  5. "Library entry name (residue name) (MK1)" Enter; "One-line description (MK1 (GAFF, AM1-BCC/RESP charges, user-parameterized))" -> `Indinavir (MK1) from 1HSG, GAFF2 + AM1-BCC charges`;
     "Parameter set name (user_MK1)" -> `MK1_AM1BCC`; "Reference(s) / citation (optional; separate multiple with ' | ')" Enter.
- What the user sees: "Parameters saved — Saved to your library. …/.proprep/forcefield_params/specialized_residues/small_molecules/MK1 set: MK1_AM1BCC.
  It will now appear in the Force Field Explorer and Topology Generator in any project." Then "Imported parameters need a transformer … No detected redox sites in this session …
  To add it later: load the structure, run the Redox Site Detector, then Redox Site Preparer -> 'Create transformer (interactive PDB editor)'. With nothing to rename, type 'save' straight away and accept the pass-through. Link it to: small_molecules/MK1".
- Result: metadata.json (description, source "user", single_state/default, residue_name MK1, files frcmod/lib, prerequisites leaprc.gaff2 detected automatically) + mk1_bcc.frcmod + mk1_bcc.lib.
- No atom-type questions were asked because GAFF2 types are already known. For new types, the wizard asks "Element for atom type …" and "Hybridization for …" (source).

### Scenario 2: Import the same set again (name collision)   — RUN
- Same answers again -> "Set 'MK1_AM1BCC' already exists for small_molecules/MK1 [single_state/default]. 1. Save as a new version — … MK1_AM1BCC_v2
  2. Overwrite the existing set 3. Cancel — leave the library untouched -- the default" -> `1` -> "set: MK1_AM1BCC_v2".
- Problems: (a) v2 became the default set and v1 was demoted, without being told. (b) Both sets point to the **same file names in the same folder**
  (mk1_bcc.frcmod/.lib), so "keep the existing set" does not keep the old files if the new ones differ. (c) The v2 description reverted to the
  auto text "AM1-BCC/RESP charges", because I accepted the default. The default description cannot tell which charge method was used.

### Scenario 3: Find and reuse it later   — RUN
- In the Explorer, category 12 now shows two MK1 rows (2b Scenario 10). Loading one gives residue MK1 with GAFF2 types and AM1-BCC charges.
- Use in Topology Generator: SOURCE-ONLY here (other scope). It needs the transformer described above.

### Scenario 4: Deposit offered at the end of a parameterizer   — SOURCE-ONLY
- After Force Field Integration: panel "Reuse these parameters? Save these parameters to your personal force-field library? … available in every future ProPrep project"
  -> "Add to your library (n)". If yes, the same entry/description/set/reference questions follow. Metal sites also get redox state, spin state and formal charge.
  Not reached because every live run stopped at the Gaussian checkpoint.

### Options reference (library)
| Prompt | Default | Meaning |
|---|---|---|
| frcmod / lib / prep | file browser (`<path>` allowed); prep `skip` | prep is converted to lib |
| Category | 1 | 1 small_molecules, 2 modified_aa, 3 metal_sites |
| Library entry name | unit name from the lib | folder = molecule; must match tLEaP unit name (case-sensitive) |
| One-line description | "<name> (GAFF, AM1-BCC/RESP charges, user-parameterized)" | free text; write the real method |
| Parameter set name | user_<name> (metal: user_<name>_RESP) | one entry can hold several sets |
| Reference(s) | empty | separate with ` | ` |
| Collision | 3 Cancel | 1 new version `_v2` (becomes default), 2 overwrite |

---

## Notes for the manual writer

**Jargon to explain:** force field, atom type, frcmod (extra parameters) vs lib/off (residue template: atoms, types, charges) vs prep/mol2;
leaprc; GAFF/GAFF2; RESP vs AM1-BCC charges; ESP (electrostatic potential); Hessian/Seminario method; PES scan; multiplicity (2S+1, e.g. high-spin Fe(III)=6, Zn(II)=1);
formal charge of a model (metal charge + charged ligands, e.g. Fe³⁺+4 Cys⁻ = −1); MCPB small model (for bonds/angles) vs large model (for charges);
bonded vs nonbonded metal model; CYM/CYX/HID/HIE; capped tripeptide ACE-X-NME; CMAP; SCRF (implicit solvent).

**Recommended beginner path:**
1. Load the structure, open **2b**, load ff19SB + OPC (+ the matching add-on, e.g. Phosphorylated AA (ff19SB)), then run **7 Type Structure Atoms**.
   Anything typed there needs no parameterization (SEP was already covered).
2. Check **category 12** for a ready-made cofactor/metal set (hemes, Fe-S, flavins, NAD(P), biopterin, Zn(Cys)4).
3. Only then use **2c** for what is left. Expect to stop at the Gaussian step and continue later (save & quit, then resume).
4. Deposit finished parameters with **2c → 2** so later projects can reuse them.

**Blockers / missing software:**
- **Gaussian is mandatory for every charge-deriving path** (small molecule RESP, modified AA, metal MCPB). There is no AM1-BCC fallback in the app.
  Users without Gaussian can make AM1-BCC parameters with antechamber outside ProPrep (2c Scenario 10) and import them (library Scenario 1).
- CREST/xtb (dihedral refinement option 2) is not installed. MODELLER is unavailable, so the metal checklist's completeness step can only cap termini.

**Things that looked like bugs (describe, don't fix):**
1. Explorer "Browse Atom Types" lists ff19SB CMAP numbers (e.g. `-0.00006`) as ~100 fake atom types (1362 / 1270 "types").
2. Metal workflow, step 12: the small/large QM models always came out empty ("Total residues: 0", .gjf files with no atoms, "Suggested charge from PDB: 0").
   Seen in 1ZNI and 1IRO. prepared_structure.pdb loses chain IDs and renumbers residues, which looks like the cause.
3. Metal workflow, step 10 when redox sites with bonds exist (the path the app itself recommends): tLEaP bond lines use raw PDB residue numbers without chains.
   Disulfides were deduplicated wrongly (6→3, one self-bond), and metal–ligand bonds point to a metal already removed. Result: fatal tLEaP error. Removing the bonds and re-running did not clear it.
4. Step 12 asks the user for a formal charge on **every atom of ordinary ligand residues** ("withheld cluster residue"): 10–24 prompts.
5. Chloride ions: "Unknown" in analysis, "Organic small molecules" in triage. Supplying AMBER's atomic_ions.lib maps CL to the first unit (Ag+) → residue renamed AG → tLEaP fatal.
   Practical advice: discard crystallization ions in step 1.
6. Small molecule: "Run anyway" on step 5 without step 4 gives "sequence item 2: expected str instance, NoneType found" plus a traceback. Normal "Gaussian output not found" is also shown with a Python traceback.
7. 2c option 4 "View parameterization status" says "No parameterizations found" while a saved half-finished job exists.
8. Identical ligands (3× EDO) run one after another in the same folder. The second offers to "resume" the first, and answering `n` deletes it.
9. Help text for metal sites still says "Coming soon".
10. Library: "new version" makes v2 the default and reuses the same file names. User entries appear under "Redox Cofactors (built-in)" and are reported as "bundled".
11. Defaults that contradict the ★ marks: Explorer default FF selection is ff14SB + TIP3P. Insulin's disulfide cysteines were recommended as CYS. A Zn-bound His was recommended as HIE.
12. Mapping table: CSD → ALA (CSD is a cysteine derivative).
13. Cosmetic: "Loaded:" blank when only cofactor sets are loaded; "press  to return" (missing key letter) in menu suggestions.

**Common mistakes to warn about:**
- Not activating the ProPrep conda env (a stray `AMBERHOME` makes ProPrep use another AmberTools).
- Typing `q` at an Explorer filter prompt (use `b`).
- Accepting the suggested QM charge (0) for metal models; work it out by hand (metal + deprotonated Cys/Asp/Glu, etc.). Multiplicity must match the metal's spin state.
- Running the de-novo modified-AA Gaussian inputs without first building the modification onto the parent residue.
- Choosing a metal that sits on a crystal symmetry axis (1ZNI Zn B:31 had no ligands in the asymmetric unit). Build the biological assembly first, or pick a fully coordinated site.
- File browsers differ: the import wizard accepts a typed path, but the metal step 6 lib browser and the Explorer custom loader only navigate by number and `..`.

**Useful discovery:** a crashed or closed session can be replayed with `proprep --full-menu --resume-session <project>/proprep_session_*.json`,
then "1. Continue from end". It re-enters every recorded answer ("Choose option [REPLAY: r]") and resumes recording ("[Replay complete - now recording new interactions]").
The workspace (loaded structures, redox sites) is **not** restored by just reopening the project folder.

**Time notes:** AM1-BCC on 92-atom MK1 (outside app) took about 4 min. PROPKA + protonation for 1ZNI took about 1 min. The tLEaP builds took seconds.
A metal run needs 60–120 typed answers before the Gaussian step, most of them in the redox-site bond editor and the per-atom charge prompts.
