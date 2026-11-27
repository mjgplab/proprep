# User-facing issues — writer w5 (Chapters 12, 13, 19)

Format: **Symptom (what the user sees)** | Why | What to do | Chapter label

## Redox Site Detector (2a) and Redox Site Preparer (4d)

- **ACE, GOL, ACT (or EDO) listed as "organic cofactor" in the inventory** | Detection treats every non-standard residue as a possible cofactor, including caps and crystallisation additives | Deselect them at "Select additional centers to include" (type only the wanted numbers) | ch:redox
- **A "bond" of 8.21 Å is created without any warning** | Target atom numbers are re-sorted by distance for every source atom; a number from an earlier list points to another atom. No length check | Type atom names (e.g. `SG`) instead of numbers; to fix, type the source atom again and choose "2 Replace one existing bond" | ch:redox
- **The same atom appears twice in search results ("SD S (2.28Å), SD S (2.32Å)")** | Each nearby site atom counts it again | Harmless; ignore | ch:redox
- **4d shows "Compatible Transformers: None" for a single heme site** | The single-site path in 2a never shows the transformer's Bond Requirements; the 6 intra-residue bonds (C2A–CAA, C3D–CAD, CA–CB of each ligand) are missing | In 4d type `why 1`, then add the missing bonds in 2a → 4 → `edit 1` → 4 → `add`; for several sites use templates in 2a, which print the Bond Requirements | ch:redox
- **Editor drops back to the 2a menu with "Error handling menu option 'edit': invalid literal for int() with base 10"; all edits lost** | Any non-number at a ResID prompt aborts the editor; Chain/Resname are not validated, so a whole line typed at "Chain" shifts the later answers | One answer per prompt; `save` often | ch:redox
- **Editor `cancel` warns "You have unsaved changes" when nothing was changed** | Cosmetic | Answer `y` | ch:redox
- **Export "1,5" exports nothing ("Invalid choices: 5")** | One bad number rejects the whole list | Use only 1, 2, 3 | ch:redox
- **Import JSON: "Successfully imported 0 redox site(s)" then "No sites found in JSON file"** | The file browser also lists ProPrep session JSON files | Pick `<PDB>_redox_sites.json` or `..._redox_sites_updated.json` | ch:redox
- **Menu hints read "...or press  to return to the main menu"; 4d prompt shows "Choice  (a)"; batch help shows "Selection Syntax: :,,..."; transformer summary "state  -> HXO"** | Text in square/angle brackets is lost by the screen markup | Type `m` for main menu; the Apply choices are a/e/r/c; the selection syntax is `<site>:<state>,...` | ch:redox
- **Residue numbers in 4d progress lines and enhanced_residue_mapping.txt do not match the final PDB (e.g. HMO 105 vs 106; "A_HEC_1299" for heme 1294)** | 4d renumbers residues; reports use intermediate IDs | Trust `final_transformed_structure.pdb` and `*_redox_sites_updated.json` | ch:redox
- **Propionate residues written between crystal waters (HMO 106, HOH 107, PRD 108, HOH 109, PRD 110) although "No ID conflicts detected"** | 4d places new residues by number | tLEaP copes; check the QM/MM HIGH layer (see 5c below) | ch:redox
- **All 4d options greyed out: "Already generated microstate(s)"** | Only one generation per session | Start a new session in a new project folder | ch:redox
- **4HHB: "Required center count (expected 1, got 4)"** | Several hemes grouped into one site; every heme transformer needs one center per site | Group one heme per site | ch:redox
- **Fe–S: "'FD1' is not a valid option for site 1"** | Broken-symmetry code must match the chosen redox state (FO*, FR*, FD*) | Use a code of the right family; ProPrep then asks per option | ch:redox
- **Batch microstate code shows "OL" for no_transformation sites** | The code is only a label | Ignore | ch:redox
- **"Using RCSB PDB PDB file for transformation" printed twice** | Cosmetic | Ignore | ch:redox
- **tLEaP/ParmEd run from another AmberTools installation (e.g. a different conda environment)** | ProPrep uses whatever tleap/AMBERHOME comes first | Start ProPrep from the activated ProPrep environment; do not set AMBERHOME by hand | ch:redox, ch:forcefields

## Force Field Explorer (2b)

- **Default force-field selection is ff14SB + TIP3P although the same screen stars ff19SB and OPC and calls TIP3P "not recommended"** | Default not aligned with the recommendations | Type the numbers you want (ff19SB + OPC) | ch:forcefields
- **Out-of-range numbers or text at force-field selection are ignored silently** | No validation message | Check the "Loading:" line | ch:forcefields
- **Browse Atom Types: first pages full of "types" like -0.00006** | ff19SB CMAP grid numbers are read as atom types | Always type a filter | ch:forcefields
- **Typing `q` at a filter prompt lists residues containing "q" (e.g. QSPCFWBOX) instead of leaving** | `q` is treated as a filter text | Type `b` to go back | ch:forcefields
- **Searching `ZN` or `SZ ZN` finds nothing in the Zn(Cys)4 set** | Search is by atom type (ZM), not atom name (ZN) | Search by type | ch:forcefields
- **"Loaded:" line empty when only cofactor sets are loaded** | Cosmetic | Ignore | ch:forcefields
- **User library entries listed under "Redox Cofactors (built-in)" and reported as "bundled"; built-in row numbers shift** | User library shares category 12 | Read the row names, not the numbers | ch:forcefields
- **Custom file loader (category 13) cannot take a typed path** | That browser only navigates by number and `..` | Navigate step by step, or copy the files into the project folder | ch:forcefields

## Force Field Parameterizer (2c) and library

- **Small molecule step 4: "✗ Step 4 failed: Gaussian output not found: mk1.log" followed by a Python Traceback** | Gaussian has not been run; the traceback is printed but ProPrep has not crashed | `s` (save & quit), run `g16 mk1.gjf`, then 2c → 3 → resume. Without Gaussian, use the AM1-BCC route outside ProPrep and import | ch:forcefields
- **No AM1-BCC (or any non-Gaussian) charge option in 2c** | Every charge route needs Gaussian | antechamber/parmchk2/tleap outside ProPrep, then 2c → 2 import (Section sec:forcefields-am1bcc) | ch:forcefields
- **"Run anyway" on step 5 before step 4: "sequence item 2: expected str instance, NoneType found" + traceback** | Step 5 needs step 4's mol2 | Do not force steps; finish step 4 first | ch:forcefields
- **2c option 4 says "No parameterizations found in workspace" although a job was saved** | Status does not read saved checklist states | Use option 3 and resume | ch:forcefields
- **Several copies of one ligand (3× EDO): second copy offers to "resume" the first copy's state; answering `n` deletes it** | Copies share one folder | Parameterize one copy only | ch:forcefields
- **Modified amino acid (SEP) Gaussian inputs are plain serine with charge −2** | The de-novo route starts from the parent residue; user must build the modification | Edit the .gjf files (keep atom order) before running Gaussian; first check the Explorer — SEP is already covered by phosaa19SB | ch:forcefields
- **Mapping table maps CSD to ALA** | CSD (3-sulfinoalanine) is a cysteine derivative | Check/edit the mapping (2c → 1 → 4 → edit_aa) | ch:forcefields
- **Library "Save as a new version" makes v2 the default silently; both versions point to the same file names** | Versioning reuses files | Rename new files before importing a second version | ch:forcefields
- **Default library description says "AM1-BCC/RESP charges"** | Generic default text | Type the real method | ch:forcefields
- **"Imported parameters need a transformer" after import** | The Topology Generator reaches library entries through a transformer | Load structure, run 2a, 4d → 1 → `save` to make a pass-through linked to the entry (not tested) | ch:forcefields
- **2c Help says metal steps are "Coming soon"** | Help text out of date | Ignore; the 15-step checklist is current | ch:forcefields
- **Metal site step 12: "Total residues: 0", "Wrote PDB file: small.pdb (0 atoms, charge=0.00)"; .gjf files have no atoms** | Likely prepared_structure.pdb loses chain IDs and renumbers residues so no site residues match (observation) | No workaround in 1.23.0; use a ready-made set (zinc_cys4, hemes, Fe-S) where possible | ch:forcefields
- **Metal site step 10: "Fatal Error! bond: Argument #2 is of type String" when 2a sites with bonds exist** | Bond lines without chain IDs; disulfides deduplicated wrongly (one self-bond); metal bonds to a removed metal | Do not run 2a beforehand; let step 11 run detection. Re-running 2a without bonds did not help | ch:forcefields
- **Chloride ions classed "Unknown"/organic; supplying atomic_ions.lib renames CL to AG and tLEaP fails** | The first unit of a multi-unit library (Ag+) is used | Discard crystallisation ions in step 1 | ch:forcefields
- **Step 12 asks for a formal charge for every ligand atom (10–24 prompts)** | Ordinary ligands treated as "withheld cluster residue" | Answer 0 (−1 for thiolate S) | ch:forcefields
- **Suggested model charge always 0** | Charge suggestion unreliable | Work it out: metal + charged ligands (Fe³⁺ + 4 Cys⁻ = −1) | ch:forcefields
- **Zn-bound His (via NE2) recommended as HIE; disulfide Cys recommended as CYS** | Protonation recommendation does not consider metal binding/disulfides | Override: HID for NE2-bound His, CYX for disulfides, CYM for metal-bound Cys | ch:forcefields
- **Zn on a crystal symmetry axis: "No nearby atoms found"** | Ligands are in symmetry copies | Choose a fully coordinated metal or build the biological assembly | ch:forcefields
- **Some file browsers accept a typed path, others do not (metal step 6, Explorer category 13)** | Different browser implementations | Copy files into the project folder | ch:forcefields

## QM/MM Preparator (5c)

- **5c: "Frame extraction failed: No prmtop file selected"** | No topology yet; the main menu only mentions redox sites | Build topology in 5a first | ch:qmmm
- **ONIOM HIGH layer contains a water (WAT109) and misses a propionate (PRD106) for 1HRC** | Site residues mapped to topology residues by number; 4d placed propionates between waters | At "Add or remove whole residues?" type `-109`, `106`; always check the list by name | ch:qmmm
- **Trim list numbers residues one higher than the add/remove prompt, and offers propionates (CA/CB names) for trimming** | Inconsistent numbering; name-based CA/CB detection | Do not trim propionates; check by name | ch:qmmm
- **3-layer ONIOM: the HMO–PRD bond crossing the boundary gets no link atom and no warning** | Boundary detection misses it when the HIGH list is wrong | Fix the HIGH list first | ch:qmmm
- **ScaleCharge digits not validated (135 accepted)** | No check of the "decreasing digits" rule | Use 500 (default) unless you know otherwise | ch:qmmm
- **Re-running step 3 with another boundary treatment leaves the old ScaleCharge in the route line** | Keyword stored and not cleared | Check route preview; delete the old keyword from the .com with a text editor | ch:qmmm
- **Charge/multiplicity line "7 2 -1 1" (real doublet, model singlet) for ferric heme** | Multiplicity default 1, no parity check; step 3 printout "7 1 -1 1" differs from file | Re-run step 4 with multiplicity 2, then steps 5–7 | ch:qmmm
- **Diagnostic-file message mentions a Test keyword but the file uses OnlyInputFiles** | Message out of date | Ignore | ch:qmmm
- **Declining "Resume from saved state?" still reuses old functional, basis, job, cores, memory, ScaleCharge, keywords** | Settings kept in configuration | Use 5c option 3 (reset) first, or check every value | ch:qmmm
- **Manual boundary charges: "Total charge (7.2400) deviates from nearest integer (7)" but "Validation passed"** | Validation only warns | Keep the total charge whole | ch:qmmm
- **MEDIUM layer method cannot be chosen (fixed HF/3-21G); 3-layer charge line "7 2 0 1 0 1 0 2" looks inconsistent** | Not configurable; values not checked | Have an ONIOM expert check/edit the .com | ch:qmmm
- **ORCA: "orca_mm not found on PATH", step 2 fails with traceback; step 8 "[Errno 2] No such file or directory: ''"; no files written** | ORCA's orca_mm is needed at preparation time | Install ORCA and give its folder | ch:qmmm
- **ORCA: non-number (e.g. `abc`) at a number prompt: "✗ Step 4 failed: could not convert string to float" + traceback** | No re-ask | Re-run the step | ch:qmmm
- **Trajectory with several frames: "Frame extraction failed: cpptraj produced no output frames"; next single-frame run picks a stale PDB** | cpptraj output naming not matched | Extract one frame at a time, in a fresh session | ch:qmmm
- **Frame syntax shown as "start🔚increment"** | ":end:" turned into an emoji | Read as start:end:increment | ch:qmmm
