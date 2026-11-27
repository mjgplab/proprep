# User-facing problems mentioned by w3 (Ch. 8 Tutorial, 14 Topology, 15 Running MD, 16 Analysis)

Format: Symptom | Why | What to do | Chapter label

## Environment
- **tLEaP/sander messages show a path in another conda env (e.g. `.../envs/<other>/bin/teLeap`)** | ProPrep keeps an `AMBERHOME` that is already set, silently using another Amber installation | Start the ProPrep env in a fresh Terminal; if `echo $AMBERHOME` points elsewhere, `unset AMBERHOME` before starting ProPrep | ch:tutorial (sec:tutorial-before), ch:topology, ch:installation

## Topology Generator (5a)
- **Back at the Topology Generator menu after "Please run RedoxSite detection and try again." — no `_tleap.in` written** | Default answer to "Do you want to continue without bond definitions? [y/n] (n)" is n, which aborts | Run 5a → 2 again (settings are remembered) and type y | ch:tutorial, ch:topology
- **"✓ tLEaP completed successfully" followed by "Errors = 42"; "Results: 1/1 files processed successfully" plus "Some expected files are missing"; menu says "✓ Topology files generated (1 prmtop files)"** | Success messages are printed even when tLEaP fails; status table counts a 0-byte prmtop as "✓ exists" | Trust "Errors =" and the file size (`ls -l *.prmtop`); read leap.log | ch:topology (sec:topology-problems)
- **A good `.prmtop` becomes 0 bytes after building a second variant (e.g. implicit solvent)** | Two tLEaP inputs write the same output names (`saveamberparm`); a failed run truncates the file; no warning | Give each variant its own names on the `saveamberparm` line when editing the template | ch:topology (sec:topology-implicit)
- **"tLEaP failed (exit code: 31)" with only "Unknown residue: HOH ... -no luck" on screen** | The FATAL reason ("Atom ... does not have a type") is only in leap.log; crystal waters in an implicit-solvent build | Remove waters with 4a before implicit-solvent builds; `grep -n "FATAL\|Error" leap.log` | ch:topology
- **Implicit-solvent topology lacks GB radii** | ProPrep does not add `set default PBradii` for implicit solvent | Add the line when editing the template | ch:topology
- **"Could not find bond parameter for atom types: SH - SH", Errors = 42/60, "ParmEd validation failed: Could not identify file format" (Python traceback)** | Disulfide cysteines still named CYS (4e does not rename to CYX); CONECT records make tLEaP bond SG atoms | Rename to CYX (2a + 4d with disulfide option) AND add `bond mol.X.SG mol.Y.SG` lines in the bond editor | ch:topology (sec:topology-disulfide), ch:redox
- **CYX build reports success but ParmEd prints "Detected CYX residue with a Sulfur atom not bonded to another CYX Sulfur! Did you forget to make the disulfide bond?"; tLEaP warns "Illegal CONECT record"** | CYX names without bond lines → topology silently without disulfides | Add the bond lines (5a → 1), regenerate input (2) and topology (4) | ch:topology
- **Bond editor saves nonsense bonds (e.g. a bond literally called `n`); typed text vanishes after adding a bond** | No validation of bond text; "Press Enter to continue" swallows input; q saves (no discard) | Check the list before saving | ch:topology
- **Heme protein (1HRC) straight into 5a: "tLEaP failed (exit code: 31)", Errors = 159, "Unknown residue: HEC"; no warning beforehand** | Heme not prepared; no parameters | Prepare via 2a → 4d (heme parameter sets) first | ch:topology (sec:topology-special), ch:redox
- **5a option 3 (all microstates) returns to the menu with no files and no message** | Batch of disulfide-only microstates needs no force-field files; generator aborts silently | Use option 2 on each microstate PDB | ch:topology
- **5a option 2 blocked: "Batch microstates already generated"** | Batch metadata from 4d present | Use option 3, or work in a separate project | ch:topology
- **"Previous settings: protein: {'name': 'ff19SB', 'leaprc': ...}" — raw computer-style dictionary on screen** | Display of stored settings | Read the last "Solvation: ..." line; press Enter to keep | ch:topology
- **No warning when ff19SB is combined with TIP3P; TIP3P line reads both "Recommended" and "not recommended" with the constant-pH FF** | Labels only | Check the "Recommended" tag | ch:topology (sec:topology-water)
- **DNA/RNA/carbohydrate/small-molecule force fields loaded for a protein-only system** | Defaults are real force fields, not None | Type 1 for each one not needed | ch:tutorial, ch:topology
- **`1UBQ.pdb` changed; `1UBQ.pdb.ter_backup` appears** | TER fixer edits the PDB in place (backup kept) | Use the backup for the original | ch:tutorial, ch:topology
- **cpin file contains only LYS/TYR; ASP/GLU/HIS missing** | Topology not built with Constant pH (ff10) FF / residues not renamed AS4/GL4; no warning | Build with protein FF 8 and rename via 4e constant-pH option | ch:topology (sec:topology-cpin)
- **After making a cpin, 5b uses `_cpin.prmtop` and asks "Enable Constant pH MD? (y)"** | cpin generation switches the default structure pair | Add normal pair (2), remove cpin pair (4); answer n | ch:topology, ch:runningmd
- **Loaded prmtop/rst7 pair: 5a options 2–6 stay blocked** | Pair usable only by 5b and QM/MM | Build in 5a if PB/cpin needed | ch:topology
- **File browser jumps back to the project folder after picking the topology** | Browser does not stay in the topology's folder | Type the full path of the rst7 | ch:topology
- **PB Titrate: Enter at "Quit without saving? (n)" stays in the checklist** | Default is no | Type q then y | ch:topology

## MD Manager (5b)
- **"Error starting simulation: [Errno 2] No such file or directory: 'pmemd'", "Simulations started: 0"** | Default engine is pmemd, not installed with ProPrep; "recommended" picks pmemd.cuda even with 0 GPUs; bulk list doesn't mark missing engines | Redo 5b → 1, hardware step 3 → 1 (sander) or 2 → sander.MPI; Execute's queue menu cannot change the engine | ch:tutorial (sec:tutorial-engine), ch:runningmd (sec:runningmd-engine)
- **Review table shows original lengths (1000 cycles / 20 ps / 300 ps) after shortening** | Display not updated | Ignore; files contain edits (check simulation.mdin) | ch:tutorial, ch:runningmd
- **Mdin saved without imin=1 (or other lines) after an editor slip — no warning** | Edits not validated | Read the MDIN Parameters list after each edit | ch:tutorial, ch:runningmd
- **Workflow preview shows `pmemd` commands and `step_N.rst` names; `files` overrides break the chain** | Preview is cosmetic | Leave `files` alone; removing a step renumbers automatically | ch:runningmd
- **`info` in the workflow editor asks for a new name instead of showing info; preview's y/n prompt rejects other commands** | Command naming | — | ch:runningmd
- **"Received '', which is not one of the available options"** | Some prompts (e.g. "Enter choice:") have no default | Type the letter | ch:tutorial, ch:runningmd
- **"Load a saved protocol" always says "No saved protocols found"** | Saved protocols not listed | Import the protocol JSON | ch:runningmd
- **Deleted custom protocols still in `~/.proprep/md_workflows/user/custom/`; "? steps" shown; imported mdin saved to project `md_templates/user/custom/` not the reported location** | Library bookkeeping | Delete files by hand | ch:runningmd
- **Wizard ("Create a step from scratch") minimization contains ntb=2, ntp=1, barostat=2; sander notes "mcbarint is greater than the number of steps"** | Pressure-control defaults applied to minimization | Answer pressure questions carefully or use a predefined stage | ch:runningmd
- **DISANG import prints "Error importing restraints: 'MDRestraintModule' object has no attribute 'get_from_workspace'"** | Error message although import worked | Check with option 3 (Display) and continue | ch:runningmd
- **DISANG interactive entry: "Residue A:1 not found" / "Invalid format"** | Topology-extracted PDB has no chain IDs | Write a `.disang` file and import it | ch:runningmd
- **DISANG importer doesn't find a restraint file** | Only `*.disang` in project folder listed | Rename to `.disang`, put in project folder | ch:runningmd
- **Restraint mask generator reports "residues 1-102" for a 76-residue protein** | Ions counted as residues | Use protein residue numbers; numbering is consecutive | ch:runningmd
- **SLURM scripts request a GPU for sander; constant-pH step2.sh uses the non-cpin topology; step2.sh has a literal `--dependency=afterok:$STEP1_JOBID`** | Script generator defaults/bugs | Edit scripts before submitting; use submit_workflow.sh | ch:runningmd (sec:runningmd-cluster)
- **Two batch folders after "Write MDIN only" then Execute** | Execution makes a fresh run folder | Use the newer folder for results | ch:runningmd
- **Each stage folder's input coordinates called `<name>.rst7` though they come from the previous stage** | Naming only; chaining correct | — | ch:runningmd
- **Monitor "Show timing information": "No timing information available"** | Seen just after a stage finished | — | ch:runningmd
- **Run plans menu says "save one from Step 4"; run plans are saved at Step 5 (Hardware)** | Wording | Use Hardware step `s` | ch:runningmd

## Analysis
- **Typing 10, 20 or any list containing a 0 prints "Returning to MD Manager menu..."; Angle (10) and Density maps (20) cannot run** | Selection test treats any '0' character as Exit | Leave 10 and 20 out | ch:tutorial, ch:analysis
- **"Failed to setup action" + Python traceback, thrown back to MD Manager menu (trajectory lost)** | Invalid atom mask (e.g. CB of glycine) in dihedral analysis | Restart analysis with valid atoms | ch:analysis
- **Water RDF accepts an empty solute mask ("Solute:" blank)** | No validation | Type a mask such as `:1-76` | ch:analysis
- **SASA default region (C-alpha) gives meaningless value** | Default region | Use region 4 | ch:analysis
- **Heating stage detected as "NVT Equilibration/Production"** | Detection by ntb/ntp only | Type 2 at "Confirm simulation type" | ch:analysis
- **Exports overwrite the same `<stage>_trajectory_analysis.csv/json`** | Fixed file names per stage | Last export has everything; copy files to keep versions | ch:analysis
- **No picture files of plots** | Plots are text-only on screen | Plot from the CSV | ch:analysis
- **`find` in 5b → 4 → 2 browser accepts only one file ("2,1" rejected)** | Single-select prompt | Use `browse` at "Select simulation" (add/done) | ch:analysis
- **`proprep --analysis` asks menu mode and "Save this as your default preference?"; exiting with 0 lands in main menu** | Startup questions; exit target | Add `--full-menu`; answer n | ch:analysis
- **Log lines with dates and INFO/WARNING appear on screen** | Internal logging leaks | Ignore | ch:analysis
- **Search list shows a `0.0 KB` 1UBQ.prmtop** | Leftover of failed tLEaP run | Pick a normal-size topology (e.g. in `shared/`) | ch:analysis
