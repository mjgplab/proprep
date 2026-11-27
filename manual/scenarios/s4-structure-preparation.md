# Structure preparation (Full Menu section 4: 4a, 4b, 4c, 4e, 4f)

Scope id: s4. ProPrep 1.23.0, started with `proprep --full-menu`. The
structures were always loaded first with 1a -> 1 (Load protein structure) -> 1 (RCSB) -> 1 (PDB ID) -> Enter (asymmetric unit) -> 6 (Done).
4d (Redox Site Preparer) and 4g (Membrane Builder) are outside this scope.

Test structures and project folders (all under `manual/runs/s4/`):

| Folder | PDB | Why | Session JSON |
|---|---|---|---|
| filter1/hsg | 1HSG (HIV protease dimer + inhibitor MK1 + waters) | filter, protonation, orientation, order effects | `filter1/hsg/proprep_session_20261121_122925.json` |
| filter2/zni | 1ZNI (insulin, 4 chains, Zn2+, Cl-) | ions in the filter | `filter2/zni/proprep_session_20261121_132807.json` |
| mut1, mut2, mut3/ubq | 1UBQ (ubiquitin, clean) | mutator, mutation application, protonation without PROPKA | `mut2/ubq/proprep_session_20261121_130405.json`, `mut3/ubq/proprep_session_20261121_132217.json` |
| fix1/rh1 | 2RH1 (beta2-adrenergic receptor/T4-lysozyme chimera) | REMARK 465 missing residues, REMARK 470 missing atoms | `fix1/rh1/proprep_session_20261121_130920.json` |
| alt1/nir | 3NIR (crambin, 0.48 A, alternate locations, H atoms) | alternate locations | `alt1/nir/proprep_session_20261121_131822.json` |
| order1/rh1 | 2RH1 | order: filter, then fix, then protonate | `order1/rh1/proprep_session_20261121_133201.json` |

Transcripts are in `manual/runs/_sessions/<session>/transcript.txt`. Sessions used: s4-filter-1, s4-filter-2, s4-mut-1, s4-mut-2, s4-mut-3, s4-fix-1, s4-alt-1, s4-order-1.

Transcript note: the driver writes `transcript.txt` one line at a time, so a prompt that does not end in a newline (for example "Choose option (r):") only shows up in the transcript once the next answer has been typed. The full byte stream, prompts included, is in `raw.log` in the same folder.

---

## Menu status marks (applies to all of section 4)

- **Main menu:** `✓` means the tool can be used now. `○` means it is blocked, and a yellow `⚠` line says why, for example "⚠ Needs a loaded structure". The main menu never shows that a tool has been *done*. After filtering, protonating and orienting, 4a, 4e and 4f still show plain `✓`.
- **Inside a tool's own menu:** `●` means the step has been done in this session, `✓` means it is available, and `○` plus `⚠` means it is blocked, for example "⚠ Need to filter structure first". Source: `utils/enhanced_menu.py`, where COMPLETED = ●, AVAILABLE/READY = ✓, BLOCKED = ○.
- Before any structure is loaded, 4a, 4b, 4e and 4f say "⚠ Needs a loaded structure". 4c says "⚠ Needs an experimental structure (RCSB or local PDB; not AlphaFold)". After a load they all turn `✓`.
- Each tool menu has a "→ Suggestion:" line. Several of these contain the text "or press  to return to the main menu", with an empty gap where the key name should be (a cosmetic bug seen in every tool).

## Which structure each tool uses (the order effect)

ProPrep keeps every intermediate structure in the "workspace" and each tool picks one by its own rules. This is the most important thing for a user to understand about section 4.

| Tool | What it picks (observed; source in brackets) |
|---|---|
| 4a PDB Filter | Only one structure: uses it. Several: **asks** in an "Available Structures in Workspace" table, for example 1 Oriented, 2 Protonation-Updated, 3 Filtered, 4 RCSB PDB. |
| 4b Mutator (sequence display) | **Asks** when there are several (for example 1 Completeness-Repaired, 2 RCSB PDB). |
| 4c Structure Fixer | Uses the Filtered structure if one exists ("Using filtered structure"), otherwise RCSB or local PDB. Never AlphaFold. (`_get_priority_structure`) |
| 4e Protonation | Picks automatically: Redox-transformed > preprocessing > **Repaired** > **Filtered** > H-stripped > RCSB > Local > AlphaFold. It **never** uses an Oriented or an earlier Protonation-updated file. It printed "Using filtered structure: filtered_structure.pdb" and "Using repaired structure: .../repaired_structure.pdb". |
| 4f Orientator | Picks automatically: **Oriented** > Protonation-updated > Redox-transformed > Repaired > Filtered > ... > RCSB. It printed "Using Protonation-Updated PDB file: structure_with_md_names.pdb", and on later runs "Using Oriented PDB file: ..._oriented.pdb". |

Consequences seen in real runs:
- **Filter, then Protonate:** works. 4e used `filtered_structure.pdb`.
- **Filter, then Fix:** works. 4c said "Using filtered structure" and still found the 58 REMARK 465 missing residues of 2RH1. The header information is taken from the download, not from the filtered file, which has no header lines.
- **Fix, then Protonate:** works, but the repaired file is **renumbered from 1** (see 4c), so PROPKA tables and residue-name overrides use the new numbers. For 2RH1, 4e listed residues up to "LYS 429"/"ARG 430", while the deposited numbering runs 29-1161.
- **Protonate, then Orient:** works. 4f took the protonation-updated file and wrote `structure_with_md_names_oriented.pdb`.
- **Orient, then re-filter and re-protonate (pitfall, RUN):** after a second filter (chain B only) and a new protonation, 4f still offered the **old** oriented file, which has two chains and the old residue names, because "Oriented" outranks everything. 4e, for its part, ignored the oriented coordinates. Orientation should be the **last** preparation step, and after any change you need to orient again.
- Filtering a second time overwrote `filtered_structure.pdb`, and setting MD names again overwrote `structure_with_md_names.pdb`. Neither asked first.

---

## Feature: 4a. PDB Filter

**What it is for (plain English):** keep only the parts of the PDB file you want to simulate (chains, protein, waters, ions, ligands/cofactors) and throw the rest away.
**Needs before you start:** a structure loaded with 1a.
**Where to find it:** Full Menu `4a`. Workflow Mode: "fix" stage (section 4 STRUCTURE PREPARATION). No command-line flag.

Sub-menu: `1. Filter PDB structure`, `2. Show filtering status`, `3. Export filter selection`. Options 2 and 3 stay `○` "⚠ Need to filter structure first" until a filter has been applied.

### Scenario 1: Opening the filter with no structure loaded — RUN
- Steps: main menu `4a`. The menu shows "⚠ Needs a loaded structure" and "→ Suggestion: Needs a loaded structure. Load one via the Structure Loader." Typing `1` printed "⚠ Needs a loaded structure" and redrew the menu. Typing `1a` inside the sub-menu gave "Received '1a', which is not one of the available options". Main-menu keys do not work inside a sub-menu: use `m` first.
- Transcript: runs/_sessions/s4-filter-1/transcript.txt

### Scenario 2: Full filter of 1HSG, keeping protein and discarding waters and the ligand — RUN
- Starting point: 1HSG loaded (RCSB).
- Steps:
  1. "Enter your choice (1):" -> `1` (Filter PDB structure). The screen says "Using RCSB PDB PDB file: .../1HSG.pdb". On first use it downloads `components.cif.gz` (the chemical component dictionary), then prints a **Chain Interface Analysis**: a buried-surface-area matrix (A↔B 3610 Å²), a colour legend, a "Chain Topology" box and a "Connectivity" list.
  2. "Detect redox-active sites before filtering? [y/n] (y):" -> `n` (see Scenario 6 for `y`).
  3. List "1. Chain A - 143 residues - Interfaces with: B (3610Å²) / 2. Chain B - 183 residues ...".
  4. "Would you like to view the structure in the interactive 3D viewer? [y/n] (n):" -> `?` gave "Received '?'; please enter y or n". Then `n`.
  5. "Select chains (comma-separated indices, or 'all') (all):" -> `all`.
  6. A **Chain Composition** table per chain: Component Type (Protein / Water / the ligand's chemical name), residue names, counts, and a Redox Sites column.
  7. "Do you want to proceed with filtering? [y/n] (y):" -> `y`.
  8. For each chain and component in turn: "Chain A - Protein [r] Retain entire component [s] Select specific residues [d] Discard entire component Choose option (r):". Answers given: A Protein `r`, A Water `d`, B Protein `r`, B ligand (shown as "Chain B - ;N-[2(R)-HYDROXY-1(S)-INDANYL]-5-[(2(S)-TERTIARY", the MK1 chemical name) `d`, B Water `d`. An invalid key (`x`) gave "Received 'x', which is not one of the available options".
  9. A "Filter Selections" box (Chain A: Amino_acid 99 residues; Chain B: Amino_acid 99 residues), then "Do you want to apply these filters? [y/n] (y):" -> `y`.
- What the user sees: "Filtered structure automatically saved to: filtered_structure.pdb" and "Structure filtering completed successfully." The sub-menu now shows `1. ● Filter PDB structure`, `2. ✓`, `3. ✓`.
- Files: `filter1/hsg/filtered_structure.pdb` (protein only, 1514 atoms, ATOM/TER/END records only, **no header/REMARK/SEQRES lines**). The first download also left `1HSG.pdb`.
- Transcript: runs/_sessions/s4-filter-1/transcript.txt

### Scenario 3: Keeping only some waters (water analysis) — RUN
- In Scenario 2's first attempt, at "Chain A - Water" -> `s`. A "Water Analysis Overview" panel lists 8 analysis methods, followed by a menu:
  `1 Distance to nearest metal ion, 2 Hydrogen bond analysis, 3 B-factor/thermal information, 4 Water burial: solvent-accessible area and enclosure, 5 Interface proximity, 6 Multi-radius atom-count profiling, 7 Directional atom-count analysis (8 sectors), 8 Water network analysis (clusters and chains), 9 Return with no waters selected and proceed with component filtering`
  Prompt: "Select metrics to display (comma-separated or single choice) (6):".
- `2,4` -> "Received '2,4', which is not one of the available options". `2 4` was rejected in the same way. **Bug:** the prompt says comma-separated is allowed, but only one number is accepted (`pdb_filter_worker.py` line ~2301 limits the choices to "1"-"9"). Single `4` worked.
- "Would you like to modify any parameters? [y/n] (n):" -> `y` shows "1. Probe radius (1.4 Å) 2. Occluding atom types (protein,hetero,metal) 3. Done with changes" -> `3`, then `n`.
- A results table "Water Burial (44 waters: 1 enclosed, 0 buried, 0 clashing)" with Residue #, SASA, Covered %, Access (bulk/enclosed), Nearest atom and Category (Clash/Enclosed/Buried/Exposed), plus a legend. HOH 324 was the only "Enclosed" water (next to ARG87).
- "Analysis Options: 1 Adjust analysis parameters and re-run, 2 Return to metrics selection menu, 3 Proceed with water selection, 4 Cancel and return to component menu" -> `3`.
- "Select waters by residue numbers (all):" (accepts `215,847,102`, ranges `100-150`, `all`, `none`). These are **PDB residue numbers**, not table row numbers.
  - `324,999` -> "Invalid residue numbers: [999] Available: 305 to 575 Please try again".
  - `324,305-310` -> "Invalid residue numbers: [306, 308, 310]". A range is rejected if any number inside it is not a water.
  - `324,305,307` -> "Keep 3 waters, discard 41 waters. Keeping: HOH305, HOH307, HOH324 Confirm selection? [y/n] (y):" -> `y`.
- Choosing `9` at the metrics menu (session s4-filter-1, chain B water) keeps **no** waters: the summary showed "Water: 0 residues".

### Scenario 4: "Select specific residues" for protein, with range or out-of-range input (crash) — RUN, FAILED
- At "Chain B - Protein" -> `s`. A table lists every residue with a "Select" index (1..99), Record Type, Residue Name and Residue Number. Prompt: "Select residues (comma-separated indices, 'all', or 'none') (all):".
- `1-50` -> **crash**: "Error executing PDB Filter.filter_pdb_structure: invalid literal for int() with base 10: '1-50'", followed by a full Python traceback. The program then drops back to the PDB Filter menu, and **every answer given so far in this filter run is lost**.
- `200` (beyond the list) -> **crash**: "Error executing PDB Filter.filter_pdb_structure: list index out of range".
- Only plain comma lists of the **Select index** work, for example `1,2,3`. Note that this screen uses table indices while the water screen uses residue numbers.

### Scenario 5: Declining "apply these filters" — RUN
- At "Do you want to apply these filters? [y/n] (y):" -> `n`. The program goes back to the filter menu without any message. Nothing is written, and the previous `filtered_structure.pdb` is kept.

### Scenario 6: Detecting redox sites first, plus protected components and a chain subset — RUN
- After an earlier filter, option `1` showed the "Available Structures in Workspace" table (Oriented / Protonation-Updated / Filtered / RCSB PDB) with the prompt "Select structure to use (1):" -> `4` (the original RCSB file).
- "Detect redox-active sites before filtering? [y/n] (y):" -> `y`, then "Import redox sites from previous JSON file? [y/n] (n):" -> `n`. This runs the full Redox Site Detector dialogue, which is covered by another scope: metals to scan, coordinating elements, radius, disulfide cutoff, non-standard residues y, redox amino acids n, center selection and grouping, refinement. Defaults were accepted, with `finish` at "Select residues to add:", `back` and then `5` (Continue). Result: "✓ Found 1 redox-active site(s)" (MK1 B:902 flagged as an "organic cofactor"). The chain list then shows "Chain B ... - [✓ 1 redox center(s)]".
- On a later filter run in the same session the redox question was **not asked again**; the stored sites were reused.
- "Select chains" -> `3` -> "Invalid chain index: 3. Must be between 1 and 2 / Invalid input: 3. Please enter comma-separated numbers or 'all'". Then `2` (chain B only).
- A component that contains a redox site gets a **different menu**: "⚠️ This component contains redox site residues: • Site site_1 (): MK1902 / [r] Retain entire component / [c] Keep redox site residues / [s] Select specific residues / [o] Override redox protection and discard (⚠️ breaks redox sites)". `c` -> "Keeping 1 redox site residues".
- Result: `filtered_structure.pdb` was overwritten (chain B protein + MK1, no waters) without asking.

### Scenario 7: Ions (1ZNI) — RUN
- 1ZNI: 4 chains. The interface analysis showed A↔B 1282, B↔D 1280 and C↔D 1619 Å². Chains `1,2` selected.
- Ions appear as their own components with chemical names: "CHLORIDE ION" (CL ×3) and "ZINC ION" (ZN ×2) in chain B. Answers: A Protein r, A Water d, B Protein r, B ZINC ION `s` -> table of ZN 31 and ZN 32 -> `1` (keep the first Zn only), B CHLORIDE ION `d`, B Water `r`.
- Result: `filter2/zni/filtered_structure.pdb` with chains A and B, 1 ZN, 20 HOH and no CL.
- The filter does **not** remove alternate locations: 10 altloc atoms (for example HIS B10 CG A/B) were still in the filtered file. Altlocs are handled in 4c (see the limitation there).

### Scenario 8: Show filtering status and export — RUN
- `2` -> table "Filtered Structure Status" (Chains 2, Residues (Standard) 198), "Original Structure Summary" (326 residues), "Filtered Structure Summary" (198), "Removed Components: Chain A: 44 residues removed, Chain B: 84 residues removed, Hetero groups: MK1, HOH".
- `3` -> "Enter output filename for statistics (filter_statistics.json):" -> Enter -> "Filter statistics exported to: filter_statistics.json". The JSON holds the original and filtered composition, what was removed, and the kept residue numbers per chain. File: `filter1/hsg/filter_statistics.json`.

### Options reference (4a)
| Prompt | Choices | Meaning |
|---|---|---|
| Select structure to use | numbers in the table | only shown when the workspace holds several structures |
| Detect redox-active sites before filtering? | y (default) / n | y runs the Redox Site Detector first, so cofactors and metals are protected from deletion |
| Import redox sites from previous JSON file? | y / n (default) | reuse a saved detection |
| View the structure in the interactive 3D viewer? | y / n (default) | opens the web viewer at http://localhost:8765/viewer |
| Select chains | `1,2`, `all` (default) | chain indices from the list |
| Proceed with filtering? | y (default) / n | |
| Component (normal) | r retain (default) / s select specific / d discard | per chain, per component type |
| Component containing a redox site | r / c keep only the redox-site residues / s / o override protection and discard | |
| Select residues (non-water) | `1,3,5` table indices, `all`, `none` | ranges and out-of-range indices crash the tool |
| Water metrics | 1-8 single choice, 9 = keep no waters | |
| Water selection | PDB residue numbers, ranges, `all`, `none` | every number must exist |
| Apply these filters? | y (default) / n | n discards the choices silently |
| Export filename | default `filter_statistics.json` | |

---

## Feature: 4b. Amino Acid Mutator

**What it is for (plain English):** write down point mutations (for example Lys48 to Arg), or swaps to modified or non-natural residues (for example Ser to phosphoserine). The mutator only **records** them. They are actually applied later by 4c Structure Fixer.
**Needs before you start:** a loaded structure. Standard mutations also need a working **MODELLER licence**. Non-standard ones do not, but their residue must exist in your force field, because tLEaP builds the missing atoms.
**Where to find it:** Full Menu `4b`. Workflow Mode "fix" stage. No CLI flag.

Sub-menu: `1. Display protein sequences with numbering`, `2. Specify amino acid mutations` (blocked: "⚠ Need to display sequences first"), `3. Show current mutation status` (blocked: "⚠ Need to specify mutations first").

### Scenario 1: Specifying before displaying the sequence — RUN
- `2` -> "⚠ Need to display sequences first" and the menu redraws. You must run `1` first.

### Scenario 2: Display sequences — RUN
- `1` -> "Using RCSB PDB structure object", then "Chain A:" with a ruler every 10 residues: "1-60 MQIFVKTLTG..." and "61-76 IQKESTLHLVLRLRGG". "Found 1 protein chain(s) with 76 total residues".
- Confusing: the ruler restarts at "10" on each line (the "10" above line 61-76 means residue 70).
- When a repaired structure also exists, it asks "Select structure to use (1):" with a table (1 Completeness-Repaired, 2 RCSB PDB).

### Scenario 3: Valid standard and non-standard mutations, plus invalid inputs — RUN
- `2` shows the format help: "Format: CHAIN:POSITION:FROM->TO", with examples for standard ("MODELLER will model the change") and non-standard mutations ("A:154:PHE->CNF ... You'll select which atoms to keep ... TLEaP will build missing atoms"). Commands: `done`, `clear`, `help`. The 20 standard codes are listed, and "Any other 3-letter code will be treated as non-standard".
- `help` -> "Mutation Help" panel: the two workflows, notes ("FROM residue must match", "Chain IDs are case-sensitive", "enter them one at a time") and common non-standard residues (CNF, SEP, TPO, PTR, MLY).
- Invalid inputs and messages (all re-prompt; nothing crashes):
  - `K48R` -> "Invalid mutation format: Invalid format. Use CHAIN:POSITION:FROM->TO (e.g., A:123:PHE->CNF)"
  - `A:48:ALA->GLY` -> "Invalid mutation format: Position 48 is LYS, not ALA"
  - `A:200:ALA->GLY` -> "Invalid mutation format: Position 200 not found in chain A"
  - `B:48:LYS->ARG` -> "Invalid mutation format: Chain B not found in structure"
  - `A:48:LYS->XX` -> invalid format (codes must be 3 letters)
  - second mutation at a position already used (`A:48:LYS->GLU`, even `A:48:LYS->LYS`) -> "Mutation already specified for chain A position 48"
- `a:48:lys->arg` (lower case) -> "✓ Added standard mutation: Chain A, Position 48, LYS->ARG". It was accepted even though the help says chain IDs are case-sensitive.
- `A:63:LYS->ALA` -> added as standard.
- Non-standard `A:65:SER->SEP` -> table "Atoms in SER 65" (index, name, element, coordinates) and "Preset Options: backbone (N, CA, C, O) / bb+cb / all / custom / cancel". `custom` -> "Enter atom indices (comma-separated) or names" -> `N,CA,C,O,CB,OG` -> "✓ Added non-standard mutation ... Keeping 6 atoms".
- Non-standard `A:4:PHE->CNF` with `bb+cb` -> "Keeping 5 atoms: N, CA, C, O, CB". With `backbone` -> "Keeping 4 atoms: N, CA, C, O".
- `clear` -> "All mutations cleared." **No confirmation is asked.**
- `done` -> summary tables "Standard Mutations" (with a change-type column, for example "positive→hydrophobic") and "Non-Standard Mutations" (atoms to keep), then "✓ Automatically stored mutations in workspace! Mutations will be applied when you run Structure Fixer → Process structure. Note: Non-standard mutations will be applied after MODELLER (or without it)".
- `3` (status) -> "Stored ... mutations (awaiting Structure Fixer)" tables.
- Transcript: runs/_sessions/s4-mut-2/transcript.txt (s4-mut-1 also shows the invalid-input tests; that session's program died mid-run, see Notes).

### Scenario 4: Applying the mutations (in 4c, without MODELLER) — RUN
- 4c menu: "2. ✓ Apply repairs/mutations/caps to structure ⚠ MODELLER off: capping only", "→ Suggestion: Found 2 pending mutations, but MODELLER is required to apply them".
- After analysis, `2` -> "MODELLER not available — gap filling and mutations are disabled. Terminal capping (ACE/NME) is still available". When the user then chose caps (`both`), the **non-standard** mutation was applied as well: "→ Mutating A:4 PHE → CNF / Mapped residue: A:4 → A:5 / Removed 7 atoms / Renamed residue to CNF ✓". The **standard** mutations (K48R, K63A) were **not** applied, with no further message, and they stay "ready" in 4b status.
- Output `mut2/ubq/repaired_structure.pdb`: ACE is residue 1, so every residue number shifts by +1 (Phe4 became CNF 5, Lys48 is still LYS 49). Waters are kept.
- **Bug (RUN, session s4-mut-3):** with one non-standard mutation pending and the cap question answered `s` (skip), ProPrep printed "No terminal caps selected." and wrote **nothing**, so the mutation was silently not applied. The same analysis summary listed "Pending Non-Standard Mutations 1" and then said "✓ No structural issues found! ... • No pending mutations". To get a non-standard mutation applied without MODELLER, you currently have to add at least one terminal cap.
- After the repair, 4b option 1 on the repaired structure showed "2-62 MQIVKTL..." with the CNF residue **missing from the one-letter sequence** (75 residues instead of 76), and ACE/NME not shown.

### Options reference (4b)
| Prompt | Input | Meaning |
|---|---|---|
| Enter mutation | `CHAIN:POS:FROM->TO` | standard if TO is one of the 20 amino acids (needs MODELLER), otherwise non-standard |
| | `done` / `clear` / `help` | finish / delete all pending (no confirmation) / help |
| Select preset or 'custom' (non-standard) | backbone / bb+cb / all / custom / cancel | which atoms of the old residue to keep; tLEaP rebuilds the rest |
| Enter selection (custom) | indices `1,2,3` or names `N,CA,C,O,CB` | |

---

## Feature: 4c. Structure Fixer

**What it is for (plain English):** checks an experimental structure for missing residues (unresolved loops or termini), missing atoms, alternate conformations and pending mutations, and repairs what it can. Real gap and loop repair and standard mutations need MODELLER. Without MODELLER it can only add ACE/NME terminal caps and apply non-standard mutations.
**Needs before you start:** an RCSB or local PDB file (not AlphaFold). A MODELLER licence key for repairs.
**Where to find it:** Full Menu `4c`. Workflow Mode "fix" stage.

MODELLER status here: at every start-up ProPrep prints "MODELLER not found. Structure repair (missing loops/residues) will be unavailable. To install: conda install -c salilab modeller ... set your license key: export KEY_MODELLER='your_license_key_here' Or edit the config.py file shown during conda installation." MODELLER **is** installed here and only the key is invalid, so "not found / To install" is misleading.
**How to add the key (from source `utils/modeller_license.py`; NOT done here):** either set `KEY_MODELLER` in the environment before starting ProPrep, or put the key on the first non-comment line of `~/.proprep/modeller_key` (for example `mkdir -p ~/.proprep && echo 'YOUR_KEY' > ~/.proprep/modeller_key`). The environment variable wins over the file, and a runtime key wins over the key baked in at install time. The key is free for academic use from https://salilab.org/modeller/registration.html.

Sub-menu: `1. Analyze structure completeness`, `2. Apply repairs/mutations/caps to structure` (with "⚠ MODELLER off: capping only"), `3. View structural issues report`, `4. Export results to JSON` (3 and 4 are blocked until analysis has run). Option 2 is also blocked until analysis has run ("⚠ Need to analyze structure first"), except when mutations are pending (seen in s4-mut-2).

### Scenario 1: Analyze 2RH1 with REMARK 465 / REMARK 470 / AltLoc — RUN
- `1` -> "Using RCSB PDB structure object". Three questions follow:
  - "Missing Residues Detection: 1. REMARK 465 records 2. SEQRES comparison 3. FASTA comparison 4. Sequence gaps 5. Skip Choose (1):" -> `1`
  - "Missing Atoms Detection: 1. REMARK 470 records 2. Template comparison 3. Skip Choose (1):" -> `1`
  - "Alternate Locations Detection: 1. AltLoc identifiers 2. Occupancy values 3. Skip Choose (1):" -> `1`
- What the user sees: a "Structure Completeness Summary" table: Missing Residues remark_465 **58** (chain A), Missing Atoms remark_470 1, Alternate Locations 0, Pending mutations 0/0. Then a colour-coded "Chain A Sequence" (green = missing, cyan = missing atoms, yellow = altlocs, magenta/orange = mutations, '-' = not in sequence).
- "Show detailed reports? [y/n] (n):" -> `y` -> "Missing Residues - remark_465" (positions -6..28 at the N-terminus and 343..365 at the C-terminus) and "Missing Atoms - remark_470".
- **Bug:** the REMARK 470 table showed "ASP 29 | Atom: A", which is the chain letter. The file says `ASP A 29 CG OD1 OD2`.
- The engineered deletion in the chimera (residues 231-262 replaced by T4-lysozyme 1002-1161) is shown as a long run of '-' between 230 and 263 and between ~365 and 1001, because the sequence view spans the whole numbering range. These are not missing residues.
- Suggestion line: "Found 58 missing residues, 1 missing atom. View details with option 3 or export with option 4 (MODELLER needed for repairs)".
- Time: a few seconds.

### Scenario 2: SEQRES comparison plus template comparison plus occupancy (2RH1) — RUN
- Choices `2`, `2`, `2`. "Fetching templates for 29 residue types..." took about 1-2 minutes (it downloads residue templates).
- Result: seqres_comparison 58 (matches REMARK 465), template_comparison **463** missing atoms (chains A, B), occupancy_value 0.
- The 463 consist of the real ASP29 CG/OD1/OD2, atoms missing from the lipids/ligands (12P 416, PLM 415 O1, GLC 2 O1), and about 440 **false "OXT" entries, one for every residue** (bug, see Scenario 4).
- **Display bug:** the coloured sequence placed the C-terminal segment "RRSSLKAYGNGYSSNGNTGEQSG" at positions ~471-490 instead of 343-365.

### Scenario 3: FASTA comparison and sequence gaps (2RH1) — RUN
- FASTA (`3`, with atoms and altlocs skipped): "Missing Residues fasta_comparison 0", then "✓ No structural issues found! Structure is complete with: No missing residues, No missing atoms, No alternate locations...". This is wrong: the atom and altloc checks were skipped, and REMARK 465 lists 58 missing residues. Likely cause (from source): the RCSB FASTA header for a one-chain entry reads "|Chain A|", but the parser only accepts "|Chains ...|", so no FASTA sequence is found and the result is silently 0.
- Sequence gaps (`4`): "Missing Residues sequence_gap **691**". Every numbering jump is counted, including the engineered 231-262 deletion and the jump from 365 to 1002 for T4-lysozyme, and the sequence view fills them with 'X'. This is misleading for chimeras and constructs with non-continuous numbering.

### Scenario 4: Template comparison on a complete protein (1UBQ) — RUN
- Choices `1` (REMARK 465), `2` (template), `1`. "Fetching templates for 19 residue types..." then Missing Atoms template_comparison **75**. Details: every residue from MET 1 to GLY 75 is listed as missing "OXT". **Bug:** OXT exists only on the C-terminal residue, so these are false positives. REMARK 470 (option 1) is the safer choice for RCSB files.

### Scenario 5: Alternate locations (3NIR) — RUN, limitation
- Loading 3NIR first asked "Hydrogen Atom Removal ... This structure contains 315 hydrogen atoms ... Remove hydrogen atoms? [y/n] (y):" -> `y` -> "H-stripped structure saved to: 3NIR_noH.pdb Removed 315 hydrogen atoms". This is part of 1a, but note that 18 H atoms (on altloc A of PRO22) remained in `3NIR_noH.pdb`. 4c then said "Using RCSB PDB structure object" (the original download, not the H-stripped file).
- Analysis (1,1,1): Alternate Locations altloc_identifier **32** (13 amino acids, 2 ethanols, 17 waters). Details table "Residue | Position | Alternates | Occupancies", for example "TYR 29 A, B, C A:0.42 B:0.30 C:0.28", "THR 2 A:0.82 B:0.26".
- `2` (Apply) -> straight to "Terminal Capping" -> `c` -> "Inserted NME ... ✓ Structure repair completed! Output: repaired_structure.pdb".
- **Limitation/bug:** without MODELLER **no alternate-location picker appears**, and `alt1/nir/repaired_structure.pdb` still contains **313 altloc atoms** (A and B copies). The source has an "Alternate Location Selection" picker ("doesn't require MODELLER" according to its comment), but the no-MODELLER branch returns before reaching it (`structure_completeness.py` `_unified_repair_workflow`). So without a MODELLER key, 4c cannot choose altlocs. **SOURCE-ONLY**, with MODELLER: the picker shows each residue's alternates and occupancies, asks which letter to keep, and fills partial alternates from the best other alternate.

### Scenario 6: Apply repairs without MODELLER (caps only) — RUN
- `2` -> "MODELLER not available — gap filling and mutations are disabled. Terminal capping (ACE/NME) is still available; install MODELLER for full repairs." (printed twice), then "═══ Terminal Capping (ACE/NME) ═══ ACE caps the N-terminus, NME caps the C-terminus. Cap atoms get placeholder geometry that tLEaP finalizes".
- Per protein chain: "Chain A — N-term: ASP29, C-term: TYR1161 / n N-terminus only / c C-terminus only / both / s Skip this chain (default) / Cap which terminus of chain A? (s):".
  - `s` -> "No terminal caps selected." and nothing is written.
  - `n` -> "Inserted ACE at A:28 (before 29) ✓ Caps added successfully ✓ Structure repair completed! Output: repaired_structure.pdb".
  - `both` (1UBQ) -> ACE + NME.
  - `c` (3NIR) -> NME only.
- Files: `repaired_structure.pdb` (+ an empty `modeller_work/` folder). In the repaired 2RH1 file the protein is **renumbered from 1** (ACE 1, ASP 2 ... TYR 506), and the original numbers 29-1161 are lost. All HETATM groups and waters are kept. Missing side-chain atoms (ASP29 CG/OD1/OD2) are **not** rebuilt. Missing loops are not built either, so the ICL3 junction and the missing termini remain as chain breaks. (Missing heavy atoms are normally added later by tLEaP; that is the Topology Generator scope.)
- What CAN be fixed without MODELLER: terminal ACE/NME caps and non-standard mutations (only when at least one cap is also added; see 4b Scenario 4). What CANNOT: missing residues/loops, missing side-chain atoms (in 4c), standard mutations, and altloc choice (blocked by the early return described above).

### Scenario 7: View report and export — RUN
- `3` reprints the detailed tables.
- `4` -> "Output filename (completeness_results.json):" -> Enter. Writes `mut2/ubq/completeness_results.json`, which lists missing residues and atoms by method and chain (including the false OXT entries).

### Options reference (4c)
| Prompt | Choices | Meaning / advice |
|---|---|---|
| Missing residues | 1 REMARK 465 (default) / 2 SEQRES / 3 FASTA / 4 numbering gaps / 5 skip | 1 or 2 are reliable for RCSB files; 3 returned 0 for a single-chain entry; 4 over-counts engineered constructs |
| Missing atoms | 1 REMARK 470 (default) / 2 template comparison / 3 skip | 2 downloads templates (slow) and reports OXT on every residue |
| Alternate locations | 1 AltLoc identifiers (default) / 2 occupancy values / 3 skip | |
| Show detailed reports? | y / n (default) | |
| Cap which terminus | n / c / both / s (default) | per chain, asked when MODELLER is off |
| Output filename (export) | default `completeness_results.json` | |

---

## Feature: 4e. Protonation State Analyzer

**What it is for (plain English):** predicts how acidic and basic side chains (Asp, Glu, His, Cys, Tyr, Lys, Arg, termini) are charged at your pH, using PROPKA or textbook pKa values. It then writes the Amber residue names (for example ASH, GLH, HID/HIE/HIP, CYM, LYN) into a new PDB file for MD, or the titratable names for constant-pH MD (AS4, GL4, HIP).
**Needs before you start:** a loaded structure. Filter or fix first if you want those changes included (see the order table above). PROPKA (propka3) is installed.
**Where to find it:** Full Menu `4e`. Workflow Mode "fix" stage.

Sub-menu: `1. Analyze protonation states (PROPKA)`, `2. Calculate ideal charge capacitance`, `3. Set MD residue names`, `4. View analysis results`, `5. Export analysis results` (2-5 are blocked: "⚠ Need to analyze protonation states first").

### Scenario 1: PROPKA analysis of the filtered 1HSG dimer at pH 7 — RUN
- `1` -> "Using filtered structure: filtered_structure.pdb", "No residues excluded from protonation analysis".
- "Enter pH value for analysis (7.0):" -> `abc` -> "Invalid pH value, using default of 7.0". It does **not** ask again.
- A "Protonation Probability Threshold" explanation (0.5 default = more likely state; >0.5 conservative; <0.5 liberal; "consider constant pH MD" for residues near their pKa). "Enter protonation probability threshold (0-1) (0.5):" -> `20` -> "Invalid threshold, using default of 0.5".
- "Run PROPKA for improved pKa prediction? [y/n] (y):" -> `y`. Spinners "Running PROPKA..." for a few seconds.
- Table "Protonation States at pH 7.0 (PROPKA pKa)": Index, Chain, Residue, pKa, ΔpKa (shift from textbook), Desolv, Protonated (True/False), Probability (fraction protonated), Charge. 48 titratable groups. "Net Charge: 5.13". Catalytic aspartates: ASP A25 pKa 5.07 (deprotonated), ASP B25 pKa 9.06 (protonated, ASH).
- "Enter indices to visualize (e.g., 1,3,5-7), or Enter to skip ():" -> `6` -> "ASP A:25 (predicted pKa 5.07, Δ +1.42 vs textbook) Buried 100 %; desolvation +2.99 (raises pKa)" and a per-residue rationale table: Contribution / Type (sidechain, backbone, coulombic) / Partner / Effect ("-0.89 sidechain ASP 25 B lowers (↓)", ...). It also opens the web viewer ("✓ Structure viewer launched at http://localhost:8765/viewer ... magenta = residue, green = lowers pKa, orange = raises pKa").
- "Would you like to generate a pH titration profile? [y/n] (n):" -> `y`. Min pH (0), max (14), step (0.5) all left at their defaults. Prints a net-charge-vs-pH table, an ASCII curve and "Isoelectric point (pI): 9.37".
- "Would you like to do interactive pH/charge lookups?" -> `y`: `5.5` -> "At pH 5.50: charge = 7.95"; `0e` -> "At charge 0.00: pH = 9.37"; `done`.
- "Save pH profile to file? [y/n] (y):" -> `y`, filename default `ph_profile.csv` -> "Saved pH profile to ph_profile.csv".
- Files: `filter1/hsg/ph_profile.csv`, `filter1/hsg/propka_results/filtered_structure_chains_A_B_pH7.0.pdb` and `.pka` (raw PROPKA output).

### Scenario 2: Set MD residue names, fixed pH, with manual overrides and His tautomers — RUN
- `3` -> "Is this for constant pH simulation? [y/n] (n):" -> `n`.
- "Set all residues to standard protonation states (ASP/GLU deprotonated, CYS/LYS/ARG/TYR protonated, HIS->HIE, heme PRN->PRD)? [y/n] (n):" -> `n` (use PROPKA).
- One table per residue type (Chain, ResID, pKa, State at pH, Probability, Recommended), followed by the possible names and "Use recommended names for all ASP residues? (y):".
  - ASP -> `n` -> syntax help: "Single name: ASH (applies to all) / Overrides: 94,119:ASP 122:ASH (applies to listed resids in every chain) / Chain-specific: A:94:ASP / Press Enter to use all recommendations". Typed `A:25:ASH B:25:ASP` -> "A:25 → ASH, B:25 → ASP". This swaps which catalytic Asp is protonated.
  - CYS (CYS = SH neutral, CYM = S⁻) -> `n` -> `XYZ` -> "Skipping invalid syntax: XYZ (expected resid:NAME)". It **moves on without re-asking**, so the recommendations are used.
  - GLU (GLU/GLH) -> `n` -> `21:LYS` -> "Invalid name 'LYS' - valid choices: GLH, GLU". It also moves on without re-asking.
  - HIS: "HIP - Fully protonated (both ND1 and NE2) / HIE - Epsilon nitrogen protonated (NE2) / HID - Delta nitrogen protonated (ND1) / Note: HIE is generally more common than HID". Recommended HIE for both. `n` -> `A:69:HIP B:69:HID`.
  - LYS (LYS/LYN), TYR: Enter (accept). ARG is never asked; it is always ARG. The TYR screen lists no alternative names (in constant-pH mode TYM is mentioned).
- "Set MD names for 44 residues / Summary of assignments: ARG 8, ASH 1, ASP 7, CYS 4, GLU 8, HID 1, HIP 1, LYS 12, TYR 2 / Update structure with these residue names? [y/n] (y):" -> `y` -> "Using analyzed structure for MD residue name update: filtered_structure.pdb / Enter output filename for updated structure (structure_with_md_names.pdb):" -> Enter -> "Modified 44 residue names in PDB file ... Saved updated structure to structure_with_md_names.pdb / ✓ Updated 3 residue names in protonation_results".
- File: `filter1/hsg/structure_with_md_names.pdb`, which contains ASH A25, HIP A69 and HID B69. Only residue **names** change; no hydrogens are added (tLEaP adds them later).
- Confusing: option 4 (View) afterwards shows the new names (ASH 25, HIP 69, HID 69), but the Protonated, Probability and Charge columns and "Net Charge: 5.13" still reflect PROPKA, not the override. For example "ASH 25 ... False ... -0.99".

### Scenario 3: Constant-pH naming (chain B only) — RUN
- `3` -> `y` (constant pH) -> "Note: Excluding 6 non-titratable residues from constant pH simulation. Excluded residue types: ARG, C_TERM, N_TERM / Note: Terminal residues are excluded".
- For each type: a table, then "Select ASP residues to TITRATE (will use AS4 name): Enter indices (e.g., 1,3 or 1-3), 'all', or 'none'". Then fixed names for the rest ("Use recommended names for all non-titratable ASP residues? (y)").
  - ASP `1` -> "B:25 → AS4 (titratable)". GLU `none`. HIS `all` -> "B:69 → HIP (titratable)". CYS, TYR, LYS: Enter (none).
- "Constant pH titration data saved to workspace / Set names for 18 residues / Update structure...? (y)". The file name was left at the default and it **overwrote** `structure_with_md_names.pdb` without asking. The file contains AS4 and HIP residues.
- Order effect seen: with chain A removed by the filter, PROPKA gave ASP B25 a pKa of **3.8** (deprotonated), against 9.06 in the dimer.

### Scenario 4: Textbook pKa (no PROPKA), pH 4.5, threshold 0.7 (1UBQ) — RUN
- "Using RCSB PDB structure: .../1UBQ.pdb", "No filter selections found - analyzing all chains: A". pH `4.5`, threshold `0.7`, PROPKA `n`.
- Table "Protonation States at pH 4.5 (Textbook pKa)" with no ΔpKa/Desolv columns: GLU 4.25, ASP 3.65, HIS 6.00, LYS 10.53, ARG 12.48, TYR 10.07, N_TERM 8.00, C_TERM 3.10. Net Charge 3.79.
- Visualize `4` -> "No PROPKA determinants (pKa near its model value, or default value used)".
- Capacitance (option 2) also reports "pKa source: Textbook values (Lehninger - not adjusted for structural environment)".

### Scenario 5: Standard states for everything (1UBQ) — RUN
- `3` -> `n` -> "Set all residues to standard protonation states ...?" -> `y` -> "Assigned standard protonation states (PROPKA recommendations ignored). Set names for 24 residues" -> `y` -> filename `ubq_standard.pdb`. File: `mut3/ubq/ubq_standard.pdb`.

### Scenario 6: Ideal charge capacitance — RUN
- `2` -> explanation panel (total vs specific capacitance, maximum 0.25 per residue at pH = pKa) -> "Enter pH for capacitance calculation (7.0):". A table per residue type: Count, pKa Range, Avg pKa, Total C, C/Residue. Then "Summary & Key Insights" (1HSG: total 0.6209, specific 0.0129, "Significant pKa shifts from textbook values: CYS (Δ=2.5)").
- "Would you like to generate a capacitance pH profile?" -> `y` (default range) -> ASCII plot "Peak capacitance: 2.5587 at pH 4.00" (1UBQ) -> lookups `n` -> "Would you like to save the capacitance profile? [y/n]:" (no default shown) -> `y` -> `capacitance_profile.csv` (`mut3/ubq/`).

### Scenario 7: Export — RUN
- `5` -> "Select output format (text):". The valid choices (`text`, `csv`) are **not shown**. `zzz` -> "Received 'zzz', which is not one of the available options", and the message still does not list them.
- `csv` -> `protonation_analysis.csv` (Residue_Key, Chain, Number, Type, pKa, Protonated, Probability, Charge).
- `text` (default) -> `protonation_analysis.text` (note the unusual extension), which has a header (pH, net charge, threshold) and a tab-separated table that includes the MD names.

### Scenario 8: Redox-site exclusion prompt — RUN
- After redox sites had been detected (4a Scenario 6), `1` showed a "Redox Sites Detected" table (site_1, chain B, MK1 902) and "Enter indices to EXCLUDE from protonation analysis (e.g., 1,3,5-7 or 'all') ():" with a note that metal-coordinating residues get strongly shifted pKa values. Enter = include all.
- **SOURCE-ONLY:** if redox "microstates" exist (made by 4d Redox Site Preparer), 4e first offers to analyse several microstates in a batch with one set of pH/threshold/PROPKA settings. Not reachable in this scope.

### Options reference (4e)
| Prompt | Choices | Meaning |
|---|---|---|
| pH value | number (default 7.0) | invalid input silently becomes 7.0 |
| Probability threshold | 0-1 (default 0.5) | invalid input silently becomes 0.5 |
| Run PROPKA? | y (default) / n | n = textbook pKa values |
| Indices to visualize | `1,3,5-7` / Enter | prints the PROPKA rationale table and opens the web viewer |
| pH titration profile | y / n (default); min 0, max 14, step 0.5 | net charge curve, pI, lookups (`7.0` or `-1.5e`), CSV |
| Constant pH simulation? | y / n (default) | y = choose residues to titrate (AS4, GL4, HIP, ...) |
| Set all to standard states? | y / n (default) | y ignores PROPKA (ASP/GLU⁻, HIS→HIE, ...) |
| Use recommended names for X? | y (default) / n | n = type one name for all, `resid:NAME`, or `CHAIN:resid:NAME`; invalid entries are skipped without re-asking |
| Output filename | default `structure_with_md_names.pdb` | overwrites without asking |
| Export format | `text` (default) / `csv` | the choices are not displayed |

---

## Feature: 4f. Structure Orientator

**What it is for (plain English):** rotates and moves the whole structure so it lines up with the x/y/z axes (useful for boxes, membranes, pulling or electric-field runs). The coordinates change; the chemistry does not.
**Needs before you start:** a loaded structure.
**Where to find it:** Full Menu `4f`. Workflow Mode "fix" stage.

Sub-menu: `1. Orient structure along Cartesian axes` (● once done). Each run prints "Using <type> PDB file: <name>", "Input structure", "Source: <workspace key>", opens the viewer with reference axes, and then lists:
`1 Align principal axes to Cartesian axes (PCA), 2 Align vector between two residues to an axis, 3 Align vector between two atoms to an axis, 4 Align plane (from residues) normal to Z, 5 Center structure at origin, 6 Apply custom rotation (Euler angles), 7 Cancel`

Every successful run prints "Applying transformation... ✓ Transformed 1,514 atoms ✓ Saved: <input>_oriented.pdb / Transformation saved: <input>_oriented.orientation.json / ✓ Registered in workspace as 'oriented_pdb_file'". The JSON records the method, rotation matrix, translation, parameters and a timestamp. Because the next run takes the previous oriented file as input, the file names grow (`structure_with_md_names_oriented_oriented_..._oriented.pdb`, 6 levels deep in this test). If the output name already exists, the source says it asks "... already exists. Overwrite?" and otherwise adds `_1`, `_2` (SOURCE-ONLY; not triggered).

### Scenario 1: PCA alignment — RUN
- `1` -> "Principal moments of inertia: Longest axis (smallest moment) λ = 21322.2, Middle 42699.6, Shortest 50600.1" -> "Align LONGEST axis to [X/Y/Z] (Z):" -> Enter (Z) -> "Align MIDDLE axis to [X/Y] (X):" -> `Z` -> "Received 'Z', which is not one of the available options" -> `Y` -> "Shortest axis will align to X" -> "Center structure at origin? [y/n] (y):" -> `y`.
- Files: `filter1/hsg/structure_with_md_names_oriented.pdb` + `.orientation.json`.

### Scenario 2: Two residues to an axis — RUN
- `2` -> "First residue (will be at negative end of axis): Chain ID (A): Residue number:" -> A, 25. "Second residue (will be at positive end of axis)" -> B, 25. "Vector length: 6.04 Å / Align this vector to axis [X/Y/Z] (Z):" -> `X` -> centre? `n`.
- Error case: second residue B 500 -> "Residue B:500 not found." and back to the menu.

### Scenario 3: Two atoms to an axis — RUN
- `3` -> A 25, atom "(e.g., CA, FE, ZN) (CA):" `CG`. Second atom B 25, with the atom name defaulting to the first one (CG) -> "Atom-atom distance: 5.17 Å" -> axis `Y` -> centre `y`.
- Error case: atom `XX` -> "Atom B:25:XX not found." and back to the menu.

### Scenario 4: Plane normal to Z — RUN
- `4` -> "Select residues that define a plane (e.g., membrane-facing residues) / format: chain:resnum, comma-separated / Example: A:10, A:20, A:30, A:40".
  - `A:10, A:20` -> "Need at least 3 residues to define a plane." and back to the menu.
  - `A:10, A:20, B:10, B:20` -> "Plane fitted to 4 residues / Normal vector: (-0.987, -0.022, 0.159)" -> centre `y`.

### Scenario 5: Center at origin — RUN
- `5` -> "Use center of mass (vs geometric center)? [y/n] (y):" -> `n` -> "Current center: (-0.00, -0.00, 0.00)" (already centred) and the file is still written.

### Scenario 6: Custom Euler rotation — RUN
- `6` -> "Enter Euler angles in degrees (ZYX convention): Rotation around Z (alpha) (0.0):" `90`, Y (beta) `abc` -> "Please enter a number" -> `0`, X (gamma) `45` -> "Center structure at origin before rotation? [y/n] (y):" `n`.

### Scenario 7: Cancel — RUN
- `7` -> "Cancelled." No file is written.

### Options reference (4f)
| Method | Prompts | Notes |
|---|---|---|
| 1 PCA | longest axis X/Y/Z (Z), middle axis of the remaining two, centre y/n | shortest axis takes the last one |
| 2 Residues | chain/resnum ×2, axis X/Y/Z (Z), centre | first residue at the negative end |
| 3 Atoms | chain/resnum/atom ×2, axis, centre | second atom name defaults to the first |
| 4 Plane | ≥3 `chain:resnum` entries, centre | plane normal goes to Z (for membranes) |
| 5 Centre | centre of mass (y) or geometric (n) | |
| 6 Euler | Z, Y, X angles in degrees, centre before rotating | ZYX convention |
| 7 Cancel | | |

---

## Notes for the manual writer

**Jargon to explain:** asymmetric unit, chain interface / buried surface area (BSA), alternate locations (altlocs) and occupancy, REMARK 465/470, SEQRES, OXT, ACE/NME caps, MODELLER, pKa / ΔpKa / desolvation, protonation probability threshold, Amber residue names (ASH, GLH, HID/HIE/HIP, CYM, CYX, LYN, TYM; constant-pH AS4/GL4), constant-pH MD, isoelectric point, charge capacitance, principal axes, Euler angles.

**Recommended beginner path:** 1a load -> 4a filter (keep protein plus the cofactors and metals you need; usually discard waters) -> 4c analyze with REMARK 465 / REMARK 470 / AltLoc identifiers, then apply caps if wanted -> 4e PROPKA at your pH, accept the recommendations, set MD names -> 4f only if you need a particular orientation, and do it **last**.

**Common mistakes / surprises**
1. The main menu tick `✓` means "available", not "done". Only the `●` inside a tool's own menu shows completion.
2. Each tool picks its input structure by its own hidden priority (see the table at the top). Orientation has top priority in 4f but is ignored by 4e. Re-doing an early step after orienting leaves a stale oriented file that 4f, and possibly later tools, will keep using.
3. 4c without MODELLER: the menu and messages say "capping only". In practice non-standard mutations are applied too, but **only if you add a cap**. Standard mutations, loops, missing side chains and **altloc selection** are not done. The altloc picker exists but cannot be reached without MODELLER.
4. 4c renumbers the protein from 1 when it writes `repaired_structure.pdb` (for example 2RH1 29-1161 becomes 1-506, and ACE shifts 1UBQ by +1). Residue numbers used later in 4e overrides, 4f residue picks or literature comparisons will not match the deposited numbering.
5. Overwrites without asking: `filtered_structure.pdb`, `structure_with_md_names.pdb`. "clear" in the mutator wipes all pending mutations without confirmation.
6. Invalid pH or threshold is silently replaced by the default. Invalid residue-name overrides are skipped without re-asking.
7. Typing a main-menu key (for example `1a`) inside a sub-menu is rejected; use `m` first.
8. Many steps open the web viewer at http://localhost:8765/viewer. Nothing in the terminal waits for it.

**Things that looked like bugs (describe, don't fix)**
- 4a: entering a range (`1-50`) or an out-of-range index at "Select residues" crashes the filter with a traceback, and all choices are lost.
- 4a: the water metrics prompt says "comma-separated" but only accepts one number.
- 4c: template comparison reports OXT missing on every residue (75 for 1UBQ, about 440 for 2RH1).
- 4c: the REMARK 470 detail table shows the chain letter ("A") as the missing atom.
- 4c: FASTA comparison returns 0 for single-chain entries (header "Chain A" vs the parser's "Chains"), and then claims "No structural issues found ... No missing atoms/altlocs" although those checks were skipped.
- 4c: the sequence display misplaces the C-terminal SEQRES segment for 2RH1. "Sequence gaps" counts engineered numbering jumps as 691 missing residues.
- 4c: the summary can show "Pending Non-Standard Mutations 1" and "No pending mutations" at the same time. Skipping caps silently drops a pending non-standard mutation.
- 4c: the no-MODELLER early return skips alternate-location selection, so the output still contains A/B altlocs.
- 4b: lower-case chain IDs are accepted although the help says they are case-sensitive. After a non-standard mutation the sequence display omits that residue.
- 4e: after manual overrides, View/Export still show the PROPKA Protonated/Charge values and net charge, not the overridden states. Export format choices are not listed. The text export gets a `.text` extension.
- 4f: repeated runs chain the `_oriented` suffix (6 levels deep in this test).
- Start-up message "MODELLER not found ... To install: conda install" appears even though MODELLER is installed and only the licence key is missing.
- 1a (other scope, seen in passing): "Remove hydrogen atoms" on 3NIR left 18 H atoms belonging to altloc A of PRO22.
- Help screen workflows (from source `menu_commands.py`) still use old numbering ("3a.Filter → 3c.Fix → 3d.Protonate → 5b.tLEaP"), while the current menu uses 4a/4c/4e/5a. SOURCE-ONLY; the `h` screen was not opened.

**Missing software / licence impact:** MODELLER licence missing, so 4c cannot do loop or residue repair, missing-atom repair, standard mutations or altloc picking. Nothing else in 4a/4b/4e/4f needed unavailable software. PROPKA worked (a few seconds for about 200 residues).

**Tooling note:** session s4-mut-1's ProPrep process vanished during a long chained shell call. This looks like a driver/shell-timeout problem, not a ProPrep crash, because nothing was printed. The same steps were then repeated successfully in s4-mut-2.
