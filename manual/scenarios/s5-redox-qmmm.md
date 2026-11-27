# Redox-active sites and QM/MM preparation (scope s5)

ProPrep 1.23.0, `--full-menu`, driven via `manual/tools/drive.py` on 2026-11-21.
All run folders are under `manual/runs/s5/`. Transcripts are under `manual/runs/_sessions/s5-*/transcript.txt`.

Structures used:

| PDB | What it is | Why it was used | Project folder |
|---|---|---|---|
| 1HRC | horse cytochrome c, one c-type heme with His18/Met80 axial ligands | single-heme path, then all the way to Gaussian/ORCA QM/MM | `runs/s5/hrc/cytc/` |
| 1UBQ | ubiquitin, no cofactors | "no redox groups" path, redox amino acids, search options | `runs/s5/ubq/ubq/` |
| 1OFW | nine-heme cytochrome c, 2 chains, 18 bis-His c-type hemes | multi-heme templates, JSON import, batch microstates | `runs/s5/c3/c3/` |
| 1FLV | flavodoxin, FMN | flavin transformer | `runs/s5/cof/flv/` |
| 1CLF | 2[4Fe-4S] ferredoxin | iron-sulfur transformer, templates | `runs/s5/cof/clf/` |
| 1PLC | plastocyanin, Cu (structure has H atoms) | copper (no transformer exists) | `runs/s5/cof/plc/` |
| 4HHB | hemoglobin, 4 chains, 4 b-type HEM | multi-chain hemes, one big site | `runs/s5/cof/hhb/` |

(3UCP was not tried. 1OFW covered the multiheme case.)

---

## Feature: 2a. Redox Site Detector

**What it is for (plain English):** It finds the parts of a protein that can gain or lose electrons: metal ions, metal cofactors such as hemes and iron-sulfur clusters, organic cofactors such as flavins, disulfide bonds and (if asked) Tyr/Trp/Phe/Met. It then helps you group them into "redox sites", add the protein residues that hold each one (the ligands), and record which atoms are bonded. Later steps (4d Redox Site Preparer, 5a Topology Generator, 5c QM/MM) read these site definitions.

**Needs before you start:** a loaded structure (1a), unless you only want to import a saved sites JSON. Nothing external is needed.

**Where to find it:** Full menu `2a` (section "2. REDOX SITE DETECTION"). Workflow mode lists it under the "detect" stage. There is no command-line flag.

The 2a menu:
```
1. Import redox sites from JSON file
2. Detect redox-active sites            (⚠ "Load a structure first" if none loaded)
3. View detected site summaries         (needs sites)
4. Edit detected sites                  (needs sites)
```
When detection starts, a "How RedoxSite Detection Works" panel lists 6 phases: Configuration (you), Inventory Scan (ProPrep), Selection & Grouping (you), Templates (optional, you), Site Refinement (you), Review & Export (you).

### Scenario 1: Single heme protein (1HRC), default scan, distance search, define bonds, export — RUN
- Starting point: new project `cytc`. 1a → 1 → 1 (RCSB) → 1 (PDB ID) → `1HRC` → 1 (asymmetric unit) → 6 → `m`. Loader warning: "Chain A: First residue is ACE #0" (non-standard numbering).
- Steps:
  1. `2a`, then `2` (Detect). Message: "Using RCSB PDB Structure object (workspace key: rcsb_structure)". The 6-phase overview panel is shown.
  2. Configuration (six questions, Enter = default):
     - "Metal elements to scan (comma-separated list, 'none', or 'all') (all)" → Enter. A table of recognised metal groups is shown first.
     - "Coordinating non-metal elements ... (all)" → Enter.
     - "Search radius around each metal (Å) (range: 0.1-20.0, default: 4.0)" → Enter.
     - "Detect CYS SG–SG pairs within this distance as candidate disulfides (Å) (range: 1.0-5.0, default: 2.5)" → Enter.
     - "Non-standard residues (HEM, FAD, etc.)? (y)" → Enter.
     - "Redox-active amino acids (TYR, TRP, PHE, MET)? (n)" → Enter.
     - A "Final Configuration" table follows.
  3. Inventory: "Found 2 potential redox centers": HEC A:105 (organometallic cofactor; nearby MET 80 SD 2.3 Å, HIS 18 NE2 2.5 Å …) and **ACE A:0 (organic cofactor)**. ACE is the acetyl N-terminal cap, not a cofactor (a false positive).
  4. Disulfides: "⚠ No SG–SG pairs were within the 2.50 Å auto-detection cutoff." A CYS list and a pairwise SG–SG distance table are shown, then "Specify disulfide bonds as index pairs (e.g. '1-2, 3-4') or 'none' (none):" → Enter.
  5. Center selection: the table renumbers the centers (ACE becomes 1, HEC becomes 2). Options: `all`, `none`, `view`, ranges, lists. "Select additional centers to include (all):" → `2` (HEC only).
  6. Grouping: "Enter site assignments (format: site:centers site:centers)", Enter = auto (one site per center) → Enter. "Auto-assigning individual sites: 1:1".
  7. No template question, because there is only one site.
  8. Refinement for site_1: "1. Distance cutoff / 2. Residue count cutoff / 3. Finish refinement" → `1`. Boundary definition "1 Geometric center / 2 Minimum distance from any site atoms / 3 custom atoms" → `2`. Distance method "1 Fixed / 2 Adaptive" → `1`. Radius `2.5`. Atom filtering "Use same atom filtering for search? (y)" → Enter.
  9. Search results: MET A:80, HIS A:18, CYS A:14, CYS A:17. "Select residues to add:" → `all`. "Added 30 atoms from 4 complete residues to site."
  10. BOND DEFINITION: "1. Define bonds between residue pairs / 2. Skip / v. View" → `1`. Pairs → `105-80 105-18 105-14 105-17`.
  11. For each source atom, a table of target atoms sorted by distance appears. Typed `FE` → `1,11` (NE2, SD). Then `CAB` → SG of Cys14, and `CAC` → SG of Cys17 (see Scenario 2). Then `done`.
  12. A site summary lists 4 bonds: Fe–NE2 (coordinate 2.04 Å), Fe–SD (coordinate 2.32 Å), CAB–SG (covalent 1.75 Å), CAC–SG (covalent 1.86 Å).
  13. "Search for additional residues to add to this site? (y)" → `n`. "Final site contains 1 center(s), 73 atoms, 4 bonds."
  14. Site Review Options: 1 View summary, 2 Edit, 3 Manually re-run, 4 Export, 5 Continue. `1` → `all` shows the table. `4` → "Choose export format (comma-separated list)" → `1,2,3`. `5` → "✓ Successfully detected 1 redox-active site(s)".
- Result / files (in `runs/s5/hrc/cytc/`):
  - `1HRC_redox_sites.json`: full site data, re-importable (2a option 1).
  - `1HRC_redox_sites.pdb`: all site atoms in one PDB.
  - `1HRC_site_1.pdb`: one PDB per site.
- Problems:
  - ACE flagged as a cofactor.
  - The search results list the same atom twice ("SD S (2.28Å), SD S (2.32Å)" counted as "2 atoms within range"), because each site atom that is near it counts it again.
  - The menu suggestion reads "...or press  to return to the main menu". The `[m]` disappears from the text on every menu (a markup bug).
  - The single-site path never shows the transformer's bond requirements. Those are only shown when a template/site type is set (multi-site runs). See Scenario 4 and the 4d section: the bonds made here were **not** enough for the 4d transformer.
- Session log: `runs/_sessions/s5-hrc/transcript.txt`; session file `runs/s5/hrc/cytc/proprep_session_20261121_122846.json`.

### Scenario 2: Wrong bond target picked, then corrected — RUN
- Starting point: Scenario 1 step 11.
- Steps: source `CAB`, then typed `28` thinking it was Cys14 SG. The target list is re-sorted by distance for every source atom, so 28 was **CYS A:17 SG at 8.21 Å**. ProPrep said "Created bond: CAB C ↔ SG S (8.21Å)" with **no warning** about the unrealistic length.
- Fix: typed `CAB` again → "Source atom CAB already has 1 bond(s)" with options "1 Delete existing bonds and create new one / 2 Replace one existing bond / 3 Add another bond (keep existing) / 4 Cancel" → `2`. Then target `SG` → "Multiple atoms named 'SG' found: 1. [1] SG in CYS A:14 (1.75Å) 2. [28] SG in CYS A:17 (8.21Å)" → `1`.
- Lesson for users: type the atom **name** (e.g. `SG`), not a number, or check the distance column. ProPrep does not reject an 8 Å "bond".

### Scenario 3: Review/edit menu (validate, list bonds, metadata, cancel) — RUN
- From Site Review Options, `2` opens the "RedoxSite Editor – Site Selection" (also reachable later as 2a option 4). Its commands are `edit <N>`, `remove <N>` (also `1-5`, `all`), `validate`, `save`, `cancel`.
- `validate` → "✓ Site validation passed - no issues found".
- `edit 1` → operations `1 centers / 2 residues / 3 atoms / 4 bonds / 5 metadata`, plus `summary`, `validate`, `back`.
- `4` (bonds) → `list` shows separate "Coordinate Bonds" and "Covalent Bonds" tables. The other bond commands are `add`, `remove`, `reclassify`, `back`.
- `5` (metadata) shows Site Type N/A, Detection Method, Boundary Definition. "Press Enter to keep current value".
- `cancel` → "⚠ You have unsaved changes", "Discard all changes and exit? (n)" → `y` → "No changes made". (It reported unsaved changes even though only Enter had been pressed.)

### Scenario 4: Adding the extra bonds a transformer needs, via the editor — RUN
- Why: 4d said no specific transformer fitted 1HRC (see 4d Scenario 1). `why 1` showed that the His/Met c-type heme transformer needs 6 bonds that are **inside** residues in the crystal structure: HEC C2A–CAA, HEC C3D–CAD, and CA–CB of each Cys, the His and the Met.
- Steps: `2a` → `4` → `edit 1` → `4` (bonds) → `add`, six times. Each `add` asks: Atom 1 Chain, Resname, ResID, Atom name; Atom 2 the same four; "Distance between atoms: 1.56Å"; bond type 1 covalent / 2 coordinate / 3 disulfide / 4 metal-metal → `1`; "Add covalent bond (1.56Å)? (y)" → `y`. Example: `A HEC 105 C2A` / `A HEC 105 CAA`.
- `list` showed 2 coordinate and 8 covalent bonds (the 6 new ones marked "intraresidue"). Then `back`, `back`, `save` → "✓ Changes saved to workspace".
- Result: 4d then listed `heme_his_met_axial_c_type` as compatible.

### Scenario 5: Editor ResID prompt given a non-number — FAILED (bug)
- Steps: in the bond `add` dialog, the Chain prompt accepted the whole line `A HEC 105 C3D A HEC 105 CAD` as a chain name (no check). Then `1` went into Resname and `y` into ResID.
- What the user sees: a log line "ERROR ... Error handling menu option 'edit': invalid literal for int() with base 10: 'y'", and the program drops straight back to the 2a menu. **All unsaved editor changes are lost.**
- The same `int(...)` pattern is used for every ResID prompt in the editor (add/remove/modify atoms, residues, bonds; source code).

### Scenario 6: Protein with no redox groups (1UBQ), defaults — RUN
- 1a load 1UBQ, then `2a` → `2`, all six defaults.
- What the user sees: "Found 0 potential redox centers", "No redox centers found in structure.", "No redox sites detected". Back to the 2a menu, with options 3/4 still blocked.
- Session: `runs/s5/ubq/ubq/proprep_session_20261121_130802.json`.

### Scenario 7: 1UBQ with redox amino acids, invalid answers, grouping errors — RUN
- Inputs to configuration:
  - Metals `XX`: accepted silently, with no check that it is an element.
  - Non-metals: `back` returns to the metals question ("(Type 'back' or 'b' to return to previous question)" is shown on questions 2–6). Then metals `none`, non-metals `N,O,S`.
  - Radius `25` → "Value must be at most 20.0", then `4.0`. Disulfide `y` → "Please enter a valid number", then `2.5`.
  - Non-standard `y`, redox amino acids `y`.
- Inventory: MET A:1 SD, PHE A:4 CD1, PHE A:45 CD1, TYR A:59 CZ. "No CYS/CYX residues found in structure."
- Center selection `1-99` → "Invalid selection. Please try again." Then `all`.
- Grouping errors:
  - `1:1 2:2,3 3:5` → "Invalid assignment format or values."
  - `1:1 2:2,3` → "Validation error: Centers not assigned to any site: 4".
  - `1:1 2:2,3 3:4` → accepted. Site 2 is a "multi-center site" (two Phe).
- Because there are now 3 sites, the **Template-Based Automation** explainer appears: templates are your own tags for detection recipes; transformers are pre-coded force-field renaming; if a template name matches a transformer name it is auto-applied in 4d. "Use template-based automation for systematic redox sites? (n)" → `n`.
- Site 1 refinement: `3` Finish → "Need at least 2 residues to define inter-residue bonds." (Bond definition is skipped automatically for one-residue sites.)

### Scenario 8: Residue-count search, custom boundary atoms, search filters, "finish" — RUN (1UBQ site 2)
- `2` (Residue count cutoff).
  - "Enter residue types to find (comma-separated, e.g., HIS,CYS,ASP) (HIS,CYS)" → `LYS,GLU`.
  - "Enter number of each type to find (comma-separated, e.g., 2,3 for 2 types)" → `1,1`.
  - Boundary `3` (custom atoms): an atom table appears. Typing `CZ` gives "Multiple atoms named 'CZ' found" → `1,2`.
  - "Use same atom filtering for search?" → `n`. Then "Include metals" → `none`, "Include non-metals" → `all`.
- Result: "Find 1 LYS, 1 GLU residues" → LYS A:48, GLU A:64. At the residue prompt, `finish` = "Finish refinement now (add nothing more, keep the site as-is)". The bond prompt follows; `2` skipped it.
- Note: the current filter is displayed as Python lists, e.g. `Include metals: []` (meaning *all*) or `['NONE']`. This is cryptic for non-programmers.

### Scenario 9: Adaptive radius, "0 = adjust", adding a water, manual re-run, export error — RUN (1UBQ site 3)
- Distance → geometric center → `2` Adaptive. Min radius `3`, max `8`, step `0.5`, target coordination `4` → "Adaptive distance from geometric center, Search radius used: 4.50Å, Total atoms found: 2".
  - The target of 4 seems to have been met by atoms of the Tyr itself, which are then hidden from the results, so only 2 atoms are shown. This is confusing.
- At "Select residues to add" typed `0` → "Adjust search parameters and repeat analysis" (back to boundary). Then min-distance, fixed `3.5` → 8 residues including HOH A:100. Typed `6 3` (space-separated list) → added HOH and ASP 58. Bonds `2` (skip).
- Site Review `3` (Manually re-run detection on specific sites) → `1` → "This will replace 1 site(s) with manually refined versions. Continue? (y)" → refinement restarts from the bare center → `3`.
- Export with `1,5` → "Invalid choices: 5. Valid options are 1, 2, 3" (nothing exported). `1` then worked → `1UBQ_redox_sites.json`.
- 2a option `3` (View detected site summaries) prints a "Site Summary" box per site, then "Press Enter to continue...".

### Scenario 10: Multiheme cytochrome (1OFW) with templates — RUN
- Inventory: 18 HEC (9 per chain, each with two His) plus 4 GOL and 1 ACT flagged as "organic cofactor". Glycerol and acetate are crystallisation additives, i.e. false positives. The SG–SG table lists 586 more pairs as hidden.
- First attempt with hemes `1-4` was lost when the terminal session died. That was a driver/environment problem, not a ProPrep message. It was redone as session `s5-c3b` with hemes `1-2`.
- Steps:
  - Grouping Enter (auto). Templates `y`. A table of all transformers by family appears (Hemes / Iron-sulfur / Organic cofactors / Metal sites / Utility).
  - "Enter site type assignments" → `1-2:"heme_bis_his_c_type"`. Output: "Transformer Mappings: heme_bis_his_c_type → heme_bis_his_c_type".
  - The template is configured on site_1, with a yellow "Caution: Distance cutoffs may not be satisfied at every site... Residue count is generally more reliable for templates." The default is now `2` (count). Typed `2`, `HIS,CYS`, `2,2`, boundary `2`, Enter → HIS 37/51, CYS 47/50 → `all`.
  - Bond step: **"Bond Requirements for heme_bis_his_c_type"** is printed, listing exactly the 6 intra-residue bonds needed (HEC C2A↔CAA, C3D↔CAD, CYS CA↔CB ×2, HIS CA↔CB ×2).
  - Same-residue pairs are accepted: `1294-1294 37-37 51-51 47-47 50-50`. Then `C2A`→`CAA`, `C3D`→`CAD`, `done`, and `CA`→`CB`,`done` for each of the four protein residues.
  - "Apply these settings to all other heme_bis_his_c_type sites (1 remaining)? (y)" → `y`. "Applied template bond: …" ×6 for HEC A:1295. A results table follows: Site 1 ✓ 75 atoms 6 bonds; Site 2 ✓ 75 atoms 6 bonds.
  - Export `1` → `1OFW_redox_sites.json`.
- Sessions: `runs/s5/c3/c3/proprep_session_20261121_131538.json` (lost run), `..._132318.json` (template run).

### Scenario 11: Import sites from JSON — RUN
- New session, local file `1OFW.pdb` loaded via 1a → 1 → 4. Then `2a` → `1`: a file browser lists the JSON files in the project folder, **including ProPrep's own session files**.
- Picking a session file (`proprep_session_….json`) gives "Successfully imported 0 redox site(s) from …" immediately followed by "No sites found in JSON file". These two messages contradict each other.
- Picking `1OFW_redox_sites.json` gives "✓ Successfully imported 2 redox site(s)" and "Restored transformer mappings for 1 site type(s)". 4d then auto-assigned the transformer (4d Scenario 4).

### Scenario 12: Flavin (1FLV), iron-sulfur (1CLF), copper (1PLC), hemoglobin (4HHB) — RUN
- **1FLV:** 1 organic cofactor FMN A:170. "Only one CYS/CYX residue is present -- a disulfide bond needs two, so there are no pairs to specify." (The disulfide question is skipped.) Finished without adding residues. 4d offered `flavin_fmn`.
- **1CLF:** 2 SF4 clusters, each listed as one "Organometallic Cofactor" with 4 Cys SG nearby (3.7–3.9 Å from the cluster reference point).
  - Templates `y` with `1-2:"iron_sulfur_4fe4s"`. Count search `CYS` `4`. Bond requirements shown: "4 coordinate Fe-SG bonds … ANY:FE1 ↔ CYS:SG … FE4".
  - Defined with pairs `56-47 56-11 56-14 56-8`, then `FE1`→`1`, … (the nearest atom was always the right SG). The template was auto-applied to SF4 A:57.
- **1PLC:** the file contains 704 H atoms. During loading: "Remove hydrogen atoms? (y)" (writes `1PLC_noH.pdb`). 2a then asks "Select structure to use (1)" between H-Stripped and RCSB.
  - Inventory: "Metal Ions: CU A:100" with His37 ND1 1.9, His87 ND1 2.1, Cys84 SG 2.1, Met92 SD 2.8 Å.
  - Fixed 3.0 Å min-distance search found all four. 4d: no specific transformer (copper needs parameters from elsewhere, e.g. 2c Force Field Parameterizer; outside this scope).
- **4HHB:** 4 HEM, each with only the proximal His NE2 at 2.6–3.1 Å. Two Cys pairs at 9.3 Å. Grouped all four as one site with `1:1-4` ("organometallic site (HEM)").
  - 4d: no transformer fits. `why 1` → "Required center count (expected 1, got 4)". The 4d table shows only "Chain A, HEM 142" as the location of the 4-center site.

### Scenario 13: No structure loaded — RUN
- New project, `2a`: option 2 shows "⚠ Load a structure first". Choosing `2` anyway prints "⚠ Load a structure first". Option 1 (import JSON) is always available.
- `4d` can be entered but shows "⚠ Needs detected redox sites — run Redox Site Detector first". 5c is marked ○ in the main menu with the same warning.

### Scenario 14: 3D viewer hooks (`view`, `v`) — SOURCE-ONLY
- Not run (it would open a browser viewer). Per the source, `view` in center selection and residue selection and `v` in bond definition open the NGL web viewer. Sites, inventories and boundaries are also highlighted automatically when a viewer is open.

### Scenario 15: Editor operations not exercised (centers/residues/atoms add/remove/replace, `remove <N>`, `reclassify`) — SOURCE-ONLY
- From `structure_prep/redox_site_editor.py`:
  - centers: list/remove/modify (chain, resname, resid, type);
  - residues: add (chain, resid, insertion code), remove, replace;
  - atoms: list with filter (chain/resid/element), remove, modify (chain, resname, resid, name, element);
  - bonds: add/remove/reclassify (e.g. `coordinate 1`);
  - metadata: site type.

### Options reference — 2a
| Prompt | Choices | Meaning |
|---|---|---|
| Metal elements to scan | `all` (default), `none`, list e.g. `FE,CU` | which metal elements become centers (no validation of typed symbols) |
| Coordinating non-metal elements | `all`, `none`, list | which atoms count as ligands around metals |
| Search radius around each metal | 0.1–20 Å (4.0) | radius for the "nearby atoms" in the inventory |
| Disulfide cutoff | 1.0–5.0 Å (2.5) | SG–SG distance for automatic disulfides |
| Non-standard residues? | y (default) / n / back | include HEM, HEC, FMN, SF4… (also catches ACE, GOL, ACT) |
| Redox-active amino acids? | n (default) / y / back | include TYR, TRP, PHE, MET |
| Specify disulfide bonds | `none` or index pairs `1-2, 3-4` | manual disulfides; only shown when ≥2 Cys. If some were auto-found: "Specify additional disulfides? (n)" |
| Select additional centers | `all`, `none` (disulfides only), `view`, `1-10`, `1,3,5-8` | which centers go forward |
| Site assignments | Enter = one site each; `1:1 2:2`, `1:1,2`, `1:1-4`, mixed | group centers into sites; every center must be assigned |
| Use template-based automation? | n (default) / y | only when >1 site. Then `1-18:"type"`, `19,21:"type", 20:"other"`. A name equal to a transformer name auto-maps it |
| Choose search constraint | 1 Distance, 2 Residue count, 3 Finish | default 1 (2 in template mode) |
| Boundary definition | 1 geometric centre, 2 min distance to any site atom, 3 custom atoms | where distances are measured from |
| Distance method | 1 Fixed (radius, default 5.0 Å), 2 Adaptive (min 3.0, max 8.0, step 0.5, target coordination 6) | |
| Count search | residue types (HIS,CYS), counts (2) | N closest residues of each type |
| Use same atom filtering for search? | y / n → metals, non-metals | filters for this search only |
| Select residues to add | numbers/ranges/`all`, `none`/`0` = adjust and repeat, `finish`, `view` | whole residues are added |
| Define bonds | 1 define pairs, 2 skip, v view | pairs by residue ID `202-80`, `MN202-GLU80`, `A:202-A:203`, row numbers, same-residue `105-105` |
| Source / target atom | number, atom name, `done` / `cancel`; ambiguous names → pick list | target numbers are **re-sorted by distance** for each source |
| Source atom already bonded | 1 delete & new, 2 replace one, 3 add another, 4 cancel (default) | |
| Search for additional residues? | y (default) / n | another search iteration |
| Apply these settings to all other <type> sites? | y / n | template mode only |
| Site Review Options | 1 view, 2 edit, 3 re-run selected sites, 4 export, 5 continue | |
| Export format | `1` JSON, `2` one PDB, `3` one PDB per site, or a list `1,2,3` | files `<PDB>_redox_sites.json/.pdb`, `<PDB>_site_N.pdb` |

---

## Feature: 4d. Redox Site Preparer

**What it is for (plain English):** AMBER force fields only know standard amino acids. A heme, flavin or iron-sulfur cluster, and the protein residues bound to it, need special residue and atom names that match a parameter library. This tool "transforms" each redox site:
- it renames residues and atoms;
- it splits residues (e.g. moves the side chains of the heme's Cys/His/Met into the heme unit, and makes the propionates separate residues);
- it chooses the library for the oxidation state you pick (oxidised/reduced, etc.).

The result is a PDB that tLEaP (5a) can build.

**Needs before you start:** detected or imported redox sites (2a). For a specific transformer to be offered, the site must contain the right residues **and** the right bonds (see Scenarios 1–2).

**Where to find it:** Full menu `4d`. Workflow mode: "fix" stage.

Menu:
```
1. Create transformer (interactive PDB editor)
2. Generate single redox microstate
3. Generate all redox microstates (batch processing)
```
After any generation, all three options are blocked ("⚠ Already generated microstate(s)"). To redo, start a new session/project.

Overview panel ("How RedoxSite Transformation Works"): 1 Transformer Assignment (you); 2 Residue ID Mapping, 3 Component Matching, 4 Apply Transformations, 5 Update Workspace (all done by ProPrep).

Built-in transformers (from the 2a template list and the source):
- **Hemes:** `heme_bis_his_c_type`, `heme_bis_his_b_type`, `heme_cys_axial_b_type`, `heme_his_met_axial_c_type`
- **Iron-sulfur:** `iron_sulfur_4fe4s`
- **Organic cofactors:** `flavin_fmn`, `flavin_fad`, `nicotinamide_nadp`, `nicotinamide_nad`, `pterin_biopterin`
- **Metal sites:** `zinc_cys4`
- **Utility:** `no_transformation`, `disulfide`

User-made transformers are stored in `~/.proprep/transformers`.

### Scenario 1: 1HRC straight after 2a Scenario 1 — no specific transformer; `why`; cancel — RUN
- `4d` → `2`.
  - "✅ Captured original states: 1 sites, 74 atoms". "Using RCSB PDB PDB file for transformation: …/1HRC.pdb" (printed twice, with a doubled "PDB PDB").
  - Table: Site 1, HEC, "2×CYS, 1×HIS, 1×MET", Compatible Transformers **None**.
  - "1 site(s) have no specific compatible transformer; detailed diagnostics written to transformer_compatibility_report.txt."
- "Enter transformer assignments:" → `why 1`. This lists every transformer and the requirements it failed. For `heme_his_met_axial_c_type`: "12/18 requirements met … HEM-variant bonds (inter-residue post-transformation): 2 thioether-Cys + axial-His + axial-Met side-chain splices + 2 propionate carve-outs (expected 6, got 0)".
- Enter on an empty assignment → "No assignments provided" → back to the 4d menu (this is how you cancel).
- File: `transformer_compatibility_report.txt` (the same diagnostics as `why`).

### Scenario 2: 1HRC single microstate with `heme_his_met_axial_c_type` — RUN
- Starting point: the 6 extra bonds had been added (2a Scenario 4). `4d` → `2`: Compatible → `heme_his_met_axial_c_type`. List: "1. heme_his_met_axial_c_type 2. no_transformation". Typed `1` (with one site the bare index works).
- Step 2, Configure Parameters. For each parameter you may type a value, batch syntax `1:value`, or Enter for per-option prompts ("Sites for 'oxidized' (numbers/ranges, empty=skip)").
  - "Iron oxidation state (oxidized/reduced)" → `oxidized`.
  - "Electronic spin state … Setting 'low_spin' for all sites" (fixed, no question). Note shown: "Ferric is S=1/2 (low-spin doublet), ferrous is S=0".
  - "Propionate/titratable-site pH treatment (fixed_pH/constant_pH)": fixed_pH = static PRP/PRD propionates with ff14SB/ff19SB; constant_pH = titratable PRN from AMBER constant-(pH,E) libraries (needs leaprc.constph + leaprc.conste, ff10 backbone). Enter → per-option → `1` for fixed_pH.
  - "Protonation state of A-ring propionate (deprotonated/protonated)" → `1:deprotonated`. D-ring → `1:deprotonated`.
- Configuration Summary, then "(a)pply | (e)dit step | (r)estart | (c)ancel". The prompt shows as "Choice  (a)" because the `[a/e/r/c]` part is eaten by markup.
- Tried `e` → "Which step to edit? (1) Assign Transformers (2) Configure Parameters" → `2`:
  - Typing `ferric` → "Skipping batch entry 'ferric': invalid site range (expected N or N-M)". It then falls back to per-option prompts.
  - `constant_pH` → "No options to set for parameter 'protonation_propionate_a' in this configuration. Skipping."
  - Edited back to oxidized/fixed_pH/deprotonated ×2, then `a`.
- "Residue ID Mapping: No ID conflicts detected". Component Matching Summary table: Bring Cys A 14, Center A 105, Cring Cys A 17, His A 18, Met A 80, Prop A A 106, Prop D A 107.
- "Review specific sites? (no)" → `1` shows the matched components → "Is this component matching correct? (y)" → `y`.
- The Transformation Plan has 12 operations, e.g.:
  - "Move ring-B thioether Cys side chain (CB/SG/HB2/HB3) to heme as CBC1/SGC1…"
  - "Reduce axial His to HIO backbone stub"
  - "Create PRD residue from propionate A group"
  - "Rename HEC to HMO"
  It is applied with progress lines. The final box shows HEC split into HMO 105 + PRD 106 + PRD 107, and MET/HIS/CYS reduced to MEO/HIO/CYO backbone stubs.
- Validation: "✅ Workspace validation passed • 1 centers validated • 73 atoms validated • 10 bonds validated".
- Files (`runs/s5/hrc/cytc/`):
  - `final_transformed_structure.pdb`: the transformed structure. Residues are **renumbered sequentially** (ACE 0→1, so HMO ends up as 106, CYO 15/18, HIO 19, MEO 81).
  - `enhanced_residue_mapping.txt`: old → new residue/atom names. It uses the *pre-renumbering* IDs (105, 106, 107), so it does not match the final PDB numbering.
  - `final_transformed_structure_redox_sites_updated.json`: the sites with IDs matching the transformed PDB.
  - Per-site plan/log JSON and TXT are written to a temporary folder, e.g. `/var/folders/.../redox_transformation_4k_9omx7/transformation_reports/`, which is not in the project folder.
- Problem found downstream: the propionate PRD residues were placed between crystal waters (PDB order HMO 106, HOH 107, PRD 108, HOH 109, PRD 110), even though "No ID conflicts detected" was printed. tLEaP coped, but 5c then mis-mapped residues (5c Scenario 3).
- Prerequisite for 5c (done here, belongs to scope 5a): `5a` → `2`.
  - Choices: ff14SB, no modified AA, no DNA/RNA/carbohydrate/lipid/GAFF, TIP3P (listed as "TIP3P (Recommended — matches ff14SB) 3-point, traditional (not recommended)", which contradicts itself), no divalent ions.
  - Forcefield set: "1. Ferric His/Met-axial c-type heme (HMO; cyt c2) — CM5". Options 2–4 are Fe-localized Δq, pooled RESP and single-conformer RESP.
  - Implicit solvent `2`, don't edit the template, default file name.
  - Then `4` → run → `final_transformed_structure.prmtop/.rst7`, charge +7, ParmEd validation passed.
  - tLEaP was run from `/opt/homebrew/Caskroom/miniforge/base/envs/etanalyze_env/bin/teLeap`, i.e. whatever `tleap` is first on the PATH, not necessarily the ProPrep env.

### Scenario 3: 1UBQ, three amino-acid sites, batch with `no_transformation` — RUN
- `4d` → `3` "Batch Redox Microstate Generation".
  - Compatibility table (all `no_transformation (1.00)`). "Enter assignment (1-3:"no_transformation")" → Enter, then preview, "Accept these assignments? (y)" → `y`.
- Step 3, state selection:
  - The table lists "1: unchanged/unchanged ✓★" per site.
  - The help reads "Selection Syntax: :,,... :,,..." because the `<site>:<state>` placeholders are eaten by markup. Examples: `1:1,2 2:1,2 3:1`, `all:1,2`, `single:1`, `single:2[1-23]`, `pairs:1`, `pairs:1[5-10]`, `single:1 pairs:1`.
  - Typed `all:1` → 1 microstate coded **"OLOLOL"**. The legend says R/O = reduced/oxidized and L/H/I = spin, so "unchanged" is shown as "Oxidized/Low", which is misleading.
- Step 5: "1. Generate all of them 2. Pick a subset" → `1` → `transformed_microstate_001.pdb`, `..._redox_sites_updated.json`, `microstate_metadata.json`. Message: "Use this file for tLEaP input generation (option 3 in Topology Generator menu)".

### Scenario 4: 1OFW, two bis-His c-type hemes, batch with subset and per-state pH treatment — RUN
- After the JSON import (2a Scenario 11), `4d` → `3`:
  - "Found transformer mappings from detection phase. Attempting automatic assignment... Site 1 (heme_bis_his_c_type): heme_bis_his_c_type ✓".
  - "Heme redox treatment: fixed_E (each microstate builds its redox state into the topology; the constant-E HEH treatment is a single-state option)."
- States per site: "1: reduced/low_spin, 2: oxidized/low_spin". `all:1,2` → 4 microstates RLRL, RLOL, OLRL, OLOL. Step 5 `2` (subset) → `1,4`.
- Step 6, Propionate/pH treatment per redox state:
  - Note shown: "leave a prompt empty to take the default (constant_pH / titratable PRN)".
  - Oxidized: `1-2:fixed_pH`, A-ring `1-2:deprotonated`, D-ring `1:deprotonated,2:protonated` (per-site values work).
  - Reduced: Enter everywhere → constant_pH (defaults).
- Output: "Generating 1/2: RLRL … ✓ A_HEC_1299 … ✓ A_HEC_1302". Files: `transformed_microstate_001.pdb`, `_002.pdb`, each `_redox_sites_updated.json`, `microstate_metadata.json`.
- Confusing: the progress lines name the hemes `A_HEC_1299` and `A_HEC_1302`, but the selected hemes were 1294 and 1295. Coordinates confirm the right hemes were transformed (Fe of 1294 → `HCO A 590`, 1295 → `HCO A 596`). All hetero groups are renumbered sequentially (HEC 587–597, PRD/PRP placed between them), while protein residue numbers stay unchanged.
- Session: `runs/s5/c3/c3/proprep_session_20261121_132916.json`.

### Scenario 5: FMN (1FLV) — RUN
- `4d` → `2` → `flavin_fmn` (`1`).
- "Redox/protonation state of the cofactor": oxidized, semiquinone_neutral (FMNH•), semiquinone_anionic (FMN•⁻), hydroquinone_anionic (FMNH⁻), hydroquinone_neutral (FMNH₂). Typed `semiquinone_neutral`. Spin is "set by the redox state; not user-configurable".
- Plan: "Rename FMN → FMS". Files `final_transformed_structure.pdb`, `enhanced_residue_mapping.txt`, `…_redox_sites_updated.json` in `runs/s5/cof/flv/`.

### Scenario 6: Two [4Fe-4S] clusters (1CLF) with broken-symmetry choice — RUN
- Auto-assigned from templates ("✓ Successfully auto-assigned transformers for all 2 sites").
- "Cluster oxidation state":
  - oxidized = [Fe4S4]^3+ (HiPIP oxidized, S=1/2)
  - reduced_hipip = [Fe4S4]^2+ (broken-symmetry singlet)
  - reduced_fd = [Fe4S4]^1+ (S=1/2)
  Typed `1:reduced_hipip,2:reduced_fd`.
- "Broken-symmetry parameter set. Valid codes per redox state: oxidized=FO1-FO6, reduced_hipip=FR1-FR3, reduced_fd=FD1-FD6". Also listed: FR0/FD0 = "Averaged parameters across all BS guesses".
  - Typed `1:FD1,2:FD0` → "'FD1' is not a valid option for site 1 (valid: ['FR1', 'FR2', 'FR3', 'FR0']); skipping." Then the per-option prompt gave site 1 `FR1`.
- Plan: SF4 56 → FR1, CYS 8/11/14/47 → C25/C26/C27/C28 (per-Fe ligand names). SF4 57 → FD0, its Cys → CXD.

### Scenario 7: No specific transformer (1PLC copper, 4HHB four hemes as one site) — RUN
- Only `no_transformation` is offered, and the diagnostics file is written. For 4HHB, `why 1` shows every heme transformer needs exactly 1 center per site. For a Cu site there is no transformer; parameters would have to come from 2c (MCPB.py, outside this scope).

### Scenario 8: Create transformer (interactive PDB editor) — RUN (editing), SOURCE-ONLY (save)
- `4d` → `1`. With several sites: "2 detected redox site(s): 1. site_1 (75 atoms) 2. site_2" → `1`.
- A table of the site's residues and atom names appears (chain, resid, resname, record, atoms). Prompt "edit".
- Commands run, with results:
  - `rename_res A 1294 HXO` → ok
  - `movenew A 1294 CAA,CBA,CGA,O1A,O2A PRX 1` → new residue A 1295 PRX
  - `rename_atom A 37 CB CBX` → ok
  - `bogus` → "! Unknown command 'bogus'. Type 'help'."
  - `summary` → numbered operation list
  - `drop` → "Reverted last operation."
  - `lib` → "No library is linked yet …; the link is chosen at 'save'."
  - `param redox_state oxidized,reduced` → "Declared parameter … (reference: oxidized)"
  - `vary 1 HXO,HXR HXO` → "! Expected 1 state value(s) ." (Correct form: `refstate redox_state=oxidized`, then `vary 1 reduced HXR`.)
  - `quit` → "Discard this transformer? (n)" → `y` → "Discarded."
- The `summary` after `vary` printed "state  -> HXO / state  -> HXR" with the state label missing (markup).
- `save` was not run, because it writes to the shared `~/.proprep/transformers/`. Per source, it asks for a template name (warns if one exists: "Replace it?") and lets you pick a deposited force-field library (numbered list or path) and the redox/spin leaf.
- Full command list (from `help`): `rename_res`, `rename_atom`, `id`, `chain`, `icode`, `hetatm`, `atom`, `movenew`, `move`, `param`, `refstate`, `vary`, `lib`, `drop`, `show`, `summary`, `save`, `quit`. Names are case-sensitive as typed.
- Note: HEC was shown with record "ATOM" in the table although it is HETATM in the PDB.

### Scenario 9: Other branches — SOURCE-ONLY
- Verbosity choice appears only for ≥20 sites: (f)ull / (c)ompact (default) / (m)inimal.
- `r` restart and `c` cancel at the Apply prompt.
- "Is this component matching correct?" → `n` → "1 Skip this site / 2 Cancel all transformations".
- In batch mode: `single:`/`pairs:` systematic selections. Step 5 asks for confirmation when the microstate count is large.

### Options reference — 4d
| Prompt | Choices | Meaning |
|---|---|---|
| Enter transformer assignments (single) | `1:"name"`, bare index or name (1 site), `why <site#>`, Enter = cancel | which transformer for each site |
| Enter assignment (batch) | `1-3:"name", 4-6:"other"` (suggested default shown) | same, for batch |
| Parameter (each) | a value; batch `1:value, 2-3:value`; Enter = per-option prompts ("Sites for 'X'") | sites left unassigned get the default |
| Heme redox_state | oxidized / reduced | e.g. HMO/HMR (His/Met), HCO for oxidized bis-His (only HCO was seen in the output files) |
| spin_state | fixed per transformer (heme low_spin); Fe4S4 BS codes FO1-6/FR1-3/FD1-6 or FR0/FD0/(FO0) averaged | |
| ph_treatment | fixed_pH (PRP/PRD, ff14SB/ff19SB) / constant_pH (PRN, ff10 + constph/conste) | single path: no default was shown; batch path default is constant_pH |
| protonation_propionate_a / _d | deprotonated / protonated | fixed_pH only |
| FMN/FAD redox_state | oxidized, semiquinone_neutral/anionic, hydroquinone_anionic/neutral | renames FMN → e.g. FMS |
| Apply prompt | a apply, e edit step (1 assign / 2 parameters), r restart, c cancel | |
| Review specific sites? | no (default), `all`, `1,2`, `site_1` | show matched atoms and confirm |
| Batch state selection | `1:1,2 2:1`, `all:1,2`, `single:N[X-Y]`, `pairs:N[X-Y]` | builds microstate combinations |
| Step 5 | 1 all, 2 subset (`1,3,5`, `1-4,7`) | which microstates to write |

---

## Feature: 5c. QM/MM Preparator

**What it is for (plain English):** It writes input files for a combined quantum/classical calculation. The redox site (e.g. the heme and its ligands) is treated with quantum chemistry (QM), and the rest of the protein with the AMBER force field (MM). Two programs are supported:
- **Gaussian**, using the ONIOM method: ProPrep writes a `.com` file.
- **ORCA**, using its QMMM module: ProPrep writes a `.inp` file, a PDB and a run script.

ProPrep only writes the files. Running them needs Gaussian 16 or ORCA, which you must install and license yourself.

**Needs before you start:**
- redox sites (2a, or imported JSON), and
- an AMBER topology + coordinates (`.prmtop` + `.rst7` from 5a), or a trajectory.

The main menu only says "Needs detected redox sites"; the topology requirement is not mentioned until you run it.

The ORCA route also needs ORCA's `orca_mm` program at preparation time. Without it, ORCA inputs **cannot be written**.

**Where to find it:** Full menu `5c`. Workflow mode: "simulate" stage.

Menu: `1. Run QM/MM preparation workflow`, `2. View current QM/MM configuration`, `3. Reset configuration to defaults`.

### Scenario 1: No topology yet — RUN
- `5c` → `1` → "Structure Loading" panel. A file browser opens on the project folder (no prmtop there) → `q`. Message: "Frame extraction failed: No prmtop file selected".
- `2` before any run → "No QM/MM workflow has been run yet." `3` → "QM/MM configuration reset."

### Scenario 2: Coordinates from rst7 (1HRC) — RUN
- `1` → "Found topology in workspace: final_transformed_structure.prmtop. Use this topology? (Enter to accept, or 'b' to browse)" → Enter.
- "Coordinate source: [r]st7/inpcrd or [t]rajectory (r)" → `r` → "Found coordinates in workspace…" → Enter.
- "Generated PDB: final_transformed_structure.pdb" (written in `qmmm_inputs/frames/`). Frame Extraction Summary table.
- "Select QM/MM interface: 1. Gaussian (ONIOM) 2. ORCA".

### Scenario 3: Gaussian ONIOM, 2 layers, full run on 1HRC — RUN
The ONIOM checklist has 7 steps, with the commands `[n]ext`, `[v]iew structure`, `[s]ave & quit`, `[q]uit`, `[#]` run/re-run a step, `[#s]` skip an optional step. State is saved in `oniom_inputs/workflow_state.json`.

1. **Select RedoxSite:** table of sites → `1` ("e.g. 1, 1-3, 1,3,5, or all"; several sites are merged into one QM region). "✓ Selected RedoxSite 1 with 73 atoms".
2. **Configure ONIOM Layers:**
   - Panel explaining HIGH (QM), MEDIUM (lower-level QM, 3-layer only) and LOW (MM).
   - "Number of ONIOM layers (2 or 3) (2)" → `2`.
   - "✓ Automatically added 7 residues from RedoxSite to HIGH layer: CYO14, CYO17, HIO18, MEO80, HMO105, PRD107, **WAT109**".
   - **BUG:** a crystal water was put in the QM region and one propionate (PRD106) was left out. The site's residue numbers (PDB numbering, where the PRDs are 108/110 with waters in between) are converted to topology positions by "number − 1", but the topology has a different residue order (PRD right after HMO). The user must check this list.
   - Workaround: "Add or remove whole residues? (n)" → `y`. Then `-109` → "✓ Removed WAT109", `106` → "✓ Added PRD106 (9 atoms)". `abc` → "Invalid input. Enter a residue number." `9999` → "Residue 9999 not found in topology". Enter to finish → "HIGH layer: 7 residues".
   - "Trim any residues to sidechain only (CA-CB cut)? (n)" → `y` lists "[1] PRD107 [2] PRD108". Here the numbering is +1 compared with the add/remove prompt just before (PRD106/107), so the two prompts in the same step number residues differently. The propionates are offered as trimmable because they have CA/CB atom names, which would be chemically wrong. Enter = skip.
   - "Freeze atoms beyond distance from QM region (Å, 0 for none) (0)" → `12`. Summary: "2-layer ONIOM, 7 HIGH residues, freeze beyond 12.0Å".
3. **Configure MM Settings & Atom Typing:**
   - Reads types and charges from the prmtop: "✓ Extracted 2121 atoms, 39 unique atom types".
   - HIGH 129 atoms, then "✓ Found 6 boundary atom(s)", "✓ Expanded HIGH by 12 atom(s) (formyl/amide caps …) 3 formamide cap(s) (Cα–C cut…) 3 primary-amide cap(s) (N–Cα cut…)".
   - Charge tables for the HIGH layer (per residue, running sum, total −1.1126) and the LOW layer (per residue type, +8.1126). Boundary table: 6 LOW Cα link hosts.
   - Writes `debug_atom_charges.txt` in the project folder.
   - "Boundary charge treatment options: 1 No adjustment, 2 Gaussian ScaleCharge keyword (recommended for EE), 3 Charge-shift redistribution, 4 Manual per-atom override" → `2`. A help panel explains the ScaleCharge digits (0 / 0.2 / … / 1.0 per bond distance; "digits must be monotonically decreasing from left to right"; default 500).
   - Typed `135` → "✓ Will add 'ScaleCharge=135'". **This was accepted although it breaks the stated rule** (no validation).
   - Then "Suggested QM region charge: -1 … Total system (real) charge: 7. Gaussian charge/multiplicity line: 7 1 -1 1". "View ONIOM layers in 3D structure viewer? (n)" → `n`.
   - Re-ran step 3 (`3` → "Re-run it? (n)" → `y`) and picked `3` charge-shift: "LYS13 CA: -0.2400 → 0.0000 (redistributed to CB, N, HA)". **However, ScaleCharge=135 from the first run stayed in the route line**, so both treatments were applied (see the file below).
4. **Configure QM Settings:** context panel (HIGH 141 atoms, total −1.1126, suggested −1).
   - "QM functional (e.g., B3LYP, M06-2X, wB97X-D) (B3LYP)" → `B3LYP`.
   - "QM basis set (e.g., 6-31G*, 6-311G(d,p), def2-SVP) (6-31G*)" → `def2-SVP`.
   - "QM region charge (HIGH layer total: -1.11) (-1)" → Enter.
   - Multiplicity guide (1 singlet, 2 doublet, 3 triplet) → `1`. **Wrong** for low-spin Fe(III), but ProPrep did not check it (see Scenario 4).
   - A yellow "Open-shell / BS-DFT advisory" panel always appears. It says ProPrep adds no broken-symmetry set-up and points to https://gaussian.com/afc/ and Guess=Read.
5. **Configure Job Settings (optional):**
   - "Job type [Opt/SP/Freq/Opt Freq] (Opt)" → `SP`; processors (4) → `8`; memory GB (8) → `16`.
   - Route preview: "#P ONIOM(B3LYP/def2-SVP:AMBER=SoftOnly)=(EmbedCharge,ScaleCharge=135) SP Geom=Connect".
   - "Extra keywords" → `SCF=XQC`.
6. **Prepare & Validate:** 6 sub-steps.
   - "Froze 359 atoms beyond 12.0Å", "✓ Placed 6 link atoms", "Built connectivity table: 2024 bonds".
   - Parameters extracted from the prmtop (78 bond, 187 angle, 331 dihedral, 39 VDW types).
   - "✓ Validation passed 0 warnings, 5 info messages".
7. **Write Gaussian Input File:** "1. Write production input file / 2. Write diagnostic file (ONIOM=OnlyInputFiles)" → `1`. "Output filename (oniom_1.com)" → Enter. "Summary report filename (oniom_summary_1.txt)" → Enter.
- Files (`runs/s5/hrc/cytc/oniom_inputs/`):
  - `oniom_1.com` (199 KB): Gaussian input. It contains `%mem`, `%nprocshared`, the route line `#P ONIOM(B3LYP/def2-SVP:AMBER=SoftOnly)=(EmbedCharge,ScaleCharge=135) SP Geom=Connect SCF=XQC`, and the charge/mult line `7 2 -1 1`. Each atom line has Element-Type-Charge, a freeze flag (0 active, −1 frozen), x y z, layer H/L, and link-atom info such as `H-HX-0.0 222`. These are followed by the connectivity, AMBER parameters (`AmbTrs`, `HrmBnd1`, `VDW`, `NonBon`).
  - `oniom_1.txt` and `oniom_summary_1.txt`: identical plain-language summaries (layers, frozen counts, method, HIGH residues, link atoms, validation). The first is written automatically, and then ProPrep asks for a name for the second.
  - `workflow_state.json`: resumable checklist state.
  - The `.com` also shows Cα link hosts with charge 0.000 (from the charge shift) together with ScaleCharge=135 (stale).
- **What the user needs to run it:** Gaussian 16 (licensed), e.g. `g16 < oniom_1.com > oniom_1.log`, with the memory/CPU stated in the file. ProPrep does not submit or run it.

### Scenario 4: Fixing multiplicity; diagnostic file — RUN
- The model region has 643 electrons, so it must have even multiplicity. ProPrep wrote real-system multiplicity 2 but model multiplicity 1 (`7 2 -1 1`), and Gaussian would reject this combination.
- Step 3's printout had said "Gaussian charge/multiplicity line: 7 1 -1 1", which is different from what was written to the file.
- Re-ran step 4 with multiplicity `2`. Re-running step 4 resets steps 5–7 to pending.
- Step 6 then `7` → `2` (diagnostic) → `oniom_1_diagnostic.com`. Its route line adds `OnlyInputFiles` and the charge line is `7 2 -1 2`.
- The message "Stop before the expensive SCF calculation (Test keyword)" does not match the file: there is no `Test` keyword in it.

### Scenario 5: Save & quit, resume, view configuration — RUN
- `s` → "State saved. You can resume later."
- Running Gaussian again later: a "Resume Available" panel ("Progress: 6/7 steps completed. Resume from saved state? (y)"). Answered `n` to start fresh.
  - **Even after declining to resume**, the earlier settings were silently reused: B3LYP/def2-SVP, SP, 8 cores, 16 GB, ScaleCharge=135 and SCF=XQC all appeared in the next file.
- Menu `2` after an ORCA run shows an "ORCA QM/MM Configuration" panel (method, charge/mult, QM/active atoms, embedding, scheme, job, rigid water).

### Scenario 6: 3-layer ONIOM with manual boundary charges — RUN
- Layers `3`. Kept the (wrong) automatic HIGH list (WAT included, PRD106 missing) to see what happens.
- MEDIUM residues `28`, `58` (THR, "MEDIUM layer: 2 residues"). Freeze `0`.
- MM step: still "6 boundary atoms". **The HMO–PRD106 covalent bond that now crosses the QM/MM boundary was not detected**, so it gets no link atom and no warning. Suggested QM charge 0.
- Boundary treatment `4` (manual): `0.0` for LYS13 CA → "✓ -0.2400 → +0.0000"; `x` → "Invalid number, keeping original". Summary: "Total System Charge: 7.2400 ⚠ Charge may be unbalanced". Validation then warned "Total charge (7.2400) deviates from nearest integer (7)" but still "✓ Validation passed".
- There is **no prompt for the MEDIUM-layer method**. It is fixed at HF/3-21G (route `ONIOM(B3LYP/def2-SVP:HF/3-21G:AMBER=SoftOnly)`).
- The charge line written was `7 2 0 1 0 1 0 2` (4 pairs). The intermediate-system multiplicity of 1 next to a model multiplicity of 2 looks inconsistent; a QM expert should check it.
- Step 5 skipped with `5s`. Files: `oniom_3layer.com`, `oniom_3layer.txt`, `oniom_3layer_summary.txt`.

### Scenario 7: ORCA path without ORCA installed — RUN / FAILED at writing
ORCA checklist steps: 1 Verify AMBER Topology, 2 Convert Force Field (orca_mm), 3 Select QM Region, 4 Define Active Region, 5 Configure QM Settings, 6 QM/MM Coupling (optional), 7 Job Settings (optional), 8 Preview & Write. State is saved in `orca_inputs/workflow_state.json`.

- Step 1: "Topology from upstream: 2121 atoms, 232 residues".
- Step 2: an explainer of `orca_mm -convff -AMBER`, then "orca_mm not found on PATH." "Path to ORCA installation directory (containing orca_mm)" → Enter → "✗ Step 2 failed: orca_mm not found in . Please ensure ORCA is installed…", followed by a **full Python traceback** on screen.
- Step 3, Select QM Region: site `1` → "Mapping RedoxSite atoms to prmtop indices (by residue ID)…" → "Warning: 5 atom(s) could not be mapped: CA:PRD110:A …". It is the same numbering problem as ONIOM, but here it is visible.
  - QM region table: CYO 14, CYO 17, HIO 18, MEO 80, HMO 105, PRD 107 (120 atoms).
  - Boundary panel ("Good cuts: Cα—C(=O), N—Cα, Cα—Cβ; Bad cuts: amide N—C(=O)"). For 7 boundary bonds, each with options:
    - amide bonds: `[1] Extend QM: add C,O (or N,H) of the neighbour`, `[2] Include entire residue`, `[0] Keep as-is`. Accepted the default `1` ×6.
    - Bond 7 "C2A(HMO105) — CA(PRD106) Other — Review": `[1] Include entire PRD106 in QM (9 atoms)` → `1`. **This fixes the missing propionate.**
  - Final: 6 boundary bonds, 141 QM atoms, "Suggested QM charge: -1", plus the warning "partial charge sum deviates from integer by 0.113. The charge alteration scheme will redistribute charges".
- Step 4: "Distance cutoff for active MM atoms (Å) (6.0)". Typed `abc` → "✗ Step 4 failed: could not convert string to float: 'abc'" with a traceback (no re-ask). Re-ran with `6` → "QM 141 + MM within 6.0 Å 613 = 754 active, 1367 frozen".
- Step 5: "QM method (e.g., B3LYP, PBE0, wB97X-D, DLPNO-CCSD(T))" → `PBE0`; basis (def2-SVP) → Enter; "Dispersion correction (D3BJ, D3, D4, or leave empty for none) (D3BJ)" → Enter; charge (−1) → Enter; multiplicity (1) → `2`. The multiplicity default is again 1.
- Step 6: "Embedding: [e]lectrostatic or [m]echanical (e)" → `e`. "Charge alteration scheme (CS, RCD, Z0, Z1, Z2, Z3) (CS)" → `RCD` → "Tune Scale_RCD parameter? (default 0.5)" → `y` → `0.5`. "Customize link atom (H) distances? (default: C=1.09, N=0.99, O=0.98 Å)" → `y` → Enter ×3. "Modify link atom double-counting settings?" → `y` → DeleteLADoubleCounting (y), DeleteLABondDoubleCounting (y).
- Step 7:
  - Job type "(Opt, EnGrad, NumFreq, SP)" → `Opt`; processors → `8`; memory per core (4000 MB) → Enter; extra keywords → `TightSCF RIJCOSX`; print level (0–4, default 1) → `1`.
  - "Use rigid MM water constraints? (n)" → `y`.
  - "Compute frozen-frozen non-bonded interactions? (y)" → `y` (the default is y below 10,000 atoms).
  - "Generate a QMMMSetup verification input (no actual calculation)? (n)" → `n`.
  - "Use PDB occupancy/B-factor columns instead of explicit atom lists? (n)" → `n`.
  - "Configure MM cutoffs and dielectric constant?" → `y`: Coulomb 12, LJ inner 10, LJ outer 12, dielectric 1.0 (defaults).
- Step 8: "Step orca-8 requires these steps first: orca-2. Run anyway? (may fail)" → `y` → "✗ Step 8 failed: [Errno 2] No such file or directory: ''" with a traceback. **No ORCA files are produced without orca_mm.**
- `q` → "Quit without saving? (Progress will be lost)" → `y`.

### Scenario 8: ORCA files when orca_mm is present — SOURCE-ONLY
From `orca_prep/orca_qmmm_preparator.py` and `orca_writer.py`, in `orca_inputs/`:
- `<prmtop>.prms`: ORCA force field from `orca_mm -convff -AMBER`.
- `<prmtop>_qmmm.inp`: `! QMMM [QMMMSetup] <method> <basis> <disp> <keywords> <job>`, `%pal nprocs`, `%maxcore`, and a `%qmmm` block (QMAtoms/ActiveAtoms ranges or Use_*_InfoFromPDB, ORCAFFFilename, Charge_Total, Mult_Total, Embedding, ChargeAlteration + Scale_*, Dist_*_HLA, double-counting flags, Rigid_MM_Water, Do_NB_For_Fixed_Fixed, PrintLevel). Optional `%mm` cutoffs/dielectric. Finally `* pdbfile <charge> <mult> <pdb>`.
- `<prmtop>_qmmm.pdb`: occupancy 1.0 marks QM atoms, B-factor 1.0 marks active atoms.
- `run_<prmtop>_qmmm.sh`: sets the ORCA path, checks that the files exist, runs `nohup orca … > .out &`, writes a `.pid`.
- For several frames: one set per frame plus `run_all_<prmtop>_qmmm.sh`.
- The final panel tells the user to `cd orca_inputs; ./run_…sh; tail -f …out` and to check the `.out`, `.xyz` and `.gbw` files.
- Possible issue (from source, worth an expert check): `Charge_Total`/`Mult_Total` are filled with the **QM-region** charge and multiplicity.

### Scenario 9: Coordinates from a trajectory — RUN (single frame works, multi-frame FAILED)
- To test, I made a fake 4-frame trajectory with cpptraj from the rst7 (`runs/s5/hrc/cytc/traj/fake_md.nc`; no MD was run).
- `t` → "Select trajectory file(s)" (file browser; a path can be typed) → `traj/fake_md.nc` → "Add another trajectory segment? (n)" → `n` → "Total frames: 4" → "Center solute in box (autoimage)? (y)" → `n`.
- Frame panel: the syntax shows as "start🔚increment" because `:end:` is turned into an emoji. Example rows such as "1:4:1 — Last 10 frames" are wrong for a 4-frame file.
- `1:4:2` (two frames) → "Running cpptraj..." → "**Frame extraction failed: cpptraj produced no output frames**". cpptraj had written `…_frame.rst7.1` and `…_frame.rst7.3` (and a single 2-model PDB), but ProPrep looks for `*_frame*.rst7`. So multi-frame extraction always fails.
- `4:4:1` (one frame) → "Extracted 1 frame(s)". The summary then listed rst7 `…_frame4.rst7` but pdb **`…_frame.pdb`**, the stale 2-model PDB left by the failed run, instead of `…_frame4.pdb`.
- When there are several frames, ProPrep would then ask which frames get inputs (not reached).

### Scenario 10: Other branches — SOURCE-ONLY
- `[v]iew structure` and "View ONIOM layers in 3D structure viewer" (browser).
- Actually trimming a residue to its side chain (moves CB+ to HIGH, cut at CA–CB).
- ScaleCharge default 500; boundary option 1 (no adjustment).
- Mechanical embedding; CS with Scale_CS tuning (default 0.06); Z0–Z3.
- QMMMSetup-only input; PDB-column atom selection.
- ORCA multi-frame output.
- "Use this topology? … 'b' to browse".

### Options reference — 5c
| Prompt | Choices | Meaning |
|---|---|---|
| Use this topology? | Enter accept / `b` browse | prmtop from workspace |
| Coordinate source | `r` rst7/inpcrd (default), `t` trajectory | |
| Trajectory | files, "Add another segment?", autoimage y/n, frames `start:end:increment` | multi-frame currently broken |
| Select interface | 1 Gaussian (ONIOM), 2 ORCA | |
| Checklist actions | n next, v view, s save & quit, q quit, `#` run/re-run, `#s` skip optional, d done | re-running a step resets later steps |
| ONIOM layers | 2 (default) or 3 | 3 adds a MEDIUM layer (fixed HF/3-21G) |
| Add/remove whole residues | number to add, `-number` to remove, Enter to finish | numbers = topology numbering as listed |
| Trim to sidechain | list numbers, `all`, Enter skip | QM/MM cut at CA–CB |
| Freeze distance | Å, 0 = none | frozen atoms get −1 in the .com |
| Boundary charge treatment | 1 none, 2 ScaleCharge (digits, default 500), 3 charge shift, 4 manual | |
| QM functional / basis | free text (B3LYP, 6-31G*) | passed straight into the route line |
| QM charge / multiplicity | integers (suggested charge; multiplicity default 1) | no parity check for the model region |
| Job type | Opt, SP, Freq, Opt Freq | Gaussian |
| Processors / Memory | integers (4 / 8 GB) | `%nprocshared`, `%mem` |
| Extra keywords | free text | appended; ScaleCharge/OnlyInputFiles go inside ONIOM() |
| Write | 1 production .com, 2 diagnostic (OnlyInputFiles) | plus a summary .txt |
| ORCA ORCA dir | path containing `orca_mm` | required |
| ORCA boundary bond | 1 extend, 2 whole residue, 0 keep | per boundary bond |
| ORCA active cutoff | Å (6.0) | non-numbers make the step fail |
| ORCA method/basis/dispersion | B3LYP / def2-SVP / D3BJ (or D3, D4, empty) | |
| ORCA embedding | e electrostatic, m mechanical | |
| ORCA charge scheme | CS (default), RCD, Z0–Z3; Scale_CS 0.06, Scale_RCD 0.5 | |
| ORCA link H distances | C 1.09, N 0.99, O 0.98 Å | |
| ORCA job | Opt, EnGrad, NumFreq, SP; procs; MB/core 4000; keywords; print level 0–4 | |
| ORCA extras | rigid MM water, frozen-frozen NB, QMMMSetup-only, PDB columns, MM cutoffs (12/10/12 Å, ε 1.0) | |

---

## Notes for the manual writer

**Jargon to explain:**
- redox site/center; c-type vs b-type heme (thioether-linked vs non-covalent); axial ligand; propionate.
- transformer (a ProPrep recipe that renames atoms for a force-field library) vs template (your tag for a repeatable detection recipe).
- microstate (one combination of oxidation states across all sites) and its 2-letter-per-site code (O/R + L/H/I).
- fixed_pH vs constant_pH; broken-symmetry (BS) state of Fe-S clusters.
- QM/MM, ONIOM layers (HIGH/MEDIUM/LOW), link atom, electrostatic vs mechanical embedding, ScaleCharge, charge shift, multiplicity, active vs frozen atoms.

**Recommended beginner path for a heme protein → QM/MM:**
1. 1a load.
2. 2a detect with defaults. With more than one site, say **yes to templates** and give the transformer name as the site type: this is the only way ProPrep shows you which bonds the transformer needs. Use residue-count search (HIS,CYS 2,2; add MET for His/Met hemes), add all residues, and define exactly the bonds listed under "Bond Requirements" (same-residue pairs like `105-105`).
3. Export JSON.
4. 4d option 2, choose the transformer, pick the oxidation state, fixed_pH, deprotonated propionates.
5. 5a topology.
6. 5c Gaussian.
7. **Check the HIGH-layer residue list by name.**

For a single-site protein, the bond requirements are never shown in 2a. Use 4d `why 1` to see them, then add the bonds in 2a → Edit → bonds → add.

**Important behaviours/bugs (describe, not fixed):**
1. 5c maps site residues to topology residues by number. When the transformed PDB has waters between the heme and its propionates (as with 1HRC), ONIOM silently puts a water in the QM region and leaves a propionate out. In ONIOM 3-layer mode, the cut HMO–PRD bond is then not even treated as a boundary. ORCA warns about unmapped atoms and offers to include the residue, which works around it.
2. In ONIOM, the multiplicity default is 1 and the model-region multiplicity is never checked against the electron count. ProPrep computed the real-system multiplicity (2) itself, so the file contained inconsistent values (`7 2 -1 1`) until it was fixed by hand.
3. ScaleCharge digits are not validated (`135` accepted). Re-running the MM step does not clear an earlier ScaleCharge. Declining "Resume" still reuses earlier settings.
4. ORCA inputs cannot be generated at all without `orca_mm` (part of ORCA). Invalid numbers in ORCA steps abort the step with a Python traceback.
5. Multi-frame trajectory extraction in 5c always fails (cpptraj output naming), and a stale PDB is then picked for single frames.
6. In the 2a editor, a non-number at any ResID prompt crashes out of the editor and loses all unsaved edits. Chain/residue names are not validated.
7. The bond tool accepts absurd bond lengths (8.2 Å) without warning, and target numbers are re-sorted by distance for every source atom.
8. False-positive "cofactors": ACE caps, glycerol (GOL), acetate (ACT). Deselect them at center selection.
9. 4d renumbers residues (sequentially for hetero groups or whole chains), and the progress lines and mapping report use intermediate IDs (e.g. `A_HEC_1299` for heme 1294, mapping report 105 vs PDB 106). Tell users to rely on the `_redox_sites_updated.json` and the final PDB, not on the printed IDs.
10. After one generation, 4d blocks all of its options for the rest of the session.
11. Cosmetic markup losses:
    - "press  to return" (missing `[m]`) on every menu
    - "Choice  (a)"
    - "Selection Syntax: :,,..."
    - "start🔚increment"
    - "state  -> HXO" in the transformer summary
    - "Using RCSB PDB PDB file" printed twice
    - "Successfully imported 0 redox site(s)" followed by "No sites found"
    - "unsaved changes" when nothing changed
    - The diagnostic-file message mentions a Test keyword that is not in the file.
12. Help → Common Workflows (source, `menu_commands.py`) still shows old keys ("4a.RedoxSite … 5d.QM/MM"). The current keys are 4d and 5c.
13. Environment note: tLEaP/teLeap and parmed ran from `etanalyze_env` (first on PATH), not from the ProPrep env. The manual should tell users to activate the ProPrep environment so the right AmberTools are used.
14. Redox-site JSON import lists ProPrep session JSON files as candidates. Choose the `<PDB>_redox_sites.json` or `_redox_sites_updated.json` file.

**Software needed beyond ProPrep:** Gaussian 16 to run `.com` files; ORCA (including `orca_mm`) to *prepare* and run ORCA inputs. A Cu (or other unsupported metal) site needs parameters from 2c/MCPB.py (other scope).

**Tally:** 34 scenario blocks. **RUN 26, SOURCE-ONLY 5, FAILED 3** (2a S5 editor crash; 5c S7 ORCA writing blocked; 5c S9 multi-frame trajectory). 4d S8 counts as RUN; only its save step was source-only. Many RUN blocks contain several sub-paths and failing sub-steps, described in place (e.g. invalid inputs, 5c S7 step 4 invalid number).
