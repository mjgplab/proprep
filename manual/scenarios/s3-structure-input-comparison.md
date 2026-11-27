# Getting a structure in, and comparing structures (scope s3)

Recorded with ProPrep 1.23.0, started as `proprep --full-menu`, on 2026-11-21.
Sessions (transcripts under `manual/runs/_sessions/<name>/transcript.txt`):

| Session | Project folder | What it covered | Replay file |
|---|---|---|---|
| `s3-main` | `runs/s3/main/s3load/` | 1a: every source, searches, errors, metadata, local files, the "view together" viewer | `proprep_session_20261121_122922.json` (the session ended when ProPrep crashed in 1b, see Scenario 1b-4) |
| `s3-asm` | `runs/s3/asm/asm/` | 1b on 1HHO, switching structures, 3a BLAST | `proprep_session_20261121_124513.json` |
| `s3-align` | `runs/s3/align/align/` | 3b sequence, TM-align and manual alignment; ua viewer | `proprep_session_20261121_125230.json` |
| `s3-misc` | `runs/s3/misc/misc/` | search filters, single download, electron density | `proprep_session_20261121_131727.json` |
| `s3-offline` | `runs/s3/offline/` | `--pdbid` on the command line, and "no internet" errors (faked with a dead proxy) | `proprep_session_20261121_131031.json` |

Notes on reading the screens below:
- A prompt shows its default in brackets, e.g. `Select source (1):`. Pressing Enter accepts it.
- Every sub-menu ends with `[b]ack · [m]main · [x]exit`. `x` asks `Are you sure you want to exit? [y/n] (n)`.

---

## Feature: 1a. Structure Loader

**What it is for (plain English):** This is how you get a protein structure into ProPrep. It can download an experimental structure from the Protein Data Bank (RCSB), a predicted structure from the AlphaFold database, or an AlphaFold model with ligands copied in from AlphaFill. It can also open a PDB file already on your computer. Every later tool works on what you load here.

**Needs before you start:** An internet connection for RCSB, UniProt, AlphaFold and AlphaFill. Nothing is needed for a local file.

**Where to find it:** Main menu `1a` (Full Menu). In Workflow Mode it is the first step of each workflow ("1a.Load → ..."). On the command line, `--pdbid <ID>` downloads a PDB entry at start-up and `--pdbfile <path>` loads a local file (see Scenario 1a-20).

The Structure Loader menu:
```
══ Structure Loader Menu ══
  1. ✓ Load protein structure
  2. ✓ Load AMBER topology & coordinate files
  3. ○ View structure metadata
       ⚠ Load a structure first
  4. ✓ Load a ligand file or library (SDF/mol2/SMILES) for docking
→ Suggestion: Find a structure by name or PDB ID across RCSB, AlphaFold, and AlphaFill, or load a local file (option 1)
```
Option 1 opens the source menu, which keeps coming back after each load so you can load several structures in a row:
```
═══ Load Structure - Select Source ═══
1. Query Protein Data Bank (RCSB)
2. Query AlphaFold Database
3. Query AlphaFill Database (AlphaFold + ligands/cofactors)
4. Load local PDB file
5. View local files together (spatial comparison, view-only)
6. Done (return to the Structure Loader menu)
```
RCSB (source 1) then asks for a method:
```
1. Enter PDB ID directly
2. Search by entry title keyword
3. Search by gene/protein (UniProt)
4. Search by ligand (code or name)
5. Back to source selection
```

### Scenario 1a-1: Metadata before anything is loaded — RUN
- Starting point: new project `s3load`, Structure Loader menu.
- Steps: `Enter your choice (1):` -> typed `3`.
- What the user sees: `⚠ Load a structure first`, and the menu again.
- Session log: runs/_sessions/s3-main/transcript.txt

### Scenario 1a-2: Download by PDB ID, asymmetric unit (4HHB) — RUN
- Starting point: Structure Loader menu.
- Steps:
  1. `Enter your choice (1):` -> `1` (Load protein structure)
  2. `Select source (1):` -> `1` (RCSB)
  3. `Select method (1):` -> `1` (Enter PDB ID directly)
  4. `Enter PDB ID (e.g., 1ABC, 7XYZ):` -> `4hhb` (lower case works; ProPrep changes it to upper case)
  5. `Structure type: 1. Asymmetric unit (default - as deposited) 2. Biological assembly (functional oligomeric form)` / `Select structure type (1):` -> `1`
- What the user sees: `Downloading 4HHB...`, `✓ Downloaded 4HHB`, then a short summary: PDB ID, Classification (`OXYGEN TRANSPORT`), deposition date, title, experimental method (`X-RAY DIFFRACTION`), resolution (`1.74 Å`), and the primary publication (title, authors, journal, PubMed and DOI links), ending with `Abstract and free full text, if any: View structure metadata`. Then back to the source menu.
- Result / files: `runs/s3/main/s3load/4HHB.pdb`, the entry as downloaded.
- Time: about 2 seconds.
- No chain or component selection is offered at load time. Choosing chains, ligands and waters is done later in 4a (PDB Filter). The 3b aligner has its own chain choice.
- Session log: runs/_sessions/s3-main/transcript.txt

### Scenario 1a-3: Recent entry with hydrogens and odd numbering (9ZZZ) — RUN
- I typed `9ZZZ` expecting an invalid ID, but it is a real 2026 cryo-EM entry (a nucleosome with PARP2). This shows two automatic checks.
- What the user sees:
  - A **Residue Numbering Warning** box: `The following chains do not start with residue number 1: Chain A: First residue is LYS #13 ...`, with the possible reasons (missing N-terminal residues, leader-sequence numbering, reference numbering).
  - **Hydrogen Atom Removal**: `This structure contains 14334 hydrogen atoms. Hydrogen naming in PDB structures is often inconsistent with force field expectations. It is recommended to remove them and let tLEaP add them back ...` / `Remove hydrogen atoms? [y/n] (y):` -> `y`
  - `H-stripped structure saved to: 9ZZZ_noH.pdb`, `Removed 14334 hydrogen atoms`
- Result / files: `9ZZZ.pdb` (original) and `9ZZZ_noH.pdb` (the same structure without hydrogens; no header records, coordinates only).
- Problems:
  - The answer is remembered for the whole session. Later loads that contain hydrogens are stripped without asking, with the message `(using previous selection)`. The 1AKK batch download in Scenario 1a-9 was not stripped at all (batch downloads skip this step).
  - **Important:** after this step the `_noH` file takes priority over every structure you load later, for every downstream tool. See Scenario 1b-4 and the notes for the writer.
- Session log: runs/_sessions/s3-main/transcript.txt

### Scenario 1a-4: PDB ID that does not exist (XQ99) — RUN
- Steps: `1` > `1` > `XQ99` > `1`.
- What the user sees (log lines, then a red line):
  ```
  ERROR - PDB XQ99 not found (404). May be a large structure only available in mmCIF format.
  WARNING - RCSB only provides mmCIF format for structures with >99,999 atoms.
  WARNING - ProPrep's metal site transformation currently requires PDB format.
  Failed to download PDB XQ99
  ```
  Back to the source menu. The ID format is not checked first: ProPrep asks asymmetric unit or assembly before it tries the download.
- Problems: the message suggests the entry may be too large, but here the ID simply does not exist. The writer should explain both meanings. Very large entries (ribosomes, viruses) that exist only as mmCIF cannot be loaded from RCSB.

### Scenario 1a-5: Malformed ID ("hello world") and empty ID — RUN
- `hello world` -> `Downloading HELLO WORLD...`, `ERROR - Error downloading PDB file HELLO WORLD: URL can't contain control characters. '/download/HELLO WORLD.pdb' (found at least ' ')`, `Failed to download PDB HELLO WORLD`. The message is technical, but nothing breaks.
- Enter with nothing typed -> `No PDB ID entered`, back to the source menu.

### Scenario 1a-6: Download a biological assembly directly (1HHO) — RUN
- Why 1HHO: its asymmetric unit (what the crystallographers deposited) holds only one alpha and one beta chain. The working hemoglobin is a tetramer of four chains.
- Steps: `1` > `1` > `1HHO` > `Select structure type (1):` -> `2` (Biological assembly)
- What the user sees: `Checking available biological assemblies for 1HHO...`, `Found 1 biological assembly`, `Downloading biological assembly 1...`, `✓ Downloaded biological assembly 1 for 1HHO`, the summary, then `Assembly contains 4 chains, 4,792 atoms`. When an entry has several assemblies, ProPrep lists them as `1. Assembly 1` and so on, and asks `Select assembly number (1)` (source; 1HHO has one). With no assembly at all it prints `No biological assemblies available for <ID>` / `Falling back to asymmetric unit...`.
- Result / files: `1HHO_assembly1.pdb`.
- Problems:
  - The summary showed no **Resolution** line, because the assembly file has no resolution record. The asymmetric-unit download shows `Resolution: 2.1 Å`.
  - The RCSB assembly file stores the tetramer as **two MODELs, each with chains A and B**. It does not contain four distinct chains A to D, and it keeps the waters (221 HOH). The 1b generator (Scenario 1b-1) instead writes chains A, B, C, D in one model. Downstream tools that read only the first model would see only two chains. The writer may want to recommend 1b for MD work. (I did not test how later tools treat the two-model file.)

### Scenario 1a-7: Search by title keyword, with organism filter, paging, metal filter — RUN
- Steps:
  1. `Select method (1):` -> `2` (Search by entry title keyword)
  2. `Enter search keyword (e.g., "carbonic anhydrase", "SOD1", "hemoglobin"):` -> `cytochrome c`
  3. `Filter by organism? (optional, press Enter to skip) ():` -> `Equus caballus`
  4. `35 structures match` / `Maximum number of structures to list, most relevant first (range: 1-10000, default: 100) (100):` -> Enter
- What the user sees: `Found 35 structures` and a table, 20 per page: `#`, `PDB ID`, `Title`, `Method` (X-ray/NMR/Cryo-EM), `Res (Å)` (`N/A` for NMR), `Metals` (e.g. `HEC (Fe)`, `HEC (Fe)  ZN`), followed by a legend explaining the Metals column. Then:
  ```
  Options:
    • Select: single (1), multiple (1,3,5), range (1-3), or 'all'
    • Type 'next' to see next page
    • Type 'filter' to refine search
    • Type 'back' to return to the Structure Loader menu
  Selection:
  ```
  The table is sorted best resolution first, with NMR entries last.
- `next` shows page 2, numbered **1-15 again** (numbers restart on each page and refer to that page only). `prev` goes back.
- `filter` opens `═══ Refine Search ═══ 1. Maximum resolution 2. Experimental method (X-ray/NMR/Cryo-EM) 3. Must contain metal 4. Clear all filters 5. Back to results`. I chose `3`, then `Metal type (e.g., Zn, Fe, Cu) or leave blank for any metal:` -> `Zn`. Result: `Filtered to 4 structures containing Zn, as an ion or within a cofactor` (1WEJ, 5IY5, 1U75, 1M60).
- `filter` > `4` (Clear all filters) -> `To clear filters, start a new search`. The list stays filtered.
- Invalid selection: `abc` -> `Invalid selection: 'abc'. Use numbers (1-4), ranges (1-3), or 'all'`. `9` gives the same message.
- `back` prints `Back to the Structure Loader menu` but actually goes back to the **source** menu (cosmetic).
- Problem: a search with several words is very broad. `horse heart cytochrome c` (no organism) gave `10104 structures match`, although the top results were relevant. Tell users to add the organism filter.

### Scenario 1a-8: Search, resolution and method filters, then download one result (crambin) — RUN
- Session `s3-misc` (project `runs/s3/misc/misc/`). Keyword `crambin`, no organism, default 100 results.
- `filter` > `1` (Maximum resolution) > `abc` -> `Invalid resolution value` (the list is unchanged).
- `filter` > `2` (Experimental method) > `1. X-ray Diffraction 2. Solution NMR 3. Cryo-EM` -> `2` -> `Filtered to 9 structures`. `Selection:` -> `1` -> `Selected: 1YVA - NMR solution structure of crambin in DPC micelles`, `Method: SOLUTION NMR`, `Download 1YVA? [y/n] (y):` -> `n` -> `Download skipped`, back to the source menu.
- New search `crambin` > `filter` > `1` > `1.0` -> `Filtered to 12 structures with resolution ≤ 1.0 Å`. `16` -> `Invalid selection: '16'. Use numbers (1-12) ...`. `1` -> `Selected: 3NIR - Crystal structure of small protein crambin at 0.48 A resolution`, `Resolution: 0.48 Å`, `Method: X-RAY DIFFRACTION`, `Download 3NIR? (y)` -> `y` -> the normal summary, then `This structure contains 315 hydrogen atoms ... Remove hydrogen atoms?` -> `n`.
- Result / files: `runs/s3/misc/misc/3NIR.pdb`.
- A single download from search results always fetches the asymmetric unit. There is no assembly choice on this route.

### Scenario 1a-9: Search, then download several results at once (1HRC + 1AKK) — RUN
- Steps: keyword `horse heart cytochrome c`, organism blank, maximum `10`, `Selection:` -> `3,9`
- What the user sees: a summary table `Selected 2 structures: 1HRC ... 1.90 Å / 1AKK ... N/A`, then `Download all 2 structures? [y/n] (y):` -> `y`, a progress bar, `Successfully downloaded 2 structure(s): ✓ 1HRC ✓ 1AKK`.
- Result / files: `1HRC.pdb`, `1AKK.pdb`.
- Problems: a batch download shows no summary or publication, gives no hydrogen-removal prompt (1AKK is an NMR model with 878 H that stayed in) and no numbering check. Only the **last** one (1AKK) becomes "the" RCSB structure that later tools use. The others are kept in a list that only some tools look at.

### Scenario 1a-10: Search by gene/protein through UniProt (SOD1, human) — RUN
- Steps: `Select method` -> `3`; `Enter gene name, protein name, or UniProt ID (e.g., "SOD1", "P00441"):` -> `SOD1`; `Filter by organism? (optional, e.g., human, mouse) ():` -> `human`.
- What the user sees: `═══ UniProt Entry Found ═══` with accession P00441, entry name, protein (Superoxide dismutase [Cu-Zn]), gene, organism, length 154, mass, function (shortened), cofactors (`Cu cation, Zn(2+)`), disulfide bonds, diseases. Then `View detailed UniProt metadata? [y/n] (y)` -> `y` shows sections IDENTITY, SEQUENCE, FUNCTION, COFACTORS, DISULFIDE BONDS, BINDING SITES, POST-TRANSLATIONAL MODIFICATIONS, DISEASE ASSOCIATIONS and DATABASE CROSS-REFERENCES (`PDB: 156`, ...). Then `This entry has 156 associated PDB structures` / `Select a PDB structure to download? [y/n] (y)` -> `y`.
- A table `PDB Structures for Superoxide dismutase [Cu-Zn] (Page 1/8)` follows, with `#`, `PDB ID`, `Method`, `Resolution`, `Chains` (e.g. `A/F=2-154`, i.e. which chains and residues the entry covers), sorted by resolution, 20 per page, with the options `Enter number (1-20)`, `next`, `back`.
- `Selection:` -> `4` (1MFM). It downloaded **immediately**: no "Download?" confirmation and no asymmetric-unit/assembly question. The usual summary followed.
- When several UniProt entries match, a numbered table is shown with `'b' to go back` (source).
- Result / files: `1MFM.pdb`.

### Scenario 1a-11: Search by ligand (tetrahydrobiopterin) — RUN
- Steps: `Select method` -> `4`. The screen explains: `Every ligand in the PDB has a code (HEM, H4B, ...). Enter a code, or words from the ligand's name to look the code up. Names match on whole words: 'biopterin' does not find 'tetrahydrobiopterin'.`
  1. `Enter ligand code or name (e.g., "HEM", "H4B", "tetrahydrobiopterin"):` -> `tetrahydrobiopterin` -> table `Ligands matching your entry` (`H4B` 736 structures, `BHS`, `H2B`, `4AB`), with columns Code, Name, Formula, Found by, Similarity, Structures.
  2. An explanation that one compound can have several codes (redox and protonation states, isomers). `Ligand code to find related ligands for (press Enter to skip) ():` -> `H4B` -> table `Ligands related to H4B`: same formula (BHS), `same atoms, -2 H` (HBI = 7,8-dihydrobiopterin, H2B), `-4 H`, and `similar structure` with a fingerprint similarity score. Footnote: `Similarity is RCSB's fingerprint score (1 = identical fingerprints). ... read the names.`
  3. `Ligand codes to search for, separated by commas (or 'back') ():` -> `BHS, XYZQ` -> `Not a PDB ligand code: XYZQ` (asked again) -> `BHS,HBI` -> `BHS: 1 structures`, `HBI: 28 structures`, `29 structures contain BHS or HBI`.
  4. `Maximum number of structures to list, best resolution first (... default: 100)` -> `5` -> `Listing the 5 best-resolved of 29 structures`, a warning that depositors' ligand codes are not always right, then the results table with two extra columns, **Ligand** and **RSCC** (how well the ligand fits the electron density, 1 = perfect; explained under the table).
  5. `back`.
- Notes: when you type a name rather than a code, the "Ligand codes to search for" prompt has an empty default. When you typed an exact code, that code is the default.

### Scenario 1a-12: AlphaFold by UniProt ID — bad ID, then P69905 — RUN
- Steps: `Select source` -> `2`; `1. Retrieve by UniProt ID 2. Search by gene name 3. Back` -> `1`; `Enter UniProt ID (e.g., P00520, Q5VSL9):`
  - `P99999X` -> `Error checking AlphaFold availability for P99999X: HTTP Error 400: Bad Request`, `No AlphaFold structure found for P99999X`, `Note: Not all proteins have AlphaFold predictions yet`.
  - `p69905` -> `✓ Retrieved AF-P69905-F1`, `Protein: Hemoglobin subunit alpha`, `Organism: Homo sapiens`, `File: AF-P69905-F1-model_v6.cif`, `Confidence: 98.0 (Very High)`, then `View AlphaFold metadata? [y/n] (y)` -> `y`.
- Metadata screen: AlphaFold ID, UniProt ID, protein, organism, **Length**, Model Quality (mean pLDDT 98.0, Very High), Confidence Distribution (Very High >90 / Confident 70-90 / Low 50-70 / Very Low <50), and a per-residue confidence map: the sequence in rows of 60, coloured by confidence (colour does not survive in the transcript).
- Result / files: `AF-P69905-F1-model_v6.cif`. This is **mmCIF**, not PDB.
- Problems / likely bugs:
  - **Length says 1077 residues, but the protein has 142 residues.** 1077 is the number of atoms. The confidence distribution is also counted per atom (`1069 residues (99.3%)` + `8 residues`), so the percentages are atom-weighted.
  - The downloaded file is `.cif`. Several file pickers in ProPrep list only `.pdb` files and reject `.cif` (`... is not a PDB file`), so this file cannot be picked by hand in the 3b aligner (Scenario 3b-3) or in 1a source 4.

### Scenario 1a-13: AlphaFold by gene name (HBB) — RUN
- Steps: `2` > `2` > `Enter gene name (e.g., SRC, TP53, BRCA1):` -> `HBB` > `Enter organism (optional, e.g., human, mouse) :` -> Enter.
- What the user sees: two INFO log lines, then `Found 10 match(es):` with a table (UniProt ID, Protein Name, Organism, Reviewed ✓): human P68871, cow, rat, mouse... `Select entry (1-10) or 'b' to go back (1):` -> `1` -> `Retrieving AlphaFold structure for P68871...`, `✓ Retrieved AF-P68871-F1`, `View AlphaFold metadata? (y)` -> `n`.
- Result / files: `AF-P68871-F1-model_v6.cif`. It **replaces** the earlier AlphaFold structure (P69905) as the workspace's AlphaFold structure. ProPrep keeps one AlphaFold structure at a time.
- Problem: metadata viewed later (Scenario 1a-18) shows `UniProt ID: N/A`, `Protein: Unknown`, `Organism: Unknown` for a structure fetched this way, and again atoms counted as residues (`Length: 1131 residues`).

### Scenario 1a-14: AlphaFill by UniProt ID (P69905) — RUN, wrong result
- Steps: `Select source` -> `3`. The screen says `AlphaFill provides AlphaFold structures enriched with ligands, cofactors, and metal ions`, then `1. Retrieve by UniProt ID 2. Search by gene name 3. Back` -> `1`; `Enter UniProt ID (e.g., P00441, Q5VSL9):` -> `P69905`.
- What the user sees: `✓ Retrieved AF-P69905-F1`, `File: ./alphafill_P69905.cif`, **`No transplanted ligands/cofactors found`**. Metadata: `Length: 168 residues`, `No transplanted ligands/cofactors found`, `This AlphaFill structure does not have enrichment data`.
- **Bug:** the file does contain transplanted ligands. It has 438 HETATM lines, including heme (HEM, 172 atoms), CMO, OXY, NO, HNI and others, and AlphaFill's own JSON for P69905 lists 22 hits. ProPrep reads the wrong field of AlphaFill's reply, so it always reports "none". The "168 residues" count is 142 amino acids plus the transplanted ligand groups. **Users should not believe the "no ligands" message.** The file itself is fine.
- Gene-name path (option 2): same flow as AlphaFold option 2 (source, not run separately).

### Scenario 1a-15: Local PDB file: browser, bad paths, wrong file types — RUN
- Setup: I downloaded `1CRN.pdb` and `1UBQ.cif` with curl into `runs/s3/localfiles/` and created `junk.txt`.
- Steps: `Select source` -> `4` -> `Load Local PDB File` and a file browser:
  ```
  Current directory: .../s3load
    [ 1] 📁 .. (parent)
    [ 2] 📄 1AKK.pdb  (160.2 KB) ...
  Commands:
    N          - open dir N, or select file N
    ..         - go up one level
    <path>     - jump to a file or directory path
    q          - cancel
  Selection (q):
  ```
  - `/nonexistent/foo.pdb` -> `Path not found: /nonexistent/foo.pdb` (browser stays open).
  - `/…/runs/s3/localfiles` (a folder) -> browser moves there. It lists only `1CRN.pdb`, hiding `.cif` and `.txt`.
  - `junk.txt` (a bare file name) -> `Enter a single number, '..', or 'q'.` Bare names are not accepted; type a full path, or `./name`.
  - full path to `junk.txt` -> `junk.txt is not a PDB file.`; full path to `1UBQ.cif` -> `1UBQ.cif is not a PDB file.`
  - `q` -> asks `Enter path to PDB file:` (a manual-entry fallback) -> `/tmp/doesnotexist.pdb` -> `File not found or invalid path`, back to the source menu.
  - Then `4` > the full path to `1CRN.pdb` -> `✓ Loaded structure from local file`, `File: 1CRN.pdb`.
- Result / files: ProPrep **copies** the chosen file into the project folder (`s3load/1CRN.pdb`) and works on the copy.
- No hydrogen prompt appeared for 1CRN: it has no hydrogens. A file with hydrogens gets the same prompt as Scenario 1a-3, unless that question was already answered in the session.
- If the file name contains `AF-` or `alphafold`, ProPrep treats it as an AlphaFold model rather than an experimental structure (source).

### Scenario 1a-16: mmCIF file through the manual-path fallback — RUN (loads, but only half works)
- Steps: `4` > `q` > `Enter path to PDB file:` -> full path to `1UBQ.cif`.
- What the user sees: `Note: BioPython parser encountered an issue with this file.` / `Storing file path only - structure parsing will be done by individual modules.` / `✓ Loaded PDB file (metadata extracted, BioPython parsing skipped)`.
- So the browser refuses a `.cif`, but the fallback accepts it, and the result is a "local structure" with no parsed atoms. Tell users to load PDB-format files. I did not test which later tools fail on it.

### Scenario 1a-17: View local files together (spatial comparison) — RUN
- Steps: `Select source` -> `5`. The screen explains: `Pick two or more local PDB files; they open overlaid in the viewer with their real relative coordinates. Enter '1,3,5' or '1-5,7,9' ... View-only — this does not load them for preparation.` The browser now also accepts `N,M  N-M` (several files). `2,5` -> `Selected: 1AKK.pdb`, `Selected: 1HRC.pdb`, `Added (1)`, `Added (2)`; `Add another file to view? [y/n] (n)` -> `n` -> `Opening 2 structures overlaid in the viewer...`, `✓ Structure viewer launched at http://localhost:8766/viewer (2 structures)`.
- The structures are shown in their **original** coordinates. They are not superimposed, so two unrelated files can appear far apart. With fewer than two files it says `Need at least two files for a spatial comparison. Nothing opened.` (source).
- `curl http://localhost:8766/viewer` -> HTTP 200, the HTML page "ProPrep Interactive Structure Viewer".
- The port was 8766, because another ProPrep (another agent) already had 8765. ProPrep tries 8765 and then the next free port.

### Scenario 1a-18: View structure metadata (several structures loaded) — RUN
- Steps: `Source menu` -> `6` (Done). The menu now reads `1. ● Load structure from PDB/AlphaFold/AlphaFill [All sources loaded: RCSB, Local, AlphaFold, AlphaFill] ✓`. `3` -> `Multiple structures available: 1. RCSB PDB structure 2. Local PDB structure 3. AlphaFold structure 4. AlphaFill structure` / `View metadata for (1):` -> `1`.
- Opening metadata also **launches the viewer** automatically (`✓ Structure viewer launched at ... (1 structure)`).
- The RCSB structure shown was **1MFM**, the most recently downloaded RCSB entry, not 4HHB.
- PDB METADATA VIEWER sections (each tested where marked):
  `1 Summary` (tested) · `2 Safety & Status Records` · `3 Header Records` · `4 Secondary Structure` · `5 Sequence Information` (tested) · `6 Database References` · `7 Citations` (tested) · `8 Hetero Compounds` (tested) · `9 Crystallographic Information` · `10 Site Information` · `11 Remark Records` · `12 View All` · `13 Export Metadata` (tested) · `14 Help` · `15 Exit`
  - **Summary** gives the textual composition description: `Compound: MOL_ID: 1; MOLECULE: PROTEIN (COPPER,ZINC SUPEROXIDE DISMUTASE); CHAIN: A; ... MUTATION: YES`, `Source: ... ORGANISM_SCIENTIFIC: HOMO SAPIENS; ... EXPRESSION_SYSTEM: ESCHERICHIA COLI ...`, keywords, crystal cell, revision history, counts (`Hetero Compounds: 4 unique types`, `Sites: 15 defined`, helices/sheets/disulfides/links), `Chain A: 153 residues`, and which REMARKs exist (`Including biological assembly information (REMARK 350)`). The quaternary-structure wording (e.g. "TETRAMERIC") appears in 1b, from REMARK 350.
  - Summary bug (cosmetic): revision history reads `1 EMARK` instead of `REMARK` (first letter cut off).
  - **Hetero Compounds**: per ligand code, the name, formula, occurrences and location (`CD CADMIUM ION 9(CD 2+) ... Chain A 156, 1 atoms`, `CU`, `ZN`, `HOH 283(H2 O)`).
  - **Sequence Information**: the three-letter sequence per chain, 10 per row.
  - **Citations**: authors, title, reference, PubMed and DOI links, then `Abstract and full-text availability:`. Here Europe PMC timed out after 15 s: `Europe PMC could not be reached (... Read timed out. (read timeout=15))`. Not fatal.
  - **Export Metadata**: `1. Export as JSON 2. Export as plain text summary 3. Return` -> `1` -> `Enter output filename (1MFM_metadata.json):` -> Enter -> `Exported metadata to 1MFM_metadata.json` (file in the project folder).
- Choice `2` (Local) showed the same viewer for `1CRN.pdb`. Choice `3` (AlphaFold) showed the AlphaFold metadata screen (see the bugs in 1a-12/1a-13). Choice `4` shows the AlphaFill screen.
- Cosmetic bug: the suggestion line reads `... or press  to return to the main menu`. The `[m]` key name is missing (it was swallowed by text formatting). It appears in several menus.

### Scenario 1a-19: Loading a second structure / switching structures — RUN
- Loading the same kind again **replaces** the previous one: a new RCSB download replaces the old RCSB structure (4HHB → 1HHO → 1MFM), and a new AlphaFold download replaces the old one. Different kinds (RCSB, Local, AlphaFold, AlphaFill) are kept side by side, and the menu lists which are loaded.
- **Which structure later tools use** is decided by a fixed priority, not by what you loaded last. Files that ProPrep made itself (the hydrogen-stripped `_noH` file, an assembly from 1b, aligned files from 3b, filtered/repaired/protonated files) all outrank a fresh download. Two runs showed this:
  - In `s3-main`, after loading 9ZZZ (H-stripped) and then 4HHB, 1HHO, 1HRC, 1AKK, 1MFM, 1CRN..., tool 1b said `Using H-Stripped PDB file: .../9ZZZ_noH.pdb`.
  - In `s3-asm`, after generating 1HHO's assembly in 1b and then loading 4HHB, 1b again used `1HHO_biological_assembly.pdb`.
- Some tools (3a BLAST, ua Viewer) show a table `Available Structures in Workspace` and let you choose. 1b does not. Workspace Options (`w`) has no "choose active structure". Its only remedy is `5. Reset workspace` or starting a new project.

### Scenario 1a-20: Command-line `--pdbid` and "no internet" — RUN
- Started as `proprep --full-menu --pdbid 1CRN` with every web request forced to fail (a dead proxy setting), to imitate a computer with no internet.
- With `--pdbid`, **no project-name question** is asked. ProPrep uses the current folder.
- What the user sees: `Downloading 1CRN...`, `ERROR - Error downloading PDB file 1CRN: <urlopen error [Errno 61] Connection refused>`, `Failed to download PDB 1CRN`, then a contradictory log line `INFO - Downloaded and loaded PDB: 1CRN`. The menu then shows nothing loaded.
- Keyword search offline: `The RCSB search service could not be reached: HTTPSConnectionPool(host='search.rcsb.org' ...)` / `A PDB search needs an internet connection.` This is clear.
- AlphaFold offline: `Error checking AlphaFold availability for P69905: <urlopen error [Errno 61] Connection refused>`, `No AlphaFold structure found for P69905`, `Note: Not all proteins have AlphaFold predictions yet`. This is **misleading**: the real cause is no internet. AlphaFill offline is the same (`No AlphaFill structure found for P69905`).
- `--pdbfile <path>` (source): copies the file into the current folder and loads it like source 4. A missing file prints `File not found: <path>`.

### Scenario 1a-21: Structure Loader options 2 and 4 — SOURCE-ONLY (other scopes)
- `2. Load AMBER topology & coordinate files`: asks for a `.parm7/.prmtop`, then `r` (restart/coordinates) or `t` (trajectory), then the files. Needed to view an MD trajectory in ua.
- `4. Load a ligand file or library (SDF/mol2/SMILES) for docking`: picks a file, lists its molecules (using RDKit), and stores it for 5d Docking.
- I left both to the MD and docking scopes.

### Options reference (1a)
| Prompt | Choice | What it does |
|---|---|---|
| Structure Loader menu | 1 | Open the source menu (below) |
| | 2 | Load AMBER parm7 + rst7/trajectory |
| | 3 | View metadata (blocked until something is loaded) |
| | 4 | Load a ligand file for docking |
| Select source | 1 RCSB / 2 AlphaFold / 3 AlphaFill / 4 Local file / 5 View local files together / 6 Done | |
| RCSB method | 1 PDB ID / 2 title keyword / 3 gene or protein via UniProt / 4 ligand / 5 back | |
| Structure type (after PDB ID) | 1 Asymmetric unit (default) / 2 Biological assembly | 2 downloads RCSB's assembly file `<ID>_assembly<N>.pdb`. If there are several assemblies it asks `Select assembly number`; with none it falls back to the asymmetric unit |
| Search keyword | free text | Searches entry titles |
| Filter by organism | free text or Enter | Restricts by organism |
| Maximum number of structures to list | 1–10000, default 100 | How many results to fetch |
| Results `Selection:` | `3`, `1,3,5`, `1-3`, `all`, `next`, `prev`, `filter`, `back` | Numbers refer to the current page |
| Refine Search | 1 max resolution (Å) / 2 method: 1 X-ray, 2 NMR, 3 Cryo-EM / 3 contains metal (element or blank = any) / 4 clear (only says "start a new search") / 5 back | |
| Download <ID>? / Download all N? | y/n (y) | Confirms a download from search results |
| UniProt search | gene, protein name or accession; organism optional | |
| View detailed UniProt metadata? | y/n (y) | |
| Select a PDB structure to download? | y/n (y) | Opens the UniProt-linked PDB list: number, `next`, `back` |
| Ligand search | code or name, then optional "related ligands" code, then the codes to search (comma-separated) or `back`, then max results | |
| AlphaFold / AlphaFill | 1 UniProt ID / 2 gene name (+ optional organism, then pick 1..N or `b`) / 3 back | |
| View AlphaFold/AlphaFill metadata? | y/n (y) | |
| Remove hydrogen atoms? | y/n (y) | y saves `<name>_noH.pdb`. The answer is reused for the rest of the session |
| File browser | `N`, `..`, full path or `./name`, `q` | `q` leads to the "Enter path to PDB file" fallback. Only `.pdb` files are listed |
| Add another file to view? | y/n (default y until 2 files are chosen) | Source 5 only |
| View metadata for | 1 RCSB / 2 Local / 3 AlphaFold / 4 AlphaFill | Only the loaded kinds are listed |
| Metadata viewer | 1–15 (see 1a-18) | 13 exports JSON or text |

---

## Feature: 1b. Biological Assembly Generator

**What it is for (plain English):** A crystal structure file often holds only part of the real molecule (the "asymmetric unit"), or a different number of copies. This tool rebuilds the working ("biological") form, for example the full hemoglobin tetramer, using the symmetry instructions stored in the PDB file (REMARK 350).

**Needs before you start:** A loaded structure (1a) whose file still has its REMARK 350 header lines, i.e. the original RCSB PDB file. A `_noH` file, an AlphaFold model or a file that ProPrep re-wrote has no REMARK 350.

**Where to find it:** Main menu `1b` > `1. Generate biological assembly from asymmetric unit`. Before anything is loaded the main menu shows `1b. ○ ... ⚠ Needs a loaded structure`.

### Scenario 1b-1: 1HHO asymmetric unit (2 chains) to tetramer (4 chains) — RUN
- Starting point: new project `asm`; loaded 1HHO with 1a > RCSB > PDB ID > asymmetric unit. No hydrogens, so no `_noH` file.
- Steps:
  1. Main menu `1b`, then `1`.
  2. Screen: `Using RCSB PDB PDB file: .../1HHO.pdb`, `Input structure: 1HHO.pdb`, `Source: rcsb_pdb_file`, `Found 1 biological assembly definition(s) in REMARK 350`, `1. TETRAMERIC (author determined) (software determined) - 2 transform(s), chains: A, B`, then the options `1-N Select assembly from REMARK 350`, `m Enter manual transformation matrix`, `c Cancel`.
  3. `Select option (1):` -> `1`
  4. `Include ligands/cofactors (HETATM)? [y/n] (y):` -> Enter (yes)
  5. `Include water molecules? [y/n] (n):` -> Enter (no)
- What the user sees: `Generating biological assembly...`, `Assembly generated successfully!`, `Original: 2 chains, 2,287 atoms`, `Assembly: 4 chains (A, B, C, D), 4,574 atoms`, `Transforms applied: 2`, `Saved: .../1HHO_biological_assembly.pdb`, `Registered in workspace as 'biological_assembly_pdb_file'`.
- Result / files: `runs/s3/asm/asm/1HHO_biological_assembly.pdb`: one model, chains A, B, C, D (1119/1168/1119/1168 atom lines), with 4 hemes (172 HEM atoms), 4 O2 and phosphate, and no water. The header records are kept, so the file still has REMARK 350.
- Time: under a second.
- From then on, this assembly is the structure that later tools use.

### Scenario 1b-2: Running 1b again on the assembly — RUN (trap)
- Choosing 1b again picked `1HHO_biological_assembly.pdb` as input and offered the same `TETRAMERIC ... chains: A, B` assembly. Accepting would apply the symmetry a second time to the already-built tetramer. I cancelled. Warn users to run 1b only once.

### Scenario 1b-3: Invalid assembly number, and cancel — RUN
- `Select option (1):` -> `7` -> `Invalid assembly number: 7`, back to the 1b menu. The prompt has no choice list, so it does not ask again.
- `c` -> `Cancelled.`

### Scenario 1b-4: Manual transformation matrix (`m`) — FAILED (crash)
- Starting point (s3-main): 1b chose `9ZZZ_noH.pdb` automatically (see 1a-19). It reported `No REMARK 350 records found in PDB file.` and offered only `m` or `c`.
- Steps: `Select option (m):` -> `m`; `Chains to transform (comma-separated, or 'all') (all):` -> `A`; the screen printed `Enter transformation matrix (3x3 rotation matrix):`.
- What the user sees:
  ```
  Error executing Show Generate biological assembly from asymmetric unit menu: name 'self' is not defined
  Error executing Show main menu: name 'self' is not defined
  ERROR - main.py:1192 - Unhandled exception: name 'self' is not defined
  Error: name 'self' is not defined
  Run with --debug flag for more details
  [Session recording stopped]
  ```
  **ProPrep quits completely.** Files already written stay on disk, and the session JSON up to that point is kept.
- Cause (from the source): the matrix-entry routine refers to the program object, but it is written outside it. This is a programming bug. Manual matrices cannot be used in 1.23.0.
- Workaround: none inside 1b. Download the assembly from RCSB instead (1a, Structure type 2).

### Scenario 1b-5: Loading another structure after 1b — RUN
- After 1b on 1HHO, I loaded 4HHB in 1a and opened 1b again. It still said `Using Biological Assembly PDB file: .../1HHO_biological_assembly.pdb`. The new download is ignored (see 1a-19).

### Scenario 1b-6: More than 62 chains — SOURCE-ONLY
- If chains × transforms is more than 62, ProPrep warns `Assembly will have N chains, but PDB format only supports 62 unique chain IDs.` and asks `Continue anyway (some chain IDs will be duplicated)? (n)`.

### Options reference (1b)
| Prompt | Choice | Effect |
|---|---|---|
| Select option | `1..N` | Use assembly N from REMARK 350 (default 1 if any exist) |
| | `m` | Type rotation matrix rows and a translation by hand (**crashes ProPrep in 1.23.0**) |
| | `c` | Cancel |
| Chains to transform | `all` or `A,B` | Manual mode only |
| Include ligands/cofactors (HETATM)? | y (default) / n | Keep hemes, ions, ligands in the copies |
| Include water molecules? | n (default) / y | Keep crystal waters |
| Continue anyway (>62 chains)? | n (default) / y | |
| Output | | `<input>_biological_assembly.pdb` next to the input |

---

## Feature: 3a. Homology Searcher (BLAST)

**What it is for (plain English):** It takes the amino-acid sequence of your protein and searches public databases (NCBI BLAST) for related proteins, for example to find a better template or an AlphaFold model of a relative.

**Needs before you start:** An internet connection. The search runs on NCBI's web BLAST server, so local BLAST+ (`blastp`) is not needed. A loaded structure is optional: you can paste a sequence or use a FASTA file. MODELLER, with a licence key, is needed only for option 5.

**Where to find it:** Main menu `3a`.

Menu:
```
══ Homology Searcher Menu ══
  1. ✓ Run BLAST homology search
  2. ○ View BLAST results               ⚠ Need to run BLAST first
  3. ○ Retrieve AlphaFold structures for hits   ⚠ Need to run BLAST first
  4. ○ Export alignments and annotations        ⚠ Need to run BLAST first
  5. ○ Build homology model (MODELLER)  ⚠ MODELLER not installed
→ Suggestion: Start by running a BLAST search (option 1) ...
```

### Scenario 3a-1: blastp of 4HHB chain A against Swiss-Prot through NCBI — RUN (see result below)
- Starting point: s3-asm, with 1HHO assembly and 4HHB loaded.
- Steps:
  1. `3a` > `1`
  2. A table `Available Structures in Workspace` (`1 Biological Assembly .../1HHO_biological_assembly.pdb`, `2 RCSB PDB .../4HHB.pdb`) / `Select structure to use (1):` -> `2` -> `✓ Selected RCSB PDB`, `Using PDB ID 4HHB extracted from filename`
  3. `===== BLAST Interactive Search Tool =====` / `Sequence input methods: 1. Enter sequence directly 2. Load from FASTA file 3. Use sequence from loaded structure (4HHB)` / `Choose input method (1):` -> `3`
  4. `Available chains in structure 4HHB: Chain A: 141 residues, Chain B: 146, Chain C: 141, Chain D: 146` / `Select chain (A):` -> `A` -> `Using sequence from chain A (141 residues)`
  5. `BLAST Program Options: 1. blastn 2. blastp 3. blastx 4. tblastn 5. tblastx` / `Select BLAST program (2):` -> `2`
  6. `Database Options: 1. nr 2. refseq_protein 3. swissprot 4. pdb` with `Note: Options 1 and 3 allow retrieval of AlphaFold structures for homologs.` / `Select database (1):` -> `3`
  7. `Enter E-value threshold (10.0):` -> Enter
  8. `Enter maximum number of hits to return (50):` -> `20`
  9. `Scoring Matrix Options: 1. BLOSUM62 (default) 2. BLOSUM45 3. BLOSUM80 4. BLOSUM90`, each with a "best for" line / `Select scoring matrix (1):` -> `1`
  10. Summary: `Program: blastp, Database: swissprot, Sequence length: 141, E-value threshold: 10.0, Maximum hits: 20` / `Run BLAST with these parameters? [y/n] (y):` -> `y`
- What the user sees: `This may take several minutes depending on sequence length and server load.` and a spinner, `Submitting BLAST search to NCBI server...`. After 10 minutes a Biopython warning appeared: `BLAST request G3XE7CNX016 is taking longer than 10 minutes, consider re-issuing it`. Checking NCBI directly about 25 minutes in, the request was still `Status=WAITING` (NCBI's queue).
- Result (completed later, recorded by the coordinator): the search finished after 13741 s (about 3 hours 50 minutes; ProPrep printed "BLAST search completed in 13741s"). ProPrep showed a table of 20 hits (rank, description, score, E-value, identity, length): all hemoglobin alpha chains, top hit `sp|P69905.2` human hemoglobin alpha, E-value 2.28e-99, 100% identity; hit 20 still 95.7%. It then asked `View detailed alignment for specific hits? [y/n] (y)` and printed hit #1 in full: score 726 bits, `Identities: 141/141 (100.0%)`, query 1-141 vs subject 2-142, and the alignment in 60-letter blocks (Query / match line / Sbjct).
- Option 3 "Retrieve AlphaFold structures for hits": shows the hits with UniProt IDs, then `Selection Methods: 1. Top N hits 2. By E-value threshold 3. By identity % threshold 4. Manual selection (enter ranks) 5. All hits` -> `1`, N = `2` -> progress bar, `✓ Successfully retrieved 2 AlphaFold structures` (files `homologs/AF-P69905-F1-model_v6.cif`, `homologs/AF-P01923-F1-model_v6.cif`). Then `Would you like to select a homolog for downstream processing? [y/n] (y)` -> Enter, a table of the retrieved models, `Select homolog (enter number) (1):` -> `2` -> `✓ Selected P01923 for downstream processing`. Menu now shows `3. ● Retrieve AlphaFold structures for hits [2 retrieved] ✓`.
- Option 4 "Export alignments and annotations": `Enter output filename for BLAST results (blast_export.xml):` -> Enter -> `BLAST results exported to: blast_export.xml` (NCBI BLAST XML), then `Export summary table as well? [y/n] (y)`. Note this prompt only accepts y/n (typing `x` gives "please enter y or n").
- Files: `runs/s3/asm/asm/blast_export.xml`, `runs/s3/asm/asm/homologs/*.cif`.
- Problems:
  - There is no progress information, no time limit and no cancel option. Ctrl-C is the only way out, and NCBI's public server can be slow. Tell users to expect long waits and to choose a small database (`swissprot` or `pdb`, not `nr`).
  - The warning showed a Biopython file from another conda environment (`envs/etanalyze_env/...`), so this computer's Python path mixes environments. That is not a ProPrep problem, but a support person might find it useful.
- Session log: runs/_sessions/s3-asm/transcript.txt

### Scenario 3a-2: Options 2–4 (view, AlphaFold for hits, export) — RUN (recorded in 3a-1 above)
### Scenario 3a-3: Build homology model (option 5) — SOURCE-ONLY / blocked
- The menu shows `⚠ MODELLER not installed` (MODELLER is installed, but its licence key is invalid, so its import fails). It cannot run without a MODELLER academic licence key (`KEY_MODELLER`). At start-up ProPrep prints how to install it and set the key.

### Options reference (3a)
| Prompt | Choices | Notes |
|---|---|---|
| Select structure to use | numbered workspace structures | Only when more than one is loaded |
| Sequence input | 1 type/paste sequence / 2 FASTA file / 3 chain of the loaded structure | |
| Select chain | chain letter | |
| BLAST program | 1 blastn, 2 blastp (default), 3 blastx, 4 tblastn, 5 tblastx | Protein work: blastp |
| Database | protein: 1 nr (default), 2 refseq_protein, 3 swissprot, 4 pdb. Nucleotide programs get a nucleotide list (…, 2 refseq_rna, 3 refseq_genomic) | AlphaFold retrieval for hits works only with nr or swissprot |
| E-value threshold | number, default 10.0 | Invalid input -> "Invalid E-value. Using default." |
| Maximum hits | integer, default 50 | |
| Scoring matrix (protein programs) | 1 BLOSUM62, 2 BLOSUM45, 3 BLOSUM80, 4 BLOSUM90 | |
| MegaBLAST? (blastn only) | y (default) / n | |
| Run BLAST with these parameters? | y/n | |

---

## Feature: 3b. Structure Aligner

**What it is for (plain English):** It lays two or more structures on top of each other as closely as possible and reports how different they are (RMSD in Å; lower means more alike). Use it, for example, to compare a crystal structure with an NMR structure or with an AlphaFold prediction, or to put several structures into the same frame.

**Needs before you start:** At least one loaded structure (the menu item is blocked otherwise). You add the second and further structures inside the tool. Structure-based (fold) alignment needs `tmtools`, which is installed.

**Where to find it:** Main menu `3b`. Menu: `1. Align multiple PDB structures`, `2. Align structures using redox site residues` (needs 2a Redox Site Detector first), `3. View alignment results`.

### Scenario 3b-1: Sequence-based alignment, 1HRC (crystal) vs 1AKK (NMR), backbone — RUN
- Starting point: project `align`; loaded 1HRC, then 1AKK (answered `n` to hydrogen removal), both from RCSB.
- Steps:
  1. `3b` > `1`. Panel `Structure Alignment Tool — Progressive superimposition with cumulative residue addition`.
  2. Step 1: `✓ Selected RCSB PDB: .../1AKK.pdb`, `✓ Added 1 workspace structure(s)`. Only the latest RCSB structure is added automatically. The screen lists `Available local structure files: • 1HRC.pdb` and the options `1. Browse for local structure file (PDB/CIF) 2. Download PDB from database 3. Done`.
  3. `Select option (1):` -> `1`; browser -> `3` (1HRC.pdb) -> `✓ Added`. `Select option (3):` -> `3` (Done).
  4. For each file: `Structure contains 1 model(s), 1 chain(s), 1745 atoms` / `Use entire structure from ...1AKK.pdb? [y/n] (y):` -> Enter (and the same for 1HRC, 993 atoms).
  5. Step 2: advice (choose an experimental reference if you want its ligands transferred, AlphaFold for full coverage) / `Enter reference structure index (1):` -> `2` (1HRC). Note: `The transformation matrix will be calculated using the selected alignment residues, but will be applied to ALL atoms in each target structure`.
  6. Step 3: `Use all common residues for alignment? [y/n] (y):` -> `y` -> `Extracted 106 residue(s) from reference structure` (104 amino acids + ACE cap + heme).
  7. Step 4: `[1] Sequence-based (automatic) [2] Structure-based (TM-align) [3] Manual specification` / `Select method (1):` -> `1`.
  8. Step 5: `[1] Full residue (all atoms) [2] Backbone only (N, CA, C, O) [3] Alpha carbon only (CA)` / `Enter mode (2):` -> Enter (backbone).
  9. The screen shows `Note: 2 non-protein residues excluded from sequence alignment`, `Chain A reference sequence: 104 residues → 1AKK chain A: 104 residues, Alignment: 104/104 identical (100.0%)`.
  10. HETATM mapping explanation / `Interactively map non-protein residues using spatial proximity? [y/n] (n):` -> `y`. For `A:ACE:0`: `No matching ACE residues found`, options `0 Skip`, `1 Add new ACE at aligned position` -> `0`. For `A:HEC:105`: `Found 1 candidate(s)`, options `0 Skip`, `1 A:HEC:105 0.30 Å`, `2 Add new HEC at aligned position` -> `1` -> `✓ Mapped to A:HEC:105`.
  11. Step 6: `Save intermediate alignment steps for debugging? [y/n] (n):` -> `n`.
- What the user sees: a progress bar with the running RMSD; `⚠ Warning: No atoms found for A:HEC:105 in reference, skipping` (in backbone mode, heme has no N/CA/C/O atoms, so mapping it was pointless); an ASCII chart `RMSD Progress: 1AKK` (RMSD as residues are added one by one); `Final RMSD: 1.679 Å`; `Saving Aligned Structures... ✓ 1AKK_aligned.pdb ✓ 1HRC_aligned.pdb`; `Reference aligned structure saved to workspace`. Then a table `Select Target Structure for Downstream Processing` / `Select target structure (or 'n' to skip) (1):` -> `1` -> `✓ 1AKK_aligned.pdb is now available for downstream modules`, `✓ Structure alignment complete!`.
- `3. View alignment results` prints a table, one row per step: residues 1, 1-2, ... 1-105 with the RMSD after each.
- Result / files: `runs/s3/align/align/1AKK_aligned.pdb` and `1HRC_aligned.pdb` (both in the common frame). A crystal vs NMR backbone RMSD of about 1.7 Å is reasonable.
- Time: a few seconds.

### Scenario 3b-2: Structure-based (TM-align), 4HHB chain A vs AlphaFold P69905 — RUN (with a workaround)
- Steps:
  1. Loaded 4HHB (RCSB) and AlphaFold P69905 in 1a. Then `3b` > `1`.
  2. `Found 2 aligned structure(s) from previous alignment ... Include these aligned structures for further alignment? [y/n] (y):` -> `n`. After a "no", ProPrep does **not** offer the other loaded structures. `Current structures: 0`.
  3. Browse -> `6` (4HHB.pdb). Browse again: the AlphaFold download `AF-P69905-F1-model_v6.cif` is **not listed**, and typing its path gives `AF-P69905-F1-model_v6.cif is not a PDB file.` This is despite the menu text "PDB/CIF". Workaround: I downloaded the PDB-format AlphaFold file (`AF-P69905-F1-model_v6.pdb`) into the project folder with curl, then picked it (`7`).
  4. `Done` -> 4HHB: `Use entire structure? (y)` -> `n` -> `Model has 4 chains 1. Chain A (198 residues) ...` (counts include waters) and a warning to select all chains wanted in the output / `Select chains (all):` -> `1` -> `✓ Added selected portion (1168 atoms)`. AF: `Use entire structure? (y)` -> Enter.
  5. Reference -> `2` (AlphaFold). `Use all common residues? (y)` -> `y` -> `Extracted 142 residue(s)`. Strategy -> `2` (TM-align).
- What the user sees: `Performing TM-align structure-based alignment...`, `Aligning: 4HHB.pdb`, `Aligning all chains combined`, `Reference: 142 CA atoms`, `Target: 141 CA atoms`, `TM-score (normalized by target): 0.9943`, `TM-score (normalized by reference): 0.9874`, `RMSD: 0.335 Å`. Then `Saving Aligned Structures...` (printed twice), `✓ 4HHB_aligned.pdb`, `✓ AF-P69905-F1-model_v6_aligned.pdb`, and the target-selection table (`n` = skip).
- There is no atom-mode step and no progressive chart. TM-score runs from 0 to 1, and above about 0.5 means the same fold. 0.99 means practically identical.
- If either structure has several chains, ProPrep first asks `[1] Align all chains combined [2] Align on specific chain` (source).
- Result / files: `4HHB_aligned.pdb`, `AF-P69905-F1-model_v6_aligned.pdb`.
- **Bug:** `3. View alignment results` after a TM-align run gives `Error executing Structure Aligner._view_alignment_results: 'step'` and a Python traceback. ProPrep survives and returns to the menu.

### Scenario 3b-3: Manual mapping, residue range, CA only, and the in-tool download — RUN
- Steps: `3b` > `1`; previous aligned? `n`; option `2. Download PDB from database` -> `Enter PDB ID to download (e.g., 1CLL) ():` -> `1hrc` -> `File 1HRC.pdb already exists, using it`; `2` -> `ZZ99` -> `✗ PDB ZZ99 not found (may be a large structure requiring mmCIF format)`; `2` -> `1AKK`; `3` Done; entire structures; reference `1` (1HRC).
- `Use all common residues? (y)` -> `n` -> `Residue Selection Options: [1] Specify residue range (e.g., A:50-150) [2] Manual residue-by-residue entry` -> `1`; `Enter residue range(s) ():` -> `A:1-40,X:5` -> `✓ Added 40 residues from A:1-40`, `✗ Error parsing 'X:5': Invalid format. Use chain:start-end` (the good part is kept, and there is no retry).
- Strategy `3` (Manual: same residue numbers in both), atom mode `3` (CA), intermediates `n`.
- Result: `Final RMSD: 1.991 Å` on residues 1–40, CA only. Files `1HRC_aligned_2.pdb` and `1AKK_aligned_2.pdb`. A second run adds `_2` to the names instead of overwriting. Target selection `n` -> `✓ Structure alignment complete!`.
- Manual residue-by-residue entry (option 2, source): `chain:resname:resid`, or `chain:resname:resid <- idx:chain:resname:resid, ...` to pair different numbering; an empty line finishes.

### Scenario 3b-4: Fewer than two structures — RUN
- `3b` > `1` > previous aligned `n` > `3` (Done) with 0 structures -> `At least 2 structures required for alignment` / `Error: At least 2 structures required for alignment`, back to the menu.

### Scenario 3b-5: Align on redox sites without detected sites — RUN
- After any alignment, the menu shows option 2 as available (●). Choosing it gives `No redox sites detected in workspace.` / `Please run redox site detection first (PDB Filter module).` The text names the wrong module: redox detection is **2a Redox Site Detector**. Before an alignment, the menu correctly shows `⚠ Detect redox sites first`.

### Options reference (3b)
| Prompt | Choices | Effect |
|---|---|---|
| Include these aligned structures for further alignment? | y (default) / n | y = start from the previous run's outputs. n = start empty; loaded structures are then not offered |
| Add Structures | 1 browse (only `.pdb` listed) / 2 download by PDB ID / 3 done | |
| Use entire structure from X? | y / n | n -> choose model (if several) and chains: `1`, `1,2`, `1-3`, `all` |
| Reference structure index | 1..N | Others are moved onto it |
| Use all common residues? | y / n | n -> 1 range `A:50-150,B:20-80` or 2 residue-by-residue |
| Mapping strategy | 1 sequence-based (homologues) / 2 TM-align (same fold, any sequence) / 3 manual (identical numbering) | |
| Atom mode (1 and 3 only) | 1 all atoms / 2 backbone N,CA,C,O (default) / 3 CA | |
| Map non-protein residues interactively? | n (default) / y | For each ligand/ion in the reference: skip, pick a nearby match, or add a copy at the aligned position |
| TM-align chain mode (multi-chain only) | 1 all chains combined / 2 one chain | |
| Save intermediate steps? | n / y | |
| Select target structure for downstream | number or `n` | Which aligned file later tools use |
| Outputs | | `<name>_aligned.pdb`, then `_aligned_2.pdb` and so on for later runs |

---

## Feature: ua. Structure Viewer

**What it is for (plain English):** It shows your structures in 3D in a web browser, with ProPrep's findings (redox sites, protonation states, filter choices, alignments) marked on them.

**Needs before you start:** A loaded structure, or a loaded topology plus trajectory for option 2. A web browser on the same computer. ProPrep runs a small web server on your machine (localhost only).

**Where to find it:** Main menu `ua` (Full Menu). The viewer also opens by itself from 1a's "View local files together", and whenever you open structure metadata.

Menu:
```
══ Structure Viewer Menu ══
  1. ✓ Launch interactive viewer
  2. ○ View the loaded MD trajectory   ⚠ Load a topology + trajectory in the Structure Loader first
  3. ✓ Show available annotations
  4. ✓ Save the scene shown in the open viewer
  5. ✓ Load a saved scene
  6. ✓ Launch viewer with electron density (X-ray structures from the PDB)
```

### Scenario ua-1: Launch the viewer with two aligned structures — RUN
- Steps: `ua` > `1` -> panel `Interactive Structure Viewer — Structure Selection ... Structures will be displayed together (superimposed if aligned).` and the table `Available Structures in Workspace` (Aligned Target 1AKK_aligned.pdb, Aligned Reference 1HRC_aligned_2.pdb, RCSB 4HHB.pdb, AlphaFold AF-P69905...) / `Select structure(s) (1):` (number, range `1-3`, list `1,3,5` or `all`) -> `1,2`.
- What the user sees: `✓ Selected 2 structures`, `✓ Structure viewer launched at http://localhost:8765/viewer (2 structures)`. That one line is the only instruction. ProPrep also tries to open the page in your default web browser by itself. When no graphical browser is available (e.g. over SSH) it does not, and the source has a hint for an `ssh -L` tunnel.
- Check with curl: `http://localhost:8765/viewer` -> HTTP 200, 163,588 bytes, title "ProPrep Interactive Structure Viewer". `/config` returned JSON with `"name": "1AKK_aligned"`, representations "Protein" (cartoon, blue) and "Ligands" (ball+stick). `/version` returned `{"version": 1, "instance": ...}`.
- Problem: the "Aligned Target" (1AKK_aligned.pdb, from the first run) and the "Aligned Reference" (1HRC_aligned_2.pdb, from the third run) came from **different** alignment runs, because I skipped target selection in later runs. So the pair on offer was not a matching pair. Advise users to pick the target at the end of each alignment, or to use fresh projects.

### Scenario ua-2: Show available annotations — RUN
- `3` -> `No annotations found in workspace`, followed by the list of possible annotation types (Redox Sites, Titratable Residues, Filter Selections, Bond Classifications, Aligned Structures). It said "none" even though aligned structures were present.

### Scenario ua-3: Save and load a scene — RUN
- `4` -> `Scene name (scene_20261121_130439):` -> Enter -> `✓ Scene saved: .../scene_20261121_130439.scene.json (2 structure(s), 10 representation(s))`.
- `5` -> `Saved scenes: 1. scene_20261121_130439.scene.json 2 structure(s), 2026-11-21T13:04:49` / `Select scene (1):` -> `1` -> `✓ Scene applied to the open viewer`.
- File: `runs/s3/align/align/scene_20261121_130439.scene.json`.

### Scenario ua-4: Electron density for 4HHB — RUN (no map available)
- `6` -> an explanation panel (2mFo-DFc = the density the model was built into; mFo-DFc = the difference map, "Trust the difference map more"; maps only for X-ray entries in their deposited frame, not for superposed structures or NMR) -> the structure table -> `3` (4HHB) -> `Density margin around the model (Angstroms) (min: 0.0, default: 5.0) (5.0):` -> Enter.
- What the user sees: `4HHB.pdb: fetching the maps of 4HHB from PDBe...`, `PDBe has no maps for 4HHB. They exist for X-ray entries deposited with their structure factors, not for NMR or cryo-EM entries`, `No electron density could be shown; launching the viewer without it`, `✓ Structure viewer launched at http://localhost:8766/viewer (1 structure)`.
- This opened a **second** port (8766) in the same ProPrep. 8765 stayed open too.

### Scenario ua-4b: Electron density for 3NIR (0.48 Å crambin) — RUN
- Session `s3-misc`, with only 3NIR loaded (hydrogens kept). `ua` > `6` -> the explanation panel. With one structure loaded, it was selected automatically (`✓ Selected RCSB PDB: .../3NIR.pdb`) -> `Density margin ... (5.0):` -> Enter.
- What the user sees: `3NIR.pdb: fetching the maps of 3NIR from PDBe...`, `✓ 3NIR.pdb: density of 3NIR (2mFo-DFc sigma 1.033, mFo-DFc difference sigma 0.263). Its atoms sit at +10.17 sigma of the 2mFo-DFc map on average`, `✓ Structure viewer launched at http://localhost:8765/viewer (1 structure)`.
- Result / files: `3nir_2fofc.ccp4` and `3nir_fofc.ccp4` (5.5 MB each, as downloaded from PDBe) and `3nir_2fofc_around_model.ccp4` / `3nir_fofc_around_model.ccp4` (22.9 MB each, the maps cut around the model with a 5 Å margin). For this ultra-high-resolution entry the cut maps are **bigger** than the originals, so expect large files. They take about 10 s to fetch.
- If the atoms are not in the density (e.g. a structure moved by 3b), the source prints `no density. Its atoms sit at ... sigma ... a model in its own density sits near +3 sigma` and shows no map.

### Scenario ua-5: Stopping the viewer — RUN
- There is **no "stop viewer" menu item**. The viewer keeps running until ProPrep exits. `x` -> `Are you sure you want to exit? [y/n] (n):` -> `y` -> the citation panel. Afterwards, `curl http://localhost:8765/viewer` and `:8766` both failed (connection refused): the viewer stopped with ProPrep.
- Note: just pressing Enter at the exit question means "no" and returns to the main menu.

### Scenario ua-6: View MD trajectory (option 2) — SOURCE-ONLY
- Blocked until 1a option 2 has loaded a topology and trajectory (MD scope).

### Options reference (ua)
| Option | Effect |
|---|---|
| 1 Launch | Pick structures from the workspace table (`1`, `1-3`, `1,3`, `all`), serve them at `http://localhost:<port>/viewer` (port 8765, or the next free one), and open the browser if one is available |
| 2 Trajectory | Show a loaded parm7 + trajectory |
| 3 Annotations | List the workspace annotations available to the viewer |
| 4 Save scene | Saves `<name>.scene.json` (structures, representations) in the project folder |
| 5 Load scene | Re-applies a saved scene to the open viewer |
| 6 Density | As 1, plus PDBe 2mFo-DFc / mFo-DFc maps; asks for a margin in Å (default 5.0) |

---

## Notes for the manual writer

**Jargon to explain:** asymmetric unit vs biological assembly; REMARK 350; mmCIF vs PDB format (ProPrep needs PDB format for most things); pLDDT (AlphaFold confidence, 0–100); AlphaFill "transplants" (ligands copied from similar experimental structures); RSCC (ligand fit to density); E-value; BLOSUM matrices; RMSD; TM-score; backbone vs CA vs all-atom alignment; reference vs target; 2mFo-DFc / mFo-DFc; UniProt accession.

**Recommended beginner path:**
1. Start a fresh project for each protein you prepare.
2. 1a > 1 > 1 > 1 > type the PDB ID > `1` (asymmetric unit), or `2` if you want the functional oligomer. Say `y` to hydrogen removal.
3. If the entry's asymmetric unit is not the working molecule, run **1b once**, choosing the REMARK 350 assembly.
4. To compare with AlphaFold, use 3b with TM-align, after getting the AlphaFold model as a `.pdb` file (see the bugs below).

**Common mistakes / surprises:**
- Loading a new structure does not make it the one later tools use if ProPrep has already made a file of its own (`_noH`, an assembly, aligned files). It silently keeps using the older file. Check the `Using ... file:` line printed at the top of each tool. Use a new project, or `w` > `5. Reset workspace`, to change structures.
- Within one kind, the last download wins (the last RCSB entry, the last AlphaFold model).
- Batch downloads from search results skip the hydrogen question and the summaries.
- On the second page of search results the numbers start at 1 again.
- The file browsers list only `.pdb`. AlphaFold and AlphaFill downloads are `.cif`, and so cannot be re-picked in the browsers.
- With no internet, the AlphaFold/AlphaFill message wrongly suggests the protein has no prediction.
- Prompts that show no choice list (1b "Select option") do not ask again after a wrong answer. They leave the tool.
- BLAST through NCBI can take a very long time, with no cancel option except Ctrl-C.

**Things that look like bugs (describe, don't fix):**
1. 1b manual matrix (`m`) crashes the whole program: `name 'self' is not defined`.
2. AlphaFill always reports "No transplanted ligands/cofactors found", even when the file has heme and others (P69905: 22 hits, 438 ligand atoms).
3. The AlphaFold metadata "Length" and the confidence distribution count atoms, not residues (P69905: "1077 residues" for 142 aa). After a gene-name search, the UniProt ID, protein and organism show as N/A/Unknown.
4. 3b "View alignment results" after a TM-align run crashes the action (`'step'` KeyError; the program survives).
5. 3b redox-site option points to the "PDB Filter module" instead of 2a, and shows as available after any alignment.
6. 3b "Browse for local structure file (PDB/CIF)" rejects `.cif`. Declining the previous aligned files also hides the loaded workspace structures.
7. A `--pdbid` download failure is followed by the log line "Downloaded and loaded PDB".
8. Cosmetic: "press  to return to the main menu" (missing `[m]`), "1 EMARK" in revision history, "Saving Aligned Structures..." printed twice, search `back` saying "Structure Loader menu" while it returns to the source menu, "Clear all filters" not clearing.
9. The RCSB assembly download (1a type 2) gives a two-MODEL file with repeated chain IDs and waters, whereas 1b gives chains A–D in one model. The two routes produce different files for the same assembly.
10. ua "Show available annotations" says "none" even with aligned structures present. ua option 6 opens a second server port.
