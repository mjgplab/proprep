# User-facing problems — writer w4 (Chapters 9, 10, 11, 20)

Format: **Symptom (what the user sees)** | Why | What to do | Chapter label

## 1a Structure Loader / 1b Biological Assembly Generator
- **A tool says "Using H-Stripped PDB file: …_noH.pdb" (or uses an old assembly) after a new structure was loaded** | Each tool picks its input by a fixed priority; files ProPrep made (_noH, assembly, aligned, filtered, repaired, protonated) outrank a fresh download | Read the "Using … file:" line; use one project per protein; Workspace Options "Reset workspace" (not tested) | ch:structure
- **"PDB XQ99 not found (404). May be a large structure only available in mmCIF format."** | The code does not exist, or the entry is mmCIF-only (>99,999 atoms) | Check the code on RCSB; very large entries cannot be loaded | ch:structure
- **"URL can't contain control characters"** | Text with spaces typed as a PDB ID | Type a valid 4-character code | ch:structure
- **"No AlphaFold structure found … Not all proteins have AlphaFold predictions yet" when offline** | Network failure reported as "no prediction" | Check the internet connection | ch:structure
- **AlphaFold metadata says "Length: 1077 residues" for a 142-residue protein; confidence percentages per atom** | Atoms counted as residues | Ignore the count; read pLDDT | ch:structure
- **After AlphaFold gene-name search, metadata shows "UniProt ID: N/A", "Organism: Unknown"** | Metadata not stored for that route | Use retrieval by UniProt ID | ch:structure
- **AlphaFill says "No transplanted ligands/cofactors found" although the file has heme etc.** | Wrong field of AlphaFill's reply read | Do not believe the message; inspect the file in the viewer | ch:structure
- **AlphaFold/AlphaFill .cif files not listed / "X.cif is not a PDB file." in file browsers (1a, 3b)** | Browsers accept only .pdb | Download the PDB-format model from the AlphaFold DB website | ch:structure, ch:compare
- **A .cif loaded via the "Enter path to PDB file" fallback: "BioPython parsing skipped"** | Only path stored, not parsed | Convert to PDB format | ch:structure
- **Bare file name in browser: "Enter a single number, '..', or 'q'."** | Browser needs number, full path or ./name | Type the full path | ch:structure
- **Batch download from search results: no summary, no hydrogen question; only last entry used** | Batch route skips these steps | Download the entry to prepare on its own | ch:structure
- **Hydrogen-removal answer reused silently "(using previous selection)"** | Answer stored for the session | Be aware; new project to answer differently | ch:structure
- **3NIR_noH.pdb still contains 18 H atoms (altloc A of Pro22)** | Hydrogen stripping misses some altloc atoms | Check the file; tLEaP may complain | ch:structure
- **Search "Clear all filters" only says "To clear filters, start a new search"; "back" says Structure Loader menu but returns to source menu** | Cosmetic | Start a new search | ch:structure
- **Page 2 of search results numbered from 1 again** | Numbers are per page | Use the number shown on the current page | ch:structure
- **`--pdbid` failure followed by "Downloaded and loaded PDB: 1CRN"** | Contradictory log line | Trust "Failed to download"; check connection | ch:structure
- **RCSB assembly download (type 2) is a two-MODEL file with repeated chain IDs and waters; no resolution line** | RCSB's own assembly file format | Prefer 1b for simulations | ch:structure
- **1b offers the same assembly again after it has been built** | 1b now uses the assembly file as input | Run 1b once; cancel with c | ch:structure
- **1b "Invalid assembly number: 7" returns to menu without re-asking** | No retry | Choose 1b again | ch:structure
- **ProPrep quits after 1b "m" (manual matrix): "name 'self' is not defined"** | Programming error in 1.23.0 | Do not use m; download the assembly with structure type 2 | ch:structure
- **1b: "No REMARK 350 records found in PDB file."** | 1b is reading a header-less file (e.g. _noH) | New project; download assembly directly | ch:structure
- **Metadata Summary revision history reads "1 EMARK"; menu suggestion "press  to return" (missing [m])** | Cosmetic | None needed | ch:structure

## 4a PDB Filter
- **"Error executing PDB Filter.filter_pdb_structure: invalid literal for int() with base 10: '1-50'" or "list index out of range"; all filter answers lost** | Ranges and out-of-range numbers not handled at "Select residues" | Type single table numbers separated by commas; restart the filter | ch:cleaning
- **Water metrics: "Received '2,4', which is not one of the available options" though prompt says comma-separated** | Only one choice accepted | Type one number | ch:cleaning
- **Water selection: "Invalid residue numbers: [306, 308, 310]" for a range** | Every number in a range must be a water | List the residue numbers | ch:cleaning
- **"n" at "apply these filters?" returns silently, nothing saved** | No message | Run the filter again | ch:cleaning
- **filtered_structure.pdb overwritten without asking** | No overwrite check | Copy the file first if you need it | ch:cleaning
- **Filtered file still contains alternate locations** | Filter does not handle altlocs | See 4c (needs MODELLER) | ch:cleaning
- **"Received '1a', which is not one of the available options" inside a tool** | Main-menu codes not accepted in sub-menus | Type m first | ch:cleaning
- **Main menu shows ✓ for tools already done** | ✓ means "available"; ● inside the tool's menu means done | Look inside the tool menu | ch:cleaning

## 4c Structure Fixer
- **Start-up "MODELLER not found … To install: conda install" although MODELLER is installed** | Message shown when the license key is missing/invalid | Set KEY_MODELLER or ~/.proprep/modeller_key | ch:cleaning
- **"MODELLER off: capping only"; loops, missing side chains, standard mutations, altloc choice not done** | No valid MODELLER key | Get a key; tLEaP adds missing side-chain atoms later | ch:cleaning
- **Repaired file renumbered from 1 (2RH1: 29–1161 → 1–506; 1UBQ +1 with ACE)** | Fixer renumbers on writing | Note the offset for later residue numbers | ch:cleaning, ch:protonation
- **Template comparison reports missing OXT on every residue (75 for 1UBQ, ~440 for 2RH1)** | False positives | Use REMARK 470 | ch:cleaning
- **REMARK 470 detail shows "ASP 29 | A" (chain letter as atom)** | Display error | Read the PDB file's REMARK 470 | ch:cleaning
- **FASTA comparison: 0 missing residues and "No structural issues found!" for 2RH1** | Single-chain FASTA header not parsed; skipped checks reported as clean | Use REMARK 465 or SEQRES | ch:cleaning
- **Sequence gaps: 691 missing residues for 2RH1** | Engineered numbering jumps counted | Use REMARK 465 for constructs | ch:cleaning
- **Long runs of "-" in the sequence view of a chimera; C-terminal SEQRES segment misplaced** | Display spans the numbering range | Not missing residues; ignore | ch:cleaning
- **Repaired file still has A/B altloc atoms (3NIR: 313)** | No-MODELLER branch skips the altloc picker | Needs MODELLER in 1.23.0 | ch:cleaning
- **"No terminal caps selected." and nothing written; pending non-standard mutation silently dropped; summary says "No pending mutations"** | Without MODELLER, mutations only applied with a cap | Add at least one cap | ch:cleaning
- **Template comparison takes 1–2 minutes ("Fetching templates…")** | Downloads templates | Wait, or use REMARK 470 | ch:cleaning

## 4b Amino Acid Mutator
- **Standard mutations (e.g. K48R) never applied, no message at repair time** | Need MODELLER | Get a MODELLER key | ch:cleaning
- **"clear" deletes all pending mutations without confirmation** | No confirmation | Type carefully | ch:cleaning
- **Option 2 blocked "Need to display sequences first"** | Sequence display required first | Run option 1 | ch:cleaning
- **Sequence ruler restarts at 10 on each line** | Display style | The ruler counts from the line start | ch:cleaning
- **After a non-standard mutation, sequence display omits that residue (75 of 76)** | Display error | Ignore; the file is correct | ch:cleaning
- **Lower-case chain accepted although help says chain IDs are case-sensitive** | Inconsistent help | None | ch:cleaning

## 4e Protonation State Analyzer
- **"Invalid pH value, using default of 7.0" / "Invalid threshold, using default of 0.5" without re-asking** | Silent default | Run option 1 again with the right value | ch:protonation
- **"Skipping invalid syntax: XYZ" / "Invalid name 'LYS' - valid choices: GLH, GLU" and ProPrep moves on** | Invalid overrides skipped, no retry | Run option 3 again | ch:protonation
- **View/Export after overrides still show PROPKA's Protonated/Charge values and net charge** | Only names updated in results | The written PDB file is correct | ch:protonation
- **structure_with_md_names.pdb overwritten without asking (also by constant-pH naming)** | No overwrite check | Type a new file name | ch:protonation
- **Export "Select output format (text)" does not list choices; "Received 'zzz'…"; text export named .text** | Choices hidden | Type text or csv | ch:protonation
- **Protonation results change with what was filtered (ASP B25 pKa 9.06 in dimer, 3.8 in chain B alone)** | PROPKA uses the structure as given | Filter carefully before protonating | ch:protonation
- **4e ignores the oriented file** | Priority list excludes oriented files | Orient last | ch:protonation

## 4f Structure Orientator
- **4f keeps using an old oriented file after re-filtering/re-protonating** | Oriented file has top priority in 4f | Orient last; new project | ch:protonation
- **File names grow: …_oriented_oriented_…_oriented.pdb (6 levels)** | Each run takes the last oriented file as input | Note that rotations add up | ch:protonation
- **"Received 'Z', which is not one of the available options" for the middle axis** | Only the remaining two axes allowed | Choose X or Y | ch:protonation
- **"Need at least 3 residues to define a plane." / "Residue B:500 not found." / "Atom B:25:XX not found."** | Input checks; return to menu | Retype | ch:protonation

## 3a Homology Searcher
- **BLAST runs for hours with no progress; "BLAST request … is taking longer than 10 minutes, consider re-issuing it"** | NCBI public queue; no timeout or cancel (our run: 13741 s ≈ 3 h 50 min) | Wait, or Ctrl+C; use swissprot/pdb and few hits | ch:compare
- **BiopythonWarning path shows another conda environment (etanalyze_env)** | Mixed Python path on the test computer, not ProPrep | Support note only | ch:compare
- **"Export summary table as well?" rejects other input ("please enter y or n")** | y/n only | Type y or n | ch:compare
- **Option 5 blocked "MODELLER not installed"** | No MODELLER key | Get a key | ch:compare
- **AlphaFold retrieval for hits only after nr/swissprot searches** | Needs UniProt IDs | Search swissprot | ch:compare

## 3b Structure Aligner
- **"Error executing Structure Aligner._view_alignment_results: 'step'" after TM-align** | View not written for TM-align runs | Ignore; open aligned files in the viewer | ch:compare
- **Browser says "PDB/CIF" but rejects .cif** | .pdb only | Use PDB-format files | ch:compare
- **After declining previous aligned files: "Current structures: 0", loaded structures not offered** | Workspace structures skipped | Add them with the browser | ch:compare
- **"At least 2 structures required for alignment"** | Fewer than two structures | Add a second | ch:compare
- **Option 2 says "Please run redox site detection first (PDB Filter module)"** | Wrong module named; real tool is 2a | Run 2a first | ch:compare
- **Residue range error "Invalid format. Use chain:start-end" skipped without retry** | No re-ask | Run again | ch:compare
- **"Saving Aligned Structures..." printed twice** | Cosmetic | None | ch:compare
- **Selected aligned target outranks later downloads** | Workspace priority | Answer n at target selection if only comparing | ch:compare

## ua Structure Viewer / --pdbview
- **Viewer pair from different alignment runs does not overlap** | Table shows latest target and reference separately | Select target at the end of each alignment; new project | ch:compare
- **"No annotations found in workspace" although aligned structures exist** | Annotation listing incomplete | None | ch:compare
- **Viewer opens on 8766/8767 instead of 8765; electron density opens a second port** | Next free port used | Use the printed address | ch:compare
- **Viewer address stops working** | Viewer stops when ProPrep exits; no stop option | Restart from ua | ch:compare
- **"PDBe has no maps for 4HHB"** | No deposited structure factors / not X-ray | Viewer opens without density | ch:compare
- **Density maps 22.9 MB each for 3NIR** | Maps cut around model with margin | Expect large files | ch:compare
- **`--pdbview` still asks the menu-mode question** | Mode not given | Add --full-menu or --workflow | ch:compare

## Note for the coordinator
- Scenario s3 3a-1 says the BLAST search finished "after roughly 40 minutes", but the transcript prints "BLAST search completed in 13741s" (about 3 h 50 min; started ~12:48, results 16:38). The chapter uses the transcript figure.
- Style file: `\begin{screen}` without the optional argument swallows the first screen line. Always write `\begin{screen}[ProPrep]` (as 00-howto does). Emoji 📁/📄 (file browser) and U+FE0F are missing from the screen font; ✓ ○ ● ⚠ are missing from the sans font used by `\prompt{}` (use `\texttt{✓}` in body text). Screens inside `steps` are about 5 characters narrower than the 86-character limit and clip silently (no overfull warning).
