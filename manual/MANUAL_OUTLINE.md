# ProPrep User Manual — proposed structure (for ProPrep 1.23.0)

Audience: scientists who are not programmers. Plain English, short sentences,
every technical term explained the first time it appears and again in the glossary.
Source material: `manual/scenarios/s0…s8` (recorded runs) and `manual/runs/` (evidence).

---

## Front matter
- Title page, version box ("Written for ProPrep 1.23.0, November 2026")
- **How to read this manual**: what you type vs what you see; the four box types (below)
- **Quick start on one page**: install → start → load 1UBQ → build a system → run a short simulation

## Part I — Before you begin
1. **What ProPrep is (and is not)** — prepares protein structures for Amber simulations; does the
   tedious work, leaves the scientific decisions to you; text based; records every answer.
2. **Requirements** — computer, disk space, internet, browser; table *"What you need for what"*
   (everyone: ProPrep + AmberTools; optional: MODELLER key, Gaussian, ORCA, pmemd/GPU). *(s0)*
3. **Installation** — Option A one-file installer (Mac / Windows via WSL / Linux, recommended);
   Option B conda script; adding the MODELLER key; checking it works; updating; removing;
   problems, including "another Amber installation is active" (AMBERHOME). *(s0, INSTALL.md)*

## Part II — Finding your way around
4. **Starting and stopping ProPrep** — in the browser (`proprep-web`) or in the terminal (`proprep`);
   the first-launch questions; naming a project; Workflow Mode vs Full Menu; exiting;
   coming back to a project later. *(s1, s2)*
5. **Menus and prompts** — menu keys like `1a`; what ✓ ○ ⚠ mean (✓ = available, not "done");
   defaults in brackets; y/n questions; going back (`b`, `m`); `undo`; Help; Preferences. *(s1)*
6. **Your project folder** — what is created where; the common files and what they are;
   *which file each tool picks up* (why step order matters, orient last). *(s3, s4)*

## Part III — What ProPrep can do (the scope)
7. **The big picture** — one diagram of the stages:
   Get a structure → Clean & check → Protonate → Special sites (redox/metals/ligands) →
   Force fields → Build the system → Simulate → Analyse; plus Compare, View, Dock, Membrane, QM/MM.
   **"Which chapters do I need?"** table by goal:
   | Your goal | Chapters |
   |---|---|
   | Simulate an ordinary protein in water | 8, 9, 10, 11, 14, 15, 16 |
   | Protein with heme / metal / redox centre | + 12, 13 |
   | Protein with a bound drug or ligand | + 13, 18 |
   | Membrane protein | + 17 |
   | QM/MM calculation | + 12, 19 |
   | Same preparation for many proteins | + 22 |

## Part IV — Step by step (how to get each part working)
Every chapter uses the same template:
**What it does · When you need it · Before you start · Step by step (the recommended path) ·
What you get (files) · Other choices · If something goes wrong.**

8. **Tutorial: ubiquitin from download to analysed simulation** (≈30 min, all verified end to end) *(s7)*
9. **Getting a structure** — by PDB ID, by searching RCSB, AlphaFold / AlphaFill, your own file;
   biological assemblies (1a, 1b) *(s3)*
10. **Cleaning and repairing** — removing waters/ions/ligands (4a), fixing missing parts (4c, needs
    MODELLER key), mutating residues (4b) *(s4)*
11. **Protonation states and orientation** — PROPKA at your pH, accepting/overriding (4e); orienting (4f) *(s4)*
12. **Redox-active sites** — detecting hemes, Fe–S, flavins, copper (2a); preparing them with the right
    names and oxidation states (4d); worked example cytochrome c *(s5)*
13. **Force fields and new parameters** — checking what is already covered (2b); making parameters for
    ligands, modified residues and metals (2c, needs Gaussian); the personal parameter library *(s6)*
14. **Building the simulation system (topology)** — force field, water, box, salt; disulfides; hemes (5a) *(s7)*
15. **Running a simulation** — protocols, shortening them, choosing the engine (sander),
    restraints, running locally or on a cluster, monitoring (5b) *(s7)*
16. **Analysing a simulation** — RMSD, RMSF, energies, distances, etc. (5b analysis, `--analysis`) *(s7)*
17. **Membrane proteins** — building a lipid bilayer around a protein, or a membrane alone (4g) *(s8)*
18. **Docking a ligand** — receptor and ligand preparation, search box, scoring, reading results (5d) *(s8)*
19. **QM/MM input files** — Gaussian ONIOM and ORCA (5c, advanced, needs Gaussian/ORCA) *(s5)*
20. **Comparing and viewing structures** — BLAST search (3a), alignment (3b), the 3D viewer (ua, `--pdbview`) *(s3)*

## Part V — Saving time
21. **Sessions: ProPrep remembers every answer** — the session file, resuming, the session editor,
    `undo`, replaying as a demonstration *(s2)*
22. **Templates and batch processing** — turning one session into a recipe for many proteins *(s2)*
23. **Command-line shortcuts** — the flags a user might want, in one table *(s1, s2)*

## Part VI — When things go wrong
24. **Troubleshooting by symptom** — "it says MODELLER not found", "topology says ✓ but the file is empty",
    "docking fails with a template error", "the membrane build was discarded", … *(all)*
25. **Known limitations in version 1.23.0** — each with a workaround *(all "looked like bugs" lists)*
26. **Getting help and sending feedback** (warning: Enter = submit) *(s1)*

## Appendices
- A. **Glossary** (≈80 terms, gathered from the "Jargon to explain" lists)
- B. **Menu map** — every key in Full Menu and Workflow Mode
- C. **Which features need which extra software** (MODELLER, Gaussian, ORCA, pmemd, GPU)
- D. **File types** — .pdb .cif .prmtop .rst7 .frcmod .lib .mol2 .pdbqt .com/.gjf .json …
- E. **Amber residue names** — HID/HIE/HIP, ASH, GLH, CYX, CYM, LYN, AS4/GL4 …

---

## Writing conventions
- Instructions use **Full Menu keys** (`4e`), which never change; a margin note gives the Workflow Mode route.
- Box types: **Tip** · **Watch out** (traps and known bugs, with the workaround) ·
  **Needs extra software** (MODELLER / Gaussian / ORCA) · **Behind the scenes** (optional science).
- Screens: short verbatim excerpts in a "screen" box; key screens as figures rendered from the
  recorded terminal output (Rich → SVG → PDF) and headless-Firefox captures of `proprep-web`.
- Worked examples use the proteins actually tested: 1UBQ, 1HRC, 1HSG, 1BL8, 4HHB.
- Bugs go in the manual only as user-facing "Watch out" boxes + Ch. 25; the full technical
  bug list is a separate document for the developers.
