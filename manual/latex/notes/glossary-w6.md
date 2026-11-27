# Glossary terms — writer w6 (Ch. 17 Membrane, 18 Docking, 21 Sessions, 22 Batch)

Format: `- **term** — plain definition.`

## Chapter 17 — Membrane proteins

- **lipid bilayer** — The two-layer sheet of fat-like molecules that forms a cell membrane. Simulations of membrane proteins build a small patch of one.
- **lipids** — The fat-like molecules a membrane is made of. Each has a water-liking head and oily tails.
- **leaflet** — One of the two layers of a bilayer. An "asymmetric" membrane has different lipids in each leaflet.
- **POPC** — A common neutral lipid of animal cell membranes (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine). The usual first choice for a simple membrane. Cholesterol is called CHL1.
- **orientation (of a membrane protein)** — Turning and shifting the protein so its water-repelling band lies in the plane of the membrane, before lipids are packed around it.
- **MEMEMBED** — A program that works out the orientation automatically. ProPrep's default.
- **PPM3** — Another orientation program, the method used to build the OPM database. Slower, and it lets ProPrep compare its result with OPM's.
- **OPM database** — The Orientations of Proteins in Membranes database: a public collection of membrane proteins already placed in a membrane. ProPrep can compare its own placement with OPM's (needs internet; reports only).
- **packmol-memgen** — The external program (part of AmberTools) that ProPrep calls to pack lipids, water and ions around the protein.
- **box padding** — The width of lipid kept between the protein and the side of the box (packmol-memgen's `--dist`). Default 15 Å.
- **water layer** — The thickness of bulk water above and below the membrane (`--dist_wat`). Default 17.5 Å.
- **molar ratio** — How many of one lipid are used for each of the others, e.g. POPC:CHL1 at 3:1.
- **loops (packing)** — Rounds of rearranging molecules. packmol-memgen does up to 20 per component and then up to 100 with everything together.
- **objective (packing)** — The number packmol lowers while packing; it measures how badly molecules still overlap. It rises and falls during a run; "lowest so far" is the value to watch.
- **force field** — The set of rules and numbers describing how atoms push and pull on each other. A membrane system needs one for the protein, one for the lipids and a matching water model. (Also used in Ch. 13.)

## Chapter 18 — Docking a ligand

- **docking** — Predicting where a small molecule sits in a protein's pocket, and how tightly, by trying many positions and scoring each one.
- **receptor** — The protein a molecule is docked into.
- **ligand** — The small molecule being docked.
- **redocking** — Taking a structure solved with its molecule in place, removing the molecule and docking it back. The test of whether the settings reproduce a known answer.
- **RMSD** — Root mean square deviation: the average distance between matching atoms of two structures. Under about 2 Å from the crystal pose means the docking reproduced reality.
- **search box** — The part of the protein that is searched. Can be set around the ligand, around named residues, or by typing the centre and edge lengths.
- **scoring function** — The formula that turns a pose into a number in kcal/mol. ProPrep offers vina, vinardo and ad4. Scores from different functions cannot be compared.
- **exhaustiveness** — How hard the search works. Higher is slower and more thorough; the default is 8.
- **PDBQT** — The file format AutoDock needs. ProPrep writes it for you from the protein and the molecule.
- **SMILES** — A way of writing a molecule as one line of text, e.g. `CC(=O)Nc1ccc(O)cc1` for paracetamol.
- **chain break** — A place where the protein chain continues in reality but the atoms are missing from the file. Such an end is left neutral; a real chain end is charged.
- **OXT** — The extra oxygen atom of a charged C-terminus. Many PDB files lack it, which makes docking fail unless the end is treated as a chain break.
- **campaign** — One receptor, box and settings run against every molecule of a library file, with a results table and a best pose per molecule.

## Chapter 21 — Sessions

- **session** — The file in which ProPrep writes down every answer you type, in order. Recording is automatic.

## Chapter 22 — Templates and batch

- **template** — A recorded session with some answers marked as blanks ("variables"), so the same steps can be run on many proteins.
- **batch processing** — Running a template over a list of proteins without typing.

## Terms used by other chapters that these chapters rely on

(Defined elsewhere; listed so the glossary stays consistent.)

- replay / hybrid mode — replaying recorded answers, then recording new ones live (Ch. 21 explains in words, no `\term{}` used).
- replay diverged — ProPrep was asked a question the recording did not expect, usually because the new protein differs (Ch. 22).
