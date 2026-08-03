# Iron-Sulfur 4Fe-4S Cluster Force Field Parameters for ProPrep

AMBER-compatible force field parameters for the [Fe4S4]^n+ cluster in three
oxidation states, packaged as `.frcmod` + `.lib` files for the ProPrep
iron_sulfur/4fe4s library.

**Multi-cluster compatible**: each parameter set uses a unique state-letter
prefix in atom-type names AND unique residue names for the 4 cysteinate
ligands (C01-C60). Two clusters in different oxidation/spin states can
coexist in a single AMBER topology without any naming conflicts.

Derived from B3LYP/6-31G(d,p) + GD3 + C-PCM (diethyl ether, eps=4) DFT
calculations followed by MCPB.py Seminario method for bond/angle force
constants and 2-stage RESP charge fitting on a large protein-cap model
(ACE-CYS-NME) with ff14SB CYM backbone constraints. Fe MK ESP radii
chosen as the IOD-weighted average of Li-Merz Fe^2+/Fe^3+ Rmin/2 values
(varies by oxidation state; see `methods.md`).

## Naming scheme

### Cluster residues

3-letter codes encode redox state and BS spin assignment:

| Code | Description | Cluster charge |
|------|-------------|----------------|
| **F** | Iron-Sulfur cluster identifier | |
| **O** | HiPIP oxidized state | [Fe4S4]^3+ |
| **R** | HiPIP reduced / Fd oxidized | [Fe4S4]^2+ |
| **D** | Fd reduced | [Fe4S4]^1+ |
| **1-6** | BS guess index (within state) | |

| Residue | State | BS partition | Spin variant | State letter | CM residues |
|---------|-------|--------------|--------------|---------------|-------------|
| FO1 | HiPIP_ox | {Fe1,Fe2}/{Fe3,Fe4} | a | A | C01, C02, C03, C04 |
| FO2 | HiPIP_ox | {Fe1,Fe2}/{Fe3,Fe4} | b | B | C05, C06, C07, C08 |
| FO3 | HiPIP_ox | {Fe1,Fe3}/{Fe2,Fe4} | a | D | C09, C10, C11, C12 |
| FO4 | HiPIP_ox | {Fe1,Fe3}/{Fe2,Fe4} | b | E | C13, C14, C15, C16 |
| FO5 | HiPIP_ox | {Fe2,Fe3}/{Fe1,Fe4} | a | G | C17, C18, C19, C20 |
| FO6 | HiPIP_ox | {Fe2,Fe3}/{Fe1,Fe4} | b | J | C21, C22, C23, C24 |
| FR1 | HiPIP_red | {Fe1,Fe2}/{Fe3,Fe4} | (single) | L | C25, C26, C27, C28 |
| FR2 | HiPIP_red | {Fe1,Fe3}/{Fe2,Fe4} | (single) | Q | C29, C30, C31, C32 |
| FR3 | HiPIP_red | {Fe2,Fe3}/{Fe1,Fe4} | (single) | R | C33, C34, C35, C36 |
| FD1 | Fd_red | {Fe1,Fe2}/{Fe3,Fe4} | a | T | C37, C38, C39, C40 |
| FD2 | Fd_red | {Fe1,Fe2}/{Fe3,Fe4} | b | U | C41, C42, C43, C44 |
| FD3 | Fd_red | {Fe1,Fe3}/{Fe2,Fe4} | a | V | C45, C46, C47, C48 |
| FD4 | Fd_red | {Fe1,Fe3}/{Fe2,Fe4} | b | W | C49, C50, C51, C52 |
| FD5 | Fd_red | {Fe2,Fe3}/{Fe1,Fe4} | a | X | C53, C54, C55, C56 |
| FD6 | Fd_red | {Fe2,Fe3}/{Fe1,Fe4} | b | Z | C57, C58, C59, C60 |

The total system charge of (cluster + 4 CMs) sums correctly to -1 / -2 /
-3 for the three states.

### Atom types (per state)

Each parameter set defines 12 unique atom types with the state letter as
the first character:

| Position char | Atom |
|--------------|------|
| 1, 2, 3, 4 | Fe1, Fe2, Fe3, Fe4 |
| 5, 6, 7, 8 | bridging S1, S2, S3, S4 |
| A, B, C, D | Cys SG of CM1, CM2, CM3, CM4 |

Example for FO1 (state letter A): atom types A1-A4, A5-A8, AA-AD.

The state letters (A, B, D, E, G, J, L, Q, R, T, U, V, W, X, Z) were
chosen to avoid collisions with old MCPB types (M*/Y*) and with common
AMBER atom-type prefixes for elements (C, H, N, O, P, S).

## Library structure

Each `.lib` file contains 5 AMBER residue units:

1. **Cluster residue** (one of FO1-FO6 / FR1-FR3 / FD1-FD6):
   - 8 atoms: 4 Fe (atom names FE1-FE4, types `<L>1`-`<L>4`) and 4
     bridging S (atom names S1-S4, types `<L>5`-`<L>8`)
   - 12 internal bonds (cubane connectivity)
   - No external bonds defined here; Fe-SG bonds to coordinating
     cysteines must be added by the topology builder

2. **4 modified cysteinate residues** (C01-C04 for FO1, C05-C08 for FO2, ...):
   - Same 10 atoms as standard CYM (N, H, CA, HA, CB, HB2, HB3, SG, C, O)
   - Backbone atoms: ff14SB CYM canonical charges (RESP-constrained)
   - Side-chain atoms: DFT-derived RESP charges
   - SG atom type: `<L>A`, `<L>B`, `<L>C`, `<L>D` (one per Cys position)

The 4 CM residues, together with the cluster residue, sum to the correct
total system charge for that redox state. Individual CM totals are NOT -1
(unlike standard CYM).

## ProPrep transformer logic

### Cys-Fe matching rule (invariant across all 15 parameter sets)

Standard PDB `SF4` representation has cluster atoms named FE1, FE2, FE3,
FE4, S1, S2, S3, S4 in a single residue. Use these atom names directly:

**Rule:** For each FE`<n>` in the input cluster, find the closest
cysteine SG. That cysteine becomes C`<idx>` where `idx` is determined by
the cluster's residue type and the n-th Cys position:

| Cluster residue | Cys 1 | Cys 2 | Cys 3 | Cys 4 |
|-----------------|-------|-------|-------|-------|
| FO1 | C01 | C02 | C03 | C04 |
| FO2 | C05 | C06 | C07 | C08 |
| ... | ... | ... | ... | ... |
| FD6 | C57 | C58 | C59 | C60 |

### Transformer steps

1. Detect a 4Fe-4S cluster in the input PDB (residue type SF4 or
   equivalent component code)
2. Read FE1, FE2, FE3, FE4 atom names from the cluster
3. For each FE`<n>`, find the bonded cysteine SG (closest, ~2.3 A)
4. Apply user-selected redox+spin state (e.g., FO1):
   - Rename cluster residue from SF4 -> the chosen FO1/FR2/etc. code
   - Atom types in the cluster: assign according to position (FE1 -> `<L>1`,
     S1 -> `<L>5`, etc.)
5. Rename the 4 cysteines to the 4 corresponding C`<idx>` residues
   (CM-of-Fe1 -> first C of the set, etc.); update SG atom types
6. Add 4 inter-residue Fe-SG bonds

### Loading order in tleap

The `.frcmod` must be loaded via `loadamberparams` BEFORE the `.lib` is
loaded; otherwise the `<L>1`-`<L>D` atom types will be undefined.

## Files

```
proprep_iron_sulfur_4fe4s/
├── README.md (this file)
├── oxidized/
│   ├── FO1/{FO1.frcmod, FO1.lib, DESCRIPTION.txt}
│   └── ... (FO2-FO6 same structure)
├── reduced_hipip/
│   └── FR1-FR3/...
└── reduced_fd/
    └── FD1-FD6/...
```

## Verification

All 15 parameter sets verified:
- Total system charge (cluster + 4 CMs) correct to within 0.001 e
- All 12 expected atom types present in each .frcmod
- All 5 expected units present in each .lib
- Cluster residue: 8 atoms, 12 internal bonds (correct cubane)
- CM residues: 10 atoms each, ff14SB-compatible backbone

## Multi-cluster usage example

A protein with two 4Fe-4S clusters in different oxidation/spin states (e.g.,
oxidized FO1 and reduced HiPIP FR2) can be set up by loading both parameter
files:

```tleap
loadamberparams FO1.frcmod  # types A1-AD
loadamberparams FR2.frcmod  # types Q1-QD (no overlap)
loadoff FO1.lib             # units FO1, C01-C04
loadoff FR2.lib             # units FR2, C29-C32 (no overlap)
```

All atom types and residue names are disjoint between the two parameter
sets, so they coexist without conflict.
