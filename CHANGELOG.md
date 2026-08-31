# Changelog

All notable changes to ProPrep are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and ProPrep adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Versions correspond to the `proprep` package on the
[mjgplab](https://anaconda.org/mjgplab/proprep) conda channel and to tagged
releases of the source repository. Builds of the same version are reserved
for recipe-level fixes (dependency repins, packaging corrections) that do
not touch the source.

## [Unreleased]

## [1.17.0] — 2026-08-31

### Changed

- PDB Filter water analysis, metric 4: burial is now computed from geometry
  alone instead of an atom count normalised by an uncalibrated constant.
  Each water oxygen gets its Lee-Richards accessible area against the
  protein and hetero atoms (1.4 Å probe, Bondi radii; other waters never
  occlude), a bulk/enclosed flag from a flood fill of probe-accessible
  space, and its nearest heavy atom. Categories are physical statements:
  Clash (nearest C/N/O/S/P atom under 2.2 Å, the wwPDB close-contact
  criterion), Enclosed (no path to bulk solvent), Buried (0 Å² but
  bulk-connected), Exposed. The count-based `burial_max_expected` and the
  `40 Å² × (1 − burial)` "estimated SASA" are gone. On 3WL2 the old score
  put every enclosed cavity water at 46-68%, and nothing above 80%. The
  multi-radius and directional profiles (metrics 6 and 7) keep the count and
  now have their own parameter entries.

- The water burial table shows a Covered column (100 × (1 − SASA / the
  area of a free water)) and prints every category rule with its cutoff,
  in the order the rules are applied, both in the method panel and in a
  legend under the table; the panel also explains that SASA is measured on
  the contact sphere while enclosure is looked for one grid cell beyond it,
  which is how a water can be untouchable yet bulk-connected. The unit
  symbol had been U+0172 "Ų" (U with ogonek) and now reads Å².

- Metrics 6 and 7 (multi-radius and directional atom counts) state the
  conventions behind their labels: the "saturation" radius is the first
  step with less than 10% growth, the directional "pattern" compares the
  sector-count range with the mean at 0.5× and 1.5×, and the compass glyphs
  are quarters of that water's own largest sector. Each label carries its
  rule, and both tables gained a legend saying so; the menu names them
  atom-count profiles.

- The structure fixer's alternate-location picker says which viewer colour
  each alternate is drawn in ("■ red in viewer") beside its occupancy, and
  the viewer's representation list is labelled "Alt A (red), occ 0.60"
  instead of "Altloc A", so the colours on screen can be matched to the
  numbers in the prompt.

- Water analysis viewer halos and the "Ordered" category. The viewer reps
  were labelled `Water Cat Hbond` / `Water Cat Ordered`; they now say the
  rule and its cutoff ("Waters: ≥ 3 H-bond partners (≤ 3.5 Å)"). With the
  burial metric selected the halos follow the burial categories (clash,
  enclosed, SASA 0, covered ≥ 90 %, covered 50–90 %) instead of the
  recommendation categories, which hid a buried water behind "Highly
  connected". "Ordered" was B-factor < 30 Å²; B-factors scale with
  resolution and refinement, so it is now B below the median of this
  structure's protein heavy atoms, and the B-factor table shows that ratio
  and the median instead of 20/40/60 Å² bins.

- Water burial treats a water as one 1.4 Å sphere whether it is the probe
  or a crystallographic water. The crystallographic water's oxygen had been
  given Bondi's 1.52 Å, the radius of an oxygen atom inside a molecule, so
  the same species had two sizes and two waters "touched" at 2.92 Å where
  real ones touch at 2.80 Å. Contact distance is now 2.80 Å, the isolated
  water area 98.5 Å², and the enclosure reach 3.30 Å; protein, heteroatom,
  and metal occluders keep Bondi radii.

- Water hydrogen-bond partners are ranked by distance alone. The previous
  ranking multiplied an "angle quality" (bins at 150/120/90° scored
  1.0/0.8/0.6/0.3) by a distance bin and then by 1/distance again. The
  angle it measured was the one at the partner atom, X–A···O, whose ideal
  value depends on hybridisation (about 120–160° at a carbonyl oxygen, about
  110° at an sp³ oxygen or amine), so "linear is best" rated a textbook
  Ser OG contact "acceptable" and an on-axis approach "excellent"; water
  partners never got an angle at all. None of the numbers had a source.
  Candidates within the 3.5 Å heavy-atom cutoff are now listed closest
  first and truncated to the configured maximum; the unused
  `hbond_angle_cutoff` parameter is gone.

- The water classification cascade is gone. Every analysed water used to be
  given one label (Metal-coordinating, Highly connected, Interface, Buried,
  Ordered, or Bulk solvent) by a first-match-wins rule chain, and the viewer
  halos followed it, so a metal-bound water that was also enclosed showed
  only as metal-bound. That is a single score under another name and runs
  against the rest of the analysis, which reports each metric on its own.
  The viewer now draws one highlight group per fact each displayed metric
  establishes (metal within cutoff; 1, 2, 3, or 4 H-bond partners; B below
  the protein median; the burial categories and coverage bands; interface),
  each labelled with its rule and cutoff, and a water may appear in several.

- Water burial counts metal ions as occluders by default. Metal ions are
  their own residue class in the burial code, and the default set was
  `protein,hetero`, so a coordinated ion neither shielded its water nor
  blocked the enclosure flood fill. Default is now `protein,hetero,metal`.

### Fixed

- Coordinated waters are part of the transformer connectivity fingerprint.
  The Weisfeiler-Lehman fingerprint that lets a reused transformer tell
  same-name residues apart excluded waters entirely, so water renames (say
  HOH to MW1 and MW2 on a di-metal site) could never auto-resolve and always
  fell to the manual prompt. Every emitter now bakes two fingerprints per
  role, the anhydrous one exactly as before and a hydrated one with waters
  as nodes; matching tries the hydrated bijection first and falls back to
  the anhydrous one, so transformers emitted before this change, and reuse
  sites whose water-metal bonds are left undefined (a restraint instead),
  match exactly as they did, while water roles resolve automatically when
  the bonds are defined on both sides.

- Structure Viewer: colouring by "Chain ID" coloured a multi-chain
  structure almost uniformly red. The option used NGL's `chainid` scheme,
  which keys on the parser's internal chain record (a new one at every TER
  and every polymer-to-HETATM break), so 9YUQ's 16 chains became 112
  records and the protein segments, first in the file, all sat at the red
  end of the scale. The option now uses `chainname`, the chain letter;
  saved scenes carrying `chainid` are mapped.

- Replaying a session recorded from `proprep --pdbid X` or `--pdbfile F`
  now reloads that structure first. The argument-driven load happens in
  `main()` with no prompts, so such a session holds no loader interactions,
  only `metadata.pdb_id`/`pdb_file`; replay printed that and then started
  at the main menu with no structure, and the first module diverged. A
  `pdb_file` recorded on another machine is found by basename in the
  project directory; if it is absent the replay says which file to copy in.

- Constant-pH/redox namelist variables (icnstph/solvph/ntcnstph, icnste/
  solve/ntcnste, saltcon) reach every production `simulation.mdin`. Only
  the live runner injected them; the batch directory writer, the
  standalone writer, and the three SLURM writers staged the cpin and passed
  `-cpin -cpout -cprestrt` but wrote the mdin untouched, so sander would
  have run plain MD with a cpin on the command line. An imported .mdin is
  now read, given its configured restraints and the titration namelist,
  and written, instead of being copied verbatim.

- The MD Manager's recommended engine assignment numbered NPT steps across
  the whole queue, so with several structures only the first structure's
  density equilibration was treated as early NPT (CPU) and every other
  structure's went to the GPU with a fixed PME grid. NPT steps are now
  numbered within each structure's own workflow.

- Every generated microstate tLEaP script starts with `logFile
  <microstate>_leap.log`, so each build writes its own log instead of all
  of them (and the info passes before them) appending to one `leap.log`
  in the working directory. Parallel builds copy that log back beside the
  inputs; the info passes log to scratch files that are removed with the
  script. The log parser and the post-run message display read the
  script's own log.

- The titration-file step accepts several topologies at once (a number, a
  comma list, a range, or `all`, which is the default). Residue and
  initial-state choices are made on the first topology and reused for any
  later one whose titratable residues match it exactly; a topology that
  differs is prompted on its own. One `titration_configs` entry is kept per
  topology alongside the legacy single `titration_config`.

- The MD Manager's constant-pH/redox offer looks up each structure's own
  titration files by prmtop name and lists every structure it applies to
  before asking once. Previously it read the single most recent config, so
  a set of five microstates was offered only the fifth's cpin.

- Batch microstate generation offered the constant-E HEH heme library next
  to the HCO library in the Topology Generator (and only HCR for the
  reduced state). The batch path recorded no redox treatment, so the set
  picker had nothing to filter on; choosing HEH loaded a library that
  defines no HCO, tLEaP built the hemes as untyped, chargeless atoms, and
  every microstate reported the same net charge. Batch microstates
  enumerate redox states explicitly, one topology each, so the treatment is
  fixed_E by construction; it is now recorded as such and the picker offers
  only the fixed_E sets. The microstate info pass also refuses to continue
  when tLEaP reports an unknown residue, instead of printing a charge that
  omits it.

- The per-water FreeSASA value that drove the "Buried" recommendation was
  always 0.0: `Result.atomArea()` takes an atom index and was given selector
  strings, and FreeSASA's PDB reader drops HETATM records so the water was
  not in the structure being queried. Every water not otherwise categorised
  was recommended as "Buried". The recommendation now uses the burial above.

### Added

- `proprep-library`: snapshot, reset and restore the user library
  (`~/.proprep/forcefield_params` and `~/.proprep/transformers`) without
  ever deleting anything. `snapshot [name] -m note` copies the library to
  `~/.proprep/library_snapshots/<name>/` with a manifest; `reset` writes a
  `before_reset_<stamp>` snapshot and then empties the library so a
  practice run starts clean; `restore <name>` writes a `before_restore_`
  snapshot and brings a snapshot back; `list` and `status` report. Settings,
  keys, templates, profiles and web sessions are outside the library and
  untouched. Restart ProPrep after a reset or restore.

- The Structure Viewer can save and reload a scene. **Save Scene** in the
  viewer's View Controls (or *Save the scene shown in the open viewer* in
  the Structure Viewer menu) writes `<name>.scene.json` into the project
  directory: every representation as shown (including ones added by hand
  in the sidebar), the camera, the camera type, and the background, with
  structure paths relative to the project so it travels with it. *Load a
  saved scene* lists the files and restores the view; a scene whose
  structure files are missing is refused with the paths named. Hand-added
  representations also now survive an annotation refresh from the CLI,
  which used to discard them.

- The Topology Generator's cpin step now also generates `cein` and `cpein`
  files, so a structure prepared for constant-redox MD can actually titrate.
  It reads the topology and generates one file per titratable family:
  `cpinutil.py` for the pH-titratable residues, `ceinutil.py` for the redox
  heme `HEH`, and both files when both families are present. A `constant_E`
  heme with titratable `PRN` propionates therefore yields a cpin *and* a cein,
  which are complementary and are passed to sander together, following Amber
  tutorial 33. Generating only one family is still offered, for holding the
  heme at a fixed oxidation state or the propionates at fixed protonation.
  `cpeinutil.py` is used only for a residue that titrates in proton and
  electron on the same site, which neither of the other two will accept. The
  MD Manager carries every generated file through to the run scripts, emitting
  each one's flag triple and chaining each one's restarts independently, and
  sets `icnstph`/`solvph`/`ntcnstph` and `icnste`/`solve`/`ntcnste` together in
  one `&cntrl` when both apply. Constant-pH runs are unchanged.

- The initial oxidation state of a redox residue can be set from a target
  potential, the redox counterpart of setting protonation states from a target
  pH. A potential above the residue's standard potential starts it oxidized,
  below starts it reduced; `HEH` has Eo = -0.203 V. In a combined run the
  propionates still follow pH against their pKa independently.

- The bis-histidine c-type heme can be prepared for constant-redox MD. A new
  `redox_treatment` choice sits alongside the existing propionate pH treatment
  and offers `constant_E` beside the existing `fixed_E`. Under `constant_E` the
  heme residue is named `HEH`, which is what AMBER's `conste.lib` calls it and
  what `ceinutil.py` matches on, so the heme reaches the cein file with both its
  ferric and ferrous charge vectors. The oxidation state is then chosen when the
  cein is generated rather than built into the topology, so one library serves
  both states, exactly as AMBER ships a single `HEH` unit. The redox-state
  question is skipped when `constant_E` is chosen, since there is no state to
  commit to at build time. The two treatments cross with the two propionate
  treatments, so the cofactor now offers four sets: `HEH` with titratable `PRN`
  for combined constant-pH/redox runs, `HEH` with static `PRP`/`PRD` for
  redox-only runs, and the two existing `HCO`/`HCR` sets unchanged.

- Self-contained installers. `constructor/` builds `ProPrep-X.Y.Z-<OS>-<arch>.sh`
  installers that carry the whole environment: ProPrep, AmberTools, MODELLER,
  `reduce`, `pdb2pqr` and every other dependency. Attendees and new users run
  one file and need no conda, no channel solve and no network during the
  install. MODELLER comes with its placeholder license; users supply their own
  key at run time (below). Built for Apple-silicon and Intel Macs and for
  Linux x86-64 (which is also the Windows route through WSL2); see
  `INSTALL.md`.

- The MODELLER license key can be supplied at run time. ProPrep reads
  `$KEY_MODELLER` or `~/.proprep/modeller_key` before the first `import modeller`
  and answers the `modeller.config` import itself, so no distributed artifact
  ever contains a key and a wrong key is fixed by editing a file rather than
  reinstalling. Outside a bundle this is a no-op when neither source exists,
  so keys baked in by `install_proprep.sh` at install time keep working; a
  runtime key wins over a baked one.

- Web shell: sessions now belong to *seats* that outlive their websockets.
  A dropped connection, a page reload or a laptop going to sleep no longer
  kills the ProPrep process; the page reconnects with backoff and receives a
  replay of the recent output, then carries on in the same session. A second
  tab opened on the same seat displaces the first (one connection per seat, as
  before).

- Web shell hosted mode. `proprep-web --seats N [--seats-dir DIR]
  [--public-url URL]` serves N independent seats from one process, each with
  its own working directory and a random link token persisted in `seats.json`
  so links survive restarts. Opening `/?seat=<token>` sets a cookie, after
  which the terminal, the docked structure viewer and the new **Download
  project** button (a zip of the seat's directory) are all seat-scoped; a
  missing or wrong token is refused. Hosted mode never auto-shuts down. The
  default single-user launch is unchanged apart from reloads resuming instead
  of restarting.

- `proprep-web` is now also available as a standalone executable next to
  `proprep` when the PyInstaller bundle is built.
### Changed

- Depositing parameters and saving a transformer ask for three names in a row
  — library entry, parameter set, template — whose wording does not say which
  level of the library each one sits at. The import wizard now shows entry
  (a directory, one per molecule or site) against set (a key inside it, one per
  derivation, several per entry) with the path each will occupy, and the
  transformer creator says the template is the rename recipe rather than the
  entry holding the parameters, naming the file it writes.

- Web shell: the structure viewer announces itself to the shell over the
  loopback bind address with a per-seat path, never through the public
  address the browser used. That callback used to take the browser's `Host`
  header and broke behind any reverse proxy, TLS terminator or tunnel; it is
  now also refused from non-loopback callers.

- Structure repair via MODELLER on a headless host (a plain `ssh` session) no
  longer launches a text browser against the viewer page. `webbrowser` falls
  back to `lynx`/`w3m` whenever `TERM` is set, which took over the terminal or
  segfaulted on the WebGL page. Nothing is opened without a display or when
  the chosen browser is text-mode (an explicit `BROWSER` still wins), and the
  launch message prints the `ssh -L` line that makes the URL reachable from
  the user's own machine. `proprep-web` gets the same guard.

- Fe4S4 parameter sets: the Cys C-beta/S-gamma bond is now derived from each
  solution's own Hessian instead of the ff14SB default that MCPB.py inserts
  (about 1.9x too stiff and 0.05 A too short; the derived values are
  bracketed by two independent published derivations). Two of the fifteen
  broken-symmetry solutions (FO4, FR2) were re-derived from unconstrained
  Hessians so all fifteen follow the same procedure; the averaged FO0/FR0/FD0
  sets are recomputed accordingly.

- Release tooling. `docs/RELEASE_PROCEDURE.md` documents every release
  endpoint in order (source tag, conda channel, public snapshot repository and
  its GitHub release, Zenodo, constructor installers). `make check-version`
  now checks all seven version-carrying files, including
  `update_proprep_in_ambertools.sh` (which had sat at 1.14.0 through two
  releases) and `constructor/construct.yaml`; `make snapshot` exports the
  public snapshot reproducibly.
### Fixed

- The titration step refuses `igb` and `intdiel` combinations that have no
  reference energies instead of writing a file full of nulls. `ceinutil.py`
  accepts any of igb 1/2/5/7/8 and intdiel 1/2, but the only redox residue
  that exists has energies for igb 2/5/7/8 in implicit solvent, 2/5/7 in
  explicit, and none at all for intdiel 2. The unsupported combinations were
  accepted silently and produced meaningless output.

- Residues that a titration utility would drop are now reported before the
  file is written. All three utilities filter chain termini by comparing a
  residue's atom count against ParmEd's definition and skip any mismatch with
  no message and a zero exit status. That is the intent for a terminal
  aspartate, but a cofactor is a different matter: `HEH` is one 87-atom
  residue spanning the porphyrin and side-chain fragments of two cysteines and
  two histidines, and a build one atom adrift simply stops titrating. Worse,
  once any residue is dropped the `-states` list no longer matches the residue
  count, so every hand-picked initial state is silently discarded in favour of
  defaults. The count is now checked against ParmEd up front and again against
  the generated file.

- Explicit-solvent constant-redox runs get the radii-corrected topology they
  need. `ceinutil.py` imports ParmEd's `changeRadii` and `change` actions and
  never calls them; it has no `-op` flag at all, unlike `cpinutil.py` and
  `cpeinutil.py`. When a cpin is generated alongside the cein, `cpinutil.py`
  writes that topology and both files share it. For a redox-only run ProPrep
  applies the same two actions directly and records the equivalent `parmed`
  input in the generated script.

- The titration step no longer asks for a tLEaP input file. It used that file
  only to look up two paths, a leftover from when residue numbers had to be
  mapped from the original PDB onto the tleap-renumbered topology; scanning
  the topology directly removed that need. The topology is now selected
  straight from the workspace.

- A constant-pH bis-His c-type heme topology no longer hides the heme from
  AMBER's redox tooling. The library named the residue `HCO` or `HCR`, and both
  `ceinutil.py` and `cpinutil.py` select residues purely by the name recorded in
  the topology, so the heme was passed over without comment and never titrated.
  The parameters themselves were always AMBER's: the shipped unit is identical
  to `conste.lib`'s `HEH` atom for atom, in name, type and charge. Only the name
  differed. Choose the `constant_E` treatment to emit `HEH`; the fixed-redox
  sets keep `HCO`/`HCR`, which is correct for them because their oxidation state
  really is fixed.
- The bis-His c-type heme library now carries the four angle and improper terms
  the titratable `PRN` propionate needs. `PRN` types both propionate oxygens as
  `OH` and gives each up to two ghost hydrogens, an arrangement no standard
  residue has, so it needs `HO-OH-HO`, `OH-C-OH` and two impropers that live in
  `frcmod.constph`. The heme library's own header said it carried them and its
  methodology said no constant-pH force field need be loaded separately, but
  neither was true: building `PRN` against the heme library plus `leaprc.conste`
  alone failed with three missing-parameter errors. The terms are copied
  verbatim from `frcmod.constph`, and since none of them exist in `parm10` and
  none can form on an ordinary hydroxyl, adding them leaves every other residue
  untouched.

- A workflow step that pauses at a checkpoint was recorded as *completed*. The
  step stops early to wait for a calculation run outside ProPrep and has
  produced none of its artifacts, but the green tick satisfied the next step's
  dependency and let the run walk past it — the failure then surfaced two steps
  later, naming the step that consumed the missing file rather than the one
  that never wrote it. A checkpoint is now `in_progress`: `[n]ext` offers it
  again, and a step depending on it warns before running.
- The modified-amino-acid RESP step reported "Missing ESP or AC file from
  previous steps" without saying which, immediately after step 7 had announced
  the AC file. It now names the artifact, the filename it expected, the
  directory it searched, and the reason step 6 may have produced nothing.
- The step-5 RESP selection was held only in memory, so a resumed session
  re-selected from scratch — silently changing which conformers feed the fit
  and demanding ESP single points the user had not run. It is now recorded in
  the workspace, keyed by residue, and reused on resume.
- An imported entry recorded its residue name from the `.lib` FILENAME rather
  than the unit inside it. `saveoff` names an entry after the tLEaP variable it
  was given, so a `GDP.lib` can hold unit `gdp` — and tLEaP matches the unit,
  case-sensitively. The unit name is now read from the library, and a
  divergence from the filename is reported.
- Depositing a set could destroy another set's parameter files. Sets are
  metadata keys, not directories, so every set in one redox/spin state shares a
  folder, and the copy was an unguarded `shutil.copy2`: two sets whose files
  shared a basename overwrote each other while both metadata leaves still
  named the file. A colliding file is now written alongside (unless it is
  byte-identical, where sharing it is what the other set already has), and a
  rolled-back deposit restores anything it overwrote instead of deleting it.
- Saving a transformer silently replaced any existing one of the same name —
  a plain `open(path, "w")`, with sanitizing that collapses `GDP ff94` and
  `gdp-ff94` onto the same file. The name is now checked before the
  force-field questions, showing what is already there.
- When Amber's parameter files could not be found, the import wizard listed
  every atom type a frcmod declares as new — standard types like `CT` and `OS`
  included — with nothing saying why. Widening is the safe direction, but it
  now says the check could not be made and where it looked, since the usual
  cause is running ProPrep outside the environment AmberTools is installed in.

- The library-collision prompt asked "How would you like to proceed?" and
  listed nothing to proceed with. `options_map` is metadata for the session
  recorder — it never renders, and passing it suppresses Rich's inline choice
  list — so the three options existed only in the source. They are printed now,
  and the version-bump option names the set it would actually write. The
  collision message also lost its `[redox/spin]` coordinates to Rich's markup
  parser, which read the brackets as a style tag; exception text is escaped
  before printing.

- The interactive transformer creator could not address any residue in its own
  table. BioPython reports a residue with no insertion code as a *space*, and
  that space reached the editor's structure intact, while every command looks a
  residue up with an empty insertion code; the comparison is exact, so
  `rename_atom A 302 O5' O5*` answered `No residue A/302 in the site.` for a
  residue printed one line above. Insertion codes are now normalized where the
  structure is built and where it is queried.
- New residue and atom names typed into the transformer creator were
  force-uppercased, so a library whose tLEaP unit is lowercase could not be
  matched by a rename — the editor could not perform the edit it exists to
  perform, and the library had to be hand-edited instead. Names are now
  recorded exactly as typed. Names being matched *on* may still be typed in any
  case, and the operation records the structure's own spelling so replay
  against a PDB still matches.

- The PyInstaller bundle had drifted from the tree: `scipy` was excluded but
  is a runtime dependency of pb_titrate and trajectory analysis, the
  pb_titrate model compounds and the cluster profiles were not bundled, and
  `qmmm_prep` was missing from the hidden imports. The bundle also carried the
  build machine's `modeller.config`, license key included; it is now excluded.

- ProPKA pKas for heme propionates were dropped entirely. ProPKA's `.pka`
  file identifies ligand groups by residue name and atom name without a
  residue number, so the summary parser skipped every `PRN`. The groups are now
  recovered through ProPKA's Python API, which keeps each group's owning atom,
  and merged additively (values parsed from the file win). Group type alone is
  not the key: a C-terminal Lys carries a side-chain and a terminus group on
  one residue and the terminus pKa was overwriting the side-chain value.

- ProPKA pKas for residues numbered 1000 and above were dropped. The summary
  column is fixed-width and loses the space between name and number at four
  digits; on a six-protomer cytochrome filament this silently lost the whole
  sixth protomer.

- Several steps reported success they did not have: a checkpoint paused for an
  external calculation was recorded as completed and the failure surfaced two
  steps later naming the wrong step; the RESP step named two missing files
  without saying which; a resumed modified-amino-acid session silently
  re-selected conformers; an imported library entry took its residue name from
  the file name rather than the unit inside it; depositing a set could
  overwrite another set's files, and rollback then deleted them; saving a
  transformer replaced an existing one of the same name without asking; and
  the import wizard listed every atom type as new when Amber's parameter files
  could not be found, without saying so.

- Re-importing a parameter set that already exists asked how to proceed
  without showing the choices (the options were passed only to the session
  recorder, which also suppressed the inline `[1/2/3]`), and Rich swallowed
  the bracketed name of the colliding set. The options are printed, the
  version-bump choice names the set it would write, and exception text is
  escaped.

- The transformer creator's interactive editor could neither address nor
  rename residues of an imported library: BioPython reports a missing
  insertion code as a space, which never matched the empty string the
  commands supply, and new names were force-uppercased, so the one rename
  that would match a lowercase tLEaP unit was refused. Insertion codes are
  normalized at the boundary and names are recorded verbatim.

- MCPB RESP input: elements were inferred from the first letter of the atom
  name when PDB columns 77-78 were absent, so `ZN` in MCPB.py's own 65-column
  `*_large.pdb` read as element `Z` and aborted. Column justification,
  monatomic-ion records and an allowlist of metal/halide names are tried
  first. Stage 2 also refit every atom under its stronger restraint, damping
  metal and ligating-atom charges toward zero; only the equivalenced atoms and
  the heavy atom each hydrogen group hangs off are now freed.

- Metal-site parameterization: a ligand parameterized standalone and again as
  part of the site (renamed in the process) left both registrations in the
  preprocessing lists, so tLEaP loaded two definitions of one ligand. Step 4
  now drops the superseded halves and de-duplicates by real path.

- The MoCo molybdenum ESP radius is pinned; tests that compare against
  MCPB.py's MK radius table are skipped where AmberTools is absent instead of
  reading as defects.

- PDB Filter water analysis: choosing metric 4 (proximity-based burial)
  crashed with `'PDBFilterWorker' object has no attribute 'parameters'`.
  The burial table read its radius, atom types, and weighting from the
  filter worker instead of from the `WaterAnalyzer` that owns them; it is
  now handed the analyzer like the metal and interface tables. Present
  since the module was written, so burial analysis had never worked.

## [1.16.0] — 2026-08-19

### Added

- Externally obtained parameters can be imported and used. The import wizard
  browses for the `.frcmod`, `.lib`/`.off` and optional `.prep` with the same
  file browser the rest of ProPrep uses, shows what each parameter category
  means and where it will be stored, and reads any atom types the frcmod
  declares so they reach tLEaP as `addAtomTypes` entries. Generated parameters
  cannot introduce new types, which is why the wizard never asked; imported
  ones can.
- A transformer can be saved that only binds a force-field library, with no
  renaming. Renaming and parameter-binding are separate jobs: MCPB output needs
  both, while a cofactor already named as its library names it needs only the
  binding. There was no way to express the second, so an imported cofactor's
  library was unreachable. The library and its redox/spin state are chosen from
  lists of what exists rather than typed as a path.
- Deposits record which force fields their parameters require, derived from the
  atom types in the files. The Topology Generator already enforced declared
  prerequisites, but almost nothing declared any.

- Pure inorganic metal clusters can be given hydrogens before parameterization.
  Hydrogen addition covered protein and organic residues only, and `reduce` has
  no chemistry for a Mo-S-O or Fe-S core, so a cofactor whose resting state
  carries a hydroxo — `Mo(=O)(=S)(OH)` in a molybdenum cofactor — reached the QM
  model as a bare oxo with the wrong electron count and charge, and there was no
  way to correct it. Editing the Gaussian input by hand is not an alternative:
  it and the model PDB are matched by index, so an atom added to one shifts the
  Seminario force-constant indices and leaves the deposited residue template
  without the hydrogen its charges were fitted with. Offered for every cluster,
  defaulting to no.
- The structure viewer focuses on a cluster while its hydrogen-addition prompt
  is on screen, and re-reads the file once a hydrogen is added so the new atom
  is visible. `ViewerCoordinator.refresh_structure()` re-serves the current
  path, which `show_structure` deliberately will not do — it treats a repeated
  path as a no-op, correct until the file is edited in place. It never starts a
  viewer or opens a tab.
- A hydrogen added to an atom with only one bond is now placed at an angle to
  that bond, rotated to the least crowded side, instead of directly opposite it
  — a hydroxo, thiol or amine hydrogen is bent, never linear. The opening angle
  can be set at the prompt.

### Changed

- The cofactor prerequisites panel reports what the selected parameter sets
  declare instead of describing their chemistry. It asserted that every
  cofactor needs a protein force field and inferred GAFF2 from residue names,
  while explaining the requirement in terms of a ribitol tail and which bond
  would fail — none of which is knowable about an arbitrary parameter set.
- MCPB step results are stored per site. A single shared record meant `step_1`
  belonged to whichever site ran last, which cross-wired one site's atom-type
  fingerprint and RESP charge constraint onto another.

- The session editor table labels its status column `Status` instead of `St`.

- The parameterization prompt selects **sites**, and each selected site goes to
  the parameterizer its own category implies. It no longer groups residues:
  that was for a modified amino acid covalently bound to a cofactor, which is
  now expressed by defining the pair as a site in the Redox Site Detector, so
  the combine/separate and category-conflict prompts are gone. A mixed
  selection is ordinary rather than a conflict.
- Selecting metal sites now parameterizes those sites and no others. MCPB
  previously re-detected and parameterized every metal site in the structure
  regardless of the selection, so a structure with two equivalent Fe2S2
  clusters could not have one derived and the other served by the reuse
  transformer. Sites left out are named, with a note that every metal site
  needs parameters before the topology build will succeed.

### Fixed

- A Gaussian output that describes a different model than the input beside it
  is refused rather than fitted. Re-running the atom-typing step rebuilds the
  models but Gaussian is run by hand, so a leftover log gives the Hessian or
  ESP of a superseded model — and the result still looks like force constants
  or charges. Compared by content: atom count, elements, then coordinates.
- A withheld metal cluster's own bonds get force constants. Fe-S inside an
  Fe2S2, Mo-S/Mo-O inside a Mo cofactor and the O-H of a hydroxo reached
  neither the coordination bond list nor the prmtop, so nothing derived a
  parameter for them and tLEaP reported one missing for each.
- Metal covalent radii are used when perceiving a cluster's bonds. Molybdenum
  fell to a carbon-like default, so a Mo cofactor was deposited with two of its
  four bonds; Fe-S was passing only because it sat just under that accidental
  cutoff. Metal-metal pairs are excluded — a cluster's metals are bridged
  through their ligands.
- An inferred hydrogen is typed for the atom it is bonded to. Amber names a
  hydrogen after its neighbour, so a metal-bound hydroxo proton typed `H` took
  the amide hydrogen's van der Waals radius where the hydroxyl convention `HO`
  is zero.
- RESP scaffolding residues are constrained to their own charge rather than to
  zero. Keeping a real ARG as a gap bridge left its +1 nowhere to go but the
  neighbouring coordinating residues, one of which came out positive.
- Charges are found for metal-site atoms again. The lookup is keyed by
  coordinate tuple, and BioPython's float32 does not compare equal to the
  float64 a resumed session restores, so every site atom summed as zero and a
  −2 site was proposed as 0.
- Atom types are read from a library once. The section test also matched
  `atomspertinfo`, whose rows repeat every name, so an 84-atom residue was
  reported as 168.
- A supplied library is matched on heavy atoms when the structure has no
  hydrogens, instead of demanding a hand-written mapping for a library that
  fits.
- Caps stay in the chain when the PDB is reordered for tLEaP. Caps bracketing
  an unfilled internal gap were moved to opposite ends of the chain, leaving a
  70 Å peptide bond and the gap they had guarded open.
- The PDB written back from a topology keeps the topology's atom names. They
  were translated to PDB v3 conventions, so a library using older names — `O1P`
  rather than `OP1`, say — built once and then failed against the file it had
  just produced.
- Triage categories survive a resume. They lived only on the instance, so a
  resumed run found no organic residues, no clusters and no isolated metals,
  and ticked those steps off as complete.
- Resuming a pending small-molecule or modified-amino-acid parameterization
  works; both dispatched to methods that were never written.
- Force Field Integration names a deposited frcmod for its library entry rather
  than the working directory's `site_N`.

- Implicit solvation for the large model is now asked for rather than
  inherited. The prompt sat entirely under an anion check, so an anionic large
  model silently adopted the small model's SCRF (announced only after the fact)
  while a neutral or cationic one dropped it with no message at all, leaving the
  two calculations at different levels of theory for no stated reason. It is now
  offered in every case, defaulting to yes whenever the small model was solvated,
  and declining says plainly that the models will differ.
- The Merz-Kollman ESP radius is taken from MCPB.py's own `vdwRadiiDict2023`
  (Smith et al. JCTC 2023, 19, 2064) rather than from the force-field IOD
  parameters. They are different quantities — the MK value decides how close to
  the nucleus ESP grid points may fall, the IOD value is a Lennard-Jones term —
  and the IOD tables stop at tetravalent, so a Mo(VI) cofactor could not resolve
  one at all and fell back to a generic 1.5 A with no cited source. Force-field
  parameters are unchanged. `tools/patch_readradii.py` adds the block to
  `large_resp.gjf` files generated before this, so they need not be rebuilt.
- Checklist state stores numbers as numbers. The serializer recognised only
  Python `int`/`float`, and a numpy scalar is neither — `np.float32` fails
  `isinstance(v, float)` and has no `__dict__` — so every coordinate read from
  BioPython was written as `{"__type__": "str", "value": "-46.078"}`. On resume,
  metal reinsertion then handed PDBIO a dict where it wanted a number and the
  step failed with "must be real number, not dict". numpy scalars now serialize
  as plain numbers, and `MetalInfo.from_dict` unwraps and coerces so a state
  file written before this still resumes.
- Resuming a metal-site run from saved checklist state no longer fails at
  structure recombination with `tLEaP error: 'site_id'`. Checklist state stores
  objects wrapped as `{"__type__", "value"}` and uses the dataclass field names,
  while the exporters write a flat dict with three fields spelled differently
  (`coordinates` for `coords`, and so on). `dict_to_redox_site`, the documented
  normalizer for exactly this case, understood only the exported shape; it now
  accepts both.
- A checklist step that reports a failure without raising is recorded as failed
  rather than completed. Structure recombination printed "tLEaP failed" and was
  ticked off, so the following step ran and died on the prepared structure it
  had never produced, hiding the real cause. Handlers signal this by returning
  `success: False`.
- RESP is constrained to the charge its own ESP was computed with. Per-site
  results are restored from a single shared workspace key, so in a multi-site
  run the charge came from whichever site ran last: a -1 Fe2S2 model was fitted
  against a -3 constraint and RESP spread the missing two electrons over the
  site, reaching +5.6 on a metal. It does not fail on a mismatched total — it
  fits a different molecule. Charge and multiplicity now come from the site's
  own Gaussian log (falling back to its input file), which cannot be another
  site's, and a disagreement with the stored step-1 value is reported.
- The large model's Gaussian input carries the metal van der Waals radii again.
  They were gathered by looking each `redox_site.centers` entry up in the atom
  assignments, which only works for a lone metal ion, where the center is the
  metal atom. For an organometallic cofactor or a pure cluster the center
  describes the residue instead, so nothing matched and the file was written
  with a bare `Pop=MK` and no radii block — the Merz-Kollman ESP then sampled
  the metals with Gaussian's defaults. The radii now come from the metal atoms
  themselves, and a site that resolves none says so instead of quietly falling
  back.
- The large model's ESP calculation now defaults to the same level of theory as
  the small model rather than `HF/6-31G*`. The MCPB.py tutorial states it
  directly: "we used the B3LYP/6-31G* level of theory to perform the
  calculations for both the small and large models". `HF/6-31G*` is the generic
  RESP convention for organics, and it silently differed from the functional
  just chosen for the small model.
- A hydrogen inside a metal cluster is no longer typed as a metal ligand. Every
  non-metal atom of a pure cluster was given a unique `Y*` type on the grounds
  that such a residue has no non-core atoms — true until a cluster could carry
  a hydroxo hydrogen, which bonds the oxygen rather than the metal.

- Session replay is strict again: a recorded answer is used only when it is the
  next unconsumed interaction and matches the prompt exactly. The forward scan
  that replaced it could not tell a recorded prompt the run skips apart from the
  same question asked at a different point in the workflow, and leapt for the
  latter — a hydrogen-editor prompt now asked during step 8 matched its
  recording from step 12 sixty-five interactions later, consuming three
  checklist decisions on the way, after which replay ran step 13 while the
  checklist sat at step 9. A prompt that does not match now announces the
  divergence once, falls through to live input, and leaves the position alone,
  so replay resynchronises as soon as the recorded question comes round again.
- Integer and decimal prompts are replayed by question rather than by position.
  `IntPrompt` and `FloatPrompt` are not subclasses of `Prompt`, so patching
  `Prompt.ask` never covered them: they fell through to the built-in input
  interception, where Rich has already printed the question itself, and every
  numeric answer was recorded with an empty prompt string. Any numeric question
  would then take the next numeric answer in the file, so adding one silently
  shifted the rest — a newly added atom-type offset prompt was answered with a
  Gaussian memory value recorded for something else. Existing logs record these
  answers anonymously and cannot be matched; those prompts now ask rather than
  supply a number from an unrelated question.
- A leftover `workflow_state.json` no longer blocks replay of a session log.

  The checklist asked "Resume from saved state?" whenever the file was present,
  a question the recorded run never faced and the log therefore cannot answer.
  Replay now starts such a workflow fresh instead of asking. Underneath, an
  unmatched prompt no longer consumes the rest of the log: the replayer scans
  forward to tolerate a recorded prompt the current run does not ask, but it
  now rewinds when nothing matches, so a single unexpected question costs only
  itself rather than ending the replay and dropping every later prompt to live
  input.
- Force Field Integration (checklist step 15) deposits one library entry per
  metal site instead of merging every site into a single entry. A merged entry
  was keyed by one site type, so it could not be reused on a structure carrying
  only one of the sites, and the reuse transformer it emitted matched neither
  site: transformer matching requires every residue in the rename table to be
  present in a *single* site, which a table spanning two sites never satisfies.
  Site type, redox state and spin state are now asked per site (the previous
  site's answers carry forward as defaults), and reusing one identity for two
  sites is refused rather than silently overwriting the first. Residue naming
  is unchanged: it still runs across all sites so names cannot collide.

- Metal-site parameterization now asks for the formal charge of a withheld
  cluster's non-metal core atoms (the bridging sulfides of an Fe-S cluster, the
  S/O core of a Mo cofactor), not only of its metals. A pure inorganic cluster
  is withheld from the force-field pass as a whole residue, so none of its atoms
  arrive with a charge; the core atoms were left unset and counted as zero, so
  the suggested QM charge came out short by their formal charge. An Fe2S2 site
  with four cysteinates proposed +2 where `[Fe2S2(SCys)4]2-` is -2. Correcting
  that by inflating a metal's charge was worse than it appeared, because the
  metal's formal charge is also the van der Waals radius key and is stored in
  the deposited library; a charge the radius database cannot resolve is now
  flagged as it is entered.

- The M\*/Y\* atom-type numbering no longer restarts at M1/Y1 when metal sites
  are parameterized in more than one pass, which would collide with types
  already deposited once every site's files load into one tLEaP session. The
  starting point comes from the workspace within a session and from earlier
  runs' fingerprint files in a fresh one; what was found is shown and can be
  overridden. Re-running a site offers to reuse the names it had before, so
  correcting one site replaces its entry instead of stranding its old types.
- Metal-site nuclearity counts metal atoms rather than metal-bearing residues,
  so an Fe2S2 cluster reads as binuclear instead of mononuclear.
- Metal-site and small-molecule counts in the menu suggestion and the status
  view count sites rather than residues; a site's coordinating ligand no longer
  inflates the total. The option numbers those messages cite now match the
  menu, which had gained an entry without them being updated.

## [1.15.0] — 2026-08-02

### Added

- Modified amino acid parameterization gained a from-structure route (Route B):
  a structure analyzer, a redesigned ten-step workflow with an explicit
  conformer-selection step, consistent residue naming, and resume resilience
  at every step. Torsional sampling is shared with Route A, which now
  auto-generates its sidechain scan rather than requiring one to be specified.
- Route B can optionally run a relaxed dihedral scan as a sampling mode and
  refit the scanned torsion with paramfit, so a rotatable bond that GAFF
  describes poorly can be corrected without leaving the workflow.
- Conjugate naming checks the force-field library for an existing residue of
  the same name before writing, and offers to reuse it instead of silently
  producing a second definition.
- The MD wizard offers the mdin keywords its own help text already referred
  to, adds `baro_stochastic` and `ninterface`, and completes four features
  that were previously only half-exposed. Coverage was swept against the full
  Amber manual by parameter index.
- Batch replay checks up front that every run supplies each variable the
  template needs, and fails loudly on divergence rather than blocking on
  stdin. The summary now distinguishes runs that failed from runs never
  attempted, and writes a retry list naming the runs still owed.
- The web shell tees each session to a plain-text transcript.
- The structure viewer can measure dihedrals, supports per-representation
  opacity, and defaults to an orthographic camera.
- The membrane builder offers solvate-only directly from the lipid menu.
- Force-field preparation deposits finished parameters from every
  parameterizer's final step, so the deposited library no longer depends on
  which route produced the parameters.
- An in-place updater, `update_proprep_in_ambertools.sh`, refreshes ProPrep
  inside an existing conda AmberTools environment without a full reinstall.
- Workspace inventory gained compact `--abbrev` output that collapses
  duplicate labels.

### Changed

- The MD manager treats the restraint manager as authoritative for imported
  mdin files, so restraints defined on import are no longer overwritten by
  the file's own restraint block.
- Heme transformers require their match criteria to be met exactly and clamp
  the bond credit awarded, so a partially matching site is no longer claimed
  by a transformer that does not fit it.
- The Seminario refinement scope reads "by-analogy" rather than "flagged",
  which describes what the option actually does.
- The membrane builder's "Advanced" geometry menu is now "Specialized
  Geometry".
- Component-type display names in the PDB filter come from a single source,
  so the triage table and the filter menu can no longer disagree.
- Checklist step numbers and section headings remain legible on light
  terminal backgrounds.
- The viewer's measurement pick marker is a fixed radius rather than scaling
  with the structure.

### Fixed

- Heme HMO (ferric) and HMR (ferrous) parameter sets had lost the trans
  pyrrole N-Fe-N angles; both are restored.
- The generated mdin carried an inverted comment for `ntmin`, and several MD
  wizard parameter descriptions were inaccurate. The help text was audited
  against the Amber manual across all 89 parameters and the advisory fields,
  which could not be kept correct, were dropped.
- Force-field parameter analysis no longer crashes on string-valued
  `*_structure` keys.
- Route B step 9 recovers the antechamber AC file on resume, and step 10
  reconstructs RedoxSite objects before syncing, so a resumed run no longer
  fails where a continuous one succeeds.
- The protonation summary surfaces desolvation, and the viewer no longer
  renders a stale structure after an edit.
- The installer no longer defaults past the orphan purge when run through
  `curl | bash`, self-heals on the update path, and forces a standalone
  ProPrep to win the `ambertools-dac` file clobber.
- Workspace inventory covers aliased receivers, wrappers, and comments, and
  excludes legacy metallo files that cannot be imported.

## [1.14.0] — 2026-07-09

### Added

- Triage category for pure inorganic metal clusters (for example 2Fe-2S,
  4Fe-4S, and Mo-S cores). A multi-atom residue that contains a metal but no
  carbon has no organic scaffold to hand to the small-molecule parameterizer,
  so it is now withheld from the standard force-field pass as a whole residue
  and its parameters are owned entirely by the metal-site (MCPB) workflow,
  instead of being misrouted to the organometallic-fragment path.
- Metal-site parameterization now assigns unique atom types to a cluster's
  internal coordinating atoms (bridging sulfides, a Mo-S-O core) rather than
  collapsing them onto a shared element type, so distinct metal-ligand bonds
  keep independent bonded parameters.
- Structure completeness now holds all resolved atoms fixed while MODELLER
  rebuilds a gap (loop refinement in a fixed environment), so gap filling
  cannot drag a metal-coordinating residue out of coordination. This makes the
  step safe to run before or after metal-site parameterization.
- Unfilled internal gaps can be closed with a TER record instead of a MODELLER
  fill or ACE/NME caps. Declining the repair plan, choosing TER on a
  single-residue gap, or running without MODELLER installed now inserts a TER
  record after the residue preceding each internal break, so tLEaP treats it as
  a chain end rather than building a long bond across the gap.
- The structure completeness step is now optional and can be skipped from the
  checklist with a `<num>s` command.
- Small-model single-residue gap filling offers a glycine-versus-actual
  bridging choice, matching the large model, so an incidental charged residue
  pulled into a gap between two coordinating ligands does not perturb the QM
  charge or Hessian.

### Changed

- The metal-cluster checklist step and its triage-table row are ordered before
  isolated metal ions.
- Single-residue gaps prompt to Fill or insert a TER record instead of being
  filled silently; internal single-residue gaps are never capped.
- Interactive repair prompts now list their choices explicitly (they were
  suppressed whenever human-readable option labels were supplied).
- Structure-alignment menus display and accept 1-based indices for reference
  structures and redox sites.

### Fixed

- Capping (rather than filling) an internal gap no longer crashes on undefined
  bracket-residue names; the gap-flanking residues are resolved through the
  residue mapper or a structure scan.
- Withheld-cluster bridging atoms with no assigned charge no longer trigger a
  NoneType error while writing the metal-site PDB.
- Cluster metals are no longer re-prompted for their element symbol, and the
  metal van der Waals lookup falls back to an element-based match when the
  residue name is not the element.
- The metal-site summary counts coordinating metal centers instead of the
  number of detected sites.

## [1.13.0] — 2026-07-08

### Added

- Transferable transformer framework. Redox/specialized-site transformers can
  now be authored as data-driven JSON specs instead of generated Python: roles
  are bound by element and atom connectivity, a split action
  (`move_to_new_residue`) is expressed through the id-mapper contract, and a
  JSON loader registers them. An interactive, table-driven creator is the
  user-facing way to author one, replacing the code generator, and
  user-authored transformers in `~/.proprep/transformers` are discovered and
  listed by the redox detector when it assigns site-type names.
- Membrane builder force-field parity with the Topology Generator for
  redox-active systems. The membrane force-field walkthrough reuses the
  Topology Generator's redox-site parameter-set selection and metal-ligand
  bond directives (stored in the shared workspace for reuse downstream),
  surfaces cofactor force-field prerequisites, and wires user-supplied solute
  frcmod/lib files through to the tLEaP topology stage.
- Modified-amino-acid from-structure multi-conformer parameterization. A
  disordered conjugated residue (e.g. a flavin in two crystallographic
  altlocs) is parameterized as a joint multi-conformer RESP fit — one capped
  model and QM optimization/ESP per altloc — instead of collapsing to a single
  conformer.
- Small-molecule parameterizer auto-emits a reuse transformer on deposit.
- Feedback utility encrypts the attached session context to the maintainer's
  key, so context can be shared in public issues while staying private by
  default.

### Changed

- 12-6-4 (Li/Merz) ion compatibility gate. ff19SB is now allowed — verified
  that AmberTools ships its renamed atom types in the polarizability file and
  that ParmEd applies the correction cleanly. ff15ipq stays gated, now for the
  correct stated reason: a modeling decision (the 12-6-4 sets were calibrated
  on ff14SB-family charges, not the implicitly-polarized IPQ charge model),
  not a technical failure.
- Menu and editor color legibility on both white and dark backgrounds (session
  editor tables, protonation-determinant headers, command-hint brackets).

### Fixed

- Fe4S4 averaged (default) RESP charge sets. Replaced a non-standard
  single-stage fit with the standard two-stage procedure (restoring the
  HB2=HB3 methylene equivalence) and re-converged a broken-symmetry SCF state
  that had collapsed and produced invalid charges. Bonded/LJ parameters
  unaffected.
- Session resume no longer crashes on templates with unfilled variables.

## [1.12.0] — 2026-07-04

### Added

- Route B modified-amino-acid from-structure parameterization. A modified
  amino acid can be parameterized directly from a capped model extracted from
  the input structure (QM optimization, ESP, and RESP charge fitting) rather
  than requiring a pre-built fragment.

## [1.11.0] — 2026-07-03

### Added

- Metal-site parameter **reuse without re-running MCPB**. A finished
  metal-site parameterization now deposits its full recipe into the user
  library at `~/.proprep` — the metal-bonded frcmod, each coordinating
  organic ligand's own parmchk2 GAFF frcmod, the OFF libraries, and an
  auto-emitted rename transformer that encodes the residue- and atom-name
  remaps the parameterization applied. Pointing the Redox Site Preparer at
  another instance of the same site replays that transformer and loads the
  deposited parameters directly, so an equivalent site elsewhere in the
  structure (or in a new structure) is parameterized by lookup instead of a
  fresh MCPB run.
- Auto-emitted reusable rename transformers. The parameterizer writes a
  data-only RedoxSite transformer describing the residue-name, atom-name, and
  membership edits it made. Matching is tiered (unique target, shared target,
  then a Weisfeiler-Lehman connectivity fingerprint on donor elements to
  disambiguate look-alike sites), self-diagnoses ambiguity, and falls back to
  the existing interactive manager when it genuinely cannot decide.
  User-authored transformers in `~/.proprep/transformers` are now discovered
  and registered alongside the built-ins.
- User force-field library promotion. Parameterized small molecules, modified
  amino acids, and metal sites can be promoted into `~/.proprep` (merge-safe
  upsert with round-trip validation and rollback) so users extend the built-in
  residue library rather than re-parameterizing.
- Multi-metal-site parameterization. A protein with two or more independent
  metal sites now parameterizes every site in one pass; metal (`M*`) and
  ligand (`Y*`) atom-type numbering is offset per site so shared-tLEaP-session
  types no longer collide, and Step 4 merges the fingerprints of all sites
  instead of only the first.
- Restrained metal–water (and other restrained) ligands. Metal-coordinating
  waters can be kept as flexible restrained contacts: a new bond `treatment`
  field round-trips through the redox-site model, DISANG restraints are
  pre-seeded for those contacts, and restrained-ligand terms are excluded from
  the bonded site frcmod.
- Per-redox-state propionate protonation in batch microstate generation. The
  batch path now prompts for heme-propionate pH treatment per (site, redox/spin
  state), matching the single-microstate path instead of silently defaulting.

### Changed

- Rigorous, calibrated MODELLER assessment for both gap-filling and homology
  modeling, replacing the earlier uncalibrated verdicts.
- Batch session replay now records the shell mode it was captured in
  (`proprep` vs `proprep-web`) and recreates that environment on replay, so
  mode-gated prompts no longer desynchronize; adds `--web`/`--no-web`
  overrides.
- The Topology Generator recommends a water model matching the selected
  protein force field and a divalent-ion set matching the selected water model.
- Redox bond and site prompts accept residue IDs (not only table row numbers)
  and auto-launch the 3D viewer at the bond-definition step; a multi-member
  redox site is presented as a single unit in the FF Parameterizer.
- Codebase-wide legibility sweep: Rich `dim` styling replaced with `grey50`,
  and named NGL representation colors resolved to hex so viewer representations
  render.

### Fixed

- MCPB frcmod finalizer no longer drops angle parameters whose atom types are a
  single character, and the metal element is emitted as `Mn` rather than `MN`.
- MCPB topology completeness: full intra-residue mol2 connectivity is derived
  by distance (fixing unbonded sidechains from truncated small models),
  ACE/NME cap-boundary bonded terms are enumerated in code, and partially
  `Y*`-substituted dihedral/improper terms are inherited.
- The metal-bond emitter honors a contact's restrained `treatment` instead of
  always emitting a rigid bond.
- Topology Generator no longer inherits or prompts for the provisional implicit
  solvent used during its own preprocessing tleap runs.
- MD restraint files (`DISANG`/`DUMPAVE`/`LISTOUT`) are emitted relative to the
  simulation directory to stay under the nmropt 80-character buffer, and
  constant-pH MD is now actually applied in the live setup path.
- Coordinated-water H–H bond artifact removed; capping a small model no longer
  severs an adjacent peptide bond.
- Dead code removed: the unused `_setup_workflows` implementation (~1060 lines),
  `structure_completeness_backup.py`, and other stale paths.

## [1.10.0] — 2026-06-26

### Added

- Session recordings now capture rich context for interactive prompts across
  the app. Roughly fifty prompts that previously recorded as bare `input()`
  lines now carry their module name, a human-readable description, and (for
  menu-style choices) an option-label map, so recorded interactions are
  legible and editable in the session editor rather than only replayable.
  Covers metal-site/ligand/antechamber parameterization, redox bond and site
  editors, alt-loc selection, MD restraint and AMBER workflow managers,
  biological assembly, and more. Classes that ran prompts without a processor
  reference now receive one from their caller so their prompts record context
  too.

### Changed

- Under `proprep-web`, the PDB Filter no longer asks "Would you like to view
  the structure in the interactive 3D viewer?" before chain selection. The
  3D viewer is already docked and open in the web shell, so the structure is
  shown there passively instead of prompting. The plain terminal still shows
  the opt-in prompt, since no viewer exists there yet.

### Fixed

- The `[exp]` experimental tag on the Redox Site Preparer's custom-transformer
  option now renders. It was being swallowed by the Rich markup parser as an
  unknown style span. The tag (and the module-level `[exp]` tags) now use
  `dark_orange3`, which stays legible on both white and black backgrounds —
  the previous bold yellow scored ~1:1 against white and was effectively
  invisible.

## [1.9.2] — 2026-06-23

### Changed

- Force Field Explorer selection is grouped, compact, and enriched. The curated
  force-field catalog is now shared with the Topology Generator (single source
  of truth in `proprep.forcefield_params.forcefield_catalog`), so the Explorer's
  scanned list shows friendly names, descriptions, a ★ recommended mark, and a
  ⚠ add-on flag, while still exposing every leaprc in AMBERHOME.

### Fixed

- Phosphorylated/modified amino-acid add-ons (`phosaa*`, `phosfb*`, `*_modAA`,
  `mimetic`, `fluorine`) are no longer grouped under "Protein". They have their
  own "Modified Amino Acids" category with a "needs a base protein FF" note, so
  selecting one alone (which leaves the whole protein untyped) is no longer an
  easy mistake.
- Force Field Explorer's smart default matched `ff14sb` as a substring and so
  defaulted to loading both `ff14SB` and `ff14SBonlysc`; it now matches the
  exact stem. Step-1 category hints/default are computed dynamically instead of
  assuming fixed category indices.

## [1.9.1] — 2026-06-23

### Fixed

- Disulfide bonds confirmed in the Redox Site Detector are now actually written
  to the generated tLEaP script. Three gaps are closed: (1) selecting "Finish"
  in site refinement now routes through the define/skip-bond prompt, so a
  disulfide site — which needs no further searching — can define its SG-SG bond;
  (2) interactively-defined SG-SG bonds are correctly classified as `disulfide`
  rather than `covalent` (the classifier was called without atom names); and
  (3) the single-state tLEaP writer now emits the `disulfide` bond category,
  which it previously omitted, silently dropping every disulfide bond.
- Bond editor: `disulfide` is now a first-class category across the filter
  (`d`), the category labels, and the add-bond menu. Also corrected the
  `metal_metal`-vs-`metal-metal` category-key mismatch that silently dropped
  metal-metal bonds from the written script and mislabeled them in the editor.

## [1.9.0] — 2026-06-22

### Added

- Per-module unmet-prerequisite notes in the menus: an unavailable (○) option
  now shows *why* it is blocked (e.g. "Needs a loaded structure") instead of a
  bare circle, and the module submenu banner shows the friendly reason rather
  than raw workspace keys.
- Module submenus reformatted to match the top-level menu: leading status glyph,
  bold-blue name, and an indented ⚠ reason line; distinct ●/✓/○/→ states.
- Session editor: delete/undo a range through the end of the session
  (`d 25-end` / `d 25-`).
- Alt-loc resolution viewer enhancement: a faint environment-shell overlay
  (residues within 5 Å) is drawn around each alternate-location residue.

### Changed

- Menu styling reworked to stay legible on both light and dark terminals
  (bold-blue names, foreground-following descriptions, dark-orange warnings,
  grey50 future stages), with two-line entries and bold-blue section headers.
- ONIOM full-residue QM/MM backbone boundary now uses a single-bond cut per
  side (formamide / primary-amide caps, flanking Cα left in the MM layer)
  instead of the acetamide/N-methylamide scheme. Cap parameters are the
  standard scaled hydrogen link atom; the formamide carbonyl C–H, absent from
  ff14SB, is supplied from GAFF.
- Menu suggestions now tell the user to "press [m] to return to the main menu"
  rather than the ambiguous "continue to next module".
- Session editor indices are now 1-based (matching the displayed numbers).
- Several module descriptions reworded; Structure Loader replaces PDB Loader in
  the full menu; the experimental ([exp]) tag removed from the now-tested
  ONIOM/ORCA/QM-MM/membrane modules.

### Fixed

- Menu availability: action options that require a loaded structure no longer
  show as available and then fail on selection — they show ○ with a reason.
- Menu completion: output-producing options now flip to ● once their result is
  written to the workspace (protonation, structure repair, alignment, topology
  load, restraint manager, force-field parameterizer).
- Redox-site restraint atom mask is now scoped to the redox-site residues
  (`(:res)&(@atoms)`); a bare name mask previously restrained same-named atoms
  (e.g. every Cα/Cβ) across the whole protein.
- The ONIOM next-step suggestion named the wrong module (PDB Filter → RedoxSite
  Detector) for detecting redox sites.
- Rich's auto-highlighter no longer recolors bracketed navigation labels
  (`[b]ack`) and syntax examples in help text.

## [1.8.1] — 2026-06-04

### Fixed

- Terminal `NME` caps now omit the methyl carbon, letting tLEaP build it from
  the loaded force field's template. Amber force fields disagree on the `NME`
  methyl carbon's atom name (`CH3` vs `C`); emitting only the amide `N` and `H`
  lets tLEaP rebuild the methyl group under the name the selected protein force
  field (ff14SB/ff19SB) expects, instead of leaving an explicitly placed `CH3`
  that may fail to match the template.

## [1.8.0] — 2026-06-04

### Added

- **Automatic ACE/NME terminal capping.** Structure completion can now cap
  chain termini with ACE (N-terminus) and NME (C-terminus) residues, so a
  fragment or isolated domain can be simulated without artificial charged
  termini. Capping works on MODELLER-renamed chains, strips the bonded
  residue's `OXT` when adding a C-terminal `NME`, and keeps cap residue
  numbers distinct from adjacent cofactors and ions.
- **PROPKA pKa-determinant inspection.** The Protonation State analyzer can
  now break a residue's predicted pKa down into its PROPKA determinants
  (desolvation, hydrogen bonds, charge–charge interactions) and highlight the
  contributing partner residues in the viewer — blue where a partner lowers
  the pKa, red where it raises it. Cofactors are folded into the calculation.

### Changed

- **The force-field-set picker is now the single source of the pH-treatment
  label.** When a cofactor ships both constant-pH and fixed-pH parameter sets,
  the Topology Generator filters the list to the treatment already chosen in
  the Redox Site Preparer, so you only choose among the remaining charge
  models. Every set now names its charge model in the title, the fixed-pH vs
  constant-pH label is applied consistently (previously only the fixed-pH sets
  announced it), and the His/Met-axial Fe-localized sensitivity set is flagged
  "not for production." Applies to all three heme leaves.

### Fixed

- Redox bond editor: when an atom name (e.g. `SG`) matches several residues,
  the "Select which '<atom>' atom(s)" prompt now accepts ranges (`1-4`),
  space/comma-separated lists, and `all` — previously a range like `1-4` was
  rejected as invalid even though the prompt displayed `1-4` as the bound.
- Terminal capping no longer emits collinear ACE/NME placeholder atoms, which
  previously triggered a PROPKA `ZeroDivisionError` on capped structures; the
  placeholder geometry is now trigonal-planar.

## [1.7.0] — 2026-06-02

### Added

- **Fixed-pH force-field sets for the c-type hemes.** The bis-histidine and
  His/Met-axial (cytochrome c2) c-type hemes now ship parallel `constpH` and
  `FixedpH` parameter sets, matching the Cys-axial b-type (NOS) heme. The
  `FixedpH` sets bundle static `PRP`/`PRD` propionates and load under modern
  ff14SB/ff19SB backbones with no `leaprc.constph`/`leaprc.conste` dependency;
  the `constpH` sets keep the titratable `PRN` propionate for constant-(pH,E)
  MD. Treatment naming is now consistent across all three heme leaves, and the
  Redox Site Preparer lets you pick each ring's protonation (A/D) under
  fixed-pH. The pH-treatment choice and its options now carry inline
  descriptions explaining the fixed-pH vs constant-pH trade-off.
- **Multi-structure spatial-comparison view.** Structure Loader option 5,
  "View local files together," opens two or more local PDB files overlaid in
  one viewer scene with their original relative coordinates preserved — for
  comparing geometries without superposition. View-only: it does not load the
  files into the preparation pipeline.

### Changed

- SLURM job scripts now default to `srun` (rather than `mpirun`) for MPI
  launches.

### Fixed

- Constant-pH heme builds now declare both `leaprc.constph` **and**
  `leaprc.conste` as prerequisites — the titratable `PRN` propionate is defined
  in `conste.lib`, not `constph.lib`, so the previous `constph`-only
  prerequisite could not resolve `PRN`.
- `--pdbview` (and local-file loading) no longer crash with `SameFileError`
  when the supplied filename is a symlink to the source file.
- MD Manager analysis browser: selecting only a topology (`.prmtop`/`.parm7`)
  file now prompts for the trajectory to pair it with, instead of dead-ending
  with "No files selected for analysis."
- Redox parameter configuration: an inapplicable (gated) option set — e.g. the
  per-ring protonation choice under constant-pH — is now shown as a calm note
  rather than a red error, and a stray debug line was removed.

## [1.6.1] — 2026-06-01

### Fixed

- MD Manager: opening **Setup and configure simulations** no longer crashes with
  `PermissionError` on a read-only install (e.g. a shared conda env). User MD
  templates now live in `~/.proprep/md_annotated_templates/user` instead of being
  created inside the package directory; built-in templates are still read from the
  package, and built-in regeneration tolerates a non-writable package dir.

## [1.6.0] — 2026-06-01

Constant-pH vs fixed-pH cofactor force-field sets, plus structure-repair and
interface-detection fixes.

### Added

- **Fixed-pH / constant-pH cofactor force-field sets.** A force-field set may now
  declare a `ph_treatment` (`constant_pH` or `fixed_pH`) and a `protonation_model`
  describing its titratable/protonation sites. Fixed-pH sets bundle static
  protomer residues (e.g. heme propionates PRD/PRP) directly in the library, so
  they load under a modern protein force field (ff14SB/ff19SB) with **no
  `leaprc.constph` dependency** — avoiding the parm10 backbone lock-in the
  constant-pH PRN model forces. The Redox Site Preparer exposes the choice
  (constant-pH vs fixed-pH, with independent per-site protomer selection), all
  output residue names are resolved from metadata by role (no hardcoded codes),
  and the Topology Generator honors per-set leaprc prerequisites. General
  infrastructure; the first consumer is the NOS Cys-axial b-type heme
  (`Guberman_H{TO,TR}_RESP_FixedpH`), shipped alongside the unchanged
  constant-pH default sets.

### Changed

- FF picker: the cofactor-prerequisite tag is reworded from "REQUIRED by your X
  selections" to "satisfies your X selections" — for an OR-group prerequisite
  (e.g. Zn(Cys)4 accepts any parm10-sourcing protein FF) several rows each
  satisfy the cofactor, so none is individually required.

### Fixed

- Structure Fixer: non-standard mutations and RedoxSite synchronization now map
  the chain id through MODELLER's chain renaming, so they work when MODELLER
  renumbers chains (e.g. a mutation on chain C that MODELLER renames to B).
- Structure Fixer: declining the repair plan returns cleanly to the menu instead
  of crashing with `'NoneType' object has no attribute 'needs_modeller'`.
- PDB Filter: SASA-based interface detection skips non-polymer (all-HETATM)
  chains — e.g. a cofactor-only chain — instead of aborting to distance-based
  detection for the entire structure.

## [1.5.0] — 2026-05-31

Coupling-aware ("effective") pKa in the PB-vs-ProPKA comparison.

### Added

- PB Titrate **effective (coupling-aware) pKa**: a new module
  (`effective_pka.py`) titrates the coupled system across a pH grid using the
  precomputed coupling matrix — Monte Carlo per pH point, no PB calls — and
  reads each site's effective pKa as the protonation-fraction 0.5 crossing.
  This is the coupled counterpart to the single-site (intrinsic) PB pKa and
  captures the neighbour-charge shift the single-site value can't see (e.g. a
  salt-bridge partner lowering a Tyr's effective pKa from ~8.0 to ~5.4).
- The "Compare pKa with ProPKA" step (pbt-cmp) now shows a **PB pKa (coupled)**
  column and computes its summary statistics against the **effective** PB pKa
  (with automatic fall-back to single-site when no coupling matrix exists). The
  comparison CSV gains effective-pKa and effective-Δ columns.

### Changed

- pbt-cmp comparison table highlights rows where PB and ProPKA **disagree on the
  protonation call** at the target pH (distinct colour from the existing
  large-|Δ| flag), and footnotes the statistics denominator (sites with a
  censored/locked effective pKa are excluded).

### Fixed

- pbt-cmp Δ column no longer shows a misleading single-site Δ for rows whose
  effective pKa is censored ('< x' / '> x'); such rows show '—' and their
  protonation call is resolved from the bound when it is decisive.

## [1.4.0] — 2026-05-31

PB Titrate coupled-solver workflow: targeted coupling, parallelization, solver
UX, and a set of correctness fixes in the redox/topology path.

### Added

- PB Titrate **targeted coupling**: when mean-field does not converge, restrict
  the coupling matrix + Monte Carlo to the unconverged ("flipper") sites plus
  their near neighbors, freezing the remainder at their converged state. Turns
  an all-sites N² coupling run into a small subgraph.
- PB Titrate **solver UX**: feasibility-aware recommendation with a graceful
  cluster-mean-field → Monte Carlo fallback when clusters are too large to
  enumerate; a per-site dominant-state + confidence table that flags ambiguous
  (near-50/50) sites; and a built-in multi-seed MC convergence check.
- PB Titrate ProPKA comparison now reports summary statistics (N, mean/median Δ,
  MAE, RMSD, Pearson r, % within 2 pH units, binary protonation-call agreement).
- Transformer assignments are persisted in the redox-sites JSON and recovered
  from disk by the Topology Generator, so they survive a ProPrep restart.

### Changed

- PB Titrate **pair-coupling PB calls now run in parallel** (the per-pair-cluster
  path was previously serial). Non-fatal pbsa exit-time-crash warnings are
  logged to a file with a per-iteration count instead of flooding the console.
  `coupling_summary` hides negligible (|W| ≈ 0) pairs.
- The `pbt-pdb` step strips hydrogens so tleap re-adds them per the rebuild
  force field (the PB topology's H come from a different FF).
- Menu gating: Redox Site Detector is available without a loaded PDB (for JSON
  import); Topology Generator accepts a loaded prmtop/rst7 pair (for PB Titrate).
- PB Titrate step order: "Write modern-FF PDB" now precedes "Persist
  Recommendations", and "Persist" no longer requires "Apply to Prmtop".

### Fixed

- Topology Generator crash (`'dict' object has no attribute 'site_id'`) when
  redox sites were loaded as dicts from a JSON or resumed session — sites are
  now normalized to RedoxSite objects.
- Mean-field captures the full flipper set across sweeps (period-2 safe) and
  merges coupled-solver assignments into the full site map instead of replacing
  it.
- Cofactor metadata: FAD/NAD/NADP `atom_types` include OD/O3 so tleap
  `addAtomTypes` covers the diphosphate.

## [1.3.0] — 2026-05-29

Usability and correctness fixes in the redox-site transformer workflow.

### Changed

- The per-site transformer-compatibility diagnostics shown during
  "Configure Site Transformers" are now written to
  `transformer_compatibility_report.txt` in the working directory instead
  of being printed in full to the screen. As the bundled transformer set
  has grown, this block could span hundreds of lines per run; the screen
  now shows a one-line pointer to the report. The on-demand
  `why <site#>` command still prints a single site's diagnosis inline.

### Fixed

- The redox-site workspace consistency check no longer reports spurious
  "Element mismatch" errors for two-letter elements. The element
  comparison is now case-insensitive, so a calcium center stored as `CA`
  is correctly matched against the PDB element symbol `Ca` (likewise for
  `Fe`, `Zn`, `Mg`, etc.). Single-letter elements were unaffected.

## [1.2.0] — 2026-05-29

Expansion of the bundled cofactor force-field library and a substantial
quality pass over the existing sets, integrating an externally-verified
parameter package. Adds the His/Met-axial c-type heme of cytochrome c2 with
four charge schemes for redox-Δq sensitivity studies; corrects silent
`ATOMIC_NUMBER` assignments that tleap had written from the atom type's
leading letter when many libraries were generated; restores the
"production sets coexist, zero conflicts" guarantee for the fragment-typed
organic cofactors after a per-cofactor QM rebuild; and rewrites every
active cofactor's methodology blurb in plain language sourced from the
package's own provenance docs.

### Added

- **His/Met-axial c-type heme** (cytochrome c2, *Blastochloris viridis*),
  ferric (HMO, +1 e) and ferrous (HMR, 0 e), low-spin, with a dedicated
  `heme_his_met_axial_c_type` transformer. Self-contained: each state's
  lib bundles the porphyrin plus the axial-His/Met side chains, the two
  thioether-Cys side chains, the corresponding backbone stubs (CYO, HIO,
  MEO), and the propionates (PRN), so only ff14SB need be loaded. Four
  charge schemes share a single bonded frcmod per redox state — CM5
  (default, density-based), pooled multi-conformer RESP (RESPmc),
  single-conformer RESP (RESPsc), and a Blumberger-style
  all-charge-on-iron limiting case (CM5_FeLocal, reduced-state endpoint
  of a redox-Δq localization sensitivity axis).
- **Cross-cofactor mean reconciliation** for shared structural-unit
  bond/angle terms across FAD/FMN/NAD/NADP after the per-cofactor
  Seminario rebuild. 45 differing terms — the ribitol/ribose–phosphate
  ester junction, the pyrophosphate bridge, and the shared
  isoalloxazine/ribose fragment terms — are set to their cross-cofactor
  mean (tagged `xcof-mean` in the frcmods); multi-cofactor topologies
  (FAD+FMN, FAD+NAD, flavocytochromes) now carry one consistent value
  per shared unit rather than the last-loaded cofactor's per-molecule
  QM value. Restores the zero-cross-cofactor-conflict guarantee.
- **`O3` terminal-phosphate type completion** in the FAD, NAD, and NADP
  frcmods. The lean rebuild forked the pyrophosphate bridge (`OD`) and
  terminals (`O3`) for bonded terms but omitted `O3`'s MASS and
  Lennard-Jones entries; these are folded in (`O3` LJ = `O2`'s, per
  Meagher 2003 *J. Comput. Chem.* **24**, 1016). No separate
  `frcmod.meagher_polyphosphate` file is required.
- **Newcomer-intelligible methodology blurbs** for every active cofactor,
  surfaced in the FF-picker `i<N>` panel. Each plain-language blurb is
  sourced from the verified package's own `METHODOLOGY.md` and per-leaf
  `PROVENANCE.json`: jargon defined in line, every term — RESP/CM5
  charges, value-neutral custom types, Seminario bond/angle from the QM
  Hessian (fragment-interior vs full-molecule-junction split),
  GAFF2/Walker/parm10-nucleic dihedrals, tleap-perceived impropers with
  the pyramidal-center drop, `O3`/`OD` phosphate forks, cross-cofactor
  averaging, validation — explained accessibly. Off-machine paths and
  external doc references stripped throughout; all provenance lives in
  the UX.

### Changed

- **`cys_axial_b_type`** rebuilt and integrated from the verified
  package: macrocycle H-name normalization for the `loadpdb`
  name-matching path, and a value-neutral atom-type fork to disjoint
  per-redox-state sets (`Ma`–`Me` ferric / `Mp`–`Mt` with `Mu`/`Mv`
  pyrrole N ferrous) that lets ferric and ferrous coexist in a
  mixed-redox topology without overwriting each other.
- **`cys_axial_b_type` spin-state directory** renamed `default →
  high_spin` (NOS-style heme is well-established high-spin in both
  redox states: ferric S=5/2, ferrous S=2); metadata key, transformer
  `spin_state` value, and directory layout updated to match the other
  heme sets' convention.
- **`bis_his_c_type`** trimmed to the verified Henriques set; the legacy
  `Guberman_HCO_RESP`/`Guberman_HCR_RESP` set, marked `is_valid:false`
  in 1.1 and with no recoverable source, is removed from the bundle.
- **Lean (delta-only) frcmods** for the flavin, nicotinamide, and pterin
  families, replacing the prior self-contained frcmods. Each lean
  frcmod contains only the parameters that the loaded base force fields
  (ff14SB/parm10 + GAFF2, plus optionally `leaprc.RNA.OL3` for
  nicotinamide adenine dihedral accuracy) do not already supply.
- **`bis_his_b_type`** parameter sets all marked `is_valid:false`
  pending re-parameterization (the reduced HBR set has a type-label
  mismatch and does not build; no recoverable source). The transformer,
  registration, and directory layout are intentionally retained so a
  re-parameterized set can drop in by replacing the lib/frcmod and
  flipping the flag.

### Fixed

- **Lib `ATOMIC_NUMBER` column** corrected across 19 bundled libs where
  tleap had silently written the element from the atom type's leading
  letter when the libs were generated — heme `FO`/`FR` → `F` = 9
  (fluorine), Fe4S4 `U*` → `U` = 92 (uranium), `Y*` → 39, `B*` → 5,
  `T*`/`A*`/`M*`/`L*` → −1. Resulting prmtops carried the wrong
  `%FLAG ATOMIC_NUMBER`, affecting downstream tools that key on element
  (GB radii, analysis, visualization).
- **`atomspertinfo` table `ptype` column** in the 15 libs that the lean
  rebuild had retyped without regenerating their `atomspertinfo`
  section (`cys_axial` + the 14 flavin/nicotinamide libs); `ptype` now
  mirrors the `atoms` table's `type` per atom, matching every other
  bundled lib. The `pelmnt = -1` convention is preserved (it is the
  standard tleap output).
- **`discover_forcefield_files` `is_valid` semantics** in
  `forcefield_params/loader.py`: when a cofactor's metadata explicitly
  defines `forcefield_sets`, the function now returns the surviving
  valid subset (even when empty) rather than falling through to
  directory-glob discovery — which ignored `is_valid` and re-exposed
  deliberately-disabled sets. Lets a cofactor be fully parked via
  `is_valid:false` while keeping its files in tree for a
  re-parameterization drop-in.

## [1.1.0] — 2026-05-24

First feature release since the AmberTools-26 snapshot was set at 1.0.0
(2026-03-10). Substantial additions across force-field bundling,
parameterization, topology generation, MD setup, and the browser-hosted
shell. The conda channel previously tracked this work by bumping the
1.0.0 build number from 1 to 42; that history is now collapsed into a
proper version bump.

### Added

- **Bundled cofactor force-field library.** Curated, internally consistent
  parameter sets for the recurring redox-cofactor inventory: flavins (FMN,
  FAD; oxidized, semiquinone, hydroquinone variants), nicotinamides
  (NAD⁺/NADH/NADP⁺/NADPH), the biopterin redox ladder (BIO, H₂B, H2Q,
  H3R, H4B, H4C), tetrahedral Zn(II)(Cys⁻)₄, cysteinate-axial b-type
  heme (NOS / cytochrome-P450 style), bis-histidine b-type and c-type
  hemes, and [4Fe-4S] clusters with broken-symmetry spin variants spanning
  the ferredoxin and high-potential redox couples. Fragment-typed naming
  convention so cofactors sharing substructure (FMN+FAD, FAD+NAD,
  pterin family) co-load without parameter-table contention.

- **FF compatibility matrix + interactive collision resolver.** Pre-flight
  compatibility matrix across all bundled parameter sets catalogs shared
  atom-type names and bonded-term disagreements. The FF picker consults
  the matrix at selection time and offers three resolution actions when
  a flagged pair is selected: pick a non-conflicting alternative,
  rename atom types in one set with `addAtomTypes` emission + lib/frcmod
  regeneration, or proceed with a logged warning.

- **PB Titrate.** Structure-specific pKa pipeline using AmberTools'
  `pbsa` solver. Eight-step interactive checklist accessible from the
  Topology Generator. Computes per-site intrinsic pKa values via the
  Bashford–Karplus four-state thermodynamic cycle and pairwise couplings
  via second-difference matrix construction; ships four solvers (mean
  field, single-site Metropolis, exact enumeration, cluster mean-field
  hybrid) with automatic recommendation. Multi-residue cofactor sites
  (heme propionates, [4Fe-4S]) handled via a Site envelope abstraction
  that keeps PB clusters integer-charge. Outputs feed the Topology
  Generator's `Use PB Titrate recommendations` option for cpin
  initial-state assignment.

- **proprep-web.** Browser-hosted shell with a docked NGL viewer next to
  the terminal pane. Detach-to-popup window, draggable splitter,
  optional install extra so users who prefer the pure-CLI experience
  carry zero extra dependency weight. Auto-port-scan from `--port`.

- **Viewer coordinator.** Process-wide singleton that drives the NGL
  viewer from intent-level methods (`show_structure`, `highlight`,
  `focus_on`, `show_sphere`, `show_bonds`). Mode-aware auto-launch:
  silent in CLI by default (push live updates if a viewer is already
  open), auto-launch in proprep-web on the first relevant event.
  Viewer hooks added across the Redox Site Detector, PDB Filter,
  Amino Acid Mutator, Structure Fixer (alt-loc migration), Redox Site
  Preparer, Protonation State Analyzer, Membrane Builder, Structure
  Orientator, MD Restraint Manager, and QM/MM Preparator.

- **Plugin infrastructure.** External packages register as ProPrep
  plugins via Python entry-points; tools are spliced into the workflow
  menu with their own stage. Workspace state hands off in-process to
  the plugin; the shared session recorder spans core + plugin modules.
  First plugin: *ETAnalyze* (manuscript in preparation).

- **Cluster profiles and run plans.** Three-object model for HPC
  submission: profiles describe a cluster's hardware palette (shareable);
  run plans bind a protocol to a cluster with per-step resource
  assignments. Bundled `tamu-aces.json` profile. SLURM generation moves
  to a dedicated Step 5 in the MD Manager.

- **Membrane Builder.** Lipid-bilayer setup integrated with the tLEaP
  topology generator (`packmol-memgen`-compatible). Full FF selection
  shared with tLEaP, anisotropic per-axis buffers, empty-bilayer
  (no-protein) support for pure-lipid simulations.

- **Topology Generator enhancements.** Per-axis buffer distances for
  rectangular boxes (DNA, membranes, anisotropic systems);
  multi-salt mixtures (`100 mM NaCl + 50 mM MgCl₂ + 5 mM CaCl₂` rather
  than NaCl only); arbitrary counter-ion selection; custom FF loading
  into the info pass; explicit-bond emission suppressed when libraries
  define `head`/`tail` atoms.

- **Final ParmEd validation pass.** Post-tLEaP `checkValidity` →
  `rediscoverMolecules` → `add12_6_4` runs automatically. Includes a
  polfile reconciliation step that resolves α conflicts when bundled
  cofactor force fields fork atom-type names that collapse into the
  same LJ slot in the assembled prmtop (e.g., heme β-pyrrole *Cp*
  sharing slot 7 with aromatic *CA*). Bonded-model metals auto-detected
  by bond count and skipped. Drops a `<stem>_parmed_replay.py` script
  next to the corrected topology so users can rerun the pass by hand.
  A new info panel orients users before the validation cascade fires.

- **MCPB cross-residue RESP equivalence.** At the RESP-fitting step of
  metal-site parameterization, users can declare cross-residue
  charge-equivalence groups using `<group>:<residue-list>` syntax
  (e.g., `1:1-4` to equate the four Cys ligands of a Zn(Cys)₄ site).
  Union-find merging composes within-residue (CH2/CH3) and
  cross-residue equivalences.

- **Site Transformers.** New transformers for the
  bis-histidine-ligated b-type heme (non-covalent ligation), the
  cysteinate-axial b-type heme (cytochrome-P450 style), the
  [4Fe-4S] cubane cluster (with three oxidation states for the
  ferredoxin and HiPIP couples and per-state broken-symmetry guess
  variants), and the tetrahedral Zn(II)(Cys⁻)₄ structural site. The
  Transformation Manager also acquired an ambiguous-role resolver
  (uniform mechanism for sites whose constituent residues cannot be
  matched purely from geometry — e.g., the front/rear axial histidines
  of a bis-his b-type heme).

- **Unified QM/MM Preparator.** ONIOM and ORCA preparators share frame
  extraction, layer assignment, and the redox-site picker. ONIOM
  output uses dedicated `HX` link-atom typing with explicit
  bond/angle/VDW parameter emission (after Gascón's t2ONIOM
  convention), and produces acetamide / N-methylamide caps rather
  than formamide when cutting at the Cα–C or N–Cα bonds.

- **Multi-structure MD setup.** Configure several structures in one
  MD-Manager session with per-structure or shared protocols.

- **Force Field Explorer module.** Browse all loaded FF parameters by
  atom type, residue name, or shared key. Custom FF file browser and
  `.prep`-file support. Bundled redox-param browsing + multi-select.

- **Session-replay backups.** Replay logs are backed up before in-place
  rewrites; select-N actions remap by filename so logs survive
  workspace reshuffles.

- **Conda recipe and `install.sh`.** ProPrep installable via
  `conda install mjgplab::proprep` against the `mjgplab` channel.
  Required dependency channels: `conda-forge`, `salilab` (MODELLER),
  `dacase` (`ambertools-dac`). The install script offers a clean
  upgrade path.

### Changed

- **Tool renames.** `TLeap` → `tLEaP` throughout the codebase, menus,
  and documentation, matching AMBER's canonical capitalization.

- **Basic Equilibration Protocol.** The single 4 ns NPT step split into
  two: a 1 ns CPU density-convergence phase (`02a_npt_density.mdin`)
  followed by a 3 ns GPU relaxation phase (`02b_npt_relaxation.mdin`).
  Was a workflow bottleneck on ~500k-atom systems because the GPU PME
  grid cannot reorganize while the box shrinks during density
  convergence; the split routes the first NPT to CPU and the second
  to GPU automatically.

- **SLURM script layout.** Adopts `stepN/` per-step subdirectories,
  matching the batch path. Scripts use `$SLURM_SUBMIT_DIR` for portable
  `cd` so bundles stay valid after rsync to a cluster.

- **SLURM custom commands split.** `custom_commands` replaced with
  separate `pre_commands` (emitted before the AMBER call) and
  `post_commands` (emitted after, regardless of AMBER exit code, with
  the exit code preserved for SLURM dependency chaining). Both
  surfaces in the cluster-profile JSON and the MD-Manager profile
  preview.

- **Hydrogen addition.** Always routes through tLEaP rather than
  `reduce`. `reduce` was unreliable on inputs with rare heteroatoms;
  tLEaP is the same engine used at topology generation.

- **ONIOM internals.** Refactored from coordinate-tuple keys to `parmed`
  integer indices. Atom types, partial charges, and FF parameters now
  read directly from the assembled prmtop rather than reconstructed
  from `AMBERHOME` files.

- **Cofactor force-field naming convention.** All bundled sets now use
  `{author}_{residue}_{charge_model}.{lib,frcmod}` so contributors can
  slot new sets into the appropriate cofactor directory and register
  them via a `metadata.json` declaring atom types, prerequisites, and
  citations. The transformer mechanism and compatibility matrix pick
  up new sets automatically.

- **FF prerequisite schema.** Cofactor metadata's `prerequisites` field
  supports OR'd `satisfied_by` lists. A cofactor that needs
  "any parm10-compatible protein backbone types" no longer marks each
  compatible protein FF as REQUIRED separately; the FF picker
  displays compatible standard FFs and surfaces the bundled
  parameter-set methodology before commit.

- **Viewer launch policy.** All CLI auto-launches gated behind explicit
  user action; the previous behavior (workflow waypoints firing
  `force=True` and popping browser tabs unbidden) is gone. proprep-web
  unaffected — its iframe still auto-fires via the `PROPREP_WEB_SHELL`
  short-circuit in the coordinator.

- **Topology Generator skip rule for membranes.** Empty-bilayer setup
  now skips STEP 2 reordering / TER-fix; protein-membrane systems
  treated identically (revisit pending when the first redox-protein-in-
  bilayer workflow comes up).

- **Conda dependency surface.** `freesasa` and `pdb2pqr` moved from
  pip to conda-forge; `tmtools` becomes an optional `tm-align` extra;
  MODELLER bundled with runtime license key from
  `KEY_MODELLER` env var.

### Fixed

- **tLEaP topology generation for MCPB metal sites.** Multi-site
  proteins no longer hit atom-type conflicts from independent
  per-site parameterizations.

- **`add12_6_4` "Zn0" crash on bonded-model metals.** Detection
  switched from atom-type-name heuristic (which keyed on `Zn0`,
  `Zn1`, ...) to bond-count check; bonded metals are now correctly
  skipped because their LJ already comes from the bundled frcmod.

- **Polfile shared-LJ-slot α conflicts.** ParmEd's `params1264`
  indexes by LJ slot; bundled forks that collapse into one slot
  (e.g., heme β-pyrrole *Cp* with aromatic *CA*) used to crash
  validation. The reconciliation pass overrides losing types to the
  canonical α and infers α for forked types not present in
  AmberTools' default polfile.

- **Multi-microstate MD-Manager handoff.** Workspace state for grouped
  microstates no longer leaks between siblings.

- **MD-Manager "Apply recommended assignments" crash.** Tuple-unpack
  bug on single-structure queues — the single-group branch was
  flat-unpacking 4 names from a nested 2-tuple — now matches the
  multi-group nested-unpack pattern.

- **Orphaned MD template.** Single-step `basic/02_npt_equilibration.mdin`
  removed; was a leftover from the Basic Equilibration split commit
  that produced a spurious "Basic (1 step): Metalloprotein NPT
  Equilibration" entry in the Template Catalog.

- **Rich-prompt edge cases.** `None` returns from interrupted prompts
  no longer crash callers; Unicode dashes accepted in range inputs;
  unescaped markup in two prompts fixed.

- **MDIN editor.** Launch uses the user's `$EDITOR` (falling back to
  `vi`/`vim` when default `nano` crashes on remotes); range selections
  and quoted-value parsing fixed.

- **Multi-chain protonation overrides.** Override dictionary keys now
  carry chain ID so `ASP15:A` and `ASP15:B` no longer collide.

- **MODELLER bundled errors.** C-level error output fully suppressed
  with `ctypes fflush`. Bundled data files alongside ProPrep so
  MODELLER works inside the conda package.

- **BioPython 1.86 compatibility.** Pinned `<1.86` until upstream PDB
  format-string regression is resolved (would otherwise produce
  79-char PDB lines).

- **Numerous cofactor-bundle correctness fixes.** Henriques c-type
  heme migrated to canonical conste residue layout; *Cp* atom-type
  fork closes Yang 2016 CB–CB override; axial His ligand identification
  via `site.atoms` instead of coordinate-derived bonds; cofactor
  transformer input matching restricted to PDB-CCD canonical
  resnames; Fe4S4 averaged libs declare head/tail and correct SG
  atom type.

### Removed

- **`InteractiveStructureViewer` bypass instantiations.** Three call
  sites (workflow_checklist, pdb_filter, small_molecule_parameterizer)
  no longer create phantom viewer instances; all route through the
  coordinator singleton so annotations from earlier tools remain
  visible.

- **Quality-assessment thresholds with fabricated numbers.** Removed
  arbitrary cutoffs from user-facing reports; quality is shown but
  not gated on invented heuristics.

---

## [1.0.0] — 2026-03-10

Initial release tagged for inclusion in AmberTools 26. Source state of
ProPrep at the version-string set point; see the ProPrep manuscript
draft for a tool-by-tool description of the included functionality.
