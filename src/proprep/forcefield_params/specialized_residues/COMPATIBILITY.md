# Bundled FF Compatibility Matrix

`compatibility.json` records every conflicting pair of bundled FF sets — same atom-type name or same bonded-term key, with different parameter values. Loading two conflicting sets into one tleap session silently lets the last-loaded definition win, which corrupts the resulting prmtop.

See `docs/ff_collision_plan.md` for the design rationale.

## Regenerating

After adding or changing any bundled FF set:

```
python -m proprep.ff_compat.build_matrix --bundled-only
```

This rewrites `compatibility.json` in summary form (counts + shared types only, no per-entry values; the entry detail can be regenerated locally with `--detailed`). The summary is small enough to commit and to diff in code review.

To include user-derived FFs at `~/.proprep/forcefield_params/`, drop `--bundled-only`.

## What's a conflict, and what isn't

| Pattern | In matrix? | Comment |
|---|---|---|
| Two sets share an atom-type with **identical** values | No | Harmless duplication (very common across ff14SB-derived sets). |
| Two sets share an atom-type with **different** values | Yes | Real conflict — load order determines outcome. |
| Same residue name in two libs, same atom-name→type mapping, **same** charges | No | Harmless duplication. |
| Same residue name in two libs, same atom-name→type mapping, **different** charges | Not yet | Currently tracked separately; in scope for follow-up. |
| Same atom-type name, both declared in `addAtomTypes` with same element/hybridization | No | Identical declaration; tleap accepts. |

## Categories of expected conflicts

Most pairs in the current matrix fall into one of these classes — they ARE conflicts at the parameter level, but ProPrep's UI prevents the user from ever creating a topology that includes both.

- **Production sets across residues**: **zero conflicts by design** for fragment-typed cofactors (FAD/FMN, NAD/NADP, biopterin). The library's fragment-typed atoms (`Aa`..`Ay` for flavin, `Qa`..`Qt` for pterin, `Ka`..`Kt` for dihydronic NAD, `Da`..` Dy`/`Ga`..`Gy`/etc. for FAD/FMN states) are chosen so that cross-residue same-fragment-type bonds are bit-identical. Consolidator at build time enforces this. This is the multi-cofactor safety guarantee.
  - *Shared inter-cofactor structural units (cross-cofactor averaging, 2026-05).* The lean (delta-only) rebuild re-derived bond/angle force constants per-cofactor by Seminario, which made the **shared** structural units disagree slightly between cofactors — the ribitol→phosphate ester junction (`OS-c3`, `P-OS-c3`, `OS-c3-c3`, `OS-c3-h1`; FAD/FMN), the pyrophosphate bridge (`OD-P`, `P-OD-P`, `O3-P-OD`, `OD-P-OS`; FAD/NAD/NADP), and the shared isoalloxazine/ribose fragment terms (`Av-Aw-c3`, `CT-CT-Ka`, …). None are native to any base FF (e.g. mixed-case `OS-c3` ≠ gaff2's lowercase `c3-os` ≠ amber `CT-OS`), so they clobber nothing outside the cofactors — the only effect was cross-cofactor disagreement. Each such term is set to the **cross-cofactor mean** (tagged `xcof-mean` in the frcmods); the spread is within force-field accuracy (~1–6% K, ~1–2° angle) and no redox-active term is touched (those use distinct per-state types). This restores bit-identical shared units → zero cross-cofactor conflicts, with no new atom types and libs untouched.
- **Alternative redox states of the same cofactor** (e.g., `oxidized` vs `semiquinone_anionic` for a flavin). The picker assigns one redox state per detected site.
- **Alternative spin states of the same Fe4S4 redox state** (e.g., FO1 vs FO2). The Fe4S4 transformer uses one spin variant per detected cluster. The averaged sets (FO0/FR0/FD0) are alternatives at the same level.
- **Atom-type-namespace overlap between Fe4S4 state letters and c-heme MCPB types** (e.g., FD1's `T1`/`T2`/`T3` for Fe atoms vs c-heme's `T1`-`T5` for histidine N atoms). Distinct chemistries reusing the same 2-letter codes. Caught by the FF-compat resolver — the user picks which side to rename if both are loaded together.
- **His/Met-axial c-heme (`HMO`/`HMR`) vs the other heme/Fe4S4 sets**: the cytochrome-c2 His/Met c-heme shares the conste macrocycle types with the bis-His c-heme (one benign `MASS` declaration differs), the standard heme types with the Cys-axial b-heme (8 `BOND`/`ANGL`/`MASS` entries differ — the b-heme macrocycle is Yang-2016-derived, not conste), and the forked `FO`/`FR` Fe types with the Fe4S4 clusters (one `NONB` entry — the same overlap the Cys-axial set already has). A heme site resolves to exactly one set, so these co-load only in a multi-heme protein that mixes ligation types; the FF-compat resolver handles that with rename-on-collision, exactly as for the Fe4S4/c-heme overlap above. The `FO`/`FR` (per-redox-state Fe) and `Cp` (β-pyrrole) forks keep the His/Met set's own ferric and ferrous parameters collision-free in a mixed-state build.

## Historical note: Set 2 removal (2026-05)

The library previously shipped two parameter sets per cofactor: Set 1 (production, multi-cofactor safe) and Set 2 (whole-residue Hessian, single-cofactor research mode). Set 2 was removed entirely because it re-derived bonds and angles for every atom in the molecule — including shared standard types like `CT`, `OS`, `OH`, `N*`. By design it would silently overwrite parm10/gaff2 globals, and to be safe would have required forking every atom (which defeats the point of using standard types). Set 1 was renamed to the cofactor's residue name (e.g., the FAD oxidized set is now named `FAD`, the FAS semiquinone-neutral set is named `FAS`, etc.) to drop the now-redundant "Set 1" label.

## Verifying after a new set lands

A new bundled FF set must not introduce any NEW conflict pairs that weren't in the prior committed matrix. CI hook (TODO) will diff the regenerated matrix against the committed one and fail the build on any new pair.
