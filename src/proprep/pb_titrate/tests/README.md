# pb_titrate test suite

Tests are organized into two tiers based on what they need to run.

## Tier 1: no external calculator required

Run directly in any ProPrep environment.

| Test | What it covers |
|------|----------------|
| `envelope_retention/` | `Site` abstraction; cluster cutout invariants (envelope retention at tight radii, all-or-nothing rule for partial overlaps); `discover_sites` with and without `detected_redox_sites`. |
| `integration/` | Workflow surface (8 steps + handlers), Topology Generator menu wiring (`pb_titrate` slotted between `generate_topology` and `generate_cpin`), enhanced-menu status logic, cpin option 3 state-assignment logic with mocked prompts. |

```bash
cd src/proprep/pb_titrate/tests/envelope_retention && python run_test.py
cd src/proprep/pb_titrate/tests/integration         && python run_test.py
```

## Tier 2: requires AmberTools `pbsa`

These run the actual PB pKa pipeline. Need `pbsa` from AmberTools on `PATH`.

| Test | System | Time | Validates |
|------|--------|------|-----------|
| `bpti_asp3/` | BPTI ASP-3 single site | ~30s | Bashford 4-state cycle, εin sweep, NMR comparison (~3.0–3.4) |
| `bpti_e2e/` | BPTI 3-site full chain | ~30s | Mean-field iterate → state map → apply → prmtop |
| `bpti_multi/` | BPTI 3-site mean-field + MC + enumerate | ~minutes | Three solvers agree on the same data |
| `bpti_tyr10/` | BPTI TYR-10 | ~30s | TYR chemistry (different from carboxylates) |
| `hewl_his15/` | HEWL HIS-15 | ~minutes | 3-state HIS titration on a real protein |
| `multistate/` | ACE-AS4-NME, ACE-HIP-NME pure-solution | ~minutes | Multi-state Boltzmann (P=0.50/0.50 at pKa, HID:HIE=4:1) |
| `rnase_a/` | RNase A active-site, HIS network | ~minutes | Coupled HIS pair + active-site cluster |

```bash
cd src/proprep/pb_titrate/tests/bpti_asp3 && python run_test.py
# etc.
```

## Tier 3: requires user-specific data paths

These reference structures on the user's workhorse system.

| Test | Path it needs |
|------|---------------|
| `bundle_smoke/` | `/workhorse/9YUQ/SecondPass/MS001/transformed_microstate_001_fixed.prmtop` |
| `rebalance_unit/` | same bundle prmtop |

These will FileNotFoundError on any other machine. They are kept for the
9YUQ multi-heme demonstration but are not part of the portable test
suite. End-to-end demonstration of the `pb_titrate` pipeline on the 9YUQ
bundle (envelope-aware sites for hemes, integer cluster charge, cpin
recommendations propagating into the cpin command) is the qualitative
acceptance test for the integration; quantitative correctness is
inherited from titrate's own validation cascade (Tier 2).

## Adding a test

For Tier 1 (no PBSA): follow the patterns in `envelope_retention/` and
`integration/`. Use `MagicMock` and `unittest.mock.patch` to stub out
PBSA / parmed loads when they're not relevant to the assertion.

For Tier 2 (PBSA): mirror an existing per-system test. Each directory
holds its own fixture prmtop+rst7 (`leaprc.constph`-built) and a
`run_test.py` that loads them, runs PB, and asserts on the result.
