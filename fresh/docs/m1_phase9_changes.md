# M1 Phase 9 — Changes

Phase 9 brings Tasks 4-9 (M2-spec layer shipped on `main` at 65de454) into alignment with M1 canonical outputs from Phases 1-8. Closes nothing new in M1 §23 — M1 closure is already 15/15 after Phase 8.

## Approach

**Additive, not replacement.** Phase 9 does NOT modify Task 4-9 module logic; the 79 prior Task 4-9 tests are unchanged. Instead it adds an alignment helper that records the cross-reference between M1 outputs and Task 4-9 modules.

See `fresh/docs/m1_phase9_tasks4to9_realignment.md` for rationale.

## Files created

```text
fresh/src/egfr_myo1d/validation/m1_alignment.py
fresh/tests/test_m1_phase9_tasks_realignment.py
fresh/docs/m1_phase9_tasks4to9_realignment.md
fresh/docs/m1_phase9_changes.md
```

## Files modified

```text
(no Task 4-9 module logic modifications in Phase 9 by design)
```

## Files deleted

None.

## Public API additions

```python
# validation/m1_alignment.py
M1AlignmentReport (dataclass)
M1_TO_TASK_CONSUMPTION  # dict mapping artifact key -> (rel_path, consumer_task_ids)
PER_STATE_RECEPTOR_GLOBS
RECEPTOR_TASK_CONSUMERS
record_m1_alignment(ctx) -> M1AlignmentReport
```

## CLI surface additions

None in Phase 9. The alignment record is a programmatic artifact (test helper), not a user-facing CLI command. Any future surface change should be driven by M2 actual execution.

## Acceptance closure

- M1 §23 acceptance: unchanged at 15/15 closed (set in Phase 8).
- Phase 9 enforces the additive policy in `manifest/m1_alignment.json` via the
  `approach` field set to `"additive_phase9_no_task_module_logic_changed"`.

## Verification

- 11 new Phase 9 tests pass:
  - 8 alignment-helper unit tests
  - 3 end-to-end coexistence tests (M1 + Task 4 / Task 6 in one run dir;
    legacy Task 4 outputs preserved)
- Total suite: **254 passing** (98 prior + 16 P1 + 27 P3 + 26 P4 + 18 P5 +
  23 P6 + 21 P7 + 14 P8 + 11 P9).
- Old workflow files unmodified.
- Task 4-9 module files have zero diff in Phase 9 commit.

## What is intentionally NOT in this phase

- Code-level changes to Task 4-9 modules (deferred to M2 runner work)
- Output path changes for any existing Task 4-9 emission
- Modifications to Task 4-9 test files
- New CLI subcommands
- M2 actual execution

## Next: post-M1 / Milestone 2

```text
M2.1  PPI input generation (real PyRosetta-ready inputs from M1 normalized outputs)
M2.2  PyRosetta adapter (smoke -> mini -> production scale ladder)
M2.3  PPI pose QC + MYO1D artifact filtering
M2.4  Symmetry-aware consensus patch
M2.5  ATP-site reference
M2.6  fpocket pocket discovery
M2.7  Membrane / dimer / PPI pocket gates
M2.8  Milestone 2 aggregation
```

When M2.1+ runners are implemented, they will consume M1 normalized outputs
natively. At that point, Task 4-9 spec modules will be reviewed for whether
they should be (a) refactored to share logic with the M2 runners, (b) kept
as audit-only spec layer, or (c) wholly replaced. That decision is a M2
implementation concern, not a Phase 9 concern.
