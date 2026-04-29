# M1 Phase 5 — Changes

Closes M1 §23 #12 (state-aware `membrane_frame.json` schema implemented and populated from coordinates) per handoff §15 and v1.0 plan §16 M1 Task 6 (membrane frame portion).

## Files created

```text
fresh/src/egfr_myo1d/model/membrane_frame.py
fresh/tests/test_m1_phase5_membrane_frame_generation.py
fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_dimer_with_TM_JM.pdb
fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_dimer_kinase_only.pdb
fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_3gt8_raw_kinase.pdb
fresh/docs/m1_phase5_membrane_frame_generation.md
fresh/docs/m1_phase5_changes.md
```

## Files modified

```text
fresh/src/egfr_myo1d/cli.py    # added compute-membrane-frame subparser + handler
```

## Files deleted

None.

## Public API additions

```python
# model/membrane_frame.py
TM_JM_RESIDUE_START = 634
TM_JM_RESIDUE_END = 674
PRIMARY_STATES, REFERENCE_CONTROL_STATES, ALL_STATES
MEMBRANE_FRAME_QC_COLUMNS  # 10 columns
COORDINATE_CONVENTION, FRAME_SOURCE_POLICY

StateMembraneFrame (dataclass)
compute_membrane_frame(state_full_frame_pdb, plus10_full_frame_pdb, state_id, profile) -> StateMembraneFrame
write_state_aware_membrane_frame_json(ctx, frames) -> Path
write_membrane_frame_qc_csv(ctx, frames) -> Path
run_membrane_frame_computation(ctx, state_ids=None, full_frame_source=None, profile="codex_dev") -> (frames, overall_status)
```

## CLI surface additions

```bash
python -m egfr_myo1d.cli compute-membrane-frame --run-id RUN \
    [--state EGFR_160-185|EGFR_170-200|3GT8_raw|all] \
    [--full-frame-source PATH] \
    [--profile codex_dev|hpc_strict] \
    [--mode smoke_env|smoke_input]
```

Total CLI subparsers after Phase 5: 15 (was 14 after Phase 4).

## Acceptance closure

- M1 §23 #12 closed: `manifest/membrane_frame.json` is state-aware (3 states), vectors computed from C-α coords via SVD principal-axis (no hardcoded fallback), 3GT8_raw marked `crystallographic_reference_control`, missing-source path produces clean WARN/FAIL with null vectors. Anti-hardcoded-vector test enforces no `[0, 0, 1]` / `[1, 0, 0]` / `[0, 0, -1]` literals in computation paths.

## Verification

- 18 new Phase 5 tests pass:
  - 8 per-state `compute_membrane_frame` cases (TM/JM, x-axis, centroids, plus10 fallback, 3GT8_raw, missing source codex_dev/hpc_strict, kinase-only no-TM)
  - 2 writer schema tests
  - 4 orchestrator tests (all-states, write containment, unknown state, invalid profile)
  - 1 anti-hardcoded-vector regression test
  - 3 CLI tests (help, subcommand help, path-traversal rejection)
- Total suite: 185 passing (98 prior + 16 Phase 1 + 27 Phase 3 + 26 Phase 4 + 18 Phase 5).
- Old workflow files unmodified.

## Out of scope (next phases)

- Phase 6: PBS generator
- Phase 7: ligand manifest shell
- Phase 8: prepare-inputs orchestrator + M1 integration test
- Phase 9: Tasks 4-9 schema realignment
- Alignment-derived 3GT8_raw frame derivation (deferred)
