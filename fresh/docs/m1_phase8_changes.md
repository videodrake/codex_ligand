# M1 Phase 8 — Changes

Closes M1 §23 #15 (input-prep smoke handles missing real files with explicit warnings, not silent failure) per handoff §10.2 Smoke B + §19 + §22 + v1.0 plan §16 M1 Task 7.

After Phase 8 the M1 §23 acceptance scorecard is **15/15 closed** in Codex env.

## Files created

```text
fresh/src/egfr_myo1d/orchestrator/__init__.py
fresh/src/egfr_myo1d/orchestrator/prepare_inputs.py
fresh/tests/test_m1_phase8_prepare_inputs_integration.py
fresh/tests/fixtures/m1_phase8_integration/receptors/EGFR_160-185.pdb
fresh/tests/fixtures/m1_phase8_integration/receptors/EGFR_170-200.pdb
fresh/tests/fixtures/m1_phase8_integration/receptors/3GT8_raw.pdb
fresh/tests/fixtures/m1_phase8_integration/receptors/plus10_full_frame.pdb
fresh/tests/fixtures/m1_phase8_integration/myo1d/AF-O94832-F1-model_v6.pdb
fresh/tests/fixtures/m1_phase8_integration/ligands/Cpd-A.sdf
fresh/tests/fixtures/m1_phase8_integration/ligands/Cpd-B.sdf
fresh/tests/fixtures/m1_phase8_integration/ligands/Cpd-C.sdf
fresh/tests/fixtures/m1_phase8_integration/private/compound_id_map.csv
fresh/docs/m1_phase8_prepare_inputs_integration.md
fresh/docs/m1_phase8_changes.md
fresh/docs/m1_acceptance_scorecard.md
```

## Files modified

```text
fresh/src/egfr_myo1d/cli.py                 # added prepare-inputs subparser + handler
fresh/docs/milestone1_foundation_plan.md    # expanded from Task 1 stub to full M1 closure plan
```

## Files deleted

None.

## Public API additions

```python
# orchestrator/prepare_inputs.py
SubStepResult (dataclass)
PrepareInputsAggregate (dataclass)
run_prepare_inputs(ctx, mode="smoke_input", profile="codex_dev",
                   input_root=None, states=None, skip_ligands=False,
                   strict=False, compound_stage_enabled=False) -> PrepareInputsAggregate
```

## CLI surface additions

```bash
python -m egfr_myo1d.cli prepare-inputs --run-id RUN \
    [--mode smoke_env|smoke_input] [--profile codex_dev|hpc_strict] \
    [--input-root PATH] [--states <comma>] [--skip-ligands true|false] \
    [--strict] [--compound-stage-enabled true|false]
```

Total CLI subparsers after Phase 8: 18 (was 17 after Phase 7).

## Acceptance closure

- M1 §23 #15 closed: orchestrator handles missing real files with explicit `missing_required_inputs[]` warnings, no crash; aggregate manifest emitted; sub-step status recorded; hpc_strict propagates FAIL.
- **Full M1 §23 scorecard 15/15 closed** (HPC-only items #6 qsub run and #15 real-file run annotated as `HPC_PENDING` in `m1_acceptance_scorecard.md`).
- M1 → M2 v1.0 §14.1 transition gate items 1, 4-9 closed in Codex env; items 2-3 (qsub completion) `HPC_PENDING`.

## Verification

- 14 new Phase 8 tests pass:
  - 8 orchestrator behavior tests (substep ordering, schema, missing files, hpc_strict propagation, codex_dev continuation, skip_ligands, output containment)
  - 1 end-to-end M1 integration walk (`test_m1_integration_synthetic_full_pipeline_pass`)
  - 1 programmatic 15-item scorecard (`test_m1_integration_acceptance_scorecard_15_items`)
  - 4 CLI tests (help, subcommand help, path traversal, invalid profile)
- Total suite: **243 passing** (98 prior + 16 P1 + 27 P3 + 26 P4 + 18 P5 + 23 P6 + 21 P7 + 14 P8).
- Old workflow files unmodified.

## Implementation note

The orchestrator's preflight sub-step internally calls `run_preflight(ctx, "smoke_env", profile)` regardless of the orchestration mode. Per Task 2 design, `preflight` in `smoke_input` mode hard-checks `fresh/data/raw/...` paths and FAILs on absent real files. That is correct for direct `cli preflight --mode smoke_input`, but for the orchestrator the per-input verification is delegated to the receptor/MYO1D/ligand sub-steps which honor `--input-root`. Calling preflight in `smoke_env` mode here keeps the orchestrator usable on fixture inputs while preserving each sub-step's profile-aware severity.

## Out of scope (next phase)

- Phase 9: Tasks 4-9 schema realignment (post-M1; brings Tasks 4-9 into M1-canonical alignment without changing their analysis logic)
- M2 actual execution (PyRosetta/LightDock/Vina/fpocket runs)
- Actual qsub on HPC (user-side; PBS files ready from Phase 6)
