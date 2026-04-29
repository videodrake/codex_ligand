# M1 Phase 9 — Tasks 4-9 Schema Realignment

Phase 9 brings the M2-spec Tasks 4-9 layer (already shipped on `main` at commit 65de454) into alignment with the M1 canonical outputs produced by Phases 1-8. Phase 9 closes nothing new in the M1 §23 scorecard; M1 closure is already 15/15 after Phase 8.

## Approach: additive, not replacement

Phase 9 adds an **alignment record**, not a code-level replacement. Specifically:

- Task 4-9 module logic is **unchanged**. Their existing 79 passing tests are untouched.
- `validation/m1_alignment.py` adds a single helper, `record_m1_alignment(ctx)`, that scans `ctx.run_dir` for M1 canonical artifacts and writes `manifest/m1_alignment.json`. The manifest captures which M1 outputs are present and which Task 4-9 modules would naturally consume them in a future M2 actual-execution phase.
- New tests verify M1 outputs and Task 4-9 outputs coexist cleanly in a single run directory.

### Why additive?

1. **Preserve passing tests.** The 79 prior Task 4-9 tests assert specific output paths and behaviors. Replacing them speculatively would risk regressions.
2. **Avoid speculative refactor.** The Task 4-9 modules contain analysis logic that was carefully written during the original Codex pass: V924R reporting, terminal-artifact detection, active-face annotation, ATP-overlap mask, pose-acceptance policy, pocket-prioritization scoring. These should not be touched without a concrete M2 driver.
3. **Defer concrete schema bindings to M2.** When the actual M2 PyRosetta/LightDock/fpocket runners are implemented, they will be built to consume M1 normalized outputs natively. The exact CSV/JSON schema bindings should be made there, not pre-emptively.

The bigger code-level realignment (e.g. Task 4 stopping its pass-through `prepared/egfr/egfr_receptor_normalized.pdb` emission and instead referencing `normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb`) will happen organically during the M2 actual-execution implementation, when the runner code makes its concrete demands clear.

## What Phase 9 records

`manifest/m1_alignment.json` (per run):

```json
{
  "run_id": "...",
  "m1_outputs_present": {
    "membrane_frame_json": true,
    "myo1d_normalized_pdb": true,
    "myo1d_construct_qc_csv": true,
    "ligand_manifest_qc_csv": true,
    "normalized/receptors/EGFR_160-185_full_frame_explicit_AB.pdb": true,
    "normalized/receptors/EGFR_160-185_dockable_669_1014_explicit_AB.pdb": true,
    "normalized/receptors/EGFR_160-185_runtime_offset_receptor_only.pdb": true,
    "qc/EGFR_160-185_receptor_mapping.csv": true,
    "qc/EGFR_160-185_receptor_normalization_audit.csv": true,
    ...  // similar entries for EGFR_170-200, 3GT8_raw
  },
  "m1_canonical_paths": {
    "membrane_frame_json": "manifest/membrane_frame.json",
    ...
  },
  "task_consumption_map": {
    "task4_prepare_ppi_inputs": [...],
    "task6_plan_ppi_sampling": [...],
    "task7_summarize_ppi_consensus": [...],
    "task8_plan_pocket_discovery": [...],
    "task9_prioritize_pocket_candidates": [...]
  },
  "aligned_count": ...,
  "expected_count": ...,
  "status": "PASS|WARN",
  "notes": [...],
  "approach": "additive_phase9_no_task_module_logic_changed",
  "score_bonus_allowed": false,
  "timestamp": "..."
}
```

## Module additions

```text
fresh/src/egfr_myo1d/validation/m1_alignment.py    new (Phase 9)
```

Public API:

```python
M1AlignmentReport (dataclass)
M1_TO_TASK_CONSUMPTION  # dict mapping M1 artifact key -> (rel_path, consumer_task_ids)
PER_STATE_RECEPTOR_GLOBS
RECEPTOR_TASK_CONSUMERS
record_m1_alignment(ctx) -> M1AlignmentReport
```

## Tasks 4-9 modules (unchanged, but documented)

| Task | Module | Will consume M1 output (future M2) |
|---|---|---|
| Task 4 — `prepare-ppi-inputs` | `validation/prepared_inputs.py` | `normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb`, `normalized/myo1d/MYO1D_955_1006.pdb`, `manifest/membrane_frame.json` |
| Task 5 — `validate-real-inputs` | `validation/real_inputs.py` | All M1 normalized outputs as readiness check targets |
| Task 6 — `plan-ppi-sampling` | `validation/ppi_sampling_plan.py` | `normalized/receptors/<state>_runtime_offset_receptor_only.pdb`, MYO1D normalized |
| Task 7 — `summarize-ppi-consensus` | `validation/ppi_consensus.py` | `qc/<state>_receptor_mapping.csv` (for residue identity validation) |
| Task 8 — `plan-pocket-discovery` | `validation/pocket_discovery.py` | `normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb` |
| Task 9 — `prioritize-pocket-candidates` | `validation/pocket_candidate_prioritization.py` | `qc/<state>_receptor_mapping.csv` (for protomer_id resolution) |

Each Task 4-9 doc has been updated with a brief "M1 dependency" note linking to this realignment doc and `m1_alignment.json`.

## End-to-end coexistence

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"

# 1. M1 layer
python -m egfr_myo1d.cli init-run --mode smoke_input --run-id m1_realign_smoke
python -m egfr_myo1d.cli prepare-inputs \
    --run-id m1_realign_smoke \
    --input-root fresh/tests/fixtures/m1_phase8_integration

# 2. Task 4-9 layer (existing CLI commands; unchanged)
python -m egfr_myo1d.cli prepare-ppi-inputs \
    --run-id m1_realign_smoke \
    --input-root fresh/tests/fixtures/task4_inputs
python -m egfr_myo1d.cli plan-ppi-sampling \
    --run-id m1_realign_smoke \
    --input-root fresh/tests/fixtures/task4_inputs \
    --contract fresh/tests/fixtures/task4_inputs/ppi_input_contract.json
# ... and similarly summarize-ppi-consensus, plan-pocket-discovery,
# prioritize-pocket-candidates

# Both layers' outputs coexist in fresh/runs/m1_realign_smoke/.
# m1_alignment.json (if recorded) captures the consumption map.
```

`record_m1_alignment(ctx)` is a Python helper invoked from tests; no CLI subcommand is added in Phase 9 because the alignment record is more useful as a programmatic artifact than as a user-facing command.

## What's NOT in Phase 9

- No code-level changes to Task 4-9 modules
- No fixture path changes in Task 4-9 tests
- No CLI surface additions
- No changes to existing output paths
- No M2 actual execution work (PyRosetta/LightDock/Vina/fpocket runners)

## Tests added

11 new tests in `fresh/tests/test_m1_phase9_tasks_realignment.py`:

- 8 alignment-helper tests (status, manifest schema, per-state receptor recording, task consumption map, run-dir containment, phase status, empty-run-dir handling)
- 3 end-to-end coexistence tests (M1 + Task 4 in one run dir, M1 + Task 4 + Task 6 with alignment record, legacy Task 4 emissions preserved)

All 254 tests pass: 98 prior + 16 P1 + 27 P3 + 26 P4 + 18 P5 + 23 P6 + 21 P7 + 14 P8 + 11 P9.

## Out of scope (post-M1, post-Phase 9)

- M2 actual execution: PyRosetta PPI docking, fpocket pocket discovery, AutoDock Vina compound docking. These will be built as separate M2 phases per `egfr_myo1d_overall_implementation_plan_milestones_1_3_v1_0.md` §M2.1 — §M2.8.
- Any code-level changes to Task 4-9 modules to natively consume M1 outputs (deferred to the M2 runner implementations).
- Real EGFR/MYO1D/ligand PDB/SDF placement (user-side; placeholders ready under `fresh/data/raw/`).
- Actual `qsub` execution on HPC (user-side).
