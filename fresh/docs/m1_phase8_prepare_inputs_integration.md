# M1 Phase 8 — prepare-inputs Orchestrator + M1 Integration

Closes **M1 §23 #15** (input-prep smoke handles missing real files with explicit warnings, not silent failure) per handoff §10.2 Smoke B + §19 + §22.

After Phase 8 the M1 §23 acceptance scorecard is **15/15 closed** in Codex env (HPC-only validation steps annotated `HPC_PENDING`).

## What it does

`fresh/src/egfr_myo1d/orchestrator/prepare_inputs.py` runs five M1 sub-steps under one `RunContext`, in order:

```text
1. preflight (env-only check; receptor/MYO1D/ligand presence is verified by the
                downstream sub-steps, not preflight, so the orchestration mode
                can be smoke_input without false-failing here)
2. prepare-receptor (per state from receptor_states.yaml: EGFR_160-185,
                      EGFR_170-200, 3GT8_raw)
3. compute-membrane-frame (state-aware membrane_frame.json + qc CSV)
4. prepare-myo1d (slice MYO1D source to 955-1006, emit normalized PDB + QC CSV)
5. manifest-ligands (Cpd-A/B/C public-only manifest; can be skipped via
                      --skip-ligands true)
```

Behavior:

- **`codex_dev` profile**: a sub-step FAIL is recorded; orchestration continues to the next sub-step. Aggregate `status` is FAIL if any sub-step FAILed.
- **`hpc_strict` profile** (or `--strict`): a sub-step FAIL stops the orchestrator immediately.
- **Missing real input files** in `smoke_input` are reported under `missing_required_inputs[]` without crashing the orchestrator.
- Aggregate manifest (`prepare_inputs_aggregate_manifest.json`) and a Markdown summary report (`prepare_inputs_summary.md`) are emitted under the run directory.

The `test_m1_integration_acceptance_scorecard_15_items` test programmatically verifies all 15 §23 items have artifact evidence after the synthetic integration walk.

## CLI

```bash
python -m egfr_myo1d.cli prepare-inputs \
    --run-id RUN \
    [--mode smoke_env|smoke_input] \
    [--profile codex_dev|hpc_strict] \
    [--input-root PATH] \
    [--states EGFR_160-185,EGFR_170-200,3GT8_raw] \
    [--skip-ligands true|false] \
    [--strict] \
    [--compound-stage-enabled true|false]
```

`--input-root` expects the layout:

```text
<input-root>/
├── receptors/
│   ├── EGFR_160-185.pdb
│   ├── EGFR_170-200.pdb
│   ├── 3GT8_raw.pdb
│   └── plus10_full_frame.pdb
├── myo1d/
│   └── AF-O94832-F1-model_v6.pdb
├── ligands/
│   ├── Cpd-A.sdf
│   ├── Cpd-B.sdf
│   └── Cpd-C.sdf
└── private/
    └── compound_id_map.csv
```

Defaults (when `--input-root` is omitted) come from `fresh/configs/{paths,receptor_states}.yaml`, pointing into `fresh/data/raw/`.

## Module additions

```text
fresh/src/egfr_myo1d/orchestrator/__init__.py
fresh/src/egfr_myo1d/orchestrator/prepare_inputs.py
```

Public API:

```python
SubStepResult (dataclass)
PrepareInputsAggregate (dataclass)
run_prepare_inputs(ctx, mode="smoke_input", profile="codex_dev",
                   input_root=None, states=None, skip_ligands=False,
                   strict=False, compound_stage_enabled=False) -> PrepareInputsAggregate
```

## Outputs

```text
fresh/runs/<run_id>/manifest/prepare_inputs_aggregate_manifest.json
fresh/runs/<run_id>/reports/prepare_inputs_summary.md
fresh/runs/<run_id>/logs/phase_status.jsonl    (appended)
fresh/runs/<run_id>/logs/master.log            (appended)

# Plus all sub-step outputs (manifest/, qc/, normalized/, scripts/) per
# Phases 1, 3, 4, 5, 7
```

## Severity / status

Aggregate status:

| Aggregate `status` | Conditions |
| --- | --- |
| `PASS` | every sub-step PASSes |
| `PASS_WITH_WARNINGS` | one or more WARN sub-steps; zero FAIL |
| `FAIL` | one or more sub-step FAIL |

Sub-step `SKIPPED` is used only for the ligand step when `--skip-ligands true`.

## Behavior policy (handoff §10.2 + §22)

- All sub-step outputs land under `fresh/runs/<run_id>/`. The aggregate manifest references each sub-step's manifest path.
- Path traversal `run_id` is rejected by `RunContext` validation.
- Missing source files in `codex_dev` produce `missing_required_inputs[<id>]` entries without crashing.
- `hpc_strict` (or `--strict`) propagates the first FAIL immediately.
- The orchestrator's preflight call uses `smoke_env` mode internally (env-only check); receptor/MYO1D/ligand presence is verified by the downstream sub-steps which honor `--input-root`.
- `score_bonus_allowed = false` is recorded in the aggregate manifest and every sub-step's manifest (no compound bonus from MYO1D key residues).

## Reusable fixture

`fresh/tests/fixtures/m1_phase8_integration/` is a composed layout that exercises the full pipeline:

```text
receptors/EGFR_160-185.pdb         (= mini explicit AB)
receptors/EGFR_170-200.pdb         (= same)
receptors/3GT8_raw.pdb             (= same)
receptors/plus10_full_frame.pdb    (= synthetic dimer with TM/JM)
myo1d/AF-O94832-F1-model_v6.pdb    (= MYO1D 955-1006 valid fixture)
ligands/Cpd-A.sdf, Cpd-B.sdf, Cpd-C.sdf
private/compound_id_map.csv       (synthetic INTERNAL_TEST_PLACEHOLDER_*)
```

## What is intentionally not in this phase

- Tasks 4-9 schema realignment (Phase 9)
- M2 actual execution (PyRosetta/LightDock/Vina/fpocket runs)
- Real EGFR/MYO1D/ligand placement (user-side, optional)
- Actual `qsub` run on HPC (user-side via Phase 6 scripts)

## End-to-end M1 closure walk

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"

python -m egfr_myo1d.cli init-run --mode smoke_input --run-id m1_closure_smoke
python -m egfr_myo1d.cli preflight --run-id m1_closure_smoke --mode smoke_env --profile codex_dev
python -m egfr_myo1d.cli prepare-inputs \
    --run-id m1_closure_smoke \
    --mode smoke_input \
    --profile codex_dev \
    --input-root fresh/tests/fixtures/m1_phase8_integration
python -m egfr_myo1d.cli prepare-pbs \
    --run-id m1_closure_smoke \
    --job-name m1_closure_smoke_env \
    --mode smoke_env --node node04
python -m egfr_myo1d.cli cleanup --run-id m1_closure_smoke --mode test --dry-run false
python -m egfr_myo1d.cli status --run-id m1_closure_smoke
```

After this walk all 15 §23 items have artifact evidence. See `fresh/docs/m1_acceptance_scorecard.md`.
