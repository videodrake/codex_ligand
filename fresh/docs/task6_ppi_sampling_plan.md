# Task 6 PPI Sampling Plan

Task 6 adds a deterministic planning layer for future EGFR-MYO1D PPI sampling.
It converts a validated Task 4-style `ppi_input_contract.json` into auditable
job specifications, a pose-acceptance policy, manifests, and QC tables.

This stage is specification-only. It does not run PyRosetta, LightDock, Vina,
fpocket, P2Rank, qsub, docking, pocket discovery, compound docking, scoring, or
candidate nomination.

## CLI

Bash:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.cli init-run --mode smoke_input --run-id test_task6_local
python -m egfr_myo1d.cli plan-ppi-sampling \
  --run-id test_task6_local \
  --mode smoke_input \
  --profile codex_dev \
  --input-root fresh/tests/fixtures/task4_inputs \
  --contract fresh/tests/fixtures/task4_inputs/ppi_input_contract.json
```

PowerShell:

```powershell
$env:PYTHONPATH = "$PWD\fresh\src;$env:PYTHONPATH"
python -m egfr_myo1d.cli init-run --mode smoke_input --run-id test_task6_local
python -m egfr_myo1d.cli plan-ppi-sampling `
  --run-id test_task6_local `
  --mode smoke_input `
  --profile codex_dev `
  --input-root fresh/tests/fixtures/task4_inputs `
  --contract fresh/tests/fixtures/task4_inputs/ppi_input_contract.json
```

## Outputs

Task 6 writes only under `fresh/runs/<run_id>/`:

- `prepared/ppi/ppi_sampling_plan.json`
- `prepared/ppi/ppi_job_specs.jsonl`
- `prepared/ppi/ppi_job_specs.csv`
- `prepared/ppi/pose_acceptance_policy.json`
- `manifest/ppi_sampling_plan_manifest.json`
- `manifest/ppi_sampling_plan_report.json`
- `qc/ppi_sampling_plan_audit.csv`
- `qc/ppi_pose_qc_policy_audit.csv`
- `qc/ppi_sampling_blockers.csv`
- `reports/task6_ppi_sampling_plan_summary.md`

## Job Specs

Each job spec is a future engine placeholder with
`execution_mode=spec_only` and `execution_allowed=false`. The smoke plan uses:

- methods: `pyrosetta_global_ppi`, `lightdock_gso_ppi`
- protomers: validated EGFR chains `A` and `B`
- primary seeds: `0`, `1`
- comparator seeds: `0`

The primary construct is `MYO1D_sheet8_9_12_core_955_1001`. The
`MYO1D_ext_beta_meander_955_1006_tail_masked` construct is emitted only as a
`noise_monitor` comparator. The `962-1006` terminal-artifact fixture receives no
PPI jobs and is listed in blockers/quarantine QC output.

## Pose Acceptance Policy

`pose_acceptance_policy.json` records deterministic future pose classes:

- `accepted_active_8_9_supported_12`
- `accepted_active_8_9_only`
- `review_sheet12_dominant`
- `reject_tail_dominant_artifact`
- `reject_flipped_or_back_face`
- `reject_membrane_proximal`
- `reject_atp_site_overlap`
- `reject_chain_or_residue_reset`
- `reject_quarantined_input`

These thresholds are future QC rules, not current scoring. Active-face residues
`961-964` and `968-972`, sheet-12 support residues `993-997`, and tail/noise
residues `998-1006` are annotations only. `score_bonus_allowed` is always false
and `key_residue_bonus_weight` is always `0.0`.

## Strictness

`codex_dev` permits explicit fixture quarantines and comparator/noise-monitor
warnings so local smoke tests can run without production inputs or engines.
`hpc_strict` blocks production-primary use of ambiguous single-chain receptors,
V924R-like receptors, invalid membrane frames, and MYO1D terminal-artifact
primary constructs. Blocked inputs produce zero planned production jobs.

## Non-Goals

Task 6 does not implement docking, receptor mutation repair, receptor
renumbering, pocket discovery, compound docking, scoring, PBS/qsub generation,
or cleanup deletion.


---

## M1 dependency (added by Phase 9 alignment)

This task module was originally designed to run before the M1 foundation
modules (cleanup, receptor normalize, membrane frame, MYO1D construct, ligand
manifest) existed. After M1 Phases 1-8 the canonical input artifacts now live
under `fresh/runs/<run_id>/normalized/`, `fresh/runs/<run_id>/manifest/membrane_frame.json`,
and `fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv`.

Phase 9 takes an **additive** approach: this module's logic and output paths
are unchanged. The alignment between M1 outputs and this task's consumption
points is recorded by `validation/m1_alignment.py`'s `record_m1_alignment(ctx)`
helper, which writes `manifest/m1_alignment.json` listing which M1 artifacts
are present and which Task 4-9 modules would naturally consume them in a future
M2 actual-execution phase.

See `fresh/docs/m1_phase9_tasks4to9_realignment.md` for the full rationale.
