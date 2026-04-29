# Task 8 Pocket Discovery Planning

Task 8 adds a deterministic pocket-selection planning and intake layer. It uses
Task 7 `ppi_consensus_patch.csv` evidence to define guarded EGFR pocket-planning
records for future pocket discovery and focused compound docking. These rows are
not detected pocket records.

Task 8 does not run Vina, PyRosetta, LightDock, fpocket, P2Rank, AlphaFold,
Boltz, Chai, qsub, PBS, sbatch, pocket discovery, docking, ligand scoring, or
candidate nomination. It does not mutate EGFR and does not rewrite biological
residue numbering.

## CLI

Bash:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id test_task8_local
python -m egfr_myo1d.cli plan-pocket-discovery \
  --run-id test_task8_local \
  --mode smoke_env \
  --profile codex_dev \
  --input-root fresh/tests/fixtures/task8_pocket_planning
```

PowerShell:

```powershell
$env:PYTHONPATH = "$PWD\fresh\src;$env:PYTHONPATH"
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id test_task8_local
python -m egfr_myo1d.cli plan-pocket-discovery `
  --run-id test_task8_local `
  --mode smoke_env `
  --profile codex_dev `
  --input-root fresh/tests/fixtures/task8_pocket_planning
```

## Inputs

Default input is `ppi_consensus_patch.csv` under `--input-root`. The table must
match the Task 7 consensus patch schema and preserve:

- `ppi_patch_id`
- `receptor_id`
- `receptor_state`
- `protomer_id`
- EGFR residue numbering in `egfr_consensus_residues`
- Task 7 evidence class and warning fractions

## Outputs

All outputs are written under `fresh/runs/<run_id>/`:

- `pockets/egfr_myo1d_ppi_adjacent_pockets.csv`
- `pockets/egfr_myo1d_ppi_guided_pocket_plan_records.csv`
- `pockets/pocket_discovery_plan.json`
- `qc/task7_consensus_schema_audit.csv`
- `qc/pocket_selection_audit.csv`
- `qc/atp_overlap_audit.csv`
- `qc/membrane_accessibility_audit.csv`
- `qc/dimer_accessibility_audit.csv`
- `manifest/pocket_discovery_manifest.json`
- `reports/pocket_discovery_summary.md`

## Guardrails

Accepted Task 7 PPI patches are limited to cautious evidence classes such as
`CONVERGENT_PATCH` and `BROAD_PATCH`, and only when tail/terminal, ATP-overlap,
and membrane-proximal fractions are not blocking. If no accepted PPI patch is
available, Task 8 emits `NO_GO` with zero planned docking jobs.

ATP-overlap and membrane-proximal evidence remain visible in audits and are
blockers or strong warnings for PPI-disruptive objectives. Protomer identity and
EGFR residue numbering are copied from the Task 7 consensus table.

For backward compatibility, `egfr_myo1d_ppi_adjacent_pockets.csv` is still
written. The clearer alias `egfr_myo1d_ppi_guided_pocket_plan_records.csv`
contains the same rows and should be preferred in later tasks.

Planning rows include explicit semantics:

- `record_semantics = ppi_guided_pocket_plan`
- `pocket_detector_runtime_executed = false`
- `detected_pocket_record = false`
- `compound_docking_or_scoring_executed = false`
- `candidate_nomination_executed = false`

Task 8 creates pocket-selection planning records only. It prepares a handoff to
future pocket discovery, but does not claim real pocket discovery was run.

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
