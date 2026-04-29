# Task 9 Pocket Candidate Prioritization

Task 9 ingests provided detector-style EGFR pocket records and prioritizes them
against Task 8 PPI-guided pocket-planning evidence. It is an intake and QC layer
only.

Task 9 does not run Vina, PyRosetta, LightDock, fpocket, P2Rank, qsub/PBS,
compound docking, compound scoring, ligand preparation, EGFR mutation repair, or
candidate nomination.

## CLI

Bash:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id test_task9_local
python -m egfr_myo1d.cli prioritize-pocket-candidates \
  --run-id test_task9_local \
  --mode smoke_env \
  --profile codex_dev \
  --input-root fresh/tests/fixtures/task9_pocket_candidates \
  --pocket-plan fresh/tests/fixtures/task9_pocket_candidates/pocket_discovery_plan.json \
  --candidate-pockets fresh/tests/fixtures/task9_pocket_candidates/pocket_detector_candidates.csv
```

PowerShell:

```powershell
$env:PYTHONPATH = "$PWD\fresh\src;$env:PYTHONPATH"
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id test_task9_local
python -m egfr_myo1d.cli prioritize-pocket-candidates `
  --run-id test_task9_local `
  --mode smoke_env `
  --profile codex_dev `
  --input-root fresh/tests/fixtures/task9_pocket_candidates `
  --pocket-plan fresh/tests/fixtures/task9_pocket_candidates/pocket_discovery_plan.json `
  --candidate-pockets fresh/tests/fixtures/task9_pocket_candidates/pocket_detector_candidates.csv
```

## Inputs

Task 9 consumes:

- a Task 8 `pocket_discovery_plan.json`
- a provided detector-style `pocket_detector_candidates.csv`
- optional Task 7 `ppi_consensus_patch.csv` provenance

Detector-style records are supplied by fixtures or user input. Task 9 records
their provenance but does not execute detector runtimes.

## Prioritization

PPI relationship uses residue overlap and nearest-residue distance before
centroid distance. Centroid-only proximity is not enough to carry a pocket
forward.

ATP-overlap pockets are blockers for this EGFR-MYO1D PPI-disruptive objective.
Membrane-blocked pockets and dimer-buried or protomer-ambiguous pockets are also
blocked or quarantined. Detector druggability is detector evidence only; it is
not compound docking evidence.

Priority classes are:

- `CARRY_FORWARD_STRONG`
- `CARRY_FORWARD_WEAK`
- `EXPLORATORY_ONLY`
- `BLOCKED_ATP`
- `BLOCKED_MEMBRANE`
- `BLOCKED_DIMER`
- `BLOCKED_SCHEMA`

No final compound or pocket candidate is nominated in Task 9.

## Outputs

All outputs are written under `fresh/runs/<run_id>/`:

- `pockets/egfr_myo1d_prioritized_pocket_candidates.csv`
- `pockets/egfr_myo1d_pocket_candidate_families.csv`
- `pockets/pocket_candidate_prioritization.json`
- `qc/pocket_candidate_schema_audit.csv`
- `qc/pocket_candidate_input_audit.csv`
- `qc/ppi_pocket_proximity_audit.csv`
- `qc/atp_exclusion_audit.csv`
- `qc/membrane_compatibility_audit.csv`
- `qc/dimer_accessibility_audit.csv`
- `qc/pocket_candidate_blockers.csv`
- `manifest/pocket_candidate_prioritization_manifest.json`
- `reports/task9_pocket_candidate_prioritization_summary.md`

The JSON and manifest explicitly state:

- `planned_compound_docking_jobs = 0`
- `compound_docking_or_scoring_executed = false`
- `candidate_nomination_executed = false`
- `pocket_detector_runtime_executed = false`

