# Task 2 Run Context, Manifests, Logging, and Preflight

Task 2 adds the operational foundation for fresh EGFR-MYO1D runs. It does not run docking, pocket discovery, receptor normalization, MYO1D slicing, qsub submission, cleanup deletion, scoring, or candidate nomination.

## Run Directory Layout

Each run is created under `fresh/runs/<run_id>/`:

```text
fresh/runs/<run_id>/
├── manifest/
│   ├── run_manifest.json
│   ├── input_manifest.json
│   ├── environment_report.json
│   └── git_snapshot.json
├── logs/
│   ├── master.log
│   ├── phase_status.jsonl
│   ├── job_status.jsonl
│   ├── jobs/
│   └── errors/
│       ├── error_summary.txt
│       └── failed_jobs.csv
├── qc/
├── reports/
├── scratch/
└── tmp/
```

`run_id` values must not contain `..`, `/`, or `\`. Outputs are restricted to the active run directory.

## PYTHONPATH

Bash:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
```

PowerShell:

```powershell
$env:PYTHONPATH = "$PWD\fresh\src;$env:PYTHONPATH"
```

## Commands

```bash
python -m egfr_myo1d.cli --help
python -m egfr_myo1d.cli version
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id test_task2_local
python -m egfr_myo1d.cli preflight --run-id test_task2_local --mode smoke_env --profile codex_dev
python -m egfr_myo1d.cli status --run-id test_task2_local
```

`init-run` creates the run tree, baseline logs, `run_manifest.json`, `input_manifest.json`, `git_snapshot.json`, and a placeholder `environment_report.json`.

`preflight` writes the real `environment_report.json` and appends a phase status record.

`status` prints manifest presence, last phase status, WARN/FAIL counts, and log paths.

## Modes

`smoke_env` checks the environment and records expected input presence, but missing receptor, MYO1D, and ligand files are warnings.

`smoke_input` is stricter for future input-prep runs: missing receptor and MYO1D files are failures, while ligand files remain warnings unless compound staging is explicitly enabled.

## Profiles

`codex_dev` is the default for local Codex smoke checks. Missing heavy scientific/HPC tools are warnings:

```text
pyrosetta, rdkit, BioPython, vina, fpocket, obabel, qsub, qstat, qdel
```

`hpc_strict` is for later HPC validation. Missing required scientific/HPC tools are failures. P2Rank/prank remains optional and warns when missing.

`numpy` and `pandas` are core Python dependencies for preflight reporting and are failures if missing.

## Logs and Manifests

Manifests live in `fresh/runs/<run_id>/manifest/`.

Central logs live in `fresh/runs/<run_id>/logs/`:

- `master.log`
- `phase_status.jsonl`
- `job_status.jsonl`
- `logs/errors/error_summary.txt`
- `logs/errors/failed_jobs.csv`

Task 2 initializes job-status logging helpers but does not submit jobs.
