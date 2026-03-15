# Config Guide

Last updated: 2026-03-12

This document explains the meaning of the configuration files in this directory. It is not the operator run sequence and it is not the command reference. Use [../docs/AI_START_HERE.md](../docs/AI_START_HERE.md) for onboarding order, [../docs/runbook.md](../docs/runbook.md) for execution procedure, and [../docs/manual_execution.md](../docs/manual_execution.md) for exact commands.

## Current Config Baseline

Read the files here with the current project baseline in mind:

| Topic | Current config interpretation |
|------|------|
| Receptor ensemble | Exactly three receptor states are the active baseline |
| Routine ligand workflow | Vina-centered |
| Phase 1 primary evidence | PyRosetta |
| Phase 1 secondary validation | LightDock |
| AFM fields | Legacy compatibility only unless explicitly re-enabled |
| Runtime environment | Production and pre-qsub lanes reuse the shared `pyrosetta` conda environment |
| Worker policy | `max_workers = 16` is the routine safe bound |

## What Lives In This Directory

| File group | Purpose |
|------|------|
| `example-project.yaml` | Main project-level YAML config for the current routine baseline |
| `full_test.yaml` | Auxiliary test-oriented YAML config, not the default baseline |
| `phase1/*.ini` | Phase 1 PyRosetta configs (18 files: 3 test + 15 production) |
| `run_pre_qsub_checks.pbs` | Scheduler wrapper for the lightweight precheck lane |
| `run_production.pbs` | Scheduler wrapper for the routine production lane |
| `run_production_fresh.pbs` | Fresh production rerun wrapper when prior production outputs must be discarded |
| `run_lightdock.pbs` | Scheduler wrapper for Phase 1 LightDock validation |
| `run_lightdock_test.pbs` | Scheduler wrapper for Phase 1 LightDock test submission |
| `run_full_test.pbs` | Auxiliary full-test submission helper |

## Config Surfaces And Their Roles

Three config surfaces coexist in the current repository.

| Surface | Scope | What it controls |
|------|------|------|
| YAML | Project-wide routine configuration | receptor list, ligand list, Vina settings, postprocess settings, worker policy, output root, registered PyRosetta result directories |
| INI | PyRosetta Phase 1 job configuration | receptor/partner metadata, construct metadata, chain mapping, run-specific output naming, PyRosetta job options |
| PBS | Submission-time wrapper configuration | environment activation, scheduler resources, precheck guard behavior, mode forwarding into production scripts |

Do not assume that the whole repository already shares one unified config schema.

## `example-project.yaml` Semantics

`example-project.yaml` is the main current project config and should be treated as the baseline semantic example for the routine workflow.

### Top-Level Meaning

| Field | Meaning |
|------|------|
| `project_name` | Names the output namespace for the routine project run |
| `output_root` | Root directory under which the project output tree is created |
| `mode` | Current Vina operating mode for the routine lane |
| `max_workers` | Parallel worker ceiling; treat `16` as the routine safe bound |
| `experimental` | Reserved area for non-baseline experiments; `null` means inactive |

### `receptors`

Each receptor entry defines one active receptor state for the routine baseline.

Expected meaning per entry:

- `id`: stable state identifier that should remain visible in downstream outputs
- `pdb`: structural input path used for receptor-side structural context
- `pdbqt`: docking-ready receptor path expected by Vina-oriented lanes
- `chain`: receptor chain identifier
- `source_type`: provenance label for how the receptor state was derived

Current baseline expectation:

- the receptor set remains exactly `3GT8_raw`, `3GT8_cl38_48`, and `3GT8_cl85_100`

### `ligands`

Each ligand entry registers one ligand across the routine ligand workflow.

Expected meaning per entry:

- `id`: stable ligand identifier
- `sdf`: chemistry input path
- `pdbqt`: docking-ready ligand path expected by Vina-oriented lanes

### `vina`

The `vina` section controls docking-search behavior for the routine ligand lane.

Representative semantics:

- `mode`: docking search style
- `exhaustiveness`: search depth
- `n_poses`: maximum poses retained per job
- `min_box` and `padding`: search box sizing controls

### `postprocess`

The `postprocess` section controls which downstream ligand analysis steps are enabled and how pocket logic is parameterized.

Representative semantics:

- parse and contact toggles determine whether pose parsing and residue-contact extraction run
- cluster and merge thresholds shape how raw contacts become pocket candidates
- comparison settings determine whether cross-state pocket comparison is produced
- report-oriented flags determine which summaries are emitted in the routine lane

Interpretation rule:

- these flags describe routine ligand postprocessing behavior; they do not define the advanced scientific Phase 2 to Phase 4 workflow by themselves

### `bootstrap`

The `bootstrap` section defines optional resampling behavior for robustness-style summaries.

Representative semantics:

- number of replicates
- sampling fraction
- random seed

### `ppi`

The `ppi` section is where the routine project config references receptor-side evidence sources.

Current meaning:

- `pyrosetta_result_dirs` registers existing PyRosetta result locations by receptor state and partner
- `afm_models` and `afm_settings` remain for legacy compatibility only

Important baseline rule:

- if `ppi.afm_models` is `null`, AFM is inactive in the routine baseline and should not be treated as a required evidence source

## `phase1/*.ini` Semantics

The Phase 1 INI files are the current PyRosetta job configs, located under `config/phase1/`. They do not replace the YAML project config; they support the still-separate PyRosetta execution surface.

Layout: 3 test configs + 15 production configs (5 seeds × 3 receptor states).

These files typically carry:

- receptor and partner identifiers
- construct metadata
- receptor and partner chain mapping
- numbering-system assumptions
- output directory naming metadata
- run-specific PyRosetta job settings

Interpretation rule:

- use INI files to understand one PyRosetta job's metadata and run shape, not to infer the whole project-level ligand workflow

## PBS Wrapper Semantics

The PBS files are submission wrappers, not the scientific source of truth.

### `run_pre_qsub_checks.pbs`

Semantic role:

- activates the baseline environment
- runs the lightweight validation lane before heavy production submission
- produces the precheck pass marker used by guarded production submission

### `run_production.pbs`

Semantic role:

- runs the routine production lane after environment setup
- forwards selected run modes into `run_production.py`
- expects the precheck guard unless explicitly bypassed

Important caution:

- production stage numbering is operational and should not be confused with the scientific Phase 1 to Phase 4 plan

### `run_production_fresh.pbs`

Semantic role:

- requests a clean rerun of the production lane when prior production outputs should no longer be trusted

### `run_lightdock.pbs` and `run_lightdock_test.pbs`

Semantic role:

- wrap the Phase 1 LightDock secondary validation entry points
- allow selection of specific receptor states at submission time

## Environment Semantics

The config and PBS files assume the shared `pyrosetta` conda environment for both pre-qsub and production lanes.

Current interpretation:

- a separate test-only conda environment is not the baseline
- environment setup should be understood as shared operational infrastructure, not as a per-lane divergence

## What Not To Infer From Config Alone

- Do not assume a field's presence means that lane is active in the current baseline.
- Do not treat AFM keys as proof that AFM is part of the routine workflow.
- Do not assume the production wrapper exposes the full scientific Phase 1 -> 2 -> 3 -> 4 architecture by default.
- Do not infer artifact meaning from config names alone; use the output docs for that.

## Use These Docs Next

- [../docs/manual_execution.md](../docs/manual_execution.md): exact commands and submission examples
- [../docs/runbook.md](../docs/runbook.md): operator sequence and stop/go rules
- [../docs/current_pipeline_status.md](../docs/current_pipeline_status.md): short summary of the current baseline
- [../docs/data_inventory.md](../docs/data_inventory.md): current input and output locations
