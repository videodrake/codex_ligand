# Milestone 1 Foundation Plan

Milestone 1 is the foundation layer of the EGFR-MYO1D fresh workflow. After completion of the M1 rework on branch `claude/task10`, the M1 §23 acceptance scorecard is **15/15 closed** (HPC-only validation steps annotated `HPC_PENDING`).

This document supersedes the original Task 1 stub.

## Purpose

M1 establishes the operational foundation required before any docking/pocket discovery / scoring work. It exists to prevent the historical failure modes catalogued in `codex_overall_project_context_v0_1.md` §13:

1. EGFR dimer stored as duplicate chain X and misread as monomer
2. PyRosetta pose manipulation resetting residue numbering
3. HETATM caps (ACE/NME) dropped by ATOM-only parsers
4. MYO1D terminal artifact (`962-`-start)
5. MYO1D beta-sheet face-flip
6. ATP pocket dominance in scoring
7. Naïve PPI/ligand centroid comparison
8. Logs and intermediate outputs scattered

M1 directly closes failure modes 1, 2, 3, 8, and embeds policies for 4-7 (which are enforced in M2 sampling/QC).

## Source-of-truth specs

```text
milestone1_foundation_codex_handoff_v0_5.md      §4 module tree
                                                   §10 PBS generation
                                                   §14 receptor normalization
                                                   §15 membrane frame
                                                   §16 MYO1D construct QC
                                                   §17 ligand manifest
                                                   §18 cleanup
                                                   §19 CLI commands
                                                   §22 HPC validation boundary
                                                   §23 acceptance scorecard (15 items)

egfr_myo1d_overall_implementation_plan_milestones_1_3_v1_0.md
                                                   §16 M1 task sequence
                                                   §14.1 M1→M2 transition gate
```

## Module tree (post-M1, matches handoff §4)

```text
fresh/src/egfr_myo1d/
├── analysis/                  (Tasks 7-9, KEEP)
├── core/
│   ├── run_context.py
│   ├── manifest.py
│   ├── logging_utils.py
│   └── cleanup.py             ✅ Phase 1
├── hpc/
│   ├── __init__.py
│   └── pbs.py                 ✅ Phase 6
├── io/
│   ├── hashing.py
│   └── residue_mapping.py     ✅ Phase 4
├── ligand/
│   ├── __init__.py
│   └── manifest.py            ✅ Phase 7
├── model/
│   ├── __init__.py
│   ├── receptor_normalize.py  ✅ Phase 4
│   ├── receptor_qc.py         ✅ Phase 4
│   └── membrane_frame.py      ✅ Phase 5
├── myo1d/
│   ├── __init__.py            ✅ Phase 2
│   ├── construct.py           ✅ Phase 2 + Phase 3
│   ├── pdb_writer.py          ✅ Phase 2
│   └── qc.py                  ✅ Phase 3
├── orchestrator/
│   ├── __init__.py            ✅ Phase 8
│   └── prepare_inputs.py      ✅ Phase 8
├── planning/                  (Task 6, KEEP)
├── preparation/               (EGFR-side masks/restraints; MYO1D moved out in P2)
├── structure/                 (Task 3)
├── validation/                (Tasks 3-9 + preflight)
├── cli.py
└── __init__.py
```

`fresh/scripts/` placeholders (`cleanup_run.py`, `generate_pbs.py`, `submit_smoke_env.sh`, `submit_smoke_input.sh`) replaced with real implementations in Phases 1 and 6.

## CLI command reference (post-M1)

```bash
# Foundation (Tasks 1-2; existing)
python -m egfr_myo1d.cli version
python -m egfr_myo1d.cli init-run    --mode smoke_env --run-id RUN
python -m egfr_myo1d.cli preflight   --run-id RUN --mode smoke_env|smoke_input --profile codex_dev|hpc_strict
python -m egfr_myo1d.cli status      --run-id RUN

# Tasks 3-9 (existing M2-spec layer)
python -m egfr_myo1d.cli validate-structures            --run-id RUN ...
python -m egfr_myo1d.cli prepare-ppi-inputs             --run-id RUN ...
python -m egfr_myo1d.cli validate-real-inputs           --run-id RUN ...
python -m egfr_myo1d.cli plan-ppi-sampling              --run-id RUN ...
python -m egfr_myo1d.cli summarize-ppi-consensus        --run-id RUN ...
python -m egfr_myo1d.cli plan-pocket-discovery          --run-id RUN ...
python -m egfr_myo1d.cli prioritize-pocket-candidates   --run-id RUN ...

# M1 completion (Phases 1, 3-8)
python -m egfr_myo1d.cli cleanup                  --run-id RUN --mode test|production [--dry-run true|false]   (P1)
python -m egfr_myo1d.cli prepare-myo1d            --run-id RUN --source PATH [--construct 955-1006]            (P3)
python -m egfr_myo1d.cli prepare-receptor         --run-id RUN --state STATE --source PATH                     (P4)
python -m egfr_myo1d.cli compute-membrane-frame   --run-id RUN [--state all] [--full-frame-source PATH]        (P5)
python -m egfr_myo1d.cli prepare-pbs              --run-id RUN --job-name NAME --mode <m> [--node NODE]        (P6)
python -m egfr_myo1d.cli manifest-ligands         --run-id RUN [--ligands-dir PATH]                            (P7)
python -m egfr_myo1d.cli prepare-inputs           --run-id RUN [--input-root PATH]                             (P8)
```

Total: 18 subcommands.

## Run output schema (post-M1)

```text
fresh/runs/<run_id>/
├── manifest/
│   ├── run_manifest.json
│   ├── input_manifest.json
│   ├── environment_report.json
│   ├── git_snapshot.json
│   ├── cleanup_report.json                                       (Phase 1)
│   ├── myo1d_construct_manifest.json                             (Phase 3)
│   ├── <state>_receptor_manifest.json                            (Phase 4 per state)
│   ├── membrane_frame.json                                       (Phase 5)
│   ├── ligand_manifest_report.json                               (Phase 7)
│   ├── prepare_inputs_aggregate_manifest.json                    (Phase 8)
│   └── <task>_manifest.json                                       (existing M2 spec)
├── normalized/
│   ├── receptors/
│   │   ├── <state>_full_frame_explicit_AB.pdb                    (Phase 4)
│   │   ├── <state>_dockable_669_1014_explicit_AB.pdb             (Phase 4 primary states)
│   │   └── <state>_runtime_offset_receptor_only.pdb              (Phase 4 +1000 offset)
│   └── myo1d/
│       └── MYO1D_955_1006.pdb                                    (Phase 3)
├── qc/
│   ├── myo1d_construct_qc.csv                                    (Phase 3)
│   ├── <state>_receptor_mapping.csv                              (Phase 4)
│   ├── <state>_receptor_normalization_audit.csv                  (Phase 4)
│   ├── membrane_frame_qc.csv                                     (Phase 5)
│   ├── ligand_manifest_qc.csv                                    (Phase 7)
│   └── <task>_*.csv                                               (existing M2 spec)
├── reports/
│   └── prepare_inputs_summary.md                                 (Phase 8)
├── scripts/
│   └── <job_name>.pbs                                            (Phase 6)
├── logs/
│   ├── master.log
│   ├── phase_status.jsonl
│   ├── job_status.jsonl
│   ├── jobs/
│   └── errors/
│       ├── error_summary.txt
│       └── failed_jobs.csv
├── scratch/
└── tmp/
```

## Acceptance scorecard

See `fresh/docs/m1_acceptance_scorecard.md` for the per-item table. **15/15 closed** in Codex env. Items #6 (qsub run) and #15 (smoke_input on real files) have additional HPC-side validation steps documented in the same scorecard.

## HPC-side validation (user)

Two scripts emit a concrete PBS file under the run dir and print the `qsub` command for manual submission:

```bash
bash fresh/scripts/submit_smoke_env.sh   [<run_id>] [<node>]
bash fresh/scripts/submit_smoke_input.sh [<run_id>] [<node>]
```

Neither script auto-calls qsub. After the user submits the generated PBS file with `qsub`, the standard `qstat` / `cli status` cycle confirms completion.

## Transition to Milestone 2

Once Phase 9 (Tasks 4-9 schema realignment) lands, the workflow can move to Milestone 2 actual execution per `egfr_myo1d_overall_implementation_plan_milestones_1_3_v1_0.md`:

```text
M2.1 PPI input generation (real PyRosetta-ready inputs from M1 normalized outputs)
M2.2 PyRosetta adapter (smoke -> mini -> production scale ladder)
M2.3 PPI pose QC + MYO1D artifact filtering
M2.4 Symmetry-aware consensus patch
M2.5 ATP-site reference
M2.6 fpocket pocket discovery
M2.7 Membrane / dimer / PPI pocket gates
M2.8 Milestone 2 aggregation
```

The CLI surface added in Phases 1-8 provides the runtime infrastructure for M2.1: receptor normalize, membrane frame, MYO1D construct, and the prepare-inputs orchestrator all emit M1 canonical outputs that M2.1 will consume.

## Documentation index

```text
fresh/docs/
├── milestone1_foundation_plan.md          (this document)
├── m1_acceptance_scorecard.md             (final §23 closure status)
├── m1_completion_rework_handoff.md        (Codex agent handoff)
├── m1_phase1_cleanup_manager.md
├── m1_phase1_changes.md
├── m1_phase2_myo1d_relocation.md
├── m1_phase2_changes.md
├── m1_phase3_myo1d_construct_qc.md
├── m1_phase3_changes.md
├── m1_phase4_receptor_normalization.md
├── m1_phase4_changes.md
├── m1_phase5_membrane_frame_generation.md
├── m1_phase5_changes.md
├── m1_phase6_pbs_generator.md
├── m1_phase6_changes.md
├── m1_phase7_ligand_manifest.md
├── m1_phase7_changes.md
├── m1_phase8_prepare_inputs_integration.md
├── m1_phase8_changes.md
└── prompts/                               (per-phase Codex-style prompt + checklist; 18 docs)
```
