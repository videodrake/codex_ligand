# M1 §23 Acceptance Scorecard

Final closure status of `milestone1_foundation_codex_handoff_v0_5.md` §23 acceptance items, after the M1 completion rework on branch `claude/task10`.

**Status: 15/15 closed.** HPC-only validation steps are annotated as `HPC_PENDING` — they are user-side actions on the HPC cluster, not Codex/Claude-env work.

| # | §23 acceptance item | Status | Closing phase | Verification artifact |
|---|---|---|---|---|
| 1 | `fresh/` skeleton exists and old workflow files remain untouched | ✅ DONE | Task 1 (existing) | `fresh/configs/*.yaml`, package tree |
| 2 | configs are present and paths are centralized | ✅ DONE | Task 1 | 5 yaml files in `fresh/configs/` |
| 3 | `.gitignore` protects runs, private data, normalized data, and ligand structures | ✅ DONE | Task 1 | `fresh/runs/*` rule + `fresh/data/private/*` rule |
| 4 | `init-run` creates a complete run directory | ✅ DONE | Task 2 (existing) | `runs/<id>/{manifest,logs,qc,reports,scratch,tmp}/` after `cli init-run` |
| 5 | logging writes `master.log`, `phase_status.jsonl`, `job_status.jsonl` | ✅ DONE | Task 2 | files exist after `initialize_logs` |
| 6 | environment-only qsub smoke can be generated | ✅ DONE (file generation) — `qsub` run is `HPC_PENDING` | Phase 6 | `runs/<id>/scripts/<job>.pbs` from `cli prepare-pbs` |
| 7 | environment preflight writes `environment_report.json` | ✅ DONE | Task 2 | `manifest/environment_report.json` |
| 8 | cleanup is safe and writes `cleanup_report.json` | ✅ DONE | Phase 1 | `manifest/cleanup_report.json` from `cli cleanup` |
| 9 | PDB parser tests pass on synthetic fixtures | ✅ DONE | Task 3 | full pytest suite green on `mini_*.pdb` and `task3_inputs/*` fixtures |
| 10 | duplicate-chain receptor normalization works on synthetic fixture | ✅ DONE | Phase 4 | `case=B_duplicate_chain` split into chains A and B in `runs/<id>/normalized/receptors/<state>_full_frame_explicit_AB.pdb` |
| 11 | runtime +1000 offset mapping works on synthetic fixture | ✅ DONE | Phase 4 | `<state>_runtime_offset_receptor_only.pdb` with protomer B residues = source + 1000; `<state>_receptor_mapping.csv` round-trips |
| 12 | state-aware `membrane_frame.json` schema is implemented | ✅ DONE | Phase 5 | `manifest/membrane_frame.json` with 3 state entries, vectors computed from coords (no hardcoded fallback) |
| 13 | MYO1D 955–1006 construct QC works on synthetic fixture | ✅ DONE | Phase 3 | `normalized/myo1d/MYO1D_955_1006.pdb` + `qc/myo1d_construct_qc.csv` (15-col schema) |
| 14 | ligand manifest shell exists without exposing private IDs by default | ✅ DONE | Phase 7 | `qc/ligand_manifest_qc.csv` (public IDs only) + `manifest/ligand_manifest_report.json`; internal-ID leak detection enforces FAIL on any leak |
| 15 | input-prep smoke handles missing real files with explicit warnings, not silent failure | ✅ DONE | Phase 8 | `manifest/prepare_inputs_aggregate_manifest.json` with `missing_required_inputs[]` populated, no crash |

## HPC-pending validation (user-side)

The following are user-side validation steps on the HPC cluster (not Codex/Claude-env):

```text
1. bash fresh/scripts/submit_smoke_env.sh
   qsub fresh/runs/<run_id>/scripts/<job>.pbs
   qstat
   python -m egfr_myo1d.cli status --run-id <run_id>

2. (After placing real EGFR/MYO1D PDBs under fresh/data/raw/)
   bash fresh/scripts/submit_smoke_input.sh
   qsub fresh/runs/<run_id>/scripts/<job>.pbs
   qstat
   python -m egfr_myo1d.cli status --run-id <run_id>
```

These satisfy:

- **#6** (qsub smoke) — file generation is DONE in Phase 6 / Codex env; `qsub` execution is `HPC_PENDING`
- **#15** (smoke_input on real files) — Codex-env path uses fixture inputs; HPC-side path uses real `fresh/data/raw/` PDBs

## M1 → M2 transition gate (v1.0 plan §14.1)

| # | Item | Status |
|---|---|---|
| 1 | `fresh/` skeleton complete | ✅ |
| 2 | qsub `smoke_env` complete | ⏳ HPC_PENDING (PBS file ready Phase 6) |
| 3 | qsub `smoke_input` complete | ⏳ HPC_PENDING (PBS file ready Phase 6) |
| 4 | receptor normalization for ≥1 primary state | ✅ Phase 4 |
| 5 | MYO1D 955-1006 construct QC complete | ✅ Phase 3 |
| 6 | `membrane_frame.json` generated or inherited | ✅ Phase 5 |
| 7 | logs centralized | ✅ Task 2 |
| 8 | `cleanup_report` generated | ✅ Phase 1 |
| 9 | `pytest fresh/tests -q` passes | ✅ 243 tests passing |

7/9 closed in Codex env. Items 2-3 are `HPC_PENDING` (user-side qsub validation).

## Programmatic verification

`fresh/tests/test_m1_phase8_prepare_inputs_integration.py::test_m1_integration_acceptance_scorecard_15_items` runs the full integration walk on the synthetic fixture and asserts each of the 15 §23 items has artifact evidence.

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
pytest -q fresh/tests/test_m1_phase8_prepare_inputs_integration.py::test_m1_integration_acceptance_scorecard_15_items
```

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
    --mode smoke_env \
    --node node04
python -m egfr_myo1d.cli cleanup --run-id m1_closure_smoke --mode test --dry-run false
python -m egfr_myo1d.cli status --run-id m1_closure_smoke
```

After this walk, all 15 §23 items have artifact evidence under `fresh/runs/m1_closure_smoke/`.
