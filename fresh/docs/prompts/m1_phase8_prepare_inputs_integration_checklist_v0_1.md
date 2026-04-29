# M1 Phase 8 Acceptance Checklist v0.1 — prepare-inputs Orchestrator + M1 Integration

Use this after the implementer applies M1 Phase 8.

## 1. Pre-Phase state preserved

```text
Old workflow files unchanged.
Phases 1-7 module logic unchanged (only cli.py has the new orchestrator subparser).
Phase 9 Tasks 4-9 realignment NOT YET DONE.
```

## 2. CLI registered

```bash
python -m egfr_myo1d.cli --help | grep prepare-inputs
python -m egfr_myo1d.cli prepare-inputs --help
```

Help text must include `--run-id`, `--mode`, `--profile`, `--input-root`, `--states`, `--skip-ligands`, `--strict`.

The full subcommand list should now include 17 entries: `version`, `init-run`, `preflight`, `status`, `validate-structures`, `prepare-ppi-inputs`, `validate-real-inputs`, `plan-ppi-sampling`, `summarize-ppi-consensus`, `plan-pocket-discovery`, `prioritize-pocket-candidates`, `cleanup`, `prepare-myo1d`, `prepare-receptor`, `compute-membrane-frame`, `prepare-pbs`, `manifest-ligands`, `prepare-inputs`.

(Note: that's 18 by my count above; the exact intended count is 17 from "11 existing + 6 new + prepare-inputs as the 7th new" — implementer should verify and document the actual count.)

## 3. Full M1 closure walk on integration fixture

```bash
python -m egfr_myo1d.cli init-run --mode smoke_input --run-id m1_closure_smoke
python -m egfr_myo1d.cli preflight --run-id m1_closure_smoke --mode smoke_input --profile codex_dev
python -m egfr_myo1d.cli prepare-inputs --run-id m1_closure_smoke --mode smoke_input --profile codex_dev --input-root fresh/tests/fixtures/m1_phase8_integration
python -m egfr_myo1d.cli prepare-pbs --run-id m1_closure_smoke --job-name m1_closure_test --mode smoke_env --node node04
python -m egfr_myo1d.cli cleanup --run-id m1_closure_smoke --mode test --dry-run false
python -m egfr_myo1d.cli status --run-id m1_closure_smoke
```

Expected: every command exits 0; aggregate manifest reports PASS or PASS_WITH_WARNINGS.

## 4. Aggregate manifest

```bash
cat fresh/runs/m1_closure_smoke/manifest/prepare_inputs_aggregate_manifest.json | python -m json.tool
```

Expected fields:

```text
- run_id, mode, profile, input_root
- states_processed: list of 3 states
- sub_steps: list of 5+ entries (preflight, prepare-receptor:STATE for each, compute-membrane-frame, prepare-myo1d, manifest-ligands)
- each sub_step has name, status (PASS/WARN/FAIL), manifest_path, warnings
- missing_required_inputs (may be empty for fixture run)
- warnings: list
- blockers: list (empty for healthy fixture)
- status: PASS or PASS_WITH_WARNINGS
- timestamp
```

## 5. Summary report

```bash
cat fresh/runs/m1_closure_smoke/reports/prepare_inputs_summary.md
```

Expected: human-readable Markdown summary with sub-step status table, warnings, links to manifests/QCs.

## 6. Required outputs from sub-steps (M1 §23 closure verification)

```bash
test -f fresh/runs/m1_closure_smoke/manifest/run_manifest.json                                           # Task 2
test -f fresh/runs/m1_closure_smoke/manifest/environment_report.json                                     # Task 2
test -f fresh/runs/m1_closure_smoke/normalized/myo1d/MYO1D_955_1006.pdb                                  # Phase 3
test -f fresh/runs/m1_closure_smoke/qc/myo1d_construct_qc.csv                                            # Phase 3
test -f fresh/runs/m1_closure_smoke/normalized/receptors/EGFR_160-185_dockable_669_1014_explicit_AB.pdb  # Phase 4
test -f fresh/runs/m1_closure_smoke/normalized/receptors/EGFR_160-185_runtime_offset_receptor_only.pdb   # Phase 4
test -f fresh/runs/m1_closure_smoke/qc/EGFR_160-185_receptor_mapping.csv                                 # Phase 4
test -f fresh/runs/m1_closure_smoke/manifest/membrane_frame.json                                         # Phase 5
test -f fresh/runs/m1_closure_smoke/qc/membrane_frame_qc.csv                                             # Phase 5
test -f fresh/runs/m1_closure_smoke/scripts/m1_closure_test.pbs                                          # Phase 6
test -f fresh/runs/m1_closure_smoke/qc/ligand_manifest_qc.csv                                            # Phase 7
test -f fresh/runs/m1_closure_smoke/manifest/cleanup_report.json                                         # Phase 1
```

All must exist.

## 7. Acceptance scorecard

```bash
cat fresh/docs/m1_acceptance_scorecard.md
```

Expected: table with 15 rows (#1-15 from handoff §23). Each row: ItemNumber, Description, Status (DONE/HPC_PENDING/PARTIAL/MISSING), Phase, Verification, Artifact path.

Programmatic check:

```bash
pytest -q fresh/tests/test_m1_phase8_prepare_inputs_integration.py::test_m1_integration_acceptance_scorecard_15_items
```

Must pass; all 15 items DONE or HPC_PENDING (none MISSING/PARTIAL).

## 8. Missing real files (smoke_input + codex_dev)

```bash
python -m egfr_myo1d.cli init-run --mode smoke_input --run-id m1_missing_files
python -m egfr_myo1d.cli prepare-inputs --run-id m1_missing_files --mode smoke_input --profile codex_dev
```

Expected (M1 §23 #15):

```text
- exit code: 0 (PASS_WITH_WARNINGS or 1 only if hpc_strict)
- aggregate manifest.missing_required_inputs is non-empty
- warnings clearly identify which files are missing
- NO crash, NO traceback
- normalized/ outputs may be missing per-state (recorded as missing in aggregate, not crashed)
- manifest/cleanup_report.json may be absent (cleanup not run in this command); independent cleanup invocation should still work
```

## 9. hpc_strict failure propagation

```bash
python -m egfr_myo1d.cli init-run --mode smoke_input --run-id m1_strict_missing
python -m egfr_myo1d.cli prepare-inputs --run-id m1_strict_missing --mode smoke_input --profile hpc_strict --input-root /tmp/no_such_input_root
```

Expected: exit 1; aggregate status FAIL; sub_steps show first failing step then stop or record subsequent skipped.

## 10. --skip-ligands

```bash
python -m egfr_myo1d.cli prepare-inputs --run-id m1_no_ligands --mode smoke_input --profile codex_dev --input-root fresh/tests/fixtures/m1_phase8_integration --skip-ligands true
```

Expected: aggregate sub_steps does NOT include manifest-ligands; ligand_manifest_qc.csv NOT written.

## 11. Status command after M1 closure walk

```bash
python -m egfr_myo1d.cli status --run-id m1_closure_smoke
```

Expected stdout includes:

```text
- run_id
- run_dir
- last phase status (one of the prepare-inputs sub-steps or aggregate)
- WARN/FAIL counts
- master.log path
- aggregate manifest path
```

## 12. Path traversal

```bash
python -m egfr_myo1d.cli prepare-inputs --run-id ../bad_run --mode smoke_input
```

Nonzero exit, no outside writes.

## 13. Tests

```bash
pytest -q fresh/tests/test_m1_phase8_prepare_inputs_integration.py
pytest -q fresh/tests
```

Expected: ≥10 phase 8 tests pass; total suite (98 prior + Phases 1-7 new + Phase 8 new) all pass.

## 14. Old workflow protection

```bash
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

Empty.

## 15. milestone1_foundation_plan.md updated

```bash
wc -l fresh/docs/milestone1_foundation_plan.md
head -10 fresh/docs/milestone1_foundation_plan.md
```

Expected: significantly more than the prior 25 lines; sections covering Purpose, Acceptance scorecard, Module tree, CLI reference, Run output schema, HPC-pending, Transition to M2.

## 16. M1→M2 transition gate (v1.0 §14.1)

Per item, verify:

```text
1. fresh/ skeleton complete                       — Phase 1+ artifacts confirm
2. qsub smoke_env complete                        — HPC_PENDING (PBS file ready Phase 6)
3. qsub smoke_input complete                      — HPC_PENDING (PBS file ready Phase 6)
4. receptor normalization for ≥1 primary state    — Phase 4 + integration fixture confirms
5. MYO1D 955-1006 construct QC                    — Phase 3 + integration fixture confirms
6. membrane_frame.json generated/inherited        — Phase 5 + integration fixture confirms
7. logs centralized                               — Task 2 + integration test confirms
8. cleanup_report generated                       — Phase 1 + integration walk confirms
9. pytest fresh/tests -q pass                     — full suite passes
```

Items 1, 4-9 DONE in Codex env. Items 2, 3 HPC_PENDING.

## 17. What must not be in this phase

```text
- Tasks 4-9 schema realignment (Phase 9)
- M2 docking work
- modifying Phases 1-7 module logic
- HPC qsub execution
- modifying old workflow files
```

## 18. Phase 8 accepted if

```text
- prepare-inputs CLI subcommand registered.
- Orchestrator runs preflight + prepare-receptor (per state) + compute-membrane-frame + prepare-myo1d + manifest-ligands in order under one RunContext.
- Aggregate manifest emitted with sub_steps detail.
- Summary report emitted in Markdown.
- Missing real files in smoke_input/codex_dev produce missing_required_inputs warnings, no crash.
- hpc_strict propagates FAIL correctly.
- --skip-ligands works.
- M1 integration smoke (synthetic fixture) PASSes; all expected M1 outputs present after walk.
- M1 acceptance scorecard 15/15 (HPC_PENDING annotated; no MISSING/PARTIAL).
- M1→M2 transition gate items 1, 4-9 closed in Codex env; 2-3 HPC_PENDING.
- milestone1_foundation_plan.md expanded from stub.
- m1_acceptance_scorecard.md created.
- ≥8 phase tests pass; full suite passes.
- M1 §23 #15 closed.
- Old workflow files unmodified.
```

## 19. Implementer final response must include

```text
M1 Phase 8 status: PASS / PASS WITH WARNINGS / FAIL
M1 closure scorecard: 15/15 (HPC_PENDING annotated)
Files created
Files modified
Commands run and results
Test summary (prior + Phases 1-7 + Phase 8 new = total)
M1 integration smoke walk evidence (commands + exit codes + key file checks)
Aggregate manifest sample
Acceptance scorecard contents
M1→M2 transition gate status per item
HPC_PENDING items documented
Old workflow protection
Known limitations: HPC qsub user-side; real input files user-side
Next: Phase 9 Tasks 4-9 schema realignment
```
