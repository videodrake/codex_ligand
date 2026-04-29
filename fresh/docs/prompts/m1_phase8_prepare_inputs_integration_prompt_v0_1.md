# Claude M1 Phase 8 Prompt — prepare-inputs Orchestrator + M1 Integration Smoke Test v0.1

Branch `claude/task10`. Phases 1-7 complete. This is **M1 Phase 8** — implements the `prepare-inputs` orchestrator CLI and M1 integration smoke test per `milestone1_foundation_codex_handoff_v0_5.md` §10.2 Smoke B, §19, §22, §23, closing M1 §23 #15 and validating items #1-15 end-to-end.

## 1. Project context

Phases 1-7 created the M1 foundation modules (cleanup, MYO1D construct + QC, receptor normalize, membrane frame, PBS gen, ligand manifest). Each has its own CLI subcommand. Phase 8 ties them together with a single `prepare-inputs` orchestrator that runs them in sequence, plus an integration smoke test that exercises the full M1 acceptance scorecard.

## 2. Absolute rules

Do not modify the old workflow. Do not modify Phases 1-7 module logic — Phase 8 is purely orchestration + integration testing. Maintain Py2/3 syntax compatibility.

The orchestrator must respect each sub-step's profile/severity policy; it should NOT downgrade FAILs to WARNs.

## 3. Scope

In scope:
- Add `prepare-inputs` CLI subcommand (orchestrator) to `cli.py`
- Update `fresh/docs/milestone1_foundation_plan.md` from Task 1 stub (~25 lines) to full M1 closure documentation (~150-300 lines)
- Add `fresh/docs/m1_acceptance_scorecard.md` — final §23 closure scorecard
- Tests under `fresh/tests/test_m1_phase8_prepare_inputs_integration.py` (≥8 tests)
- Fixtures under `fresh/tests/fixtures/m1_phase8_integration/` — combined fixture set for end-to-end run
- Docs `fresh/docs/m1_phase8_prepare_inputs_integration.md` and `m1_phase8_changes.md`

Out of scope:
- Modifying any Phase 1-7 module
- Tasks 4-9 realignment (Phase 9)
- M2 work
- Real qsub on HPC

## 4. Required CLI behavior

```bash
python -m egfr_myo1d.cli prepare-inputs \
  --run-id RUN \
  --mode smoke_input \
  [--profile codex_dev|hpc_strict] \
  [--input-root fresh/data/raw] \
  [--states EGFR_160-185,EGFR_170-200,3GT8_raw]   # comma-sep; default: all per receptor_states.yaml
  [--skip-ligands true|false]                       # default false
  [--strict]
```

Behavior — runs sub-steps in this order:

```text
1. preflight (mode=smoke_input, profile inherited)
2. for each state: prepare-receptor (state, source from receptor_states.yaml)
3. compute-membrane-frame (state=all, full-frame source from receptor_states.yaml)
4. prepare-myo1d (source from fresh_run.yaml)
5. manifest-ligands (unless --skip-ligands true)
```

Each sub-step:
- Runs in the same RunContext (single run_id)
- Appends its phase status entry
- If a sub-step FAILs and `--profile=hpc_strict`: orchestrator stops and exits 1
- If a sub-step FAILs and `--profile=codex_dev`: orchestrator records FAIL and continues to next sub-step (so user sees full picture). Final exit is 1 if any sub-step FAILed.
- Missing real input files (e.g., `fresh/data/raw/receptors/EGFR_160-185.pdb`): produces `missing_required_inputs` in manifest, NOT crash.

Aggregated output:
- `fresh/runs/<run_id>/manifest/prepare_inputs_aggregate_manifest.json`
- `fresh/runs/<run_id>/reports/prepare_inputs_summary.md`

## 5. Files to create / modify

Create:

```text
fresh/tests/test_m1_phase8_prepare_inputs_integration.py
fresh/tests/fixtures/m1_phase8_integration/receptors/EGFR_160-185_synthetic.pdb
fresh/tests/fixtures/m1_phase8_integration/receptors/EGFR_170-200_synthetic.pdb
fresh/tests/fixtures/m1_phase8_integration/receptors/3GT8_raw_synthetic.pdb
fresh/tests/fixtures/m1_phase8_integration/receptors/plus10_full_frame_synthetic.pdb
fresh/tests/fixtures/m1_phase8_integration/myo1d/AF-O94832-F1-model_v6_synthetic.pdb
fresh/tests/fixtures/m1_phase8_integration/ligands/Cpd-A.sdf
fresh/tests/fixtures/m1_phase8_integration/ligands/Cpd-B.sdf
fresh/tests/fixtures/m1_phase8_integration/ligands/Cpd-C.sdf
fresh/tests/fixtures/m1_phase8_integration/private/compound_id_map.csv     # synthetic placeholder internals
fresh/tests/fixtures/m1_phase8_integration/test_input_root.yaml             # local override pointing all paths into the fixture tree
fresh/docs/m1_phase8_prepare_inputs_integration.md
fresh/docs/m1_phase8_changes.md
fresh/docs/m1_acceptance_scorecard.md
```

(Reuse fixtures from Phases 3-7 where possible — copy, do not duplicate logic.)

Modify:

```text
fresh/src/egfr_myo1d/cli.py                       # add prepare-inputs orchestrator subparser + handler
fresh/docs/milestone1_foundation_plan.md          # replace Task 1 stub with full M1 closure plan
```

## 6. Public API

The orchestrator can be a `_cmd_prepare_inputs` handler in cli.py that calls the sub-step Python functions directly (not via subprocess). This avoids spawning new processes and keeps a single RunContext.

```python
def _cmd_prepare_inputs(args):
    ctx = RunContext.create(args.run_id, args.mode) if not exists else RunContext.for_existing(args.run_id)
    initialize_logs(ctx)

    aggregate = {"sub_steps": [], "status": "PASS", "warnings": [], "blockers": []}

    # Step 1: preflight
    preflight_result = run_preflight(ctx, args.mode, args.profile)
    aggregate["sub_steps"].append({"name": "preflight", "status": preflight_result.status})
    if preflight_result.status == "FAIL" and args.profile == "hpc_strict":
        return _finalize_aggregate(ctx, aggregate, exit_code=1)

    # Step 2: prepare-receptor per state
    for state in resolve_states(args.states):
        source = resolve_state_source(state, args.input_root)
        result = normalize_receptor(ctx, source, state, args.profile, args.strict)
        aggregate["sub_steps"].append({"name": "prepare-receptor:" + state, "status": result.status})
        if result.status == "FAIL" and args.profile == "hpc_strict":
            return _finalize_aggregate(ctx, aggregate, exit_code=1)

    # Step 3: compute-membrane-frame
    frame_result = compute_all_membrane_frames(ctx, args.input_root, args.profile)
    aggregate["sub_steps"].append({"name": "compute-membrane-frame", "status": frame_result.status})

    # Step 4: prepare-myo1d
    myo1d_result = run_myo1d_qc(ctx, ...)
    aggregate["sub_steps"].append({"name": "prepare-myo1d", "status": myo1d_result.status})

    # Step 5: manifest-ligands
    if not args.skip_ligands:
        ligand_result = build_ligand_manifest(ctx, ...)
        aggregate["sub_steps"].append({"name": "manifest-ligands", "status": ligand_result.status})

    return _finalize_aggregate(ctx, aggregate)
```

Aggregate manifest schema:

```json
{
  "run_id": "...",
  "mode": "smoke_input",
  "profile": "codex_dev",
  "input_root": "...",
  "states_processed": [...],
  "sub_steps": [
    {"name": "preflight", "status": "PASS|WARN|FAIL", "manifest_path": "...", "warnings": [...]},
    {"name": "prepare-receptor:EGFR_160-185", "status": "...", ...},
    ...
  ],
  "missing_required_inputs": [...],
  "warnings": [...],
  "blockers": [...],
  "status": "PASS|PASS_WITH_WARNINGS|FAIL",
  "timestamp": "..."
}
```

## 7. Required output files

After `prepare-inputs` runs (in addition to outputs from each sub-step):

```text
fresh/runs/<run_id>/manifest/prepare_inputs_aggregate_manifest.json
fresh/runs/<run_id>/reports/prepare_inputs_summary.md
fresh/runs/<run_id>/logs/phase_status.jsonl                    # appended (one entry per sub-step + one aggregate)
fresh/runs/<run_id>/logs/master.log                            # appended
```

## 8. M1 integration smoke test fixture

The fixture at `fresh/tests/fixtures/m1_phase8_integration/` must contain enough material for one full pass of `prepare-inputs` in `--profile codex_dev` to PASS_WITH_WARNINGS or PASS:

```text
- 3 synthetic receptor PDBs (EGFR_160-185, EGFR_170-200, 3GT8_raw) — explicit AB chains, residues spanning 634-1014 (or close enough for crop), at least one with V924R-like marker, at least one with TM/JM 634-674 residues for membrane frame
- 1 synthetic plus10_full_frame PDB with TM/JM
- 1 synthetic MYO1D PDB with residues 955-1006 + key residues
- 3 minimal valid SDF placeholders for Cpd-A/B/C
- 1 synthetic compound_id_map.csv with placeholder internals
- 1 test_input_root.yaml pointing all paths into the fixture tree
```

The fixtures may be small but must be syntactically valid PDBs/SDFs.

## 9. Tests required (≥8)

```text
test_prepare_inputs_runs_all_substeps_in_order
test_prepare_inputs_handles_missing_real_files_with_explicit_warnings_in_smoke_input
test_prepare_inputs_aggregate_phase_status_records_each_substep
test_prepare_inputs_fails_cleanly_when_a_substep_fails_in_hpc_strict
test_prepare_inputs_continues_when_a_substep_fails_in_codex_dev_but_aggregate_status_is_FAIL
test_prepare_inputs_outputs_under_run_dir_only
test_prepare_inputs_skip_ligands_flag_skips_manifest_step
test_m1_integration_synthetic_full_pipeline_pass
test_m1_integration_acceptance_scorecard_15_items
test_cli_help_includes_prepare_inputs
```

(10 tests; ≥8 required.)

`test_m1_integration_synthetic_full_pipeline_pass` runs (programmatically, not via subprocess):

```python
def test_m1_integration_synthetic_full_pipeline_pass(tmp_path):
    run_id = unique_run_id("m1_int")
    fixture_root = REPO_ROOT / "fresh" / "tests" / "fixtures" / "m1_phase8_integration"
    # init-run
    ctx = init_run(run_id, mode="smoke_input")
    # preflight
    preflight_result = run_preflight(ctx, "smoke_input", "codex_dev")
    assert preflight_result.status in ("PASS", "WARN")
    # prepare-inputs
    aggregate = run_prepare_inputs(ctx, mode="smoke_input", profile="codex_dev", input_root=fixture_root)
    assert aggregate["status"] in ("PASS", "PASS_WITH_WARNINGS")
    # cleanup test mode
    cleanup_result = run_cleanup(ctx, mode="test", dry_run=False, profile="codex_dev")
    assert cleanup_result.status in ("PASS", "WARN")
    # Verify all expected M1 outputs present
    assert (ctx.run_dir / "normalized" / "myo1d" / "MYO1D_955_1006.pdb").is_file()
    assert (ctx.run_dir / "normalized" / "receptors" / "EGFR_160-185_dockable_669_1014_explicit_AB.pdb").is_file()
    assert (ctx.run_dir / "normalized" / "receptors" / "EGFR_160-185_runtime_offset_receptor_only.pdb").is_file()
    assert (ctx.run_dir / "manifest" / "membrane_frame.json").is_file()
    assert (ctx.run_dir / "qc" / "ligand_manifest_qc.csv").is_file()
    assert (ctx.run_dir / "manifest" / "cleanup_report.json").is_file()
```

`test_m1_integration_acceptance_scorecard_15_items`: programmatically verify each of M1 §23 #1-15 by checking file existence, manifest fields, or running `cli status` and inspecting the result. Mark HPC-only items (#6 qsub run) as `HPC_PENDING` rather than failed.

## 10. Acceptance commands (M1 closure walk)

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"

# Targeted phase test
pytest -q fresh/tests/test_m1_phase8_prepare_inputs_integration.py

# Full suite — must pass with prior 98 + Phases 1-7 new + Phase 8 new tests
pytest -q fresh/tests

# Full M1 acceptance walk on the integration fixture
python -m egfr_myo1d.cli init-run --mode smoke_input --run-id m1_closure_smoke
python -m egfr_myo1d.cli preflight --run-id m1_closure_smoke --mode smoke_input --profile codex_dev
python -m egfr_myo1d.cli prepare-inputs --run-id m1_closure_smoke --mode smoke_input --profile codex_dev --input-root fresh/tests/fixtures/m1_phase8_integration
python -m egfr_myo1d.cli prepare-pbs --run-id m1_closure_smoke --job-name m1_closure_test --mode smoke_env --node node04
python -m egfr_myo1d.cli cleanup --run-id m1_closure_smoke --mode test --dry-run false
python -m egfr_myo1d.cli status --run-id m1_closure_smoke

# Acceptance scorecard inspection
cat fresh/runs/m1_closure_smoke/manifest/prepare_inputs_aggregate_manifest.json
cat fresh/runs/m1_closure_smoke/reports/prepare_inputs_summary.md

# Path traversal
python -m egfr_myo1d.cli prepare-inputs --run-id ../bad_run --mode smoke_input --profile codex_dev

# Missing real files behavior (no fixture override; should produce missing_required_inputs)
python -m egfr_myo1d.cli init-run --mode smoke_input --run-id m1_missing_files
python -m egfr_myo1d.cli prepare-inputs --run-id m1_missing_files --mode smoke_input --profile codex_dev

# Old workflow protection
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

## 11. Documentation requirements

`fresh/docs/milestone1_foundation_plan.md` (replace stub):

Sections:
- Purpose of M1
- M1 acceptance scorecard mapping (link to m1_acceptance_scorecard.md)
- Module tree (post-M1) reference
- CLI command reference (all 17 subparsers including new ones)
- Run output schema reference
- HPC-pending items (qsub validation)
- Transition to M2

`fresh/docs/m1_acceptance_scorecard.md`:

Per M1 §23 #1-15, table with columns: ItemNumber, Description, Status (DONE/HPC_PENDING/PARTIAL/MISSING), Phase that closed it, Verification command, Output artifact path.

For HPC_PENDING items, document the user-side step needed.

## 12. Final response format

```text
M1 Phase 8 status: PASS / PASS WITH WARNINGS / FAIL
M1 closure status: 15/15 items closed (HPC_PENDING annotated)

Files created:
- ...

Files modified:
- cli.py (added prepare-inputs orchestrator)
- fresh/docs/milestone1_foundation_plan.md (expanded from stub)
- fresh/docs/m1_acceptance_scorecard.md (new)

Commands run and results

Test summary:
- prior: 98
- Phases 1-7 new: ~75
- Phase 8 new: ~10
- total: ~183 passing

M1 integration scorecard (#1-15):
1. fresh/ skeleton — DONE (Phase 1 onwards retained)
2. configs — DONE
3. .gitignore — DONE
4. init-run — DONE (Task 2)
5. logging — DONE (Task 2)
6. qsub smoke generation — DONE (Phase 6); qsub run HPC_PENDING
7. preflight — DONE (Task 2)
8. cleanup — DONE (Phase 1)
9. PDB parser — DONE (Task 3)
10. duplicate-chain normalization — DONE (Phase 4)
11. +1000 runtime offset — DONE (Phase 4)
12. state-aware membrane_frame.json — DONE (Phase 5)
13. MYO1D 955-1006 construct QC — DONE (Phase 3)
14. ligand manifest shell — DONE (Phase 7)
15. input-prep smoke handles missing real files — DONE (Phase 8)

M1 → M2 transition gate (v1.0 §14.1):
1. fresh/ skeleton complete — yes
2. qsub smoke_env complete — HPC_PENDING (file ready, user runs qsub)
3. qsub smoke_input complete — HPC_PENDING
4. receptor normalization for ≥1 primary state — DONE (synthetic fixture)
5. MYO1D 955-1006 construct QC — DONE
6. membrane_frame.json generated/inherited — DONE
7. logs centralized — DONE
8. cleanup_report generated — DONE
9. pytest fresh/tests -q pass — DONE

Old workflow protection: empty diff
Acceptance closure: M1 §23 #15 closed; full M1 acceptance scorecard 15/15 (HPC_PENDING annotated for #6 qsub run, #15 real-file run)

Known limitations:
- HPC user runs `bash fresh/scripts/submit_smoke_env.sh` and `bash fresh/scripts/submit_smoke_input.sh` for actual qsub validation
- Real EGFR/MYO1D PDB files placed by user under fresh/data/raw/ (not in repo); fixture tree exercises the codex_dev path
- Phase 9 (Tasks 4-9 schema realignment) is the next phase, post-M1-closure
```
