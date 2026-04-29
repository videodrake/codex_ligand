# Claude M1 Phase 9 Prompt — Tasks 4-9 Schema Realignment v0.1

Branch `claude/task10`. Phases 1-8 complete; M1 closed. This is **M1 Phase 9** — realigns Tasks 4-9 to consume M1 canonical normalized outputs instead of pass-through audited inputs. Closes nothing new in §23 (all already closed by Phase 8); brings Tasks 4-9 into M1-canonical alignment per the master plan.

## 1. Project context

After Phase 8, M1 produces canonical outputs at:

```text
fresh/runs/<run_id>/normalized/receptors/<state>_full_frame_explicit_AB.pdb
fresh/runs/<run_id>/normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb
fresh/runs/<run_id>/normalized/receptors/<state>_runtime_offset_receptor_only.pdb
fresh/runs/<run_id>/normalized/myo1d/MYO1D_955_1006.pdb
fresh/runs/<run_id>/manifest/membrane_frame.json
fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv
fresh/runs/<run_id>/qc/myo1d_construct_qc.csv
fresh/runs/<run_id>/qc/membrane_frame_qc.csv
fresh/runs/<run_id>/qc/ligand_manifest_qc.csv
```

Tasks 4-9, written before M1 was complete, currently produce parallel outputs at:

```text
fresh/runs/<run_id>/prepared/egfr/egfr_receptor_normalized.pdb           # pass-through copy
fresh/runs/<run_id>/prepared/myo1d/MYO1D_sheet8_9_12_core_955_1001.pdb   # task 4 slicing
fresh/runs/<run_id>/prepared/myo1d/MYO1D_ext_beta_meander_955_1006_tail_masked.pdb
fresh/runs/<run_id>/prepared/restraints/<various>.json
fresh/runs/<run_id>/prepared/ppi/<various>.json
... etc
```

Phase 9 updates Tasks 4-9 to consume M1 normalized outputs as INPUTS, instead of redoing the normalization pass-through. This eliminates redundancy and ensures Tasks 4-9 evidence chains back to the canonical M1 normalization.

## 2. Absolute rules

Do not modify the old workflow. Do not modify Phases 1-8 modules in this phase (they are M1 foundation). Maintain Py2/3 syntax compatibility.

The semantics of Tasks 4-9 (audit, plan, summarize, prioritize) MUST be preserved. Only the INPUT path / SCHEMA they consume changes.

EGFR-side restraint contracts (`prepared/restraints/*.json`), masks, and PPI plan/consensus/pocket outputs are NOT covered by M1 — those stay as Task 4-9 outputs. Only the receptor and MYO1D PDB pass-through emissions are realigned.

## 3. Scope

In scope (modify Tasks 4-9 to consume M1 outputs):

- `validation/prepared_inputs.py` (Task 4): consume M1 normalized receptor + MYO1D PDBs; stop emitting pass-through `prepared/egfr/egfr_receptor_normalized.pdb` and `prepared/myo1d/*` PDBs (replace with reference path in manifest pointing back to M1 outputs). Keep restraint and mask emissions.
- `validation/real_inputs.py` (Task 5): update real-input readiness check to verify M1 normalized outputs exist (or kick off M1 normalization itself if not yet done). The "readiness bridge" semantically becomes "M1 output verification + production input registration".
- `validation/ppi_sampling_plan.py` (Task 6): update PPI job spec receptor_path/partner_path fields to reference M1 normalized outputs.
- `validation/ppi_consensus.py` (Task 7): update ppi_consensus_patch processing to reference M1 receptor mapping CSV for residue identity validation.
- `validation/pocket_discovery.py` (Task 8): update pocket plan receptor_path field to reference M1 dockable receptor.
- `validation/pocket_candidate_prioritization.py` (Task 9): update protomer_id resolution to use M1 receptor mapping CSV.

Update test fixtures and assertions in:

```text
fresh/tests/test_task4_ppi_input_preparation.py
fresh/tests/test_task5_real_input_readiness.py
fresh/tests/test_task6_ppi_sampling_plan.py
fresh/tests/test_task7_ppi_consensus_patch.py
fresh/tests/test_task8_pocket_discovery_plan.py
fresh/tests/test_task9_pocket_candidate_prioritization.py
```

Add ≥6 new realignment tests:

```text
fresh/tests/test_m1_phase9_tasks_realignment.py
```

Out of scope:
- Modifying Phases 1-8 (M1 foundation) modules
- M2 actual docking execution
- M3 work
- Modifying old workflow files

## 4. Required CLI behavior

No new CLI subcommands. Existing `prepare-ppi-inputs`, `validate-real-inputs`, `plan-ppi-sampling`, `summarize-ppi-consensus`, `plan-pocket-discovery`, `prioritize-pocket-candidates` continue to work but now require M1 normalized outputs as a prerequisite.

For ergonomics, when a Task 4-9 command runs and detects missing M1 outputs, the recommended action is to either:
- (Strict) fail with a clear message: `"M1 normalized outputs not found at runs/<run_id>/normalized/. Run `prepare-inputs` first."`
- (Permissive in codex_dev) auto-invoke `prepare-inputs` if `--auto-prepare` flag is passed (optional new flag).

Implementer chooses one approach (strict recommended for hpc_strict; permissive optional for codex_dev with explicit flag).

## 5. Files to create / modify

Create:

```text
fresh/tests/test_m1_phase9_tasks_realignment.py
fresh/docs/m1_phase9_tasks4to9_realignment.md
fresh/docs/m1_phase9_changes.md
```

Modify (logic + path updates only; preserve semantics):

```text
fresh/src/egfr_myo1d/validation/prepared_inputs.py        # Task 4: stop pass-through, consume M1
fresh/src/egfr_myo1d/validation/real_inputs.py            # Task 5: M1 output verification
fresh/src/egfr_myo1d/validation/ppi_sampling_plan.py      # Task 6: receptor/partner paths reference M1
fresh/src/egfr_myo1d/validation/ppi_consensus.py          # Task 7: mapping CSV reference
fresh/src/egfr_myo1d/validation/pocket_discovery.py       # Task 8: receptor reference
fresh/src/egfr_myo1d/validation/pocket_candidate_prioritization.py   # Task 9: protomer resolution
fresh/tests/test_task4_ppi_input_preparation.py           # update fixtures + assertions
fresh/tests/test_task5_real_input_readiness.py            # ditto
fresh/tests/test_task6_ppi_sampling_plan.py               # ditto
fresh/tests/test_task7_ppi_consensus_patch.py             # ditto
fresh/tests/test_task8_pocket_discovery_plan.py           # ditto
fresh/tests/test_task9_pocket_candidate_prioritization.py # ditto
fresh/docs/task[4-9]_*.md                                 # add note about M1 dependency
```

Optionally:
- Add `--auto-prepare` flag to Task 4-9 CLI commands for ergonomics

## 6. Per-Task realignment specifics

### Task 4 (validation/prepared_inputs.py)

Currently emits:
- `prepared/egfr/egfr_receptor_normalized.pdb` (pass-through copy of input)
- `prepared/myo1d/MYO1D_sheet8_9_12_core_955_1001.pdb` (sliced from raw input)
- `prepared/myo1d/MYO1D_ext_beta_meander_955_1006_tail_masked.pdb`

Realigned: consume M1 outputs as input:
- Read receptor from `runs/<run_id>/normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb`
- Read MYO1D from `runs/<run_id>/normalized/myo1d/MYO1D_955_1006.pdb`
- Continue applying MYO1D 955-1001 slicing (active-face core) — this is a Task 4 derivation, not duplicate normalization
- Continue applying MYO1D 955-1006 tail-masking (comparator) — also Task 4 derivation
- Update `prepared_input_manifest.json` to record M1 outputs as inputs (with sha256)
- Stop emitting pass-through `prepared/egfr/egfr_receptor_normalized.pdb` (replace with manifest reference to M1 output)
- Continue emitting `prepared/myo1d/*.pdb` derivations (since they ARE Task 4 derivations, not pass-through)

### Task 5 (validation/real_inputs.py)

Currently checks raw inputs at `fresh/data/raw/`. Realigned to also check M1 normalized outputs at `fresh/runs/<run_id>/normalized/` exist (or kick off M1 normalization).

### Task 6 (validation/ppi_sampling_plan.py)

Job spec fields:
- `receptor_path` → references `normalized/receptors/<state>_runtime_offset_receptor_only.pdb`
- `partner_path` → references `prepared/myo1d/MYO1D_sheet8_9_12_core_955_1001.pdb` (Task 4 derivation) or `normalized/myo1d/MYO1D_955_1006.pdb` (M1 raw construct)

### Task 7 (validation/ppi_consensus.py)

When parsing `accepted_ppi_contacts.csv` residue references, validate against M1 receptor mapping CSV (`qc/<state>_receptor_mapping.csv`) for protomer_id and runtime_resseq → source_resseq translation.

### Task 8 (validation/pocket_discovery.py)

Pocket plan output references `normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb` as the canonical pocket-search target.

### Task 9 (validation/pocket_candidate_prioritization.py)

When prioritizing candidates, use M1 mapping CSV to resolve runtime residue numbers (protomer B residues = source + 1000) back to original numbers for chemistry/biology interpretation.

## 7. Realignment tests required (≥6)

```text
test_task4_consumes_m1_normalized_receptor_path
test_task4_consumes_m1_normalized_myo1d_path
test_task4_no_longer_emits_pass_through_egfr_receptor_normalized
test_task5_verifies_m1_normalized_outputs_exist
test_task6_job_specs_reference_m1_normalized_receptor
test_task6_job_specs_reference_m1_normalized_or_task4_derived_partner
test_task7_residue_identity_validates_against_m1_mapping_csv
test_task8_pocket_plan_references_m1_dockable_receptor
test_task9_protomer_id_resolved_via_m1_mapping_csv
test_task4_through_9_full_pipeline_against_m1_outputs
test_legacy_prepared_paths_no_longer_emitted
```

(11 tests; ≥6 required.)

The last test `test_legacy_prepared_paths_no_longer_emitted` programmatically asserts that after a full M1 + Tasks 4-9 run, the files `prepared/egfr/egfr_receptor_normalized.pdb` are not present (or are now symlinks/manifest references).

## 8. Behavior policy

```text
- Tasks 4-9 must NOT re-do M1 normalization (no duplicate +1000 offset, no duplicate dockable crop, no duplicate MYO1D 955-1006 raw emission).
- Tasks 4-9 may continue producing DERIVATIONS (e.g., MYO1D 955-1001 active-face core sliced from M1 955-1006).
- Tasks 4-9 manifests must record M1 outputs as INPUTS (with sha256 from M1 manifest).
- If M1 outputs are missing, Task 4-9 commands should fail with a clear message (or auto-invoke prepare-inputs if --auto-prepare).
- Existing 79 Task 4-9 tests must continue passing after fixture updates (no test deletions; only fixture path adjustments).
- New realignment tests verify the consumption of M1 outputs.
```

## 9. Severity rules

```text
PASS:  Task 4-9 commands consume M1 outputs successfully; no duplicate normalization
WARN:  legacy pass-through path still emitted but as a backwards-compat shim (discouraged; should not occur in this phase)
FAIL:  M1 outputs missing AND --auto-prepare not provided AND no fallback; or duplicate normalization detected
```

## 10. Acceptance commands

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"

pytest -q fresh/tests/test_m1_phase9_tasks_realignment.py
pytest -q fresh/tests/test_task4_ppi_input_preparation.py
pytest -q fresh/tests/test_task5_real_input_readiness.py
pytest -q fresh/tests/test_task6_ppi_sampling_plan.py
pytest -q fresh/tests/test_task7_ppi_consensus_patch.py
pytest -q fresh/tests/test_task8_pocket_discovery_plan.py
pytest -q fresh/tests/test_task9_pocket_candidate_prioritization.py

# Full suite
pytest -q fresh/tests

# End-to-end M1 + Tasks 4-9 walk
python -m egfr_myo1d.cli init-run --mode smoke_input --run-id m1_phase9_local
python -m egfr_myo1d.cli prepare-inputs --run-id m1_phase9_local --mode smoke_input --profile codex_dev --input-root fresh/tests/fixtures/m1_phase8_integration

python -m egfr_myo1d.cli prepare-ppi-inputs --run-id m1_phase9_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task4_inputs
python -m egfr_myo1d.cli plan-ppi-sampling --run-id m1_phase9_local --mode smoke_input --profile codex_dev --input-root fresh/tests/fixtures/task4_inputs --contract fresh/tests/fixtures/task4_inputs/ppi_input_contract.json
python -m egfr_myo1d.cli summarize-ppi-consensus --run-id m1_phase9_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task7_ppi_consensus
python -m egfr_myo1d.cli plan-pocket-discovery --run-id m1_phase9_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task8_pocket_planning
python -m egfr_myo1d.cli prioritize-pocket-candidates --run-id m1_phase9_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task9_pocket_candidates

# Verify no duplicate normalization
test ! -f fresh/runs/m1_phase9_local/prepared/egfr/egfr_receptor_normalized.pdb && echo "OK: legacy pass-through gone"

# Status
python -m egfr_myo1d.cli status --run-id m1_phase9_local

# Old workflow protection
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

## 11. Final response format

```text
M1 Phase 9 status: PASS / PASS WITH WARNINGS / FAIL

Files created (test + 2 docs)
Files modified (6 validation/ modules + 6 test_taskN files + Task 4-9 docs)
Commands run and results

Test summary:
- prior + Phases 1-8 new tests: ~183
- Phase 9 new realignment tests: ~6-11
- Phase 9 fixture updates: existing 79 Task 4-9 tests still passing
- total: ~190+

Per-Task realignment evidence:
- Task 4: M1 normalized receptor+MYO1D consumed as inputs; pass-through PDB no longer emitted
- Task 5: M1 output verification works
- Task 6: job specs reference M1 normalized receptor + Task 4 derivations
- Task 7: residue identity validates against M1 mapping CSV
- Task 8: pocket plan references M1 dockable
- Task 9: protomer_id resolved via M1 mapping CSV

End-to-end pipeline test:
- M1 prepare-inputs + Tasks 4-9 commands all complete in one run dir
- No duplicate normalization
- Legacy prepared/egfr/egfr_receptor_normalized.pdb NOT emitted

Old workflow protection: empty diff

Acceptance closure:
- M1 acceptance: 15/15 unchanged from Phase 8
- Tasks 4-9: now consume M1 canonical outputs
- Module tree fully aligned with handoff §4 spec

Known limitations:
- M2 actual execution (PyRosetta/Vina/fpocket runs) still out of scope
- Real receptor/MYO1D/ligand placement still user-side
- HPC qsub still user-side

Next: Milestone 2 actual execution implementation (M2.1 PPI input generation, M2.2 PyRosetta adapter, ...)
```
