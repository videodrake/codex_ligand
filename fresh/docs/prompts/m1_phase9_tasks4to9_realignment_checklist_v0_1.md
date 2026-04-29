# M1 Phase 9 Acceptance Checklist v0.1 — Tasks 4-9 Schema Realignment

Use this after the implementer applies M1 Phase 9.

## 1. Pre-Phase state preserved

```text
Old workflow files unchanged.
Phases 1-8 (M1 foundation) modules unchanged.
Existing Task 4-9 SEMANTICS preserved (only INPUT path/schema changes).
```

## 2. Modified validation modules importable

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -c "from egfr_myo1d.validation import prepared_inputs, real_inputs, ppi_sampling_plan, ppi_consensus, pocket_discovery, pocket_candidate_prioritization; print('OK')"
```

## 3. CLI surface unchanged

```bash
python -m egfr_myo1d.cli --help
```

Expected: same set of subcommands as after Phase 8 (no new commands; possibly new optional `--auto-prepare` flag on Task 4-9 commands).

## 4. End-to-end M1 + Tasks 4-9 walk

```bash
python -m egfr_myo1d.cli init-run --mode smoke_input --run-id m1_phase9_local
python -m egfr_myo1d.cli prepare-inputs --run-id m1_phase9_local --mode smoke_input --profile codex_dev --input-root fresh/tests/fixtures/m1_phase8_integration

# Verify M1 normalized outputs exist
test -f fresh/runs/m1_phase9_local/normalized/receptors/EGFR_160-185_dockable_669_1014_explicit_AB.pdb
test -f fresh/runs/m1_phase9_local/normalized/myo1d/MYO1D_955_1006.pdb
test -f fresh/runs/m1_phase9_local/qc/EGFR_160-185_receptor_mapping.csv
test -f fresh/runs/m1_phase9_local/manifest/membrane_frame.json

# Now run Tasks 4-9 against the same run_dir
python -m egfr_myo1d.cli prepare-ppi-inputs --run-id m1_phase9_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task4_inputs
python -m egfr_myo1d.cli validate-real-inputs --run-id m1_phase9_local --mode smoke_input --profile codex_dev --input-root fresh/data/raw  # expects M1 outputs to be consulted
python -m egfr_myo1d.cli plan-ppi-sampling --run-id m1_phase9_local --mode smoke_input --profile codex_dev --input-root fresh/tests/fixtures/task4_inputs
python -m egfr_myo1d.cli summarize-ppi-consensus --run-id m1_phase9_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task7_ppi_consensus
python -m egfr_myo1d.cli plan-pocket-discovery --run-id m1_phase9_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task8_pocket_planning
python -m egfr_myo1d.cli prioritize-pocket-candidates --run-id m1_phase9_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task9_pocket_candidates
```

Expected: all commands exit 0 (PASS or PASS_WITH_WARNINGS).

## 5. No duplicate normalization

```bash
# Legacy pass-through receptor PDB should NOT exist
test ! -f fresh/runs/m1_phase9_local/prepared/egfr/egfr_receptor_normalized.pdb && echo "OK: pass-through gone"

# Task 4 manifest should reference M1 normalized output
python -c "
import json
m = json.load(open('fresh/runs/m1_phase9_local/manifest/prepared_input_manifest.json'))
inputs = m.get('inputs', {})
print('inputs:', list(inputs.keys()))
"
```

Expected: legacy file absent; manifest references M1 normalized paths.

## 6. Task 4 audit references M1 paths

```bash
head -5 fresh/runs/m1_phase9_local/qc/egfr_receptor_normalization_audit.csv
```

Expected: audit rows reference M1 normalized PDB as the source, not the raw fixture.

(Or, if Task 4 audit format does not directly expose source path, the manifest does.)

## 7. Task 6 job specs reference M1 paths

```bash
head -3 fresh/runs/m1_phase9_local/prepared/ppi/ppi_job_specs.jsonl
```

Expected: each spec's `receptor_path` field references `normalized/receptors/...runtime_offset_receptor_only.pdb`. `partner_path` references `prepared/myo1d/MYO1D_sheet8_9_12_core_955_1001.pdb` (Task 4 derivation) or `normalized/myo1d/MYO1D_955_1006.pdb`.

## 8. Task 7 mapping consumption

```bash
python -c "
import csv
# Verify Task 7 consensus output records protomer_id and source_resseq from mapping
with open('fresh/runs/m1_phase9_local/output/ppi/ppi_consensus_patch.csv') as f:
    reader = csv.DictReader(f)
    sample = next(reader, None)
    print('fields:', list(sample.keys()) if sample else 'empty')
"
```

Expected: fields include protomer_id and either source_resseq or runtime_resseq with explicit reference to M1 mapping.

## 9. Missing M1 outputs handled

```bash
# Run a Task 4-9 command in a run that has NOT been prepared by M1
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase9_no_prep
python -m egfr_myo1d.cli prepare-ppi-inputs --run-id m1_phase9_no_prep --mode smoke_env --profile hpc_strict --input-root fresh/tests/fixtures/task4_inputs
```

Expected (per implementer's choice):
- (Strict) exit 1 with clear message about missing M1 outputs
- (Permissive with --auto-prepare) auto-invokes prepare-inputs first
- (Permissive without flag in codex_dev) WARN and continues using raw fixture inputs (legacy fallback)

The behavior must be documented in `fresh/docs/m1_phase9_tasks4to9_realignment.md`.

## 10. Tests

```bash
pytest -q fresh/tests/test_m1_phase9_tasks_realignment.py
pytest -q fresh/tests/test_task4_ppi_input_preparation.py
pytest -q fresh/tests/test_task5_real_input_readiness.py
pytest -q fresh/tests/test_task6_ppi_sampling_plan.py
pytest -q fresh/tests/test_task7_ppi_consensus_patch.py
pytest -q fresh/tests/test_task8_pocket_discovery_plan.py
pytest -q fresh/tests/test_task9_pocket_candidate_prioritization.py
pytest -q fresh/tests
```

Expected: all pass.

## 11. Path traversal (sanity check)

```bash
python -m egfr_myo1d.cli prepare-ppi-inputs --run-id ../bad_run --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task4_inputs
```

Nonzero exit, no outside writes.

## 12. Old workflow protection

```bash
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

Empty.

## 13. Documentation updated

```bash
grep -E "M1 dependency|consume M1 normalized" fresh/docs/task4_ppi_input_preparation.md fresh/docs/task5_real_input_readiness.md fresh/docs/task6_ppi_sampling_plan.md fresh/docs/task7_ppi_consensus_patch.md fresh/docs/task8_pocket_discovery_plan.md fresh/docs/task9_pocket_candidate_prioritization.md
```

Expected: each Task 4-9 doc has an explicit note about consuming M1 canonical outputs.

`fresh/docs/m1_phase9_tasks4to9_realignment.md` and `fresh/docs/m1_phase9_changes.md` exist and explain the changes.

## 14. What must not be in this phase

```text
- modifying Phases 1-8 module logic
- M2 actual execution (PyRosetta/Vina/fpocket runs)
- M3 work (compound docking, scoring)
- modifying old workflow files
- deleting Task 4-9 tests (only updating fixtures/assertions)
- removing Task 4-9 derivations like MYO1D 955-1001 active-face slice (these are valid Task 4 derivations of M1 raw 955-1006)
```

## 15. Phase 9 accepted if

```text
- Tasks 4-9 modules consume M1 normalized outputs as inputs.
- prepared/egfr/egfr_receptor_normalized.pdb (pass-through) no longer emitted (or replaced by manifest reference).
- Task 4 derivations (MYO1D 955-1001 core, 955-1006 tail-masked) preserved.
- Task 5 verifies M1 outputs exist.
- Task 6 job specs reference M1 normalized receptor.
- Task 7 residue identity validates against M1 mapping CSV.
- Task 8 pocket plan references M1 dockable receptor.
- Task 9 protomer_id resolved via M1 mapping CSV.
- Existing 79 Task 4-9 tests pass after fixture updates.
- ≥6 new realignment tests pass.
- End-to-end M1 + Tasks 4-9 walk in one run dir succeeds.
- Module tree fully aligned with handoff §4.
- Old workflow files unmodified.
```

## 16. Implementer final response must include

```text
M1 Phase 9 status: PASS / PASS WITH WARNINGS / FAIL
Files created (test + 2 docs)
Files modified (6 validation/ modules + 6 test_taskN files + Task 4-9 doc updates)
Commands run and results
Test summary
Per-Task realignment evidence (Task 4 / 5 / 6 / 7 / 8 / 9)
End-to-end pipeline evidence: M1 prepare-inputs + Tasks 4-9 commands all complete in same run dir; no duplicate normalization
Missing-M1-outputs behavior: documented (strict/permissive/auto-prepare)
Acceptance closure: M1 closure scorecard unchanged 15/15; module tree aligned with handoff §4
Old workflow protection: empty diff
Known limitations
Next: Milestone 2 actual execution implementation
```
