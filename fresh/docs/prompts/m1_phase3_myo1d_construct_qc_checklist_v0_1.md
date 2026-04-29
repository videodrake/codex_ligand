# M1 Phase 3 Acceptance Checklist v0.1 — MYO1D Construct + QC

Use this after the implementer applies M1 Phase 3.

## 1. Confirm pre-Phase state preserved

```text
Old workflow files unchanged.
Phase 1 cleanup (core/cleanup.py, cleanup CLI, scripts/cleanup_run.py) unchanged.
Phase 2 relocation (myo1d/construct.py, myo1d/pdb_writer.py) intact and Phase 3 only EXTENDS construct.py.
preparation/{masks.py, restraints.py} unchanged.
```

## 2. New module + extension

```bash
test -f fresh/src/egfr_myo1d/myo1d/qc.py
python -c "from egfr_myo1d.myo1d.construct import slice_myo1d_construct, emit_myo1d_construct_pdb; print('OK')"
python -c "from egfr_myo1d.myo1d.qc import run_myo1d_qc; print('OK')"
```

## 3. CLI registered

```bash
python -m egfr_myo1d.cli --help | grep prepare-myo1d
python -m egfr_myo1d.cli prepare-myo1d --help
```

Expected: `prepare-myo1d` listed in main help; subparser help shows `--run-id`, `--source`, `--construct`, `--profile`, `--mode`.

## 4. Primary construct fixture (955-1006, codex_dev)

```bash
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase3_local
python -m egfr_myo1d.cli prepare-myo1d \
  --run-id m1_phase3_local \
  --source fresh/tests/fixtures/m1_phase3_myo1d/myo1d_955_1006_valid.pdb \
  --construct 955-1006 \
  --profile codex_dev
```

Expected:
- exit 0
- file `fresh/runs/m1_phase3_local/normalized/myo1d/MYO1D_955_1006.pdb` exists
- file `fresh/runs/m1_phase3_local/qc/myo1d_construct_qc.csv` exists with column header line matching prompt §7
- file `fresh/runs/m1_phase3_local/manifest/myo1d_construct_manifest.json` exists with source + output sha256
- output PDB residues span 955-1006 with original numbering preserved
- key_sheet8_present, key_sheet9_present, key_sheet12_present all true
- status: PASS

## 5. Comparator construct fixture (955-1001)

```bash
python -m egfr_myo1d.cli prepare-myo1d \
  --run-id m1_phase3_local \
  --source fresh/tests/fixtures/m1_phase3_myo1d/myo1d_955_1001_short.pdb \
  --construct 955-1001 \
  --profile codex_dev
```

Expected:
- exit 0
- output PDB spans 955-1001
- c_watch_present is false or partial (1002-1006 absent by definition of construct)
- status: PASS or WARN with explanatory note

## 6. Terminal-artifact negative fixture (962-1006, codex_dev)

```bash
python -m egfr_myo1d.cli prepare-myo1d \
  --run-id m1_phase3_local \
  --source fresh/tests/fixtures/m1_phase3_myo1d/myo1d_962_1006_terminal_bad.pdb \
  --construct 962-1006 \
  --profile codex_dev
```

Expected:
- exit 0 (WARN does not fail)
- status: WARN
- warnings include `myo1d_962_start_terminal_artifact`
- output PDB written under run_dir, but downstream consumers are advised by the warning

## 7. Terminal-artifact negative fixture (962-1006, hpc_strict)

```bash
python -m egfr_myo1d.cli prepare-myo1d \
  --run-id m1_phase3_local \
  --source fresh/tests/fixtures/m1_phase3_myo1d/myo1d_962_1006_terminal_bad.pdb \
  --construct 962-1006 \
  --profile hpc_strict
```

Expected:
- exit 1 (FAIL in strict)
- status: FAIL
- no normalized PDB emitted (or emitted but tagged blocked)
- manifest records the blocker

## 8. Cap-HETATM preservation fixture

```bash
python -m egfr_myo1d.cli prepare-myo1d \
  --run-id m1_phase3_local \
  --source fresh/tests/fixtures/m1_phase3_myo1d/myo1d_with_ace_nme_caps.pdb \
  --construct 955-1006 \
  --profile codex_dev
```

Expected:
- ACE / NME records in output PDB
- ace_nme_caps_present: true in QC
- standard AA written as HETATM in source (e.g., ILE1000) is preserved as biological residue, not dropped

## 9. Missing source in codex_dev

```bash
python -m egfr_myo1d.cli prepare-myo1d \
  --run-id m1_phase3_local \
  --source fresh/data/raw/myo1d/does_not_exist.pdb \
  --construct 955-1006 \
  --profile codex_dev
```

Expected: exit 0, status WARN, manifest lists missing_required_inputs.

## 10. Missing source in hpc_strict

```bash
python -m egfr_myo1d.cli prepare-myo1d \
  --run-id m1_phase3_local \
  --source fresh/data/raw/myo1d/does_not_exist.pdb \
  --construct 955-1006 \
  --profile hpc_strict
```

Expected: exit 1, status FAIL.

## 11. Score bonus enforcement

Inspect `fresh/runs/m1_phase3_local/manifest/myo1d_construct_manifest.json` or the QC report. The implementation must read `key_residue_bonus_weight` from `fresh/configs/gates.yaml` and assert/record it as `0.0`. No hardcoding inside the module.

```bash
grep -n "key_residue_bonus_weight" fresh/src/egfr_myo1d/myo1d/qc.py
```

Expected: any reference uses the gates.yaml key, not a literal `0.0` written into Python.

## 12. Path traversal safety

```bash
python -m egfr_myo1d.cli prepare-myo1d --run-id ../bad_run --source fresh/tests/fixtures/m1_phase3_myo1d/myo1d_955_1006_valid.pdb --construct 955-1006
```

Expected: nonzero exit, no outside directory created.

## 13. Tests

```bash
pytest -q fresh/tests/test_m1_phase3_myo1d_construct_qc.py
pytest -q fresh/tests
```

Required tests pass (≥10). Existing 98 + Phase 1 (8) + Phase 3 new tests all pass.

## 14. Old workflow protection

```bash
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

Expected: empty.

## 15. What must not be in this phase

```text
- receptor normalization (Phase 4)
- membrane frame computation (Phase 5)
- ligand work (Phase 7)
- integration orchestrator (Phase 8)
- Task 4-9 schema realignment (Phase 9)
- any docking, scoring, candidate
```

## 16. Phase 3 accepted if

```text
- myo1d/qc.py exists, run_myo1d_qc public.
- myo1d/construct.py extended with slice_*, emit_* without breaking Task 4 callers.
- prepare-myo1d CLI subcommand registered.
- normalized/myo1d/MYO1D_955_1006.pdb emitted under run_dir.
- qc/myo1d_construct_qc.csv columns match handoff §16.
- manifest/myo1d_construct_manifest.json includes sha256.
- 955-1006 codex_dev = PASS, 955-1001 comparator = PASS/WARN, 962-1006 codex_dev = WARN, 962-1006 hpc_strict = FAIL.
- ACE/NME caps preserved.
- key_residue_bonus_weight read from gates.yaml.
- ≥10 phase tests pass; existing 98 + Phase 1 (8) tests still pass.
- M1 §23 #13 closed.
- Old workflow files untouched.
```

## 17. Implementer final response must include

```text
M1 Phase 3 status: PASS / PASS WITH WARNINGS / FAIL
Files created
Files modified
Commands run and results
Test summary (prior + Phase 1 + Phase 3 new)
Output artifact verification
Negative-regression behavior (962-1006, missing source) profile-aware results
Cap preservation evidence
Score bonus enforcement source (gates.yaml read)
Old workflow protection
Acceptance closure: M1 §23 #13 closed
Known limitations
```
