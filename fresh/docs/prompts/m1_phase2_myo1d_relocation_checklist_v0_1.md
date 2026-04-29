# M1 Phase 2 Acceptance Checklist v0.1 — MYO1D Module Relocation

Use this after the implementer applies M1 Phase 2.

## 1. Confirm pre-Phase state preserved

Old workflow files/directories must not be modified:

```text
run_production.py
main.py
egfr_pipeline/**
config/**
docs/runbook.md
output/**
results_export/**
```

Tasks 1-9 logic must not be modified — only import statements in `validation/prepared_inputs.py` and `validation/real_inputs.py` should differ.

## 2. New module locations exist

```bash
test -f fresh/src/egfr_myo1d/myo1d/__init__.py
test -f fresh/src/egfr_myo1d/myo1d/construct.py
test -f fresh/src/egfr_myo1d/myo1d/pdb_writer.py
```

Expected: all three files exist.

## 3. Old module locations are gone

```bash
test ! -f fresh/src/egfr_myo1d/preparation/constructs.py
test ! -f fresh/src/egfr_myo1d/preparation/pdb_writer.py
```

Expected: neither file exists.

`fresh/src/egfr_myo1d/preparation/` should still contain `__init__.py`, `masks.py`, `restraints.py` (these are EGFR-side and remain).

## 4. Import audit

```bash
grep -rn "preparation.constructs\|preparation.pdb_writer" fresh/src fresh/tests
```

Expected: empty output (no remaining references to old paths).

```bash
grep -rn "from egfr_myo1d.myo1d.construct" fresh/src
grep -rn "from egfr_myo1d.myo1d.pdb_writer" fresh/src
```

Expected: matches in `validation/prepared_inputs.py` and `validation/real_inputs.py` (and possibly `myo1d/__init__.py` if it re-exports).

## 5. New modules importable

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -c "from egfr_myo1d.myo1d import construct, pdb_writer; print('OK', construct.__name__, pdb_writer.__name__)"
```

Expected:

```text
OK egfr_myo1d.myo1d.construct egfr_myo1d.myo1d.pdb_writer
```

## 6. Logic preservation (byte/symbol comparison)

Compared to the pre-relocation state, the moved files must contain the same:

```text
- function names and signatures
- module-level constants (CAP_RESNAMES, STANDARD_AA, TASK3_WARNING_CLASSES if present)
- public symbols (anything imported by validation/prepared_inputs.py or real_inputs.py)
- docstrings
- behavior on the existing test fixtures
```

The only allowed difference is the module's own internal self-imports (e.g., if `pdb_writer.py` imported a sibling, the path may need updating).

## 7. Full pytest must pass

```bash
pytest -q fresh/tests
```

Expected: all 98 prior tests pass, zero new failures, zero new warnings.

If a Task 4 test (test_task4_ppi_input_preparation.py) or Task 5 test (test_task5_real_input_readiness.py) fails, the failure is almost certainly a missed import update.

## 8. CLI smoke

```bash
python -m egfr_myo1d.cli --help
python -m egfr_myo1d.cli prepare-ppi-inputs --help
python -m egfr_myo1d.cli validate-real-inputs --help
```

Expected: all three exit 0 with help text.

## 9. End-to-end round-trip on existing fixtures

```bash
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase2_local
python -m egfr_myo1d.cli prepare-ppi-inputs --run-id m1_phase2_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task4_inputs
python -m egfr_myo1d.cli status --run-id m1_phase2_local
```

Expected: identical outputs to pre-relocation behavior. Files emitted under `fresh/runs/m1_phase2_local/prepared/` match the Task 4 schema unchanged.

## 10. Old workflow protection

```bash
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

Expected: empty output.

## 11. What must not be in this phase

```text
- new MYO1D functions in construct.py (Phase 3)
- new myo1d/qc.py module (Phase 3)
- new prepare-myo1d CLI subcommand (Phase 3)
- any logic change to moved files
- any modification to preparation/masks.py or preparation/restraints.py (EGFR-side, untouched)
- any new tests (regression net is the existing 98 tests)
```

## 12. Phase 2 accepted if

```text
- myo1d/__init__.py, myo1d/construct.py, myo1d/pdb_writer.py exist.
- preparation/constructs.py and preparation/pdb_writer.py are deleted.
- imports in validation/prepared_inputs.py and validation/real_inputs.py updated.
- no remaining references to preparation.constructs or preparation.pdb_writer.
- moved files have identical logic, docstrings, signatures, constants.
- all 98 existing tests pass.
- existing CLI commands work unchanged.
- end-to-end Task 4 fixture run produces identical output paths/files.
- old workflow files unmodified.
```

## 13. Implementer final response must include

```text
M1 Phase 2 status: PASS / FAIL
Files created (3 src + 2 docs)
Files modified (2 validation/ imports + optional preparation/__init__.py)
Files deleted (2 preparation/ files)
Test result: 98 prior tests still pass, no new tests
Import audit: grep returns zero matches for old paths
Round-trip verification: Task 4 fixture run identical to pre-relocation
Old workflow protection: empty diff
Acceptance closure: Phase 3 unblocked
Known limitations: no new functionality; structural move only
```
