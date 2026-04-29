# M1 Phase 7 Acceptance Checklist v0.1 — Ligand Manifest Shell

Use this after the implementer applies M1 Phase 7.

## 1. Pre-Phase state preserved

```text
Old workflow files unchanged.
Phases 1-6 outputs/modules unchanged.
fresh/data/private/compound_id_map.csv (real one, if present) NOT modified.
```

## 2. New module importable

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -c "from egfr_myo1d.ligand.manifest import build_ligand_manifest, LIGAND_MANIFEST_QC_COLUMNS; print('OK', len(LIGAND_MANIFEST_QC_COLUMNS))"
```

## 3. CLI registered

```bash
python -m egfr_myo1d.cli --help | grep manifest-ligands
python -m egfr_myo1d.cli manifest-ligands --help
```

Help text must include `--run-id`, `--ligands-dir`, `--private-mapping`, `--profile`, `--mode`, `--compound-stage-enabled`.

## 4. All ligands present + private mapping (codex_dev)

```bash
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase7_local
python -m egfr_myo1d.cli manifest-ligands \
  --run-id m1_phase7_local \
  --ligands-dir fresh/tests/fixtures/m1_phase7_ligand \
  --private-mapping fresh/tests/fixtures/m1_phase7_ligand/compound_id_map.csv \
  --profile codex_dev
```

Expected:

```text
- exit 0
- fresh/runs/m1_phase7_local/qc/ligand_manifest_qc.csv exists
- fresh/runs/m1_phase7_local/manifest/ligand_manifest_report.json exists
- report.status == "PASS"
- report.present_count == 3
- report.internal_ids_leaked_into_outputs == false
```

## 5. QC CSV schema

```bash
head -1 fresh/runs/m1_phase7_local/qc/ligand_manifest_qc.csv
```

Expected header (exact column names per prompt §7):

```text
public_id,expected_path,exists,format,sha256_or_missing,size_bytes_or_missing,status,notes
```

Each row should have `public_id` ∈ {Cpd-A, Cpd-B, Cpd-C}; sha256 for present files; status PASS or WARN.

## 6. No internal-ID leak

```bash
# Synthetic fixture uses placeholder INTERNAL_TEST_* internals
grep -E "INTERNAL_TEST_|173940|97806|VAX" fresh/runs/m1_phase7_local/qc/ligand_manifest_qc.csv \
  fresh/runs/m1_phase7_local/manifest/ligand_manifest_report.json && echo "FAIL: leak" || echo "OK: no leak"
```

Expected: "OK: no leak". (Even synthetic INTERNAL_TEST_ identifiers must not appear in run outputs.)

## 7. Missing ligands in codex_dev

```bash
python -m egfr_myo1d.cli manifest-ligands \
  --run-id m1_phase7_local \
  --ligands-dir /tmp/does_not_exist \
  --profile codex_dev
```

Expected:

```text
- exit 0 (WARN)
- report.status == "WARN"
- report.missing_count == 3
- per-ligand status == "WARN"
- notes describe missing path
```

## 8. Missing ligands in hpc_strict + compound stage = FAIL

```bash
python -m egfr_myo1d.cli manifest-ligands \
  --run-id m1_phase7_local \
  --ligands-dir /tmp/does_not_exist \
  --profile hpc_strict \
  --compound-stage-enabled true
```

Expected:

```text
- exit 1 (FAIL)
- report.status == "FAIL"
```

## 9. Missing ligands in hpc_strict + no compound stage = WARN

```bash
python -m egfr_myo1d.cli manifest-ligands \
  --run-id m1_phase7_local \
  --ligands-dir /tmp/does_not_exist \
  --profile hpc_strict \
  --compound-stage-enabled false
```

Expected:

```text
- exit 0 (WARN)
- report.status == "WARN"
```

## 10. Private mapping absent

```bash
python -m egfr_myo1d.cli manifest-ligands \
  --run-id m1_phase7_local \
  --ligands-dir fresh/tests/fixtures/m1_phase7_ligand \
  --private-mapping /tmp/no_such_mapping.csv \
  --profile codex_dev
```

Expected:

```text
- exit 0 (WARN)
- report.private_mapping_present == false
- warnings include "private_mapping_not_present_or_unreadable"
- still no internal IDs in outputs (because mapping wasn't read)
```

## 11. .gitignore protects private mapping

```bash
grep -E "fresh/data/private|compound_id_map" .gitignore
```

Expected: at least one rule covering `fresh/data/private/*` or `fresh/data/private/compound_id_map.csv`.

```bash
# Try staging the private mapping path; should be ignored
mkdir -p fresh/data/private
echo "public_id,internal_id,notes" > fresh/data/private/compound_id_map.csv
git status --short fresh/data/private/compound_id_map.csv
rm fresh/data/private/compound_id_map.csv
```

Expected: empty status output (file ignored).

## 12. Path traversal

```bash
python -m egfr_myo1d.cli manifest-ligands --run-id ../bad_run
```

Nonzero exit, no outside writes.

## 13. Tests

```bash
pytest -q fresh/tests/test_m1_phase7_ligand_manifest.py
pytest -q fresh/tests
```

## 14. Old workflow protection

```bash
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

Empty.

## 15. What must not be in this phase

```text
- ligand structure preparation (PDBQT, RDKit, Open Babel) — that's M3.2
- ligand docking — that's M3.3+
- modifying real fresh/data/private/compound_id_map.csv
- exposing real internal compound IDs anywhere
- modifying old workflow files
```

## 16. Phase 7 accepted if

```text
- ligand/__init__.py and ligand/manifest.py created.
- manifest-ligands CLI subcommand registered.
- 3 fixture SDFs (Cpd-A/B/C) and synthetic compound_id_map.csv created.
- QC CSV columns match prompt §7.
- Public-IDs-only enforced; zero internal ID leaks (programmatically verified).
- Private mapping schema validation works; missing mapping warns.
- Profile/stage matrix correct (4 combinations per §8).
- .gitignore protects private mapping path.
- ≥8 phase tests pass; existing tests pass.
- M1 §23 #14 closed.
- Old workflow files unmodified.
```

## 17. Implementer final response must include

```text
M1 Phase 7 status: PASS / PASS WITH WARNINGS / FAIL
Files created
Files modified (cli.py, .gitignore)
Commands run and results
Test summary
Public-ID-only enforcement evidence
Profile/stage matrix evidence
.gitignore validation
Acceptance closure: M1 §23 #14 closed
Old workflow protection
Known limitations: no PDBQT prep, no docking
```
