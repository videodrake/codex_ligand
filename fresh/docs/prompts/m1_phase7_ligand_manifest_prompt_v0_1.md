# Claude M1 Phase 7 Prompt — Ligand Manifest Shell v0.1

Branch `claude/task10`. Phases 1-6 complete. This is **M1 Phase 7** — implements the ligand manifest/confidentiality shell per `milestone1_foundation_codex_handoff_v0_5.md` §17, closing M1 §23 #14.

## 1. Project context

The project uses three confidential MYO1D activity-associated chemical anchors with public IDs `Cpd-A`, `Cpd-B`, `Cpd-C` and internal IDs (e.g., 173940, 97806, VAX-C12_0). Public outputs must use only public IDs; internal mapping lives in `fresh/data/private/compound_id_map.csv` (gitignored).

M1 does NOT run ligand docking. This phase only creates the manifest/shell so that:
- Future Vina/PyRosetta runs (M3) can hash and locate ligand structures.
- Public outputs never leak internal IDs.
- Missing ligand files are reported cleanly without crashing the pipeline.

## 2. Absolute rules

Do not modify the old workflow. Maintain Py2/3 syntax compatibility.

Source-of-truth values from `fresh/configs/fresh_run.yaml`:

```yaml
ligands:
  public_ids: [Cpd-A, Cpd-B, Cpd-C]
  confidentiality: public_outputs_use_public_ids_only
```

`fresh/configs/paths.yaml`:

```yaml
raw_ligands: fresh/data/raw/ligands
private_data: fresh/data/private
```

No hardcoding inside the module.

## 3. Scope

In scope:
- Create `fresh/src/egfr_myo1d/ligand/__init__.py`
- Create `fresh/src/egfr_myo1d/ligand/manifest.py`
- Add `manifest-ligands` CLI subcommand
- Tests under `fresh/tests/test_m1_phase7_ligand_manifest.py` (≥8 tests)
- Fixtures under `fresh/tests/fixtures/m1_phase7_ligand/` — minimal valid SDF placeholders for Cpd-A/B/C and a synthetic private mapping
- Docs `fresh/docs/m1_phase7_ligand_manifest.md` and `m1_phase7_changes.md`

Out of scope:
- Ligand structure preparation (PDBQT, RDKit, Open Babel — that's M3.2)
- Ligand docking (M3.3+)
- Modifying real `fresh/data/private/compound_id_map.csv` (it's gitignored; user owns it)
- Modifying old workflow files

## 4. Required CLI behavior

```bash
python -m egfr_myo1d.cli manifest-ligands \
  --run-id RUN \
  [--ligands-dir PATH]                 # default fresh/data/raw/ligands
  [--private-mapping PATH]              # default fresh/data/private/compound_id_map.csv
  [--profile codex_dev|hpc_strict]
  [--mode smoke_env|smoke_input]
  [--compound-stage-enabled true|false]  # default false; controls whether missing ligands FAIL or WARN
```

Behavior:
- Read public IDs from fresh_run.yaml (Cpd-A, Cpd-B, Cpd-C)
- For each public ID, look for `<ligands-dir>/<public_id>.sdf` (also accept .mol, .mol2, .pdb)
- Hash present files with sha256
- Read private mapping if present (read-only); validate schema
- Write manifest CSV under runs/<run_id>/qc/
- Public outputs: ONLY public_ids. Internal IDs MUST NOT appear in any file under runs/<run_id>/.

## 5. Files to create / modify

Create:

```text
fresh/src/egfr_myo1d/ligand/__init__.py
fresh/src/egfr_myo1d/ligand/manifest.py
fresh/tests/test_m1_phase7_ligand_manifest.py
fresh/tests/fixtures/m1_phase7_ligand/Cpd-A.sdf
fresh/tests/fixtures/m1_phase7_ligand/Cpd-B.sdf
fresh/tests/fixtures/m1_phase7_ligand/Cpd-C.sdf
fresh/tests/fixtures/m1_phase7_ligand/compound_id_map.csv      # synthetic, NOT gitignored as it's a test fixture with placeholder IDs
fresh/docs/m1_phase7_ligand_manifest.md
fresh/docs/m1_phase7_changes.md
```

The fixture SDF files must be minimal valid SDF (e.g., 5 atoms of methane-like H/C; not real chemistry). The fixture compound_id_map.csv uses fake placeholder internals like `INTERNAL_TEST_173940` — it must NOT contain real internal IDs.

Modify:

```text
fresh/src/egfr_myo1d/cli.py             # add manifest-ligands subparser + handler
.gitignore                              # verify fresh/data/private/compound_id_map.csv pattern; add if missing
```

## 6. Public API

`ligand/manifest.py`:

```python
def build_ligand_manifest(ctx, public_ids, ligands_dir, private_mapping_path, profile, compound_stage_enabled):
    # type: (RunContext, list[str], Path, Path | None, str, bool) -> LigandManifest
    """
    For each public_id in public_ids:
      Look for <ligands_dir>/<public_id>.{sdf,mol,mol2,pdb}
      If present: record sha256, format, file size
      If absent: record exists=false; severity per profile/compound_stage_enabled
    Read private_mapping_path if present; validate schema (public_id,internal_id,notes).
    Write qc/ligand_manifest_qc.csv (public IDs only).
    Verify no internal IDs appear in any output file.
    Append phase status. Return LigandManifest report.
    """

LIGAND_MANIFEST_QC_COLUMNS = [
    "public_id", "expected_path", "exists", "format",
    "sha256_or_missing", "size_bytes_or_missing", "status", "notes"
]
```

`LigandManifest` dataclass:

```text
public_ids: list[str]
present_count: int
missing_count: int
ligand_files: list[dict]               # public_id -> {path, exists, sha256, format, ...}
private_mapping_present: bool
private_mapping_schema_valid: bool | None
internal_ids_leaked: bool              # MUST be false
status: "PASS" | "WARN" | "FAIL"
warnings: list[str]
output_qc_csv: Path
```

## 7. Required output files

After `manifest-ligands` runs:

```text
fresh/runs/<run_id>/qc/ligand_manifest_qc.csv
fresh/runs/<run_id>/manifest/ligand_manifest_report.json
fresh/runs/<run_id>/logs/phase_status.jsonl              # appended
fresh/runs/<run_id>/logs/master.log                      # appended
```

QC CSV columns:

```csv
public_id,expected_path,exists,format,sha256_or_missing,size_bytes_or_missing,status,notes
```

Report JSON:

```json
{
  "run_id": "...",
  "public_ids": ["Cpd-A", "Cpd-B", "Cpd-C"],
  "ligands_dir": "...",
  "private_mapping_present": true,
  "private_mapping_schema_valid": true,
  "internal_ids_leaked_into_outputs": false,
  "ligands": [
    {"public_id": "Cpd-A", "expected_path": "...", "exists": true, "format": "sdf", "sha256": "...", "size_bytes": 1234, "status": "PASS"},
    ...
  ],
  "compound_stage_enabled": false,
  "status": "PASS",
  "warnings": []
}
```

## 8. Behavior policy

```text
- Always use public IDs in all output files. Never write internal IDs to runs/<run_id>/.
- Verify internal_ids_leaked=false by scanning the QC CSV and report JSON for any private mapping internal_id strings; FAIL if any leak.
- Private mapping schema: 3 columns (public_id, internal_id, notes). Validate column names. If schema invalid, status=WARN, message recorded.
- Missing ligand files:
  - codex_dev + compound_stage_enabled=false: status=WARN, OK to proceed
  - codex_dev + compound_stage_enabled=true:  status=WARN
  - hpc_strict + compound_stage_enabled=false: status=WARN
  - hpc_strict + compound_stage_enabled=true:  status=FAIL
- Accept multiple ligand formats: sdf, mol, mol2, pdb. The "format" field reflects the actual file extension found. If multiple formats exist for the same public_id, prefer .sdf, then .mol, then .mol2, then .pdb; report alternatives in notes.
- Do NOT modify ligand files in any way (no parsing beyond hash + format detection).
- private_mapping path may be absent: status=WARN with message "private_mapping_not_present_or_unreadable".
```

## 9. Severity rules

```text
PASS:  all public IDs found, private mapping present and valid, zero internal-ID leaks
WARN:  some ligand files missing while compound_stage_enabled=false; private mapping absent or schema invalid; multiple formats per public_id
FAIL:  any internal ID leaked into outputs; ligand files missing while hpc_strict + compound_stage_enabled=true; path-traversal run_id; write outside run_dir
```

## 10. Tests required (≥8)

```text
test_ligand_manifest_uses_only_public_ids_in_qc_csv
test_ligand_internal_ids_never_appear_in_public_outputs
test_ligand_files_present_records_sha256
test_ligand_files_missing_reported_cleanly_in_codex_dev
test_ligand_files_missing_in_hpc_strict_with_compound_stage_enabled_fails
test_ligand_files_missing_in_hpc_strict_without_compound_stage_enabled_warns
test_private_mapping_schema_validated_when_present
test_private_mapping_absent_warns
test_internal_id_leak_detection_fails_status
test_ligand_manifest_writes_under_run_dir_only
test_path_traversal_run_id_rejected
test_cli_help_includes_manifest_ligands
```

(12 tests; ≥8 required.)

## 11. Acceptance commands

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"

pytest -q fresh/tests/test_m1_phase7_ligand_manifest.py
pytest -q fresh/tests

python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase7_local

# All ligands present + private mapping present
python -m egfr_myo1d.cli manifest-ligands --run-id m1_phase7_local --ligands-dir fresh/tests/fixtures/m1_phase7_ligand --private-mapping fresh/tests/fixtures/m1_phase7_ligand/compound_id_map.csv --profile codex_dev

# Inspect QC CSV
head -5 fresh/runs/m1_phase7_local/qc/ligand_manifest_qc.csv
grep -E "INTERNAL|173940|97806|VAX" fresh/runs/m1_phase7_local/qc/ligand_manifest_qc.csv && echo "LEAK!" || echo "No internal IDs leaked"

# Missing ligands in codex_dev (no compound stage)
python -m egfr_myo1d.cli manifest-ligands --run-id m1_phase7_local --ligands-dir /tmp/no_such_dir --profile codex_dev

# Missing ligands in hpc_strict + compound stage = FAIL
python -m egfr_myo1d.cli manifest-ligands --run-id m1_phase7_local --ligands-dir /tmp/no_such_dir --profile hpc_strict --compound-stage-enabled true || echo "Expected FAIL"

# Path traversal
python -m egfr_myo1d.cli manifest-ligands --run-id ../bad_run

# Old workflow protection
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/

# Verify gitignore
grep -E "fresh/data/private|compound_id_map" .gitignore
```

## 12. Final response format

```text
M1 Phase 7 status: PASS / PASS WITH WARNINGS / FAIL
Files created (incl. fixture SDFs and synthetic compound_id_map.csv)
Files modified (cli.py, .gitignore if needed)
Commands run and results
Test summary
Public-ID-only enforcement: confirmed (no internal ID leaks)
Profile/stage-aware severity matrix verified:
- codex_dev / stage=false / missing: WARN
- codex_dev / stage=true / missing:  WARN
- hpc_strict / stage=false / missing: WARN
- hpc_strict / stage=true / missing:  FAIL
SDF file hashing: sha256 recorded
Private mapping schema validation
Acceptance closure: M1 §23 #14 closed
Old workflow protection: empty diff
Known limitations:
- No PDBQT prep (M3.2)
- No actual ligand docking (M3.3+)
- No fresh compound library expansion (post-M3)
```
