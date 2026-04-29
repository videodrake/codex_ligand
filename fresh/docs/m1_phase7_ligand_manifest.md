# M1 Phase 7 — Ligand Manifest Shell

Closes M1 §23 #14 (ligand manifest shell exists without exposing private IDs by default) per `milestone1_foundation_codex_handoff_v0_5.md` §17.

## What it does

For each public ligand ID resolved from `fresh/configs/fresh_run.yaml` (`Cpd-A`, `Cpd-B`, `Cpd-C`):

1. Search `<ligands_dir>` for a file named `<public_id>.{sdf|mol|mol2|pdb}` (priority sdf > mol > mol2 > pdb).
2. If present, record `sha256` + `format` + `size_bytes` + alternative formats.
3. If absent, severity is determined by the profile/stage matrix below.

The private mapping at `<private_mapping_path>` (default `fresh/data/private/compound_id_map.csv`) is read **only** for internal-ID leak detection. Its content (internal_ids) is never copied into any run output. The 3-column schema `(public_id, internal_id, notes)` is validated.

After writing the QC CSV and report JSON, the module re-reads those files and scans for any internal_id substring. If any internal_id appears in any output file, status is escalated to `FAIL` and the leak is recorded.

## CLI

```bash
python -m egfr_myo1d.cli manifest-ligands \
    --run-id RUN \
    [--ligands-dir PATH] \
    [--private-mapping PATH] \
    [--profile codex_dev|hpc_strict] \
    [--mode smoke_env|smoke_input] \
    [--compound-stage-enabled true|false]
```

Defaults are resolved from `fresh/configs/paths.yaml`:
- `--ligands-dir` → `raw_ligands` (`fresh/data/raw/ligands`)
- `--private-mapping` → `<private_data>/compound_id_map.csv` (`fresh/data/private/compound_id_map.csv`)

## Module additions

```text
fresh/src/egfr_myo1d/ligand/__init__.py    new
fresh/src/egfr_myo1d/ligand/manifest.py    new
```

Public API:

```python
SUPPORTED_FORMATS = ("sdf", "mol", "mol2", "pdb")
PRIVATE_MAPPING_REQUIRED_COLUMNS = ("public_id", "internal_id", "notes")
LIGAND_MANIFEST_QC_COLUMNS = [...]   # 8 columns

LigandFileRecord (dataclass)
LigandManifest   (dataclass)
load_public_ids(ctx) -> list[str]
load_default_paths(ctx) -> (raw_ligands, private_mapping)
build_ligand_manifest(ctx, public_ids=None, ligands_dir=None,
                       private_mapping_path=None, profile="codex_dev",
                       compound_stage_enabled=False) -> LigandManifest
```

## Severity matrix (handoff §17)

| profile | compound_stage_enabled | missing files | status |
| --- | --- | --- | --- |
| `codex_dev` | `false` | yes | `WARN` |
| `codex_dev` | `true`  | yes | `WARN` |
| `hpc_strict`| `false` | yes | `WARN` |
| `hpc_strict`| `true`  | yes | `FAIL` |

Other escalators:
- Internal-ID leak detected in outputs → `FAIL` (no profile/stage gate)
- Private mapping schema invalid → `WARN` (informational)
- Private mapping absent → `WARN`

## Outputs

```text
fresh/runs/<run_id>/qc/ligand_manifest_qc.csv
fresh/runs/<run_id>/manifest/ligand_manifest_report.json
fresh/runs/<run_id>/logs/phase_status.jsonl    (appended)
fresh/runs/<run_id>/logs/master.log            (appended)
```

QC CSV columns:

```csv
public_id,expected_path,exists,format,sha256_or_missing,size_bytes_or_missing,status,notes
```

Manifest JSON content (excerpt):

```json
{
  "public_ids": ["Cpd-A", "Cpd-B", "Cpd-C"],
  "ligands_dir": "...",
  "private_mapping_path": "...",
  "private_mapping_present": true,
  "private_mapping_schema_valid": true,
  "internal_ids_leaked_into_outputs": false,
  "leaked_internal_ids": [],
  "compound_stage_enabled": false,
  "ligands": [
    {"public_id": "Cpd-A", "exists": true, "format": "sdf",
     "sha256": "...", "size_bytes": 412, "alternative_formats": [],
     "status": "PASS", "notes": ""}
  ],
  "present_count": 3,
  "missing_count": 0,
  "status": "PASS",
  "warnings": [],
  "profile": "codex_dev",
  "score_bonus_allowed": false,
  "timestamp": "..."
}
```

## Behavior policy (handoff §17)

- **Public IDs only** in any run output. The private mapping is read once into memory for leak detection and never serialized.
- **Internal-ID leak detection** scans the QC CSV and report JSON after writing. Any substring match against any internal_id from the private mapping triggers FAIL and the leak is recorded in `leaked_internal_ids`.
- M1 does **not** run ligand docking — this is shell only. PDBQT prep, RDKit/OpenBabel conversion, and Vina docking are M3.
- Private mapping at `fresh/data/private/compound_id_map.csv` is gitignored. The .gitignore rule was set in Task 1; Phase 7 includes a regression test.

## Reusable fixtures

```text
fresh/tests/fixtures/m1_phase7_ligand/Cpd-A.sdf       (3-atom placeholder)
fresh/tests/fixtures/m1_phase7_ligand/Cpd-B.sdf       (4-atom placeholder)
fresh/tests/fixtures/m1_phase7_ligand/Cpd-C.sdf       (5-atom placeholder)
fresh/tests/fixtures/m1_phase7_ligand/compound_id_map.csv (synthetic INTERNAL_TEST_PLACEHOLDER_*)
```

The synthetic compound_id_map.csv contains placeholder internal IDs (`INTERNAL_TEST_PLACEHOLDER_A/B/C`). Tests assert these never appear in any run output.

## What is intentionally not in this phase

- Ligand structure preparation (PDBQT generation, RDKit parsing, OpenBabel conversion) — that's M3.2
- Ligand docking (Vina, etc.) — that's M3.3+
- Modifying the real `fresh/data/private/compound_id_map.csv` (user owns it)
- Exposing real internal compound IDs anywhere
- Modifying old workflow files
