# Claude M1 Phase 4 Prompt — Receptor Normalization v0.1

Branch `claude/task10`. Phases 1-3 complete. This is **M1 Phase 4** — implements EGFR receptor normalization per `milestone1_foundation_codex_handoff_v0_5.md` §14, closing M1 §23 #10 (explicit A/B normalization) and #11 (+1000 runtime offset).

## 1. Project context

The fresh workflow expects EGFR dimer receptor inputs in three states:

```text
EGFR_160-185      primary membrane-validated state
EGFR_170-200      primary membrane-validated state
3GT8_raw          crystallographic reference/control
```

Real receptor PDBs may arrive in three shapes:

```text
Case A: explicit A/B chains
Case B: same-chain duplicate (e.g., both copies stored as chain X with residue-number reset)
Case C: single-chain monomer (NOT VALID for dimer-only main analysis)
```

Currently `validation/prepared_inputs.py` (Task 4) audits the receptor (chain count, duplicate atoms, V924R-like mutation) but does NOT actually normalize: no chain split, no 669-1014 dockable crop, no +1000 runtime offset, no mapping CSV.

Phase 4 implements the actual normalization producing canonical M1 outputs.

## 2. Absolute rules

Do not modify the old workflow. Do not modify Tasks 1-9 outputs paths in this phase (Phase 9 reconciles). Maintain Py2/3 syntax compatibility.

Source-of-truth values from `fresh/configs/gates.yaml`:

```yaml
receptor:
  default_dockable_crop_start: 669
  default_dockable_crop_end: 1014
  excluded_tm_core_default: "634-668"
  runtime_offset_second_protomer: 1000
```

No hardcoding inside the module.

## 3. Scope

In scope:
- Create `fresh/src/egfr_myo1d/model/__init__.py`
- Create `fresh/src/egfr_myo1d/model/receptor_normalize.py`
- Create `fresh/src/egfr_myo1d/model/receptor_qc.py`
- Create `fresh/src/egfr_myo1d/io/residue_mapping.py`
- Add `prepare-receptor` CLI subcommand
- Tests under `fresh/tests/test_m1_phase4_receptor_normalization.py` (≥12 tests)
- Fixtures under `fresh/tests/fixtures/m1_phase4_receptor/` (or reuse mini_*.pdb where applicable)
- Docs `fresh/docs/m1_phase4_receptor_normalization.md` and `m1_phase4_changes.md`

Out of scope:
- Membrane frame generation (Phase 5)
- MYO1D work (Phase 3 already covered)
- Modifying Task 4 prepared_inputs.py output paths (Phase 9)
- M2 docking work
- V924R mutation repair (must remain WARN-only, never silently mutated)

## 4. Required CLI behavior

```bash
python -m egfr_myo1d.cli prepare-receptor \
  --run-id RUN \
  --state EGFR_160-185|EGFR_170-200|3GT8_raw \
  --source PATH/to/receptor.pdb \
  [--profile codex_dev|hpc_strict] \
  [--mode smoke_env|smoke_input] \
  [--strict]
```

Behavior:
- `--state`: one of the three known states; determines role mapping (primary vs reference_control)
- `--source`: receptor PDB path (real or fixture)
- Process exit: 0 PASS/WARN, 1 FAIL
- Stdout: short summary with state, role, observed_chains, normalized_paths, qc_status

## 5. Files to create / modify

Create:

```text
fresh/src/egfr_myo1d/model/__init__.py
fresh/src/egfr_myo1d/model/receptor_normalize.py
fresh/src/egfr_myo1d/model/receptor_qc.py
fresh/src/egfr_myo1d/io/residue_mapping.py
fresh/tests/test_m1_phase4_receptor_normalization.py
fresh/tests/fixtures/m1_phase4_receptor/explicit_AB_dimer.pdb
fresh/tests/fixtures/m1_phase4_receptor/duplicate_chain_X_dimer.pdb
fresh/tests/fixtures/m1_phase4_receptor/single_chain_monomer.pdb
fresh/tests/fixtures/m1_phase4_receptor/v924r_warn.pdb
fresh/tests/fixtures/m1_phase4_receptor/atom_hetatm_mixed.pdb
fresh/docs/m1_phase4_receptor_normalization.md
fresh/docs/m1_phase4_changes.md
```

(Reuse `fresh/tests/fixtures/{mini_explicit_AB.pdb, mini_duplicate_chain_X.pdb, task3_inputs/egfr_*.pdb}` where appropriate.)

Modify:

```text
fresh/src/egfr_myo1d/cli.py   # add prepare-receptor subparser + handler
```

## 6. Public API

`model/receptor_normalize.py`:

```python
def normalize_receptor(ctx, source_pdb, state_id, profile, strict=False):
    # type: (RunContext, Path, str, str, bool) -> NormalizedReceptor
    """
    Parse receptor PDB.
    Detect Case A/B/C.
    Split into protomer A and B if Case B.
    Refuse to promote Case C to dimer.
    Apply 669-1014 dockable crop.
    Apply +1000 runtime offset to protomer B.
    Write normalized PDBs and mapping CSV under ctx.run_dir.
    Append phase status. Return NormalizedReceptor report.
    """
```

`NormalizedReceptor` dataclass:

```text
state_id: str
role: "primary_membrane_validated_state" | "crystallographic_reference_control_not_primary_membrane_state"
source_file: str
case: "A_explicit_AB" | "B_duplicate_chain" | "C_monomer"
observed_chains: list[str]
protomer_count: int
warnings: list[dict]            # severity, classification, code, message
full_frame_pdb: Path | None     # normalized full frame, explicit A/B
dockable_pdb: Path | None       # 669-1014 cropped
runtime_offset_pdb: Path | None # protomer B residues + 1000
mapping_csv: Path
audit_csv: Path
manifest_json: Path
v924r_warn: bool                # true if ARG924 (or equivalent) detected
v924r_handled: "warn_only_not_mutated"
status: "PASS" | "WARN" | "FAIL"
```

`io/residue_mapping.py`:

```python
def write_residue_mapping(path, rows, ctx):
    # rows: list of dicts with keys per spec column

def read_residue_mapping(path):
    # returns list of dicts

MAPPING_CSV_COLUMNS = [
    "state", "source_file", "protomer_id",
    "source_chain", "source_resseq", "source_icode", "source_resname",
    "runtime_chain", "runtime_resseq", "atom_count", "role"
]
```

`model/receptor_qc.py`:

```python
def audit_receptor_residues(structure, state_id, dockable_crop, expected_chains, warn_mutations):
    # returns list of audit rows: residue_number, chain, in_dockable_crop, missing, warn_mutation_hit, ...
```

## 7. Required output files (handoff §14.2)

After `prepare-receptor` runs:

```text
# For EGFR_160-185 / EGFR_170-200:
fresh/runs/<run_id>/normalized/receptors/<state>_full_frame_explicit_AB.pdb
fresh/runs/<run_id>/normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb
fresh/runs/<run_id>/normalized/receptors/<state>_runtime_offset_receptor_only.pdb
fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv
fresh/runs/<run_id>/qc/<state>_receptor_normalization_audit.csv
fresh/runs/<run_id>/manifest/<state>_receptor_manifest.json

# For 3GT8_raw:
fresh/runs/<run_id>/normalized/receptors/3GT8_raw_explicit_AB.pdb
fresh/runs/<run_id>/normalized/receptors/3GT8_raw_runtime_offset_receptor_only.pdb
fresh/runs/<run_id>/qc/3GT8_raw_receptor_mapping.csv
fresh/runs/<run_id>/manifest/3GT8_raw_receptor_manifest.json
```

(3GT8_raw is a kinase-domain-only crystal; the dockable crop output may be omitted or mark `requires_membrane_gate: false` per receptor_states.yaml.)

Mapping CSV columns:

```csv
state,source_file,protomer_id,source_chain,source_resseq,source_icode,source_resname,runtime_chain,runtime_resseq,atom_count,role
```

Audit CSV columns (suggested, mirror Task 4 audit shape):

```csv
state,fixture_role,chain_id,protomer_id,residue_number,insertion_code,residue_name,record_type,biological,is_cap,is_warn_mutation,warning_classification,production_policy,notes
```

Manifest JSON includes: state, role, source_sha256, output_pdb_sha256s, mapping_csv_sha256, normalization_case, warnings, status.

## 8. Behavior policy

```text
- Original residue numbers preserved on protomer A.
- Protomer B receives runtime residue number = source_resseq + 1000.
- Insertion codes preserved.
- ATOM and HETATM records both preserved.
- Standard AAs written as HETATM are kept as biological residues.
- ACE/NME or other cap HETATMs preserved (passed through to dockable PDB).
- Lipid HETATM (POPC, POPS, etc.): NOT written into dockable_*.pdb (receptor-only).
  Full-frame PDB may include them if requested by future flag (this phase: drop them from dockable, retain in full_frame if input has them).
- Water/ion (HOH, NA, CL): dropped from dockable receptor.
- Ligand HETATMs in receptor input: dropped from dockable receptor (this phase does not write ligand-receptor complexes).
- V924R (ARG at residue 924 where VAL is expected): WARN only, never mutated. v924r_warn=true in report.
- Case B detection: duplicate atom identifiers OR residue-number reset within same chain.
  Splitting strategy: first half of residues → protomer A, second half → protomer B (or use blank line / TER record if present).
- Case C detection: single chain, no duplicate identifiers, no residue reset, residue count consistent with monomer.
  Status: WARN in codex_dev (recorded as not-valid-for-dimer-analysis), FAIL in hpc_strict.
- Missing source: WARN in codex_dev, FAIL in hpc_strict.
- 3GT8_raw: dockable crop optional; mark role as reference_control_not_primary_membrane_state in manifest.
```

## 9. Severity rules

```text
PASS:  Case A, all expected chains present, V924R absent, all required residues in dockable crop, no resets
WARN:  Case A with V924R; Case B successfully split; missing residues in 669-1014 (recorded); HETATM caps present; codex_dev with missing source
FAIL:  Case C in hpc_strict; malformed PDB; unwriteable output; path-traversal run_id; explicit chain split impossible (e.g., 3 chains); hpc_strict with missing source
```

## 10. Tests required (≥12)

```text
test_explicit_AB_passthrough_preserves_chains_and_numbers
test_duplicate_chain_X_detected_and_split_into_A_B
test_monomer_only_input_warns_in_codex_dev_fails_in_hpc_strict
test_dockable_crop_excludes_634_668_keeps_669_1014
test_dockable_crop_records_missing_ranges_in_audit
test_runtime_offset_protomer_B_plus_1000_only
test_runtime_offset_protomer_A_unchanged
test_mapping_csv_round_trip_residue_identity
test_mapping_csv_columns_match_spec
test_3gt8_raw_marked_reference_control_not_primary
test_v924r_residue_warned_not_mutated
test_atom_and_hetatm_records_preserved
test_lipid_hetatm_dropped_from_dockable
test_ace_nme_caps_preserved_in_dockable
test_normalization_writes_under_run_dir_only
test_path_traversal_run_id_rejected
test_cli_help_includes_prepare_receptor
```

(17 tests above; ≥12 required.)

## 11. Acceptance commands

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"

pytest -q fresh/tests/test_m1_phase4_receptor_normalization.py
pytest -q fresh/tests

# Explicit AB Case A
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase4_local
python -m egfr_myo1d.cli prepare-receptor --run-id m1_phase4_local --state EGFR_160-185 --source fresh/tests/fixtures/m1_phase4_receptor/explicit_AB_dimer.pdb --profile codex_dev

# Duplicate chain X Case B
python -m egfr_myo1d.cli prepare-receptor --run-id m1_phase4_local --state EGFR_170-200 --source fresh/tests/fixtures/m1_phase4_receptor/duplicate_chain_X_dimer.pdb --profile codex_dev

# Monomer Case C in hpc_strict
python -m egfr_myo1d.cli prepare-receptor --run-id m1_phase4_local --state EGFR_160-185 --source fresh/tests/fixtures/m1_phase4_receptor/single_chain_monomer.pdb --profile hpc_strict || echo "Expected FAIL"

# 3GT8_raw reference control
python -m egfr_myo1d.cli prepare-receptor --run-id m1_phase4_local --state 3GT8_raw --source fresh/tests/fixtures/m1_phase4_receptor/explicit_AB_dimer.pdb --profile codex_dev

# Path traversal
python -m egfr_myo1d.cli prepare-receptor --run-id ../bad_run --state EGFR_160-185 --source fresh/tests/fixtures/m1_phase4_receptor/explicit_AB_dimer.pdb

# Old workflow protection
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

## 12. Final response format

```text
M1 Phase 4 status: PASS / PASS WITH WARNINGS / FAIL
Files created
Files modified
Commands run and results
Test summary (prior + Phase 1-3 + Phase 4 new)
Output verification per state:
- explicit_AB output path exists and chains preserved
- duplicate_chain_X split into A and B
- runtime offset: protomer A unchanged, protomer B +1000
- mapping CSV round-trips
- 3GT8_raw marked reference_control
V924R behavior: WARN, not mutated
Lipid HETATM behavior: dropped from dockable
Cap HETATM behavior: preserved
Acceptance closure: M1 §23 #10 and #11 closed
Old workflow protection: empty diff
Known limitations: no membrane frame yet (Phase 5)
```
