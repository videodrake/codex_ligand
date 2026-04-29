# Claude M1 Phase 3 Prompt — MYO1D Construct Slicing + QC v0.1

Branch `claude/task10`. Phases 1-2 complete (cleanup + relocation). This is **M1 Phase 3** — implements MYO1D construct preparation per `milestone1_foundation_codex_handoff_v0_5.md` §16, closing M1 §23 #13.

## 1. Project context

Phase 2 relocated existing MYO1D logic to `fresh/src/egfr_myo1d/myo1d/{construct.py, pdb_writer.py}`. Those files currently contain Task 4-flavored helpers (residue annotation, terminal-artifact detection, validate_active_face_presence). They do NOT yet emit a canonical M1 MYO1D construct PDB or write the M1 QC CSV.

Phase 3 extends `myo1d/construct.py` with M1-spec emit functions and adds a new `myo1d/qc.py` for the QC CSV. A new CLI subcommand `prepare-myo1d` orchestrates them.

## 2. Absolute rules

Do not modify the old workflow:

```text
run_production.py
main.py
egfr_pipeline/**
config/**
docs/runbook.md
output/**
results_export/**
```

Maintain Python 2.7.11 / 3.9 syntax compatibility.

Source-of-truth values must come from `fresh/configs/gates.yaml` (or `fresh_run.yaml`):

```text
construct: "955-1006"
key_residues:
  sheet8:  "961-964"
  sheet9:  "968-972"
  sheet12: "993-997"
watch_residues:
  n_terminal: "955-957"
  c_terminal: "1001-1006"
key_residue_bonus_weight: 0.0
```

No hardcoding of these values inside the module.

## 3. Scope

In scope:
- Extend `myo1d/construct.py` with: `slice_myo1d_construct(structure, start, end)`, `emit_myo1d_construct_pdb(ctx, structure, output_path, residue_range)`
- Create `myo1d/qc.py` with: `run_myo1d_qc(ctx, source_pdb, construct_range, key_residues, watch_residues, profile) -> Myo1dQcReport`
- Add `prepare-myo1d` CLI subcommand to `cli.py`
- Tests under `fresh/tests/test_m1_phase3_myo1d_construct_qc.py` (≥10 tests)
- Fixtures under `fresh/tests/fixtures/m1_phase3_myo1d/` (slicing variants)
- Docs `fresh/docs/m1_phase3_myo1d_construct_qc.md` and `fresh/docs/m1_phase3_changes.md`

Out of scope:
- Modifying existing `myo1d/construct.py` Task 4 functions (they continue serving validation/prepared_inputs.py)
- Receptor normalization (Phase 4)
- Membrane frame computation (Phase 5)
- Any docking, scoring, ligand work
- Phase 9 schema realignment of Task 4 to consume M1 outputs

## 4. Required CLI behavior

```bash
python -m egfr_myo1d.cli prepare-myo1d \
  --run-id RUN \
  --source PATH/to/myo1d.pdb \
  [--construct 955-1006] \
  [--profile codex_dev|hpc_strict] \
  [--mode smoke_env|smoke_input]
```

Behavior:
- `--source`: path to source MYO1D PDB (e.g., `fresh/data/raw/myo1d/AF-O94832-F1-model_v6.pdb` or fixture)
- `--construct`: residue range string `start-end`. Default `955-1006` from gates.yaml. Range `962-1006` triggers terminal-artifact warning (or FAIL in hpc_strict).
- `--profile`: codex_dev tolerates missing source (WARN); hpc_strict treats missing source as FAIL.
- `--mode`: same as other CLI commands.
- Process exit: 0 PASS/WARN, 1 FAIL.

## 5. Files to create / modify

Create:

```text
fresh/src/egfr_myo1d/myo1d/qc.py
fresh/tests/test_m1_phase3_myo1d_construct_qc.py
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_955_1006_valid.pdb       # primary construct fixture
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_962_1006_terminal_bad.pdb # negative regression
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_955_1001_short.pdb        # comparator construct
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_with_ace_nme_caps.pdb     # cap preservation test
fresh/docs/m1_phase3_myo1d_construct_qc.md
fresh/docs/m1_phase3_changes.md
```

(Reuse `fresh/tests/fixtures/task3_inputs/myo1d_*.pdb` where appropriate; only add new fixtures for cases not already present.)

Modify:

```text
fresh/src/egfr_myo1d/myo1d/construct.py    # extend with slice_*, emit_* functions
fresh/src/egfr_myo1d/cli.py                 # add prepare-myo1d subparser + handler
```

Optionally modify:

```text
fresh/configs/gates.yaml                    # only if a key value is missing; do NOT change existing values
```

## 6. Public API

`myo1d/construct.py` (added functions):

```python
def slice_myo1d_construct(structure, start, end):
    # type: (PDBStructure, int, int) -> PDBStructure
    """Return a new PDBStructure containing residues whose source_resseq is in [start, end].
    Preserves chain_id, residue numbering, insertion_code, ATOM and HETATM records.
    Does NOT renumber.
    """

def emit_myo1d_construct_pdb(ctx, structure, output_path, residue_range):
    # type: (RunContext, PDBStructure, Path, tuple) -> None
    """Write structure to output_path under ctx.run_dir.
    Refuses to write outside run_dir (uses ctx.require_within_run_dir).
    Preserves ACE/NME/cap HETATM records.
    """
```

`myo1d/qc.py` (new module):

```python
def run_myo1d_qc(ctx, source_pdb, construct_range, key_residues, watch_residues, profile):
    # type: (RunContext, Path, tuple, dict, dict, str) -> Myo1dQcReport
    """Parse source MYO1D PDB, slice to construct_range, write normalized PDB,
    write QC CSV, append phase status, return report."""
```

`Myo1dQcReport` dataclass:

```text
construct_id: str
source_file: str
chain_id: str
start_residue: int
end_residue: int
n_residues: int
missing_residues: list[int]
key_sheet8_present: bool
key_sheet9_present: bool
key_sheet12_present: bool
n_watch_present: bool
c_watch_present: bool
ace_nme_caps_present: bool
status: "PASS" | "WARN" | "FAIL"
warnings: list[str]
output_pdb: Path | None
output_qc_csv: Path
```

## 7. Required output files

After `prepare-myo1d` runs:

```text
fresh/runs/<run_id>/normalized/myo1d/MYO1D_955_1006.pdb     # the normalized construct PDB
fresh/runs/<run_id>/qc/myo1d_construct_qc.csv               # one row per QC criterion or one combined row
fresh/runs/<run_id>/manifest/myo1d_construct_manifest.json  # metadata + sha256 of source + sha256 of output
fresh/runs/<run_id>/logs/phase_status.jsonl                 # appended
fresh/runs/<run_id>/logs/master.log                         # appended
```

QC CSV columns (handoff §16):

```csv
construct_id,source_file,chain_id,start_residue,end_residue,n_residues,missing_residues,key_sheet8_present,key_sheet9_present,key_sheet12_present,n_watch_present,c_watch_present,ace_nme_caps_present,status,warnings
```

`missing_residues` and `warnings` may be `;`-separated strings inside a single CSV cell, or use a sub-CSV column convention consistent with existing audits.

## 8. Behavior policy

```text
- Source PDB residue numbering is preserved (no renumber).
- Construct range default 955-1006; values come from gates.yaml.
- Construct 962-1006 produces WARN in codex_dev (terminal artifact), FAIL in hpc_strict.
- Construct 955-1001 is allowed as comparator/short variant; produces PASS with `c_watch_present=false`.
- Active-face residues 961-964 (sheet8) + 968-972 (sheet9): annotation/QC only. No score bonus.
- sheet12 support residues 993-997: annotation/QC only.
- N-terminal watch 955-957, C-terminal watch 1001-1006: counted but not used to promote/reject.
- ACE/NME/other cap HETATM: preserved in output PDB; reported in `ace_nme_caps_present`.
- Standard AA written as HETATM (e.g., ILE1000 in capped variant): kept as biological residue, not dropped.
- Source missing in codex_dev: WARN with status="WARN", missing_required_inputs recorded in manifest.
- Source missing in hpc_strict: FAIL.
```

## 9. Severity rules

```text
PASS:  source parsed, construct sliced, all required key residues present, no terminal artifact
WARN:  source missing in codex_dev, comparator (955-1001) accepted, capped HETATM present, watch residues partially present
FAIL:  source missing in hpc_strict, malformed source, residue range invalid (e.g., end<start), terminal artifact in hpc_strict, write attempt outside run_dir
```

## 10. Tests required (≥10)

```text
test_myo1d_955_1006_construct_emitted_with_original_numbering
test_myo1d_key_residues_present_assertions
test_myo1d_watch_residues_recorded_not_promoted
test_myo1d_962_start_terminal_artifact_warned_in_codex_dev
test_myo1d_962_start_terminal_artifact_fails_in_hpc_strict
test_myo1d_ace_nme_cap_preserved_in_output_pdb
test_myo1d_qc_csv_columns_match_spec
test_myo1d_missing_source_reports_cleanly_in_codex_dev
test_myo1d_score_bonus_zero_enforced_via_gates_yaml
test_cli_help_includes_prepare_myo1d
test_prepare_myo1d_writes_under_run_dir_only
test_myo1d_955_1001_comparator_pass_with_c_watch_partial
```

(13 tests above; ≥10 required.)

Test conventions: match `test_task9_*.py` style; local helpers; pytest tmp_path; one CLI help via subprocess.

## 11. Acceptance commands

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"

# Targeted phase tests
pytest -q fresh/tests/test_m1_phase3_myo1d_construct_qc.py

# Full suite (must include prior 98 + Phase 1 (8) + Phase 3 new)
pytest -q fresh/tests

# CLI smoke (codex_dev with valid fixture)
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase3_local
python -m egfr_myo1d.cli prepare-myo1d --run-id m1_phase3_local --source fresh/tests/fixtures/m1_phase3_myo1d/myo1d_955_1006_valid.pdb --construct 955-1006 --profile codex_dev

# CLI smoke (negative regression with terminal artifact, expect WARN exit 0 in codex_dev)
python -m egfr_myo1d.cli prepare-myo1d --run-id m1_phase3_local --source fresh/tests/fixtures/m1_phase3_myo1d/myo1d_962_1006_terminal_bad.pdb --construct 962-1006 --profile codex_dev

# Same fixture with hpc_strict, expect FAIL exit 1
python -m egfr_myo1d.cli prepare-myo1d --run-id m1_phase3_local --source fresh/tests/fixtures/m1_phase3_myo1d/myo1d_962_1006_terminal_bad.pdb --construct 962-1006 --profile hpc_strict || echo "Expected nonzero exit"

# Path traversal
python -m egfr_myo1d.cli prepare-myo1d --run-id ../bad_run --source fresh/tests/fixtures/m1_phase3_myo1d/myo1d_955_1006_valid.pdb --construct 955-1006

# Old workflow protection
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

## 12. Final response format

```text
M1 Phase 3 status: PASS / PASS WITH WARNINGS / FAIL
Files created
Files modified
Commands run and results
Test summary (prior + Phase 1 + Phase 3 new = total)
Output artifact verification:
- normalized/myo1d/MYO1D_955_1006.pdb exists with original numbering
- qc/myo1d_construct_qc.csv columns match §7
- manifest/myo1d_construct_manifest.json includes source + output sha256
Negative-regression behavior:
- 962-1006 in codex_dev = WARN; in hpc_strict = FAIL
- 955-1001 in either profile = PASS with c_watch partial
Cap preservation:
- ACE/NME records preserved in output PDB
Score bonus:
- key_residue_bonus_weight=0.0 enforced from gates.yaml
Old workflow protection: empty diff
Acceptance closure: M1 §23 #13 closed
Known limitations / not implemented by design:
- No receptor work (Phase 4)
- No membrane frame (Phase 5)
- No M2 work
```
