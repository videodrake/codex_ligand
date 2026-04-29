# M1 Phase 3 — MYO1D Construct + QC

Closes M1 §23 #13 (MYO1D 955-1006 construct QC) per `milestone1_foundation_codex_handoff_v0_5.md` §16.

## What it does

1. Reads source-of-truth MYO1D values from `fresh/configs/gates.yaml`:
   - `myo1d.construct: "955-1006"` (default range)
   - `myo1d.key_residues.{sheet8,sheet9,sheet12}` (annotation only)
   - `myo1d.watch_residues.{n_terminal,c_terminal}` (annotation only)
   - `myo1d.key_residue_bonus_weight: 0.0` (enforced; `RuntimeError` if non-zero)
2. Slices the source MYO1D PDB to the requested residue range, preserving original residue numbering, chain identity, insertion codes, and ACE/NME caps.
3. Emits the canonical M1 normalized PDB at `fresh/runs/<run_id>/normalized/myo1d/MYO1D_<start>_<end>.pdb`.
4. Computes the QC fields per handoff §16.3 and writes:
   - `fresh/runs/<run_id>/qc/myo1d_construct_qc.csv`
   - `fresh/runs/<run_id>/manifest/myo1d_construct_manifest.json`
5. Appends a `prepare-myo1d` phase status entry.

## CLI

```bash
python -m egfr_myo1d.cli prepare-myo1d \
    --run-id RUN \
    --source PATH/to/myo1d.pdb \
    [--construct 955-1006] \
    [--profile codex_dev|hpc_strict] \
    [--mode smoke_env|smoke_input]
```

`--construct` defaults to the value in `gates.yaml` (`955-1006`).

## Module additions

### `myo1d/construct.py` (extended in Phase 3)

```python
def slice_myo1d_construct(structure, start, end, include_caps=True): ...
def emit_myo1d_construct_pdb(ctx, structure, output_path): ...
```

### `myo1d/qc.py` (new in Phase 3)

```python
def run_myo1d_qc(ctx, source_pdb, construct_range=None, profile="codex_dev") -> Myo1dQcReport: ...
def parse_construct_range(text_or_tuple) -> (start, end): ...
def expand_residue_set(text) -> [int]: ...
def load_myo1d_gate(ctx) -> dict: ...

MYO1D_CONSTRUCT_QC_COLUMNS = [
    "construct_id", "source_file", "chain_id",
    "start_residue", "end_residue", "n_residues", "missing_residues",
    "key_sheet8_present", "key_sheet9_present", "key_sheet12_present",
    "n_watch_present", "c_watch_present", "ace_nme_caps_present",
    "status", "warnings",
]
```

## Severity

| Status | Conditions |
| --- | --- |
| `PASS` | source parsed, all sheet 8/9/12 key residues present, no terminal artifact, no missing key residues |
| `WARN` | source missing in `codex_dev`; `962`-start terminal artifact in `codex_dev`; capped HETATMs preserved (informational); missing residues in range; key residues incomplete |
| `FAIL` | source missing in `hpc_strict`; `962`-start in `hpc_strict`; parse error; invalid construct range; `key_residue_bonus_weight != 0` in gates.yaml; write attempt outside `run_dir` |

## Outputs

```text
fresh/runs/<run_id>/normalized/myo1d/MYO1D_<start>_<end>.pdb
fresh/runs/<run_id>/qc/myo1d_construct_qc.csv
fresh/runs/<run_id>/manifest/myo1d_construct_manifest.json
fresh/runs/<run_id>/logs/phase_status.jsonl                       (appended)
fresh/runs/<run_id>/logs/master.log                               (appended)
```

`myo1d_construct_qc.csv` has one row with all 15 columns above.

`myo1d_construct_manifest.json` records source + output sha256, residue range, presence flags, status, warnings, key_residue_bonus_weight, and `score_bonus_allowed: false`.

## Behavior policy (handoff §16.4 + §3 of phase prompt)

- Original source residue numbers preserved on emit (no renumber).
- Sheet 8/9/12 key residues are annotation/QC only. No score bonus. `key_residue_bonus_weight` is read from `gates.yaml` and asserted to equal `0.0` at every invocation.
- ACE/NME (or other caps) preserved; standard amino acids written as `HETATM` are kept as biological residues.
- `962`-start construct produces WARN in `codex_dev`; FAIL in `hpc_strict` — must not be promoted as a production partner.
- `955-1001` short construct passes when sheet 8/9/12 are present. `c_watch` (1001-1006) only partially fits in this short range; the missing 1002-1006 residues are listed in `missing_residues` but do not trigger FAIL.

## Reusable fixtures

```text
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_955_1006_valid.pdb         (29 residues 955-1006, no caps)
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_955_1001_short.pdb         (24 residues 955-1001)
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_962_1006_terminal_bad.pdb  (18 residues 962-1006, terminal artifact)
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_with_ace_nme_caps.pdb      (ACE+NME + 19 residues 955-1001)
```

## What is intentionally not in this phase

- Receptor normalization (Phase 4)
- Membrane frame computation (Phase 5)
- PBS generation (Phase 6)
- Ligand manifest (Phase 7)
- Integration orchestrator (Phase 8)
- Tasks 4-9 schema realignment (Phase 9)
- M2 docking, scoring, candidate work
