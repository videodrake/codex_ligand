# Task 3 Structural Input Contracts

Task 3 adds lightweight structural input contracts and QC reporting for the fresh EGFR-MYO1D workflow. It validates whether receptor, MYO1D, and membrane-frame inputs are structurally and metadata-ready enough for later tasks.

It does not validate biological binding and does not run docking, pocket discovery, receptor normalization, MYO1D production slicing, qsub/PBS, cleanup deletion, scoring, or candidate nomination.

## CLI Smoke

Bash:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id test_task3_local
python -m egfr_myo1d.cli validate-structures --run-id test_task3_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task3_inputs
python -m egfr_myo1d.cli status --run-id test_task3_local
```

PowerShell:

```powershell
$env:PYTHONPATH = "$PWD\fresh\src;$env:PYTHONPATH"
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id test_task3_local
python -m egfr_myo1d.cli validate-structures --run-id test_task3_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task3_inputs
python -m egfr_myo1d.cli status --run-id test_task3_local
```

## Outputs

Task 3 writes only under `fresh/runs/<run_id>/`:

```text
manifest/structure_input_manifest.json
manifest/structure_qc_report.json
qc/chain_protomer_audit.csv
qc/residue_mapping_audit.csv
qc/membrane_frame_report.json
qc/myo1d_annotation_report.json
logs/master.log
logs/phase_status.jsonl
logs/job_status.jsonl
```

## EGFR Metadata

EGFR receptor contracts describe expected dimer/protomer metadata:

```json
{
  "receptor_id": "EGFR_valid_dimer_AB",
  "path": "egfr_valid_dimer_AB.pdb",
  "expected_assembly": "dimer",
  "expected_protomers": [
    {"protomer_id": "A", "chain_ids": ["A"], "required_ranges": [[699, 700]]},
    {"protomer_id": "B", "chain_ids": ["B"], "required_ranges": [[699, 700]]}
  ],
  "allow_single_chain_dimer_with_mapping": true
}
```

The validator preserves chain IDs, residue numbers, insertion codes, ATOM/HETATM records, and residue identities. It reports, but does not fix:

- single-chain dimer ambiguity
- duplicate atom identities from same-chain duplicate dimers
- missing required residue ranges
- residue-number reset risk
- large residue-number jumps
- HETATM caps or modified residues
- 3GT8 V924R-like mutation warnings

V924R is reported as a warning so it is not silently treated as normal WT evidence. Task 3 does not mutate or normalize the receptor.

## MYO1D Annotation

Task 3 records the MYO1D beta-meander annotation:

```text
primary_construct = MYO1D_sheet8_9_12_core_955_1001
primary_range = 955-1001
comparator_construct = MYO1D_ext_beta_meander_955_1006_tail_masked
comparator_range = 955-1006
structural_buffer = 955-960
primary_active_face = 961-964,968-972
support_region = 993-997
short_c_terminal_cap = 998-1001
extended_tail_noise_zone = 998-1006
key_residue_bonus_weight = 0.0
key_residue_contact_is_annotation_only = true
```

Residues 961-964 and 968-972 are required active-face annotations. Residues 993-997 are sheet-12 support annotations. Residues 998-1006 are tracked as terminal/tail monitoring regions when present because they can dominate contact accounting without necessarily representing the intended beta-meander face.

The key-residue annotation is QC-only. Task 3 assigns no score bonus.

## Membrane Frame

Membrane-frame JSON validation checks that `membrane_normal` and `dimer_axis` are numeric length-3 nonzero vectors. The QC report stores normalized vectors but does not overwrite user input or infer a final biological frame from fixtures.

Warnings are produced when:

- membrane normal and dimer axis are nearly parallel
- protomer centroids are missing, invalid, or identical

Zero-norm vectors are failures.

## Expected codex_dev Warnings

`codex_dev` / `smoke_env` is intentionally forgiving about future scientific inputs. Expected warnings include:

- missing real production structures
- single-chain dimer ambiguity
- V924R-like mutation fixture
- MYO1D extended tail/noise-zone fixtures
- missing membrane-frame metadata when no contract fixture is supplied

Malformed files and unsafe paths still fail.

