# Task 5 Real Input Readiness

Task 5 adds a strict readiness bridge between synthetic fixtures and real uploaded EGFR/MYO1D/membrane-frame inputs. It consumes the Task 4 `ppi_input_contract.json` style and writes auditable PASS/WARN/FAIL reports without preparing docking jobs or modifying scientific structures.

Task 5 does not run docking, PyRosetta, LightDock, Vina, fpocket, P2Rank, scoring, candidate nomination, qsub/PBS generation, cleanup deletion, mutation repair, residue renumbering, or silent receptor normalization.

## Command Examples

Bash:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.cli validate-real-inputs \
  --run-id real_readiness_001 \
  --mode smoke_input \
  --profile hpc_strict \
  --input-root fresh/data/raw \
  --contract fresh/data/raw/ppi_input_contract.json
```

PowerShell:

```powershell
$env:PYTHONPATH = "$PWD\fresh\src;$env:PYTHONPATH"
python -m egfr_myo1d.cli validate-real-inputs `
  --run-id real_readiness_001 `
  --mode smoke_input `
  --profile hpc_strict `
  --input-root fresh/data/raw `
  --contract fresh/data/raw/ppi_input_contract.json
```

## Required Real Input Layout

Recommended public-safe layout:

```text
fresh/data/raw/
├── ppi_input_contract.json
├── receptors/
│   └── <real_egfr_dimer>.pdb
├── myo1d/
│   └── <real_myo1d_model>.pdb
└── membrane_frame.json
```

Confidential ligand structures and private compound mappings should remain out of git-tracked public files. Task 5 does not require ligand structures.

## Contract Shape

The real contract follows Task 4 fields:

```json
{
  "task3_closeout": {
    "pass_with_warnings_allowed_only_for_explicit_fixtures": true,
    "production_same_warning_policy": "fail_or_quarantine",
    "fixture_warning_whitelist": []
  },
  "receptor": {
    "id": "EGFR_real_dimer",
    "path": "receptors/EGFR_real_dimer.pdb",
    "expected_chains": ["A", "B"],
    "required_protomer_count": 2,
    "biological_residue_range": [699, 1007],
    "warn_mutations": [
      {"residue_number": 924, "expected": "VAL", "warn_if": "ARG", "reason": "3GT8 V924R-like artifact"}
    ],
    "non_target_masks": {
      "atp_site_residues": [745, 855],
      "membrane_proximal_label": "computed_from_membrane_frame"
    }
  },
  "myo1d": {
    "primary_construct": {
      "id": "MYO1D_sheet8_9_12_core_955_1001",
      "path": "myo1d/AF-O94832-F1-model_v6.pdb",
      "biological_residue_range": [955, 1001]
    },
    "comparator_construct": {
      "id": "MYO1D_ext_beta_meander_955_1006_tail_masked",
      "path": "myo1d/AF-O94832-F1-model_v6.pdb",
      "biological_residue_range": [955, 1006]
    },
    "active_face": [961, 962, 963, 964, 968, 969, 970, 971, 972],
    "sheet12_support": [993, 994, 995, 996, 997],
    "structural_buffer": [955, 956, 957, 958, 959, 960],
    "contact_monitoring_cap": [998, 999, 1000, 1001, 1002, 1003, 1004, 1005, 1006],
    "key_residue_contact_is_annotation_only": true
  },
  "membrane_frame": {
    "path": "membrane_frame.json"
  }
}
```

## Warning vs Blocker Policy

In `hpc_strict` or strict real-input mode:

- ambiguous single-chain receptor dimers fail or quarantine
- missing explicit A/B protomer chains fail or quarantine
- duplicate atom identities fail or quarantine
- V924R-like ARG924 is reported and not repaired silently
- MYO1D 962-start terminal artifacts fail or quarantine
- MYO1D 955-1006 is comparator/noise-monitor only, not primary
- active-face residues are annotation-only and `score_bonus_allowed` remains `false`

In `codex_dev`, explicitly allowlisted fixture warnings may remain warnings for regression fixtures. Real production inputs should not rely on fixture allowlists.

## Outputs

Task 5 writes under `fresh/runs/<run_id>/`:

```text
manifest/real_input_manifest.json
manifest/real_input_readiness_report.json
qc/real_egfr_readiness_audit.csv
qc/real_myo1d_readiness_audit.csv
qc/real_membrane_frame_audit.csv
qc/real_input_blockers.csv
```

The report includes input paths, SHA256 values, PASS/WARN/FAIL counts, fixture-warning count, production-blocker count, quarantine count, final verdict, blocker list, future mask serialization previews, and a `not_implemented` list.

## Relationship To Task 4

Task 5 consumes Task 4 policy concepts but does not emit prepared PDB packs. It validates whether real inputs are ready to be passed into a future preparation or sampling step. It keeps Task 4 warning classes explicit so synthetic regression warnings cannot become silent production acceptance.

## Intentionally Not Implemented

- PPI docking
- PyRosetta docking or relaxation
- LightDock
- Vina
- fpocket/P2Rank pocket discovery
- compound docking
- scoring or candidate nomination
- qsub/PBS generation
- cleanup deletion
- EGFR mutation repair
- receptor remodeling
- residue renumbering
- silent receptor normalization


---

## M1 dependency (added by Phase 9 alignment)

This task module was originally designed to run before the M1 foundation
modules (cleanup, receptor normalize, membrane frame, MYO1D construct, ligand
manifest) existed. After M1 Phases 1-8 the canonical input artifacts now live
under `fresh/runs/<run_id>/normalized/`, `fresh/runs/<run_id>/manifest/membrane_frame.json`,
and `fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv`.

Phase 9 takes an **additive** approach: this module's logic and output paths
are unchanged. The alignment between M1 outputs and this task's consumption
points is recorded by `validation/m1_alignment.py`'s `record_m1_alignment(ctx)`
helper, which writes `manifest/m1_alignment.json` listing which M1 artifacts
are present and which Task 4-9 modules would naturally consume them in a future
M2 actual-execution phase.

See `fresh/docs/m1_phase9_tasks4to9_realignment.md` for the full rationale.
