# M2.1 PPI Input Generation

M2.1 consumes canonical M1 outputs and emits staged EGFR-MYO1D PPI input packs for later PyRosetta and LightDock phases.

This phase does not run docking, relaxation, LightDock, Vina, fpocket, P2Rank, compound scoring, candidate nomination, qsub/PBS submission, or cleanup deletion.

## Inputs

Required per primary membrane-validated EGFR state:

```text
fresh/runs/<run_id>/normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb
fresh/runs/<run_id>/normalized/receptors/<state>_runtime_offset_receptor_only.pdb
fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv
fresh/runs/<run_id>/manifest/membrane_frame.json
fresh/runs/<run_id>/normalized/myo1d/MYO1D_955_1006.pdb
fresh/runs/<run_id>/qc/myo1d_construct_qc.csv
```

The default state set is read from `fresh/configs/receptor_states.yaml` and excludes `3GT8_raw`, which remains a reference/control state.

## Outputs

```text
fresh/runs/<run_id>/prepared/m2_1_ppi_inputs/<state>/receptor/
fresh/runs/<run_id>/prepared/m2_1_ppi_inputs/<state>/myo1d/
fresh/runs/<run_id>/prepared/m2_1_ppi_inputs/<state>/specs/
fresh/runs/<run_id>/qc/m2_1_ppi_input_qc.csv
fresh/runs/<run_id>/manifest/m2_1_ppi_input_manifest.json
fresh/runs/<run_id>/reports/m2_1_ppi_input_generation.md
```

The PyRosetta spec references the M1 runtime-offset receptor plus mapping CSV. The LightDock spec references the M1 dockable explicit A/B receptor plus mapping CSV. Both specs set `execution_allowed=false`, `docking_executed=false`, and `engine_command=null`.

## MYO1D Policy

The M1 `MYO1D_955_1006.pdb` artifact is consumed as the comparator/noise-monitor source. M2.1 derives `MYO1D_955_1001_primary.pdb` as the production-oriented PPI partner while preserving MYO1D residue numbering. The 955-1006 comparator remains marked `noise_monitor_only`.

Sheet 8/9 active-face and sheet 12 support residues are carried as QC evidence only. `score_bonus_allowed=false` and `key_residue_bonus_weight=0.0` are enforced in every spec.

## Guardrails

- EGFR chain/protomer identity remains explicit A/B.
- EGFR source residue identity is preserved through the M1 mapping CSV.
- Protomer B runtime offset must remain `source_resseq + 1000`.
- V924R-like artifacts are never repaired or silently treated as WT.
- ATP-site overlap remains blocking/non-target for later PPI-disruptive pocket objectives.
- Membrane-proximal regions remain flag/block material for later pose QC.
- Confidential ligand structures and internal compound IDs are not consumed by M2.1.

## CLI

```bash
python -m egfr_myo1d.cli generate-m2-ppi-inputs \
  --run-id <run_id> \
  --mode smoke_input \
  --profile codex_dev
```

Use `--states EGFR_160-185,EGFR_170-200` to limit state generation explicitly.

