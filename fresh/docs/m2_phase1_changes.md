# M2 Phase 1 Changes

Implements M2.1 PPI input generation as a fresh-only staging layer.

## Files Created

```text
fresh/src/egfr_myo1d/m2/__init__.py
fresh/src/egfr_myo1d/m2/ppi_inputs.py
fresh/tests/test_m2_phase1_ppi_input_generation.py
fresh/docs/m2_phase1_ppi_input_generation.md
fresh/docs/m2_phase1_changes.md
```

## Files Modified

```text
fresh/src/egfr_myo1d/cli.py
```

## CLI Addition

```bash
python -m egfr_myo1d.cli generate-m2-ppi-inputs --run-id RUN [--states <comma-list>]
```

## Behavior

- Consumes M1 normalized EGFR/MYO1D/membrane artifacts from `fresh/runs/<run_id>/`.
- Emits input packs/spec JSON under `fresh/runs/<run_id>/prepared/m2_1_ppi_inputs/`.
- Derives MYO1D 955-1001 as the primary production-oriented partner from M1 `MYO1D_955_1006.pdb`.
- Keeps MYO1D 955-1006 as comparator/noise-monitor only.
- Writes `manifest/m2_1_ppi_input_manifest.json`, `qc/m2_1_ppi_input_qc.csv`, and `reports/m2_1_ppi_input_generation.md`.
- Records `execution_allowed=false`, `docking_executed=false`, `score_bonus_allowed=false`, and `ligand_or_compound_inputs_used=false`.

## Not Implemented By Design

- PyRosetta docking or relaxation
- LightDock execution
- Vina, fpocket, or P2Rank runtime
- Compound docking, scoring, or candidate nomination
- qsub/PBS/sbatch submission
- Cleanup deletion
- Receptor mutation repair or renumbering

