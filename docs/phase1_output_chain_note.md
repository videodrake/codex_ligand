# Phase 1 Output Chain Note

This note summarizes the current Phase 1 output chain as it exists in the
repository today.

It is meant for operators and follow-up coding work. It does not replace the
PRD or task documents.

## Current Phase 1 Flow

The current Phase 1 output chain is:

1. PyRosetta input validation and run metadata
2. PyRosetta docking and decoy score export
3. Interface residue extraction
4. Orientation filtering
5. Cluster consensus
6. Cross-state comparison
7. LightDock secondary validation
8. Phase 1 review report and Phase 2 handoff

## Current Primary Evidence Outputs

Primary evidence is the PyRosetta path.

Important files currently produced in or under `output/phase1_ppi/<state>/`:
- `phase1_input_validation_report.json`
- `phase1_input_validation_summary.md`
- `pyrosetta_run_metadata.json`
- `pyrosetta_decoy_scores.csv`
- `pyrosetta_interface_models.csv`
- `pyrosetta_interface_residue_table.csv`
- `orientation_filter_log.csv`
- `ppi_cluster_summary.csv`
- `ppi_hotspot_residues.csv`
- `ppi_interface_patch_table.csv`

Current metadata carried across the main PyRosetta downstream chain includes:
- `receptor_id`
- `construct_type`
- `orientation_validation_status`

## Current Cross-State Outputs

At the Phase 1 output root, the current cross-state comparison layer produces:
- `ppi_patch_cross_state_comparison.csv`
- `ppi_patch_state_robustness.csv`
- `phase1_interface_comparison_report.md`

The comparison and robustness layer now preserves:
- `construct_type`
- `orientation_validation_status`

Cross-state interpretation should still defer to the Phase 1 input validation
outputs if there is any numbering or chain mismatch warning upstream.

## Current Secondary Evidence Outputs

Secondary evidence is the LightDock path.

Important files currently produced in or under `output/phase1_ppi/<state>/lightdock/`:
- `lightdock_run_metadata.json`
- `run_lightdock_<state>.sh`
- `lightdock_interface_support_table.csv`
- `lightdock_model_summary.csv`

Cross-method output is written to:
- `output/phase1_ppi/<state>/cross_method_convergence.csv`

Current LightDock notes:
- LightDock remains secondary evidence only.
- `construct_type` is propagated through the LightDock support path.
- `orientation_validation_status` is currently present but uses the honest
  placeholder value `not_available` in raw LightDock support outputs.
- LightDock-only residues should not be treated as primary Phase 2 patch input
  without additional support.
- AFM is not part of the current active secondary-validation baseline.

## Current Final Phase 1 Outputs

The current Phase 1 review and handoff layer produces:
- `phase1_interface_report.md`
- `phase1_downstream_patch_reference.csv`

The downstream patch reference currently preserves:
- `construct_type`
- `orientation_validation_status`
- robustness labels
- method agreement
- confidence

## Practical Reading Order

If you need to inspect a single Phase 1 run quickly, the current practical
reading order is:

1. `phase1_input_validation_summary.md`
2. `pyrosetta_run_metadata.json`
3. `ppi_cluster_summary.csv`
4. `ppi_patch_state_robustness.csv`
5. `cross_method_convergence.csv`
6. `phase1_interface_report.md`

## Related Notes

For more focused implementation details, see:
- `docs/phase1_pyrosetta_execution_note.md`
- `docs/phase1_ppi_handoff_note.md`
- `docs/phase1_orientation_filter_note.md`
- `docs/phase1_lightdock_validation_note.md`
