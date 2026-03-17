# Phase 1 PPI Handoff Note

This note captures the current repository baseline for the legacy
PyRosetta-postprocess handoff tables under `output/egfr_myo1d_vina/ppi/`.

## Scope and relationship to sibling docs

- This file documents legacy handoff tables: `ppi_pyrosetta_residues.csv` and `ppi_pyrosetta_summary.csv`.
- For the structured Phase 1 branch under `output/phase1_ppi/`, read `docs/phase1_output_chain_note.md` first.
- For phase-wide artifact meaning and priority, read `docs/output_artifact_map.md`.

## Legacy postprocess CSVs (routine baseline)

- `ppi_pyrosetta_residues.csv`
- `ppi_pyrosetta_summary.csv`

## Legacy residue-table baseline

`ppi_pyrosetta_residues.csv` currently includes:

- `receptor_id`
- `partner_id`
- `source`
- `chain`
- `residue_id`
- `residue_num`
- `residue_name`
- `lobe_label`
- `construct_type`
- `orientation_validation_status`
- `frequency_final_ranking`
- `frequency_cluster_summary`
- `n_models_final_ranking`
- `occupancy`
- `mean_interface_delta_e`
- `best_interface_delta_e`

## Legacy summary-table baseline

`ppi_pyrosetta_summary.csv` currently includes:

- `receptor_id`
- `partner_id`
- `source`
- `construct_type`
- `orientation_validation_status`
- `n_final_models`
- `n_clusters`
- `n_interface_residues`
- `n_nlobe_interface_residues`
- `n_clobe_interface_residues`
- `top_residues`
- `best_dg`
- `mean_dg`
- `best_dsasa`

## Orientation field contract across both paths

`orientation_validation_status` is present in both legacy and structured Phase 1 contracts so downstream readers can use one field name.

Interpretation by source:

- Legacy postprocess path (`output/egfr_myo1d_vina/ppi/`): value is usually `not_available`.
- Structured Phase 1 path (`output/phase1_ppi/`): value is expected to reflect orientation filtering output (`pass`, `fail`, or `ambiguous`) when the orientation module is run.

This split is intentional and should not be treated as a contradiction.
