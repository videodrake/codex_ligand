# Phase 1 PPI Handoff Note

This note captures the current repository baseline for the PyRosetta-to-downstream residue handoff.

## Current CSVs

- `ppi_pyrosetta_residues.csv`
- `ppi_pyrosetta_summary.csv`

## Current residue-table baseline

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

## Current summary-table baseline

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

## Important note

`orientation_validation_status` is present as a handoff field, but in the current postprocess path it is usually populated as `not_available`.

This is intentional. The orientation-aware filtering stage has not yet been wired into the current PyRosetta postprocess pipeline, so the repository now preserves the field contract without pretending the evidence already exists.
