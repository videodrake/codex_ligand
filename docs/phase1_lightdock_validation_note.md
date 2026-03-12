# Phase 1 LightDock Validation Note

This note records the current implementation baseline for Phase 1 Task Group 1.4.
It is intentionally narrow: the goal is to document what the repository
currently exports, not to restate the full PRD.

## Current Role

LightDock is treated as secondary evidence only.
It supports Phase 1 receptor-side patch review, but it does not replace
PyRosetta as the primary evidence source for downstream handoff.

Current evidence hierarchy:
- Primary: PyRosetta interface consensus and cross-state robustness
- Secondary: LightDock interface support and cross-method convergence
- Legacy optional only: AlphaFold-Multimer parser path if explicitly re-enabled later

## Current Output Files

Per receptor state, the LightDock path is currently:
- `output/phase1_ppi/<state>/lightdock/`

Current files generated or expected in that directory:
- `lightdock_run_metadata.json`
- `run_lightdock_<state>.sh`
- `lightdock_interface_support_table.csv`
- `lightdock_model_summary.csv`

Cross-method output is written to:
- `output/phase1_ppi/<state>/cross_method_convergence.csv`

## Current CSV Baselines

`lightdock_interface_support_table.csv` currently includes:
- `model_id`
- `receptor_id`
- `construct_type`
- `orientation_validation_status`
- `swarm_id`
- `pose_rank`
- `scoring_value`
- `chain`
- `residue_id`
- `residue_num`
- `residue_name`
- `lobe_label`
- `source`

`lightdock_model_summary.csv` currently includes:
- `model_id`
- `receptor_id`
- `construct_type`
- `orientation_validation_status`
- `swarm_id`
- `pose_rank`
- `scoring_value`
- `n_receptor_interface_residues`
- `n_partner_interface_residues`
- `n_nlobe_interface_residues`
- `n_clobe_interface_residues`
- `receptor_interface_residues`
- `partner_interface_residues`
- `source`

`cross_method_convergence.csv` currently includes:
- `receptor_id`
- `construct_type`
- `orientation_validation_status`
- `chain`
- `residue_id`
- `residue_num`
- `residue_name`
- `lobe_label`
- `in_pyrosetta`
- `in_lightdock`
- `pyrosetta_max_occupancy`
- `lightdock_frequency`
- `convergence_class`
- `method_agreement`

## Metadata Behavior

Current metadata behavior is:
- `construct_type` is propagated from `lightdock_run_metadata.json` into
  LightDock support tables.
- `orientation_validation_status` is currently emitted as `not_available` in
  raw LightDock support tables.
- `construct_type` and `orientation_validation_status` are propagated into
  `cross_method_convergence.csv`.
- For PyRosetta-supported residues, `orientation_validation_status` comes from
  `ppi_interface_patch_table.csv`.
- For LightDock-only residues, `orientation_validation_status` is currently
  recorded as `not_available`.

## Current Limitations

The current implementation still has a few honest limits:
- LightDock extraction does not yet apply an orientation-aware filter that is
  equivalent to the PyRosetta orientation filter path.
- LightDock raw support tables now expose `orientation_validation_status`, but
  the current value is a placeholder (`not_available`) until an equivalent
  LightDock orientation filter path exists.
- `cross_method_convergence.csv` is a comparison layer, not a replacement for
  PyRosetta patch definition.

## Operational Note

If a residue appears as `lightdock_only`, it should be interpreted as
secondary method-specific support unless later promoted by additional evidence.
It should not be treated as primary Phase 2 patch input on its own.
