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

## Execution Workflow

The LightDock module supports the following workflow:

```
--setup       Generate run scripts and metadata (all 3 states)
--run         Execute generated bash scripts (requires LightDock on PATH)
--extract     Parse LightDock output PDBs + apply orientation filter
--convergence Cross-method comparison (PyRosetta vs LightDock)
--note        Generate dynamic results summary
--all         setup + extract + convergence + note (no execution)
--check       Verify LightDock availability on PATH
```

PBS scripts:
- `config/run_lightdock.pbs` — production (400 swarms, 200 glowworms, 100 steps)
- `config/run_lightdock_test.pbs` — test (50 swarms, 50 glowworms, 50 steps)

### Generated bash script workflow

The generated `run_lightdock_<state>.sh` performs:
1. Pre-flight availability check (all 5 LightDock executables)
2. PDB splitting (chain A/B, preserving TER records)
3. `lightdock3_setup.py` (swarm placement)
4. `lightdock3.py` (optimization)
5. `lgd_generate_conformations.py` (per-swarm pose generation)
6. `lgd_rank.py` (global ranking → `rank_by_scoring.list`)
7. `lgd_cluster_bsas.py` (clustering)
8. Completion marker (`.lightdock_complete`)

## Orientation Filter

LightDock poses now support orientation-aware filtering via a PDB-based
implementation (`compute_orientation_score_from_pdb`) that mirrors the
PyRosetta orientation filter algorithm without requiring PyRosetta:

- Parse CA + CB atoms from PDB files (pure text parsing)
- PCA via numpy SVD to compute sheet-plane normal
- Multi-probe CA→CB consensus to orient normal toward active face
- Dot product with receptor-direction vector → pass/fail/ambiguous

This is the same algorithm as `orientation_filter.py:compute_orientation_score()`
but operates on raw PDB files instead of PyRosetta Pose objects.

The `orientation_validation_status` column in LightDock output tables now
contains actual values (`pass`, `fail`, `ambiguous`, `insufficient_data`)
instead of the previous placeholder `not_available`.

Convergence analysis (`_load_lightdock_residues`) respects this filter:
only `pass` or `not_available` models contribute to residue frequency
calculations, matching PyRosetta's orientation-validated consensus.

## Current Output Files

Per receptor state, the LightDock path is:
- `output/phase1_ppi/<state>/lightdock/`

Files generated or expected in that directory:
- `lightdock_run_metadata.json`
- `run_lightdock_<state>.sh`
- `lightdock_interface_support_table.csv`
- `lightdock_model_summary.csv`
- `.lightdock_complete` (marker after successful execution)

Cross-method output:
- `output/phase1_ppi/<state>/cross_method_convergence.csv`
- `output/phase1_ppi/<state>/cross_method_convergence_summary.json`

Results summary:
- `output/phase1_ppi/phase1_lightdock_validation_results.md`

## Current CSV Baselines

`lightdock_interface_support_table.csv` columns:
- `model_id`, `receptor_id`, `construct_type`
- `orientation_validation_status`, `orientation_score`
- `swarm_id`, `pose_rank`, `scoring_value`
- `chain`, `residue_id`, `residue_num`, `residue_name`, `lobe_label`
- `source`

`lightdock_model_summary.csv` columns:
- `model_id`, `receptor_id`, `construct_type`
- `orientation_validation_status`, `orientation_score`
- `swarm_id`, `pose_rank`, `scoring_value`
- `n_receptor_interface_residues`, `n_partner_interface_residues`
- `n_nlobe_interface_residues`, `n_clobe_interface_residues`
- `receptor_interface_residues`, `partner_interface_residues`
- `source`

`cross_method_convergence.csv` columns:
- `receptor_id`, `construct_type`, `orientation_validation_status`
- `chain`, `residue_id`, `residue_num`, `residue_name`, `lobe_label`
- `in_pyrosetta`, `in_lightdock`
- `pyrosetta_max_occupancy`, `lightdock_frequency`
- `convergence_class`, `method_agreement`

`cross_method_convergence_summary.json` fields:
- `receptor_id`, `n_convergent`, `n_pyrosetta_only`, `n_lightdock_only`
- `n_total`, `jaccard_overlap`, `jaccard_nlobe`, `jaccard_clobe`

## Score Comparison Note

PyRosetta uses REU (dG_separated), LightDock uses DFIRE2 scoring (fastdfire).
These scoring functions are not directly comparable. Cross-method validation
relies on residue-level Jaccard overlap, which is unit-agnostic.

## Metadata Behavior

- `construct_type` is propagated from `lightdock_run_metadata.json` into
  LightDock support tables.
- `orientation_validation_status` is now computed per-pose from the PDB-based
  orientation filter (previously hardcoded as `not_available`).
- `orientation_score` (dot product, -1.0 to +1.0) is stored alongside the
  classification.
- For convergence analysis, LightDock residues are filtered by orientation
  status (only `pass` and `not_available` models contribute).

## Remaining Limitations

- LightDock chain reassignment: If LightDock renumbers residues in output
  PDBs, the orientation filter may fail to find active-face residues.
  Mitigation: falls back to `insufficient_data` classification.
- `cross_method_convergence.csv` is a comparison layer, not a replacement
  for PyRosetta patch definition.

## Operational Note

If a residue appears as `lightdock_only`, it should be interpreted as
secondary method-specific support unless later promoted by additional evidence.
It should not be treated as primary Phase 2 patch input on its own.
