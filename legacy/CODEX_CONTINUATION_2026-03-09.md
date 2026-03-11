# Codex Continuation Note

Read this file first when resuming work in another environment.

## Current State

This repository was refactored through Task Group 4 for the Vina-centered workflow.

Implemented scope:
- Task Group 1: project-level config support in `run_docking.py`
- Task Group 2: receptor-sequential / ligand-parallel dispatch with configurable `max_workers`
- Task Group 3: pose parsing and receptor contact extraction
- Task Group 4: receptor-level pocket clustering and pocket summarization
- Task Group 5: cross-receptor pocket comparison (`compare_pockets.py`)
- Task Group 6: PPI output standardization (`extract_ppi_residues.py`)
- Task Group 7: reporting and manual review exports (`generate_report.py`)
- Task Group 8: validation, regression checks, and handoff readiness (`validate_outputs.py`)

All planned Task Groups (0-8) are complete. Remaining:
- production Vina runtime validation on the real Ubuntu server

## Files Added Or Updated

Core updated:
- `run_docking.py`
- `README.md`

Task Group 1-4 modules added:
- `config/example-project.yaml`
- `parse_vina_results.py`
- `extract_contacts.py`
- `cluster_pockets.py`
- `summarize_pockets.py`

Task Group 5 module added:
- `compare_pockets.py`

Task Group 6 module added:
- `extract_ppi_residues.py`
- Mock PyRosetta data: `smoke_test/mock_pyrosetta/3GT8_raw/`

Task Group 7 module added:
- `generate_report.py`

Task Group 8 module added:
- `validate_outputs.py`
- `smoke_test/config_real_pdbs.yaml` (validation config with real PDB paths)

Local smoke-test assets added:
- `smoke_test/config.yaml`
- `smoke_test/input/*`
- `smoke_test/output/smoke_vina/*`

## Important Execution Notes

- The VMware workspace is not the production server.
- `max_workers=16` is kept as the intended server-side operating default.
- Do not tune performance from this local environment.
- Treat this environment as a functional validation workspace only.

## Real Input Files Currently Present

The user placed real source files under:
- `smoke_test/input/original/3gt8_dimer_anp.pdb`
- `smoke_test/input/original/-10_38-48.pdb`
- `smoke_test/input/original/-10_85-100.pdb`
- `smoke_test/input/original/173940_ligand.sdf`
- `smoke_test/input/original/97806_ligand.sdf`
- `smoke_test/input/original/VAX-C12_0_ligand.sdf`

These are real source files, not the earlier dummy structures.

## Environment Setup Completed Here

Local virtualenv created:
- `.venv`

Verified installed packages:
- `pyyaml`
- `numpy`
- `pandas`
- `matplotlib`
- `rdkit`

Notes:
- Python venv was created with Python 3.9 because the default Python 3.13 environment had `ensurepip` issues in this workspace.
- `vina` Python package did not install cleanly here.
- No system `vina` executable was detected in this workspace.

## Current Data Contracts

### `vina_pose_table.csv`
- `receptor_id`
- `ligand_id`
- `pose_rank`
- `affinity`
- `rmsd_lb`
- `rmsd_ub`
- `centroid_x`
- `centroid_y`
- `centroid_z`
- `raw_pose_file`
- `pocket_id`
- `contact_residues`
- `n_contact_residues`

### `vina_pocket_table.csv`
- `receptor_id`
- `pocket_id`
- `centroid_x`
- `centroid_y`
- `centroid_z`
- `n_pose`
- `n_ligand`
- `best_affinity`
- `mean_affinity`
- `union_contact_residues`
- `top_residues`

### `vina_drug_pocket_map.csv`
- `receptor_id`
- `ligand_id`
- `dominant_pocket_id`
- `dominant_pocket_pose_count`
- `dominant_pocket_fraction`
- `best_affinity`
- `best_pose_rank`
- `top_pose_residues`
- `alternative_pockets`
- `is_multimodal_binding`

### `vina_pocket_comparison.csv` (Task Group 5)
- `receptor_a`, `pocket_a`, `receptor_b`, `pocket_b`
- `centroid_dist` — pocket centroid distance in Angstrom
- `residue_jaccard` — Jaccard similarity of union_contact_residues (normalized)
- `residue_overlap_coeff` — overlap coefficient (lenient when pocket sizes differ)
- `shared_residues`, `n_shared_residues`
- `residues_only_a`, `residues_only_b`
- `n_residues_a`, `n_residues_b`
- `shared_ligands`, `n_shared_ligands`, `n_ligands_a`, `n_ligands_b`
- `affinity_a`, `affinity_b`, `n_pose_a`, `n_pose_b`
- `same_patch_candidate` — auxiliary flag (centroid < 8 A AND jaccard >= 0.3 OR overlap >= 0.5)

### `ppi_pyrosetta_residues.csv` (Task Group 6)
- `receptor_id`, `source` (= "pyrosetta_ppi")
- `residue_id` — normalized (chain stripped, CHARMM names converted)
- `residue_num`
- `frequency_final_ranking` — how many final models contact this residue
- `frequency_cluster_summary` — how many cluster representatives contact it
- `n_models_final_ranking`
- `occupancy` — frequency / n_models
- `mean_interface_delta_e`, `best_interface_delta_e` — from InterfaceEnergies CSVs if available

### `ppi_pyrosetta_summary.csv` (Task Group 6)
- `receptor_id`, `source`, `n_final_models`, `n_clusters`, `n_interface_residues`
- `top_residues`, `best_dg`, `mean_dg`, `best_dsasa`

### `ppi_afm_residues.csv` (Task Group 6, stub)
- `receptor_id`, `source` (= "alphafold_multimer")
- `residue_id`, `residue_num`, `min_ca_distance`

## Current Rules

- Receptors are limited to exactly three IDs in project config:
  - `3GT8_raw`
  - `3GT8_cl38_48`
  - `3GT8_cl85_100`
- Residue strings are stored as `A:MET971` style.
- Contact cutoff and pocket clustering cutoff are separate settings.
- Clustering is receptor-wide, not ligand-local.
- Pocket assignment is deterministic via sorted greedy assignment.

## Recommended First Step After Moving To Ubuntu

1. Read this file.
2. Read `README.md`.
3. Inspect `smoke_test/input/original/` or move those real inputs into the preferred server-side input layout.
4. Rebuild or verify the Python environment on Ubuntu.
5. Install or expose a working Vina runtime.
6. Update the project config to point at the real receptor/ligand preparation paths.
7. Run a small end-to-end functional test before any high-load batch.

## Recommended Next Work Item

All Task Groups (0-8) are complete. Deployment steps:

1. Wire real input files into project config on Ubuntu
2. Establish SDF -> PDBQT preparation path (obabel or mk_prepare_ligand)
3. Install/expose Vina runtime
4. Run full end-to-end with real data: Vina docking → TG3-7 postprocess → TG8 validation
5. Point ppi.pyrosetta_result_dirs at actual PyRosetta output directories
6. Point ppi.afm_models at actual AlphaFold-Multimer model PDBs (when available)
7. Run `python validate_outputs.py --config <config> --repo-root .` to verify everything

Residue numbering note: residue numbers are consistent across receptors (699-1007 overlap), chain IDs differ (A/B vs X), CHARMM names (HSD→HIS) normalized in compare_pockets.py and extract_ppi_residues.py. validate_outputs.py detects and warns about chain ID differences.

## Resume Prompt

Use this prompt in the next Codex environment:

`Read README.md and CODEX_CONTINUATION_2026-03-09.md first. All Task Groups 0-8 are complete. Deploy on Ubuntu: wire real inputs, set up Vina runtime, run end-to-end, then validate with validate_outputs.py.`
