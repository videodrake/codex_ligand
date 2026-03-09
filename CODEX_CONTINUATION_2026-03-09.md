# Codex Continuation Note

Read this file first when resuming work in another environment.

## Current State

This repository was refactored through Task Group 4 for the Vina-centered workflow.

Implemented scope:
- Task Group 1: project-level config support in `run_docking.py`
- Task Group 2: receptor-sequential / ligand-parallel dispatch with configurable `max_workers`
- Task Group 3: pose parsing and receptor contact extraction
- Task Group 4: receptor-level pocket clustering and pocket summarization

Not implemented yet:
- Task Group 5: cross-receptor pocket comparison
- report generation beyond existing README notes
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

Before Task Group 5:
- wire the real input files into a project config
- establish the SDF -> PDBQT preparation path on Ubuntu
- verify pose parsing, contacts, clustering, and summary outputs on a small real sample
- verify residue numbering consistency across the three receptor PDBs

After that, begin Task Group 5:
- `compare_pockets.py`
- raw cross-receptor metrics only first
- no premature report layer

## Resume Prompt

Use this prompt in the next Codex environment:

`Read README.md and CODEX_CONTINUATION_2026-03-09.md first, then continue from Task Group 4 completion state. Do not redo Tasks 1-4. First wire the real receptor/ligand inputs from smoke_test/input/original or their moved Ubuntu equivalents into a project config, set up the actual Vina/SDF->PDBQT execution path on Ubuntu, run a small functional validation, and only then proceed toward Task Group 5.`
