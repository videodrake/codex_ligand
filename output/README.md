# Output Documentation Index

This folder is the output root for the `codex_ligand` workspace. Use this README to find the correct project root and reading order quickly.

## Folder Purpose

- `egfr_myo1d_vina/`: Routine baseline project outputs from the Vina-centered workflow. Canonical runtime outputs remain here, and the additive step folders live here as derived interpretation views when the production lane refreshes them.
- `phase1_ppi/`: Phase 1 PPI evidence reports.
- `phase2_pockets/`: Phase 2 pocket proposal and druggability reports.
- `phase3_docking/`: Phase 3 diversity-aware docking reports.
- `phase4_perturbation/`: Phase 4 perturbation ranking and prioritization reports.

## First Checks For New Runs

1. Open `egfr_myo1d_vina/step_index.md` when it exists. That is the first human-readable entry point for a completed production run.
2. If you need canonical files directly, open `egfr_myo1d_vina/vina_pocket_table.csv`, `egfr_myo1d_vina/valid_sites.csv`, and `egfr_myo1d_vina/project_report.txt`.
3. Open the relevant phase folder README below for the phase-separated scientific branches.
4. Cross-check artifact meaning in `../docs/output_artifact_map.md`.

## Canonical Vs Derived

- Canonical runtime outputs remain under `output/{project}/` and stay the source of truth.
- `step1_vina_raw/` through `step7_validate/` are derived interpretation views.
- Large raw PyRosetta directories are referenced from manifests and index files rather than duplicated into the step view.

## Phase README Files

- [phase1_ppi/README.md](phase1_ppi/README.md)
- [phase2_pockets/README.md](phase2_pockets/README.md)
- [phase3_docking/README.md](phase3_docking/README.md)
- [phase4_perturbation/README.md](phase4_perturbation/README.md)

## Cross-Reference

- [../docs/output_artifact_map.md](../docs/output_artifact_map.md): Canonical artifact semantics and interpretation priority.
- [../docs/current_pipeline_status.md](../docs/current_pipeline_status.md): Current baseline summary before interpreting outputs.
