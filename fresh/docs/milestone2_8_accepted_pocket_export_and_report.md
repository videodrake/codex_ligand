# M2.8 Accepted Pocket Export and Milestone 2 Report

## Label and Aliases

User-facing label:

```text
M2.8 - Accepted pocket export and Milestone 2 report
```

Documented aliases:

- `M2.8 Milestone 2 final aggregation`
- `M2 Task 8: Milestone 2 aggregation and report`
- `M2 Prompt 9: accepted pocket export and M2 report`
- `M2-T12 Accepted pocket export and Milestone 2 report`

The implemented scope matches the documented M2.8 scope: final Milestone 2
aggregation, accepted pocket export metadata for Milestone 3, and the M2 report.

## Purpose

M2.8 turns M2.4-M2.7 evidence into a clean handoff package for Milestone 3.
It does not start Milestone 3 and does not prepare ligands, PDBQT files, Vina
jobs, compound scores, or candidate tiers.

M2.8 preserves these scientific guardrails:

- ATP-rejected pockets are not promoted.
- 3GT8_raw-only pockets remain `reference_only` and are not primary M3 targets.
- Soft scores do not rescue hard-gate failures.
- Chain, protomer, state, and residue identity are retained in export files.

## Required Inputs

Production export requires:

- M2.4 PPI consensus patch:
  - `fresh/runs/<run_id>/phase1_ppi/consensus/ppi_consensus_patch_merged.csv`
  - or `fresh/runs/<run_id>/phase1_ppi/tables/ppi_consensus_patch.csv`
  - or `fresh/runs/<run_id>/output/ppi/ppi_consensus_patch.csv`
- M2.5 ATP reference:
  - `fresh/runs/<run_id>/phase2_pockets/atp_reference/atp_site_reference.csv`
- M2.6 merged pocket families:
  - `fresh/runs/<run_id>/phase2_pockets/merged/pocket_candidates_merged.csv`
- M2.7 gate QC and pocket classifications:
  - `fresh/runs/<run_id>/phase2_pockets/gated/pocket_gate_qc.csv`
  - `fresh/runs/<run_id>/phase2_pockets/gated/accepted_pocket_families.csv`
  - `fresh/runs/<run_id>/phase2_pockets/gated/rejected_pocket_families.csv`
- M1 receptor mapping for exported states:
  - `fresh/runs/<run_id>/qc/<state_id>_receptor_mapping.csv`
- M1 membrane frame:
  - `fresh/runs/<run_id>/manifest/membrane_frame.json`

The command also accepts explicit input paths, but every input must resolve
inside `fresh/` or the run directory.

## Transition Gates

M2.8 writes:

```text
fresh/runs/<run_id>/phase2_pockets/final/m2_transition_gate_qc.csv
fresh/runs/<run_id>/phase2_pockets/final/m2_transition_gate_status.json
```

Transition gates:

```text
G1 PPI pose QC complete
G2 accepted/rejected PPI pose reasons recorded
G3 PPI consensus patch generated
G4 ATP reference generated
G5 raw pocket discovery complete
G6 pocket hard gates applied
G7 accepted pocket family table generated
G8 if no accepted pockets, M3 docking is blocked pending manual review
```

Allowed statuses:

```text
PASS
PASS_WITH_WARNINGS
FAIL
NOT_EVALUABLE
```

In production mode, missing required evidence is `FAIL`. In synthetic/smoke
mode, missing evidence produces explicit warnings and schema-valid empty
outputs where possible. Synthetic fixture outputs are always marked as
non-production evidence and do not enable M3 docking.

## Export Classes

M2.8 exports only these tiers:

```text
primary
secondary_review
reference_only
```

Rules:

- `primary`: `accepted_primary_pocket`, all hard gates pass, primary state
  recurrence is sufficient, and the family is not reference-only.
- `secondary_review`: `accepted_secondary_cryptic`, all hard gates pass, at
  least one primary state supports the family, and manual review is required.
- `reference_only`: 3GT8_raw/reference-only evidence. These rows are retained
  for traceability but are not primary Milestone 3 targets.

Rejected gate classes are not exported as primary:

```text
atp_reject
ppi_low_relevance_reject
membrane_geometry_reject
dimer_buried_reject
mapping_reject
dimer_origin_reject
not_evaluable
manual_review
```

## Box Metadata Policy

M2.8 writes box metadata because M2-T12 documents
`accepted_pocket_boxes.csv`. These boxes are metadata only:

- no receptor PDBQT is generated
- no ligand PDBQT is generated
- no docking is run

Default policy when no more detailed geometry is available:

```text
box_center = accepted family center
box_size = 12.0 x 12.0 x 12.0 angstrom
box_margin = 4.0 angstrom
box_source = family_center_minimum_box_m2_8_metadata_only
```

Box QC fails if the center is missing/non-numeric, any size is non-positive,
state/protomer identity is missing, a 3GT8_raw row is marked primary, or an
ATP-failed non-reference row is exported.

## Output Paths

Final aggregation outputs:

```text
fresh/runs/<run_id>/phase2_pockets/final/milestone2_summary.md
fresh/runs/<run_id>/phase2_pockets/final/ppi_to_pocket_evidence_table.csv
fresh/runs/<run_id>/phase2_pockets/final/accepted_pocket_families_for_compound_docking.csv
fresh/runs/<run_id>/phase2_pockets/final/m2_transition_gate_qc.csv
fresh/runs/<run_id>/phase2_pockets/final/m2_transition_gate_status.json
fresh/runs/<run_id>/phase2_pockets/final/m2_final_status.json
```

M3 handoff package outputs:

```text
fresh/runs/<run_id>/phase2_pockets/export_for_m3/accepted_pockets_for_m3.csv
fresh/runs/<run_id>/phase2_pockets/export_for_m3/accepted_pocket_boxes.csv
fresh/runs/<run_id>/phase2_pockets/export_for_m3/accepted_pocket_residue_sets.json
fresh/runs/<run_id>/phase2_pockets/export_for_m3/accepted_pocket_receptor_state_map.csv
fresh/runs/<run_id>/phase2_pockets/export_for_m3/pocket_gate_qc.csv
fresh/runs/<run_id>/phase2_pockets/export_for_m3/ppi_consensus_patch.csv
fresh/runs/<run_id>/phase2_pockets/export_for_m3/m2_final_report.md
fresh/runs/<run_id>/phase2_pockets/export_for_m3/m2_export_manifest.json
fresh/runs/<run_id>/phase2_pockets/export_for_m3/m2_cleanup_report.json
```

Compatibility report:

```text
fresh/runs/<run_id>/reports/milestone2_summary.md
```

M2.8 does not create:

```text
fresh/runs/<run_id>/phase3_compounds/
fresh/runs/<run_id>/phase4_scoring/
```

## Cleanup

Cleanup is report-only in M2.8. No files are deleted in production. The cleanup
report preserves manifests, logs, QC, PPI evidence, ATP reference, merged
pockets, gated outputs, final outputs, export package, and reports.

## CLI

```bash
python -m egfr_myo1d.cli export-m2-results --run-id <RUN_ID>
```

Useful smoke/test option:

```bash
python -m egfr_myo1d.cli export-m2-results --run-id <RUN_ID> --synthetic-fixture true
```

The command does not require PyRosetta, fpocket, P2Rank, Vina, LightDock, or
qsub/sbatch.

## PASS / WARN / FAIL

- `PASS`: all transition gates pass and at least one primary/secondary export
  row is available.
- `PASS_WITH_WARNINGS`: smoke/synthetic evidence is incomplete, no primary or
  secondary pocket is exportable, reference-only evidence is present, or manual
  review is required.
- `FAIL`: production-required PPI, ATP, merged pocket, gate QC, accepted/rejected
  classification, membrane-frame, or exported-state mapping evidence is missing,
  or any input/output path escapes allowed directories.

If no accepted primary or secondary pocket exists, schema-correct export files
are still written and the report states that Milestone 3 docking should not
proceed without manual review or gate-threshold review.

## Non-goals

M2.8 does not run or create:

- Milestone 3 directories
- Vina docking
- ligand preparation
- receptor PDBQT preparation
- PyRosetta docking or relaxation
- LightDock
- P2Rank
- compound scoring
- qsub/PBS/sbatch production submission
- cleanup deletion
- compound candidate nomination

## Next Step

```text
Milestone 3 / M3 Codex Task 0 - readiness audit and M3 directory skeleton
```
