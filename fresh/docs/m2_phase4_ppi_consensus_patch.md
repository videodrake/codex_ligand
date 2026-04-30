# M2.4 Symmetry-aware PPI Consensus Patch

## Scope

M2.4 builds a deterministic EGFR-side PPI consensus patch table from the M2.3
restored contact table. The milestone label follows `fresh/docs/milestone1_foundation_plan.md`:

```text
M2.4 Symmetry-aware consensus patch
```

The more detailed M2 task plan splits this area into `M2-T4` pose QC and
`M2-T5` consensus patch builder. For this implementation, the fixed milestone
label is M2.4 and the implemented scope is the consensus patch builder layer.

## Inputs

Default inputs under `fresh/runs/<run_id>/`:

- `phase1_ppi/tables/ppi_pose_contacts.csv`
- `phase1_ppi/pyrosetta_adapter/pyrosetta_job_manifest.jsonl`
- `manifest/m2_1_ppi_input_manifest.json`

`ppi_pose_contacts.csv` is expected to come from M2.3 mapping restoration.
The M2.2 job manifest is used to find the M1 receptor mapping CSV and re-check
runtime residue identity before consensus rows are written. It also provides the
planned seed/model universe used as the support denominator.

## Outputs

Primary output:

- `output/ppi/ppi_consensus_patch.csv`

M2 table-compatible mirror and support tables:

- `phase1_ppi/tables/ppi_consensus_patch.csv`
- `phase1_ppi/tables/ppi_consensus_residues.csv`
- `phase1_ppi/tables/ppi_patch_state_support.csv`
- `phase1_ppi/tables/ppi_patch_symmetry_support.csv`
- `qc/m2_4_ppi_consensus_patch_qc.csv`
- `manifest/m2_4_ppi_consensus_patch_manifest.json`
- `reports/m2_4_ppi_consensus_patch.md`

## Consensus Schema

`ppi_consensus_patch.csv` columns:

```csv
ppi_patch_id,receptor_id,state_id,state_role,egfr_chain_id,egfr_protomer_id,egfr_residue_number,egfr_residue_name,egfr_contact_count,seed_support_count,pose_support_count,support_fraction,egfr_contact_centroid_x,egfr_contact_centroid_y,egfr_contact_centroid_z,myo1d_contact_residues,myo1d_sheet_support,myo1d_active_face_score,sheet8_support_score,sheet9_support_score,sheet12_support_score,membrane_side_class,evidence_status,warnings
```

Rows are chain- and protomer-resolved. Chain A and chain B are not merged into a
single row. Symmetry support is reported separately in
`ppi_patch_symmetry_support.csv`.

Only accepted pose-QC contact rows are emitted as consensus rows. Rejected and
pending/blank pose-QC rows are excluded from consensus construction; pending
evidence is reported as `PASS_WITH_WARNINGS` with an empty consensus table if no
accepted contact evidence remains.

`support_fraction` is conservative: it is
`pose_support_count / planned_pose_count_for_state_protomer`, where the planned
pose count comes from the M2.2 job manifest `models_per_seed` records. If that
denominator is unavailable, M2.4 falls back to contact-bearing accepted evidence
and records a warning.

## PASS / WARN / FAIL

- `PASS`: restored contact evidence exists, residue mapping re-checks pass, and at least one consensus patch row is written.
- `PASS_WITH_WARNINGS`: M2.3 contact table exists but has no accepted contact evidence, or contacts lack final pose-QC acceptance labels.
- `FAIL`: required input tables are missing, paths escape the run directory, a contact cannot map back to the M1 receptor mapping, or residue-number drift is detected.

## Dry-run Limitation

In a dry-run-only local environment, M2.3 may have no real contact rows. M2.4
then writes a schema-valid empty `output/ppi/ppi_consensus_patch.csv` and reports
`PASS_WITH_WARNINGS`. This is not biological evidence and must not be used to
advance pocket discovery.

## Non-goals

M2.4 does not:

- import or run PyRosetta
- run docking or relaxation
- run LightDock, Vina, fpocket, or P2Rank
- run pocket discovery
- score compounds
- nominate candidates
- submit qsub/PBS/sbatch jobs
- delete cleanup candidates

## Next Milestone

The next recommended milestone is M2.5 ATP-site reference, but only after real
PPI contact evidence and pose-QC acceptance are available or the dry-run warning
is explicitly accepted for documentation-only progression.
