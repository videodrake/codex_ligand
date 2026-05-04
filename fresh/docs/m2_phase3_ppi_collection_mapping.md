# M2.3 PPI Pose Collection and Mapping Restoration

## Scope

M2.3 starts the PPI-side output interpretation layer without running any
scientific engine. It consumes the M2.2 PyRosetta dry-run job manifest and
optionally a run-local raw contact CSV, then restores EGFR runtime residue
numbers through the M1 receptor mapping CSV.

Current implementation is intentionally staged:

- Collects M2.2 dry-run job status into `phase1_ppi/tables/ppi_raw_pose_table.csv`.
- Writes `phase1_ppi/tables/ppi_pose_contacts.csv` with the M2 plan contact schema.
- Writes `phase1_ppi/tables/ppi_mapping_restoration_qc.csv`.
- Writes `phase1_ppi/tables/unmapped_contacts.csv` for contacts that cannot be restored.
- Writes `manifest/m2_3_ppi_collection_manifest.json`.
- Writes `reports/m2_3_ppi_collection_mapping.md`.

## Inputs

Required:

- `fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/pyrosetta_job_manifest.jsonl`

Optional:

- A raw contact CSV under the same `fresh/runs/<run_id>/` directory.

The optional raw contact CSV is for mapping restoration only. It is not a
docking result claim, not a pose acceptance table, and not a consensus patch.

## Mapping Policy

EGFR contact residue identity is restored from:

```text
runtime_resseq -> source_resseq/source_icode/source_resname
```

via the M1 receptor mapping CSV recorded in the M2.2 job manifest. The current
M1 mapping schema does not yet have a separate UniProt column, so
`egfr_uniprot_residue` is populated from the restored source residue and the QC
table records that policy explicitly.

Every restored contact keeps a `protomer_id`. A B-protomer runtime residue such
as source +1000 must map back to the original residue number before downstream
QC or consensus can use it.

## Non-Goals

M2.3 does not:

- import PyRosetta
- run docking or relaxation
- run LightDock, Vina, fpocket, or P2Rank
- run scheduler submission
- parse engine pose PDBs as scientific output
- classify accepted/rejected poses
- score compounds or nominate candidates
- delete cleanup candidates

## Expected Statuses

- `PASS`: job manifest exists and supplied raw contacts restore cleanly.
- `PASS_WITH_WARNINGS`: job manifest exists but no raw contact table was supplied, or non-blocking warnings exist.
- `FAIL`: M2.2 job manifest is missing, raw contact schema is invalid, or unmapped contact fraction exceeds the configured threshold.
