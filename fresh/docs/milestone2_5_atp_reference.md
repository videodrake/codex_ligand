# M2.5 ATP-site Reference

## Milestone Identity

User-facing label:

```text
M2.5 ATP-site reference
```

Documented aliases:

- `M2 Prompt 6: ATP reference builder`
- `M2-T7 ATP reference builder`
- `M2 Task 5: ATP reference generation`

The purpose is to create an ATP-site reference that later pocket/pose gates can
use to block ATP-cleft or ATP-migrated results from being promoted as
EGFR-MYO1D PPI-disruptive non-ATP candidates.

## Inputs

Default inputs are discovered from:

- `fresh/configs/atp_reference.yaml`
- `fresh/configs/receptor_states.yaml`
- `fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv`
- `fresh/runs/<run_id>/normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb`
- optional `3GT8_raw` ligand-bearing reference PDB under the run directory or `fresh/data/raw/receptors/`
- optional M2.4 `fresh/runs/<run_id>/output/ppi/ppi_consensus_patch.csv`

The M2.4 consensus patch is not required to build the ATP reference. If it is
missing or empty, M2.5 records a warning because ATP reference generation is a
hard-gate setup step, not a pocket discovery step.

## Reference Modes

M2.5 supports:

- `ligand_bearing_reference`: detects configured ATP-like ligand residue names
  from a ligand-bearing PDB and collects nearby receptor residues. Neighbor
  residue identity is restored through the M1 receptor mapping CSV before it is
  written as UniProt/source numbering.
- `configured_residue_set`: uses `fresh/configs/atp_reference.yaml`.
- `synthetic_fixture`: used by tests to mark synthetic PDB/mapping evidence.
- `missing_or_unresolved`: represented by `FAIL` status when neither ligand nor
  configured residues exist.

Allowed ATP-like ligand residue names are configured in
`fresh/configs/atp_reference.yaml`; arbitrary HETATM ligands are not treated as
ATP. A configured fallback residue set must meet its
`minimum_residue_count` unless ligand-derived mapped evidence is available.

## Outputs

Canonical ATP reference outputs:

- `fresh/runs/<run_id>/phase2_pockets/atp_reference/atp_site_reference.csv`
- `fresh/runs/<run_id>/phase2_pockets/atp_reference/atp_site_residue_mapping.csv`
- `fresh/runs/<run_id>/phase2_pockets/atp_reference/atp_site_centroid_by_state.csv`

Additional M2 compatibility/reporting outputs:

- `fresh/runs/<run_id>/phase2_pockets/atp_reference/atp_site_mapping_qc.csv`
- `fresh/runs/<run_id>/phase2_pockets/atp_reference/atp_reference_status.json`
- `fresh/runs/<run_id>/phase2_pockets/tables/atp_site_reference.csv`
- `fresh/runs/<run_id>/manifest/m2_5_atp_reference_manifest.json`
- `fresh/runs/<run_id>/reports/m2_5_atp_reference.md`

The canonical `atp_reference/` files are source-of-truth for later M2.6/M2.7.
The `phase2_pockets/tables/` copy is a reporting compatibility copy.

## Schemas

`atp_site_reference.csv`:

```csv
atp_site_id,reference_source_id,reference_source_path,reference_mode,ligand_resname,ligand_chain_id,ligand_residue_number,uniprot_residue_number,residue_name,atp_site_region,evidence_source,evidence_status,warnings
```

`atp_site_residue_mapping.csv`:

```csv
state_id,receptor_id,egfr_chain_id,egfr_protomer_id,uniprot_residue_number,reference_residue_number,receptor_residue_number,reference_residue_name,receptor_residue_name,mapping_status,residue_name_match,ca_x,ca_y,ca_z,heavy_atom_count,evidence_status,warnings
```

`atp_site_centroid_by_state.csv`:

```csv
state_id,receptor_id,egfr_chain_id,egfr_protomer_id,centroid_method,atp_centroid_x,atp_centroid_y,atp_centroid_z,residue_count,atom_count,reference_radius_angstrom,evidence_status,warnings
```

## ATP Gate Policy

Configured defaults:

```text
atp_overlap_fraction_reject_threshold = 0.15
atp_centroid_distance_reject_threshold_angstrom = 10.0
```

Future gate rule:

```text
future non_atp_pass = false if:
    atp_overlap_fraction >= threshold
    OR atp_centroid_distance <= threshold
```

M2.5 records this policy but does not apply it to pocket candidates.

## PASS / WARN / FAIL

- `PASS`: ATP reference rows, primary-state mappings, and primary-state
  centroids are available without warnings.
- `PASS_WITH_WARNINGS`: reference rows and primary-state centroids exist but
  ligand evidence is absent, M2.4 runtime contact evidence is absent, some
  mapped residues are missing, or residue-name mismatches are detected.
- `FAIL`: no ligand-bearing reference and no configured residue set are
  available, the config is missing, input paths escape `fresh/` / the run dir,
  no primary-state ATP centroid can be built, or the configured fallback residue
  set is too sparse for the hard gate.

## 3GT8_raw Handling

`3GT8_raw` may contribute ligand/reference evidence or frame-transfer
comparison evidence. It is always reference/control evidence and must not be
counted as a primary membrane-validated receptor state.

## Non-goals

M2.5 does not run PyRosetta, LightDock, Vina, fpocket, P2Rank, qsub/PBS/sbatch,
cleanup deletion, pocket discovery, compound scoring, or candidate nomination.

## Next Step

Next milestone:

```text
M2.6 fpocket adapter and pocket normalization
```
