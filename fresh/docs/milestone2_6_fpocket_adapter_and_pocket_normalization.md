# M2.6 fpocket Adapter and Pocket Normalization

## Milestone Identity

User-facing label:

```text
M2.6 fpocket adapter and pocket normalization
```

Documented aliases:

- `M2.6 Pocket discovery`
- `M2 Prompt 7: fpocket adapter and pocket normalization`
- `M2-T8 fpocket pocket discovery adapter`
- `M2-T9 Pocket normalization and pocket family merge`
- `M2 Task 6: fpocket adapter and raw pocket parser`

The documented scope is raw pocket discovery and normalization only. M2.6
creates reproducible, state/protomer-aware raw EGFR dimer pocket candidates for
later M2.7 gates. It does not decide accepted pockets.

## Inputs

Default inputs are discovered from:

- `fresh/configs/pocket.yaml`
- `fresh/configs/receptor_states.yaml`
- `fresh/runs/<run_id>/normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb`
- `fresh/runs/<run_id>/normalized/receptors/<state>_dockable_reference_explicit_AB.pdb`
- `fresh/runs/<run_id>/normalized/receptors/<state>_full_frame_explicit_AB.pdb`
- `fresh/runs/<run_id>/normalized/receptors/<state>_runtime_offset_receptor_only.pdb`
- `fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv`
- optional M2.4 `fresh/runs/<run_id>/output/ppi/ppi_consensus_patch.csv`
- optional M2.5 `fresh/runs/<run_id>/phase2_pockets/atp_reference/atp_site_reference.csv`
- optional parser-only fpocket output root under `fresh/` or the run dir

M2.4 consensus patch evidence is required for production interpretation. If it
is missing in production mode, M2.6 fails before invoking fpocket. Parser-only
and synthetic fixture modes may proceed with a warning.

M2.5 ATP reference is required before M2.7 can apply the ATP hard gate. M2.6
does not apply the ATP gate and only warns if the ATP reference is absent.

## fpocket Tool Policy

Primary pocket tool:

```text
fpocket
```

Optional later:

```text
P2Rank, only if path is configured
```

M2.6 implements fpocket only. P2Rank is not implemented or run.

## Receptor Staging

Production fpocket execution stages each receptor PDB into:

```text
fresh/runs/<run_id>/phase2_pockets/raw/<state>_fpocket_pockets/
```

fpocket is run with `subprocess.run([...], shell=False)` from that run-local
directory, so it cannot write into `fresh/data/normalized/` or legacy workflow
paths. Parser-only mode copies supplied fpocket output directories into the same
run-local state directory before parsing.

## Modes

- `production`: stage receptor PDBs, require fpocket binary, run fpocket, parse
  fpocket output.
- `parser_only`: parse existing fpocket output directories; does not require
  fpocket to be installed.
- `synthetic_fixture`: parser-only test provenance label. Synthetic fixtures are
  test data only and are not production evidence.

## Outputs

Canonical M2.6 outputs:

- `fresh/runs/<run_id>/phase2_pockets/raw/<state>_fpocket_raw.csv`
- `fresh/runs/<run_id>/phase2_pockets/raw/<state>_fpocket_pockets/`
- `fresh/runs/<run_id>/phase2_pockets/merged/pocket_candidates_raw.csv`
- `fresh/runs/<run_id>/phase2_pockets/merged/pocket_candidates_merged.csv`

Additional M2 compatibility/reporting outputs:

- `fresh/runs/<run_id>/phase2_pockets/merged/pocket_candidate_members.csv`
- `fresh/runs/<run_id>/phase2_pockets/tables/pocket_candidates_merged.csv`
- `fresh/runs/<run_id>/phase2_pockets/raw/fpocket_discovery_status.json`
- `fresh/runs/<run_id>/phase2_pockets/merged/pocket_normalization_status.json`
- `fresh/runs/<run_id>/manifest/m2_6_fpocket_adapter_manifest.json`
- `fresh/runs/<run_id>/reports/m2_6_fpocket_adapter_and_pocket_normalization.md`
- `fresh/runs/<run_id>/qc/m2_6_fpocket_adapter_qc.csv`

The canonical `raw/` and `merged/` files are the source of truth for M2.7.
The `phase2_pockets/tables/` copy is compatibility/reporting output only.

M2.6 does not write:

- `fresh/runs/<run_id>/phase2_pockets/tables/pocket_gate_qc.csv`
- `fresh/runs/<run_id>/phase2_pockets/tables/pocket_rejects.csv`
- `fresh/runs/<run_id>/phase2_pockets/tables/accepted_pockets_for_m3.csv`
- `fresh/runs/<run_id>/phase2_pockets/export_for_m3/accepted_pocket_boxes.csv`

## Schemas

Per-state raw fpocket table:

```csv
state_id,state_role,receptor_id,receptor_path,source_tool,source_tool_version,source_pocket_id,source_pocket_dir,fpocket_score,fpocket_druggability_score,fpocket_volume,num_alpha_spheres,num_pocket_atoms,centroid_x,centroid_y,centroid_z,raw_chain_ids,raw_residue_ids,raw_residue_names,egfr_chain_ids,egfr_protomer_ids,mapped_uniprot_residues,mapping_status,parser_status,evidence_status,warnings
```

`pocket_candidates_raw.csv`:

```csv
candidate_pocket_id,state_id,state_role,receptor_id,receptor_path,source_tool,source_pocket_id,egfr_chain_ids,egfr_protomer_ids,multi_protomer,pocket_center_x,pocket_center_y,pocket_center_z,pocket_extent_x,pocket_extent_y,pocket_extent_z,pocket_volume,fpocket_score,fpocket_druggability_score,pocket_residue_count,receptor_residues,mapped_uniprot_residues,mapping_status,raw_candidate_status,evidence_status,warnings
```

`raw_candidate_status` is always:

```text
raw_ungated
```

`pocket_candidates_merged.csv`:

```csv
pocket_family_id,representative_candidate_pocket_id,representative_state_id,representative_state_role,representative_protomer_id,member_candidate_pocket_ids,member_state_ids,member_state_roles,member_protomer_ids,primary_state_count,reference_state_count,raw_candidate_count,family_center_x,family_center_y,family_center_z,centroid_spread_angstrom,residue_union_uniprot,residue_intersection_uniprot,max_pairwise_residue_jaccard,median_pairwise_residue_jaccard,merge_method,merge_threshold_centroid_angstrom,merge_threshold_residue_jaccard,merge_status,evidence_status,warnings
```

`merge_status` is always raw and ungated:

```text
raw_family_ungated
```

## Merge Policy

Defaults are configured in `fresh/configs/pocket.yaml`:

```text
centroid_merge_threshold_angstrom = 6.0
residue_jaccard_merge_threshold = 0.30
merge_rule = centroid_distance <= threshold OR residue_jaccard >= threshold
```

Pocket IDs alone never drive family merge. The merge keeps member pocket IDs,
state roles, and protomer IDs explicit.

## Chain, Protomer, and Mapping Safeguards

Pocket residues are parsed by chain and receptor residue number. If an M1
mapping CSV is available, it is used as the source of truth for UniProt/source
residue numbers and protomer IDs. Runtime-number mapping, missing mappings, and
residue-name mismatches are recorded in row warnings rather than silently
dropping the pocket.

Chain A and Chain B are not collapsed. Multi-protomer pockets are marked
explicitly.

## 3GT8_raw Handling

`3GT8_raw` may be parsed and normalized as reference-only raw evidence. It is
reported with `state_role=reference`. It does not count toward
`primary_state_count`, and a `3GT8_raw`-only family is marked as not promotable
by M2.6.

## PASS / WARN / FAIL

- `PASS`: raw pocket candidates are parsed and normalized without warnings.
- `PASS_WITH_WARNINGS`: parser-only smoke execution lacks M2.4/M2.5 evidence,
  raw pockets are absent, mappings are partial, residue names drift, or ATP
  reference is missing for later M2.7.
- `FAIL`: production mode lacks M2.4 consensus evidence, fpocket is required
  but unavailable or fails, production receptor inputs are missing, or any
  input/output path escapes `fresh/` or the run dir.

## Non-goals

M2.6 does not run PyRosetta docking/relaxation, LightDock, P2Rank, Vina,
qsub/PBS/sbatch production, cleanup deletion, compound docking, compound
scoring, candidate nomination, M2.7 gates, or accepted pocket export.

## Next Step

Next milestone:

```text
M2.7 Membrane/dimer/PPI pocket gates
```
