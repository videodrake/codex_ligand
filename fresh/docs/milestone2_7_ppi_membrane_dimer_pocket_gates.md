# M2.7 PPI/membrane/dimer Pocket Gates

## Milestone Identity

User-facing label:

```text
M2.7 PPI/membrane/dimer pocket gates
```

Documented aliases:

- `M2.7 Membrane/dimer/PPI pocket gates`
- `M2 Prompt 8: PPI/membrane/dimer gates`
- `M2 Task 7: membrane/dimer/PPI pocket gates`
- `M2-T10 PPI-pocket relationship classifier`
- `M2-T11 Membrane/lateral/dimer accessibility gates`

The documented scope is hard-gate classification for M2.6 raw pocket families.
M2.7 does not perform M2.8 final aggregation/export packaging and does not
start compound docking.

## Inputs

Default inputs are discovered from:

- `fresh/configs/pocket.yaml`
- `fresh/configs/receptor_states.yaml`
- `fresh/runs/<run_id>/output/ppi/ppi_consensus_patch.csv`
- `fresh/runs/<run_id>/phase1_ppi/tables/ppi_consensus_patch.csv`
- `fresh/runs/<run_id>/phase2_pockets/atp_reference/atp_site_reference.csv`
- `fresh/runs/<run_id>/phase2_pockets/atp_reference/atp_site_centroid_by_state.csv`
- `fresh/runs/<run_id>/phase2_pockets/merged/pocket_candidates_raw.csv`
- `fresh/runs/<run_id>/phase2_pockets/merged/pocket_candidates_merged.csv`
- `fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv`
- `fresh/runs/<run_id>/normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb`
- `fresh/runs/<run_id>/manifest/membrane_frame.json`

Production M2.7 gating fails if M2.4 PPI consensus, M2.5 ATP reference, M2.6
merged pocket families, receptor mapping/coordinates, or membrane frame are not
available. Synthetic/smoke fixtures may run with explicit warnings and are not
production evidence.

## Hard Gates

M2.7 applies hard gates separately from soft priority scores:

```text
G0 dimer_origin_pass
G1 non_atp_pass
G2 ppi_relationship_pass
G3 lower_lateral_pass
G4 dimer_accessibility_pass
G5 primary_state_support_or_review_flag
```

A hard-gate-failed pocket family cannot be rescued by a soft score.

## ATP Gate

Configured defaults in `fresh/configs/pocket.yaml`:

```text
atp_overlap_fraction_reject_threshold = 0.15
atp_centroid_distance_reject_angstrom = 10.0
```

Reject rule:

```text
non_atp_pass = false if:
    atp_overlap_fraction >= threshold
    OR atp_centroid_distance <= threshold
```

No pocket is promoted in production without a valid ATP reference.

## PPI Relationship Classes

M2.7 writes one relationship row per pocket family:

- `orthosteric`: pocket residues overlap PPI patch residues.
- `rim`: no overlap, but residue/atom distance is near the PPI patch boundary.
- `allosteric_near`: no direct overlap, but residue distance supports a nearby
  same lower/lateral neighborhood.
- `low_relevance`: no overlap or near-distance support.
- `mapping_ambiguous`: mapping cannot support reliable classification.
- `not_evaluable`: required evidence is missing.

Centroid distance alone is not sufficient to classify a pocket as PPI-relevant.

## Membrane Geometry

M2.7 reads `membrane_frame.json` for the membrane normal, origin/JM anchor, and
protomer centroids. It does not hard-code a fallback membrane normal.

The lateral score is:

```text
u_out = projected outward vector from dimer center to protomer centroid
r_pocket = projected vector from protomer centroid to pocket centroid
lateral_score = dot(normalize(u_out), normalize(r_pocket))
```

Default pass:

```text
lateral_score >= 0.15
```

If lower/cytosolic direction is ambiguous, the lower gate fails with a warning.

## Dimer Accessibility

M2.7 detects inter-protomer interface residues from dimer receptor heavy atoms.
Default interface cutoff:

```text
interface_heavy_atom_cutoff_angstrom = 5.0
```

Default reject:

```text
dimer_accessibility_pass = false if interface_overlap_fraction >= 0.25
```

SASA is not computed in this increment:

```text
sasa_status = not_computed_m2_7_optional
```

## 3GT8_raw Handling

`3GT8_raw` may contribute reference evidence. It is not counted as a primary
membrane-validated receptor state. A `3GT8_raw`-only family is marked
`reference_only`, never `accepted_primary_pocket`.

## Gate Classes

- `accepted_primary_pocket`
- `accepted_secondary_cryptic`
- `reference_only`
- `atp_reject`
- `membrane_geometry_reject`
- `dimer_buried_reject`
- `ppi_low_relevance_reject`
- `dimer_origin_reject`
- `mapping_reject`
- `not_evaluable`
- `manual_review`

Default primary support thresholds:

```text
required_primary_state_count_for_primary = 2
allow_secondary_if_primary_state_count = 1
```

## Outputs

Canonical gated outputs:

- `fresh/runs/<run_id>/phase2_pockets/gated/pocket_gate_qc.csv`
- `fresh/runs/<run_id>/phase2_pockets/gated/accepted_pocket_families.csv`
- `fresh/runs/<run_id>/phase2_pockets/gated/rejected_pocket_families.csv`
- `fresh/runs/<run_id>/phase2_pockets/gated/pocket_ppi_relationship.csv`
- `fresh/runs/<run_id>/phase2_pockets/gated/pocket_membrane_geometry.csv`
- `fresh/runs/<run_id>/phase2_pockets/gated/pocket_dimer_accessibility.csv`
- `fresh/runs/<run_id>/phase2_pockets/gated/pocket_gate_status.json`

Compatibility/reporting outputs:

- `fresh/runs/<run_id>/phase2_pockets/tables/pocket_gate_qc.csv`
- `fresh/runs/<run_id>/phase2_pockets/tables/pocket_rejects.csv`
- `fresh/runs/<run_id>/phase2_pockets/tables/accepted_pockets_for_m3.csv`
- `fresh/runs/<run_id>/phase2_pockets/tables/pocket_ppi_relationship.csv`
- `fresh/runs/<run_id>/phase2_pockets/tables/pocket_membrane_geometry.csv`
- `fresh/runs/<run_id>/phase2_pockets/tables/pocket_dimer_accessibility.csv`
- `fresh/runs/<run_id>/phase2_pockets/tables/pocket_gate_status.json`

`accepted_pockets_for_m3.csv` is a gate-result compatibility table only. M2.7
does not create the M2.8 export package under `phase2_pockets/export_for_m3/`.

## PASS / WARN / FAIL

- `PASS`: families are gated without warnings.
- `PASS_WITH_WARNINGS`: synthetic/smoke inputs are missing, mappings are
  partial, SASA is not computed, or reference-only evidence is present.
- `FAIL`: production-required PPI, ATP, pocket-family, membrane-frame, or
  receptor mapping/coordinate inputs are missing, or an input path escapes
  `fresh/` / the run dir.

## Non-goals

M2.7 does not run Vina, ligand preparation, PyRosetta docking/relaxation,
LightDock, P2Rank, qsub/PBS/sbatch production, cleanup deletion, compound
scoring, candidate nomination, M2.8 final export/reporting, or Milestone 3.

## Next Step

Next milestone:

```text
M2.8 / M2 Prompt 9 — accepted pocket export and Milestone 2 report
```
