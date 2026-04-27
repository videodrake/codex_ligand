# EGFR-MYO1D membrane-aware PPI-to-pocket-to-compound workflow

## 1. Final objective

The final scientific objective is not simply to dock MYO1D to EGFR. The objective is:

> Identify the EGFR-side PPI surface used by MYO1D, find membrane-compatible druggable pockets near that surface in the EGFR dimer context, and nominate compounds that can occupy those pockets to perturb the EGFR-MYO1D interaction.

This changes the workflow logic. The pipeline must be PPI-first and dimer-first:

1. Build biologically valid EGFR kinase-domain dimers.
2. Dock MYO1D only against those dimers.
3. Retain only poses compatible with the MYO1D active face and the plasma-membrane geometry.
4. Derive an EGFR-side consensus PPI patch from repeated seeds/states.
5. Search for pockets on dimeric EGFR, not on isolated monomer EGFR.
6. Keep only pockets that are near the PPI patch, non-ATP, druggable, and accessible from the lower/lateral membrane-proximal side.
7. Use focused docking to nominate compounds for the surviving pockets.

## 2. Key biological constraints

### EGFR must be treated as a dimer

EGFR kinase activation is structurally coupled to asymmetric kinase-domain dimerization and juxtamembrane/TM geometry. Therefore, monomer-only pocket discovery is not sufficient for the main paper claim. Monomer results can be used only as technical controls or historical comparison.

Main-run receptor rule:

- Allowed: EGFR kinase-domain dimer.
- Disallowed for main conclusion: isolated monomer EGFR.
- Required metadata: protomer identity, dimer source, chain mapping, residue offset mapping, dimer-interface geometry.

### The MYO1D-binding site must be membrane-compatible

MYO1D is a membrane-associated myosin, and the reported EGFR-MYO1D biology describes EGFR family kinase-domain anchoring at the plasma membrane. Therefore, the acceptable EGFR-side PPI patch should not be an arbitrary solvent-exposed kinase surface. It should be compatible with a membrane-proximal cytoplasmic encounter.

Operational interpretation:

- The valid PPI patch should lie on the lower/lateral side of the EGFR dimer.
- "Lower" means closer to the juxtamembrane/TM side of EGFR after a membrane reference frame is assigned.
- "Lateral" means outward-facing on the side of the dimer, not buried in the central protomer-protomer interface and not on a top/distal face that MYO1D is unlikely to approach from the membrane.
- The pocket must be accessible to a cytosolic small molecule and must not require penetration through the membrane, the dimer interface, or the kinase ATP cleft.

## 3. Required coordinate frame

The current kinase-domain-only inputs do not contain a real membrane. We must therefore construct a reproducible virtual membrane frame before classifying "side" and "below".

Recommended references:

- EGFR juxtamembrane-kinase asymmetric dimer reference: PDB `3GOP`.
- EGFR TM-JM membrane reference: PDB `2M20`.
- Existing inactive kinase dimer state: PDB `3GT8` / local `3GT8_raw`.

Frame construction:

1. Align each EGFR kinase protomer to a reference EGFR JM-kinase structure.
2. Transfer the approximate JM/TM anchor direction from the reference.
3. Define a membrane normal vector `membrane_normal`.
4. Define the kinase-side cytoplasmic coordinate as the direction away from the membrane plane.
5. Define a dimer axis from protomer-1 centroid to protomer-2 centroid.
6. Define lower/lateral pocket coordinates by projecting pocket centroids onto:
   - the membrane-normal axis,
   - the dimer axis,
   - and the outward radial vector from the dimer center.

Required output:

```text
output/phase0_geometry/membrane_frame.json
output/phase0_geometry/dimer_geometry_qc.csv
```

Minimum `membrane_frame.json` fields:

```json
{
  "receptor_id": "3GT8_raw_dimer",
  "reference_jm_kinase_pdb": "3GOP",
  "reference_tm_jm_pdb": "2M20",
  "membrane_normal": [0.0, 0.0, 1.0],
  "membrane_plane_point": [0.0, 0.0, 0.0],
  "dimer_axis": [1.0, 0.0, 0.0],
  "protomer_a_centroid": [0.0, 0.0, 0.0],
  "protomer_b_centroid": [0.0, 0.0, 0.0],
  "coordinate_convention": "positive membrane_normal points from cytosol toward membrane"
}
```

The numeric values above are placeholders. The actual values must be computed and stored per receptor state.

## 4. MYO1D construct and terminal-noise control

MYO1D beta sheets 8, 9, and 12 must be included. The preferred decision is:

- Candidate main construct: `955-1001`.
- Conservative comparator: `955-1006` with tail masking.

Binding evidence zones:

- Sheet 8/9 active face: `961-964,968-972`.
- Sheet 12 support zone: `993-997`.
- N-terminal buffer: `955-960`.
- C-terminal cap/noise zone: `998-1001` for `955-1001`, or `998-1006` for `955-1006`.

Pipeline annotation:

```ini
[Constraints]
key_residues_B = 961-964,968-972,993-997
key_residue_bonus_weight = 0.0

[ExperimentalData]
critical_residues_B = 961-964,968-972,993-997
non_binding_residues_B = 998-1006
```

`key_residue_bonus_weight = 0.0` is important. The first rerun should record key-residue contact ratios without changing pose ranking. Any positive bonus weight changes model selection and needs explicit approval before a paper-grade production run.

## 5. PPI docking workflow

### Phase 1A. Dimer receptor preparation

For every receptor state:

- Prepare EGFR dimer receptor.
- Preserve protomer mapping:
  - protomer A: original numbering, e.g. `699-1007`.
  - protomer B: offset numbering, e.g. `1699-2007`, or equivalent mapping CSV.
- Record whether the dimer is experimental, template-derived, or MD-derived.

Important decision:

- If MD-derived states are monomer-only, they must be embedded into a common dimer geometry or replaced with true dimer state representatives.
- A 3-state dimer conclusion is not defensible unless every state has a documented dimer construction method.

### Phase 1B. MYO1D PPI docking

Run PyRosetta docking against dimer receptors only.

Primary filters:

- MYO1D active-face orientation pass.
- MYO1D beta sheet 8/9 contact present.
- Sheet 12 support tracked separately.
- Tail-dominant contacts rejected or quarantined.
- ATP-site overlap excluded.

New membrane-aware PPI filter:

- MYO1D centroid must be on the cytosolic side of the virtual membrane plane.
- MYO1D must not penetrate the virtual membrane plane.
- The EGFR contact centroid must be lower/lateral relative to the dimer frame.
- Contacts buried inside the central protomer-protomer interface are rejected unless explicitly classified as biologically accessible.

Required output:

```text
output/phase1_ppi/ppi_pose_membrane_qc.csv
output/phase1_ppi/ppi_consensus_patch_dimer_aware.csv
```

Suggested columns for `ppi_pose_membrane_qc.csv`:

- `receptor_id`
- `seed`
- `pose_id`
- `orientation_class`
- `myo1d_contact_class`
- `tail_contact_fraction`
- `egfr_contact_centroid_x/y/z`
- `membrane_z_percentile`
- `lower_side_class`
- `lateral_side_class`
- `dimer_interface_overlap_fraction`
- `membrane_compatible`
- `reject_reason`

## 6. Pocket discovery workflow

Pocket discovery must run on the dimer receptor used for PPI docking.

Required tools:

- `fpocket`: geometric pocket detection.
- `P2Rank`: machine-learning pocket proposal.
- Optional for publication strength: MDpocket or trajectory-based pocket occupancy if trajectories are available.

Hard gates before compound docking:

| Gate | Rule | Reason |
|---|---|---|
| G0 dimer gate | pocket generated on EGFR dimer receptor | monomer pockets may be artificial |
| G1 ATP exclusion | ATP-site pockets cannot be final PPI-modulator candidates | avoids rediscovering kinase inhibitors |
| G2 PPI proximity | pocket is orthosteric/rim/allosteric relative to consensus PPI patch | target must perturb EGFR-MYO1D |
| G3 membrane-side gate | pocket is lower/lateral in the dimer membrane frame | matches membrane-associated biology |
| G4 dimer-accessibility gate | pocket is not buried in the central dimer interface | compound must be physically accessible |
| G5 state support | pocket appears across states or is a reproducible cryptic pocket | avoids single-frame artifacts |

Suggested new output:

```text
output/phase2_pockets/pocket_membrane_geometry.csv
output/phase2_pockets/pocket_dimer_accessibility.csv
output/phase2_pockets/phase3_candidate_pocket_reference_membrane_aware.csv
```

Suggested columns for `pocket_membrane_geometry.csv`:

- `receptor_id`
- `pocket_id`
- `protomer_assignment`
- `centroid_x/y/z`
- `membrane_z_A`
- `membrane_z_percentile`
- `lower_gate_pass`
- `lateral_gate_pass`
- `outward_radial_score`
- `dimer_interface_distance_A`
- `dimer_interface_overlap_fraction`
- `membrane_geometry_class`
- `membrane_gate_reason`

Pocket geometry classes:

- `lower_lateral_accessible`: preferred class.
- `lower_central_interface`: biologically interesting but low compound accessibility.
- `upper_or_distal`: reject for main PPI-modulator claim.
- `membrane_occluded`: reject.
- `ambiguous_geometry`: report but do not prioritize.

## 7. Pocket ranking logic

Do not let Vina affinity alone dominate the result. The ATP cleft and other deep pockets will often score well but may be irrelevant to EGFR-MYO1D disruption.

Recommended priority order:

1. Pass all hard gates.
2. Strong PPI relationship: rim or allosteric within the lower/lateral PPI neighborhood.
3. Good druggability: volume, enclosure, hydrophobic/hydrogen-bond balance, multi-tool support.
4. State robustness or reproducible cryptic opening.
5. Focused docking support from chemically diverse compounds.
6. Pose mechanism: compound points toward or stabilizes the MYO1D-disruptive surface.

Proposed scoring axes, pending approval before implementation:

| Axis | Meaning |
|---|---|
| A1 PPI relevance | overlap/distance to dimer-aware EGFR-MYO1D patch |
| A2 membrane-dimer compatibility | lower/lateral/accessibility geometry |
| A3 druggability | fpocket/P2Rank/volume/shape/chemistry |
| A4 compound docking support | Vina affinity, convergence, ligand diversity, pose consistency |
| A5 robustness | state recurrence, seed recurrence, method agreement |

These axes should be documented even if the current code keeps the old weights. Changing numeric weights in production requires explicit approval.

## 8. Compound derivation strategy

The existing three ligands are useful as probe compounds for pocket characterization, but they are not sufficient to claim compound discovery.

Recommended compound funnel:

### Stage C1. Probe docking

- Use the existing three diverse ligands only to test whether a pocket can accept small molecules.
- Output: pocket-level ligandability, not final compounds.

### Stage C2. Fragment/enriched-library docking

Use a real compound library:

- fragment-like molecules for shallow rim pockets,
- PPI-oriented libraries for protein-surface cavities,
- filtered by PAINS/reactive groups,
- diverse by scaffold and physicochemical properties.

Dock only to pockets that passed the dimer/membrane/PPI gates.

### Stage C3. Ensemble focused docking

For each surviving pocket:

- Dock across all receptor states where the pocket exists.
- Use repeated seeds.
- Keep pose clusters, not single best poses.
- Prefer compounds that retain the same interaction mode across states.

### Stage C4. Rescoring and mechanism filter

A compound is a candidate only if it satisfies:

- binds in the lower/lateral dimer-compatible pocket,
- avoids the ATP site,
- contacts residues geometrically coupled to the MYO1D patch,
- has stable pose convergence,
- does not rely on clashes with the receptor or impossible membrane penetration,
- ideally points substituents toward the MYO1D approach path.

### Stage C5. MD/relaxation validation

For the top compounds:

- local minimization or induced-fit relaxation,
- short MD stability test,
- contact persistence analysis,
- pocket hydration/occlusion check,
- optional MM/GBSA or equivalent rescoring.

Final deliverable:

```text
output/final_compounds/membrane_aware_compound_shortlist.csv
```

Suggested columns:

- `compound_id`
- `smiles`
- `source_library`
- `target_pocket_id`
- `receptor_state_support`
- `best_affinity`
- `median_affinity`
- `pose_cluster_occupancy`
- `ppi_disruption_residue_contacts`
- `membrane_geometry_pass`
- `atp_site_flag`
- `admet_flags`
- `selection_tier`
- `mechanistic_rationale`

## 9. What changes relative to the current workflow

Current advanced pipeline:

```text
PPI docking -> PPI patch -> fpocket/P2Rank -> PPI relationship -> focused Vina -> perturbation score
```

Required membrane-aware pipeline:

```text
dimer receptor + membrane frame
    -> MYO1D active-face PPI docking
    -> membrane-compatible dimer-side PPI patch
    -> dimer pocket proposal
    -> PPI + membrane + dimer accessibility gates
    -> focused docking against surviving pockets
    -> compound shortlist
```

The existing `PKT07` and `PKT34` style conclusions should be treated as hypotheses until rerun through the dimer/membrane-aware gates. They should not be final manuscript claims unless:

- the pockets are reproduced on dimer receptors,
- they lie on the lower/lateral membrane-compatible side,
- they remain close to the dimer-aware MYO1D patch,
- and focused compound docking remains stable.

## 10. Implementation targets in this repository

Recommended new or modified components:

- New: `egfr_pipeline/phase0/membrane_frame.py`
  - builds `membrane_frame.json`.
- New: `egfr_pipeline/phase1/membrane_pose_filter.py`
  - classifies MYO1D poses by membrane compatibility.
- Modify: `egfr_pipeline/phase2/patch_relationship.py`
  - use dimer-aware patch residue IDs and protomer mapping.
- New: `egfr_pipeline/phase2/membrane_geometry.py`
  - classifies each pocket as lower/lateral/central/upper.
- Modify: `egfr_pipeline/phase2/phase3_export.py`
  - export only pockets passing dimer, PPI, ATP, and membrane gates as primary/secondary candidates.
- Modify: `egfr_pipeline/phase3/run_diverse_docking.py`
  - add support for approved compound libraries beyond the three probe ligands.
- New: `egfr_pipeline/final_compounds/compound_shortlist.py`
  - integrates docking, membrane geometry, PPI relevance, and compound filters.

Important compatibility note:

- Do not edit `egfr_pipeline/paths.py` unless absolutely necessary.
- Do not change scoring weights without explicit approval.
- Preserve existing CSV columns; add new columns rather than renaming old ones.

## 11. Paper-grade decision criteria

A final pocket can be claimed as an EGFR-MYO1D modulator pocket only if it satisfies all of the following:

1. It is found on EGFR dimer receptor input.
2. It is non-ATP.
3. It is linked to the orientation-filtered MYO1D PPI patch.
4. It is lower/lateral and membrane-compatible.
5. It is compound-accessible in the dimer geometry.
6. It has reproducible pocket evidence across states or a documented cryptic-state mechanism.
7. It supports focused docking by multiple chemically diverse compounds or fragments.
8. Its proposed mechanism is PPI perturbation, not generic kinase inhibition.

A final compound can be nominated only if it satisfies all of the following:

1. It docks reproducibly into a final accepted pocket.
2. It avoids ATP-site binding.
3. Its pose is compatible with membrane/dimer access.
4. It interacts with residues coupled to the MYO1D patch.
5. It is chemically plausible after PAINS/reactivity and basic property filtering.
6. It remains stable after local relaxation or short MD validation.

## 12. References used for design rationale

- Ko et al., 2019, MYO1D binds the EGFR family kinase domain and anchors them to the plasma membrane: https://www.nature.com/articles/s41388-019-0954-8
- Ko et al., 2021 correction: https://www.nature.com/articles/s41388-021-01675-y
- RCSB `3GOP`, EGFR juxtamembrane and kinase domains: https://www.rcsb.org/structure/3GOP
- RCSB `2M20`, EGFR TM-JM segment in bicelles: https://www.rcsb.org/structure/2m20
- RCSB `3GT8`, inactive EGFR kinase-domain structure: https://www.rcsb.org/structure/3GT8
- Jura et al., 2009, juxtamembrane activation mechanism: https://doi.org/10.1016/j.cell.2009.04.025
