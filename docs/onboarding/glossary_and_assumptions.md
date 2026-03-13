# Glossary And Assumptions

This document defines the working vocabulary for the EGFR-MYO1D repository as it exists today. It is not a generic structural-biology glossary. It is a project-specific interpretation guide for new GPTs so that the same term is not read differently across old planning docs, current outputs, and code.

## Reading Rule

- If a term has both a project-specific meaning and a generic scientific meaning, prefer the project-specific meaning in this document.
- If an older document uses the same word differently, trust current code, current outputs, and the current onboarding docs first.

## 1. Core Project Terms

| Term | Meaning in this repository | Current note |
|------|------|------|
| `EGFR-MYO1D pipeline` | A state-comparison research pipeline combining a Vina-centered ligand layer with a Phase 1 receptor-side PPI evidence layer | Current routine baseline is Vina-centered, not fully PPI-first end to end |
| `receptor state` | One of the three fixed EGFR structural states used for comparison | Exactly `3GT8_raw`, `3GT8_cl38_48`, `3GT8_cl85_100` |
| `state comparison` | Comparing evidence across those three receptor states rather than treating the receptor as a single static object | Applies to Vina pockets and Phase 1 PPI residues |
| `baseline` | The current default operational interpretation of the repo | PyRosetta primary Phase 1, LightDock active secondary, AFM legacy optional, `max_workers = 16` |

## 2. Fixed Biological Entities

| Term | Meaning in this repository | Current note |
|------|------|------|
| `3GT8_raw` | Raw 3GT8-derived receptor state | Current active receptor state |
| `3GT8_cl38_48` | MD cluster representative from the 38-48 ns window | Current active receptor state |
| `3GT8_cl85_100` | MD cluster representative from the 85-100 ns window | Current active receptor state |
| `extended beta-meander` | Current primary Phase 1 MYO1D partner construct spanning residues `955-1006` | Replaces older truncated `962-1006` framing |
| `TH1 domain` | Larger MYO1D source/domain context retained in repo | Current Phase 1 docs treat it as plausibility context, not the primary search input |
| `AFM` | AlphaFold-Multimer-derived interface support path | Legacy optional only unless explicitly re-enabled |

## 3. Evidence Vocabulary

| Term | Meaning in this repository | Current note |
|------|------|------|
| `primary evidence` | The method or artifact type that should be trusted first for the relevant question | For ligand questions: Vina outputs. For Phase 1 PPI patch questions: PyRosetta |
| `secondary evidence` | Independent supporting evidence that can strengthen or weaken confidence but should not replace the primary method by default | LightDock is the active Phase 1 secondary validation path |
| `legacy evidence` | Older or optional evidence retained for reference but not part of the default baseline | AFM parser outputs and older fragment-based pilot results |
| `method agreement` | Whether two independent evidence paths support the same residue/site | Common values include `both`, `pyrosetta_only`, `lightdock_only`, `none` |
| `confidence` | Repository-specific downstream confidence label, not absolute biological truth | Often `high`, `medium`, or `low` in Phase 1 and Phase 4 outputs |

Practical interpretation:

- `primary` does not mean infallible.
- `secondary` means supportive and independent, not disposable.
- `legacy` means preserved but not baseline.

## 4. Structural And Output Terms

| Term | Meaning in this repository | Where it appears |
|------|------|------|
| `pose` | One docked ligand placement from a Vina run | `vina_pose_table.csv` |
| `pocket` | A grouped ligandable region inferred from multiple poses or structure-based pocket proposals | `vina_pocket_table.csv`, `candidate_pockets.csv` |
| `patch` | The receptor-side PPI interface region inferred from Phase 1 structural evidence | `ppi_interface_patch_table.csv`, `phase1_downstream_patch_reference.csv` |
| `hotspot residue` | A receptor residue repeatedly supported within a patch or cluster summary | `ppi_hotspot_residues.csv`, `phase1_downstream_patch_reference.csv` |
| `consensus` | A summary built from multiple models, clusters, poses, or states rather than a single structure | Used in cluster consensus, cross-state robustness, and consensus site outputs |
| `robustness` | Persistence of a residue or site across receptor states, not just within one run | `ppi_patch_state_robustness.csv`, `candidate_pocket_state_classes.csv` |
| `handoff file` | A structured artifact meant to be consumed by the next phase | `phase1_downstream_patch_reference.csv`, `phase3_candidate_pocket_reference.csv`, `phase4_docking_evidence_reference.csv` |

Important distinctions:

- A `patch` is a Phase 1 PPI concept.
- A `pocket` is a ligandability/docking concept.
- A `pose` is a single docking hypothesis, not a site-level conclusion.

## 5. Phase 1 Terms

| Term | Meaning in this repository | Current note |
|------|------|------|
| `full_kinase_domain` | Current structured Phase 1 receptor construct, covering N-lobe and C-lobe | Current Phase 1 baseline construct type |
| `legacy_clobe_fragment` | Older fragment-based receptor context used in pilot-style PPI work | Historical reference only |
| `orientation-aware filtering` | Mandatory Phase 1 quality control step intended to reject face-flipped or biologically implausible partner orientations | Part of the current Phase 1 design and output chain |
| `cross-method convergence` | Residue-level agreement between PyRosetta and LightDock | `cross_method_convergence.csv` |
| `Phase 1 downstream patch reference` | The structured patch handoff exported for later phases | `phase1_downstream_patch_reference.csv` |

## 6. Phase 2 Relationship Classes

These classes describe how a candidate pocket relates to the Phase 1 patch.

| Term | Meaning in this repository | Current note |
|------|------|------|
| `orthosteric_candidate` | Pocket directly overlaps Phase 1 hotspot residues | Highest direct PPI relevance class in Phase 2 relationship output |
| `rim_candidate` | Pocket partially overlaps the patch or borders it closely | Near-interface modulator candidate |
| `allosteric_candidate` | Pocket is spatially near the patch but does not directly overlap patch hotspots | Mechanistically indirect candidate |
| `low_relevance_candidate` | Pocket is structurally distant or not meaningfully connected to the patch | Ligandable maybe, but not currently PPI-relevant |

Relationship-class caution:

- These are Phase 2 structural relation labels.
- They are not the same thing as the final Phase 4 mechanistic class names.

## 7. Phase 4 Mechanistic Classes

These classes are the advanced final-interpretation labels used in the Phase 4 stack.

| Term | Meaning in this repository | Current note |
|------|------|------|
| `orthosteric_disruptor_candidate` | Site directly overlaps the MYO1D attachment patch and is interpreted as a direct disruption candidate | Phase 4 final mechanistic class |
| `interface_rim_modulator_candidate` | Site sits at the interface rim and may modulate attachment indirectly from the boundary | Phase 4 final mechanistic class |
| `allosteric_modulator_candidate` | Site is spatially distinct but interpreted as a plausible indirect perturbation route | Phase 4 final mechanistic class |
| `ligandable_but_ppi_irrelevant_candidate` | Site may be druggable but lacks a mechanistic link to MYO1D attachment disruption | Phase 4 final mechanistic class |
| `uncertain_mechanism_candidate` | Evidence is incomplete or contradictory, so a stronger mechanism label is withheld | Honest low-certainty fallback class |

## 8. Druggability And State Terms

| Term | Meaning in this repository | Current note |
|------|------|------|
| `druggability_confidence` | Confidence that a pocket is chemically tractable based on proposal evidence | Usually `high`, `medium`, or `low` |
| `overall_druggability_tier` | Tiered summary of druggability plus relevance support | Common values: `tier_1`, `tier_2`, `tier_3` |
| `state_robust` | Pocket pattern matches across all states with comparable support | Phase 2 cross-state class |
| `state_shifted` | Pocket appears related across states but with shifted placement or matching behavior | Phase 2 cross-state class |
| `state_specific_pocket` | Pocket currently appears in only one state or only one state has usable matching data | Common current value in Phase 2/4 outputs |
| `uncertain` | Cross-state evidence is insufficient to assign a stronger state class | Conservative fallback |

Important distinction:

- `robustness` in Phase 1 usually refers to residue persistence across receptor states.
- `state_class` in Phase 2/4 usually refers to pocket-level cross-state behavior.

## 9. Numbering And Chain-Mapping Assumptions

These are working assumptions that a new GPT should use unless a user asks for a different numbering frame explicitly.

### Numbering assumptions

| Topic | Current assumption |
|------|------|
| Main receptor numbering | Use PDB-consistent numbering for current structured Phase 1 outputs |
| Current structured Phase 1 receptor range | `699-1007` |
| N-lobe/C-lobe boundary | Residue `838` in current PDB-consistent numbering |
| UniProt relation | Current docs mention approximate UniProt = PDB + 24 mapping, but current outputs are primarily expressed in PDB-consistent numbering |
| Partner numbering | Extended beta-meander is treated as `955-1006` in current Phase 1 metadata |

### Chain assumptions

| Topic | Current assumption |
|------|------|
| Structured Phase 1 receptor chain | Receptor is normalized to chain `A` |
| Structured Phase 1 partner chain | Partner is normalized to chain `B` in docking-pair metadata |
| Legacy prepared dimer assets | Older dimer-prepared files may include remapped chain/residue spaces and offset handling |
| Legacy offset pattern | Older prepared dimer assets may use a chain-B renumbering offset of `+1000`, producing values like `1701-2007` before restoration |

Practical rule:

- For current structured Phase 1 interpretation, assume chain `A` is the receptor and use the current normalized numbering.
- If you are reading older prepared dimer outputs, check mapping CSV or restoration notes before comparing residue IDs directly.

## 10. `orientation_validation_status` Meaning

This field is easy to misread, so use the meanings below.

| Value | Meaning in this repository | Current note |
|------|------|------|
| `orientation_validated` | The relevant Phase 1 evidence survived the orientation-aware filter and can be treated as orientation-checked | Current desired positive status for structured Phase 1 PyRosetta-derived patch evidence |
| `not_available` | No equivalent orientation validation signal is available for that artifact | Common honest placeholder for LightDock raw support tables and some legacy-linked paths |
| `mixed_orientation_status` | Aggregated view contains mixed upstream orientation statuses | Summary/report-only mixed label, not a direct per-model judgment |

Interpretation rule:

- `not_available` does not mean "failed".
- It means the orientation check is absent or not directly propagated for that artifact.

## 11. Output-Reading Assumptions

| Term | Meaning in this repository | Current note |
|------|------|------|
| `routine baseline outputs` | The outputs most directly tied to the default Vina-centered operational flow | `valid_sites.csv`, `cross_method_agreement.csv`, `project_report.txt`, Vina tables |
| `phase-separated outputs` | Outputs tied to the newer Phase 1-4 scientific workflow modules | `output/phase1_ppi/` through `output/phase4_perturbation/` |
| `pointer stub file` | A small file that points to another payload path instead of containing the actual dataset body | Seen under `output/egfr_myo1d_vina/` root |

## 12. Safe Default Assumptions For New GPTs

Start with these assumptions unless the user explicitly redirects you.

- The repository currently has one active receptor ensemble of exactly three states.
- Vina is still the current center of gravity for routine ligand evidence.
- PyRosetta is the primary Phase 1 structural evidence layer.
- LightDock is the active secondary Phase 1 validation layer.
- AFM is legacy optional only.
- Current structured Phase 1 outputs prefer PDB-consistent numbering with receptor chain `A`.
- `orientation_validation_status = not_available` is an honesty marker, not an automatic failure label.
