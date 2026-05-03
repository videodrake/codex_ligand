# M3 ext20k mini publication interpretation notes

이 문서는 `m3_ext20k_clean3_20260503_023548` 결과를 논문 Methods/Results/Discussion에 반영하기 위한 해석 노트입니다. 핵심 원칙은 이 결과를 **validated inhibitor**, **proven binder**, **drug candidate**로 쓰지 않고, **computational compound-pocket hypothesis**로 표현하는 것입니다.

## Run Summary

```text
M2 source run: prod_ppi_ext20k_20260502_124707
M3 docking run: m3_ext20k_clean3_20260503_023548
Final status: phase3_compounds PASS, M3-T11 report and handoff completed
Profile: mini focused docking
```

Final M3 stage summary:

```text
M3-T6 Vina collection: PASS, expected_jobs=36, completed=36, parsed_pose_rows=186
M3-T7 pose attribution: PASS, poses_classified=186, hard_gate_pass=134, rejected=52
M3-T8 pose clustering: PASS, eligible=134, clusters=78, converged=6
M3-T9 anchor convergence: PASS, support_rows=7, anchor_rows=9, evidence_ready=True
M3-T10 evidence tiering: PASS, evidence_rows=11, candidate_rows=11, tier1=2, tier2=0
M3-T11 report cleanup: PASS, reports=11, cleanup_mode=execute, planned=0, deleted=0
```

## Overall Study Logic

이 workflow는 EGFR dimer 상태에서 MYO1D 접촉 가능성이 높은 PPI/membrane-proximal pocket을 찾고, 공개 compound identifier 3개(`Cpd-A`, `Cpd-B`, `Cpd-C`)가 해당 pocket들에서 반복 docking 및 geometry-based filtering을 통해 일관된 computational support를 보이는지 평가한 것입니다.

분석 흐름은 다음과 같습니다.

```text
M2 PPI ensemble
  -> pocket discovery
  -> ATP/PPI/membrane/dimer hard gating
  -> accepted pocket boxes for M3
  -> focused Vina mini docking
  -> pose attribution
  -> pose clustering
  -> anchor convergence
  -> evidence tiering
  -> report and handoff
```

## M2 Pocket Discovery Interpretation

M2에서는 EGFR-MYO1D PPI model ensemble에서 pocket 후보를 찾고, 다음 hard gate를 적용했습니다.

- ATP site overlap exclusion
- PPI relevance
- membrane lower-lateral geometry
- dimer accessibility
- mapping/origin traceability

M2 결과:

```text
pocket families evaluated: 24
accepted primary pocket: 1
accepted secondary cryptic pockets: 3
exported for M3: 4 pockets
```

대표 primary pocket:

```text
pocket_family_id: m2_6_pocket_family_0015
states: EGFR_160-185, EGFR_170-200
residue union: 752;754;755;758;759;761;762;860;861;862;863
classification: non-ATP, PPI-relevant, lower-lateral accessible, dimer-accessible
```

논문식 표현 예시:

> M2 pocket gating identified one recurrent primary pocket and three secondary cryptic pocket families that passed non-ATP, PPI-relevance, membrane-geometry, and dimer-accessibility gates. The primary pocket family contained EGFR residues 752, 754, 755, 758, 759, 761, 762, and 860-863.

## M3 Docking Setup

M3에서는 M2에서 export된 accepted pocket boxes를 대상으로 mini focused docking을 수행했습니다.

```text
compounds: Cpd-A, Cpd-B, Cpd-C
accepted docking boxes: 4
receptor states: EGFR_160-185, EGFR_170-200
planned jobs: 36
completed jobs: 36
parsed pose rows: 186
```

36개 job 구성:

```text
3 compounds x 4 pocket boxes x 3 repeats = 36
```

논문식 표현 예시:

> The M3 mini focused docking stage evaluated three public compound identifiers against four M2-derived EGFR pocket boxes using three repeats per compound-pocket combination, producing 36 completed docking jobs and 186 parsed pose rows.

## Pose Attribution

Docking pose는 Vina affinity만으로 우선순위화하지 않았습니다. 각 pose는 M2-derived pocket framework에 다시 매핑되었고, 다음 hard gate로 분류되었습니다.

- pocket retention
- ATP-site migration absence
- membrane penetration absence
- dimer-interface clash absence
- PPI-proximal relationship confirmation
- receptor/pocket traceability

결과:

```text
poses classified: 186
hard-gate pass: 134
rejected: 52
```

논문식 표현 예시:

> Docked poses were not prioritized by affinity alone. Instead, poses were attributed to the M2-derived pocket framework and filtered using geometric hard gates, including pocket retention, ATP-site migration, membrane penetration, dimer-interface clash, and PPI-proximal relationship checks.

## Pose Clustering

Hard-gate pass pose 134개를 clustering하여 반복 docking에서 유사 pose가 재현되는지 확인했습니다.

```text
eligible pose rows: 134
clusters written: 78
converged clusters: 6
```

Mini relaxed clustering configuration:

```text
centroid_distance_cutoff_A: 2.0
heavy_atom_rmsd_cutoff_A: 2.0
cluster_fraction_min: 0.20
cluster_size_min: 2
pocket_retention_fraction_min: 0.95
atp_migration_fraction_max: 0.0
membrane_penetration_fraction_max: 0.0
dimer_interface_clash_fraction_max: 0.0
```

해석:

- 78개 cluster 중 6개가 convergence 기준을 통과했습니다.
- 많은 pose는 scattered/singleton으로 남았고, 일부 compound-pocket 조합만 반복적인 support를 보였습니다.
- 이 단계는 affinity ranking이 아니라 pose recurrence와 geometry consistency를 보는 단계입니다.

## Anchor Convergence

Pose clustering 결과를 compound-pocket-state 단위로 통합하여 anchor convergence를 평가했습니다.

```text
compound_pocket_support_rows: 7
anchor_rows: 9
anchors_allowed_for_evidence_integration: 3
evidence_ready: True
```

해석:

- 3개 anchor가 evidence integration으로 넘어갈 만큼 충분한 computational support를 보였습니다.
- multi-compound primary-state support가 관찰되었습니다.
- reference-only support는 최종 promotion에 사용하지 않았습니다.

논문식 표현 예시:

> Anchor convergence analysis aggregated converged pose clusters into compound-pocket support rows and anchor-level evidence. Seven compound-pocket support rows and nine anchor rows were generated, with three anchors allowed for downstream evidence integration.

## Evidence Tiering

최종 evidence integration에서는 M2 pocket gate, M3 pose attribution, clustering, anchor convergence를 통합해 candidate hypothesis tier를 부여했습니다.

```text
evidence rows: 11
candidate rows: 11
Tier 1: 2
Tier 2: 0
report_ready: True
```

중요한 원칙:

- Vina affinity는 descriptive metadata로 기록되었습니다.
- Vina affinity는 Tier assignment, candidate promotion, best-compound selection에 사용하지 않았습니다.
- ATP-confounded evidence, reference-only evidence, broad-scan-only evidence는 Tier 1로 승격하지 않았습니다.
- Tier 1은 experimental validation이 아니라 가장 강한 computational hypothesis입니다.

논문식 표현 예시:

> Final evidence tiering produced 11 compound-pocket hypothesis rows, including two Tier 1 computational hypotheses. Tier assignment was based on hard-gate compliance and convergence evidence rather than Vina affinity-based ranking.

## Tier Interpretation

```text
Tier 1
  M2 pocket gate, M3 pose attribution, clustering, anchor convergence, and evidence tiering을 모두 통과한 최상위 computational hypothesis.

Tier 2/3
  일부 computational support는 있으나 Tier 1보다 state, compound, cluster, or convergence support가 약한 가설.

Reject
  ATP migration, membrane/dimer conflict, missing mapping, pocket-retention failure, or other hard-gate/quarantine reason으로 제외된 항목.
```

권장 표현:

```text
computational hypothesis
prioritized compound-pocket hypothesis
docking-supported hypothesis
PPI-proximal pocket engagement hypothesis
candidate hypothesis for follow-up validation
```

피해야 할 표현:

```text
validated inhibitor
proven binder
drug candidate
clinically relevant compound
confirmed EGFR-MYO1D PPI inhibitor
```

## Suggested Results Paragraph

> The M3 mini focused docking workflow completed successfully using four M2-derived EGFR pocket boxes and three public compound identifiers. All 36 planned mini docking jobs completed, yielding 186 parsed docking pose rows. Pose attribution classified all poses against the M2 pocket framework, with 134 poses passing hard geometric filters and 52 poses rejected. Subsequent clustering assigned the 134 eligible poses into 78 clusters, of which six met the relaxed mini convergence criteria. Anchor convergence analysis produced seven compound-pocket support rows and nine anchor rows, with three anchors allowed for downstream evidence integration. Final evidence tiering generated 11 compound-pocket hypothesis rows, including two Tier 1 computational hypotheses. Importantly, Vina affinity values were retained as descriptive metadata only and were not used for candidate promotion or best-compound selection.

## Suggested Methods Paragraph

> M3 focused docking was performed on M2-derived non-ATP, PPI-relevant EGFR pocket families. Receptor boxes were prepared from accepted M2 pocket exports, and three public compound identifiers were docked against four accepted pocket boxes with three repeats per compound-pocket combination. Docking outputs were parsed into raw pose records and subjected to pose attribution against the M2 pocket framework. Poses were filtered using hard geometric gates for pocket retention, ATP-site migration, membrane penetration, dimer-interface clash, and PPI-proximal relationship. Passing poses were clustered using a hybrid centroid/RMSD-based approach, followed by compound-pocket support aggregation, anchor convergence analysis, and final evidence tiering. Affinity scores were recorded but were not used for candidate promotion, ranking, or Tier 1 assignment.

## Suggested Limitations Paragraph

> These results should be interpreted as computational hypotheses rather than experimentally validated binding or inhibitory activity. The mini docking profile was designed for focused prioritization and workflow validation, not exhaustive chemical screening. The Tier 1 hypotheses require further validation through expanded docking, orthogonal structural modeling, biochemical binding assays, and EGFR-MYO1D PPI disruption assays.

## Suggested Figure/Table Plan

Recommended figures:

1. Workflow schematic: M2 PPI ensemble -> pocket gating -> M3 docking -> pose attribution -> clustering -> anchor convergence -> evidence tiering.
2. M2 accepted pocket map: highlight `m2_6_pocket_family_0015` and secondary cryptic pockets.
3. M3 funnel plot: 36 jobs -> 186 poses -> 134 hard-gate pass poses -> 78 clusters -> 6 converged clusters -> 3 evidence-ready anchors -> 2 Tier 1 hypotheses.
4. Evidence heatmap: compounds vs pocket families with tier/evidence support.

Recommended tables:

1. Accepted pocket families from M2.
2. Docking and pose attribution summary.
3. Cluster convergence summary.
4. Anchor convergence summary.
5. Final candidate hypothesis tier table.

## Files to Cite Internally

```text
fresh/runs/m3_ext20k_clean3_20260503_023548/report/milestone3_summary.md
fresh/runs/m3_ext20k_clean3_20260503_023548/report/final_candidate_summary.md
fresh/runs/m3_ext20k_clean3_20260503_023548/report/reviewer_risk_notes.md
fresh/runs/m3_ext20k_clean3_20260503_023548/phase3_compounds/tables/final_m3_candidate_hypotheses.csv
fresh/runs/m3_ext20k_clean3_20260503_023548/phase3_compounds/tables/pocket_compound_evidence_table.csv
fresh/runs/m3_ext20k_clean3_20260503_023548/phase3_compounds/tables/compound_anchor_convergence.csv
fresh/runs/m3_ext20k_clean3_20260503_023548/phase3_compounds/qc/final_candidate_gate_qc.json
fresh/runs/m3_ext20k_clean3_20260503_023548/phase3_compounds/qc/m3_report_qc.json
```

## Minimal Final Summary

```text
M3 mini focused docking pipeline completed successfully.
Final status: phase3_compounds PASS, M3-T11 report and handoff completed.
Run ID: m3_ext20k_clean3_20260503_023548
M2 source: prod_ppi_ext20k_20260502_124707
Evidence tiering: PASS, evidence_rows=11, candidate_rows=11, tier1=2
Report cleanup: PASS, reports=11, cleanup_mode=execute, planned=0, deleted=0
```
