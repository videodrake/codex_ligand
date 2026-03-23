# AC Coverage Checklist — Final Verification

> PRD Must-Have (F-1~F-4) + Should-Have (F-5~F-6) + Success Criteria (SC-1~SC-6)

## Must-Have Features (F-1 ~ F-4)

### F-1: 실험 데이터 통합

| AC | 설명 | 구현 파일 | 테스트 | 상태 |
|---|---|---|---|---|
| AC-1.1 | 리간드 3종 Tanimoto < 0.4 다양성 | scripts/assess_ligand_diversity.py | test_experimental_data_group1 | PASS |
| AC-1.2 | ATP site 37잔기 region_definitions.py 통합 | egfr_pipeline/region_definitions.py | test_setup_group0 | PASS |
| AC-1.3 | Ko et al. sheet 일관성 체크 | run_production.py | test_experimental_data_group1 | PASS |
| AC-1.4 | ATP 실험 사실 문서 반영 | docs/manual_vina.md, CLAUDE.md | test_experimental_data_group1 | PASS |
| AC-1.5 | Ko sheet 접촉 정보 pocket_summary | egfr_pipeline/vina/pocket_summary.py | test_experimental_data_group1 | PASS |

### F-2: Vina 편향 정량화

| AC | 설명 | 구현 파일 | 테스트 | 상태 |
|---|---|---|---|---|
| AC-2.1 | 포즈 영역별 분포 분석 | egfr_pipeline/vina/pose_region_classifier.py | test_vina_bias_group2 | PASS |
| AC-2.2 | Cross-receptor 잔기 비교 + validate 통합 | egfr_pipeline/validate.py | test_vina_bias_group2 | PASS |
| AC-2.3 | Bootstrap-Verdict 연동 | egfr_pipeline/verdict.py | test_vina_bias_group2 | PASS |
| AC-2.4 | C-lobe surface affinity 분석 | scripts/analyze_affinity_distribution.py | test_vina_bias_group2 | PASS |
| AC-2.5 | Vina scoring 편향 문서 | docs/manual_vina.md | test_vina_bias_group2 | PASS |

### F-3: PPI 강건성 검증

| AC | 설명 | 구현 파일 | 테스트 | 상태 |
|---|---|---|---|---|
| AC-3.1 | Threshold sweep 6단계 | egfr_pipeline/phase1/orientation_filter.py | test_ppi_group3 | PASS |
| AC-3.2 | Ambiguous model report | egfr_pipeline/phase1/orientation_filter.py | test_ppi_group3 | PASS |
| AC-3.3 | Fragment pilot 3범위 | scripts/pilot_fragment_range.py | test_ppi_group3 | PASS |
| AC-3.4 | Sheet 12 sensitivity | egfr_pipeline/phase1/orientation_filter.py | test_ppi_group3 | PASS |
| AC-3.5 | Handoff quality guards | run_production.py | test_ppi_group3 | PASS |
| AC-3.6 | Orthosteric+rim >80% WARNING | egfr_pipeline/phase2/patch_relationship.py | test_ppi_group3 | PASS |
| AC-3.7 | Conformational selection candidate | egfr_pipeline/phase1/compare_states.py | test_ppi_group3 | PASS |
| AC-3.8 | n_valid_models + 2x imbalance | egfr_pipeline/phase1/compare_states.py | test_ppi_group3 | PASS |
| AC-3.9 | Concordance strong/weak_both | egfr_pipeline/phase1/lightdock_validation.py | test_ppi_group3 | PASS |

### F-4: Verdict 메트릭 교정

| AC | 설명 | 구현 파일 | 테스트 | 상태 |
|---|---|---|---|---|
| AC-4.1 | 6개 가중치 시뮬레이션 | scripts/verdict_weight_sensitivity.py | test_verdict_group4 | PASS |
| AC-4.2 | pocket_depth_A + offset 보정 | egfr_pipeline/vina/pocket_summary.py, scripts/centroid_offset_analysis.py | test_verdict_group4 | PASS |
| AC-4.3 | Cross-receptor 차등 (20/14/0) | egfr_pipeline/verdict.py | test_verdict_group4 | PASS |
| AC-4.4 | PPI 거리 분포 분석 | scripts/verdict_weight_sensitivity.py | test_verdict_group4 | PASS |
| AC-4.5 | Phase 4 A3 축 문서 | docs/phase4_A3_axis_specification.md | test_verdict_group4 | PASS |

## Should-Have Features (F-5 ~ F-6)

### F-5: 워크플로우 비교 및 문서

| AC | 설명 | 구현 파일 | 테스트 | 상태 |
|---|---|---|---|---|
| AC-5.1 | Workflow A↔B 비교 모듈 | egfr_pipeline/workflow_comparison.py | test_workflow_group5 | PASS |
| AC-5.2 | 불일치 해석 가이드 | docs/workflow_comparison_guide.md | test_workflow_group5 | PASS |
| AC-5.3 | Allosteric candidate 태그 + report | egfr_pipeline/verdict.py, report.py | test_workflow_group5 | PASS |
| AC-5.4 | 방법론 한계 5섹션 | docs/methodology_limitations.md | test_workflow_group5 | PASS |
| AC-5.5 | 면책 조항 | verdict.py, report.py, phase4/review_report.py | test_workflow_group5 | PASS |
| AC-5.6 | Output 경로 가이드 + DEPRECATED | docs/output_path_guide.md, paths.py | test_workflow_group5 | PASS |

### F-6: 인프라 (Optional)

| AC | 설명 | 구현 파일 | 테스트 | 상태 |
|---|---|---|---|---|
| AC-6.1 | 모듈 분리 분석 | docs/module_separation_analysis.md | test_infra_group6 | PASS |
| AC-6.2 | Precheck Workflow B 도구 체크 | scripts/run_pre_qsub_checks.sh | test_infra_group6 | PASS |
| AC-6.3 | Phase 3 병렬화 확인 | (이미 구현, 문서화) | test_infra_group6 | PASS |

## Success Criteria

| SC | 설명 | 상태 | 근거 |
|---|---|---|---|
| SC-1 | Must-Have F-1~F-4 AC 전체 통과 | PASS | 위 체크리스트 23/23 |
| SC-2 | 기존 파이프라인 회귀 없이 실행 | PASS | VERDICT_FIELDS 45개 baseline 유지, import 체인 정상 |
| SC-3 | valid_sites.csv에 신규 컬럼 추가 | PASS | is_atp_site, exclusion_reason, bootstrap_confidence, contacts_sheet_8_9, allosteric_candidate |
| SC-4 | 핸드오프 가드 edge case | PASS | test_ppi_group3: 0 hotspot FAIL, 0 pocket FAIL, <5 poses FAIL |
| SC-5 | Orientation/sheet12 sensitivity 문서화 | PASS | sweep 6단계 + sheet12 2행 출력 구현 |
| SC-6 | Verdict 가중치 민감도 + 근거 문서 | PASS | 6조합 시뮬레이션 + cross_coverage 차등 근거 |

## 미충족 AC

없음.
