# Before/After Comparison — EGFR-MYO1D Pipeline Refactoring

> Group 0 baseline (2026-03-19) vs. 최종 상태 (2026-03-24)

## 1. Baseline 상태 (Group 0)

- **테스트:** 292 passed
- **valid_sites.csv:** 미존재 (output/ 비어 있음, 서버 미실행)
- **포켓 수 baseline:** N/A (프로덕션 결과 없음)
- **VERDICT_FIELDS:** 45 컬럼

## 2. 최종 상태

- **테스트:** 529 passed (+237)
- **valid_sites.csv:** 미존재 (서버 미실행 — 동일 상태)
- **VERDICT_FIELDS:** 46 컬럼 (+1)

## 3. 컬럼 변화

### 기존 컬럼 보존 (45개 전체 유지)
모든 baseline 컬럼이 삭제/이름변경 없이 유지됨.

### 신규 컬럼 (+1)
| 컬럼 | 그룹 | 설명 |
|------|------|------|
| `allosteric_candidate` | Group 5 | Vina ≥ 35 AND PPI ≤ 5 → True |

### pocket_table.csv 신규 컬럼
| 컬럼 | 그룹 | 설명 |
|------|------|------|
| `pocket_depth_A` | Group 4 | Centroid → 최근접 비접촉 표면 Cα 거리 |

## 4. 로직 변경

### Verdict scoring (verdict.py)
| 변경 | 그룹 | 영향 |
|------|------|------|
| cross_coverage 차등 (20/14/0) | Group 4 | 2/3 states: 20→14 pts, 1/3 states: 10→0 pts |
| allosteric_candidate 태깅 | Group 5 | 신규 boolean 컬럼 |

### Handoff quality guards (run_production.py)
| 가드 | 그룹 | 조건 |
|------|------|------|
| Phase 1→2 hotspot=0 | Group 3 | FAIL (FileNotFoundError) |
| Phase 2→3 pockets=0 | Group 3 | FAIL |
| Phase 2→3 all irrelevant | Group 3 | WARNING |
| Phase 3→4 poses<5 | Group 3 | FAIL |
| Phase 3→4 ligand=1 | Group 3 | WARNING |

### Report sections (report.py)
| 섹션 | 그룹 | 내용 |
|------|------|------|
| DISCLAIMER (최상단) | Group 5 | 면책 조항 |
| 4.7 Allosteric Candidates | Group 5 | Allosteric 후보 목록 |

## 5. 신규 모듈/스크립트

| 파일 | 그룹 | 유형 |
|------|------|------|
| `egfr_pipeline/region_definitions.py` | Group 0 | 상수 중앙화 |
| `egfr_pipeline/vina/pose_region_classifier.py` | Group 2 | 포즈 영역 분류 |
| `egfr_pipeline/workflow_comparison.py` | Group 5 | A↔B 비교 |
| `scripts/assess_ligand_diversity.py` | Group 1 | 리간드 다양성 |
| `scripts/analyze_pose_distribution.py` | Group 2 | 포즈 분포 |
| `scripts/analyze_affinity_distribution.py` | Group 2 | Affinity 분포 |
| `scripts/verdict_weight_sensitivity.py` | Group 4 | 가중치 민감도 |
| `scripts/centroid_offset_analysis.py` | Group 4 | Centroid 오프셋 |
| `scripts/pilot_fragment_range.py` | Group 3 | Fragment pilot |

## 6. 신규 문서

| 문서 | 그룹 |
|------|------|
| `docs/phase4_A3_axis_specification.md` | Group 4 |
| `docs/methodology_limitations.md` | Group 5 |
| `docs/workflow_comparison_guide.md` | Group 5 |
| `docs/workflow_comparison_design.md` | Group 5 |
| `docs/output_path_guide.md` | Group 5 |
| `docs/module_separation_analysis.md` | Group 6 |
| `docs/manual_vina.md` (수정) | Group 2 |

## 7. Verdict 가중치

| 항목 | Before | After | 근거 |
|------|--------|-------|------|
| Vina/PPI/Cross (with PPI) | 50/20/30 | 50/20/30 | 유지 (민감도 분석 결과 robust) |
| cross_coverage_3of3 | 20 | 20 | 유지 |
| cross_coverage_2of3 | 20 | 14 | 차등 도입 (AC-4.3) |
| cross_coverage_1of3 | 10 | 0 | coverage via support only |

가중치 변경이 있으나 전체 배분(50/20/30)은 유지. 민감도 시뮬레이션에서 6개 조합 테스트 완료.

## 8. 테스트 증가

| 그룹 | 테스트 파일 | 건수 |
|------|------------|------|
| Group 2 | test_vina_bias_group2.py | 39 |
| Group 3 | test_ppi_group3.py | 60 |
| Group 4 | test_verdict_group4.py | 43 |
| Group 5 | test_workflow_group5.py | 39 |
| Group 6 | test_infra_group6.py | 22 |
| 기존 수정 | test_adv_pipeline_validation, test_lightdock_validation | +4 |
| **합계** | | **+207** |
