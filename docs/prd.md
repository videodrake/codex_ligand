# PRD: EGFR-MYO1D 파이프라인 Major Refactoring

## 1. Introduction & Overview

EGFR-MYO1D 결합 부위 탐색 파이프라인은 소분자 도킹(Vina)과 단백질-단백질 도킹(PyRosetta)의 두 독립 증거원을 운용하여, 약물 3종이 EGFR C-lobe 표면 어디에 결합하여 MYO1D를 해리시키는지를 계산적으로 탐색하는 통합 시스템이다. 현재 파이프라인은 프로덕션 수준의 자동화와 HPC 오케스트레이션을 갖추고 있으나, 심층 분석 보고서에서 37개의 과학적·기술적 개선 사항(FIX-01~37)이 식별되었다.

본 프로젝트는 이 37개 수정 사항을 6개 계획안(A~F)으로 묶어 체계적으로 구현한다. 핵심 목표는 세 가지이다:

1. **실험 데이터 반영**: 확보된 실험적 사실(ATP 결합 유지, Ko et al. alanine substitution)을 파이프라인 로직에 직접 통합
2. **편향 정량화 및 보정**: Vina blind docking의 구조적 편향과 PPI 예측의 파라미터 민감도를 측정하고 보정
3. **결과 통합 및 해석 체계 구축**: 두 워크플로우(A: 독립 병렬 탐색, B: PPI-guided 정밀 탐색)의 결과를 체계적으로 비교하는 모듈 구현

이 프로젝트는 기존 복잡한 코드를 수정하는 작업이며, 기존 프로덕션 결과와의 호환성을 유지하면서 과학적 신뢰도를 향상시키는 것이 핵심이다. 코드 변경은 대부분 기존 모듈에 기능을 추가하거나 검증 로직을 강화하는 형태이며, 재도킹(re-docking) 없이 기존 결과 데이터를 재분석하는 방식으로 진행된다.

## 2. Terminology

| 용어 | 정의 |
|------|------|
| Workflow A | Vina blind docking + PyRosetta PPI를 독립 병렬로 수행하고 Verdict로 통합하는 워크플로우 |
| Workflow B | PPI-first → Pocket Analysis → Focused Docking → Perturbation Scoring 순차 구조의 정밀 워크플로우 |
| Verdict | Workflow A의 최종 3축 스코어링 시스템 (Vina 50 + PPI 20 + Cross-Receptor 30) |
| Hotspot | PPI 도킹에서 반복적으로 관찰되는 interface 잔기. Phase 1에서 식별 |
| Orientation filter | MYO1D beta-meander의 active face가 EGFR 쪽을 향하는지 판정하는 PPI 후처리 필터 |
| Handoff | Workflow B의 Phase 간 데이터를 전달하는 CSV 파일 (3개: patch_reference, candidate_pocket, docking_evidence) |
| Active face | Ko et al. 실험으로 확인된 MYO1D의 EGFR 결합면 (sheet 8: 961-964, sheet 9: 968-972) |
| Bootstrap | Vina 포즈의 무작위 리샘플링을 통한 포켓 안정성 통계 검증 |
| 바이브코딩 | 기존 코드를 읽고 분석하여 판단을 내린 후 구현하는 작업 유형 (단순 자동화와 구분) |

## 3. User Personas

### 계산생물학 연구자 (Primary)
- **Who**: EGFR-MYO1D 프로젝트의 파이프라인을 운영하는 계산생물학자
- **Goal**: 파이프라인의 과학적 신뢰도를 높이고, 두 워크플로우의 결과를 체계적으로 비교하여 실험적 검증 우선순위를 결정
- **Pain Point**: 실험 데이터가 파이프라인에 반영되지 않아 결과 해석 시 수동 필터링이 필요하고, 두 워크플로우 결과를 수동으로 대조해야 함. Verdict의 가중치 근거가 불명확하여 결과 신뢰도에 대한 확신이 부족
- **Tech Comfort**: Python 능숙, HPC/PBS 운용 경험, 구조생물학 도메인 지식 보유

## 4. Feature Requirements

### Feature F-1: 실험 데이터 통합 및 즉시 검증
**Priority:** Must-Have
**Description:** 확보된 실험적 사실(리간드 파일, ATP 결합 유지, Ko et al. alanine substitution)을 파이프라인 로직과 문서에 직접 반영한다. 리간드 3종의 화학적 다양성을 검증하여 cross-chemical consensus의 유효성을 확인하고, ATP binding site 포즈를 자동으로 false positive 처리하며, Ko et al. 데이터를 핸드오프 검증과 포켓 해석에 활용한다.
**User Story:** As a 계산생물학 연구자, I want 실험적으로 확인된 사실이 파이프라인에 자동 반영되기를, so that 결과 해석 시 수동 필터링이 불필요하고 실험-계산 정합성이 보장된다.

**Acceptance Criteria:**
- [ ] AC-1.1: 리간드 3종(173940, 97806, VAX-C12_0)의 쌍별 Tanimoto similarity, MW, LogP, TPSA, HBD/HBA가 `ligand_diversity_assessment.csv`로 출력되며, 다양성 판정(< 0.4: 다양 / 0.4-0.7: 중간 / > 0.7: 유사)이 자동 수행된다
- [ ] AC-1.2: pocket_table.csv에 `is_atp_site` 컬럼이 추가되며, `region_definitions.py`의 ATP_SITE_RESIDUES(37개)와 접촉 잔기의 겹침 비율 > 50%인 포켓이 `True`로 표기된다. valid_sites.csv에 `exclusion_reason` 컬럼이 추가되어 해당 포켓이 "ATP_site_experimental"로 태깅된다
- [ ] AC-1.3: `_validate_adv_handoff()` 또는 Phase 1 핸드오프 생성 시, hotspot 잔기에 sheet 8/9 active face(961-964, 968-972) 잔기 3개 미만이면 FAIL이 발생하고 Workflow B가 중단된다. sheet 10/11 잔기가 hotspot에 포함되면 WARNING이 출력된다
- [ ] AC-1.4: research_overview_full.md와 CLAUDE.md에 "ATP 결합 유지 + 활성 소실" 실험 사실이 명시되어 있으며, Vina 결과 해석 가이드에 ATP site 포켓 배제 근거가 기술된다
- [ ] AC-1.5: pocket_table.csv에 `contacts_sheet_8_9`, `contacts_sheet_10_11`, `contacts_sheet_12` 컬럼이 추가되며, 각 포켓이 MYO1D beta-meander의 어떤 sheet와 접촉하는지 정량적으로 기록된다

**Edge Cases:**
- EC-1.1: 리간드 SDF 파일이 RDKit으로 파싱되지 않을 경우 → 에러 메시지와 함께 중단, 수동 구조 확인 요청
- EC-1.2: 포켓의 접촉 잔기가 ATP site와 C-lobe surface 영역에 걸쳐 있을 경우 (비율 40-60%) → `is_atp_site_borderline = True` 추가 태그, 수동 검토 권고
- EC-1.3: Ko et al. 체크에서 hotspot 잔기에 active face가 3개 이상 포함되지만 전부 낮은 occupancy인 경우 → WARNING("Active face 잔기 포함되나 occupancy 낮음")

**Dependencies:**
- Depends on: 없음 (즉시 시작 가능)
- Blocked by: 없음. 리간드 SDF 3종 및 ATP site 핵심 잔기(Lys745, Thr790, Gln791, Met793, Asp855)는 이미 확정됨

---

### Feature F-2: Vina Blind Docking 편향 정량화 및 보정
**Priority:** Must-Have
**Description:** Vina blind docking의 구조적 편향(강한 포켓이 약한 포켓을 압도)을 정량적으로 측정하고, 그 결과를 Verdict 점수에 반영한다. 포즈의 영역별 분포를 추적하여 C-lobe surface 포켓의 탐색 충분성을 확인하고, cross-receptor 잔기 번호 일관성을 강화하며, bootstrap 결과를 Verdict에 자동 연동한다.
**User Story:** As a 계산생물학 연구자, I want Vina blind docking의 편향 크기를 정량적으로 파악하기를, so that Workflow A 결과의 C-lobe surface 포켓 판정을 어디까지 신뢰할 수 있는지 근거에 기반하여 판단할 수 있다.

**Acceptance Criteria:**
- [ ] AC-2.1: `vina_pose_distribution_by_region.csv`가 생성되며, 각 포즈의 접촉 잔기를 `region_definitions.py` 기준으로 n_lobe/atp_site/c_lobe_surface/c_lobe_core로 분류한 영역별 포즈 수·비율·평균 affinity가 receptor_id × ligand_id별로 기록된다. C-lobe surface 포즈가 전체의 10% 미만이면 WARNING이 출력된다
- [ ] AC-2.2: `residue_alignment_check.csv`가 생성되며, 3개 receptor state 쌍별로 동일 잔기 번호에서 아미노산 타입이 다른 경우가 탐지된다. 불일치 발견 시 WARNING과 구체적 잔기 번호·구간이 보고된다. validate.py의 기존 (8.3) 체크와 통합된다
- [ ] AC-2.3: bootstrap이 실행된 경우 pocket_table.csv에 `bootstrap_confidence` 컬럼(high/medium/low)이 추가되며, Verdict의 vina_stability_pts 계산에 pocket_exists_frac이 반영된다(frac < 0.5 → 0점, 0.5-0.8 → 절반, > 0.8 → 만점). bootstrap 미실행 시 기존 로직이 유지된다(후방 호환)
- [ ] AC-2.4: C-lobe surface 포켓의 affinity 분포(25/50/75/90 percentile)가 분석되며, 현재 임계값(-8.0/-6.5/-5.0)이 이 분포에서 적절한 차별력을 갖는지 평가 결과가 문서화된다
- [ ] AC-2.5: 결과 해석 가이드에 Vina scoring function의 소수성 과대평가 편향과 C-lobe surface 포켓의 affinity 해석 지침이 기술된다

**Edge Cases:**
- EC-2.1: C-lobe surface 포즈가 0개인 경우 → WARNING("Vina blind docking에서 C-lobe surface에 포즈 미도달"), Workflow B의 focused docking 결과에 의존하라는 가이드 출력
- EC-2.2: 잔기 번호 불일치가 10개 이상 발견된 경우 → FAIL로 격상, PDB 파일 수동 확인 요청
- EC-2.3: bootstrap pocket_exists_frac가 모든 포켓에서 > 0.8인 경우 → bootstrap 연동이 기존 결과를 변화시키지 않음을 확인, 문서에 기록

**Dependencies:**
- Depends on: F-1 (ATP site 잔기 정의, region_definitions.py 필요)
- Blocked by: 없음

---

### Feature F-3: PPI Branch 강건성 검증 및 안전망
**Priority:** Must-Have
**Description:** PPI 도킹 결과의 파라미터 민감도를 검증하고, Workflow B 핸드오프의 데이터 품질 가드를 강화한다. Orientation filter threshold, fragment 범위, sheet 12 working assumption의 세 가지 핵심 파라미터에 대한 sensitivity 분석을 수행하고, 핸드오프 파일의 edge case 방어를 추가하며, hotspot occupancy 정규화와 cross-method agreement 정량화를 개선한다.
**User Story:** As a 계산생물학 연구자, I want PPI 예측 결과가 파라미터 선택에 robust한지 검증하고, Workflow B에서 오류가 전파되는 것을 자동으로 방어하기를, so that PPI 기반 결과를 신뢰하고 Workflow B를 안전하게 운영할 수 있다.

**Acceptance Criteria:**
- [ ] AC-3.1: orientation_filter.py에서 threshold [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]에 대한 sweep이 수행되며, 각 threshold별 pass/fail/ambiguous 비율과 hotspot 잔기 목록이 `orientation_threshold_sensitivity.csv`로 출력된다. threshold 간 hotspot robustness 점수가 계산된다
- [ ] AC-3.2: ambiguous 모델의 state/seed별 비율, 평균 dG_separated(pass 대비), unique interface 잔기 목록이 `orientation_ambiguous_report.csv`로 출력되며, Phase 1 리뷰 리포트에 "Ambiguous Models Summary" 섹션이 추가된다
- [ ] AC-3.3: 3가지 fragment 범위(945-1006, 955-1006, 955-1015)로 pilot PPI docking을 실행하기 위한 config YAML 3개와 PBS 스크립트가 생성되며, pilot 결과의 hotspot Jaccard similarity 비교 스크립트가 포함된다
- [ ] AC-3.4: active face 정의를 sheet 8+9 only(현재)와 sheet 8+9+12 두 설정으로 orientation filter를 재실행한 결과가 `sheet12_sensitivity.csv`로 출력되며, 기존 scored_all_models.csv를 재활용하여 재도킹 없이 비교된다
- [ ] AC-3.5: 3개 핸드오프 파일에 데이터 품질 가드가 추가된다. Phase 1: hotspot 잔기 0개이면 FAIL. Phase 2: 후보 포켓 0개이면 FAIL, 전부 irrelevant이면 WARNING. Phase 3: 유효 포즈 5개 미만이면 FAIL, 리간드 1종뿐이면 WARNING. 기존 `_validate_adv_handoff()`와 통합된다
- [ ] AC-3.6: Phase 2 분류 후 orthosteric+rim 비율이 80%를 초과하면 WARNING("PPI 패치 과대추정 가능성")이 출력된다
- [ ] AC-3.7: state-specific 분류이면서 해당 state에서 dG_separated 상위 10% 이내인 잔기에 `conformational_selection_candidate = True` 태그가 추가된다
- [ ] AC-3.8: ppi_hotspot_residues.csv에 `n_valid_models` 컬럼이 추가되며, seed/state 간 orientation-valid 모델 수가 2배 이상 차이나면 WARNING이 발생한다
- [ ] AC-3.9: cross_method_convergence.csv에 `pyrosetta_occupancy`, `lightdock_occupancy`, `method_concordance_score` 컬럼이 추가되며, "both" 분류 내에서 concordance > 0.5는 "strong_both", < 0.5는 "weak_both"로 세분화된다

**Edge Cases:**
- EC-3.1: 모든 threshold에서 pass 비율이 동일한 경우 → "Threshold에 무관한 결과, 현재 값 유지"로 판정
- EC-3.2: Pilot fragment 범위 3종의 hotspot Jaccard가 모두 < 0.3인 경우 → "Fragment 범위 선택이 결과에 결정적 영향, 구조적 검토 필요" WARNING
- EC-3.3: 핸드오프 품질 가드에서 기존 `_validate_adv_handoff()`의 반환 패턴(FAIL/WARNING/PASS)과 충돌하는 경우 → 기존 패턴에 맞춰 조정, 충돌 시 기존 방식 우선

**Dependencies:**
- Depends on: F-1 (Ko et al. 잔기 번호, sheet 정의)
- Blocked by: Fragment pilot (C-3.3)은 서버 작업 필요

---

### Feature F-4: Verdict 메트릭 체계 검증 및 교정
**Priority:** Must-Have
**Description:** Verdict 3축 스코어링(Vina 50/PPI 20/Cross-Receptor 30)의 가중치 민감도를 시뮬레이션하여 근거 기반 파라미터를 확정한다. Centroid 의미론적 오프셋의 포켓 깊이별 보정, cross-receptor 2/3 vs 3/3 차등 점수, PPI spatial 임계값 실효 범위를 검토하고, Phase 4의 A3(Perturbation relevance) 축 정의를 명확화한다.
**User Story:** As a 계산생물학 연구자, I want Verdict 가중치와 임계값이 정량적 근거에 기반하여 설정되기를, so that 포켓 판정(STRONG/MODERATE/WEAK)의 신뢰도가 높아지고, PPI 관련 포켓이 적절히 평가된다.

**Acceptance Criteria:**
- [ ] AC-4.1: 6개 가중치 조합(현재 50/20/30, PPI 강화 40/30/30, PPI 최대 35/35/30, Vina 약화 30/40/30, 균등 33/33/34 포함)에서 포켓별 총점과 STRONG/MODERATE/WEAK 판정이 재계산되며, 조합 간 판정 변화 포켓 목록이 `verdict_weight_sensitivity.csv`로 출력된다. PPI 가중치 증가로 새로 STRONG에 진입하는 포켓이 식별된다
- [ ] AC-4.2: pocket_table.csv에 `pocket_depth_A` 컬럼이 추가되며, centroid 오프셋 보정(corrected_distance = raw_distance - alpha × pocket_depth, alpha 0.5~1.0 시뮬레이션) 전후의 ppi_spatial_pts 변화가 비교된다
- [ ] AC-4.3: cross_receptor_pts에 2/3 vs 3/3 차등이 도입된다(예: 1/3 → 10-15점, 2/3 → 22-24점, 3/3 → 30점). 차등 전후 STRONG 경계(55점)에 걸친 포켓의 순위 변화가 문서화된다
- [ ] AC-4.4: Vina 포켓 centroid와 PPI cluster centroid 간 실제 거리 분포가 분석되며, 현재 임계값(8/15/25Å)이 이 분포를 적절히 분할하는지 평가된다. 하나의 구간에 포켓 70% 이상이 몰리면 구간 세분화가 권장된다
- [ ] AC-4.5: Phase 4의 A3(Perturbation relevance) 축의 계산 로직(입력 데이터, 임계값, 가중치, 출력 범위)이 코드에서 추출되어 `phase4_A3_axis_specification.md`로 문서화되며, PIPELINE_ARCHITECTURE_REPORT.md가 업데이트된다

**Edge Cases:**
- EC-4.1: 6개 가중치 조합에서 STRONG 포켓이 완전히 동일한 경우 → "가중치에 robust, 현재 배분 유지" 판정
- EC-4.2: 포켓 깊이 보정이 5개 이상 포켓의 판정을 변경하는 경우 → 보정 적용 전 개별 포켓의 생물학적 의미를 확인해야 함
- EC-4.3: A3 축의 코드 구현이 문서와 불일치하는 경우 → 불일치 내용을 명시하고, 코드 기준으로 문서를 수정

**Dependencies:**
- Depends on: F-1 (리간드 다양성 결과), F-2 (affinity 분포, 편향 크기), F-3 (PPI robustness 확인)
- Blocked by: F-2의 affinity 임계값 결정(B-4.2)이 완료되어야 Verdict 전체 조정 가능

---

### Feature F-5: 워크플로우 비교 및 문서 체계 정비
**Priority:** Should-Have
**Description:** Workflow A(valid_sites.csv)와 Workflow B(phase4_final_review_table.csv)의 결과를 포켓 단위로 체계적으로 비교하는 모듈을 구현하고, allosteric 메커니즘 해석 프레임워크를 추가하며, 방법론적 한계와 결과 면책 조항을 정비한다. 이를 통해 두 워크플로우의 설계 의도(상호 보완)가 운영 수준에서 실현된다.
**User Story:** As a 계산생물학 연구자, I want 두 워크플로우의 결과를 자동으로 비교하여 합의 포켓과 불일치 포켓을 식별하기를, so that 실험적 검증 우선순위를 체계적으로 결정하고, 결과의 한계를 명확히 소통할 수 있다.

**Acceptance Criteria:**
- [ ] AC-5.1: `workflow_comparison.py` 모듈이 구현되며, valid_sites.csv와 phase4_final_review_table.csv를 포켓 centroid 거리 + 잔기 Jaccard로 매칭하여 4가지(Consensus/A-only/B-only/Conflict)로 분류한 `workflow_comparison.csv`가 생성된다
- [ ] AC-5.2: 불일치 시나리오별(A=STRONG+B=irrelevant, B=상위+A=WEAK, 둘 다 무관심) 해석 가이드가 `workflow_comparison_guide.md`로 작성된다
- [ ] AC-5.3: Verdict에서 Vina 축 ≥ 35점 AND PPI 축 ≤ 5점인 포켓에 `allosteric_candidate = True` 태그가 추가되며, project_report.txt에 "Allosteric 후보 포켓" 섹션이 신설된다
- [ ] AC-5.4: 방법론적 한계(rigid-body, LightDock 독립성, 입력 구조 편향, solvent 효과, Vina scoring 편향) 5개 섹션이 `methodology_limitations.md`로 통합되며, CLAUDE.md와 research_overview_full.md에서 참조된다
- [ ] AC-5.5: valid_sites.csv, project_report.txt, phase4_final_review_table.csv의 최상단에 "계산적 예측이며 실험적 검증 필요" 면책 조항이 추가된다
- [ ] AC-5.6: `output_path_guide.md`가 작성되며 Workflow A/B/비교/Legacy 경로가 명확히 기술된다. paths.py의 legacy 함수에 deprecation 주석이 추가된다

**Edge Cases:**
- EC-5.1: 두 워크플로우의 포켓 centroid가 모두 8Å 이내에 매칭되지 않는 경우 → "미매칭" 포켓으로 분류, 각 워크플로우에서 독립적으로 해석
- EC-5.2: allosteric 후보가 0개인 경우 → 정상적 결과로 간주, 보고서에 "allosteric 후보 미식별" 기록
- EC-5.3: valid_sites.csv에 주석 행(#으로 시작)을 추가하면 downstream 파서가 실패하는 경우 → 주석 행 대신 별도 metadata 파일로 면책 조항 분리

**Dependencies:**
- Depends on: F-4 (Verdict 가중치 확정 후 최종 결과 기반 비교)
- Blocked by: Workflow B 프로덕션 결과가 존재해야 비교 가능

---

### Feature F-6: 인프라 및 장기 개선
**Priority:** Nice-to-Have
**Description:** 파이프라인의 장기 유지보수성과 운영 편의성을 개선한다. 대규모 단일 파일(pipeline_manager.py 4079줄, vina_executor.py 2813줄)의 모듈 분리를 검토하고, Workflow B의 외부 도구(fpocket, P2Rank) 설치 검증을 precheck에 포함하며, Phase 3 라운드 내 병렬화를 검토한다.
**User Story:** As a 계산생물학 연구자, I want 파이프라인의 코드 구조와 실행 환경 검증이 개선되기를, so that 유지보수가 용이하고 환경 설정 오류를 사전에 방지할 수 있다.

**Acceptance Criteria:**
- [ ] AC-6.1: pipeline_manager.py와 vina_executor.py의 내부 함수/클래스 경계가 분석되며, 분리가 필요한 경우(단일 함수 500줄 이상, 순환 참조 존재) 구체적 분리 계획이 문서화된다. 분리 시 기존 import 경로와 호환되는 wrapper가 제공된다
- [ ] AC-6.2: `--check-workflow-b` 또는 `--profile advanced` 플래그로 fpocket/P2Rank/LightDock 가용성을 검증하는 기능이 precheck에 추가된다. run_advanced_pipeline.pbs의 Phase 2 시작 전 자동 검증이 연결된다
- [ ] AC-6.3: Phase 3의 각 라운드 내 open 포켓 간 도킹이 현재 순차인지 병렬인지 확인되며, 순차인 경우 PBS job array 또는 ThreadPoolExecutor 기반 라운드 내 병렬화가 구현된다

**Edge Cases:**
- EC-6.1: 모듈 분리 후 기존 테스트가 실패하는 경우 → 분리 취소, 현재 구조 유지
- EC-6.2: fpocket은 설치되어 있으나 P2Rank이 없는 경우 → 개별 도구별 가용 여부를 보고, P2Rank 없이도 Phase 2가 동작 가능한지 확인
- EC-6.3: 라운드 내 병렬화가 이미 구현되어 있는 경우 → 문서화만 수행

**Dependencies:**
- Depends on: F-1~F-5 완료 후 여유가 있을 때
- Blocked by: 없음

## 5. User Flows

### Flow 1: 실험 데이터 통합 후 Workflow A 재분석
1. 연구자가 리간드 다양성 스크립트를 실행하여 3종의 화학적 다양성을 확인
2. 시스템이 `ligand_diversity_assessment.csv`를 생성하고 다양성 판정을 출력
3. 연구자가 파이프라인을 재실행(`--from 5`, Phase 5 Verdict부터)
4. 시스템이 ATP site 포켓을 자동 배제하고 Ko et al. 일관성 체크를 수행
5. 시스템이 영역별 포즈 분포, sheet 접촉 정보가 포함된 valid_sites.csv를 생성
→ 결과: 실험 데이터가 반영된 Verdict 판정, ATP site false positive 자동 제거

### Flow 2: PPI 파라미터 민감도 검증
1. 연구자가 orientation threshold sweep을 실행
2. 시스템이 6개 threshold별 결과를 `orientation_threshold_sensitivity.csv`로 출력
3. 연구자(또는 바이브코딩)가 hotspot 안정성을 확인하고 threshold를 확정
4. 연구자가 sheet 12 sensitivity 분석을 실행
5. 시스템이 sheet 8+9 vs sheet 8+9+12 결과를 비교하여 working assumption을 검증
→ 결과: 파라미터 robustness 확인, 필요 시 프로덕션 파라미터 조정

### Flow 3: 두 워크플로우 결과 비교
1. 연구자가 workflow_comparison.py를 실행
2. 시스템이 valid_sites.csv와 phase4_final_review_table.csv를 매칭
3. 시스템이 Consensus/A-only/B-only/Conflict 4가지로 분류한 비교 결과를 출력
4. 연구자가 Consensus 포켓을 최우선 실험 후보로 선정
5. A-only 포켓에 대해 allosteric 가능성을 검토
→ 결과: 실험적 검증 우선순위 목록, 불일치 해석

### Flow 4: Verdict 메트릭 교정
1. 연구자가 가중치 민감도 시뮬레이션을 실행
2. 시스템이 6개 조합별 포켓 순위와 판정 변화를 출력
3. 바이브코딩이 PPI 승격 포켓을 식별하고 권장안을 제시
4. 연구자가 PyMOL 확인 후 최종 가중치를 확정
5. 시스템이 확정된 가중치로 프로덕션 재실행
→ 결과: 근거 기반 Verdict 파라미터, 재실행된 valid_sites.csv

## 6. Non-Goals

- **AlphaFold-Multimer 통합**: 진정한 독립 검증을 위한 장기 과제이나, 현재 파이프라인의 rigid-body docking 프레임워크 내에서는 구현하지 않음
- **MM-PBSA/GBSA 자유에너지 재스코어링**: 상위 후보에 대한 정밀 에너지 계산은 본 프로젝트의 범위 밖
- **실험적 검증(HDX-MS, XL-MS, SPR/ITC)**: 계산 결과의 실험적 확인은 별도 프로젝트
- **새로운 도킹 실행**: 모든 분석은 기존 프로덕션 결과를 재활용. 재도킹은 fragment 범위 변경(C-3.4) 또는 Verdict 가중치 변경(D-1.3) 시에만 부분적으로 수행
- **전장 TH1 도킹**: beta-meander fragment(955-1006)에서 전장 MYO1D TH1 domain으로의 확장은 현재 수행하지 않음
- **GUI/웹 인터페이스**: 모든 작업은 CLI + PBS 스크립트 기반

## 7. Technical Considerations

- **언어 및 환경**: Python 3.x, HPC 클러스터 + PBS 스케줄러, RDKit/PyRosetta/Vina 환경 구축 완료
- **기존 코드 호환성**: 모든 수정은 기존 함수 시그니처와 반환값을 유지해야 하며, 새 컬럼/파일 추가 방식으로 구현. 기존 valid_sites.csv를 파싱하는 downstream 코드가 깨지지 않아야 함
- **재도킹 회피**: 대부분의 분석은 기존 pose_table.csv, scored_all_models.csv를 입력으로 사용. 4,500 Vina 포즈와 300K PPI 모델을 재생성하지 않음
- **확정된 상수**: ATP site 37개 잔기, Ko et al. sheet 잔기 번호(8:961-964, 9:968-972, 10:975-984, 11:985-991, 12:998-1004), SASA 기반 5영역 309잔기 분류가 이미 확정됨
- **바이브코딩 11건**: 기존 코드를 읽고 분석하여 판단을 내린 후 구현하는 작업. 자동화 코드와 달리 기존 로직의 이해가 선행되어야 함
- **사람 개입 항목**: PyMOL 시각 확인, 서버 pilot 실행, 과학적 최종 판단은 사람이 수행. 3건 완료(A-1.1, A-2.1, A-4.2), 1건(B-1.1) 완료됨

## 8. MVP Success Criteria

- [ ] SC-1: 모든 Must-Have Feature(F-1~F-4)의 AC가 통과한다
- [ ] SC-2: 기존 프로덕션 파이프라인(`run_production.py`)이 새 코드로 회귀 없이 실행된다
- [ ] SC-3: valid_sites.csv에 `is_atp_site`, `exclusion_reason`, `bootstrap_confidence`, `contacts_sheet_8_9` 컬럼이 추가되어 있다
- [ ] SC-4: 핸드오프 품질 가드(AC-3.5)가 빈 CSV, hotspot 0개, 유효 포즈 부족 케이스에서 FAIL을 반환한다
- [ ] SC-5: orientation threshold sensitivity와 sheet 12 sensitivity 분석 결과가 문서화되어 있다
- [ ] SC-6: Verdict 가중치 민감도 시뮬레이션 결과가 존재하고, 최종 가중치 선택 근거가 문서화되어 있다
