# Task List — EGFR-MYO1D 파이프라인 Major Refactoring

## Dependency Order

```
Group 0 (Setup) ──→ Group 1 (F-1: 실험 데이터)
                         │
                         ├──→ Group 2 (F-2: Vina 편향) ──┐
                         │                                │
                         └──→ Group 3 (F-3: PPI 강건성) ──┤
                                                          │
                                                          v
                                              Group 4 (F-4: Verdict 교정)
                                                          │
                                                          v
                                              Group 5 (F-5: 비교 & 문서)
                                                          │
                                                          v
                                              Group 6 (F-6: 인프라) (Optional)
                                                          │
                                                          v
                                              Group 7 (E2E 통합)
                                                          │
                                                          v
                                              Group 8 (로깅/관찰성)

Group G (Dimer 수정) — 완료 ✅
  G-1 → G-5 → G-2 → G-6
                ↓
          G-3 → G-4
                ↓
               G-7

Group H (상태 오염 방지) — 독립 실행 가능, 우선순위 높음
  H-6 (PHASE_CORE_OUTPUTS 확장)
    → H-4 (_validate_lane_outputs, 모든 lane type)
      → H-1 (marker 쓰기 전 output 검증)
      → H-2 (skip 시 stale marker 감지)
  H-5 (seed marker 무결성) — 독립
  H-3 (PBS 로그) — 독립
```

## Project Structure

```
egfr_pipeline/
├── verdict.py                      [수정]
├── report.py                       [수정]
├── validate.py                     [수정]
├── region_definitions.py           [신규]
├── workflow_comparison.py          [신규]
├── paths.py                        [수정]
├── vina/
│   ├── pocket_summary.py           [수정]
│   ├── pocket_stability.py         [수정]
│   ├── vina_executor.py            [읽기] F-6 분리 검토 대상
│   └── pose_region_classifier.py   [신규]
├── phase1/
│   ├── orientation_filter.py       [수정]
│   ├── cluster_consensus.py        [수정]
│   ├── compare_states.py           [수정]
│   ├── review_report.py            [수정]
│   └── lightdock_validation.py     [수정]
├── phase2/                         [수정]
├── phase3/
│   └── run_diverse_docking.py      [수정]
├── phase4/                         [수정]
└── pyrosetta_docking/pipeline_manager.py [수정]
scripts/
├── assess_ligand_diversity.py      [신규]
├── analyze_pose_distribution.py    [신규]
├── analyze_affinity_distribution.py [신규]
├── verdict_weight_sensitivity.py   [신규]
├── centroid_offset_analysis.py     [신규]
└── pilot_fragment_range.py         [신규]
run_production.py                   [수정]
config/run_pre_qsub_checks.pbs     [수정]
config/run_advanced_pipeline.pbs    [수정]
docs/                               [신규 다수]
```

---

### 0. Project Setup 🟢

- [x] **0.1** — 기존 코드베이스 구조 파악 및 변경 대상 식별
  - **What:** 기존 파이프라인의 디렉토리 트리를 확인하고, 수정 대상 파일(verdict.py, report.py, validate.py, pocket_summary.py 등)의 현재 상태를 파악한다. 각 파일의 주요 함수/클래스 목록과 줄 수를 기록한다.
  - **Files:** `egfr_pipeline/` 전체 (읽기 전용 탐색)
  - **AC:** 셋업
  - **Test:** 디렉토리 트리 출력이 tech-stack의 Project Structure와 일치하며, 수정 대상 파일이 모두 존재

- [x] **0.2** — 확정된 상수 및 기존 산출물 확인
  - **What:** 이미 완료된 사람 개입 항목(A-1.1, A-2.1, A-4.2, B-1.1)의 산출물을 확인한다. ATP site 37잔기(region_definitions.py), Ko et al. sheet 번호, SDF 3종 존재, SASA 5영역 분류가 사용 가능한 상태인지 검증한다.
  - **Files:** `input/` (SDF 파일), `region_definitions.py` (B-1.1 산출물), 마스터 계획 상수 목록
  - **AC:** 셋업
  - **Test:** ATP_SITE_RESIDUES 37개 잔기 set, KO_SHEET_8~12 잔기 범위, SDF 3종 파일이 모두 로드 가능

- [x] **0.3** — 기존 테스트 스위트 실행 및 baseline 기록
  - **What:** 기존 파이프라인의 테스트를 실행하여 현재 상태에서 모두 통과하는지 확인한다. 현재 프로덕션 valid_sites.csv의 STRONG/MODERATE/WEAK 포켓 수를 baseline으로 기록한다.
  - **Files:** 기존 테스트 파일, `output/` 내 valid_sites.csv
  - **AC:** 셋업
  - **Test:** 기존 테스트 전체 통과, baseline 포켓 판정 수 기록 완료

- [x] **0.4** — region_definitions.py를 파이프라인 내부로 통합
  - **What:** B-1.1에서 생성된 `region_definitions.py`를 `egfr_pipeline/region_definitions.py`로 복사/이동하고, 다운스트림 모듈에서 import 가능하도록 한다. `REGION_LOOKUP`, `get_region()`, `ATP_SITE_RESIDUES` 등의 함수/상수가 정상 동작하는지 확인한다.
  - **Files:** `region_definitions.py` → `egfr_pipeline/region_definitions.py`
  - **AC:** 셋업
  - **Test:** `from egfr_pipeline.region_definitions import get_region, ATP_SITE_RESIDUES` 성공, `get_region(745)` → `"atp_site"`, `len(ATP_SITE_RESIDUES)` == 37

- [x] **0.T** — Tests for Project Setup
  - **What:** Setup 결과 검증: (1) 모든 수정 대상 파일 존재 확인, (2) 확정 상수 로드 테스트, (3) 기존 테스트 스위트 통과, (4) region_definitions import 테스트
  - **Files:** 테스트 파일
  - **AC:** 셋업
  - **Test:** pytest로 setup 검증 테스트 전체 통과

**Definition of Done:** 기존 코드베이스의 수정 대상이 모두 식별되고, 확정된 상수가 import 가능하며, 기존 테스트가 회귀 없이 통과한다.

---

### 1. 실험 데이터 통합 및 즉시 검증 (F-1) 🟡

- [x] **1.1** — 리간드 3종 화학 다양성 분석 스크립트
  - **What:** RDKit으로 3종 리간드의 SDF를 로드하여 쌍별 Morgan fingerprint(radius=2, nBits=2048) Tanimoto similarity, 개별 physicochemical descriptors(MW, LogP, TPSA, HBD, HBA, RotatableBonds, RingCount, FractionCSP3)를 계산한다. 다양성 판정(< 0.4: 다양 / 0.4-0.7: 중간 / > 0.7: 유사)을 자동 수행한다.
  - **Files:** `scripts/assess_ligand_diversity.py` [신규], `output/ligand_diversity_assessment.csv` [신규]
  - **AC:** AC-1.1
  - **Test:** 스크립트 실행 시 CSV 생성, 3개 쌍의 Tanimoto 값이 0-1 범위, 판정 문자열("diverse"/"moderate"/"similar") 포함

- [x] **1.2** — ATP site false positive 파이프라인 통합
  - **What:** (A-2.2) pocket_summary.py에서 포켓의 union_contact_residues 중 ATP_SITE_RESIDUES 비율 > 50% → `is_atp_site = True` 플래그 추가. (A-2.3) verdict.py에서 `is_atp_site = True` 포켓에 `exclusion_reason = "ATP_site_experimental"` 태그 부여, valid_sites.csv에 컬럼 추가. (A-2.4) report.py에서 project_report.txt에 "실험적 배제 포켓" 섹션 추가.
  - **Files:** `egfr_pipeline/vina/pocket_summary.py` [수정], `egfr_pipeline/verdict.py` [수정], `egfr_pipeline/report.py` [수정]
  - **AC:** AC-1.2
  - **Test:** ATP site 잔기가 주로 접촉하는 포켓에 `is_atp_site=True` 설정, valid_sites.csv에 `exclusion_reason` 컬럼 존재, report에 배제 섹션 포함

- [x] **1.3** — ATP 실험 사실 핵심 문서 반영
  - **What:** (A-3.1) research_overview_full.md 섹션 1.2에 실험적 사실 3("ATP 결합 유지 + 활성 소실") 추가. CLAUDE.md에 "실험적 근거" 요약 추가. Vina 결과 해석 가이드에 ATP site 배제 근거 명시.
  - **Files:** `docs/research_overview_full.md` [수정], `CLAUDE.md` [수정], 결과 해석 가이드 [수정]
  - **AC:** AC-1.4
  - **Test:** 각 문서에서 "ATP 결합 유지" 또는 "ATP binding maintained" 문구 검색 가능

- [x] **1.4** — Ko et al. 일관성 체크 로직 구현
  - **What:** (A-4.1) `_validate_adv_handoff()` 또는 Phase 1 핸드오프 생성 시 Ko et al. 체크 추가: hotspot에 sheet 8/9 active face(961-964, 968-972) 잔기 3개 미만 → FAIL, sheet 10/11(975-991) 잔기가 hotspot → WARNING, sheet 12(998-1004) 잔기 감지 → INFO.
  - **Files:** `run_production.py` [수정] (`_validate_adv_handoff` 함수)
  - **AC:** AC-1.3
  - **Test:** 테스트 CSV(active face 2개만 포함) → FAIL 반환, 정상 CSV(active face 5개) → PASS 반환, sheet 10 포함 CSV → WARNING 반환

- [x] **1.5** — 포켓별 Ko et al. sheet 접촉 정보 추가
  - **What:** (A-5.1) pocket_table.csv에 `contacts_sheet_8_9`, `contacts_sheet_10_11`, `contacts_sheet_12` 컬럼 추가. 각 Vina 포켓의 MYO1D side 접촉 잔기를 KO_SHEET_8|9, KO_SHEET_10|11, KO_SHEET_12 기준으로 분류. (A-5.2) 해석 가이드에 sheet별 의미 추가.
  - **Files:** `egfr_pipeline/vina/pocket_summary.py` 또는 `egfr_pipeline/verdict.py` [수정]
  - **AC:** AC-1.5
  - **Test:** pocket_table.csv에 3개 신규 컬럼 존재, 정수값 기록, sheet 8/9 잔기와 접촉하는 포켓에서 `contacts_sheet_8_9 > 0`

- [x] **1.T** — Tests for 실험 데이터 통합
  - **What:** (1) ligand_diversity_assessment.csv 생성 및 포맷 검증, (2) is_atp_site 판정 로직 단위 테스트 (ATP site 잔기 비율 50% 경계), (3) Ko et al. 체크 3가지 시나리오(PASS/WARNING/FAIL), (4) sheet 접촉 컬럼 정합성, (5) 기존 valid_sites.csv 파서가 신규 컬럼으로 깨지지 않는지 후방 호환 테스트
  - **Files:** 테스트 파일
  - **AC:** AC-1.1 ~ AC-1.5
  - **Test:** pytest 전체 통과

**Definition of Done:** 리간드 다양성이 정량적으로 평가되고, ATP site 포켓이 자동 배제되며, Ko et al. 일관성 체크가 핸드오프 검증에 통합되어 실험 데이터가 파이프라인 로직에 직접 반영된다.

---

### 2. Vina Blind Docking 편향 정량화 (F-2) 🟡

- [x] **2.1** — 포즈 영역별 분포 분석 및 보고
  - **What:** (B-1.2) `region_definitions.py`의 `get_region()`을 사용하여 vina_pose_table.csv 각 포즈의 contact_residues를 n_lobe/atp_site/c_lobe_surface/c_lobe_core로 분류. 영역별 포즈 수·비율·평균 affinity를 receptor_id × ligand_id별로 `vina_pose_distribution_by_region.csv` 출력. C-lobe surface < 10% → WARNING. (B-1.4) report.py에 분포 요약 섹션 추가.
  - **Files:** `scripts/analyze_pose_distribution.py` 또는 `egfr_pipeline/vina/pose_region_classifier.py` [신규], `egfr_pipeline/report.py` [수정]
  - **AC:** AC-2.1
  - **Test:** CSV 생성, 영역별 fraction 합계 ≈ 1.0, report에 분포 섹션 포함

- [x] **2.2** — Cross-receptor 잔기 번호 구간별 비교 + validate.py 통합
  - **What:** (B-2.1) 3개 receptor state PDB에서 (잔기번호, 아미노산) 매핑 추출, state 쌍별 불일치 탐지, `residue_alignment_check.csv` 출력. (B-2.2 🔧바이브코딩) validate.py의 (8.3) 기존 체크 코드를 읽고, 구간별 비교 로직을 적절히 통합. 기존 종료 코드 체계(0/1/2) 호환.
  - **Files:** `egfr_pipeline/validate.py` [수정]
  - **AC:** AC-2.2
  - **Test:** 인위적으로 loop 영역 1개 잔기 아미노산을 변경한 테스트 PDB → WARNING 발생 확인

- [x] **2.3** — Bootstrap 결과와 Verdict 자동 연동
  - **What:** (B-3.1) pocket_stability.py의 pocket_exists_frac을 pocket_table.csv에 병합. verdict.py에서 vina_stability_pts에 frac 반영 (< 0.5 → 0점, 0.5-0.8 → 절반, > 0.8 → 만점). `bootstrap_confidence` 컬럼 추가. (B-3.2 🔧바이브코딩) verdict.py의 현재 vina_stability_pts 계산 로직을 확인하고, bootstrap 연동 전후 기존 프로덕션 결과(valid_sites.csv) 비교. 판정 변화 3개 이상이면 임계값 재조정.
  - **Files:** `egfr_pipeline/vina/pocket_stability.py` [수정], `egfr_pipeline/verdict.py` [수정]
  - **AC:** AC-2.3
  - **Test:** bootstrap 미실행 프로젝트에서 `bootstrap_confidence = "not_assessed"`, 기존 점수 변화 없음 확인. bootstrap 실행 프로젝트에서 frac < 0.5 포켓의 stability_pts = 0 확인

- [x] **2.4** — C-lobe surface affinity 분포 분석 및 임계값 검토
  - **What:** (B-4.1) C-lobe surface 포켓(is_atp_site=False)의 affinity 분포 추출, 25/50/75/90 percentile 계산, 현재 임계값(-8.0/-6.5/-5.0) 위치 확인. (B-4.2 🔧바이브코딩) 70% 이상이 하나의 구간에 몰리면 세분화 제안, 고르면 유지. 조정 시 전후 Verdict 순위 변화 비교. 결정을 문서화.
  - **Files:** `scripts/analyze_affinity_distribution.py` [신규], `egfr_pipeline/verdict.py` [수정 가능]
  - **AC:** AC-2.4
  - **Test:** 분포 분석 결과에 percentile 값 포함, 차별력 판정("sufficient"/"insufficient") 출력

- [x] **2.5** — Vina scoring 편향 문서화
  - **What:** (B-5.1) 결과 해석 가이드에 Vina scoring function의 소수성 과대평가 편향, C-lobe surface 포켓의 affinity 해석 지침(-5~-7 kcal/mol도 표면 포켓으로서 의미 있음) 기술.
  - **Files:** 결과 해석 가이드 문서 [수정]
  - **AC:** AC-2.5
  - **Test:** 문서에 "소수성 과대평가" 또는 "hydrophobic bias" 관련 기술 존재

- [x] **2.T** — Tests for Vina 편향 정량화
  - **What:** (1) 포즈 분포 CSV 포맷 및 영역별 fraction 합계 검증, (2) 잔기 매핑 비교의 불일치 감지 단위 테스트, (3) bootstrap 연동의 후방 호환 테스트, (4) 기존 valid_sites.csv와 비교하여 의도치 않은 판정 뒤집힘 없음 확인, (5) 영역 분류 일관성 (동일 잔기 → 동일 영역)
  - **Files:** 테스트 파일
  - **AC:** AC-2.1 ~ AC-2.5
  - **Test:** pytest 전체 통과

**Definition of Done:** Vina blind docking의 편향이 정량적으로 파악되어 영역별 포즈 분포가 보고되고, bootstrap-Verdict 연동이 자동화되며, C-lobe surface 포켓에 적합한 affinity 임계값이 검토되어 문서화된다.

---

### 3. PPI Branch 강건성 검증 및 안전망 (F-3) 🔴

- [x] **3.1** — Orientation filter threshold 민감도 분석
  - **What:** (C-1.1) orientation_filter.py에 threshold sweep 모드 추가 [0.05~0.30], 각 threshold별 pass/fail/ambiguous 비율 + hotspot 잔기 → `orientation_threshold_sensitivity.csv`. (C-1.2) threshold 간 hotspot 변화 비교, threshold_robustness 점수 계산. (C-1.3 🔧바이브코딩) pass 비율 민감도, sheet 8/9 안정성, ambiguous 비율 기준으로 최적 threshold 판정 초안 생성.
  - **Files:** `egfr_pipeline/phase1/orientation_filter.py` [수정]
  - **AC:** AC-3.1
  - **Test:** sweep 결과 CSV에 6개 threshold 행, hotspot_residues 컬럼이 세미콜론 구분 잔기 목록

- [x] **3.2** — Ambiguous 모델 모니터링 및 보고
  - **What:** (C-2.1) state/seed별 ambiguous 비율, ambiguous vs pass 평균 dG_separated, ambiguous-only unique 잔기 → `orientation_ambiguous_report.csv`. (C-2.2) Phase 1 리뷰 리포트에 "Ambiguous Models Summary" 섹션 추가. unique 잔기 3개 이상이면 "추가 조사 권장" 주석.
  - **Files:** `egfr_pipeline/phase1/orientation_filter.py` [수정], `egfr_pipeline/phase1/review_report.py` [수정]
  - **AC:** AC-3.2
  - **Test:** ambiguous report CSV 생성, unique 잔기 목록이 비어 있지 않으면 리포트에 주석 포함

- [x] **3.3** — Fragment 범위 sensitivity pilot 준비
  - **What:** (C-3.2) 3가지 범위(945-1006, 955-1006, 955-1015)로 config YAML 3개 + PBS 스크립트 생성. 결과 비교 스크립트(hotspot Jaccard, centroid 거리, orientation pass 비율). 서버 실행(C-3.3)은 사람이 수행. (C-3.4 🔧바이브코딩) pilot 결과 Jaccard 기준 범위 확정 판정 초안.
  - **Files:** `scripts/pilot_fragment_range.py` [신규], config YAML 3개 [신규]
  - **AC:** AC-3.3
  - **Test:** 3개 config YAML 생성, PBS 스크립트 문법 검증, 비교 스크립트에 Jaccard 계산 함수 포함

- [x] **3.4** — Sheet 12 포함/제외 sensitivity 분석
  - **What:** (C-4.1) 설정 A(sheet 8+9 only)와 설정 B(sheet 8+9+12) 두 active_face로 orientation filter 재실행(기존 scored_all_models.csv 재활용, 재도킹 불필요). `sheet12_sensitivity.csv` 출력. (C-4.2 🔧바이브코딩) hotspot 겹침 비율 계산 + 판정 초안(≥80% → 유지, 50-80% → 검토, <50% → 병행).
  - **Files:** `egfr_pipeline/phase1/orientation_filter.py` [수정]
  - **AC:** AC-3.4
  - **Test:** sheet12_sensitivity.csv에 2행(설정 A/B), hotspot_overlap 컬럼 존재

- [x] **3.5** — 핸드오프 데이터 품질 가드 + 기존 검증 통합
  - **What:** (C-5.1) 3개 핸드오프에 데이터 품질 체크: Phase 1 hotspot 0개 → FAIL, Phase 2 포켓 0개 → FAIL / 전부 irrelevant → WARNING, Phase 3 유효 포즈 5개 미만 → FAIL / 리간드 1종 → WARNING. (C-5.2 🔧바이브코딩) `_validate_adv_handoff()` 읽고 삽입 위치 결정, 기존 FAIL/WARNING/PASS 패턴 호환, 빈 CSV/hotspot 0개/정상 3케이스 테스트.
  - **Files:** `run_production.py` [수정] (`_validate_adv_handoff`)
  - **AC:** AC-3.5
  - **Test:** 빈 CSV → FAIL, hotspot 0개 CSV → FAIL, 정상 CSV → PASS, 전부 irrelevant → WARNING

- [x] **3.6** — PPI 모니터링 확장 (패치 감지, conf selection, occupancy, cross-method)
  - **What:** (C-6.1) Phase 2 분류 후 orthosteric+rim > 80% → WARNING. (C-7.1) state-specific이면서 dG_separated 상위 10% → `conformational_selection_candidate = True`. (C-8.1) hotspot CSV에 `n_valid_models` 컬럼, seed/state 간 2배 차이 → WARNING. (C-9.1) cross_method_convergence에 pyrosetta/lightdock occupancy + concordance_score 추가, both 세분화.
  - **Files:** `egfr_pipeline/phase2/` [수정], `egfr_pipeline/phase1/compare_states.py` [수정], `egfr_pipeline/phase1/cluster_consensus.py` [수정], `egfr_pipeline/phase1/lightdock_validation.py` [수정]
  - **AC:** AC-3.6, AC-3.7, AC-3.8, AC-3.9
  - **Test:** 각 WARNING 조건에 대한 단위 테스트, conformational_selection_candidate 태그 존재, concordance_score 범위 0-1

- [x] **3.T** — Tests for PPI 강건성 검증
  - **What:** (1) threshold sweep 출력 포맷 검증, (2) ambiguous report 정합성, (3) pilot 스크립트 생성 검증, (4) sheet12 sensitivity 2행 출력, (5) 핸드오프 품질 가드 edge case (빈 CSV, 0개 hotspot, 전부 irrelevant), (6) 기존 _validate_adv_handoff 후방 호환
  - **Files:** 테스트 파일
  - **AC:** AC-3.1 ~ AC-3.9
  - **Test:** pytest 전체 통과

**Definition of Done:** Orientation threshold, fragment 범위, sheet 12의 sensitivity가 검증되고, 핸드오프 데이터 품질 가드가 활성화되어, PPI 분석의 파라미터 robustness가 확인되고 Workflow B의 오류 전파가 방어된다.

---

### 4. Verdict 메트릭 체계 검증 및 교정 (F-4) 🟡

- [x] **4.1** — Verdict 가중치 민감도 시뮬레이션
  - **What:** (D-1.1) 기존 valid_sites.csv를 입력으로 6개 가중치 조합(50/20/30, 40/30/30, 35/35/30, 30/40/30, 33/33/34 + 1개 추가)에서 포켓별 총점·판정 재계산. 조합 간 판정 변화 목록 + PPI 승격 포켓 식별 → `verdict_weight_sensitivity.csv`. (D-1.2 🔧바이브코딩→사람) 안정·민감·승격 포켓 분석 + 권장안 → `verdict_weight_sensitivity_report.md`.
  - **Files:** `scripts/verdict_weight_sensitivity.py` [신규]
  - **AC:** AC-4.1
  - **Test:** 6개 조합 결과 포함 CSV, 판정 변화 포켓 목록 비어 있지 않음 (또는 "변화 없음" 기록)

- [x] **4.2** — 포켓 깊이 지표 + centroid 오프셋 보정 시뮬레이션
  - **What:** (D-2.1) pocket_summary.py에 `pocket_depth_A` 계산 추가 (centroid → 가장 가까운 비접촉 수용체 표면 Cα 거리). pocket_table.csv에 컬럼 추가. (D-2.2) 보정 거리(raw - alpha × depth, alpha 0.5~1.0)에서 ppi_spatial_pts 변화 비교.
  - **Files:** `egfr_pipeline/vina/pocket_summary.py` [수정], `scripts/centroid_offset_analysis.py` [신규]
  - **AC:** AC-4.2
  - **Test:** pocket_depth_A 값이 양수, 보정 전후 점수 비교 테이블 생성

- [x] **4.3** — Cross-receptor pts 2/3 vs 3/3 차등 도입
  - **What:** (D-3.1) verdict.py의 cross_receptor_pts에 차등 구현: 1/3 → 10-15점, 2/3 → 22-24점, 3/3 → 30점 (변경안 2개 시뮬레이션). 차등 전후 STRONG 경계(55점)에 걸친 포켓 순위 변화 확인.
  - **Files:** `egfr_pipeline/verdict.py` [수정]
  - **AC:** AC-4.3
  - **Test:** 2/3 포켓과 3/3 포켓의 점수가 다름, 기존 동점이었던 포켓 쌍에서 차등 발생

- [x] **4.4** — PPI spatial pts 실효 범위 분석
  - **What:** (D-4.1) Vina-PPI centroid 간 거리 분포 추출, 현재 임계값(8/15/25Å)의 분포 분할 적절성 확인. (D-4.2 🔧바이브코딩→사람) 구간별 포켓 집중도 분석, 차별력 부족 시 새 임계값 제안, centroid 오프셋(3-5Å) 감안한 실질 범위 보고.
  - **Files:** `scripts/verdict_weight_sensitivity.py`에 통합 또는 별도 스크립트
  - **AC:** AC-4.4
  - **Test:** 거리 분포 히스토그램 또는 percentile 테이블 생성, 구간별 포켓 수 기록

- [x] **4.5** — Phase 4 A3 축 정의 추출 및 문서화
  - **What:** (D-5.1 🔧바이브코딩) `egfr_pipeline/phase4/` 모듈을 읽어 A3(Perturbation relevance) 계산 로직 추적: 입력 데이터, 메커니즘 분류→점수 변환, sub-score 가중치, 출력 범위 추출. (D-5.2) `phase4_A3_axis_specification.md` 작성, PIPELINE_ARCHITECTURE_REPORT.md 업데이트.
  - **Files:** `egfr_pipeline/phase4/` [읽기], `docs/phase4_A3_axis_specification.md` [신규]
  - **AC:** AC-4.5
  - **Test:** A3 스펙 문서에 입력/임계값/가중치/출력범위 4가지 항목 명시

- [x] **4.T** — Tests for Verdict 메트릭 교정
  - **What:** (1) 가중치 시뮬레이션 결과 6개 조합 검증, (2) pocket_depth_A 양수 검증, (3) 차등 점수 로직 단위 테스트 (n_receptors=1,2,3별 점수), (4) 거리 분포 통계 유효성, (5) 기존 valid_sites.csv와 비교하여 변경 사항 추적
  - **Files:** 테스트 파일
  - **AC:** AC-4.1 ~ AC-4.5
  - **Test:** pytest 전체 통과

**Definition of Done:** Verdict의 3축 가중치와 세부 임계값이 정량적 근거에 기반하여 검증되고, 필요한 조정이 문서화되며, Phase 4 A3 축의 정의가 명확히 기록된다.

---

### 5. 워크플로우 비교 및 문서 체계 정비 (F-5) 🟡

- [x] **5.1** — Workflow A↔B 비교 모듈 설계
  - **What:** (E-1.1 🔧바이브코딩→사람) valid_sites.csv와 phase4_final_review_table.csv의 컬럼 구조 확인. 매칭 방법 옵션(centroid 거리 + Jaccard vs centroid만) 비교. ATP site 포켓 제외 여부, irrelevant 포켓 포함 여부, "상위" 기준 등 설계 결정을 정리하여 `workflow_comparison_design.md` 출력.
  - **Files:** `docs/workflow_comparison_design.md` [신규]
  - **AC:** AC-5.1 (설계 파트)
  - **Test:** 설계 문서에 매칭 방법, 분류 기준, 제외 조건이 명시

- [x] **5.2** — Workflow A↔B 비교 모듈 구현
  - **What:** (E-1.2) `workflow_comparison.py` 구현. 5.1에서 확정된 매칭 방법으로 두 CSV의 포켓을 매칭하여 4가지(Consensus/A-only/B-only/Conflict)로 분류. A-only에 "allosteric 후보?" 플래그, B-only에 "blind docking 편향?" 플래그. 출력: `workflow_comparison.csv` + `workflow_comparison_report.md`.
  - **Files:** `egfr_pipeline/workflow_comparison.py` [신규]
  - **AC:** AC-5.1 (구현 파트)
  - **Test:** 테스트 입력(valid_sites + phase4_table)으로 4가지 분류 결과 생성, CSV 포맷 검증

- [x] **5.3** — Allosteric 후보 분류 및 보고
  - **What:** (E-3.1) verdict.py에서 Vina 축 ≥ 35 AND PPI 축 ≤ 5인 포켓에 `allosteric_candidate = True`. valid_sites.csv에 컬럼 추가. (E-3.2) report.py의 project_report.txt에 "Allosteric 후보 포켓" 섹션 신설.
  - **Files:** `egfr_pipeline/verdict.py` [수정], `egfr_pipeline/report.py` [수정]
  - **AC:** AC-5.3
  - **Test:** Vina 축 높고 PPI 축 낮은 포켓에 태그 설정, report에 섹션 존재

- [x] **5.4** — 방법론적 한계 통합 문서 + 해석 가이드
  - **What:** (E-4.1) 5개 섹션(rigid-body, LightDock, 입력 구조, solvent, Vina scoring) → `methodology_limitations.md`. cross_method_convergence에 independence_level 메타데이터 추가. (E-4.2) CLAUDE.md, research_overview_full.md에서 참조. (E-2.1) 불일치 시나리오별 해석 가이드 → `workflow_comparison_guide.md`.
  - **Files:** `docs/methodology_limitations.md` [신규], `docs/workflow_comparison_guide.md` [신규], `CLAUDE.md` [수정]
  - **AC:** AC-5.2, AC-5.4
  - **Test:** 두 문서 존재, 5개 한계 섹션 포함, 3가지 불일치 시나리오(A=STRONG B=irrelevant 등) 기술

- [x] **5.5** — 결과 면책 조항 추가
  - **What:** (E-5.1) valid_sites.csv, project_report.txt, phase4_final_review_table.csv 최상단에 "계산적 예측이며 실험적 검증 필요" 면책 조항 추가. downstream 파서와 충돌 시 별도 metadata 파일로 분리.
  - **Files:** `egfr_pipeline/report.py` [수정], `egfr_pipeline/verdict.py` [수정]
  - **AC:** AC-5.5
  - **Test:** 생성된 파일 최상단에 "DISCLAIMER" 또는 "면책" 문구 포함

- [x] **5.6** — Output 경로 가이드 + Legacy deprecation
  - **What:** (E-6.1) Workflow A/B/비교/Legacy 경로를 명확히 기술하는 `output_path_guide.md` 작성. (E-6.2) paths.py의 legacy 함수에 `# DEPRECATED` 주석 추가.
  - **Files:** `docs/output_path_guide.md` [신규], `egfr_pipeline/paths.py` [수정]
  - **AC:** AC-5.6
  - **Test:** 경로 가이드에 4개 섹션(Workflow A/B/비교/Legacy), paths.py에 DEPRECATED 주석 존재

- [x] **5.T** — Tests for 워크플로우 비교 및 문서
  - **What:** (1) workflow_comparison.py 출력 포맷 검증 (4가지 분류), (2) allosteric 태그 로직 단위 테스트, (3) 면책 조항이 파일 최상단에 위치, (4) 문서 파일 존재 및 필수 섹션 포함, (5) paths.py legacy 함수 deprecation 주석 존재
  - **Files:** 테스트 파일
  - **AC:** AC-5.1 ~ AC-5.6
  - **Test:** pytest 전체 통과

**Definition of Done:** 두 워크플로우의 결과가 체계적으로 비교되어 Consensus/A-only/B-only/Conflict로 분류되고, allosteric 후보가 식별되며, 방법론적 한계와 면책 조항이 정비되어 결과의 소통 품질이 확보된다.

---

### 6. 인프라 및 장기 개선 (F-6) (Optional) 🟢

- [x] **6.1** — 대규모 파일 모듈 분리 검토
  - **What:** (F-1.1) pipeline_manager.py(4079줄)와 vina_executor.py(2813줄)의 내부 함수/클래스 경계를 분석. 단일 함수 500줄 이상, 순환 참조 여부 확인. 분리 필요 시 구체적 분리 계획 문서화. 불필요 시 "현재 구조 유지" 결정 + 근거 기록.
  - **Files:** `egfr_pipeline/pyrosetta_docking/pipeline_manager.py` [읽기], `egfr_pipeline/vina/vina_executor.py` [읽기]
  - **AC:** AC-6.1
  - **Test:** 분석 결과 문서 존재, 결정(분리/유지) 근거 포함

- [x] **6.2** — fpocket/P2Rank/LightDock 설치 검증 precheck
  - **What:** (F-2.1) precheck에 `--check-workflow-b` 플래그 추가, fpocket/P2Rank/LightDock 가용성 체크. (F-2.2) run_advanced_pipeline.pbs의 Phase 2 시작 전 자동 검증 연결.
  - **Files:** `config/run_pre_qsub_checks.pbs` 또는 `scripts/run_pre_qsub_checks.sh` [수정], `config/run_advanced_pipeline.pbs` [수정]
  - **AC:** AC-6.2
  - **Test:** `--check-workflow-b` 플래그로 실행 시 fpocket/P2Rank 경로 체크 출력

- [x] **6.3** — Phase 3 라운드 내 병렬화 검토
  - **What:** (F-3.1) `run_diverse_docking.py`에서 각 라운드 내 open 포켓 도킹이 순차인지 병렬인지 확인. (F-3.2) 순차인 경우 PBS job array 또는 ThreadPoolExecutor 기반 라운드 내 병렬화 구현.
  - **Files:** `egfr_pipeline/phase3/run_diverse_docking.py` [읽기/수정]
  - **AC:** AC-6.3
  - **Test:** 이미 병렬이면 확인 문서화, 병렬화 추가 시 순차 대비 실행 시간 단축 확인

- [x] **6.T** — Tests for 인프라 개선
  - **What:** (1) precheck 플래그 동작 확인, (2) 모듈 분리 시 기존 import 호환 테스트, (3) 병렬화 시 결과 일관성 확인
  - **Files:** 테스트 파일
  - **AC:** AC-6.1 ~ AC-6.3
  - **Test:** pytest 전체 통과

**Definition of Done:** 대규모 파일의 분리 여부가 결정되고, Workflow B 외부 도구 검증이 precheck에 통합되며, Phase 3 라운드 내 병렬화 수준이 확인/개선된다.

---

### 7. E2E 통합 및 최종 검증 🔴

- [x] **7.1** — 전체 파이프라인 회귀 테스트
  - **What:** 모든 수정이 반영된 코드로 `run_production.py`를 실행하여 기존 기능이 정상 동작하는지 확인. 새로 추가된 컬럼(is_atp_site, exclusion_reason, bootstrap_confidence, contacts_sheet_*, allosteric_candidate)이 valid_sites.csv에 포함되고, 기존 컬럼은 변경되지 않음을 확인.
  - **Files:** `run_production.py`, 전체 output
  - **AC:** SC-2 (기존 파이프라인 회귀 없이 실행)
  - **Test:** 프로덕션 실행 완료, 기존 컬럼 값이 baseline과 일치 (또는 의도된 변경만 존재)

- [x] **7.2** — Before/After Verdict 결과 비교
  - **What:** Group 0에서 기록한 baseline과 최종 결과를 비교. STRONG/MODERATE/WEAK 포켓 수 변화, ATP site 배제 포켓 수, allosteric 후보 수를 정리. Verdict 가중치 변경이 있었다면 해당 영향도 분석. 비교 결과를 `before_after_comparison.md`로 출력.
  - **Files:** baseline 기록, 최종 valid_sites.csv
  - **AC:** SC-3, SC-6
  - **Test:** 비교 문서에 포켓 수 변화 테이블 포함

- [x] **7.3** — AC 커버리지 검증 및 최종 문서
  - **What:** PRD의 모든 Must-Have AC(AC-1.1~AC-4.5)가 충족되었는지 체크리스트 확인. Should-Have(AC-5.1~AC-5.6) 충족 여부 기록. 미충족 AC가 있으면 사유 기록. 최종 문서 일관성 검토 (문서 간 상호 참조 정합).
  - **Files:** `docs/prd.md`, 전체 산출물
  - **AC:** SC-1
  - **Test:** Must-Have AC 전체 통과, 미충족 사유 없음

- [x] **7.T** — E2E 통합 테스트
  - **What:** (1) 전체 테스트 스위트 실행 (기존 + 신규), (2) 핸드오프 품질 가드 edge case 통합 테스트, (3) valid_sites.csv 파서 호환성, (4) workflow_comparison.py 실 데이터 실행, (5) 문서 파일 전체 존재 확인
  - **Files:** 전체 테스트 파일
  - **AC:** SC-1 ~ SC-6
  - **Test:** 모든 테스트 통과, SC 체크리스트 전체 확인

**Definition of Done:** 모든 Must-Have AC가 충족되고, 기존 파이프라인이 회귀 없이 동작하며, 최종 결과가 baseline과 비교 문서화되어 프로젝트가 완결된다.

---

### 8. 파이프라인 로깅 및 관찰성 개선 (F-4 Logging) 🔴

- [x] **8.1** — `organize_outputs()` signature mismatch 수정 (F-4.4)
  - **What:** run_production.py에서 `repo_root=REPO_ROOT` kwarg 전달하나 output_organizer.py가 미수용 → TypeError. 호출부 또는 시그니처 수정.
  - **Files:** `run_production.py`, `egfr_pipeline/output_organizer.py`
  - **Test:** organize_outputs 호출 시 TypeError 미발생

- [x] **8.2** — Derived step view 파일 permission 수정 (F-4.7)
  - **What:** step_view.py의 `_atomic_write_json/text()`가 NamedTemporaryFile 기본 0600 permission 사용 → 0644로 수정. _staged_step_dir 디렉토리도 0755.
  - **Files:** `egfr_pipeline/step_view.py`
  - **Test:** 생성 파일 permission 0644, 디렉토리 0755

- [x] **8.3** — Phase skip 조건에 내용 검증 추가 (F-4.1)
  - **What:** `_canonical_output_valid()` 함수 추가: 파일 존재 + 최소 크기(100B) + CSV header+data 확인. Phase 4~6 skip 조건에 적용.
  - **Files:** `run_production.py`
  - **Test:** 빈 파일 → skip 안됨, 정상 파일 → skip됨

- [x] **8.4** — 후처리 실패 시 exit code 반영 (F-4.3)
  - **What:** organize_outputs/step view 실패를 `post_run_warnings` 리스트에 수집, 최종 요약 후 `sys.exit(1)`.
  - **Files:** `run_production.py`
  - **Test:** 후처리 warning 존재 시 exit code 1

- [x] **8.5** — 완료 요약에 derived step health 표시 (F-4.2)
  - **What:** `_print_completion_summary()`에 step view generated/not_generated 카운트 추가.
  - **Files:** `run_production.py`
  - **Test:** 요약 출력에 "step view" 문자열 포함

- [x] **8.6** — Phase 완료 시 핵심 output sanity check 로깅 (F-4.5)
  - **What:** Phase 4~7 완료 후 핵심 output 파일 존재/크기 확인 → [✓]/[△]/[✗] 터미널 출력.
  - **Files:** `run_production.py`
  - **Test:** sanity check 함수 존재, 빈 파일에 [△] 출력

- [x] **8.T** — Tests for 로깅 및 관찰성
  - **Files:** 테스트 파일
  - **Test:** pytest 전체 통과

---

### G. Dimer 입력 인식 및 전처리 수정 🔴

> **출처**: Production 결과 PyMOL/ChimeraX 검증에서 발견
> **긴급도**: 높음 — 현재 파이프라인이 dimer 입력을 monomer로 오인하고 있음

**Dependency:**
```
G-1 (진단) → G-5 (PPI_TARGETS 확인) → G-2 (메타데이터 정정) → G-6 (문서 정정)
                                         ↓
                               G-3 (전처리 로직) → G-4 (excluded residues)
                                         ↓
                                    G-7 (자동 감지)
```

- [x] **G-1** — 입력 PDB 구조 진단 (dimer 여부 확인)
  - **What:** 입력 PDB 3종(`3GT8_raw.pdb`, `EGFR_160-185.pdb`, `EGFR_170-200.pdb`)의 chain 구성, 잔기 번호 범위, atom 수를 확인한다. Phase 1 전처리 후 도킹 입력 PDB와 실제 도킹 결과 PDB의 chain 구성도 확인하여, 각 단계에서 dimer/monomer 상태를 파악한다.
  - **Files:** `input/receptors/*.pdb`, `output/phase1_ppi/runtime_inputs/docking_*.pdb`, `phase2_ppi_docking/*/prod_seed0/docking_*/cluster_results/*.pdb`
  - **AC:** 진단
  - **Test:** 각 입력 PDB별 chain ID, 잔기 범위, CA atom 수 기록 완료. Dimer/monomer 판정 근거 문서화.

- [x] **G-2** — `prepare_inputs.py` RECEPTOR_STATES 메타데이터 정정
  - **What:** G-1 진단 결과에 따라 RECEPTOR_STATES dict의 `source_chain`, `description` 필드를 실제 구조에 맞게 수정한다. `is_dimer` 필드를 추가하여 dimer 여부를 명시한다. 3GT8_raw도 PDB 3GT8 자체가 asymmetric dimer이므로 확인 후 정정.
  - **Files:** `egfr_pipeline/phase1/prepare_inputs.py`
  - **AC:** 코드 수정
  - **Test:** RECEPTOR_STATES의 description이 실제 구조와 일치. is_dimer 필드가 진단 결과와 일치.

- [x] **G-3** — Dimer-aware 전처리 로직 구현
  - **What:** G-1 진단 결과에 따라 `prepare_receptor()` 함수가 dimer 구조를 올바르게 처리하도록 수정한다. 핵심 결정: dimer를 도킹 입력으로 사용하는 것이 의도된 설계라면 메타데이터만 수정, 아니라면 전처리 로직 수정. Chain X 내 두 monomer 분리 또는 dimer 유지 전략 결정.
  - **Files:** `egfr_pipeline/phase1/prepare_inputs.py`
  - **AC:** 코드 수정
  - **Test:** 전처리 후 PDB가 의도된 구조(dimer 또는 monomer)와 일치. 기존 테스트 회귀 없음.

- [x] **G-4** — `generate_configs.py` excluded_residues 업데이트
  - **What:** Dimer 구조라면 두 번째 monomer의 membrane-proximal 잔기도 excluded_residues에 포함해야 하는지 검토한다. `prepare_dimer_pdb.py`의 `MEMBRANE_PROXIMAL_B`가 이미 이를 처리하는지 확인. 필요시 EXCLUDED_RESIDUES_A 업데이트.
  - **Files:** `egfr_pipeline/phase1/generate_configs.py`, `egfr_pipeline/ppi/prepare_dimer_pdb.py`
  - **AC:** 코드 수정
  - **Test:** Excluded residues가 dimer 구조에서 과학적으로 적절한 범위를 커버.

- [x] **G-5** — `run_production.py` PPI_TARGETS 입력 경로 검증
  - **What:** PPI_TARGETS 리스트가 실제로 어떤 config/input PDB를 가리키는지 확인하고, 해당 입력이 dimer/monomer 어느 것인지 명시적으로 기록한다.
  - **Files:** `run_production.py`
  - **AC:** 코드 수정
  - **Test:** PPI_TARGETS의 각 입력 경로가 올바른 PDB를 참조하며, dimer/monomer 상태가 주석으로 명시됨.

- [x] **G-6** — Report/메타데이터 monomer→dimer 용어 정정
  - **What:** "single chain", "monomer" 등 잘못된 표기를 dimer 반영으로 수정한다. 대상: `prepare_inputs.py`, `generate_configs.py`, `phase1_pyrosetta_execution_note.md`, `phase1_input_validation_report.md`, 보고서 등.
  - **Files:** `egfr_pipeline/phase1/prepare_inputs.py`, `egfr_pipeline/phase1/generate_configs.py`, `docs/phase1_pyrosetta_execution_note.md`, `docs/phase1_input_validation_report.md`
  - **AC:** 문서 수정
  - **Test:** grep으로 "single chain" 또는 "monomer receptor" 등 잘못된 표기가 남아있지 않음.

- [x] **G-7** — Phase 1 input validation에 dimer 감지 체크 추가
  - **What:** 입력 PDB를 읽을 때 자동으로 dimer 여부를 감지하는 `detect_dimer()` 함수를 구현한다. 판정 기준: (1) 같은 chain에서 동일 잔기 번호 2번 이상 등장, (2) CA atom 수가 kinase domain 1개(~309) 대비 비정상적으로 많음(>400). Phase 1 validation report에 결과 포함.
  - **Files:** `egfr_pipeline/phase1/prepare_inputs.py` 또는 `egfr_pipeline/validate.py`
  - **AC:** 코드 수정
  - **Test:** monomer PDB → `is_dimer: False`, dimer PDB → `is_dimer: True`. 기존 테스트 회귀 없음.

- [x] **G-T** — Tests for Dimer 수정
  - **What:** G-1~G-7 전체에 대한 테스트. detect_dimer 단위 테스트, RECEPTOR_STATES 메타데이터 일관성 테스트, excluded_residues 범위 테스트.
  - **Files:** 테스트 파일
  - **Test:** pytest 전체 통과

**Definition of Done:** 입력 PDB의 dimer/monomer 상태가 정확히 진단되고, 코드/문서/메타데이터가 실제 구조와 일치하며, 향후 입력에 대해 자동 감지가 동작한다.

**Definition of Done:** 빈 파일 skip 방지, 후처리 실패 시 exit code 반영, 완료 요약에 step health 표시, permission 수정으로 파이프라인 관찰성이 확보된다.

---

### H. 상태 오염 방지 — Lane Marker 검증 및 파이프라인 강건성 🔴

> **출처**: Production 실행 시 빈 lane marker가 전체 Workflow B를 무력화한 사건
> **긴급도**: 높음 — 테스트↔프로덕션 상태 오염이 반복 발생
>
> **이미 수정된 항목 (Group 8, commit 7e2661c + 87711bb):**
> - NamedTemporaryFile permission → 0644, 디렉토리 → 0755
> - `organize_outputs()` signature mismatch 해결
> - Phase 6 빈 파일/placeholder skip 방지 (`check_phase6()` placeholder 감지)
> - WARN → exit code 반영 (`sys.exit(2)` for post_run_warnings)
> - Phase 4-7 `_sanity_check_phase_outputs()` 추가

**근본 원인:** `_run_lane()`이 phase 함수 반환 후 핸드오프 파일 존재를 확인하지 않고 `lane_complete.json`을 생성. 이전 테스트 실행에서 빈 마커가 남아 실제 실행을 skip.

**기존 코드 흐름 (수정 전):**
```
_run_lane():
  L1612: lane_is_complete() → 마커만 확인, 산출물 안 봄    ← H-2 수정 지점
  L1622: reset_lane_completion()
  L1643: write_lane_manifest(status="running")
  L1645: phase 함수 실행
    └→ _validate_adv_handoff() ← INPUT 검증 (이미 존재, 변경 불필요)
    └→ 실제 phase 로직
  L1694: write_lane_completion(status="completed")         ← H-1 수정 지점
         (현재: 산출물 검증 없이 바로 "completed" 마커 생성)
```

**검증 레이어 구분:**
- `_validate_adv_handoff()`: **입력(INPUT)** 검증 — "이 lane을 실행하려면 선행 lane의 산출물이 있어야 한다" (이미 구현됨)
- `_validate_lane_outputs()`: **출력(OUTPUT)** 검증 — "이 lane이 완료됐으면 자기 산출물이 있어야 한다" (신규)

**Dependency:**
```
H-6 (PHASE_CORE_OUTPUTS 확장)
  → H-4 (_validate_lane_outputs, 모든 lane type)
    → H-1 (marker 쓰기 전 output 검증)
    → H-2 (skip 시 stale marker 감지)
H-5 (seed marker 무결성) — 독립
H-3 (PBS 로그) — 독립
H-T (통합 테스트)
```

- [x] **H-6** — `PHASE_CORE_OUTPUTS`를 Phase 1-3으로 확장 + `_sanity_check_phase_outputs()` 범위 확대
  - **What:** 현재 `PHASE_CORE_OUTPUTS`(L311-326)가 Phase 4-7만 정의하고, `_sanity_check_phase_outputs()`(L329-354)도 `phase_dirs`에 Phase 4-7만 있어 Phase 1-3은 early return된다. Phase 1(vina-cpu: `*_vina_results/` 내 `.pdbqt` 파일), Phase 2(ppi: `final_ranking.csv` per seed), Phase 3(ppi-post: `ppi_pyrosetta_residues.csv`, `ppi_pyrosetta_summary.csv`)의 핵심 산출물을 `PHASE_CORE_OUTPUTS`에 추가하고, `phase_dirs`에 Phase 1-3 경로를 추가한다. 검증 기준: 파일 존재 + 최소 크기(100B) + CSV는 header+1행 이상.
  - **Files:** `run_production.py` (`PHASE_CORE_OUTPUTS` dict, `_sanity_check_phase_outputs()`)
  - **AC:** 코드 수정
  - **Test:** Phase 1 완료 후 pdbqt 없으면 `[△]` 경고. Phase 2 ranking 비어있으면 `[✗]`. Phase 3 CSV 없으면 `[✗]`.

- [x] **H-4** — `_validate_lane_outputs()` 함수 및 전체 lane별 산출물 매핑 테이블 정의
  - **What:** 모든 lane type에 대해 기대 산출물 파일 목록을 dict으로 정의한다. H-6에서 확장된 `PHASE_CORE_OUTPUTS`를 기반으로 lane→phase 매핑을 구성한다. `_validate_adv_handoff()`가 "다음 lane의 입력" 관점이라면, 이 함수는 "현재 lane의 출력" 관점으로 검증한다. 매핑이 정의되지 않은 lane은 검증을 건너뛴다(기존 동작 유지).
  - **대상 lane:**
    - `vina-cpu` → Phase 1 산출물 (pdbqt 파일)
    - `ppi` → Phase 2 산출물 (final_ranking.csv, state/seed별)
    - `ppi-post` → Phase 3 산출물 (ppi_pyrosetta_*.csv)
    - `vina-post` → Phase 4 산출물 (이미 PHASE_CORE_OUTPUTS에 있음)
    - `finalize` → Phase 5-7 산출물 (이미 있음)
    - `adv-phase1` → `phase1_downstream_patch_reference.csv`
    - `adv-phase2` → `phase3_candidate_pocket_reference.csv`
    - `adv-phase3-setup` → `phase3_docking_job_table.csv`
    - `adv-phase3-execute` → `phase3_round_log.csv` (append)
    - `adv-phase3-post` → `phase4_docking_evidence_reference.csv`
    - `adv-phase4` → Phase 4 scoring 산출물
  - **Files:** `run_production.py`
  - **AC:** 코드 수정
  - **Test:** 모든 lane type에 대해 산출물 목록이 정의되어 있고, 파일 존재+비어있지 않으면 `(True, [])`, 미존재 또는 빈 파일이면 `(False, [missing_files])` 반환. 매핑 미정의 lane은 `(True, [])` 반환.

- [x] **H-1** — `_run_lane()` 완료 마커 생성 전 출력(OUTPUT) 검증 추가
  - **What:** `_run_lane()`에서 `write_lane_completion(status="completed")` 호출 직전(L1694)에, `_validate_lane_outputs(lane)`를 호출하여 해당 lane의 핵심 산출물이 실제로 존재하고 비어있지 않은지 검증한다. 검증 실패 시 status="failed"로 마커를 쓰고 `RuntimeError("Lane {lane} completed without valid outputs: {missing}")` raise. 이 검증은 `_validate_adv_handoff()`(입력 검증)와는 별개의 **출력 검증**이며, 양쪽이 같은 파일 목록을 양방향으로 체크하여 빈 마커 문제를 이중 방어한다. 수정 전: phase 함수 반환 → 즉시 "completed" 마커. 수정 후: phase 함수 반환 → 산출물 검증 → 통과 시에만 "completed" 마커.
  - **Files:** `run_production.py` (`_run_lane()` 함수, L1694 부근)
  - **AC:** 코드 수정
  - **Test:** phase 함수가 정상 반환했지만 핸드오프 파일이 없으면 status="failed" 마커 + exception. 핸드오프 파일이 있으면 정상 "completed" 마커.

- [x] **H-2** — Lane skip 시 stale marker 감지
  - **What:** `_run_lane()`의 skip 분기(L1612-1620)에서, `lane_is_complete()`가 True를 반환하더라도 `_validate_lane_outputs(lane)`로 실제 산출물 존재를 재확인한다. 마커만 있고 산출물이 없으면 stale marker로 판정 → `reset_lane_completion()` 호출 → `[WARN] Stale marker detected for {lane}: marker exists but outputs missing. Re-executing.` 출력 → skip하지 않고 lane 재실행으로 fall through. 매핑 미정의 lane은 기존 동작 유지(마커만으로 skip). 수정 전: 마커 status=="completed"이면 무조건 skip. 수정 후: 마커 + 산출물 둘 다 있어야 skip.
  - **Files:** `run_production.py` (`_run_lane()` 함수, L1612 부근)
  - **AC:** 코드 수정
  - **Test:** (1) 빈 마커 + 산출물 없음 → skip 안 하고 재실행. (2) 정상 마커 + 산출물 있음 → skip. (3) 매핑 미정의 lane → 마커만으로 skip(기존 동작).

- [x] **H-5** — `seed_complete.json` 무결성 검증 강화
  - **What:** `_seed_is_complete()`(L403-423)의 legacy migration 경로에 검증을 추가한다. 현재: `final_ranking.csv` 존재하면 즉시 마커 자동 생성. 수정 후: (1) ranking.csv가 최소 100B 이상이고 header+1행 이상인지 체크, (2) 검증 실패 시 마커 생성하지 않고 False 반환, (3) 자동 생성 시 marker에 `"migration": true, "ranking_rows": N, "ranking_size_bytes": N` 메타데이터 포함. 기존 프로덕션 실행에는 이미 marker가 있으므로 migration 경로를 타지 않아 하위호환성 유지.
  - **Files:** `run_production.py` (`_seed_is_complete()`, `_write_seed_completion()`)
  - **AC:** 코드 수정
  - **Test:** 빈 ranking.csv → seed 미완료. header만 있는 ranking.csv → seed 미완료. 정상 ranking.csv(header+1행 이상) → seed 완료 + migration 메타데이터 포함.

- [x] **H-3** — PBS 의존성 체인 실패 시 job ID 기록 및 로그 보존
  - **What:** `run_advanced_pipeline.pbs`에서 제출된 모든 job ID를 `output/workflow_a/logs/pbs_job_chain.log`에 기록한다. 형식: `YYYY-MM-DD HH:MM:SS | lane_name | job_id | depends_on`. 기존 로그가 있으면 append. `log_job()` 헬퍼 함수를 스크립트 내에 정의하여 각 qsub 직후 호출. 후속 디버깅 시 `qstat -f <job_id>` 또는 `tracejob <job_id>`로 실패 원인 추적 가능.
  - **Files:** `config/run_advanced_pipeline.pbs`
  - **AC:** 인프라 수정
  - **Test:** 파이프라인 제출 후 `pbs_job_chain.log`가 생성되고, 각 lane의 job ID/시각/의존이 기록됨.

- [x] **H-T** — Tests for 상태 오염 방지
  - **What:** 아래 시나리오를 테스트:
    1. `_validate_lane_outputs()` — 모든 lane type에 대해 산출물 매핑 정확성
    2. `_validate_lane_outputs()` — 빈 파일 → 검증 실패, 정상 파일 → 검증 성공
    3. `_validate_lane_outputs()` — 매핑 미정의 lane → True 반환 (기존 동작 유지)
    4. stale marker 시나리오: 마커 있고 산출물 없음 → skip 안 하고 재실행
    5. output 검증 시나리오: phase 성공했지만 산출물 없음 → "failed" 마커 + exception
    6. seed_complete.json — 빈 ranking → 미완료, header-only → 미완료, 정상 → 완료+migration
    7. PHASE_CORE_OUTPUTS — Phase 1-3 항목 존재 확인
    8. (regression) Phase 6 빈 report → skip 안 함 (Group 8에서 수정)
    9. (regression) 후처리 실패 → exit code ≠ 0 (Group 8에서 수정)
    10. (regression) step_view 파일 permission 0644 (Group 8에서 수정)
  - **Files:** `tests/test_state_contamination_group_h.py`
  - **Test:** pytest 전체 통과

**Definition of Done:** 마커가 산출물 존재를 보증하고(H-1), stale marker가 자동 감지되며(H-2), seed marker 무결성이 검증되고(H-5), Phase 1-7 전체에 sanity check가 있으며(H-6), 전체 lane type의 산출물 매핑이 정의되고(H-4), PBS 체인 실패 시 job ID 추적이 가능하다(H-3).
