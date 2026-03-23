# CONTEXT.md — 결정 배경 및 작업 기억

> 프로젝트 전체에 걸쳐 누적됩니다. 새 세션은 이 파일로 맥락을 복원합니다.

## 아키텍처 결정

| 결정 | 이유 | 대안 | 날짜 |
|------|------|------|------|
| CSV 확장 전략: 기존 컬럼 유지 + 신규 컬럼만 추가 | downstream 파서(PyMOL 스크립트, 보고서 생성기 등)가 기존 컬럼에 의존 | 새 스키마로 전면 교체 → 호환성 파괴 위험 | 2026-03-19 |
| 재도킹 회피: 기존 4,500 Vina + 300K PPI 모델 재활용 | 도킹 비용(CPU 수백 시간) 대비 분석만으로 목표 달성 가능 | 전체 재실행 → 비용 과다, fragment 범위 변경 시에만 부분 재실행 | 2026-03-19 |
| region_definitions.py 중앙 모듈: 모든 영역/잔기 상수를 한 곳에서 관리 | ATP site, Ko et al., SASA 영역 분류가 여러 모듈에서 참조됨 | 각 모듈에 하드코딩 → 불일치 위험 | 2026-03-19 |
| 핸드오프 검증 강화: 스키마 + 데이터 품질 2단계 | 빈 CSV나 hotspot 0개가 조용히 통과하여 의미 없는 결과 생성 방지 | 현재 스키마 체크만 유지 → edge case 취약 | 2026-03-19 |
| 순차 의존 실행: A → B,C(병렬) → D → E → F | 리간드 다양성(A)과 편향 크기(B)가 메트릭 조정(D)의 근거 | 모든 계획안 병렬 → 메트릭 조정 근거 부재 | 2026-03-19 |
| 바이브코딩 11건 필수: 기존 코드 분석 후 판단하여 구현 | 기존 로직의 암묵적 계약(함수 시그니처, 데이터 형태)을 파악해야 올바른 통합 가능 | 문서만 보고 구현 → 기존 코드와 충돌 위험 | 2026-03-19 |

## 알려진 리스크

| 리스크 | 영향 | 완화책 | 상태 |
|--------|------|--------|------|
| PPI 오류 전파: Phase 1 hotspot이 잘못되면 Workflow B 전체 왜곡 | Workflow B 결과 신뢰 불가 | AC-1.3: Ko et al. 일관성 체크로 조기 경보 | open |
| 리간드 유사성: 3종이 화학적으로 유사하면 consensus 무의미 | Vina consensus 메트릭 가치 저하 | AC-1.1(Task 1.1)에서 즉시 확인, 결과에 따라 가중치 조정 | open |
| Fragment 범위 변경: pilot에서 955-1006이 부적합하면 PPI 재실행 필요 | 서버 자원 수백 CPU시간, 일정 지연 | AC-3.3: 소규모 pilot(2K 모델)으로 사전 확인 | open |
| Verdict 가중치 변경: 프로덕션 재실행 필요 | Phase 5(Verdict)부터 재실행, valid_sites.csv 재생성 | `--from 5` 옵션으로 최소 범위 재실행 | open |
| Bootstrap 연동의 기존 결과 변동: 연동 후 STRONG/MODERATE 판정 뒤집힘 | 기존 결과와 비교 불가, 해석 혼란 | AC-2.3: 연동 전후 비교, 변화 3개 이상이면 임계값 재조정 | open |
| valid_sites.csv 주석 행: # 주석이 downstream 파서를 깨뜨릴 가능성 | 면책 조항 추가 실패 | EC-5.3: 충돌 시 별도 metadata 파일로 분리 | open |

## 그룹 노트

### Group 0 — Project Setup (2026-03-19)

**계획:**
- 0.1: 12개 수정 대상 파일 존재·구조 확인 (탐색 완료)
- 0.2: SDF 3종 존재 확인, ATP_SITE 37개 확정, Ko sheet 코드값 확인
- 0.3: pytest 실행, baseline 포켓 수 N/A (output/ 전체 비어 있음)
- 0.4: egfr_pipeline/region_definitions.py 신규 작성 (B-1.1 미완, 문헌 근사치)
- 0.T: Setup 검증 테스트

**기술 결정:**
- ATP_SITE_RESIDUES 37개 구성: P-loop(718-723), β3(726), αC-helix(743-745,762), β4(766), β5 linker(777), hinge(788-796), post-hinge(797,800), catalytic(831-837), pre-DFG(844), A-loop(854-858)
- 309개 전체 분류: SASA 미계산 → 구조 문헌 기반 근사치, B-1.1 FreeSASA 분석 시 업데이트 필요
- Ko et al. sheet: orientation_filter.py 코드값 (sheet_8: 961-964, sheet_9: 968-972, sheet_10: 977-980, sheet_11: 985-988, sheet_12: 993-997)
- SDF RDKit 로드: Group 0에서는 파일 존재만 확인, RDKit 의존 테스트는 Task 1.1
- Baseline: 292 tests 전체 통과 (2026-03-19). valid_sites.csv 미존재 → 포켓 수 baseline N/A

### Group 1 — 실험 데이터 통합 F-1 (2026-03-19)

**순서:** 1.1 → 1.2 → 1.5 → 1.3 → 1.4 → 1.T
**의존:** 1.2와 1.5가 pocket_summary.py 수정 (순차). 나머지 독립적.

---

**수정 대상 파일 (탐색 결과):**

| 파일 | 줄 수 | 주요 함수/클래스 |
|------|-------|-----------------|
| verdict.py | 1808 | score_pocket(), generate_verdict(), _check_membrane_overlap() |
| report.py | 655 | generate_report(), format_verdict_section() |
| validate.py | 805 | run_validation(), check_residue_numbering() |
| paths.py | 250 | workflow_a/b 경로 함수들, legacy 함수들 |
| vina/pocket_summary.py | 227 | summarize_pose_rows(), summarize_from_config() |
| vina/pocket_stability.py | 406 | run_bootstrap_replicates(), compute_bootstrap_stats() |
| phase1/orientation_filter.py | 772 | process_state_orientation(), SHEET_*_RESIDUES |
| phase1/cluster_consensus.py | 820 | compute_cluster_consensus(), NLOBE_CLOBE_BOUNDARY |
| phase1/compare_states.py | 534 | compare_across_states() |
| phase1/review_report.py | 668 | generate_review_report() |
| phase1/lightdock_validation.py | 1314 | compute_cross_method_convergence() |
| pyrosetta_docking/pipeline_manager.py | 4156 | PipelineManager class |

### Group 2 — Vina Blind Docking 편향 정량화 F-2 (2026-03-20)

**수정된 실행 순서:** 2.1, 2.2, 2.4 (독립) → 2.3 (축소) → 2.5 → 2.T
**2.3 범위 축소:** Bootstrap-Verdict 연동 코어 로직이 이미 verdict.py에 구현됨(score_pocket 4단계 stability_pts, bootstrap CSV 로드·병합·감쇠 전부 동작) — bootstrap_confidence 카테고리 컬럼 추가 + 임계값 정합성 확인만 수행.

**기술 결정:**
- pose_region_classifier.py 별도 모듈: pocket_summary.py와 관심사 분리 (포즈 수준 영역 분석 vs 포켓 수준 요약)
- 영역 판정: 접촉 잔기 과반(>50%) 기준, 과반 미달 시 "mixed", 잔기 없으면 "unknown"
- C-lobe surface < 10% → WARNING (AC-2.1), 0% → EC-2.1 특별 경고
- 2.4는 2.1과 독립: is_atp_site 플래그(Group 1 구현)로 필터링, pose region 분류 불필요
- Affinity 임계값: 현재 -8.0/-6.5 2단계만 존재 (tasks.md의 -5.0은 코드에 없음)
- validate.py 통합 시 기존 parse_pdb_residue_identity() 재활용, exit code 0/1/2 준수

**완료 (2026-03-22):**
- 2.5: manual_vina.md 해석 가이드에 소수성 과대평가 편향, C-lobe surface affinity 지침, ATP site 배제 근거, Ko sheet 접촉 해석 추가
- 2.T: test_vina_bias_group2.py 39건 — AC-2.1~2.5 커버, EC-2.1(0% C-lobe), EC-2.2(≥10 mismatch FAIL), EC-2.3(all high bootstrap), 경계값(stability 0.40/0.60/0.80), 잘못된 입력, CHARMM HIS 정규화
- Group 2 전체 완료: 365 passed (ligand_diversity 2건 제외 — output CSV 미존재, Group 2 무관)

### Group 3 — PPI Branch 강건성 검증 F-3 (2026-03-22)

**실행 순서:** 3.1 → 3.2 → 3.4 (orientation_filter.py 집중) → 3.3 (독립) → 3.5 (run_production.py) → 3.6 (4파일) → 3.T

**기술 결정:**
- 3.1 threshold sweep: orientation_log.csv의 기존 dot product 재분류 + interface_models.csv에서 PASS 모델 잔기 집계 → PyRosetta 불필요
- 3.4 sheet 12 sensitivity: compute_orientation_score()에 active_face_residues 선택 인자 추가. 실행은 서버 필요(PCA 재계산)
- 3.5 핸드오프 가드: 기존 FileNotFoundError 패턴 유지 (EC-3.3)
- 3.6 concordance_score = min(occ_a, occ_b) / max(occ_a, occ_b). both → strong_both(>0.5) / weak_both(≤0.5)

**완료 (2026-03-23):**
- 3.T: test_ppi_group3.py 60건 — AC-3.1~3.9 전체 커버 + EC-3.1~3.3 + adversarial 20건 (malformed scores, boundary values, empty CSV, missing files, NaN handling)
- Group 3 전체 완료: 425 passed (ligand_diversity 2건 제외 — output CSV 미존재, Group 3 무관)
- 컬럼명 의도적 편차: PRD 축약명 대신 구현에서 더 정확한 명명 사용 (n_orientation_valid_models, pyrosetta_max_occupancy, lightdock_frequency, concordance_score). 기능은 PRD 의도 100% 일치.

### Group 4 — Verdict 메트릭 체계 검증 및 교정 F-4 (2026-03-23)

**실행 순서:** 4.5 (문서) → 4.1 → 4.3 (verdict 시뮬레이션) → 4.4 (거리 분포) → 4.2 (pocket_depth) → 4.T

**기술 결정:**
- output/ 디렉토리에 실제 파이프라인 결과 없음 (서버 미실행) → 합성 데이터 + score_pocket() 직접 호출로 what-if 시뮬레이션
- 4.5 먼저: phase4/ 코드 읽기 전용, 의존 없음. A3 축 정의 추출 + 문서화
- 4.1 가중치 시뮬레이션: score_pocket() 래핑하여 6개 조합 what-if. 실제 valid_sites.csv 없이도 동작
- 4.3 cross_receptor_pts 차등: 시뮬레이션 스크립트에서 what-if (verdict.py 직접 수정 아님)
- 4.2 pocket_depth: pocket_summary.py에 함수 추가, PDB 존재 시에만 실행. 오프셋 시뮬레이션은 합성 데이터
- Phase 4 A3 축: 이미 perturbation_scoring.py + score_framework.py에 완전 구현됨 (A1=0.30, A2=0.25, A3=0.30, A4=0.15)

**완료 (2026-03-23):**
- 4.5: phase4_A3_axis_specification.md — A3 계산로직 4항목 추출 (입력/임계값/가중치/출력범위)
- 4.1: verdict_weight_sensitivity.py — 6개 가중치 조합 시뮬레이션, 합성 데이터 + score_pocket 래핑
- 4.3: verdict.py cross_coverage 차등 (20/14/0), 시뮬레이션 3안 비교. 1/3 coverage=0 (PRD 10 대비 의도적 하향)
- 4.2: pocket_depth_A 컬럼 + compute_pocket_depth() + centroid_offset_analysis.py alpha 0.5-1.0
- 4.4: PPI 거리 분포 분석 — 4구간, 70% 집중 경고, percentile, 오프셋 실질 범위
- 4.T: test_verdict_group4.py 43건 — AC-4.1~4.5 + EC-4.1~4.3 + adversarial 12건
- Group 4 전체 완료: 468 passed (ligand_diversity 2건 제외)

### Group 5 — 워크플로우 비교 및 문서 체계 정비 F-5 (2026-03-23)

**실행 순서:** 5.1 → 5.3 → 5.2 → 5.5 → 5.4 → 5.6 → 5.T

**기술 결정:**
- 포켓 매칭: valid_sites.csv에 centroid 미포함 → vina_pocket_table.csv + candidate_pockets.csv에서 centroid 로딩. centroid 거리 < 8Å + 잔기 Jaccard ≥ 0.3
- CSV 면책: csv.DictReader가 # 주석 미지원 (21+ 곳 사용) → valid_sites_disclaimer.md 별도 파일로 분리 (EC-5.3)
- Allosteric 태깅: verdict.py에서 vina_quality_score ≥ 35 AND ppi_proximity_score ≤ 5 → allosteric_candidate = True. Phase 4의 allosteric_modulator_candidate와 병행
- 5.3을 5.2보다 먼저: allosteric_candidate 컬럼이 비교 모듈의 A-only 플래그에 필요
- legacy paths.py: 2개 함수 (legacy_project_root, legacy_phase1_ppi_dir) — 외부 호출 없음, DEPRECATED 주석만

**완료 (2026-03-24):**
- 5.1: workflow_comparison_design.md — 매칭 방법, 4분류, 제외 조건 설계
- 5.3: verdict.py allosteric_candidate (Vina≥35 AND PPI≤5) + report.py 섹션 4.7
- 5.2: workflow_comparison.py — centroid+Jaccard 매칭, 4분류, CSV+보고서
- 5.5: 면책 조항 — report.txt 최상단 DISCLAIMER + valid_sites_disclaimer.md + phase4_review_disclaimer.md (EC-5.3: 별도 파일)
- 5.4: methodology_limitations.md 5섹션 + workflow_comparison_guide.md 3시나리오 + CLAUDE.md/research_overview 참조
- 5.6: output_path_guide.md 4섹션 + paths.py DEPRECATED 주석 2개
- 5.T: test_workflow_group5.py 39건 — AC-5.1~5.6 + EC-5.1~5.3 + adversarial 9건
- Group 5 전체 완료: 507 passed (ligand_diversity 2건 제외)

### Group 6 — 인프라 및 장기 개선 F-6 (2026-03-24)

**실행 순서:** 6.1 → 6.3 → 6.2 → 6.T

**기술 결정:**
- 6.1: pipeline_manager.py (4156줄, 순환 참조 없음, 500줄+ 함수 없음), vina_executor.py (2792줄, 43함수, 순환 참조 없음) → "현재 구조 유지" 결정 예상
- 6.3: Phase 3 ThreadPoolExecutor 이미 구현 (run_diverse_docking.py line 304-308) → EC-6.3 문서화만
- 6.2: precheck에 fpocket/P2Rank/LightDock 체크 추가 + run_advanced_pipeline.pbs 연결

**완료 (2026-03-24):**
- 6.1: module_separation_analysis.md — "현재 구조 유지" 결정. pipeline_manager 500줄+ 2개(보고서), vina_executor 없음, 순환 참조 없음
- 6.3: Phase 3 ThreadPoolExecutor 확인 → EC-6.3 문서화만 (이미 구현)
- 6.2: run_pre_qsub_checks.sh에 CHECK_WORKFLOW_B 플래그 추가 (fpocket critical, P2Rank/LightDock warning) + run_advanced_pipeline.pbs 연결
- 6.T: test_infra_group6.py 17건
- Group 6 전체 완료: 524 passed

## 발견된 이슈

(범위 밖 버그, 기술 부채, 추후 개선 사항을 여기에 기록)
- [ ] vina_executor.py(2813줄)가 tech-stack에 [읽기]로만 기록됨 — Group 6에서 분리 검토 시 구조 파악 필요 — 발견: Stage 2 교차검증
- [ ] verdict.py `_MEMBRANE_RESIDUES`(line 968)에 잔기 747 누락 — generate_configs.py와 prepare_inputs.py에는 747 포함 — 발견: Group 0 탐색
- [ ] NLOBE_CLOBE_BOUNDARY=838이 7개 파일에 중복 정의 — 향후 region_definitions.py에서 중앙화 대상
- [ ] SHEET_DEFINITIONS가 prepare_inputs.py와 orientation_filter.py 두 곳에 중복 — 향후 region_definitions.py에서 중앙화 대상
