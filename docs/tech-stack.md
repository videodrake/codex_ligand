# Tech Stack: EGFR-MYO1D 파이프라인 Major Refactoring

## Language & Runtime

- **Python 3.x** — 기존 파이프라인 전체가 Python으로 구현되어 있음. 모든 수정 및 신규 모듈도 Python으로 작성
- **HPC 환경** — PBS (Portable Batch System) 기반 작업 스케줄링, qsub를 통한 도킹 작업 제출

## Framework & Libraries

| 카테고리 | 선택 | 이유 | 대안 검토 |
|----------|------|------|-----------|
| 화학정보학 | RDKit | 리간드 다양성 분석(Morgan fingerprint, Tanimoto, physicochemical descriptors). 이미 환경에 설치됨 | OpenBabel — fingerprint 지원이 제한적 |
| 소분자 도킹 | AutoDock Vina | 기존 blind/focused docking 결과 재활용. exhaustiveness 384, n_poses 500 설정 | Glide — 라이선스 비용, 기존 결과 호환 불가 |
| PPI 도킹 | PyRosetta (ref2015) | 기존 300K 모델 결과 재활용. InterfaceAnalyzerMover 기반 정밀 스코어링 | HADDOCK — 기존 파이프라인과 통합 불가 |
| PPI 검증 | LightDock (DFIRE2) | PyRosetta와 독립적 scoring function으로 cross-method validation | ClusPro — 서버 기반이라 자동화 어려움 |
| 표면적 계산 | FreeSASA | SASA 기반 영역 분류(n_lobe/atp_site/c_lobe_surface/c_lobe_core). B-1.1에서 사용 완료 | NACCESS — Python 바인딩 없음 |
| 구조 파싱 | BioPython (PDB) | PDB 파일에서 잔기 번호, 아미노산 타입, 좌표 추출 | MDAnalysis — 과도한 기능, 이 용도에 불필요 |
| 데이터 처리 | pandas, numpy | CSV 기반 데이터 파이프라인의 핵심. 기존 코드 전체에서 사용 | polars — 기존 코드 호환성 유지 위해 pandas 계속 사용 |
| 포켓 탐지 | fpocket, P2Rank | Workflow B Phase 2에서 포켓 제안. 기존 설치 전제 | SiteMap — 라이선스 비용 |

## External Tools & Services

| 도구 | 용도 | 필수 여부 |
|------|------|-----------|
| AutoDock Vina CLI | 소분자 도킹 실행 (blind/focused) | Workflow A/B 모두 필수 |
| PBS/qsub | HPC 작업 스케줄링, lane별 병렬 제출 | 프로덕션 실행 필수 |
| fpocket | Workflow B Phase 2 포켓 탐지 | Workflow B 필수 |
| P2Rank | Workflow B Phase 2 포켓 탐지 (fpocket과 보완) | Workflow B 필수 |
| PyMOL | 구조 시각화, 사람 검토 시 사용 | 자동화 불필요, 사람 개입 시만 |

## Project Structure

기존 파이프라인 코드에 대한 수정([수정])과 신규 파일([신규]) 목록. 기존 구조를 유지하면서 기능을 추가하는 방식이다.

```
egfr_pipeline/                                  [기존]
├── verdict.py                                  [수정] F-1,2,4,5: ATP 태깅, bootstrap 연동, 차등 점수, allosteric 분류
├── report.py                                   [수정] F-1,2,5: 배제 섹션, 분포 요약, allosteric 섹션, 면책 조항
├── validate.py                                 [수정] F-2: 구간별 잔기 매핑 비교 추가
├── region_definitions.py                       [신규] F-2: SASA 기반 영역 분류 모듈 (B-1.1 산출물 이동)
├── workflow_comparison.py                      [신규] F-5: Workflow A↔B 포켓 비교 모듈
├── paths.py                                    [수정] F-5: legacy 함수 deprecation 주석
├── vina/
│   ├── pocket_summary.py                       [수정] F-1,4: ATP site 판정, pocket_depth, sheet 접촉
│   ├── pocket_stability.py                     [수정] F-2: bootstrap 결과 pocket_table 병합
│   ├── vina_executor.py                        [읽기] F-6: 모듈 분리 검토 대상 (2813줄)
│   └── pose_region_classifier.py               [신규] F-2: 포즈 영역별 분류
├── phase1/
│   ├── orientation_filter.py                   [수정] F-3: threshold sweep, ambiguous 분석, sheet 12 설정 변이
│   ├── cluster_consensus.py                    [수정] F-3: n_valid_models 기록
│   ├── compare_states.py                       [수정] F-3: conformational_selection_candidate 태그
│   ├── review_report.py                        [수정] F-3,5: ambiguous 섹션, LightDock 독립성 한계
│   └── lightdock_validation.py                 [수정] F-3: 정량적 concordance 확장
├── phase2/                                     [수정] F-3: 분류 비율 모니터링 통계
├── phase3/                                     [수정] F-6: 라운드 내 병렬화 (선택)
├── phase4/                                     [수정] F-4: A3 축 로직 문서화
└── pyrosetta_docking/
    └── pipeline_manager.py                     [수정] F-5: rigid-body 한계 문서 (코드 내 주석)

scripts/
├── assess_ligand_diversity.py                  [신규] F-1: RDKit 리간드 다양성 분석
├── analyze_pose_distribution.py                [신규] F-2: 포즈 영역별 분포 분석
├── analyze_affinity_distribution.py            [신규] F-2: C-lobe affinity 분포 분석
├── verdict_weight_sensitivity.py               [신규] F-4: 가중치 민감도 시뮬레이션
├── centroid_offset_analysis.py                 [신규] F-4: centroid 오프셋 보정 시뮬레이션
└── pilot_fragment_range.py                     [신규] F-3: fragment 범위 pilot config/PBS 생성

run_production.py                               [수정] F-1,3: _validate_adv_handoff 강화 (Ko 체크, 품질 가드)

config/
├── run_pre_qsub_checks.pbs                     [수정] F-6: --check-workflow-b 플래그 추가
└── run_advanced_pipeline.pbs                   [수정] F-6: Phase 2 전 도구 자동 검증

docs/
├── prd.md                                      [신규] 이 PRD
├── tasks.md                                    [신규] 태스크 목록
├── tech-stack.md                               [신규] 이 문서
├── methodology_limitations.md                  [신규] F-5: 방법론적 한계 통합 (5개 섹션)
├── workflow_comparison_guide.md                [신규] F-5: 불일치 해석 가이드
├── output_path_guide.md                        [신규] F-5: 출력 경로 가이드
└── phase1_notes.md                             [수정] F-3: threshold/범위/sheet12 결정 근거 추가

output/ (실행 시 생성)
├── ligand_diversity_assessment.csv             [신규] F-1: A-1.2 산출물
├── vina_pose_distribution_by_region.csv        [신규] F-2: B-1.2 산출물
├── residue_alignment_check.csv                 [신규] F-2: B-2.1 산출물
├── orientation_threshold_sensitivity.csv       [신규] F-3: C-1.1 산출물
├── orientation_ambiguous_report.csv            [신규] F-3: C-2.1 산출물
├── sheet12_sensitivity.csv                     [신규] F-3: C-4.1 산출물
├── verdict_weight_sensitivity.csv              [신규] F-4: D-1.1 산출물
└── workflow_comparison/
    ├── workflow_comparison.csv                 [신규] F-5: E-1.2 산출물
    └── workflow_comparison_report.md           [신규] F-5: E-1.2 산출물
```

## Development Environment

- **테스트**: pytest 기반 단위 테스트 + 프로덕션 회귀 테스트 (기존 valid_sites.csv 비교)
- **린터/포매터**: 기존 프로젝트 설정 따름
- **버전 관리**: Git, 각 계획안(A~F)을 별도 브랜치 또는 커밋 그룹으로 관리
- **CI**: PBS 기반 프로덕션 실행이 CI 역할. `run_pre_qsub_checks.pbs`가 사전 환경 검증 수행
- **데이터 관리**: 프로덕션 결과(pose_table.csv, scored_all_models.csv 등)는 기존 output/ 디렉토리에 존재. 신규 분석 결과는 동일 디렉토리 내 별도 파일로 생성

## Key Constraints

- **후방 호환**: 기존 valid_sites.csv, pocket_table.csv의 기존 컬럼은 변경하지 않고 새 컬럼만 추가. 기존 파서가 깨지지 않아야 함
- **재도킹 최소화**: 4,500 Vina 포즈와 300K PPI 모델 재생성은 원칙적으로 하지 않음. Fragment 범위 변경(C-3.4)이나 Verdict 가중치 변경(D-1.3) 시에만 부분 재실행
- **사람 개입 의존**: 11건의 바이브코딩 필수 항목과 추가 사람 판단 항목이 있어, 완전 자동화는 불가능. 바이브코딩 결과물은 사람 검토를 거쳐야 프로덕션에 반영
- **서버 자원**: Fragment pilot(C-3.3)은 서버에서 3개 PBS 작업 제출 필요 (각 ~수 시간). 프로덕션 재실행 시에도 서버 자원 필요
- **확정된 상수 불변**: ATP site 37잔기, Ko et al. sheet 번호, SASA 5영역 분류는 이미 확정되어 변경하지 않음
- **프로덕션 재실행 시점**: Verdict 가중치/임계값 변경 시 `python run_production.py --from 5` (Phase 5부터). PPI 재실행은 fragment 범위 변경 시에만
