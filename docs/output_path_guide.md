# Output Path Guide

> `egfr_pipeline/paths.py`가 단일 진실 소스(single source of truth)이다.
> 이 문서는 경로 구조를 사람이 읽을 수 있도록 정리한 것이다.

## Workflow A — Standard Production

Vina blind docking + PPI blind docking → 통합 Verdict.

```
output/workflow_a/
├── phase1_vina_docking/          Vina blind docking 결과
│   └── {receptor_id}/            수용체별 .pdbqt 파일
├── phase2_ppi_docking/           PPI global docking (PyRosetta)
│   ├── {state}/prod_seed{n}/     상태×시드별 도킹 결과
│   └── runtime_inputs/           준비된 입력 파일
├── phase3_ppi_postprocess/       PPI 후처리 (잔기 추출, 스코어 표준화)
│   └── restored_runs/{receptor}/{partner}/
│                                postprocess 복원 산출물(restored PDB/CSV)
├── phase4_vina_postprocess/      Vina 후처리
│   ├── vina_pose_table.csv       전체 포즈 테이블
│   ├── vina_pocket_table.csv     포켓 요약 (centroid, affinity, depth)
│   ├── vina_drug_pocket_map.csv  리간드-포켓 매핑
│   ├── vina_pocket_comparison.csv  교차 수용체 비교
│   └── vina_pocket_bootstrap.csv   부트스트랩 안정성
├── phase5_verdict/               최종 판정
│   ├── valid_sites.csv           3축 스코어링 결과
│   ├── valid_sites_disclaimer.md 면책 조항
│   └── cross_method_agreement.csv
├── phase6_report/                종합 보고서
│   └── project_report.txt
├── phase7_validation/            출력 검증
└── logs/                         실행 로그
```

**진입점:** `run_production.py --lane {phase}`
**경로 함수:** `paths.wa_phase{N}_*(config)`

### Step-view overlay (optional)

Step-view를 활성화하면 요약/인덱스 폴더가 추가로 생성된다.

- 기본(root_mode 미지정): `output/workflow_a/step_views/{project_name}/step{N}_*`
- `step_output_view.root_mode = legacy_project_root`:
  `output/{project_name}/step{N}_*`
- 새 기본 모드에서 legacy 위치(`output/{project_name}`)가 이미 있으면
  `MOVED_TO_WORKFLOW_A_STEP_VIEWS.txt` 안내 파일이 남을 수 있다.

## Workflow B — PPI-First Advanced Pipeline

PPI 인터페이스 분석 → 포켓 제안 → Focused docking → Perturbation scoring.

```
output/workflow_b/
├── phase1_ppi_analysis/          PPI 인터페이스 매핑
│   ├── {state}/                  수용체 상태별
│   │   ├── orientation_filter_log.csv
│   │   ├── ppi_hotspot_residues.csv
│   │   ├── ppi_interface_patch_table.csv
│   │   └── lightdock/            LightDock 검증
│   ├── ppi_patch_cross_state_comparison.csv
│   ├── ppi_patch_state_robustness.csv
│   ├── phase1_downstream_patch_reference.csv
│   └── phase1_interface_report.md
├── phase2_pocket_analysis/       포켓 분석 및 분류
│   ├── candidate_pockets.csv
│   ├── pocket_patch_relationship.csv
│   └── phase3_candidate_pocket_reference.csv
├── phase3_focused_docking/       타겟 Vina docking
│   └── phase4_docking_evidence_reference.csv
└── phase4_scoring/               Perturbation relevance scoring
    ├── phase4_evidence_normalized.csv
    ├── phase4_axis_scores.csv
    ├── perturbation_candidate_table.csv
    ├── phase4_final_review_table.csv
    ├── phase4_review_disclaimer.md
    └── phase4_ranking_method_note.md
```

**진입점:** `run_production.py --lane adv-phase{N}`
**경로 함수:** `paths.wb_phase{N}_*(config)`

## Historical Reference Data (Non-Runtime)

`pilot_data_reference.csv` 같은 historical 참조 데이터는 **실행 입력(runtime input)이 아니다**.

- historical 항목은 `legacy_reference.*` 또는 `historical_reference.*` 네임스페이스 키로만 기록한다.
- 실행 파이프라인에서 쓰는 일반 키(`results_dir`, `runtime_inputs`)는 historical 레코드에서 직접 재사용하지 않는다.
  - 불가피하게 동일 키를 써야 한다면 `historical_only=true`를 반드시 포함해야 한다.
- downstream consumer는 다음 조건을 만족하는 항목만 historical로 읽어야 한다.
  - `historical_only == true AND status == historical_reference_only`

## Legacy PPI Fallback Compatibility Mode

Step-view Step 3(`step3_ppi_postprocess`)는 기본적으로 legacy 경로
`output/phase1_ppi/phase1_interface_report.md`를 읽지 않는다.

- 기본값: `compat.allow_legacy_ppi_paths = false`
- 명시적으로 `true`를 설정한 경우에만 legacy Phase 1 report fallback 허용
- legacy 경로를 실제로 사용하면 `step_manifest.json`의 `warning_codes`에
  `LEGACY_PPI_PATH_USED` 코드가 기록됨

```json
{
  "compat": {
    "allow_legacy_ppi_paths": true
  }
}
```

## Comparison Output

Workflow A↔B 비교 결과 (Group 5).

```
output/
├── workflow_comparison.csv        4분류 비교 (Consensus/A-only/B-only/Conflict)
├── workflow_comparison_report.md  비교 보고서
├── verdict_weight_sensitivity.csv 가중치 민감도 시뮬레이션
└── centroid_offset_correction.csv centroid 오프셋 보정 시뮬레이션
```

## Legacy (DEPRECATED)

이전 flat 레이아웃. 새 코드에서 사용 금지.

| 이전 경로 | 대체 | paths.py 함수 |
|-----------|------|---------------|
| `output/{project_name}/` | `output/workflow_a/` | `legacy_project_root()` → `workflow_a_root()` |
| `output/phase1_ppi/` | `output/workflow_b/phase1_ppi_analysis/` | `legacy_phase1_ppi_dir()` → `wb_phase1_ppi_analysis()` |

`paths.py`의 legacy 함수에는 `# DEPRECATED` 주석이 추가되어 있다. 기존 코드에서 이 함수를 호출하는 곳은 없다.
