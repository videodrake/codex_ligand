# Architecture Overview

> 현재 유효한 진입점, 데이터 흐름, 스키마 계약을 한곳에 정리한 문서.
> 최종 업데이트: 2026-03-10

---

## 유효 진입점

| 진입점 | 용도 | 실행 환경 |
|--------|------|----------|
| `python main.py` | 통합 CLI (인터랙티브 메뉴 + 서브커맨드) | 로컬/서버 |
| `python main.py vina` | Vina 리간드 도킹 | 로컬 (Vina 필요) |
| `python main.py postprocess` | Vina 후처리 체인 | 로컬 |
| `python main.py verdict` | 사이트 판정 (evidence classification) | 로컬 |
| `python main.py report` | 종합 보고서 생성 | 로컬 |
| `python main.py validate` | 출력 검증 | 로컬 |
| `python main.py full` | 전체 파이프라인 (vina→postprocess→verdict→report→validate) | 로컬 |
| `python main.py pyrosetta` | PyRosetta PPI 도킹 (서브메뉴) | HPC (PBS) |
| `python -m egfr_pipeline.vina.bootstrap config.yaml` | 부트스트랩 안정성 분석 (독립 CLI) | 로컬 |
| `python -m egfr_pipeline.vina.sweep config.yaml` | 커트오프 감도 스윕 (독립 CLI) | 로컬 |
| `qsub config/run_ppi_test.pbs` | PyRosetta PBS 배치 (테스트) | HPC |
| `qsub config/run_ppi_prod.pbs` | PyRosetta PBS 배치 (프로덕션) | HPC |
| `.venv/bin/pytest tests/ -v` | 테스트 스위트 (46 tests) | 로컬 |

### Deprecated 진입점 (실행 금지)

| 경로 | 대체 | 비고 |
|------|------|------|
| `legacy/run_docking.py` | `main.py vina` | v1.0 인터랙티브 |
| `legacy/pipeline_manager.py` | `main.py pyrosetta` | v1.0 PyRosetta |
| `legacy/*.py` (전부) | `main.py` 서브커맨드 | `legacy/README.md` 매핑 참조 |

---

## 데이터 흐름

```
                        ┌──────────────────────────────────────────┐
                        │            config (YAML/JSON)            │
                        └────────────────┬─────────────────────────┘
                                         │
                    ┌────────────────────┬┴──────────────────────┐
                    ▼                    ▼                       ▼
            ┌──────────────┐   ┌─────────────────┐   ┌──────────────────┐
            │  Vina Docking │   │ PyRosetta PPI   │   │ AlphaFold-Multi  │
            │  (dock.py)    │   │ (pipeline_mgr)  │   │ (afm_extract)    │
            └──────┬───────┘   └────────┬────────┘   └────────┬─────────┘
                   │                    │                      │
                   ▼                    ▼                      ▼
         vina_pose_table.csv    ppi_pyrosetta_         ppi_afm_
                   │            residues.csv            residues.csv
                   │            ppi_pyrosetta_               │
                   │            summary.csv                  │
                   ▼                    │                     │
    ┌──────────────────────┐            │                     │
    │  Postprocess Chain   │            │                     │
    │  parse → contacts    │            │                     │
    │  → cluster → summarize│           │                     │
    │  → compare           │            │                     │
    └──────────┬───────────┘            │                     │
               │                        │                     │
               ▼                        │                     │
    vina_pocket_table.csv               │                     │
    vina_drug_pocket_map.csv            │                     │
    vina_pocket_comparison.csv          │                     │
               │                        │                     │
               ├────────────────────────┤                     │
               │         Bootstrap (optional)                 │
               │         → vina_pocket_bootstrap.csv          │
               │                        │                     │
               ▼                        ▼                     ▼
    ┌──────────────────────────────────────────────────────────┐
    │                    Verdict Module                         │
    │  (Vina quality + PPI proximity + Cross-receptor)         │
    │  + Experimental priors (optional)                        │
    │  + Bootstrap stability (optional)                        │
    └──────────────────────┬───────────────────────────────────┘
                           │
                           ▼
              cross_method_agreement.csv
              valid_sites.csv
              vina_consensus_sites.csv
                           │
                           ▼
    ┌──────────────────────────────────┐
    │           Report                  │
    │  project_report.txt               │
    │  combined_residue_evidence.csv    │
    └──────────────────────────────────┘
                           │
                           ▼
    ┌──────────────────────────────────┐
    │          Validate                 │
    │  스키마 일치 / ID 일관성 /        │
    │  잔기 번호 / handoff 준비 확인    │
    └──────────────────────────────────┘
```

---

## 모듈 의존 관계

```
main.py (CLI/메뉴만 — 비즈니스 로직 없음)
  └─ egfr_pipeline/
       ├─ config.py          ← 모든 모듈이 사용
       ├─ residue_utils.py   ← compare, verdict가 사용
       │
       ├─ vina/
       │   ├─ dock.py        → (외부 Vina 실행)
       │   ├─ parse_poses.py → vina_pose_table.csv
       │   ├─ contacts.py    → vina_pose_table.csv (contact_residues 추가)
       │   ├─ cluster.py     → vina_pose_table.csv (pocket_id 추가)
       │   ├─ summarize.py   → vina_pocket_table.csv, vina_drug_pocket_map.csv
       │   ├─ compare.py     → vina_pocket_comparison.csv
       │   ├─ bootstrap.py   → vina_pocket_bootstrap.csv
       │   └─ sweep.py       (분석 도구, 출력 파일 없음)
       │
       ├─ ppi/
       │   ├─ pyrosetta_extract.py → ppi_pyrosetta_residues/summary.csv
       │   ├─ afm_extract.py       → ppi_afm_residues.csv (stub)
       │   ├─ postprocess_ppi.py   (chain 원복 자동화)
       │   ├─ prepare_dimer_pdb.py (dimer PDB 준비)
       │   └─ submit.py            (PBS qsub 제출)
       │
       ├─ verdict.py → cross_method_agreement.csv, valid_sites.csv,
       │               vina_consensus_sites.csv
       ├─ report.py  → project_report.txt, combined_residue_evidence.csv
       └─ validate.py (검증만, 출력 파일 없음)
```

---

## CSV 스키마 계약

모든 CSV 스키마는 `egfr_pipeline/validate.py`의 `EXPECTED_SCHEMAS` dict에 정의.
각 모듈의 `*_FIELDS` 상수와 1:1 대응하며, 테스트(`tests/test_pipeline.py::TestSchemaConsistency`)가 동기화를 보장.

### Core Outputs (반드시 존재)

| 파일 | 생성 모듈 | 컬럼 수 |
|------|----------|---------|
| `vina_pose_table.csv` | parse_poses → contacts → cluster | 13 |
| `vina_pocket_table.csv` | summarize | 14 |
| `vina_drug_pocket_map.csv` | summarize | 10 |

### Optional Outputs

| 파일 | 생성 모듈 | 컬럼 수 | 조건 |
|------|----------|---------|------|
| `vina_pocket_comparison.csv` | compare | 23 | 2+ receptors |
| `vina_pocket_bootstrap.csv` | bootstrap | 10 | bootstrap 실행 시 |
| `cross_method_agreement.csv` | verdict | 16 | verdict 실행 시 |
| `valid_sites.csv` | verdict | 23 | verdict 실행 시 |
| `vina_consensus_sites.csv` | verdict | 11 | cross-receptor 매치 존재 시 |
| `ppi_pyrosetta_residues.csv` | pyrosetta_extract | 10 | PPI 데이터 존재 시 |
| `combined_residue_evidence.csv` | report | 8 | report 실행 시 |

### 스키마 변경 절차

1. 모듈의 `*_FIELDS` 리스트 수정
2. `validate.py`의 `EXPECTED_SCHEMAS` 동일하게 수정
3. `tests/test_pipeline.py` 실행 → `TestSchemaConsistency` 자동 검증
4. 기존 smoke_test 출력 갱신 (필요 시)

---

## 설정 체계

두 독립 파이프라인이 각자의 config 형식 사용:

| 파이프라인 | 형식 | 로더 | 실행 환경 |
|-----------|------|------|----------|
| Vina 리간드 도킹 | YAML (또는 JSON) | `egfr_pipeline/config.py` | 로컬 |
| PyRosetta PPI | INI | `configparser` (pipeline_manager.py 내장) | HPC PBS |

공유 설정값 없음 — 통일 불필요. 상세: `config/README.md`.

---

## 테스트

```bash
.venv/bin/pytest tests/ -v              # 전체 46개 (0.5초)
.venv/bin/pytest tests/ -k bootstrap    # Bootstrap만
.venv/bin/pytest tests/ -k e2e          # End-to-end
.venv/bin/pytest tests/ -k schema       # 스키마 일치
```

테스트는 합성 fixture 사용 (실제 Vina/PyRosetta 불필요).
JSON config로 pyyaml 의존성 우회 (YAML 테스트는 별도 클래스).

---

## 핵심 불변 규칙

1. **main.py에 비즈니스 로직 금지** — import → call 패턴만
2. **CSV 스키마는 validate.py가 진실의 원천** — 모듈 FIELDS와 반드시 동기화
3. **Verdict 100점 총점 불변** — experimental priors는 정보 태그만 (점수 변경 X)
4. **PyRosetta 구버전 호환** — analysis.py의 try/except/hasattr 절대 제거 금지
5. **legacy/ 실행 금지** — 참조 전용, legacy/README.md에 매핑 문서
