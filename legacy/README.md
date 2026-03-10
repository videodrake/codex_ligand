# Legacy Scripts — 구버전 파일 매핑

> 이 디렉토리의 스크립트들은 참조용으로 보관되어 있습니다.
> 모든 기능은 `main.py` 통합 CLI를 통해 접근하세요.

## 파일 매핑 (구버전 -> 신버전)

### Vina 도킹 파이프라인

| 구버전 (legacy/) | 신버전 | 설명 |
|-------------------|--------|------|
| `run_docking.py` | `egfr_pipeline/vina/dock.py` (코어) + `main.py` (CLI 진입점) | Vina 도킹 실행 |
| `parse_vina_results.py` | `egfr_pipeline/vina/parse_poses.py` | 도킹 결과 파싱 |
| `extract_contacts.py` | `egfr_pipeline/vina/contacts.py` | 접촉 잔기 추출 |
| `cluster_pockets.py` | `egfr_pipeline/vina/cluster.py` | 포켓 클러스터링 |
| `summarize_pockets.py` | `egfr_pipeline/vina/summarize.py` | 포켓 요약 |
| `compare_pockets.py` | `egfr_pipeline/vina/compare.py` | 교차 receptor 비교 |

### PPI 도킹 파이프라인

| 구버전 (legacy/) | 신버전 | 설명 |
|-------------------|--------|------|
| `pipeline_manager.py` | `egfr_pipeline/pyrosetta_docking/pipeline_manager.py` | PyRosetta PPI 오케스트레이터 |
| `docking.py` | `egfr_pipeline/pyrosetta_docking/docking.py` | PyRosetta 워커 (Relax, Docking, Refinement) |
| `analysis.py` | `egfr_pipeline/pyrosetta_docking/analysis.py` | 스코어링, RMSD, Interface 분석 |
| `common.py` | `egfr_pipeline/pyrosetta_docking/common.py` | PyRosetta 유틸리티 |
| `extract_ppi_residues.py` | `egfr_pipeline/ppi/pyrosetta_extract.py` | PPI 잔기 추출 표준화 |

### 보고서 및 검증

| 구버전 (legacy/) | 신버전 | 설명 |
|-------------------|--------|------|
| `generate_report.py` | `egfr_pipeline/report.py` | 종합 보고서 생성 |
| `validate_outputs.py` | `egfr_pipeline/validate.py` | 출력 검증 |

### Config 및 PBS 파일

| 구버전 (legacy/) | 신버전 | 설명 |
|-------------------|--------|------|
| `config.ini` | `config/ppi_test_*.ini` | PyRosetta PPI 설정 (1K 모델, 테스트) |
| `config_10k.ini` | - | 10K 모델 설정 (사용 중단) |
| `config_20k.ini` | `config/ppi_prod_*.ini` | PyRosetta PPI 설정 (프로덕션) |
| `config_100k.ini` | - | 100K 모델 설정 (사용 중단) |
| `config_1M.ini` | - | 1M 모델 설정 (사용 중단) |
| `run_v1.pbs` | `config/run_ppi_prod.pbs` | PBS 프로덕션 스크립트 |
| `run_v2_test.pbs` | `config/run_ppi_test.pbs` | PBS 테스트 스크립트 |

### 유틸리티/호환성 점검 스크립트

| 구버전 (legacy/) | 신버전 | 설명 |
|-------------------|--------|------|
| `check_filter_v2_compat.py` | - | v2.0 필터 호환성 점검 (배포 전 1회성) |
| `check_improvements_compat.py` | - | v1.0 개선사항 호환성 점검 |
| `check_clustering_compat.py` | - | 클러스터링 호환성 점검 |
| `check_clustering_compat_v2.py` | - | 클러스터링 v2 호환성 점검 |
| `check_availability.py` | - | 의존성 확인 |
| `check_syntax.py` | - | 문법 확인 |
| `analyze.py` | - | 단독 분석 스크립트 |
| `analyze_ligand_contacts.py` | - | 리간드 접촉 분석 |

## 신규 모듈 (legacy에 대응 없음)

| 신버전 | 설명 |
|--------|------|
| `main.py` | 통합 CLI 진입점 (인터랙티브 메뉴 + 서브커맨드) |
| `egfr_pipeline/verdict.py` | 사이트 판정 (STRONG/MODERATE/WEAK) |
| `egfr_pipeline/ppi/postprocess_ppi.py` | PPI 후처리 자동화 (chain 원복 + 잔기 추출) |
| `egfr_pipeline/ppi/prepare_dimer_pdb.py` | Dimer PDB 준비 + chain 원복 |
| `egfr_pipeline/ppi/afm_extract.py` | AlphaFold-Multimer 파서 (stub) |
| `egfr_pipeline/config.py` | YAML config 로드 유틸리티 |
| `config/example-project.yaml` | 프로젝트 설정 (YAML) |

## 실행 방법 대응

```bash
# 구버전
python run_docking.py                          # 인터랙티브
python run_docking.py --mode blind             # CLI blind docking
python run_docking.py --config output/.../config.yaml  # config 재실행

# 신버전
python main.py                                 # 인터랙티브 (메뉴 [1] = Vina)
python main.py vina -c config/example-project.yaml     # CLI Vina docking
python main.py full -c config/example-project.yaml     # 전체 파이프라인
```
