# config/ — 설정 파일 가이드

## 어떤 config가 어떤 파이프라인용?

| 파일 | 파이프라인 | 형식 | 용도 |
|------|-----------|------|------|
| `example-project.yaml` | Vina 리간드 도킹 | YAML | 프로젝트 config 템플릿 |
| `ppi_test_TH1.ini` | PyRosetta PPI | INI | TH1 테스트 (1K 모델, ~3-4시간) |
| `ppi_test_beta_meander.ini` | PyRosetta PPI | INI | Beta-meander 테스트 (1K 모델) |
| `ppi_prod_TH1.ini` | PyRosetta PPI | INI | TH1 프로덕션 (20K 모델, ~24-36시간) |
| `ppi_prod_beta_meander.ini` | PyRosetta PPI | INI | Beta-meander 프로덕션 (20K 모델) |
| `run_ppi_test.pbs` | PyRosetta PPI | PBS | HPC 테스트 배치 스크립트 |
| `run_ppi_prod.pbs` | PyRosetta PPI | PBS | HPC 프로덕션 배치 스크립트 |

## 왜 형식이 다른가?

- **Vina (YAML)**: 최신 통합 파이프라인. `egfr_pipeline/config.py`가 로딩 (JSON도 가능).
- **PyRosetta (INI)**: v1.0부터 사용된 레거시 형식. `pipeline_manager.py`가 `configparser`로 직접 로딩.
- 두 파이프라인은 **공유 설정값이 없음** — 별도 실행 환경 (Vina=로컬, PyRosetta=HPC PBS).

## 사용법

```bash
# Vina — YAML config 지정
python main.py vina --config config/example-project.yaml
python main.py postprocess --config config/example-project.yaml

# PyRosetta — INI config 지정 (인터랙티브 메뉴에서 선택)
python main.py pyrosetta
# 또는 직접:
python -m egfr_pipeline.pyrosetta_docking.pipeline_manager config/ppi_test_TH1.ini input_PDB/C-lobe_TH1.pdb

# PBS 배치
qsub config/run_ppi_test.pbs
qsub -v CONFIG_FILE=config/ppi_prod_TH1.ini config/run_ppi_prod.pbs
```

## Config 섹션 참조

### Vina YAML (`example-project.yaml`)
- `receptors` / `ligands`: 입력 구조 목록
- `vina`: 도킹 파라미터 (exhaustiveness, n_poses, box)
- `postprocess`: 후처리 단계 on/off 및 cutoff
- `bootstrap`: 포켓 안정성 분석 파라미터
- `experimental`: 실험 데이터 (known binding residues)
- `ppi`: PyRosetta/AFM 결과 경로 (교차 검증용)

### PyRosetta INI (`ppi_*.ini`)
- `[Path]` / `[System]` / `[Docking]`: 입력 및 실행 환경
- `[FilterStage1]` / `[FilterStage2]`: v2.0 2-Pass 필터
- `[MiniRefinement]`: Stage 1→2 사이 인터페이스 리패킹
- `[Constraints]`: 제외 구역 / 핵심 잔기 (실험 데이터)
- `[ExperimentalData]`: 알려진 결합 영역 (Jura 2009)
