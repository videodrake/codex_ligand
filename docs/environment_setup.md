# 환경 설정 가이드

이 문서는 프로젝트를 처음 사용하는 사람을 위한 서버 환경 설치 및 설정 가이드이다.

목표:

- 프로젝트가 기대하는 환경을 파악
- 필수 도구와 선택 도구를 구분
- `PATH`에 있어야 하는 외부 바이너리를 확인
- 대규모 작업 제출 전에 설정을 검증

이 레포지토리는 프로덕션 환경을 자동으로 생성하거나 관리하지 **않는다**.
공유 서버 환경은 사용자가 직접 관리해야 한다.

## 1. 서버 환경 기본 사항

서버 환경의 공식 기준:

- 공유 conda 환경 이름: `pyrosetta`
- 프로덕션과 pre-qsub 검증이 동일한 환경 사용
- LightDock이 현재 Phase 1 보조 검증 경로로 활성화
- AlphaFold-Multimer는 일상 기준에 포함되지 않음
- 레포지토리 스크립트는 별도의 테스트 환경을 생성하지 않음

처음부터 시작하는 경우, **하나의** `pyrosetta` 환경을 준비하고
필요한 도구를 해당 환경에 설치하거나 서버 `PATH`에서 사용 가능하게 한다.

## 2. 기본 Python 패키지

pre-qsub 검증을 포함한 모든 작업에 필수인 패키지:

- `pytest`
- `pyyaml`
- `numpy`
- `pandas`
- `scipy`
- `matplotlib`

패키지 참조 목록은 [requirements-test.txt](../requirements-test.txt)에도 저장되어 있다.

이 패키지들은 다음 스크립트에서 검사된다:

- `scripts/setup_test_env.sh`
- `scripts/run_pre_qsub_checks.sh`
- `docs/pre_qsub_test_line.md`

### 수동 설치 방법

환경을 먼저 활성화한다:

```bash
source ~/.bashrc
conda activate pyrosetta
```

그런 다음 서버에서 사용 가능한 패키지 소스를 사용하여 수동으로 설치한다:

```bash
conda install pytest pyyaml numpy pandas scipy matplotlib
```

서버가 오프라인인 경우, 사이트의 내부 미러, 사전 로드된 채널, 또는
관리자 승인 패키지 경로를 사용한다.

## 3. 파이프라인 단계별 도구 설치

### 3.1 Profile A: Pre-qsub 검증 전용

레포지토리 상태 검증만 필요한 경우:

- 공유 `pyrosetta` conda 환경
- 위의 기본 Python 패키지

이 프로파일로 실행 가능한 것:

- `bash scripts/setup_test_env.sh`
- `bash scripts/run_pre_qsub_checks.sh`
- `qsub config/run_pre_qsub_checks.pbs`

전체 프로덕션 도킹에는 **불충분**하다.

### 3.2 Profile B: 현재 활성 기준

현재 활성 워크플로우를 실행하려면 다음을 설치하거나 제공한다.

#### A. PyRosetta

필요 대상:

- Phase 1 receptor-side PPI mapping
- `qsub config/run_lightdock.pbs`
- `qsub config/run_lightdock_test.pbs`
- `qsub config/run_production.pbs` (Phase 2/3 PPI 단계 포함 시)

참고:

- PyRosetta는 이 레포지토리에서 설치되지 않음
- `pyrosetta` 환경에 직접 설치해야 함
- 라이선스와 설치 방법은 로컬 연구실/서버 정책에 따름

#### B. AutoDock Vina Python API

필요 대상:

- `python main.py vina`
- `run_production.py` 내 Vina 도킹

필요 Python 패키지:

- `vina`

#### C. 리간드/수용체 PDBQT 준비 도구

Vina는 `.pdbqt` 파일을 사용한다. 레포지토리는 사전 빌드된 `.pdbqt` 파일을
재사용할 수 있지만, 누락 시 자동 준비를 시도한다.

`egfr_pipeline/vina/vina_executor.py`의 현재 폴백 순서:

수용체 준비:

- `prepare_receptor` (ADFR)
- `prepare_receptor4.py` via `pythonsh` (MGLTools)
- `obabel`

리간드 준비:

- `mk_prepare_ligand.py` (Meeko)
- `prepare_ligand` (ADFR)
- `prepare_ligand4.py` via `pythonsh` (MGLTools)
- `obabel`

실용적으로, 최소한 다음을 설치한다:

- 수용체 준비 도구 하나
- 리간드 준비 도구 하나

가장 안정적인 조합:

- Meeko + OpenBabel
- ADFR suite + OpenBabel

입력 `.pdbqt` 파일이 이미 준비되어 유효한 경우, 이 도구들은 일상 재실행에
덜 중요하다.

#### D. RDKit

용도:

- 리간드 파일 처리
- SDF 기반 워크플로우
- 일부 로컬 분자 준비 경로

패키지:

- `rdkit`

현재 레포지토리에서 강력히 권장되며, 일부 워크플로우는 없이도 동작할 수 있다.

#### E. LightDock

필요 대상:

- 활성 Phase 1 보조 검증
- `egfr_pipeline.phase1.lightdock_validation`
- `run_lightdock_<state>.sh` 생성 및 실행

레포지토리가 기대하는 LightDock 커맨드라인 도구:

- `lightdock3_setup.py`
- `lightdock3.py`
- `lgd_generate_conformations.py`
- `lgd_cluster_bsas.py`

`egfr_pipeline/phase1/lightdock_validation.py`에 문서화된 일반적인 설치 경로:

- `pip install lightdock3`
- `conda install -c bioconda lightdock`

중요한 실용 규칙:

- Python 패키지가 존재하는 것만으로는 불충분
- LightDock 커맨드라인 실행 파일이 `PATH`에서 호출 가능해야 함

### 3.3 Profile C: 확장 연구 도구

최소 pre-qsub 경로에는 불필요하지만, 확장 또는 후속 단계 작업에 필요하다.

#### A. 포켓 제안 도구

Phase 2 포켓 제안 모듈이 사용하는 도구:

- `fpocket`
- `P2Rank`
- 선택적으로 FTMap 등의 핫스팟 증거

현재 상태:

- 코드가 파싱 및 설정 생성을 지원
- 이 도구들은 기본 CLI가 아닌 서버 측에서 실행될 수 있음

전체 Phase 2 브랜치를 사용할 계획이라면 설치:

- `fpocket`
- `P2Rank`

#### B. MD 스택

MD는 초기 온보딩 대상이 아니며 현재 일상 프로덕션 경로에서 자동 실행되지 않지만,
광범위한 연구 파이프라인의 일부이다.

MD 관련 작업에 필요한 도구:

- `gromacs` / `gmx`
- `MDAnalysis`
- 선택적으로 `gmx_MMPBSA` 또는 사이트 표준 동등물

## 4. 최소 설치 매트릭스

| 컴포넌트 | Pre-qsub 전용 | 현재 활성 기준 | 확장 연구 |
|---|---|---|---|
| `pytest`, `pyyaml`, `numpy`, `pandas`, `scipy`, `matplotlib` | 필수 | 필수 | 필수 |
| `vina` Python 패키지 | 불필요 | 필수 | 필수 |
| PDBQT 준비 도구 (`mk_prepare_ligand.py`, `prepare_receptor`, `obabel` 등) | 불필요 | 권장, 종종 필수 | 권장 |
| `rdkit` | 선택 | 권장 | 권장 |
| PyRosetta | 불필요 | 필수 | 필수 |
| LightDock | 불필요 | Phase 1 검증에 필수 | Phase 1 수렴 사용 시 필수 |
| `fpocket` / `P2Rank` | 불필요 | 선택 | Phase 2 포켓 제안에 필수 |
| `gromacs` / `MDAnalysis` | 불필요 | 선택 | MD 안정성 작업에 필수 |

## 5. 레포지토리 스크립트 동작 방식

- `scripts/setup_test_env.sh`
  - 환경을 생성하지 않음
  - `pyrosetta`에 필수 패키지가 있는지만 검증

- `scripts/run_pre_qsub_checks.sh`
  - `pyrosetta`를 활성화
  - 경량 pre-qsub 검증 레인 실행

- `config/run_pre_qsub_checks.pbs`
  - PBS를 통해 동일한 검증 레인 실행
  - `output/pre_qsub_status/last_pass.json`에 기록
  - 환경 검사 또는 테스트 실행 실패 시 `status: failed` 기록

## 6. 첫 실행 검증 체크리스트

환경이 준비되었다고 판단되면, 다음 순서로 확인한다.

### 6.1 공유 환경 활성화

```bash
source ~/.bashrc
conda activate pyrosetta
```

### 6.2 기본 레포지토리 상태 확인

```bash
cd ~/codex_ligand
bash scripts/setup_test_env.sh
bash scripts/run_pre_qsub_checks.sh
```

또는 PBS를 통해:

```bash
qsub config/run_pre_qsub_checks.pbs
```

### 6.3 메인 CLI 확인

```bash
python main.py --help
python main.py validate --help
```

### 6.4 LightDock 가용성 확인

최소한 다음 명령어가 "command not found" 없이 해석되어야 한다:

```bash
which lightdock3_setup.py
which lightdock3.py
which lgd_generate_conformations.py
which lgd_cluster_bsas.py
```

레포지토리 모듈 접근 가능 여부도 확인 가능:

```bash
python -m egfr_pipeline.phase1.lightdock_validation --help
```

### 6.5 PyRosetta 가용성 확인

사이트 정책이 허용하는 경우:

```bash
python -c "import pyrosetta; print('PyRosetta import OK')"
```

### 6.6 Vina 스택 확인

리간드 도킹을 실행할 계획인 경우:

```bash
python -c "import vina; print('vina import OK')"
which obabel
```

Meeko 또는 ADFR 준비를 사용하는 경우:

```bash
which mk_prepare_ligand.py
which prepare_receptor
which prepare_ligand
```

## 7. 이 레포지토리가 하지 않는 것

이 레포지토리는 현재 다음을 **하지 않는다**:

- conda 환경 자동 생성
- PyRosetta 자동 설치
- LightDock 자동 설치
- 서버 측 외부 도구 자동 설치
- 패키지 미러 또는 오프라인 패키지 소스 선택

이러한 부분은 연구실/서버 정책에 따라 관리해야 한다.

## 8. 일반적인 첫 실행 실패 유형

대부분의 초기 실패는 다음 중 하나에 해당한다:

- `pyrosetta` 환경이 활성화되어 있지만 필수 Python 패키지 누락
- LightDock 패키지가 설치되었지만 커맨드라인 도구가 `PATH`에 없음
- PDBQT 준비 도구가 없어서 Vina 입력 변환 실패
- 사용자가 AFM이 여전히 일상 기준의 일부라고 가정
- 사용자가 MD가 기본 프로덕션 경로의 일부라고 가정 (현재는 다운스트림 게이트)

## 9. 하나의 실용 규칙

첫 성공적인 설정을 위해 다음이 참인지 확인한다:

- 공유 conda env `pyrosetta` 활성화
- 기본 Python 패키지 설치
- PyRosetta 동작
- `vina` 동작
- 최소한 하나의 PDBQT 준비 경로 동작
- LightDock 실행 파일이 `PATH`에 있음
- 대규모 제출 전에 pre-qsub 검사 통과

## 10. 환경 설정 후 권장 읽기 순서

환경이 준비되면 다음을 읽는다:

1. `docs/AI_START_HERE.md`
2. `docs/current_pipeline_status.md`
3. `docs/runbook.md`
4. `docs/data_flow_guide.md`
5. `docs/phase1_lightdock_validation_note.md`
6. `docs/pre_qsub_test_line.md`
