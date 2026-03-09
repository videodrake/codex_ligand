# 🚀 Project Specification: PyRosetta PPI Docking Pipeline (V1.0 Stable)

## 1. 프로젝트 개요

이 프로젝트는 **PyRosetta**를 기반으로 한 **고성능 단백질-단백질 상호작용(PPI) 도킹 및 분석 파이프라인**입니다.
단일 PDB 파일을 입력받아 구조 이완(Relax), 전역 도킹(Global Docking), 클러스터링(Clustering), 정밀화(Refinement)를 거쳐 최적의 결합 구조를 예측합니다.

* **상태:** V1.0 Stable (모든 주요 버그 수정 완료)
* **환경:** Linux Cluster, Python 3.x, PyRosetta (구버전/신버전 호환성 확보)
* **핵심 기술:** Multiprocessing(병렬 처리), Custom Interface Analysis, Greedy Clustering

## 2. 파일 구조 및 역할

### 📁 `pipeline_manager.py` (Orchestrator)

* **역할:** 전체 파이프라인의 흐름 제어, 멀티프로세싱 관리, 파일 입출력, 설정 로드.
* **주요 로직:**
* `config.ini` 파싱 및 유효성 검사.
* **Step 1:** `FastRelax` 실행 (결과는 `relaxed_cache`에 캐싱하여 재사용).
* **Step 2:** Global Blind Docking (수천~수만 개 모델 생성).
* **Step 2.5:** Scoring & Filtering (평균 - N*표준편차 컷오프). **중요: `pool.imap`과 `zip`을 사용하여 입력 구조와 점수 결과의 1:1 매핑 무결성 보장.**
* **Step 3:** RMSD 기반 Greedy Clustering (대표 구조 선발).
* **Step 4:** Local Refinement (대표 구조 주변 미세 조정).
* **Step 5:** 최종 Scoring 및 CSV/PDB 출력.
* **Step 7:** PyMol 스크립트 및 에너지 Funnel Plot 생성.



### 📁 `docking.py` (Worker: Sampling)

* **역할:** 실제 Rosetta 프로토콜을 수행하여 구조를 변형(Sampling)함.
* **주요 함수:**
* `run_relax_task`: `FastRelax` 수행. 딕셔너리 형태 반환(`status`, `pdb_data`).
* `run_global_docking_task`: `RigidBodyPerturbMover`(360도 회전, 100A 이동) + `DockMCMProtocol` 수행.
* `run_refinement_task`: `RigidBodyPerturbMover`(작은 범위) + `MinMover`(Sidechain/Backbone/Jump 최소화).


* **주의사항:** 모든 함수는 예외 발생 시 `traceback`을 포함한 에러 딕셔너리를 반환함.

### 📁 `analysis.py` (Worker: Scoring)

* **역할:** 생성된 구조의 에너지 계산, RMSD 측정, 상호작용 잔기 분석.
* **호환성 방어 코드 (Safe Mode):**
* **`residue_energies`:** 구버전 PyRosetta에서 해당 속성이 없을 경우, 에러 없이 총점(Total Score)만 기록하고 상세 점수는 0 처리.
* **`set_calc_sc`:** 최신 버전에서는 `apply()`에 통합되었으므로, 명시적 호출 제거 및 `get_sc_value` 안전 조회.
* **`inject_bfactors`:** B-factor 주입 시, PDBInfo의 엄격한 타입 체크를 피하기 위해 **모든 원자(Atom)를 순회하며 값을 주입**하는 방식 사용.
* **Vector Math:** `xyzVector` 객체 간 연산 오류 방지를 위해 `float` 값으로 변환 후 계산.



### 📁 `common.py` (Utility)

* **역할:** PyRosetta 초기화 및 데이터 변환.
* **특징:**
* **Idempotency:** `pyrosetta.init()` 중복 호출 시 `RuntimeError`를 무시하여 안전성 확보.
* **Validation:** 빈 문자열이나 잘못된 PDB 데이터 입력 시 `None`을 반환하여 파이프라인 중단 방지.



## 3. 핵심 알고리즘 및 설정 (`config_10k.ini` 기준)

1. **System:** 32 CPU Cores.
2. **Docking:** Global Blind Docking 10,000 모델 생성.
3. **Filter:** `sigma = -1.0` (상위 약 16% 수준 통과), 최소 500개 생존 보장.
4. **Cluster:** RMSD Threshold 4.0Å. 상위 50개 클러스터 그룹 선정.
5. **Refinement:** 선정된 각 클러스터 대표 구조당 50번씩 추가 정밀 도킹 (0.1Å / 2.0도 섭동).
6. **Result:** 최종 상위 20개 모델 저장.

## 4. 해결된 주요 이슈 (Context Memory)

이 프로젝트를 수정하거나 확장할 때 **절대 건드리지 말아야 할(또는 주의해야 할) 부분**입니다.

1. **PyRosetta 버전 호환성:**
* 서버의 PyRosetta 버전이 `residue_energies` 속성이나 `set_calc_sc` 메서드를 지원하지 않을 수 있음.
* 따라서 `analysis.py`에는 `try-except`와 `hasattr` 체크가 필수적임. 함부로 최신 API로 변경 금지.


2. **Multiprocessing Index Mapping:**
* 병렬 처리 시 순서가 섞이거나 에러 발생 시 인덱스가 밀리는 현상이 있었음.
* 현재 `pool.imap`의 결과를 `zip(inputs, results)`로 묶어서 처리하는 방식으로 해결함. 이 구조 유지 필수.


3. **네트워크 차단 환경:**
* 클러스터가 인터넷(DNS)에 연결되지 않음. `conda update` 불가능.
* 따라서 코드는 **현재 설치된 구버전 라이브러리에서도 돌아가도록 작성됨 (Polyfill/Fallback 전략 사용).**



## 5. 실행 가이드

1. **환경 준비:** `conda activate pyrosetta`
2. **설정 파일:** `config_10k.ini` 준비 (위 내용 참조).
3. **실행:**
```bash
python pipeline_manager.py config_10k.ini

```


4. **출력:** `[PDB파일명]/final_result/` 폴더에 `Rank01_...pdb` 및 CSV 파일 생성됨.

---

### 🤖 다음 AI를 위한 프롬프트

> "위의 **'Project Specification'** 내용을 바탕으로 현재 PyRosetta 도킹 파이프라인의 구조와 제약 사항(특히 버전 호환성 및 네트워크 차단 문제)을 완벽히 이해해.
> 현재 코드는 V1.0 Stable 상태이며, `pipeline_manager.py`, `docking.py`, `analysis.py`, `common.py` 4개의 파일로 구성되어 있어.
> **작업 목표:** [이 프로젝트가 안정적으로 구동 될 수 있도록, 모든 분석 지표들을 출력할 수 잇도록 하는것]"