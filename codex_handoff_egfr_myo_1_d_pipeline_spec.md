# 코덱스 인수인계 문서

## EGFR C-lobe 표적 탐색 및 PPI/리간드 통합 분석 파이프라인 리팩터링 명세서

---

# 0. 문서의 목적

이 문서는 현재 진행 중인 **EGFR–MYO1D 연구 프로젝트**에 사용되는 세 가지 계산 프로그램

1. **PyRosetta 기반 global docking 프로그램**
2. **로컬 AlphaFold-Multimer 실행 프로그램**
3. **AutoDock Vina 실행 프로그램**

을 **정리, 수정, 통합, 표준화**하기 위해 작성된 **코덱스용 인수인계 문서**다.

이 문서는 단순한 기능 추가 요청서가 아니라, 현재 연구의 방향, 계산의 우선순위, 출력 데이터 구조, 파일 체계, 해석 목적, 후속 비교 분석까지 모두 포함한 **개발 명세서 + 분석 설계서 + 리팩터링 가이드**다.

코덱스는 이 문서를 바탕으로 아래 목표를 달성해야 한다.

- 기존에 따로 돌던 프로그램들을 **하나의 일관된 파이프라인 관점**으로 재구성한다.
- 모든 결과가 사람이 보기 좋은 수준을 넘어서, **후속 구조 해석과 비교 분석에 바로 들어갈 수 있는 표준화된 출력**으로 나오도록 한다.
- 현재 연구의 핵심 질문인 **EGFR kinase domain C-lobe 상의 MYO1D 관련 결합면 및 소분자 민감 포켓 후보를 구조 상태별로 비교 가능하게 한다.**
- 기존 보고서의 잔기나 결론을 절대 기준으로 삼지 않고, **새롭게 생산되는 구조/도킹/포켓 데이터가 더 높은 우선순위**를 가지는 방향으로 분석 체계를 바꾼다.
- 최종적으로는 사용자가 각 구조와 각 리간드에 대해 얻은 결과를 **표, CSV, JSON, pocket summary, receptor 간 비교표**로 쉽게 확인할 수 있어야 한다.

이 문서는 가능한 한 누락 없이 상세하게 작성되며, 코덱스는 이 문서를 기준으로 모듈 설계, CLI 설계, 출력 포맷 설계, 코드 리팩터링, 통합 실행 스크립트 작성, 테스트 코드 작성까지 수행해야 한다.

---

# 1. 연구 배경과 현재 계산의 역할

## 1.1 연구의 생물학적 배경

이 프로젝트의 중심은 **EGFR 비활성 상태와 MYO1D 결합 축**이다.

기존 많은 EGFR 연구는 활성화 이후의 신호전달, kinase 활성, 리간드 결합 후 내재화, 돌연변이 활성화, ATP pocket 억제제에 집중해 왔다. 그러나 이 프로젝트는 그 이전 단계인 **비활성 EGFR이 세포막에서 어떤 구조 상태로 존재하며, 그 상태에서 MYO1D와 어떤 방식으로 연결될 수 있는가**를 구조적으로 파악하는 데 목적이 있다.

특히 현 시점에서 가장 중요한 계산 과제는 다음과 같다.

- EGFR kinase domain, 특히 **C-lobe 표면**에 대해
- MYO1D 결합과 관련될 가능성이 있는 **PPI 결합면 후보**를 정리하고
- 동시에 소분자가 결합 가능한 **대체 포켓 후보**를 찾고
- 이 포켓 및 잔기 패턴이 구조 상태에 따라 얼마나 재현되는지 비교하는 것

즉, 계산의 역할은 단순 docking 점수 비교가 아니라 다음과 같이 정의된다.

1. **구조 상태별 receptor 표면의 차이 파악**
2. **MYO1D 관련 결합면을 추적할 수 있는 PPI 후보 표면 분석**
3. **소분자 민감 포켓의 구조 상태별 재현성 평가**
4. **같은 pocket/patch가 서로 다른 receptor 상태에서 반복되는지 비교**
5. **후속 해석을 위한 residue-level 데이터셋 구축**

이 연구는 현재 wet 실험 확장이 어려운 상황을 전제로 하므로, 계산 결과는 반드시 **재현성, 비교 가능성, 표준화**가 높아야 한다.

---

## 1.2 현재 연구에서 고정된 전제

아래 내용은 현재 사용자와 지도교수 간 합의 또는 연구 방향상 고정된 전제로 간주한다.

### 1.2.1 receptor 구조는 현재 아래 3개를 중심으로 비교한다.

1. **3GT8 원본 구조**
2. **MD 진행 후 38–48 구간 클러스터 대표 구조**
3. **MD 진행 후 85–100 구간 클러스터 대표 구조**

즉, 본 파이프라인의 기본 receptor ensemble은 위 3개다.

이 세 구조는 동일한 EGFR kinase domain 맥락에서 비교되며, residue numbering과 chain 처리의 일관성이 필수다.

### 1.2.2 현재 중요한 목표는 “완벽한 단일 결론”이 아니라 “정량적으로 비교 가능한 결과물”이다.

과거 분석에서는 기존 C02/C04/C07 같은 site naming, 과거 blind docking 결과, 일부 기존 보고서 내 잔기 해석 등이 존재하지만, 앞으로는 그것들을 **후순위 참고사항**으로 두며, **새롭게 산출되는 데이터가 우선**이다.

즉, 소프트웨어는 과거 이름 체계에 종속되지 않아야 한다.

### 1.2.3 포켓 비교의 기본 단위는 receptor마다 따로 정의된다.

사용자는 현재 Vina 결과를 바탕으로 **pose centroid 간 4 Å 이하를 같은 포켓으로 간주하여 clustering**하고 있다. 이 기준은 일단 유지할 수 있으나, 이 결과는 어디까지나 **포켓 분류의 1차 작업 기준**일 뿐이다.

따라서 소프트웨어는 반드시 다음을 모두 남겨야 한다.

- 원시 pose centroid 좌표
- pose별 affinity
- pose별 contact residue
- pocket clustering 결과
- pocket centroid
- pocket residue union 및 residue frequency

즉, 최종 pocket label만 저장하는 방식은 불충분하다.

### 1.2.4 공통 잔기 탐색은 조건부 분석이다.

포켓 4개가 서로 완전히 다른 위치일 수도 있다. 이 경우 전체 포켓의 공통 잔기를 찾는 것은 무의미할 수 있다. 따라서 파이프라인은 다음을 구분할 수 있어야 한다.

- 같은 큰 표면 patch의 sub-pocket인지
- 일부 겹치는 pocket인지
- 완전히 독립적인 pocket인지

이를 위해 receptor 간 pocket 비교표와 residue overlap 계산 기능이 반드시 필요하다.

### 1.2.5 계산 결과는 사람이 바로 판단할 수 있도록 “보이는 출력”이 중요하다.

현재 사용자는 여러 프로그램을 개별적으로 돌려왔고, 결과도 산발적으로 존재할 가능성이 높다. 따라서 각 프로그램은 다음 기준을 만족해야 한다.

- 실행이 자동화되어야 한다.
- 로그가 명확해야 한다.
- 실패 시 어느 단계에서 실패했는지 보여야 한다.
- 결과가 CSV/JSON/요약 텍스트/이미지로 동시에 남아야 한다.
- 사람이 포켓별 residue와 약물 분포를 직관적으로 볼 수 있어야 한다.

---

# 2. 현재 존재하는 프로그램과 그 재정의된 역할

본 프로젝트에는 현재 최소 세 종류의 프로그램이 존재한다.

1. PyRosetta 기반 global docking
2. 로컬 AlphaFold-Multimer
3. AutoDock Vina 실행 프로그램

이 세 프로그램은 각각 개별 목적이 있었겠지만, 앞으로는 아래처럼 역할을 재정의한다.

---

## 2.1 PyRosetta global docking 프로그램

### 현재 역할(추정)

- receptor–ligand가 아니라, **protein–protein global docking**을 수행
- 주로 EGFR kinase domain과 MYO1D 일부 domain/segment의 상대 배향 후보 탐색에 사용
- 대규모 decoy 생성 및 scoring/filtering용

### 앞으로의 역할

PyRosetta 프로그램은 단순 decoy 생성기가 아니라 아래 기능을 수행해야 한다.

1. **입력된 receptor–partner 조합에 대해 reproducible한 global docking 수행**
2. 각 decoy 또는 cluster 결과를 후속 분석 가능한 표준화된 형식으로 출력
3. interface residue, interface score, dSASA, shape complementarity, packstat, buried unsat HB 등을 **구조화된 CSV/JSON으로 저장**
4. decoy clustering 후 cluster representative를 자동 추출
5. 대표 pose에 대해 interface residue 목록과 residue frequency를 자동 계산
6. Vina pocket 결과와 later stage에서 비교할 수 있도록 receptor residue patch를 뽑아내는 기반 자료 제공

### 반드시 수정할 점

- 현재 프로그램이 단순 silent file 또는 score file만 남긴다면, 사람이 읽을 수 있는 **summary table**을 추가해야 한다.
- interface residue 계산이 내부적으로만 존재한다면 외부 파일로 내보내야 한다.
- docking decoy가 많을 경우 전체를 저장하되, top N / cluster representative / per-cluster residue consensus까지 자동 생성해야 한다.

### 핵심 출력

- `global_docking_decoy_scores.csv`
- `global_docking_cluster_summary.csv`
- `global_docking_interface_residues.csv`
- `global_docking_top_models/`
- `global_docking_report.md`

---

## 2.2 로컬 AlphaFold-Multimer 프로그램

### 현재 역할(추정)

- EGFR와 MYO1D 조합 또는 일부 도메인 조합에 대해 복합체 예측 수행
- 배향 후보 또는 접촉 잔기 후보 확인용

### 앞으로의 역할

AlphaFold-Multimer는 절대 정답 생성기가 아니라 **보조 구조 가설 생성기**로 사용한다. 따라서 소프트웨어는 다음을 만족해야 한다.

1. 입력 sequence 또는 domain fragment 정의를 명확하게 저장한다.
2. 실행 파라미터, model preset, recycles, seed 등을 기록한다.
3. 결과 모델마다 pLDDT, pTM/ipTM, PAE 등 주요 지표를 표로 저장한다.
4. 복합체 접촉 잔기와 interface residue를 구조 파일로부터 재추출한다.
5. 구조 신뢰도가 낮은 경우에도 어떤 residue patch가 반복되는지 비교할 수 있게 한다.
6. PyRosetta docking 결과와 비교할 수 있도록 residue-level output을 생성한다.

### 반드시 수정할 점

- 결과를 단순 PDB와 ranking\_debug 정도만 남기지 말고, **모델별 요약 CSV와 interface residue table**을 남겨야 한다.
- ipTM 또는 interface confidence가 낮더라도, residue contact map 추출 기능이 필요하다.
- 여러 입력 조합을 batch로 돌릴 수 있어야 하며, 실패한 job과 성공한 job이 구분되어야 한다.

### 핵심 출력

- `afm_run_metadata.json`
- `afm_model_summary.csv`
- `afm_interface_residue_table.csv`
- `afm_selected_models/`
- `afm_report.md`

---

## 2.3 AutoDock Vina 프로그램

### 현재 역할(추정)

- receptor 3개 상태에 대해 여러 후보약물 docking
- pose centroid를 기준으로 pocket clustering
- C-lobe 상의 포켓 후보 탐색

### 앞으로의 역할

Vina 프로그램은 이번 리팩터링의 중심이다. 단순 docking launcher가 아니라 아래 기능을 갖는 **포켓 데이터 생성 엔진**이 되어야 한다.

1. receptor 3개 각각에 대해 후보약물 세트를 batch docking한다.
2. 모든 ligand pose의 affinity와 좌표를 구조화된 표로 추출한다.
3. 각 pose의 centroid를 계산한다.
4. 사용자가 설정한 기준(현재 기본값 4 Å)으로 pose들을 pocket으로 cluster한다.
5. 각 pose에 대해 receptor contact residue를 계산한다.
6. pocket별 residue union, residue frequency, ligand diversity를 계산한다.
7. receptor 간 pocket overlap 분석을 자동화한다.
8. top pose, top pocket, ligand-to-pocket mapping을 보기 쉽게 출력한다.

### 반드시 수정할 점

- 현재 Vina 결과를 사람이 직접 열어보거나 PyMOL에서 수동 확인해야 하는 상태라면, **기계가 읽는 데이터 + 사람이 읽는 요약표**를 동시에 생성해야 한다.
- 단순 affinity 기준 정렬만으로는 부족하며, pocket 중심의 요약이 필요하다.
- receptor별, ligand별, pocket별 결과를 각각 다른 관점에서 볼 수 있어야 한다.

### 핵심 출력

- `vina_pose_table.csv`
- `vina_pocket_table.csv`
- `vina_drug_pocket_map.csv`
- `vina_pocket_overlap_table.csv`
- `vina_top_poses/`
- `vina_report.md`

---

# 3. 전체 파이프라인 아키텍처의 목표

세 프로그램은 각자 돌고 끝나는 것이 아니라, 아래와 같은 통합 흐름을 가지도록 설계한다.

## 3.1 통합 분석의 기본 흐름

1. receptor 준비
2. receptor metadata 저장
3. PyRosetta global docking 수행(선택적)
4. AlphaFold-Multimer 수행(선택적)
5. AutoDock Vina 수행
6. Vina pose 파싱
7. contact residue 계산
8. pose centroid 기반 pocket clustering
9. pocket 요약 생성
10. receptor 간 pocket 비교
11. 기존 PPI 결과와의 참조 비교(선택적)
12. 최종 summary report 생성

즉, 세 프로그램은 앞으로 하나의 상위 orchestration 레이어 아래 묶일 수 있어야 한다.

---

## 3.2 사용자 관점의 핵심 요구사항

사용자는 복잡한 내부 모듈 구조보다 다음을 원한다.

- receptor 세트와 ligand 세트를 지정하면 한 번에 계산이 돌아갈 것
- 결과가 폴더 구조상 정돈되어 있을 것
- 구조별로 어떤 pocket이 나왔는지 바로 보일 것
- 특정 ligand가 어느 pocket에 들어갔는지 바로 보일 것
- receptor 상태가 바뀌어도 같은 patch인지 비교할 수 있을 것
- 나중에 ChatGPT나 코덱스, 혹은 다른 분석 도구에 다시 넘기기 쉬운 CSV와 JSON이 생성될 것

따라서 파이프라인은 “연구자 혼자 보는 코드”가 아니라 **추후 반복 사용 가능한 반제품 소프트웨어** 수준으로 정리되어야 한다.

---

# 4. 결과 해석을 위해 필수인 데이터 종류

이 절은 가장 중요하다. 코덱스는 “무슨 계산을 돌릴 것인가”보다 **무슨 데이터를 표준화해서 남겨야 하는가**에 집중해야 한다.

## 4.1 receptor 메타데이터

각 receptor 구조에 대해 반드시 아래를 저장한다.

- receptor\_id
- source\_type (`raw_3GT8`, `md_cluster_38_48`, `md_cluster_85_100`)
- original\_file\_path
- cleaned\_file\_path
- pdbqt\_file\_path
- chain\_used
- residue\_range
- receptor\_preparation\_method
- protonation / hydrogen treatment 정보
- remarks

예시:

| receptor\_id    | source\_type         | chain | residue\_range | pdbqt\_file                          |
| --------------- | -------------------- | ----- | -------------- | ------------------------------------ |
| 3GT8\_raw       | raw\_3GT8            | A     | 672-998        | receptors/3GT8\_raw/3GT8\_raw\.pdbqt |
| 3GT8\_cl38\_48  | md\_cluster\_38\_48  | A     | 672-998        | receptors/cluster\_38\_48/rep.pdbqt  |
| 3GT8\_cl85\_100 | md\_cluster\_85\_100 | A     | 672-998        | receptors/cluster\_85\_100/rep.pdbqt |

---

## 4.2 ligand 메타데이터

각 ligand에 대해 반드시 아래를 저장한다.

- ligand\_id
- display\_name
- input\_file\_path
- molecular\_weight (선택)
- formal\_charge (선택)
- source\_or\_series\_name (선택)
- experiment\_annotation (강/중/약 등)
- notes

특히 가능하다면 다음 annotation을 받아서 저장하도록 한다.

- MYO1D binding inhibition strength
- cell death / cytotoxicity annotation
- priority group

이 정보는 계산 결과 해석에 매우 유용하다.

---

## 4.3 Vina pose 원시 데이터

각 ligand × receptor 조합에 대해 모든 pose에 대해 아래를 저장한다.

- receptor\_id
- ligand\_id
- pose\_rank
- affinity
- rmsd\_lb (있으면)
- rmsd\_ub (있으면)
- centroid\_x
- centroid\_y
- centroid\_z
- contact\_residues
- n\_contact\_residues
- pocket\_id

이 데이터를 담는 파일이 `vina_pose_table.csv`다.

이 파일은 모든 후속 분석의 가장 중요한 원시표다.

---

## 4.4 pocket 요약 데이터

각 receptor에서 형성된 pocket마다 아래를 계산해 저장한다.

- receptor\_id
- pocket\_id
- centroid\_x
- centroid\_y
- centroid\_z
- n\_pose
- n\_ligand
- ligand\_ids
- union\_contact\_residues
- residue\_frequency\_table 또는 직렬화된 요약
- top\_residues
- representative\_pose
- mean\_affinity
- best\_affinity
- notes

이 데이터를 담는 파일이 `vina_pocket_table.csv`다.

이 파일은 “이 receptor에서 어떤 pocket들이 나왔는가”를 보는 핵심 표다.

---

## 4.5 ligand-to-pocket 매핑 데이터

각 ligand가 각 receptor에서 주로 어느 pocket에 들어갔는지를 정리한다.

- ligand\_id
- receptor\_id
- dominant\_pocket\_id
- best\_affinity
- dominant\_pocket\_pose\_count
- dominant\_pocket\_fraction
- top\_pose\_residues
- alternative\_pockets (있으면)

이 데이터를 담는 파일이 `vina_drug_pocket_map.csv`다.

이 파일은 “특정 약물이 어느 pocket으로 유도되는가”를 보기 위한 표다.

---

## 4.6 receptor 간 pocket 비교 데이터

이 데이터는 매우 중요하다. receptor 3개 간 pocket 재현성을 보기 위해 필요하다.

각 receptor A의 pocket i와 receptor B의 pocket j에 대해 다음을 계산한다.

- receptor\_a
- pocket\_a
- receptor\_b
- pocket\_b
- centroid\_distance
- residue\_overlap\_n
- residue\_union\_n
- jaccard\_index
- shared\_ligands
- same\_patch\_candidate (heuristic boolean)
- notes

이 파일은 `vina_pocket_overlap_table.csv`로 저장한다.

---

# 5. 개발자가 반드시 이해해야 할 핵심 분석 철학

이 절은 중요하다. 코덱스는 단순 자동화 코드만 짜면 안 되고, 왜 이 출력이 필요한지 이해해야 한다.

## 5.1 포켓은 절대값이 아니라 구조 상태 의존적 가설이다.

동일한 ligand라도 3GT8 원본, 38–48 cluster 대표, 85–100 cluster 대표 구조에서 서로 다른 pocket으로 갈 수 있다. 이는 오류가 아니라 **구조 상태 의존성의 신호**일 수 있다.

따라서 파이프라인은 결과를 한 pocket으로 억지 통합하지 말고, 먼저 **구조별 pocket landscape**를 충실히 남겨야 한다.

## 5.2 공통 잔기 탐색은 전체 포켓에 일괄 적용하는 기본 기능이 아니다.

포켓이 서로 다른 위치라면 공통 잔기를 찾는 것이 무의미할 수 있다. 따라서 코드 수준에서 다음을 지원해야 한다.

- pocket 간 centroid 거리 비교
- residue overlap 비교
- 일부 pocket만 묶어서 patch로 정의하는 기능
- 전체가 아니라 사용자 지정 subset pocket 비교

## 5.3 점수 하나로 결론 내리면 안 된다.

affinity, n\_pose, residue frequency, ligand diversity, receptor 간 재현성은 서로 다른 의미를 가진다. 따라서 summary report는 다중 기준을 나열해야 한다.

예를 들어 특정 pocket이

- affinity는 좋지만 ligand 다양성은 낮을 수 있고
- 반대로 affinity는 보통이지만 여러 receptor 상태에서 반복될 수 있다.

따라서 보고서는 단일 rank만 출력하지 말고, 여러 지표를 병렬 표시해야 한다.

## 5.4 “보이는 데이터”가 중요하다.

사용자는 계산을 직접 설계/판단해야 하므로, CSV만 있으면 충분하지 않다. 최소한 다음이 필요하다.

- pocket별 top residue bar chart 또는 frequency 요약
- receptor별 pocket 개요 텍스트 리포트
- ligand별 dominant pocket 요약
- pocket overlap 요약

그래프는 선택사항이지만 텍스트 리포트는 필수다.

---

# 6. 구체적 리팩터링 요구사항

이 절에서는 코덱스가 각 프로그램을 어떻게 다듬어야 하는지 명시한다.

## 6.1 공통 요구사항

세 프로그램 모두 아래 규칙을 따라야 한다.

### 6.1.1 설정 파일 기반 실행

가능하면 YAML 또는 JSON 설정 파일을 입력으로 받는다.

예:

```yaml
project_name: egfr_myo1d_clobe
receptors:
  - id: 3GT8_raw
    pdb: receptors/3GT8_raw/3GT8_raw.pdb
    pdbqt: receptors/3GT8_raw/3GT8_raw.pdbqt
  - id: 3GT8_cl38_48
    pdb: receptors/cluster_38_48/rep.pdb
    pdbqt: receptors/cluster_38_48/rep.pdbqt
  - id: 3GT8_cl85_100
    pdb: receptors/cluster_85_100/rep.pdb
    pdbqt: receptors/cluster_85_100/rep.pdbqt
ligands:
  - id: drugA
    pdbqt: ligands/drugA.pdbqt
  - id: drugB
    pdbqt: ligands/drugB.pdbqt
vina:
  center: [x, y, z]
  size: [sx, sy, sz]
  exhaustiveness: 32
  num_modes: 20
  energy_range: 6
pocket_clustering:
  centroid_cutoff_angstrom: 4.0
contact_cutoff_angstrom: 4.0
```

### 6.1.2 모든 단계에 대해 로그 저장

- stdout/stderr 로깅
- 단계별 성공/실패 상태
- 입력 파일 존재 여부
- 결과 파일 저장 위치

### 6.1.3 출력 파일명 규칙 통일

모든 산출물은 project root 아래 `parsed/`, `reports/`, `logs/`, `figures/` 등에 정리한다.

### 6.1.4 dry-run 모드 지원

실제 실행 대신 입력만 점검하고 어떤 작업이 수행될지 보여주는 모드가 필요하다.

### 6.1.5 재현성 확보

- random seed 기록
- 사용한 소프트웨어 버전 기록
- 실행 파라미터 저장

---

## 6.2 PyRosetta 프로그램 리팩터링 요구사항

### 6.2.1 입력 명세

입력은 최소 다음을 받는다.

- receptor structure file
- partner structure file
- chain selection
- number of decoys
- docking mode
- score function 선택
- output directory

### 6.2.2 기능 요구사항

- decoy generation
- score extraction
- interface residue extraction
- clustering
- top model saving
- cluster representative saving
- per-cluster residue frequency 계산

### 6.2.3 출력 요구사항

#### 필수 CSV

- `pyrosetta_decoy_scores.csv`
- `pyrosetta_cluster_summary.csv`
- `pyrosetta_interface_residue_table.csv`

#### 필수 구조 파일

- top N models
- cluster representative models

#### 필수 요약 리포트

- docking total count
- best score models
- cluster count
- interface residue hotspots

### 6.2.4 추가 요구사항

- 기존 C02/C04 같은 naming을 하드코딩하지 말 것
- interface residue를 residue index와 residue name으로 표준화할 것
- later comparison을 위해 receptor-side residues와 partner-side residues를 분리해서 기록할 것

---

## 6.3 AlphaFold-Multimer 프로그램 리팩터링 요구사항

### 6.3.1 입력 명세

- sequence FASTA 또는 domain fragment definition
- chain labels
- output dir
- model preset
- num recycles
- random seeds

### 6.3.2 기능 요구사항

- batch run 지원
- metadata 저장
- model ranking 저장
- interface residue extraction
- PAE / ipTM / pLDDT 정리

### 6.3.3 출력 요구사항

- `afm_model_summary.csv`
- `afm_interface_residue_table.csv`
- selected PDB models
- `afm_run_metadata.json`
- `afm_report.md`

### 6.3.4 추가 요구사항

- 각 모델에 대해 receptor-side contact residues를 뽑을 것
- residue contact가 일정 threshold 이상 반복되는 모델만 추려볼 수 있는 옵션 제공
- low confidence result라도 버리지 말고 별도 flag를 붙일 것

---

## 6.4 AutoDock Vina 프로그램 리팩터링 요구사항

이 절이 가장 중요하다.

### 6.4.1 입력 명세

- receptor list
- ligand list
- docking box center/size
- exhaustiveness
- num\_modes
- energy\_range
- contact cutoff
- pocket centroid clustering cutoff (기본 4.0 Å)

### 6.4.2 기능 요구사항

#### A. 배치 실행

각 receptor × ligand 조합에 대해 자동 실행

#### B. 출력 파싱

각 output PDBQT에서

- mode rank
- affinity
- 원자 좌표
- centroid 추출

#### C. contact residue 계산

각 pose에 대해 receptor residue와의 최소 heavy atom 거리 기반으로 contact residue 계산

#### D. pocket clustering

같은 receptor 내 모든 pose centroid를 모아 사용자 cutoff 기준으로 clustering

#### E. pocket summary 생성

pocket별

- centroid
- pose 수
- ligand 수
- residue union
- residue frequency
- top residues
- best affinity
- mean affinity 정리

#### F. ligand-to-pocket mapping 생성

ligand별 dominant pocket 계산

#### G. receptor 간 pocket overlap 분석

각 receptor pair에 대해 pocket끼리 비교

### 6.4.3 출력 요구사항

#### 필수 CSV

- `vina_pose_table.csv`
- `vina_pocket_table.csv`
- `vina_drug_pocket_map.csv`
- `vina_pocket_overlap_table.csv`

#### 필수 JSON

- `vina_run_metadata.json`
- `vina_pocket_table.json`

#### 필수 리포트

- `vina_report.md`

#### 선택 이미지

- pocket별 대표 pose 이미지 생성 스크립트 훅

### 6.4.4 특별 요구사항

- ligand pose clustering cutoff는 설정 가능해야 하나 기본은 4.0 Å
- contact cutoff도 설정 가능하게 할 것
- receptor residue numbering은 원본 PDB와 일치하도록 유지할 것
- centroid 계산과 contact residue 계산은 별도 함수로 모듈화할 것
- pocket assignment는 deterministic해야 한다

---

# 7. 추천 디렉터리 구조

아래 구조를 기본 템플릿으로 제안한다.

```text
project_root/
├── config/
│   └── project.yaml
├── receptors/
│   ├── 3GT8_raw/
│   │   ├── 3GT8_raw.pdb
│   │   └── 3GT8_raw.pdbqt
│   ├── cluster_38_48/
│   │   ├── rep_38_48.pdb
│   │   └── rep_38_48.pdbqt
│   └── cluster_85_100/
│       ├── rep_85_100.pdb
│       └── rep_85_100.pdbqt
├── ligands/
│   ├── drugA.pdbqt
│   ├── drugB.pdbqt
│   └── ...
├── runs/
│   ├── pyrosetta/
│   ├── afm/
│   └── vina/
├── parsed/
│   ├── vina_pose_table.csv
│   ├── vina_pocket_table.csv
│   ├── vina_drug_pocket_map.csv
│   ├── vina_pocket_overlap_table.csv
│   ├── pyrosetta_decoy_scores.csv
│   ├── pyrosetta_cluster_summary.csv
│   ├── pyrosetta_interface_residue_table.csv
│   ├── afm_model_summary.csv
│   └── afm_interface_residue_table.csv
├── reports/
│   ├── vina_report.md
│   ├── pyrosetta_report.md
│   ├── afm_report.md
│   └── integrated_report.md
├── figures/
│   ├── pocket_maps/
│   ├── overlap_heatmaps/
│   └── ligand_assignments/
├── logs/
└── scripts/
```

---

# 8. Vina 결과 파싱 세부 명세

코덱스가 가장 먼저 안정화해야 하는 부분이다.

## 8.1 Vina output 파싱 목표

각 `*_out.pdbqt` 파일에서 모든 pose를 읽어 다음 데이터를 산출한다.

- pose\_rank
- affinity
- ligand atom coordinates
- centroid
- pose별 contact residues

### 파싱 포인트

- `MODEL` / `ENDMDL`
- `REMARK VINA RESULT:`
- ligand atom coordinate lines

## 8.2 centroid 계산 규칙

기본은 **center of geometry**를 사용한다. 질량 중심이 아니라 좌표 평균이다.

### 이유

- 구현 단순성
- ligand 간 비교에 충분함
- pocket clustering 목적상 적절함

## 8.3 contact residue 계산 규칙

기본 기준:

- ligand heavy atom과 receptor residue heavy atom 사이 최소 거리 ≤ `contact_cutoff` Å
- 기본값 4.0 Å

선택 옵션:

- 4.5 Å도 지원 가능
- hydrophobic / polar / hbond 등 추가 분류는 옵션으로 추가 가능

## 8.4 pocket clustering 규칙

기본 규칙:

- 같은 receptor에서 나온 모든 pose centroid를 모음
- centroid 간 거리 cutoff ≤ 4.0 Å면 같은 cluster로 묶음
- clustering 알고리즘은 DBSCAN 또는 연결요소 기반 단순 그래프 clustering 가능
- `min_samples=1` 허용

### 주의

- pocket clustering 기준과 contact cutoff는 완전히 분리된 옵션이어야 한다.
- 둘을 같은 값으로 두더라도 내부적으로는 다른 변수여야 한다.

---

# 9. PyRosetta 결과 파싱 세부 명세

## 9.1 interface residue 추출

각 decoy 또는 cluster representative에 대해 receptor-side interface residue와 partner-side interface residue를 구분해 저장한다.

필수 컬럼 예시:

- decoy\_id
- cluster\_id
- score
- receptor\_interface\_residues
- partner\_interface\_residues
- n\_receptor\_interface\_residues
- n\_partner\_interface\_residues
- dSASA
- sc\_value
- packstat
- buried\_unsat\_hbonds

## 9.2 cluster summary

cluster 단위로 다음을 정리한다.

- cluster\_id
- n\_members
- best\_score
- representative\_model
- receptor residue frequency
- partner residue frequency

### 중요한 점

이 결과는 Vina 결과와 직접 결합하지 않아도 된다. 하지만 later stage에서 receptor residue patch 비교에 사용할 수 있도록 residue 표준화는 필수다.

---

# 10. AlphaFold-Multimer 결과 파싱 세부 명세

## 10.1 모델별 요약

- model\_name
- rank
- ipTM
- pTM
- mean\_pLDDT
- interface\_contact\_count
- receptor\_contact\_residues
- partner\_contact\_residues

## 10.2 interface residue 추출 기준

기본적으로 receptor-side residue와 partner-side residue의 heavy atom 최소 거리 cutoff를 사용한다.

- 기본 cutoff 4.0 Å 또는 5.0 Å
- 설정 파일에서 조절 가능

## 10.3 신뢰도 지표의 역할

- 높은 ipTM만 남기는 hard filter 기능은 optional
- 기본은 모든 모델을 저장하되 flag를 다는 방식

---

# 11. 통합 오케스트레이션 요구사항

코덱스는 개별 프로그램만 고치는 것이 아니라, 가능하다면 상위 실행기 하나를 만들어야 한다.

예시 CLI:

```bash
python run_pipeline.py --config config/project.yaml --steps vina,parse,report
python run_pipeline.py --config config/project.yaml --steps pyrosetta,report
python run_pipeline.py --config config/project.yaml --steps afm,parse,report
python run_pipeline.py --config config/project.yaml --steps all
```

이 상위 실행기는 아래를 해야 한다.

- config 로드
- 입력 파일 검사
- 각 스텝 호출
- 로그 저장
- 실패 시 종료 코드 반환
- 성공 시 리포트 생성

---

# 12. 리포트 설계 요구사항

모든 프로그램은 markdown report를 생성해야 한다. 연구자는 CSV만 보는 것이 아니라, 결과를 빠르게 요약해 읽을 수 있어야 한다.

## 12.1 Vina report 필수 내용

1. receptor 목록
2. ligand 목록
3. receptor별 pocket 수
4. receptor별 top pocket 요약
5. ligand별 dominant pocket 요약
6. receptor 간 pocket overlap 상위 후보
7. 해석시 유의사항

## 12.2 PyRosetta report 필수 내용

1. docking 설정
2. decoy 수
3. cluster 수
4. top cluster summary
5. interface residue hotspot 요약
6. score distribution 요약

## 12.3 AFM report 필수 내용

1. run 설정
2. 모델 ranking
3. confidence 요약
4. interface residue 반복 패턴 요약
5. low confidence 주의사항

## 12.4 integrated report 필수 내용

1. receptor ensemble 설명
2. Vina pocket landscape 요약
3. PyRosetta/AFM의 receptor-side residue 후보 참고 요약
4. receptor 상태별 pocket 재현성 요약
5. 앞으로 사람이 직접 볼 우선 pocket 후보 목록

---

# 13. 사람이 보기 쉬운 출력 형식 요구사항

사용자는 결과를 “보기 쉽게” 받고 싶어 한다. 따라서 코덱스는 단순 CSV 생성에 그치지 말고 다음을 고려해야 한다.

## 13.1 긴 residue list 가독성

CSV 한 셀에 residue list가 너무 길어질 수 있다. 이를 위해 다음 두 방식을 병행한다.

- CSV: `MET971;HIS972;LEU973` 같이 세미콜론 구분 문자열
- JSON: residue frequency를 key-value 형태로 저장

## 13.2 pocket summary용 간단한 정렬 버전

예시:

| receptor  | pocket | ligands | best\_aff | top residues           |
| --------- | ------ | ------- | --------- | ---------------------- |
| 3GT8\_raw | P1     | 4       | -7.3      | MET971, HIS972, LEU973 |

이런 표가 markdown report에 자동 삽입되면 좋다.

## 13.3 receptor 간 overlap heatmap

선택사항이지만 가능하면 생성한다.

- pocket A vs pocket B Jaccard index heatmap
- receptor pair별 table 이미지

## 13.4 ligand별 pocket assignment 요약표

예시:

| ligand | 3GT8\_raw | 38\_48 | 85\_100 |
| ------ | --------- | ------ | ------- |
| drugA  | P1        | P2     | P2      |
| drugB  | P3        | P1     | P4      |

이 표는 매우 중요하다.

---

# 14. 지금 당장 필요한 “데이터를 뽑아내는 실험”이란 무엇인가

이 문맥에서 “실험”은 wet lab 실험이 아니라 **계산 결과를 표준화된 데이터로 추출하는 작업**을 의미한다.

코덱스는 다음 데이터 추출 실험을 가능하게 해야 한다.

## 14.1 receptor metadata extraction

- receptor 구조 파일을 읽고 metadata 저장

## 14.2 Vina pose extraction experiment

- output PDBQT에서 모든 pose와 affinity 추출

## 14.3 contact residue extraction experiment

- pose별 receptor contact residue 추출

## 14.4 pocket clustering experiment

- pose centroid 기반 clustering 수행

## 14.5 pocket summary experiment

- pocket별 residue union / frequency 계산

## 14.6 receptor overlap experiment

- receptor 간 pocket 비교

## 14.7 PyRosetta interface extraction experiment

- top decoy/cluster representative의 interface residue 추출

## 14.8 AFM interface extraction experiment

- 모델별 interface residue 추출

즉, 지금 가장 필요한 것은 새로운 생물학 실험이 아니라 **기존 계산 산출물을 비교 가능한 데이터셋으로 바꾸는 계산 실험**이다.

---

# 15. 코덱스가 구현해야 하는 스크립트 목록

최소 구현 단위는 아래처럼 나눈다.

## 15.1 공통 유틸

- `io_utils.py`
- `path_utils.py`
- `logging_utils.py`
- `config_utils.py`

## 15.2 receptor/ligand metadata

- `collect_metadata.py`

## 15.3 Vina 관련

- `run_vina_batch.py`
- `parse_vina_results.py`
- `extract_contacts.py`
- `cluster_pockets.py`
- `summarize_pockets.py`
- `compare_pockets.py`
- `build_vina_report.py`

## 15.4 PyRosetta 관련

- `run_pyrosetta_global_docking.py`
- `parse_pyrosetta_scores.py`
- `extract_pyrosetta_interface.py`
- `cluster_pyrosetta_models.py`
- `build_pyrosetta_report.py`

## 15.5 AFM 관련

- `run_afm_batch.py`
- `parse_afm_outputs.py`
- `extract_afm_interface.py`
- `build_afm_report.py`

## 15.6 통합

- `run_pipeline.py`
- `build_integrated_report.py`

---

# 16. Vina 핵심 스크립트의 상세 요구사항

## 16.1 `run_vina_batch.py`

### 입력

- config file
- receptor list
- ligand list
- optional subset filters

### 기능

- receptor × ligand 조합 반복
- vina command 생성
- 실행
- stdout/stderr 저장
- 실패 시 기록

### 출력

- raw output PDBQT
- run log
- metadata json

---

## 16.2 `parse_vina_results.py`

### 목표

각 output PDBQT에서 pose별 정보 추출

### 출력 컬럼

- receptor\_id
- ligand\_id
- pose\_rank
- affinity
- centroid\_x
- centroid\_y
- centroid\_z

### 추가 요구

- rmsd 정보가 있으면 포함
- output이 비어 있거나 손상된 경우 graceful handling

---

## 16.3 `extract_contacts.py`

### 목표

pose별 contact residue 추출

### 입력

- receptor PDB
- pose coordinates
- contact cutoff

### 출력

- pose별 residue list
- n\_contact\_residues

### 추가 옵션

- heavy atom only 여부
- per-residue min distance 저장 여부

---

## 16.4 `cluster_pockets.py`

### 목표

동일 receptor 내 모든 pose centroid clustering

### 입력

- pose table
- centroid cutoff

### 출력

- pocket\_id assignment
- pocket centroid
- cluster member list

### 주의

- ligand별로 따로 clustering하면 안 되고 receptor별 전체 pose를 대상으로 해야 한다.

---

## 16.5 `summarize_pockets.py`

### 목표

pocket별 요약 생성

### 출력

- n\_pose
- n\_ligand
- union residues
- residue frequency
- top residues
- mean/best affinity

### 추가

- dominant ligands
- representative pose 선정 규칙 정의

---

## 16.6 `compare_pockets.py`

### 목표

receptor 간 pocket overlap 비교

### 입력

- pocket table
- optional pose table

### 출력

- centroid distance
- residue overlap
- jaccard
- shared ligands

### 선택 기능

- same\_patch\_candidate flag
- threshold 기반 match 추천

---

# 17. PyRosetta 핵심 스크립트의 상세 요구사항

## 17.1 `run_pyrosetta_global_docking.py`

- config 기반 실행
- seed 기록
- decoy PDB 또는 silent file 저장
- score file 저장

## 17.2 `parse_pyrosetta_scores.py`

- score 파일 파싱
- top N, best by score, best by interface score 정리

## 17.3 `extract_pyrosetta_interface.py`

- 구조 파일 기반 interface residue 추출
- receptor-side / partner-side 분리
- residue frequency table 생성

## 17.4 `cluster_pyrosetta_models.py`

- decoy clustering
- representative model 선택
- cluster summary 생성

---

# 18. AFM 핵심 스크립트의 상세 요구사항

## 18.1 `run_afm_batch.py`

- sequence pair batch 실행
- metadata 기록
- run failure logging

## 18.2 `parse_afm_outputs.py`

- ranking/debug JSON 파싱
- pLDDT, ipTM, pTM 추출
- model summary 생성

## 18.3 `extract_afm_interface.py`

- 모델별 interface residue 추출
- receptor-side residue frequency 생성

---

# 19. 우선순위 계획

코덱스는 모든 기능을 한 번에 완벽하게 만들려 하지 말고, 아래 순서로 구현한다.

## 1차 우선순위

1. Vina batch 실행 정리
2. Vina output 파싱
3. contact residue 추출
4. pocket clustering
5. pocket summary CSV 출력

## 2차 우선순위

6. receptor 간 pocket overlap 비교
7. ligand-to-pocket map 출력
8. markdown report 생성

## 3차 우선순위

9. PyRosetta interface output 정리
10. AFM interface output 정리
11. integrated report 생성

## 4차 우선순위

12. 시각화 도구 추가
13. heatmap / summary figure 추가
14. advanced rescoring hooks 추가

---

# 20. 절대로 놓치면 안 되는 설계 원칙

1. 기존 보고서 내용은 절대 기준이 아니다.
2. 새 계산 결과가 항상 우선이다.
3. receptor 3개는 각각 동등하게 취급한다.
4. 포켓 비교 전에 pocket location과 residue set을 남겨야 한다.
5. “공통 잔기”는 자동 기본 출력이 아니라 조건부 해석 대상이다.
6. 출력은 반드시 표준화되어야 한다.
7. 사람이 읽기 쉬운 리포트와 기계가 읽기 쉬운 CSV/JSON을 둘 다 남겨야 한다.
8. residue numbering 일관성이 최우선이다.
9. 모든 CLI는 실패 시 명확한 에러 메시지를 출력해야 한다.
10. 모든 실행 파라미터와 버전 정보를 기록해야 한다.

---

# 21. 코덱스에게 직접 전달할 개발 지시문

아래는 코덱스에게 전달할 직접 지시문으로 사용해도 좋다.

## 코덱스 개발 지시문

당신은 EGFR C-lobe 표면의 PPI/리간드 후보를 구조 상태별로 비교하는 계산 파이프라인을 리팩터링해야 한다. 현재 세 종류의 기존 프로그램이 있다. PyRosetta global docking, 로컬 AlphaFold-Multimer, AutoDock Vina 실행 프로그램이다. 이 프로그램들은 현재 각자 따로 돌아가고 결과도 사람이 재활용하기 어렵다. 이를 하나의 일관된 분석 체계로 개편해야 한다.

가장 우선순위가 높은 것은 AutoDock Vina 프로그램 정리다. receptor는 현재 정확히 세 개다. 3GT8 원본 구조, MD 후 38–48 구간 클러스터 대표 구조, MD 후 85–100 구간 클러스터 대표 구조다. ligand 세트를 이 세 receptor 각각에 docking하고, output pose들을 centroid 기준 4 Å cutoff로 pocket clustering한다. 그러나 최종 pocket ID만 남기는 방식은 불충분하다. 각 pose의 affinity, centroid 좌표, contact residues를 먼저 추출하고, 그 다음 pocket assignment, pocket centroid, union residues, residue frequency, ligand diversity를 계산해야 한다.

결과는 최소한 다음 CSV 네 개를 생성해야 한다. `vina_pose_table.csv`, `vina_pocket_table.csv`, `vina_drug_pocket_map.csv`, `vina_pocket_overlap_table.csv`. 첫 번째는 pose 원시표다. 두 번째는 receptor별 pocket summary다. 세 번째는 ligand가 receptor별로 어느 pocket에 주로 들어갔는지 보여준다. 네 번째는 receptor 간 pocket overlap과 centroid distance, residue overlap, Jaccard index를 정리한다.

PyRosetta global docking 프로그램도 수정해야 한다. 지금은 대량 decoy와 score가 존재할 가능성이 높지만, 이후 비교 분석에 바로 쓰기 어렵다. decoy score summary, cluster summary, interface residue table을 CSV로 내보내야 한다. receptor-side residue와 partner-side residue는 반드시 분리 기록해야 한다. 기존 C02/C04/C07 같은 과거 site 이름에 의존하지 마라.

AlphaFold-Multimer도 마찬가지다. 단순 PDB 생성만으로 끝내지 말고 model summary, ipTM/pTM/pLDDT, interface residue table을 만들어야 한다. low-confidence 모델이라도 완전히 버리지 말고 flag를 달아 남겨라.

전체 프로젝트는 설정 파일 기반으로 실행되어야 하고, 입력 파일 존재 검사, 로그 기록, dry-run 모드, 실행 파라미터 저장, 버전 정보 저장이 가능해야 한다. 폴더 구조는 receptors, ligands, runs, parsed, reports, logs, figures, scripts로 정리한다. CSV 외에 markdown report도 생성해야 한다. 사용자는 사람이 읽기 쉬운 결과를 원한다.

중요한 점은, 이 프로젝트에서 과거 보고서의 잔기나 pocket 해석은 절대 기준이 아니라는 것이다. 앞으로 새로 생성되는 데이터가 더 높은 우선순위를 가진다. 또한 pocket이 여러 개일 때 무조건 공통 잔기를 찾는 것은 잘못될 수 있다. 따라서 먼저 pocket 간 위치 관계와 residue overlap을 계산할 수 있어야 하며, 그 다음에야 같은 patch인지 아닌지 해석할 수 있다. 즉, 공통 잔기 탐색은 기본 기능이 아니라 조건부 비교 기능이다.

우선순위는 다음과 같다. 1차로 Vina batch 실행, Vina parsing, contact extraction, pocket clustering, pocket summary를 완성한다. 2차로 receptor 간 pocket overlap, ligand-to-pocket map, report 생성까지 완성한다. 3차로 PyRosetta와 AFM 결과 파싱을 정리한다. 4차로 integrated report와 시각화 기능을 추가한다.

코드는 모듈화되어야 하고, 향후 사용자가 미니콘다 환경에서 쉽게 실행할 수 있어야 한다. 출력은 반드시 사람이 보기 좋은 형태와 후처리 가능한 형태 둘 다 제공해야 한다.

---

# 22. 코덱스가 생성해야 할 테스트 시나리오

1. receptor 1개 + ligand 2개만 넣었을 때 정상 실행되는지
2. receptor 3개 + ligand 여러 개에서 batch가 정상 반복되는지
3. output PDBQT가 비어 있을 때 오류 처리가 되는지
4. residue numbering이 다른 receptor가 들어오면 경고를 띄우는지
5. contact cutoff와 pocket cutoff를 바꿨을 때 결과가 반영되는지
6. receptor 간 pocket overlap table이 생성되는지
7. report 파일이 빈 내용 없이 생성되는지

---

# 23. 구현 언어 및 패키지 권장

- Python 중심 구현
- `pathlib`, `argparse`, `subprocess`, `csv`, `json`, `logging`
- `pandas`
- `numpy`
- `scikit-learn` 또는 직접 clustering 구현
- `MDAnalysis` 또는 `Biopython` for structure parsing
- 선택적으로 `matplotlib` for figure generation

중요한 것은 의존성이 과도하지 않아야 하고, 설치 문서가 함께 제공되어야 한다는 점이다.

---

# 24. 마무리 요약

이 프로젝트의 핵심은 계산 프로그램을 더 많이 만드는 것이 아니라, **기존 계산 결과를 반복 가능하고 비교 가능한 구조 데이터셋으로 정리하는 것**이다. 세 receptor 상태(3GT8 원본, 38–48 클러스터 대표, 85–100 클러스터 대표)를 중심으로, ligand docking 결과를 pose-level, pocket-level, receptor-comparison-level로 정리해야 한다. PyRosetta와 AlphaFold-Multimer는 보조적 구조 가설 제공 도구로 남기되, 이 역시 residue-level output을 표준화해야 한다.

코덱스는 이 문서를 기준으로 기존 코드베이스를 정리하고, 결과가 사람이 보기 쉽고 후속 해석에 즉시 사용 가능한 수준으로 리팩터링해야 한다. 가장 먼저 완성해야 할 것은 Vina 결과 표준화와 pocket 비교 체계다. 그 위에 PyRosetta/AFM residue output을 얹고, 마지막에 통합 리포트를 만드는 방향으로 개발을 진행한다.

이 문서를 구현 지침의 단일 기준 문서로 사용하라.

---

# 25. 세부 데이터 계약(Data Contract)

이 절은 코덱스가 실제 구현 시 가장 중요하게 지켜야 할 **입력/출력 계약**을 기술한다. 단순한 컬럼 이름 예시가 아니라, 각 파일이 어떤 의미를 가지며 어떤 타입을 가져야 하는지 명시한다.

## 25.1 `receptor_metadata.csv`

필수 컬럼:

- `receptor_id`: 문자열, 전역 고유값
- `source_type`: `raw_3GT8`, `md_cluster_38_48`, `md_cluster_85_100` 중 하나
- `pdb_path`: 절대 또는 프로젝트 루트 기준 상대경로
- `pdbqt_path`: 동일
- `chain_id`: 문자열
- `residue_start`: 정수
- `residue_end`: 정수
- `n_residues`: 정수
- `preparation_notes`: 문자열
- `numbering_checksum`: 문자열 또는 해시

### 목적

이 파일은 모든 후속 산출물이 어떤 receptor를 참조하는지 정의하는 표준 기준 파일이다.

### 주의

- `receptor_id`는 절대 바뀌지 않아야 한다.
- residue numbering이 다르면 `numbering_checksum` 비교에서 드러나야 한다.

---

## 25.2 `ligand_metadata.csv`

필수 컬럼:

- `ligand_id`
- `display_name`
- `pdbqt_path`
- `annotation_myo1d_inhibition` (선택값: `strong`, `medium`, `weak`, `unknown`)
- `annotation_cell_viability` (선택값: `high_killing`, `moderate`, `low`, `unknown`)
- `priority_tier` (`tier1`, `tier2`, `tier3`, `unknown`)
- `notes`

### 목적

실험 annotation이 있는 ligand와 없는 ligand를 later analysis에서 구분할 수 있게 한다.

---

## 25.3 `vina_pose_table.csv`

필수 컬럼 상세:

- `receptor_id`: 문자열
- `ligand_id`: 문자열
- `pose_rank`: 정수
- `affinity`: 실수
- `rmsd_lb`: 실수 또는 결측
- `rmsd_ub`: 실수 또는 결측
- `centroid_x`: 실수
- `centroid_y`: 실수
- `centroid_z`: 실수
- `pocket_id`: 문자열
- `n_contact_residues`: 정수
- `contact_residues`: 세미콜론 구분 문자열
- `contact_residue_min_distances`: JSON 문자열 또는 선택적 별도 파일 참조
- `raw_pose_file`: 문자열

### 권장 추가 컬럼

- `cluster_label_before_normalization`
- `pose_energy_rank_within_ligand`
- `pose_energy_rank_within_receptor`

### 목적

모든 후속 pocket 분석의 원천 테이블이다.

---

## 25.4 `vina_pocket_table.csv`

필수 컬럼 상세:

- `receptor_id`
- `pocket_id`
- `centroid_x`
- `centroid_y`
- `centroid_z`
- `n_pose`
- `n_ligand`
- `ligand_ids`
- `best_affinity`
- `mean_affinity`
- `median_affinity`
- `union_contact_residues`
- `top_residues`
- `residue_frequency_json`
- `representative_pose`
- `representative_ligand`
- `notes`

### 목적

사람이 receptor별 pocket landscape를 빠르게 이해하는 요약 표

---

## 25.5 `vina_drug_pocket_map.csv`

필수 컬럼 상세:

- `ligand_id`
- `receptor_id`
- `dominant_pocket_id`
- `dominant_pocket_pose_count`
- `dominant_pocket_fraction`
- `best_affinity`
- `best_pose_rank`
- `top_pose_residues`
- `alternative_pockets`
- `is_multimodal_binding` (boolean)

### 목적

특정 ligand가 한 receptor 상태에서 어떤 pocket으로 수렴하는지 요약

---

## 25.6 `vina_pocket_overlap_table.csv`

필수 컬럼 상세:

- `receptor_a`
- `pocket_a`
- `receptor_b`
- `pocket_b`
- `centroid_distance`
- `residue_overlap_n`
- `residue_union_n`
- `jaccard_index`
- `shared_ligands`
- `shared_ligands_n`
- `same_patch_candidate`
- `matching_confidence` (`high`, `medium`, `low`)
- `notes`

### 목적

receptor 상태 간 pocket이 같은 구조 patch인지 아닌지 판단하는 1차 비교표

---

# 26. pocket 비교 로직의 상세 설계

이 절은 코덱스가 반드시 정확히 이해해야 한다. 사용자는 이미 “포켓이 4개면 위치가 다를 수 있는데 공통잔기를 찾는 것이 왜 의미가 있느냐”는 핵심 질문을 제기했다. 따라서 소프트웨어는 무조건 공통잔기를 뽑는 것이 아니라, 먼저 pocket 간 관계를 분류할 수 있어야 한다.

## 26.1 pocket 관계 분류 레벨

각 pocket pair는 아래 셋 중 하나로 분류할 수 있어야 한다.

### Type A: 같은 patch의 변형

- centroid distance가 작음
- residue overlap이 큼
- 일부 ligand가 공통적으로 배정됨
- top residues가 상당 부분 겹침

### Type B: 부분 overlap

- centroid distance는 중간
- residue overlap은 제한적
- 같은 큰 면의 edge pocket일 가능성

### Type C: 독립 pocket

- centroid distance 큼
- residue overlap 낮음 또는 0
- shared ligand도 낮음

이 분류는 hard-coded rule로 시작해도 되지만, 출력에서는 항상 raw metric이 먼저 보이고, 판정은 보조적이어야 한다.

---

## 26.2 same\_patch\_candidate 추천 규칙(heuristic)

초기 구현에서는 다음과 같은 휴리스틱을 사용할 수 있다.

- `centroid_distance <= 6.0 Å` AND `jaccard_index >= 0.25` → `same_patch_candidate = True`
- `centroid_distance <= 8.0 Å` AND `shared_ligands_n >= 2` → `same_patch_candidate = True`
- 그 외는 `False`

이 규칙은 프로젝트 초기 기준이며, 반드시 config에서 조정 가능해야 한다.

---

## 26.3 왜 raw metric이 중요하냐

사용자는 각 포켓이 정말 같은 면인지 직접 판단하고 싶어 한다. 따라서 코드는 `same_patch_candidate=True` 같은 이진 플래그를 제공하더라도, 다음 raw 값들을 항상 같이 출력해야 한다.

- centroid distance
- overlap residue count
- jaccard index
- shared ligand count
- 대표 residue 목록

즉, 소프트웨어가 결론을 대신 내려서는 안 된다. 소프트웨어는 **판단 재료를 구조화해서 제공**해야 한다.

---

# 27. contact residue 계산의 상세 명세

## 27.1 기본 원칙

contact residue 계산은 이 프로젝트에서 residue-level 해석의 기반이 되므로 반드시 deterministic하고 재현 가능해야 한다.

### 기본 규칙

- receptor residue heavy atom과 ligand heavy atom 간 최소 거리 사용
- 수소는 기본적으로 제외
- cutoff 기본값은 4.0 Å
- cutoff는 config로 조절 가능

## 27.2 출력 방식

pose 단위로 다음 두 가지를 모두 남긴다.

1. `contact_residues`: residue 문자열 목록
2. `contact_residue_min_distances`: residue별 최소 거리 매핑

예:

```json
{
  "MET971": 3.42,
  "HIS972": 3.81,
  "LEU973": 3.66
}
```

### 이유

단순 residue 존재 여부만으로는 이후 ranking이 어렵다. 최소 거리 정보가 있으면 나중에 tighter filter를 적용할 수 있다.

---

## 27.3 residue string 포맷 통일

반드시 하나의 포맷만 사용한다.

권장 포맷:

- `RESNAME+RESID`, 예: `MET971`, `HIS972`

chain이 여러 개인 경우:

- `A:MET971`

프로젝트 전체에서 같은 포맷을 유지한다.

---

# 28. representative pose 선정 규칙

각 pocket에 대해 대표 pose 하나를 자동 선정해야 한다. 이는 report와 후속 PyMOL 확인에 중요하다.

## 28.1 기본 규칙

다음 우선순위를 추천한다.

1. 해당 pocket 내 가장 낮은 affinity pose
2. affinity가 같다면 contact residue 수가 더 많은 pose
3. 그래도 같다면 pose\_rank가 더 높은 pose

## 28.2 저장 정보

- representative\_pose\_rank
- representative\_ligand\_id
- representative\_pose\_file\_path

## 28.3 주의

대표 pose는 pocket의 “가장 진실한 구조”가 아니라 시각화용 대표값이므로, pocket summary에서는 top residues와 residue frequency가 더 우선이다.

---

# 29. report 문서의 세부 항목 설계

코덱스는 markdown report를 만들 때 아래 항목을 필수로 포함해야 한다.

## 29.1 `vina_report.md`

### 1절. 실행 개요

- 날짜
- receptor 수
- ligand 수
- Vina 파라미터
- contact cutoff
- pocket clustering cutoff

### 2절. receptor별 pocket 수 요약

예:

- 3GT8\_raw: 4 pockets
- 3GT8\_cl38\_48: 5 pockets
- 3GT8\_cl85\_100: 3 pockets

### 3절. receptor별 pocket summary 표

- pocket ID
- n\_pose
- n\_ligand
- best affinity
- top residues

### 4절. ligand별 dominant pocket 표

- ligand ID
- 각 receptor에서의 dominant pocket
- best affinity

### 5절. receptor 간 pocket overlap 상위 후보

- top 10 same\_patch\_candidate pairs

### 6절. 주의사항

- 공통 잔기 자동 해석의 한계
- centroid cutoff는 pocket 정의의 편의 기준일 뿐이라는 설명
- receptor 구조 상태 차이로 인한 pocket 이동 가능성

---

## 29.2 `integrated_report.md`

이 리포트는 Vina, PyRosetta, AFM 결과를 통합하는 요약 문서다.

### 포함 요소

1. 프로젝트 개요
2. receptor ensemble 설명
3. Vina pocket landscape 요약
4. PyRosetta receptor-side interface residue hotspot 요약
5. AFM receptor-side contact residue 요약
6. 서로 교차 등장하는 residue patch 후보 요약
7. 다음 수동 검토 우선순위 제안

### 주의

이 문서는 자동 생성되더라도, 최종 결론문이 아니라 **검토용 통합 브리핑 문서**여야 한다.

---

# 30. PyRosetta와 AFM의 역할을 과대평가하지 않기 위한 주석 규칙

현재 사용자는 가장 우선적으로 보고 싶은 것은 Vina 기반 pocket 분석이다. PyRosetta와 AFM은 구조 참고자료이자 residue patch 해석의 보조축이다. 따라서 report 문구도 다음 원칙을 따라야 한다.

## 30.1 금지할 표현

- “AFM이 정답 구조를 제시했다” 같은 표현
- “PyRosetta 결과가 pocket을 확정한다” 같은 표현

## 30.2 권장 표현

- “reference interface residues”
- “supporting receptor-side residue patch candidates”
- “auxiliary comparison set”
- “complementary structural hypothesis output”

즉, 이 프로그램들은 보조적 구조 가설 제공 도구이지, 메인 포켓 판정 엔진이 아니다.

---

# 31. 수동 검토를 위한 파일 자동 생성 요구사항

사용자는 최종적으로 PyMOL이나 Chimera 등에서 몇 개 구조를 직접 확인할 가능성이 높다. 따라서 코덱스는 수동 검토용 파일도 자동으로 만들어야 한다.

## 31.1 receptor별 top pocket pose export

각 receptor에서 pocket별 대표 pose를 별도 디렉터리에 export

예:

```text
exports/manual_review/
├── 3GT8_raw/
│   ├── P1_rep_drugA_pose1.pdbqt
│   ├── P2_rep_drugC_pose2.pdbqt
│   └── ...
├── 3GT8_cl38_48/
└── 3GT8_cl85_100/
```

## 31.2 PyMOL selection snippet 자동 생성

가능하면 report와 함께 다음 형식의 selection snippet도 생성한다.

```python
select pocket_P1_res, resi 971+972+973+974
show sticks, pocket_P1_res
```

### 목적

사용자가 residue patch를 수동 시각화하기 쉽게 함

---

# 32. configuration 설계 상세

`project.yaml`은 반드시 확장 가능하게 설계한다.

권장 필드:

```yaml
project_name: egfr_myo1d_pipeline
output_root: ./project_root
random_seed: 20260308
receptors:
  - id: 3GT8_raw
    pdb: receptors/3GT8_raw/3GT8_raw.pdb
    pdbqt: receptors/3GT8_raw/3GT8_raw.pdbqt
    chain: A
    source_type: raw_3GT8
  - id: 3GT8_cl38_48
    pdb: receptors/cluster_38_48/rep_38_48.pdb
    pdbqt: receptors/cluster_38_48/rep_38_48.pdbqt
    chain: A
    source_type: md_cluster_38_48
  - id: 3GT8_cl85_100
    pdb: receptors/cluster_85_100/rep_85_100.pdb
    pdbqt: receptors/cluster_85_100/rep_85_100.pdbqt
    chain: A
    source_type: md_cluster_85_100
ligands:
  - id: drugA
    pdbqt: ligands/drugA.pdbqt
    annotation_myo1d_inhibition: strong
    annotation_cell_viability: high_killing
    priority_tier: tier1
  - id: drugB
    pdbqt: ligands/drugB.pdbqt
    annotation_myo1d_inhibition: medium
    annotation_cell_viability: moderate
    priority_tier: tier2
vina:
  center: [0.0, 0.0, 0.0]
  size: [30.0, 30.0, 30.0]
  exhaustiveness: 32
  num_modes: 20
  energy_range: 6
contacts:
  cutoff: 4.0
  heavy_atom_only: true
pocket_clustering:
  method: dbscan
  centroid_cutoff: 4.0
  min_samples: 1
pocket_overlap:
  same_patch_centroid_cutoff: 6.0
  same_patch_jaccard_cutoff: 0.25
reports:
  export_markdown: true
  export_json: true
  export_manual_review_files: true
```

---

# 33. 코드 품질 요구사항

코덱스는 빠르게 동작하는 스크립트만 만들지 말고, 장기 유지보수 가능한 코드로 정리해야 한다.

## 33.1 함수 책임 분리

나쁜 예:

- 하나의 함수가 Vina 실행, parsing, contact 계산, pocket clustering, report 생성까지 다 하는 것

좋은 예:

- 실행
- 파싱
- contact 계산
- clustering
- summary 생성
- report 생성 각 단계가 분리된 함수/모듈로 존재

## 33.2 타입 힌트

가능한 모든 공개 함수에 타입 힌트 추가

## 33.3 예외 처리

- 파일 없음
- 비어 있는 output
- malformed PDBQT
- residue numbering mismatch
- config key 누락 이런 경우 명확한 예외 메시지와 로그를 남길 것

## 33.4 테스트 가능성

핵심 함수는 순수 함수에 가깝게 설계해 단위 테스트가 가능해야 한다.

---

# 34. acceptance criteria

코덱스는 아래 조건을 만족해야 완료로 본다.

## 34.1 최소 기능 완료 기준

1. receptor 3개와 ligand 여러 개를 넣어 Vina batch가 정상 실행된다.
2. `vina_pose_table.csv`가 생성된다.
3. `vina_pocket_table.csv`가 생성된다.
4. `vina_drug_pocket_map.csv`가 생성된다.
5. `vina_pocket_overlap_table.csv`가 생성된다.
6. `vina_report.md`가 생성된다.
7. residue numbering mismatch 시 경고 또는 오류가 뜬다.
8. config 기반 실행이 가능하다.

## 34.2 품질 기준

1. 같은 입력으로 재실행 시 같은 pocket assignment가 나온다.
2. output 파일 경로가 예측 가능하다.
3. 로그가 사람이 읽을 수 있다.
4. CSV가 UTF-8로 저장된다.
5. report가 빈 섹션 없이 생성된다.

## 34.3 확장성 기준

1. ligand 수가 늘어나도 batch로 처리 가능하다.
2. receptor 구조가 3개보다 많아져도 config만 바꾸면 처리 가능하다.
3. contact cutoff와 centroid cutoff를 config에서 쉽게 조정할 수 있다.

---

# 35. migration plan

기존 코드베이스를 바로 폐기하지 말고 아래 순서로 이행한다.

## 35.1 1단계

기존 Vina 실행 스크립트 보존 + 새 parser/summarizer 추가

## 35.2 2단계

기존 출력 포맷과 새 출력 포맷 비교 검증

## 35.3 3단계

새 output 체계가 안정되면 기존 수동 정리 방식 제거

## 35.4 4단계

PyRosetta와 AFM도 같은 방식으로 output 표준화

---

# 36. 사용자가 코덱스에게 붙여넣기 좋은 작업 분할 요청 예시

## 작업 묶음 A: Vina parser 구현

“기존 AutoDock Vina 실행 스크립트는 유지하고, 결과 PDBQT를 파싱해서 pose별 affinity, centroid, receptor contact residues를 추출하는 `parse_vina_results.py`와 `extract_contacts.py`를 새로 만들어줘. receptor는 3GT8\_raw, 3GT8\_cl38\_48, 3GT8\_cl85\_100 세 개를 config로 받게 하고, output은 `vina_pose_table.csv`로 저장되게 해줘.”

## 작업 묶음 B: pocket clustering 구현

“`vina_pose_table.csv`를 읽어서 receptor별 centroid 기반 4.0 Å clustering을 수행하고 `vina_pocket_table.csv`와 `vina_drug_pocket_map.csv`를 생성하는 스크립트를 만들어줘. pocket별 union residues와 residue frequency를 꼭 포함해줘.”

## 작업 묶음 C: receptor 간 pocket 비교 구현

“`vina_pocket_table.csv`를 바탕으로 receptor 간 pocket overlap을 비교하는 `compare_pockets.py`를 만들어줘. centroid distance, residue overlap, jaccard index, shared ligands를 계산하고 `vina_pocket_overlap_table.csv`로 저장해줘.”

## 작업 묶음 D: report 생성 구현

“위 CSV 파일들을 읽어서 receptor별 pocket summary, ligand별 dominant pocket, receptor 간 overlap top candidates를 포함한 `vina_report.md`를 생성하는 스크립트를 만들어줘.”

## 작업 묶음 E: 통합 실행기 구현

“config 파일 하나를 입력받아 `run_vina_batch.py`, `parse_vina_results.py`, `extract_contacts.py`, `cluster_pockets.py`, `summarize_pockets.py`, `compare_pockets.py`, `build_vina_report.py`를 순서대로 호출하는 `run_pipeline.py`를 만들어줘.”

---

# 37. 개발 시 피해야 할 흔한 실수

1. receptor별가 아니라 ligand별로 pocket clustering하는 것
2. pocket ID를 receptor 전역이 아니라 전체 프로젝트 전역으로 섞는 것
3. affinity만 보고 pocket를 대표하는 것
4. contact residue를 centroid 근처 residue로 대충 대체하는 것
5. residue numbering mismatch를 무시하는 것
6. report 없이 CSV만 남기는 것
7. 기존 보고서의 old site 이름을 결과 파일명에 하드코딩하는 것
8. pocket overlap을 boolean만 저장하고 raw metric을 버리는 것
9. 수동 검토용 representative pose export를 생략하는 것
10. config 없이 상수값을 코드에 박아넣는 것

---

# 38. 향후 확장 가능성 메모

현재는 아래 기능이 후순위지만, 구조를 잘 잡아두면 이후 쉽게 붙일 수 있다.

1. MM/GBSA rescoring hook
2. short MD validation hook
3. PLIP interaction typing
4. residue heatmap figure generation
5. pocket graph visualization
6. web UI 또는 notebook summary
7. residue-level ranking score 계산기
8. 약물 효과 annotation과 pocket 패턴 연계 분석

이 기능들은 지금 당장 필수는 아니지만, 코드 구조는 나중에 붙이기 쉽게 열어두는 것이 좋다.

---

# 39. 최종 개발 목표의 다시 한 번 명확한 정의

이 프로젝트의 목표는 “도킹을 더 많이 돌리는 것”이 아니다. 목표는 아래와 같다.

- 세 receptor 상태(3GT8 원본, 38–48 대표, 85–100 대표)에 대해
- 여러 ligand docking 결과를
- pose, pocket, receptor-overlap 수준으로 구조화하고
- 사람이 보기 쉬운 리포트와 재분석 가능한 표준 출력으로 정리하며
- PyRosetta/AFM 결과를 residue-level 참고 자료로 함께 관리하는 것

즉, 이 프로젝트는 **자동 도킹 실행기 개발**이 아니라 **연구용 구조 데이터 인프라 구축**에 가깝다.

---

# 40. 코덱스에게 전달할 최종 마감 문구

이 문서는 단순한 아이디어 메모가 아니라, 현재 진행 중인 EGFR–MYO1D 연구의 계산 파이프라인을 실제로 재정비하기 위한 명세서다. 구현 시 가장 중요한 것은 새로 생성되는 데이터를 우선시하고, 결과를 비교 가능한 표준 형식으로 출력하며, 사용자가 pocket 위치와 residue patch를 직접 판단할 수 있도록 raw metric을 충분히 보존하는 것이다. 기존 보고서나 과거 pocket 이름은 절대 기준이 아니다. 지금부터 만들어지는 표준화된 output이 이후 연구 해석의 기준이 된다. 따라서 코드의 목적은 정답을 미리 가정하는 것이 아니라, **판단 가능한 구조 데이터셋을 안정적으로 생산하는 것**이다.

이 문서를 기준으로 작업을 분할하고, 먼저 Vina 계열 파서와 pocket summary 체계를 완성한 뒤, PyRosetta와 AFM residue output을 같은 철학으로 정리하라. 최종 목표는 사용자가 receptor 상태별 pocket landscape와 ligand 배정 패턴을 한눈에 읽을 수 있도록 만드는 것이다.

