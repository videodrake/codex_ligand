# EGFR-MYO1D PPT 중학생용 쉬운 설명 텍스트/데이터

이 문서는 PPT에 들어간 말을 더 쉽게 고치기 위한 참고 문서입니다. 발표에서는 "무엇을 했는지", "무엇을 아직 안 했는지", "다음에 무엇을 하면 되는지"가 먼저 보이게 설명하면 됩니다.

## 공통 설명

- EGFR은 세포막에 있는 단백질입니다.
- MYO1D가 EGFR에 붙는지 보려면, 먼저 EGFR의 모양을 정해야 합니다.
- 이번 발표는 MYO1D 결과 발표가 아니라, MYO1D 분석을 시작하기 전 EGFR 모델과 대표 구조를 준비한 진행 보고입니다.
- 데이터가 완벽하지 않아도, "여기까지 했다"와 "아직 확인할 부분이 있다"를 분명히 말하면 됩니다.

## 1장. 표지

### 슬라이드에 들어간 말

- 제목: EGFR 모델을 먼저 만들고 확인했습니다. 대표 구조도 정했습니다
- 한 줄 설명: MYO1D가 EGFR에 붙는지 보려면, 먼저 EGFR의 모양과 방향을 정해야 합니다.
- 범위: 이번 발표는 EGFR 대표 구조를 정하고 MYO1D 분석으로 넘어가는 준비 상황입니다.

### 데이터/자료

- 표지 이미지는 발표용 개념 그림입니다.
- 파일: `reports/assets/slide01_cover_background.png`

## 2장. 왜 EGFR 모델부터 정하나?

### 쉬운 설명

EGFR 모양이 바뀌면 MYO1D가 붙는 위치도 달라 보일 수 있습니다. 그래서 MYO1D를 붙여보기 전에 기준이 될 EGFR 모델을 먼저 정했습니다.

### 발표할 내용

- EGFR의 여러 부분을 이어 붙여 기본 모델을 만들었습니다.
- +10도 방향 모델이 현재 지표에서 가장 좋아 보였습니다.
- Cluster 1/2 대표 구조를 다음 입력 후보로 정했습니다.

### 아직 다음 단계

- MYO1D를 붙여본 최종 결과는 아직 아닙니다.
- 붙는 자리, 작은 홈, 화합물 docking은 다음 단계입니다.

## 3장. EGFR 모델을 만든 방법

### 쉬운 설명

EGFR은 한 번에 완성된 구조 자료가 부족해서, 여러 구조 자료를 이어 붙이고 빈 곳은 계산으로 채웠습니다.

### 데이터/자료

- 막을 지나는 부분: `2M0B`
- 막 근처 연결부: `2M20`
- 세포 안쪽 효소부: `3GT8`
- 빈 구간 채우기: `MODELLER`
- 참고 파일: `1.align/*`, `01_numbering_provenance/template_model_mapping.csv`

### 주의할 점

이 모델은 실험으로 바로 찍은 완성 구조가 아니라, 기존 자료와 계산을 합쳐 만든 모델입니다.

## 4장. +10도 방향 모델 선택

### 쉬운 설명

EGFR 안쪽 큰 부분을 여러 각도로 돌려 보면서, 어느 방향이 더 안정적으로 보이는지 비교했습니다. 현재 가진 지표에서는 +10도 모델이 가장 좋은 후보로 보였습니다.

### 발표 숫자

- 움직임 크기: `2.68 Å`
- 전체 크기: `28.31 Å`
- 두 가닥 연결: `11.49`

### 데이터/자료

- `02_orientation_metrics/orientation_comparison_table.csv`
- `reports/assets/slide04_rmsd_by_orientation.png`
- `reports/assets/slide04_rg_by_orientation.png`
- `reports/assets/slide04_hbond_by_orientation.png`
- `reports/assets/slide04_orientation_score_summary.png`

## 5장. 200 ns 움직임 기록 확인

### 쉬운 설명

trajectory는 단백질이 시간에 따라 어떻게 움직였는지 저장한 기록입니다. 이번에는 +10도 EGFR 모델이 200 ns 동안 움직인 기록을 확인했습니다.

### 발표 숫자

- 원자 수: `49,854`
- 저장 장면 수: `2,001`
- 확인한 시간: `0-200 ns`
- 검사 결과: `OK`

### 데이터/자료

- `gmx_check_EGFR_plus10_step7_200ns_nw.txt`
- `EGFR_plus10_step7_production_nw.gro`
- `EGFR_plus10_step7_200ns_nw.xtc`
- `reports/assets/slide05_kd_membrane_distance_200ns.png`
- `reports/assets/slide05_interchain_contacts_200ns.png`
- `reports/assets/slide05_tm_kd_tilt_200ns.png`

## 6장. EGFR 구조 그림

### 쉬운 설명

실제 PDB 좌표를 이용해서 EGFR 두 가닥이 어떻게 놓여 있는지 간단히 그렸습니다. 논문용 고품질 그림은 아니지만, 발표에서 현재 진행 상황을 설명하기에는 충분합니다.

### 발표 숫자

- 막과의 거리: `25.66 Å`
- 두 가닥 접촉: `505.28`

### 데이터/자료

- `03_receptor_pack/receptor/EGFR_plus10_final_model_EGFR_rot10_best.pdb`
- `reports/assets/slide6_plus10_side_view.png`
- `reports/assets/slide6_plus10_top_view.png`
- 생성 스크립트: `reports/ppt_src/make_structure_overview.py`

## 7장. 대표 구조 선택

### 쉬운 설명

clustering은 비슷한 구조끼리 묶어서, 그중 가장 대표적인 구조를 고르는 과정입니다. 이번에는 Cluster 1과 Cluster 2가 나오는 구간을 다음 MYO1D docking 입력으로 사용하기로 했습니다.

### 발표 숫자

- 선택 구간: `70-100 ns`
- cutoff: `0.173 nm`
- Cluster 1: `121/301장`, `40.2%`, 대표 시간 `79.3 ns`
- Cluster 2: `68/301장`, `22.6%`, 대표 시간 `93.5 ns`

### 데이터/자료

- `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173/clustering_summary.txt`
- `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173/cluster_population.png`
- `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173/cluster_timeline.png`
- C1 구조 파일: `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173/EGFR_170-200.pdb`
- C2 구조 파일: `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173/cluster_2_representative.pdb`

### 발표 표현

쉽게 말하면, EGFR 대표 구조는 C1 하나만 보지 않고 C1과 C2 두 모양을 다음 MYO1D docking 입력 후보로 가져갑니다.

## 8장. 완료한 일과 남은 일

### 지금 말할 수 있음

- EGFR 기본 모델을 만들었습니다.
- +10도 방향 모델을 현재 후보로 골랐습니다.
- 200 ns 움직임 기록 파일을 확인했습니다.
- Cluster 1/2 대표 구조를 다음 입력 후보로 정했습니다.

### 아직 확인 필요

- C1/C2 구조 파일을 MYO1D docking 입력 형식으로 정리합니다.
- 대표 구조별로 MYO1D docking 결과를 비교해야 합니다.
- 더 예쁜 구조 그림은 나중에 PyMOL/ChimeraX로 만들 수 있습니다.
- 일부 연결부 제작 과정은 추가 확인이 필요합니다.

### 다음 단계

- MYO1D를 EGFR에 붙여보는 docking
- MYO1D가 닿는 EGFR 위치 찾기
- 약물이 들어갈 만한 작은 홈 찾기
- 대표 구조별 docking 결과를 비교한 뒤 화합물 docking

## 9장. 다음 결정

### 교수님께 확인할 질문

다음 분석에서는 C1과 C2 대표 구조에 MYO1D를 붙여보고, 어느 위치가 반복해서 나오는지 비교합니다.

- Cluster 1은 가장 많이 나온 대표 모양입니다.
- Cluster 2는 두 번째로 많이 나온 대표 모양입니다.
- 두 구조를 모두 docking input 후보로 사용합니다.

### 다음 작업 순서

1. C1/C2 구조 준비
2. MYO1D 붙여보기
3. 닿는 자리 정리하기
4. 작은 홈 찾기
5. 화합물 넣어보기

## 현재 PPT 파일

- C1/C2 선택 반영 PPTX: `reports/ppt_output/EGFR_MYO1D_receptor_modeling_clustering_deck_cluster12_selected.pptx`
- 진행본 PPTX: `reports/ppt_output/EGFR_MYO1D_progress_until_slide9.pptx`
- 미리보기 PNG 폴더: `reports/ppt_output/previews/`
