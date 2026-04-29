# EGFR-MYO1D PPT 슬라이드별 텍스트와 데이터 자료

목적: 교수님께 “현재 어디까지 진행했는지”를 발표하기 위한 슬라이드별 원고/자료 정리이다. 데이터가 최종 논문 수준으로 완벽하다는 의미가 아니라, 현재 폴더 안에서 확인되는 진행 상황과 사용 가능한 근거를 발표용으로 묶은 것이다.

## Slide 1. Title

### 슬라이드 텍스트

- Membrane-compatible EGFR dimer receptor model 구축 및 clustering
- 후속 MYO1D PPI 분석 전, EGFR receptor input을 +10° model과 final 200 ns MD validation 기준으로 정리
- Scope: PPI docking, pocket discovery, compound docking은 이번 발표에서 결과로 다루지 않음

### 사용 자료

- `reports/assets/slide01_cover_background.png`
- `reports/EGFR_modeling_clustering_PPT_data_package.md`

### 발표 멘트

이번 발표는 MYO1D docking 결과가 아니라, 그 전에 필요한 EGFR receptor input을 어디까지 정리했는지 보고하는 것이다. 현재는 receptor modeling, orientation scan, final 200 ns validation, clustering 검토 단계까지 진행했다.

## Slide 2. Why receptor model first?

### 슬라이드 텍스트

- MYO1D PPI docking 전에 EGFR receptor frame을 먼저 고정해야 downstream 해석이 흔들리지 않는다.
- 현재 deck에서 다루는 확정 범위:
  - template assembly provenance
  - +10° orientation 선택 근거
  - final no-water 200 ns trajectory identity
- 결과로 주장하지 않는 범위:
  - MYO1D docking result
  - EGFR-side PPI patch / pocket discovery
  - compound docking

### 사용 자료

- `reports/EGFR_modeling_clustering_PPT_data_package.md`

### 발표 멘트

후속 PPI 분석을 하려면 먼저 EGFR receptor frame이 고정되어야 한다. 그래서 이번에는 docking 결과보다 receptor 준비 상태를 먼저 보고하고, downstream 단계는 다음 작업으로 둔다.

## Slide 3. Composite model construction

### 슬라이드 텍스트

- TM-JM-KD inactive dimer model assembly
- 2M0B TM dimer, 2M20 JM-A reference, 3GT8 inactive KD dimer, MODELLER gap filling을 조합
- 구성 요소:
  - TM dimer: 2M0B, local script retained 634-670
  - JM-A reference: 2M20, membrane-proximal JM-A orientation
  - inactive KD dimer: 3GT8, symmetric inactive kinase-domain dimer
  - modeled gaps: JM-B 684-704 and A-loop 862-876

### 사용 자료

- `1.align/2m0b.cif`
- `1.align/2m20.cif`
- `1.align/3gt8.cif`
- `1.align/assemble_egfr.py`
- `1.align/assemble_tmjma_kd.py`
- `1.align/fill_tmjma_gaps.py`
- `01_numbering_provenance/template_model_mapping.csv`
- `1.align/final_models/EGFR_rot10_best.pdb`

### 발표 멘트

EGFR receptor는 단일 template이 아니라 TM, JM, KD 영역을 조합해서 만들었다. gap은 MODELLER로 채웠고, 이후 orientation scan을 위해 여러 rotational variant를 만들었다.

## Slide 4. Rotational orientation scan

### 슬라이드 텍스트

- Five TM-KD rotations screened; +10° ranked best
- +10° selected metric snapshot:
  - RMSD mean: 2.68 Å ± 0.39
  - Rg mean: 28.31 Å ± 0.12
  - Inter-chain H-bond: 11.49 ± 3.11 where available
- Available parsed metrics 기준으로 +10°를 working receptor orientation으로 선택

### 사용 자료

- `02_orientation_metrics/orientation_comparison_table.csv`
- `02_orientation_metrics/orientation_metric_summary.md`
- `reports/assets/slide04_rmsd_by_orientation.png`
- `reports/assets/slide04_rg_by_orientation.png`
- `reports/assets/slide04_hbond_by_orientation.png`
- `reports/assets/slide04_orientation_score_summary.png`

### 발표 멘트

-20, -10, 0, +10, +20 다섯 orientation을 비교했다. 완전히 모든 지표가 다 갖춰진 것은 아니지만, 현재 파싱된 RMSD/Rg/H-bond/score 기준으로 +10°가 가장 좋아 보여서 working receptor로 선택했다.

## Slide 5. MD validation result

### 슬라이드 텍스트

- Final +10 no-water 200 ns trajectory is identified
- Final pair:
  - `EGFR_plus10_step7_production_nw.gro`
  - `EGFR_plus10_step7_200ns_nw.xtc`
- gmx_check:
  - 49,854 atoms
  - 2,001 coordinate frames
  - 100 ps output interval
  - 0-200 ns
- analysis status: OK

### 사용 자료

- `03_receptor_pack/receptor/EGFR_plus10_step7_production_nw.gro`
- `03_receptor_pack/trajectory/EGFR_plus10_step7_200ns_nw.xtc`
- `03_receptor_pack/metadata/gmx_check_EGFR_plus10_step7_200ns_nw.txt`
- `reports/assets/slide05_kd_membrane_distance_200ns.png`
- `reports/assets/slide05_interchain_contacts_200ns.png`
- `reports/assets/slide05_tm_kd_tilt_200ns.png`

### 발표 멘트

최종적으로 사용할 +10 no-water trajectory pair를 정리했고, gmx_check에서 200 ns trajectory identity를 확인했다. 이 trajectory에서 KD-membrane distance, interchain contacts, TM-KD tilt 같은 지표를 계산했다.

## Slide 6. Structural interpretation

### 슬라이드 텍스트

- +10° receptor frame is the downstream input candidate
- 현재 단계에서는 membrane-compatible inactive dimer frame으로 해석
- 150-200 ns KD distance: 25.66 Å ± 0.86
- 150-200 ns contacts: 505.28 ± 38.10 heavy contacts within 4.5 Å
- Caveat: GRO chain ID가 없어 protomer는 residue order로 split

### 사용 자료

- `03_receptor_pack/receptor/EGFR_plus10_final_model_EGFR_rot10_best.pdb`
- `03_receptor_pack/receptor/EGFR_plus10_step7_production_nw.gro`
- `03_receptor_pack/metadata/membrane_stability_plus10_last50ns_150_200ns_mdanalysis.csv`
- `reports/assets/slide6_plus10_side_view.png`
- `reports/assets/slide6_plus10_top_view.png`
- 생성 스크립트: `reports/ppt_src/make_structure_overview.py`

### 발표 멘트

전문 구조 렌더는 아니지만, 실제 PDB CA 좌표를 이용해 side/top overview를 만들었다. 발표에서는 이 이미지를 “현재 구조 해석용 overview”로만 사용하고, 필요하면 이후 PyMOL/ChimeraX로 publication-quality figure를 만들면 된다.

## Slide 7. Trajectory clustering

### 슬라이드 텍스트

- Older representative exists; final-window reps need confirmation
- 80-100 ns clustering:
  - C1 = 150/201 frames
  - C1 population = 74.6%
  - representative time = 84.0 ns
- 160-185 ns C1: missing / not confirmed
- 170-200 ns C1/C2: missing / not confirmed
- `EGFR_160-185.pdb`, `EGFR_170-200.pdb`: file exists but interval provenance ambiguous

### 사용 자료

- `4.+10_clustering/cluster_results/80-100ns_cut0.15/clustering_summary.txt`
- `4.+10_clustering/cluster_results/80-100ns_cut0.15/cluster_population.png`
- `4.+10_clustering/cluster_results/80-100ns_cut0.15/cluster_timeline.png`
- `reports/assets/slide07_cluster_population_80_100.png`
- `reports/assets/slide07_cluster_timeline_80_100.png`
- `reports/clustering_file_inventory.txt`
- `reports/clustering_160_200_search_hits.txt`

### 발표 멘트

clustering은 진행되어 있고 80-100 ns 대표 구조는 확인된다. 다만 교수님께 요청받은 160-185 또는 170-200 ns representative는 파일명만으로 확정하기에는 provenance가 애매해서, 현재는 needs confirmation으로 보고하는 것이 안전하다.

## Slide 8. Current confirmed status

### 슬라이드 텍스트

- Confirmed now:
  - EGFR TM-JM-KD inactive dimer composite model constructed
  - +10° current receptor input candidate
  - final no-water trajectory identity confirmed
  - final-pair membrane-stability metrics available
- Needs confirmation:
  - final-window clustering representative
  - ambiguous representative PDB filenames
  - publication-quality structure render
  - JM-A coordinate-generation provenance
- Future work only:
  - MYO1D docking
  - EGFR-side PPI patch
  - PPI-adjacent pocket discovery
  - compound docking

### 사용 자료

- `reports/EGFR_modeling_clustering_PPT_data_package.md`
- `reports/다음할일.md`

### 발표 멘트

이 장에서는 교수님께 현재 확정해서 말할 수 있는 것과 아직 확인이 필요한 것, 그리고 다음 단계로 남겨둘 것을 한 번에 보여준다.

## Slide 9. Next decision

### 슬라이드 텍스트

- Confirm receptor representative, then start MYO1D PPI work
- Feedback request:
  - Which receptor state should become the MYO1D docking input?
  - Pre-MD final +10 PDB vs final 200 ns coordinate vs confirmed cluster representative
  - Use 150-200 ns, 160-200 ns, or another final-window clustering interval?
  - Redefine positive-patch residue labels using final topology residue identity
- Sequential next work:
  1. Confirm representative receptor state
  2. Run/finalize MYO1D docking
  3. Extract EGFR-side PPI patch
  4. Identify PPI-adjacent pocket
  5. Compound docking

### 사용 자료

- `reports/다음할일.md`
- `reports/EGFR_modeling_clustering_PPT_data_package.md`

### 발표 멘트

교수님께 받을 핵심 결정은 후속 docking input으로 어떤 receptor state를 사용할지이다. 이 결정 후 MYO1D docking, PPI patch, pocket, compound 순서로 넘어가면 된다.

## 최종 산출물 경로

- PPTX: `reports/ppt_output/EGFR_MYO1D_receptor_modeling_clustering_deck.pptx`
- PNG preview: `reports/ppt_output/previews/`
- PPT build script: `reports/ppt_src/build_egfr_deck.mjs`
- Structure overview image script: `reports/ppt_src/make_structure_overview.py`
