# EGFR-MYO1D PPT 이미지/페이지 제작 계획

목표: 교수님 발표용 9장 PPT를 한 장씩 제작한다. 과학적 근거가 필요한 구조/trajectory/metric 이미지는 실제 데이터 경로를 사용하고, GPT Image 2 계열 생성 이미지는 표지, 개념 배경, 비데이터성 보조 시각화에만 사용한다.

## 공통 원칙

- 실제 결과를 보여주는 이미지는 AI 생성 이미지로 대체하지 않는다.
- PDB/GRO/XTC/CSV/XVG/log 기반 claim은 반드시 해당 경로를 slide note 또는 보고서 provenance에 남긴다.
- `3.rotate_EGFR/10/full_200ns.xtc`는 final evidence로 사용하지 않는다.
- 160-185 ns / 170-200 ns clustering representative는 현재 `ambiguous / needs confirmation`이다.
- MYO1D docking, pocket, compound docking은 future work로만 언급한다.
- 생성 이미지는 “분자 구조 결과”가 아니라 “conceptual cover/background”로만 사용한다.

## Deck Style

- Audience: 교수님께 진행상황 보고
- Tone: technical, sober, high-confidence but caveated
- Visual style: clean scientific presentation, white/ink background with cyan/green/orange accents
- Avoid: 과장된 신약개발/compound 느낌, fake molecular result, generic consulting card grid
- Main visual language: open typography, real figures, native tables, restrained generated scientific background

## Slide 1. Title

### Message

EGFR-MYO1D PPI 분석 전, membrane-compatible EGFR dimer receptor input을 +10° model과 final MD validation 기준으로 정리했다.

### Primary visual

GPT Image 2 generated cover background.

### Image prompt draft

Use case: scientific-educational
Asset type: PowerPoint title slide background
Primary request: premium scientific cover image for an EGFR membrane receptor modeling presentation, showing an abstract membrane bilayer horizon with a subtle dimeric receptor silhouette, no specific molecular structure claims
Style/medium: polished biomedical editorial illustration, realistic but abstract, not a literal protein model
Composition/framing: 16:9 wide, large negative space on left for Korean title text, receptor silhouette on right, membrane plane running horizontally
Lighting/mood: calm, rigorous, deep white-to-soft-cyan scientific lighting
Color palette: white, charcoal, cyan, muted green, restrained orange accents
Constraints: no labels, no readable text, no logos, no drug molecules, no docking pose, no fake detailed protein structure, no watermark

### Actual data needed

- none for generated background

### Existing data references

- `reports/EGFR_modeling_clustering_PPT_data_package.md`
- `03_receptor_pack/receptor/EGFR_plus10_final_model_EGFR_rot10_best.pdb` for optional later structure render, not for generated cover

### Status

Ready to generate cover image.

## Slide 2. Why receptor model first?

### Message

후속 MYO1D PPI docking 전에 EGFR receptor frame을 먼저 확정해야 한다.

### Primary visual

Native PPT schematic: receptor input 확정 -> MYO1D docking -> PPI patch -> pocket -> compound docking.

### GPT Image 2 use

Optional small conceptual background only. Not required.

### Actual data needed

- none

### Existing data references

- `reports/EGFR_modeling_clustering_PPT_data_package.md`

### Missing data to add later

- none

## Slide 3. Composite model construction

### Message

2M0B TM dimer, 2M20 JM-A reference, 3GT8 inactive KD dimer, MODELLER gap filling으로 EGFR TM-JM-KD composite model을 만들었다.

### Primary visual

Template/workflow table + optional real structure render.

### GPT Image 2 use

No AI-generated molecular result. Use native table/diagram.

### Actual data needed

- `1.align/2m0b.cif`
- `1.align/2m20.cif`
- `1.align/3gt8.cif`
- `1.align/assemble_egfr.py`
- `1.align/assemble_tmjma_kd.py`
- `1.align/fill_tmjma_gaps.py`
- `1.align/final_models/EGFR_rot10_best.pdb`
- `01_numbering_provenance/template_model_mapping.csv`

### Missing data to add later

- If a publication-quality structure render is desired: rendered PNG from PyMOL/ChimeraX using `03_receptor_pack/figures/EGFR_plus10_domain_coloring.pml` or `.cxc`.

## Slide 4. Rotational orientation scan

### Message

Five TM-KD rotational variants were screened; +10° is top-ranked among available parsed metrics.

### Primary visual

Real orientation metric plots and a concise native comparison table.

### GPT Image 2 use

Not needed.

### Actual data needed

- `02_orientation_metrics/orientation_comparison_table.csv`
- `02_orientation_metrics/orientation_metric_summary.md`
- `02_orientation_metrics/figures/rmsd_by_orientation.png`
- `02_orientation_metrics/figures/rg_by_orientation.png`
- `02_orientation_metrics/figures/interchain_hbond_by_orientation.png`
- `02_orientation_metrics/figures/orientation_score_summary.png`

### Missing data to add later

- Optional: one composite image comparing five structural variants side-by-side.

## Slide 5. MD validation result: +10° selected

### Message

Final +10 no-water trajectory identity is confirmed by gmx_check, and final-pair 200 ns metrics are available.

### Primary visual

Large validation number strip + one or two final trajectory metric plots.

### GPT Image 2 use

Not needed.

### Actual data needed

- `03_receptor_pack/metadata/gmx_check_EGFR_plus10_step7_200ns_nw.txt`
- `03_receptor_pack/metadata/membrane_stability_plus10_200ns_mdanalysis.csv`
- `03_receptor_pack/metadata/membrane_stability_plus10_last50ns_150_200ns_mdanalysis.csv`
- `03_receptor_pack/figures/trajectory_metrics/plus10_200ns_kd_membrane_abs_distance_A.png`
- `03_receptor_pack/figures/trajectory_metrics/plus10_200ns_interchain_heavy_contacts_4p5A.png`
- `03_receptor_pack/figures/trajectory_metrics/plus10_200ns_tm_kd_tilt_angle_deg.png`

### Missing data to add later

- none for current draft

## Slide 6. Structural interpretation of +10° model

### Message

+10° model preserves a membrane-compatible inactive dimer frame suitable as downstream receptor input candidate.

### Primary visual

Real structural render: side view/top view with membrane plane.

### GPT Image 2 use

No fake molecular result. Use actual structure render only.

### Actual data needed

- `03_receptor_pack/receptor/EGFR_plus10_final_model_EGFR_rot10_best.pdb`
- `03_receptor_pack/receptor/EGFR_plus10_step7_production_nw.gro`
- `03_receptor_pack/figures/EGFR_plus10_domain_coloring.pml`
- `03_receptor_pack/figures/EGFR_plus10_domain_coloring.cxc`

### Missing data to add later

- `reports/assets/slide6_plus10_side_view.png`
- `reports/assets/slide6_plus10_top_view.png`

If these are not available when making the first deck draft, use a placeholder panel with exact data path and render instructions.

## Slide 7. Trajectory clustering

### Message

Older clustering representative exists, but requested 160-185 / 170-200 ns representatives are not yet confirmed.

### Primary visual

Cluster population/timeline figure from older confirmed 80-100 ns output, plus a “needs confirmation” status table.

### GPT Image 2 use

Not needed.

### Actual data needed

- `4.+10_clustering/cluster_results/80-100ns_cut0.15/clustering_summary.txt`
- `4.+10_clustering/cluster_results/80-100ns_cut0.15/cluster_population.png`
- `4.+10_clustering/cluster_results/80-100ns_cut0.15/cluster_timeline.png`
- `03_receptor_pack/receptor/EGFR_plus10_cluster_rep_80-100ns_cut0.15_cluster1.pdb`

### Ambiguous files

- `3.rotate_EGFR/10/cluster_results/60-85ns_cut0.16/EGFR_160-185.pdb`
- `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173/EGFR_170-200.pdb`

### Missing data to add later

- Confirmed 160-185 ns C1 representative PDB and summary
- Confirmed 170-200 ns C1/C2 representative PDBs and summary

## Slide 8. Current confirmed status

### Message

확정된 receptor-modeling/MD validation 범위와 아직 미확정인 downstream work를 분리한다.

### Primary visual

Native table: Confirmed / Needs confirmation / Future work.

### GPT Image 2 use

Not needed.

### Actual data needed

- `reports/EGFR_modeling_clustering_PPT_data_package.md`
- `reports/다음할일.md`

### Missing data to add later

- none

## Slide 9. Next steps and feedback request

### Message

교수님께 final receptor state와 clustering window 사용 여부를 확인받고 다음 단계로 넘어간다.

### Primary visual

Native roadmap: receptor representative confirmation -> MYO1D docking -> EGFR-side PPI patch -> pocket -> compound docking.

### GPT Image 2 use

Optional subtle closing background, not necessary.

### Actual data needed

- none beyond report

### Missing data to add later

- 교수님 피드백 후 updated decision log

## Immediate build order

1. Generate Slide 1 cover background with GPT Image 2 / built-in image generation.
2. Save generated cover under `reports/assets/`.
3. Build Slide 1 as editable PPT page with the generated background and title text.
4. Render/check Slide 1.
5. Continue Slide 2 using native diagram.
6. Use real figures from the listed paths for Slide 4, 5, and 7.
7. Use placeholders with exact data paths for Slide 6 if structure renders are not yet available.
