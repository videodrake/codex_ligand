# EGFR-MYO1D receptor-modeling/clustering PPT data package

Generated: 2026-04-29
Dataset root: `C:\Users\dudwn\OneDrive\문서\EGFR_MYO1D\analysis_v0_5`

## 2026-04-30 slide update

사용자 결정에 따라 PPT 수정본에서는 EGFR 대표 구조를 `Cluster 1/2` 기반으로 정리했다.

- 선택 구간: `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173`
- Cluster 1: 121/301 frames, 40.2%, representative time 79.3 ns
- Cluster 2: 68/301 frames, 22.6%, representative time 93.5 ns
- C1 구조 파일: `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173/EGFR_170-200.pdb`
- C2 구조 파일: `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173/cluster_2_representative.pdb`
- 수정 PPT: `reports/ppt_output/EGFR_MYO1D_receptor_modeling_clustering_deck_cluster12_selected.pptx`

아래 기존 audit 메모는 2026-04-29 당시 파일명/summary 불일치를 조심하라는 기록으로 남겨둔다. 2026-04-30 PPT에서는 사용자의 대표 구조 결정에 맞춰 C1/C2 구간 사용으로 문구를 갱신했다.

## Executive status

이번 보고서는 교수님 발표용 PPT를 만들기 전, 현재 폴더에서 확인되는 EGFR receptor model 구축, orientation scan, 200 ns MD validation, clustering 근거를 정리한 데이터 패키지이다. 새 MD simulation은 수행하지 않았고, 파일/로그/CSV/기존 문서에서 확인되는 내용만 사용했다.

핵심 결론은 다음처럼 제한해서 말하는 것이 안전하다.

- Final receptor input으로는 `03_receptor_pack/receptor/EGFR_plus10_step7_production_nw.gro`와 `03_receptor_pack/trajectory/EGFR_plus10_step7_200ns_nw.xtc` 조합을 우선 사용한다.
- `gmx_check_EGFR_plus10_step7_200ns_nw.txt`에서 final no-water trajectory는 49,854 atoms, 2,001 coordinate frames, 100 ps timestep, 0-200,000 ps 범위로 확인된다.
- +10° orientation 선택은 기존 parsed orientation metrics에서 +10이 가장 favorable한 것으로 정리할 수 있다. 단, 이 orientation 비교 metric은 주로 50-100 ns window 기반이며, 모든 membrane/contact metric이 orientation별로 완비된 것은 아니다.
- Final 200 ns pair에 대한 새 membrane-stability rerun output은 `analysis_status=OK`이며, `NOT_AVAILABLE` 값이 없는 OK rows로 확인된다.
- 요청된 160-185 ns / 170-200 ns clustering 결과는 현재 파일명만으로 확정하면 안 된다. `EGFR_160-185.pdb`, `EGFR_170-200.pdb` 파일은 존재하지만, 주변 summary는 각각 60-85 ns, 70-100 ns clustering으로 기록되어 있어 `ambiguous / needs confirmation`으로 둔다.

## Provenance and generated inventories

작업 시작 후 생성한 provenance 파일:

- `reports/EGFR_modeling_clustering_reference_map.md`
- `reports/rotate_EGFR_file_inventory.txt`
- `reports/clustering_file_inventory.txt`
- `reports/clustering_160_200_search_hits.txt`
- `reports/figure_file_inventory.txt`

`reports/clustering_160_200_search_hits.txt`는 넓은 텍스트 검색 결과라 매우 크다. 실제 slide 판단에는 개별 summary 파일과 inventory를 우선 사용했다.

## A. 전체 데이터 인벤토리

### Final receptor / selected trajectory

| Category | File path | Role | Slide use | Status | Source type | Caveat |
|---|---|---|---|---|---|---|
| final receptor coordinate | `03_receptor_pack/receptor/EGFR_plus10_step7_production_nw.gro` | final +10 no-water coordinate/reference | Slide 5-6 receptor frame, final status | found | raw file | final coordinate is GRO, not PDB |
| final receptor PDB | `03_receptor_pack/receptor/EGFR_plus10_final_model_EGFR_rot10_best.pdb` | folder 1 final +10 model PDB | Slide 3-4 construction/orientation visual | found | raw file | pre-MD final model, not trajectory-derived representative |
| cluster representative PDB | `03_receptor_pack/receptor/EGFR_plus10_cluster_rep_80-100ns_cut0.15_cluster1.pdb` | existing 80-100 ns dominant cluster representative | Slide 7 supporting older clustering | found | raw file/documentation | not the requested 160-185 or 170-200 ns interval |
| CHARMM-GUI starting system | `03_receptor_pack/receptor/EGFR_plus10_charmm_step5_input.pdb` / `.gro` | membrane system starting structure | Slide 5/appendix setup provenance | found | raw file | starting system, not final MD coordinate |
| final 200 ns trajectory | `03_receptor_pack/trajectory/EGFR_plus10_step7_200ns_nw.xtc` | final no-water 200 ns +10 trajectory | Slide 5-7 validation/clustering source | found | raw file/log | validated by gmx_check log |
| canonical selected trajectory | `03_receptor_pack/trajectory/EGFR_plus10_selected.xtc` | canonical copy of final selected trajectory | provenance only | found | raw file/documentation | regular file, not symlink; target metadata not verifiable |
| last-50-ns trajectory | `03_receptor_pack/trajectory/EGFR_plus10_last50ns_150_200ns_nw.xtc` | last 50 ns no-water trajectory | possible final-window analysis source | found | raw file | no separate uploaded gmx_check log found |
| final validation log | `03_receptor_pack/metadata/gmx_check_EGFR_plus10_step7_200ns_nw.txt` | trajectory identity validation | Slide 5 evidence | found | log | confirms 49,854 atoms / 2,001 frames / 100 ps / 200 ns |

### Composite model construction files

| Category | File path | Role | Slide use | Status | Source type | Caveat |
|---|---|---|---|---|---|---|
| TM template | `1.align/2m0b.cif` | EGFR TM dimer template | Slide 3 template table | found | raw file/documentation | local script retained 634-670 segment |
| JM-A reference | `1.align/2m20.cif` | JM-A / membrane-proximal orientation reference | Slide 3 template table | found | raw file/documentation | exact coordinate-generation path for JM-A remains incomplete |
| inactive KD template | `1.align/3gt8.cif` | inactive symmetric kinase-domain dimer template | Slide 3 template table | found | raw file/documentation | PDB numbering and UniProt-like model numbering must not be merged blindly |
| assembly script | `1.align/assemble_egfr.py` | early TM-JM-KD assembly/provenance | Slide 3 workflow | found | script | source evidence, not necessarily full execution log |
| TM-JM-A/KD assembly script | `1.align/assemble_tmjma_kd.py` | TM-JM-A + KD assembly path | Slide 3 workflow | found | script | chain swap logic exists; avoid over-describing exact swap state |
| gap filling script | `1.align/fill_tmjma_gaps.py` | MODELLER gap filling for JM-B/A-loop | Slide 3 workflow | found | script | modeled gap, not template-derived coordinates |
| variant analysis | `1.align/variant_analysis.csv` | rotational geometry summary | Slide 4 orientation scan | found | CSV | initial structural metrics, not MD stability proof |
| final rotational PDBs | `1.align/final_models/EGFR_rot_20_best.pdb`, `EGFR_rot_10_best.pdb`, `EGFR_rot00_best.pdb`, `EGFR_rot10_best.pdb`, `EGFR_rot20_best.pdb` | five rotational variants | Slide 4 visual comparison | found | raw file | use +10 as selected, others as scan inputs |
| MODELLER outputs | `1.align/rot_outputs/rot*/all_models.json`, `alignment.ali`, `run_modeller.py` | model generation provenance | appendix/provenance | found | raw/script/log-like | detailed model score interpretation not extracted here |

### Rotational orientation scan / MD outputs

| Category | File path | Role | Slide use | Status | Source type | Caveat |
|---|---|---|---|---|---|---|
| orientation comparison table | `02_orientation_metrics/orientation_comparison_table.csv` | parsed RMSD/Rg/H-bond comparison across orientations | Slide 4-5 metric table | found | CSV | mostly 50-100 ns; membrane/contact metrics missing for some orientations |
| orientation summary | `02_orientation_metrics/orientation_metric_summary.md` | narrative summary and ranking | Slide 4-5 rationale | found | documentation | states +10 top-ranked among available parsed metrics, with caveats |
| selected orientation rationale | `03_receptor_pack/metadata/selected_orientation_rationale.md` | packaged +10 selection rationale | Slide 5 current status | found | documentation | older clustering support is 80-100 ns, not requested 160-200 ns |
| orientation figures | `02_orientation_metrics/figures/rmsd_by_orientation.png`, `rg_by_orientation.png`, `interchain_hbond_by_orientation.png`, `orientation_score_summary.png` | ready plots for scan comparison | Slide 4-5 figures | found | figure | KD-membrane distance and positive-patch figures not created |
| +10 processed outputs | `3.rotate_EGFR/10/plot/stable_p10_200ns_*.xvg` | parsed +10 stable series used in orientation comparison | Slide 5 supporting plots | found | raw XVG | filename says 200 ns, parsed x-axis in Task 2 was ~100 ns |
| -10 reference outputs | `3.rotate_EGFR/10/plot/collapsed_m10_100ns_*.xvg` | -10 reference/collapsed comparator series | Slide 5 comparison | found | raw XVG | stored under `3.rotate_EGFR/10/plot`; orientation label must be read from filenames/context |
| final metrics figures | `03_receptor_pack/figures/trajectory_metrics/plus10_200ns_*.png` | final 200 ns metric plots from final pair | Slide 5-6 final +10 validation visuals | found | computed from existing trajectory | output from existing trajectory, not new simulation |

### MD trajectory / topology / setup

| Category | File path | Role | Slide use | Status | Source type | Caveat |
|---|---|---|---|---|---|---|
| GROMACS topology | `03_receptor_pack/topology/topol.top` and `toppar/*.itp` | final topology package | appendix/provenance | found | raw file | exact lipid/ion counts can be extracted from topology if needed |
| TPR | `03_receptor_pack/topology/step5_production.tpr` | supporting topology/run input | appendix/provenance | found | raw file | not directly inspected with GROMACS in this run |
| CHARMM-GUI input configs | `2.charmm-gui/charmm*/input.config.dat` | bilayer setup and forcefield info | Slide 5 setup note | found | raw config | sign mapping from folder labels inferred unless separately documented |
| production MDP | `2.charmm-gui/charmm*/gromacs/step7_production.mdp` | MD production parameters | appendix/provenance | found | raw config | `dt=0.002 ps`, `nsteps=50000000`, output every 50,000 steps, v-rescale, C-rescale, LINCS/h-bonds |
| CHARMM-GUI GROMACS README | `2.charmm-gui/charmm*/gromacs/README` | setup instructions/provenance | appendix | found | documentation | not scientific validation |

### Final membrane-stability / validation metrics

| Category | File path | Role | Slide use | Status | Source type | Caveat |
|---|---|---|---|---|---|---|
| final 200 ns membrane summary | `03_receptor_pack/metadata/membrane_stability_plus10_200ns_mdanalysis.csv` | full-window trajectory-derived metrics | Slide 5-6 validation | found | computed from existing trajectory | `analysis_status=OK`; selection split inferred because GRO lacks chain IDs |
| final last50 summary | `03_receptor_pack/metadata/membrane_stability_plus10_last50ns_150_200ns_mdanalysis.csv` | 150-200 ns final-window metrics | Slide 5-6 stability note | found | computed from existing trajectory | derived from final 200 ns trajectory; no separate gmx_check for last50 file |
| final timeseries | `03_receptor_pack/metadata/membrane_stability_timeseries_plus10_200ns.csv` | full-window metric timeseries | possible plot source | found | computed from existing trajectory | large CSV; use existing PNGs where possible |
| final analysis log | `03_receptor_pack/metadata/membrane_stability_analysis_log.txt` | method, selections, validation pass | provenance/caveats | found | log | positive-patch labels 689/692/713/715 are VAL/LEU/LYS/ILE, not all lysines |
| computed metric manifest | `03_receptor_pack/metadata/computed_metric_manifest.csv` | SHA256 and generated output manifest | provenance | found | documentation/manifest | paths in manifest may reflect original generation machine |
| older partial step5 metrics | `03_receptor_pack/metadata/run_step5_clustering_20260429/*` | archived/provenance-only | do not use in main slide metrics | found | documentation/provenance | contains partial/NOT_AVAILABLE outputs; not final quantitative evidence |

Selected final 200 ns membrane-stability values suitable for cautious use:

| Metric | 0-200 ns mean ± SD | 150-200 ns mean ± SD | Source |
|---|---:|---:|---|
| KD-membrane absolute distance | 25.63 ± 1.23 Å | 25.66 ± 0.86 Å | `membrane_stability_*_mdanalysis.csv` |
| Chain A/B KD-membrane distance difference | 0.49 ± 0.40 Å | 0.40 ± 0.26 Å | same |
| TM-KD tilt angle | 176.39 ± 1.86° | 174.85 ± 1.00° | same |
| KD dimer out-of-plane tilt | 0.70 ± 0.58° | 0.57 ± 0.37° | same |
| Interchain heavy contacts within 4.5 Å | 507.17 ± 43.48 | 505.28 ± 38.10 | same |

### Clustering output

| Category | File path | Role | Slide use | Status | Source type | Caveat |
|---|---|---|---|---|---|---|
| packaged older cluster summary | `03_receptor_pack/metadata/clustering/clustering_summary_80-100ns_cut0.15.txt` | dominant 80-100 ns cluster info | Slide 7 supporting older clustering | found | documentation/log | not requested 160-185/170-200 ns |
| packaged older representative | `03_receptor_pack/receptor/EGFR_plus10_cluster_rep_80-100ns_cut0.15_cluster1.pdb` | selected cluster 1 representative | Slide 7 visual candidate | found | raw file | representative time 84.0 ns according to summary |
| `4.+10_clustering` summaries | `4.+10_clustering/cluster_results/70-85ns_cut0.15/*`, `80-100ns_cut0.15/*` | older clustering outputs | Slide 7 caveat/supporting evidence | found | raw/log/figure | 70-85 and 80-100 ns only |
| apparent 160-185 PDB | `3.rotate_EGFR/10/cluster_results/60-85ns_cut0.16/EGFR_160-185.pdb` | possible structure file by filename | Slide 7: do not use as confirmed interval | ambiguous | raw file/provenance | folder summary says 60-85 ns; audit marks source_provenance |
| apparent 170-200 PDB | `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173/EGFR_170-200.pdb` | possible structure file by filename | Slide 7: do not use as confirmed interval | ambiguous | raw file/provenance | folder summary says 70-100 ns; audit marks source_provenance |
| requested 160-185 ns C1 representative | not found as confirmed clustering result | requested interval representative | Slide 7 missing row | missing / not found | missing | only ambiguous filename found under 60-85 ns result folder |
| requested 170-200 ns C1/C2 representatives | not found as confirmed clustering result | requested interval representatives | Slide 7 missing row | missing / not found | missing | only ambiguous filename found under 70-100 ns result folder |

Confirmed older clustering summary for possible PPT backup:

| Interval | Cutoff | Cluster | Frames | Population | Representative frame/time | Source |
|---|---:|---:|---:|---:|---|---|
| 70-85 ns | 0.15 nm | C1 | 138/151 | 91.4% | frame 749, 74.9 ns | `4.+10_clustering/cluster_results/70-85ns_cut0.15/clustering_summary.txt` |
| 80-100 ns | 0.15 nm | C1 | 150/201 | 74.6% | frame 840, 84.0 ns | `4.+10_clustering/cluster_results/80-100ns_cut0.15/clustering_summary.txt` |
| 80-100 ns | 0.15 nm | C2 | 38/201 | 18.9% | frame 978, 97.8 ns | same |

### Figure / image candidates

| Slide use | File path | Status | Caption draft | Caveat |
|---|---|---|---|---|
| orientation scan comparison | `02_orientation_metrics/figures/rmsd_by_orientation.png` | found | Orientation scan RMSD comparison across -20°, -10°, 0°, +10°, +20°. | parsed available metrics only |
| orientation score | `02_orientation_metrics/figures/orientation_score_summary.png` | found | Available-metric ranking supporting +10° selection. | not a complete biological stability score |
| Rg comparison | `02_orientation_metrics/figures/rg_by_orientation.png` | found | Radius of gyration comparison by orientation. | 50-100 ns parsed window |
| H-bond comparison | `02_orientation_metrics/figures/interchain_hbond_by_orientation.png` | found | Inter-chain H-bond comparison where available. | missing for -20/0/+20 |
| final 200 ns KD-membrane distance | `03_receptor_pack/figures/trajectory_metrics/plus10_200ns_kd_membrane_abs_distance_A.png` | found | Final +10 200 ns KD-to-cytosolic-leaflet distance. | computed from existing trajectory |
| final 200 ns interchain contacts | `03_receptor_pack/figures/trajectory_metrics/plus10_200ns_interchain_heavy_contacts_4p5A.png` | found | Inter-protomer heavy-atom contacts across final +10 200 ns trajectory. | contacts, not H-bonds |
| final 200 ns TM-KD tilt | `03_receptor_pack/figures/trajectory_metrics/plus10_200ns_tm_kd_tilt_angle_deg.png` | found | TM-KD vector angle relative to membrane normal. | selection definitions from analysis log |
| cluster population | `4.+10_clustering/cluster_results/80-100ns_cut0.15/cluster_population.png` | found | 80-100 ns cluster population distribution. | older/supporting clustering only |
| structural coloring scripts | `03_receptor_pack/figures/EGFR_plus10_domain_coloring.pml`, `.cxc` | found | Suggested PyMOL/ChimeraX domain coloring for +10 receptor. | script present; rendered image not generated in this run |

If new structure renders are needed, use the existing scripts rather than inventing images. Suggested coloring: chain A blue/cyan, chain B orange/magenta, TM dark/slate, JM gold, KD lighter chain colors, membrane plane gray transparent, +10 selected model highlighted in green.

### Command / script / log files

| File path | Role | Slide/provenance use | Status | Caveat |
|---|---|---|---|---|
| `scripts/audit_dataset.py` | dataset audit and gmx_check parser | provenance | found | standard-library audit; does not run new MD |
| `03_receptor_pack/scripts/hpc_plus10_mdanalysis.py` | final-pair MDAnalysis rerun script | provenance for final metrics | found | uses existing trajectory only |
| `03_receptor_pack/metadata/membrane_stability_analysis_log.txt` | analysis run log | final metric validation | found | includes selection caveats |
| `4.+10_clustering/clustering.py` | GROMOS clustering script | clustering method provenance | found | older clustering input is step5 production trajectory |
| `1.align/rot_outputs/rot+10/run_modeller.py` | +10 MODELLER run script | construction provenance | found | execution state should not be over-claimed |
| `1.align/rot_outputs/rot+10/alignment.ali` | MODELLER alignment | construction provenance | found | useful for appendix only |

## B. 슬라이드별 내용 초안

### Slide 1. Title

Title: EGFR-MYO1D PPI 분석을 위한 membrane-compatible EGFR dimer receptor model 구축 및 clustering

One-line summary: 후속 MYO1D docking/PPI 분석에 앞서, EGFR inactive TM-JM-KD dimer receptor input을 +10° orientation과 final no-water MD trajectory 기준으로 정리했다.

Scope: receptor modeling → TM-KD orientation scan → 200 ns MD validation → clustering representative 검토. PPI docking, pocket discovery, compound docking 결과는 이번 발표 범위에서 제외한다.

### Slide 2. Why receptor model first?

- MYO1D PPI 분석을 시작하기 전 EGFR 쪽 receptor frame이 고정되어야 docking pose와 interface 해석이 흔들리지 않는다.
- EGFR은 membrane-proximal TM/JM/KD geometry가 중요하므로 monomer 또는 isolated KD만으로는 후속 PPI context를 충분히 반영하기 어렵다.
- 이번 발표는 receptor input 확정 단계까지이며, MYO1D docking consensus, PPI-adjacent pocket, compound docking은 future work로만 둔다.

### Slide 3. Composite model construction

| Component | Source | Role | Current status | Caveat |
|---|---|---|---|---|
| TM dimer | 2M0B | EGFR transmembrane dimer anchor | found | local script uses 634-670 segment |
| JM-A reference | 2M20 | membrane-proximal JM-A orientation | found | exact coordinate-generation path remains `needs confirmation` |
| inactive KD dimer | 3GT8 | inactive symmetric kinase-domain dimer | found | PDB numbering vs model/UniProt-like numbering caveat |
| gap filling | MODELLER | JM-B 684-704 and A-loop 862-876 modeled gaps | found | modeled, not template-derived |

Workflow draft: 2M0B TM dimer alignment → 2M20 JM-A orientation reference → 3GT8 inactive symmetric KD dimer placement → MODELLER gap filling → five TM-KD rotational variants → CHARMM-GUI membrane setup.

Figure suggestion: use `EGFR_plus10_final_model_EGFR_rot10_best.pdb` with domain coloring script from `03_receptor_pack/figures/EGFR_plus10_domain_coloring.pml` or `.cxc`.

### Slide 4. Rotational orientation scan

- Five variants: -20°, -10°, 0°, +10°, +20°.
- Purpose: TM/JM anchor is fixed near membrane, but KD orientation relative to membrane and dimer interface can vary; scan selects a receptor frame that is geometrically and dynamically more suitable.
- Available scan metrics: RMSD, Rg, inter-chain H-bonds where available, clustering availability, visual/membrane compatibility.
- Missing for scan-level quantitative use: orientation-wide KD-to-lipid-plane distance, positive patch contacts, some inter-chain H-bonds.

| Orientation | RMSD mean ± SD Å | Rg mean ± SD Å | H-bond mean ± SD | Cluster rep | Quality flag |
|---|---:|---:|---:|---|---|
| -20° | 4.62 ± 0.50 | 32.51 ± 0.28 | not available | no | moderate RMSD available |
| -10° | 5.00 ± 0.25 | 28.85 ± 0.15 | 6.96 ± 2.16 | yes | moderate RMSD available |
| 0° | 8.50 ± 0.36 | 32.48 ± 0.15 | not available | no | high RMSD available |
| +10° | 2.68 ± 0.39 | 28.31 ± 0.12 | 11.49 ± 3.11 | yes | low RMSD available |
| +20° | 7.54 ± 0.51 | 33.46 ± 0.17 | not available | no | high RMSD available |

Suggested figures: `rmsd_by_orientation.png`, `rg_by_orientation.png`, `interchain_hbond_by_orientation.png`, `orientation_score_summary.png`.

### Slide 5. MD validation result: +10° selected

Main message: +10°는 available parsed orientation metrics에서 가장 favorable하며, final +10 no-water trajectory는 gmx_check log로 200 ns identity가 확인된다.

Evidence bullets:

- Final +10 trajectory log: 49,854 atoms, 2,001 frames, 100 ps timestep, last frame 200,000 ps.
- Final membrane-stability rerun uses `EGFR_plus10_step7_production_nw.gro` + `EGFR_plus10_step7_200ns_nw.xtc` and reports `analysis_status=OK`.
- 150-200 ns window에서 KD-membrane absolute distance는 25.66 ± 0.86 Å, chain A/B distance difference는 0.40 ± 0.26 Å, interchain heavy contacts는 505.28 ± 38.10이다.
- `3.rotate_EGFR/10/full_200ns.xtc`는 final evidence로 사용하지 않는다. excluded trajectory이며 atom-count/timestep red flags가 문서화되어 있다.

Conclusion sentence: 현재 데이터 기준으로는 +10° EGFR dimer를 후속 MYO1D PPI 분석을 위한 receptor input의 우선 후보로 사용하는 것이 가장 합리적이며, 이 주장은 receptor modeling/MD validation 범위 안에서만 말한다.

### Slide 6. Structural interpretation of +10° model

- +10° model은 symmetric inactive TM-JM-KD dimer context를 유지하는 receptor frame으로 정리되어 있다.
- Final 200 ns analysis에서 KD-membrane distance와 interchain heavy contacts가 trajectory-derived metric으로 계산되어, receptor frame을 후속 분석 input으로 검토할 근거가 생겼다.
- 단, final topology의 chain ID가 `.gro`에 없어서 protomer A/B는 residue-order split으로 추정했다. 이 방법은 report caveat로 유지한다.
- Positive-patch label 689/692/713/715는 final topology에서 VAL689, LEU692, LYS713, ILE715이므로 모두 lysine처럼 표현하면 안 된다.

Figure suggestion: side view/top view structural render from `EGFR_plus10_step7_production_nw.gro` or `EGFR_plus10_final_model_EGFR_rot10_best.pdb`; membrane plane 표시; chain/domain coloring script 사용.

### Slide 7. Trajectory clustering

Purpose: 단일 snapshot이 아니라 MD ensemble에서 representative receptor state를 확보하기 위한 단계이다.

Current confirmed and missing items:

| Requested item | Current finding | Status | PPT handling |
|---|---|---|---|
| 160-185 ns C1 representative | confirmed clustering summary not found | missing / not found | main claim으로 쓰지 않음 |
| 170-200 ns C1 representative | confirmed clustering summary not found | missing / not found | main claim으로 쓰지 않음 |
| 170-200 ns C2 representative | confirmed clustering summary not found | missing / not found | main claim으로 쓰지 않음 |
| `EGFR_160-185.pdb` | file exists under `3.rotate_EGFR/10/cluster_results/60-85ns_cut0.16/` | ambiguous | file name만 언급 가능, interval claim 금지 |
| `EGFR_170-200.pdb` | file exists under `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173/` | ambiguous | file name만 언급 가능, interval claim 금지 |
| 80-100 ns cluster 1 representative | cluster 1 = 150/201 frames, 74.6%, representative time 84.0 ns | found, supporting older clustering | backup figure/evidence로 사용 가능 |

Recommended PPT choice: 지금은 `EGFR_plus10_cluster_rep_80-100ns_cut0.15_cluster1.pdb`를 “supporting older clustering representative”로만 보여주고, 교수님께 final 160-200 ns window에서 representative를 새로 확정/계산해도 되는지 확인하는 구성이 안전하다.

### Slide 8. Current confirmed status

| Item | Can we say it now? | Evidence | Wording for professor |
|---|---|---|---|
| EGFR symmetric inactive TM-JM-KD dimer composite model was built | yes, with construction caveats | `1.align/`, final PDBs, template mapping | EGFR TM/JM/KD composite receptor model을 구성했다 |
| +10° selected receptor is current working receptor input | yes, with metric caveats | orientation metrics, selected orientation rationale | available parsed metrics 기준 +10°를 우선 receptor input으로 선택했다 |
| final 200 ns no-water trajectory exists and passes uploaded gmx_check identity check | yes | final gmx_check log | final +10 no-water trajectory는 200 ns / 2,001 frames로 확인된다 |
| final 200 ns membrane-stability metrics were generated from existing trajectory | yes | MDAnalysis output CSV/log | 기존 trajectory에서 final-pair metric을 계산했고 OK rows로 확인된다 |
| 160-185/170-200 ns clustering representatives are secured | no | no confirmed summary found | 해당 window representative는 needs confirmation으로 남긴다 |
| MYO1D docking/PPI/pocket/compound results are final | no | not in current scope | 이번 발표에서는 future work로만 언급한다 |

### Slide 9. Next steps and feedback request

- 교수님께 확인할 핵심 질문: final +10 trajectory의 150-200 ns 또는 160-200 ns window에서 clustering representative를 새로 확정해서 후속 docking input으로 사용해도 되는가?
- 다음 단계 후보:
  1. MYO1D TH1 beta-meander docking
  2. EGFR-side consensus PPI patch extraction
  3. PPI-adjacent non-ATP pocket analysis
  4. compound docking
- 위 단계들은 이번 발표에서 결과처럼 말하지 않고, receptor input 확정 이후의 future work로만 언급한다.

## C. 각 슬라이드별 필요한 figure 목록

| Slide | Figure title | File path | Caption draft | If missing |
|---|---|---|---|---|
| 1 | Selected EGFR +10 receptor cover visual | `03_receptor_pack/receptor/EGFR_plus10_final_model_EGFR_rot10_best.pdb` + `.pml/.cxc` scripts | Membrane-compatible EGFR inactive TM-JM-KD dimer receptor input for MYO1D analysis. | Render in PyMOL/ChimeraX with chain/domain coloring |
| 2 | Scope boundary schematic | none required | Receptor modeling and MD validation are reported before PPI/pocket/compound stages. | Native PPT diagram recommended |
| 3 | Composite model construction workflow | `1.align/2m0b.cif`, `2m20.cif`, `3gt8.cif`, final +10 PDB | TM, JM-A, inactive KD, and modeled gap components used to build the receptor. | Use simple template table and one structural render |
| 4 | Orientation scan metric summary | `02_orientation_metrics/figures/orientation_score_summary.png` | Available parsed metrics rank +10° highest among tested TM-KD rotations. | Use table if figure is not presentation-ready |
| 4 | RMSD/Rg/H-bond by orientation | `rmsd_by_orientation.png`, `rg_by_orientation.png`, `interchain_hbond_by_orientation.png` | +10° shows lower RMSD and higher available H-bond count than -10° reference. | Keep caveat for missing H-bond data |
| 5 | Final trajectory validation / metrics | `03_receptor_pack/figures/trajectory_metrics/plus10_200ns_kd_membrane_abs_distance_A.png` | Final +10 trajectory-derived KD-membrane distance over 200 ns. | Plot from timeseries CSV if figure needs restyling |
| 5 | Final interchain contact metric | `plus10_200ns_interchain_heavy_contacts_4p5A.png` | Inter-protomer heavy contacts remain measurable across final 200 ns trajectory. | Explain as heavy contacts, not H-bonds |
| 6 | +10 model side/top view | final model PDB/GRO + coloring scripts | Structural interpretation of +10 receptor geometry in membrane context. | Generate with PyMOL/ChimeraX; show membrane plane and normal |
| 7 | Older clustering population | `4.+10_clustering/cluster_results/80-100ns_cut0.15/cluster_population.png` | Supporting older 80-100 ns clustering; dominant C1 is 74.6%. | Mark as supporting older clustering, not requested 160-200 |
| 8 | Confirmed vs needs confirmation table | no raster required | Current confirmed receptor-modeling status and unresolved items. | Native PPT table recommended |
| 9 | Next-step flow | no raster required | Receptor input confirmation precedes MYO1D docking, PPI patch extraction, pocket analysis, and compound docking. | Native PPT flow diagram recommended |

Suggested structure visualization command drafts:

```text
PyMOL: load 03_receptor_pack/receptor/EGFR_plus10_final_model_EGFR_rot10_best.pdb
PyMOL: run 03_receptor_pack/figures/EGFR_plus10_domain_coloring.pml
ChimeraX: open 03_receptor_pack/receptor/EGFR_plus10_final_model_EGFR_rot10_best.pdb
ChimeraX: open 03_receptor_pack/figures/EGFR_plus10_domain_coloring.cxc
```

## D. 발표 멘트

### Slide 1 speaker note

오늘은 EGFR-MYO1D 프로젝트에서 MYO1D docking이나 pocket 분석 결과가 아니라, 그 전에 필요한 EGFR receptor input을 어떻게 정리했는지 보고드리겠습니다. 현재 범위는 EGFR inactive TM-JM-KD dimer model 구축, TM-KD orientation scan, +10° receptor 선택, final 200 ns no-water trajectory validation, 그리고 clustering representative 검토까지입니다.

### Slide 2 speaker note

후속 MYO1D PPI docking을 해석하려면 먼저 EGFR 쪽 receptor frame이 고정되어야 합니다. 특히 EGFR은 membrane-proximal TM/JM orientation과 kinase domain dimer context가 중요하기 때문에 isolated KD나 monomer만으로는 후속 interface 해석이 흔들릴 수 있습니다. 그래서 이번 발표에서는 receptor model과 MD validation까지만 정리하고, pocket이나 compound 관련 내용은 결과로 말하지 않겠습니다.

### Slide 3 speaker note

Composite model은 2M0B의 TM dimer, 2M20의 JM-A reference, 3GT8의 inactive symmetric kinase-domain dimer를 조합하고, JM-B와 일부 loop/gap은 MODELLER로 보완하는 방식으로 구성되어 있습니다. 다만 2M20 기반 JM-A 좌표 생성 경로는 완전히 독립적인 run log로 확인된 것은 아니어서 caveat로 남겨두겠습니다. 최종 모델 residue numbering은 chain A 634-1014, chain B 634-1010 범위로 확인됩니다.

### Slide 4 speaker note

TM-KD orientation은 -20°, -10°, 0°, +10°, +20° 다섯 가지로 scan했습니다. 사용 가능한 parsed metric 기준으로 +10°가 RMSD가 낮고, Rg가 안정적이며, inter-chain H-bond도 상대적으로 높게 나타났습니다. 단, 모든 orientation에서 membrane-distance나 contact metric이 완비된 것은 아니므로, 이것은 available metric 기반의 선택 근거로 표현하는 것이 맞습니다.

### Slide 5 speaker note

최종 +10° pair는 `EGFR_plus10_step7_production_nw.gro`와 `EGFR_plus10_step7_200ns_nw.xtc`를 기준으로 정리했습니다. uploaded gmx_check log에서는 49,854 atoms, 2,001 frames, 100 ps timestep, 200 ns 범위가 확인됩니다. 이 trajectory에서 계산한 final membrane-stability output도 OK 상태이며, 따라서 현재 단계에서는 +10° receptor를 후속 분석 input 후보로 쓰는 근거를 확보했다고 말할 수 있습니다.

### Slide 6 speaker note

구조적으로는 +10° model이 symmetric inactive dimer context와 membrane-proximal frame을 유지하는 receptor input으로 정리되어 있습니다. final 200 ns trajectory에서는 KD-membrane distance, TM-KD tilt, interchain heavy contact 같은 지표를 계산해 receptor geometry를 점검했습니다. 다만 chain ID가 없는 GRO에서 protomer를 residue order로 나눈 점과, 일부 positive-patch label residue가 실제로 lysine이 아니라는 점은 해석 caveat로 유지하겠습니다.

### Slide 7 speaker note

Clustering은 단일 snapshot이 아니라 MD에서 대표 receptor state를 확보하기 위한 단계입니다. 현재 확인된 cluster output 중에서는 80-100 ns 구간의 dominant cluster representative가 있고, cluster 1이 150/201 frames, 74.6%로 보고됩니다. 반면 이번에 확인 요청된 160-185 ns와 170-200 ns representative는 파일명만으로는 확정할 수 없고, 주변 summary가 60-85 ns 또는 70-100 ns로 되어 있어 needs confirmation으로 두겠습니다.

### Slide 8 speaker note

현재 확정해서 말씀드릴 수 있는 것은 EGFR TM-JM-KD composite model을 만들었고, available metric 기준으로 +10° receptor를 우선 선택했으며, final no-water 200 ns trajectory identity를 log로 확인했다는 점입니다. 아직 확정했다고 말하면 안 되는 것은 MYO1D docking 결과, PPI pocket, compound docking 결과, 그리고 160-200 ns clustering representative의 확정입니다. 이 구분을 명확히 해서 다음 단계로 넘어가고자 합니다.

### Slide 9 speaker note

교수님께 확인드리고 싶은 핵심은 final +10 trajectory의 후반부, 예를 들면 150-200 ns 또는 160-200 ns window에서 clustering representative를 새로 확정해서 후속 docking input으로 사용해도 되는지입니다. 승인이 되면 다음 단계는 MYO1D TH1 beta-meander docking, EGFR-side PPI patch extraction, PPI-adjacent pocket analysis, compound docking 순서로 진행하겠습니다. 다만 오늘 발표에서는 이 단계들을 future work로만 언급하겠습니다.

## E. 최종 summary

### 1. 현재 확정했다고 말할 수 있는 claim

- EGFR-MYO1D 후속 분석을 위한 EGFR TM-JM-KD dimer receptor model을 구성했다.
- Available parsed orientation metrics 기준으로 +10° model을 우선 receptor input으로 선택했다.
- Final +10 no-water trajectory `EGFR_plus10_step7_200ns_nw.xtc`는 uploaded gmx_check log에서 49,854 atoms, 2,001 frames, 100 ps timestep, 0-200 ns로 확인된다.
- Existing final trajectory에서 membrane-stability metrics를 계산한 output이 있고, OK rows에는 `NOT_AVAILABLE` 값이 없다.

### 2. 아직 말하면 안 되는 claim

- MYO1D docking consensus가 확정됐다.
- PPI pocket 또는 PPI-adjacent non-ATP pocket을 확정했다.
- compound 후보 또는 drug candidate를 제시했다.
- 요청된 160-185 ns / 170-200 ns clustering representative가 확정됐다.
- Positive-patch 689/692/713/715를 모두 lysine 기반 patch로 해석했다.

### 3. 교수님께 질문해야 할 point

- Final +10 trajectory의 150-200 ns 또는 160-200 ns window에서 representative clustering을 새로 확정해도 되는가?
- 후속 MYO1D docking input으로 pre-MD final model, final 200 ns coordinate/GRO-derived structure, 또는 confirmed cluster representative 중 무엇을 우선 사용할 것인가?
- Positive-patch residue labeling은 final topology residue identity 기준으로 재정의할 것인가?

## Commands/checks run in this pass

- Created reference map and inventories under `reports/`.
- Read top-level docs: `README.md`, `DATASET_CARD.md`, `MANIFEST.tsv`, `VALIDATION.md`, `EXCLUDE_FROM_DATASET.md`, `ADD_TO_GIT.md`.
- Checked final gmx_check log for atoms/frame/time identity.
- Read final membrane-stability summary CSVs and analysis log.
- Read construction/template mapping and numbering provenance.
- Read orientation comparison/summary and selected-orientation rationale.
- Read clustering summaries for 70-85 ns and 80-100 ns outputs, plus ambiguous `EGFR_160-185.pdb` / `EGFR_170-200.pdb` headers.

## Output paths

- Reference map: `reports/EGFR_modeling_clustering_reference_map.md`
- PPT data package report: `reports/EGFR_modeling_clustering_PPT_data_package.md`
