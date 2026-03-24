# Phase 1 참고 노트

> Phase 1 PPI-first Interface Mapping 관련 기술 노트를 통합한 문서.
> 원본 파일: `phase1_pyrosetta_execution_note.md`, `phase1_sampling_rationale.md`,
> `phase1_ppi_handoff_note.md`, `phase1_lightdock_validation_note.md`,
> `phase1_output_chain_note.md`, `phase1_orientation_filter_note.md`,
> `orientation_filter_design.md` (모두 `docs/archive/`로 이동됨)

---

## 1. 실행 (Execution)

### 1.1 파이프라인 구성 요소

| Component | File | Role |
|-----------|------|------|
| Orchestrator | `egfr_pipeline/pyrosetta_docking/pipeline_manager.py` | 7-step pipeline execution |
| Docking worker | `egfr_pipeline/pyrosetta_docking/movers.py` | Relax, Global Docking, Refinement |
| Scoring worker | `egfr_pipeline/pyrosetta_docking/scoring.py` | InterfaceAnalyzer metrics |
| Utilities | `egfr_pipeline/pyrosetta_docking/pyrosetta_init.py` | PyRosetta init, Pose I/O |

### 1.2 7-Step Pipeline

1. **Relax** -- FastRelax (ref2015), cached in `relaxed_cache/`
2. **Global Docking** -- RigidBodyPerturbMover(360deg, 100A) -> SlideIntoContact -> DockMCMProtocol
3. **Fast Scoring & Filtering** -- v2.0 2-pass or v1.0 single-pass
4. **Full Scoring** -- Complete InterfaceAnalyzer metrics + L_RMSD
5. **Clustering** -- CoM pre-filter + L_RMSD greedy clustering
6. **Selection & Save** -- Round-robin diversity + L_RMSD deduplication
7. **Visualization & Report** -- PyMOL scripts + validation report

### 1.3 Full-Kinase-Domain 호환성

기존 파이프라인은 full-kinase-domain 입력을 코드 수정 없이 처리할 수 있다:

- 임의의 2-chain PDB 입력 (chain A = receptor, chain B = partner)
- FoldTree setup: generic `"A_B"` partner definition
- Scoring, filtering, clustering 모두 잔기 수 독립적
- Excluded residues는 config에서 파싱 (하드코딩 아님)
- Auto-threshold clustering은 chain size에 따라 자동 조정

**Phase 1을 위한 pipeline_manager.py, movers.py, scoring.py, pyrosetta_init.py 수정 불필요.**

### 1.4 Config 변경사항 (Legacy -> Phase 1)

| Setting | Legacy (dimer) | Phase 1 (dimer, +1000 offset) | Why |
|---------|---------------|-------------------|-----|
| input_pdb | `EGFR_dimer_*.pdb` | `docking_*_ext_beta_meander.pdb` | Dimer receptor (+1000 offset) |
| excluded_residues_A | Includes `1713-1720,...` | Includes `1713-1720,...` | 양쪽 monomer membrane-proximal |
| total_global_models | 50,000 (single seed) | 20,000 x 5 seeds | Multi-seed strategy |
| random_seed | auto | Deterministic per seed | Reproducibility |
| n_cpus | 32 | 16 | Shared HPC safety |

변경되지 않은 설정: filter thresholds (v2.0 2-pass), MiniRefinement, clustering parameters (auto-adaptive), refinement protocol, output format, ExperimentalData section.

### 1.5 Run Metadata Schema

모든 Phase 1 실행은 `pyrosetta_run_metadata.json`을 생성한다:

```json
{
    "receptor_id": "3GT8_raw",
    "partner_id": "extended_beta_meander_955_1006",
    "construct_type": "dimer_offset",
    "config_file": "config/phase1/phase1_prod_3GT8_raw_seed0.ini",
    "input_pdb": "input/PPI/phase1/docking_3GT8_raw_ext_beta_meander.pdb",
    "total_global_models": 20000,
    "n_cpus": 16,
    "random_seed": "42",
    "seed_index": 0,
    "is_production": true,
    "filter_version": "v2.0",
    "phase": "Phase 1: PPI-first Interface Mapping",
    "task_group": "TG 1.1: PyRosetta Global Docking Standardization"
}
```

**현재 구현 참고**: 위 JSON은 Phase 1 v2 목표 상태를 반영한다. 현재 구현은 `partner_id`로 `MYO1D_beta_meander` 또는 `MYO1D_TH1` 같은 값을 사용하며, `partner_construct`, `receptor_chain_ids`, `partner_chain_ids`, `numbering_system`, `n_cpus_requested`/`n_cpus_used`, `run_label`, `run_status`, `input_validation_status` 등 추가 필드를 포함한다.

### 1.6 Key Metadata Fields for Downstream Traceability

| Field | Purpose | Used by |
|-------|---------|---------|
| `receptor_id` | Cross-state comparison | TG 1.5 |
| `partner_id` | Partner identification | All downstream |
| `construct_type` | Distinguish from legacy | TG 1.7 pilot comparison |
| `seed_index` | Multi-seed consolidation | Score standardization |
| `filter_version` | Audit trail | Reproducibility |

### 1.7 실행 명령어

```bash
# Config 생성
python -m egfr_pipeline.phase1.generate_configs

# Codex 워크스페이스 dry-run
python -m egfr_pipeline.phase1.launch_docking --test --dry-run

# 서버 테스트 (단일 state, 1K models, ~2-4시간)
conda activate pyrosetta
python -m egfr_pipeline.phase1.launch_docking --test --state 3GT8_raw

# 서버 프로덕션 (전체 state, 전체 seed)
python -m egfr_pipeline.phase1.launch_docking --production

# 개별 seed 병렬 실행
for state in 3GT8_raw EGFR_160-185 EGFR_170-200; do
    for seed in 0 1 2 3 4; do
        python -m egfr_pipeline.phase1.launch_docking \
            --config config/phase1/phase1_prod_${state}_seed${seed}.ini &
    done
done

# 점수 표준화 (도킹 완료 후)
python -m egfr_pipeline.phase1.standardize_scores
```

---

## 2. 샘플링 (Sampling Rationale)

### 2.1 Production Budget

- 20,000 models per seed
- 5 seeds per receptor state
- 100,000 total models per receptor state
- 300,000 total models across the 3 receptor states

### 2.2 Multi-Seed 설계 근거

`20k x 5 seeds`를 단일 `100k` 실행보다 선호하는 이유:

1. **충분한 global-search coverage**: Rosetta global docking은 수렴을 위해 10,000-100,000 trajectories가 필요하다. 20,000은 허용 범위 내이지만 low-to-middle에 해당한다.
2. **명시적 수렴 및 robustness 평가**: 5개 독립 stochastic restart로 인터페이스 패치의 재현성을 평가할 수 있다. 단일 seed에서만 발견되는 패치 vs 여러 seed에서 반복 발견되는 패치를 구분 가능.

### 2.3 해석 정책

| Budget | 용도 |
|--------|------|
| 20k x 1 seed | Smoke test, pipeline validation, first-pass interface discovery |
| 20k x 3 seeds (60k/state) | Intermediate production checkpoint |
| 20k x 5 seeds (100k/state) | Final production, manuscript-supported interpretation |

### 2.4 시스템 크기

현재 Phase 1 setup: receptor ~309 residues + extended beta-meander partner ~52 residues = ~361 residues. Rosetta tutorial에서 global docking은 비교적 작은 복합체에 적합하다고 명시하며, 이 시스템은 해당 범위 내.

### 2.5 Conformational Flexibility 한계

Rosetta global docking은 fixed-backbone을 가정한다. 20k -> 100k로 증가해도 backbone mismatch나 partner flexibility 문제는 해결되지 않는다. Seed-to-seed agreement가 100k에서도 낮다면 limiting factor는 conformational modeling.

### 2.6 Manuscript-Ready Wording

**Methods:**

> For Phase 1 protein-protein docking, we performed blind global docking with PyRosetta/RosettaDock using 20,000 trajectories per random seed and five independent seeds for each receptor state, yielding 100,000 trajectories per state. This sampling level was chosen because Rosetta global docking typically requires on the order of 10,000-100,000 decoys for convergence, with 100,000 decoys representing a commonly recommended upper production target for global searches. We distributed this budget across independent seeds rather than a single run in order to evaluate reproducibility of recovered interface patches across stochastic restarts.

**Results:**

> Interface patches that were recovered across multiple independent seeds were interpreted as more robust than patches observed in only a single seed, even when the total number of sampled decoys was comparable. This multi-seed design allowed us to distinguish reproducible docking features from seed-specific stochastic outcomes.

**Limitations:**

> Increasing the total number of rigid-body docking trajectories improves search coverage but does not fully address limitations associated with backbone mismatch or partner flexibility. Accordingly, lack of convergence at a given site was interpreted cautiously, particularly for a flexible beta-meander-like partner, and was not treated as definitive evidence of absence of binding.

### 2.7 References

1. Rosetta protein-protein docking tutorial: https://docs.rosettacommons.org/demos/latest/tutorials/Protein-Protein-Docking/Protein-Protein-Docking
2. RosettaDock application documentation: https://docs.rosettacommons.org/docs/latest/application_documentation/docking/docking-protocol
3. Marze NA, Roy Burman SS, Sheffler W, Gray JJ. Efficient flexible backbone protein-protein docking for challenging targets. Bioinformatics. 2018. https://pmc.ncbi.nlm.nih.gov/articles/PMC6184633/
4. Alam N, Goldstein O, Xia B, Porter KA, Kozakov D, Schueler-Furman O. High-resolution global peptide-protein docking using fragments-based PIPER-FlexPepDock. PLoS Comput Biol. 2017. https://pmc.ncbi.nlm.nih.gov/articles/PMC5760072/

---

## 3. 인터페이스 추출 및 점수 표준화 (Interface Extraction & Score Standardization)

### 3.1 pyrosetta_decoy_scores.csv Schema

| Column | Type | Source | Notes |
|--------|------|--------|-------|
| decoy_id | str | File_PDB | PDB filename |
| receptor_id | str | metadata | e.g., "3GT8_raw" |
| partner_id | str | metadata | "extended_beta_meander_955_1006" |
| construct_type | str | metadata | "dimer_offset" |
| seed_index | int | metadata | 0-4 for production |
| run_type | str | metadata | "test" or "production" |
| Rank | int | final_ranking.csv | 1-20 within seed |
| cluster_id | str | Parent -> "C01" | Cluster identifier |
| total_score | float | Total_Score | Full Rosetta score (REU) |
| I_sc | float | dG_separated | Interface score proxy |
| dG_separated | float | InterfaceAnalyzer | Binding energy (REU) |
| dSASA | float | InterfaceAnalyzer | Buried surface (A^2) |
| dG_density | float | derived | dG/dSASA x 100 |
| sc_value | float | InterfaceAnalyzer | Shape complementarity (0-1) |
| packstat | float | InterfaceAnalyzer | Packing density (0-1) |
| delta_unsatHbonds | int | InterfaceAnalyzer | Unsatisfied H-bonds |
| nres_int | int | InterfaceAnalyzer | Interface residue count |
| hbonds_int | int | InterfaceAnalyzer | Interface H-bond count |
| L_RMSD | float | CalphaSuperimpose | vs relaxed reference (A) |
| L_RMSD_best | float | CalphaSuperimpose | vs best-dG model (A) |
| Binding_Residues_A | str | ContactAnalysis | Receptor interface residues |
| Binding_Residues_B | str | ContactAnalysis | Partner interface residues |
| key_contact_ratio | float | ContactAnalysis | Key residue contact fraction |
| source_file | str | -- | Provenance tracking |

**현재 구현 참고**: 현재 export에는 위 목표 스키마 외에 `receptor_construct`, `partner_construct`, `receptor_chain_ids`, `partner_chain_ids`, `center_x/y/z` 등의 실용적 필드가 추가로 포함되어 있다.

### 3.2 I_sc Note

Task document는 `I_sc`를 primary ranking metric으로 지정한다. 현재 파이프라인에서는 `InterfaceAnalyzerMover`의 `dG_separated`가 이 역할을 한다. True `I_sc` (`InterfaceScoreCalculator`)와 `dG_separated`는 rigid-body docking에서 높은 상관관계를 보인다. 필요시 `scoring.py`에 true `I_sc` 추출을 추가할 수 있다.

---

## 4. Orientation 필터 (Orientation-Aware Filtering)

### 4.1 문제 정의

MYO1D beta-meander는 flat, thin beta-sheet ribbon (5 consecutive beta-strands, sheets 8-12)이다. Global blind docking에서 두 가지 방향으로 착지할 수 있다:

- **Correct**: sheets 8/9의 active face가 receptor와 접촉
- **Flipped**: back face가 receptor와 접촉

Contact-count threshold만으로는 이 구분이 불가능하다. Beta-sheet에서 side chain은 평면 위아래로 교대 배치되며, sheet 두께(~5-7 A)가 접촉 cutoff(8-10 A)보다 얇아 양 면의 잔기가 모두 "접촉"으로 판정될 수 있다. Backbone hydrogen bond도 양 면에서 형성되어 비특이적이지만 에너지적으로 유리한 접촉을 만든다.

### 4.2 Algorithm: Dual-Vector Orientation Test

**Step 1: Active-face CA 좌표 수집**

Sheet 8 (VAL961, VAL962, ASN963, VAL964) + Sheet 9 (VAL968, GLN969, CYS970, SER971, LEU972)의 CA 좌표 추출. 총 9 residues가 active face를 정의한다.

**Step 2: PCA로 sheet-plane normal 계산**

9개 CA 위치에 PCA 적용. 2개의 큰 principal component가 sheet plane을 span하고, 가장 작은 component가 plane normal이 된다. PCA는 단순 cross product보다 robust하며, sheet이 약간 curved/twisted된 경우에도 작동한다.

**Step 3: Multi-probe CA->CB consensus로 normal 방향 결정**

PCA normal은 부호가 모호하다. 이를 해결하기 위해:
- Active-face residues의 CA->CB vector 사용 (default: VAL962, VAL964, SER971 -- sheets 8, 9 모두 커버)
- 각 probe에 대해 dot(normal, CA->CB) 부호 확인
- Majority vote로 normal 방향 결정 (다수가 flip 필요 표시 시 flip)
- 단일 probe보다 rotamer-dependent noise에 강건
- Valid CA/CB가 1개뿐이면 single-probe(VAL962)로 fallback

**Step 4: Local receptor-direction vector 계산**

전체 receptor centroid 대신 localized centroid 사용:
- Active-face CA에서 10 A 이내의 receptor CA atoms 수집
- 이들의 centroid를 receptor contact centroid로 사용
- Sheet 8/9 centroid -> receptor contact centroid 방향 벡터 생성

이 접근은 실제 접촉 기하학을 반영하며, full kinase domain의 먼 global centroid에 영향받지 않는다.

**Step 5: Dot product 분류**

```
orientation_score = dot(active_face_normal, receptor_direction)

if |score| < 0.15:     -> ambiguous (edge-on contact)
elif score > 0:         -> pass (active face toward receptor)
else:                   -> fail (active face away from receptor)
```

### 4.3 Sheet Residue 정의

Ko et al. 실험 데이터 기반:

| Sheet | Residues | Role | Experimental Evidence |
|-------|----------|------|----------------------|
| 8 | 961, 962, 963, 964 | Active face (primary) | Ala substitution -> function abolished |
| 9 | 968, 969, 970, 971, 972 | Active face (primary) | Ala substitution -> function abolished |
| 10 | ~977-980 | Neutral | Ala substitution -> WT-level function |
| 11 | ~985-988 | Neutral | Ala substitution -> WT-level function |
| 12 | ~993-997 | Structural support (role TBD) | Ala substitution -> function abolished |

Sheet 10-12 boundaries는 approximate. Sheet 12는 기능적으로 essential이지만 direct contact인지 structural support인지 미확정. 현재 working assumption: structural support. Active-face 정의에는 sheets 8, 9만 사용하고 sheet 12는 모니터링만 수행.

### 4.4 구현

파일: `egfr_pipeline/phase1/orientation_filter.py`

Core functions:
- `compute_orientation_score()` / `compute_orientation_score_from_pdb()` (PDB-based, PyRosetta 불필요)
- `compute_sheet_plane_normal()`
- `orient_normal_to_active_face()`
- `compute_receptor_interface_centroid()`
- `process_state_orientation()`
- `merge_orientation_into_models()`
- `validate_pilot_structures()`

PDB-based 구현(`compute_orientation_score_from_pdb`)은 LightDock 포즈에도 적용 가능하며, PyRosetta 없이 순수 텍스트 파싱 + numpy SVD로 동일 알고리즘을 실행한다.

### 4.5 출력 스키마

`orientation_filter_log.csv`:

| Column | Type | Description |
|--------|------|-------------|
| filename | str | PDB filename |
| orientation_score | float | Dot product (-1 to +1) |
| orientation_class | str | pass / fail / ambiguous / error codes |
| n_active_face_ca | int | Number of active-face CA found |
| n_back_face_ca | int | Number of back-face CA found |
| sheet_centroid | str | [x, y, z] of sheet 8/9 centroid |
| normal_vector | str | [x, y, z] of active-face normal |
| receptor_direction | str | [x, y, z] toward receptor contact centroid |
| error | str | Error message if any |

### 4.6 Pipeline 내 위치

```
TG 1.1 PyRosetta docking
  -> TG 1.2 Interface extraction (all models)
  -> TG 1.2A Orientation filter
  -> TG 1.3 Consensus (pass models only)
  -> TG 1.5 Multi-state comparison
```

Downstream 사용:
- `orientation_class == "pass"`: consensus building 진입
- `orientation_class == "ambiguous"`: manual review 대상
- `orientation_class == "fail"`: consensus 제외, raw data 보존

### 4.7 실행 명령어

```bash
conda activate pyrosetta
python -m egfr_pipeline.phase1.orientation_filter --merge
python -m egfr_pipeline.phase1.orientation_filter --state 3GT8_raw --merge

# Pilot validation
python -m egfr_pipeline.phase1.orientation_filter \
  --pilot_dir /path/to/pilot/final_result/ \
  --pilot_out output/workflow_a/phase2_ppi_docking/orientation_filter_pilot_validation.csv
```

### 4.8 검증 전략

1. **Retroactive pilot validation**: 기존 5개 유효 pilot 구조 (C02_M01, C02_M03, C04_M01, C04_M02, C07_M03)에 적용. 수동으로 "biologically reasonable"로 판정된 구조는 PASS 예상. FAIL 시 threshold가 너무 strict하거나 실제 face-flip인 것.
2. **Synthetic face-flip test**: 알려진 PASS 구조의 beta-meander를 180도 회전하여 FAIL 확인.
3. **Edge case documentation**: AMBIGUOUS 구조를 PyMOL에서 시각적 검사, ambiguous band threshold 보정에 활용.
4. **Statistical validation**: pass/fail/ambiguous 분포 보고, orientation class vs energy metrics (dG, I_sc) cross-tabulation, PASS 모델의 평균 인터페이스 품질이 FAIL보다 높은지 검증.

### 4.9 한계

- Multi-probe consensus는 3개 probe 중 2개 이상의 valid CA/CB가 필요 (core beta-sheet residues이므로 missing atoms 가능성 낮음)
- Rigid-body assumption: beta-sheet geometry가 도킹 중 유지된다고 가정 (standard RosettaDock에서는 타당)
- Ambiguous band threshold (0.15)는 초기값이며 pilot validation 후 보정 필요
- Sheet 12 역할 미확정: direct contact face로 밝혀지면 active-face definition 확장 필요

---

## 5. LightDock 검증 (LightDock Validation)

### 5.1 역할

LightDock는 secondary evidence only. Phase 1 receptor-side patch review를 지원하지만 PyRosetta를 대체하지 않는다.

Evidence hierarchy:
- **Primary**: PyRosetta interface consensus + cross-state robustness
- **Secondary**: LightDock interface support + cross-method convergence
- **Legacy optional**: AlphaFold-Multimer parser (re-enable 시에만)

### 5.2 실행 워크플로우

```
--setup       Generate run scripts and metadata (all 3 states)
--run         Execute generated bash scripts (requires LightDock on PATH)
--extract     Parse LightDock output PDBs + apply orientation filter
--convergence Cross-method comparison (PyRosetta vs LightDock)
--note        Generate dynamic results summary
--all         setup + extract + convergence + note (no execution)
--check       Verify LightDock availability on PATH
```

PBS scripts:
- `config/run_lightdock.pbs` -- production (400 swarms, 200 glowworms, 100 steps)
- `config/run_lightdock_test.pbs` -- test (50 swarms, 50 glowworms, 50 steps)

Generated `run_lightdock_<state>.sh` 수행 단계:
1. Pre-flight availability check (5 LightDock executables)
2. PDB splitting (chain A/B, preserving TER records)
3. `lightdock3_setup.py` (swarm placement)
4. `lightdock3.py` (optimization)
5. `lgd_generate_conformations.py` (per-swarm pose generation)
6. `lgd_rank.py` (global ranking -> `rank_by_scoring.list`)
7. `lgd_cluster_bsas.py` (clustering)
8. Completion marker (`.lightdock_complete`)

### 5.3 Orientation Filter 적용

LightDock 포즈에도 PDB-based orientation filter 적용:
- `compute_orientation_score_from_pdb()`: pure text parsing + numpy SVD
- PyRosetta 없이 동일 알고리즘 실행
- `orientation_validation_status`: `pass`, `fail`, `ambiguous`, `insufficient_data`
- Convergence analysis에서 `pass` 또는 `not_available` 모델만 residue frequency에 반영

### 5.4 출력 파일

Per receptor state (`output/workflow_a/phase2_ppi_docking/<state>/lightdock/`):
- `lightdock_run_metadata.json`
- `run_lightdock_<state>.sh`
- `lightdock_interface_support_table.csv`
- `lightdock_model_summary.csv`
- `.lightdock_complete` (marker)

Cross-method output:
- `output/workflow_a/phase2_ppi_docking/<state>/cross_method_convergence.csv`
- `output/workflow_a/phase2_ppi_docking/<state>/cross_method_convergence_summary.json`

Results summary:
- `output/workflow_a/phase2_ppi_docking/phase1_lightdock_validation_results.md`

### 5.5 CSV Baselines

**`lightdock_interface_support_table.csv`**: `model_id`, `receptor_id`, `construct_type`, `orientation_validation_status`, `orientation_score`, `swarm_id`, `pose_rank`, `scoring_value`, `chain`, `residue_id`, `residue_num`, `residue_name`, `lobe_label`, `source`

**`lightdock_model_summary.csv`**: `model_id`, `receptor_id`, `construct_type`, `orientation_validation_status`, `orientation_score`, `swarm_id`, `pose_rank`, `scoring_value`, `n_receptor_interface_residues`, `n_partner_interface_residues`, `n_nlobe_interface_residues`, `n_clobe_interface_residues`, `receptor_interface_residues`, `partner_interface_residues`, `source`

**`cross_method_convergence.csv`**: `receptor_id`, `construct_type`, `orientation_validation_status`, `chain`, `residue_id`, `residue_num`, `residue_name`, `lobe_label`, `in_pyrosetta`, `in_lightdock`, `pyrosetta_max_occupancy`, `lightdock_frequency`, `convergence_class`, `method_agreement`

**`cross_method_convergence_summary.json`**: `receptor_id`, `n_convergent`, `n_pyrosetta_only`, `n_lightdock_only`, `n_total`, `jaccard_overlap`, `jaccard_nlobe`, `jaccard_clobe`

### 5.6 Score Comparison Note

PyRosetta는 REU (dG_separated), LightDock는 DFIRE2 scoring (fastdfire) 사용. 직접 비교 불가. Cross-method validation은 residue-level Jaccard overlap으로 수행 (unit-agnostic).

### 5.7 한계

- LightDock chain reassignment: output PDB에서 잔기 번호가 재할당될 경우 orientation filter가 active-face residues를 찾지 못할 수 있음. `insufficient_data` classification으로 fallback.
- `cross_method_convergence.csv`는 비교 레이어이며 PyRosetta patch definition을 대체하지 않는다.
- `lightdock_only` 잔기는 추가 증거 없이 primary Phase 2 patch input으로 사용 불가.

---

## 6. 출력 구조 (Output Chain)

### 6.1 Phase 1 전체 흐름

1. PyRosetta input validation + run metadata
2. PyRosetta docking + decoy score export
3. Interface residue extraction
4. Orientation filtering
5. Cluster consensus
6. Cross-state comparison
7. LightDock secondary validation
8. Phase 1 review report + Phase 2 handoff

### 6.2 디렉토리 구조

```
output/workflow_a/phase2_ppi_docking/
  <state>/                              # e.g., 3GT8_raw, EGFR_160-185, EGFR_170-200
    test_seed0/  |  prod_seed0/ ... prod_seed4/
      pyrosetta_run_metadata.json
      docking_<input_stem>/
        final_ranking.csv
        cluster_results/
        final_result/
        ...
    pyrosetta_decoy_scores.csv          # consolidated across seeds
```

**현재 구현 참고**: 현재 naming rule은 `<input_stem>__<receptor_id>__<partner_construct>__<construct_type>__<run_label>`으로, test/prod 실행이 동일 input PDB를 공유할 때 덮어쓰기를 방지한다.

### 6.3 Primary Evidence 출력 (`output/workflow_a/phase2_ppi_docking/<state>/`)

- `phase1_input_validation_report.json`
- `phase1_input_validation_summary.md`
- `pyrosetta_run_metadata.json`
- `pyrosetta_decoy_scores.csv`
- `pyrosetta_interface_models.csv`
- `pyrosetta_interface_residue_table.csv`
- `orientation_filter_log.csv`
- `ppi_cluster_summary.csv`
- `ppi_hotspot_residues.csv`
- `ppi_interface_patch_table.csv`

전 과정에 전파되는 metadata: `receptor_id`, `construct_type`, `orientation_validation_status`

### 6.4 Cross-State 출력 (Phase 1 output root)

- `ppi_patch_cross_state_comparison.csv`
- `ppi_patch_state_robustness.csv`
- `phase1_interface_comparison_report.md`

`construct_type`과 `orientation_validation_status`가 보존됨. Numbering 또는 chain mismatch warning이 있을 경우 input validation 출력을 우선 참조.

### 6.5 Final Phase 1 출력

- `phase1_interface_report.md`
- `phase1_downstream_patch_reference.csv`

Downstream patch reference에 포함: `construct_type`, `orientation_validation_status`, robustness labels, method agreement, confidence. Orientation filtering을 건너뛴 경우 `not_available`로 fallback.

### 6.6 실무 열람 순서

1. `phase1_input_validation_summary.md`
2. `pyrosetta_run_metadata.json`
3. `ppi_cluster_summary.csv`
4. `ppi_patch_state_robustness.csv`
5. `cross_method_convergence.csv`
6. `phase1_interface_report.md`

---

## 7. 핸드오프 (Handoff)

### 7.1 Legacy Postprocess CSVs

Legacy 경로 (`output/workflow_a/phase3_ppi_postprocess/`):

**`ppi_pyrosetta_residues.csv`**: `receptor_id`, `partner_id`, `source`, `chain`, `residue_id`, `residue_num`, `residue_name`, `lobe_label`, `construct_type`, `orientation_validation_status`, `frequency_final_ranking`, `frequency_cluster_summary`, `n_models_final_ranking`, `occupancy`, `mean_interface_delta_e`, `best_interface_delta_e`

**`ppi_pyrosetta_summary.csv`**: `receptor_id`, `partner_id`, `source`, `construct_type`, `orientation_validation_status`, `n_final_models`, `n_clusters`, `n_interface_residues`, `n_nlobe_interface_residues`, `n_clobe_interface_residues`, `top_residues`, `best_dg`, `mean_dg`, `best_dsasa`

### 7.2 `orientation_validation_status` 필드 계약

이 필드는 legacy와 structured Phase 1 양쪽에 모두 존재하며, downstream reader가 단일 필드명을 사용할 수 있도록 한다.

| Source | Value |
|--------|-------|
| Legacy postprocess (`output/workflow_a/phase3_ppi_postprocess/`) | Usually `not_available` |
| Structured Phase 1 (`output/workflow_a/phase2_ppi_docking/`) | `pass`, `fail`, `ambiguous`, `insufficient_data` |

이 차이는 의도적이며 모순이 아님.

### 7.3 Validation Checklist

- [x] Existing pipeline handles full-kinase-domain inputs without modification
- [x] receptor_id, partner_id, construct_type recorded in run metadata
- [x] Config files generated for all 3 states x (1 test + 5 production seeds)
- [x] Output directories separate receptor states and seed indices
- [x] Score standardization utility consolidates multi-seed outputs
- [x] Compute scaling documented with multi-seed recommendation
- [x] All production runs designated as server-side only
