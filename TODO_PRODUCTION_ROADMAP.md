# EGFR-MYO1D PPI 도킹 — 프로덕션 로드맵

Last updated: 2026-04-03

다이머 상태 키나아제 도메인을 사용한 최종 프로덕션 실행을 위한
단계별 체크리스트입니다. 각 Phase는 이전 Phase의 완료를 전제합니다.

---

## Phase 0: 다이머 입력 준비 (프로덕션 실행 전)

### 0-1. 다이머 PDB 준비

- [ ] EGFR 키나아제 다이머 PDB 확보 (3GT8 비대칭 단위 또는 생물학적 단위)
- [ ] Chain A = 수용체, Chain B = 파트너 확인 (또는 A/B = 다이머 두 사슬)
- [ ] 다이머 인터페이스 잔기 확인 → `excluded_residues_a`에 추가
  - 다이머 접촉면은 MYO1D 결합에 사용할 수 없으므로 제외 필수
  - `egfr_pipeline/ppi/prepare_dimer_pdb.py`의 +1000 offset 처리 확인
- [ ] 누락 잔기/루프 확인 → Rosetta로 모델링 또는 제외
- [ ] 양성자화 상태 검토 (EGFR 활성 사이트 HIS 잔기들)
- [ ] MD 클러스터 대표 구조 2-3개도 다이머 상태로 준비

### 0-2. Config 업데이트

- [ ] `config/phase1/phase1_prod_*.ini` 15개 파일 검토:
  ```ini
  [Path]
  input_pdb_name = <다이머 PDB 경로>

  [Constraints]
  # 기존 모노머 제외 잔기 + 다이머 인터페이스 잔기 추가
  excluded_residues_a = 709-720,724-731,...,<다이머 인터페이스 잔기>

  # Ko et al. 실험 데이터 반영 (현재 비어 있음!)
  key_residues_b = 961,962,963,964,968,969,970,971,972
  key_residue_bonus_weight = 2.0
  ```
- [ ] `known_binding_region_a`의 출처 문서화 (현재 confidence=low)
- [ ] `config/example-project.yaml`에 실험 잔기 설정:
  ```yaml
  experimental:
    source: "Ko_et_al_mutagenesis"
    known_binding_residues: [961, 962, 963, 964, 968, 969, 970, 971, 972]
    known_non_binding_residues: [975, 976, 977, 978]  # sheet 10-11 (WT-level)
  ```

### 0-3. Precheck

- [ ] `qsub config/run_pre_qsub_checks.pbs` 실행
- [ ] `output/precheck/last_pass.json` → `"status": "passed"` 확인

---

## Phase 1: 프로덕션 PPI 도킹 (HPC)

### 1-1. 듀얼 시드 제출 (32코어 전체 활용)

```bash
# 3 states × 5 seeds = 9 PBS 잡 (듀얼 시드 전략)
for state in 3GT8_raw EGFR_160-185 EGFR_170-200; do
  qsub -v STATE=$state,SEED_A=0,SEED_B=1 config/run_ppi_dual_seed.pbs
  qsub -v STATE=$state,SEED_A=2,SEED_B=3 config/run_ppi_dual_seed.pbs
  qsub -v STATE=$state,SEED_A=4 config/run_ppi_state_seed.pbs
done
```

- [ ] 9개 잡 제출 완료
- [ ] `qstat` 으로 진행 확인
- [ ] 각 시드의 `PROGRESS.log` tail -f 모니터링

### 1-2. 완료 확인

- [ ] 15개 시드 모두 `seed_complete.json` 존재 확인:
  ```bash
  ls output/workflow_a/phase2_ppi_docking/*/prod_seed*/seed_complete.json | wc -l
  # 기대값: 15
  ```
- [ ] 에러 확인:
  ```bash
  grep 'ERROR\|CRITICAL' output/workflow_a/phase2_ppi_docking/*/prod_seed*/logs/workers.log
  ```
- [ ] 각 시드의 `final_ranking.csv` 존재 확인

---

## Phase 2: 즉시 검증 (기존 데이터로, 재도킹 불필요)

### 2-1. 임계값 민감도 분석

각 시드의 필터 임계값을 ±20% 변동시켜 결과 안정성을 확인합니다.

```bash
# 각 시드에 대해 실행
for state in 3GT8_raw EGFR_160-185 EGFR_170-200; do
  for s in 0 1 2 3 4; do
    python -m egfr_pipeline.posthoc_analysis sensitivity \
      --scored output/workflow_a/phase2_ppi_docking/$state/prod_seed$s/scored_all_models.csv
  done
done
```

- [ ] 전체 시드 실행 완료
- [ ] `sensitivity_analysis.csv` 검토
- [ ] **핵심 확인**: `top20_overlap_frac ≥ 0.75`이면 필터가 강건함
- [ ] 만약 `UNSTABLE` → 해당 임계값 재검토 필요

### 2-2. 시드 간 수렴 분석

5개 시드가 동일한 결합 사이트로 수렴하는지 확인합니다.

```bash
for state in 3GT8_raw EGFR_160-185 EGFR_170-200; do
  python -m egfr_pipeline.posthoc_analysis convergence \
    --state_dir output/workflow_a/phase2_ppi_docking/$state/
done
```

- [ ] 3개 상태 모두 실행 완료
- [ ] `convergence_analysis.csv` 검토
- [ ] **핵심 확인**: `convergence_verdict = CONVERGED` 또는 `PARTIALLY_CONVERGED`
- [ ] `NOT_CONVERGED`이면 → 디코이 수를 50K로 증가하거나 시드 추가 고려

### 2-3. 엔트로피 보정 적용

ref2015의 엔탈피 전용 dG에 결합 엔트로피 페널티를 추가합니다.

```bash
for state in 3GT8_raw EGFR_160-185 EGFR_170-200; do
  for s in 0 1 2 3 4; do
    python -m egfr_pipeline.posthoc_analysis entropy \
      --ranking output/workflow_a/phase2_ppi_docking/$state/prod_seed$s/final_ranking.csv
  done
done
```

- [ ] `final_ranking_entropy_corrected.csv` 생성 확인
- [ ] 보정 후 순위 변동 확인 (큰 변동 = 원래 순위가 약한 결합에 의존)

### 2-4. Verdict 재실행 (4축 스코어링 적용)

```bash
python main.py -c config/example-project.yaml verdict
```

- [ ] `valid_sites.csv`에 `exp_score` 컬럼 확인
- [ ] `decision_trace`에 `EXP[...]` 섹션 확인
- [ ] STRONG/MODERATE/WEAK 분류 검토

---

## Phase 3: 추가 검증 (HPC 필요)

### 3-1. 다중 스코어링 합의 (재도킹 아님, 재스코어링)

기존 PDB를 ref2015 + beta_nov16 + ref2015_cart 3개 함수로 재스코어링합니다.

```bash
for state in 3GT8_raw EGFR_160-185 EGFR_170-200; do
  for s in 0 1 2 3 4; do
    qsub -v STATE=$state,SEED=$s config/run_consensus_scoring.pbs
  done
done
```

- [ ] 15개 잡 제출 및 완료
- [ ] `consensus_scores.csv` 검토
- [ ] **핵심 확인**: `consensus_hit=yes` 비율 ≥ 30%이면 스코어링 간 합의 양호
- [ ] consensus_hit 모델들이 최종 랭킹 상위에 있는지 교차 확인

### 3-2. 디코이 농축 검증 (추가 도킹 필요)

결과가 무작위 대비 통계적으로 유의미한지 검증합니다.

```bash
# Step 1: scrambled PDB 생성 (개발 환경에서 가능)
python -m egfr_pipeline.decoy_enrichment generate \
    --partner_pdb input/PPI/phase1/MYO1D_TH1.pdb \
    --output_dir output/decoy_controls/scrambled_pdbs/ \
    --n_scrambles 5

# Step 2: HPC에서 컨트롤 도킹
for state in 3GT8_raw EGFR_160-185 EGFR_170-200; do
  qsub -v STATE=$state config/run_decoy_docking.pbs
done

# Step 3: 통계 분석 (도킹 완료 후)
for state in 3GT8_raw EGFR_160-185 EGFR_170-200; do
  python -m egfr_pipeline.decoy_enrichment analyse \
    --real_scores output/workflow_a/phase2_ppi_docking/$state/prod_seed0/scored_all_models.csv \
    --control_dir output/decoy_controls/$state/ \
    --output output/decoy_controls/$state/enrichment_report.csv \
    --state $state --seed 0
done
```

- [ ] scrambled PDB 5개 생성 확인
- [ ] 3개 state 컨트롤 도킹 완료
- [ ] `enrichment_report.csv` 검토
- [ ] **핵심 확인**: `enrichment_verdict = significant` (Z < -2, p < 0.05)
- [ ] `no_enrichment`이면 → 결과 해석에 심각한 주의 필요

---

## Phase 4: LightDock 교차 검증 (선택적, 권장)

### 4-1. LightDock 실행

```bash
qsub config/run_lightdock.pbs
```

- [ ] 3개 상태 모두 완료
- [ ] `lightdock_interface_support_table.csv` 생성 확인

### 4-2. Convergence 분석 및 verdict 통합

```bash
# Convergence CSV 생성
for state in 3GT8_raw EGFR_160-185 EGFR_170-200; do
  python -m egfr_pipeline.phase1.lightdock_validation --convergence --state $state
done

# verdict 통합을 위해 convergence CSV를 ppi_postprocess에 복사
cp output/workflow_b/phase1_ppi_analysis/*/lightdock/lightdock_convergence.csv \
   output/workflow_a/phase3_ppi_postprocess/lightdock_convergence.csv

# verdict 재실행
python main.py -c config/example-project.yaml verdict
```

- [ ] `reason_tags`에 `lightdock_validated` 태그 확인
- [ ] PyRosetta와 LightDock이 동일 사이트를 발견하면 → 높은 신뢰도

---

## Phase 5: Workflow A 후처리 및 통합

### 5-1. PPI 후처리

```bash
python main.py -c config/example-project.yaml ppi-postprocess
```

- [ ] `ppi_pyrosetta_residues.csv` 생성
- [ ] `ppi_pyrosetta_summary.csv` 생성

### 5-2. Vina Blind Docking (병렬 실행 가능)

```bash
qsub config/run_vina_cpu.pbs
# 완료 후:
python main.py -c config/example-project.yaml postprocess
```

- [ ] Vina 포즈 테이블 생성
- [ ] 포켓 클러스터링 완료
- [ ] bootstrap 안정성 분석 완료

### 5-3. 최종 판정

```bash
python main.py -c config/example-project.yaml verdict
python main.py -c config/example-project.yaml report
python main.py -c config/example-project.yaml validate
```

- [ ] `valid_sites.csv` — STRONG 사이트 목록 확인
- [ ] `project_report.txt` — 종합 리포트 확인
- [ ] `validation_status.json` — 모든 체크 PASS 확인

---

## Phase 6: 결과 해석 및 논문 준비

### 6-1. 구조 분석 (PyMOL)

```bash
# 서버에서 로컬로 복사
scp -r user@node05:~/codex_ligand/output/workflow_a/phase2_ppi_docking/ ~/Desktop/ppi_results/

# PyMOL 시각화
pymol 1_OVERVIEW_Clusters.pml          # 전체 사이트 분포
pymol final_result/2_DETAIL_C01.pml    # 개별 사이트
pymol final_result/view_results.pml    # 최종 랭킹
```

- [ ] 상위 클러스터의 생물학적 위치 확인:
  - C-lobe 표면인가? (기대되는 MYO1D 결합 위치)
  - 막 근접면에 있지 않은가?
  - 다이머 인터페이스와 겹치지 않는가?
- [ ] Ko et al. sheet 8-9 잔기가 인터페이스에 포함되는지 확인
- [ ] 3개 상태 간 일관된 사이트 확인

### 6-2. 논문 Methods 섹션 체크리스트

- [ ] 에너지 함수: "ref2015 (Rosetta 에너지 함수)" 명시
- [ ] 샘플링: "20,000 decoys × 5 seeds × 3 receptor states = 300,000 total" 명시
- [ ] 필터링: Stage 1/2 임계값과 그 근거 명시
- [ ] 민감도 분석 결과 인용 ("top-20 overlap ≥ X% under ±20% threshold variation")
- [ ] 수렴 분석 결과 인용 ("N consensus sites across 5 seeds")
- [ ] 다이머 상태 사용 근거 명시
- [ ] 엔트로피 보정 방법 설명
- [ ] 제한사항:
  - "ref2015는 PPI 전용이 아닌 단백질 접힘 최적화 에너지 함수"
  - "3개 수용체 상태는 제한적 컨포메이션 앙상블"
  - "rigid-body 도킹으로 파트너 유연성 미반영"
  - "결과는 후보 인터페이스이며, 실험적 검증 필요"

### 6-3. Supplementary 데이터

- [ ] `sensitivity_analysis.csv` → Table S1
- [ ] `convergence_analysis.csv` → Table S2
- [ ] `consensus_scores.csv` (상위 모델) → Table S3
- [ ] `enrichment_report.csv` → Table S4
- [ ] `energy_funnel.png` (대표 시드) → Figure S1
- [ ] `1_OVERVIEW_Clusters.pml` 캡처 → Figure S2

---

## Phase 7: 실험 검증 계획 (연구실)

이 파이프라인의 결과는 **"계산 기반 후보 인터페이스"**입니다.
확정적 결론을 위해 다음 실험이 필요합니다:

### 7-1. 우선순위 1 (필수)

- [ ] **SPR/BLI**: EGFR-MYO1D 결합 친화도 (KD) 측정
- [ ] **알라닌 스캐닝**: 상위 10-20개 예측 EGFR 인터페이스 잔기에 대해
  - 각 잔기를 Ala로 치환
  - MYO1D 결합 변화 측정
  - 예측 vs 실험 일치율 = 파이프라인 검증

### 7-2. 우선순위 2 (강력 권장)

- [ ] **XL-MS (교차결합 질량분석)**: 실제 접촉 토폴로지 결정
- [ ] **HDX-MS**: 결합 시 backbone dynamics 변화 매핑

### 7-3. 우선순위 3 (이상적)

- [ ] **Cryo-EM**: EGFR-MYO1D 복합체 구조 직접 결정
- [ ] **NMR CSP**: 원자 수준 접촉 맵

---

## Quick Reference: 자주 쓰는 명령

```bash
# === 모니터링 ===
tail -f output/workflow_a/phase2_ppi_docking/3GT8_raw/prod_seed0/PROGRESS.log
qstat -u $USER
grep 'ERROR' output/workflow_a/phase2_ppi_docking/*/prod_seed*/logs/workers.log

# === 분석 (개발 환경 / 서버 모두 가능) ===
python -m egfr_pipeline.posthoc_analysis sensitivity --scored <scored_csv>
python -m egfr_pipeline.posthoc_analysis convergence --state_dir <state_dir>
python -m egfr_pipeline.posthoc_analysis entropy --ranking <ranking_csv>

# === verdict (개발 환경 / 서버) ===
python main.py -c config/example-project.yaml verdict
python main.py -c config/example-project.yaml report
python main.py -c config/example-project.yaml validate

# === HPC 제출 ===
qsub -v STATE=3GT8_raw,SEED_A=0,SEED_B=1 config/run_ppi_dual_seed.pbs
qsub -v STATE=3GT8_raw,SEED=0 config/run_consensus_scoring.pbs
qsub -v STATE=3GT8_raw config/run_decoy_docking.pbs
```

---

## 결정 트리: 결과에 따른 분기

```
Phase 2-1 민감도 분석
  ├─ top20_overlap ≥ 75% → 필터 강건, Phase 2-2로 진행
  └─ top20_overlap < 50% → 임계값 재조정 필요 (scored_all_models.csv 재분석)

Phase 2-2 수렴 분석
  ├─ CONVERGED → 샘플링 충분, Phase 3으로 진행
  ├─ PARTIALLY_CONVERGED → 결과 사용 가능하나 논문에 한계 명시
  └─ NOT_CONVERGED → 디코이 50K로 증가 또는 시드 추가 후 Phase 1 재실행

Phase 3-2 디코이 농축
  ├─ significant → 결과 통계적으로 유의미
  ├─ suggestive → 컨트롤 추가 (n_scrambles 10으로)
  └─ no_enrichment → 결과 신뢰 불가 — 근본적 재검토 필요
```
