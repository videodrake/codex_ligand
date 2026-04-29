# M1 Completion Rework — Handoff to Codex Agent

작성일: 2026-04-29
대상: Codex (이전에 Tasks 1-9를 구현한 에이전트)
상태: 진행 중 (Phases 0-4 완료, Phases 5-9 남음)
브랜치: `claude/task10` (`main` 65de454에서 분기)

---

## 0. TL;DR

- 너의 작업 (commits up to 65de454, Tasks 1-9)은 그대로 보존돼. 삭제 없음.
- 다만 M1 foundation 일부와 M2 spec 작업이 섞여서 진행됐고, 그 결과 M1 §23 acceptance scorecard가 ~50% (6 DONE / 4 PARTIAL / 5 MISSING)에 머물렀어.
- 사용자 결정: M1을 spec대로 깔끔히 마무리한 다음 M2 실행에 들어가자.
- 새 브랜치 `claude/task10`에서 9개 phase로 세분화해서 진행 중. **Phases 0-4 완료 (5 commits, +9788 lines, 167 tests passing).** Phases 5-9 남음.
- 너의 Tasks 1-9 코드와 테스트는 유지되며, Phase 9에서 input path/schema만 M1 canonical로 정렬할 예정.

---

## 1. 프로젝트 컨텍스트

EGFR-MYO1D fresh workflow. membrane-compatible PPI-to-pocket-to-compound discovery.

두 개의 source-of-truth 문서:

```text
milestone1_foundation_codex_handoff_v0_5.md
egfr_myo1d_overall_implementation_plan_milestones_1_3_v1_0.md
```

repo: `https://github.com/videodrake/codex_ligand`

---

## 2. 무엇이 아쉬웠는지 (솔직하게)

너의 commit history는 많은 걸 해냈어:

```text
Task 1: fresh/ skeleton ✓
Task 2: run context + manifest + logging + preflight ✓
Task 3: structure input contracts QC ✓
Task 4: PPI input preparation + restraints + masks ✓ (이건 사실 M2 spec 영역)
Task 5: real-input readiness bridge ✓ (M2)
Task 6: PPI sampling plan + pose-QC scaffold ✓ (M2)
Task 7: PPI consensus patch summary ✓ (M2)
Task 8: pocket discovery plan ✓ (M2)
Task 9: pocket candidate prioritization ✓ (M2)
```

근데 v1.0 plan §16의 M1 task 시퀀스와 비교하면 두 가지 drift가 있었어.

### 2.1 건너뛴 M1 sub-task

```text
M1 Task 3 (PBS generator + qsub smoke harness)        — fresh/scripts/generate_pbs.py가 placeholder만
M1 Task 4 (cleanup manager)                            — fresh/scripts/cleanup_run.py가 "no files deleted" 메시지만 출력
M1 Task 5 (PDB parser + receptor normalization)        — parser는 있지만 +1000 offset, 669-1014 dockable crop, mapping CSV가 없음
M1 Task 6 (membrane frame + MYO1D construct prep)      — validator는 있지만 generator가 없음
M1 Task 7 (M1 integration test)                        — 없음
```

### 2.2 출력 path/schema가 spec과 다름

Task 4가 만든 출력:

```text
fresh/runs/<run_id>/prepared/egfr/egfr_receptor_normalized.pdb
fresh/runs/<run_id>/prepared/myo1d/MYO1D_sheet8_9_12_core_955_1001.pdb
fresh/runs/<run_id>/prepared/myo1d/MYO1D_ext_beta_meander_955_1006_tail_masked.pdb
```

handoff §14.2가 명시하는 M1 canonical 출력:

```text
fresh/runs/<run_id>/normalized/receptors/<state>_full_frame_explicit_AB.pdb
fresh/runs/<run_id>/normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb
fresh/runs/<run_id>/normalized/receptors/<state>_runtime_offset_receptor_only.pdb
fresh/runs/<run_id>/normalized/myo1d/MYO1D_955_1006.pdb
fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv
fresh/runs/<run_id>/manifest/membrane_frame.json   (state-aware schema)
```

특히 빠진 것:

- 어떤 receptor PDB 변형도 **+1000 runtime offset**을 protomer B에 적용하지 않음 (PyRosetta pose manipulation에서 residue-number reset을 막는 source-of-truth)
- mapping CSV (source residue identity ↔ runtime residue identity 추적용) 없음
- `membrane_frame.json` **생성** 로직 없음 (validation만 있음). 즉 좌표에서 membrane normal이나 X 다이머 축을 계산하는 모듈이 없었음
- M1 raw `MYO1D_955_1006.pdb` 출력 없음 (Task 4가 derivation으로 955-1001 / 955-1006 tail-masked만 생성)
- ligand manifest shell 없음

### 2.3 결과

M1 §23 acceptance scorecard:

```text
6 DONE / 4 PARTIAL / 5 MISSING out of 15 (~50%)

DONE:    #1 skeleton, #2 configs, #3 .gitignore, #4 init-run, #5 logging, #7 preflight, #9 PDB parser
PARTIAL: #10 chain normalization, #12 membrane_frame schema, #13 MYO1D QC, #15 input-prep smoke
MISSING: #6 qsub smoke gen, #8 cleanup, #11 +1000 offset, #14 ligand manifest, M1 Task 7 integration
```

그리고 Tasks 5-9 (M2 spec 작업)는 위 missing input artifact들을 묵시적으로 가정하고 있었어. 예: prepare-ppi-inputs는 normalized/receptors가 있다고 가정하지만, 실제로는 없었음.

이게 의도적인 trade-off였다면 OK였겠지만, v1.0 plan §16의 sequencing과 어긋난다는 점, 그리고 M1→M2 transition gate (§14.1)의 절반 이상이 비어 있다는 점에서 구조적 정합성 문제로 판단했어.

---

## 3. 사용자 결정 (verbatim)

```text
"마일스톤 1을 완성시키기위해 세분화하고 체계화해서 계획 세워서구현해보자"
"이미 만든 파일들중에서 문제가 있는것들은 과감하게 변경하거나 지워도 돼"
"너무 기존코드를 유지하려고 하지 않아도댐"
"전체 계획대로 완성되었으면 좋겠어"
"시간과 자원은 충분함"
```

해석:

- M1을 spec대로 마무리. 패치 위에 패치 금지.
- 정합성 안 맞는 기존 파일은 옮기거나 지워도 OK.
- 시간/자원 충분 → 최대 quality.
- 결과적으로 M1 closure + Tasks 4-9 schema realignment 둘 다 한 브랜치에서 처리.

---

## 4. 새 작업 계획

마스터 플랜: `C:\Users\admin\.claude\plans\1-enchanted-pumpkin.md` (executor 로컬)

브랜치: `claude/task10` (origin/main 65de454에서 분기, Phase 0-4 진행 commits 추가)

10단계 구성 (Phase 0 + 1-9):

```text
Phase 0:  9개 phase prompt + 9개 acceptance checklist 사전 작성 (codex_taskN_*.md 컨벤션 매칭)
Phase 1:  core/cleanup.py + cleanup CLI                              → §23 #8
Phase 2:  preparation/{constructs,pdb_writer}.py → myo1d/ 이동 (logic 무변경)
Phase 3:  myo1d/qc.py 신규 + myo1d/construct.py 확장 + prepare-myo1d → §23 #13
Phase 4:  model/{receptor_normalize,receptor_qc}.py +
          io/residue_mapping.py + prepare-receptor                   → §23 #10, #11
Phase 5:  model/membrane_frame.py + compute-membrane-frame           → §23 #12
Phase 6:  hpc/pbs.py + scripts/{generate_pbs,submit_smoke_*} 실구현
          + prepare-pbs                                               → §23 #6
Phase 7:  ligand/manifest.py + manifest-ligands                      → §23 #14
Phase 8:  prepare-inputs orchestrator + M1 integration smoke test    → §23 #15
Phase 9:  Tasks 4-9 schema realignment (M1 canonical 입력 소비)        → §23에 없음
```

Phase 8 후 M1 §23 scorecard 15/15 (HPC-only 항목 #6 qsub run, #15 real-file run은 사용자 측 검증으로 표기).

---

## 5. 지금까지 완료된 것 (claude/task10 brach commits)

```text
7c14aa2  Phase 0: 18 prompt+checklist docs                          (4509 lines)
f6e31f0  Phase 1: cleanup manager                    (16 tests,      934 lines)
95d27da  Phase 2: MYO1D module relocation            (regression net, 128 lines)
a735dd4  Phase 3: MYO1D construct + QC               (27 tests,     1206 lines)
6efc422  Phase 4: receptor normalization +1000 offset (26 tests,     1605 lines)
```

총: 167 tests passing (98 prior + 69 신규). 11/15 acceptance items closed.

### 5.1 모듈 트리 진척 (handoff §4 spec 대비, HEAD 6efc422 기준)

```text
fresh/src/egfr_myo1d/
├── core/
│   ├── run_context.py            ✓
│   ├── manifest.py               ✓
│   ├── logging_utils.py          ✓
│   └── cleanup.py                ✓ Phase 1
├── hpc/
│   └── pbs.py                    ⏳ Phase 6
├── io/
│   ├── hashing.py                ✓
│   └── residue_mapping.py        ✓ Phase 4
├── ligand/
│   └── manifest.py               ⏳ Phase 7
├── model/
│   ├── receptor_normalize.py     ✓ Phase 4
│   ├── receptor_qc.py            ✓ Phase 4
│   └── membrane_frame.py         ⏳ Phase 5
├── myo1d/
│   ├── construct.py              ✓ Phase 2(이동) + Phase 3(확장)
│   ├── pdb_writer.py             ✓ Phase 2(이동)
│   └── qc.py                     ✓ Phase 3
├── structure/                    ✓ (Task 3, 그대로 유지)
├── preparation/                  ✓ (EGFR-side masks, restraints만 남음)
├── planning/                     ✓ (Task 6, 그대로)
├── analysis/                     ✓ (Tasks 7-9, 그대로)
├── validation/                   ✓ (Phase 9에서 input path만 갱신 예정)
├── cli.py                        ✓ (subparser 4개 추가됨: cleanup, prepare-myo1d, prepare-receptor)
└── __init__.py                   ✓
```

### 5.2 CLI surface 추가됨

```bash
python -m egfr_myo1d.cli cleanup --run-id RUN --mode test|production [--dry-run true|false]
python -m egfr_myo1d.cli prepare-myo1d --run-id RUN --source PATH [--construct 955-1006] [--profile ...]
python -m egfr_myo1d.cli prepare-receptor --run-id RUN --state EGFR_160-185|EGFR_170-200|3GT8_raw --source PATH [--profile ...]
```

기존 11개 subparser는 그대로. 총 14개.

### 5.3 §23 acceptance scorecard 현재 상태

```text
#1  fresh/ skeleton                       ✓ 기존
#2  configs                               ✓ 기존
#3  .gitignore                            ✓ 기존
#4  init-run                              ✓ 기존
#5  logging                               ✓ 기존
#6  qsub smoke generated                  ⏳ Phase 6
#7  preflight                             ✓ 기존
#8  cleanup safe + report                 ✓ Phase 1
#9  PDB parser tests                      ✓ 기존
#10 explicit A/B normalization            ✓ Phase 4
#11 +1000 runtime offset                  ✓ Phase 4
#12 state-aware membrane_frame.json       ⏳ Phase 5
#13 MYO1D 955-1006 construct QC           ✓ Phase 3
#14 ligand manifest shell                 ⏳ Phase 7
#15 input-prep smoke handles missing      ⏳ Phase 8
```

11/15 closed. 4/15 남음.

---

## 6. 남은 phase (5-9) — Codex가 이어받을 경우 참고

각 phase는 사전에 prompt + checklist 문서가 작성되어 있어. 그대로 사용 가능.

```text
fresh/docs/prompts/m1_phase5_membrane_frame_generation_prompt_v0_1.md
fresh/docs/prompts/m1_phase5_membrane_frame_generation_checklist_v0_1.md
fresh/docs/prompts/m1_phase6_pbs_generator_prompt_v0_1.md
fresh/docs/prompts/m1_phase6_pbs_generator_checklist_v0_1.md
fresh/docs/prompts/m1_phase7_ligand_manifest_prompt_v0_1.md
fresh/docs/prompts/m1_phase7_ligand_manifest_checklist_v0_1.md
fresh/docs/prompts/m1_phase8_prepare_inputs_integration_prompt_v0_1.md
fresh/docs/prompts/m1_phase8_prepare_inputs_integration_checklist_v0_1.md
fresh/docs/prompts/m1_phase9_tasks4to9_realignment_prompt_v0_1.md
fresh/docs/prompts/m1_phase9_tasks4to9_realignment_checklist_v0_1.md
```

### 6.1 Phase 5 — membrane frame generation (M1 §23 #12)

handoff §15. **좌표에서 계산** (hardcode 금지). Z+ = C2 axis = membrane normal (extracellular → cytosolic), X+ = chain A → B. 출력:

```text
fresh/runs/<run_id>/manifest/membrane_frame.json   (state-aware: EGFR_160-185, EGFR_170-200, 3GT8_raw)
fresh/runs/<run_id>/qc/membrane_frame_qc.csv
```

3GT8_raw는 `crystallographic_reference_control_not_primary` 표기, vector 필드는 null 또는 alignment-derived only (이번 phase 범위 밖). plus10_full_frame fallback. 둘 다 없으면 `status=missing_frame_source`, vector 발명 금지.

`[0,0,1]`이나 `[1,0,0]` 같은 리터럴이 module source에 fallback 결과로 등장하면 안 됨 (docstring/example 텍스트 안에서만 허용).

### 6.2 Phase 6 — PBS generator (M1 §23 #6)

handoff §10. `hpc/pbs.py` + `prepare-pbs` CLI. 절대 경로 stdout/stderr (`#PBS -o /ABS/REPO/...`), PYTHONPATH=fresh/src export, 5개 thread limit (OMP/OPENBLAS/MKL/NUMEXPR/VECLIB) = 1, conda activate pyrosetta. Mode → ppn 매핑 hpc.yaml (smoke=4, mini=16, scaling=32, prod=32).

`fresh/scripts/{generate_pbs.py, submit_smoke_env.sh, submit_smoke_input.sh}` placeholder를 **실제 emit + qsub-ready** 스크립트로 교체. 자동으로 qsub은 호출 **금지** (사용자가 HPC에서 직접 실행).

골든 PBS 파일을 `fresh/tests/fixtures/m1_phase6_pbs/` 아래 두고 diff 테스트.

### 6.3 Phase 7 — ligand manifest shell (M1 §23 #14)

handoff §17. `ligand/manifest.py` + `manifest-ligands` CLI. **public ID만 (Cpd-A/B/C)**, 내부 ID는 `fresh/data/private/compound_id_map.csv` (gitignored)에만. 출력에 internal ID 누설 시 FAIL. 4가지 profile/stage 매트릭스:

```text
codex_dev / stage=false / missing  → WARN
codex_dev / stage=true  / missing  → WARN
hpc_strict / stage=false / missing → WARN
hpc_strict / stage=true  / missing → FAIL
```

ligand 도킹은 M1 범위 밖. 이 phase는 manifest shell만.

### 6.4 Phase 8 — prepare-inputs orchestrator + M1 integration test (M1 §23 #15)

handoff §10.2 Smoke B + §19 + §22. 단일 CLI가 한 run_id 안에서 sub-step을 순차 실행:

```text
preflight
  → prepare-receptor (per state)
  → compute-membrane-frame
  → prepare-myo1d
  → manifest-ligands
```

Aggregate manifest + Markdown summary 출력. 실제 입력 파일 없으면 `missing_required_inputs` 보고하되 crash 금지. hpc_strict에서 sub-step FAIL 시 즉시 종료.

추가 산출물:

```text
fresh/docs/m1_acceptance_scorecard.md          (15개 항목 표; HPC_PENDING 표기)
fresh/docs/milestone1_foundation_plan.md       (Task 1 stub에서 본격 M1 closure 문서로 확장)
```

테스트 `test_m1_integration_acceptance_scorecard_15_items`: 프로그래매틱하게 §23 1-15 항목을 검사 (DONE 또는 HPC_PENDING).

### 6.5 Phase 9 — Tasks 4-9 schema realignment (post-M1)

§23에 새로 닫는 항목 없음. Tasks 4-9가 M1 canonical 출력을 INPUT으로 소비하도록 갱신:

| Task | 현재 동작 | Phase 9 후 |
|---|---|---|
| Task 4 (`prepared_inputs.py`) | `prepared/egfr/egfr_receptor_normalized.pdb`를 raw input의 pass-through 복사로 emit | M1의 `normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb`를 INPUT으로 참조; pass-through emit 제거. `prepared/myo1d/MYO1D_*` 파생 출력은 유지 (Task 4 derivation) |
| Task 5 (`real_inputs.py`) | raw input만 검사 | M1 normalized 출력 존재도 검증 |
| Task 6 (`ppi_sampling_plan.py`) | job spec receptor_path가 raw input 가리킴 | M1 `runtime_offset_receptor_only.pdb` 참조 |
| Task 7 (`ppi_consensus.py`) | residue identity 자체 검증 | M1 mapping CSV (`<state>_receptor_mapping.csv`)로 protomer/runtime resseq 검증 |
| Task 8 (`pocket_discovery.py`) | receptor_path raw | M1 dockable_*.pdb 참조 |
| Task 9 (`pocket_candidate_prioritization.py`) | 자체 protomer 추론 | M1 mapping CSV로 protomer_id resolve |

기존 79개 Task 4-9 테스트는 fixture path 조정만으로 그대로 통과해야 함. Logic 변경 금지. 추가로 ≥6개 realignment 회귀 테스트.

---

## 7. 새로 도입된 컨벤션

```text
브랜치: claude/task10 (Claude 작업; 기존 codex/<topic> · claude/<topic> 패턴과 동일)
phase prompt:    fresh/docs/prompts/m1_phaseN_<topic>_prompt_v0_1.md
phase checklist: fresh/docs/prompts/m1_phaseN_<topic>_checklist_v0_1.md
phase user doc:  fresh/docs/m1_phaseN_<topic>.md
phase changelog: fresh/docs/m1_phaseN_changes.md
test 파일:       fresh/tests/test_m1_phaseN_<topic>.py
fixture:         fresh/tests/fixtures/m1_phaseN_<topic>/
```

threshold/policy 값은 모두 `fresh/configs/gates.yaml` source-of-truth:

```text
key_residue_bonus_weight: 0.0   (런타임에 != 0이면 RuntimeError)
dockable_crop_default: "669-1014"
excluded_tm_core_default: "634-668"
runtime_offset_second_protomer: 1000
myo1d.construct: "955-1006"
myo1d.key_residues.{sheet8,sheet9,sheet12}: ...
myo1d.watch_residues.{n_terminal,c_terminal}: ...
cleanup.{test_policy, allow_delete_outside_run_dir, production_cleanup_default}: ...
```

Python 스타일: Python 3.7+ 가정, `from __future__ import annotations` 사용, modern type hint는 annotation에만 (런타임 `dict[...]` 사용 회피 — commit 40336dd 패턴 유지).

---

## 8. 너의 기존 작업 (Tasks 1-9) — 어떻게 다뤘는지

### 8.1 그대로 유지된 모듈

```text
core/run_context.py, core/manifest.py, core/logging_utils.py
io/hashing.py
structure/{pdb_parser, contracts, geometry, myo1d_annotation}.py    (Task 3 전부)
validation/preflight.py
validation/structure_inputs.py                                       (Task 3)
preparation/{masks.py, restraints.py}                                (EGFR-side, MYO1D 아님)
planning/{ppi_sampling.py, pose_qc_policy.py}                        (Task 6)
analysis/{ppi_consensus, pocket_selection, pocket_candidate_prioritization}.py  (Tasks 7-9)
validation/ppi_consensus.py, pocket_discovery.py, pocket_candidate_prioritization.py (Tasks 7-9)
모든 fresh/tests/test_taskN_*.py 파일                                (Tasks 1-9 테스트)
모든 fresh/tests/fixtures/taskN_*/ 디렉토리
```

### 8.2 이동된 모듈 (logic 무변경, Phase 2)

```text
preparation/constructs.py    →  myo1d/construct.py
preparation/pdb_writer.py    →  myo1d/pdb_writer.py
```

`git mv`로 history 유지. import path만 4개 파일에서 갱신:

```text
validation/prepared_inputs.py
validation/real_inputs.py
myo1d/construct.py (self-import)
tests/test_task4_ppi_input_preparation.py
```

### 8.3 Phase 9에서 path/schema만 갱신 예정

```text
validation/prepared_inputs.py        (Task 4)
validation/real_inputs.py            (Task 5)
validation/ppi_sampling_plan.py      (Task 6)
validation/ppi_consensus.py          (Task 7)
validation/pocket_discovery.py       (Task 8)
validation/pocket_candidate_prioritization.py  (Task 9)
```

각 모듈의 **분석 로직** (V924R 보고, terminal artifact 감지, active-face 주석, ATP-overlap mask, pose acceptance policy, pocket prioritization scoring 등)은 보존. 입력 path만 raw fixture에서 M1 normalized 출력으로 변경. 기존 79개 테스트는 fixture path만 조정.

### 8.4 Placeholder 교체

```text
fresh/scripts/cleanup_run.py    Phase 1에서 실제 구현으로 교체 (이전 placeholder는 "no files deleted"만 출력)
fresh/scripts/generate_pbs.py   Phase 6에서 교체 예정
fresh/scripts/submit_smoke_env.sh   Phase 6에서 교체 예정
fresh/scripts/submit_smoke_input.sh Phase 6에서 교체 예정
```

---

## 9. 사용자가 명시적으로 원치 않은 것

```text
- 패치 위에 패치 (정합성 안 맞는 상태에서 위에 덮어쓰기)
- Tasks 1-9 분석 로직 삭제
- 3GT8_raw를 primary membrane-validated state로 승격
- V924R을 silent하게 WT (VAL)로 mutate
- MYO1D key residue에 score bonus 추가 (key_residue_bonus_weight는 0.0 고정)
- old workflow 파일 수정: run_production.py, main.py, egfr_pipeline/, config/, docs/runbook.md, output/, results_export/
- 962-start MYO1D 구조를 production partner로 승격 (negative regression fixture로만 유지)
- 3GT8_raw에만 지지된 pocket을 final candidate로 promote
```

---

## 10. HPC-pending (사용자 측 검증)

다음은 Codex 환경이나 Claude 환경에서 검증 불가 — 사용자가 HPC에서 직접 실행:

```text
- bash fresh/scripts/submit_smoke_env.sh        (qsub smoke A)
- bash fresh/scripts/submit_smoke_input.sh      (qsub smoke B)
- 실제 EGFR_160-185.pdb / EGFR_170-200.pdb / 3GT8_raw.pdb / plus10_full_frame.pdb 배치
- 실제 AF-O94832-F1-model_v6.pdb (MYO1D) 배치
- 실제 Cpd-A/B/C SDF + compound_id_map.csv 배치 (private)
- conda env pyrosetta (Python 3.9.25) 실 검증
```

이 항목들은 acceptance scorecard에서 `HPC_PENDING`으로 표기될 예정.

---

## 11. 만약 네가 이어받는다면

1. 이 문서를 먼저 읽기.
2. 현재 브랜치 확인: `git checkout claude/task10`. HEAD는 `6efc422` (Phase 4 commit).
3. Phase 5 prompt 읽기: `fresh/docs/prompts/m1_phase5_membrane_frame_generation_prompt_v0_1.md`. 그 안에 self-contained instruction이 다 들어있음.
4. Phase 5 acceptance criteria: `fresh/docs/prompts/m1_phase5_membrane_frame_generation_checklist_v0_1.md`.
5. 구현 → 테스트 → docs → commit → 다음 phase. 패턴은 Phases 1-4와 동일.
6. **하지 말 것**: Phases 0-4를 다시 하지 말 것. Phase 2에서 이동된 `myo1d/{construct.py, pdb_writer.py}`를 `preparation/`으로 되돌리지 말 것. 새 모듈 (`core/cleanup.py`, `model/*.py`, `myo1d/qc.py`, `io/residue_mapping.py`) 삭제하지 말 것.
7. **테스트 회귀 보장**: 매 phase 끝에 `pytest -q fresh/tests` 통과 확인. 현재 167 passing.

각 phase 끝에서 commit message 컨벤션:

```text
M1 Phase N: <topic> (closes M1 §23 #X[, #Y])

<3-5 줄 요약>

CLI surface:
- python -m egfr_myo1d.cli <new-command> ...

Tests: <count> new (<breakdown>). All <total> tests pass.

Co-Authored-By: <yourself>
```

---

## 12. 참고

```text
저장소:        https://github.com/videodrake/codex_ligand
브랜치:         claude/task10
부모 commit:    65de454 (origin/main; 너의 Tasks 1-9 마지막)
이 문서:        fresh/docs/m1_completion_rework_handoff.md
phase prompt:   fresh/docs/prompts/m1_phaseN_*_prompt_v0_1.md
phase checklist: fresh/docs/prompts/m1_phaseN_*_checklist_v0_1.md
master spec:    milestone1_foundation_codex_handoff_v0_5.md (root)
overall plan:   egfr_myo1d_overall_implementation_plan_milestones_1_3_v1_0.md (root)
```

피드백은 phase prompt 문서 안의 `mandatory end-of-task TODO and self-test block` 컨벤션을 따라줘 — 실제 명령 실행 결과 + 테스트 카운트 + old workflow diff 확인을 명시.

이상.
