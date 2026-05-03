# M3 ext20k mini 결과 확인 가이드

이 문서는 완료된 M3 mini focused docking run 결과를 온라인/GitHub에서 확인하거나, HPC에서 직접 열어볼 때 쓰는 안내서입니다.

논문 Methods/Results/Discussion에 넣을 전체 해석은 [`m3_ext20k_clean3_publication_interpretation.md`](m3_ext20k_clean3_publication_interpretation.md)를 참고하세요.

## Run IDs

```bash
M2_RUN_ID=prod_ppi_ext20k_20260502_124707
M3_RUN_ID=m3_ext20k_clean3_20260503_023548
RUN=fresh/runs/$M3_RUN_ID
```

## 최종 상태

완료된 run의 최종 상태는 다음과 같습니다.

```text
last_phase=phase3_compounds PASS: M3-T11 report and handoff completed
```

핵심 M3 결과:

```text
M3-T6 Vina collection: PASS, expected_jobs=36, completed=36, parsed_pose_rows=186
M3-T7 pose attribution: PASS, poses_classified=186, hard_gate_pass=134, rejected=52
M3-T8 pose clustering: PASS, eligible=134, clusters=78, converged=6
M3-T9 anchor convergence: PASS, support_rows=7, anchor_rows=9, evidence_ready=True
M3-T10 evidence tiering: PASS, evidence_rows=11, candidate_rows=11, tier1=2, tier2=0
M3-T11 report cleanup: PASS, reports=11, cleanup_mode=execute, planned=0, deleted=0
```

`status` 명령에서 보이는 `phase_WARN_count`, `phase_FAIL_count`는 수정 전 재시도 이력이 로그에 남아 카운트된 것입니다. 현재 최종 phase는 `PASS`입니다.

## 1. 사람이 먼저 읽을 최종 보고서

아래 파일 4개부터 보면 됩니다.

```bash
less $RUN/report/milestone3_summary.md
less $RUN/report/final_candidate_summary.md
less $RUN/report/reviewer_risk_notes.md
less $RUN/report/cleanup_report.json
```

각 파일의 의미:

- `milestone3_summary.md`: 전체 M3 요약과 최종 상태
- `final_candidate_summary.md`: 최종 Tier 1/2/3 후보 가설 요약
- `reviewer_risk_notes.md`: 해석할 때 주의할 점과 리스크
- `cleanup_report.json`: cleanup 실행 요약

## 2. Tier 1 후보 가설만 보기

최상위 Tier 1 후보 가설만 보기 위한 명령어입니다.

```bash
python - <<'PY'
from pathlib import Path
import csv

run=Path("fresh/runs/m3_ext20k_clean3_20260503_023548")
p=run/"phase3_compounds/tables/final_m3_candidate_hypotheses.csv"

rows=list(csv.DictReader(open(p)))
print("total candidates:", len(rows))
print()

for r in rows:
    if r.get("candidate_tier") == "Tier 1":
        print("== TIER 1 ==")
        print("candidate:", r.get("candidate_hypothesis_id"))
        print("compound:", r.get("compound_public_id"))
        print("pocket:", r.get("pocket_family_id"))
        print("score:", r.get("candidate_priority_score"))
        print("mechanism:", r.get("dominant_mechanism_class"))
        print("states:", r.get("primary_state_ids_supported"))
        print("reason:", r.get("candidate_accept_reason"))
        print("statement:", r.get("candidate_hypothesis_statement"))
        print()
PY
```

## 3. 전체 Tier 분포 보기

전체 후보 row를 compact하게 확인합니다.

```bash
python - <<'PY'
from pathlib import Path
import csv
from collections import Counter

run=Path("fresh/runs/m3_ext20k_clean3_20260503_023548")
rows=list(csv.DictReader(open(run/"phase3_compounds/tables/final_m3_candidate_hypotheses.csv")))

print(Counter(r["candidate_tier"] for r in rows))
print()

for r in rows:
    print(
        r.get("candidate_tier"),
        r.get("compound_public_id"),
        r.get("pocket_family_id"),
        "score="+r.get("candidate_priority_score",""),
        "reason="+r.get("candidate_accept_reason", r.get("candidate_reject_reason","")),
    )
PY
```

## 4. 근거 체인 확인

후보가 왜 올라왔는지, 또는 왜 제외됐는지 보려면 아래 표들을 순서대로 보면 됩니다.

```bash
less $RUN/phase3_compounds/tables/compound_pocket_support.csv
less $RUN/phase3_compounds/tables/compound_anchor_convergence.csv
less $RUN/phase3_compounds/tables/pocket_compound_evidence_table.csv
less $RUN/phase3_compounds/tables/final_m3_candidate_hypotheses.csv
less $RUN/phase3_compounds/tables/rejected_candidate_reasons.csv
```

표의 의미:

- `compound_pocket_support.csv`: compound-pocket별 pose/cluster support
- `compound_anchor_convergence.csv`: compound/state에 걸친 anchor convergence
- `pocket_compound_evidence_table.csv`: tiering에 사용된 통합 근거 row
- `final_m3_candidate_hypotheses.csv`: 최종 후보 가설과 tier
- `rejected_candidate_reasons.csv`: hard gate 또는 quarantine 사유

## 5. QC 상태 한 번에 보기

주요 QC 파일들의 status, blocker, warning만 빠르게 봅니다.

```bash
python - <<'PY'
from pathlib import Path
import json

run=Path("fresh/runs/m3_ext20k_clean3_20260503_023548")
for f in [
    "phase3_compounds/qc/docking_completion_qc.json",
    "phase3_compounds/qc/pose_attribution_qc.json",
    "phase3_compounds/qc/pose_clustering_qc.json",
    "phase3_compounds/qc/anchor_convergence_qc.json",
    "phase3_compounds/qc/final_candidate_gate_qc.json",
    "phase3_compounds/qc/m3_report_qc.json",
]:
    q=json.loads((run/f).read_text())
    print("\n==", f, "==")
    print("status:", q.get("overall_status"))
    print("blockers:", q.get("blockers"))
    print("warnings:", q.get("warnings"))
PY
```

## 해석 기준

- Tier 1: M2 pocket gate, M3 pose attribution, clustering, anchor convergence, evidence tiering을 모두 통과한 최상위 computational hypothesis입니다.
- Tier 2/3: 근거는 있으나 Tier 1보다 state/compound/cluster support가 약하거나 덜 반복적인 가설입니다.
- Reject: ATP migration, membrane/dimer conflict, mapping/retention fail 등 hard gate 또는 quarantine 사유로 제외된 항목입니다.
- 중요한 점: 이 결과는 검증된 inhibitor, binder, drug candidate, clinical claim이 아닙니다. 다음 실험 또는 확장 docking에서 우선 검토할 computational hypothesis입니다.

## 최종 보고용 짧은 요약

```text
M3 mini focused docking pipeline completed successfully.
Final status: phase3_compounds PASS, M3-T11 report and handoff completed.
Run ID: m3_ext20k_clean3_20260503_023548
M2 source: prod_ppi_ext20k_20260502_124707
Evidence tiering: PASS, evidence_rows=11, candidate_rows=11, tier1=2
Report cleanup: PASS, reports=11, cleanup_mode=execute, planned=0, deleted=0
```
