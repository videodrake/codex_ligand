# Nightly Incremental Improvement Automation

Last updated: 2026-03-17

이 문서는 이 저장소를 매일 새벽마다 아주 작은 단위로 안전하게 개선해 나가기 위한 자동화 설계를 설명한다.

기존 [nightly_review_automation.md](nightly_review_automation.md)가 "깊은 리뷰 번들 생성"에 초점을 둔다면, 이 문서는 그 다음 단계인 "리뷰 -> 미세 수정 -> 검증 -> 커밋 -> 브랜치 푸시"까지 포함하는 운영 루프를 정의한다.

## 목표

자동화의 목표는 다음 한 문장으로 요약할 수 있다.

"매일 새벽, 저장소 전체를 넓게 훑되 한 번에 하나의 작은 변경만 안전하게 적용하고, 검증을 통과한 경우에만 전용 브랜치로 푸시한다."

이 자동화는 대규모 리팩터링이나 HPC 실험 실행기가 아니라, 점진적 품질 개선기이다.

## 왜 이 저장소에 맞는가

이 저장소는 이미 자동화 친화적인 구조를 갖고 있다.

- `scripts/nightly_review.py`가 리뷰 단위를 세분화해서 번들을 생성한다.
- `run_production.py`와 `egfr_pipeline/runtime.py`가 운영 레이어를 분리한다.
- `egfr_pipeline/step_view.py`와 `egfr_pipeline/validate.py`가 출력 계약과 회귀를 강하게 통제한다.
- `tests/`가 orchestration, step view, verdict, phase별 계약을 꽤 넓게 커버한다.

반대로, 자동화가 조심해야 할 부분도 분명하다.

- 실제 도킹/PPI 계산은 HPC `qsub` 환경을 전제로 하므로 새벽 자동화가 로컬에서 실행하면 안 된다.
- `output/`의 canonical 산출물은 해석 대상이지 자동 수정 대상이 아니다.
- scientific threshold, residue numbering, output schema는 작은 수정이라도 파급력이 크다.

즉, 이 저장소에서의 자동화는 "무거운 계산 실행"이 아니라 "코드/문서/테스트/운영 계약의 미세 보수"에 맞춰져야 한다.

## 자동화 이름

권장 이름: `Nightly Micro-Improvement`

권장 실행 시각: 사용자 로컬 시간대의 새벽 1회

## 한 번의 실행에서 해야 하는 일

한 번의 자동화 실행은 아래 둘 중 하나만 성공으로 본다.

1. 작은 변경 1건을 만들고, 관련 검증을 통과한 뒤, 전용 브랜치에 커밋/푸시한다.
2. 안전한 변경 후보가 없다고 판단하고, 이유와 다음 우선 후보를 리뷰 리포트에 남긴 뒤 종료한다.

중요한 점은 "매일 무조건 많이 바꾸는 것"이 아니라 "매일 하나씩이라도 축적 가능한 변경만 남기는 것"이다.

## 브랜치 전략

자동화는 장기 실행 브랜치 하나를 사용한다.

- 브랜치명: `codex/nightly-incremental`
- 수정 시작 전 선행 조건: 항상 최신 `main` 코드를 먼저 가져온다.
- 원칙: 새벽 실행마다 최대 1개 논리 변경 = 최대 1개 커밋
- 동기화 방식: 실행 시작 시 `origin/main` 기준으로 rebase 또는 fast-forward 정렬
- 푸시 방식: 검증 통과 후에만 `origin/codex/nightly-incremental`로 푸시

이 전략을 쓰는 이유는 다음과 같다.

- 하루치 변경을 시간 순으로 추적하기 쉽다.
- 실패 시 범위를 하루 단위로 되돌리기 쉽다.
- "지금 자동화가 어디까지 손댔는가"를 하나의 브랜치에서 보기 쉽다.

## 실행 루프

### 1. Preflight

자동화는 먼저 저장소를 안전 상태로 맞춘다.

- `origin/main` 최신 상태 fetch
- 로컬 `main`을 최신 `origin/main` 기준으로 fast-forward
- 자동화 브랜치를 checkout
- 자동화 브랜치를 최신 `main` 기준으로 rebase
- 충돌 발생 시 변경 없이 종료

중요한 원칙은 "수정 전에 main부터 최신화"이다. 자동화 브랜치가 최신 `main`을 흡수하기 전에는 review unit 선택이나 파일 수정을 시작하지 않는다.

이 단계에서 충돌이 나면 그날은 수정 작업을 하지 않는다. 새벽 자동화는 충돌 해결기가 아니라 유지보수 보조기이기 때문이다.

### 2. 리뷰 번들 생성

기존 도구를 그대로 사용한다.

```bash
python scripts/nightly_review.py --label nightly-auto
```

자동화는 최소한 아래 파일을 읽는다.

- `review_manifest.json`
- `review_checklist.md`
- `review_prompt.md`
- 직전 `nightly_review_report.md`가 있으면 함께 참고

### 3. 대상 선택

한 번의 실행에서는 review unit 하나만 선택한다.

우선순위는 아래 순서를 권장한다.

1. 최근 변경된 파일과 직접 연결된 `critical` unit
2. 테스트는 있는데 문서/계약/예외처리가 약한 unit
3. step view, validation, runtime guard처럼 운영 리스크를 낮추는 unit
4. README/runbook/doc drift 같은 저위험 문서 unit
5. rotation이 오래 밀린 unit

선택 로직은 가능하면 `scripts/nightly_review.py`의 `review_units()` 정의를 그대로 따른다. 별도 우선순위 체계를 이중으로 만들면 drift가 생기기 쉽다.

### 4. 변경 타입 결정

자동화는 후보를 아래 네 가지 타입으로만 분류한다.

### A. 문서 정렬

예시:

- README와 runbook 간 명령어 drift 수정
- docs 인덱스에 누락된 문서 연결
- verdict/report 표현에서 과도한 확정 표현 완화

가장 안전한 시작 타입이다.

### B. 테스트 보강

예시:

- 이미 존재하는 동작을 보호하는 regression test 추가
- step view / runtime / nightly review의 경계 조건 테스트 추가

코드 수정 전의 안전망을 까는 용도다.

### C. 저위험 코드 수정

예시:

- null/empty input 방어
- 잘못된 경로 메시지나 상태 보고 보완
- 문서와 어긋난 작은 동작 회귀 수정

반드시 관련 테스트가 있거나 함께 추가되어야 한다.

### D. 구조 리팩터링

허용은 가능하지만 새벽 자동화의 기본값으로는 비권장이다.

선택 조건:

- 동작 변화가 없어야 한다.
- diff가 작아야 한다.
- 관련 테스트를 즉시 돌릴 수 있어야 한다.

### 5. 가드레일 확인

아래 조건 중 하나라도 걸리면 그 실행은 "수정 없이 리포트만 남김"으로 종료한다.

- `qsub`, 도킹, PyRosetta, LightDock 같은 무거운 계산 실행이 필요한 경우
- output schema 변경이 필요한데 관련 테스트/문서까지 한 번에 손대기 어려운 경우
- residue numbering, scientific threshold, verdict weight처럼 과학적 의미가 큰 상수 수정
- 한 실행에서 2개 이상의 review unit을 건드려야 해결되는 문제
- 변경 범위가 커져서 사실상 기능 개발이 되는 경우

권장 소규모 제한은 다음과 같다.

- 최대 1 review unit
- 최대 3개 핵심 파일 수정
- 가능하면 150라인 내외의 순증/순감
- 코드 수정 시 최소 1개 관련 테스트 실행

### 6. 검증

자동화는 "전체 테스트" 대신 "가장 작은 관련 테스트"를 우선한다.

우선순위:

1. `review_manifest.json` 또는 `review_units()`가 지목한 테스트
2. 수정 파일명과 직접 연결된 단일 pytest 파일
3. orchestration 변경이면 관련 smoke test 1개 추가

예시:

- `step_view.py` 수정: `pytest -m "unit and not slow" tests/unit -q` + 필요 시 `pytest -m integration tests/integration -q`
- `run_production.py` 수정: `pytest -m "unit and not slow" tests/unit -q` + `pytest -m integration tests/integration -q`
- `scripts/nightly_review.py` 수정: `pytest -m "unit and not slow" tests/unit -q`

테스트가 없고 문서만 바꿨다면 테스트를 생략할 수 있지만, 그 이유를 리포트에 남겨야 한다.

### 7. 리포트 기록

각 실행은 아래 산출물을 남긴다.

- review bundle
- 변경 diff
- 실행한 테스트와 결과
- `nightly_review_report.md`의 갱신 내용
- 성공 시 커밋 1개와 원격 브랜치 푸시

리포트에는 최소한 아래 항목이 있어야 한다.

- 이번에 고른 unit
- 왜 이 unit을 골랐는지
- 무엇을 바꿨는지
- 어떤 테스트를 돌렸는지
- 남은 리스크가 무엇인지
- 다음 밤에 보기 좋은 후보가 무엇인지

### 8. Commit and Push

커밋은 한 줄로 의도를 드러내는 쪽이 좋다.

예시:

- `docs: add nightly micro-improvement automation design`
- `tests: cover step view stale manifest recovery`
- `runtime: guard empty scratch workspace path`

#### 푸시 규칙

자동화는 아래 조건을 모두 만족할 때만 푸시한다.

- 워킹트리에 의도한 변경만 존재
- 관련 테스트가 통과했거나 문서-only 변경임이 명확
- 커밋 메시지가 변경 단위를 설명

푸시 후에는 리포트에 원격 브랜치명을 명시한다.

## 첫 도입 단계

이 자동화는 한 번에 full power로 켜기보다 단계적으로 여는 것이 안전하다.

### Stage 0. Review-only

첫 3회 정도는 리뷰 번들 생성과 우선순위 선택만 한다. 실제 수정은 하지 않는다.

목적:

- 어떤 unit이 자주 걸리는지 파악
- 너무 넓은 범위를 잡는 경향이 있는지 확인

### Stage 1. Docs and tests only

그 다음 1주 정도는 문서 정렬과 테스트 보강만 허용한다.

목적:

- 자동화가 저장소 구조를 안정적으로 따라가는지 검증
- 코드를 바꾸지 않고도 축적 가능한 개선을 먼저 확보

### Stage 2. Low-risk code fixes

이후에만 저위험 코드 수정을 허용한다.

조건:

- 관련 테스트가 이미 있거나 함께 추가 가능
- scientific semantics를 바꾸지 않음
- 운영 계약을 깨지 않음

## 이 저장소에서 좋은 초기 대상

새벽 자동화의 초기 대상은 아래 영역이 좋다.

- README / `docs/README.md` / runbook 간 문서 드리프트 정렬
- `scripts/nightly_review.py`와 관련 문서의 계약 일치
- `step_view.py`의 요약/리커버리 문구와 테스트 정렬
- `validate.py`의 경고 메시지 또는 스키마 가드 보강
- `run_production.py`의 status / lane 메시지 개선과 테스트 보강

반대로 초기에는 피하는 게 좋은 영역은 아래와 같다.

- Vina/PPI scoring threshold 조정
- Phase 2~4 scientific scoring 규칙 변경
- 실제 HPC 제출 스크립트의 동작 변경
- 입력 구조물이나 canonical output 파일 수정

## 제안 자동화 프롬프트

아래 프롬프트를 Codex 자동화의 초안으로 사용할 수 있다.

```text
Review this repository and make at most one safe micro-improvement.

Workflow:
1. Fetch origin/main and update local main before doing any review or edits.
2. Rebase the automation branch onto the latest main snapshot.
3. Generate the nightly review bundle with python scripts/nightly_review.py --label nightly-auto.
4. Read the review manifest and checklist.
5. Select exactly one review unit, preferring changed critical files first, then low-risk operational or documentation gaps.
6. Make only a small change that is safe to verify locally. Avoid qsub, heavy docking, PyRosetta execution, and canonical output edits.
7. Run the smallest relevant pytest target when executable code changes. If you only changed docs, say why tests were skipped.
8. Overwrite the nightly review report with findings, change summary, tests run, and residual risks.
9. Commit the change to the automation branch and push only if verification passed.

Constraints:
- Touch only one logical area.
- Prefer docs/tests/guardrails before scientific behavior changes.
- Do not change verdict thresholds, residue numbering logic, or output schemas unless the change is purely corrective and covered by tests.
```

## 성공 기준

이 자동화가 잘 설계되었다고 볼 수 있는 기준은 다음과 같다.

- 며칠이 지나면 브랜치 히스토리에 작은 단위의 읽기 쉬운 커밋이 쌓인다.
- 실패한 날에도 "왜 안 바꿨는지"가 리포트에 남는다.
- step view, validation, docs, test guard처럼 운영 안전성이 먼저 좋아진다.
- scientific core를 함부로 건드리지 않으면서도 저장소 품질이 꾸준히 개선된다.

## 다음 구현 순서

이 문서는 설계 문서이므로, 실제 구현은 아래 순서를 권장한다.

1. 새벽 자동화 프롬프트를 Codex automation으로 등록
2. 첫 주는 review-only 또는 docs/tests only로 제한
3. 브랜치 `codex/nightly-incremental` 운영 시작
4. 1주 후 리포트를 보고 low-risk code fix 허용 범위를 조정
