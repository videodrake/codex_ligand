# Module Separation Analysis (AC-6.1)

> pipeline_manager.py와 vina_executor.py의 내부 구조를 분석하고 분리 필요성을 판단한다.

## 1. pipeline_manager.py (4,156줄)

### 구조
- **클래스 1개:** PipelineManager (line 67-4156, 4,090줄)
- **메서드 37개**, 그 중 11개가 100줄 이상

### 500줄 이상 함수
| 함수 | 줄 수 | 역할 |
|------|-------|------|
| `generate_validation_report()` | 1,036줄 (1218-2253) | PyRosetta 도킹 검증 보고서 생성 |
| `_hotspot_section()` (nested) | 540줄 (1714-2253) | 보고서 hotspot 섹션 (generate_validation_report 내부) |

### 순환 참조
- 없음. 내부 모듈(scoring, movers, pyrosetta_init)과 표준 라이브러리만 import.

### 외부 의존
- main.py, run_production.py, launch_docking.py에서 PipelineManager를 import.

### 분리 후보
1. `generate_validation_report()` + `generate_html_dashboard()` → `pyrosetta_docking/reporting.py` (1,257줄)
2. `_filter_v1()` + `_filter_v2()` → `pyrosetta_docking/model_filter.py` (663줄)
3. `_step4_clustering()` + `_step5_selection_and_save()` → `pyrosetta_docking/clustering.py` (509줄)

## 2. vina_executor.py (2,792줄)

### 구조
- **클래스 없음** (순수 절차형)
- **함수 55개**, 그 중 7개가 100줄 이상

### 500줄 이상 함수
- **없음.** 최대 함수: `interactive_mode()` 378줄.

### 함수 크기 분포
| 구간 | 함수 수 |
|------|---------|
| 0-20줄 | 30 |
| 21-50줄 | 14 |
| 51-100줄 | 4 |
| 101-400줄 | 7 |
| 400줄+ | 0 |

### 순환 참조
- 없음. vina_executor → parse_poses (단방향), parse_poses → derive_docking_seed (역방향이지만 함수 내부 conditional import로 안전).

### 외부 의존
- main.py, run_diverse_docking.py (Phase 3)에서 import.

## 3. 결정: 현재 구조 유지

### 근거

**pipeline_manager.py:**
- 500줄 이상 함수 2개 존재하나, 둘 다 보고서 생성 로직 (`generate_validation_report` + 내부 `_hotspot_section`)
- 보고서 분리는 기술적으로 가능하나, PipelineManager 내부 상태(self.run_results, self.config)에 밀접하게 의존
- 분리 시 10+ 개의 self 참조를 매개변수로 전환해야 하며, 기존 import 경로(main.py, run_production.py)에 wrapper 필요
- **위험 대비 이점이 낮음**: 현재 코드가 정상 동작하고, 보고서 로직은 변경 빈도가 낮음

**vina_executor.py:**
- 500줄 이상 함수 없음 → AC-6.1 분리 기준 미충족
- 55개 함수가 잘 분리되어 있고, 대부분 20줄 이하
- 이미 모듈형 구조 (prepare → run → parse → summarize)

**순환 참조:** 두 파일 모두 없음.

### 향후 검토 조건
다음 중 하나라도 발생하면 분리 재검토:
1. `generate_validation_report()`에 새 기능 추가 요청 시 → reporting.py 분리
2. pipeline_manager.py가 5,000줄 초과 시
3. PipelineManager 클래스에 새 단계(step 6+) 추가 시

## 4. 요약

| 파일 | 줄 수 | 500줄+ 함수 | 순환 참조 | 결정 |
|------|-------|-------------|-----------|------|
| pipeline_manager.py | 4,156 | 2개 (보고서) | 없음 | **유지** — 위험 대비 이점 낮음 |
| vina_executor.py | 2,792 | 0개 | 없음 | **유지** — 이미 모듈형 |
