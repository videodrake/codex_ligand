## Context Summary
- Project: EGFR–MYO1D structural analysis pipeline for standardized docking, pocket, and residue-level comparison across three receptor states
- Current Phase: Phase 2 (Full Task Breakdown)
- Status: Full task list drafted from approved brief, PRD outline, PRD, and task outline
- Phase Inputs: `brief-egfr-myo1d-pipeline.md`, `prd-outline-egfr-myo1d-pipeline.md`, `prd-egfr-myo1d-pipeline.md`, `tasks-outline-egfr-myo1d-pipeline.md`
- Key Constraints: Existing GitHub codebase must be reused; practical shared-server ceiling is 16 usable CPU cores; tasks must be dependency-ordered, independently testable, and Codex-friendly; legacy residue/site labels are reference-only
- Output Rule: English first, Korean translation after the English task list

---

# Task Breakdown
## EGFR–MYO1D Pipeline

This document breaks the approved PRD into implementation-facing task groups. The purpose of this file is to support incremental repository improvement by Codex or a human developer without losing scientific context or dependency order.

Each task group below includes:
- Objective
- Priority
- Zone
- Main tasks
- Subtasks
- Test tasks
- Dependencies
- Deliverables

The order below is intentional and should be followed unless a later review explicitly changes dependencies.

---

# Task Group 0: Project Setup and Repository Baseline
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Establish a clean baseline understanding of the current GitHub repository, define the working project context, standardize documentation placement, and set up the minimum execution conventions required before feature implementation begins.

## Main Tasks

### 0.1 Inspect the current repository structure
- Identify existing Vina-related scripts
- Identify existing PyRosetta-related scripts
- Identify existing AlphaFold-Multimer-related scripts
- Identify current config patterns, if any
- Identify current output directories and naming patterns
- Identify obvious duplication or dead code regions

### 0.2 Create or update repository-level project context files
- Add the current handoff/specification document to the repository root or docs area
- Ensure the approved brief and PRD documents are stored in a stable location
- Define where future generated reports, parsed tables, and logs should live

### 0.3 Define baseline run assumptions
- Record supported operating environment assumptions
- Record the practical 16-core usable CPU ceiling
- Record that the project is file-based and session-recoverable
- Record that receptor numbering consistency must be preserved whenever possible

### 0.4 Establish output naming and identity conventions
- Define stable receptor IDs
- Define stable ligand IDs
- Define run identifiers or timestamp conventions if needed
- Define rules that prevent receptor-state output mixing

## Subtasks
- Review repository root files and folders
- Inventory runnable scripts and their inputs/outputs
- Inventory existing helper modules
- Draft a baseline repository map for later development
- Draft a naming convention table for receptors, ligands, and outputs

## Test Tasks
- Verify that all core scripts can be located and classified by role
- Verify that the three receptor states can be named consistently in documentation
- Verify that a future developer can understand where to place config, parsed data, and reports

## Dependencies
None.

## Deliverables
- Repository inventory note
- Baseline naming convention note
- Stable placement of brief/PRD/task documents inside repo
- Initial repository map for developers

---

# Task Group 1: Structured Input and Run Management
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Create a consistent project-level input layer for receptor definitions, ligand definitions, metadata, and shared run settings so that all later outputs remain traceable and comparable.

## Main Tasks

### 1.1 Define a project-level configuration format
- Represent receptor definitions in a structured format
- Represent ligand definitions in a structured format
- Represent shared run settings in a structured format
- Allow future extension without breaking existing runs

### 1.2 Implement receptor metadata handling
- Register receptor ID
- Register receptor source type
- Register receptor file paths
- Register chain and relevant notes
- Preserve receptor identity for downstream outputs

### 1.3 Implement ligand metadata handling
- Register ligand ID
- Register ligand file path
- Support optional scientific annotations
- Preserve ligand identity for downstream outputs

### 1.4 Implement configuration validation
- Detect missing receptor files
- Detect missing ligand files
- Detect missing required config fields
- Detect malformed inputs before execution begins

### 1.5 Implement run metadata persistence
- Record project-level run settings
- Record the date or run context if appropriate
- Preserve run parameters for later session recovery

## Subtasks
- Draft a minimal configuration schema
- Draft receptor metadata schema
- Draft ligand metadata schema
- Implement validation checks for required fields
- Implement serialization of run metadata into stable files

## Test Tasks
- Validate that the three receptor states can all be declared in one config
- Validate that multiple ligands can be declared and read correctly
- Validate that invalid file paths trigger clear failures
- Validate that missing required fields trigger clear failures
- Validate that downstream outputs can reference receptor and ligand IDs without ambiguity

## Dependencies
Depends on Task Group 0.

## Deliverables
- Project config template
- Receptor metadata table or export
- Ligand metadata table or export
- Input validation behavior specification
- Run metadata output file

---

# Task Group 2: Parallel Batch Docking Execution
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Implement receptor × ligand docking batch execution with configurable parallel workers, safe operation on a shared server, and stable output structure that remains traceable across runs.

## Main Tasks

### 2.1 Refactor or wrap the existing Vina batch runner
- Reuse the existing GitHub codebase rather than replacing it from scratch
- Support receptor × ligand job matrix execution
- Preserve receptor and ligand identities in outputs

### 2.2 Add configurable parallel execution support
- Add user-defined max worker count
- Use the practical 16-core ceiling as the expected routine upper bound
- Ensure worker configuration is not hard-coded
- Prevent silent oversubscription assumptions

### 2.3 Standardize run logging
- Record per-job success/failure visibility
- Record input receptor and ligand identities per job
- Make failures inspectable instead of silent

### 2.4 Standardize raw output placement
- Ensure receptor-specific outputs are separated
- Ensure ligand-specific outputs are traceable
- Avoid output file collisions during parallel execution

### 2.5 Preserve sequential/parallel equivalence
- Ensure output structure remains logically identical whether run sequentially or in parallel

## Subtasks
- Review current Vina execution entry points
- Identify where concurrency can be introduced safely
- Add worker-limit parameter handling
- Add job queue or dispatch logic
- Add output path isolation rules
- Add job-level logging structure

## Test Tasks
- Run a minimal sequential batch and confirm output placement
- Run the same minimal batch with more than one worker and confirm output equivalence
- Run a larger batch with worker limit set below 16 and confirm correct throttling behavior
- Confirm that a failed job is visible in logs and does not silently disappear
- Confirm that receptor outputs are not mixed under parallel execution

## Dependencies
Depends on Task Group 1.

## Deliverables
- Refactored batch docking runner
- Worker-limit support
- Job logging behavior
- Stable raw output directory convention
- Parallel execution usage documentation

---

# Task Group 3: Pose Parsing and Contact Extraction
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Transform raw Vina outputs into a standardized pose-level dataset with receptor identity, ligand identity, pose metrics, centroid coordinates, and receptor contact residue information.

## Main Tasks

### 3.1 Parse raw docking output files
- Extract pose rank
- Extract affinity or equivalent pose score
- Extract pose-level coordinates needed for centroid calculation
- Associate each pose with receptor and ligand identity

### 3.2 Compute pose centroids
- Generate reproducible pose centroid coordinates
- Store centroid data in a reusable parsed format

### 3.3 Extract receptor contact residues for each pose
- Implement reproducible contact extraction logic
- Preserve residue identity in a standard format
- Preserve raw contact evidence for later review

### 3.4 Build a standardized pose-level output table
- Include receptor ID
- Include ligand ID
- Include pose rank
- Include pose score
- Include centroid coordinates
- Include contact residue summary

### 3.5 Preserve traceability back to raw output
- Retain references to original raw result files where appropriate
- Ensure that parsed outputs do not become detached from source files

## Subtasks
- Review current Vina output format assumptions in existing code
- Define pose-level field set
- Implement centroid extraction function
- Implement residue-contact extraction function
- Implement pose table generation logic
- Define standard residue string formatting

## Test Tasks
- Confirm that a known Vina result produces the expected number of parsed poses
- Confirm that each pose is assigned the correct receptor and ligand identity
- Confirm that centroid values are stored for every parsed pose
- Confirm that contact residues are extracted consistently from repeated runs
- Confirm that parsed pose data can be reviewed without reopening the raw docking file

## Dependencies
Depends on Task Group 2.

## Deliverables
- Standardized pose parsing logic
- Pose-level parsed output table
- Contact residue extraction logic
- Pose centroid extraction logic
- Parsing validation examples

---

# Task Group 4: Pocket Clustering and Pocket Summary Generation
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Convert receptor-local pose collections into receptor-local pockets using the configured centroid rule, then summarize those pockets in a form that preserves residue-level evidence and ligand distribution.

## Main Tasks

### 4.1 Implement receptor-local pocket assignment
- Group poses into pockets within each receptor state
- Use the configured centroid-based grouping rule
- Keep pocket assignment reproducible

### 4.2 Compute pocket-level summary fields
- Pose count per pocket
- Ligand count per pocket
- Pocket centroid
- Best and summary pose metrics
- Residue union and/or residue frequency

### 4.3 Select representative poses for pockets
- Choose a representative pose per pocket for later review
- Preserve the identity of the representative ligand and pose

### 4.4 Build a standardized pocket table
- Represent each receptor-local pocket as a stable record
- Preserve enough detail for later cross-receptor comparison

### 4.5 Build a ligand-to-pocket mapping summary
- Identify dominant pocket per ligand per receptor
- Preserve alternative pocket information where relevant

## Subtasks
- Define pocket ID naming rules within each receptor
- Implement clustering from pose centroids
- Implement pocket centroid generation
- Implement residue frequency summary logic
- Implement representative pose selection rule
- Implement ligand-to-pocket mapping logic

## Test Tasks
- Confirm that pocket assignment is stable when the same parsed pose table is reprocessed
- Confirm that pocket summaries include pose count and ligand count correctly
- Confirm that representative-pose selection is deterministic
- Confirm that ligand-to-pocket mapping can distinguish dominant from non-dominant pockets
- Confirm that receptor-local pockets remain clearly separated across receptor states

## Dependencies
Depends on Task Group 3.

## Deliverables
- Pocket assignment logic
- Pocket summary table
- Ligand-to-pocket mapping table
- Representative-pose selection behavior
- Pocket-level validation examples

---

# Task Group 5: Cross-Receptor Pocket Comparison
**Priority:** Must-Have  
**Zone:** 🟡 Yellow

## Objective
Compare pockets across receptor states using raw comparison metrics that help determine whether two pockets may represent the same patch, a partial overlap, or distinct pockets.

## Main Tasks

### 5.1 Define pocket comparison metrics
- Compare pocket location across receptor states
- Compare residue overlap across receptor states
- Compare ligand evidence across receptor states

### 5.2 Build pocket-to-pocket comparison logic
- Compare every relevant pocket pair across receptor states
- Preserve raw metrics rather than only a binary classification

### 5.3 Support same-patch candidate interpretation without over-automation
- Allow comparison outputs to support human interpretation
- Avoid forcing a single black-box conclusion without raw evidence

### 5.4 Build a standardized cross-receptor comparison table
- Store comparison metrics in a stable reusable format
- Make it possible to inspect comparison results without raw docking files

## Subtasks
- Define comparison field set
- Implement centroid-distance comparison
- Implement residue-overlap comparison
- Implement ligand-overlap or ligand-sharing comparison if applicable
- Define optional same-patch candidate flags as auxiliary labels only

## Test Tasks
- Confirm that comparison results can be generated for all receptor-state pairs
- Confirm that raw comparison metrics are preserved in output
- Confirm that obviously distinct pockets remain distinguishable in the comparison output
- Confirm that the output can support manual scientific review without requiring raw file reinspection

## Dependencies
Depends on Task Group 4.

## Deliverables
- Cross-receptor pocket comparison logic
- Comparison output table
- Optional same-patch candidate support behavior
- Pocket comparison review examples

---

# Task Group 6: Supporting PPI Output Standardization
**Priority:** Should-Have  
**Zone:** 🟡 Yellow

## Objective
Standardize residue-level outputs from PyRosetta global docking and AlphaFold-Multimer so they can be reviewed as auxiliary receptor-side evidence alongside Vina-derived pocket summaries.

## Main Tasks

### 6.1 Inspect the current PyRosetta output structure
- Identify score output forms
- Identify model or decoy outputs
- Identify whether interface residues are already extractable

### 6.2 Standardize PyRosetta residue-side outputs
- Extract receptor-side interface residue information
- Preserve partner-side information separately if useful
- Build reviewable summary outputs

### 6.3 Inspect the current AlphaFold-Multimer output structure
- Identify available model summary outputs
- Identify available confidence summaries
- Identify whether receptor-side contact residues are extractable

### 6.4 Standardize AlphaFold-Multimer residue-side outputs
- Extract receptor-side contact/interface residues
- Preserve model-level context needed for interpretation
- Build reviewable summary outputs

### 6.5 Keep PPI outputs explicitly auxiliary
- Ensure PPI outputs are clearly labeled as supporting evidence rather than pocket truth definitions

## Subtasks
- Inventory current PyRosetta parsing possibilities
- Inventory current AFM parsing possibilities
- Define receptor-side residue summary format
- Define model/cluster summary format for both modules
- Add stable output naming for PPI summaries

## Test Tasks
- Confirm that at least one PyRosetta result can be turned into a receptor-side residue summary
- Confirm that at least one AFM result can be turned into a receptor-side residue summary
- Confirm that the resulting outputs are readable without manual raw-file reconstruction
- Confirm that Vina and PPI outputs remain conceptually separated in naming and documentation

## Dependencies
Depends on Task Group 0 for repository inspection and later report usage in Task Group 7.

## Deliverables
- PyRosetta residue-summary output format
- AFM residue-summary output format
- Parsing or extraction notes for both modules
- Auxiliary-evidence labeling convention

---

# Task Group 7: Reporting and Manual Review Exports
**Priority:** Should-Have  
**Zone:** 🟡 Yellow

## Objective
Generate readable summaries and manual-review helper outputs that allow the researcher and collaborators to inspect pocket and residue evidence without reconstructing the workflow from raw files.

## Main Tasks

### 7.1 Build receptor-level and project-level summary reports
- Summarize receptor-local pockets
- Summarize ligand-to-pocket mappings
- Summarize cross-receptor pocket comparison highlights

### 7.2 Integrate auxiliary PPI summaries into reports
- Include PyRosetta residue-side summaries when available
- Include AFM residue-side summaries when available
- Preserve the distinction between primary and auxiliary evidence

### 7.3 Export manual-review helper artifacts
- Export representative poses for visual inspection where useful
- Export helper summaries that reduce manual lookup effort

### 7.4 Preserve handoff readability
- Ensure reports support future Codex or collaborator use without full scientific rebriefing

## Subtasks
- Define minimum report sections
- Define project-level summary table layout
- Define ligand summary table layout
- Define pocket comparison highlight section
- Define optional manual-review export rules

## Test Tasks
- Confirm that a report can be generated from the core Vina outputs alone
- Confirm that auxiliary PPI outputs can be attached when available
- Confirm that the report is readable without opening raw pose files
- Confirm that the report does not collapse raw evidence into unsupported conclusions

## Dependencies
Depends on Task Groups 4, 5, and optionally 6.

## Deliverables
- Project summary report
- Pocket summary section
- Ligand-to-pocket summary section
- Cross-receptor comparison highlights
- Optional manual-review export set

---

# Task Group 8: Validation, Regression Checks, and Handoff Readiness
**Priority:** Must-Have  
**Zone:** 🟡 Yellow

## Objective
Add validation steps, consistency checks, regression protections, and repository-readiness tasks so that the pipeline can be safely resumed, reviewed, and extended by future sessions or coding agents.

## Main Tasks

### 8.1 Validate output consistency
- Confirm that core outputs exist after expected runs
- Confirm that receptor IDs and ligand IDs remain consistent across output files
- Confirm that parsed and summarized outputs remain connected to raw sources

### 8.2 Add regression-oriented checks
- Ensure that key outputs do not silently change structure across refactors
- Ensure that naming conventions remain stable
- Ensure that parallel execution does not change logical output content

### 8.3 Check residue numbering consistency assumptions
- Detect or flag receptor numbering mismatches where possible
- Preserve warnings where direct comparison may be unsafe

### 8.4 Improve handoff readiness
- Ensure the repository contains the approved brief, PRD, and task documents
- Ensure future developers can locate generated outputs
- Ensure project context is recoverable without prior chat history

## Subtasks
- Define core output checklist
- Define consistency checks for receptor/ligand identity
- Define minimal regression test targets
- Define residue numbering mismatch warning behavior
- Define handoff-readiness checklist

## Test Tasks
- Confirm that all required output tables are present after a standard run
- Confirm that future sessions can identify the latest relevant docs and outputs
- Confirm that a change in one module does not silently break downstream file expectations
- Confirm that numbering mismatch conditions are surfaced to the user when possible

## Dependencies
Depends on all previous task groups.

## Deliverables
- Output validation checklist
- Regression check targets
- Numbering-consistency warning behavior
- Handoff-readiness checklist
- Repository continuation note

---

# Suggested Initial Codex Execution Order

If these tasks are handed to Codex incrementally, the recommended first implementation batch is:

1. **Task Group 0** — repository inspection and baseline conventions  
2. **Task Group 1** — structured input and run management  
3. **Task Group 2** — parallel Vina batch execution  
4. **Task Group 3** — pose parsing and contact extraction  
5. **Task Group 4** — pocket clustering and summary generation

This first batch is enough to create the core Vina-centered MVP evidence layer.

The recommended second implementation batch is:

6. **Task Group 5** — cross-receptor pocket comparison  
7. **Task Group 7** — reporting and manual review exports

The recommended third batch is:

8. **Task Group 6** — supporting PPI output standardization  
9. **Task Group 8** — validation, regression checks, and handoff readiness

This ordering reflects the actual research priority: Vina-based pocket standardization first, comparison second, auxiliary evidence standardization after the core flow is stable.

---

# Korean Translation (한글 번역)

## EGFR–MYO1D Pipeline 작업 분해 문서

이 문서는 승인된 PRD를 실제 구현 가능한 task group으로 분해한 것이다. 목적은 Codex 또는 개발자가 과학적 맥락과 의존성 순서를 잃지 않고 저장소를 점진적으로 개선할 수 있게 하는 것이다.

각 task group은 다음을 포함한다.
- 목표
- 우선순위
- Zone
- 메인 task
- 서브태스크
- 테스트 태스크
- 의존성
- 예상 산출물

아래 순서는 의도된 것이며, 이후 검토에서 명시적으로 바꾸지 않는 한 그대로 따라야 한다.

---

# Task Group 0: Project Setup and Repository Baseline
**우선순위:** Must-Have  
**Zone:** 🟢 Green

## 목표
현재 GitHub 저장소를 명확히 파악하고, 프로젝트 문맥을 정리하고, 문서 위치를 표준화하고, feature 구현 전에 필요한 최소 실행 규칙을 세운다.

## 메인 태스크

### 0.1 현재 저장소 구조 점검
- 기존 Vina 관련 스크립트 식별
- 기존 PyRosetta 관련 스크립트 식별
- 기존 AlphaFold-Multimer 관련 스크립트 식별
- 현재 config 패턴 유무 확인
- 현재 output 디렉터리와 naming 패턴 확인
- 중복 코드나 dead code 후보 확인

### 0.2 저장소 수준 프로젝트 컨텍스트 파일 정리
- 현재 handoff/specification 문서를 repo root 또는 docs 영역에 추가
- 승인된 brief와 PRD 문서가 안정된 위치에 저장되도록 정리
- 향후 report, parsed table, logs 저장 위치 정의

### 0.3 기본 실행 가정 정리
- 지원 운영 환경 가정 기록
- 실사용 가능 16코어 ceiling 기록
- 프로젝트가 file-based이며 session-recoverable이라는 점 기록
- receptor numbering consistency가 가능한 한 유지되어야 함을 기록

### 0.4 output naming 및 identity 규칙 정리
- 안정적인 receptor ID 정의
- 안정적인 ligand ID 정의
- 필요 시 run identifier 또는 timestamp 규칙 정의
- receptor 상태 간 output 섞임 방지 규칙 정의

## 서브태스크
- repo root 파일/폴더 검토
- 실행 가능한 스크립트와 입출력 인벤토리 작성
- 기존 helper module 인벤토리 작성
- 이후 개발용 baseline repository map 초안 작성
- receptor, ligand, output naming convention 표 작성

## 테스트 태스크
- 핵심 스크립트를 모두 찾고 역할별 분류가 가능한지 확인
- 세 receptor 상태를 문서에서 일관되게 이름 붙일 수 있는지 확인
- 미래 개발자가 config, parsed data, report 위치를 이해할 수 있는지 확인

## 의존성
없음.

## 예상 산출물
- 저장소 인벤토리 노트
- baseline naming convention 노트
- repo 내부에 brief/PRD/task 문서 안정 배치
- 개발자용 초기 저장소 맵

---

# Task Group 1: Structured Input and Run Management
**우선순위:** Must-Have  
**Zone:** 🟢 Green

## 목표
receptor 정의, ligand 정의, metadata, shared run setting을 위한 일관된 프로젝트 입력 계층을 만들어 모든 후속 output이 추적 가능하고 비교 가능하도록 한다.

## 메인 태스크

### 1.1 프로젝트 수준 config 형식 정의
- receptor 정의를 구조화 형식으로 표현
- ligand 정의를 구조화 형식으로 표현
- shared run setting을 구조화 형식으로 표현
- 향후 확장이 기존 실행을 깨지 않도록 설계

### 1.2 receptor metadata 처리 구현
- receptor ID 등록
- receptor source type 등록
- receptor file path 등록
- chain 및 관련 note 등록
- downstream output에 receptor identity 보존

### 1.3 ligand metadata 처리 구현
- ligand ID 등록
- ligand file path 등록
- 선택적 scientific annotation 지원
- downstream output에 ligand identity 보존

### 1.4 config validation 구현
- receptor file 누락 감지
- ligand file 누락 감지
- 필수 config field 누락 감지
- 실행 전 malformed input 감지

### 1.5 run metadata 저장 구현
- 프로젝트 수준 실행 설정 기록
- 필요 시 날짜 또는 run context 기록
- 이후 session recovery를 위한 run parameter 저장

## 서브태스크
- 최소 config schema 초안 작성
- receptor metadata schema 초안 작성
- ligand metadata schema 초안 작성
- 필수 field validation 구현
- run metadata를 stable file로 직렬화하는 기능 구현

## 테스트 태스크
- 세 receptor 상태를 한 config에 선언할 수 있는지 검증
- 여러 ligand가 정상적으로 선언되고 읽히는지 검증
- 잘못된 file path가 명확한 실패를 내는지 검증
- 필수 field 누락 시 명확한 실패를 내는지 검증
- downstream output이 receptor/ligand ID를 모호함 없이 참조하는지 검증

## 의존성
Task Group 0에 의존.

## 예상 산출물
- 프로젝트 config template
- receptor metadata table/export
- ligand metadata table/export
- input validation 동작 명세
- run metadata output file

---

# Task Group 2: Parallel Batch Docking Execution
**우선순위:** Must-Have  
**Zone:** 🟢 Green

## 목표
receptor × ligand docking batch 실행을 구현하고, 설정 가능한 병렬 worker와 공유 서버 환경에서의 안전한 실행을 지원하며, 결과 구조가 항상 추적 가능하도록 한다.

## 메인 태스크

### 2.1 기존 Vina batch runner 리팩터링 또는 래핑
- 기존 GitHub 코드베이스 재사용
- receptor × ligand job matrix 실행 지원
- output에 receptor/ligand identity 보존

### 2.2 설정 가능한 병렬 실행 지원 추가
- 사용자 정의 max worker count 추가
- 실질적 16코어 ceiling을 routine upper bound로 사용
- worker 설정 하드코딩 금지
- silent oversubscription 가정 금지

### 2.3 run logging 표준화
- job별 성공/실패 기록
- 각 job의 receptor/ligand identity 기록
- 실패가 숨겨지지 않도록 구성

### 2.4 raw output 배치 표준화
- receptor별 output 분리
- ligand별 output 추적 가능성 확보
- 병렬 실행 중 output file collision 방지

### 2.5 순차/병렬 동등성 유지
- sequential/parallel 여부와 관계없이 논리적 output 구조 동일 유지

## 서브태스크
- 현재 Vina 실행 진입점 검토
- concurrency 도입 가능 지점 확인
- worker-limit parameter 처리 추가
- job queue 또는 dispatch logic 추가
- output path isolation 규칙 추가
- job-level logging 구조 추가

## 테스트 태스크
- 최소 sequential batch 실행 후 output placement 확인
- 같은 batch를 worker > 1로 실행 후 output equivalence 확인
- worker limit를 16 이하로 설정한 큰 batch 실행 후 throttling 동작 확인
- 실패 job이 로그에 보이고 숨겨지지 않는지 확인
- 병렬 실행 중 receptor output이 섞이지 않는지 확인

## 의존성
Task Group 1에 의존.

## 예상 산출물
- 리팩터링된 batch docking runner
- worker-limit 지원
- job logging 동작
- 안정적인 raw output 디렉터리 규칙
- 병렬 실행 사용 문서

---

# Task Group 3: Pose Parsing and Contact Extraction
**우선순위:** Must-Have  
**Zone:** 🟢 Green

## 목표
raw Vina output을 receptor identity, ligand identity, pose metric, centroid coordinate, receptor contact residue를 포함하는 표준 pose-level 데이터셋으로 변환한다.

## 메인 태스크

### 3.1 raw docking output parsing
- pose rank 추출
- affinity 또는 해당 pose score 추출
- centroid 계산을 위한 pose-level coordinate 추출
- pose와 receptor/ligand identity 연결

### 3.2 pose centroid 계산
- 재현 가능한 centroid coordinate 생성
- centroid를 재사용 가능한 parsed format으로 저장

### 3.3 receptor contact residue 추출
- 재현 가능한 contact extraction logic 구현
- residue identity를 표준 형식으로 보존
- 이후 검토를 위한 raw contact evidence 보존

### 3.4 표준 pose-level output table 생성
- receptor ID 포함
- ligand ID 포함
- pose rank 포함
- pose score 포함
- centroid coordinate 포함
- contact residue summary 포함

### 3.5 raw output과의 traceability 보존
- 필요 시 원본 raw result file 참조 정보 유지
- parsed output이 source file와 분리되지 않도록 함

## 서브태스크
- 기존 code의 Vina output format 가정 검토
- pose-level field set 정의
- centroid extraction function 구현
- residue-contact extraction function 구현
- pose table generation logic 구현
- standard residue string formatting 정의

## 테스트 태스크
- 알려진 Vina 결과에서 기대한 수의 parsed pose가 나오는지 확인
- 각 pose가 올바른 receptor/ligand identity를 가지는지 확인
- 모든 pose에 centroid 값이 저장되는지 확인
- repeated run에서 contact residue가 일관되게 추출되는지 확인
- raw docking file을 다시 열지 않고 parsed pose data를 검토할 수 있는지 확인

## 의존성
Task Group 2에 의존.

## 예상 산출물
- 표준 pose parsing logic
- pose-level parsed output table
- contact residue extraction logic
- pose centroid extraction logic
- parsing validation 예시

---

# Task Group 4: Pocket Clustering and Pocket Summary Generation
**우선순위:** Must-Have  
**Zone:** 🟢 Green

## 목표
receptor-local pose 묶음을 설정된 centroid rule로 pocket으로 변환하고, residue-level evidence와 ligand 분포를 보존하는 pocket summary를 생성한다.

## 메인 태스크

### 4.1 receptor-local pocket assignment 구현
- 각 receptor 상태 내부에서 pose를 pocket으로 그룹화
- 설정된 centroid-based grouping rule 사용
- pocket assignment 재현성 유지

### 4.2 pocket-level summary field 계산
- pocket별 pose count
- pocket별 ligand count
- pocket centroid
- best 및 summary pose metric
- residue union 및/또는 residue frequency

### 4.3 representative pose 선택
- 이후 검토를 위한 pocket 대표 pose 선택
- representative ligand 및 pose identity 보존

### 4.4 표준 pocket table 생성
- 각 receptor-local pocket을 안정적인 record로 표현
- 이후 cross-receptor comparison에 충분한 디테일 유지

### 4.5 ligand-to-pocket mapping summary 생성
- ligand별 dominant pocket 식별
- relevant한 경우 alternative pocket 정보 보존

## 서브태스크
- receptor별 pocket ID naming rule 정의
- pose centroid 기반 clustering 구현
- pocket centroid 생성 구현
- residue frequency summary logic 구현
- representative pose selection rule 구현
- ligand-to-pocket mapping logic 구현

## 테스트 태스크
- 같은 parsed pose table을 다시 처리해도 pocket assignment가 안정적인지 확인
- pocket summary가 pose count와 ligand count를 올바르게 포함하는지 확인
- representative-pose selection이 deterministic한지 확인
- ligand-to-pocket mapping이 dominant와 non-dominant pocket을 구분하는지 확인
- receptor 상태 간 receptor-local pocket이 명확히 분리되는지 확인

## 의존성
Task Group 3에 의존.

## 예상 산출물
- pocket assignment logic
- pocket summary table
- ligand-to-pocket mapping table
- representative-pose selection 동작
- pocket-level validation 예시

---

# Task Group 5: Cross-Receptor Pocket Comparison
**우선순위:** Must-Have  
**Zone:** 🟡 Yellow

## 목표
위치, residue overlap, shared ligand 근거를 사용해 receptor 상태 간 pocket을 비교하고, 같은 patch/부분 overlap/독립 pocket 여부를 사람이 판단할 수 있게 raw metric을 제공한다.

## 메인 태스크

### 5.1 pocket comparison metric 정의
- receptor 상태 간 pocket location 비교
- receptor 상태 간 residue overlap 비교
- receptor 상태 간 ligand evidence 비교

### 5.2 pocket-to-pocket comparison logic 구현
- 관련 pocket pair를 모두 비교
- binary classification이 아니라 raw metric 보존

### 5.3 same-patch candidate 보조 해석 지원
- comparison output이 human interpretation을 돕도록 설계
- raw evidence 없이 black-box conclusion을 강요하지 않음

### 5.4 표준 cross-receptor comparison table 생성
- comparison metric을 안정적인 재사용 형식으로 저장
- raw docking file 없이도 비교 결과를 검토할 수 있게 함

## 서브태스크
- comparison field set 정의
- centroid-distance comparison 구현
- residue-overlap comparison 구현
- 적용 가능하면 ligand-overlap/shared ligand comparison 구현
- same-patch candidate flag는 auxiliary label로만 정의

## 테스트 태스크
- 모든 receptor-state pair에 대해 comparison result 생성 가능 여부 확인
- raw comparison metric이 output에 보존되는지 확인
- 명백히 다른 pocket이 comparison output에서 구분되는지 확인
- raw file 재확인 없이 scientific review가 가능한지 확인

## 의존성
Task Group 4에 의존.

## 예상 산출물
- cross-receptor pocket comparison logic
- comparison output table
- optional same-patch candidate support behavior
- pocket comparison review examples

---

# Task Group 6: Supporting PPI Output Standardization
**우선순위:** Should-Have  
**Zone:** 🟡 Yellow

## 목표
PyRosetta global docking과 AlphaFold-Multimer의 residue-level output을 표준화하여, Vina 기반 pocket summary와 나란히 볼 수 있는 auxiliary receptor-side evidence로 만든다.

## 메인 태스크

### 6.1 현재 PyRosetta output 구조 점검
- score output 형식 확인
- model/decoy output 형식 확인
- interface residue 추출 가능성 확인

### 6.2 PyRosetta residue-side output 표준화
- receptor-side interface residue 정보 추출
- 필요 시 partner-side 정보 분리 저장
- 검토 가능한 summary output 생성

### 6.3 현재 AFM output 구조 점검
- model summary output 확인
- confidence summary 확인
- receptor-side contact residue 추출 가능성 확인

### 6.4 AFM residue-side output 표준화
- receptor-side contact/interface residue 추출
- 해석에 필요한 model-level context 보존
- 검토 가능한 summary output 생성

### 6.5 PPI output의 auxiliary 성격 명시
- PPI output이 pocket truth definition이 아니라 supporting evidence로 명확히 표기되도록 함

## 서브태스크
- 현재 PyRosetta parsing 가능성 인벤토리 작성
- 현재 AFM parsing 가능성 인벤토리 작성
- receptor-side residue summary format 정의
- 두 모듈용 model/cluster summary format 정의
- stable output naming 추가

## 테스트 태스크
- 최소 하나의 PyRosetta 결과가 receptor-side residue summary로 변환되는지 확인
- 최소 하나의 AFM 결과가 receptor-side residue summary로 변환되는지 확인
- resulting output이 raw file 재구성 없이 읽을 수 있는지 확인
- Vina와 PPI output naming/documentation이 개념적으로 분리되는지 확인

## 의존성
Task Group 0에 의존하며, 이후 Task Group 7 report에 연결됨.

## 예상 산출물
- PyRosetta residue-summary output format
- AFM residue-summary output format
- 두 모듈용 parsing/extraction note
- auxiliary-evidence labeling convention

---

# Task Group 7: Reporting and Manual Review Exports
**우선순위:** Should-Have  
**Zone:** 🟡 Yellow

## 목표
연구자와 협업자가 raw file을 다시 재구성하지 않고도 pocket 및 residue evidence를 읽을 수 있는 report와 manual-review 보조 파일을 생성한다.

## 메인 태스크

### 7.1 receptor-level 및 project-level summary report 생성
- receptor-local pocket 요약
- ligand-to-pocket mapping 요약
- cross-receptor pocket comparison highlight 요약

### 7.2 auxiliary PPI summary를 report에 통합
- 가능할 때 PyRosetta residue-side summary 포함
- 가능할 때 AFM residue-side summary 포함
- primary evidence와 auxiliary evidence 구분 유지

### 7.3 manual-review helper artifact export
- 필요 시 visual inspection용 representative pose export
- 수동 lookup effort를 줄이는 helper summary export

### 7.4 handoff readability 유지
- future Codex/collaborator가 full scientific rebrief 없이도 사용할 수 있게 report 구성

## 서브태스크
- 최소 report section 정의
- project-level summary table layout 정의
- ligand summary table layout 정의
- pocket comparison highlight section 정의
- optional manual-review export 규칙 정의

## 테스트 태스크
- core Vina output만으로 report 생성 가능한지 확인
- available한 경우 auxiliary PPI output이 붙는지 확인
- raw pose file을 열지 않고 report를 읽을 수 있는지 확인
- report가 raw evidence 없이 unsupported conclusion을 강요하지 않는지 확인

## 의존성
Task Groups 4, 5, 선택적으로 6에 의존.

## 예상 산출물
- project summary report
- pocket summary section
- ligand-to-pocket summary section
- cross-receptor comparison highlights
- optional manual-review export set

---

# Task Group 8: Validation, Regression Checks, and Handoff Readiness
**우선순위:** Must-Have  
**Zone:** 🟡 Yellow

## 목표
validation 단계, consistency check, regression protection, repository readiness task를 추가해, 이후 세션과 코딩 에이전트가 안전하게 작업을 이어갈 수 있게 한다.

## 메인 태스크

### 8.1 output consistency 검증
- expected run 후 core output 존재 여부 확인
- output file 간 receptor ID / ligand ID 일관성 확인
- parsed/summarized output이 raw source와 연결되어 있는지 확인

### 8.2 regression-oriented check 추가
- refactor 후 key output structure가 조용히 바뀌지 않도록 보호
- naming convention 안정성 확인
- 병렬 실행이 논리적 output content를 바꾸지 않는지 확인

### 8.3 residue numbering consistency 가정 점검
- 가능한 경우 receptor numbering mismatch 감지 또는 flag 처리
- direct comparison이 unsafe할 수 있는 경우 warning 보존

### 8.4 handoff readiness 개선
- repo에 승인된 brief, PRD, task 문서 존재 여부 확인
- future developer가 generated output 위치를 찾을 수 있는지 확인
- prior chat history 없이도 project context recovery가 가능한지 확인

## 서브태스크
- core output checklist 정의
- receptor/ligand identity consistency check 정의
- 최소 regression test target 정의
- residue numbering mismatch warning behavior 정의
- handoff-readiness checklist 정의

## 테스트 태스크
- standard run 후 required output table 존재 확인
- future session이 latest docs와 outputs를 식별할 수 있는지 확인
- 한 모듈 수정이 downstream file expectation을 조용히 깨지 않는지 확인
- numbering mismatch 상황이 가능한 경우 사용자에게 surface되는지 확인

## 의존성
모든 이전 task group에 의존.

## 예상 산출물
- output validation checklist
- regression check targets
- numbering-consistency warning behavior
- handoff-readiness checklist
- repository continuation note

---

# Codex에게 넘길 때 추천하는 초기 실행 순서

Codex에게 이 task를 점진적으로 넘긴다면, 추천하는 첫 번째 구현 배치는 다음과 같다.

1. **Task Group 0** — repository inspection and baseline conventions  
2. **Task Group 1** — structured input and run management  
3. **Task Group 2** — parallel Vina batch execution  
4. **Task Group 3** — pose parsing and contact extraction  
5. **Task Group 4** — pocket clustering and summary generation

이 첫 배치만으로도 core Vina-centered MVP evidence layer를 만들 수 있다.

두 번째 구현 배치는 다음을 추천한다.

6. **Task Group 5** — cross-receptor pocket comparison  
7. **Task Group 7** — reporting and manual review exports

세 번째 배치는 다음을 추천한다.

8. **Task Group 6** — supporting PPI output standardization  
9. **Task Group 8** — validation, regression checks, and handoff readiness

이 순서는 실제 연구 우선순위를 반영한다. 즉, 먼저 Vina 기반 pocket standardization을 만들고, 그 다음 comparison, 마지막으로 auxiliary evidence standardization을 붙이는 구조다.

