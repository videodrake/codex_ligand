# Test Suite Triage

This repository now uses two practical test lanes instead of treating every test as equally blocking.

## Keep as blocking pre-qsub tests

These tests directly protect code paths that would waste compute time if broken:

- phase 2 core logic and schema tests
- phase 3 core logic and schema tests
- phase 4 core scoring and schema tests
- CLI smoke tests
- validation smoke tests

These remain in the pre-qsub lane because they guard:

- schema contracts
- data handoff between phases
- deterministic job construction
- scoring behavior
- entrypoint safety

## Keep, but move out of pre-qsub blocking lane

The following tests are still useful, but they are not the best gate before spending cluster resources:

- report-generation tests
- review-table formatting tests
- presentation-summary tests

They are marked with `@pytest.mark.reporting`.

Why they stay in the repository:

- they still protect user-facing artifacts
- they still catch drift in summaries and presentation outputs

Why they are excluded from the pre-qsub lane:

- they do not usually predict whether a heavy compute submission will fail
- they are lower priority than schema, handoff, and execution safety

## Current decision

No existing test was deleted in this cleanup.

Reason:

- the large synthetic phase tests are still doing valuable contract work
- the lower-value group is small enough to keep
- classification is safer than deletion at this stage

## Commands

Pre-qsub lane:

```bash
bash scripts/run_pre_qsub_checks.sh
```

Full regression lane:

```bash
python -m pytest tests -q
```

Reporting-only lane:

```bash
python -m pytest -m reporting tests -q
```
