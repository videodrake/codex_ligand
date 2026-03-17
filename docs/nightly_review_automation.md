# Nightly Review Automation

Last updated: 2026-03-16

This document explains the repository-side tooling that supports the nightly deep-review automation.

For the patch-producing follow-up loop that turns review output into a single daily micro-change, see [nightly_incremental_improvement_automation.md](nightly_incremental_improvement_automation.md).

## Purpose

`scripts/nightly_review.py` builds a review bundle that is meant to be consumed by an automated Codex review run.

Each bundle captures:

- project config context
- git change summary
- granular review units across runtime, Vina, PPI, Phase 1-4, config, inputs, outputs, tests, and docs
- related test and documentation pointers for each unit
- a checklist, prompt packet, and report template

## Manual Usage

Run from the repository root:

```bash
python scripts/nightly_review.py
```

Optional flags:

```bash
python scripts/nightly_review.py --label nightly-manual
python scripts/nightly_review.py --output-root output/review_automation
python scripts/nightly_review.py -c config/example-project.yaml
```

## Generated Files

Each run creates `output/review_automation/<label-or-timestamp>/` with:

- `review_manifest.json`: machine-readable inventory and git/context snapshot
- `review_checklist.md`: review order plus per-unit scope, tests, docs, and focus points
- `review_prompt.md`: reviewer instructions for the automation run
- `nightly_review_report.md`: report target that the automation should overwrite with findings

The latest bundle pointer is also written to:

- `output/review_automation/latest.json`

## Review Scope

The bundle covers these major surfaces:

- entry points and orchestration
- runtime core utilities
- workflow scripts and derived step views
- routine Vina execution and postprocess layers
- verdict, report, and validation integration
- PPI support and PyRosetta core
- MD analysis
- scientific Phase 1 through Phase 4 modules
- config, input, output, test, and documentation contracts

## Automation Expectations

The downstream automation should:

1. generate the latest bundle with `python scripts/nightly_review.py`
2. read `review_manifest.json` and `review_checklist.md`
3. review changed files first, then all critical units
4. run the smallest relevant pytest targets when executable code is involved
5. overwrite `nightly_review_report.md` with severity-ordered findings
6. state residual risks or test gaps even when no findings are discovered

## Scheduling Note

The Codex automation for this tool is intended to run hourly during the night review window in the user's local time zone.
