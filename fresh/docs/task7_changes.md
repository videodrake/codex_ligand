# Task 7 Changes

Created:

- `fresh/src/egfr_myo1d/analysis/__init__.py`
- `fresh/src/egfr_myo1d/analysis/ppi_consensus.py`
- `fresh/src/egfr_myo1d/validation/ppi_consensus.py`
- `fresh/tests/fixtures/task7_ppi_consensus/accepted_ppi_contacts.csv`
- `fresh/tests/test_task7_ppi_consensus_patch.py`
- `fresh/docs/task7_ppi_consensus_patch.md`
- `fresh/docs/task7_changes.md`

Modified:

- `fresh/src/egfr_myo1d/cli.py`

Summary:

- Added `summarize-ppi-consensus`.
- Added schema validation for supplied PPI contact tables.
- Added EGFR and MYO1D residue-list parsing audits.
- Added active-face, sheet-12, tail/noise, ATP-overlap, membrane-proximal, and
  convergence audits.
- Added deterministic EGFR-side consensus patch CSV generation.
- Preserved Task 6 spec-only posture and no-score-bonus policy.

Validation:

- `python -m py_compile fresh/src/egfr_myo1d/analysis/ppi_consensus.py fresh/src/egfr_myo1d/validation/ppi_consensus.py fresh/tests/test_task7_ppi_consensus_patch.py`
- `pytest -q fresh/tests/test_task7_ppi_consensus_patch.py`
- `pytest -q fresh/tests`

