# M2.8 Changes

Implemented M2.8 accepted pocket export and Milestone 2 reporting under
`fresh/`.

## Added

- `fresh/src/egfr_myo1d/pocket/m2_export.py`
- CLI command:

```bash
python -m egfr_myo1d.cli export-m2-results --run-id <RUN_ID>
```

- Focused tests:
  - `fresh/tests/test_m2_phase8_export_report.py`
- Documentation:
  - `fresh/docs/milestone2_8_accepted_pocket_export_and_report.md`

## Outputs

M2.8 writes final aggregation outputs under:

```text
fresh/runs/<run_id>/phase2_pockets/final/
```

and the M3 handoff package under:

```text
fresh/runs/<run_id>/phase2_pockets/export_for_m3/
```

It also writes a compatibility report:

```text
fresh/runs/<run_id>/reports/milestone2_summary.md
```

## Guardrails

- ATP rejects are not exported.
- Hard-gate failures are not rescued by soft score.
- 3GT8_raw-only pockets are kept as `reference_only`.
- Cleanup is report-only and deletes nothing.
- No Vina, ligand prep, receptor PDBQT prep, PyRosetta docking/relaxation,
  LightDock, P2Rank, qsub/sbatch production, compound scoring, or candidate
  nomination is run.
