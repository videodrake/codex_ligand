# M2.3 Change Notes

Added an auditable M2.3 scaffold under `fresh/`:

- `fresh/src/egfr_myo1d/ppi/collect_ppi_outputs.py`
- `fresh/src/egfr_myo1d/ppi/restore_residue_mapping.py`
- CLI command `collect-m2-ppi-outputs`
- `fresh/tests/test_m2_phase3_ppi_collection_mapping.py`
- `fresh/docs/m2_phase3_ppi_collection_mapping.md`

The implementation collects M2.2 dry-run status rows and restores optional
run-local raw contacts through M1 receptor mapping. It writes only under
`fresh/runs/<run_id>/`.

Scientific execution remains out of scope. No docking, relaxation, LightDock,
Vina, fpocket/P2Rank, compound scoring, candidate nomination, scheduler
submission, or cleanup deletion is performed.
