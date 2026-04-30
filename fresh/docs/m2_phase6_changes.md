# M2.6 Changes

Added fpocket parser/adapter and raw pocket normalization under `fresh/`:

- `fresh/configs/pocket.yaml`
- `fresh/src/egfr_myo1d/pocket/fpocket_adapter.py`
- `fresh/src/egfr_myo1d/pocket/normalization.py`
- `fresh/src/egfr_myo1d/pocket/family_merge.py`
- CLI command `run-m2-fpocket-discovery`
- `fresh/tests/test_m2_phase6_fpocket_discovery.py`
- `fresh/docs/milestone2_6_fpocket_adapter_and_pocket_normalization.md`

The command writes canonical raw and merged M2.6 outputs under:

- `fresh/runs/<run_id>/phase2_pockets/raw/`
- `fresh/runs/<run_id>/phase2_pockets/merged/`

It supports parser-only/synthetic fixture parsing without fpocket installed and
production fpocket execution only on run-local staged receptor PDBs.

M2.6 leaves all candidates as `raw_ungated`. It does not apply ATP/PPI/membrane
/dimer gates and does not create accepted M3 pocket exports.
