# M2.7 Changes

Added PPI/membrane/dimer hard gates under `fresh/`:

- `fresh/src/egfr_myo1d/pocket/ppi_relationship.py`
- `fresh/src/egfr_myo1d/pocket/membrane_gate.py`
- `fresh/src/egfr_myo1d/pocket/dimer_accessibility.py`
- `fresh/src/egfr_myo1d/pocket/pocket_gate_apply.py`
- CLI command `gate-m2-pockets`
- `fresh/tests/test_m2_phase7_pocket_gates.py`
- `fresh/docs/milestone2_7_ppi_membrane_dimer_pocket_gates.md`

The gate layer reads M2.4 PPI patch, M2.5 ATP reference, M2.6 pocket family
tables, M1 receptor mappings, receptor PDBs, and membrane frame evidence.

Canonical outputs are written under:

- `fresh/runs/<run_id>/phase2_pockets/gated/`

Compatibility copies are written under:

- `fresh/runs/<run_id>/phase2_pockets/tables/`

All accepted rows remain gate-result evidence only. M2.7 does not create the
M2.8 `export_for_m3/` package, does not run docking, and does not nominate
compound candidates.
