# M2.4 Change Notes

Added the M2.4 consensus patch layer under `fresh/`:

- `fresh/src/egfr_myo1d/ppi/consensus_patch.py`
- `fresh/src/egfr_myo1d/build_ppi_consensus_patch.py`
- CLI command `build-m2-ppi-consensus-patch`
- `fresh/tests/test_m2_phase4_ppi_consensus_patch.py`
- `fresh/docs/m2_phase4_ppi_consensus_patch.md`

The builder consumes M2.3 restored contact evidence and writes:

- `fresh/runs/<run_id>/output/ppi/ppi_consensus_patch.csv`
- `fresh/runs/<run_id>/phase1_ppi/tables/ppi_consensus_patch.csv`
- support, QC, manifest, and report files under the same run directory

The implementation re-checks M1 receptor mapping before consensus rows are
emitted, keeps A/B chain and protomer identity explicit, and reports residue
drift as a failure.

Review hardening:

- Pending or blank pose-QC contact evidence is no longer emitted as consensus.
- `support_fraction` uses the M2.2 planned seed/model denominator instead of
  only contact-bearing poses.
- Existing chain/protomer IDs are compared against the M1 mapping and drift is
  reported as a failure.

Scientific execution remains out of scope: no docking, relaxation, pocket
discovery, scheduler submission, compound scoring, cleanup deletion, or
candidate nomination is performed.
