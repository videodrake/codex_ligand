# M2.5 Change Notes

Added ATP reference generation under `fresh/`:

- `fresh/configs/atp_reference.yaml`
- `fresh/src/egfr_myo1d/pocket/atp_reference.py`
- CLI command `build-m2-atp-reference`
- `fresh/tests/test_m2_phase5_atp_reference.py`
- `fresh/docs/milestone2_5_atp_reference.md`

The builder writes canonical ATP reference outputs under:

- `fresh/runs/<run_id>/phase2_pockets/atp_reference/`

It also writes a compatibility copy:

- `fresh/runs/<run_id>/phase2_pockets/tables/atp_site_reference.csv`

The implementation supports ligand-bearing reference PDBs, configured residue
sets from `fresh/configs/atp_reference.yaml`, and synthetic test fixtures. It
preserves chain/protomer identity and reports missing mappings or residue-name
mismatches as warnings.

Review hardening:

- M2.5 now fails when no primary membrane-validated state ATP centroid can be
  built.
- Ligand-neighbor residues are restored through M1 receptor mapping before they
  are written as UniProt/source residue numbers.
- The configured fallback ATP-cleft set was expanded and now declares a minimum
  residue count so sparse fallback configs cannot silently advance as a hard
  gate.

No fpocket/P2Rank, PyRosetta docking/relaxation, Vina, LightDock, qsub/sbatch,
cleanup deletion, compound scoring, pocket discovery, or candidate nomination is
performed.
