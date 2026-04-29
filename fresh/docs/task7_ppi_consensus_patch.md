# Task 7 PPI Consensus Patch

Task 7 adds a pure-Python post-processing layer for supplied EGFR-MYO1D PPI
contact records. It reads accepted or future accepted-pose contact tables,
validates their schema, preserves EGFR protomer/chain/residue evidence, tracks
MYO1D annotation classes, and writes guarded EGFR-side consensus patch evidence.

Task 7 is post-processing/specification only. It does not run docking,
PyRosetta, LightDock, Vina, fpocket, P2Rank, AlphaFold, Boltz, Chai, qsub, PBS,
sbatch, pocket discovery, ligand docking, scoring, or candidate nomination.

## CLI

Bash:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id test_task7_local
python -m egfr_myo1d.cli summarize-ppi-consensus \
  --run-id test_task7_local \
  --mode smoke_env \
  --profile codex_dev \
  --input-root fresh/tests/fixtures/task7_ppi_consensus
```

PowerShell:

```powershell
$env:PYTHONPATH = "$PWD\fresh\src;$env:PYTHONPATH"
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id test_task7_local
python -m egfr_myo1d.cli summarize-ppi-consensus `
  --run-id test_task7_local `
  --mode smoke_env `
  --profile codex_dev `
  --input-root fresh/tests/fixtures/task7_ppi_consensus
```

## Input

The default input table is `accepted_ppi_contacts.csv` under `--input-root`.
Required columns include receptor metadata, method/seed/replicate, pose ID,
EGFR contact residues, EGFR contact centroid, MYO1D contact residues, active-face
contacts, sheet-12 support contacts, tail/noise contacts, ATP and membrane flags,
terminal-artifact flags, and pose acceptance class.

EGFR residues must preserve chain identity, for example `A:742;A:746;A:923`.
MYO1D residues may be listed as residue numbers, for example `961;962;968;993`.
Malformed residue IDs are reported; Task 7 does not silently repair them.

## Outputs

All outputs are written under `fresh/runs/<run_id>/`:

- `manifest/ppi_consensus_input_manifest.json`
- `manifest/ppi_consensus_qc_report.json`
- `qc/ppi_consensus_schema_audit.csv`
- `qc/ppi_residue_parsing_audit.csv`
- `qc/ppi_tail_noise_audit.csv`
- `qc/ppi_active_face_audit.csv`
- `qc/ppi_membrane_atp_audit.csv`
- `qc/ppi_convergence_audit.csv`
- `output/ppi/ppi_consensus_patch.csv`
- `reports/ppi_consensus_report.txt`

`ppi_consensus_patch.csv` is a consensus evidence table. It is not a compound
candidate table and does not nominate ligands or pockets.

## Guardrails

MYO1D active-face residues `961-964` and `968-972`, sheet-12 support residues
`993-997`, and tail/noise residues `998-1006` are QC annotations only. They do
not create score bonuses. Tail/noise and terminal-artifact evidence remains
visible and can classify a patch as `TAIL_DOMINATED`.

ATP-overlap and membrane-proximal flags are propagated to QC reports and the
consensus patch CSV. Dimer/protomer context remains explicit through
`protomer_id` and chain-preserved EGFR residue IDs.

Evidence classes are cautious computational labels:

- `CONVERGENT_PATCH`
- `BROAD_PATCH`
- `DISPERSED`
- `TAIL_DOMINATED`
- `ATP_OVERLAPPING`
- `MEMBRANE_PROXIMAL`
- `INSUFFICIENT_EVIDENCE`
- `SCHEMA_INVALID`

These classes prepare the handoff to future pocket discovery and require visual
and scientific review before interpretation.

