# Fresh EGFR-MYO1D Workflow

This is the fresh dimer/membrane-aware EGFR-MYO1D workflow area.

Milestone 1 is foundation only. No docking, pocket discovery, scoring, final ranking, or candidate nomination is implemented in Task 1.

All run outputs must go under:

```text
fresh/runs/<run_id>/
```

Do not write fresh workflow outputs into the old `output/` or `results_export/` areas.

Ligand structures are confidential by default and ignored by git. Public-facing files should use only:

```text
Cpd-A
Cpd-B
Cpd-C
```

Package commands are run from the repository root with `PYTHONPATH`:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.cli --help
python -m egfr_myo1d.cli version
```

Future milestones will add runtime context, manifests, logging, receptor normalization, membrane-frame handling, MYO1D preparation, and smoke harnesses. Task 1 intentionally creates only the clean skeleton and guardrails.
