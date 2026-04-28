# Milestone 1 Foundation Plan

Milestone 1 creates the clean `fresh/` workflow foundation for the EGFR-MYO1D PPI-to-pocket-to-compound project.

Scope for Task 1:

- Create the `fresh/` directory skeleton.
- Add configuration placeholders.
- Add documentation placeholders and do-not-touch guardrails.
- Add small synthetic fixtures for future parser and normalization tests.
- Add a minimal importable package and CLI.
- Add `.gitignore` protections for fresh runtime outputs and confidential ligand/private data.

Out of scope for Task 1:

- PyRosetta production docking.
- AutoDock Vina docking.
- fpocket or P2Rank production pocket discovery.
- PPI pose scoring.
- Compound docking or final scoring.
- Candidate nomination.
- Real receptor normalization, membrane-frame computation, or MYO1D slicing.

All fresh run outputs must be written under `fresh/runs/<run_id>/` in later tasks.
