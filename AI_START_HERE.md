# AI Start Here

This repository contains legacy code and legacy documents from earlier EGFR-MYO1D
pipeline attempts. They are not all authoritative. Every AI agent must start
here before reading other files.

## Current Authority

The active implementation target is the fresh workflow under `fresh/`.

Read these files in this order:

1. `AI_START_HERE.md`
2. `fresh/docs/AGENT_DOC_ROUTING.md`
3. `fresh/docs/project_knowledge_final_clean_2026-04-27.md`
4. `fresh/docs/egfr_myo1d_final_research_overview_2026-04-27.md`
5. `ppi_surface_tool_installation_and_integration_handoff_v1_0.md`
6. `fresh/docs/hpc_tool_environment_status.md`
7. `fresh/docs/optional_ai_tool_runtime_status.md`
8. `fresh/configs/tool_envs.yaml`
9. The relevant `fresh/configs/*.yaml`, `fresh/docs/m*_*.md`, source files, and
   tests for the task you are changing.

If two documents conflict, use this precedence:

1. The latest user instruction in the current conversation.
2. Machine-readable `fresh/configs/*.yaml` for runtime settings and gates.
3. Fresh workflow docs under `fresh/docs/`.
4. Root handoff docs explicitly listed above.
5. Source code and tests.
6. Legacy materials only as non-authoritative implementation examples.

## Do Not Use As Controlling Specs

Do not use these as scientific, workflow, environment, or milestone authority:

- `codex_ligand*` folders
- root `docs/`
- root `config/`
- root `egfr_pipeline/`
- root `scripts/`
- root `tests/`
- `results_export/` and `results_export.tar.gz`
- `analysis_v0_5/`
- `analysis_v0_4_folder1_4/`
- analysis zip files
- `node_modules/`

Old code may be inspected only for reusable implementation patterns, and only
when the active fresh workflow has no local pattern. Old documents must not be
used to decide project scope, scientific claims, tool requirements, or output
contracts.

## Active HPC Reality

The verified HPC layout is split by environment:

- `pyrosetta`: core workflow orchestration, PyRosetta, RDKit, BioPython, numpy,
  pandas, Vina, fpocket, mdpocket, Open Babel, PBS tools, and GROMACS.
- `ppi_surface`: isolated pyKVFinder environment.
- `p2rank_java11`: isolated Java 11 and P2Rank environment.

Do not modify system Java. P2Rank needs Java 11 and must be run from the
isolated `p2rank_java11` env.

Use `fresh/docs/hpc_tool_environment_status.md` for human-readable verification
notes, `fresh/docs/optional_ai_tool_runtime_status.md` for optional AI tool
runtime details, and `fresh/configs/tool_envs.yaml` for machine-readable adapter
policy.

## Output and Git Rules

All generated fresh workflow outputs must stay under:

```text
fresh/runs/<run_id>/
```

Do not commit large local analysis folders, generated run outputs, private
ligands, external tool clones, or downloaded data bundles.
