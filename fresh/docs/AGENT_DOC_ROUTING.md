# Agent Document Routing

This file tells implementation agents which documents are authoritative and how
to resolve conflicts. It exists because this repository intentionally carries
legacy code for reference while the active workflow is being rebuilt under
`fresh/`.

## First Rule

Use the fresh workflow as the active product. Legacy folders and old documents
are not controlling specifications.

## Authoritative Inputs

| Need | Read |
| --- | --- |
| Project goal and single-source scientific context | `fresh/docs/project_knowledge_final_clean_2026-04-27.md` |
| Detailed receptor, MYO1D, PPI, pocket, and compound rationale | `fresh/docs/egfr_myo1d_final_research_overview_2026-04-27.md` |
| PPI-surface tool installation and Milestone 3 integration plan | `ppi_surface_tool_installation_and_integration_handoff_v1_0.md` |
| Verified HPC environment state | `fresh/docs/hpc_tool_environment_status.md` |
| Machine-readable tool environment policy | `fresh/configs/tool_envs.yaml` |
| Tool preflight registry | `fresh/configs/tool_registry.yaml` |
| Receptor, pocket, ATP-reference, path, gate, and run settings | `fresh/configs/*.yaml` |
| Milestone 1 and Milestone 2 implementation contracts | `fresh/docs/m1_*.md`, `fresh/docs/m2_*.md`, `fresh/docs/milestone2_*.md` |
| Current behavior | `fresh/src/`, `fresh/tests/` |

## Conflict Rules

Use this order when two sources disagree:

1. Latest explicit user instruction.
2. `fresh/configs/*.yaml` for executable configuration.
3. `fresh/docs/` documents listed in this routing file.
4. Root handoff files explicitly listed in this routing file.
5. Source code and tests.
6. Legacy code as pattern reference only.

If prose and config disagree about environment names, executable paths, or tool
activation, update the prose or code to match `fresh/configs/tool_envs.yaml`
unless the user gives a newer instruction.

## Non-Authoritative Areas

Do not use the following as controlling project knowledge:

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

Permitted use of old material:

- Inspect old source code only for small implementation patterns.
- Do not import old science text, milestone definitions, environment claims, or
  output contracts.
- Do not copy old generated outputs into `fresh/`.

## Adapter Implementation Rules

Milestone 3 tool adapters must preserve explicit environment boundaries:

- Run orchestration and core tools from `pyrosetta`.
- Run pyKVFinder from `ppi_surface`.
- Run P2Rank from `p2rank_java11`.
- Keep InDeep, PeSTo, MaSIF, PocketMiner, and PASSer optional.
- Record activated environment, command line, stdout, stderr, return code, and
  output manifest path for every external tool call.
- Never silently fall back to system Java or global system tools.

## Output Rules

All fresh workflow outputs belong under:

```text
fresh/runs/<run_id>/
```

Generated data, local analysis folders, private ligand files, and external tool
clones must remain untracked.
