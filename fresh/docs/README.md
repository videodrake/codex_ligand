# Fresh Workflow Docs

This directory contains the authoritative docs for the active fresh
EGFR-MYO1D workflow. Legacy root docs, `codex_ligand*` folders, old result
exports, and analysis data folders are non-authoritative; use them only for
narrow code-pattern reference when no fresh workflow pattern exists.

Start here:

1. `../../AI_START_HERE.md`
2. `AGENT_DOC_ROUTING.md`
3. `project_knowledge_final_clean_2026-04-27.md`
4. `egfr_myo1d_final_research_overview_2026-04-27.md`
5. `hpc_tool_environment_status.md`
6. `optional_ai_tool_runtime_status.md`
7. `../configs/tool_envs.yaml`

## Current Operational Notes

- [M2 PyRosetta Input Sanitation Note](m2_pyrosetta_input_sanitation.md): documents why M2.2 rewrites AB_C inputs into Rosetta-friendly protein-only PDBs, how that compares with the legacy docking inputs, how the 3-node/32-core real PBS plan is generated, and how to debug `fill_missing_atoms` failures on HPC.
- [PPI Surface Tool Integration Handoff](ppi_surface_tool_installation_and_integration_handoff_v1_0.md): controlled install, preflight, and evidence-integration plan for fpocket/mdpocket, pyKVFinder, mini-FTMap, InDeep, PeSTo, MaSIF, PocketMiner, PASSer, and related M2/M3 surface-zone gates.
