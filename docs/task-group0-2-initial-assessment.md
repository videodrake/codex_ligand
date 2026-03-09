# EGFR–MYO1D Pipeline: Task Group 0~2 Initial Assessment

## 1) Repository map (current state)

- Root
  - `.gitkeep`
  - `docs/`
    - `task-group0-2-initial-assessment.md` (this file)

## 2) Vina entrypoint/config/output structure scan

### Requested files not found
The following requested documents were not present in the repository checkout, so no spec-driven implementation analysis was possible:

1. `README.md`
2. `docs/project-context.md`
3. `docs/brief-egfr-myo1d-pipeline.md`
4. `docs/prd-egfr-myo1d-pipeline.md`
5. `docs/tasks-egfr-myo1d-pipeline.md`
6. `CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md`
7. `docs/runbook.md`

### Code scan result
- No pipeline source code, scripts, configs, or Vina-related files were found in the current branch.
- No existing entrypoints were found (e.g., `main.py`, `cli.py`, `run_*.sh`, `Snakefile`, `nextflow.config`, etc.).
- No Vina config templates or output tree conventions were found.

## 3) Minimal refactoring plan (blocked until source/docs are available)

When the intended codebase and docs are restored, execute the following low-risk sequence:

1. **Inventory current execution paths**
   - Locate all Vina invocation points and wrappers.
   - Build a call graph from top-level CLI/script to Vina subprocess call.

2. **Normalize receptor set as explicit inputs (exactly 3)**
   - Define one typed receptor registry with:
     - `3GT8 raw`
     - `MD cluster representative 38–48`
     - `MD cluster representative 85–100`
   - Avoid hardcoding legacy residue/site labels from reports.

3. **Extract runtime config layer**
   - Add/standardize a config model for:
     - receptor list
     - ligand set
     - docking box parameters
     - CPU allocation default `16` (despite 32-core host)
     - output root and run-id

4. **Standardize output schema**
   - Ensure per-run deterministic layout (e.g., per receptor / per ligand / scores).
   - Capture provenance metadata so newer computation can supersede older evidence.

5. **Introduce compatibility wrapper**
   - Preserve existing CLI/entrypoint behavior while routing through new config layer.

## 4) File-by-file change plan (initial target list, pending real files)

Because source files are missing, this is a placeholder plan by likely concern:

- `README.md`
  - Add canonical receptor triad and 16-core operational default.

- `docs/runbook.md`
  - Document execution profile for 16-core parallelism and result precedence policy.

- `pipeline entrypoint` (actual path TBD)
  - Refactor to use centralized config object and receptor registry.

- `vina config builder` (actual path TBD)
  - Ensure receptor-specific config generation from runtime inputs.

- `output writer` (actual path TBD)
  - Enforce deterministic output tree + evidence timestamp/provenance.

## 5) Assumptions requiring confirmation

1. The current branch appears to be an empty scaffold and may not be the intended working tree.
2. The 7 requested docs may exist in another branch or remote but are absent here.
3. The existing "GitHub codebase" to be reused is not present in local checkout.
4. Expected tech stack (Python/Snakemake/Nextflow/Bash) is unknown until files are available.
5. Any receptor/site naming from prior reports should be treated as non-authoritative unless recomputed.

## 6) Smallest first patch proposal

**Proposed first low-risk code patch once files are available:**

- Add a single source of truth for receptor identifiers and CPU default in a small config module.
- Wire only the top-level CLI to consume these values (no scoring/output algorithm changes).
- Keep old flags backward-compatible and map legacy receptor aliases to the new canonical triad.

This limits blast radius while satisfying key constraints:
- exactly 3 receptors,
- 16-core execution baseline,
- avoidance of hardcoded legacy residue/site assumptions,
- enabling fresh calculations as higher-priority evidence.
