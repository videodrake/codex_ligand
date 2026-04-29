"""M1 Phase 8: prepare-inputs orchestrator (handoff §10.2 Smoke B + §19 + §22).

Sequences the M1 input-preparation sub-steps under one ``RunContext``:

    preflight
        -> prepare-receptor (per state from receptor_states.yaml)
        -> compute-membrane-frame
        -> prepare-myo1d
        -> manifest-ligands

Behavior:
- ``codex_dev`` profile: a sub-step FAIL is recorded and the orchestrator
  continues to the next sub-step. Aggregate ``status`` is FAIL if any
  sub-step FAILed.
- ``hpc_strict`` profile (or ``strict=True``): a sub-step FAIL stops the
  orchestrator immediately.
- Missing real input files in ``smoke_input`` are reported under
  ``missing_required_inputs`` without crashing.
- Aggregate manifest (``prepare_inputs_aggregate_manifest.json``) and
  Markdown summary report are emitted under the run directory.
- Ligand step can be skipped with ``skip_ligands=True``.

Closes M1 §23 #15.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from egfr_myo1d.core.logging_utils import append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError


@dataclass
class SubStepResult:
    name: str
    status: str  # PASS / WARN / FAIL / SKIPPED
    manifest_path: str | None = None
    warnings: list[str] = field(default_factory=list)
    fail_reason: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class PrepareInputsAggregate:
    run_id: str
    mode: str
    profile: str
    input_root: str
    states_processed: list[str]
    sub_steps: list[SubStepResult] = field(default_factory=list)
    missing_required_inputs: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    status: str = "PASS"
    timestamp: str = ""
    aggregate_manifest_path: Path | None = None
    summary_report_path: Path | None = None


# ---------------------------------------------------------------------------
# Path resolution helpers
# ---------------------------------------------------------------------------


def _resolve_state_paths(
    ctx: RunContext, input_root: Path | None
) -> tuple[dict[str, Path], Path | None]:
    """Resolve per-state PDB paths and the plus10_full_frame fallback.

    Order of resolution:
    1. If ``input_root/receptors/<state>.pdb`` exists, use it.
    2. Else if ``receptor_states.yaml.input_files[<state>]`` exists, use it
       relative to repo root.
    3. Else None (the receptor sub-step will then mark the state missing).
    """
    config = load_yaml_config(
        ctx.repo_root / "fresh" / "configs" / "receptor_states.yaml"
    )
    input_files = config.get("input_files", {})
    state_paths: dict[str, Path] = {}
    states = (
        config.get("primary_membrane_validated_states", [])
        + config.get("reference_control_states", [])
    )
    for state_id in states:
        if input_root is not None:
            override = input_root / "receptors" / "{0}.pdb".format(state_id)
            if override.is_file():
                state_paths[state_id] = override
                continue
        rel = input_files.get(state_id)
        if rel:
            candidate = ctx.repo_root / rel
            state_paths[state_id] = candidate if candidate.is_file() else candidate
    plus10 = None
    if input_root is not None:
        cand = input_root / "receptors" / "plus10_full_frame.pdb"
        if cand.is_file():
            plus10 = cand
    if plus10 is None:
        rel_plus10 = input_files.get("plus10_full_frame")
        if rel_plus10:
            cand = ctx.repo_root / rel_plus10
            plus10 = cand if cand.is_file() else None
    return state_paths, plus10


def _resolve_myo1d_source(ctx: RunContext, input_root: Path | None) -> Path:
    if input_root is not None:
        candidate = input_root / "myo1d" / "AF-O94832-F1-model_v6.pdb"
        if candidate.is_file():
            return candidate
    return ctx.repo_root / "fresh" / "data" / "raw" / "myo1d" / "AF-O94832-F1-model_v6.pdb"


def _resolve_ligands_dir(ctx: RunContext, input_root: Path | None) -> Path:
    if input_root is not None:
        candidate = input_root / "ligands"
        if candidate.is_dir():
            return candidate
    paths_cfg = load_yaml_config(ctx.repo_root / "fresh" / "configs" / "paths.yaml")
    return ctx.repo_root / paths_cfg.get("raw_ligands", "fresh/data/raw/ligands")


def _resolve_private_mapping(ctx: RunContext, input_root: Path | None) -> Path:
    if input_root is not None:
        candidate = input_root / "private" / "compound_id_map.csv"
        if candidate.is_file():
            return candidate
    paths_cfg = load_yaml_config(ctx.repo_root / "fresh" / "configs" / "paths.yaml")
    return (
        ctx.repo_root
        / paths_cfg.get("private_data", "fresh/data/private")
        / "compound_id_map.csv"
    )


def _resolve_state_list(ctx: RunContext, states: Iterable[str] | None) -> list[str]:
    if states is not None:
        return list(states)
    config = load_yaml_config(
        ctx.repo_root / "fresh" / "configs" / "receptor_states.yaml"
    )
    return list(config.get("primary_membrane_validated_states", [])) + list(
        config.get("reference_control_states", [])
    )


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------


def _write_aggregate_manifest(
    ctx: RunContext, aggregate: PrepareInputsAggregate
) -> Path:
    payload: dict[str, Any] = {
        "run_id": aggregate.run_id,
        "mode": aggregate.mode,
        "profile": aggregate.profile,
        "input_root": aggregate.input_root,
        "states_processed": aggregate.states_processed,
        "sub_steps": [
            {
                "name": step.name,
                "status": step.status,
                "manifest_path": step.manifest_path,
                "warnings": step.warnings,
                "fail_reason": step.fail_reason,
                "extra": step.extra,
            }
            for step in aggregate.sub_steps
        ],
        "missing_required_inputs": aggregate.missing_required_inputs,
        "warnings": aggregate.warnings,
        "blockers": aggregate.blockers,
        "status": aggregate.status,
        "timestamp": aggregate.timestamp,
        "score_bonus_allowed": False,
    }
    target = ctx.manifest_dir / "prepare_inputs_aggregate_manifest.json"
    write_json(target, payload, ctx)
    return target


def _write_summary_report(ctx: RunContext, aggregate: PrepareInputsAggregate) -> Path:
    target = ctx.reports_dir / "prepare_inputs_summary.md"
    safe = ctx.require_within_run_dir(target)
    safe.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("# prepare-inputs summary")
    lines.append("")
    lines.append("- run_id: `{0}`".format(aggregate.run_id))
    lines.append("- mode: `{0}`".format(aggregate.mode))
    lines.append("- profile: `{0}`".format(aggregate.profile))
    lines.append("- input_root: `{0}`".format(aggregate.input_root))
    lines.append("- status: **{0}**".format(aggregate.status))
    lines.append("- timestamp: `{0}`".format(aggregate.timestamp))
    lines.append("")
    lines.append("## Sub-steps")
    lines.append("")
    lines.append("| Step | Status | Warnings | Manifest |")
    lines.append("| --- | --- | --- | --- |")
    for step in aggregate.sub_steps:
        warnings_text = ";".join(step.warnings) if step.warnings else ""
        manifest_text = step.manifest_path or ""
        lines.append(
            "| {0} | {1} | {2} | {3} |".format(
                step.name, step.status, warnings_text, manifest_text
            )
        )
    lines.append("")
    if aggregate.missing_required_inputs:
        lines.append("## Missing required inputs")
        lines.append("")
        for item in aggregate.missing_required_inputs:
            lines.append("- `{0}`".format(item))
        lines.append("")
    if aggregate.blockers:
        lines.append("## Blockers")
        lines.append("")
        for item in aggregate.blockers:
            lines.append("- {0}".format(item))
        lines.append("")
    if aggregate.warnings:
        lines.append("## Warnings")
        lines.append("")
        for item in aggregate.warnings:
            lines.append("- {0}".format(item))
        lines.append("")
    safe.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return safe


# ---------------------------------------------------------------------------
# Sub-step runners (each catches exceptions and returns a SubStepResult)
# ---------------------------------------------------------------------------


def _run_preflight_step(
    ctx: RunContext, mode: str, profile: str
) -> SubStepResult:
    try:
        from egfr_myo1d.validation.preflight import run_preflight

        # The orchestrator runs preflight as an env-only check. Per-input
        # presence is verified by the downstream receptor / myo1d / ligand
        # sub-steps (which honor --input-root). Calling preflight in
        # smoke_input mode here would FAIL on missing default raw/ inputs even
        # when the orchestrator has been pointed at a fixture --input-root.
        report = run_preflight(ctx, "smoke_env", profile)
        return SubStepResult(
            name="preflight",
            status=report.get("status", "PASS"),
            manifest_path=str(ctx.manifest_dir / "environment_report.json"),
            warnings=[],
            extra={
                "preflight_mode": "smoke_env",
                "orchestration_mode": mode,
                "counts": report.get("counts", {}),
            },
        )
    except Exception as exc:
        return SubStepResult(
            name="preflight",
            status="FAIL",
            manifest_path=None,
            warnings=[],
            fail_reason="preflight_exception: {0}".format(exc),
        )


def _run_receptor_step(
    ctx: RunContext,
    state_id: str,
    source: Path | None,
    profile: str,
    strict: bool,
) -> tuple[SubStepResult, list[str]]:
    missing: list[str] = []
    if source is None or not source.is_file():
        missing.append("receptor:{0}".format(state_id))
        return (
            SubStepResult(
                name="prepare-receptor:{0}".format(state_id),
                status="FAIL" if profile == "hpc_strict" or strict else "WARN",
                manifest_path=None,
                warnings=["receptor_source_missing"],
                fail_reason=None,
                extra={"state_id": state_id, "source": str(source) if source else None},
            ),
            missing,
        )
    try:
        from egfr_myo1d.model.receptor_normalize import normalize_receptor

        receptor = normalize_receptor(
            ctx, source, state_id, profile=profile, strict=strict
        )
        return (
            SubStepResult(
                name="prepare-receptor:{0}".format(state_id),
                status=receptor.status,
                manifest_path=str(receptor.manifest_json) if receptor.manifest_json else None,
                warnings=list(receptor.warnings),
                extra={
                    "state_id": state_id,
                    "case": receptor.case,
                    "v924r_warn": receptor.v924r_warn,
                    "n_residues_a": receptor.n_residues_a,
                    "n_residues_b": receptor.n_residues_b,
                },
            ),
            missing,
        )
    except Exception as exc:
        return (
            SubStepResult(
                name="prepare-receptor:{0}".format(state_id),
                status="FAIL",
                manifest_path=None,
                warnings=[],
                fail_reason="prepare_receptor_exception: {0}".format(exc),
            ),
            missing,
        )


def _run_membrane_frame_step(
    ctx: RunContext,
    state_paths: dict[str, Path],
    plus10: Path | None,
    state_ids: Iterable[str],
    profile: str,
) -> SubStepResult:
    try:
        from egfr_myo1d.model.membrane_frame import (
            run_membrane_frame_computation,
        )

        # Choose a non-None source as the unified --full-frame-source.
        # Prefer plus10, else any state path that exists.
        source: Path | None = plus10
        if source is None:
            for sid in state_ids:
                p = state_paths.get(sid)
                if p is not None and p.is_file():
                    source = p
                    break
        frames, overall = run_membrane_frame_computation(
            ctx,
            state_ids=list(state_ids),
            full_frame_source=source,
            profile=profile,
        )
        # Map PASS_WITH_WARNINGS to WARN for the sub-step status field
        mapped = "WARN" if overall == "PASS_WITH_WARNINGS" else overall
        return SubStepResult(
            name="compute-membrane-frame",
            status=mapped,
            manifest_path=str(ctx.manifest_dir / "membrane_frame.json"),
            warnings=[],
            extra={
                "states_in_frame": [f.state_id for f in frames],
                "overall_native": overall,
            },
        )
    except Exception as exc:
        return SubStepResult(
            name="compute-membrane-frame",
            status="FAIL",
            manifest_path=None,
            warnings=[],
            fail_reason="compute_membrane_frame_exception: {0}".format(exc),
        )


def _run_myo1d_step(
    ctx: RunContext, source: Path, profile: str
) -> tuple[SubStepResult, list[str]]:
    missing: list[str] = []
    if not source.is_file():
        missing.append("myo1d:{0}".format(source))
    try:
        from egfr_myo1d.myo1d.qc import run_myo1d_qc

        report = run_myo1d_qc(ctx, source, construct_range=None, profile=profile)
        return (
            SubStepResult(
                name="prepare-myo1d",
                status=report.status,
                manifest_path=str(report.output_manifest) if report.output_manifest else None,
                warnings=list(report.warnings),
                extra={
                    "construct_id": report.construct_id,
                    "n_residues": report.n_residues,
                    "key_sheet8_present": report.key_sheet8_present,
                    "key_sheet9_present": report.key_sheet9_present,
                    "key_sheet12_present": report.key_sheet12_present,
                },
            ),
            missing,
        )
    except Exception as exc:
        return (
            SubStepResult(
                name="prepare-myo1d",
                status="FAIL",
                manifest_path=None,
                warnings=[],
                fail_reason="prepare_myo1d_exception: {0}".format(exc),
            ),
            missing,
        )


def _run_ligand_step(
    ctx: RunContext,
    ligands_dir: Path,
    private_mapping: Path,
    profile: str,
    compound_stage_enabled: bool,
) -> SubStepResult:
    try:
        from egfr_myo1d.ligand.manifest import build_ligand_manifest

        manifest = build_ligand_manifest(
            ctx,
            ligands_dir=ligands_dir,
            private_mapping_path=private_mapping,
            profile=profile,
            compound_stage_enabled=compound_stage_enabled,
        )
        return SubStepResult(
            name="manifest-ligands",
            status=manifest.status,
            manifest_path=str(manifest.output_manifest_json)
            if manifest.output_manifest_json
            else None,
            warnings=list(manifest.warnings),
            extra={
                "present_count": manifest.present_count,
                "missing_count": manifest.missing_count,
                "internal_ids_leaked_into_outputs": manifest.internal_ids_leaked_into_outputs,
            },
        )
    except Exception as exc:
        return SubStepResult(
            name="manifest-ligands",
            status="FAIL",
            manifest_path=None,
            warnings=[],
            fail_reason="manifest_ligands_exception: {0}".format(exc),
        )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def run_prepare_inputs(
    ctx: RunContext,
    mode: str = "smoke_input",
    profile: str = "codex_dev",
    input_root: Path | None = None,
    states: Iterable[str] | None = None,
    skip_ligands: bool = False,
    strict: bool = False,
    compound_stage_enabled: bool = False,
) -> PrepareInputsAggregate:
    """Run the M1 input-preparation orchestration."""
    if profile not in {"codex_dev", "hpc_strict"}:
        raise ValueError("profile must be codex_dev or hpc_strict; got: {0}".format(profile))
    if mode not in {"smoke_env", "smoke_input"}:
        raise ValueError("mode must be smoke_env or smoke_input; got: {0}".format(mode))

    state_ids = _resolve_state_list(ctx, states)
    state_paths, plus10_path = _resolve_state_paths(ctx, input_root)
    myo1d_source = _resolve_myo1d_source(ctx, input_root)
    ligands_dir = _resolve_ligands_dir(ctx, input_root)
    private_mapping = _resolve_private_mapping(ctx, input_root)

    aggregate = PrepareInputsAggregate(
        run_id=ctx.run_id,
        mode=mode,
        profile=profile,
        input_root=str(input_root) if input_root else "",
        states_processed=list(state_ids),
        timestamp=now_iso(),
    )

    def _fold(step: SubStepResult, missing: list[str] | None = None) -> bool:
        aggregate.sub_steps.append(step)
        if missing:
            aggregate.missing_required_inputs.extend(missing)
        if step.status == "FAIL":
            aggregate.blockers.append("{0}:{1}".format(step.name, step.fail_reason or "FAIL"))
            if profile == "hpc_strict" or strict:
                return True  # signal stop
        return False

    # 1. preflight
    if _fold(_run_preflight_step(ctx, mode, profile)):
        return _finalize(ctx, aggregate)

    # 2. prepare-receptor per state
    for state_id in state_ids:
        source = state_paths.get(state_id)
        step, missing = _run_receptor_step(ctx, state_id, source, profile, strict)
        if _fold(step, missing):
            return _finalize(ctx, aggregate)

    # 3. compute-membrane-frame
    if _fold(
        _run_membrane_frame_step(ctx, state_paths, plus10_path, state_ids, profile)
    ):
        return _finalize(ctx, aggregate)

    # 4. prepare-myo1d
    myo_step, myo_missing = _run_myo1d_step(ctx, myo1d_source, profile)
    if _fold(myo_step, myo_missing):
        return _finalize(ctx, aggregate)

    # 5. manifest-ligands (optional)
    if skip_ligands:
        aggregate.sub_steps.append(
            SubStepResult(
                name="manifest-ligands", status="SKIPPED", manifest_path=None, warnings=[]
            )
        )
    else:
        ligand_step = _run_ligand_step(
            ctx, ligands_dir, private_mapping, profile, compound_stage_enabled
        )
        _fold(ligand_step)

    return _finalize(ctx, aggregate)


def _finalize(
    ctx: RunContext, aggregate: PrepareInputsAggregate
) -> PrepareInputsAggregate:
    fail_count = sum(1 for s in aggregate.sub_steps if s.status == "FAIL")
    warn_count = sum(
        1 for s in aggregate.sub_steps if s.status in {"WARN", "PASS_WITH_WARNINGS"}
    )
    if fail_count:
        aggregate.status = "FAIL"
    elif warn_count:
        aggregate.status = "PASS_WITH_WARNINGS"
    else:
        aggregate.status = "PASS"

    aggregate.aggregate_manifest_path = _write_aggregate_manifest(ctx, aggregate)
    aggregate.summary_report_path = _write_summary_report(ctx, aggregate)

    phase_status_value = "WARN" if aggregate.status == "PASS_WITH_WARNINGS" else aggregate.status
    append_phase_status(
        ctx,
        phase="prepare-inputs",
        status=phase_status_value,
        message=(
            "prepare-inputs aggregate={0} substeps={1} fails={2} warns={3} missing={4}".format(
                aggregate.status,
                len(aggregate.sub_steps),
                fail_count,
                warn_count,
                len(aggregate.missing_required_inputs),
            )
        ),
        details={
            "mode": aggregate.mode,
            "profile": aggregate.profile,
            "input_root": aggregate.input_root,
            "states_processed": aggregate.states_processed,
            "missing_required_inputs": aggregate.missing_required_inputs,
            "blockers": aggregate.blockers,
            "aggregate_manifest": str(aggregate.aggregate_manifest_path),
            "summary_report": str(aggregate.summary_report_path),
            "score_bonus_allowed": False,
        },
    )
    return aggregate
