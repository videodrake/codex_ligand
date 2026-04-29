"""M1 ↔ Tasks 4-9 alignment helper (Phase 9).

After M1 (Phases 1-8) emits canonical outputs into ``fresh/runs/<run_id>/``,
the existing Tasks 4-9 modules can be invoked in the SAME run directory and
their outputs naturally coexist. This module records the cross-reference so
downstream consumers (M2 actual execution, future tooling) can locate the
M1 normalized artifacts that should feed Task 4-9 logic.

Approach
--------
Phase 9 takes an **additive** approach: Tasks 4-9 module logic is unchanged;
their existing 79 tests are unmodified. ``record_m1_alignment(ctx)`` walks
``ctx.run_dir`` for M1 canonical artifacts and writes a single
``manifest/m1_alignment.json`` mapping that ties M1 outputs to the Task 4-9
inputs they would feed. Code-level realignment (e.g. Task 4 stopping
emission of its pass-through ``prepared/egfr/egfr_receptor_normalized.pdb``)
is deferred to actual M2 execution work, where the concrete schema needs
will drive the changes.

Why additive?
~~~~~~~~~~~~~
- Preserves the 79 passing Task 4-9 tests without churn.
- Avoids speculative refactoring of analysis logic (V924R reporting,
  terminal-artifact detection, active-face annotation, ATP-overlap mask,
  pose-acceptance policy, pocket-prioritization scoring) that was carefully
  written earlier and is referenced by docs.
- The actual M2 PyRosetta/Vina/fpocket runners (future M2 phases) will be
  built to consume M1 normalized outputs natively; that's where the
  concrete schema bindings should be made, not pre-emptively here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_phase_status
from egfr_myo1d.core.manifest import now_iso, write_json
from egfr_myo1d.core.run_context import RunContext


@dataclass
class M1AlignmentReport:
    run_id: str
    m1_outputs_present: dict[str, bool] = field(default_factory=dict)
    m1_canonical_paths: dict[str, str] = field(default_factory=dict)
    task_consumption_map: dict[str, list[str]] = field(default_factory=dict)
    aligned_count: int = 0
    expected_count: int = 0
    status: str = "PASS"
    notes: list[str] = field(default_factory=list)


# Mapping: M1 artifact key -> (relative path under run_dir, list of Task IDs that
# would naturally consume this artifact in a future M2 execution layer).
M1_TO_TASK_CONSUMPTION = {
    "membrane_frame_json": (
        "manifest/membrane_frame.json",
        ["task4_prepare_ppi_inputs", "task6_plan_ppi_sampling"],
    ),
    "myo1d_normalized_pdb": (
        "normalized/myo1d/MYO1D_955_1006.pdb",
        ["task4_prepare_ppi_inputs", "task6_plan_ppi_sampling"],
    ),
    "myo1d_construct_qc_csv": (
        "qc/myo1d_construct_qc.csv",
        ["task4_prepare_ppi_inputs"],
    ),
    "ligand_manifest_qc_csv": (
        "qc/ligand_manifest_qc.csv",
        ["task9_prioritize_pocket_candidates"],
    ),
    # Receptor outputs are state-aware (multiple per run); enumerated dynamically
    # by glob below. The "task_consumption_map" entry keys those scans.
}

PER_STATE_RECEPTOR_GLOBS = (
    ("normalized/receptors", "*_full_frame_explicit_AB.pdb"),
    ("normalized/receptors", "*_dockable_669_1014_explicit_AB.pdb"),
    ("normalized/receptors", "*_runtime_offset_receptor_only.pdb"),
    ("qc", "*_receptor_mapping.csv"),
    ("qc", "*_receptor_normalization_audit.csv"),
)

RECEPTOR_TASK_CONSUMERS = [
    "task4_prepare_ppi_inputs",
    "task6_plan_ppi_sampling",
    "task7_summarize_ppi_consensus",
    "task8_plan_pocket_discovery",
    "task9_prioritize_pocket_candidates",
]


def _safe_rel(path: Path, ctx: RunContext) -> str:
    try:
        return path.resolve().relative_to(ctx.run_dir.resolve()).as_posix()
    except ValueError:
        return str(path)


def record_m1_alignment(ctx: RunContext) -> M1AlignmentReport:
    """Scan ``ctx.run_dir`` for M1 canonical artifacts, write
    ``manifest/m1_alignment.json``, and return the alignment report.

    The function never modifies existing artifacts; it only records the
    cross-reference between M1 outputs and the Task 4-9 modules that would
    consume them in a future M2 actual-execution phase.
    """
    report = M1AlignmentReport(run_id=ctx.run_id)

    # Scan for fixed-path M1 artifacts
    for key, (rel_path, consumers) in M1_TO_TASK_CONSUMPTION.items():
        full_path = ctx.run_dir / rel_path
        present = full_path.is_file()
        report.m1_outputs_present[key] = present
        report.expected_count += 1
        if present:
            report.aligned_count += 1
            report.m1_canonical_paths[key] = _safe_rel(full_path, ctx)
            for task_id in consumers:
                report.task_consumption_map.setdefault(task_id, []).append(key)

    # Scan for per-state receptor artifacts
    receptor_artifacts: dict[str, list[str]] = {}
    for subdir, pattern in PER_STATE_RECEPTOR_GLOBS:
        target_dir = ctx.run_dir / subdir
        if not target_dir.is_dir():
            continue
        for path in sorted(target_dir.glob(pattern)):
            key = "{0}/{1}".format(subdir, path.name)
            report.m1_outputs_present[key] = True
            report.aligned_count += 1
            report.m1_canonical_paths[key] = _safe_rel(path, ctx)
            receptor_artifacts.setdefault(subdir, []).append(path.name)

    if any(receptor_artifacts.values()):
        for task_id in RECEPTOR_TASK_CONSUMERS:
            report.task_consumption_map.setdefault(task_id, []).extend(
                "{0}/{1}".format(sub, name)
                for sub, names in receptor_artifacts.items()
                for name in names
            )

    if report.aligned_count == 0:
        report.status = "WARN"
        report.notes.append("no_m1_canonical_artifacts_found_in_run_dir")
    elif report.aligned_count < report.expected_count:
        report.status = "WARN"
        report.notes.append(
            "partial_m1_artifacts_present:{0}_of_{1}".format(
                report.aligned_count, report.expected_count
            )
        )

    payload: dict[str, Any] = {
        "run_id": report.run_id,
        "m1_outputs_present": report.m1_outputs_present,
        "m1_canonical_paths": report.m1_canonical_paths,
        "task_consumption_map": report.task_consumption_map,
        "aligned_count": report.aligned_count,
        "expected_count": report.expected_count,
        "status": report.status,
        "notes": report.notes,
        "approach": "additive_phase9_no_task_module_logic_changed",
        "score_bonus_allowed": False,
        "timestamp": now_iso(),
    }
    target = ctx.manifest_dir / "m1_alignment.json"
    write_json(target, payload, ctx)

    append_phase_status(
        ctx,
        phase="m1-alignment",
        status=report.status,
        message=(
            "m1-alignment status={0} aligned={1}/{2}".format(
                report.status, report.aligned_count, report.expected_count
            )
        ),
        details={
            "approach": "additive",
            "m1_alignment_manifest": _safe_rel(target, ctx),
            "consumer_tasks": list(report.task_consumption_map.keys()),
            "score_bonus_allowed": False,
        },
    )
    return report
