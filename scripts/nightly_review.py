#!/usr/bin/env python3
"""Build a detailed nightly review bundle for this repository."""

from __future__ import annotations

import argparse
import ast
import json
import subprocess
import warnings
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Iterable, List, Optional, Sequence, Tuple

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "example-project.yaml"
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "output" / "review_automation"

PRIORITY_ORDER = {"critical": 0, "high": 1, "medium": 2, "context": 3}
IGNORED_SCAN_PARTS = {".git", ".pytest_cache", "__pycache__", ".venv", "output/review_automation"}


@dataclass(frozen=True)
class ReviewUnit:
    unit_id: str
    title: str
    area: str
    priority: str
    path_patterns: Tuple[str, ...]
    focus_points: Tuple[str, ...]
    test_patterns: Tuple[str, ...] = ()
    doc_patterns: Tuple[str, ...] = ()


def unit(
    unit_id: str,
    title: str,
    area: str,
    priority: str,
    paths: Sequence[str],
    focus: Sequence[str],
    tests: Sequence[str] = (),
    docs: Sequence[str] = (),
) -> ReviewUnit:
    return ReviewUnit(
        unit_id=unit_id,
        title=title,
        area=area,
        priority=priority,
        path_patterns=tuple(paths),
        focus_points=tuple(focus),
        test_patterns=tuple(tests),
        doc_patterns=tuple(docs),
    )


def review_units() -> List[ReviewUnit]:
    return [
        unit(
            "control_plane",
            "Entry Points and Orchestration",
            "runtime",
            "critical",
            ("main.py", "run_production.py", "run_full_test.py"),
            (
                "Verify CLI dispatch, failure propagation, and operator-facing flow.",
                "Check skip/resume/force semantics and operational phase numbering.",
                "Review status output and recovery behavior for incomplete runs.",
            ),
            tests=("tests/test_smoke_cli.py", "tests/test_run_production.py", "tests/test_validation_smoke.py"),
            docs=("README.md", "docs/runbook.md", "docs/manual_execution.md", "docs/architecture.md"),
        ),
        unit(
            "runtime_core",
            "Runtime Core and Shared Utilities",
            "runtime",
            "critical",
            ("egfr_pipeline/config.py", "egfr_pipeline/runtime_support.py", "egfr_pipeline/residue_utils.py"),
            (
                "Check config loading, path normalization, and repo-root assumptions.",
                "Review residue-numbering helpers and shared metadata handling.",
                "Look for duplicated runtime logic or fragile helper behavior.",
            ),
            tests=("tests/test_vina_project_config.py", "tests/test_validation_smoke.py"),
            docs=("config/README.md", "docs/architecture.md"),
        ),
        unit(
            "workflow_utilities",
            "Workflow Utility Scripts",
            "ops",
            "high",
            ("scripts/precheck_guard.py", "scripts/rebuild_step_views.py", "scripts/reset_production_outputs.py", "scripts/run_pre_qsub_checks.sh", "scripts/setup_test_env.sh"),
            (
                "Check guardrails around destructive actions and scheduler gating.",
                "Review local-shell versus scheduler portability assumptions.",
                "Confirm helper scripts preserve canonical outputs and expected state.",
            ),
            tests=("tests/test_precheck_guard.py", "tests/test_rebuild_step_views.py"),
            docs=("docs/pre_qsub_test_line.md", "config/run_pre_qsub_checks.pbs"),
        ),
        unit(
            "step_output_layer",
            "Derived Step Output Layer",
            "outputs",
            "critical",
            ("egfr_pipeline/output_steps.py",),
            (
                "Verify canonical outputs remain the source of truth and step views stay additive.",
                "Check manifest freshness, stale-step handling, and regeneration safety.",
                "Look for drift between output docs, manifest fields, and code behavior.",
            ),
            tests=("tests/test_output_steps.py", "tests/test_step_e2e.py", "tests/test_step_modes.py", "tests/test_rebuild_step_views.py"),
            docs=("docs/architecture.md", "docs/output_artifact_map.md", "output/README.md"),
        ),
        unit(
            "vina_execution",
            "Routine Vina Execution",
            "vina",
            "critical",
            ("egfr_pipeline/vina/dock.py", "egfr_pipeline/vina/bootstrap.py", "egfr_pipeline/vina/sweep.py"),
            (
                "Check docking invocation, worker usage, and runtime-prepared input assumptions.",
                "Review reproducibility controls such as seeds and search parameters.",
                "Look for silent failures around missing binaries or malformed inputs.",
            ),
            tests=("tests/test_vina_runtime.py", "tests/test_vina_project_config.py"),
            docs=("docs/manual_vina.md", "config/example-project.yaml"),
        ),
        unit(
            "vina_postprocess",
            "Routine Vina Parsing, Contacts, and Clustering",
            "vina",
            "critical",
            ("egfr_pipeline/vina/parse_poses.py", "egfr_pipeline/vina/contacts.py", "egfr_pipeline/vina/cluster.py", "egfr_pipeline/vina/summarize.py", "egfr_pipeline/vina/compare.py"),
            (
                "Check malformed pose handling, contact-distance logic, and null safety.",
                "Review pocket merge heuristics and cross-state comparison contracts.",
                "Look for unstable thresholds that could distort downstream verdicts.",
            ),
            tests=("tests/test_vina_postprocess.py", "tests/test_verdict.py", "tests/test_report_outputs.py"),
            docs=("docs/output_artifact_map.md", "docs/data_flow_guide.md"),
        ),
        unit(
            "verdict_reporting",
            "Verdict, Reporting, and Validation",
            "integration",
            "critical",
            ("egfr_pipeline/verdict.py", "egfr_pipeline/report.py", "egfr_pipeline/validate.py"),
            (
                "Check scoring semantics, evidence weighting, and false-confidence risks.",
                "Review user-facing reporting for overstated claims or stale references.",
                "Confirm validation coverage matches the current output contract.",
            ),
            tests=("tests/test_verdict.py", "tests/test_report_outputs.py", "tests/test_validation_smoke.py"),
            docs=("docs/output_artifact_map.md", "docs/runbook.md", "docs/architecture.md"),
        ),
        unit(
            "ppi_support",
            "PPI Extraction and Submission Support",
            "ppi",
            "high",
            ("egfr_pipeline/ppi/*.py",),
            (
                "Check residue extraction, project-root exports, and metadata preservation.",
                "Review separation between active PyRosetta paths and legacy AFM compatibility.",
                "Look for mismatches between prepared inputs, extraction outputs, and scheduler helpers.",
            ),
            tests=("tests/test_postprocess_ppi.py", "tests/test_pyrosetta_extract.py", "tests/test_phase1_smoke.py"),
            docs=("docs/manual_pyrosetta.md", "docs/data_inventory.md"),
        ),
        unit(
            "pyrosetta_core",
            "PyRosetta Docking Core",
            "ppi",
            "critical",
            ("egfr_pipeline/pyrosetta_docking/*.py",),
            (
                "Review heavy-compute orchestration, score bookkeeping, and metadata contracts.",
                "Check cluster and final ranking outputs consumed by downstream layers.",
                "Look for partial-run behaviors that can silently corrupt downstream assumptions.",
            ),
            tests=("tests/test_pyrosetta_metadata.py", "tests/test_pyrosetta_extract.py", "tests/test_phase1_smoke.py"),
            docs=("docs/manual_pyrosetta.md", "docs/architecture.md"),
        ),
        unit(
            "md_surface",
            "MD Analysis Surface",
            "md",
            "medium",
            ("egfr_pipeline/md/*.py", "scripts/md_interface_analysis.py"),
            (
                "Check submenu entrypoints and standalone execution assumptions.",
                "Review output schema compatibility with the rest of the pipeline.",
                "Look for lightly-tested code living outside the main package tree.",
            ),
            tests=("tests/test_md_interface_analysis.py",),
            docs=("docs/manual_execution.md",),
        ),
        unit(
            "phase1_preparation",
            "Scientific Phase 1 Preparation",
            "phase1",
            "high",
            ("egfr_pipeline/phase1/generate_configs.py", "egfr_pipeline/phase1/prepare_inputs.py", "egfr_pipeline/phase1/launch_docking.py"),
            (
                "Check runtime input preparation against current construct metadata.",
                "Review generated config correctness and per-state output routing.",
                "Look for hidden assumptions about monomer/dimer constructs or seed layout.",
            ),
            tests=("tests/test_phase1_smoke.py",),
            docs=("docs/phase1_pyrosetta_execution_note.md", "docs/phase1_output_chain_note.md"),
        ),
        unit(
            "phase1_analysis",
            "Scientific Phase 1 Filtering, Validation, and Review",
            "phase1",
            "critical",
            ("egfr_pipeline/phase1/extract_interface.py", "egfr_pipeline/phase1/orientation_filter.py", "egfr_pipeline/phase1/standardize_scores.py", "egfr_pipeline/phase1/cluster_consensus.py", "egfr_pipeline/phase1/compare_states.py", "egfr_pipeline/phase1/pilot_comparison.py", "egfr_pipeline/phase1/lightdock_validation.py", "egfr_pipeline/phase1/review_report.py"),
            (
                "Check interface extraction, orientation filter semantics, and consensus logic.",
                "Review cross-state robustness tables and LightDock secondary validation handling.",
                "Look for report claims that outrun the underlying evidence or test coverage.",
            ),
            tests=("tests/test_phase1_orientation_filter.py", "tests/test_cluster_consensus.py", "tests/test_compare_states.py", "tests/test_lightdock_validation.py", "tests/test_lightdock_orientation.py", "tests/test_review_report.py"),
            docs=("docs/orientation_filter_design.md", "docs/phase1_lightdock_validation_note.md", "docs/phase1_sampling_rationale.md"),
        ),
        unit(
            "phase2_generation",
            "Scientific Phase 2 Pocket Generation",
            "phase2",
            "high",
            ("egfr_pipeline/phase2/patch_ingestion.py", "egfr_pipeline/phase2/pocket_proposal.py", "egfr_pipeline/phase2/pocket_merge.py"),
            (
                "Check Phase 1 handoff ingestion and normalized proposal schema.",
                "Review merge rules for overlapping pockets and provenance tracking.",
                "Look for brittle behavior when state coverage is partial or sparse.",
            ),
            tests=("tests/test_patch_ingestion_md_gate.py", "tests/test_pocket_proposal_fpocket.py", "tests/test_phase2.py"),
            docs=("docs/data_flow_guide.md", "docs/phase1_ppi_handoff_note.md"),
        ),
        unit(
            "phase2_export",
            "Scientific Phase 2 Scoring and Phase 3 Export",
            "phase2",
            "critical",
            ("egfr_pipeline/phase2/patch_relationship.py", "egfr_pipeline/phase2/druggability_confidence.py", "egfr_pipeline/phase2/cross_state_alignment.py", "egfr_pipeline/phase2/phase3_export.py", "egfr_pipeline/phase2/rerun_cascade.py", "egfr_pipeline/phase2/review_report.py"),
            (
                "Review relationship-class rules and hotspot-distance metrics.",
                "Check druggability confidence semantics and cross-state classification output.",
                "Verify Phase 3 export completeness and field naming stability.",
            ),
            tests=("tests/test_patch_relationship_state_scope.py", "tests/test_phase2.py"),
            docs=("docs/output_artifact_map.md", "docs/data_flow_guide.md"),
        ),
        unit(
            "phase3_surface",
            "Scientific Phase 3 Planning, Execution, and Export",
            "phase3",
            "critical",
            ("egfr_pipeline/phase3/*.py",),
            (
                "Check job construction, box definitions, budget policy, and execution bookkeeping.",
                "Review pose attribution, diversity metrics, and pending-vs-complete evidence states.",
                "Verify Phase 4 handoff fields remain stable and internally consistent.",
            ),
            tests=("tests/test_run_diverse_docking_execution.py", "tests/test_phase3.py"),
            docs=("docs/output_artifact_map.md", "docs/data_flow_guide.md"),
        ),
        unit(
            "phase4_surface",
            "Scientific Phase 4 Scoring and Reporting",
            "phase4",
            "critical",
            ("egfr_pipeline/phase4/*.py",),
            (
                "Review evidence ingestion, axis scoring, confidence modifiers, and mechanistic classification.",
                "Check final ranking tables and review/report outputs for contract drift.",
                "Look for narrative outputs that blur exploratory and robust evidence.",
            ),
            tests=("tests/test_phase4.py", "tests/test_phase4_affinity_cap.py"),
            docs=("docs/output_artifact_map.md", "docs/current_pipeline_status.md"),
        ),
        unit(
            "config_surfaces",
            "Config Surfaces and Scheduler Wrappers",
            "config",
            "critical",
            ("config/example-project.yaml", "config/full_test.yaml", "config/phase1/*.ini", "config/*.pbs", "config/README.md"),
            (
                "Check YAML, INI, and PBS boundaries for conflicting semantics.",
                "Review active-baseline versus legacy-key ambiguity and resource defaults.",
                "Look for config fields referenced in code but undocumented in current docs.",
            ),
            tests=("tests/test_vina_project_config.py", "tests/test_run_production.py"),
            docs=("docs/first_time_environment_setup.md", "docs/runbook.md"),
        ),
        unit(
            "input_contract",
            "Input Data Contract",
            "data",
            "high",
            ("input/receptors/*.pdb", "input/ligands/*.sdf", "input/PPI/*.pdb", "input/PPI/phase1/*", "input/PPI/prepared/*"),
            (
                "Check raw, prepared, and legacy inputs for provenance drift.",
                "Review naming consistency between checked-in files and config identifiers.",
                "Look for runtime-required prepared inputs that are absent from the workspace snapshot.",
            ),
            docs=("docs/data_inventory.md", "config/example-project.yaml"),
        ),
        unit(
            "output_contract",
            "Output Contract and Artifact Map",
            "outputs",
            "critical",
            ("output/README.md",),
            (
                "Review root-level output pointer behavior and canonical-versus-derived separation.",
                "Check output docs for stale filenames or misleading reading-order guidance.",
                "Look for mismatches between output contracts and validation logic.",
            ),
            tests=("tests/test_output_steps.py", "tests/test_rebuild_step_views.py", "tests/test_step_e2e.py", "tests/test_validation_smoke.py"),
            docs=("docs/output_artifact_map.md", "docs/data_inventory.md", "docs/architecture.md"),
        ),
        unit(
            "test_surface",
            "Test Surface and Regression Coverage",
            "quality",
            "critical",
            ("tests/*.py", "pytest.ini", "requirements-test.txt"),
            (
                "Check which critical runtime paths lack direct regression coverage.",
                "Review smoke-versus-full-test boundaries and marker semantics.",
                "Look for tests that encode stale behavior or miss current failure modes.",
            ),
            docs=("docs/test_suite_triage.md",),
        ),
        unit(
            "docs_surface",
            "Documentation and Onboarding Surface",
            "docs",
            "high",
            ("README.md", "CLAUDE.md", "docs/*.md", "docs/archive/*.md"),
            (
                "Check onboarding order and architecture docs for drift against current code.",
                "Review archive and active docs for conflicting instructions or stale terminology.",
                "Look for operator steps that do not match actual commands, outputs, or schedules.",
            ),
        ),
    ]


def _relative_to_repo(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def _is_ignored(relative_path: str) -> bool:
    normalized = relative_path.replace("\\", "/")
    return any(
        normalized == part or normalized.startswith(f"{part}/")
        for part in IGNORED_SCAN_PARTS
    )


def list_repo_files(repo_root: Path) -> List[Path]:
    files: List[Path] = []
    for path in repo_root.rglob("*"):
        if not path.is_file():
            continue
        rel_path = _relative_to_repo(path, repo_root)
        if _is_ignored(rel_path):
            continue
        files.append(path)
    return files


def _matches_pattern(relative_path: str, pattern: str) -> bool:
    return PurePosixPath(relative_path).match(pattern)


def _resolve_patterns(
    all_files: Sequence[Path],
    repo_root: Path,
    patterns: Sequence[str],
) -> Tuple[List[str], List[str]]:
    all_relative = [_relative_to_repo(path, repo_root) for path in all_files]
    resolved: List[str] = []
    missing: List[str] = []

    for pattern in patterns:
        matches = [rel for rel in all_relative if _matches_pattern(rel, pattern)]
        if matches:
            resolved.extend(matches)
            continue
        exact = repo_root / pattern
        if exact.exists() and exact.is_file():
            resolved.append(_relative_to_repo(exact, repo_root))
        else:
            missing.append(pattern)

    return sorted(set(resolved)), sorted(set(missing))


def _run_command(args: Sequence[str], repo_root: Path) -> Optional[str]:
    try:
        result = subprocess.run(
            args,
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def collect_git_summary(repo_root: Path) -> dict:
    branch = _run_command(("git", "branch", "--show-current"), repo_root)
    head = _run_command(("git", "rev-parse", "HEAD"), repo_root)
    if branch is None and head is None:
        return {
            "available": False,
            "branch": "",
            "head": "",
            "changed_files": [],
            "staged_files": [],
            "unstaged_files": [],
            "untracked_files": [],
            "recent_commits": [],
        }

    status_output = _run_command(("git", "status", "--short"), repo_root) or ""
    staged_files: List[str] = []
    unstaged_files: List[str] = []
    untracked_files: List[str] = []

    for raw_line in status_output.splitlines():
        if not raw_line.strip():
            continue
        code = raw_line[:2]
        path = raw_line[2:].strip()
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        path = path.replace("\\", "/")
        if code == "??":
            untracked_files.append(path)
            continue
        if code[0] != " ":
            staged_files.append(path)
        if code[1] != " ":
            unstaged_files.append(path)

    commit_output = _run_command(
        ("git", "log", "-5", "--date=iso-strict", "--pretty=format:%H%x1f%ad%x1f%s"),
        repo_root,
    ) or ""
    recent_commits = []
    for line in commit_output.splitlines():
        parts = line.split("\x1f", 2)
        if len(parts) != 3:
            continue
        sha, committed_at, subject = parts
        recent_commits.append(
            {"sha": sha, "committed_at": committed_at, "subject": subject}
        )

    changed_files = sorted(set(staged_files + unstaged_files + untracked_files))
    return {
        "available": True,
        "branch": branch or "",
        "head": head or "",
        "changed_files": changed_files,
        "staged_files": sorted(set(staged_files)),
        "unstaged_files": sorted(set(unstaged_files)),
        "untracked_files": sorted(set(untracked_files)),
        "recent_commits": recent_commits,
    }


def load_project_config_summary(repo_root: Path, config_path: Path) -> dict:
    if yaml is None or not config_path.exists():
        return {"path": str(config_path), "loaded": False}

    payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    receptors = payload.get("receptors") or []
    ligands = payload.get("ligands") or []
    return {
        "path": _relative_to_repo(config_path, repo_root),
        "loaded": True,
        "project_name": payload.get("project_name", ""),
        "output_root": payload.get("output_root", ""),
        "mode": payload.get("mode", ""),
        "max_workers": payload.get("max_workers"),
        "receptors": [str(row.get("id", "")) for row in receptors],
        "ligands": [str(row.get("id", "")) for row in ligands],
    }


def analyze_python_files(paths: Iterable[str], repo_root: Path) -> dict:
    python_files = [path for path in sorted(set(paths)) if path.endswith(".py")]
    function_count = 0
    class_count = 0
    syntax_errors: List[dict] = []

    for rel_path in python_files:
        file_path = repo_root / rel_path
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", SyntaxWarning)
                tree = ast.parse(
                    file_path.read_text(encoding="utf-8"),
                    filename=str(file_path),
                )
        except SyntaxError as exc:
            syntax_errors.append(
                {"path": rel_path, "line": exc.lineno or 0, "message": exc.msg}
            )
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                function_count += 1
            elif isinstance(node, ast.ClassDef):
                class_count += 1

    return {
        "python_file_count": len(python_files),
        "function_count": function_count,
        "class_count": class_count,
        "syntax_errors": syntax_errors,
    }


def build_unit_entries(
    repo_root: Path,
    all_files: Sequence[Path],
    changed_files: Sequence[str],
) -> List[dict]:
    changed_set = set(changed_files)
    entries: List[dict] = []

    for current in review_units():
        matched_files, missing_patterns = _resolve_patterns(
            all_files, repo_root, current.path_patterns
        )
        related_tests, missing_tests = _resolve_patterns(
            all_files, repo_root, current.test_patterns
        )
        related_docs, missing_docs = _resolve_patterns(
            all_files, repo_root, current.doc_patterns
        )
        changed_matches = sorted(
            path
            for path in matched_files + related_tests + related_docs
            if path in changed_set
        )
        if matched_files and not missing_patterns:
            status = "present"
        elif matched_files:
            status = "partial"
        else:
            status = "missing"

        newest_modified_at = ""
        if matched_files:
            newest = max(
                (repo_root / rel_path for rel_path in matched_files),
                key=lambda item: item.stat().st_mtime,
            )
            newest_modified_at = datetime.fromtimestamp(
                newest.stat().st_mtime, tz=timezone.utc
            ).isoformat()

        python_stats = analyze_python_files(matched_files + related_tests, repo_root)
        review_score = PRIORITY_ORDER.get(current.priority, 9) * 100
        if status == "missing":
            review_score += 40
        elif status == "partial":
            review_score += 10
        if changed_matches:
            review_score -= 35

        entries.append(
            {
                **asdict(current),
                "status": status,
                "matched_files": matched_files,
                "missing_patterns": missing_patterns,
                "related_tests": related_tests,
                "missing_tests": missing_tests,
                "related_docs": related_docs,
                "missing_docs": missing_docs,
                "changed_files": changed_matches,
                "newest_modified_at": newest_modified_at,
                "pytest_command": "pytest " + " ".join(related_tests) if related_tests else "",
                "review_score": review_score,
                **python_stats,
            }
        )

    entries.sort(key=lambda item: (item["review_score"], item["title"].lower()))
    for index, entry in enumerate(entries, start=1):
        entry["review_rank"] = index
    return entries


def build_manifest(
    repo_root: Path,
    config_path: Path,
    label: str,
    all_files: Sequence[Path],
    git_summary: dict,
    unit_entries: Sequence[dict],
) -> dict:
    markdown_count = sum(1 for path in all_files if path.suffix.lower() == ".md")
    python_count = sum(1 for path in all_files if path.suffix.lower() == ".py")
    config_count = sum(
        1
        for path in all_files
        if path.suffix.lower() in {".yaml", ".yml", ".ini", ".toml", ".pbs"}
    )
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "label": label,
        "repo_root": str(repo_root),
        "config_summary": load_project_config_summary(repo_root, config_path),
        "inventory_summary": {
            "file_count": len(all_files),
            "python_file_count": python_count,
            "markdown_file_count": markdown_count,
            "config_file_count": config_count,
            "review_unit_count": len(unit_entries),
        },
        "git": git_summary,
        "review_units": list(unit_entries),
    }


def render_checklist(manifest: dict, bundle_dir: Path) -> str:
    config_summary = manifest["config_summary"]
    git_summary = manifest["git"]
    lines = [
        "# Nightly Review Checklist",
        "",
        f"- Generated at (UTC): {manifest['generated_at_utc']}",
        f"- Bundle directory: `{bundle_dir}`",
        f"- Project name: `{config_summary.get('project_name', 'unknown')}`",
        f"- Config path: `{config_summary.get('path', '')}`",
        f"- Active receptors: {', '.join(config_summary.get('receptors', [])) or 'unknown'}",
        f"- Active ligands: {', '.join(config_summary.get('ligands', [])) or 'unknown'}",
        "",
        "## Goals",
        "",
        "1. Review changed files first, then every critical unit, then the remaining units in order.",
        "2. Prioritize bugs, regressions, brittle interfaces, config drift, doc drift, and missing tests.",
        "3. Findings should come first and include file references whenever possible.",
        "",
        "## Git Context",
        "",
        f"- Git available: `{git_summary.get('available', False)}`",
        f"- Branch: `{git_summary.get('branch', '') or 'unknown'}`",
        f"- Head: `{git_summary.get('head', '') or 'unknown'}`",
        f"- Changed files: `{len(git_summary.get('changed_files', []))}`",
        "",
        "## Review Order",
        "",
        "| Rank | Priority | Unit | Status | Changed | Files | Tests | Docs |",
        "|------|----------|------|--------|---------|-------|-------|------|",
    ]

    for current in manifest["review_units"]:
        lines.append(
            f"| {current['review_rank']} | {current['priority']} | {current['title']} | "
            f"{current['status']} | {len(current['changed_files'])} | {len(current['matched_files'])} | "
            f"{len(current['related_tests'])} | {len(current['related_docs'])} |"
        )

    lines.append("")
    for current in manifest["review_units"]:
        lines.extend(
            [
                f"## {current['review_rank']}. {current['title']} (`{current['unit_id']}`)",
                "",
                f"- Area: `{current['area']}`",
                f"- Priority: `{current['priority']}`",
                f"- Status: `{current['status']}`",
                f"- Changed matches: `{len(current['changed_files'])}`",
                f"- Python files: `{current['python_file_count']}`",
                f"- Functions/classes: `{current['function_count']}/{current['class_count']}`",
            ]
        )
        if current["newest_modified_at"]:
            lines.append(f"- Newest file mtime (UTC): `{current['newest_modified_at']}`")
        if current["changed_files"]:
            lines.append("- Changed files:")
            lines.extend(f"  - `{path}`" for path in current["changed_files"])
        if current["matched_files"]:
            lines.append("- Scope files:")
            lines.extend(f"  - `{path}`" for path in current["matched_files"])
        if current["missing_patterns"]:
            lines.append("- Missing expected patterns:")
            lines.extend(f"  - `{pattern}`" for pattern in current["missing_patterns"])
        if current["related_tests"]:
            lines.append("- Related tests:")
            lines.extend(f"  - `{path}`" for path in current["related_tests"])
        if current["related_docs"]:
            lines.append("- Related docs:")
            lines.extend(f"  - `{path}`" for path in current["related_docs"])
        if current["focus_points"]:
            lines.append("- Focus points:")
            lines.extend(f"  - {item}" for item in current["focus_points"])
        if current["pytest_command"]:
            lines.append(f"- Suggested pytest command: `{current['pytest_command']}`")
        if current["syntax_errors"]:
            lines.append("- Syntax errors:")
            for item in current["syntax_errors"]:
                lines.append(f"  - `{item['path']}` line {item['line']}: {item['message']}")
        lines.append("")

    return "\n".join(lines)


def render_prompt_packet(bundle_dir: Path) -> str:
    return "\n".join(
        [
            "# Nightly Review Prompt Packet",
            "",
            "1. Read `review_manifest.json` and `review_checklist.md` in this bundle.",
            "2. Review changed files first, then every critical unit, then the remaining units in checklist order.",
            "3. Focus on bugs, regressions, contract drift, config drift, doc drift, and missing tests.",
            "4. Run the smallest relevant pytest targets listed in the checklist when code is implicated.",
            "5. Overwrite `nightly_review_report.md` in this bundle with the final findings.",
            "",
            "Output rules:",
            "- Findings first, ordered by severity.",
            "- Cite file paths and line numbers when possible.",
            "- State which tests you ran and which you could not run.",
            "- If there are no findings, say that explicitly and list residual risks or coverage gaps.",
            "",
            f"Bundle path: `{bundle_dir}`",
        ]
    )


def render_report_template(manifest: dict, bundle_dir: Path) -> str:
    git_summary = manifest["git"]
    changed_count = len(git_summary.get("changed_files", []))
    return "\n".join(
        [
            "# Nightly Review Report",
            "",
            f"- Generated bundle: `{bundle_dir}`",
            f"- Prepared at (UTC): {manifest['generated_at_utc']}",
            f"- Branch: `{git_summary.get('branch', '') or 'unknown'}`",
            f"- Head: `{git_summary.get('head', '') or 'unknown'}`",
            f"- Changed files in git status: `{changed_count}`",
            "",
            "## Findings",
            "",
            "- Pending automation run.",
            "",
            "## Tests",
            "",
            "- Pending automation run.",
            "",
            "## Residual Risks",
            "",
            "- Pending automation run.",
            "",
        ]
    )


def create_review_bundle(
    repo_root: Path,
    output_root: Path,
    config_path: Path,
    label: Optional[str] = None,
) -> dict:
    timestamp = datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
    run_label = label or timestamp
    bundle_dir = output_root / run_label
    bundle_dir.mkdir(parents=True, exist_ok=True)

    all_files = list_repo_files(repo_root)
    git_summary = collect_git_summary(repo_root)
    unit_entries = build_unit_entries(
        repo_root=repo_root,
        all_files=all_files,
        changed_files=git_summary.get("changed_files", []),
    )
    manifest = build_manifest(
        repo_root=repo_root,
        config_path=config_path,
        label=run_label,
        all_files=all_files,
        git_summary=git_summary,
        unit_entries=unit_entries,
    )

    manifest_path = bundle_dir / "review_manifest.json"
    checklist_path = bundle_dir / "review_checklist.md"
    prompt_path = bundle_dir / "review_prompt.md"
    report_path = bundle_dir / "nightly_review_report.md"
    latest_path = output_root / "latest.json"

    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    checklist_path.write_text(render_checklist(manifest, bundle_dir), encoding="utf-8")
    prompt_path.write_text(render_prompt_packet(bundle_dir), encoding="utf-8")
    report_path.write_text(
        render_report_template(manifest, bundle_dir),
        encoding="utf-8",
    )
    latest_payload = {
        "label": run_label,
        "generated_at_utc": manifest["generated_at_utc"],
        "bundle_dir": str(bundle_dir),
        "manifest_path": str(manifest_path),
        "checklist_path": str(checklist_path),
        "prompt_path": str(prompt_path),
        "report_path": str(report_path),
    }
    latest_path.write_text(json.dumps(latest_payload, indent=2), encoding="utf-8")

    return {
        "bundle_dir": bundle_dir,
        "manifest_path": manifest_path,
        "checklist_path": checklist_path,
        "prompt_path": prompt_path,
        "report_path": report_path,
        "latest_path": latest_path,
        "manifest": manifest,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a detailed nightly review bundle for this repository."
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Directory used to store generated review bundles.",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Project config used for review context.",
    )
    parser.add_argument(
        "--label",
        default="",
        help="Optional bundle label. Defaults to a timestamp.",
    )
    args = parser.parse_args()

    result = create_review_bundle(
        repo_root=PROJECT_ROOT,
        output_root=args.output_root.resolve(),
        config_path=args.config.resolve(),
        label=args.label or None,
    )
    print(f"Nightly review bundle: {result['bundle_dir']}")
    print(f"Manifest: {result['manifest_path']}")
    print(f"Checklist: {result['checklist_path']}")
    print(f"Prompt: {result['prompt_path']}")
    print(f"Report template: {result['report_path']}")
    print(f"Latest pointer: {result['latest_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
