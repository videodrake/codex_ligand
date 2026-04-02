"""Derived step-view helpers for production outputs.

This module builds the additive step output layer described in the planning
docs. Canonical runtime outputs under ``output/{project}`` remain the source
of truth; step folders are curated copies and indexes for interpretation.
"""

from __future__ import annotations

import configparser
import csv
import json
import re
import shutil
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple, Union
from urllib.parse import urlencode

from egfr_pipeline.config import load_config
from egfr_pipeline import paths
from egfr_pipeline.parsing_utils import safe_float, safe_int
from egfr_pipeline.pyrosetta_docking.run_metadata import build_output_root_name


@dataclass(frozen=True)
class StepSpec:
    """Configuration for a derived step folder."""

    step_number: int
    phase_number: int
    step_name: str
    folder_name: str
    purpose: str
    summary_description: str
    primary_files: Tuple[str, ...]
    required_artifacts: Tuple[str, ...]
    optional_artifacts: Tuple[str, ...] = ()
    upstream_steps: Tuple[int, ...] = ()
    next_step_reads: Tuple[str, ...] = ()


STEP_SPECS: Dict[int, StepSpec] = {
    1: StepSpec(
        step_number=1,
        phase_number=1,
        step_name="vina_raw",
        folder_name="step1_vina_raw",
        purpose="Raw blind docking pose inventory across receptor states.",
        summary_description=(
            "This step captures the raw Vina pose outputs before pocket-level "
            "interpretation."
        ),
        primary_files=("raw_pose_index.csv",),
        required_artifacts=("raw_pose_index.csv",),
        upstream_steps=(),
        next_step_reads=("step4_vina_postprocess/vina_pocket_table.csv",),
    ),
    2: StepSpec(
        step_number=2,
        phase_number=2,
        step_name="ppi_raw",
        folder_name="step2_ppi_raw",
        purpose="Raw PyRosetta docking summaries and raw run references.",
        summary_description=(
            "This step captures the raw PPI docking outputs without copying the "
            "heavy raw run directories."
        ),
        primary_files=("raw_run_paths.tsv",),
        required_artifacts=("raw_run_paths.tsv",),
        optional_artifacts=("pyrosetta_run_metadata.json",),
        upstream_steps=(),
        next_step_reads=("step3_ppi_postprocess/ppi_pyrosetta_residues.csv",),
    ),
    3: StepSpec(
        step_number=3,
        phase_number=3,
        step_name="ppi_postprocess",
        folder_name="step3_ppi_postprocess",
        purpose="Receptor-side PPI residue evidence extracted from docking.",
        summary_description=(
            "This step converts raw PPI docking into receptor residue evidence, "
            "partner-aware summaries, and provenance-preserving long-form tables."
        ),
        primary_files=("ppi_pyrosetta_residues.csv", "ppi_pyrosetta_summary.csv"),
        required_artifacts=("ppi_pyrosetta_residues.csv", "ppi_pyrosetta_summary.csv"),
        optional_artifacts=(
            "ppi_pyrosetta_residue_long.csv",
            "ppi_pyrosetta_model_table.csv",
            "phase1_interface_report.md",
        ),
        upstream_steps=(2,),
        next_step_reads=(
            "step5_verdict/cross_method_agreement.csv",
            "step6_report/combined_residue_evidence.csv",
        ),
    ),
    4: StepSpec(
        step_number=4,
        phase_number=4,
        step_name="vina_postprocess",
        folder_name="step4_vina_postprocess",
        purpose="Pocket-level interpretation of Vina docking outputs.",
        summary_description=(
            "This step turns raw Vina poses into pocket tables, ligand-pocket "
            "mappings, coverage/provenance summaries, and cross-receptor comparisons."
        ),
        primary_files=(
            "vina_pocket_table.csv",
            "vina_postprocess_coverage.csv",
            "vina_pose_table.csv",
        ),
        required_artifacts=(
            "vina_pose_table.csv",
            "vina_pocket_table.csv",
            "vina_drug_pocket_map.csv",
            "vina_pocket_comparison.csv",
            "vina_pocket_bootstrap.csv",
            "vina_postprocess_coverage.csv",
        ),
        optional_artifacts=(
            "vina_contact_distances.csv",
            "vina_pocket_residue_occupancy.csv",
            "vina_clustering_merge_log.csv",
            "vina_clustering_parameters.json",
        ),
        upstream_steps=(1,),
        next_step_reads=(
            "step5_verdict/valid_sites.csv",
            "step6_report/project_report.txt",
        ),
    ),
    5: StepSpec(
        step_number=5,
        phase_number=5,
        step_name="verdict",
        folder_name="step5_verdict",
        purpose="Reproducibility-aware evidence classification for candidate sites.",
        summary_description=(
            "This step combines blind-docking, PPI, and cross-state recurrence "
            "signals to prioritize sites for interpretation."
        ),
        primary_files=("valid_sites.csv", "cross_method_agreement.csv"),
        required_artifacts=("valid_sites.csv", "cross_method_agreement.csv"),
        optional_artifacts=("vina_consensus_sites.csv",),
        upstream_steps=(3, 4),
        next_step_reads=("step6_report/project_report.txt",),
    ),
    6: StepSpec(
        step_number=6,
        phase_number=6,
        step_name="report",
        folder_name="step6_report",
        purpose="Narrative project summary and combined residue evidence.",
        summary_description=(
            "This step is the report-first view for reading a completed run with "
            "exploratory versus recurrent evidence separated."
        ),
        primary_files=("project_report.txt", "combined_residue_evidence.csv"),
        required_artifacts=("project_report.txt", "combined_residue_evidence.csv"),
        upstream_steps=(4, 5),
        next_step_reads=("step7_validate/validation_status.json",),
    ),
    7: StepSpec(
        step_number=7,
        phase_number=7,
        step_name="validate",
        folder_name="step7_validate",
        purpose="Persisted validation status for output integrity.",
        summary_description=(
            "This step records structured validation results for the derived view."
        ),
        primary_files=("validation_status.json", "validation_summary.txt"),
        required_artifacts=("validation_status.json", "validation_summary.txt"),
        optional_artifacts=("validation_recovery_plan.json", "validation_recovery_playbook.md"),
        upstream_steps=(1, 2, 3, 4, 5, 6),
        next_step_reads=(),
    ),
}


STEP_NAME_TO_FOLDER = {
    spec.step_name: spec.folder_name for spec in STEP_SPECS.values()
}


# ---------------------------------------------------------------------------
# Phase-aware project root shim
# ---------------------------------------------------------------------------

# Artifact filenames that live in specific phase directories
_ARTIFACT_ROUTES: Dict[str, str] = {}
# Populated from STEP_SPECS: map filename → "vina_post" | "ppi_post" | "verdict" | "report"
for _spec in STEP_SPECS.values():
    _phase = _spec.phase_number
    for _name in (*_spec.required_artifacts, *_spec.optional_artifacts, *_spec.primary_files):
        if not _name:
            continue
        if _name.startswith("vina_"):
            _ARTIFACT_ROUTES[_name] = "vina_post"
        elif _name.startswith("ppi_"):
            _ARTIFACT_ROUTES[_name] = "ppi_post"
        elif _name in ("valid_sites.csv", "cross_method_agreement.csv", "vina_consensus_sites.csv"):
            _ARTIFACT_ROUTES[_name] = "verdict"
        elif _name in ("project_report.txt", "combined_residue_evidence.csv"):
            _ARTIFACT_ROUTES[_name] = "report"
# Extras not captured by STEP_SPECS
_ARTIFACT_ROUTES.update({
    "ppi_afm_residues.csv": "ppi_post",
    "ppi_afm_summary.csv": "ppi_post",
    "ppi_pyrosetta_residue_long.csv": "ppi_post",
    "ppi_pyrosetta_model_table.csv": "ppi_post",
    "vina_pocket_bootstrap.csv": "vina_post",
    "vina_contact_distances.csv": "vina_post",
    "vina_pocket_residue_occupancy.csv": "vina_post",
    "vina_clustering_merge_log.csv": "vina_post",
    "vina_clustering_parameters.json": "vina_post",
})


class _ProjectView:
    """Drop-in shim for the old flat project_root Path.

    Routes artifact file lookups to phase-specific directories in the
    new ``output/workflow_a/`` layout while keeping step directories
    and metadata files at the project directory level.
    """

    def __init__(self, project_dir: Path, abs_output_root: Path):
        self._dir = project_dir
        wa = abs_output_root / "workflow_a"
        self._routes = {
            "vina_post": wa / "phase4_vina_postprocess",
            "ppi_post": wa / "phase3_ppi_postprocess",
            "verdict": wa / "phase5_verdict",
            "report": wa / "phase6_report",
        }

    def __truediv__(self, other) -> Path:
        name = str(other)
        route_key = _ARTIFACT_ROUTES.get(name)
        if route_key:
            return self._routes[route_key] / name
        return self._dir / name

    # Delegate common Path operations to the underlying directory
    def mkdir(self, *args, **kwargs):
        self._dir.mkdir(*args, **kwargs)

    def resolve(self, strict: bool = False) -> Path:
        return self._dir.resolve(strict=strict)

    def exists(self) -> bool:
        return self._dir.exists()

    def is_dir(self) -> bool:
        return self._dir.is_dir()

    @property
    def parent(self) -> Path:
        return self._dir.parent

    @property
    def name(self) -> str:
        return self._dir.name

    def __str__(self) -> str:
        return str(self._dir)

    def __repr__(self) -> str:
        return f"_ProjectView({self._dir})"

    def __fspath__(self) -> str:
        return str(self._dir)

    def as_posix(self) -> str:
        return self._dir.as_posix()

    def iterdir(self):
        return self._dir.iterdir()

    @property
    def parts(self):
        return self._dir.parts


def _validation_artifact_phase_map() -> Dict[str, int]:
    artifact_phase: Dict[str, int] = {}
    for step_num, spec in STEP_SPECS.items():
        if step_num == 7:
            continue
        for name in (*spec.required_artifacts, *spec.optional_artifacts, *spec.primary_files):
            if name and name not in artifact_phase:
                artifact_phase[name] = spec.phase_number
    artifact_phase.setdefault("ppi_afm_residues.csv", 3)
    return artifact_phase


VALIDATION_ARTIFACT_PHASES = _validation_artifact_phase_map()

VALIDATION_CATEGORY_METADATA = {
    "missing_artifact": {
        "label": "Missing artifact",
        "priority_rank": 0,
        "action_type": "rerun_then_validate",
    },
    "schema": {
        "label": "Schema mismatch",
        "priority_rank": 0,
        "action_type": "rerun_then_validate",
    },
    "validation_error": {
        "label": "Validation runtime error",
        "priority_rank": 0,
        "action_type": "manual_then_validate",
    },
    "id_consistency": {
        "label": "ID consistency mismatch",
        "priority_rank": 1,
        "action_type": "manual_then_rerun",
    },
    "residue_consistency": {
        "label": "Residue numbering or identity review",
        "priority_rank": 1,
        "action_type": "manual_then_validate",
    },
    "handoff": {
        "label": "Repository or handoff file issue",
        "priority_rank": 1,
        "action_type": "manual_then_validate",
    },
    "empty_output": {
        "label": "Empty output",
        "priority_rank": 1,
        "action_type": "manual_then_rerun",
    },
    "coverage": {
        "label": "Coverage or parse gap",
        "priority_rank": 1,
        "action_type": "manual_then_rerun",
    },
    "coverage_gap": {
        "label": "Coverage gap",
        "priority_rank": 1,
        "action_type": "manual_then_rerun",
    },
    "schema_warning": {
        "label": "Schema warning",
        "priority_rank": 2,
        "action_type": "manual_then_validate",
    },
    "validation_check": {
        "label": "Validation follow-up",
        "priority_rank": 2,
        "action_type": "manual_then_validate",
    },
}

VALIDATION_ACTION_LABELS = {
    "rerun_then_validate": "Rerun upstream phase, then validate",
    "manual_then_rerun": "Manual review first, then rerun upstream phase",
    "manual_then_validate": "Manual review first, then rerun validation",
}

VALIDATION_PRIORITY_LABELS = {
    0: "immediate",
    1: "high",
    2: "medium",
    3: "low",
}

VALIDATION_SEVERITY_RANKS = {
    "error": 0,
    "failure": 1,
    "warning": 2,
}

OPERATIONAL_CATEGORY_METADATA = {
    "step_not_generated": {
        "label": "Step view missing",
        "priority_rank": 0,
        "action_type": "rerun_upstream",
    },
    "missing_raw_pose": {
        "label": "Missing raw pose coverage",
        "priority_rank": 0,
        "action_type": "rerun_upstream",
    },
    "missing_ppi_ranking": {
        "label": "Missing PyRosetta ranking summary",
        "priority_rank": 0,
        "action_type": "rerun_upstream",
    },
    "missing_ppi_evidence": {
        "label": "Missing PPI evidence table",
        "priority_rank": 0,
        "action_type": "rerun_upstream",
    },
    "stale_step_view": {
        "label": "Stale derived step view",
        "priority_rank": 1,
        "action_type": "refresh_step_view",
    },
    "metadata_review": {
        "label": "PyRosetta metadata review",
        "priority_rank": 1,
        "action_type": "manual_then_rerun",
    },
    "fallback_source": {
        "label": "Fallback canonical source in use",
        "priority_rank": 1,
        "action_type": "manual_review",
    },
    "optional_supporting_artifact": {
        "label": "Optional supporting artifact missing",
        "priority_rank": 2,
        "action_type": "manual_review",
    },
}

OPERATIONAL_ACTION_LABELS = {
    "rerun_upstream": "Rerun upstream phase",
    "refresh_step_view": "Refresh this derived step view",
    "manual_then_rerun": "Manual review first, then rerun upstream phase",
    "manual_review": "Manual review only",
}


def step_output_view_enabled(
    config: Optional[dict],
    *,
    cli_disabled: bool = False,
) -> bool:
    """Return whether the additive step output layer should be generated."""

    if cli_disabled:
        return False
    if not config:
        return True

    step_config = config.get("step_output_view")
    if isinstance(step_config, dict) and "enabled" in step_config:
        return bool(step_config["enabled"])

    if "step_output_view_enabled" in config:
        return bool(config["step_output_view_enabled"])

    return True


def utc_now_iso() -> str:
    """Return an ISO-8601 UTC timestamp without fractional seconds."""

    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _as_path(value: Union[Path, str]) -> Path:
    return value if isinstance(value, Path) else Path(value)


def _repo_root_for(config_path: Path, repo_root: Optional[Union[Path, str]] = None) -> Path:
    if repo_root is not None:
        return _as_path(repo_root).resolve()
    resolved = config_path.resolve()
    if resolved.parent.name == "config" and len(resolved.parents) > 1:
        return resolved.parents[1]
    return resolved.parent


def resolve_project_root(
    config_path: Union[Path, str],
    config: Optional[dict] = None,
    repo_root: Optional[Union[Path, str]] = None,
) -> Path:
    """Resolve the canonical project root from a config file."""

    config_path = _as_path(config_path)
    config_data = config or load_config(str(config_path))
    root = Path(config_data.get("output_root", "./output"))
    base = _repo_root_for(config_path, repo_root)
    if not root.is_absolute():
        root = base / root
    project_name = config_data.get("project_name")
    return root / project_name if project_name else root


def _display_path(
    path: Union[Path, str],
    relative_to: Optional[Union[Path, str]] = None,
) -> str:
    path_obj = _as_path(path)
    if not path_obj.is_absolute():
        return path_obj.as_posix()
    if relative_to is not None:
        try:
            return (
                path_obj.resolve(strict=False)
                .relative_to(_as_path(relative_to).resolve(strict=False))
                .as_posix()
            )
        except ValueError:
            pass
    return path_obj.resolve(strict=False).as_posix()


def _slugify_group_token(value: object) -> str:
    token = re.sub(r"[^a-z0-9]+", "-", str(value or "").strip().lower()).strip("-")
    return token or "item"


def _recovery_group_key(source_kind: str, *parts: object) -> str:
    tokens = [_slugify_group_token(source_kind)]
    tokens.extend(_slugify_group_token(part) for part in parts if str(part or "").strip())
    return "-".join(tokens)


def _playbook_group_link(path: str, group_key: str) -> str:
    cleaned_path = str(path or "").strip()
    cleaned_group_key = str(group_key or "").strip()
    if cleaned_path and cleaned_group_key:
        return f"{cleaned_path}#{cleaned_group_key}"
    return cleaned_path


def _read_json(path: Path) -> dict:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def _atomic_write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(
        "w",
        delete=False,
        dir=path.parent,
        encoding="utf-8",
        newline="",
    ) as handle:
        handle.write(text)
        temp_path = Path(handle.name)
    temp_path.replace(path)
    return path


def _atomic_write_json(path: Path, data: dict) -> Path:
    text = json.dumps(data, indent=2, ensure_ascii=False, sort_keys=False) + "\n"
    return _atomic_write_text(path, text)


def _step_spec(step_num: int, step_name: Optional[str] = None) -> StepSpec:
    spec = STEP_SPECS.get(step_num)
    if spec is not None:
        return spec
    if step_name is None:
        raise KeyError(f"Unknown step number: {step_num}")
    folder_name = (
        step_name
        if step_name.startswith(f"step{step_num}_")
        else f"step{step_num}_{step_name}"
    )
    return StepSpec(
        step_number=step_num,
        phase_number=step_num,
        step_name=step_name.replace(f"step{step_num}_", ""),
        folder_name=folder_name,
        purpose="Derived step output.",
        summary_description="Derived step output.",
        primary_files=(),
        required_artifacts=(),
    )


def ensure_step_dir(
    project_root: Union[Path, str],
    step_num: int,
    step_name: Optional[str] = None,
) -> Path:
    """Ensure a derived step directory exists and return it."""

    spec = _step_spec(step_num, step_name)
    step_dir = _as_path(project_root) / spec.folder_name
    step_dir.mkdir(parents=True, exist_ok=True)
    return step_dir


@contextmanager
def _staged_step_dir(project_root: Path, step_num: int) -> Iterator[Tuple[Path, Path, bool]]:
    spec = _step_spec(step_num)
    final_dir = project_root / spec.folder_name
    temp_dir = project_root / f".{spec.folder_name}.tmp"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)
    existed = final_dir.exists()
    try:
        yield temp_dir, final_dir, existed
    except Exception:
        shutil.rmtree(temp_dir, ignore_errors=True)
        raise
    else:
        if final_dir.exists():
            shutil.rmtree(final_dir)
        temp_dir.replace(final_dir)


def write_step_manifest(step_dir: Union[Path, str], manifest_data: dict) -> Path:
    """Write ``step_manifest.json`` for a step directory."""

    return _atomic_write_json(_as_path(step_dir) / "step_manifest.json", manifest_data)


def write_current_run_manifest(project_root: Union[Path, str], manifest_data: dict) -> Path:
    """Write ``current_run_manifest.json`` at the project root."""

    return _atomic_write_json(
        _as_path(project_root) / "current_run_manifest.json",
        manifest_data,
    )


def _format_summary_entries(entries: Sequence[dict]) -> List[str]:
    lines: List[str] = []
    for entry in entries:
        path_text = entry.get("path", "")
        label = f"`{path_text}`" if path_text else "artifact"
        description = entry.get("description", "").strip()
        status = entry.get("status", "").strip()
        detail = description
        if status:
            detail = f"{detail} ({status})" if detail else status
        if detail:
            lines.append(f"- {label}: {detail}")
        else:
            lines.append(f"- {label}")
    return lines


def write_step_summary(
    step_dir: Union[Path, str],
    step_num: int,
    description: str,
    key_files: Sequence[dict],
    next_step_reads: Sequence[dict],
    warnings: Optional[Sequence[str]] = None,
) -> Path:
    """Write ``summary.md`` with stable headings for a step folder."""

    spec = _step_spec(step_num)
    lines = [
        f"# Step {step_num}: {spec.step_name}",
        "",
        "## Scientific Meaning",
        description,
        "",
        "## Key Files",
    ]
    lines.extend(_format_summary_entries(key_files) or ["- No step artifacts were captured."])
    lines.extend(["", "## Next Step Reads"])
    lines.extend(
        _format_summary_entries(next_step_reads)
        or ["- No downstream read order is defined for this step."]
    )
    if warnings:
        lines.extend(["", "## Notes", *[f"- {warning}" for warning in warnings]])
    lines.append("")
    return _atomic_write_text(_as_path(step_dir) / "summary.md", "\n".join(lines))


def copy_artifact_if_exists(
    source: Union[Path, str],
    dest: Union[Path, str],
    missing_list: List[str],
    artifact_label: Optional[str] = None,
) -> bool:
    """Copy a small canonical artifact when present and record misses."""

    source_path = _as_path(source)
    if not source_path.exists():
        missing_list.append(artifact_label or source_path.name)
        return False
    dest_path = _as_path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_path, dest_path)
    return True


def write_artifact_index(
    step_dir: Union[Path, str],
    filename: str,
    rows: Iterable[dict],
    fieldnames: Sequence[str],
) -> Path:
    """Write a CSV or TSV artifact index atomically."""

    step_dir = _as_path(step_dir)
    target = step_dir / filename
    delimiter = "\t" if target.suffix.lower() == ".tsv" else ","
    target.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(
        "w",
        delete=False,
        dir=target.parent,
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter=delimiter)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
        temp_path = Path(handle.name)
    temp_path.replace(target)
    return target


def _artifact_entry(
    name: str,
    source_path: Path,
    repo_root: Path,
    required: bool,
    copied: bool,
) -> dict:
    entry = {
        "name": name,
        "required": required,
        "status": "copied" if copied else "missing",
        "canonical_path": _display_path(source_path, repo_root),
        "step_path": name if copied else "",
    }
    if source_path.exists():
        stat = source_path.stat()
        entry["source_mtime"] = datetime.fromtimestamp(
            stat.st_mtime,
            tz=timezone.utc,
        ).replace(microsecond=0).isoformat().replace("+00:00", "Z")
        entry["size_bytes"] = stat.st_size
    return entry


def _collect_base_context(
    config_path: Union[Path, str],
    repo_root: Optional[Union[Path, str]] = None,
    *,
    create_project_root: bool = True,
) -> Tuple[Path, Path, dict, Path]:
    config_path = _as_path(config_path)
    config = load_config(str(config_path))
    repo_root_path = _repo_root_for(config_path, repo_root)
    project_root = resolve_project_root(config_path, config=config, repo_root=repo_root_path)
    if create_project_root:
        project_root.mkdir(parents=True, exist_ok=True)
    return config_path, repo_root_path, config, project_root


def _config_receptor_ids(config: dict) -> List[str]:
    return [str(item.get("id", "")) for item in config.get("receptors", []) if item.get("id")]


def _config_ligand_ids(config: dict) -> List[str]:
    return [str(item.get("id", "")) for item in config.get("ligands", []) if item.get("id")]


def _manifest_status(required_missing: Sequence[str], copied_required_count: int) -> str:
    if not required_missing:
        return "complete"
    if copied_required_count > 0:
        return "partial"
    return "missing"


def _vina_mode(config: dict) -> str:
    return str(config.get("mode", (config.get("vina") or {}).get("mode", "blind")))


def _vina_exhaustiveness(config: dict) -> int:
    vina_cfg = config.get("vina") or {}
    return int(vina_cfg.get("exhaustiveness", config.get("exhaustiveness", 0) or 0))


def _vina_n_poses(config: dict) -> int:
    vina_cfg = config.get("vina") or {}
    return int(vina_cfg.get("n_poses", config.get("n_poses", 0) or 0))


def _ligand_output_name(ligand: dict) -> str:
    pdbqt_path = ligand.get("pdbqt")
    if pdbqt_path:
        return Path(str(pdbqt_path)).stem.replace("_ligand", "")
    ligand_id = ligand.get("id")
    if ligand_id:
        return str(ligand_id)
    raise KeyError("Ligand entry is missing both 'id' and 'pdbqt'.")


def _count_pdbqt_models(path: Path) -> int:
    if not path.exists():
        return 0
    model_count = 0
    has_content = False
    with open(path, encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.strip():
                has_content = True
            if line.startswith("MODEL"):
                model_count += 1
    if model_count:
        return model_count
    return 1 if has_content else 0


def _step2_target_slug(name: str) -> str:
    return name.replace("-", "_").replace(" ", "_")


def _resolve_ppi_target_dir(
    target: dict,
    repo_root: Path,
) -> Optional[Path]:
    direct_dir = target.get("docking_dir")
    if direct_dir:
        path = Path(str(direct_dir))
        return path if path.is_absolute() else repo_root / path

    config_ini = target.get("config_ini")
    input_pdb_hint = target.get("input_pdb")
    if not config_ini or not input_pdb_hint:
        return None

    config_path = repo_root / str(config_ini)
    if not config_path.exists():
        return None

    parser = configparser.ConfigParser()
    parser.read(config_path, encoding="utf-8")

    input_pdb_rel = parser.get("Path", "input_pdb_name", fallback=str(input_pdb_hint))
    input_pdb = repo_root / input_pdb_rel
    if not input_pdb.exists():
        alt_input = repo_root / str(input_pdb_hint)
        if alt_input.exists():
            input_pdb = alt_input
        else:
            return None

    root_name = build_output_root_name(parser, str(input_pdb), input_pdb.stem)
    return repo_root / root_name


def _ppi_target_ranking_path(run_dir: Path) -> Path:
    preferred = run_dir / "final_result" / "final_ranking.csv"
    if preferred.exists():
        return preferred
    fallback = run_dir / "final_ranking.csv"
    return preferred if preferred.exists() else fallback


def _ppi_target_metadata_path(run_dir: Path, explicit_path: Optional[Union[Path, str]] = None) -> Optional[Path]:
    if explicit_path:
        path = _as_path(explicit_path)
        if path.exists():
            return path
    for candidate in (
        run_dir / "pyrosetta_run_metadata.json",
        run_dir / "final_result" / "pyrosetta_run_metadata.json",
        run_dir / "restored" / "pyrosetta_run_metadata.json",
    ):
        if candidate.exists():
            return candidate
    return None


def _resolve_project_artifact_path(project_root: Path, name: str) -> Optional[Path]:
    preferred = project_root / name
    if preferred.exists():
        return preferred
    fallback = project_root / "ppi" / name
    if fallback.exists():
        return fallback
    # New layout: check per-phase directories under workflow_a/
    wa_root = project_root.parent / "workflow_a" if project_root.parent else None
    if wa_root and wa_root.is_dir():
        for phase_dir_name in (
            "phase4_vina_postprocess",
            "phase3_ppi_postprocess",
            "phase5_verdict",
            "phase6_report",
            "phase7_validation",
        ):
            candidate = wa_root / phase_dir_name / name
            if candidate.exists():
                return candidate
    return None


def _phase1_interface_report_path(repo_root: Path) -> Optional[Path]:
    # New layout
    new_path = repo_root / "output" / "workflow_b" / "phase1_ppi_analysis" / "phase1_interface_report.md"
    if new_path.exists():
        return new_path
    # Legacy fallback
    legacy = repo_root / "output" / "phase1_ppi" / "phase1_interface_report.md"
    return legacy if legacy.exists() else None


def _validation_message_groups(messages: Sequence[str]) -> Tuple[List[str], List[str]]:
    missing_files: List[str] = []
    schema_errors: List[str] = []
    for message in messages:
        lowered = message.lower()
        if "schema" in lowered or "column" in lowered:
            schema_errors.append(message)
            continue
        if "missing" in lowered or "not found" in lowered:
            missing_files.append(message)
    return missing_files, schema_errors


def _phase_title(phase_number: Optional[int]) -> str:
    if phase_number in STEP_SPECS:
        spec = STEP_SPECS[int(phase_number)]
        return f"Phase {phase_number}: {spec.step_name.replace('_', ' ').title()}"
    return "Manual review"


def _rerun_command_for_phase(phase_number: Optional[int]) -> str:
    if phase_number in STEP_SPECS and int(phase_number) < 7:
        return f"python run_production.py --from {int(phase_number)}"
    return "python run_production.py --only 7"


def _validation_category_meta(category: str) -> dict:
    return dict(
        VALIDATION_CATEGORY_METADATA.get(
            category,
            {
                "label": "Validation follow-up",
                "priority_rank": 2,
                "action_type": "manual_then_validate",
            },
        )
    )


def _validation_action_label(action_type: str) -> str:
    return VALIDATION_ACTION_LABELS.get(action_type, "Manual review first, then rerun validation")


def _validation_priority_label(priority_rank: int) -> str:
    return VALIDATION_PRIORITY_LABELS.get(int(priority_rank), "low")


def _extract_validation_artifacts(message: str) -> List[str]:
    matches = re.findall(
        r"[A-Za-z0-9_./-]+\.(?:csv|txt|json|md|ya?ml|ini|pdb|pdbqt|py)",
        message,
    )
    artifacts: List[str] = []
    seen = set()
    for match in matches:
        name = Path(match).name
        if name not in seen:
            seen.add(name)
            artifacts.append(name)
    return artifacts


def _validation_issue_category(message: str) -> str:
    lowered = message.lower()
    if "schema" in lowered or "missing expected columns" in lowered or "could not read header" in lowered:
        return "schema"
    if "extra columns" in lowered:
        return "schema_warning"
    if "required doc missing" in lowered or "optional doc missing" in lowered or "module missing" in lowered:
        return "handoff"
    if "exists but is empty" in lowered:
        return "empty_output"
    if "missing" in lowered or "not found" in lowered:
        return "missing_artifact"
    if "unexpected receptor ids" in lowered or "unexpected ligand ids" in lowered:
        return "id_consistency"
    if "missing receptor ids" in lowered:
        return "coverage_gap"
    if "raw pose" in lowered or "raw_pose_file" in lowered or "coverage" in lowered:
        return "coverage"
    if (
        "overlap" in lowered
        or "identity mismatches" in lowered
        or "numbering offset" in lowered
        or "mutation" in lowered
        or "receptor pdb" in lowered
    ):
        return "residue_consistency"
    if "validation crashed" in lowered or "exception" in lowered:
        return "validation_error"
    return "validation_check"


def _validation_issue_phase(message: str) -> Tuple[Optional[int], Optional[str]]:
    for artifact_name in _extract_validation_artifacts(message):
        phase_number = VALIDATION_ARTIFACT_PHASES.get(artifact_name)
        if phase_number is not None:
            return phase_number, artifact_name

    lowered = message.lower()
    if "raw pose" in lowered or "raw_pose_file" in lowered:
        return 1, None
    if "vina postprocess coverage" in lowered:
        return 4, "vina_postprocess_coverage.csv"
    if "validation crashed" in lowered or "exception" in lowered:
        return 7, None
    return None, None


def _validation_resolution_hint(
    *,
    message: str,
    artifact_name: Optional[str],
    phase_number: Optional[int],
    category: str,
) -> str:
    if category == "validation_error":
        return "Fix the validation exception itself, then rerun Phase 7."
    if phase_number in STEP_SPECS and int(phase_number) < 7:
        phase_label = _phase_title(phase_number)
        if artifact_name:
            return f"Rebuild `{artifact_name}` starting from {phase_label} so downstream outputs stay in sync."
        return f"Inspect {phase_label} outputs, correct the issue, and rerun validation."
    if category == "handoff":
        return "Review the referenced repository file manually, then rerun Phase 7 if you changed anything."
    return "Inspect the referenced input or validation finding manually, then rerun Phase 7."


def _validation_checklist_for_issue(
    *,
    category: str,
    artifact_name: Optional[str],
    phase_number: Optional[int],
    recommended_command: str,
) -> List[str]:
    artifact_label = f"`{artifact_name}`" if artifact_name else "the affected output"
    phase_label = _phase_title(phase_number)
    if category == "missing_artifact":
        return [
            f"Confirm whether {artifact_label} is missing from the canonical project root, not only from the step view.",
            f"Review {phase_label} inputs and manifest data to see why the artifact was not produced.",
            f"Rerun with `{recommended_command}`, then rerun validation.",
        ]
    if category == "schema":
        return [
            f"Compare the columns in {artifact_label} against the schema expected by `egfr_pipeline/validate.py`.",
            f"Inspect the writer behind {phase_label} and confirm all required fields are still emitted.",
            f"Rerun with `{recommended_command}`, then rerun validation.",
        ]
    if category == "schema_warning":
        return [
            f"Check whether the extra columns in {artifact_label} are intentional and backward-compatible.",
            "If the schema change is expected, decide whether the validator should be updated to accept it.",
            "Rerun Phase 7 after that decision is made.",
        ]
    if category == "empty_output":
        return [
            f"Open {artifact_label} and confirm whether zero rows are expected for this project.",
            f"Review filters, thresholds, or coverage inputs upstream of {phase_label}.",
            f"If the file should not be empty, rerun with `{recommended_command}` before rerunning validation.",
        ]
    if category == "id_consistency":
        return [
            f"Compare receptor and ligand IDs in {artifact_label} against the project config.",
            "Check for renamed receptors, stale aliases, or mixed inputs from another run.",
            f"After correcting the IDs, rerun with `{recommended_command}` and then rerun validation.",
        ]
    if category in {"coverage", "coverage_gap"}:
        return [
            "Inspect the expected receptor/ligand pairs and the raw pose or coverage inputs that feed this check.",
            f"Confirm whether {artifact_label} reflects all expected parsed records.",
            f"If coverage is incomplete, rerun with `{recommended_command}` and validate again.",
        ]
    if category == "residue_consistency":
        return [
            "Inspect receptor PDB numbering, chain assignments, and mutation notes used in cross-receptor comparison.",
            "Decide whether cross-receptor interpretation is still safe, or whether the inputs need to be corrected first.",
            "After documenting or fixing the issue, rerun Phase 7. If you changed receptor inputs, rerun the affected upstream phases manually first.",
        ]
    if category == "handoff":
        return [
            "Restore the missing repository file or documentation artifact referenced by validation.",
            "Confirm the file lives in the expected repository path and is readable from this workspace.",
            "Rerun Phase 7 after the repository artifact has been restored.",
        ]
    if category == "validation_error":
        return [
            "Inspect the console traceback or run telemetry to identify where validation crashed.",
            "Fix the validation code path or malformed input that triggered the exception.",
            "Rerun Phase 7 after the validation exception is resolved.",
        ]
    return [
        "Review the validation message in context and identify the source input or artifact.",
        "Decide whether the issue requires a content fix, a rerun, or a validation-only follow-up.",
        "Rerun Phase 7 after the issue is resolved or documented.",
    ]


def _validation_issue_entry(message: str, severity: str) -> dict:
    phase_number, artifact_name = _validation_issue_phase(message)
    if severity == "error" and phase_number is None:
        phase_number = 7
    category = "validation_error" if severity == "error" else _validation_issue_category(message)
    meta = _validation_category_meta(category)
    action_type = str(meta["action_type"])
    recommended_command = _rerun_command_for_phase(phase_number)
    requires_manual_review = action_type != "rerun_then_validate" or phase_number is None
    return {
        "message": message,
        "severity": severity,
        "category": category,
        "category_label": str(meta["label"]),
        "artifact_name": artifact_name or "",
        "phase_number": phase_number,
        "phase_label": _phase_title(phase_number),
        "recommended_command": recommended_command,
        "requires_manual_review": requires_manual_review,
        "action_type": action_type,
        "action_label": _validation_action_label(action_type),
        "priority_rank": int(meta["priority_rank"]),
        "priority_label": _validation_priority_label(int(meta["priority_rank"])),
        "resolution_hint": _validation_resolution_hint(
            message=message,
            artifact_name=artifact_name,
            phase_number=phase_number,
            category=category,
        ),
        "checklist": _validation_checklist_for_issue(
            category=category,
            artifact_name=artifact_name,
            phase_number=phase_number,
            recommended_command=recommended_command,
        ),
    }


def _sorted_validation_issues(issues: Sequence[dict]) -> List[dict]:
    return sorted(
        issues,
        key=lambda issue: (
            int(issue.get("priority_rank", 2)),
            VALIDATION_SEVERITY_RANKS.get(str(issue.get("severity", "warning")), 9),
            99 if issue.get("phase_number") is None else int(issue.get("phase_number")),
            str(issue.get("message", "")),
        ),
    )


def _validation_action_groups(issues: Sequence[dict]) -> List[dict]:
    groups: Dict[Tuple[str, Optional[int], str], dict] = {}
    for issue in _sorted_validation_issues(issues):
        key = (
            str(issue.get("category", "validation_check")),
            issue.get("phase_number"),
            str(issue.get("recommended_command", "")),
        )
        group = groups.get(key)
        if group is None:
            group = {
                "category": issue["category"],
                "category_label": issue.get("category_label", issue["category"]),
                "phase_number": issue.get("phase_number"),
                "phase_label": issue.get("phase_label", "Manual review"),
                "recommended_command": issue.get("recommended_command", ""),
                "group_key": _recovery_group_key(
                    "validation",
                    issue.get("category", "validation_check"),
                    issue.get("phase_number"),
                    issue.get("recommended_command", ""),
                ),
                "action_type": issue.get("action_type", "manual_then_validate"),
                "action_label": issue.get("action_label", "Manual review first, then rerun validation"),
                "priority_rank": int(issue.get("priority_rank", 2)),
                "priority_label": issue.get("priority_label", "medium"),
                "manual_review_required": bool(issue.get("requires_manual_review")),
                "checklist": list(issue.get("checklist") or []),
                "issue_count": 0,
                "artifacts": [],
                "findings": [],
            }
            groups[key] = group
        group["issue_count"] += 1
        artifact_name = str(issue.get("artifact_name", "")).strip()
        if artifact_name and artifact_name not in group["artifacts"]:
            group["artifacts"].append(artifact_name)
        group["findings"].append(
            {
                "message": issue["message"],
                "severity": issue["severity"],
                "resolution_hint": issue["resolution_hint"],
            }
        )
    return sorted(
        groups.values(),
        key=lambda group: (
            int(group.get("priority_rank", 2)),
            99 if group.get("phase_number") is None else int(group.get("phase_number")),
            str(group.get("category_label", "")),
        ),
    )


def build_validation_recovery_plan(
    project_root: Union[Path, str],
    validation_status: dict,
) -> dict:
    project_root = _as_path(project_root)
    issues: List[dict] = []
    for message in validation_status.get("failure_messages") or []:
        issues.append(_validation_issue_entry(str(message), "failure"))
    for message in validation_status.get("warnings") or []:
        issues.append(_validation_issue_entry(str(message), "warning"))
    if validation_status.get("error_message"):
        issues.append(_validation_issue_entry(str(validation_status["error_message"]), "error"))
    issues = _sorted_validation_issues(issues)
    action_groups = _validation_action_groups(issues)

    upstream_phases = sorted(
        {
            int(issue["phase_number"])
            for issue in issues
            if issue.get("phase_number") in STEP_SPECS and int(issue["phase_number"]) < 7
        }
    )
    primary_phase = upstream_phases[0] if upstream_phases else (7 if issues else None)
    primary_command = _rerun_command_for_phase(primary_phase) if primary_phase is not None else ""
    primary_reason = ""
    if primary_phase is not None:
        primary_issue = next(
            (issue for issue in issues if issue.get("phase_number") == primary_phase),
            None,
        )
        if primary_issue is not None:
            primary_reason = str(primary_issue.get("resolution_hint", ""))

    rerun_commands: List[dict] = []
    seen_commands = set()
    for phase_number in upstream_phases:
        command = _rerun_command_for_phase(phase_number)
        if command in seen_commands:
            continue
        seen_commands.add(command)
        rerun_commands.append(
            {
                "phase_number": phase_number,
                "phase_label": _phase_title(phase_number),
                "command": command,
                "reason": f"Refresh outputs starting from {_phase_title(phase_number)}.",
            }
        )

    validation_only_command = "python run_production.py --only 7"
    if issues or validation_status.get("status") == "error":
        rerun_commands.append(
            {
                "phase_number": 7,
                "phase_label": _phase_title(7),
                "command": validation_only_command,
                "reason": "Use this after manual fixes or after upstream reruns to confirm validation is clean.",
            }
        )

    triage_counts = {
        "rerun_then_validate": len(
            [group for group in action_groups if group.get("action_type") == "rerun_then_validate"]
        ),
        "manual_then_rerun": len(
            [group for group in action_groups if group.get("action_type") == "manual_then_rerun"]
        ),
        "manual_then_validate": len(
            [group for group in action_groups if group.get("action_type") == "manual_then_validate"]
        ),
    }

    return {
        "status": str(validation_status.get("status", "unknown")),
        "generated_at": utc_now_iso(),
        "project_root": _display_path(project_root, project_root.parent.parent),
        "recommended_phase_number": primary_phase,
        "recommended_phase_label": _phase_title(primary_phase),
        "recommended_command": primary_command,
        "recommended_reason": primary_reason or "Review the findings below before rerunning validation.",
        "validation_only_command": validation_only_command,
        "manual_review_required": any(bool(issue.get("requires_manual_review")) for issue in issues),
        "issues": issues,
        "action_groups": action_groups,
        "triage": triage_counts,
        "rerun_commands": rerun_commands,
        "summary": {
            "issue_count": len(issues),
            "failure_count": safe_int(validation_status.get("failure_count"), 0),
            "warning_count": safe_int(validation_status.get("warning_count"), 0),
            "upstream_issue_count": len(
                [issue for issue in issues if issue.get("phase_number") in STEP_SPECS and int(issue["phase_number"]) < 7]
            ),
            "manual_review_count": len([issue for issue in issues if issue.get("requires_manual_review")]),
            "action_group_count": len(action_groups),
        },
    }


def write_validation_recovery_outputs(
    step_dir: Union[Path, str],
    recovery_plan: dict,
) -> Tuple[Path, Path]:
    """Write Step 7 recovery artifacts for validation follow-up."""

    step_dir = _as_path(step_dir)
    json_path = _atomic_write_json(step_dir / "validation_recovery_plan.json", recovery_plan)

    lines = [
        "# Validation Recovery Playbook",
        "",
        "## Status",
        f"- Validation status: `{recovery_plan['status']}`",
        f"- Finding count: `{recovery_plan['summary']['issue_count']}`",
        f"- Failure count: `{recovery_plan['summary']['failure_count']}`",
        f"- Warning count: `{recovery_plan['summary']['warning_count']}`",
    ]
    if recovery_plan.get("recommended_command"):
        lines.append(f"- Safest repair command: `{recovery_plan['recommended_command']}`")
        lines.append(f"- Earliest affected phase: `{recovery_plan['recommended_phase_label']}`")
    lines.append(f"- Validation-only rerun: `{recovery_plan['validation_only_command']}`")
    lines.append(
        f"- Manual review required: `{'yes' if recovery_plan.get('manual_review_required') else 'no'}`"
    )
    lines.append(f"- Action groups: `{recovery_plan['summary'].get('action_group_count', 0)}`")

    lines.extend(
        [
            "",
            "## Triage Summary",
            (
                f"- `{recovery_plan['triage'].get('rerun_then_validate', 0)}` group(s) can be repaired by rerunning an upstream phase."
            ),
            (
                f"- `{recovery_plan['triage'].get('manual_then_rerun', 0)}` group(s) need manual review before an upstream rerun."
            ),
            (
                f"- `{recovery_plan['triage'].get('manual_then_validate', 0)}` group(s) need manual review before a validation-only rerun."
            ),
            "",
            "## Why This Command",
            f"- {recovery_plan.get('recommended_reason') or 'No repair command is required.'}",
            "",
            "## Action Groups",
        ]
    )
    if recovery_plan["action_groups"]:
        for index, group in enumerate(recovery_plan["action_groups"], start=1):
            lines.extend(
                [
                    f"<a id=\"{group['group_key']}\"></a>",
                    f"### Action Group {index}: {group['category_label']}",
                    f"- Deep link: `{_playbook_group_link('validation_recovery_playbook.md', group['group_key'])}`",
                    f"- Priority: `{group['priority_label']}`",
                    f"- Action path: `{group['action_label']}`",
                    f"- Likely source phase: `{group['phase_label']}`",
                    f"- Suggested command: `{group['recommended_command']}`",
                    f"- Finding count: `{group['issue_count']}`",
                    "",
                ]
            )
            if group["artifacts"]:
                lines.append("- Related artifacts: " + ", ".join(f"`{name}`" for name in group["artifacts"]))
            lines.append("- Checklist:")
            lines.extend(f"  - {item}" for item in group["checklist"])
            lines.append("- Findings:")
            lines.extend(
                f"  - [{finding['severity'].upper()}] {finding['message']}"
                for finding in group["findings"]
            )
            lines.append("")
    else:
        lines.extend(
            [
                "- No grouped recovery actions were recorded.",
                "",
            ]
        )

    lines.append("## Individual Findings")
    if recovery_plan["issues"]:
        for index, issue in enumerate(recovery_plan["issues"], start=1):
            lines.extend(
                [
                    f"### Finding {index}",
                    f"- Message: {issue['message']}",
                    f"- Severity: `{issue['severity']}`",
                    f"- Category: `{issue['category_label']}`",
                    f"- Priority: `{issue['priority_label']}`",
                    f"- Action path: `{issue['action_label']}`",
                    f"- Likely source phase: `{issue['phase_label']}`",
                    f"- Suggested command: `{issue['recommended_command']}`",
                    f"- Guidance: {issue['resolution_hint']}",
                    "",
                ]
            )
    else:
        lines.extend(
            [
                "- No actionable validation findings were recorded.",
                "",
            ]
        )

    lines.append("## Rerun Options")
    if recovery_plan["rerun_commands"]:
        for item in recovery_plan["rerun_commands"]:
            lines.append(f"- `{item['command']}`: {item['reason']}")
    else:
        lines.append("- No rerun command is required based on the current validation record.")

    markdown_path = _atomic_write_text(
        step_dir / "validation_recovery_playbook.md",
        "\n".join(lines).rstrip() + "\n",
    )
    return json_path, markdown_path


def _validation_payload(
    project_root: Path,
    validation_result: Optional[object] = None,
    error_message: Optional[str] = None,
) -> Tuple[dict, List[str]]:
    warnings = list(getattr(validation_result, "warnings", []) or [])
    failures = list(getattr(validation_result, "failures", []) or [])
    passes = list(getattr(validation_result, "passes", []) or [])

    missing_files, schema_errors = _validation_message_groups([*failures, *warnings])
    if error_message:
        status = "error"
    elif failures:
        status = "failed"
    else:
        status = "passed"

    payload = {
        "status": status,
        "timestamp": utc_now_iso(),
        "project_root": _display_path(project_root, project_root.parent.parent),
        "missing_files": missing_files,
        "schema_errors": schema_errors,
        "warnings": warnings,
        "failure_messages": failures,
        "validated_steps": [1, 2, 3, 4, 5, 6, 7],
        "pass_count": len(passes),
        "warning_count": len(warnings),
        "failure_count": len(failures),
    }
    if error_message:
        payload["error_message"] = error_message
    return payload, failures


def write_validation_outputs(
    step_dir: Union[Path, str],
    validation_status: dict,
    summary_lines: Sequence[str],
) -> Tuple[Path, Path]:
    """Write Step 7 validation artifacts."""

    step_dir = _as_path(step_dir)
    status_path = _atomic_write_json(step_dir / "validation_status.json", validation_status)
    summary_path = _atomic_write_text(
        step_dir / "validation_summary.txt",
        "\n".join(summary_lines).rstrip() + "\n",
    )
    return status_path, summary_path


def _build_step_manifest(
    *,
    spec: StepSpec,
    config_path: Path,
    repo_root: Path,
    project_root: Path,
    config: dict,
    artifact_entries: Sequence[dict],
    missing_files: Sequence[str],
    warnings: Sequence[str],
    notes: str,
    regenerated: bool,
) -> dict:
    copied_required_count = sum(
        1 for entry in artifact_entries if entry["required"] and entry["status"] == "copied"
    )
    manifest = {
        "step_number": spec.step_number,
        "step_name": spec.step_name,
        "phase_number": spec.phase_number,
        "generated_at": utc_now_iso(),
        "project_name": config.get("project_name", project_root.name),
        "project_root": _display_path(project_root, repo_root),
        "source_config": _display_path(config_path, repo_root),
        "receptor_ids": _config_receptor_ids(config),
        "ligand_ids": _config_ligand_ids(config),
        "upstream_steps": list(spec.upstream_steps),
        "artifact_paths": [entry["step_path"] for entry in artifact_entries if entry["step_path"]],
        "notes": notes,
        "status": _manifest_status(missing_files, copied_required_count),
        "missing_files": list(missing_files),
        "warnings": list(warnings),
        "source_artifacts": list(artifact_entries),
    }
    if regenerated:
        manifest["regenerated_at"] = manifest["generated_at"]
    return manifest


def _record_copy_step(
    *,
    step_num: int,
    config_path: Union[Path, str],
    repo_root: Optional[Union[Path, str]] = None,
) -> Path:
    config_path, repo_root_path, config, project_root = _collect_base_context(
        config_path,
        repo_root=repo_root,
    )
    spec = _step_spec(step_num)
    missing_required: List[str] = []
    warnings: List[str] = []
    with _staged_step_dir(project_root, step_num) as (temp_dir, _, existed):
        artifact_entries: List[dict] = []

        for name in spec.required_artifacts:
            source = project_root / name
            copied = copy_artifact_if_exists(source, temp_dir / name, missing_required, name)
            artifact_entries.append(
                _artifact_entry(
                    name=name,
                    source_path=source,
                    repo_root=repo_root_path,
                    required=True,
                    copied=copied,
                )
            )

        for name in spec.optional_artifacts:
            source = project_root / name
            copied = copy_artifact_if_exists(source, temp_dir / name, [], name)
            artifact_entries.append(
                _artifact_entry(
                    name=name,
                    source_path=source,
                    repo_root=repo_root_path,
                    required=False,
                    copied=copied,
                )
            )
            if not copied:
                warnings.append(f"Optional artifact missing: {name}")

        notes = ""
        if missing_required:
            notes = "Some canonical artifacts were missing during step view generation."

        manifest = _build_step_manifest(
            spec=spec,
            config_path=config_path,
            repo_root=repo_root_path,
            project_root=project_root,
            config=config,
            artifact_entries=artifact_entries,
            missing_files=missing_required,
            warnings=warnings,
            notes=notes,
            regenerated=existed,
        )
        write_step_manifest(temp_dir, manifest)

        key_files = []
        for name in spec.primary_files:
            matching = next((entry for entry in artifact_entries if entry["name"] == name), None)
            if matching and matching["status"] == "copied":
                description = "Inspect this first." if name == spec.primary_files[0] else "Supporting artifact."
                key_files.append({"path": name, "description": description})
            else:
                key_files.append({"path": name, "description": "Expected artifact.", "status": "missing"})

        next_reads = [
            {"path": path, "description": "Downstream interpretation path."}
            for path in spec.next_step_reads
        ]
        if spec.step_number == 4:
            key_files[0]["description"] = "Pocket cluster view. Inspect this first."
            if len(key_files) > 1 and key_files[1]["path"] == "vina_postprocess_coverage.csv":
                key_files[1]["description"] = "Coverage/provenance matrix for expected receptor-ligand raw pose files."
        elif spec.step_number == 5:
            key_files[0]["description"] = "Final site prioritization. Inspect this first."
        elif spec.step_number == 6:
            key_files[0]["description"] = "Narrative report. Inspect this first."

        write_step_summary(
            temp_dir,
            step_num=step_num,
            description=spec.summary_description,
            key_files=key_files,
            next_step_reads=next_reads,
            warnings=[
                *warnings,
                *(
                    [f"Missing required artifacts: {', '.join(missing_required)}"]
                    if missing_required
                    else []
                ),
            ],
        )

    return project_root / spec.folder_name


def record_step1_outputs(
    config_path: Union[Path, str],
    repo_root: Optional[Union[Path, str]] = None,
) -> Path:
    """Capture Phase 1 canonical Vina raw outputs in ``step1_vina_raw``."""

    config_path, repo_root_path, config, project_root = _collect_base_context(
        config_path,
        repo_root=repo_root,
    )
    spec = _step_spec(1)
    mode = _vina_mode(config)
    exhaustiveness = _vina_exhaustiveness(config)
    n_poses = _vina_n_poses(config)
    missing_pairs: List[str] = []
    warnings: List[str] = []
    expected_pairs = 0
    found_pairs = 0

    with _staged_step_dir(project_root, 1) as (temp_dir, _, existed):
        rows: List[dict] = []
        source_artifacts: List[dict] = []

        for receptor in config.get("receptors", []):
            receptor_id = str(receptor.get("id", ""))
            if not receptor_id:
                continue
            for ligand in config.get("ligands", []):
                ligand_id = str(ligand.get("id", ""))
                ligand_name = _ligand_output_name(ligand)
                pose_path = project_root / receptor_id / f"{ligand_name}_{mode}.pdbqt"
                relative_pose_path = _display_path(pose_path, project_root)
                model_count = _count_pdbqt_models(pose_path)
                exists = pose_path.exists()
                expected_pairs += 1
                if exists:
                    found_pairs += 1
                else:
                    missing_pairs.append(relative_pose_path)

                rows.append(
                    {
                        "receptor_id": receptor_id,
                        "ligand_id": ligand_id,
                        "raw_pose_file": relative_pose_path if exists else "",
                        "n_models": model_count,
                        "source": "canonical_output" if exists else "missing",
                        "docking_mode": mode,
                        "exhaustiveness": exhaustiveness,
                        "n_poses": n_poses,
                    }
                )
                source_artifacts.append(
                    {
                        "name": f"{receptor_id}/{ligand_name}_{mode}.pdbqt",
                        "required": True,
                        "status": "indexed" if exists else "missing",
                        "canonical_path": relative_pose_path,
                        "step_path": "raw_pose_index.csv",
                        "n_models": model_count,
                    }
                )

        write_artifact_index(
            temp_dir,
            "raw_pose_index.csv",
            rows=rows,
            fieldnames=[
                "receptor_id",
                "ligand_id",
                "raw_pose_file",
                "n_models",
                "source",
                "docking_mode",
                "exhaustiveness",
                "n_poses",
            ],
        )

        if missing_pairs:
            warnings.append(
                f"Missing canonical raw pose files for {len(missing_pairs)} expected docking pair(s)."
            )
        notes = (
            f"Indexed {found_pairs} of {expected_pairs} expected receptor/ligand pose files."
        )
        manifest = _build_step_manifest(
            spec=spec,
            config_path=config_path,
            repo_root=repo_root_path,
            project_root=project_root,
            config=config,
            artifact_entries=[
                {
                    "name": "raw_pose_index.csv",
                    "required": True,
                    "status": "generated",
                    "canonical_path": "",
                    "step_path": "raw_pose_index.csv",
                }
            ],
            missing_files=missing_pairs,
            warnings=warnings,
            notes=notes,
            regenerated=existed,
        )
        if missing_pairs and found_pairs == 0:
            manifest["status"] = "missing"
        elif missing_pairs:
            manifest["status"] = "partial"
        else:
            manifest["status"] = "complete"
        manifest["source_artifacts"] = source_artifacts
        write_step_manifest(temp_dir, manifest)

        summary_description = (
            f"{spec.summary_description} Indexed {found_pairs} of {expected_pairs} expected "
            f"docking pairs without duplicating the raw .pdbqt files."
        )
        write_step_summary(
            temp_dir,
            step_num=1,
            description=summary_description,
            key_files=[
                {
                    "path": "raw_pose_index.csv",
                    "description": "Pose inventory with project-relative canonical raw pose paths. Inspect this first.",
                }
            ],
            next_step_reads=[
                {
                    "path": "step4_vina_postprocess/vina_pocket_table.csv",
                    "description": "Pocket interpretation built from these raw poses.",
                }
            ],
            warnings=[
                *warnings,
                f"Expected receptor/ligand pairs: {expected_pairs}.",
                f"Found canonical raw pose files: {found_pairs}.",
            ],
        )

    return project_root / spec.folder_name


def record_step2_outputs(
    config_path: Union[Path, str],
    repo_root: Optional[Union[Path, str]] = None,
    ppi_targets: Optional[Sequence[dict]] = None,
) -> Path:
    """Capture Phase 2 PyRosetta raw summaries in ``step2_ppi_raw``."""

    config_path, repo_root_path, config, project_root = _collect_base_context(
        config_path,
        repo_root=repo_root,
    )
    spec = _step_spec(2)
    targets = list(ppi_targets or [])
    missing_required: List[str] = []
    warnings: List[str] = []
    raw_run_rows: List[dict] = []
    metadata_records: Dict[str, dict] = {}

    with _staged_step_dir(project_root, 2) as (temp_dir, _, existed):
        artifact_entries: List[dict] = []

        for target in targets:
            target_name = str(target.get("name", "")).strip() or "unknown"
            target_slug = _step2_target_slug(target_name)
            step_ranking_name = f"{target_slug}_final_ranking.csv"
            run_dir = _resolve_ppi_target_dir(target, repo_root_path)
            ranking_path = _ppi_target_ranking_path(run_dir) if run_dir else Path()
            ranking_exists = bool(run_dir) and ranking_path.exists()

            if ranking_exists:
                copy_artifact_if_exists(ranking_path, temp_dir / step_ranking_name, [], step_ranking_name)
            else:
                missing_required.append(step_ranking_name)

            artifact_entries.append(
                {
                    "name": step_ranking_name,
                    "required": True,
                    "status": "copied" if ranking_exists else "missing",
                    "canonical_path": _display_path(ranking_path, repo_root_path) if run_dir else "",
                    "step_path": step_ranking_name if ranking_exists else "",
                }
            )

            metadata_path = _ppi_target_metadata_path(run_dir, target.get("metadata_path")) if run_dir else None
            if metadata_path:
                try:
                    metadata_records[target_name] = {
                        "source_path": _display_path(metadata_path, repo_root_path),
                        "data": _read_json(metadata_path),
                    }
                except (json.JSONDecodeError, OSError) as exc:
                    warnings.append(f"Unable to read metadata JSON for {target_name}: {exc}")
            else:
                warnings.append(f"Optional metadata JSON missing for {target_name}.")

            raw_run_rows.append(
                {
                    "target_name": target_name,
                    "raw_run_dir": _display_path(run_dir, repo_root_path) if run_dir else "",
                    "final_ranking_csv": _display_path(ranking_path, repo_root_path) if ranking_exists else "",
                    "metadata_json": _display_path(metadata_path, repo_root_path) if metadata_path else "",
                    "source": "canonical_output" if ranking_exists else "missing",
                }
            )

            if run_dir:
                warnings.append(
                    f"Raw run directory intentionally referenced, not duplicated: {_display_path(run_dir, repo_root_path)}"
                )

        write_artifact_index(
            temp_dir,
            "raw_run_paths.tsv",
            rows=raw_run_rows,
            fieldnames=[
                "target_name",
                "raw_run_dir",
                "final_ranking_csv",
                "metadata_json",
                "source",
            ],
        )
        artifact_entries.append(
            {
                "name": "raw_run_paths.tsv",
                "required": True,
                "status": "generated",
                "canonical_path": "",
                "step_path": "raw_run_paths.tsv",
            }
        )

        if metadata_records:
            _atomic_write_json(
                temp_dir / "pyrosetta_run_metadata.json",
                metadata_records,
            )
            artifact_entries.append(
                {
                    "name": "pyrosetta_run_metadata.json",
                    "required": False,
                    "status": "generated",
                    "canonical_path": "",
                    "step_path": "pyrosetta_run_metadata.json",
                }
            )

        notes = "Raw PyRosetta run directories are referenced by path and are not copied into the step view."
        manifest = _build_step_manifest(
            spec=spec,
            config_path=config_path,
            repo_root=repo_root_path,
            project_root=project_root,
            config=config,
            artifact_entries=artifact_entries,
            missing_files=missing_required,
            warnings=warnings,
            notes=notes,
            regenerated=existed,
        )
        write_step_manifest(temp_dir, manifest)

        key_files = [
            {
                "path": "raw_run_paths.tsv",
                "description": "Path index back to the canonical raw PyRosetta run directories. Inspect this first.",
            },
        ]
        for target in targets:
            target_name = str(target.get("name", "")).strip() or "unknown"
            target_slug = _step2_target_slug(target_name)
            step_ranking_name = f"{target_slug}_final_ranking.csv"
            key_files.append(
                {
                    "path": step_ranking_name,
                    "description": f"{target_name} ranking summary.",
                    "status": "missing" if step_ranking_name in missing_required else "",
                }
            )
        if metadata_records:
            key_files.append(
                {
                    "path": "pyrosetta_run_metadata.json",
                    "description": "Available metadata JSON aggregated by target.",
                }
            )

        write_step_summary(
            temp_dir,
            step_num=2,
            description=(
                f"{spec.summary_description} Ranking summaries are copied for quick review, "
                "while the heavy PyRosetta run directories stay in their canonical locations."
            ),
            key_files=key_files,
            next_step_reads=[
                {
                    "path": "step3_ppi_postprocess/ppi_pyrosetta_residues.csv",
                    "description": "Residue evidence extracted from these raw runs.",
                }
            ],
            warnings=warnings,
        )

    return project_root / spec.folder_name


def record_step3_outputs(
    config_path: Union[Path, str],
    repo_root: Optional[Union[Path, str]] = None,
) -> Path:
    """Capture Phase 3 PPI postprocess outputs in ``step3_ppi_postprocess``."""

    config_path, repo_root_path, config, project_root = _collect_base_context(
        config_path,
        repo_root=repo_root,
    )
    spec = _step_spec(3)
    missing_required: List[str] = []
    warnings: List[str] = []

    with _staged_step_dir(project_root, 3) as (temp_dir, _, existed):
        artifact_entries: List[dict] = []

        for name in ("ppi_pyrosetta_residues.csv", "ppi_pyrosetta_summary.csv"):
            source = _resolve_project_artifact_path(project_root, name)
            copied = False
            if source is not None:
                copy_artifact_if_exists(source, temp_dir / name, [], name)
                copied = True
                if source.parent.name == "ppi":
                    warnings.append(
                        f"Using fallback canonical source for {name}: {_display_path(source, repo_root_path)}"
                    )
            else:
                missing_required.append(name)

            artifact_entries.append(
                {
                    "name": name,
                    "required": True,
                    "status": "copied" if copied else "missing",
                    "canonical_path": _display_path(source, repo_root_path) if source else "",
                    "step_path": name if copied else "",
                }
            )

        for name in ("ppi_pyrosetta_residue_long.csv", "ppi_pyrosetta_model_table.csv"):
            source = _resolve_project_artifact_path(project_root, name)
            copied = False
            if source is not None:
                copy_artifact_if_exists(source, temp_dir / name, [], name)
                copied = True
                if source.parent.name == "ppi":
                    warnings.append(
                        f"Using fallback canonical source for {name}: {_display_path(source, repo_root_path)}"
                    )
            else:
                warnings.append(f"Optional artifact missing: {name}")

            artifact_entries.append(
                {
                    "name": name,
                    "required": False,
                    "status": "copied" if copied else "missing",
                    "canonical_path": _display_path(source, repo_root_path) if source else "",
                    "step_path": name if copied else "",
                }
            )

        interface_report = _phase1_interface_report_path(repo_root_path)
        if interface_report is not None:
            copy_artifact_if_exists(
                interface_report,
                temp_dir / "phase1_interface_report.md",
                [],
                "phase1_interface_report.md",
            )
            artifact_entries.append(
                {
                    "name": "phase1_interface_report.md",
                    "required": False,
                    "status": "copied",
                    "canonical_path": _display_path(interface_report, repo_root_path),
                    "step_path": "phase1_interface_report.md",
                }
            )
        else:
            warnings.append("Optional Phase 1 interface report is not available.")

        notes = (
            "Step 3 reflects the aggregated residue evidence used by verdict/report plus "
            "optional long-form provenance tables for run/seed tracing."
        )
        manifest = _build_step_manifest(
            spec=spec,
            config_path=config_path,
            repo_root=repo_root_path,
            project_root=project_root,
            config=config,
            artifact_entries=artifact_entries,
            missing_files=missing_required,
            warnings=warnings,
            notes=notes,
            regenerated=existed,
        )
        write_step_manifest(temp_dir, manifest)

        key_files = [
            {
                "path": "ppi_pyrosetta_residues.csv",
                "description": "Receptor-side interface residue evidence. Inspect this first.",
                "status": "missing" if "ppi_pyrosetta_residues.csv" in missing_required else "",
            },
            {
                "path": "ppi_pyrosetta_summary.csv",
                "description": "Aggregated PPI summary across registered runs.",
                "status": "missing" if "ppi_pyrosetta_summary.csv" in missing_required else "",
            },
        ]
        for optional_name, description in (
            (
                "ppi_pyrosetta_residue_long.csv",
                "Per-model residue provenance table for seed/run reproducibility checks.",
            ),
            (
                "ppi_pyrosetta_model_table.csv",
                "Per-model score table with run/seed metadata and interface counts.",
            ),
        ):
            if (temp_dir / optional_name).exists():
                key_files.append(
                    {
                        "path": optional_name,
                        "description": description,
                        "status": "",
                    }
                )
        if interface_report is not None:
            key_files.append(
                {
                    "path": "phase1_interface_report.md",
                    "description": "Optional Phase 1 human-readable interface report.",
                }
            )

        write_step_summary(
            temp_dir,
            step_num=3,
            description=(
                f"{spec.summary_description} This step copies the compact residue evidence tables "
                "that downstream verdict and report steps read from the canonical project root."
            ),
            key_files=key_files,
            next_step_reads=[
                {
                    "path": "step5_verdict/cross_method_agreement.csv",
                    "description": "Cross-method agreement calculated from this residue evidence.",
                },
                {
                    "path": "step6_report/combined_residue_evidence.csv",
                    "description": "Combined residue evidence table using these PPI summaries.",
                },
            ],
            warnings=warnings,
        )

    return project_root / spec.folder_name


def record_step4_outputs(
    config_path: Union[Path, str],
    repo_root: Optional[Union[Path, str]] = None,
) -> Path:
    """Capture Phase 4 canonical outputs in ``step4_vina_postprocess``."""

    return _record_copy_step(step_num=4, config_path=config_path, repo_root=repo_root)


def record_step5_outputs(
    config_path: Union[Path, str],
    repo_root: Optional[Union[Path, str]] = None,
) -> Path:
    """Capture Phase 5 canonical outputs in ``step5_verdict``."""

    return _record_copy_step(step_num=5, config_path=config_path, repo_root=repo_root)


def record_step6_outputs(
    config_path: Union[Path, str],
    repo_root: Optional[Union[Path, str]] = None,
) -> Path:
    """Capture Phase 6 canonical outputs in ``step6_report``."""

    return _record_copy_step(step_num=6, config_path=config_path, repo_root=repo_root)


def record_step7_outputs(
    config_path: Union[Path, str],
    repo_root: Optional[Union[Path, str]] = None,
    validation_result: Optional[object] = None,
    error_message: Optional[str] = None,
) -> Path:
    """Persist validation results in ``step7_validate``."""

    config_path, repo_root_path, config, project_root = _collect_base_context(
        config_path,
        repo_root=repo_root,
    )
    spec = _step_spec(7)
    validation_status, failures = _validation_payload(
        project_root,
        validation_result=validation_result,
        error_message=error_message,
    )
    recovery_plan = build_validation_recovery_plan(project_root, validation_status)

    summary_lines = [
        f"Overall status: {validation_status['status']}",
        "",
        "Missing artifacts:",
    ]
    if validation_status["missing_files"]:
        summary_lines.extend(f"- {item}" for item in validation_status["missing_files"])
    else:
        summary_lines.append("- None")

    summary_lines.extend(["", "Schema mismatches:"])
    if validation_status["schema_errors"]:
        summary_lines.extend(f"- {item}" for item in validation_status["schema_errors"])
    else:
        summary_lines.append("- None")

    summary_lines.extend(["", "Warning summary:"])
    if validation_status["warnings"]:
        summary_lines.extend(f"- {item}" for item in validation_status["warnings"])
    else:
        summary_lines.append("- None")

    if failures:
        summary_lines.extend(["", "Validation failures:"])
        summary_lines.extend(f"- {item}" for item in failures)

    summary_lines.extend(["", "Next action:"])
    if validation_status["status"] == "error":
        summary_lines.append("- Fix the validation exception and rerun Phase 7.")
    elif validation_status["status"] == "failed":
        summary_lines.append("- Resolve the failed validation checks before interpreting downstream results.")
    else:
        summary_lines.append("- Review step_index.md and proceed with report/verdict interpretation.")
    if recovery_plan.get("recommended_command"):
        summary_lines.append(
            f"- Safest repair command: {recovery_plan['recommended_command']} ({recovery_plan['recommended_phase_label']})."
        )
    summary_lines.append("- Open validation_recovery_playbook.md for phase-by-phase recovery guidance.")
    summary_lines.append(f"- Validation-only rerun: {recovery_plan['validation_only_command']}.")

    with _staged_step_dir(project_root, 7) as (temp_dir, _, existed):
        status_path, summary_path = write_validation_outputs(
            temp_dir,
            validation_status,
            summary_lines,
        )
        recovery_json_path, recovery_md_path = write_validation_recovery_outputs(
            temp_dir,
            recovery_plan,
        )

        artifact_entries = [
            {
                "name": "validation_status.json",
                "required": True,
                "status": "generated",
                "canonical_path": "",
                "step_path": status_path.name,
            },
            {
                "name": "validation_summary.txt",
                "required": True,
                "status": "generated",
                "canonical_path": "",
                "step_path": summary_path.name,
            },
            {
                "name": "validation_recovery_plan.json",
                "required": False,
                "status": "generated",
                "canonical_path": "",
                "step_path": recovery_json_path.name,
            },
            {
                "name": "validation_recovery_playbook.md",
                "required": False,
                "status": "generated",
                "canonical_path": "",
                "step_path": recovery_md_path.name,
            },
        ]

        notes = "Validation artifacts are derived from egfr_pipeline.validate and do not replace canonical outputs."
        manifest = _build_step_manifest(
            spec=spec,
            config_path=config_path,
            repo_root=repo_root_path,
            project_root=project_root,
            config=config,
            artifact_entries=artifact_entries,
            missing_files=validation_status["missing_files"],
            warnings=validation_status["warnings"],
            notes=notes,
            regenerated=existed,
        )
        manifest["status"] = validation_status["status"]
        if error_message:
            manifest["error_message"] = error_message
        write_step_manifest(temp_dir, manifest)

        write_step_summary(
            temp_dir,
            step_num=7,
            description=(
                f"{spec.summary_description} The persisted validation files summarize "
                "the integrity state of the canonical project outputs at the time Phase 7 ran."
            ),
            key_files=[
                {
                    "path": "validation_status.json",
                    "description": "Structured validation result with missing/schema/warning lists. Inspect this first.",
                },
                {
                    "path": "validation_summary.txt",
                    "description": "Human-readable validation summary and next action.",
                },
                {
                    "path": "validation_recovery_playbook.md",
                    "description": "Phase-by-phase recovery guide with safest rerun commands.",
                },
                {
                    "path": "validation_recovery_plan.json",
                    "description": "Structured recovery plan for tooling or downstream dashboards.",
                },
            ],
            next_step_reads=[
                {
                    "path": "step_index.md",
                    "description": "Project-level read order reflecting the latest validation state.",
                }
            ],
            warnings=validation_status["warnings"],
        )

    return project_root / spec.folder_name


def _step_status_map(project_root: Path) -> Dict[str, str]:
    statuses: Dict[str, str] = {}
    for step_num, spec in STEP_SPECS.items():
        manifest_path = project_root / spec.folder_name / "step_manifest.json"
        if manifest_path.exists():
            try:
                statuses[f"step{step_num}"] = _read_json(manifest_path).get("status", "complete")
            except (json.JSONDecodeError, OSError):
                statuses[f"step{step_num}"] = "error"
        elif (project_root / spec.folder_name).exists():
            statuses[f"step{step_num}"] = "incomplete"
        else:
            statuses[f"step{step_num}"] = "not_generated"
    return statuses


def clear_step_views(
    config_path: Union[Path, str],
    *,
    repo_root: Optional[Union[Path, str]] = None,
    step_numbers: Optional[Sequence[int]] = None,
    include_root_files: bool = False,
) -> List[Path]:
    """Remove selected derived step folders without touching canonical outputs."""

    _, _, _, project_root = _collect_base_context(
        config_path,
        repo_root=repo_root,
        create_project_root=False,
    )
    targets = sorted(
        {
            int(step_num)
            for step_num in (step_numbers or STEP_SPECS.keys())
            if int(step_num) in STEP_SPECS
        }
    )
    removed_paths: List[Path] = []

    for step_num in targets:
        step_dir = project_root / STEP_SPECS[step_num].folder_name
        if not step_dir.exists():
            continue
        if step_dir.is_dir():
            shutil.rmtree(step_dir)
        else:
            step_dir.unlink()
        removed_paths.append(step_dir)

    if include_root_files:
        for name in (
            "current_run_manifest.json",
            "step_index.md",
            "run_status.json",
            "run_overview.md",
            "run_overview.html",
            "report_digest.md",
            "operational_recovery_playbook.md",
            "operational_recovery_plan.json",
        ):
            path = project_root / name
            if not path.exists():
                continue
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
            removed_paths.append(path)

    return removed_paths


def _pyrosetta_raw_run_paths_from_config(config: dict) -> List[str]:
    raw_paths: List[str] = []
    ppi_config = config.get("ppi") or {}
    result_dirs = ppi_config.get("pyrosetta_result_dirs") or {}
    for entries in result_dirs.values():
        if not entries:
            continue
        for entry in entries:
            path = entry.get("path")
            if path and path not in raw_paths:
                raw_paths.append(str(path))
    return raw_paths


def build_current_run_manifest(
    config_path: Union[Path, str],
    *,
    repo_root: Optional[Union[Path, str]] = None,
    execution_mode: str = "full",
    fresh_run: bool = False,
    stale_steps: Optional[Sequence[int]] = None,
    ppi_config_paths: Optional[Sequence[str]] = None,
    pyrosetta_raw_run_paths: Optional[Sequence[str]] = None,
) -> dict:
    """Build the root-level run manifest for the current derived view state."""

    config_path, repo_root_path, config, project_root = _collect_base_context(
        config_path,
        repo_root=repo_root,
        create_project_root=False,
    )
    raw_run_paths = list(pyrosetta_raw_run_paths or _pyrosetta_raw_run_paths_from_config(config))
    display_raw_paths = [
        _display_path(Path(path), repo_root_path) if Path(path).is_absolute() else Path(path).as_posix()
        for path in raw_run_paths
    ]
    step_status = _step_status_map(project_root)
    stale_step_numbers = sorted(
        {
            int(step_num)
            for step_num in (stale_steps or [])
            if int(step_num) in STEP_SPECS
        }
    )
    for step_num in stale_step_numbers:
        step_status[f"step{step_num}"] = "stale"

    return {
        "project_name": config.get("project_name", project_root.name),
        "generated_at": utc_now_iso(),
        "project_root": _display_path(project_root, repo_root_path),
        "receptors": _config_receptor_ids(config),
        "ligands": _config_ligand_ids(config),
        "vina_config_path": _display_path(config_path, repo_root_path),
        "ppi_config_paths": list(ppi_config_paths or []),
        "pyrosetta_raw_run_paths": display_raw_paths,
        "step_status": step_status,
        "fresh_run": bool(fresh_run),
        "stale_steps": stale_step_numbers,
        "execution_mode": execution_mode,
    }


def _load_current_run_manifest(project_root: Path) -> Optional[dict]:
    path = project_root / "current_run_manifest.json"
    if not path.exists():
        return None
    try:
        return _read_json(path)
    except (json.JSONDecodeError, OSError):
        return None


def _load_run_status(project_root: Path) -> Optional[dict]:
    path = project_root / "run_status.json"
    if not path.exists():
        return None
    try:
        return _read_json(path)
    except (json.JSONDecodeError, OSError):
        return None


def write_run_status(project_root: Union[Path, str], status_data: dict) -> Path:
    """Write ``run_status.json`` at the project root."""

    return _atomic_write_json(_as_path(project_root) / "run_status.json", status_data)


def _step_primary_paths(spec: StepSpec) -> str:
    return ", ".join(f"`{spec.folder_name}/{name}`" for name in spec.primary_files) or "-"


def _load_csv_rows(path: Path) -> List[dict]:
    if not path.exists():
        return []
    try:
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
        with open(path, newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter=delimiter))
    except OSError:
        return []



def _status_label(status: str) -> str:
    return str(status or "unknown").replace("_", " ").strip() or "unknown"


def _load_step_manifest(project_root: Path, step_num: int) -> Optional[dict]:
    path = project_root / STEP_SPECS[step_num].folder_name / "step_manifest.json"
    if not path.exists():
        return None
    try:
        return _read_json(path)
    except (json.JSONDecodeError, OSError):
        return None


def _operational_category_meta(category: str) -> dict:
    return dict(
        OPERATIONAL_CATEGORY_METADATA.get(
            category,
            {
                "label": "Operational follow-up",
                "priority_rank": 2,
                "action_type": "manual_review",
            },
        )
    )


def _operational_action_label(action_type: str) -> str:
    return OPERATIONAL_ACTION_LABELS.get(action_type, "Manual review only")


def _operational_command_for_step(step_number: int, action_type: str) -> str:
    if action_type == "refresh_step_view":
        return f"python run_production.py --only {step_number}"
    if action_type in {"rerun_upstream", "manual_then_rerun"}:
        return f"python run_production.py --from {step_number}"
    return ""


def _operational_checklist_for_issue(
    *,
    category: str,
    step_number: int,
    artifact_name: str,
    recommended_command: str,
) -> List[str]:
    phase_label = _phase_title(step_number)
    artifact_label = f"`{artifact_name}`" if artifact_name else "the affected step artifact"
    if category == "step_not_generated":
        return [
            f"Confirm whether {phase_label} was skipped entirely or whether the derived step folder was never refreshed.",
            f"Inspect canonical outputs that feed {phase_label} before rerunning.",
            f"Use `{recommended_command}` to rebuild this phase and its downstream views.",
        ]
    if category == "missing_raw_pose":
        return [
            f"Inspect the receptor/ligand pair behind {artifact_label} and confirm the canonical `.pdbqt` file exists.",
            "Check whether the docking job failed upstream or whether the file path changed after the run.",
            f"Use `{recommended_command}` to refresh raw poses and downstream interpretation.",
        ]
    if category == "missing_ppi_ranking":
        return [
            f"Open `step2_ppi_raw/raw_run_paths.tsv` and confirm which target is missing {artifact_label}.",
            "Check whether the raw PyRosetta run finished and whether `final_result/final_ranking.csv` exists.",
            f"Use `{recommended_command}` to rebuild Phase 2 and downstream PPI interpretation.",
        ]
    if category == "missing_ppi_evidence":
        return [
            f"Confirm whether {artifact_label} is absent from the canonical project root or only from the derived step view.",
            "Review whether Phase 3 postprocess completed and whether Step 2 supplied the required raw ranking inputs.",
            f"Use `{recommended_command}` to rebuild Phase 3 and downstream consumers.",
        ]
    if category == "stale_step_view":
        return [
            f"Compare the timestamp of {phase_label} against `current_run_manifest.json` to confirm the stale state.",
            "Decide whether only the derived step view is stale or whether the upstream canonical outputs changed too.",
            f"Use `{recommended_command}` to refresh this step view before mixing interpretations.",
        ]
    if category == "metadata_review":
        return [
            "Check whether the referenced PyRosetta metadata JSON is missing, unreadable, or malformed.",
            f"Confirm that the raw run directory recorded for {phase_label} still exists and is readable.",
            f"After fixing the metadata issue, use `{recommended_command}` to rebuild the affected scope.",
        ]
    if category == "fallback_source":
        return [
            f"Confirm why {artifact_label} was loaded from a fallback canonical location instead of the project root.",
            "Decide whether the fallback is intentional for this run or whether outputs from another location are being mixed in.",
            "Document the reason before relying on downstream verdict/report interpretation.",
        ]
    if category == "optional_supporting_artifact":
        return [
            f"Decide whether {artifact_label} is actually needed for this review or only for deeper provenance/debugging.",
            "If it should exist, restore or regenerate it before sharing the run externally.",
            "If it is optional for the current interpretation, document the omission and continue carefully.",
        ]
    return [
        "Inspect the step manifest and the referenced artifact in context.",
        "Decide whether the issue requires a rerun, a step refresh, or manual documentation.",
        "Reopen the overview after the issue is addressed.",
    ]


def _build_operational_issue(
    *,
    step_number: int,
    category: str,
    message: str,
    artifact_name: str = "",
) -> dict:
    meta = _operational_category_meta(category)
    action_type = str(meta["action_type"])
    recommended_command = _operational_command_for_step(step_number, action_type)
    return {
        "step_number": step_number,
        "step_label": _phase_title(step_number),
        "category": category,
        "category_label": str(meta["label"]),
        "message": message,
        "artifact_name": artifact_name,
        "priority_rank": int(meta["priority_rank"]),
        "priority_label": _validation_priority_label(int(meta["priority_rank"])),
        "action_type": action_type,
        "action_label": _operational_action_label(action_type),
        "recommended_command": recommended_command,
        "requires_manual_review": action_type in {"manual_then_rerun", "manual_review"},
        "checklist": _operational_checklist_for_issue(
            category=category,
            step_number=step_number,
            artifact_name=artifact_name,
            recommended_command=recommended_command,
        ),
    }


def _operational_issues_for_step(
    step_number: int,
    step_status: str,
    manifest: Optional[dict],
) -> List[dict]:
    issues: List[dict] = []
    if step_status in {"stale"}:
        issues.append(
            _build_operational_issue(
                step_number=step_number,
                category="stale_step_view",
                message=f"{_phase_title(step_number)} is marked stale in the current manifest.",
            )
        )
    if manifest is None:
        if step_status in {"not_generated", "incomplete"}:
            issues.append(
                _build_operational_issue(
                    step_number=step_number,
                    category="step_not_generated",
                    message=f"{_phase_title(step_number)} has not generated its derived step outputs yet.",
                )
            )
        return issues

    missing_files = list(manifest.get("missing_files") or [])
    warnings = [str(item) for item in (manifest.get("warnings") or [])]

    if step_number == 1:
        for item in missing_files:
            issues.append(
                _build_operational_issue(
                    step_number=1,
                    category="missing_raw_pose",
                    message=f"Missing raw pose file: {item}",
                    artifact_name=str(item),
                )
            )
    elif step_number == 2:
        for item in missing_files:
            issues.append(
                _build_operational_issue(
                    step_number=2,
                    category="missing_ppi_ranking",
                    message=f"Missing PyRosetta final ranking: {item}",
                    artifact_name=str(item),
                )
            )
        for warning in warnings:
            if warning.startswith("Optional metadata JSON missing for ") or warning.startswith(
                "Unable to read metadata JSON for "
            ):
                issues.append(
                    _build_operational_issue(
                        step_number=2,
                        category="metadata_review",
                        message=warning,
                    )
                )
    elif step_number == 3:
        for item in missing_files:
            issues.append(
                _build_operational_issue(
                    step_number=3,
                    category="missing_ppi_evidence",
                    message=f"Missing PPI evidence artifact: {item}",
                    artifact_name=str(item),
                )
            )
        for warning in warnings:
            if warning.startswith("Using fallback canonical source for "):
                artifact_name = warning.split("Using fallback canonical source for ", 1)[1].split(":", 1)[0].strip()
                issues.append(
                    _build_operational_issue(
                        step_number=3,
                        category="fallback_source",
                        message=warning,
                        artifact_name=artifact_name,
                    )
                )
            elif warning.startswith("Optional artifact missing: ") or warning == "Optional Phase 1 interface report is not available.":
                artifact_name = warning.split(": ", 1)[1].strip() if ": " in warning else "phase1_interface_report.md"
                issues.append(
                    _build_operational_issue(
                        step_number=3,
                        category="optional_supporting_artifact",
                        message=warning,
                        artifact_name=artifact_name,
                    )
                )
    return issues


def _sorted_operational_issues(issues: Sequence[dict]) -> List[dict]:
    return sorted(
        issues,
        key=lambda issue: (
            int(issue.get("priority_rank", 2)),
            int(issue.get("step_number", 99)),
            str(issue.get("category_label", "")),
            str(issue.get("message", "")),
        ),
    )


def _operational_action_groups(issues: Sequence[dict]) -> List[dict]:
    groups: Dict[Tuple[int, str, str], dict] = {}
    for issue in _sorted_operational_issues(issues):
        key = (
            int(issue.get("step_number", 99)),
            str(issue.get("category", "")),
            str(issue.get("recommended_command", "")),
        )
        group = groups.get(key)
        if group is None:
            group = {
                "step_number": issue["step_number"],
                "step_label": issue["step_label"],
                "category": issue["category"],
                "category_label": issue["category_label"],
                "group_key": _recovery_group_key(
                    "operational",
                    issue.get("step_number"),
                    issue.get("category", ""),
                    issue.get("recommended_command", ""),
                ),
                "priority_rank": issue["priority_rank"],
                "priority_label": issue["priority_label"],
                "action_type": issue["action_type"],
                "action_label": issue["action_label"],
                "recommended_command": issue["recommended_command"],
                "manual_review_required": issue["requires_manual_review"],
                "checklist": list(issue.get("checklist") or []),
                "issue_count": 0,
                "artifacts": [],
                "findings": [],
            }
            groups[key] = group
        group["issue_count"] += 1
        artifact_name = str(issue.get("artifact_name", "")).strip()
        if artifact_name and artifact_name not in group["artifacts"]:
            group["artifacts"].append(artifact_name)
        group["findings"].append({"message": issue["message"]})
    return sorted(
        groups.values(),
        key=lambda group: (
            int(group.get("priority_rank", 2)),
            int(group.get("step_number", 99)),
            str(group.get("category_label", "")),
        ),
    )


def build_operational_recovery_plan(
    project_root: Union[Path, str],
    *,
    current_manifest: Optional[dict] = None,
) -> dict:
    project_root = _as_path(project_root)
    current_manifest = current_manifest or _load_current_run_manifest(project_root) or {}
    step_status = current_manifest.get("step_status") or _step_status_map(project_root)
    issues: List[dict] = []
    for step_number in (1, 2, 3):
        issues.extend(
            _operational_issues_for_step(
                step_number,
                str(step_status.get(f"step{step_number}", "not_generated")),
                _load_step_manifest(project_root, step_number),
            )
        )
    issues = _sorted_operational_issues(issues)
    action_groups = _operational_action_groups(issues)

    rerun_groups = [group for group in action_groups if group.get("recommended_command")]
    primary_group = rerun_groups[0] if rerun_groups else None
    summary = {
        "issue_count": len(issues),
        "action_group_count": len(action_groups),
        "manual_review_count": len([issue for issue in issues if issue.get("requires_manual_review")]),
        "blocking_count": len([issue for issue in issues if int(issue.get("priority_rank", 2)) == 0]),
        "stale_count": len([issue for issue in issues if issue.get("category") == "stale_step_view"]),
    }
    triage = {
        "rerun_upstream": len([group for group in action_groups if group.get("action_type") == "rerun_upstream"]),
        "refresh_step_view": len([group for group in action_groups if group.get("action_type") == "refresh_step_view"]),
        "manual_then_rerun": len([group for group in action_groups if group.get("action_type") == "manual_then_rerun"]),
        "manual_review": len([group for group in action_groups if group.get("action_type") == "manual_review"]),
    }

    return {
        "generated_at": utc_now_iso(),
        "project_root": _display_path(project_root, project_root.parent.parent),
        "issues": issues,
        "action_groups": action_groups,
        "recommended_step_number": primary_group.get("step_number") if primary_group else None,
        "recommended_step_label": primary_group.get("step_label", "No rerun required") if primary_group else "No rerun required",
        "recommended_command": primary_group.get("recommended_command", "") if primary_group else "",
        "recommended_reason": (
            f"Address {primary_group['category_label']} in {primary_group['step_label']} first."
            if primary_group
            else "No upstream operational blockers are currently recorded for Steps 1-3."
        ),
        "triage": triage,
        "summary": summary,
    }


def write_operational_recovery_outputs(
    project_root: Union[Path, str],
    recovery_plan: dict,
) -> Tuple[Path, Path]:
    project_root = _as_path(project_root)
    json_path = _atomic_write_json(project_root / "operational_recovery_plan.json", recovery_plan)

    lines = [
        "# Operational Recovery Playbook",
        "",
        "## Status",
        f"- Action groups: `{recovery_plan['summary']['action_group_count']}`",
        f"- Findings: `{recovery_plan['summary']['issue_count']}`",
        f"- Blocking findings: `{recovery_plan['summary']['blocking_count']}`",
        f"- Manual review findings: `{recovery_plan['summary']['manual_review_count']}`",
        f"- Stale step findings: `{recovery_plan['summary']['stale_count']}`",
    ]
    if recovery_plan.get("recommended_command"):
        lines.append(f"- Safest repair command: `{recovery_plan['recommended_command']}`")
        lines.append(f"- Earliest affected step: `{recovery_plan['recommended_step_label']}`")
    lines.extend(
        [
            "",
            "## Triage Summary",
            f"- `{recovery_plan['triage']['rerun_upstream']}` group(s) need an upstream rerun.",
            f"- `{recovery_plan['triage']['refresh_step_view']}` group(s) only need a derived-step refresh.",
            f"- `{recovery_plan['triage']['manual_then_rerun']}` group(s) need manual review before rerun.",
            f"- `{recovery_plan['triage']['manual_review']}` group(s) are manual review only.",
            "",
            "## Why This Command",
            f"- {recovery_plan['recommended_reason']}",
            "",
            "## Action Groups",
        ]
    )
    if recovery_plan["action_groups"]:
        for index, group in enumerate(recovery_plan["action_groups"], start=1):
            lines.extend(
                [
                    f"<a id=\"{group['group_key']}\"></a>",
                    f"### Action Group {index}: {group['category_label']}",
                    f"- Deep link: `{_playbook_group_link('operational_recovery_playbook.md', group['group_key'])}`",
                    f"- Priority: `{group['priority_label']}`",
                    f"- Scope: `{group['step_label']}`",
                    f"- Action path: `{group['action_label']}`",
                    (
                        f"- Suggested command: `{group['recommended_command']}`"
                        if group.get("recommended_command")
                        else "- Suggested command: manual review only"
                    ),
                    f"- Finding count: `{group['issue_count']}`",
                ]
            )
            if group["artifacts"]:
                lines.append("- Related artifacts: " + ", ".join(f"`{item}`" for item in group["artifacts"]))
            lines.append("- Checklist:")
            lines.extend(f"  - {item}" for item in group["checklist"])
            lines.append("- Findings:")
            lines.extend(f"  - {finding['message']}" for finding in group["findings"])
            lines.append("")
    else:
        lines.extend(["- No upstream operational blockers were recorded for Steps 1-3.", ""])

    markdown_path = _atomic_write_text(
        project_root / "operational_recovery_playbook.md",
        "\n".join(lines).rstrip() + "\n",
    )
    return json_path, markdown_path


def _progress_summary(phase_states: Sequence[dict]) -> Tuple[int, int, int]:
    total = len(phase_states)
    resolved = sum(
        1
        for entry in phase_states
        if str(entry.get("status", "")) in {"completed", "skipped", "failed"}
    )
    percent = int((resolved / total) * 100) if total else 0
    return resolved, total, percent


def _step_health_summary(current_manifest: Optional[dict], project_root: Path) -> Dict[str, object]:
    step_status = (current_manifest or {}).get("step_status") or _step_status_map(project_root)
    counts: Counter = Counter(step_status.values())
    return {
        "statuses": step_status,
        "counts": dict(counts),
        "complete_count": int(counts.get("complete", 0)),
        "partial_count": int(counts.get("partial", 0)),
        "stale_count": int(counts.get("stale", 0)),
        "failed_count": int(counts.get("failed", 0) + counts.get("error", 0)),
        "not_generated_count": int(counts.get("not_generated", 0)),
    }


def _normalize_step_numbers(step_numbers: Optional[Iterable[object]]) -> List[int]:
    normalized: List[int] = []
    for item in step_numbers or []:
        step_number = safe_int(item, None)
        if step_number is None or step_number not in STEP_SPECS or step_number in normalized:
            continue
        normalized.append(step_number)
    return normalized


def _step_context_label(step_numbers: Optional[Iterable[object]]) -> str:
    normalized = _normalize_step_numbers(step_numbers)
    if not normalized:
        return ""
    if len(normalized) == 1:
        return f"Step {normalized[0]}"
    return "Steps " + ", ".join(str(step_number) for step_number in normalized)


def _step_summary_path(step_number: int) -> str:
    spec = _step_spec(step_number)
    return f"{spec.folder_name}/summary.md"


def _step_summary_links(step_numbers: Optional[Iterable[object]]) -> List[dict]:
    links: List[dict] = []
    for step_number in _normalize_step_numbers(step_numbers):
        spec = _step_spec(step_number)
        links.append(
            {
                "step_number": step_number,
                "label": f"Step {step_number}: {spec.step_name.replace('_', ' ').title()}",
                "path": _step_summary_path(step_number),
            }
        )
    return links


def _with_step_summary_links(payload: Mapping[str, object]) -> dict:
    enriched = dict(payload)
    enriched["step_numbers"] = _normalize_step_numbers(payload.get("step_numbers"))
    enriched["step_summary_links"] = _step_summary_links(enriched.get("step_numbers"))
    return enriched


def _step_summary_links_text(step_summary_links: Sequence[Mapping[str, object]]) -> str:
    return ", ".join(f"`{str(item.get('path', ''))}`" for item in step_summary_links if str(item.get("path", "")).strip())


def _overview_group_link(item: Mapping[str, object]) -> str:
    group_key = str(item.get("group_key", "")).strip()
    if not group_key:
        return "run_overview.html#recovery-radar"
    params = urlencode(
        {
            "section": "recovery-radar",
            "scope": str(item.get("source_kind", "all")),
            "priority": str(item.get("priority_label", "all")),
            "group": group_key,
        }
    )
    return f"run_overview.html#{params}"


def _overview_step_link(step_number: int) -> str:
    section = "operational-readiness" if int(step_number) <= 3 else "result-highlights"
    params = urlencode({"section": section, "step": int(step_number)})
    return f"run_overview.html#{params}"


def _step_triage_items(recovery_radar: Mapping[str, object]) -> Dict[int, List[dict]]:
    triage: Dict[int, List[dict]] = {step_num: [] for step_num in STEP_SPECS}
    radar_items = recovery_radar.get("all_items") or recovery_radar.get("items") or []
    for index, item in enumerate(radar_items, start=1):
        triage_entry = {
            "group_key": str(item.get("group_key", "")),
            "title": str(item.get("title", "Recovery group")),
            "source_kind": str(item.get("source_kind", "general")),
            "source_label": str(item.get("source_label", "General")),
            "priority_rank": safe_int(item.get("priority_rank"), 3),
            "priority_label": str(item.get("priority_label", "low")),
            "scope_label": str(item.get("scope_label", "")),
            "recommended_command": str(item.get("recommended_command", "")).strip(),
            "playbook_link": str(item.get("playbook_link", item.get("playbook_path", ""))),
            "overview_link": _overview_group_link(item),
            "step_summary_links": list(item.get("step_summary_links") or []),
            "default_order": index,
        }
        for step_number in _normalize_step_numbers(item.get("step_numbers")):
            triage.setdefault(step_number, []).append(dict(triage_entry))
    for step_number in list(triage):
        triage[step_number] = sorted(
            triage[step_number],
            key=lambda entry: (
                safe_int(entry.get("priority_rank"), 3),
                safe_int(entry.get("default_order"), 99),
                str(entry.get("title", "")),
            ),
        )
    return triage


def _step_triage_summary(triage_items: Sequence[Mapping[str, object]]) -> str:
    if not triage_items:
        return "`clear`"
    first_item = triage_items[0]
    summary = f"`{str(first_item.get('overview_link', 'run_overview.html#recovery-radar'))}`"
    overflow = len(triage_items) - 1
    if overflow > 0:
        summary += f" (+{overflow} more)"
    return summary


def _step1_input_summary(project_root: Path) -> dict:
    path = project_root / STEP_SPECS[1].folder_name / "raw_pose_index.csv"
    rows = _load_csv_rows(path)
    manifest = _load_step_manifest(project_root, 1)
    if not rows:
        return {
            "title": "Step 1 Raw Pose Health",
            "available": False,
            "path": f"{STEP_SPECS[1].folder_name}/raw_pose_index.csv",
            "headline": "Raw Vina pose coverage is not indexed yet.",
            "details": [],
            "step_numbers": [1],
        }

    expected = len(rows)
    found = sum(1 for row in rows if str(row.get("raw_pose_file", "")).strip())
    missing_pairs = [
        f"{row.get('receptor_id', '?')}/{row.get('ligand_id', '?')}"
        for row in rows
        if not str(row.get("raw_pose_file", "")).strip()
    ]
    model_counts = [safe_int(row.get("n_models"), 0) for row in rows if safe_int(row.get("n_models"), 0) > 0]
    return {
        "title": "Step 1 Raw Pose Health",
        "available": True,
        "path": f"{STEP_SPECS[1].folder_name}/raw_pose_index.csv",
        "headline": f"Indexed raw poses for {found}/{expected} receptor-ligand pair(s).",
        "details": [
            (
                f"Missing raw pose pairs: {', '.join(missing_pairs[:3])}."
                + (f" (+{len(missing_pairs) - 3} more)" if len(missing_pairs) > 3 else "")
                if missing_pairs
                else "All expected raw pose pairs are present."
            ),
            (
                f"Observed pose model count range: {min(model_counts)}-{max(model_counts)}."
                if model_counts
                else "No pose model counts were recorded."
            ),
            (
                f"Step status: {manifest.get('status', 'unknown')}."
                if manifest
                else "Step manifest is not available."
            ),
        ],
        "step_numbers": [1],
    }


def _step2_ppi_summary(project_root: Path) -> dict:
    path = project_root / STEP_SPECS[2].folder_name / "raw_run_paths.tsv"
    rows = _load_csv_rows(path)
    manifest = _load_step_manifest(project_root, 2)
    if not rows:
        return {
            "title": "Step 2 PPI Run Health",
            "available": False,
            "path": f"{STEP_SPECS[2].folder_name}/raw_run_paths.tsv",
            "headline": "Raw PyRosetta run coverage is not indexed yet.",
            "details": [],
            "step_numbers": [2],
        }

    ranking_ready = [row for row in rows if str(row.get("final_ranking_csv", "")).strip()]
    metadata_ready = [row for row in rows if str(row.get("metadata_json", "")).strip()]
    missing_targets = [
        str(row.get("target_name", "?"))
        for row in rows
        if not str(row.get("final_ranking_csv", "")).strip()
    ]
    return {
        "title": "Step 2 PPI Run Health",
        "available": True,
        "path": f"{STEP_SPECS[2].folder_name}/raw_run_paths.tsv",
        "headline": f"PyRosetta ranking summaries are available for {len(ranking_ready)}/{len(rows)} target(s).",
        "details": [
            (
                f"Targets missing final ranking files: {', '.join(missing_targets[:3])}."
                + (f" (+{len(missing_targets) - 3} more)" if len(missing_targets) > 3 else "")
                if missing_targets
                else "All indexed raw runs include a final ranking file."
            ),
            f"Targets with metadata JSON: {len(metadata_ready)}/{len(rows)}.",
            (
                f"Step status: {manifest.get('status', 'unknown')}."
                if manifest
                else "Step manifest is not available."
            ),
        ],
        "step_numbers": [2],
    }


def _step3_evidence_summary(project_root: Path) -> dict:
    path = project_root / STEP_SPECS[3].folder_name / "ppi_pyrosetta_residues.csv"
    rows = _load_csv_rows(path)
    summary_rows = _load_csv_rows(project_root / STEP_SPECS[3].folder_name / "ppi_pyrosetta_summary.csv")
    manifest = _load_step_manifest(project_root, 3)
    if not rows:
        return {
            "title": "Step 3 PPI Evidence",
            "available": False,
            "path": f"{STEP_SPECS[3].folder_name}/ppi_pyrosetta_residues.csv",
            "headline": "Aggregated PPI residue evidence is not available yet.",
            "details": [],
            "step_numbers": [3],
        }

    def _priority(row: dict) -> Tuple[float, float, float]:
        return (
            safe_float(row.get("frac_runs_supporting"), 0.0),
            safe_float(row.get("occupancy"), 0.0),
            safe_float(row.get("n_runs_supporting"), 0.0),
        )

    ranked = sorted(rows, key=_priority, reverse=True)
    top_residues = [
        f"{row.get('receptor_id', '?')}:{row.get('residue_id', '?')} occ={row.get('occupancy', 'n/a')}"
        for row in ranked[:3]
    ]
    receptors = sorted({str(row.get("receptor_id", "")).strip() for row in rows if row.get("receptor_id")})
    return {
        "title": "Step 3 PPI Evidence",
        "available": True,
        "path": f"{STEP_SPECS[3].folder_name}/ppi_pyrosetta_residues.csv",
        "headline": f"PPI residue evidence contains {len(rows)} row(s) across {len(receptors)} receptor state(s).",
        "details": [
            "Top residues: " + ", ".join(top_residues) + "." if top_residues else "Top residues are not available.",
            f"Summary rows available: {len(summary_rows)}.",
            (
                f"Step status: {manifest.get('status', 'unknown')}."
                if manifest
                else "Step manifest is not available."
            ),
        ],
        "step_numbers": [3],
    }


def _pocket_summary(project_root: Path) -> dict:
    rows = _load_csv_rows(project_root / "vina_pocket_table.csv")
    if not rows:
        return {
            "available": False,
            "path": "step4_vina_postprocess/vina_pocket_table.csv",
            "headline": "Pocket summary is not available yet.",
            "details": [],
            "step_numbers": [4],
        }

    receptors = sorted({str(row.get("receptor_id", "")).strip() for row in rows if row.get("receptor_id")})
    best_row = min(rows, key=lambda row: safe_float(row.get("best_affinity"), default=9999.0))
    multi_ligand = sum(1 for row in rows if safe_int(row.get("n_ligand"), 0) >= 2)
    coverage_rows = _load_csv_rows(project_root / "vina_postprocess_coverage.csv")
    parsed_pairs = sum(1 for row in coverage_rows if str(row.get("status", "")).strip().lower() == "parsed")
    comparison_rows = _load_csv_rows(project_root / "vina_pocket_comparison.csv")
    same_patch_candidates = sum(
        1
        for row in comparison_rows
        if str(row.get("same_patch_candidate", "")).strip().lower() == "true"
    )
    top_residues = str(best_row.get("top_residues", "")).strip()
    return {
        "available": True,
        "path": "step4_vina_postprocess/vina_pocket_table.csv",
        "headline": f"{len(rows)} pockets summarized across {len(receptors)} receptor state(s).",
        "details": [
            (
                "Best affinity pocket: "
                f"{best_row.get('receptor_id', '?')}/{best_row.get('pocket_id', '?')} "
                f"({best_row.get('best_affinity', 'n/a')} kcal/mol)."
            ),
            f"Multi-ligand pockets: {multi_ligand}.",
            (
                f"Parsed receptor-ligand pairs: {parsed_pairs}/{len(coverage_rows)}."
                if coverage_rows
                else "Coverage table is not available."
            ),
            (
                f"Cross-receptor same-patch candidates: {same_patch_candidates}."
                if comparison_rows
                else "Cross-receptor comparison table is not available."
            ),
            (f"Top residues for the strongest pocket: {top_residues}." if top_residues else "Top residues are not available."),
        ],
        "step_numbers": [4],
    }


def _verdict_summary(project_root: Path) -> dict:
    rows = _load_csv_rows(project_root / "valid_sites.csv")
    if not rows:
        return {
            "available": False,
            "path": "step5_verdict/valid_sites.csv",
            "headline": "Verdict summary is not available yet.",
            "details": [],
            "step_numbers": [5],
        }

    verdict_counts: Counter = Counter(
        str(row.get("verdict", "unknown")).strip().upper() or "UNKNOWN"
        for row in rows
    )
    ranked = sorted(
        rows,
        key=lambda row: (
            -safe_float(row.get("confidence_score"), default=-1.0),
            safe_float(row.get("best_affinity"), default=9999.0),
        ),
    )
    consensus_sites = {
        str(row.get("consensus_site_id", "")).strip()
        for row in rows
        if str(row.get("consensus_site_id", "")).strip()
    }
    top_row = ranked[0] if ranked else {}
    top_labels = [
        (
            f"{row.get('receptor_id', '?')}/{row.get('pocket_id', '?')}"
            f" [{row.get('verdict', 'UNKNOWN')}, score={row.get('confidence_score', 'n/a')}]"
        )
        for row in ranked[:3]
    ]
    return {
        "available": True,
        "path": "step5_verdict/valid_sites.csv",
        "headline": (
            f"{len(rows)} verdict row(s). "
            + ", ".join(f"{label}={count}" for label, count in sorted(verdict_counts.items()))
        ),
        "details": [
            "Top ranked sites: " + ", ".join(top_labels) + "." if top_labels else "Top ranked sites are not available.",
            f"Consensus site groups: {len(consensus_sites)}." if consensus_sites else "Consensus site groups were not assigned.",
            (
                f"Top site reasons: {top_row.get('reasons', '')}."
                if str(top_row.get("reasons", "")).strip()
                else "Top site reasons are not available."
            ),
        ],
        "step_numbers": [5],
    }


def _validation_summary_for_overview(project_root: Path) -> dict:
    path = project_root / STEP_SPECS[7].folder_name / "validation_status.json"
    if not path.exists():
        return {
            "available": False,
            "path": f"{STEP_SPECS[7].folder_name}/validation_summary.txt",
            "headline": "Validation has not been recorded yet.",
            "details": [],
            "status": "not_run",
            "recommended_command": "",
            "validation_only_command": "python run_production.py --only 7",
            "playbook_path": f"{STEP_SPECS[7].folder_name}/validation_recovery_playbook.md",
            "step_numbers": [7],
        }

    payload = _read_json(path)
    recovery_plan_path = project_root / STEP_SPECS[7].folder_name / "validation_recovery_plan.json"
    recovery_plan = _read_json(recovery_plan_path) if recovery_plan_path.exists() else {}
    status = str(payload.get("status", "unknown"))
    warnings = safe_int(payload.get("warning_count"), 0)
    failures = safe_int(payload.get("failure_count"), 0)
    details: List[str] = []
    if payload.get("missing_files"):
        details.append(f"Missing artifacts: {len(payload['missing_files'])}.")
    if payload.get("schema_errors"):
        details.append(f"Schema mismatches: {len(payload['schema_errors'])}.")
    if payload.get("error_message"):
        details.append(f"Last error: {payload['error_message']}.")
    if recovery_plan.get("recommended_command"):
        details.append(f"Safest repair command: {recovery_plan['recommended_command']}.")
    if recovery_plan.get("recommended_phase_number") in STEP_SPECS:
        details.append(f"Earliest affected phase: {recovery_plan['recommended_phase_label']}.")
    manual_review_count = safe_int((recovery_plan.get("summary") or {}).get("manual_review_count"), 0)
    action_group_count = safe_int((recovery_plan.get("summary") or {}).get("action_group_count"), 0)
    if manual_review_count:
        details.append(f"Manual review findings: {manual_review_count}.")
    if action_group_count:
        details.append(f"Recovery checklist groups: {action_group_count}.")
    if recovery_plan_path.exists():
        details.append("Recovery playbook: step7_validate/validation_recovery_playbook.md.")
    return {
        "available": True,
        "path": f"{STEP_SPECS[7].folder_name}/validation_summary.txt",
        "headline": f"Validation status: {status}. {warnings} warning(s), {failures} failure(s).",
        "details": details,
        "status": status,
        "recommended_command": str(recovery_plan.get("recommended_command", "")),
        "validation_only_command": str(
            recovery_plan.get("validation_only_command", "python run_production.py --only 7")
        ),
        "playbook_path": f"{STEP_SPECS[7].folder_name}/validation_recovery_playbook.md",
        "manual_review_count": manual_review_count,
        "action_group_count": action_group_count,
        "manual_review_required": bool(recovery_plan.get("manual_review_required", False)),
        "step_numbers": [7],
    }


def _load_validation_recovery_plan(project_root: Path) -> Optional[dict]:
    path = project_root / STEP_SPECS[7].folder_name / "validation_recovery_plan.json"
    if not path.exists():
        return None
    try:
        return _read_json(path)
    except (json.JSONDecodeError, OSError):
        return None


def _load_operational_recovery_plan(project_root: Path) -> Optional[dict]:
    path = project_root / "operational_recovery_plan.json"
    if not path.exists():
        return None
    try:
        return _read_json(path)
    except (json.JSONDecodeError, OSError):
        return None


def _operational_recovery_summary_for_overview(
    project_root: Path,
    *,
    current_manifest: Optional[dict] = None,
    recovery_plan: Optional[dict] = None,
) -> dict:
    plan = recovery_plan or build_operational_recovery_plan(project_root, current_manifest=current_manifest)
    issue_count = safe_int((plan.get("summary") or {}).get("issue_count"), 0)
    action_group_count = safe_int((plan.get("summary") or {}).get("action_group_count"), 0)
    manual_review_count = safe_int((plan.get("summary") or {}).get("manual_review_count"), 0)
    if issue_count == 0:
        return {
            "title": "Operational Recovery",
            "available": True,
            "path": "operational_recovery_playbook.md",
            "headline": "No upstream operational blockers are recorded for Steps 1-3.",
            "details": [
                "Raw pose coverage, raw PPI runs, and PPI residue evidence do not currently need recovery actions.",
                "Open operational_recovery_playbook.md for the checklist format, even when there are no blockers.",
            ],
            "issue_count": 0,
            "action_group_count": 0,
            "manual_review_count": 0,
            "recommended_command": "",
            "step_numbers": [1, 2, 3],
        }

    details: List[str] = []
    if plan.get("recommended_command"):
        details.append(f"Safest repair command: {plan['recommended_command']}.")
    details.append(f"Action groups: {action_group_count}.")
    details.append(f"Manual review findings: {manual_review_count}.")
    details.append("Recovery playbook: operational_recovery_playbook.md.")
    return {
        "title": "Operational Recovery",
        "available": True,
        "path": "operational_recovery_playbook.md",
        "headline": (
            f"{action_group_count} upstream recovery group(s) across Steps 1-3."
            + (f" {manual_review_count} finding(s) require manual review." if manual_review_count else "")
        ),
        "details": details,
        "issue_count": issue_count,
        "action_group_count": action_group_count,
        "manual_review_count": manual_review_count,
        "recommended_command": str(plan.get("recommended_command", "")),
        "step_numbers": [1, 2, 3],
    }


def _recovery_step_numbers(source_kind: str, group: Mapping[str, object]) -> List[int]:
    explicit = _normalize_step_numbers(group.get("step_numbers"))
    if explicit:
        return explicit
    if source_kind == "validation":
        phase_number = safe_int(group.get("phase_number"), None)
        if phase_number is not None and phase_number in STEP_SPECS and phase_number != 7:
            return [phase_number, 7]
        return [7]
    return _normalize_step_numbers([group.get("step_number")])


def _recovery_radar_item(
    *,
    source_kind: str,
    source_label: str,
    group: dict,
    playbook_path: str,
) -> dict:
    findings = list(group.get("findings") or [])
    messages = [str(item.get("message", "")).strip() for item in findings if str(item.get("message", "")).strip()]
    artifacts = [str(item).strip() for item in (group.get("artifacts") or []) if str(item).strip()]
    scope_label = str(group.get("phase_label") or group.get("step_label") or "Manual review")
    step_numbers = _recovery_step_numbers(source_kind, group)
    group_key = str(
        group.get("group_key")
        or _recovery_group_key(
            source_kind,
            group.get("category_label", "Recovery group"),
            scope_label,
            group.get("recommended_command", ""),
        )
    )
    return {
        "source_kind": source_kind,
        "source_label": source_label,
        "group_key": group_key,
        "title": str(group.get("category_label", "Recovery group")),
        "scope_label": scope_label,
        "priority_rank": safe_int(group.get("priority_rank"), 2),
        "priority_label": str(group.get("priority_label", "medium")),
        "action_label": str(group.get("action_label", "Manual review")),
        "recommended_command": str(group.get("recommended_command", "")).strip(),
        "manual_review_required": bool(group.get("manual_review_required")),
        "issue_count": safe_int(group.get("issue_count"), len(messages)),
        "messages": messages[:2],
        "message_overflow_count": max(0, len(messages) - 2),
        "artifacts": artifacts[:2],
        "artifact_overflow_count": max(0, len(artifacts) - 2),
        "playbook_path": playbook_path,
        "playbook_link": _playbook_group_link(playbook_path, group_key),
        "step_numbers": step_numbers,
        "step_summary_links": _step_summary_links(step_numbers),
    }


def _build_recovery_radar(
    *,
    operational_recovery_plan: Optional[dict],
    validation_recovery_plan: Optional[dict],
) -> dict:
    def _priority_counts(rows: Sequence[dict]) -> dict:
        counts = {label: 0 for label in VALIDATION_PRIORITY_LABELS.values()}
        for row in rows:
            label = str(row.get("priority_label") or "low")
            counts[label] = counts.get(label, 0) + 1
        return counts

    def _filter_view(*, filter_key: str, label: str, rows: Sequence[dict]) -> dict:
        counts = _priority_counts(rows)
        manual_review_groups = len([row for row in rows if row.get("manual_review_required")])
        if not rows:
            headline = "No groups currently match this filter."
        else:
            summary_parts = [f"{len(rows)} group(s)"]
            for priority_label in ("immediate", "high", "medium", "low"):
                count = int(counts.get(priority_label, 0))
                if count:
                    summary_parts.append(f"{count} {priority_label}")
            if manual_review_groups and filter_key != "manual":
                summary_parts.append(f"{manual_review_groups} manual review")
            headline = ", ".join(summary_parts) + "."
        return {
            "filter_key": filter_key,
            "label": label,
            "group_count": len(rows),
            "manual_review_groups": manual_review_groups,
            "priority_counts": counts,
            "headline": headline,
        }

    items: List[dict] = []
    for group in (validation_recovery_plan or {}).get("action_groups") or []:
        items.append(
            _recovery_radar_item(
                source_kind="validation",
                source_label="Validation",
                group=group,
                playbook_path=f"{STEP_SPECS[7].folder_name}/validation_recovery_playbook.md",
            )
        )
    for group in (operational_recovery_plan or {}).get("action_groups") or []:
        items.append(
            _recovery_radar_item(
                source_kind="operational",
                source_label="Operational",
                group=group,
                playbook_path="operational_recovery_playbook.md",
            )
        )

    source_rank = {"validation": 0, "operational": 1}
    items = sorted(
        items,
        key=lambda item: (
            int(item.get("priority_rank", 2)),
            source_rank.get(str(item.get("source_kind", "operational")), 9),
            str(item.get("scope_label", "")),
            str(item.get("title", "")),
        ),
    )
    visible_items = items[:6]
    filter_views = [
        _filter_view(filter_key="all", label="All Groups", rows=items),
        _filter_view(
            filter_key="validation",
            label="Validation",
            rows=[item for item in items if item.get("source_kind") == "validation"],
        ),
        _filter_view(
            filter_key="operational",
            label="Operational",
            rows=[item for item in items if item.get("source_kind") == "operational"],
        ),
        _filter_view(
            filter_key="manual",
            label="Manual Review",
            rows=[item for item in items if item.get("manual_review_required")],
        ),
    ]
    summary = {
        "total_groups": len(items),
        "visible_groups": len(visible_items),
        "overflow_count": max(0, len(items) - len(visible_items)),
        "immediate_count": len([item for item in items if int(item.get("priority_rank", 2)) == 0]),
        "manual_review_groups": len([item for item in items if item.get("manual_review_required")]),
        "validation_groups": len([item for item in items if item.get("source_kind") == "validation"]),
        "operational_groups": len([item for item in items if item.get("source_kind") == "operational"]),
        "priority_counts": _priority_counts(items),
    }
    headline = (
        "No recovery action groups need attention right now."
        if not items
        else (
            f"{summary['total_groups']} recovery action group(s): "
            f"{summary['immediate_count']} immediate, "
            f"{summary['manual_review_groups']} needing manual review."
        )
    )
    return {
        "headline": headline,
        "items": visible_items,
        "all_items": items,
        "summary": summary,
        "filter_views": filter_views,
        "available": bool(items),
    }


def _report_summary(project_root: Path) -> dict:
    report_path = project_root / "project_report.txt"
    combined_path = project_root / "combined_residue_evidence.csv"
    if not report_path.exists():
        return {
            "available": False,
            "path": "step6_report/project_report.txt",
            "headline": "Narrative report is not available yet.",
            "details": [],
            "step_numbers": [6],
        }

    try:
        raw_lines = report_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        raw_lines = []
    digest_lines: List[str] = []
    for line in raw_lines:
        stripped = line.strip()
        if not stripped:
            continue
        if set(stripped) <= {"=", "-"}:
            continue
        digest_lines.append(stripped)
        if len(digest_lines) >= 3:
            break

    combined_rows = _load_csv_rows(combined_path)
    return {
        "available": True,
        "path": "step6_report/project_report.txt",
        "headline": "Narrative report and combined residue evidence are available.",
        "details": [
            *(f"Report digest: {line}" for line in digest_lines),
            (
                f"Combined residue evidence rows: {len(combined_rows)}."
                if combined_path.exists()
                else "Combined residue evidence table is not available."
            ),
        ],
        "step_numbers": [6],
    }


def _quick_links(project_root: Path, *, expected_paths: Optional[Sequence[str]] = None) -> List[dict]:
    expected = {str(item) for item in (expected_paths or [])}
    link_specs = [
        ("run_overview.md", "Markdown overview"),
        ("run_overview.html", "HTML overview"),
        ("report_digest.md", "Condensed report digest"),
        ("operational_recovery_playbook.md", "Operational recovery playbook"),
        ("operational_recovery_plan.json", "Operational recovery plan"),
        ("step_index.md", "Step index"),
        ("run_status.json", "Execution telemetry"),
        ("current_run_manifest.json", "Derived step manifest"),
        ("step6_report/project_report.txt", "Narrative report"),
        ("step6_report/combined_residue_evidence.csv", "Combined residue evidence"),
        ("step5_verdict/valid_sites.csv", "Prioritized site verdicts"),
        ("step4_vina_postprocess/vina_pocket_table.csv", "Pocket summary table"),
        ("step7_validate/validation_summary.txt", "Validation summary"),
        ("step7_validate/validation_recovery_playbook.md", "Validation recovery playbook"),
    ]
    links: List[dict] = []
    for rel_path, description in link_specs:
        links.append(
            {
                "path": rel_path,
                "description": description,
                "available": (project_root / Path(rel_path)).exists() or rel_path in expected,
            }
        )
    return links


def _next_actions(
    *,
    run_status: Optional[dict],
    current_manifest: Optional[dict],
    validation_summary: dict,
    operational_recovery_summary: dict,
    quick_links: Sequence[dict],
) -> List[str]:
    step_status = (current_manifest or {}).get("step_status") or {}
    stale_steps = sorted(
        int(name.replace("step", ""))
        for name, status in step_status.items()
        if status == "stale" and name.startswith("step")
    )
    overall_status = str((run_status or {}).get("overall_status", "available"))
    current_phase_name = str((run_status or {}).get("current_phase_name", "")).strip()
    current_phase_number = (run_status or {}).get("current_phase_number")
    last_error = str((run_status or {}).get("last_error", "")).strip()

    if overall_status == "failed":
        rerun_hint = (
            f"Rerun from Phase {current_phase_number} with `python run_production.py --from {current_phase_number}`."
            if current_phase_number
            else "Rerun the failing phase after fixing the error."
        )
        return [
            f"Execution stopped during {current_phase_name or 'an unknown phase'}.",
            f"Last error: {last_error or 'check run_status.json for the traceback context.'}",
            rerun_hint,
        ]
    if overall_status == "running":
        return [
            f"Current phase: {current_phase_name or 'preparing run metadata'}.",
            "Refresh this overview or open `run_status.json` to monitor progress.",
            "Avoid mixing partially refreshed derived steps with older results until the run finishes.",
        ]
    if validation_summary.get("status") in {"failed", "error"}:
        playbook_path = str(
            validation_summary.get("playbook_path", "step7_validate/validation_recovery_playbook.md")
        )
        recommended_command = str(validation_summary.get("recommended_command", "")).strip()
        manual_review_count = safe_int(validation_summary.get("manual_review_count"), 0)
        return [
            (
                f"Open `{playbook_path}` first. {manual_review_count} finding(s) require manual review before rerun."
                if manual_review_count
                else f"Open `{playbook_path}` and `step7_validate/validation_summary.txt` before trusting the results."
            ),
            (
                f"Safest repair path: `{recommended_command}`."
                if recommended_command and recommended_command != "python run_production.py --only 7"
                else "Fix the reported validation issues, then rerun Phase 7."
            ),
            "After upstream fixes, rerun `python run_production.py --only 7` to confirm the outputs are healthy.",
            "Do not treat the report or verdict as final until validation passes.",
        ]
    if safe_int(operational_recovery_summary.get("issue_count"), 0) > 0:
        recommended_command = str(operational_recovery_summary.get("recommended_command", "")).strip()
        manual_review_count = safe_int(operational_recovery_summary.get("manual_review_count"), 0)
        return [
            (
                f"Open `operational_recovery_playbook.md` first. {manual_review_count} finding(s) need manual review."
                if manual_review_count
                else "Open `operational_recovery_playbook.md` before interpreting downstream results."
            ),
            (
                f"Safest repair path: `{recommended_command}`."
                if recommended_command
                else "Review the operational action groups and clear them before trusting downstream outputs."
            ),
            "Refresh the affected upstream steps before relying on Step 4-6 outputs.",
        ]
    if stale_steps:
        stale_text = ", ".join(f"Step {step_num}" for step_num in stale_steps)
        return [
            f"Rebuild stale derived views before interpretation: {stale_text}.",
            "Start with `step_index.md` to confirm which steps are still stale.",
            "Once refreshed, review the report and verdict outputs in order.",
        ]
    available_links = {item["path"] for item in quick_links if item["available"]}
    if "step6_report/project_report.txt" in available_links:
        return [
            "Start with `step6_report/project_report.txt` for the narrative summary.",
            "Then open `step5_verdict/valid_sites.csv` to inspect prioritized pockets.",
            "Use `step4_vina_postprocess/vina_pocket_table.csv` for detailed pocket evidence.",
        ]
    return [
        "Start with `step_index.md` to see which derived steps are available.",
        "Use `run_status.json` for execution history and `current_run_manifest.json` for step health.",
    ]


def _command_suggestions(
    *,
    run_status: Optional[dict],
    current_manifest: Optional[dict],
    validation_summary: dict,
    operational_recovery_summary: dict,
    quick_links: Sequence[dict],
) -> List[dict]:
    suggestions: List[dict] = []
    seen_commands = set()
    step_status = (current_manifest or {}).get("step_status") or {}
    stale_steps = sorted(
        int(name.replace("step", ""))
        for name, status in step_status.items()
        if status == "stale" and name.startswith("step")
    )
    overall_status = str((run_status or {}).get("overall_status", "available"))
    current_phase_number = (run_status or {}).get("current_phase_number")
    last_error = str((run_status or {}).get("last_error", "")).strip()
    available_links = {item["path"] for item in quick_links if item["available"]}

    def _add(label: str, command: str, reason: str) -> None:
        if command in seen_commands:
            return
        seen_commands.add(command)
        suggestions.append({"label": label, "command": command, "reason": reason})

    if overall_status == "failed" and current_phase_number:
        _add(
            "Resume from failed phase",
            f"python run_production.py --from {current_phase_number}",
            f"Use this after fixing the failure cause{': ' + last_error if last_error else ''}.",
        )
    if validation_summary.get("status") in {"failed", "error"}:
        recommended_command = str(validation_summary.get("recommended_command", "")).strip()
        if recommended_command and recommended_command != "python run_production.py --only 7":
            _add(
                "Repair validation findings",
                recommended_command,
                "Refresh the earliest affected phase and all downstream outputs before validating again.",
            )
        _add(
            "Rerun validation",
            str(validation_summary.get("validation_only_command", "python run_production.py --only 7")),
            "Recheck output integrity after fixing missing files, schema issues, or manual validation blockers.",
        )
    operational_command = str(operational_recovery_summary.get("recommended_command", "")).strip()
    if safe_int(operational_recovery_summary.get("issue_count"), 0) > 0 and operational_command:
        _add(
            "Repair upstream operational issues",
            operational_command,
            "Restore Step 1-3 readiness before relying on pocket, verdict, or report interpretation.",
        )
    if stale_steps:
        _add(
            "Refresh stale derived steps",
            f"python run_production.py --only {','.join(str(step) for step in stale_steps)}",
            "Regenerate the stale interpretation layers without rerunning unrelated phases.",
        )
    if "step5_verdict/valid_sites.csv" not in available_links and "step4_vina_postprocess/vina_pocket_table.csv" in available_links:
        _add(
            "Generate verdict and downstream views",
            "python run_production.py --only 5,6,7",
            "Build the site prioritization, report, and validation outputs from existing pocket summaries.",
        )
    if "step6_report/project_report.txt" not in available_links and "step5_verdict/valid_sites.csv" in available_links:
        _add(
            "Generate report and validation",
            "python run_production.py --only 6,7",
            "Create the narrative report and immediately validate it.",
        )
    if not suggestions:
        _add(
            "Run the standard production flow",
            "python run_production.py",
            "Use the default pipeline entry point when you want the normal end-to-end refresh.",
        )
    return suggestions


def _guided_action_focus(
    *,
    recovery_radar: dict,
    next_actions: Sequence[str],
) -> List[dict]:
    items: List[dict] = []
    for index, radar_item in enumerate(recovery_radar.get("items") or [], start=1):
        message_preview = str((radar_item.get("messages") or [""])[0]).strip()
        preview_text = (
            message_preview
            if message_preview
            else f"{radar_item.get('scope_label', 'This scope')} needs attention."
        )
        command_text = str(radar_item.get("recommended_command") or "").strip()
        if command_text:
            next_step = f"Suggested command: `{command_text}`."
        elif radar_item.get("manual_review_required"):
            next_step = "Manual review is required before any rerun."
        else:
            next_step = "Review this group before deciding whether a rerun is necessary."
        items.append(
            {
                "title": f"{radar_item['source_label']}: {radar_item['title']}",
                "summary": preview_text,
                "detail": next_step,
                "path": str(radar_item.get("playbook_link") or radar_item.get("playbook_path", "")),
                "command": command_text,
                "group_keys": [str(radar_item.get("group_key", ""))] if str(radar_item.get("group_key", "")).strip() else [],
                "source_kind": str(radar_item.get("source_kind", "general")),
                "source_label": str(radar_item.get("source_label", "General")),
                "priority_label": str(radar_item.get("priority_label", "low")),
                "priority_rank": safe_int(radar_item.get("priority_rank"), 3),
                "manual_review_required": bool(radar_item.get("manual_review_required")),
                "default_order": index,
            }
        )
    if items:
        return items

    fallback: List[dict] = []
    for index, action in enumerate(next_actions, start=1):
        fallback.append(
            {
                "title": f"Check {index}",
                "summary": str(action),
                "detail": "Default review order from the latest overview snapshot.",
                "path": "",
                "command": "",
                "group_keys": [],
                "source_kind": "general",
                "source_label": "General",
                "priority_label": "low",
                "priority_rank": 3,
                "manual_review_required": False,
                "default_order": index,
            }
        )
    return fallback


def _annotate_command_suggestions(
    *,
    command_suggestions: Sequence[dict],
    recovery_radar: dict,
    validation_summary: dict,
    operational_recovery_summary: dict,
    run_status: Optional[dict],
) -> List[dict]:
    command_matches: dict = {}
    for radar_item in recovery_radar.get("items") or []:
        command_text = str(radar_item.get("recommended_command") or "").strip()
        if command_text:
            command_matches.setdefault(command_text, []).append(radar_item)

    annotated: List[dict] = []
    failed_phase_command = ""
    current_phase_number = (run_status or {}).get("current_phase_number")
    if str((run_status or {}).get("overall_status", "")) == "failed" and current_phase_number:
        failed_phase_command = f"python run_production.py --from {current_phase_number}"
    validation_only_command = str(validation_summary.get("validation_only_command", "")).strip()
    operational_command = str(operational_recovery_summary.get("recommended_command", "")).strip()
    validation_group_keys = [
        str(item.get("group_key", "")).strip()
        for item in recovery_radar.get("items") or []
        if str(item.get("source_kind", "")) == "validation" and str(item.get("group_key", "")).strip()
    ]
    operational_group_keys = [
        str(item.get("group_key", "")).strip()
        for item in recovery_radar.get("items") or []
        if str(item.get("source_kind", "")) == "operational" and str(item.get("group_key", "")).strip()
    ]

    for index, item in enumerate(command_suggestions, start=1):
        command_text = str(item.get("command", "")).strip()
        source_kind = "general"
        source_label = "General"
        priority_rank = 3
        priority_label = "low"
        manual_review_required = False
        group_keys: List[str] = []
        matched_groups = command_matches.get(command_text) or []
        if matched_groups:
            best_group = min(
                matched_groups,
                key=lambda group: (
                    safe_int(group.get("priority_rank"), 3),
                    0 if not group.get("manual_review_required") else 1,
                ),
            )
            source_kind = str(best_group.get("source_kind", "general"))
            source_label = str(best_group.get("source_label", "General"))
            priority_rank = safe_int(best_group.get("priority_rank"), 3)
            priority_label = str(best_group.get("priority_label", _validation_priority_label(priority_rank)))
            manual_review_required = any(group.get("manual_review_required") for group in matched_groups)
            group_keys = [
                str(group.get("group_key", "")).strip()
                for group in matched_groups
                if str(group.get("group_key", "")).strip()
            ]
        elif validation_only_command and command_text == validation_only_command:
            source_kind = "validation"
            source_label = "Validation"
            priority_rank = 1 if validation_summary.get("status") in {"failed", "error"} else 2
            priority_label = _validation_priority_label(priority_rank)
            manual_review_required = safe_int(validation_summary.get("manual_review_count"), 0) > 0
            group_keys = validation_group_keys
        elif operational_command and command_text == operational_command:
            source_kind = "operational"
            source_label = "Operational"
            priority_rank = 0 if safe_int(operational_recovery_summary.get("issue_count"), 0) > 0 else 2
            priority_label = _validation_priority_label(priority_rank)
            manual_review_required = safe_int(operational_recovery_summary.get("manual_review_count"), 0) > 0
            group_keys = operational_group_keys
        elif failed_phase_command and command_text == failed_phase_command:
            priority_rank = 0
            priority_label = _validation_priority_label(priority_rank)
        elif str(item.get("label", "")).startswith("Refresh stale derived steps"):
            source_kind = "operational"
            source_label = "Operational"
            priority_rank = 1
            priority_label = _validation_priority_label(priority_rank)
            group_keys = operational_group_keys

        annotated.append(
            {
                "label": str(item.get("label", "")),
                "command": command_text,
                "reason": str(item.get("reason", "")),
                "group_keys": group_keys,
                "source_kind": source_kind,
                "source_label": source_label,
                "priority_label": priority_label,
                "priority_rank": priority_rank,
                "manual_review_required": manual_review_required,
                "default_order": index,
            }
        )
    return annotated


def _annotate_quick_links(
    *,
    quick_links: Sequence[dict],
    recovery_radar: dict,
    validation_summary: dict,
    operational_recovery_summary: dict,
    run_status: Optional[dict],
    current_manifest: Optional[dict],
) -> List[dict]:
    overall_status = str((run_status or {}).get("overall_status", "available"))
    step_status = (current_manifest or {}).get("step_status") or {}
    stale_present = any(status == "stale" for status in step_status.values())
    validation_active = validation_summary.get("status") in {"failed", "error"}
    operational_active = safe_int(operational_recovery_summary.get("issue_count"), 0) > 0
    validation_manual = safe_int(validation_summary.get("manual_review_count"), 0) > 0
    operational_manual = safe_int(operational_recovery_summary.get("manual_review_count"), 0) > 0
    validation_group_keys = [
        str(item.get("group_key", "")).strip()
        for item in recovery_radar.get("items") or []
        if str(item.get("source_kind", "")) == "validation" and str(item.get("group_key", "")).strip()
    ]
    operational_group_keys = [
        str(item.get("group_key", "")).strip()
        for item in recovery_radar.get("items") or []
        if str(item.get("source_kind", "")) == "operational" and str(item.get("group_key", "")).strip()
    ]

    annotated: List[dict] = []
    for index, item in enumerate(quick_links, start=1):
        path = str(item.get("path", ""))
        description = str(item.get("description", ""))
        available = bool(item.get("available"))
        source_kind = "overview"
        source_label = "Overview"
        priority_rank = 3
        manual_review_required = False
        reason = "Reference entry point for the latest derived outputs."
        group_keys: List[str] = []

        if path == "step7_validate/validation_recovery_playbook.md":
            source_kind = "validation"
            source_label = "Validation"
            priority_rank = 0 if validation_active else 2
            manual_review_required = validation_manual
            reason = "Primary recovery checklist when validation blocks final sign-off."
            group_keys = validation_group_keys
        elif path == "step7_validate/validation_summary.txt":
            source_kind = "validation"
            source_label = "Validation"
            priority_rank = 1 if validation_active else 2
            manual_review_required = validation_manual
            reason = "Fast validation snapshot with the safest rerun path."
            group_keys = validation_group_keys
        elif path == "operational_recovery_playbook.md":
            source_kind = "operational"
            source_label = "Operational"
            priority_rank = 0 if operational_active else 2
            manual_review_required = operational_manual
            reason = "Primary upstream recovery checklist before trusting downstream results."
            group_keys = operational_group_keys
        elif path == "operational_recovery_plan.json":
            source_kind = "operational"
            source_label = "Operational"
            priority_rank = 1 if operational_active else 2
            manual_review_required = operational_manual
            reason = "Machine-readable operational triage summary."
            group_keys = operational_group_keys
        elif path in {"run_status.json", "current_run_manifest.json", "step_index.md"}:
            source_kind = "operational" if overall_status in {"running", "failed"} or operational_active or stale_present else "overview"
            source_label = "Operational" if source_kind == "operational" else "Overview"
            priority_rank = 1 if overall_status in {"running", "failed"} or stale_present else (2 if operational_active else 3)
            reason = "Best source for live execution state, stale steps, and derived health."
            group_keys = operational_group_keys if source_kind == "operational" else []
        elif path in {"run_overview.md", "run_overview.html", "report_digest.md"}:
            source_kind = "overview"
            source_label = "Overview"
            priority_rank = 2 if not (validation_active or operational_active or overall_status == "running") else 3
            reason = "User-facing summary layer for fast scanning and triage."
        elif path in {
            "step6_report/project_report.txt",
            "step6_report/combined_residue_evidence.csv",
            "step5_verdict/valid_sites.csv",
            "step4_vina_postprocess/vina_pocket_table.csv",
        }:
            source_kind = "results"
            source_label = "Results"
            priority_rank = 1 if not (validation_active or operational_active or overall_status == "running") else 3
            reason = "Core interpretation artifact once recovery blockers are cleared."

        annotated.append(
            {
                "path": path,
                "description": description,
                "available": available,
                "group_keys": group_keys,
                "source_kind": source_kind,
                "source_label": source_label,
                "priority_rank": priority_rank,
                "priority_label": _validation_priority_label(priority_rank),
                "manual_review_required": manual_review_required,
                "reason": reason,
                "availability_rank": 0 if available else 1,
                "default_order": index,
            }
        )
    return annotated


def build_run_overview_data(
    project_root: Union[Path, str],
    *,
    current_run_manifest: Optional[dict] = None,
    run_status: Optional[dict] = None,
    operational_recovery_plan: Optional[dict] = None,
    validation_recovery_plan: Optional[dict] = None,
    expected_paths: Optional[Sequence[str]] = None,
) -> dict:
    """Build a user-facing overview payload for the project root."""

    project_root = _as_path(project_root)
    current = current_run_manifest or _load_current_run_manifest(project_root) or {}
    status_payload = run_status or _load_run_status(project_root) or {}
    operational_plan = operational_recovery_plan or _load_operational_recovery_plan(project_root)
    validation_plan = validation_recovery_plan or _load_validation_recovery_plan(project_root)
    phase_rows = list(status_payload.get("phase_states") or [])
    if not phase_rows:
        step_status = current.get("step_status") or _step_status_map(project_root)
        for step_num, spec in STEP_SPECS.items():
            phase_rows.append(
                {
                    "phase_number": step_num,
                    "phase_name": f"Phase {step_num}: {spec.step_name}",
                    "status": step_status.get(f"step{step_num}", "not_started"),
                    "started_at": "",
                    "completed_at": "",
                    "duration_seconds": None,
                    "skip_reason": "",
                    "last_error": "",
                }
            )

    resolved_count, total_count, progress_percent = _progress_summary(phase_rows)
    step_health = _step_health_summary(current, project_root)
    operational_recovery_summary = _operational_recovery_summary_for_overview(
        project_root,
        current_manifest=current or None,
        recovery_plan=operational_plan,
    )
    operational_summaries = [
        _with_step_summary_links(_step1_input_summary(project_root)),
        _with_step_summary_links(_step2_ppi_summary(project_root)),
        _with_step_summary_links(_step3_evidence_summary(project_root)),
        _with_step_summary_links(operational_recovery_summary),
    ]
    pockets = _with_step_summary_links(_pocket_summary(project_root))
    verdicts = _with_step_summary_links(_verdict_summary(project_root))
    reports = _with_step_summary_links(_report_summary(project_root))
    validation_summary = _with_step_summary_links(_validation_summary_for_overview(project_root))
    recovery_radar = _build_recovery_radar(
        operational_recovery_plan=operational_plan,
        validation_recovery_plan=validation_plan,
    )
    quick_links = _quick_links(project_root, expected_paths=expected_paths)
    next_actions = _next_actions(
        run_status=status_payload or None,
        current_manifest=current or None,
        validation_summary=validation_summary,
        operational_recovery_summary=operational_recovery_summary,
        quick_links=quick_links,
    )
    command_suggestions = _command_suggestions(
        run_status=status_payload or None,
        current_manifest=current or None,
        validation_summary=validation_summary,
        operational_recovery_summary=operational_recovery_summary,
        quick_links=quick_links,
    )
    action_focus = _guided_action_focus(
        recovery_radar=recovery_radar,
        next_actions=next_actions,
    )
    annotated_command_suggestions = _annotate_command_suggestions(
        command_suggestions=command_suggestions,
        recovery_radar=recovery_radar,
        validation_summary=validation_summary,
        operational_recovery_summary=operational_recovery_summary,
        run_status=status_payload or None,
    )
    annotated_quick_links = _annotate_quick_links(
        quick_links=quick_links,
        recovery_radar=recovery_radar,
        validation_summary=validation_summary,
        operational_recovery_summary=operational_recovery_summary,
        run_status=status_payload or None,
        current_manifest=current or None,
    )

    overall_status = str(status_payload.get("overall_status", "available"))
    current_phase_name = str(status_payload.get("current_phase_name", "")).strip()
    if overall_status == "running":
        hero_title = "Pipeline run is in progress."
        hero_body = f"The latest update is tracking {current_phase_name or 'the current phase'}."
    elif overall_status == "failed":
        hero_title = "Pipeline run stopped before completion."
        hero_body = "Review the failing phase and rerun only the affected scope."
    elif validation_summary.get("status") in {"failed", "error"}:
        hero_title = "Pipeline run finished, but validation found issues."
        hero_body = "Resolve validation findings before treating the outputs as final."
    elif overall_status == "completed":
        hero_title = "Pipeline run completed."
        hero_body = "Start with the report and verdict outputs for interpretation."
    else:
        hero_title = "Derived outputs are available for review."
        hero_body = "This overview summarizes the latest known step and execution state."

    return {
        "project_name": current.get("project_name", status_payload.get("project_name", project_root.name)),
        "overall_status": overall_status,
        "hero_title": hero_title,
        "hero_body": hero_body,
        "execution_mode": current.get("execution_mode", status_payload.get("execution_mode", "unknown")),
        "generated_at": current.get("generated_at", status_payload.get("updated_at", utc_now_iso())),
        "updated_at": status_payload.get("updated_at", current.get("generated_at", utc_now_iso())),
        "current_phase_name": current_phase_name,
        "current_phase_number": status_payload.get("current_phase_number"),
        "resolved_count": resolved_count,
        "total_count": total_count,
        "progress_percent": progress_percent,
        "phase_rows": phase_rows,
        "step_health": step_health,
        "operational_summaries": operational_summaries,
        "operational_recovery_summary": operational_recovery_summary,
        "recovery_radar": recovery_radar,
        "pocket_summary": pockets,
        "verdict_summary": verdicts,
        "report_summary": reports,
        "validation_summary": validation_summary,
        "quick_links": quick_links,
        "next_actions": next_actions,
        "command_suggestions": command_suggestions,
        "action_focus": action_focus,
        "annotated_command_suggestions": annotated_command_suggestions,
        "annotated_quick_links": annotated_quick_links,
    }


def write_run_overview_markdown(project_root: Union[Path, str], overview_data: dict) -> Path:
    """Write ``run_overview.md`` at the project root."""

    project_root = _as_path(project_root)
    recovery_radar = overview_data["recovery_radar"]
    step_health = overview_data["step_health"]
    validation_summary = overview_data["validation_summary"]
    lines = [
        "# Run Overview",
        "",
        f"## {overview_data['hero_title']}",
        overview_data["hero_body"],
        "",
        "## At a Glance",
        f"- Project name: `{overview_data['project_name']}`",
        f"- Overall execution status: `{overview_data['overall_status']}`",
        f"- Execution mode: `{overview_data['execution_mode']}`",
        f"- Progress: `{overview_data['resolved_count']}/{overview_data['total_count']}` phase(s) resolved ({overview_data['progress_percent']}%).",
        f"- Current phase: `{overview_data['current_phase_name'] or 'n/a'}`",
        f"- Last update: `{overview_data['updated_at']}`",
        (
            f"- Validation: `{validation_summary['status']}`"
            if validation_summary.get("available")
            else "- Validation: `not recorded`"
        ),
        (
            "- Derived step health: "
            f"`{step_health['complete_count']}` complete, "
            f"`{step_health['partial_count']}` partial, "
            f"`{step_health['stale_count']}` stale, "
            f"`{step_health['failed_count']}` failed/error, "
            f"`{step_health['not_generated_count']}` not generated."
        ),
        "",
        "## What To Check First",
    ]
    for idx, action in enumerate(overview_data["next_actions"], start=1):
        lines.append(f"{idx}. {action}")

    lines.extend(["", "## Recovery Radar", f"- {recovery_radar['headline']}"])
    radar_summary = recovery_radar["summary"]
    lines.append(
        "- Radar counts: "
        f"`{radar_summary['validation_groups']}` validation, "
        f"`{radar_summary['operational_groups']}` operational, "
        f"`{radar_summary['manual_review_groups']}` manual-review group(s)."
    )
    lines.extend(["", "### Filter Views"])
    for filter_view in recovery_radar["filter_views"]:
        lines.append(f"- {filter_view['label']}: {filter_view['headline']}")
    if recovery_radar["items"]:
        for index, item in enumerate(recovery_radar["items"], start=1):
            lines.extend(
                [
                    f"### Radar {index}: {item['source_label']} - {item['title']}",
                    f"- Scope: `{item['scope_label']}`",
                    f"- Priority: `{item['priority_label']}`",
                    f"- Action path: `{item['action_label']}`",
                    (
                        f"- Suggested command: `{item['recommended_command']}`"
                        if item["recommended_command"]
                        else "- Suggested command: manual review only"
                    ),
                    f"- Playbook: `{item.get('playbook_link', item['playbook_path'])}`",
                ]
            )
            if item.get("step_summary_links"):
                lines.append(
                    f"- Step summaries: {_step_summary_links_text(item['step_summary_links'])}"
                )
            if item["artifacts"]:
                artifact_text = ", ".join(f"`{name}`" for name in item["artifacts"])
                if item["artifact_overflow_count"]:
                    artifact_text += f" (+{item['artifact_overflow_count']} more)"
                lines.append(f"- Related artifacts: {artifact_text}")
            if item["messages"]:
                message_text = " | ".join(item["messages"])
                if item["message_overflow_count"]:
                    message_text += f" (+{item['message_overflow_count']} more)"
                lines.append(f"- Findings: {message_text}")
            lines.append("")
    else:
        lines.append("- No recovery action groups are active right now.")

    lines.extend(
        [
            "",
            "## Pipeline Progress",
            "| Phase | Status | Started | Finished | Notes |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for row in overview_data["phase_rows"]:
        note = str(row.get("skip_reason") or row.get("last_error") or "").replace("\n", " ").strip()
        lines.append(
            f"| {row.get('phase_number', '?')} | `{row.get('status', 'unknown')}` | "
            f"`{row.get('started_at', '') or '-'}` | `{row.get('completed_at', '') or '-'}` | "
            f"{note or '-'} |"
        )

    lines.extend(["", "## Operational Readiness"])
    for item in overview_data["operational_summaries"]:
        lines.append(f"- {item['title']}: {item['headline']}")
        if item.get("step_summary_links"):
            lines.append(f"  - Step summaries: {_step_summary_links_text(item['step_summary_links'])}")
        lines.extend(f"  - {detail}" for detail in item["details"])

    lines.extend(
        [
            "",
            "## Result Highlights",
            f"- Pocket summary: {overview_data['pocket_summary']['headline']}",
        ]
    )
    if overview_data["pocket_summary"].get("step_summary_links"):
        lines.append(
            f"  - Step summaries: {_step_summary_links_text(overview_data['pocket_summary']['step_summary_links'])}"
        )
    lines.extend(f"  - {detail}" for detail in overview_data["pocket_summary"]["details"])
    lines.append(f"- Verdict summary: {overview_data['verdict_summary']['headline']}")
    if overview_data["verdict_summary"].get("step_summary_links"):
        lines.append(
            f"  - Step summaries: {_step_summary_links_text(overview_data['verdict_summary']['step_summary_links'])}"
        )
    lines.extend(f"  - {detail}" for detail in overview_data["verdict_summary"]["details"])
    lines.append(f"- Report summary: {overview_data['report_summary']['headline']}")
    if overview_data["report_summary"].get("step_summary_links"):
        lines.append(
            f"  - Step summaries: {_step_summary_links_text(overview_data['report_summary']['step_summary_links'])}"
        )
    lines.extend(f"  - {detail}" for detail in overview_data["report_summary"]["details"])
    lines.append(f"- Validation summary: {overview_data['validation_summary']['headline']}")
    if overview_data["validation_summary"].get("step_summary_links"):
        lines.append(
            f"  - Step summaries: {_step_summary_links_text(overview_data['validation_summary']['step_summary_links'])}"
        )
    lines.extend(f"  - {detail}" for detail in overview_data["validation_summary"]["details"])

    lines.extend(["", "## Quick Links"])
    for item in overview_data["quick_links"]:
        status = "available" if item["available"] else "missing"
        lines.append(f"- `{item['path']}`: {item['description']} ({status})")
    lines.extend(["", "## Suggested Commands"])
    for item in overview_data["command_suggestions"]:
        lines.append(f"- `{item['command']}`: {item['label']} - {item['reason']}")
    lines.append("")
    return _atomic_write_text(project_root / "run_overview.md", "\n".join(lines))


def write_run_overview_html(project_root: Union[Path, str], overview_data: dict) -> Path:
    """Write ``run_overview.html`` at the project root."""

    project_root = _as_path(project_root)
    phase_rows_html = []
    for row in overview_data["phase_rows"]:
        note = str(row.get("skip_reason") or row.get("last_error") or "").replace("\n", " ").strip() or "-"
        phase_rows_html.append(
            "<tr>"
            f"<td>{escape(str(row.get('phase_number', '?')))}</td>"
            f"<td><span class=\"badge badge-{escape(str(row.get('status', 'unknown')))}\">{escape(_status_label(str(row.get('status', 'unknown'))))}</span></td>"
            f"<td>{escape(str(row.get('phase_name', '')))}</td>"
            f"<td>{escape(str(row.get('started_at', '') or '-'))}</td>"
            f"<td>{escape(str(row.get('completed_at', '') or '-'))}</td>"
            f"<td>{escape(note)}</td>"
            "</tr>"
        )

    def _step_summary_links_html(step_summary_links: Sequence[dict], *, label: str = "Step summaries") -> str:
        if not step_summary_links:
            return ""
        link_html = ", ".join(
            f"<a href=\"{escape(str(item.get('path', '')))}\"><code>{escape(str(item.get('path', '')))}</code></a>"
            for item in step_summary_links
            if str(item.get("path", "")).strip()
        )
        if not link_html:
            return ""
        return f"<p><strong>{escape(label)}:</strong> {link_html}</p>"

    quick_links_html = []
    for item in overview_data.get("annotated_quick_links", overview_data["quick_links"]):
        if item["available"]:
            link_body = f"<a href=\"{escape(item['path'])}\"><code>{escape(item['path'])}</code></a>"
        else:
            link_body = f"<code>{escape(item['path'])}</code>"
        availability_text = "Available" if item["available"] else "Missing"
        group_keys = " ".join(
            str(group_key).strip()
            for group_key in (item.get("group_keys") or [])
            if str(group_key).strip()
        )
        quick_links_html.append(
            f"<li class=\"guided-item guided-link\" "
            f"data-source-kind=\"{escape(str(item.get('source_kind', 'overview')))}\" "
            f"data-manual-review=\"{'yes' if item.get('manual_review_required') else 'no'}\" "
            f"data-priority=\"{escape(str(item.get('priority_label', 'low')))}\" "
            f"data-priority-rank=\"{escape(str(safe_int(item.get('priority_rank'), 3)))}\" "
            f"data-availability-rank=\"{escape(str(safe_int(item.get('availability_rank'), 1)))}\" "
            f"data-group-keys=\"{escape(group_keys)}\" "
            f"data-default-order=\"{escape(str(safe_int(item.get('default_order'), 0)))}\">"
            f"<div class=\"tone-row\">"
            f"<span class=\"tone tone-{escape(str(item.get('source_kind', 'overview')))}\">{escape(str(item.get('source_label', 'Overview')))}</span>"
            f"<span class=\"tone priority-{escape(str(item.get('priority_label', 'low')))}\">{escape(str(item.get('priority_label', 'low')).title())}</span>"
            f"<span class=\"tone tone-summary\">{escape(availability_text)}</span>"
            f"</div>"
            f"<p>{link_body} - {escape(str(item.get('description', '')))}</p>"
            f"<p>{escape(str(item.get('reason', '')))}</p>"
            "</li>"
        )

    card_specs = [
        ("Pocket Summary", overview_data["pocket_summary"]),
        ("Verdict Summary", overview_data["verdict_summary"]),
        ("Report Summary", overview_data["report_summary"]),
        ("Validation Summary", overview_data["validation_summary"]),
    ]
    cards_html = []
    for index, (title, payload) in enumerate(card_specs, start=1):
        details_html = "".join(f"<li>{escape(detail)}</li>" for detail in payload["details"]) or "<li>No extra details.</li>"
        step_numbers = _normalize_step_numbers(payload.get("step_numbers"))
        step_label = _step_context_label(step_numbers)
        step_numbers_attr = " ".join(str(step_number) for step_number in step_numbers)
        step_chip_html = (
            f"<div class=\"tone-row\"><span class=\"tone tone-step\">{escape(step_label)}</span></div>"
            if step_label
            else ""
        )
        step_summary_links_html = _step_summary_links_html(payload.get("step_summary_links") or [])
        cards_html.append(
            f"<section class=\"card result-card\" data-step-numbers=\"{escape(step_numbers_attr)}\" data-default-order=\"{index}\">"
            f"{step_chip_html}"
            f"<h3>{escape(title)}</h3>"
            f"<p>{escape(payload['headline'])}</p>"
            f"<ul>{details_html}</ul>"
            f"{step_summary_links_html}"
            f"<p class=\"path\"><code>{escape(payload['path'])}</code></p>"
            "</section>"
        )

    action_html = []
    for item in overview_data.get("action_focus", []):
        source_kind = escape(str(item.get("source_kind", "general")))
        manual_review = "yes" if item.get("manual_review_required") else "no"
        priority_label = escape(str(item.get("priority_label", "low")))
        priority_rank = escape(str(safe_int(item.get("priority_rank"), 3)))
        default_order = escape(str(safe_int(item.get("default_order"), 0)))
        chips = [
            f"<span class=\"tone tone-{source_kind}\">{escape(str(item.get('source_label', 'General')))}</span>",
            f"<span class=\"tone priority-{priority_label}\">{escape(str(item.get('priority_label', 'low')).title())}</span>",
        ]
        if item.get("manual_review_required"):
            chips.append("<span class=\"tone tone-manual\">Manual review</span>")
        path_html = ""
        if item.get("path"):
            path_html = (
                f"<p class=\"path\"><a href=\"{escape(str(item['path']))}\"><code>{escape(str(item['path']))}</code></a></p>"
            )
        command_line = (
            f"<p><strong>Command:</strong> <code>{escape(str(item['command']))}</code></p>"
            if item.get("command")
            else ""
        )
        group_keys = " ".join(
            str(group_key).strip()
            for group_key in (item.get("group_keys") or [])
            if str(group_key).strip()
        )
        action_html.append(
            f"<li class=\"guided-item guided-action\" "
            f"data-source-kind=\"{source_kind}\" "
            f"data-manual-review=\"{manual_review}\" "
            f"data-priority=\"{priority_label}\" "
            f"data-priority-rank=\"{priority_rank}\" "
            f"data-group-keys=\"{escape(group_keys)}\" "
            f"data-default-order=\"{default_order}\">"
            f"<div class=\"tone-row\">{''.join(chips)}</div>"
            f"<strong>{escape(str(item.get('title', 'Action')))}</strong>"
            f"<p>{escape(str(item.get('summary', '')))}</p>"
            f"<p>{escape(str(item.get('detail', '')))}</p>"
            f"{command_line}"
            f"{path_html}"
            "</li>"
        )
    operational_cards_html = []
    for index, payload in enumerate(overview_data["operational_summaries"], start=1):
        details_html = "".join(f"<li>{escape(detail)}</li>" for detail in payload["details"]) or "<li>No extra details.</li>"
        step_numbers = _normalize_step_numbers(payload.get("step_numbers"))
        step_label = _step_context_label(step_numbers)
        step_numbers_attr = " ".join(str(step_number) for step_number in step_numbers)
        step_chip_html = (
            f"<div class=\"tone-row\"><span class=\"tone tone-step\">{escape(step_label)}</span></div>"
            if step_label
            else ""
        )
        step_summary_links_html = _step_summary_links_html(payload.get("step_summary_links") or [])
        operational_cards_html.append(
            f"<section class=\"card operational-card\" data-step-numbers=\"{escape(step_numbers_attr)}\" data-default-order=\"{index}\">"
            f"{step_chip_html}"
            f"<h3>{escape(payload['title'])}</h3>"
            f"<p>{escape(payload['headline'])}</p>"
            f"<ul>{details_html}</ul>"
            f"{step_summary_links_html}"
            f"<p class=\"path\"><code>{escape(payload['path'])}</code></p>"
            "</section>"
        )
    radar_summary = overview_data["recovery_radar"]["summary"]
    radar_filter_views = {
        item["filter_key"]: item for item in overview_data["recovery_radar"].get("filter_views", [])
    }
    radar_summary_html = "".join(
        (
            f"<span class=\"tone tone-summary\">Validation {escape(str(radar_summary['validation_groups']))}</span>"
            f"<span class=\"tone tone-summary\">Operational {escape(str(radar_summary['operational_groups']))}</span>"
            f"<span class=\"tone tone-summary\">Immediate {escape(str(radar_summary['immediate_count']))}</span>"
            f"<span class=\"tone tone-summary\">Manual review {escape(str(radar_summary['manual_review_groups']))}</span>"
        )
    )
    recovery_radar_html = []
    for index, item in enumerate(overview_data["recovery_radar"]["items"], start=1):
        chips = [
            f"<span class=\"tone tone-{escape(item['source_kind'])}\">{escape(item['source_label'])}</span>",
            f"<span class=\"tone priority-{escape(item['priority_label'])}\">{escape(item['priority_label'].title())}</span>",
            f"<span class=\"tone tone-action\">{escape(item['action_label'])}</span>",
        ]
        if item["manual_review_required"]:
            chips.append("<span class=\"tone tone-manual\">Manual review</span>")
        findings_html = "".join(f"<li>{escape(message)}</li>" for message in item["messages"]) or "<li>No finding preview.</li>"
        if item["message_overflow_count"]:
            findings_html += f"<li>+{escape(str(item['message_overflow_count']))} more finding(s) in the playbook.</li>"
        artifact_html = ""
        if item["artifacts"]:
            artifact_list = ", ".join(f"<code>{escape(name)}</code>" for name in item["artifacts"])
            if item["artifact_overflow_count"]:
                artifact_list += f" (+{escape(str(item['artifact_overflow_count']))} more)"
            artifact_html = f"<p><strong>Artifacts:</strong> {artifact_list}</p>"
        command_html_line = (
            f"<code>{escape(item['recommended_command'])}</code>"
            if item["recommended_command"]
            else "manual review only"
        )
        source_kind = escape(item["source_kind"])
        manual_review = "yes" if item["manual_review_required"] else "no"
        priority_label = escape(item["priority_label"])
        group_key = escape(str(item.get("group_key", "")))
        step_numbers = _normalize_step_numbers(item.get("step_numbers"))
        step_numbers_attr = escape(" ".join(str(step_number) for step_number in step_numbers))
        step_label = _step_context_label(step_numbers)
        step_scope_html = (
            f"<p><strong>Related steps:</strong> {escape(step_label)}</p>"
            if step_label
            else ""
        )
        step_summary_links_html = _step_summary_links_html(item.get("step_summary_links") or [])
        recovery_radar_html.append(
            f"<article id=\"radar-card-{index}\" class=\"card radar-card\" "
            f"tabindex=\"0\" "
            f"role=\"button\" "
            f"aria-pressed=\"false\" "
            f"data-group-key=\"{group_key}\" "
            f"data-step-numbers=\"{step_numbers_attr}\" "
            f"data-source-kind=\"{source_kind}\" "
            f"data-manual-review=\"{manual_review}\" "
            f"data-priority=\"{priority_label}\">"
            f"<div class=\"tone-row\">{''.join(chips)}</div>"
            f"<h3>{escape(item['title'])}</h3>"
            f"<p><strong>Scope:</strong> {escape(item['scope_label'])}</p>"
            f"{step_scope_html}"
            f"{step_summary_links_html}"
            f"<p><strong>Suggested command:</strong> {command_html_line}</p>"
            f"{artifact_html}"
            f"<ul>{findings_html}</ul>"
            f"<p class=\"path\"><a href=\"{escape(item.get('playbook_link', item['playbook_path']))}\"><code>{escape(item.get('playbook_link', item['playbook_path']))}</code></a></p>"
            "</article>"
        )
    jump_links_html = "".join(
        f"<a class=\"jump-link\" href=\"#{escape(target)}\">{escape(label)}</a>"
        for target, label in (
            ("what-to-check-first", "What To Check First"),
            ("recovery-radar", "Recovery Radar"),
            ("result-highlights", "Result Highlights"),
            ("operational-readiness", "Operational Readiness"),
            ("pipeline-progress", "Pipeline Progress"),
            ("quick-links", "Quick Links"),
            ("suggested-commands", "Suggested Commands"),
        )
    )
    radar_filter_controls_html = ""
    radar_filter_script = ""
    if recovery_radar_html:
        all_groups = radar_filter_views.get("all", {}).get("group_count", len(recovery_radar_html))
        validation_groups = radar_filter_views.get("validation", {}).get("group_count", 0)
        operational_groups = radar_filter_views.get("operational", {}).get("group_count", 0)
        manual_groups = radar_filter_views.get("manual", {}).get("group_count", 0)
        priority_counts = radar_summary.get("priority_counts", {})
        radar_filter_controls_html = "".join(
            (
                '<div class="radar-filter-stack">'
                '<div class="radar-filter-group">'
                '<p class="filter-label">Scope</p>'
                '<div class="radar-filters" role="toolbar" aria-label="Recovery radar scope filters">'
                f'<button type="button" class="filter-chip is-active" data-radar-filter="all" aria-pressed="true">All Groups ({escape(str(all_groups))})</button>'
                f'<button type="button" class="filter-chip" data-radar-filter="validation" aria-pressed="false">Validation ({escape(str(validation_groups))})</button>'
                f'<button type="button" class="filter-chip" data-radar-filter="operational" aria-pressed="false">Operational ({escape(str(operational_groups))})</button>'
                f'<button type="button" class="filter-chip" data-radar-filter="manual" aria-pressed="false">Manual Review ({escape(str(manual_groups))})</button>'
                "</div>"
                "</div>"
                '<div class="radar-filter-group">'
                '<p class="filter-label">Priority</p>'
                '<div class="radar-filters" role="toolbar" aria-label="Recovery radar priority filters">'
                f'<button type="button" class="filter-chip is-active" data-radar-priority-filter="all" aria-pressed="true">All Priorities ({escape(str(all_groups))})</button>'
                f'<button type="button" class="filter-chip" data-radar-priority-filter="immediate" aria-pressed="false">Immediate ({escape(str(priority_counts.get("immediate", 0)))})</button>'
                f'<button type="button" class="filter-chip" data-radar-priority-filter="high" aria-pressed="false">High ({escape(str(priority_counts.get("high", 0)))})</button>'
                f'<button type="button" class="filter-chip" data-radar-priority-filter="medium" aria-pressed="false">Medium ({escape(str(priority_counts.get("medium", 0)))})</button>'
                f'<button type="button" class="filter-chip" data-radar-priority-filter="low" aria-pressed="false">Low ({escape(str(priority_counts.get("low", 0)))})</button>'
                "</div>"
                "</div>"
                "</div>"
                f'<p id="radar-filter-status" class="filter-status" role="status">Showing {escape(str(len(recovery_radar_html)))} recovery group(s) for all scopes and all priorities.</p>'
            )
        )
        radar_filter_script = """
  <script>
    (function () {
      const scopeButtons = Array.from(document.querySelectorAll("[data-radar-filter]"));
      const priorityButtons = Array.from(document.querySelectorAll("[data-radar-priority-filter]"));
      const sectionLinks = Array.from(document.querySelectorAll(".jump-link"));
      const cards = Array.from(document.querySelectorAll(".radar-card"));
      const radarList = document.getElementById("recovery-radar-list");
      const actionList = document.getElementById("guided-action-list");
      const actionItems = Array.from(document.querySelectorAll(".guided-action"));
      const actionStatus = document.getElementById("action-filter-status");
      const actionEmpty = document.getElementById("action-filter-empty");
      const commandList = document.getElementById("guided-command-list");
      const commandItems = Array.from(document.querySelectorAll(".guided-command"));
      const commandStatus = document.getElementById("command-filter-status");
      const commandEmpty = document.getElementById("command-filter-empty");
      const linkList = document.getElementById("prioritized-quick-links");
      const linkItems = Array.from(document.querySelectorAll(".guided-link"));
      const linkStatus = document.getElementById("quick-link-filter-status");
      const resultCardList = document.getElementById("result-card-list");
      const resultCards = Array.from(document.querySelectorAll(".result-card"));
      const resultStatus = document.getElementById("result-context-status");
      const operationalCardList = document.getElementById("operational-card-list");
      const operationalCards = Array.from(document.querySelectorAll(".operational-card"));
      const operationalStatus = document.getElementById("operational-context-status");
      const groupStatus = document.getElementById("group-focus-status");
      const clearGroupButton = document.getElementById("clear-group-focus");
      const toggleStepRadarButton = document.getElementById("toggle-step-radar-visibility");
      const status = document.getElementById("radar-filter-status");
      const knownSections = sectionLinks
        .map(function (link) { return (link.getAttribute("href") || "").replace(/^#/, ""); })
        .filter(Boolean);
      const knownGroupKeys = cards
        .map(function (card) { return card.dataset.groupKey || ""; })
        .filter(Boolean);
      const knownStepValues = Array.from(new Set(
        cards.concat(resultCards).concat(operationalCards)
          .reduce(function (values, card) {
            return values.concat(stepNumbersFor(card));
          }, [])
      ));
      const state = { scope: "all", priority: "all", section: "", group: "", step: "", stepRadarExpanded: false };
      if (!scopeButtons.length || !priorityButtons.length || !cards.length || !status) {
        return;
      }

      function normalizeSelection(buttons, attributeName, value, fallbackValue) {
        const normalized = value || fallbackValue;
        const supported = buttons.some(function (button) {
          return button.getAttribute(attributeName) === normalized;
        });
        return supported ? normalized : fallbackValue;
      }

      function parseHashState() {
        const rawHash = window.location.hash.replace(/^#/, "");
        const parsed = { scope: "all", priority: "all", section: "", group: "", step: "" };
        if (!rawHash) {
          return parsed;
        }
        if (rawHash.indexOf("=") === -1 && rawHash.indexOf("&") === -1) {
          parsed.section = rawHash;
        } else {
          const params = new URLSearchParams(rawHash);
          parsed.section = params.get("section") || "";
          parsed.scope = params.get("scope") || "all";
          parsed.priority = params.get("priority") || "all";
          parsed.group = params.get("group") || "";
          parsed.step = params.get("step") || "";
        }
        parsed.scope = normalizeSelection(scopeButtons, "data-radar-filter", parsed.scope, "all");
        parsed.priority = normalizeSelection(priorityButtons, "data-radar-priority-filter", parsed.priority, "all");
        if (knownSections.indexOf(parsed.section) === -1) {
          parsed.section = "";
        }
        if (knownGroupKeys.indexOf(parsed.group) === -1) {
          parsed.group = "";
        }
        if (knownStepValues.indexOf(Number(parsed.step)) === -1) {
          parsed.step = "";
        }
        return parsed;
      }

      function buildHashState() {
        if (state.section && state.scope === "all" && state.priority === "all" && !state.group && !state.step) {
          return "#" + state.section;
        }
        const params = new URLSearchParams();
        if (state.section) {
          params.set("section", state.section);
        }
        if (state.scope !== "all") {
          params.set("scope", state.scope);
        }
        if (state.priority !== "all") {
          params.set("priority", state.priority);
        }
        if (state.group) {
          params.set("group", state.group);
        }
        if (state.step) {
          params.set("step", state.step);
        }
        const nextState = params.toString();
        return nextState ? "#" + nextState : "";
      }

      function syncHashState(pushEntry) {
        const nextHash = buildHashState();
        const baseUrl = window.location.pathname + window.location.search;
        if (pushEntry && nextHash) {
          window.location.hash = nextHash;
          return;
        }
        if (window.history && window.history.replaceState) {
          window.history.replaceState(null, "", nextHash ? baseUrl + nextHash : baseUrl);
        } else if (nextHash) {
          window.location.hash = nextHash;
        }
      }

      function scrollToSection(sectionId) {
        if (!sectionId) {
          return;
        }
        const target = document.getElementById(sectionId);
        if (!target) {
          return;
        }
        target.scrollIntoView({ block: "start" });
      }

      function matchesScope(card, filterValue) {
        if (filterValue === "all") {
          return true;
        }
        if (filterValue === "manual") {
          return card.dataset.manualReview === "yes";
        }
        return card.dataset.sourceKind === filterValue;
      }

      function matchesPriority(card, filterValue) {
        if (filterValue === "all") {
          return true;
        }
        return card.dataset.priority === filterValue;
      }

      function matchesCurrentFilters(item) {
        return matchesScope(item, state.scope) && matchesPriority(item, state.priority);
      }

      function groupKeysFor(element) {
        const rawValue = element.dataset.groupKeys || "";
        return rawValue ? rawValue.split(/\s+/).filter(Boolean) : [];
      }

      function stepNumbersFor(element) {
        const rawValue = element.dataset.stepNumbers || "";
        return rawValue
          ? rawValue.split(/\s+/).map(function (value) { return Number(value); }).filter(function (value) { return !Number.isNaN(value); })
          : [];
      }

      function updateButtons(buttons, attributeName, activeValue) {
        buttons.forEach(function (button) {
          const active = button.getAttribute(attributeName) === activeValue;
          button.classList.toggle("is-active", active);
          button.setAttribute("aria-pressed", active ? "true" : "false");
        });
      }

      function filterLabel(filterValue, mode) {
        if (filterValue === "all") {
          return mode === "scope" ? "all scopes" : "all priorities";
        }
        if (filterValue === "manual") {
          return "manual-review scope";
        }
        return filterValue + (mode === "scope" ? " scope" : " priority");
      }

      function sortByRank(items) {
        return items.slice().sort(function (left, right) {
          const priorityDiff = Number(left.dataset.priorityRank || 3) - Number(right.dataset.priorityRank || 3);
          if (priorityDiff !== 0) {
            return priorityDiff;
          }
          const availabilityDiff = Number(left.dataset.availabilityRank || 0) - Number(right.dataset.availabilityRank || 0);
          if (availabilityDiff !== 0) {
            return availabilityDiff;
          }
          return Number(left.dataset.defaultOrder || 0) - Number(right.dataset.defaultOrder || 0);
        });
      }

      function sortByDefaultOrder(items) {
        return items.slice().sort(function (left, right) {
          return Number(left.dataset.defaultOrder || 0) - Number(right.dataset.defaultOrder || 0);
        });
      }

      function formatStepContextLabel(stepNumbers) {
        if (!stepNumbers.length) {
          return "the selected steps";
        }
        if (stepNumbers.length === 1) {
          return "Step " + stepNumbers[0];
        }
        return "Steps " + stepNumbers.join(", ");
      }

      function activeStepNumber() {
        return Number(state.step || 0);
      }

      function relatedGroupKeysForActiveStep() {
        const activeStep = activeStepNumber();
        if (!activeStep) {
          return [];
        }
        return cards
          .filter(function (card) {
            return stepNumbersFor(card).indexOf(activeStep) !== -1;
          })
          .map(function (card) { return card.dataset.groupKey || ""; })
          .filter(Boolean);
      }

      function syncContextCards(items, container, statusNode, noun, activeSteps) {
        if (!container || !statusNode || !items.length) {
          return 0;
        }
        if (!activeSteps.length) {
          sortByDefaultOrder(items).forEach(function (item) {
            item.classList.remove("is-related");
            item.classList.remove("is-dimmed");
            container.appendChild(item);
          });
          statusNode.textContent = "Showing all " + items.length + " " + noun + " in default review order.";
          return 0;
        }
        const relatedItems = [];
        const otherItems = [];
        items.forEach(function (item) {
          const isRelated = stepNumbersFor(item).some(function (stepNumber) {
            return activeSteps.indexOf(stepNumber) !== -1;
          });
          item.classList.toggle("is-related", isRelated);
          item.classList.toggle("is-dimmed", !isRelated);
          if (isRelated) {
            relatedItems.push(item);
            return;
          }
          otherItems.push(item);
        });
        sortByDefaultOrder(relatedItems).concat(sortByDefaultOrder(otherItems)).forEach(function (item) {
          container.appendChild(item);
        });
        if (relatedItems.length) {
          statusNode.textContent = "Promoting " + relatedItems.length + " " + noun + " linked to "
            + formatStepContextLabel(activeSteps) + ".";
        } else {
          statusNode.textContent = "No " + noun + " are directly linked to "
            + formatStepContextLabel(activeSteps) + "; showing the default review order.";
        }
        return relatedItems.length;
      }

      function syncCollection(items, container, emptyState, statusNode, noun) {
        if (!container || !statusNode) {
          return;
        }
        const visibleItems = [];
        const hiddenItems = [];
        items.forEach(function (item) {
          const visible = matchesCurrentFilters(item);
          item.hidden = !visible;
          if (visible) {
            visibleItems.push(item);
          } else {
            hiddenItems.push(item);
          }
        });
        sortByRank(visibleItems).concat(sortByRank(hiddenItems)).forEach(function (item) {
          container.appendChild(item);
        });
        if (emptyState) {
          emptyState.hidden = visibleItems.length !== 0;
        }
        if (visibleItems.length) {
          statusNode.textContent = "Showing " + visibleItems.length + " " + noun + " for "
            + filterLabel(state.scope, "scope") + " and "
            + filterLabel(state.priority, "priority") + ".";
        } else {
          statusNode.textContent = "No " + noun + " match "
            + filterLabel(state.scope, "scope") + " and "
            + filterLabel(state.priority, "priority") + ".";
        }
      }

      function syncPriorityCollection(items, container, statusNode, noun) {
        if (!container || !statusNode) {
          return;
        }
        const useDefaultOrder = state.scope === "all" && state.priority === "all";
        const matchingItems = [];
        const nonMatchingItems = [];
        items.forEach(function (item) {
          item.hidden = false;
          item.classList.remove("is-prioritized");
          if (useDefaultOrder || matchesCurrentFilters(item)) {
            matchingItems.push(item);
          } else {
            nonMatchingItems.push(item);
          }
        });
        const orderedItems = useDefaultOrder
          ? sortByDefaultOrder(items)
          : sortByRank(matchingItems).concat(sortByDefaultOrder(nonMatchingItems));
        orderedItems.forEach(function (item, index) {
          if (!useDefaultOrder && index < matchingItems.length) {
            item.classList.add("is-prioritized");
          }
          container.appendChild(item);
        });
        if (useDefaultOrder) {
          statusNode.textContent = "Showing all " + items.length + " " + noun + " in default review order.";
          return;
        }
        if (matchingItems.length) {
          statusNode.textContent = "Prioritizing " + matchingItems.length + " " + noun + " for "
            + filterLabel(state.scope, "scope") + " and "
            + filterLabel(state.priority, "priority") + ".";
        } else {
          statusNode.textContent = "No " + noun + " directly match "
            + filterLabel(state.scope, "scope") + " and "
            + filterLabel(state.priority, "priority") + "; showing the default order.";
        }
      }

      function applyGroupFocus() {
        const activeGroup = state.group;
        const activeStep = activeStepNumber();
        const activeStepGroupKeys = activeGroup ? [] : relatedGroupKeysForActiveStep();
        const activeCards = cards.filter(function (card) {
          if (activeGroup) {
            return card.dataset.groupKey === activeGroup;
          }
          return Boolean(activeStep) && stepNumbersFor(card).indexOf(activeStep) !== -1;
        });
        cards.forEach(function (card) {
          const isRelated = activeCards.indexOf(card) !== -1;
          card.classList.toggle("is-related", isRelated);
          card.classList.toggle("is-dimmed", Boolean(activeGroup || activeStep) && !isRelated);
          card.setAttribute("aria-pressed", isRelated ? "true" : "false");
        });

        function syncRelated(items, container) {
          let relatedCount = 0;
          const relatedVisible = [];
          const otherVisible = [];
          const hiddenItems = [];
          items.forEach(function (item) {
            const groupKeys = groupKeysFor(item);
            const isRelated = Boolean(activeGroup)
              ? groupKeys.indexOf(activeGroup) !== -1
              : activeStepGroupKeys.some(function (groupKey) { return groupKeys.indexOf(groupKey) !== -1; });
            item.classList.toggle("is-related", isRelated);
            item.classList.toggle("is-dimmed", Boolean(activeGroup || activeStep) && !isRelated);
            if (isRelated && !item.hidden) {
              relatedCount += 1;
            }
            if (item.hidden) {
              hiddenItems.push(item);
            } else if (isRelated) {
              relatedVisible.push(item);
            } else {
              otherVisible.push(item);
            }
          });
          if (container) {
            sortByDefaultOrder(relatedVisible)
              .concat(sortByDefaultOrder(otherVisible))
              .concat(sortByDefaultOrder(hiddenItems))
              .forEach(function (item) {
                container.appendChild(item);
              });
          }
          return relatedCount;
        }

        const relatedActions = syncRelated(actionItems, actionList);
        const relatedCommands = syncRelated(commandItems, commandList);
        const relatedLinks = syncRelated(linkItems, linkList);
        const activeSteps = activeCards.length ? stepNumbersFor(activeCards[0]) : (activeStep ? [activeStep] : []);
        const relatedResultCards = syncContextCards(
          resultCards,
          resultCardList,
          resultStatus,
          "result highlight card(s)",
          activeSteps
        );
        const relatedOperationalCards = syncContextCards(
          operationalCards,
          operationalCardList,
          operationalStatus,
          "operational readiness card(s)",
          activeSteps
        );

        if (clearGroupButton) {
          clearGroupButton.hidden = !activeGroup && !activeStep;
        }
        if (!groupStatus) {
          return;
        }
        if (!activeGroup && !activeStep) {
          groupStatus.textContent = "Select a recovery card to highlight matching actions, commands, links, and step cards.";
          return;
        }
        if (!activeGroup && activeStep) {
          groupStatus.textContent = "Step-linked focus for Step " + activeStep + ": "
            + activeCards.length + " recovery card(s), "
            + relatedActions + " action(s), "
            + relatedCommands + " command suggestion(s), "
            + relatedLinks + " quick link(s), "
            + relatedResultCards + " result card(s), and "
            + relatedOperationalCards + " operational card(s).";
          return;
        }
        if (!activeCards.length) {
          groupStatus.textContent = "Select a recovery card to highlight matching actions, commands, links, and step cards.";
          return;
        }
        const cardTitle = activeCards[0].querySelector("h3");
        const label = cardTitle ? cardTitle.textContent : "selected recovery group";
        groupStatus.textContent = "Linked guidance for " + label + ": "
          + relatedActions + " action(s), "
          + relatedCommands + " command suggestion(s), and "
          + relatedLinks + " quick link(s), "
          + relatedResultCards + " result card(s), and "
          + relatedOperationalCards + " operational card(s).";
      }

      function applyFilters() {
        const activeGroup = state.group;
        const activeStep = activeStepNumber();
        const visibleCards = [];
        const hiddenCards = [];
        cards.forEach(function (card) {
          const visible = matchesCurrentFilters(card);
          if (visible) {
            visibleCards.push(card);
          } else {
            hiddenCards.push(card);
          }
        });
        const relatedVisibleCards = visibleCards.filter(function (card) {
          if (activeGroup) {
            return card.dataset.groupKey === activeGroup;
          }
          return Boolean(activeStep) && stepNumbersFor(card).indexOf(activeStep) !== -1;
        });
        const unrelatedVisibleCards = visibleCards.filter(function (card) {
          return relatedVisibleCards.indexOf(card) === -1;
        });
        const collapseForStep = Boolean(activeStep) && !activeGroup && !state.stepRadarExpanded && relatedVisibleCards.length;
        visibleCards.forEach(function (card) {
          card.hidden = false;
          card.classList.remove("is-prioritized");
        });
        hiddenCards.forEach(function (card) {
          card.hidden = true;
          card.classList.remove("is-prioritized");
        });
        if (collapseForStep) {
          unrelatedVisibleCards.forEach(function (card) {
            card.hidden = true;
          });
        }
        if (radarList) {
          const orderedVisibleCards = (activeGroup || activeStep)
            ? sortByRank(relatedVisibleCards).concat(sortByRank(unrelatedVisibleCards))
            : sortByRank(visibleCards);
          orderedVisibleCards.concat(sortByRank(hiddenCards)).forEach(function (card, index) {
            if ((activeGroup || activeStep) && index < relatedVisibleCards.length) {
              card.classList.add("is-prioritized");
            }
            radarList.appendChild(card);
          });
        }
        updateButtons(scopeButtons, "data-radar-filter", state.scope);
        updateButtons(priorityButtons, "data-radar-priority-filter", state.priority);
        if (toggleStepRadarButton) {
          const canToggleStepView = Boolean(activeStep) && !activeGroup && relatedVisibleCards.length > 0;
          toggleStepRadarButton.hidden = !canToggleStepView;
          toggleStepRadarButton.textContent = collapseForStep
            ? "Show Unrelated Recovery Cards"
            : "Collapse Unrelated Recovery Cards";
        }
        if (activeStep && !activeGroup && relatedVisibleCards.length) {
          const hiddenCount = collapseForStep ? unrelatedVisibleCards.length : 0;
          status.textContent = (collapseForStep ? "Showing " : "Prioritizing ")
            + relatedVisibleCards.length + " recovery group(s) linked to Step " + activeStep
            + " for " + filterLabel(state.scope, "scope") + " and "
            + filterLabel(state.priority, "priority")
            + (hiddenCount ? "; " + hiddenCount + " unrelated group(s) hidden." : ".");
        } else {
          status.textContent = "Showing " + visibleCards.length + " recovery group(s) for "
            + filterLabel(state.scope, "scope") + " and "
            + filterLabel(state.priority, "priority") + ".";
        }
        syncCollection(actionItems, actionList, actionEmpty, actionStatus, "top action(s)");
        syncCollection(commandItems, commandList, commandEmpty, commandStatus, "command suggestion(s)");
        syncPriorityCollection(linkItems, linkList, linkStatus, "quick link(s)");
        applyGroupFocus();
      }

      scopeButtons.forEach(function (button) {
        button.addEventListener("click", function () {
          state.scope = button.dataset.radarFilter || "all";
          state.stepRadarExpanded = false;
          syncHashState(false);
          applyFilters();
        });
      });
      priorityButtons.forEach(function (button) {
        button.addEventListener("click", function () {
          state.priority = button.dataset.radarPriorityFilter || "all";
          state.stepRadarExpanded = false;
          syncHashState(false);
          applyFilters();
        });
      });
      cards.forEach(function (card) {
        function activateCardGroup(toggleIfSame) {
          const nextGroup = card.dataset.groupKey || "";
          state.group = toggleIfSame && state.group === nextGroup ? "" : nextGroup;
          state.step = "";
          syncHashState(false);
          applyFilters();
        }
        card.addEventListener("click", function () {
          activateCardGroup(true);
        });
        card.addEventListener("keydown", function (event) {
          if (event.key !== "Enter" && event.key !== " ") {
            return;
          }
          event.preventDefault();
          activateCardGroup(true);
        });
      });
      if (clearGroupButton) {
        clearGroupButton.addEventListener("click", function () {
          state.group = "";
          state.step = "";
          state.stepRadarExpanded = false;
          syncHashState(false);
          applyFilters();
        });
      }
      if (toggleStepRadarButton) {
        toggleStepRadarButton.addEventListener("click", function () {
          state.stepRadarExpanded = !state.stepRadarExpanded;
          applyFilters();
        });
      }
      sectionLinks.forEach(function (link) {
        link.addEventListener("click", function (event) {
          const targetId = (link.getAttribute("href") || "").replace(/^#/, "");
          if (!targetId) {
            return;
          }
          event.preventDefault();
          state.section = targetId;
          syncHashState(true);
        });
      });
      window.addEventListener("hashchange", function () {
        const hashState = parseHashState();
        state.scope = hashState.scope;
        state.priority = hashState.priority;
        state.section = hashState.section;
        state.group = hashState.group;
        state.step = hashState.step;
        state.stepRadarExpanded = false;
        applyFilters();
        scrollToSection(state.section);
      });
      const initialHashState = parseHashState();
      state.scope = initialHashState.scope;
      state.priority = initialHashState.priority;
      state.section = initialHashState.section;
      state.group = initialHashState.group;
      state.step = initialHashState.step;
      state.stepRadarExpanded = false;
      applyFilters();
      if (state.section) {
        scrollToSection(state.section);
      }
    }());
  </script>"""
    command_html = []
    for item in overview_data.get("annotated_command_suggestions", overview_data["command_suggestions"]):
        source_kind = escape(str(item.get("source_kind", "general")))
        manual_review = "yes" if item.get("manual_review_required") else "no"
        priority_label = escape(str(item.get("priority_label", "low")))
        priority_rank = escape(str(safe_int(item.get("priority_rank"), 3)))
        default_order = escape(str(safe_int(item.get("default_order"), 0)))
        chips = [
            f"<span class=\"tone tone-{source_kind}\">{escape(str(item.get('source_label', 'General')))}</span>",
            f"<span class=\"tone priority-{priority_label}\">{escape(str(item.get('priority_label', 'low')).title())}</span>",
        ]
        if item.get("manual_review_required"):
            chips.append("<span class=\"tone tone-manual\">Manual review</span>")
        group_keys = " ".join(
            str(group_key).strip()
            for group_key in (item.get("group_keys") or [])
            if str(group_key).strip()
        )
        command_html.append(
            f"<li class=\"guided-item guided-command\" "
            f"data-source-kind=\"{source_kind}\" "
            f"data-manual-review=\"{manual_review}\" "
            f"data-priority=\"{priority_label}\" "
            f"data-priority-rank=\"{priority_rank}\" "
            f"data-group-keys=\"{escape(group_keys)}\" "
            f"data-default-order=\"{default_order}\">"
            f"<div class=\"tone-row\">{''.join(chips)}</div>"
            f"<strong>{escape(item['label'])}</strong><br>"
            f"<code>{escape(item['command'])}</code><br>"
            f"<span>{escape(item['reason'])}</span>"
            "</li>"
        )
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Run Overview - {escape(str(overview_data['project_name']))}</title>
  <style>
    :root {{
      --bg: #f3efe7;
      --surface: #fffdf8;
      --ink: #1d1b18;
      --muted: #615b53;
      --line: #d2c7b8;
      --accent: #0f6d58;
      --warn: #9b5c00;
      --danger: #9c2f1f;
      --shadow: 0 18px 40px rgba(38, 31, 21, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: "Segoe UI", Tahoma, sans-serif;
      background: radial-gradient(circle at top, #faf7f0 0%, var(--bg) 45%, #ece4d8 100%);
      color: var(--ink);
      line-height: 1.5;
    }}
    main {{
      max-width: 1120px;
      margin: 0 auto;
      padding: 32px 20px 48px;
    }}
    header, section {{
      background: var(--surface);
      border: 1px solid var(--line);
      border-radius: 20px;
      box-shadow: var(--shadow);
      padding: 24px;
      margin-bottom: 20px;
    }}
    h1, h2, h3 {{ margin-top: 0; }}
    .lede {{ color: var(--muted); font-size: 1.05rem; }}
    .stats {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 12px;
      margin-top: 20px;
    }}
    .section-nav {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-top: 18px;
    }}
    .jump-link {{
      display: inline-flex;
      align-items: center;
      padding: 8px 12px;
      border-radius: 999px;
      border: 1px solid var(--line);
      background: #fbf6ee;
      color: #0a4d8c;
      text-decoration: none;
      font-weight: 600;
    }}
    .stat {{
      border: 1px solid var(--line);
      border-radius: 16px;
      padding: 14px;
      background: #fcfaf5;
    }}
    .stat strong {{
      display: block;
      font-size: 1.3rem;
    }}
    .progress {{
      height: 16px;
      background: #e5ddd0;
      border-radius: 999px;
      overflow: hidden;
      margin: 16px 0 8px;
    }}
    .bar {{
      height: 100%;
      width: {overview_data['progress_percent']}%;
      background: linear-gradient(90deg, #3c947d, var(--accent));
    }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      gap: 14px;
    }}
    .tone-row {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin-bottom: 12px;
    }}
    .tone {{
      display: inline-flex;
      align-items: center;
      border-radius: 999px;
      border: 1px solid var(--line);
      padding: 4px 10px;
      font-size: 0.8rem;
      font-weight: 700;
      background: #f7f1e8;
      color: var(--ink);
    }}
    .tone-summary {{
      background: #f7f1e8;
    }}
    .tone-validation {{
      background: #f9ece8;
      border-color: #d9b0a6;
      color: #8a2f1f;
    }}
    .tone-operational {{
      background: #ebf4ef;
      border-color: #9cc8b6;
      color: #0f6d58;
    }}
    .tone-action {{
      background: #f2ece3;
    }}
    .tone-step {{
      background: #eef3f9;
      border-color: #b9c9db;
      color: #244b73;
    }}
    .tone-manual {{
      background: #fff4de;
      border-color: #d8ba75;
      color: #8a5a00;
    }}
    .priority-immediate {{
      background: #f9e3df;
      border-color: #d8a19a;
      color: var(--danger);
    }}
    .priority-high {{
      background: #fff0d8;
      border-color: #d9b36b;
      color: var(--warn);
    }}
    .priority-medium, .priority-low {{
      background: #eef3ea;
      border-color: #b9c8b2;
      color: #49603f;
    }}
    .card {{
      border: 1px solid var(--line);
      border-radius: 16px;
      padding: 16px;
      background: #fff;
    }}
    .radar-card {{
      border-width: 2px;
      cursor: pointer;
    }}
    .radar-card:focus {{
      outline: 3px solid rgba(10, 77, 140, 0.28);
      outline-offset: 3px;
    }}
    .radar-empty {{
      color: var(--muted);
      margin: 0;
    }}
    .radar-filter-stack {{
      display: grid;
      gap: 10px;
      margin: 16px 0 8px;
    }}
    .radar-filter-group {{
      display: grid;
      gap: 6px;
    }}
    .filter-label {{
      margin: 0;
      color: var(--muted);
      font-size: 0.82rem;
      font-weight: 700;
      letter-spacing: 0.04em;
      text-transform: uppercase;
    }}
    .radar-filters {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }}
    .filter-chip {{
      border: 1px solid var(--line);
      border-radius: 999px;
      background: #fbf6ee;
      color: var(--ink);
      padding: 8px 12px;
      font: inherit;
      font-weight: 600;
      cursor: pointer;
    }}
    .filter-chip.is-active {{
      background: #e3f0ea;
      border-color: #8db9aa;
      color: var(--accent);
    }}
    .filter-status {{
      color: var(--muted);
      margin-top: 0;
    }}
    .guided-list {{
      margin: 0;
      padding-left: 20px;
      display: grid;
      gap: 12px;
    }}
    .guided-item {{
      border: 1px solid var(--line);
      border-radius: 16px;
      background: #fcfaf5;
      padding: 14px 16px;
    }}
    .card.is-prioritized, .guided-item.is-prioritized {{
      border-color: #9cc8b6;
      background: #fffdf8;
      box-shadow: inset 0 0 0 1px rgba(15, 109, 88, 0.08);
    }}
    .card.is-related, .guided-item.is-related {{
      border-color: #0a4d8c;
      box-shadow: 0 0 0 2px rgba(10, 77, 140, 0.16), var(--shadow);
      background: #fffefa;
    }}
    .card.is-dimmed, .guided-item.is-dimmed {{
      opacity: 0.45;
    }}
    .guided-item p {{
      margin: 8px 0 0;
    }}
    .guided-empty {{
      color: var(--muted);
      margin-bottom: 0;
    }}
    .card .path {{
      color: var(--muted);
      margin-bottom: 0;
    }}
    .badge {{
      display: inline-block;
      padding: 4px 10px;
      border-radius: 999px;
      border: 1px solid currentColor;
      font-size: 0.85rem;
      font-weight: 600;
      text-transform: capitalize;
      background: #fff;
    }}
    .badge-completed, .badge-complete {{ color: var(--accent); }}
    .badge-running {{ color: #0b4ea2; }}
    .badge-partial, .badge-stale, .badge-skipped {{ color: var(--warn); }}
    .badge-failed, .badge-error {{ color: var(--danger); }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.95rem;
    }}
    th, td {{
      border-bottom: 1px solid var(--line);
      text-align: left;
      vertical-align: top;
      padding: 10px 8px;
    }}
    th {{
      color: var(--muted);
      font-weight: 700;
    }}
    code {{
      background: #f3eee6;
      padding: 1px 5px;
      border-radius: 6px;
    }}
    a {{ color: #0a4d8c; }}
    @media (max-width: 720px) {{
      header, section {{ padding: 18px; }}
      table {{ font-size: 0.88rem; }}
      .section-nav {{ gap: 8px; }}
    }}
  </style>
</head>
<body>
  <main>
    <header>
      <h1>Run Overview</h1>
      <h2>{escape(str(overview_data['hero_title']))}</h2>
      <p class="lede">{escape(str(overview_data['hero_body']))}</p>
      <div class="progress" role="progressbar" aria-valuemin="0" aria-valuemax="100" aria-valuenow="{overview_data['progress_percent']}" aria-label="Pipeline progress">
        <div class="bar"></div>
      </div>
      <p>{escape(str(overview_data['resolved_count']))} of {escape(str(overview_data['total_count']))} phases resolved ({escape(str(overview_data['progress_percent']))}%).</p>
      <div class="stats">
        <div class="stat"><span>Overall status</span><strong>{escape(_status_label(str(overview_data['overall_status'])))}</strong></div>
        <div class="stat"><span>Execution mode</span><strong>{escape(str(overview_data['execution_mode']))}</strong></div>
        <div class="stat"><span>Current phase</span><strong>{escape(str(overview_data['current_phase_name'] or 'n/a'))}</strong></div>
        <div class="stat"><span>Last update</span><strong>{escape(str(overview_data['updated_at']))}</strong></div>
      </div>
      <nav class="section-nav" aria-label="Overview sections">{jump_links_html}</nav>
    </header>
    <section id="what-to-check-first">
      <h2>What To Check First</h2>
      <p id="action-filter-status" class="filter-status" role="status">Showing {escape(str(len(action_html)))} top action(s) in urgency order.</p>
      <ol id="guided-action-list" class="guided-list">{''.join(action_html)}</ol>
      <p id="action-filter-empty" class="guided-empty" hidden>No top actions match the active recovery radar filter.</p>
    </section>
    <section id="recovery-radar">
      <h2>Recovery Radar</h2>
      <p class="lede">{escape(str(overview_data['recovery_radar']['headline']))}</p>
      <div class="tone-row">{radar_summary_html}</div>
      <div class="tone-row">
        <button type="button" id="clear-group-focus" class="filter-chip" hidden>Clear Linked Focus</button>
        <button type="button" id="toggle-step-radar-visibility" class="filter-chip" hidden>Show Unrelated Recovery Cards</button>
      </div>
      <p id="group-focus-status" class="filter-status" role="status">Select a recovery card to highlight matching actions, commands, links, and step cards.</p>
      {radar_filter_controls_html}
      {f'<div id="recovery-radar-list" class="cards">{"".join(recovery_radar_html)}</div>' if recovery_radar_html else '<p class="radar-empty">No recovery action groups are active right now.</p>'}
    </section>
    <section id="result-highlights">
      <h2>Result Highlights</h2>
      <p id="result-context-status" class="filter-status" role="status">Showing all {escape(str(len(cards_html)))} result highlight card(s) in default review order.</p>
      <div id="result-card-list" class="cards">{''.join(cards_html)}</div>
    </section>
    <section id="operational-readiness">
      <h2>Operational Readiness</h2>
      <p id="operational-context-status" class="filter-status" role="status">Showing all {escape(str(len(operational_cards_html)))} operational readiness card(s) in default review order.</p>
      <div id="operational-card-list" class="cards">{''.join(operational_cards_html)}</div>
    </section>
    <section id="pipeline-progress">
      <h2>Pipeline Progress</h2>
      <table>
        <thead>
          <tr>
            <th scope="col">Phase</th>
            <th scope="col">Status</th>
            <th scope="col">Name</th>
            <th scope="col">Started</th>
            <th scope="col">Finished</th>
            <th scope="col">Notes</th>
          </tr>
        </thead>
        <tbody>{''.join(phase_rows_html)}</tbody>
      </table>
    </section>
    <section id="quick-links">
      <h2>Quick Links</h2>
      <p id="quick-link-filter-status" class="filter-status" role="status">Showing {escape(str(len(quick_links_html)))} quick link(s) in default review order.</p>
      <ul id="prioritized-quick-links" class="guided-list">{''.join(quick_links_html)}</ul>
    </section>
    <section id="suggested-commands">
      <h2>Suggested Commands</h2>
      <p id="command-filter-status" class="filter-status" role="status">Showing {escape(str(len(command_html)))} command suggestion(s) in urgency order.</p>
      <ul id="guided-command-list" class="guided-list">{''.join(command_html)}</ul>
      <p id="command-filter-empty" class="guided-empty" hidden>No command suggestions match the active recovery radar filter.</p>
    </section>
  </main>
{radar_filter_script}
</body>
</html>
"""
    return _atomic_write_text(project_root / "run_overview.html", html)


def write_report_digest(project_root: Union[Path, str], overview_data: dict) -> Path:
    """Write a compact `report_digest.md` for fast report-first review."""

    project_root = _as_path(project_root)
    operational_recovery = overview_data["operational_recovery_summary"]
    recovery_radar = overview_data["recovery_radar"]
    validation_summary = overview_data["validation_summary"]
    lines = [
        "# Report Digest",
        "",
        "## At a Glance",
        f"- {overview_data['hero_title']}",
        f"- Pocket summary: {overview_data['pocket_summary']['headline']}",
        f"- Verdict summary: {overview_data['verdict_summary']['headline']}",
        f"- Report summary: {overview_data['report_summary']['headline']}",
        f"- Validation summary: {overview_data['validation_summary']['headline']}",
        "",
        "## Recovery Snapshot",
        f"- {recovery_radar['headline']}",
    ]
    lines.extend(["", "### Filter Views"])
    for filter_view in recovery_radar["filter_views"]:
        lines.append(f"- {filter_view['label']}: {filter_view['headline']}")
    if recovery_radar["items"]:
        for index, item in enumerate(recovery_radar["items"], start=1):
            command_text = item["recommended_command"] or "manual review only"
            lines.extend(
                [
                    f"{index}. {item['source_label']} - {item['title']} ({item['scope_label']})",
                    f"Action: {item['action_label']}; Priority: {item['priority_label']}; Command: `{command_text}`",
                    f"Playbook: `{item.get('playbook_link', item['playbook_path'])}`",
                ]
            )
            if item.get("step_summary_links"):
                lines.append(f"Step summaries: {_step_summary_links_text(item['step_summary_links'])}")
        if recovery_radar["summary"]["overflow_count"]:
            lines.append(
                f"Additional recovery groups: {recovery_radar['summary']['overflow_count']} more in `run_overview.md` or `run_overview.html`."
            )
    else:
        lines.append("- No recovery action groups are active right now.")

    lines.extend(
        [
            "",
        "## Narrative Highlights",
        ]
    )
    report_details = overview_data["report_summary"]["details"] or ["Report digest is not available yet."]
    lines.extend(f"- {detail}" for detail in report_details)
    lines.extend(["", "## Priority Findings"])
    verdict_details = overview_data["verdict_summary"]["details"] or ["Verdict details are not available yet."]
    lines.extend(f"- {detail}" for detail in verdict_details)
    lines.extend(["", "## Operational Follow-up", f"- {operational_recovery['headline']}"])
    if operational_recovery.get("step_summary_links"):
        lines.append(f"- Step summaries: {_step_summary_links_text(operational_recovery['step_summary_links'])}")
    operational_details = operational_recovery["details"] or ["Operational recovery follow-up is not available yet."]
    lines.extend(f"- {detail}" for detail in operational_details)
    lines.extend(["", "## Validation Follow-up", f"- {validation_summary['headline']}"])
    if validation_summary.get("step_summary_links"):
        lines.append(f"- Step summaries: {_step_summary_links_text(validation_summary['step_summary_links'])}")
    validation_details = validation_summary["details"] or ["Validation follow-up is not available yet."]
    lines.extend(f"- {detail}" for detail in validation_details)
    lines.extend(["", "## Suggested Commands"])
    lines.extend(
        f"- `{item['command']}`: {item['label']} - {item['reason']}"
        for item in overview_data["command_suggestions"]
    )
    lines.extend(
        [
            "",
            "## Open Next",
            "1. `step6_report/project_report.txt`",
            "2. `step5_verdict/valid_sites.csv`",
            "3. `step4_vina_postprocess/vina_pocket_table.csv`",
            "4. `operational_recovery_playbook.md`",
            "5. `step7_validate/validation_recovery_playbook.md`",
            "6. `step7_validate/validation_summary.txt`",
            "",
        ]
    )
    return _atomic_write_text(project_root / "report_digest.md", "\n".join(lines))


def update_run_overview(
    project_root: Union[Path, str],
    *,
    current_run_manifest: Optional[dict] = None,
    run_status: Optional[dict] = None,
) -> Tuple[Path, Path]:
    """Refresh Markdown and HTML overview files for the project."""

    project_root = _as_path(project_root)
    operational_recovery_plan = build_operational_recovery_plan(
        project_root,
        current_manifest=current_run_manifest,
    )
    validation_recovery_plan = _load_validation_recovery_plan(project_root)
    write_operational_recovery_outputs(project_root, operational_recovery_plan)
    overview_data = build_run_overview_data(
        project_root,
        current_run_manifest=current_run_manifest,
        run_status=run_status,
        operational_recovery_plan=operational_recovery_plan,
        validation_recovery_plan=validation_recovery_plan,
        expected_paths=[
            "report_digest.md",
            "operational_recovery_playbook.md",
            "operational_recovery_plan.json",
        ],
    )
    markdown_path = write_run_overview_markdown(project_root, overview_data)
    html_path = write_run_overview_html(project_root, overview_data)
    write_report_digest(project_root, overview_data)
    return markdown_path, html_path


def write_step_index(project_root: Union[Path, str], index_data: dict) -> Path:
    """Write the project-level ``step_index.md`` entry point."""

    project_root = _as_path(project_root)
    run_summary = index_data["run_summary"]
    steps = index_data["steps"]
    triage_by_step = index_data["triage_by_step"]
    raw_debug_paths = index_data["raw_debug_paths"]
    notes = index_data["notes"]

    lines = [
        "# Step Output Index",
        "",
        "## Run Summary",
        f"- Project name: `{run_summary['project_name']}`",
        f"- Generated at: `{run_summary['generated_at']}`",
        f"- Execution mode: `{run_summary['execution_mode']}`",
        f"- Receptor IDs: `{', '.join(run_summary['receptor_ids']) or 'n/a'}`",
        f"- Ligand IDs: `{', '.join(run_summary['ligand_ids']) or 'n/a'}`",
        f"- Config paths: `{', '.join(run_summary['config_paths']) or 'n/a'}`",
        "",
        "## Step Overview Table",
        "| Step | Folder | Purpose | Status | Step Summary | Triage | Inspect First |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for item in steps:
        lines.append(
            f"| [{item['step_number']}]({item['overview_focus_link']}) | `{item['folder_name']}` | {item['purpose']} | "
            f"`{item['status']}` | `{item['summary_path']}` | {item['triage_summary']} | {item['primary_files']} |"
        )

    lines.extend(["", "## Step Triage"])
    if not any(triage_by_step.values()):
        lines.extend(
            [
                "- No active recovery groups are linked to any step right now.",
                "- Use `run_overview.md` or `run_overview.html` for the clean run summary.",
            ]
        )
    else:
        for item in steps:
            triage_items = triage_by_step.get(item["step_number"]) or []
            lines.extend(["", f"### Step {item['step_number']}: {_step_spec(item['step_number']).step_name}"])
            if not triage_items:
                lines.append("- No active recovery groups are linked to this step.")
                continue
            lines.append(f"- Active recovery groups: {len(triage_items)}.")
            for index, triage_item in enumerate(triage_items, start=1):
                command_text = str(triage_item.get("recommended_command", "")).strip() or "manual review only"
                lines.append(
                    f"- Group {index}: {triage_item['source_label']} - {triage_item['title']} "
                    f"(`{triage_item['priority_label']}`, {triage_item['scope_label']})."
                )
                lines.append(f"- Overview link: `{triage_item['overview_link']}`")
                lines.append(f"- Playbook: `{triage_item['playbook_link']}`")
                lines.append(f"- Suggested command: `{command_text}`")
                if triage_item.get("step_summary_links"):
                    lines.append(
                        f"- Related step summaries: {_step_summary_links_text(triage_item['step_summary_links'])}"
                    )

    lines.extend(
        [
            "",
            "## Where To Read First",
            f"1. `{STEP_SPECS[6].folder_name}/project_report.txt` - Narrative summary ({index_data['read_first'][0]}).",
            f"2. `{STEP_SPECS[5].folder_name}/valid_sites.csv` - Prioritized sites ({index_data['read_first'][1]}).",
            f"3. `{STEP_SPECS[4].folder_name}/vina_pocket_table.csv` - Pocket-level interpretation ({index_data['read_first'][2]}).",
            f"4. `{STEP_SPECS[3].folder_name}/ppi_pyrosetta_residues.csv` - Receptor residue evidence ({index_data['read_first'][3]}).",
            "",
            "## Raw Debug Paths",
            f"- Canonical project root: `{raw_debug_paths['project_root']}`",
        ]
    )
    if raw_debug_paths["pyrosetta_raw_run_paths"]:
        lines.append(
            "- PyRosetta raw run paths: "
            + ", ".join(f"`{path}`" for path in raw_debug_paths["pyrosetta_raw_run_paths"])
        )
    else:
        lines.append("- PyRosetta raw run paths: `not recorded`")
    if raw_debug_paths["receptor_docking_paths"]:
        lines.append(
            "- Receptor docking paths: "
            + ", ".join(f"`{path}`" for path in raw_debug_paths["receptor_docking_paths"])
        )
    else:
        lines.append("- Receptor docking paths: `not recorded`")

    lines.extend(["", "## Notes and Warnings"])
    if notes:
        lines.extend(f"- {note}" for note in notes)
    else:
        lines.append("- No current warnings.")
    lines.append("")

    return _atomic_write_text(project_root / "step_index.md", "\n".join(lines))


def update_step_index(
    project_root: Union[Path, str],
    *,
    current_run_manifest: Optional[dict] = None,
    notes: Optional[Sequence[str]] = None,
) -> Path:
    """Create or refresh the project-level step index."""

    project_root = _as_path(project_root)
    current = current_run_manifest or _load_current_run_manifest(project_root) or {}
    step_status = current.get("step_status") or _step_status_map(project_root)
    overview_data = build_run_overview_data(project_root, current_run_manifest=current)
    triage_by_step = _step_triage_items(overview_data["recovery_radar"])
    steps = []
    for step_num, spec in STEP_SPECS.items():
        triage_items = triage_by_step.get(step_num) or []
        steps.append(
            {
                "step_number": step_num,
                "folder_name": spec.folder_name,
                "purpose": spec.purpose,
                "status": step_status.get(f"step{step_num}", "not_generated"),
                "summary_path": _step_summary_path(step_num),
                "overview_focus_link": _overview_step_link(step_num),
                "triage_summary": _step_triage_summary(triage_items),
                "primary_files": _step_primary_paths(spec),
            }
        )

    read_first_paths = [
        project_root / STEP_SPECS[6].folder_name / "project_report.txt",
        project_root / STEP_SPECS[5].folder_name / "valid_sites.csv",
        project_root / STEP_SPECS[4].folder_name / "vina_pocket_table.csv",
        project_root / STEP_SPECS[3].folder_name / "ppi_pyrosetta_residues.csv",
    ]
    read_first = ["available" if path.exists() else "missing" for path in read_first_paths]

    receptor_docking_paths = []
    for receptor_id in current.get("receptors", []):
        receptor_path = project_root / receptor_id
        if receptor_path.exists():
            receptor_docking_paths.append(_display_path(receptor_path, project_root.parent))

    note_lines = list(notes or [])
    if current.get("execution_mode", "full") != "full":
        note_lines.append(
            "Partial rerun detected. Compare step manifest timestamps before mixing interpretations."
        )
    if any(item.get("status") == "stale" for item in steps):
        note_lines.append(
            "Some derived steps are stale for the current rerun scope. Rebuild them before mixing interpretations."
        )
    if (project_root / "current_run_manifest.json").exists():
        note_lines.append("Current run manifest: `current_run_manifest.json`.")
    note_lines.append("User-facing overview: `run_overview.md` and `run_overview.html`.")
    if (project_root / "report_digest.md").exists():
        note_lines.append("Condensed report digest: `report_digest.md`.")
    if (project_root / "operational_recovery_playbook.md").exists():
        note_lines.append("Operational recovery playbook: `operational_recovery_playbook.md`.")
    validation_path = project_root / STEP_SPECS[7].folder_name / "validation_status.json"
    if validation_path.exists():
        note_lines.append(
            f"Validation summary available at `{STEP_SPECS[7].folder_name}/validation_status.json`."
        )
    recovery_playbook_path = project_root / STEP_SPECS[7].folder_name / "validation_recovery_playbook.md"
    if recovery_playbook_path.exists():
        note_lines.append(
            f"Validation recovery playbook available at `{STEP_SPECS[7].folder_name}/validation_recovery_playbook.md`."
        )

    index_data = {
        "run_summary": {
            "project_name": current.get("project_name", project_root.name),
            "generated_at": current.get("generated_at", utc_now_iso()),
            "execution_mode": current.get("execution_mode", "full"),
            "receptor_ids": list(current.get("receptors", [])),
            "ligand_ids": list(current.get("ligands", [])),
            "config_paths": [current.get("vina_config_path", "n/a"), *current.get("ppi_config_paths", [])],
        },
        "steps": steps,
        "triage_by_step": triage_by_step,
        "read_first": read_first,
        "raw_debug_paths": {
            "project_root": current.get("project_root", _display_path(project_root, project_root.parent)),
            "pyrosetta_raw_run_paths": list(current.get("pyrosetta_raw_run_paths", [])),
            "receptor_docking_paths": receptor_docking_paths,
        },
        "notes": note_lines,
    }
    return write_step_index(project_root, index_data)


def refresh_root_step_views(
    config_path: Union[Path, str],
    *,
    repo_root: Optional[Union[Path, str]] = None,
    execution_mode: str = "full",
    fresh_run: bool = False,
    stale_steps: Optional[Sequence[int]] = None,
    ppi_config_paths: Optional[Sequence[str]] = None,
    pyrosetta_raw_run_paths: Optional[Sequence[str]] = None,
    notes: Optional[Sequence[str]] = None,
) -> Tuple[Path, Path]:
    """Refresh the root manifest and step index after step generation."""

    config_path, repo_root_path, _, project_root = _collect_base_context(
        config_path,
        repo_root=repo_root,
        create_project_root=False,
    )
    manifest = build_current_run_manifest(
        config_path,
        repo_root=repo_root_path,
        execution_mode=execution_mode,
        fresh_run=fresh_run,
        stale_steps=stale_steps,
        ppi_config_paths=ppi_config_paths,
        pyrosetta_raw_run_paths=pyrosetta_raw_run_paths,
    )
    manifest_path = write_current_run_manifest(project_root, manifest)
    index_path = update_step_index(project_root, current_run_manifest=manifest, notes=notes)
    update_run_overview(project_root, current_run_manifest=manifest)
    return manifest_path, index_path


__all__ = [
    "STEP_SPECS",
    "build_current_run_manifest",
    "clear_step_views",
    "copy_artifact_if_exists",
    "ensure_step_dir",
    "record_step1_outputs",
    "record_step2_outputs",
    "record_step3_outputs",
    "record_step4_outputs",
    "record_step5_outputs",
    "record_step6_outputs",
    "record_step7_outputs",
    "refresh_root_step_views",
    "step_output_view_enabled",
    "resolve_project_root",
    "update_run_overview",
    "update_step_index",
    "write_validation_outputs",
    "write_artifact_index",
    "write_current_run_manifest",
    "write_run_overview_html",
    "write_run_overview_markdown",
    "write_run_status",
    "write_step_index",
    "write_step_manifest",
    "write_step_summary",
]
