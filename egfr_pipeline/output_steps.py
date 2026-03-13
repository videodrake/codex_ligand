"""Derived step-view helpers for production outputs.

This module builds the additive step output layer described in the planning
docs. Canonical runtime outputs under ``output/{project}`` remain the source
of truth; step folders are curated copies and indexes for interpretation.
"""

from __future__ import annotations

import configparser
import csv
import json
import shutil
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union

from egfr_pipeline.config import load_config
from egfr_pipeline.pyrosetta_docking.metadata import build_output_root_name


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
        primary_files=("TH1_final_ranking.csv", "beta_meander_final_ranking.csv"),
        required_artifacts=(
            "TH1_final_ranking.csv",
            "beta_meander_final_ranking.csv",
            "raw_run_paths.tsv",
        ),
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
            "This step converts raw PPI docking into receptor residue evidence "
            "that can be compared with Vina pockets."
        ),
        primary_files=("ppi_pyrosetta_residues.csv", "ppi_pyrosetta_summary.csv"),
        required_artifacts=("ppi_pyrosetta_residues.csv", "ppi_pyrosetta_summary.csv"),
        optional_artifacts=("phase1_interface_report.md",),
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
            "mappings, and cross-receptor comparisons."
        ),
        primary_files=("vina_pocket_table.csv", "vina_pose_table.csv"),
        required_artifacts=(
            "vina_pose_table.csv",
            "vina_pocket_table.csv",
            "vina_drug_pocket_map.csv",
            "vina_pocket_comparison.csv",
            "vina_pocket_bootstrap.csv",
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
        purpose="Evidence strength classification for candidate sites.",
        summary_description=(
            "This step combines Vina pocket evidence with PPI support to "
            "prioritize sites for interpretation."
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
            "This step is the report-first view for reading a completed run."
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
        upstream_steps=(1, 2, 3, 4, 5, 6),
        next_step_reads=(),
    ),
}


STEP_NAME_TO_FOLDER = {
    spec.step_name: spec.folder_name for spec in STEP_SPECS.values()
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
    return None


def _phase1_interface_report_path(repo_root: Path) -> Optional[Path]:
    report_path = repo_root / "output" / "phase1_ppi" / "phase1_interface_report.md"
    return report_path if report_path.exists() else None


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
                "path": "TH1_final_ranking.csv",
                "description": "TH1 partner ranking summary.",
                "status": "missing" if "TH1_final_ranking.csv" in missing_required else "",
            },
            {
                "path": "beta_meander_final_ranking.csv",
                "description": "Beta-meander partner ranking summary.",
                "status": "missing" if "beta_meander_final_ranking.csv" in missing_required else "",
            },
            {
                "path": "raw_run_paths.tsv",
                "description": "Path index back to the canonical raw PyRosetta run directories.",
            },
        ]
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

        notes = "Step 3 reflects the residue evidence files used by the current report and verdict logic."
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

    with _staged_step_dir(project_root, 7) as (temp_dir, _, existed):
        status_path, summary_path = write_validation_outputs(
            temp_dir,
            validation_status,
            summary_lines,
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
        for name in ("current_run_manifest.json", "step_index.md"):
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


def _step_primary_paths(spec: StepSpec) -> str:
    return ", ".join(f"`{spec.folder_name}/{name}`" for name in spec.primary_files) or "-"


def write_step_index(project_root: Union[Path, str], index_data: dict) -> Path:
    """Write the project-level ``step_index.md`` entry point."""

    project_root = _as_path(project_root)
    run_summary = index_data["run_summary"]
    steps = index_data["steps"]
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
        "| Step | Folder | Purpose | Status | Inspect First |",
        "| --- | --- | --- | --- | --- |",
    ]
    for item in steps:
        lines.append(
            f"| {item['step_number']} | `{item['folder_name']}` | {item['purpose']} | "
            f"`{item['status']}` | {item['primary_files']} |"
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
    steps = []
    for step_num, spec in STEP_SPECS.items():
        steps.append(
            {
                "step_number": step_num,
                "folder_name": spec.folder_name,
                "purpose": spec.purpose,
                "status": step_status.get(f"step{step_num}", "not_generated"),
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
    validation_path = project_root / STEP_SPECS[7].folder_name / "validation_status.json"
    if validation_path.exists():
        note_lines.append(
            f"Validation summary available at `{STEP_SPECS[7].folder_name}/validation_status.json`."
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
    "update_step_index",
    "write_validation_outputs",
    "write_artifact_index",
    "write_current_run_manifest",
    "write_step_index",
    "write_step_manifest",
    "write_step_summary",
]
