#!/usr/bin/env python3
"""Rebuild additive step-view outputs from existing canonical project outputs."""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path
from typing import Callable, Dict, List


REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from egfr_pipeline.step_view import (
    record_step1_outputs,
    record_step2_outputs,
    record_step3_outputs,
    record_step4_outputs,
    record_step5_outputs,
    record_step6_outputs,
    record_step7_outputs,
    refresh_root_step_views,
)
from egfr_pipeline.validate import run_validation
from run_production import PPI_TARGETS


DEFAULT_CONFIG = REPO_ROOT / "config" / "example-project.yaml"


def _parse_steps(value: str) -> List[int]:
    steps = sorted({int(item.strip()) for item in value.split(",") if item.strip()})
    invalid = [step for step in steps if step < 1 or step > 7]
    if invalid:
        raise argparse.ArgumentTypeError(f"Invalid step numbers: {invalid}")
    return steps


def _record_step2(config_path: Path) -> Path:
    return record_step2_outputs(
        config_path,
        repo_root=REPO_ROOT,
        ppi_targets=PPI_TARGETS,
    )


STEP_BUILDERS: Dict[int, Callable[[Path], Path]] = {
    1: lambda config_path: record_step1_outputs(config_path, repo_root=REPO_ROOT),
    2: _record_step2,
    3: lambda config_path: record_step3_outputs(config_path, repo_root=REPO_ROOT),
    4: lambda config_path: record_step4_outputs(config_path, repo_root=REPO_ROOT),
    5: lambda config_path: record_step5_outputs(config_path, repo_root=REPO_ROOT),
    6: lambda config_path: record_step6_outputs(config_path, repo_root=REPO_ROOT),
}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rebuild derived step-view outputs from existing canonical outputs.",
    )
    parser.add_argument(
        "-c",
        "--config",
        default=str(DEFAULT_CONFIG),
        help="Project config used to resolve the canonical output root.",
    )
    parser.add_argument(
        "--steps",
        type=_parse_steps,
        default=[1, 2, 3, 4, 5, 6, 7],
        help="Comma-separated step numbers to rebuild. Defaults to 1,2,3,4,5,6,7.",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    selected_steps = list(args.steps)
    failures: List[str] = []

    print(f"Rebuilding step views from canonical outputs: {config_path}")
    print(f"Selected steps: {selected_steps}")
    print("")

    for step_num in selected_steps:
        if step_num == 7:
            continue
        try:
            step_dir = STEP_BUILDERS[step_num](config_path)
            print(f"[OK] step{step_num}: {step_dir}")
        except Exception as exc:  # pragma: no cover - exercised via real runs
            failures.append(f"step{step_num}: {exc}")
            print(f"[WARN] step{step_num}: {exc}")
            traceback.print_exc()

    if 7 in selected_steps:
        try:
            validation_result = run_validation(str(config_path), repo_root=str(REPO_ROOT))
            step7_dir = record_step7_outputs(
                config_path,
                repo_root=REPO_ROOT,
                validation_result=validation_result,
            )
            print(f"[OK] step7: {step7_dir}")
        except Exception as exc:  # pragma: no cover - exercised via real runs
            failures.append(f"step7: {exc}")
            step7_dir = record_step7_outputs(
                config_path,
                repo_root=REPO_ROOT,
                error_message=str(exc),
            )
            print(f"[WARN] step7 validation failed, persisted error state: {step7_dir}")
            traceback.print_exc()

    manifest_path, index_path = refresh_root_step_views(
        config_path,
        repo_root=REPO_ROOT,
        execution_mode=f"backfill:{','.join(str(step) for step in selected_steps)}",
        ppi_config_paths=[target["config_ini"] for target in PPI_TARGETS],
    )
    print("")
    print(f"Run manifest: {manifest_path}")
    print(f"Step index:   {index_path}")

    if failures:
        print("")
        print("Completed with warnings:")
        for item in failures:
            print(f"- {item}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
