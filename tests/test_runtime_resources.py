from __future__ import annotations

import json
from pathlib import Path

from egfr_pipeline.runtime import (
    lane_is_complete,
    resolve_runtime_resources,
    write_lane_completion,
)


def test_resolve_runtime_resources_prefers_pbs_and_detects_visible_gpus() -> None:
    runtime = resolve_runtime_resources(
        requested_cpus=32,
        environ={
            "PBS_NP": "24",
            "CUDA_VISIBLE_DEVICES": "0,2",
            "PBS_JOBID": "123.server",
            "HOSTNAME": "node05",
        },
    )

    assert runtime.allocated_cpus == 24
    assert runtime.effective_cpus == 24
    assert runtime.visible_gpu_ids == ["0", "2"]
    assert runtime.pbs_job_id == "123.server"
    assert runtime.host == "node05"
    assert runtime.source == "pbs_np"


def test_lane_is_complete_requires_completed_marker(tmp_path: Path) -> None:
    project_root = tmp_path / "output" / "demo"
    write_lane_completion(project_root, "vina-cpu", status="failed", engine="cpu-vina")
    assert lane_is_complete(project_root, "vina-cpu", engine="cpu-vina") is False

    write_lane_completion(project_root, "vina-cpu", status="completed", engine="cpu-vina")
    assert lane_is_complete(project_root, "vina-cpu", engine="cpu-vina") is True

    manifest_path = (
        project_root
        / "runtime"
        / "lanes"
        / "vina-cpu__cpu-vina"
        / "runtime_manifest.json"
    )
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["status"] == "completed"
    assert payload["completed_at"]
