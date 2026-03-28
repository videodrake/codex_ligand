import json
from pathlib import Path

import pytest

from egfr_pipeline.step_view import record_step2_outputs

pytestmark = pytest.mark.integration


def test_record_step2_outputs_emits_legacy_warning_codes(tmp_path: Path) -> None:
    repo_root = tmp_path
    config_path = repo_root / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "project_name": "demo_project",
                "output_root": "./output",
                "step_output_view": {"root_mode": "workflow_a_step_views"},
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    run_dir = repo_root / "output" / "workflow_a" / "phase2_ppi_docking" / "active" / "prod_seed0"
    run_dir.mkdir(parents=True)
    # Legacy-only ranking location (canonical final_result/final_ranking.csv is absent)
    (run_dir / "final_ranking.csv").write_text("rank,score\n1,-10\n", encoding="utf-8")

    # Legacy metadata location
    legacy_meta = run_dir / "restored" / "pyrosetta_run_metadata.json"
    legacy_meta.parent.mkdir(parents=True)
    legacy_meta.write_text("{}", encoding="utf-8")

    step2_dir = record_step2_outputs(
        config_path,
        repo_root=repo_root,
        ppi_targets=[{"name": "target A", "docking_dir": str(run_dir)}],
    )

    manifest = json.loads((step2_dir / "step_manifest.json").read_text(encoding="utf-8"))
    warning_codes = manifest.get("warning_codes", [])

    assert warning_codes
    assert any(item.get("code") == "LEGACY_PPI_PATH_USED" for item in warning_codes)

    # Ranking file should still be copied into step2 view
    copied_rankings = list(step2_dir.glob("*_final_ranking.csv"))
    assert copied_rankings, "Expected copied ranking artifact in step2 output"
