import json
from pathlib import Path

import pytest

import run_production
from egfr_pipeline import paths
from egfr_pipeline.verdict import generate_verdict

pytestmark = pytest.mark.unit


def _write_min_config(path: Path) -> Path:
    cfg = {
        "project_name": "demo",
        "output_root": str(path / "output"),
        "receptors": [],
        "ligands": [],
    }
    cfg_path = path / "config.json"
    cfg_path.write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding="utf-8")
    return cfg_path


def test_verdict_raises_when_phase2_exists_but_phase3_missing(tmp_path: Path) -> None:
    config_path = _write_min_config(tmp_path)
    output_root = tmp_path / "output"

    phase2_run = output_root / "workflow_a" / "phase2_ppi_docking" / "active" / "prod_seed0"
    phase2_run.mkdir(parents=True)
    (phase2_run / "final_ranking.csv").write_text("rank,score\n1,-1\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="Phase 2 PPI 도킹 결과가 존재하지만 Phase 3 후처리 CSV가 없습니다"):
        generate_verdict(str(config_path))


def test_check_verdict_prerequisites_blocks_finalize_when_phase3_missing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config_path = _write_min_config(tmp_path)
    monkeypatch.setattr(run_production, "CONFIG_PATH", config_path)

    cfg = run_production._load_config()
    phase4 = paths.wa_phase4_vina_postprocess(cfg)
    phase4.mkdir(parents=True, exist_ok=True)
    (phase4 / "vina_pocket_table.csv").write_text("receptor_id,pocket_id\n", encoding="utf-8")

    phase2_run = paths.wa_phase2_ppi_docking(cfg) / "active" / "prod_seed0"
    phase2_run.mkdir(parents=True, exist_ok=True)
    (phase2_run / "final_ranking.csv").write_text("rank,score\n1,-1\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match=r"Phase 3 \(PPI 후처리\)"):
        run_production._check_verdict_prerequisites(cfg)
