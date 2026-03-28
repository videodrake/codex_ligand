import json
from pathlib import Path

import pytest

from egfr_pipeline.ppi.postprocess_ppi import (
    _parse_run_dir_context,
    _write_legacy_migration_marker,
    _write_restored_manifest,
)

pytestmark = pytest.mark.unit


def test_parse_run_dir_context() -> None:
    run_type, seed = _parse_run_dir_context(Path("active/dock_seed7"))
    assert run_type == "dock"
    assert seed == "7"


def test_write_restored_manifest_has_expected_fields(tmp_path: Path) -> None:
    restored = tmp_path / "restored"
    restored.mkdir()
    source = tmp_path / "inactive" / "dock_seed11"
    source.parent.mkdir(parents=True)

    manifest_path = _write_restored_manifest(
        restored,
        source_docking_dir=source,
        receptor_id="3GT8_raw",
        partner_name="beta",
    )
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert payload["source_state"] == "inactive"
    assert payload["source_seed"] == "11"
    assert payload["receptor_id"] == "3GT8_raw"


def test_legacy_migration_marker_is_idempotent(tmp_path: Path) -> None:
    legacy = tmp_path / "legacy"
    legacy.mkdir()

    _write_legacy_migration_marker(legacy, tmp_path / "new_a")
    first = (legacy / "MOVED_TO.txt").read_text(encoding="utf-8")

    _write_legacy_migration_marker(legacy, tmp_path / "new_b")
    second = (legacy / "MOVED_TO.txt").read_text(encoding="utf-8")

    assert first == second
    assert "new_a" in first
