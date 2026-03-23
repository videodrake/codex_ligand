"""Group 8 tests: Pipeline logging and observability (F-4 Logging).

Covers:
1. organize_outputs() signature fix (F-4.4 / 8.1)
2. step_view permission fix (F-4.7 / 8.2)
3. Phase skip content validation (F-4.1 / 8.3)
4. Exit code on post-run failure (F-4.3 / 8.4)
5. Step health in completion summary (F-4.2 / 8.5)
6. Phase sanity check logging (F-4.5 / 8.6)
"""

import csv
import os
import tempfile
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent


# =====================================================================
# 8.1: organize_outputs() signature fix
# =====================================================================

class TestOrganizeOutputsSignature:
    """F-4.4: no repo_root kwarg in call site."""

    def test_call_site_has_no_repo_root(self):
        source = (PROJECT_ROOT / "run_production.py").read_text()
        # Should NOT pass repo_root to organize_outputs
        assert "organize_outputs(str(CONFIG_PATH))" in source
        assert "repo_root=REPO_ROOT" not in source.split("organize_outputs")[1].split(")")[0]

    def test_organize_outputs_importable(self):
        from egfr_pipeline.output_organizer import organize_outputs
        import inspect
        sig = inspect.signature(organize_outputs)
        assert "repo_root" not in sig.parameters


# =====================================================================
# 8.2: step_view permission fix
# =====================================================================

class TestStepViewPermissions:
    """F-4.7: _atomic_write_text sets 0644, _staged_step_dir sets 0755."""

    def test_atomic_write_text_sets_644(self):
        source = (PROJECT_ROOT / "egfr_pipeline" / "step_view.py").read_text()
        assert "0o644" in source or "0644" in source

    def test_staged_step_dir_sets_755(self):
        source = (PROJECT_ROOT / "egfr_pipeline" / "step_view.py").read_text()
        assert "0o755" in source or "0755" in source

    def test_atomic_write_creates_readable_file(self, tmp_path):
        from egfr_pipeline.step_view import _atomic_write_text

        test_file = tmp_path / "test_output.txt"
        _atomic_write_text(test_file, "hello world")
        assert test_file.exists()
        # Check permission: should be readable by group/other
        mode = test_file.stat().st_mode
        assert mode & 0o044 != 0  # group+other read bits set


# =====================================================================
# 8.3: Phase skip content validation
# =====================================================================

class TestPhaseSkipValidation:
    """F-4.1: check_phase6 detects placeholder reports."""

    def test_csv_has_rows_empty_file(self, tmp_path):
        from run_production import _csv_has_rows

        empty = tmp_path / "empty.csv"
        empty.write_text("")
        assert _csv_has_rows(empty) is False

    def test_csv_has_rows_header_only(self, tmp_path):
        from run_production import _csv_has_rows

        header = tmp_path / "header.csv"
        header.write_text("col1,col2\n")
        assert _csv_has_rows(header) is False

    def test_csv_has_rows_with_data(self, tmp_path):
        from run_production import _csv_has_rows

        data = tmp_path / "data.csv"
        data.write_text("col1,col2\na,b\n")
        assert _csv_has_rows(data) is True

    def test_check_phase6_placeholder_detection(self):
        """Verify check_phase6 source contains placeholder pattern detection."""
        source = (PROJECT_ROOT / "run_production.py").read_text()
        assert "placeholder_patterns" in source
        assert "No Vina pocket data available" in source


# =====================================================================
# 8.4: Exit code on failure
# =====================================================================

class TestExitCodeOnFailure:
    """F-4.3: sys.exit(1) on post-run warnings or phase failures."""

    def test_post_run_warnings_pattern_exists(self):
        source = (PROJECT_ROOT / "run_production.py").read_text()
        assert "post_run_warnings" in source
        assert "sys.exit(1)" in source

    def test_post_run_warnings_on_organize_failure(self):
        source = (PROJECT_ROOT / "run_production.py").read_text()
        assert "post_run_warnings.append" in source


# =====================================================================
# 8.5: Step health in summary
# =====================================================================

class TestStepHealthInSummary:
    """F-4.2: completion summary includes step view count."""

    def test_summary_has_step_view_section(self):
        source = (PROJECT_ROOT / "run_production.py").read_text()
        assert "Derived Step View" in source
        assert "generated" in source


# =====================================================================
# 8.6: Phase sanity check
# =====================================================================

class TestPhaseSanityCheck:
    """F-4.5: core output health check after phase completion."""

    def test_sanity_check_function_exists(self):
        from run_production import _sanity_check_phase_outputs
        assert callable(_sanity_check_phase_outputs)

    def test_phase_core_outputs_defined(self):
        from run_production import PHASE_CORE_OUTPUTS
        assert 4 in PHASE_CORE_OUTPUTS
        assert 5 in PHASE_CORE_OUTPUTS
        assert 6 in PHASE_CORE_OUTPUTS
        assert 7 in PHASE_CORE_OUTPUTS

    def test_sanity_check_called_after_completion(self):
        source = (PROJECT_ROOT / "run_production.py").read_text()
        assert "_sanity_check_phase_outputs(phase_num, config)" in source
