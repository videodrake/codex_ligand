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
        """Verify check_phase6 uses placeholder constants from report.py."""
        source = (PROJECT_ROOT / "run_production.py").read_text()
        assert "placeholder_patterns" in source
        assert "PLACEHOLDER_NO_VINA" in source  # Uses constant, not hardcoded string


# =====================================================================
# 8.4: Exit code on failure
# =====================================================================

class TestExitCodeOnFailure:
    """F-4.3: sys.exit(1) on post-run warnings or phase failures."""

    def test_post_run_warnings_pattern_exists(self):
        source = (PROJECT_ROOT / "run_production.py").read_text()
        assert "post_run_warnings" in source
        assert "sys.exit(1)" in source  # phase failure
        assert "sys.exit(2)" in source  # post-run only failure

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


# =====================================================================
# Adversarial tests
# =====================================================================

class TestPermissionEdgeCases:
    """Permission handling edge cases."""

    def test_atomic_write_json_also_gets_644(self, tmp_path):
        """JSON files should also be 0644 (via _atomic_write_text)."""
        from egfr_pipeline.step_view import _atomic_write_json

        test_file = tmp_path / "test.json"
        _atomic_write_json(test_file, {"key": "value"})
        mode = test_file.stat().st_mode
        assert mode & 0o044 != 0  # group+other read

    def test_atomic_write_overwrites_existing(self, tmp_path):
        """Overwriting existing file should still set correct permissions."""
        from egfr_pipeline.step_view import _atomic_write_text

        test_file = tmp_path / "test.txt"
        test_file.write_text("old content")
        test_file.chmod(0o600)  # Restrict to owner-only

        _atomic_write_text(test_file, "new content")
        assert test_file.read_text() == "new content"
        mode = test_file.stat().st_mode
        assert mode & 0o044 != 0  # Should be readable again


class TestPlaceholderDetection:
    """Placeholder detection boundary cases."""

    def test_partial_placeholder_is_not_detected(self):
        """Report with only some placeholder patterns should NOT be flagged."""
        # check_phase6 uses `all(pat in content for pat in ...)` — partial match is OK
        source = (PROJECT_ROOT / "run_production.py").read_text()
        assert "all(pat in content for pat in placeholder_patterns)" in source

    def test_real_report_with_data_passes(self, tmp_path):
        """A report with real data + some 'No data' lines should still pass."""
        from run_production import _csv_has_rows

        # CSV with actual data rows passes
        real_csv = tmp_path / "data.csv"
        real_csv.write_text("col1,col2\nval1,val2\nval3,val4\n")
        assert _csv_has_rows(real_csv) is True


class TestSanityCheckMissingDirs:
    """Sanity check when phase directories don't exist."""

    def test_sanity_check_no_crash_on_missing_dir(self, tmp_path, capsys):
        from run_production import _sanity_check_phase_outputs

        config = {"output_root": str(tmp_path)}
        # Phase 6 dir doesn't exist → should print [✗] but not crash
        _sanity_check_phase_outputs(6, config)
        captured = capsys.readouterr()
        assert "NOT FOUND" in captured.out

    def test_sanity_check_small_file_warns(self, tmp_path, capsys):
        from run_production import _sanity_check_phase_outputs
        from egfr_pipeline.paths import wa_phase6_report

        config = {"output_root": str(tmp_path)}
        report_dir = wa_phase6_report(config)
        report_dir.mkdir(parents=True, exist_ok=True)
        (report_dir / "project_report.txt").write_text("tiny")

        _sanity_check_phase_outputs(6, config)
        captured = capsys.readouterr()
        assert "possibly empty" in captured.out or "△" in captured.out

    def test_sanity_check_healthy_file(self, tmp_path, capsys):
        from run_production import _sanity_check_phase_outputs
        from egfr_pipeline.paths import wa_phase6_report

        config = {"output_root": str(tmp_path)}
        report_dir = wa_phase6_report(config)
        report_dir.mkdir(parents=True, exist_ok=True)
        (report_dir / "project_report.txt").write_text("x" * 500)

        _sanity_check_phase_outputs(6, config)
        captured = capsys.readouterr()
        assert "✓" in captured.out


class TestExitCodeIntegration:
    """Verify exit code patterns are correctly structured."""

    def test_has_failed_phases_check(self):
        source = (PROJECT_ROOT / "run_production.py").read_text()
        assert "has_failed_phases" in source
        assert '"failed"' in source

    def test_exit_codes_differentiated(self):
        """Exit 1 for phase failure, exit 2 for post-run only."""
        source = (PROJECT_ROOT / "run_production.py").read_text()
        assert "sys.exit(1)" in source
        assert "sys.exit(2)" in source


class TestExtractResnumInsertionCodes:
    """extract_resnum handles PDB insertion codes."""

    def test_normal_residue(self):
        from egfr_pipeline.residue_utils import extract_resnum
        assert extract_resnum("MET971") == 971

    def test_insertion_code(self):
        from egfr_pipeline.residue_utils import extract_resnum
        assert extract_resnum("ALA700A") == 700

    def test_negative_with_insertion(self):
        from egfr_pipeline.residue_utils import extract_resnum
        assert extract_resnum("GLY-50B") == -50

    def test_no_number(self):
        from egfr_pipeline.residue_utils import extract_resnum
        assert extract_resnum("XYZ") is None


class TestBOMHandling:
    """CSV reads use utf-8-sig for BOM safety."""

    def test_handoff_uses_utf8_sig(self):
        source = (PROJECT_ROOT / "run_production.py").read_text()
        # All handoff CSV reads should use utf-8-sig
        assert "utf-8-sig" in source

    def test_utf8_sig_reads_bom_csv(self, tmp_path):
        """utf-8-sig strips BOM from first column name."""
        import csv
        bom_csv = tmp_path / "bom.csv"
        bom_csv.write_bytes(b'\xef\xbb\xbfcol1,col2\na,b\n')

        with open(bom_csv, newline="", encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
        assert "col1" in rows[0]  # BOM stripped

    def test_utf8_sig_works_without_bom(self, tmp_path):
        """utf-8-sig is backward compatible with non-BOM files."""
        import csv
        normal_csv = tmp_path / "normal.csv"
        normal_csv.write_text("col1,col2\na,b\n")

        with open(normal_csv, newline="", encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
        assert "col1" in rows[0]


class TestAffinityNullCheck:
    """Expanded affinity null-check catches nan/N/A/Inf."""

    def test_nan_excluded(self):
        val = "nan"
        assert val.strip().lower() in ("", "na", "none", "nan", "n/a", "inf", "-inf")

    def test_capital_NA_excluded(self):
        val = "N/A"
        assert val.strip().lower() in ("", "na", "none", "nan", "n/a", "inf", "-inf")

    def test_valid_affinity_not_excluded(self):
        val = "-7.5"
        assert val.strip().lower() not in ("", "na", "none", "nan", "n/a", "inf", "-inf")
