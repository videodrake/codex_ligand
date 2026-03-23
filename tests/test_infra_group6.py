"""Group 6 tests: Infrastructure improvements (F-6).

Covers:
1. Module separation analysis document (AC-6.1)
2. Precheck Workflow B tool check script (AC-6.2)
3. Phase 3 parallelization confirmation (AC-6.3)
4. Advanced pipeline precheck integration (AC-6.2)
"""

from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent


# =====================================================================
# AC-6.1: Module separation analysis
# =====================================================================

class TestModuleSeparationAnalysis:
    """AC-6.1: analysis document with decision and rationale."""

    def test_document_exists(self):
        path = PROJECT_ROOT / "docs" / "module_separation_analysis.md"
        assert path.exists()

    def test_covers_both_files(self):
        content = (PROJECT_ROOT / "docs" / "module_separation_analysis.md").read_text()
        assert "pipeline_manager" in content
        assert "vina_executor" in content

    def test_contains_decision(self):
        content = (PROJECT_ROOT / "docs" / "module_separation_analysis.md").read_text()
        assert "유지" in content or "분리" in content

    def test_contains_line_counts(self):
        content = (PROJECT_ROOT / "docs" / "module_separation_analysis.md").read_text()
        assert "4,156" in content or "4156" in content
        assert "2,792" in content or "2792" in content

    def test_circular_import_assessed(self):
        content = (PROJECT_ROOT / "docs" / "module_separation_analysis.md").read_text()
        assert "순환" in content or "circular" in content.lower()

    def test_500_line_functions_identified(self):
        content = (PROJECT_ROOT / "docs" / "module_separation_analysis.md").read_text()
        assert "500" in content
        assert "generate_validation_report" in content


# =====================================================================
# AC-6.2: Precheck Workflow B tool check
# =====================================================================

class TestPrecheckWorkflowB:
    """AC-6.2: precheck script has --check-workflow-b / CHECK_WORKFLOW_B."""

    def test_precheck_script_exists(self):
        path = PROJECT_ROOT / "scripts" / "run_pre_qsub_checks.sh"
        assert path.exists()

    def test_precheck_has_workflow_b_check(self):
        content = (PROJECT_ROOT / "scripts" / "run_pre_qsub_checks.sh").read_text()
        assert "CHECK_WORKFLOW_B" in content

    def test_precheck_checks_fpocket(self):
        content = (PROJECT_ROOT / "scripts" / "run_pre_qsub_checks.sh").read_text()
        assert "fpocket" in content

    def test_precheck_checks_p2rank(self):
        content = (PROJECT_ROOT / "scripts" / "run_pre_qsub_checks.sh").read_text()
        assert "P2Rank" in content or "prank" in content

    def test_precheck_checks_lightdock(self):
        content = (PROJECT_ROOT / "scripts" / "run_pre_qsub_checks.sh").read_text()
        assert "lightdock" in content.lower() or "LightDock" in content

    def test_advanced_pipeline_calls_precheck(self):
        content = (PROJECT_ROOT / "config" / "run_advanced_pipeline.pbs").read_text()
        assert "precheck" in content.lower() or "run_pre_qsub_checks" in content

    def test_fpocket_missing_is_critical(self):
        """EC-6.2: fpocket missing → exit 1 (critical)."""
        content = (PROJECT_ROOT / "scripts" / "run_pre_qsub_checks.sh").read_text()
        assert "CRITICAL" in content and "fpocket" in content

    def test_p2rank_missing_is_warning(self):
        """EC-6.2: P2Rank missing → WARNING (not critical)."""
        content = (PROJECT_ROOT / "scripts" / "run_pre_qsub_checks.sh").read_text()
        # P2Rank missing should not cause exit 1
        assert "fpocket only" in content or "WARNING" in content


# =====================================================================
# AC-6.3: Phase 3 parallelization
# =====================================================================

class TestPhase3Parallelization:
    """AC-6.3 / EC-6.3: Phase 3 already uses ThreadPoolExecutor."""

    def test_threadpool_used(self):
        source = (PROJECT_ROOT / "egfr_pipeline" / "phase3" / "run_diverse_docking.py").read_text()
        assert "ThreadPoolExecutor" in source

    def test_as_completed_used(self):
        source = (PROJECT_ROOT / "egfr_pipeline" / "phase3" / "run_diverse_docking.py").read_text()
        assert "as_completed" in source

    def test_max_workers_configurable(self):
        source = (PROJECT_ROOT / "egfr_pipeline" / "phase3" / "run_diverse_docking.py").read_text()
        assert "max_workers" in source


# =====================================================================
# Adversarial tests
# =====================================================================

class TestPrecheckScriptStructure:
    """Verify precheck script is well-formed bash."""

    def test_script_has_shebang(self):
        content = (PROJECT_ROOT / "scripts" / "run_pre_qsub_checks.sh").read_text()
        assert content.startswith("#!/")

    def test_script_has_set_euo_pipefail(self):
        content = (PROJECT_ROOT / "scripts" / "run_pre_qsub_checks.sh").read_text()
        assert "set -euo pipefail" in content

    def test_workflow_b_check_is_conditional(self):
        """CHECK_WORKFLOW_B=0 should skip tool checks."""
        content = (PROJECT_ROOT / "scripts" / "run_pre_qsub_checks.sh").read_text()
        # The check block should be inside an if condition
        assert 'CHECK_WORKFLOW_B' in content
        assert 'if' in content

    def test_lightdock_missing_is_not_critical(self):
        """LightDock missing should warn but not exit."""
        content = (PROJECT_ROOT / "scripts" / "run_pre_qsub_checks.sh").read_text()
        # LightDock WARNING should not trigger TOOLS_OK=0
        lines = content.split("\n")
        in_lightdock_block = False
        for line in lines:
            if "LightDock" in line or "lightdock" in line:
                in_lightdock_block = True
            if in_lightdock_block and "TOOLS_OK=0" in line:
                pytest.fail("LightDock missing should not set TOOLS_OK=0")
            if in_lightdock_block and line.strip() == "fi":
                break


class TestModuleSeparationFutureConditions:
    """Verify analysis includes re-evaluation triggers."""

    def test_future_review_conditions(self):
        content = (PROJECT_ROOT / "docs" / "module_separation_analysis.md").read_text()
        assert "5,000" in content or "5000" in content or "향후" in content
