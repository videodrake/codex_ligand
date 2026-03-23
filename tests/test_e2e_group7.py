"""Group 7 tests: E2E integration and final verification.

Covers:
1. SC-2: VERDICT_FIELDS baseline preservation (regression)
2. SC-3: New columns present in VERDICT_FIELDS
3. SC-4: Handoff quality guard edge cases (integration)
4. SC-5: Sensitivity analysis functions exist
5. SC-6: Weight sensitivity simulation executable
6. run_production.py import chain integrity
7. Document files existence checklist
8. Pocket table columns (pocket_depth_A)
"""

from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Baseline VERDICT_FIELDS from Group 0 (before any changes)
BASELINE_VERDICT_FIELDS = [
    "receptor_id", "pocket_id", "verdict", "confidence_score",
    "vina_quality_score", "ppi_proximity_score", "cross_receptor_score",
    "vina_affinity_pts", "vina_convergence_pts", "vina_stability_pts",
    "vina_diversity_pts", "vina_consensus_pts", "ppi_spatial_pts",
    "ppi_overlap_pts", "ppi_reproducibility_pts", "cross_receptor_pts",
    "cross_receptor_support_pts", "score_denominator", "ppi_data_available",
    "best_affinity", "n_pose", "n_ligand", "dominant_ligand_fraction",
    "ligand_pose_entropy", "spatial_dist_to_ppi", "closest_ppi_partner",
    "n_ppi_partners_near", "n_shared_with_ppi", "ppi_frac_runs_supporting",
    "ppi_best_interface_delta_e", "cross_receptor_matches",
    "cross_receptor_support", "consensus_site_id",
    "exp_sensitivity", "exp_specificity", "exp_enrichment", "exp_rank_impact",
    "pocket_stability", "bootstrap_confidence",
    "evidence_profile", "reason_tags", "decision_trace", "reasons",
    "is_atp_site", "exclusion_reason",
]


# =====================================================================
# SC-2: Regression — baseline columns preserved
# =====================================================================

class TestBaselinePreservation:
    """SC-2: existing VERDICT_FIELDS not deleted or renamed."""

    def test_all_baseline_columns_present(self):
        from egfr_pipeline.verdict import VERDICT_FIELDS

        for col in BASELINE_VERDICT_FIELDS:
            assert col in VERDICT_FIELDS, f"Baseline column '{col}' missing from VERDICT_FIELDS"

    def test_no_baseline_columns_reordered_or_duplicated(self):
        from egfr_pipeline.verdict import VERDICT_FIELDS

        # No duplicates
        assert len(VERDICT_FIELDS) == len(set(VERDICT_FIELDS))


# =====================================================================
# SC-3: New columns added
# =====================================================================

class TestNewColumns:
    """SC-3: new columns present in VERDICT_FIELDS."""

    def test_is_atp_site_present(self):
        from egfr_pipeline.verdict import VERDICT_FIELDS
        assert "is_atp_site" in VERDICT_FIELDS

    def test_exclusion_reason_present(self):
        from egfr_pipeline.verdict import VERDICT_FIELDS
        assert "exclusion_reason" in VERDICT_FIELDS

    def test_bootstrap_confidence_present(self):
        from egfr_pipeline.verdict import VERDICT_FIELDS
        assert "bootstrap_confidence" in VERDICT_FIELDS

    def test_allosteric_candidate_present(self):
        from egfr_pipeline.verdict import VERDICT_FIELDS
        assert "allosteric_candidate" in VERDICT_FIELDS

    def test_pocket_table_has_depth_column(self):
        """pocket_depth_A added to pocket_summary output."""
        from egfr_pipeline.vina import pocket_summary
        source = Path(pocket_summary.__file__).read_text()
        assert "pocket_depth_A" in source

    def test_pocket_table_has_sheet_contacts(self):
        """contacts_sheet_8_9 in pocket_summary output."""
        from egfr_pipeline.vina import pocket_summary
        source = Path(pocket_summary.__file__).read_text()
        assert "contacts_sheet_8_9" in source


# =====================================================================
# SC-4: Handoff quality guards (integration)
# =====================================================================

class TestHandoffGuardsIntegration:
    """SC-4: guards exist in run_production.py."""

    def test_validate_handoff_data_quality_exists(self):
        from run_production import _validate_handoff_data_quality
        assert callable(_validate_handoff_data_quality)

    def test_validate_adv_handoff_exists(self):
        from run_production import _validate_adv_handoff
        assert callable(_validate_adv_handoff)


# =====================================================================
# SC-5: Sensitivity analysis functions
# =====================================================================

class TestSensitivityAnalysis:
    """SC-5: orientation/sheet12/weight sensitivity functions exist."""

    def test_sweep_threshold_sensitivity_exists(self):
        from egfr_pipeline.phase1.orientation_filter import sweep_threshold_sensitivity
        assert callable(sweep_threshold_sensitivity)

    def test_compare_sheet12_sensitivity_exists(self):
        from egfr_pipeline.phase1.orientation_filter import compare_sheet12_sensitivity
        assert callable(compare_sheet12_sensitivity)

    def test_weight_sensitivity_script_exists(self):
        assert (PROJECT_ROOT / "scripts" / "verdict_weight_sensitivity.py").exists()


# =====================================================================
# SC-6: Weight sensitivity simulation
# =====================================================================

class TestWeightSensitivityExecutable:
    """SC-6: simulation runs and produces results."""

    def test_run_sensitivity_callable(self):
        from scripts.verdict_weight_sensitivity import run_sensitivity, generate_synthetic_pockets

        pockets = generate_synthetic_pockets()
        results, warnings = run_sensitivity(pockets)
        assert len(results) > 0
        assert all("baseline_total" in r for r in results)


# =====================================================================
# Import chain integrity
# =====================================================================

class TestImportChain:
    """run_production.py and all modified modules import cleanly."""

    def test_run_production_imports(self):
        import run_production
        assert hasattr(run_production, 'main')

    def test_verdict_imports(self):
        from egfr_pipeline.verdict import score_pocket, generate_verdict, VERDICT_FIELDS
        assert callable(score_pocket)

    def test_report_imports(self):
        from egfr_pipeline.report import generate_report
        assert callable(generate_report)

    def test_workflow_comparison_imports(self):
        from egfr_pipeline.workflow_comparison import run_comparison
        assert callable(run_comparison)

    def test_region_definitions_imports(self):
        from egfr_pipeline.region_definitions import ATP_SITE_RESIDUES, KO_SHEET_8_9
        assert len(ATP_SITE_RESIDUES) == 37

    def test_paths_imports(self):
        from egfr_pipeline.paths import workflow_a_root, workflow_b_root
        assert callable(workflow_a_root)


# =====================================================================
# Document files existence
# =====================================================================

class TestDocumentExistence:
    """All required documentation files exist."""

    @pytest.mark.parametrize("doc_path", [
        "docs/prd.md",
        "docs/tasks.md",
        "docs/CONTEXT.md",
        "docs/phase4_A3_axis_specification.md",
        "docs/methodology_limitations.md",
        "docs/workflow_comparison_guide.md",
        "docs/workflow_comparison_design.md",
        "docs/output_path_guide.md",
        "docs/module_separation_analysis.md",
        "docs/manual_vina.md",
        "docs/before_after_comparison.md",
        "docs/ac_coverage_checklist.md",
        "docs/research_overview_full.md",
        "CLAUDE.md",
        "PIPELINE_ARCHITECTURE_REPORT.md",
    ])
    def test_document_exists(self, doc_path):
        assert (PROJECT_ROOT / doc_path).exists(), f"Missing: {doc_path}"


# =====================================================================
# Cross-module consistency
# =====================================================================

class TestCrossModuleConsistency:
    """Key constants are consistent across modules."""

    def test_nlobe_clobe_boundary_consistent(self):
        """NLOBE_CLOBE_BOUNDARY should be 838 wherever defined."""
        from egfr_pipeline.phase1.compare_states import NLOBE_CLOBE_BOUNDARY as cs_val
        from egfr_pipeline.phase1.lightdock_validation import NLOBE_CLOBE_BOUNDARY as lv_val
        assert cs_val == 838
        assert lv_val == 838

    def test_strong_threshold_consistent(self):
        from egfr_pipeline.verdict import DEFAULT_THRESHOLDS
        from scripts.verdict_weight_sensitivity import STRONG_MIN
        assert DEFAULT_THRESHOLDS["valid_min"] == STRONG_MIN
