#!/usr/bin/env python3
"""Tests for Phase 4: Perturbation Relevance Scoring (TG 4.0–4.7)."""

import csv
import json
from pathlib import Path
from typing import Dict, List

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write(path: Path, lines: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def _load_csv(path: Path) -> List[dict]:
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _csv_columns(path: Path) -> List[str]:
    with open(path, encoding="utf-8") as f:
        return next(csv.reader(f))


# ---------------------------------------------------------------------------
# CSV header constants (matching production schemas)
# ---------------------------------------------------------------------------

PHASE1_PATCH_HEADER = (
    "residue_id,residue_num,residue_name,chain,lobe_label,evidence_source,"
    "construct_type,orientation_validation_status,robustness_class,"
    "n_states_present,states_present,global_max_occupancy,is_hotspot_any_state,"
    "pyrosetta_evidence,lightdock_evidence,method_agreement,confidence,notes"
)

PHASE2_RELATIONSHIP_HEADER = (
    "receptor_id,pocket_id,relationship_class,hotspot_overlap_count,"
    "hotspot_overlap_fraction,centroid_distance_A,overlapping_residues,"
    "n_pocket_residues,n_patch_hotspots,classification_basis"
)

PHASE2_DRUGGABILITY_HEADER = (
    "receptor_id,pocket_id,relationship_class,best_proposal_score,"
    "best_druggability_score,normalized_proposal_score,"
    "normalized_druggability_score,druggability_confidence,n_sources,sources,"
    "multi_source_support,mean_volume_A3,n_residues,ftmap_hotspot_count,"
    "ftmap_overlap_residues,hotspot_support_available,overall_druggability_tier"
)

PHASE2_STATE_CLASS_HEADER = (
    "pocket_id,receptor_id,state_class,n_states_matched,matched_states,"
    "matched_pocket_ids,best_match_distance_A,best_match_jaccard,"
    "state_class_basis,druggability_confidence,overall_druggability_tier,"
    "relationship_class"
)

PHASE3_EVIDENCE_HEADER = (
    "receptor_id,candidate_pocket_id,ligand_id,pose_support_count,"
    "best_affinity_kcal,mean_affinity_kcal,n_seeds_completed,docking_mode,"
    "contact_residues,n_contact_residues,relationship_class,"
    "overall_druggability_tier,druggability_confidence,docking_priority,"
    "state_class,n_states_matched,pocket_status,diversity_verdict,"
    "coverage_fraction,normalized_entropy,ligand_support_strength,"
    "total_poses_in_pocket,fraction_of_receptor_poses"
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_phase4(tmp_path):
    """Full multi-pocket, multi-ligand test environment.

    Creates:
    - 3 hotspot residues (2 robust, 1 state-specific)
    - 3 pockets: PKT01 (orthosteric, tier_1), PKT02 (orthosteric, tier_2),
      PKT03 (low_relevance, tier_3)
    - 2 ligands: ligA, ligB
    - Phase 3 evidence for PKT01 and PKT02 (2 pockets x 2 ligands = 4 rows)
    """
    phase1 = tmp_path / "phase1"
    phase2 = tmp_path / "phase2"
    phase3 = tmp_path / "phase3"
    phase4 = tmp_path / "phase4"

    # Phase 1: patch reference
    _write(phase1 / "phase1_downstream_patch_reference.csv", [
        PHASE1_PATCH_HEADER,
        "ASP855,855,ASP,A,C-lobe,primary,full_kinase_domain,orientation_validated,moderate,2,3GT8_raw;3GT8_cl38,1.0,True,True,True,both,medium,",
        "GLU866,866,GLU,A,C-lobe,primary,full_kinase_domain,orientation_validated,moderate,2,3GT8_raw;3GT8_cl38,1.0,True,True,True,both,medium,",
        "LEU819,819,LEU,A,N-lobe,primary,full_kinase_domain,orientation_validated,state_specific,1,3GT8_raw,1.0,True,True,False,single,low,",
    ])

    # Phase 2: relationship
    _write(phase2 / "pocket_patch_relationship.csv", [
        PHASE2_RELATIONSHIP_HEADER,
        "3GT8_raw,3GT8_raw_PKT01,orthosteric_candidate,2,0.667,130.0,ASP855;GLU866,3,3,hotspot_overlap",
        "3GT8_raw,3GT8_raw_PKT02,orthosteric_candidate,1,0.333,125.0,LEU819,3,3,hotspot_overlap",
        "3GT8_raw,3GT8_raw_PKT03,low_relevance_candidate,0,0.000,150.0,,2,3,no_overlap",
    ])

    # Phase 2: druggability
    _write(phase2 / "druggability_proposal_summary.csv", [
        PHASE2_DRUGGABILITY_HEADER,
        "3GT8_raw,3GT8_raw_PKT01,orthosteric_candidate,0.82,0.71,1.0,1.0,high,1,fpocket,False,890.0,3,,,False,tier_1",
        "3GT8_raw,3GT8_raw_PKT02,orthosteric_candidate,0.56,0.44,0.48,0.48,medium,1,fpocket,False,560.0,3,,,False,tier_2",
        "3GT8_raw,3GT8_raw_PKT03,low_relevance_candidate,0.32,0.20,0.0,0.0,low,1,fpocket,False,310.0,2,,,False,tier_3",
    ])

    # Phase 2: state classes
    _write(phase2 / "candidate_pocket_state_classes.csv", [
        PHASE2_STATE_CLASS_HEADER,
        "3GT8_raw_PKT01,3GT8_raw,state_specific_pocket,1,,,,,single state,high,tier_1,orthosteric_candidate",
        "3GT8_raw_PKT02,3GT8_raw,state_specific_pocket,1,,,,,single state,medium,tier_2,orthosteric_candidate",
        "3GT8_raw_PKT03,3GT8_raw,state_specific_pocket,1,,,,,single state,low,tier_3,low_relevance_candidate",
    ])

    # Phase 3: docking evidence (2 pockets x 2 ligands = 4 rows)
    _write(phase3 / "phase4_docking_evidence_reference.csv", [
        PHASE3_EVIDENCE_HEADER,
        "3GT8_raw,3GT8_raw_PKT01,ligA,5,-6.5,-6.0,3,focused,ASP855;GLU866,2,orthosteric_candidate,tier_1,high,primary,state_specific_pocket,1,open,well_distributed,1.000,1.000,moderate,5,0.500",
        "3GT8_raw,3GT8_raw_PKT02,ligA,2,-5.0,-4.8,3,focused,LEU819,1,orthosteric_candidate,tier_2,medium,primary,state_specific_pocket,1,open,well_distributed,1.000,1.000,weak,2,0.500",
        "3GT8_raw,3GT8_raw_PKT01,ligB,3,-5.8,-5.5,3,focused,ASP855,1,orthosteric_candidate,tier_1,high,primary,state_specific_pocket,1,open,well_distributed,1.000,1.000,moderate,3,0.500",
        "3GT8_raw,3GT8_raw_PKT02,ligB,1,,,3,focused,,0,orthosteric_candidate,tier_2,medium,primary,state_specific_pocket,1,open,well_distributed,1.000,1.000,weak,1,0.500",
    ])

    return {
        "phase1": phase1,
        "phase2": phase2,
        "phase3": phase3,
        "phase4": phase4,
    }


@pytest.fixture
def tmp_phase4_no_phase1(tmp_path):
    """Environment with missing Phase 1 data (graceful degradation)."""
    phase1 = tmp_path / "phase1"  # intentionally empty
    phase2 = tmp_path / "phase2"
    phase3 = tmp_path / "phase3"
    phase4 = tmp_path / "phase4"

    phase1.mkdir(parents=True, exist_ok=True)

    _write(phase2 / "pocket_patch_relationship.csv", [
        PHASE2_RELATIONSHIP_HEADER,
        "3GT8_raw,3GT8_raw_PKT01,orthosteric_candidate,0,0.0,130.0,,3,0,no_patch_data",
    ])
    _write(phase2 / "druggability_proposal_summary.csv", [
        PHASE2_DRUGGABILITY_HEADER,
        "3GT8_raw,3GT8_raw_PKT01,orthosteric_candidate,0.80,0.70,1.0,1.0,high,1,fpocket,False,800.0,3,,,False,tier_1",
    ])
    _write(phase2 / "candidate_pocket_state_classes.csv", [
        PHASE2_STATE_CLASS_HEADER,
        "3GT8_raw_PKT01,3GT8_raw,state_specific_pocket,1,,,,,single state,high,tier_1,orthosteric_candidate",
    ])
    _write(phase3 / "phase4_docking_evidence_reference.csv", [
        PHASE3_EVIDENCE_HEADER,
        "3GT8_raw,3GT8_raw_PKT01,ligA,3,-6.0,-5.5,3,focused,,0,orthosteric_candidate,tier_1,high,primary,state_specific_pocket,1,open,well_distributed,1.000,1.000,moderate,3,1.0",
    ])

    return {
        "phase1": phase1,
        "phase2": phase2,
        "phase3": phase3,
        "phase4": phase4,
    }


@pytest.fixture
def tmp_phase4_multi_state(tmp_path):
    """Environment with 2 receptor states (robust pocket)."""
    phase1 = tmp_path / "phase1"
    phase2 = tmp_path / "phase2"
    phase3 = tmp_path / "phase3"
    phase4 = tmp_path / "phase4"

    _write(phase1 / "phase1_downstream_patch_reference.csv", [
        PHASE1_PATCH_HEADER,
        "ASP855,855,ASP,A,C-lobe,primary,full_kinase_domain,orientation_validated,robust,3,all,1.0,True,True,True,both,high,",
        "GLU866,866,GLU,A,C-lobe,primary,full_kinase_domain,orientation_validated,robust,3,all,1.0,True,True,True,both,high,",
    ])
    _write(phase2 / "pocket_patch_relationship.csv", [
        PHASE2_RELATIONSHIP_HEADER,
        "3GT8_raw,3GT8_raw_PKT01,orthosteric_candidate,2,1.0,10.0,ASP855;GLU866,2,2,hotspot_overlap",
    ])
    _write(phase2 / "druggability_proposal_summary.csv", [
        PHASE2_DRUGGABILITY_HEADER,
        "3GT8_raw,3GT8_raw_PKT01,orthosteric_candidate,0.90,0.85,1.0,1.0,high,2,fpocket;p2rank,True,1000.0,5,,,False,tier_1",
    ])
    _write(phase2 / "candidate_pocket_state_classes.csv", [
        PHASE2_STATE_CLASS_HEADER,
        "3GT8_raw_PKT01,3GT8_raw,robust_pocket,3,3GT8_raw;cl38;cl85,PKT01;PKT01;PKT01,2.0,0.8,matched,high,tier_1,orthosteric_candidate",
    ])
    _write(phase3 / "phase4_docking_evidence_reference.csv", [
        PHASE3_EVIDENCE_HEADER,
        "3GT8_raw,3GT8_raw_PKT01,ligA,8,-7.5,-7.0,3,focused,ASP855,1,orthosteric_candidate,tier_1,high,primary,robust_pocket,3,open,well_distributed,1.000,1.000,strong,8,1.0",
    ])

    return {
        "phase1": phase1,
        "phase2": phase2,
        "phase3": phase3,
        "phase4": phase4,
    }


# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------

def _run_phase4_pipeline(env) -> Path:
    """Run full TG 4.0 → 4.7 on test environment. Returns output dir."""
    from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
    from egfr_pipeline.phase4.score_framework import run_score_framework
    from egfr_pipeline.phase4.mechanistic_classification import run_mechanistic_classification
    from egfr_pipeline.phase4.perturbation_scoring import run_perturbation_scoring
    from egfr_pipeline.phase4.state_interpretation import run_state_interpretation
    from egfr_pipeline.phase4.review_output import run_review_output
    from egfr_pipeline.phase4.final_report import run_final_report
    from egfr_pipeline.phase4.presentation_summary import run_presentation_summary

    out = env["phase4"]
    run_evidence_ingestion(env["phase1"], env["phase2"], env["phase3"], out)
    run_score_framework(out)
    run_mechanistic_classification(out)
    run_perturbation_scoring(out)
    run_state_interpretation(out)
    run_review_output(out)
    run_final_report(out)
    run_presentation_summary(out)
    return out


# ===================================================================
# TG 4.0: Evidence Ingestion
# ===================================================================

class TestEvidenceIngestion:
    """Tests for TG 4.0: Multi-Phase Evidence Ingestion."""

    def test_outputs_created(self, tmp_phase4):
        from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
        e = tmp_phase4
        norm, val = run_evidence_ingestion(e["phase1"], e["phase2"], e["phase3"], e["phase4"])
        assert norm.exists()
        assert val.exists()

    def test_normalized_row_count(self, tmp_phase4):
        from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
        e = tmp_phase4
        run_evidence_ingestion(e["phase1"], e["phase2"], e["phase3"], e["phase4"])
        rows = _load_csv(e["phase4"] / "phase4_evidence_normalized.csv")
        # 4 Phase 3 evidence rows + 1 skipped pocket (PKT03) = 5
        assert len(rows) == 5

    def test_skipped_pocket_included(self, tmp_phase4):
        from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
        e = tmp_phase4
        run_evidence_ingestion(e["phase1"], e["phase2"], e["phase3"], e["phase4"])
        rows = _load_csv(e["phase4"] / "phase4_evidence_normalized.csv")
        pkt03 = [r for r in rows if r["pocket_id"] == "3GT8_raw_PKT03"]
        assert len(pkt03) == 1
        assert pkt03[0]["phase3_loaded"] == "False"
        assert pkt03[0]["ligand_support_strength"] == "none"

    def test_phase1_summary_populated(self, tmp_phase4):
        from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
        e = tmp_phase4
        run_evidence_ingestion(e["phase1"], e["phase2"], e["phase3"], e["phase4"])
        rows = _load_csv(e["phase4"] / "phase4_evidence_normalized.csv")
        row = rows[0]
        assert row["phase1_loaded"] == "True"
        assert int(row["n_hotspot_residues"]) == 3
        assert int(row["n_method_agreement_both"]) == 2

    def test_validation_report_pass(self, tmp_phase4):
        from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
        e = tmp_phase4
        run_evidence_ingestion(e["phase1"], e["phase2"], e["phase3"], e["phase4"])
        text = (e["phase4"] / "phase4_evidence_validation.md").read_text(encoding="utf-8")
        assert "PASS" in text or "CAUTION" in text
        assert "0 errors" in text

    def test_no_phase1_graceful(self, tmp_phase4_no_phase1):
        from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
        e = tmp_phase4_no_phase1
        norm, val = run_evidence_ingestion(e["phase1"], e["phase2"], e["phase3"], e["phase4"])
        assert norm.exists()
        rows = _load_csv(norm)
        assert rows[0]["phase1_loaded"] == "False"


# ===================================================================
# TG 4.1: Score Framework
# ===================================================================

class TestScoreFrameworkUnit:
    """Unit tests for TG 4.1 axis score computation."""

    def test_a1_perfect_overlap(self):
        from egfr_pipeline.phase4.score_framework import compute_a1_ppi_interface
        row = {
            "hotspot_overlap_fraction": "1.0",
            "mean_hotspot_confidence": "3.0",
            "n_hotspot_residues": "5",
            "n_robust_hotspots": "5",
            "n_method_agreement_both": "5",
        }
        score = compute_a1_ppi_interface(row)
        assert score > 0.9

    def test_a1_zero_overlap(self):
        from egfr_pipeline.phase4.score_framework import compute_a1_ppi_interface
        row = {
            "hotspot_overlap_fraction": "0",
            "mean_hotspot_confidence": "0",
            "n_hotspot_residues": "0",
            "n_robust_hotspots": "0",
            "n_method_agreement_both": "0",
        }
        score = compute_a1_ppi_interface(row)
        assert score == 0.0

    def test_a2_tier1_high(self):
        from egfr_pipeline.phase4.score_framework import compute_a2_druggability
        row = {
            "overall_druggability_tier": "tier_1",
            "druggability_confidence": "high",
            "best_proposal_score": "0.9",
        }
        score = compute_a2_druggability(row)
        assert score > 0.8

    def test_a2_tier3_low(self):
        from egfr_pipeline.phase4.score_framework import compute_a2_druggability
        row = {
            "overall_druggability_tier": "tier_3",
            "druggability_confidence": "low",
            "best_proposal_score": "0.1",
        }
        score = compute_a2_druggability(row)
        assert score < 0.3

    def test_a3_orthosteric_with_support(self):
        from egfr_pipeline.phase4.score_framework import compute_a3_perturbation_relevance
        row = {
            "relationship_class": "orthosteric_candidate",
            "hotspot_overlap_count": "4",
            "n_hotspot_residues": "7",
            "ligand_support_strength": "strong",
            "pose_support_count": "10",
        }
        score = compute_a3_perturbation_relevance(row)
        assert score > 0.7

    def test_a3_irrelevant_no_support(self):
        from egfr_pipeline.phase4.score_framework import compute_a3_perturbation_relevance
        row = {
            "relationship_class": "low_relevance_candidate",
            "hotspot_overlap_count": "0",
            "n_hotspot_residues": "7",
            "ligand_support_strength": "none",
            "pose_support_count": "0",
        }
        score = compute_a3_perturbation_relevance(row)
        assert score < 0.1

    def test_a4_robust_3_states(self):
        from egfr_pipeline.phase4.score_framework import compute_a4_state_robustness
        row = {
            "state_class": "robust_pocket",
            "n_states_matched": "3",
            "diversity_verdict": "well_distributed",
            "coverage_fraction": "1.0",
        }
        score = compute_a4_state_robustness(row)
        assert score > 0.8

    def test_a4_state_specific(self):
        from egfr_pipeline.phase4.score_framework import compute_a4_state_robustness
        row = {
            "state_class": "state_specific_pocket",
            "n_states_matched": "1",
            "diversity_verdict": "",
            "coverage_fraction": "0",
        }
        score = compute_a4_state_robustness(row)
        assert 0.2 < score < 0.5


class TestScoreFrameworkIntegration:
    """Integration tests for TG 4.1."""

    def test_axis_definition_table(self, tmp_phase4):
        from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
        from egfr_pipeline.phase4.score_framework import run_score_framework
        e = tmp_phase4
        run_evidence_ingestion(e["phase1"], e["phase2"], e["phase3"], e["phase4"])
        run_score_framework(e["phase4"])
        axes = _load_csv(e["phase4"] / "phase4_axis_definition_table.csv")
        assert len(axes) == 4
        weights = [float(a["weight"]) for a in axes]
        assert abs(sum(weights) - 1.0) < 0.001

    def test_axis_scores_computed(self, tmp_phase4):
        from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
        from egfr_pipeline.phase4.score_framework import run_score_framework
        e = tmp_phase4
        run_evidence_ingestion(e["phase1"], e["phase2"], e["phase3"], e["phase4"])
        run_score_framework(e["phase4"])
        scored = _load_csv(e["phase4"] / "phase4_axis_scores.csv")
        assert len(scored) == 5
        for s in scored:
            for axis in ["A1_ppi_interface", "A2_druggability",
                         "A3_perturbation_relevance", "A4_state_robustness"]:
                val = float(s[axis])
                assert 0.0 <= val <= 1.0, f"{axis} out of range: {val}"


# ===================================================================
# TG 4.2: Mechanistic Classification
# ===================================================================

class TestMechanisticClassificationUnit:
    """Unit tests for TG 4.2 classification logic."""

    def test_orthosteric_disruptor(self):
        from egfr_pipeline.phase4.mechanistic_classification import classify_candidate
        row = {
            "relationship_class": "orthosteric_candidate",
            "hotspot_overlap_count": "3",
            "hotspot_overlap_fraction": "0.5",
            "overall_druggability_tier": "tier_1",
            "ligand_support_strength": "moderate",
            "pose_support_count": "5",
        }
        cls, basis, conf = classify_candidate(row)
        assert cls == "orthosteric_disruptor_candidate"
        assert conf == "high"

    def test_orthosteric_weak_overlap_is_uncertain(self):
        from egfr_pipeline.phase4.mechanistic_classification import classify_candidate
        row = {
            "relationship_class": "orthosteric_candidate",
            "hotspot_overlap_count": "1",
            "hotspot_overlap_fraction": "0.14",
            "overall_druggability_tier": "tier_2",
            "ligand_support_strength": "weak",
            "pose_support_count": "1",
        }
        cls, basis, conf = classify_candidate(row)
        assert cls == "uncertain_mechanism_candidate"

    def test_low_relevance_classified_irrelevant(self):
        from egfr_pipeline.phase4.mechanistic_classification import classify_candidate
        row = {
            "relationship_class": "low_relevance_candidate",
            "hotspot_overlap_count": "0",
            "hotspot_overlap_fraction": "0.0",
            "overall_druggability_tier": "tier_3",
            "ligand_support_strength": "none",
            "pose_support_count": "0",
        }
        cls, basis, conf = classify_candidate(row)
        assert cls == "ligandable_but_ppi_irrelevant_candidate"
        assert conf == "high"

    def test_allosteric_tier1(self):
        from egfr_pipeline.phase4.mechanistic_classification import classify_candidate
        row = {
            "relationship_class": "allosteric_candidate",
            "hotspot_overlap_count": "0",
            "hotspot_overlap_fraction": "0.0",
            "overall_druggability_tier": "tier_1",
            "ligand_support_strength": "moderate",
            "pose_support_count": "3",
        }
        cls, basis, conf = classify_candidate(row)
        assert cls == "allosteric_modulator_candidate"
        assert conf == "medium"

    def test_rim_candidate(self):
        from egfr_pipeline.phase4.mechanistic_classification import classify_candidate
        row = {
            "relationship_class": "rim_candidate",
            "hotspot_overlap_count": "2",
            "hotspot_overlap_fraction": "0.3",
            "overall_druggability_tier": "tier_2",
            "ligand_support_strength": "weak",
            "pose_support_count": "1",
        }
        cls, basis, conf = classify_candidate(row)
        assert cls == "interface_rim_modulator_candidate"
        assert conf == "medium"


class TestMechanisticClassificationIntegration:
    """Integration tests for TG 4.2."""

    def test_classification_output(self, tmp_phase4):
        from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
        from egfr_pipeline.phase4.score_framework import run_score_framework
        from egfr_pipeline.phase4.mechanistic_classification import run_mechanistic_classification
        e = tmp_phase4
        run_evidence_ingestion(e["phase1"], e["phase2"], e["phase3"], e["phase4"])
        run_score_framework(e["phase4"])
        run_mechanistic_classification(e["phase4"])
        classes = _load_csv(e["phase4"] / "final_candidate_classes.csv")
        assert len(classes) == 5
        mech_set = set(c["mechanistic_class"] for c in classes)
        assert "orthosteric_disruptor_candidate" in mech_set


# ===================================================================
# TG 4.3: Perturbation Scoring
# ===================================================================

class TestPerturbationScoringUnit:
    """Unit tests for TG 4.3 scoring logic."""

    def test_weighted_score(self):
        from egfr_pipeline.phase4.perturbation_scoring import compute_perturbation_score
        row = {
            "A1_ppi_interface": "0.5",
            "A2_druggability": "0.8",
            "A3_perturbation_relevance": "0.6",
            "A4_state_robustness": "0.4",
        }
        score, contribs = compute_perturbation_score(row)
        expected = 0.5 * 0.30 + 0.8 * 0.25 + 0.6 * 0.30 + 0.4 * 0.15
        assert abs(score - round(expected, 4)) < 0.001

    def test_affinity_cap_irrelevant_site(self):
        from egfr_pipeline.phase4.perturbation_scoring import apply_affinity_cap
        row = {
            "mechanistic_class": "ligandable_but_ppi_irrelevant_candidate",
            "A1_ppi_interface": "0.1",
            "A2_druggability": "0.9",
            "A3_perturbation_relevance": "0.8",
            "A4_state_robustness": "0.1",
        }
        raw_score = 0.1 * 0.30 + 0.9 * 0.25 + 0.8 * 0.30 + 0.1 * 0.15
        capped = apply_affinity_cap(row, raw_score)
        assert capped < raw_score, "Cap should reduce score for irrelevant site"

    def test_affinity_cap_not_applied_orthosteric(self):
        from egfr_pipeline.phase4.perturbation_scoring import apply_affinity_cap
        row = {
            "mechanistic_class": "orthosteric_disruptor_candidate",
            "A1_ppi_interface": "0.5",
            "A2_druggability": "0.9",
            "A3_perturbation_relevance": "0.7",
            "A4_state_robustness": "0.5",
        }
        raw_score = 0.5 * 0.30 + 0.9 * 0.25 + 0.7 * 0.30 + 0.5 * 0.15
        capped = apply_affinity_cap(row, raw_score)
        assert capped == raw_score, "No cap for orthosteric sites"

    def test_ranking_order(self, tmp_phase4):
        """PKT01 should rank above PKT02, PKT03 last."""
        out = _run_phase4_pipeline(tmp_phase4)
        table = _load_csv(out / "perturbation_candidate_table.csv")
        rank1 = table[0]
        assert rank1["pocket_id"] == "3GT8_raw_PKT01"
        last = table[-1]
        assert last["pocket_id"] == "3GT8_raw_PKT03"


# ===================================================================
# TG 4.4: State Interpretation
# ===================================================================

class TestStateInterpretation:
    """Tests for TG 4.4 state robustness."""

    def test_single_state_flagged(self, tmp_phase4):
        out = _run_phase4_pipeline(tmp_phase4)
        interp = _load_csv(out / "phase4_state_interpretation.csv")
        for r in interp:
            assert r["accessibility_label"] == "single_state_only"
            assert "provisional" in r["adjusted_confidence"]

    def test_robust_pocket_upgraded(self, tmp_phase4_multi_state):
        out = _run_phase4_pipeline(tmp_phase4_multi_state)
        interp = _load_csv(out / "phase4_state_interpretation.csv")
        assert len(interp) == 1
        r = interp[0]
        assert r["state_interpretation"] == "persistent"
        assert r["accessibility_label"] == "always_accessible"
        # Robust should strengthen confidence
        assert r["adjusted_confidence"] == "high"

    def test_accessibility_note_created(self, tmp_phase4):
        out = _run_phase4_pipeline(tmp_phase4)
        note = (out / "phase4_accessibility_note.md").read_text(encoding="utf-8")
        assert "State-Robustness" in note


# ===================================================================
# TG 4.5: Review Output
# ===================================================================

@pytest.mark.reporting
class TestReviewOutput:
    """Tests for TG 4.5 review tables."""

    def test_review_table_columns(self, tmp_phase4):
        from egfr_pipeline.phase4.review_output import REVIEW_TABLE_COLUMNS
        out = _run_phase4_pipeline(tmp_phase4)
        actual = _csv_columns(out / "phase4_final_review_table.csv")
        assert actual == REVIEW_TABLE_COLUMNS

    def test_expanded_table_columns(self, tmp_phase4):
        from egfr_pipeline.phase4.review_output import EXPANDED_TABLE_COLUMNS
        out = _run_phase4_pipeline(tmp_phase4)
        actual = _csv_columns(out / "phase4_expanded_evidence_table.csv")
        assert actual == EXPANDED_TABLE_COLUMNS

    def test_expanded_has_all_phases(self, tmp_phase4):
        out = _run_phase4_pipeline(tmp_phase4)
        expanded = _load_csv(out / "phase4_expanded_evidence_table.csv")
        for r in expanded:
            if r["ligand_id"]:
                assert r["phase1_loaded"] == "True"
                assert r["phase2_loaded"] == "True"
                assert r["phase3_loaded"] == "True"


# ===================================================================
# TG 4.6: Final Report
# ===================================================================

@pytest.mark.reporting
class TestFinalReport:
    """Tests for TG 4.6 integrated report."""

    def test_report_created(self, tmp_phase4):
        out = _run_phase4_pipeline(tmp_phase4)
        report = out / "integrated_phase4_report.md"
        assert report.exists()

    def test_report_sections(self, tmp_phase4):
        out = _run_phase4_pipeline(tmp_phase4)
        text = (out / "integrated_phase4_report.md").read_text(encoding="utf-8")
        assert "Executive Summary" in text
        assert "Why Affinity Alone" in text
        assert "Scoring Framework" in text
        assert "Ranked Candidates" in text
        assert "State Robustness" in text
        assert "Validation and Caveats" in text
        assert "File Inventory" in text

    def test_report_mentions_all_pockets(self, tmp_phase4):
        out = _run_phase4_pipeline(tmp_phase4)
        text = (out / "integrated_phase4_report.md").read_text(encoding="utf-8")
        assert "3GT8_raw_PKT01" in text
        assert "3GT8_raw_PKT02" in text
        assert "3GT8_raw_PKT03" in text


# ===================================================================
# TG 4.7: Presentation Summary
# ===================================================================

@pytest.mark.reporting
class TestPresentationSummary:
    """Tests for TG 4.7 presentation layer."""

    def test_shortlist_pocket_level(self, tmp_phase4):
        out = _run_phase4_pipeline(tmp_phase4)
        shortlist = _load_csv(out / "phase4_presentation_shortlist.csv")
        # 3 pockets, deduplicated from 5 review rows
        assert len(shortlist) == 3
        pids = [s["pocket_id"] for s in shortlist]
        assert len(pids) == len(set(pids)), "Shortlist should be pocket-level"

    def test_brief_markdown(self, tmp_phase4):
        out = _run_phase4_pipeline(tmp_phase4)
        brief = (out / "phase4_top_candidates_brief.md").read_text(encoding="utf-8")
        assert "Top Candidates Brief" in brief
        assert "Ligand Coverage Matrix" in brief

    def test_shortlist_has_rationale(self, tmp_phase4):
        out = _run_phase4_pipeline(tmp_phase4)
        shortlist = _load_csv(out / "phase4_presentation_shortlist.csv")
        for s in shortlist:
            assert s["rationale"], f"Missing rationale for {s['pocket_id']}"


# ===================================================================
# Schema Consistency
# ===================================================================

class TestSchemaConsistency:
    """Verify declared schemas match actual CSV output."""

    def test_evidence_normalized_schema(self, tmp_phase4):
        from egfr_pipeline.phase4.evidence_ingestion import NORMALIZED_COLUMNS
        from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
        e = tmp_phase4
        run_evidence_ingestion(e["phase1"], e["phase2"], e["phase3"], e["phase4"])
        actual = _csv_columns(e["phase4"] / "phase4_evidence_normalized.csv")
        assert actual == NORMALIZED_COLUMNS

    def test_axis_scores_schema(self, tmp_phase4):
        from egfr_pipeline.phase4.score_framework import SCORE_COLUMNS
        out = _run_phase4_pipeline(tmp_phase4)
        actual = _csv_columns(out / "phase4_axis_scores.csv")
        assert actual == SCORE_COLUMNS

    def test_classification_schema(self, tmp_phase4):
        from egfr_pipeline.phase4.mechanistic_classification import CLASSIFICATION_COLUMNS
        out = _run_phase4_pipeline(tmp_phase4)
        actual = _csv_columns(out / "final_candidate_classes.csv")
        assert actual == CLASSIFICATION_COLUMNS

    def test_perturbation_axis_scores_schema(self, tmp_phase4):
        from egfr_pipeline.phase4.perturbation_scoring import AXIS_SCORE_COLUMNS
        out = _run_phase4_pipeline(tmp_phase4)
        actual = _csv_columns(out / "perturbation_axis_scores.csv")
        assert actual == AXIS_SCORE_COLUMNS

    def test_candidate_table_schema(self, tmp_phase4):
        from egfr_pipeline.phase4.perturbation_scoring import CANDIDATE_TABLE_COLUMNS
        out = _run_phase4_pipeline(tmp_phase4)
        actual = _csv_columns(out / "perturbation_candidate_table.csv")
        assert actual == CANDIDATE_TABLE_COLUMNS

    def test_state_interpretation_schema(self, tmp_phase4):
        from egfr_pipeline.phase4.state_interpretation import INTERPRETATION_COLUMNS
        out = _run_phase4_pipeline(tmp_phase4)
        actual = _csv_columns(out / "phase4_state_interpretation.csv")
        assert actual == INTERPRETATION_COLUMNS

    def test_shortlist_schema(self, tmp_phase4):
        from egfr_pipeline.phase4.presentation_summary import SHORTLIST_COLUMNS
        out = _run_phase4_pipeline(tmp_phase4)
        actual = _csv_columns(out / "phase4_presentation_shortlist.csv")
        assert actual == SHORTLIST_COLUMNS


# ===================================================================
# End-to-End
# ===================================================================

class TestEndToEnd:
    """Full pipeline integration tests."""

    def test_full_pipeline(self, tmp_phase4):
        """Full pipeline: 3 pockets, 2 ligands, 5 ranked entries."""
        out = _run_phase4_pipeline(tmp_phase4)
        review = _load_csv(out / "phase4_final_review_table.csv")
        assert len(review) == 5
        # PKT01 should be rank 1
        assert review[0]["pocket_id"] == "3GT8_raw_PKT01"

    def test_no_phase1_pipeline(self, tmp_phase4_no_phase1):
        """Pipeline works without Phase 1 data."""
        out = _run_phase4_pipeline(tmp_phase4_no_phase1)
        review = _load_csv(out / "phase4_final_review_table.csv")
        assert len(review) >= 1
        # A1 should be 0 without Phase 1
        assert float(review[0]["A1_ppi_interface"]) == 0.0

    def test_multi_state_pipeline(self, tmp_phase4_multi_state):
        """Pipeline with robust pocket (3 states)."""
        out = _run_phase4_pipeline(tmp_phase4_multi_state)
        review = _load_csv(out / "phase4_final_review_table.csv")
        assert len(review) == 1
        r = review[0]
        assert r["mechanistic_class"] == "orthosteric_disruptor_candidate"
        assert float(r["A4_state_robustness"]) > 0.7

    def test_idempotent(self, tmp_phase4):
        """Running pipeline twice produces same output."""
        out1 = _run_phase4_pipeline(tmp_phase4)
        scores1 = (out1 / "perturbation_candidate_table.csv").read_text(encoding="utf-8")
        out2 = _run_phase4_pipeline(tmp_phase4)
        scores2 = (out2 / "perturbation_candidate_table.csv").read_text(encoding="utf-8")
        assert scores1 == scores2

    def test_orthosteric_beats_irrelevant(self, tmp_phase4):
        """Orthosteric disruptors always rank above irrelevant sites."""
        out = _run_phase4_pipeline(tmp_phase4)
        table = _load_csv(out / "perturbation_candidate_table.csv")
        ortho = [r for r in table if "orthosteric" in r["mechanistic_class"]]
        irrel = [r for r in table if "irrelevant" in r["mechanistic_class"]]
        if ortho and irrel:
            best_ortho = min(float(r["perturbation_score"]) for r in ortho)
            worst_irrel = max(float(r["perturbation_score"]) for r in irrel)
            assert best_ortho > worst_irrel
